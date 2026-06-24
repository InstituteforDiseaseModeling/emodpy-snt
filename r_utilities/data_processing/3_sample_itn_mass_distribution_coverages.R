# from DHS survey on ITN use, calculate starting use rates right after mass distributions (starting value for distributions - 'effective coverage')
#    note that we assume that use declines following the same rate/pattern as ITN retention as calculated in Amelia's work
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
}
library(plyr); library(dplyr)  # get error messages when load plyr after have already loaded dplyr
library(reshape2)
library(raster)
library(sp)
library(pals)
library(prettyGraphs) 
library(lubridate)
library(gridExtra) 
library(geofacet)
library(tidyr)

##############################################################
# sample net retention times and parameterize lognormal
##############################################################

# ============================================================
# function: get_lognormal_mu_from_A_halflife
# get the mu parameter for the lognormal net retention time (the function used in dtk simulations) that produces a distribution with the same mean or median value (depending on user selection) as the retention-time function used in Amelia's paper
#     note that matching the median retention time will have a longer lognormal retention than matching the mean retention time
#
# Matches a lognormal distribution to the smooth-compact Loss
# retention function using two approaches:
#
#   Method 1 (median match): set lognormal median = Loss half-life
#     mu = log(halflife_day)
#
#   Method 2 (mean match): set lognormal mean = Loss mean
#     mu = log(E[T_Loss]) - sigma^2 / 2
#     where E[T_Loss] = integral of Loss(t, kappa, tau) from 0 to tau
#       (computed numerically; requires kappa and tau)
#
# Args:
#   halflife_day  : half-life of the Loss function (days); the time at
#                   which Loss(t) = 0.5
#   sigma         : SD of the log (log-scale SD) of the lognormal
#   kappa         : shape parameter of the Loss function
#   match_type    : "median" to match the half-life (Method 1, default)
#                   "mean"   to match the mean retention time (Method 2)
#
# Returns:
#   A single numeric value: the lognormal mu for the chosen match_type.
#   A summary data.frame comparing both methods is printed as a side effect.
# ============================================================

get_lognormal_mu_from_A_halflife = function(A_halflife_day=(2.2*365), itn_lognorm_sigma=0.8, kappa=20, match_type="median", verbose=FALSE) {
  match_type <- match.arg(match_type, c("median","mean"))
  
  # --- Derive tau from halflife and kappa ------------------------------
  # From Loss(h, kappa, tau) = 0.5:  h = tau * sqrt(log(2) / (kappa + log(2)))
  # Rearranging:
  tau <- A_halflife_day * sqrt((kappa + log(2)) / log(2))
  
  # --- Method 1: match lognormal median to Loss half-life --------------
  mu_m1 <- log(A_halflife_day)

  # --- Method 2 / verbose: requires numerical integration --------------
  # Only compute when needed to avoid overhead during repeated sampling.
  if (match_type == "mean" || verbose) {
    loss_fn  <- function(t) exp(kappa - kappa / (1 - (t / tau)^2))
    E_T_loss <- integrate(
      loss_fn,
      lower        = 0,
      upper        = tau * (1 - 1e-8),
      subdivisions = 500L,
      rel.tol      = 1e-8
    )$value
    mu_m2 <- log(E_T_loss) - itn_lognorm_sigma^2 / 2
  }

  # --- Verbose comparison table ----------------------------------------
  if (verbose) {
    df <- data.frame(
      method           = c("M1_median_match", "M2_mean_match"),
      mu               = c(mu_m1, mu_m2),
      lognormal_median = c(exp(mu_m1),                            exp(mu_m2)),
      lognormal_mean   = c(exp(mu_m1 + itn_lognorm_sigma^2 / 2), exp(mu_m2 + itn_lognorm_sigma^2 / 2)),
      loss_halflife    = A_halflife_day,
      loss_mean        = E_T_loss,
      tau              = tau
    )
    print(df)
  }

  if (match_type == "median") mu_m1 else mu_m2
}




# OLD VERSION
# #    this value is calculated from the net halflife in the function from Amelia's paper
# #    process: 1) use median value from distribution A to get the mean value from distribution A
# #             2) find the mu lognormal parameter such that the mean of the lognormal distribution has the same mean as distribution A
# # retention-time function (1-CDF) from Amelia's paper:
# #    A(tt) = exp(18 - 18 / (1 - (tt/tau_days)^2)), 
# #    where tau_days = A_halflife_days * (1 - 18 / (18 - log(0.5)))^(-1/2)
# # retention-time probability distribution = d/dt[1-CDF] = A'(tt)
# #    A'(tt) = 18*exp(18 - 18 / (1 - (tt/tau_days)^2)) * 2*tt / tau_days^2 / (1 - (tt/tau_days)^2)^2)
# get_lognormal_mu_from_A_halflife = function(A_halflife_days=(1.31*365), itn_lognorm_sigma=0.8){
#   # calculate tau as used in Amelia's retention function
#   tau_days = A_halflife_days * (1 - 18 / (18 - log(0.5)))^(-1/2)  
#   
#   # find the mean (expected value) of the retention time distribution A(tt)... this will later be set as the mean of the dtk lognormal retention time distribution
#   # numeric approximation to get mean
#   xx = seq(0,(4.5*A_halflife_days),0.1)
#   A_mean = weighted.mean(xx,  18*exp(18 - 18 / (1 - (xx/tau_days)^2)) * 2*xx / tau_days^2 / (1 - (xx/tau_days)^2)^2)
#   
#   # the mean of the lognormal distribution is exp(mu + sigma^2 / 2), so mu = log(A_mean) - sigma^2 / 2
#   lognormal_mu = log(A_mean) - itn_lognorm_sigma^2 / 2
#   if(is.na(lognormal_mu)) warning('PROBLEM DETEcTED: sampled lognormal mu for net retention distribution is NA. This may be because the expected value was numerically calculated using too large a maximum retention time.')
#   return(lognormal_mu)
# }



# from the distribution of median LLIN retention times, sample median retention times and convert them into mus for the lognormal distribution used in the dtk
sample_net_longevity_params = function(hbhi_dir, num_samples=50, itn_lognorm_sigma=0.8,
                                       median_llin_retention_est = 1.31,  # years (from Amelia's retention time plot, SLE=1.47(1.31,1.63), BDI=1.31(1.14,1.47), BFA=1.58(1.41,1.76))
                                       median_llin_retention_lowerCI = 1.14,  # years (from Amelia's retention time plot, SLE=1.47(1.31,1.63), BDI=1.31(1.14,1.47), BFA=1.58(1.41,1.76))
                                       kill_decay_time=1460,
                                       block_decay_time=730){
  
  # draw S samples of the median net retention time (where S = num_samples values plus the mean value) and convert into the mu value for the lognormal distribution so that the means of both distributions match
  sampled_retentions_median_years = c(median_llin_retention_est, rnorm(n=(num_samples), mean=median_llin_retention_est, sd=((median_llin_retention_est-median_llin_retention_lowerCI)/2)))
  sampled_retentions_median_days = sampled_retentions_median_years * 365
  sampled_retentions_lognormal_mu = sapply(sampled_retentions_median_days, get_lognormal_mu_from_A_halflife, itn_lognorm_sigma=itn_lognorm_sigma)
  
  # create dataframe of net discard/decay parameters
  net_discard_decay = data.frame('seed'=seq(1,(num_samples+1)), 
                                 'net_life_median_years'=sampled_retentions_median_years, 
                                 'net_life_median_days'=sampled_retentions_median_days, 
                                 'net_life_lognormal_mu'=sampled_retentions_lognormal_mu, 
                                 'net_life_lognormal_sigma'=rep(itn_lognorm_sigma,(num_samples+1)), 
                                 'kill_decay_time'=rep(kill_decay_time,(num_samples+1)),
                                 'block_decay_time'=rep(block_decay_time,(num_samples+1)))
  # save csv of ITN discard and decay parameters
  write.csv(net_discard_decay, paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'), row.names=FALSE)
}


# from the distribution of median LLIN retention times, get quantiles for median retention times and convert them into mus for the lognormal distribution used in the dtk
quantile_net_longevity_params = function(hbhi_dir, quantiles=c(0.1, 0.3, 0.5, 0.7, 0.9), itn_lognorm_sigma=0.8,
                                       median_llin_retention_est = 1.31,  # years (from Amelia's retention time plot, SLE=1.47(1.31,1.63), BDI=1.31(1.14,1.47), BFA=1.58(1.41,1.76))
                                       median_llin_retention_lowerCI = 1.14,  # years (from Amelia's retention time plot, SLE=1.47(1.31,1.63), BDI=1.31(1.14,1.47), BFA=1.58(1.41,1.76))
                                       kill_decay_time=1460,
                                       block_decay_time=730){
  
  # get values at a series of quantiles for the median net retention time and convert into the mu value for the lognormal distribution so that the means of both distributions match
  sampled_retentions_median_years = qnorm(p=quantiles, mean=median_llin_retention_est, sd=((median_llin_retention_est-median_llin_retention_lowerCI)/2))
  sampled_retentions_median_days = sampled_retentions_median_years * 365
  sampled_retentions_lognormal_mu = sapply(sampled_retentions_median_days, get_lognormal_mu_from_A_halflife, itn_lognorm_sigma=itn_lognorm_sigma)
  
  # create dataframe of net discard/decay parameters
  net_discard_decay = data.frame('seed'=seq(1,length(quantiles)), 
                                 'net_life_median_years'=sampled_retentions_median_years, 
                                 'net_life_median_days'=sampled_retentions_median_days, 
                                 'net_life_lognormal_mu'=sampled_retentions_lognormal_mu, 
                                 'net_life_lognormal_sigma'=rep(itn_lognorm_sigma,length(quantiles)), 
                                 'kill_decay_time'=rep(kill_decay_time,length(quantiles)),
                                 'block_decay_time'=rep(block_decay_time,length(quantiles)))
  # save csv of ITN discard and decay parameters
  write.csv(net_discard_decay, paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'), row.names=FALSE)
}



#######################################################
# calculate access from nets distributed per capita
#######################################################
calc_access_from_npc = function(access, npc_access_param1=0.98, npc_access_param2=3){
  return(npc_access_param1*(1-(1-access)*npc_access_param2^(-1*access)))
}



##############################################################
# when cluster-level DHS data are not available but state-level ones are, add that information in the same format so that it can be used when estimating ITN coverage
##############################################################
add_state_DHS_ITN = function(additional_ITN_DHS_years, additional_ITN_DHS_files, ds_pop_df_filename){
  admin_info = read.csv(ds_pop_df_filename)
  net_dhs_info = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_ITN_dates_and_rates.csv'))
  net_dhs_info$date = as.Date(net_dhs_info$date, tryFormats=c('%m/%d/%Y', '%m-%d-%Y', '%Y-%m-%d'))  
  net_dhs_info_original = net_dhs_info
  cur_dhs_years = unique(year(net_dhs_info$date))
  file_revised = FALSE
  
  for(yy in 1:length(additional_ITN_DHS_years)){
    new_dhs_year = additional_ITN_DHS_years[yy]
    # check this year doesn't already exist in net_dhs_info
    if(!(new_dhs_year %in% cur_dhs_years)){
      file_revised=TRUE
      # read in state-level DHS ITN data for this year
      dhs_itn_state_cur = read.csv(additional_ITN_DHS_files[yy]) %>%
        rename(date = assumed_survey_date) %>%
        dplyr::select(State, itn_u5_rate, date) %>%
        mutate(date = as.Date(date, format='%m/%d/%Y'),
               State = str_trim(State, side = "right"))  # remove spaces from end of State names
      
      # standardize state names
      dhs_itn_state_cur = standardize_state_names_in_df(target_names_df=admin_info, origin_names_df=dhs_itn_state_cur, target_names_col='State', origin_names_col='State')
      
      # create values at the LGA level
      dhs_itn_admin_cur = merge(dhs_itn_state_cur, admin_info[,c('admin_name', 'State')], all=TRUE) %>%
        rename(NOMREGION = State) %>%
        dplyr::select(-matched_name)
      
      # merge into existing net_dhs_info
      net_dhs_info = merge(net_dhs_info, dhs_itn_admin_cur, all=TRUE)
    }
  }
  if(file_revised){
    write.csv(net_dhs_info, paste0(hbhi_dir, '/estimates_from_DHS/DHS_ITN_dates_and_rates.csv'), row.names=FALSE)
    write.csv(net_dhs_info_original, paste0(hbhi_dir, '/estimates_from_DHS/DHS_ITN_dates_and_rates_priorVersion.csv'), row.names=FALSE)
  }
  return(net_dhs_info)
}





##############################################################
# get state-level fraction of ITNs that came from mass distributions
##############################################################
get_itn_cov_attrib_campaigns = function(net_dhs_info, hbhi_dir, years=NA, overlap_dhs_years=NA){
  campaign_source_all = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/itn_source_state_level_all_years.csv')) %>%
    dplyr::select(State, year, frac_campaign)
  net_dhs_info$year = year(net_dhs_info$date)
  all_survey_years = unique(net_dhs_info$year)
  
  # some DHS/MIS surveys span over multiple years. in these cases, we want the fraction of ITNs from campaign to apply to both years
  if(any(!is.na(overlap_dhs_years))){
    # create a lookup table to see the additional years for each survey
    additional_year_lookup = data.frame(years, overlap_dhs_years)
    lookup_non_na = additional_year_lookup %>% filter(!is.na(overlap_dhs_years))
    
    campaign_source_all = lookup_non_na %>%
      # keep only rows with non-missing years and overlap_dhs_years
      filter(!is.na(years), !is.na(overlap_dhs_years)) %>%
      group_map(~ {
        orig_year <- .x$years[1]
        new_year  <- .x$overlap_dhs_years[1]
        
        campaign_source_all %>%
          filter(year == orig_year) %>%
          mutate(year = new_year)
      }) %>%
      bind_rows() %>% # combine all the new rows
      bind_rows(campaign_source_all) %>%  # add to the original df
      distinct() 
  }
  
  # some survey years do not have estimates of the fraction of ITNs from campaigns. in those years, we use the average for each state when data is available
  campaign_complete = campaign_source_all %>%
    group_by(State) %>%
    # make sure every state has all years
    complete(year = all_survey_years) %>%
    # fill estimate with state mean if missing
    mutate(frac_campaign = ifelse(is.na(frac_campaign), mean(frac_campaign, na.rm = TRUE), frac_campaign)) %>%
    ungroup() %>%
    rename(NOMREGION = State)

  # merge into the DHS coverage estimates df and adjust the coverage rates to reflect coverage from campaign-distributed ITNs
  net_campaign_dhs_info = merge(net_dhs_info, campaign_complete, all.x=TRUE) %>%
    mutate(across(contains("_rate"), ~ .x * frac_campaign))
  return(net_campaign_dhs_info)
}


############################################################################################################################
############################################################################################################################
# back-calculated mass distribution coverage from observed net use at time of DHS survey
############################################################################################################################
############################################################################################################################


# - - - - - - version of functions with single seed but multiple different distribution dates depending on admin - - - - - - - #

# prepare DHS ITN dataframe
aggregate_itn_dhs_data_across_years = function(hbhi_dir, years, itn_variables, min_num_total=30, overwrite=FALSE){
  net_dhs_filename = paste0(hbhi_dir, '/estimates_from_DHS/DHS_ITN_dates_and_rates.csv')
  if(file.exists(net_dhs_filename) & !(overwrite)){
    net_dhs_info = read.csv(net_dhs_filename)
  } else{
    net_dhs_info = data.frame()
    for(yy in 1:length(years)){
      admin_sums = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_minN', min_num_total,'_', years[yy], '.csv'))[,-1]
      net_dhs_info_yy = admin_sums[, c('NOMREGION', 'NOMDEP', 'Archetype', 'mean_date', paste0(itn_variables, '_rate'))]
      if(nrow(net_dhs_info)>0){
        net_dhs_info = merge(net_dhs_info, net_dhs_info_yy, all=TRUE)
      } else {
        net_dhs_info = net_dhs_info_yy
      }
    }
    colnames(net_dhs_info)[which(colnames(net_dhs_info)=='NOMDEP')] = 'admin_name'
    colnames(net_dhs_info)[which(colnames(net_dhs_info)=='mean_date')] = 'date'
    write.csv(net_dhs_info, net_dhs_filename, row.names=FALSE)
  }
  return(net_dhs_info)
}


##############################################################
# create input file for mass distribution of LLINs
##############################################################
# from DHS output on cm rates in U5, create a csv file that will be read in by simulations. each row corresponds to an admin and specifies some period of implementation in that admin. included in the csv will be
#    - admin name
#    - seed number (determines which sampled parameter set is used)
#    - start year
#    - simday when distribution occurs
#    - U5 coverage
#    - 5-10 coverage
#    - 10-15 coverage
#    - 15-20 coverage
#    - adult coverage
#    - type of net
#    - blocking_rate
#    - mortality_rate
#    - kill_rate
#    - net_life_lognormal_mu  # for the Expiration_Period_Log_Normal_Mu parameter in the lognormal decay distribution (time before nets discarded, lost, forgotten, etc.).
#    - net_life_lognormal_sigma
# single seed and admins may have different mass distributions dates
create_itn_input_from_DHS_differentDates = function(hbhi_dir, itn_variables, itn_distributions_by_admin_filename, grid_layout_state_locations, ds_pop_df_filename, sim_start_year=2010, maximum_coverage=0.9,
                                                    seasonality_monthly_scalar,  # adjust net usage for seasonality
                                                    years, min_num_total=30, default_first_coverage=0.1, itn_variable_base='itn_u5', save_age_ratio_plots=FALSE, save_timeseries_coverage_plots=FALSE,
                                                    additional_ITN_DHS_years=c(), additional_ITN_DHS_files=c(), overlap_dhs_years=NA
){
  # get the distribution dates for each admin
  itn_distributions_by_admin = read.csv(itn_distributions_by_admin_filename) 
  itn_distributions_by_admin$date = as.Date(itn_distributions_by_admin$date)
  # get the DHS value and date in each admin (includes all ITN variables)
  net_dhs_info = aggregate_itn_dhs_data_across_years(hbhi_dir, years, itn_variables, min_num_total)
  
  # add any additional DHS years that were recorded at the state level rather than cluster level
  if((length(additional_ITN_DHS_years)>0) & is.numeric(additional_ITN_DHS_years[1])){
    net_dhs_info = add_state_DHS_ITN(additional_ITN_DHS_years=additional_ITN_DHS_years, additional_ITN_DHS_files=additional_ITN_DHS_files, ds_pop_df_filename=ds_pop_df_filename)
  }
  net_dhs_info$date = as.Date(net_dhs_info$date, tryFormats=c('%m/%d/%Y', '%m-%d-%Y', '%Y-%m-%d'))  
  
  # adjust DHS/MIS survey coverage to capture the ITNs from mass campaigns
  net_dhs_info = get_itn_cov_attrib_campaigns(net_dhs_info=net_dhs_info, hbhi_dir=hbhi_dir, years=years, overlap_dhs_years=overlap_dhs_years)

  # read in sampled retention lognormal mus and sigmas and take the first (expected) value
  net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
  net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[1]
  net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[1]
  
  
  # iterate through admins where distributions occurred, calculating the distribution coverages from the coverage at the next DHS survey year and the retention times
  all_admins = unique(itn_distributions_by_admin$admin_name)
  all_admins = all_admins[!is.na(all_admins)]
  coverage_df = data.frame()
  for(aa in 1:length(all_admins)){
    # subset dhs coverage data and itn distribution dates to this admin
    cur_net_dhs_info = net_dhs_info[net_dhs_info$admin_name == all_admins[aa],]
    cur_distributions = itn_distributions_by_admin[itn_distributions_by_admin$admin_name == all_admins[aa],]
    cur_distributions = cur_distributions[order(cur_distributions$date),]

    # assign which dhs survey should be used to ascertain coverage for each distribution
    cur_distributions$matched_dhs_date = as.Date(NA)
    cur_distributions$matched_dhs_val = NA
    for(i_dist in 1:nrow(cur_distributions)){
      if(any(cur_net_dhs_info$date>cur_distributions$date[i_dist])){
        cur_distributions$matched_dhs_date[i_dist] = cur_net_dhs_info$date[min(which(cur_net_dhs_info$date>cur_distributions$date[i_dist]))]
        cur_distributions$matched_dhs_val[i_dist] = cur_net_dhs_info[[paste0(itn_variable_base,'_rate')]][min(which(cur_net_dhs_info$date>cur_distributions$date[i_dist]))]
        if(i_dist<nrow(cur_distributions)){
          if(cur_distributions$date[i_dist+1] < cur_distributions$matched_dhs_date[i_dist]){ # if next distribution is before dhs survey, we don't know coverage from current distribution
            cur_distributions$matched_dhs_date[i_dist] = NA
            cur_distributions$matched_dhs_val[i_dist] = NA
          }
        }
      }else{
        cur_distributions$matched_dhs_date[i_dist] = NA
        cur_distributions$matched_dhs_val[i_dist] = NA
      }
    }
    # assign distribution use-coverages
    # if a distribution is associated with a dhs survey,
    #   1) adjust for the seasonal use at the time of the survey
    #   2) subtract the residual coverage from prior net distributions (assume that net distributions before the previous one have negligible contribution to net total)
    #   3) back-calculate what coverage would have been from distribution to have given rise to observed coverage
    # if the first distribution is not associated with a dhs survey, use a default coverage (default_first_coverage).
    # if a later distribution is not associated with a dhs survey, use the coverage from the prior survey.
    cur_distributions$coverage = NA
    for(i_dist in 1:nrow(cur_distributions)){
      if(is.na(cur_distributions$matched_dhs_date[i_dist])){
        if(i_dist==1){
          cur_distributions$coverage[i_dist] = default_first_coverage
        } else{
          cur_distributions$coverage[i_dist] = cur_distributions$coverage[i_dist-1]
        }
      } else{
        if(i_dist>1){
          # time between the current dhs survey and the prior distribution (i_dist-1)
          time_lapse_days_residual = as.integer(cur_distributions$matched_dhs_date[i_dist] - cur_distributions$date[i_dist-1])
          # calculate what coverage remains from the previous distribution
          residual_coverage = cur_distributions$coverage[i_dist-1] * (1-plnorm(time_lapse_days_residual, meanlog=net_life_lognormal_mu, sdlog=net_life_lognormal_sigma))
        } else{
          residual_coverage = 0
        }
        season_adjusted_dhs_coverage = min(maximum_coverage, cur_distributions$matched_dhs_val[i_dist]  / seasonality_monthly_scalar[month(cur_distributions$matched_dhs_date[i_dist])])
        # remove nets leftover from prior distribution
        # season_adjusted_dhs_coverage = max(0, season_adjusted_dhs_coverage-residual_coverage)
        # UPDATE: nets are distributed at random, not selectively to people who didn't have them before, so we don't just subtract old nets
        #     total_coverage = 1- ((1-coverage_from_residual_nets)*(1-coverage_from_new_nets))
        season_adjusted_dhs_coverage = max(0, 1-(1-season_adjusted_dhs_coverage)/(1-residual_coverage))
        
        # extrapolate back to what the distribution coverage would need to be to obtain this observed coverage
        # time between the dhs survey and the distribution
        time_lapse_days = as.integer(cur_distributions$matched_dhs_date[i_dist] - cur_distributions$date[i_dist])
        # fraction of orginally-disributed nets remaining at the time of the survey
        remaining_net_frac = 1-plnorm(time_lapse_days, meanlog=net_life_lognormal_mu, sdlog=net_life_lognormal_sigma)
        backward_projected_dist_coverage = season_adjusted_dhs_coverage / remaining_net_frac
        backward_projected_dist_coverage = min(maximum_coverage, backward_projected_dist_coverage)
        cur_distributions$coverage[i_dist] = backward_projected_dist_coverage   
      }
    }
    # add current admin into dataframe containing distribution coverage for all admins
    if(nrow(coverage_df)<1){
      coverage_df = cur_distributions
    } else{
      coverage_df = rbind(coverage_df, cur_distributions)
    }
  }
  # add in age-specific net-use rates at time of distribution, taking the mean across all dhs-years and LGAs
  itn_variable_age_scales = rep(NA,length(itn_variables))
  net_dhs_info$year = year(net_dhs_info$date)
  gg_list = list()
  gg_list_index = 1
  for(i_var in 1:length(itn_variables)){
    if(itn_variable_base != itn_variables[i_var]){
      net_dhs_info$ratio = net_dhs_info[[paste0(itn_variables[i_var],'_rate')]] / net_dhs_info[[paste0(itn_variable_base,'_rate')]]
      itn_variable_age_scales[i_var] = median(net_dhs_info$ratio, na.rm=TRUE)
      if(save_age_ratio_plots){
        gg = ggplot(data=net_dhs_info, aes(x=ratio, color=as.factor(year), fill=as.factor(year))) + 
          geom_histogram(aes(y=..density..), alpha=.2, position="identity")+
          # geom_density(alpha=.2) +
          geom_vline(xintercept=itn_variable_age_scales[i_var], size=1.2) +
          coord_cartesian(xlim=c(0,8)) +
          ggtitle(paste0(itn_variables[i_var],' / ', itn_variable_base))+
          theme_bw() + 
          theme(legend.position='none')
      gg_list[[gg_list_index]] = gg
      gg_list_index = gg_list_index + 1
      }
    } else{
      itn_variable_age_scales[i_var] = 1
    }
  }
  if(save_age_ratio_plots){
    # do.call("grid.arrange", c(gg_list, nrow=1)) 
    ggsave(filename=paste0(hbhi_dir, '/estimates_from_DHS/plots/itn_use_rate_age_ratios.png'), plot=arrangeGrob(gg_list[[1]], gg_list[[2]], gg_list[[3]], gg_list[[4]],  nrow=1), width=10, height=3, units='in')
  }
  # add columns for all of the itn age groups to coverage dataframe
  for(i_var in 1:length(itn_variables)){
    coverage_df[[itn_variables[i_var]]] = coverage_df$coverage * itn_variable_age_scales[i_var]
  }

  # add in the day of the simulation each intervention should start and the duration
  coverage_df$simday = coverage_df$date - as.Date(paste0(sim_start_year,'-01-01'))
  # add in net decay info
  coverage_df$net_life_lognormal_mu =  net_life_lognormal_mu
  coverage_df$net_life_lognormal_sigma = net_life_lognormal_sigma 
  
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
  write.csv(coverage_df, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_toPresent.csv'), row.names=FALSE)
  
  if(save_timeseries_coverage_plots){
    
    # plot timeseries of ITN coverage from mass distributions in each LGA given distribution schedule, seasonality in ITN use, and net decay parameters
    # also include dots for DHS observations
    # data format: a long dataframe with columns for LGA name, state name, archetype name, date, and coverage
    net_dhs_info$State = net_dhs_info$NOMREGION
    coverage_df = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_toPresent.csv'))
    coverage_df$date = as.Date(coverage_df$date)
    coverage_df$matched_dhs_date = as.Date(coverage_df$matched_dhs_date)
    coverage_df = coverage_df[!is.na(coverage_df$date), ]
    if(nrow(coverage_df) == 0){
      stop("No valid dates found in coverage_df$date after parsing.")
    }
    all_admins = unique(coverage_df$admin_name)
    coverage_timeseries = data.frame('admin_name'=c(), 'State'=c(), 'date'=c(), 'coverage'=c())
    first_day = as.Date('2010-01-01')
    final_day = max(coverage_df$date, na.rm=TRUE)
    # date_vector = seq(first_day, as.Date('2022-01-01'), by=1)
    date_vector = seq.Date(first_day, final_day, by='month')
    timeseries_length = length(date_vector)
    for(aa in 1:length(all_admins)){
      cur_distributions = coverage_df[coverage_df$admin_name == all_admins[aa],]
      coverage_vector = rep(0, timeseries_length)
      for(i_dist in 1:nrow(cur_distributions)){
        dist_index = which.min(abs(date_vector - cur_distributions$date[i_dist]))
        indices_updated = seq(from=dist_index, by=1, length.out=min((12*5), (timeseries_length-dist_index+1)))
        coverage_vector[indices_updated] = 1-((1-coverage_vector[indices_updated]) * (1-(cur_distributions$itn_u5[i_dist] * (1-plnorm((round(seq(from=1, by=30.4, length.out=(12*5)))), meanlog=net_life_lognormal_mu, sdlog=net_life_lognormal_sigma)))[1:length(indices_updated)]))
      }
      admin_df = data.frame('admin_name'=rep(all_admins[aa], timeseries_length), 
                            'State'=rep(coverage_df$State[coverage_df$admin_name == all_admins[aa]][1], timeseries_length),
                            'date'=date_vector, 
                            'coverage'=coverage_vector)
      coverage_timeseries = rbind(coverage_timeseries, admin_df)
    }
    # set maximum coverage at 1
    coverage_timeseries$coverage = sapply(coverage_timeseries$coverage, min, 1)
    # adjust for seasonal net usage
    coverage_timeseries$month = month(coverage_timeseries$date)
    coverage_timeseries$seasonal_adjust = seasonality_monthly_scalar[coverage_timeseries$month]
    coverage_timeseries$adjusted_coverage = coverage_timeseries$coverage * coverage_timeseries$seasonal_adjust
    
    gg=ggplot()+
      geom_line(data=coverage_timeseries, aes(x=date, y=adjusted_coverage, col=admin_name)) +
      geom_point(data=net_dhs_info, aes(x=date,y=itn_u5_rate, col=admin_name)) +
      geom_point(data=net_dhs_info, aes(x=date,y=itn_u5_rate), shape=21, col='black') +
      geom_point(data=coverage_df, aes(x=date), y=1, col='black', shape='|', size=1) +
      geom_point(data=coverage_df, aes(x=date), y=0.98, col='black', shape='V', size=1) +
      coord_cartesian(xlim=c(first_day, final_day)) +
      theme_bw()+
      theme(legend.position='none') +
      facet_wrap('State', nrow=5)
    ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_rate_timeseries_extrapolation_wSeasonAdjust.png'), plot=gg, width=18, height=15, units='in', dpi=900)
    
    gg=ggplot()+
      geom_line(data=coverage_timeseries, aes(x=date, y=coverage, col=admin_name)) +
      geom_point(data=net_dhs_info, aes(x=date,y=itn_u5_rate, col=admin_name)) +
      geom_point(data=net_dhs_info, aes(x=date,y=itn_u5_rate), shape=21, col='black') +
      coord_cartesian(xlim=c(first_day, final_day)) +
      theme_bw()+
      theme(legend.position='none') +
      facet_wrap('State', nrow=5)
    # coverage_timeseries_daily = coverage_timeseries
    # gg_daily=gg
    ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_rate_timeseries_extrapolation_noSeasonAdjust.png'), plot=gg, width=18, height=15, units='in', dpi=900)
    
    # for build slide
    ggb=ggplot()+
      geom_point(data=net_dhs_info, aes(x=date,y=itn_u5_rate, col=admin_name)) +
      geom_point(data=net_dhs_info, aes(x=date,y=itn_u5_rate), shape=21, col='black') +
      geom_point(data=coverage_df, aes(x=date), y=1, col='black', shape='|', size=1) +
      geom_point(data=coverage_df, aes(x=date), y=0.98, col='black', shape='V', size=1) +
      coord_cartesian(xlim=c(first_day, final_day)) +
      theme_bw()+
      theme(legend.position='none') +
      facet_wrap('State', nrow=5)
    ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_rate_timeseries_extrapolation_build0.png'), plot=ggb, width=18, height=15, units='in', dpi=900)
    
    
    itn_distributions_by_state_raw_filename = gsub('/[^/]*.csv','/uncorrected_mass_dist_dates_state.csv',itn_distributions_by_admin_filename)
    if(file.exists(itn_distributions_by_state_raw_filename)){
      llin_info_raw_state = read.csv(itn_distributions_by_state_raw_filename)
      net_dhs_info$code = net_dhs_info$NOMREGION
      llin_info_raw_state$code = llin_info_raw_state$State
      llin_info_raw_state$date = as.Date(llin_info_raw_state$date)

      ggc1 = ggplot(grid_layout_state_locations)+
        # geom_vline(data=llin_info_raw_state, aes(xintercept=date), color='darkgreen')+
        geom_hline(aes(yintercept=0), color='black')+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate), color='black', size=0.8)+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate, color=admin_name), size=0.6)+
        scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
        scale_x_date(guide = guide_axis(check.overlap = TRUE), breaks=as.Date(paste0(c(2012,2016,2020),'/01/01')), labels=c(2012,2016,2020)) +
        coord_cartesian(xlim=c(first_day, final_day)) +
        theme_classic()+
        theme(legend.position='none') +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') 
      ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_from_dhs.png'), plot=ggc1, width=9, height=6, units='in', dpi=1600)
      
      ggc = ggplot(grid_layout_state_locations)+
        geom_vline(data=llin_info_raw_state, aes(xintercept=date), color='darkgreen')+
        geom_hline(aes(yintercept=0), color='black')+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate), color='black', size=0.8)+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate, color=admin_name), size=0.6)+
        scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
        scale_x_date(guide = guide_axis(check.overlap = TRUE), breaks=as.Date(paste0(c(2012,2016,2020),'/01/01')), labels=c(2012,2016,2020)) +
        coord_cartesian(xlim=c(first_day, final_day)) +
        theme_classic()+
        theme(legend.position='none') +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') 
      ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_from_dhs_and_raw_dist_dates.png'), plot=ggc, width=9, height=6, units='in', dpi=1600)
      
      coverage_df$code = coverage_df$State
      ggd = ggplot(grid_layout_state_locations)+
        geom_vline(data=llin_info_raw_state, aes(xintercept=date), color='darkgreen')+
        geom_point(data=coverage_df, aes(x=date), y=1, col='blue', shape='|', size=7) +
        geom_hline(aes(yintercept=0), color='black')+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate), color='black', size=0.8)+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate, color=admin_name), size=0.6)+
        scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
        scale_x_date(guide = guide_axis(check.overlap = TRUE), breaks=as.Date(paste0(c(2012,2016,2020),'/01/01')), labels=c(2012,2016,2020)) +
        coord_cartesian(xlim=c(first_day, final_day)) +
        theme_classic()+
        theme(legend.position='none') +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') 
      ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_from_dhs_and_estimated_dist_dates.png'), plot=ggd, width=9, height=6, units='in', dpi=1600)
      
      ggdv2 = ggplot(grid_layout_state_locations)+
        # geom_vline(data=llin_info_raw_state, aes(xintercept=date), color='darkgreen')+
        geom_vline(data=coverage_df, aes(xintercept=date), color='blue',  size=0.3) +
        geom_hline(aes(yintercept=0), color='black')+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate), color='black', size=0.8)+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate, color=admin_name), size=0.6)+
        scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
        scale_x_date(guide = guide_axis(check.overlap = TRUE), breaks=as.Date(paste0(c(2012,2016,2020),'/01/01')), labels=c(2012,2016,2020)) +
        coord_cartesian(xlim=c(first_day, final_day)) +
        theme_classic()+
        theme(legend.position='none') +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') 
      ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_from_dhs_and_estimated_dist_dates_v2.png'), plot=ggdv2, width=9, height=6, units='in', dpi=1600)
      
      
      coverage_timeseries$code = coverage_timeseries$State
      gge = ggplot(grid_layout_state_locations)+
        geom_line(data=coverage_timeseries, aes(x=date, y=adjusted_coverage, col=admin_name), linewidth=0.15) +
        geom_vline(data=llin_info_raw_state, aes(xintercept=date), color='darkgreen')+
        geom_point(data=coverage_df, aes(x=date), y=1, col='blue', shape='|', size=7) +
        geom_hline(aes(yintercept=0), color='black')+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate), color='black', size=0.8)+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate, color=admin_name), size=0.6)+
        scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
        scale_x_date(guide = guide_axis(check.overlap = TRUE), breaks=as.Date(paste0(c(2012,2016,2020),'/01/01')), labels=c(2012,2016,2020)) +
        coord_cartesian(xlim=c(first_day, final_day)) +
        theme_classic()+
        theme(legend.position='none') +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') 
      ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_timeseries_from_dhs_and_estimated_dist_dates.png'), plot=gge, width=9, height=6, units='in', dpi=1600)
      
      
      ggev2 = ggplot(grid_layout_state_locations)+
        geom_line(data=coverage_timeseries, aes(x=date, y=adjusted_coverage, col=admin_name), linewidth=0.15) +
        # geom_vline(data=llin_info_raw_state, aes(xintercept=date), color='darkgreen')+
        geom_vline(data=coverage_df, aes(xintercept=date), color='blue',  size=0.3) +
        geom_hline(aes(yintercept=0), color='black')+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate), color='black', size=0.8)+
        geom_point(data=net_dhs_info, aes(x=date, y=itn_u5_rate, color=admin_name), size=0.6)+
        scale_y_continuous(guide = guide_axis(check.overlap = TRUE)) +
        scale_x_date(guide = guide_axis(check.overlap = TRUE), breaks=as.Date(paste0(c(2012,2016,2020),'/01/01')), labels=c(2012,2016,2020)) +
        coord_cartesian(xlim=c(first_day, final_day)) +
        theme_classic()+
        theme(legend.position='none') +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') 
      ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/itn_use_timeseries_from_dhs_and_estimated_dist_dates_v2.png'), plot=ggev2, width=9, height=6, units='in', dpi=1600)
      
      
    }

  }
}



# - - - - - - version of functions with multiple seeds but same distribution dates for all admin - - - - - - - #


##############################################################
# helper function that takes as input coverages at time of DHS survey for a particular age group and calculates the coverage at the time of net distribution
##############################################################
# when the entire country gets nets distributed at the same time, can use function below:
get_dist_coverage = function(coverage_df,  # coverage at the time of DHS survey; also contains parameters describing net retention
                             maximum_coverage,
                             itn_distribution_years, 
                             itn_distribution_months,
                             matching_dhs_years,  # DHS survey year following mass ITN distributions (must match order of itn_distribution_years)
                             matching_dhs_months,  # chose the median month across all dhs entries, not specific for each admin
                             dhs_rescale_multipliers,
                             itn_variables,
                             i_var){
  # iterate through distribution years, calculating the distribution coverages from the coverage at the next DHS survey year and the retention times
  for(yy in 1:length(itn_distribution_years)){
    
    # get the coverage values for relevant DHS survey year
    coverage_df_cur = coverage_df[coverage_df$year == matching_dhs_years[yy],]
    coverage_df_cur$dist_year = itn_distribution_years[yy]
    
    # adjust to account for lower seasonal ITN use
    coverage_df_cur$value = sapply((coverage_df_cur$value * dhs_rescale_multipliers[yy]), min, maximum_coverage)
    # if this isn't the first distribution, check how much of this coverage is attributable to earlier distributions and subtract
    if(yy>1){
      for(y_prev in 1:(yy-1)){
        
        # calculate the time betweeen the distribution and DHS survey
        time_lapse_days = as.integer(as.Date(paste0(matching_dhs_years[yy], '/',matching_dhs_months[yy],'/1')) - as.Date(paste0(itn_distribution_years[y_prev], '/',itn_distribution_months[y_prev],'/1')))
        
        # calculate what coverage remains from the previous distribution
        leftover_nets = dist_coverage_var[dist_coverage_var$dist_year == itn_distribution_years[y_prev],]
        leftover_nets$from_previous_dist = leftover_nets[[itn_variables[i_var]]] * (1-plnorm(time_lapse_days, meanlog=leftover_nets$net_life_lognormal_mu, sdlog=leftover_nets$net_life_lognormal_sigma))
        leftover_nets = leftover_nets[,colnames(leftover_nets) %in% c('seed','admin_name','net_life_lognormal_mu', 'net_life_lognormal_sigma','from_previous_dist')]
        
        if('seed' %in% colnames(leftover_nets)){
          coverage_df_cur = merge(coverage_df_cur, leftover_nets, by=c('seed','admin_name','net_life_lognormal_mu', 'net_life_lognormal_sigma'))
        } else{
          coverage_df_cur = merge(coverage_df_cur, leftover_nets, by=c('admin_name','net_life_lognormal_mu', 'net_life_lognormal_sigma'))
        }
        # subtract the nets that remain from the previous distribution from the coverage observed in the survey, since those weren't distributed during the current campaign
        coverage_df_cur$value = coverage_df_cur$value - coverage_df_cur$from_previous_dist
        coverage_df_cur$value = sapply(coverage_df_cur$value, max, 0)  # make sure coverage isn't less than zero
        coverage_df_cur = coverage_df_cur[,colnames(coverage_df_cur) %in% c('seed','admin_name','net_life_lognormal_mu', 'net_life_lognormal_sigma','year','value','dist_year')]  # get rid of previous distribution column
      }
    }
    
    # calculate the time betweeen the distribution and DHS survey
    time_lapse_days = as.integer(as.Date(paste0(matching_dhs_years[yy], '/',matching_dhs_months[yy],'/1')) - as.Date(paste0(itn_distribution_years[yy], '/',itn_distribution_months[yy],'/1')))
    
    # extrapolate back to what the distribution coverage would need to be to obtain this observed coverage
    # fraction of orginally-disributed nets remaining at the time of the survey
    remaining_nets = 1-plnorm(time_lapse_days, meanlog=coverage_df_cur$net_life_lognormal_mu, sdlog=coverage_df_cur$net_life_lognormal_sigma)
    # original coverage
    coverage_df_cur[[itn_variables[i_var]]] = coverage_df_cur$value / remaining_nets
    coverage_df_cur[[itn_variables[i_var]]] = sapply(coverage_df_cur[[itn_variables[i_var]]], min, maximum_coverage)
    coverage_df_cur$month = itn_distribution_months[yy]
    
    # save this year's distribution coverage to data frame
    colnames(coverage_df_cur)[which(colnames(coverage_df_cur)=='value')] = paste0('dhs_',itn_variables[i_var])
    colnames(coverage_df_cur)[which(colnames(coverage_df_cur)=='year')] = 'dhs_year'
    if(yy == 1){  # create data frame for this variable
      dist_coverage_var = coverage_df_cur#[,which(colnames(coverage_df_cur) %in% c('seed','admin_name','dist_year',itn_variables[i_var],'net_life_lognormal_mu'))]
    } else{ # add to this variable's data frame
      dist_coverage_var = rbind(dist_coverage_var, coverage_df_cur)#[,which(colnames(coverage_df_cur) %in% c('seed','admin_name','dist_year',itn_variables[i_var],'net_life_lognormal_mu'))])
    }
  }
  return(dist_coverage_var)
}


##############################################################
# create input file for mass distribution of LLINs
##############################################################
# from DHS output on cm rates in U5, create a csv file that will be read in by simulations. each row corresponds to an admin and specifies some period of implementation in that admin. included in the csv will be
#    - admin name
#    - seed number (determines which sampled parameter set is used)
#    - start year
#    - simday when distribution occurs
#    - U5 coverage
#    - 5-10 coverage
#    - 10-15 coverage
#    - 15-20 coverage
#    - adult coverage
#    - type of net
#    - blocking_rate
#    - mortality_rate
#    - kill_rate
#    - net_life_lognormal_mu  # for the Expiration_Period_Log_Normal_Mu paraeter in the lognormal decay distribution (time before nets discarded, lost, forgotten, etc.).
#    - net_life_lognormal_sigma
# multiple seeds, but assumes that all admins have mass distributions at the same time
create_itn_input_from_DHS = function(hbhi_dir, itn_variables, num_samples=50, sim_start_year=2010, maximum_coverage=0.9,
                                     itn_distribution_years = c(2010, 2011, 2014),  # since the 2010-2020 simulation starts at the first day of 2010, consider the 2009 distribution to occur on the first day of 2010 to get those nets in the simulation
                                     itn_distribution_months = c(1, 7 ,7),
                                     matching_dhs_years = c(2010,2012,2016),  # DHS survey year following mass ITN distributions (must match order of itn_distribution_years)
                                     matching_dhs_months = c(10, 12, 10),  # chose the median month across all dhs entries, not specific for each admin
                                     dhs_rescale_multipliers  # adjust net usage for seasonality
){
  
  # read in sampled retention lognormal mus and sigmas
  net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
  net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
  
  # iterate through variables
  for(i_var in 1:length(itn_variables)){
    
    # get the sampled value from the DHS survey in each admin (num_samples values plus the mean value)
    sample_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_', itn_variables[i_var],'.csv'))[,-1]
    # convert from wide to long format to make each seed its own row
    coverage_df = reshape2::melt(sample_df, id.vars=c("admin_name", "year"))
    colnames(coverage_df)[colnames(coverage_df)=='variable'] = 'seed'
    # colnames(coverage_df)[colnames(coverage_df)=='value'] = itn_variables[i_var]
    # replace "sample_" with "" in the seed column
    coverage_df$seed = gsub('sample_','',coverage_df$seed)
    # merge in the net halflife for each seed
    coverage_df = merge(coverage_df, net_discard_decay, by='seed')
    
    # iterate through distribution years, calculating the distribution coverages from the coverage at the next DHS survey year and the retention times
    dist_coverage_var = get_dist_coverage(coverage_df, maximum_coverage, itn_distribution_years, itn_distribution_months, matching_dhs_years, matching_dhs_months, 
                                          dhs_rescale_multipliers, itn_variables, i_var)
    
    # save this coverage variable to data frame
    if(i_var == 1){
      dist_coverage = dist_coverage_var
    } else{
      dist_coverage = merge(dist_coverage, dist_coverage_var, by=c('seed','admin_name','dist_year','net_life_lognormal_mu','net_life_lognormal_sigma', 'dhs_year', 'month'))
    }
  }
  
  # add in the day of the simulation each intervention should start and the duration
  dist_coverage$simday = (dist_coverage$dist_year - sim_start_year) * 365 + round(dist_coverage$month*30)
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
  write.csv(dist_coverage, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_2020.csv'), row.names=FALSE)
}














# create two intervention files: one for the burnin and one for the calibration simulation. Use only mean rates for all-age ITN use and retention function.
create_seasonality_calibration_itn_from_DHS = function(hbhi_dir, itn_variables, burn_start_year=1961, calib_start_year=2011, maximum_coverage=0.9,
                                           itn_distribution_years_calib = c(2009, 2011, 2014),  # since the 2010-2020 simulation starts at the first day of 2010, consider the 2009 distribution to occur on the first day of 2010 to get those nets in the simulation
                                           itn_distribution_months_calib = c(7, 7 ,7),
                                           matching_dhs_years = c(2010,2012,2016),  # DHS survey year following mass ITN distributions (must match order of itn_distribution_years)
                                           matching_dhs_months = c(10, 12, 10),  # chose the median month across all dhs entries, not specific for each admin
                                           dhs_rescale_multipliers){
  
  archetype_rates = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_archetype_rates.csv'))
 
  # read in sampled retention lognormal mus and sigmas and use values for seed=1 (the expected values)
  net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
  calib_net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[net_discard_decay$seed ==1]
  calib_net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[net_discard_decay$seed ==1]
  
  # iterate through variables (age groups)
  for(i_var in 1:length(itn_variables)){
    coverage_df = data.frame('value' = archetype_rates[[paste0(itn_variables[i_var],'_rate')]], 'admin_name'=archetype_rates$seasonality_archetype, 'year'=archetype_rates$year,'net_life_lognormal_mu'=rep(calib_net_life_lognormal_mu, nrow(archetype_rates)),'net_life_lognormal_sigma'=rep(calib_net_life_lognormal_sigma, nrow(archetype_rates)))
    
    
    # iterate through distribution years, calculating the distribution coverages from the coverage at the next DHS survey year and the retention times
    dist_coverage_var = get_dist_coverage(coverage_df, maximum_coverage, itn_distribution_years=itn_distribution_years_calib, itn_distribution_months=itn_distribution_months_calib, matching_dhs_years, matching_dhs_months, 
                                          dhs_rescale_multipliers, itn_variables, i_var)
    
    # save this coverage variable to data frame
    if(i_var == 1){
      dist_coverage = dist_coverage_var
    } else{
      dist_coverage = merge(dist_coverage, dist_coverage_var, by=c('admin_name', 'dist_year','net_life_lognormal_mu', 'net_life_lognormal_sigma', 'dhs_year', 'month'))
    }
  }
  
  # break into two dataframes: one for distributions occurring during the burnin and one for distributions occurring during the main calibration run
  dist_coverage_burn = dist_coverage[dist_coverage$dist_year < calib_start_year,]
  dist_coverage_calib = dist_coverage[dist_coverage$dist_year >= calib_start_year,]
  
  # add in the day of the simulation each intervention should start and the duration
  dist_coverage_burn$simday = (dist_coverage_burn$dist_year - burn_start_year) * 365 + round((dist_coverage_burn$month-1)*30+1)
  dist_coverage_calib$simday = (dist_coverage_calib$dist_year - calib_start_year) * 365 + round((dist_coverage_calib$month-1)*30+1)
  
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
  write.csv(dist_coverage_burn, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_burnin.csv'), row.names=FALSE)
  write.csv(dist_coverage_calib, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_main.csv'), row.names=FALSE)
}






# create one intervention file the calibration simulation (none in the burnin). Each seed comes from a different quantile.
create_transmission_calibration_itn_from_DHS = function(hbhi_dir, itn_variables,sim_start_year=2010, maximum_coverage=0.9,
                                                      itn_distribution_years = c(2010, 2011, 2014),  # since the 2010-2020 simulation starts at the first day of 2010, consider the 2009 distribution to occur on the first day of 2010 to get those nets in the simulation
                                                      itn_distribution_months = c(1, 7 ,7),
                                                      matching_dhs_years = c(2010,2012,2016),  # DHS survey year following mass ITN distributions (must match order of itn_distribution_years)
                                                      matching_dhs_months = c(10, 12, 10),  # chose the median month across all dhs entries, not specific for each admin
                                                      dhs_rescale_multipliers  # adjust net usage for seasonality
                                                      ){
  
  # read in sampled retention lognormal mus and sigmas
  net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'))
  net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
  
  # iterate through variables
  for(i_var in 1:length(itn_variables)){
    
    # get the sampled value from the DHS survey in each admin (num_samples values plus the mean value)
    sample_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_quantiles_', itn_variables[i_var],'.csv'))[,-1]
    # convert from wide to long format to make each seed its own row
    coverage_df = reshape2::melt(sample_df, id.vars=c("admin_name", "year"))
    colnames(coverage_df)[colnames(coverage_df)=='variable'] = 'seed'
    # colnames(coverage_df)[colnames(coverage_df)=='value'] = itn_variables[i_var]
    # replace "sample_" with "" in the seed column
    coverage_df$seed = gsub('sample_','',coverage_df$seed)
    # merge in the net halflife for each seed
    coverage_df = merge(coverage_df, net_discard_decay, by='seed')
    
    # iterate through distribution years, calculating the distribution coverages from the coverage at the next DHS survey year and the retention times
    dist_coverage_var = get_dist_coverage(coverage_df, maximum_coverage, itn_distribution_years, itn_distribution_months, matching_dhs_years, matching_dhs_months, 
                      dhs_rescale_multipliers, itn_variables, i_var)
    # save this coverage variable to data frame
    if(i_var == 1){
      dist_coverage = dist_coverage_var
    } else{
      dist_coverage = merge(dist_coverage, dist_coverage_var, by=c('seed','admin_name','dist_year','net_life_lognormal_mu','net_life_lognormal_sigma', 'dhs_year', 'month'))
    }
  }
  
  # add in the day of the simulation each intervention should start and the duration
  dist_coverage$simday = (dist_coverage$dist_year - sim_start_year) * 365 + round(dist_coverage$month*30)
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
  write.csv(dist_coverage, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_transmission_calib.csv'), row.names=FALSE)
}







plot_itn_campaign_coverage = function(hbhi_dir, itn_variables, admin_shape, maximum_coverage=0.9, 
                                     itn_distribution_years = c(2010, 2011, 2014),  # since the 2010-2020 simulation starts at the first day of 2010, consider the 2009 distribution to occur on the first day of 2010 to get those nets in the simulation
                                     itn_distribution_months = c(1, 7 ,7),
                                     matching_dhs_years = c(2010,2012,2016),  # DHS survey year following mass ITN distributions (must match order of itn_distribution_years)
                                     matching_dhs_months = c(10, 12, 10),  # chose the median month across all dhs entries, not specific for each admin
                                     admin_name_plot=NA, colors_range_0_to_1=NA, plot_sim_type='main_simulation'){
  if(plot_sim_type=='main_simulation'){
      dist_coverage = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_toPresent.csv'))
      file_string = ''
  } else if(plot_sim_type=='seasonality_calib'){
      dist_coverage_burn = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_burnin.csv'))
      dist_coverage_calib = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_main.csv'))
      dist_coverage = rbind(dist_coverage_burn, dist_coverage_calib)
      file_string = '_season_calib'
  } else if(plot_sim_type=='transmission_calib'){
    dist_coverage = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_transmission_calib.csv'))
    file_string = '_transmission_calib'
  }

  ############################################################################################
  # plot sampled values for ITN distributions (for one and all admins)
  ############################################################################################
  if(!dir.exists(paste0(hbhi_dir,'/simulation_inputs/plots'))) dir.create(paste0(hbhi_dir,'/simulation_inputs/plots'))
  par(mfrow=c(1,1))
  
  # all admins/archetypes: results of ITN coverate sampling for time of DHS, also showing extrapolation back to ITN coverage at time of campaign
  itn_ages_plot=c('u5','5_10','10_15','15_20','o20')
  itn_ages_ylab=paste0(c('U5', '5-10','10-15','15-20','>20'), ' ITN use')
  for(i_age in 1:length(itn_ages_plot)){
    colname_dhs = paste0('dhs_itn_',itn_ages_plot[i_age])
    colname_dist = paste0('itn_',itn_ages_plot[i_age])
    
    nrow_plot = round(length(unique(dist_coverage$admin_name))^(1/2))
    ncol_plot = ceiling(length(unique(dist_coverage$admin_name))/nrow_plot)
    png(paste0(hbhi_dir,'/simulation_inputs/plots/itn_extrapolation_', itn_ages_plot[i_age],file_string,'.png'), width=max(4,3*ncol_plot), height=max(4,2.5*nrow_plot), units='in',res=300)
    par(mfrow=c(nrow_plot, ncol_plot))
    for(i_admin in 1:length(unique(dist_coverage$admin_name))){
      admin_name = unique(dist_coverage$admin_name)[i_admin]
      plot(NA, ylim=c(0,1), xlim=c(min(itn_distribution_years), max(matching_dhs_years)+1), bty='L', ylab=itn_ages_ylab[i_age], xlab='year', main=admin_name)
      for(yy in 1:length(itn_distribution_years)){
        # subset to admin and first distribution year
        dist_coverage_ad = dist_coverage[intersect(which(dist_coverage$admin_name == admin_name), which(dist_coverage$dist_year==itn_distribution_years[yy])),]
        # plot lines of coverage decline between DHS and distribution
        plot_times = seq((itn_distribution_years[yy]+itn_distribution_months[yy]/12), (matching_dhs_years[yy]+matching_dhs_months[yy]/12), 0.01)
        time_between_dist_dhs = as.integer(as.Date(paste0(matching_dhs_years[yy], '/',matching_dhs_months[yy],'/1')) - as.Date(paste0(itn_distribution_years[yy], '/',itn_distribution_months[yy],'/1')))
        for(ss in 1:nrow(dist_coverage_ad)){
          lines(plot_times, sapply((dist_coverage_ad[[colname_dhs]][ss] / (1-plnorm(time_between_dist_dhs, meanlog=dist_coverage_ad$net_life_lognormal_mu[ss], sdlog=dist_coverage_ad$net_life_lognormal_sigma[ss])) * (1-plnorm((plot_times-(itn_distribution_years[yy]+itn_distribution_months[yy]/12))*365, meanlog=dist_coverage_ad$net_life_lognormal_mu[ss], sdlog=dist_coverage_ad$net_life_lognormal_sigma[ss]))), min, maximum_coverage), col='thistle')
        }
        # dhs coverage observed/sampled
        points(rep(matching_dhs_years[yy]+matching_dhs_months[yy]/12, nrow(dist_coverage_ad)),dist_coverage_ad[[colname_dhs]], col=rgb(0.2,0.7,0.3), pch=20)
        # extrapolated coverage at distribution
        points(rep(itn_distribution_years[yy]+itn_distribution_months[yy]/12, nrow(dist_coverage_ad)),dist_coverage_ad[[colname_dist]], col=rgb(0.5,0.1,0.7), pch=20)
      }
    }
    dev.off()
  }

  if(plot_sim_type == 'main_simulation'){
    # plot seasonality of net usage
    png(paste0(hbhi_dir,'/simulation_inputs/plots/itn_use_seasonality.png'), width=5, height=3, units='in',res=300)
    par(mfrow=c(1,1))
    plot(c(1,rep(2:12, each=2),13),rep(seasonality_monthly_scalar, each=2), type='l', bty='n', ylim=c(0.5,1), ylab = 'seasonality multiplier', xlab='month')
    dev.off()
    
    # plot of sampling distribution for median net retention time
    png(paste0(hbhi_dir,'/simulation_inputs/plots/itn_median_retention_distribution.png'), width=5, height=3, units='in',res=300)
    x_vals = seq(0.8,1.9,0.01)
    plot(x_vals, dnorm(x=x_vals, mean=median_llin_retention_est, sd=((median_llin_retention_est-median_llin_retention_lowerCI)/2)), type='l', ylab='probability density', xlab='net retention time', bty='L', main=c('Distribution from which median','net retention time is sampled'))
    points(x=c(median_llin_retention_lowerCI, median_llin_retention_est, median_llin_retention_est +(median_llin_retention_est-median_llin_retention_lowerCI)), y=c(0,0,0), cex=5, pch='|', col='red')
    legend('topright', c('mean and 95% CI from Amelia'), lwd=c(5), col=c('red'), bty='n')
    # legend('topright', c('sampling distribution','mean and 95%CI from Amelia'), lwd=c(1,5), col=c('black','red'), bty='n')
    dev.off()
    
    if(!is.na(admin_name_plot)){
      # single admin: results of ITN coverate sampling for time of DHS
      png(paste0(hbhi_dir,'/simulation_inputs/plots/itn_u5_',admin_name_plot,'0.png'), width=5, height=3.5, units='in',res=300)
      plot(NA, ylim=c(0,1), xlim=c(min(itn_distribution_years), max(matching_dhs_years)+1), bty='L', ylab='U5 ITN use', xlab='year', main=admin_name_plot)
      for(yy in 1:length(itn_distribution_years)){
        # subset to admin and first distribution year
        dist_coverage_ad = dist_coverage[intersect(which(dist_coverage$admin_name == admin_name_plot), which(dist_coverage$dist_year==itn_distribution_years[yy])),]
        # dhs coverage observed/sampled
        points(rep(matching_dhs_years[yy]+matching_dhs_months[yy]/12, nrow(dist_coverage_ad)),dist_coverage_ad$dhs_itn_u5, col=rgb(0.2,0.7,0.3, alpha=0.25), pch=20, cex=0.8)
      }
      dev.off()
      
      # single admin: results of ITN coverate sampling for time of DHS, also showing extrapolation back to ITN coverage at time of campaign
      png(paste0(hbhi_dir,'/simulation_inputs/plots/itn_u5_',admin_name_plot,'1.png'), width=5, height=3.5, units='in',res=300)
      plot(NA, ylim=c(0,1), xlim=c(min(itn_distribution_years), max(matching_dhs_years)+1), bty='L', ylab='U5 ITN use', xlab='year', main=admin_name_plot)
      for(yy in 1:length(itn_distribution_years)){
        # subset to admin and first distribution year
        dist_coverage_ad = dist_coverage[intersect(which(dist_coverage$admin_name == admin_name_plot), which(dist_coverage$dist_year==itn_distribution_years[yy])),]
        # plot lines of coverage decline between DHS and distribution
        plot_times = seq((itn_distribution_years[yy]+itn_distribution_months[yy]/12), (matching_dhs_years[yy]+matching_dhs_months[yy]/12), 0.01)
        time_between_dist_dhs = as.integer(as.Date(paste0(matching_dhs_years[yy], '/',matching_dhs_months[yy],'/1')) - as.Date(paste0(itn_distribution_years[yy], '/',itn_distribution_months[yy],'/1')))
        for(ss in 1:nrow(dist_coverage_ad)){
          lines(plot_times, sapply((dist_coverage_ad$dhs_itn_u5[ss] / (1-plnorm(time_between_dist_dhs, meanlog=dist_coverage_ad$net_life_lognormal_mu[ss], sdlog=dist_coverage_ad$net_life_lognormal_sigma[ss])) * (1-plnorm((plot_times-(itn_distribution_years[yy]+itn_distribution_months[yy]/12))*365, meanlog=dist_coverage_ad$net_life_lognormal_mu[ss], sdlog=dist_coverage_ad$net_life_lognormal_sigma[ss]))), min, maximum_coverage), col='thistle')
        }
        # dhs coverage observed/sampled
        points(rep(matching_dhs_years[yy]+matching_dhs_months[yy]/12, nrow(dist_coverage_ad)),dist_coverage_ad$dhs_itn_u5, col=rgb(0.2,0.7,0.3, alpha=0.25), pch=20, cex=0.8)
        # extrapolated coverage at distribution
        points(rep(itn_distribution_years[yy]+itn_distribution_months[yy]/12, nrow(dist_coverage_ad)),dist_coverage_ad$itn_u5, col=rgb(0.5,0.1,0.7, alpha=0.25), pch=20, cex=0.8)
      }
      dev.off()
    }
  }
  
  
  
}






plot_itn_campaign_coverage_map = function(hbhi_dir, itn_variables, admin_shape, 
                                          admin_name_plot=NA, colors_range_0_to_1=NA, plot_sim_type='main_simulation'){
  if(plot_sim_type=='main_simulation'){
    dist_coverage = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_toPresent.csv'))
    file_string = ''
  } else if(plot_sim_type=='seasonality_calib'){
    dist_coverage_burn = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_burnin.csv'))
    dist_coverage_calib = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_main.csv'))
    dist_coverage = rbind(dist_coverage_burn, dist_coverage_calib)
    file_string = '_season_calib'
  } else if(plot_sim_type=='transmission_calib'){
    dist_coverage = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_transmission_calib.csv'))
    file_string = '_transmission_calib'
  }
  
  ############################################################################################
  # plot sampled values for ITN distributions (for one and all admins)
  ############################################################################################
  if(!dir.exists(paste0(hbhi_dir,'/simulation_inputs/plots'))) dir.create(paste0(hbhi_dir,'/simulation_inputs/plots'))
  par(mfrow=c(1,1))
  
  # all admins/archetypes: results of ITN coverate sampling for time of DHS, also showing extrapolation back to ITN coverage at time of campaign
  itn_ages_plot=c('u5','5_10','10_15','15_20','o20')
  itn_ages_ylab=paste0(c('U5', '5-10','10-15','15-20','>20'), ' ITN use')

  itn_distribution_years = sort(unique(dist_coverage$year))
  ############################################################################################
  # plot map of ITN use at time of distribution
  ############################################################################################
  if(plot_sim_type=='main_simulation'){
    if(any(is.na(colors_range_0_to_1))){
      colors_range_0_to_1 = add.alpha(pals::parula(101), alpha=0.5)
    }
    if('seed' %in% colnames(dist_coverage)){
      dist_coverage_exp = dist_coverage[dist_coverage$seed == 1,]
    } else{
      dist_coverage_exp = dist_coverage
    }
    
    # par(mfrow=c(1,length(itn_distribution_years)))
    for(yy in 1:length(itn_distribution_years)){
      dist_coverage_yy_0 = dist_coverage_exp[dist_coverage_exp$year==itn_distribution_years[yy],]
      # add in NA coverage for admins without distributions in that year
      dist_coverage_yy_0 = merge(dist_coverage_yy_0, data.frame('admin_name'=admin_shape$NOMDEP), all=TRUE)
      reorder_admins = match(admin_shape$NOMDEP, dist_coverage_yy_0$admin_name)
      dist_coverage_yy = dist_coverage_yy_0[reorder_admins,]
      if(all(dist_coverage_yy$admin_name == admin_shape$NOMDEP)){
        png(paste0(hbhi_dir,'/simulation_inputs/plots/itn_campaign_u5_coverage_map_',itn_distribution_years[yy], file_string,'.png'), width=5, height=4, units='in',res=300)
        layout(matrix(c(1,1,1,2, 1,1,1,3),nrow=2, byrow=TRUE))
        admin_colors = colors_range_0_to_1[1+round(dist_coverage_yy$itn_u5*100)]
        plot(admin_shape, col = admin_colors, main=itn_distribution_years[yy])
        
        # legend - colorbar
        legend_image = as.raster(matrix(rev(colors_range_0_to_1[1+round(seq(0,1,length.out=20)*100)]), ncol=1))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'U5 ITN use')
        text(x=1.5, y = seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        rasterImage(legend_image, 0, 0, 1,1)
        dev.off()
      }
    }
  }
}











############################################################################################################################
############################################################################################################################
# get coverage from WHO-provided numbers of nets distributed per capita and age-specific use rates
############################################################################################################################
############################################################################################################################

create_mass_itn_input_from_npc = function(hbhi_dir, ds_pop_df_filename, itn_mass_coverage_filepath,  
                                          num_samples=50, sim_start_year=2010, season_calib_start_year_burn=1961, season_calib_start_year=2011, exclude_calib_itn_mass_after=2016, trans_calib_start_year=2010,
                                          U5_use_given_access=1, 
                                          npc_access_param1_vals = c(0.98,0.85, 0.92,0.98,0.99),  # first value is best estimate, then two lower and two upper bounds from lowest to highest
                                          npc_access_param2_vals = c(3, 1.8, 2.2, 5, 9),
                                          age_group_names=c('u5', '5_10', '10_15' ,'15_20', 'o20'),
                                          relative_use_by_age=c(1, 0.735771738, 0.684344698, 0.72043687, 1)
){
  ds_pop_df = read.csv(ds_pop_df_filename)
  itn_mass_coverage = read.csv(itn_mass_coverage_filepath)
  itn_mass_coverage$dist_date = as.Date(itn_mass_coverage$campaign_date, '%m/%d/%Y')
  
  # check that all of the mass campaign DS are in the DS_pop file
  if(!all(itn_mass_coverage$admin_name %in% ds_pop_df$admin_name)) warning("POTENTIAL ISSUE: some admin names from the mass ITN distribution file did not match the admin names used for the simulation")
  
  # read in sampled retention lognormal mus and sigmas
  net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
  net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
  
  # # create plot of nets per capita
  # # itn_mass_coverage2 = itn_mass_coverage
  # # itn_mass_coverage2 = itn_mass_coverage2[itn_mass_coverage2$dist_date>as.Date('2014-01-01'),]
  # # jitter = runif(length(itn_mass_coverage2$dist_date),-80,80)
  # # # plot(itn_mass_coverage2$dist_date + jitter, itn_mass_coverage2$nets_per_capita, pch=20, bty='L', ylim=c(0.2,0.8), ylab='nets per capita', xlab='year',cex=2, col=rgb(0.3,0.7,0.6, 0.7))
  # # plot(itn_mass_coverage2$dist_date + jitter, itn_mass_coverage2$nets_per_capita, pch=20, bty='L', ylim=c(0.2,0.8), ylab='nets per capita', xlab='year',cex=1.5, col=rgb(0.3,0.7,0.6))
  # 
  # itn_mass_coverage2$year = lubridate::year(itn_mass_coverage2$dist_date)
  # ggplot(itn_mass_coverage2, aes(x=year, y=nets_per_capita, group=dist_date)) + 
  #   # geom_jitter(color=rgb(0.3,0.7,0.6), size=3, alpha=0.4) +
  #   geom_jitter(color=rgb(0.65,0.9,0.88), size=3) +
  #   geom_boxplot(fill=NA) + 
  #   ylab('nets per capita') + 
  #   xlab('year of distribution') + 
  #   theme_classic() +
  #   theme(
  #     legend.position="none",
  #     plot.title = element_blank())
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  ######   main simulation input files   ##########
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #  
  itn_mass_coverage_main = itn_mass_coverage
  
  # add simulation start day
  itn_mass_coverage_main$simday = as.numeric(itn_mass_coverage_main$dist_date - as.Date(paste0('1/1/',sim_start_year), '%m/%d/%Y'))
  
  # expand to more seeds
  itn_mass_coverage_main_with_seed = itn_mass_coverage_main %>% slice(rep(1:n(), each=(num_samples+1)))
  itn_mass_coverage_main_with_seed$seed = rep(1:(num_samples+1), times=nrow(itn_mass_coverage_main))
  itn_mass_coverage_main = itn_mass_coverage_main_with_seed
  
  # calculate expected access given reported npc for each seed
  itn_mass_coverage_main$access = NA
  for(ff in 1:length(npc_access_param1_vals)){
    rows_ff = which(((itn_mass_coverage_main$seed) %% length(npc_access_param1_vals)) == (ff %% length(npc_access_param1_vals)))
    itn_mass_coverage_main$access[rows_ff] = calc_access_from_npc(access=itn_mass_coverage_main$nets_per_capita[rows_ff], npc_access_param1=npc_access_param1_vals[ff], npc_access_param2=npc_access_param2_vals[ff])
  }
  
  # estimate expected high-season U5 use given access
  itn_mass_coverage_main$base_use = itn_mass_coverage_main$access * U5_use_given_access
  
  # add in age-specific use given U5 use and relative use rates
  for(aa in 1:length(age_group_names)){
    itn_mass_coverage_main[[paste0('itn_', age_group_names[aa])]] = itn_mass_coverage_main$base_use * relative_use_by_age[aa]
  }
  
  # add net retention mus and sigmas
  itn_mass_coverage_main = merge(itn_mass_coverage_main, net_discard_decay, by=c('seed'), all.x=TRUE)
  
  # only keep relevant columns
  itn_mass_coverage_main$year = lubridate::year(itn_mass_coverage_main$dist_date)
  itn_mass_coverage_main = arrange(itn_mass_coverage_main, simday,admin_name,seed)
  itn_mass_coverage_main = itn_mass_coverage_main[,which(colnames(itn_mass_coverage_main) %in% c('admin_name', 'seed', 'nets_per_capita', 'dist_date', 'year', 'simday', 'net_life_lognormal_mu', 'net_life_lognormal_sigma', paste0('itn_',age_group_names)))]
  
  
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
  ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
  write.csv(itn_mass_coverage_main, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_toPresent.csv'), row.names=FALSE)
  
  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  ######   seasonality calibration input files   ##########
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - # 
  # get population-weighted average NPC in each seasonality archetype
  itn_mass_coverage_s_calib = merge(itn_mass_coverage, ds_pop_df[,c('admin_name','seasonality_archetype', 'pop_size')], by='admin_name')
  itn_mass_coverage_s_calib = itn_mass_coverage_s_calib[,-which(colnames(itn_mass_coverage_s_calib) %in% c('admin_name','campaign_date'))] %>% group_by(seasonality_archetype, dist_date) %>%
    summarise(nets_per_capita = weighted.mean(nets_per_capita, pop_size))
  colnames(itn_mass_coverage_s_calib)[colnames(itn_mass_coverage_s_calib)=='seasonality_archetype'] = 'admin_name'
  
  # calculate expected access given reported npc for each seed, using the best-guess estimate (the first value in the npc_access vectors)
  itn_mass_coverage_s_calib$access = calc_access_from_npc(access=itn_mass_coverage_s_calib$nets_per_capita, npc_access_param1=npc_access_param1_vals[1], npc_access_param2=npc_access_param2_vals[1])
  
  # estimate expected high-season U5 use given access
  itn_mass_coverage_s_calib$base_use = itn_mass_coverage_s_calib$access * U5_use_given_access
  
  # add in age-specific use given U5 use and relative use rates
  for(aa in 1:length(age_group_names)){
    itn_mass_coverage_s_calib[[paste0('itn_', age_group_names[aa])]] = itn_mass_coverage_s_calib$base_use * relative_use_by_age[aa]
  }
  
  # add best-guess net retention mu and sigma
  itn_mass_coverage_s_calib$net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[net_discard_decay$seed == 1]
  itn_mass_coverage_s_calib$net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[net_discard_decay$seed == 1]
  
  # only keep relevant columns
  itn_mass_coverage_s_calib$year = lubridate::year(itn_mass_coverage_s_calib$dist_date)
  itn_mass_coverage_s_calib = arrange(itn_mass_coverage_s_calib, year, admin_name)
  itn_mass_coverage_s_calib = itn_mass_coverage_s_calib[,which(colnames(itn_mass_coverage_s_calib) %in% c('admin_name','nets_per_capita', 'dist_date', 'year', 'simday', 'net_life_lognormal_mu', 'net_life_lognormal_sigma', paste0('itn_',age_group_names)))]
  
  # remove distributions after the last desired calibration ITN mass distribution (sometimes want to do this to avoid transient net behavior during calibration year)
  itn_mass_coverage_s_calib = itn_mass_coverage_s_calib[itn_mass_coverage_s_calib$year < exclude_calib_itn_mass_after,]
  
  # break into two dataframes: one for distributions occurring during the burnin and one for distributions occurring during the main calibration run
  itn_mass_coverage_s_calib_burn = itn_mass_coverage_s_calib[itn_mass_coverage_s_calib$year < season_calib_start_year,]
  itn_mass_coverage_s_calib_main = itn_mass_coverage_s_calib[itn_mass_coverage_s_calib$year >= season_calib_start_year,]
  
  # add simulation start day
  itn_mass_coverage_s_calib_main$simday = as.numeric(itn_mass_coverage_s_calib_main$dist_date - as.Date(paste0('1/1/', season_calib_start_year), '%m/%d/%Y'))
  itn_mass_coverage_s_calib_burn$simday = as.numeric(itn_mass_coverage_s_calib_burn$dist_date - as.Date(paste0('1/1/', season_calib_start_year_burn), '%m/%d/%Y'))
  
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  write.csv(itn_mass_coverage_s_calib_burn, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_burnin.csv'), row.names=FALSE)
  write.csv(itn_mass_coverage_s_calib_main, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_main.csv'), row.names=FALSE)
  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  ######   transmission calibration input files   ##########
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #  
  itn_mass_coverage_t_calib = itn_mass_coverage
  
  # add simulation start day
  itn_mass_coverage_t_calib$simday = as.numeric(itn_mass_coverage_t_calib$dist_date - as.Date(paste0('1/1/',trans_calib_start_year), '%m/%d/%Y'))
  
  
  # read in quantile retention lognormal mus and sigmas
  net_discard_decay_quant = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'))
  net_discard_decay_quant = net_discard_decay_quant[,which(colnames(net_discard_decay_quant) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
  
  # expand to more seeds
  itn_mass_coverage_t_calib_with_seed = itn_mass_coverage_t_calib %>% slice(rep(1:n(), each=max(net_discard_decay_quant$seed)))
  itn_mass_coverage_t_calib_with_seed$seed = rep(1:max(net_discard_decay_quant$seed), times=nrow(itn_mass_coverage_t_calib))
  itn_mass_coverage_t_calib = itn_mass_coverage_t_calib_with_seed
  
  # calculate expected access given reported npc for each seed
  # re-order npc-access params from smallest to largest for quantiles
  npc_access_param1_vals_2 = sort(npc_access_param1_vals)
  npc_access_param2_vals_2 = sort(npc_access_param2_vals)
  itn_mass_coverage_t_calib$access = NA
  for(ff in 1:length(npc_access_param1_vals_2)){
    rows_ff = which(((itn_mass_coverage_t_calib$seed) %% length(npc_access_param1_vals_2)) == (ff %% length(npc_access_param1_vals_2)))
    itn_mass_coverage_t_calib$access[rows_ff] = calc_access_from_npc(access=itn_mass_coverage_t_calib$nets_per_capita[rows_ff], npc_access_param1=npc_access_param1_vals_2[ff], npc_access_param2=npc_access_param2_vals_2[ff])
  }
  
  # estimate expected high-season U5 use given access
  itn_mass_coverage_t_calib$base_use = itn_mass_coverage_t_calib$access * U5_use_given_access
  
  # add in age-specific use given U5 use and relative use rates
  for(aa in 1:length(age_group_names)){
    itn_mass_coverage_t_calib[[paste0('itn_', age_group_names[aa])]] = itn_mass_coverage_t_calib$base_use * relative_use_by_age[aa]
  }
  
  # add net retention mus and sigmas
  itn_mass_coverage_t_calib = merge(itn_mass_coverage_t_calib, net_discard_decay_quant, by=c('seed'), all.x=TRUE)
  
  # only keep relevant columns
  itn_mass_coverage_t_calib$year = lubridate::year(itn_mass_coverage_t_calib$dist_date)
  itn_mass_coverage_t_calib = arrange(itn_mass_coverage_t_calib, simday,admin_name,seed)
  itn_mass_coverage_t_calib = itn_mass_coverage_t_calib[,which(colnames(itn_mass_coverage_t_calib) %in% c('admin_name', 'seed', 'nets_per_capita', 'dist_date', 'year', 'simday', 'net_life_lognormal_mu', 'net_life_lognormal_sigma', paste0('itn_',age_group_names)))]
  
  
  # save starting ITN use at time of distribution. still need to add kill, blocking rates before it can be consumed by simulation scripts
  write.csv(itn_mass_coverage_t_calib, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_transmission_calib.csv'), row.names=FALSE)
  
}




plot_mass_itn_input_from_npc = function(){
  itn_mass_coverage_main = read.csv(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_2020.csv'))
  # write.csv(itn_mass_coverage_s_calib_burn, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_burnin.csv'), row.names=FALSE)
  # write.csv(itn_mass_coverage_s_calib_main, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_season_calib_main.csv'), row.names=FALSE)
  # write.csv(itn_mass_coverage_t_calib, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_transmission_calib.csv'), row.names=FALSE)
  # 
  
  # subset to single seed
  itn_mass_coverage_main_s1 = itn_mass_coverage_main[itn_mass_coverage_main$seed ==1,]
  jitter = runif(nrow(itn_mass_coverage_main_s1), min=-20, max=20)
  plot(as.Date(itn_mass_coverage_main_s1$dist_date)+jitter, itn_mass_coverage_main_s1$itn_u5, type='p', pch=20, bty='L', xlab='date',ylab='U5 access from mass distribution')
}










# ############################################################################################
# # plot relative rates of ITN use in different age groups at country level
# ############################################################################################
# 
# par(mfrow=c(2,2))
# save_weighted_ratios = matrix(NA, nrow=0,ncol=4)
# save_equal_weights = matrix(NA, nrow=0,ncol=4)
# for(yy in 1:length(years)){
#   clust_sums = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_cluster_outputs_', years[yy], '.csv'))[,-1]
#   clust_sums$rel_5_10 = clust_sums$itn_5_10_rate / clust_sums$itn_u5_rate
#   clust_sums$rel_10_15 = clust_sums$itn_10_15_rate / clust_sums$itn_u5_rate
#   clust_sums$rel_15_20 = clust_sums$itn_15_20_rate / clust_sums$itn_u5_rate
#   clust_sums$rel_o20 = clust_sums$itn_o20_rate / clust_sums$itn_u5_rate
#   
#   
#   # remove any clusters with Inf
#   clust_sums_2 = clust_sums[which(clust_sums$rel_5_10<Inf),]
#   
#   # par(mfrow=c(2,2))
#   # hist(clust_sums_2$rel_5_10)
#   # hist(clust_sums_2$rel_10_15)
#   # hist(clust_sums_2$rel_15_20)
#   # hist(clust_sums_2$rel_o20)
#   
#   # unweighted mean of all cluster values
#   unweighted_cluster_means = c(
#     mean(clust_sums_2$rel_5_10, na.rm=TRUE),
#     mean(clust_sums_2$rel_10_15, na.rm=TRUE),
#     mean(clust_sums_2$rel_15_20, na.rm=TRUE),
#     mean(clust_sums_2$rel_o20, na.rm=TRUE)
#   )
#   # weighted mean of cluster values
#   weighted_cluster_means = c(
#     weighted.mean(clust_sums_2$rel_5_10, clust_sums_2$itn_weights, na.rm=TRUE),
#     weighted.mean(clust_sums_2$rel_10_15, clust_sums_2$itn_weights, na.rm=TRUE),
#     weighted.mean(clust_sums_2$rel_15_20, clust_sums_2$itn_weights, na.rm=TRUE),
#     weighted.mean(clust_sums_2$rel_o20, clust_sums_2$itn_weights, na.rm=TRUE)
#   )
#   # without aggregating or weighting by cluster, equal weight to all individuals included in survey
#   unweighted_survey_means = c(
#     (sum(clust_sums$itn_5_10_num_true) / sum(clust_sums$itn_5_10_num_total) ) / (sum(clust_sums$itn_u5_num_true) / sum(clust_sums$itn_u5_num_total)),
#     (sum(clust_sums$itn_10_15_num_true) / sum(clust_sums$itn_10_15_num_total) ) / (sum(clust_sums$itn_u5_num_true) / sum(clust_sums$itn_u5_num_total)),
#     (sum(clust_sums$itn_15_20_num_true) / sum(clust_sums$itn_15_20_num_total) ) / (sum(clust_sums$itn_u5_num_true) / sum(clust_sums$itn_u5_num_total)),
#     (sum(clust_sums$itn_o20_num_true) / sum(clust_sums$itn_o20_num_total) ) / (sum(clust_sums$itn_u5_num_true) / sum(clust_sums$itn_u5_num_total))
#   )
#   save_weighted_ratios = rbind(save_weighted_ratios, weighted_cluster_means)
#   save_equal_weights = rbind(save_equal_weights, unweighted_survey_means)
#   ratios = matrix(c(unweighted_cluster_means, weighted_cluster_means, unweighted_survey_means), byrow=TRUE, nrow=3)
#   barplot(ratios, beside=TRUE, col=c(rgb(0.3,0.3,1), rgb(0.6,0.6,1), rgb(0.9,0.9,1)), names.arg=c('5-10','10-15','15-20','>20'), xlab='age group',ylab='coverage relative to U5', main= paste0('net usage relative to U5 - ', years[yy]))
# }
# plot(NA,ylab='', xlab='', xlim=c(0,1), ylim=c(0,1),axes=FALSE)
# legend('left', c('unweighted mean of clusters', 'DHS-weighted mean of clusters', 'equal weights for all individuals'), col=c(rgb(0.3,0.3,1), rgb(0.6,0.6,1), rgb(0.9,0.9,1)), pch=15, bty='n')
# par(mfrow=c(1,1))
# 
# colMeans(save_weighted_ratios)
# colMeans(save_equal_weights)
# 
# 
# 
# 
# 
# # get DHS-household-weighted and unweighted ITN use rates in each age group
# for(yy in 1:length(years)){
#   year = years[yy]
#   DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',year,'_files_recodes_for_sims.csv'))
# 
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   # how similar are region, country weighted versus non-weighted values? uses non-cluster-aggregated values, so won't work for comparing coverage relative to U5 coverage
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   var_index = which(DHS_file_recode_df$variable == 'itn')
#   if(!is.na(DHS_file_recode_df$filename[var_index])){
#     cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
#     cur_dta = cur_dta[cur_dta[[DHS_file_recode_df$age_code[var_index]]]<=5,]
#     cur_dta$pos = NA
#     cur_dta$pos[which(cur_dta[,grep(colnames(cur_dta), pattern = DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$pos_pattern[var_index])] = 1
#     cur_dta$pos[which(cur_dta[,grep(colnames(cur_dta), pattern = DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$neg_pattern[var_index])] = 0
#     cur_dta$wt = cur_dta$hv005/1000000
#     
#     # Tabulate indicator by region
#     ddply(cur_dta,~hv024,summarise,mean=weighted.mean(pos, wt))
#     
#     # tabulate weighted mean for entire country
#     weighted.mean(cur_dta$pos, cur_dta$wt)
#     
#     # tabulate unweighted mean for entire country
#     mean(cur_dta$pos)
#   }
#   
#   
# }













