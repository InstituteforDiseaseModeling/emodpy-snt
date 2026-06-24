# process_sim_input_functions.R


library(data.table)
library(dplyr)
library(tidyverse)
library(lubridate)
# library(rgdal)
library(lubridate)




#####################################################################################################
# ================================================================================================= #
# calculate state- or LGA-level intervention  timeseries for included scenarios, based on inputs
# ================================================================================================= #
#####################################################################################################

##############
# CM
##############

# ----- admin- or state-level CM intervention coverage ----- #
get_cm_timeseries_by_state = function(cm_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(cm_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  input_df = merge(input_df, admin_info)
  
  # CM is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_cm_years = unique(input_df$year)
    for(yy in cur_cm_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  # Per-row extrapolation of duration == -1 rows up to end_year. See the
  # parallel block in get_vacc_timeseries_by_state for the reasoning -- the
  # previous version (using cm_repeated$year[1] for the loop length and a
  # global max(year) < end_year gate) under-extrapolated some admins and
  # would silently skip extrapolation when any admin's rows already reached
  # end_year.
  if(any(input_df$duration == -1)){
    cm_repeated = input_df[input_df$duration == -1, ]
    n_extra     = pmax(end_year - cm_repeated$year, 0L)
    if(any(n_extra > 0)){
      rep_idx          = rep(seq_len(nrow(cm_repeated)), n_extra)
      offset           = sequence(n_extra)
      cm_extended      = cm_repeated[rep_idx, , drop = FALSE]
      cm_extended$year   = cm_extended$year   + offset
      cm_extended$simday = cm_extended$simday + offset * 365
      input_df = rbind(input_df, cm_extended)
    }
  }
  
  if(get_state_level){
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise_all(mean) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    # get population-weighted CM coverage across admins
    input_df$multiplied_U5_cm = input_df$U5_coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied CM coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cm, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cm / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    cm_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise_all(mean) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    input_df_agg_admin = input_df
    
    # take average across simulation seeds
    cm_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  cm_agg$scenario = exp_name
  return(cm_agg)
}


# ----- national-level CM intervention coverage ----- #
get_cm_timeseries_exp = function(cm_filepath, pop_sizes, end_year, exp_name, cur_admins, min_year, plot_by_month=TRUE){
  
  input_df = read.csv(cm_filepath)
  # subset to appropriate admins
  input_df = input_df[input_df$admin_name %in% cur_admins,] 
  if(!('seed' %in% input_df)) input_df$seed = 1
  
  # CM is sometimes repeated for several years but only listed once; change to repeate the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_cm_years = unique(input_df$year)
    for(yy in cur_cm_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  # Per-row extrapolation of duration == -1 rows up to end_year. See the
  # parallel block in get_vacc_timeseries_by_state for the reasoning -- the
  # previous version (using cm_repeated$year[1] for the loop length and a
  # global max(year) < end_year gate) under-extrapolated some admins and
  # would silently skip extrapolation when any admin's rows already reached
  # end_year.
  if(any(input_df$duration == -1)){
    cm_repeated = input_df[input_df$duration == -1, ]
    n_extra     = pmax(end_year - cm_repeated$year, 0L)
    if(any(n_extra > 0)){
      rep_idx          = rep(seq_len(nrow(cm_repeated)), n_extra)
      offset           = sequence(n_extra)
      cm_extended      = cm_repeated[rep_idx, , drop = FALSE]
      cm_extended$year   = cm_extended$year   + offset
      cm_extended$simday = cm_extended$simday + offset * 365
      input_df = rbind(input_df, cm_extended)
    }
  }
  # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
  input_df <- input_df %>% group_by(year, admin_name, seed) %>%
    summarise_all(mean) %>% ungroup()
  
  input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
  
  
  # get population-weighted CM coverage across admins
  input_df = merge(input_df, pop_sizes)
  input_df$multiplied_U5_cm = input_df$U5_coverage * input_df$pop_size
  
  # get sum of population sizes and multiplied CM coverage across included admins
  input_df_agg_admin <- input_df %>% dplyr::select(year, seed, multiplied_U5_cm, pop_size) %>% group_by(year, seed) %>%
    summarise_all(sum) %>% ungroup()
  # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
  input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cm / input_df_agg_admin$pop_size
  
  
  # take average, max, and min burdens across simulation seeds
  if(plot_by_month){
    # subdivide year values and add dates
    # date dataframe
    included_years = unique(input_df_agg_admin$year)
    all_months = as.Date(paste0(rep(included_years, each=12),'-',c('01','02','03','04','05','06','07','08','09','10','11','12'), '-01' ))
    date_df = data.frame(year=rep(included_years, each=12), date=all_months)
    input_df_agg_admin_monthly = merge(input_df_agg_admin, date_df, by='year', all.x=TRUE, all.y=TRUE)
    cm_agg = as.data.frame(input_df_agg_admin_monthly) %>% dplyr::select(year, date, U5_coverage) %>%
      dplyr::group_by(year, date) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    cm_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, U5_coverage) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  
  cm_agg$scenario = exp_name
  return(cm_agg)
  
}




##############
# ANC ITN
##############

# ----- admin- or state-level ANC intervention coverage ----- #
get_itn_anc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  if('seasonality_archetype' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-seasonality_archetype)
  if('State' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-State)
  if('pop_size' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-pop_size)
  
  input_df = merge(input_df, admin_info)
  
  # ANC ITN is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_years = unique(input_df$year)
    for(yy in cur_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  # Extrapolate every duration == -1 row forward to end_year, per row. The
  # previous version used df_repeated$year[1] for the loop length and applied
  # the same offset to every row, which (1) made some admins fall short of
  # end_year (line ends at 2034 instead of 2035) when their derived year was
  # later than the first row's, and (2) was gated by a global
  # max(input_df$year) < end_year check that incorrectly skipped extrapolation
  # whenever a single admin's rescaled rows already reached end_year while
  # other admins still had a duration=-1 row at projection_start_year.
  # n_extra is computed per row so each row extends exactly to end_year and
  # rows already at/past end_year produce no extras.
  if(any(input_df$duration == -1)){
    df_repeated = input_df[input_df$duration == -1, ]
    n_extra     = pmax(end_year - df_repeated$year, 0L)
    if(any(n_extra > 0)){
      rep_idx          = rep(seq_len(nrow(df_repeated)), n_extra)
      offset           = sequence(n_extra)
      df_extended      = df_repeated[rep_idx, , drop = FALSE]
      df_extended$year   = df_extended$year   + offset
      df_extended$simday = df_extended$simday + offset * 365
      input_df = rbind(input_df, df_extended)
    }
  }
  
  if(get_state_level){
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise(coverage=mean(coverage), pop_size=mean(pop_size)) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    # get population-weighted CM coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied CM coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise(coverage=mean(coverage)) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    input_df_agg_admin = input_df
    
    # take average across simulation seeds
    input_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_agg$scenario = exp_name
  return(input_agg)
}




##############
# EPI ITN
##############

# ----- admin- or state-level AEPINC intervention coverage ----- #
get_itn_epi_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  if('seasonality_archetype' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-seasonality_archetype)
  if('State' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-State)
  if('pop_size' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-pop_size)
  
  input_df = merge(input_df, admin_info)
  
  # EPI ITN is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_years = unique(input_df$year)
    for(yy in cur_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  # Extrapolate every duration == -1 row forward to end_year, per row. The
  # previous version used df_repeated$year[1] for the loop length and applied
  # the same offset to every row, which (1) made some admins fall short of
  # end_year (line ends at 2034 instead of 2035) when their derived year was
  # later than the first row's, and (2) was gated by a global
  # max(input_df$year) < end_year check that incorrectly skipped extrapolation
  # whenever a single admin's rescaled rows already reached end_year while
  # other admins still had a duration=-1 row at projection_start_year.
  # n_extra is computed per row so each row extends exactly to end_year and
  # rows already at/past end_year produce no extras.
  if(any(input_df$duration == -1)){
    df_repeated = input_df[input_df$duration == -1, ]
    n_extra     = pmax(end_year - df_repeated$year, 0L)
    if(any(n_extra > 0)){
      rep_idx          = rep(seq_len(nrow(df_repeated)), n_extra)
      offset           = sequence(n_extra)
      df_extended      = df_repeated[rep_idx, , drop = FALSE]
      df_extended$year   = df_extended$year   + offset
      df_extended$simday = df_extended$simday + offset * 365
      input_df = rbind(input_df, df_extended)
    }
  }
  
  if(get_state_level){
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise(coverage=mean(coverage), pop_size=mean(pop_size)) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    # get population-weighted CM coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied CM coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
    input_df <- input_df %>% group_by(year, admin_name, State, seed) %>%
      summarise(coverage=mean(coverage)) %>% ungroup()
    
    # input_df = input_df[intersect(which(input_df$year >= min_year), which(input_df$year <= end_year)),]
    input_df_agg_admin = input_df
    
    # take average across simulation seeds
    input_agg = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_agg$scenario = exp_name
  return(input_agg)
}






##############
# SMC
##############
get_smc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  if('seasonality_archetype' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-seasonality_archetype)
  if('pop_size' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-pop_size)
  if('date_estimate' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-date_estimate)
  
  input_df$u5_coverage_total = input_df$coverage_high_access_U5 * input_df$high_access_U5 + input_df$coverage_low_access_U5 * (1-input_df$high_access_U5)
  
  # get mean across rounds within a year
  input_df = input_df %>% dplyr::select(-State) %>% 
    filter(year>=min_year, year<=end_year) %>%
    group_by(year, admin_name, seed) %>%
    summarise_all(mean) %>% ungroup()
  
  # create data frame with all admin-years as zeros if SMC isn't specified as given in input file
  admin_years = merge(admin_info, data.frame('year'=seq(min_year, end_year)), all=TRUE)
  input_df = merge(input_df, admin_years, by=c('admin_name','year'), all=TRUE)
  input_df$u5_coverage_total[is.na(input_df$u5_coverage_total)] = 0
  input_df$U5_coverage = input_df$u5_coverage_total
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_U5_cov = input_df$U5_coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  input_ave_df$scenario = exp_name

  return(input_ave_df)
}





##############
# SPAQ (combined SMC + IPTsc) -- coverage for ONE age band (min_age/max_age rows)
##############
# SPAQ files carry several age bands as separate rows: [0.25,5) & [5,10) from SMC,
# plus [5,12) & [12,16) from the school-based IPTsc fold-in. Select one band per call
# via band_min/band_max. Effective coverage =
#   coverage_high_access*high_access_frac + coverage_low_access*(1-high_access_frac).
# A band absent from a file (e.g. [5,12)/[12,16) in spaq_bau) yields all-zero coverage.
get_spaq_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year,
                                        band_min, band_max, get_state_level=TRUE){

  input_df = read.csv(input_filepath)
  if(!('seed' %in% colnames(input_df))) input_df$seed = 1
  if('seasonality_archetype' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-seasonality_archetype)
  if('pop_size' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-pop_size)
  if('date_estimate' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-date_estimate)

  # subset to the requested age band BEFORE any averaging (otherwise bands/rounds mix)
  input_df = input_df[input_df$min_age == band_min & input_df$max_age == band_max, , drop=FALSE]

  # effective per-round coverage for this band
  input_df$band_coverage = input_df$coverage_high_access * input_df$high_access_frac +
                            input_df$coverage_low_access * (1 - input_df$high_access_frac)

  # mean across rounds within a year
  input_df = input_df %>% dplyr::select(-State) %>%
    filter(year>=min_year, year<=end_year) %>%
    group_by(year, admin_name, seed) %>%
    summarise_all(mean) %>% ungroup()

  # fill all admin-years with 0 where this band's SPAQ isn't delivered in the input file
  admin_years = merge(admin_info, data.frame('year'=seq(min_year, end_year)), all=TRUE)
  input_df = merge(input_df, admin_years, by=c('admin_name','year'), all=TRUE)
  input_df$band_coverage[is.na(input_df$band_coverage)] = 0
  input_df$U5_coverage = input_df$band_coverage   # reuse the field name used by the aggregation below

  # calculate values within state or admin
  if(get_state_level){
    # population-weighted coverage across admins
    input_df$multiplied_U5_cov = input_df$U5_coverage * input_df$pop_size
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cov, pop_size) %>%
      group_by(year, State, seed) %>% summarise_all(sum) %>% ungroup()
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cov / input_df_agg_admin$pop_size
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  input_ave_df$scenario = exp_name

  return(input_ave_df)
}




##############
# ITN mass campaign
##############
get_itn_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1

  # get mean across rounds within a year
  input_df <- input_df %>% group_by(year, admin_name, seed) %>%
    summarise(itn_u5 = mean(itn_u5)) %>% ungroup() %>%
    filter(year<=end_year)
  
  # create data frame with all admin-years as zeros if SMC isn't specified as given in input file
  admin_years = merge(admin_info, data.frame('year'=seq(min_year, end_year)), all=TRUE)
  input_df = merge(input_df, admin_years, by=c('admin_name','year'), all=TRUE)
  input_df$itn_u5[is.na(input_df$itn_u5)] = 0
  input_df$U5_coverage = input_df$itn_u5
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_U5_cov = input_df$U5_coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_U5_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$U5_coverage = input_df_agg_admin$multiplied_U5_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, U5_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  input_ave_df$scenario = exp_name
  
  return(input_ave_df)
}






##############
# EPI vaccine
##############
get_vacc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  if('State' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-State)
  input_df = merge(input_df, admin_info)
  
  # look just at coverage for primary series
  input_df = input_df[input_df$vaccine == 'primary',]
  # true start day of administration is RTSS_day+rtss_touchpoint (due to how things are setup in campaign for birth-triggered delayed vaccine)
  input_df$simday = input_df$RTSS_day + input_df$rtss_touchpoint
  input_df$year = floor(input_df$simday/365) + min_year
  
  # EPI vaccine is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
  input_df$years_repeated = input_df$duration/365
  if(any(input_df$years_repeated>1)){
    cur_years = unique(input_df$year)
    for(yy in cur_years){
      # get first instance of this year
      cur_year = input_df[input_df$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          input_df = rbind(input_df, temp_year)
        }
      }
    }
  }
  # Extrapolate every duration == -1 row forward to end_year, per row. The
  # previous version used df_repeated$year[1] for the loop length and applied
  # the same offset to every row, which (1) made some admins fall short of
  # end_year (line ends at 2034 instead of 2035) when their derived year was
  # later than the first row's, and (2) was gated by a global
  # max(input_df$year) < end_year check that incorrectly skipped extrapolation
  # whenever a single admin's rescaled rows already reached end_year while
  # other admins still had a duration=-1 row at projection_start_year.
  # n_extra is computed per row so each row extends exactly to end_year and
  # rows already at/past end_year produce no extras.
  if(any(input_df$duration == -1)){
    df_repeated = input_df[input_df$duration == -1, ]
    n_extra     = pmax(end_year - df_repeated$year, 0L)
    if(any(n_extra > 0)){
      rep_idx          = rep(seq_len(nrow(df_repeated)), n_extra)
      offset           = sequence(n_extra)
      df_extended      = df_repeated[rep_idx, , drop = FALSE]
      df_extended$year   = df_extended$year   + offset
      df_extended$simday = df_extended$simday + offset * 365
      input_df = rbind(input_df, df_extended)
    }
  }
  
  
  # # create data frame with all admin-years as zeros before start of vaccine
  # admin_years = merge(admin_info, data.frame('year'=seq(min_year, vacc_start_year-1), coverage=0), all=TRUE)
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_ave_df$scenario = exp_name
  
  return(input_ave_df)
}





##############
# PMC
##############
get_pmc_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){
  
  input_df = read.csv(input_filepath)
  if(!('seed' %in% input_df)) input_df$seed = 1
  
  # true start day of administration is simday+pmc_touchpoints (due to how things are setup in campaign for birth-triggered delayed vaccine)
  input_df$simday = input_df$simday + input_df$pmc_touchpoints
  
  # get average coverage across touchpoints; for simday, take date of earliest touchpoint simday
  input_df <- input_df %>% group_by(admin_name, seed) %>%
    summarise(coverage = mean(coverage), simday=min(simday)) %>% 
    ungroup()  
  input_df$year = floor(input_df$simday/365) + min_year
  input_df = merge(input_df, admin_info)
  
  if('duration' %in% colnames(input_df)){
    # EPI PMC is sometimes repeated for several years but only listed once; change to repeat the appropriate number of times
    input_df$years_repeated = input_df$duration/365
    if(any(input_df$years_repeated>1)){
      cur_years = unique(input_df$year)
      for(yy in cur_years){
        # get first instance of this year
        cur_year = input_df[input_df$year==yy,]
        if(cur_year$years_repeated[1]>1){
          for(rr in 1:(cur_year$years_repeated[1] - 1)){
            temp_year = cur_year
            temp_year$year = cur_year$year + rr
            temp_year$simday = cur_year$simday + rr*365
            input_df = rbind(input_df, temp_year)
          }
        }
      }
    }
  } else{
    input_df$duration = -1  # default value
  }

  # Extrapolate every duration == -1 row forward to end_year, per row. The
  # previous version used df_repeated$year[1] for the loop length and applied
  # the same offset to every row, which (1) made some admins fall short of
  # end_year (line ends at 2034 instead of 2035) when their derived year was
  # later than the first row's, and (2) was gated by a global
  # max(input_df$year) < end_year check that incorrectly skipped extrapolation
  # whenever a single admin's rescaled rows already reached end_year while
  # other admins still had a duration=-1 row at projection_start_year.
  # n_extra is computed per row so each row extends exactly to end_year and
  # rows already at/past end_year produce no extras.
  if(any(input_df$duration == -1)){
    df_repeated = input_df[input_df$duration == -1, ]
    n_extra     = pmax(end_year - df_repeated$year, 0L)
    if(any(n_extra > 0)){
      rep_idx          = rep(seq_len(nrow(df_repeated)), n_extra)
      offset           = sequence(n_extra)
      df_extended      = df_repeated[rep_idx, , drop = FALSE]
      df_extended$year   = df_extended$year   + offset
      df_extended$simday = df_extended$simday + offset * 365
      input_df = rbind(input_df, df_extended)
    }
  }
  
  
  # # create data frame with all admin-years as zeros before start of vaccine
  # admin_years = merge(admin_info, data.frame('year'=seq(min_year, vacc_start_year-1), coverage=0), all=TRUE)
  
  # calculate values within state or admin
  if(get_state_level){
    # get population-weighted CM coverage across admins
    input_df$multiplied_cov = input_df$coverage * input_df$pop_size
    
    # get sum of population sizes and multiplied SMC coverage across included admins
    input_df_agg_admin <- input_df %>% dplyr::select(year, State, seed, multiplied_cov, pop_size) %>% group_by(year, State, seed) %>%
      summarise_all(sum) %>% ungroup()
    # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
    input_df_agg_admin$coverage = input_df_agg_admin$multiplied_cov / input_df_agg_admin$pop_size
    
    # get average within state and across seeds
    input_ave_df = as.data.frame(input_df_agg_admin) %>% dplyr::select(year, State, coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  } else{
    # take average across simulation seeds
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(coverage),
                       max_coverage = max(coverage),
                       min_coverage = min(coverage))
  }
  input_ave_df$scenario = exp_name
  
  return(input_ave_df)
}



##############
# IRS
##############
get_irs_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){

  input_df = read.csv(input_filepath)
  if(!('seed' %in% colnames(input_df))) input_df$seed = 1
  if('seasonality_archetype' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-seasonality_archetype)
  if('pop_size' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-pop_size)
  if('State' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-State)

  # filter to the simulation window
  input_df = input_df %>% filter(year >= min_year, year <= end_year)

  # admin-year skeleton so admins without IRS in a given year show coverage 0
  admin_years = merge(admin_info, data.frame('year' = seq(min_year, end_year)), all = TRUE)
  input_df = merge(input_df, admin_years, by = c('admin_name', 'year'), all = TRUE)
  input_df$effective_coverage[is.na(input_df$effective_coverage)] = 0
  input_df$IRS_coverage = input_df$effective_coverage

  if(get_state_level){
    # population-weighted state-level coverage
    input_df$multiplied_IRS_cov = input_df$IRS_coverage * input_df$pop_size
    input_df_agg <- input_df %>% dplyr::select(year, State, seed, multiplied_IRS_cov, pop_size) %>%
      dplyr::group_by(year, State, seed) %>%
      dplyr::summarise_all(sum) %>% ungroup()
    input_df_agg$IRS_coverage = input_df_agg$multiplied_IRS_cov / input_df_agg$pop_size

    input_ave_df = as.data.frame(input_df_agg) %>% dplyr::select(year, State, IRS_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(IRS_coverage),
                       max_coverage  = max(IRS_coverage),
                       min_coverage  = min(IRS_coverage))
  } else{
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, IRS_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(IRS_coverage),
                       max_coverage  = max(IRS_coverage),
                       min_coverage  = min(IRS_coverage))
  }
  input_ave_df$scenario = exp_name

  return(input_ave_df)
}



##############
# IPTp
##############
get_iptp_timeseries_by_state = function(input_filepath, admin_info, end_year, exp_name, min_year, get_state_level=TRUE){

  input_df = read.csv(input_filepath)
  if(!('seed' %in% colnames(input_df))) input_df$seed = 1
  if('State' %in% colnames(input_df)) input_df = input_df %>% dplyr::select(-State)

  # filter to the simulation window
  input_df = input_df %>% filter(year >= min_year, year <= end_year)

  # admin-year skeleton so admins without an IPTp row in a given year show 0
  admin_years = merge(admin_info, data.frame('year' = seq(min_year, end_year)), all = TRUE)
  input_df = merge(input_df, admin_years, by = c('admin_name', 'year'), all = TRUE)
  input_df$coverage[is.na(input_df$coverage)] = 0
  input_df$IPTp_coverage = input_df$coverage

  if(get_state_level){
    # population-weighted state-level coverage
    input_df$multiplied_IPTp_cov = input_df$IPTp_coverage * input_df$pop_size
    input_df_agg <- input_df %>% dplyr::select(year, State, seed, multiplied_IPTp_cov, pop_size) %>%
      dplyr::group_by(year, State, seed) %>%
      dplyr::summarise_all(sum) %>% ungroup()
    input_df_agg$IPTp_coverage = input_df_agg$multiplied_IPTp_cov / input_df_agg$pop_size

    input_ave_df = as.data.frame(input_df_agg) %>% dplyr::select(year, State, IPTp_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(IPTp_coverage),
                       max_coverage  = max(IPTp_coverage),
                       min_coverage  = min(IPTp_coverage))
  } else{
    input_ave_df = as.data.frame(input_df) %>% dplyr::select(year, IPTp_coverage, State, admin_name) %>%
      dplyr::group_by(year, State, admin_name) %>%
      dplyr::summarise(mean_coverage = mean(IPTp_coverage),
                       max_coverage  = max(IPTp_coverage),
                       min_coverage  = min(IPTp_coverage))
  }
  input_ave_df$scenario = exp_name

  return(input_ave_df)
}
