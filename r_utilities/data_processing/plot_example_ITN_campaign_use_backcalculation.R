# plot_example_ITN_campaign_use_backcalculation.R

#### setup inputs ###

# set filepaths
user = Sys.getenv("USERNAME")
user_path = file.path("C:/Users",user)
hbhi_dir = paste0(user_path, '/Gates Foundation Dropbox/Malaria Team Folder/projects/snt/Nigeria/snt_2025')

# read in sampled retention lognormal mus and sigmas and take the first (expected) value
net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[1]
net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[1]
# read in seasonality pattern
seasonality_monthly_scalar = read.csv(paste0(hbhi_dir, '/simulation_inputs/ITN_use_seasonality.csv'))$itn_use_scalar  # seasonality scalar for net use: gives rate of net use in each month relative to month with maximum use rate
maximum_coverage = 0.8
default_first_coverage = 0.1

# set example distribution dates and use rates at surveys
dist_dates = c(as.Date('2014-08-20'),as.Date('2017-05-01'),as.Date('2020-04-15'),as.Date('2023-12-20'))
matched_dhs_date = c(as.Date('2015-11-15'),as.Date('2018-11-20'),as.Date('2021-10-01'),as.Date('2024-07-15'))
matched_dhs_val = c(0.4, 0.35, 0.43, 0.55)
itn_distribution_info = data.frame('admin_name' = rep('Example', length(dist_dates)), 'State'=rep('Example', length(dist_dates)), date=dist_dates, matched_dhs_date=matched_dhs_date, matched_dhs_val=matched_dhs_val)



# assign distribution use-coverages
# if a distribution is associated with a dhs survey,
#   1) adjust for the seasonal use at the time of the survey
#   2) subtract the residual coverage from prior net distributions (assume that net distributions before the previous one have negligible contribution to net total)
#   3) back-calculate what coverage would have been from distribution to have given rise to observed coverage
# if the first distribution is not associated with a dhs survey, use a default coverage (default_first_coverage).
# if a later distribution is not associated with a dhs survey, use the coverage from the prior survey.
itn_distribution_info$coverage = NA
for(i_dist in 1:nrow(itn_distribution_info)){
  if(is.na(itn_distribution_info$matched_dhs_date[i_dist])){
    if(i_dist==1){
      itn_distribution_info$coverage[i_dist] = default_first_coverage
    } else{
      itn_distribution_info$coverage[i_dist] = itn_distribution_info$coverage[i_dist-1]
    }
  } else{
    if(i_dist>1){
      # time between the current dhs survey and the prior distribution (i_dist-1)
      time_lapse_days_residual = as.integer(itn_distribution_info$matched_dhs_date[i_dist] - itn_distribution_info$date[i_dist-1])
      # calculate what coverage remains from the previous distribution
      residual_coverage = itn_distribution_info$coverage[i_dist-1] * (1-plnorm(time_lapse_days_residual, meanlog=net_life_lognormal_mu, sdlog=net_life_lognormal_sigma))
    } else{
      residual_coverage = 0
    }
    season_adjusted_dhs_coverage = min(maximum_coverage, itn_distribution_info$matched_dhs_val[i_dist]  / seasonality_monthly_scalar[month(itn_distribution_info$matched_dhs_date[i_dist])])
    # remove nets leftover from prior distribution
    # NOTE: nets are distributed at random, not selectively to people who didn't have them before, so we don't just subtract old nets
    #     total_coverage = 1- ((1-coverage_from_residual_nets)*(1-coverage_from_new_nets))
    season_adjusted_dhs_coverage = max(0, 1-(1-season_adjusted_dhs_coverage)/(1-residual_coverage))
    
    # extrapolate back to what the distribution coverage would need to be to obtain this observed coverage
    # time between the dhs survey and the distribution
    time_lapse_days = as.integer(itn_distribution_info$matched_dhs_date[i_dist] - itn_distribution_info$date[i_dist])
    # fraction of orginally-disributed nets remaining at the time of the survey
    remaining_net_frac = 1-plnorm(time_lapse_days, meanlog=net_life_lognormal_mu, sdlog=net_life_lognormal_sigma)
    backward_projected_dist_coverage = season_adjusted_dhs_coverage / remaining_net_frac
    backward_projected_dist_coverage = min(maximum_coverage, backward_projected_dist_coverage)
    itn_distribution_info$coverage[i_dist] = backward_projected_dist_coverage   

  }
}

# plot timeseries of ITN coverage from mass distributions in each LGA given distribution schedule, seasonality in ITN use, and net decay parameters
# also include dots for DHS observations
# data format: a long dataframe with columns for LGA name, state name, archetype name, date, and coverage
coverage_df = itn_distribution_info
coverage_df$itn_u5 = coverage_df$coverage
all_admins = unique(coverage_df$admin_name)  
coverage_timeseries = data.frame('admin_name'=c(), 'State'=c(), 'date'=c(), 'coverage'=c())
first_day = min(dist_dates)
date_vector = seq.Date(first_day, as.Date('2025-01-01'), by='month')
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
  geom_vline(data=coverage_df, aes(xintercept=date), col='lightgrey', size=2) +
  geom_line(data=coverage_timeseries, aes(x=date, y=adjusted_coverage), col='purple', size=1) +
  geom_point(data=itn_distribution_info, aes(x=matched_dhs_date,y=matched_dhs_val), col='black', size=3.5) +
  ylab('estimated ITN use rate ') +
  # geom_point(data=itn_distribution_info, aes(x=matched_dhs_date,y=matched_dhs_val), shape=21, col='black') +
  # coord_cartesian(xlim=c(as.Date('2010-01-01'), as.Date('2025-01-01'))) +
  theme_classic()+
  theme(legend.position='none') +
  facet_wrap('State', nrow=5)
ggsave(filename=paste0(hbhi_dir, '/simulation_inputs/plots/EXAMPLE_itn_use_rate_timeseries_extrapolation_wSeasonAdjust.png'), plot=gg, width=18, height=15, units='in', dpi=900)
