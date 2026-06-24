# create_sim_input_cm.R

library(reshape2)
library(dplyr)


###############################
# case management
###############################
# from DHS output on cm rates in U5, create a csv file that will be read in by simulations. each row corresponds to an admin and specifies some period of implementation in that admin. included in the csv will be
#    - admin name
#    - seed number (determines which sampled parameter set is used)
#    - start year
#    - simday when these CM rates start
#    - duration of these CM rates 
#    - U5 coverage
#    - adult coverage
#    - severe coverage

# # function input examples:
# hbhi_dir = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/burundi_hbhi'
# cm_variable_name = 'cm'
# # assume that adult coverage is half that of U5
# adult_multiplier = 0.5
# # assume that severe coverage is twice that of U5
# severe_multiplier = 2
# # assume minimum severe coverage is 60%
# severe_minimum = 0.6
# # assume maximum coverage for any age/severity is 90%
# maximum_coverage = 0.9
# sim_start_year = 2010

create_cm_input_from_DHS = function(hbhi_dir, cm_variable_name='cm', sim_start_year=2010, adult_multiplier=0.5, severe_multiplier=2, severe_minimum=0.6, maximum_coverage=0.9, act_adherence_effective_multiplier=1, cm_suffix=''){
  # read in coverages for U5
  sample_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_',cm_variable_name,'.csv'))[,-1]
  # convert from wide to long format to make each seed its own row
  coverage_df = reshape2::melt(sample_df, id.vars=c("admin_name", "year"))
  colnames(coverage_df)[colnames(coverage_df)=='variable'] = 'seed'
  colnames(coverage_df)[colnames(coverage_df)=='value'] = 'U5_coverage'
  # replace "sample_" with "" in the seed column
  coverage_df$seed = gsub('sample_','',coverage_df$seed)
  if(length(unique(coverage_df$seed))<=1){
    coverage_df = coverage_df[,-which(colnames(coverage_df)=='seed')]
  }
  
  # add in coverage for adults
  coverage_df$adult_coverage = sapply((coverage_df$U5_coverage * adult_multiplier), min, maximum_coverage)
  # add severe disease coverage
  coverage_df$severe_coverage = sapply(sapply((coverage_df$U5_coverage * severe_multiplier), min, maximum_coverage), max, severe_minimum)
  
  # decrease the effective coverage to account for imperfect adherence to the full ACT regimine
  coverage_df$adult_coverage = coverage_df$adult_coverage * act_adherence_effective_multiplier
  coverage_df$U5_coverage = coverage_df$U5_coverage * act_adherence_effective_multiplier
  
  # add in the day of the simulation each intervention should start and the duration
  coverage_df$simday = (coverage_df$year - sim_start_year) * 365
  coverage_df$duration = 365  # duration in days

  # re-order rows for easier review
  coverage_df = coverage_df %>%
    arrange(admin_name, year)

  if(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent'))) dir.create(paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent'))
  write.csv(coverage_df, paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/cm_2010', cm_suffix,'_toPresent.csv'), row.names=FALSE)
}


create_season_calib_cm_input_from_DHS = function(hbhi_dir, cm_variable_name='cm', dhs_year_cm_burnin=2010, dhs_year_cm_calib=2016, adult_multiplier=0.5, severe_multiplier=2, severe_minimum=0.6, maximum_coverage=0.9){
  # use the coverage observed in the dhs_year_cm_burnin DHS survey for all burnin years and the coverage observed in dhs_year_cm_calib for all main calib years
  # read in coverages for U5
  season_arch_rates = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_archetype_rates.csv'))
  season_arch_rates = season_arch_rates[,which(colnames(season_arch_rates) %in% c('seasonality_archetype','year', paste0(cm_variable_name,'_rate')))]
  colnames(season_arch_rates)[colnames(season_arch_rates)==paste0(cm_variable_name,'_rate')] = 'U5_coverage'
  colnames(season_arch_rates)[colnames(season_arch_rates)=='seasonality_archetype'] = 'admin_name'
  
  
  # add in coverage for adults
  season_arch_rates$adult_coverage = sapply((season_arch_rates$U5_coverage * adult_multiplier), min, maximum_coverage)
  # add severe disease coverage
  season_arch_rates$severe_coverage = sapply(sapply((season_arch_rates$U5_coverage * severe_multiplier), min, maximum_coverage), max, severe_minimum)
  
  # add in the day of the simulation each intervention should start and the duration  (CM is constant for entirity of burnin and for calib, so start day=0 and duration =-1)
  #   only a single year's value will be retained for each archetype (set by dhs_year_cm_burnin and dhs_year_cm_calib)
  season_arch_rates$simday = 0
  season_arch_rates$duration = -1
  
  # extract U5 coverage for burnin and for main calibration
  burnin_coverage = season_arch_rates[season_arch_rates$year == dhs_year_cm_burnin,]
  calib_coverage = season_arch_rates[season_arch_rates$year == dhs_year_cm_calib,]
  
  if(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/interventions_calib'))) dir.create(paste0(hbhi_dir, '/simulation_inputs/interventions_calib'))
  write.csv(burnin_coverage, paste0(hbhi_dir, '/simulation_inputs/interventions_calib/cm_season_burnin.csv'))
  write.csv(calib_coverage, paste0(hbhi_dir, '/simulation_inputs/interventions_calib/cm_season_calib.csv'))
}


create_trans_calib_cm_input_from_DHS = function(hbhi_dir, cm_variable_name='cm', sim_start_year=2010, adult_multiplier=0.5, severe_multiplier=2, severe_minimum=0.6, maximum_coverage=0.9){
  # read in coverages for U5
  sample_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_quantiles_',cm_variable_name,'.csv'))[,-1]
  # convert from wide to long format to make each seed its own row
  coverage_df = reshape2::melt(sample_df, id.vars=c("admin_name", "year"))
  colnames(coverage_df)[colnames(coverage_df)=='variable'] = 'seed'
  colnames(coverage_df)[colnames(coverage_df)=='value'] = 'U5_coverage'
  # replace "sample_" with "" in the seed column
  coverage_df$seed = gsub('sample_','',coverage_df$seed)
  
  
  # add in coverage for adults
  coverage_df$adult_coverage = sapply((coverage_df$U5_coverage * adult_multiplier), min, maximum_coverage)
  # add severe disease coverage
  coverage_df$severe_coverage = sapply(sapply((coverage_df$U5_coverage * severe_multiplier), min, maximum_coverage), max, severe_minimum)
  
  # add in the day of the simulation each intervention should start and the duration
  coverage_df$simday = (coverage_df$year - sim_start_year) * 365
  coverage_df$duration = 365  # duration in days
  
  if(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/interventions_calib'))) dir.create(paste0(hbhi_dir, '/simulation_inputs/interventions_calib'))
  write.csv(coverage_df, paste0(hbhi_dir, '/simulation_inputs/interventions_calib/cm_trans_calib.csv'))
}







###############################
# IPTp
###############################
# from DHS output on cm rates in U5, create a csv file that will be read in by simulations. each row corresponds to an admin and specifies some period of implementation in that admin. included in the csv will be
#    - admin name
#    - seed number (determines which sampled parameter set is used)
#    - start year
#    - simday when these IPTp rates start
#    - duration of these IPTp rates 

# # function input examples:
# hbhi_dir = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/burundi_hbhi'
# iptp_variable_name = 'iptp'
# # assume maximum coverage for any age/severity is 90%
# maximum_coverage = 0.9
# sim_start_year = 2010

create_iptp_input_from_DHS = function(hbhi_dir, iptp_variable_name='iptp', sim_start_year=2010, maximum_coverage=0.9){
  # read in coverages for U5
  sample_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_', iptp_variable_name,'.csv'))[,-1]
  # convert from wide to long format to make each seed its own row
  coverage_df = reshape2::melt(sample_df, id.vars=c("admin_name", "year"))
  colnames(coverage_df)[colnames(coverage_df)=='variable'] = 'seed'
  colnames(coverage_df)[colnames(coverage_df)=='value'] = 'coverage'
  # replace "sample_" with "" in the seed column
  coverage_df$seed = gsub('sample_','',coverage_df$seed)
  if(length(unique(coverage_df$seed))<=1){
    coverage_df = coverage_df[,-which(colnames(coverage_df)=='seed')]
  }
  
  # add in the day of the simulation each intervention should start and the duration
  coverage_df$simday = (coverage_df$year - sim_start_year) * 365
  coverage_df$duration = 365  # duration in days
  
  if(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent'))) dir.create(paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent'))
  write.csv(coverage_df, paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/iptp_2010_toPresent.csv'))
}




