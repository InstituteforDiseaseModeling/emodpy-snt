# process_aggregate_sim_output_functions.R


library(data.table)
library(dplyr)
library(tidyverse)
library(lubridate)
# library(rgdal)
library(lubridate)


###########################################################################################################
# total simulation burden over specified interval, averaged across seeds but separated by district
###########################################################################################################


# I think the code below isn't used anymore... I'm not updating the direct versus direct+indirect mortality in it for now

# function that returns dataframe where each row is an admin and columns contain U5 and all-age total number of cases, total number of deaths, incidence, death rate, and average pfpr within the time period in each admin, averaged over all seeds
get_total_burden = function(sim_output_filepath, experiment_name, admin_pop, comparison_start_year, comparison_end_year, district_subset, cur_admins='all', overwrite_files=FALSE){
  output_filename = paste0(sim_output_filepath, '/', experiment_name, '/totalBurden_', comparison_start_year,'_', comparison_end_year,'_', district_subset,'.csv')
  if(file.exists(output_filename) & !overwrite_files){
    burden_means = read.csv(output_filename)
  } else{
    burden_df = fread(file = paste0(sim_output_filepath, '/', experiment_name, '/malariaBurden_withAdjustments.csv'))
    # subset to appropriate time period
    burden_df = burden_df[intersect(which(burden_df$year >= comparison_start_year), which(burden_df$year <= comparison_end_year)),]
    if(!(cur_admins[1] == 'all')){
      # subset to appropriate admins
      burden_df = burden_df[burden_df$admin_name %in% cur_admins,]
    }

    # PfPR
    pfpr_u5_means = burden_df[,c('admin_name', 'PfPR_U5')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()
    pfpr_all_means = burden_df[,c('admin_name', 'PfPR_MiP_adjusted')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()
    pfpr_means = merge(pfpr_u5_means, pfpr_all_means, by='admin_name')
    colnames(pfpr_means)[which(colnames(pfpr_means) == 'PfPR_U5')] = 'pfpr_u5'
    colnames(pfpr_means)[which(colnames(pfpr_means) == 'PfPR_MiP_adjusted')] = 'pfpr_all'


    # mortality
    # divide mortality by total population size in simulation; will later multiply by true population of DS
    burden_df$mortality_pp_u5 = (burden_df$total_mortality_U5_1*1 + burden_df$total_mortality_U5_2*1) / 2 / burden_df$Statistical_Population  # weighted average of mortality estimates (1/2 mort_1, 1/2 mort_2)
    burden_df$mortality_pp_all = (burden_df$total_mortality_1*1 + burden_df$total_mortality_2*1) / 2 / burden_df$Statistical_Population  # weighted average of mortality estimates (1/2 mort_1, 1/2 mort_2)
    mortality_sums = burden_df[,c('admin_name', 'Run_Number', 'mortality_pp_u5', 'mortality_pp_all')] %>% dplyr::group_by(admin_name, Run_Number) %>% dplyr::summarise_all(sum) %>% dplyr::ungroup()
    mortality_sums = mortality_sums[,c('admin_name', 'mortality_pp_u5', 'mortality_pp_all')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()
    # multiply by admin population to get total number of deaths
    mortality_sums = merge(mortality_sums, admin_pop, by='admin_name')
    mortality_sums$deaths_u5 = mortality_sums$mortality_pp_u5 * mortality_sums$pop_size
    mortality_sums$deaths_all = mortality_sums$mortality_pp_all * mortality_sums$pop_size

    # mortality rate  (number of deaths in a year in each admin divided by the U5 or all-age population size times 1000)
    # take sum of deaths within each year
    burden_df_annual =  burden_df %>% dplyr::group_by(admin_name, Run_Number, year) %>%
      dplyr::summarise(total_mortality_all_1 = sum(total_mortality_1),
                       total_mortality_all_2 = sum(total_mortality_2),
                       New_clinical_cases_all = sum(New_Clinical_Cases),
                       Pop_all = mean(Statistical_Population),
                       total_mortality_U5_1 = sum(total_mortality_U5_1),
                       total_mortality_U5_2 = sum(total_mortality_U5_2),
                       New_clinical_cases_U5 = sum(New_clinical_cases_U5),
                       Pop_U5 = mean(Pop_U5)
      ) %>%
      dplyr::ungroup()
    burden_df_annual$mortality_rate_u5 = (burden_df_annual$total_mortality_U5_1*1 + burden_df_annual$total_mortality_U5_2*1) / 2 / burden_df_annual$Pop_U5 * 1000  # weighted average of mortality estimates (1/2 mort_1, 1/2 mort_2)
    burden_df_annual$mortality_rate_all = (burden_df_annual$total_mortality_all_1*1 + burden_df_annual$total_mortality_all_2*1) / 2 / burden_df_annual$Pop_all * 1000  # weighted average of mortality estimates (1/2 mort_1, 1/2 mort_2)
    # remove any NA rows
    burden_df_annual = burden_df_annual[!is.na(burden_df_annual$mortality_rate_u5),]
    burden_df_annual = burden_df_annual[!is.na(burden_df_annual$mortality_rate_all),]
    # get average annual rate across all included years
    mortality_rate_means = burden_df_annual[,c('admin_name', 'mortality_rate_u5', 'mortality_rate_all')] %>% dplyr::group_by(admin_name) %>%
      dplyr::summarise(mortality_rate_u5 = mean(mortality_rate_u5),
                       mortality_rate_all = mean(mortality_rate_all)) %>%
      dplyr::ungroup()



    # clinical cases
    # divide number of clinical cases by total population size in simulation; will later multiply by true population of DS to get estimated number of clinical cases in that DS
    burden_df$cases_pp_u5 = burden_df$New_clinical_cases_U5 / burden_df$Statistical_Population
    burden_df$cases_pp_all = burden_df$New_Clinical_Cases / burden_df$Statistical_Population
    case_sums = burden_df[,c('admin_name', 'Run_Number', 'cases_pp_u5', 'cases_pp_all')] %>% dplyr::group_by(admin_name, Run_Number) %>% dplyr::summarise_all(sum) %>% dplyr::ungroup()
    case_sums = case_sums[,c('admin_name', 'cases_pp_u5', 'cases_pp_all')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()
    # multiply by admin population to get total number of cases
    case_sums = merge(case_sums, admin_pop, by='admin_name')
    case_sums$clinical_cases_u5 = case_sums$cases_pp_u5 * case_sums$pop_size
    case_sums$clinical_cases_all = case_sums$cases_pp_all * case_sums$pop_size

    # incidence (number of cases in a year in each admin divided by the U5 or all-age population size times 1000)
    burden_df_annual$incidence_u5 = burden_df_annual$New_clinical_cases_U5 / burden_df_annual$Pop_U5 * 1000
    burden_df_annual$incidence_all = burden_df_annual$New_clinical_cases_all / burden_df_annual$Pop_all * 1000
    # remove any NA rows
    burden_df_annual = burden_df_annual[!is.na(burden_df_annual$incidence_u5),]
    burden_df_annual = burden_df_annual[!is.na(burden_df_annual$incidence_all),]
    # get average annual rate across all included years
    incidence_means = burden_df_annual[,c('admin_name', 'incidence_u5', 'incidence_all')] %>% dplyr::group_by(admin_name) %>%
      dplyr::summarise(incidence_u5 = mean(incidence_u5),
                       incidence_all = mean(incidence_all)) %>%
      dplyr::ungroup()




    # fraction of population U5
    burden_df$frac_u5 = burden_df$Pop_U5 / burden_df$Statistical_Population
    frac_u5 = burden_df[,c('admin_name', 'frac_u5')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()

    burden_means = merge(pfpr_means, mortality_sums)#, by='admin_name')
    burden_means = merge(burden_means, case_sums)#, by=c('admin_name', 'pop_size'))
    burden_means = merge(burden_means, mortality_rate_means)#, by=c('admin_name'))
    burden_means = merge(burden_means, incidence_means)#, by=c('admin_name'))
    burden_means = merge(burden_means, frac_u5)#, by=c('admin_name'))
    burden_means$pop_size_u5 = burden_means$pop_size * burden_means$frac_u5

    write.csv(burden_means, output_filename, row.names=FALSE)
  }
  return(burden_means)
}




# 
# 
# 
# 
# 
# 
# 
# # function that returns dataframe where each row is an admin and columns contain U1 total number of cases, total number of deaths, incidence, death rate, and average pfpr within the time period in each admin, averaged over all seeds
# get_total_U1_burden = function(sim_output_filepath, experiment_name, admin_pop, comparison_start_year, comparison_end_year, district_subset, cur_admins='all', overwrite_files=FALSE){
#   output_filename = paste0(sim_output_filepath, '/', experiment_name, '/totalU1Burden_', comparison_start_year,'_', comparison_end_year,'_', district_subset,'.csv')
#   if(file.exists(output_filename) & !overwrite_files){
#     burden_means = read.csv(output_filename)
#   } else{
#     
#     burden_df = fread(file = paste0(sim_output_filepath, '/', experiment_name, '/malariaBurden_withAdjustments.csv'), check.names=TRUE)
#     # subset to appropriate time period
#     burden_df = burden_df[intersect(which(burden_df$year >= comparison_start_year), which(burden_df$year <= comparison_end_year)),]
#     if(!(cur_admins[1] == 'all')){
#       # subset to appropriate admins
#       burden_df = burden_df[burden_df$admin_name %in% cur_admins,] 
#     }
#     
#     # PfPR
#     pfpr_means = burden_df[,c('admin_name', 'PfPR_U1')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()
#     colnames(pfpr_means)[which(colnames(pfpr_means) == 'PfPR_U1')] = 'pfpr_u1'
#     
#     # number of deaths
#     # divide mortality by total population size in simulation; will later multiply by true population of DS
#     burden_df$mortality_pp_u1 = (burden_df$total_mortality_U1_1*1 + burden_df$total_mortality_U1_2*1) / 2 / burden_df$Statistical_Population  # weighted average of mortality estimates (1/2 mort_1, 1/2 mort_2)
#     mortality_sums = burden_df[,c('admin_name', 'Run_Number', 'mortality_pp_u1')] %>% dplyr::group_by(admin_name, Run_Number) %>% dplyr::summarise_all(sum) %>% dplyr::ungroup()
#     mortality_sums = mortality_sums[,c('admin_name', 'mortality_pp_u1')] %>% dplyr::group_by(admin_name) %>% dplyr::summarise_all(mean) %>% dplyr::ungroup()
#     # multiply by admin population to get total number of deaths
#     mortality_sums = merge(mortality_sums, admin_pop, by='admin_name')
#     mortality_sums$deaths_u1 = mortality_sums$mortality_pp_u1 * mortality_sums$pop_size
#     
#     # mortality rate  (number of deaths in a year in each admin divided by the U1 population size times 1000)
#     # take sum of deaths within each year
#     burden_df_annual =  burden_df %>% dplyr::group_by(admin_name, Run_Number, year) %>% 
#       dplyr::summarise(total_mortality_U1_1 = sum(total_mortality_U1_1),
#                        total_mortality_U1_2 = sum(total_mortality_U1_2),
#                        New_clinical_cases_U1 = sum(New_clinical_cases_U1),
#                        Pop_U1 = mean(Pop_U1)
#       ) %>% 
#       dplyr::ungroup()
#     burden_df_annual$mortality_rate_u1 = (burden_df_annual$total_mortality_U1_1*1 + burden_df_annual$total_mortality_U1_2*1) / 2 / burden_df_annual$Pop_U1 * 1000  # weighted average of mortality estimates (1/2 mort_1, 1/2 mort_2)
#     # remove any NA rows
#     burden_df_annual = burden_df_annual[!is.na(burden_df_annual$mortality_rate_u1),]
#     # get average annual rate across all included years
#     mortality_rate_means = burden_df_annual[,c('admin_name', 'mortality_rate_u1')] %>% dplyr::group_by(admin_name) %>% 
#       dplyr::summarise(mortality_rate_u1 = mean(mortality_rate_u1)) %>% 
#       dplyr::ungroup()
#     
#     
#     
#     # clinical cases
#     # divide number of clinical cases by total population size in simulation; will later multiply by true population of DS to get estimated number of clinical cases in that DS
#     burden_df$cases_pp_u1 = burden_df$New_clinical_cases_U1 / burden_df$Statistical_Population
#     case_sums = burden_df[,c('admin_name', 'Run_Number', 'cases_pp_u1')] %>% group_by(admin_name, Run_Number) %>% summarise_all(sum) %>% ungroup()
#     case_sums = case_sums[,c('admin_name', 'cases_pp_u1')] %>% group_by(admin_name) %>% summarise_all(mean) %>% ungroup()
#     # multiply by admin population to get total number of cases
#     case_sums = merge(case_sums, admin_pop, by='admin_name')
#     case_sums$clinical_cases_u1 = case_sums$cases_pp_u1 * case_sums$pop_size
#     
#     
#     
#     # incidence (number of cases in a year in each admin divided by the U1 population size times 1000)
#     burden_df_annual$incidence_u1 = burden_df_annual$New_clinical_cases_U1 / burden_df_annual$Pop_U1 * 1000
#     # remove any NA rows
#     burden_df_annual = burden_df_annual[!is.na(burden_df_annual$incidence_u1),]
#     # get average annual rate across all included years
#     incidence_means = burden_df_annual[,c('admin_name', 'incidence_u1')] %>% dplyr::group_by(admin_name) %>% 
#       dplyr::summarise(incidence_u1 = mean(incidence_u1)) %>% 
#       dplyr::ungroup()
#     
#     
#     # fraction of population U1
#     burden_df$frac_u1 = burden_df$Pop_U1 / burden_df$Statistical_Population
#     frac_u1 = burden_df[,c('admin_name', 'frac_u1')] %>% group_by(admin_name) %>% summarise_all(mean) %>% ungroup()
#     
#     burden_means = merge(pfpr_means, mortality_sums, by='admin_name')
#     burden_means = merge(burden_means, case_sums, by=c('admin_name', 'pop_size'))
#     burden_means = merge(burden_means, mortality_rate_means, by=c('admin_name'))
#     burden_means = merge(burden_means, incidence_means, by=c('admin_name'))
#     burden_means = merge(burden_means, frac_u1, by=c('admin_name'))
#     burden_means$pop_size_u1 = burden_means$pop_size * burden_means$frac_u1
#     
#     write.csv(burden_means, output_filename, row.names=FALSE)
#   }
#   
#   return(burden_means)
# }






###############################################################################################################################
# total simulation burden over specified interval, aggregated across all included districts and separated by seed
###############################################################################################################################

# total over a time period
get_cumulative_burden = function(sim_output_filepath, experiment_name, start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all', LLIN2y_flag=FALSE, overwrite_files=FALSE){
#'  @description get cumulative U5 and all-age burden over specified years in specified districts (for all malaria metrics, separate values for each seed)
#'  @return save and return data frame where each row is a seed and each column is the total over all included years of different burden metrics:
#'      sum of:
#'         - clinical cases (all ages)
#'         - clinical cases (U5)
#'         - deaths (all ages) - upper, lower, average parameter estimates
#'         - deaths (U5) - upper, lower, average parameter estimates
#'         - mLBW
#'         - malaria-attributable stillbirths
#'       average of annual (population weighted):
#'         - PfPR (all ages)
#'         - PfPR (U5)
#'         - incidence (all ages)
#'         - incidence (U5)
#'         - death rate (all ages)
#'         - death rate (U5)
#'         - mLBW
#'         - stillbirths


  output_filename = paste0(sim_output_filepath, '/', experiment_name, '/cumulativeBurden_', start_year, '_', end_year, '_', district_subset, '.csv')
  if(file.exists(output_filename) & !overwrite_files){
    df_aggregated = read.csv(output_filename)
  }else{
    # if we include all admins, get list of names from population size dataframe
    if(cur_admins[1] == 'all'){
      cur_admins = unique(admin_pop$admin_name)
    }
    
    cur_file = fread(file = paste0(sim_output_filepath, '/', experiment_name, '/malariaBurden_withAdjustments.csv'), check.names=TRUE)
    # filter to relevant years
    cur_file = cur_file[cur_file$year <= end_year,]
    cur_file = cur_file[cur_file$year >= start_year,]
    # merge population sizes in each admin
    df = merge(cur_file, admin_pop, by='admin_name')
    # subset to appropriate admins
    df = df[df$admin_name %in% cur_admins,]  
    
    # all age metrics - rescaled to full population
    df$positives_all_ages = df$PfPR_MiP_adjusted * df$pop_size
    df$cases_all_ages = df$New_Clinical_Cases * (df$pop_size / df$Statistical_Population)
    df$direct_deaths_1_all_ages = df$direct_mortality_nonMiP_1 * df$pop_size / df$Statistical_Population
    df$direct_deaths_2_all_ages = df$direct_mortality_nonMiP_2 * df$pop_size / df$Statistical_Population
    df$all_deaths_1_all_ages = df$total_mortality_1 * df$pop_size / df$Statistical_Population
    df$all_deaths_2_all_ages = df$total_mortality_2 * df$pop_size / df$Statistical_Population
    df$num_mLBW = df$mLBW_births * df$pop_size / df$Statistical_Population
    df$num_mStillbirths = df$MiP_stillbirths * df$pop_size / df$Statistical_Population
    # U5 metrics - rescaled to full population
    df$pop_size_U5 = df$pop_size * (df$Pop_U5 / df$Statistical_Population)  # assumes fraction of individual U5 in simulation is same as fraction in full population
    df$positives_U5 = df$PfPR_U5 * df$pop_size_U5
    df$cases_U5 = df$New_clinical_cases_U5 * df$pop_size_U5 / df$Pop_U5
    df$direct_deaths_1_U5 = df$direct_mortality_nonMiP_U5_1 * df$pop_size_U5 / df$Pop_U5
    df$direct_deaths_2_U5 = df$direct_mortality_nonMiP_U5_2 * df$pop_size_U5 / df$Pop_U5
    df$all_deaths_1_U5 = df$total_mortality_U5_1 * df$pop_size_U5 / df$Pop_U5
    df$all_deaths_2_U5 = df$total_mortality_U5_2 * df$pop_size_U5 / df$Pop_U5
    
    df_aggregated = df %>% group_by(Run_Number) %>%
      dplyr::summarize(pop_all_sum = sum(pop_size),
                       pop_U5_sum = sum(pop_size_U5),
                       cases_all_sum = sum(cases_all_ages),
                       cases_U5_sum = sum(cases_U5),
                       positives_all_sum = sum(positives_all_ages),
                       positives_U5_sum = sum(positives_U5),
                       direct_deaths_1_all_sum = sum(direct_deaths_1_all_ages),
                       direct_deaths_2_all_sum = sum(direct_deaths_2_all_ages),
                       all_deaths_1_all_sum = sum(all_deaths_1_all_ages),
                       all_deaths_2_all_sum = sum(all_deaths_2_all_ages),
                       direct_deaths_1_U5_sum = sum(direct_deaths_1_U5),
                       direct_deaths_2_U5_sum = sum(direct_deaths_2_U5),
                       all_deaths_1_U5_sum = sum(all_deaths_1_U5),
                       all_deaths_2_U5_sum = sum(all_deaths_2_U5),
                       mLBW_sum = sum(num_mLBW),
                       mStill_sum = sum(num_mStillbirths),
                       num_values_grouped = n())
    
    
    # clinical incidence (annual): sum of number of cases over all months / (pop size) / (number of years) * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months /  number of months) / (number of months / 12)  * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months / 12)  * 1000
    df_aggregated$incidence_all = (df_aggregated$cases_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$incidence_U5 = (df_aggregated$cases_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_1_all = (df_aggregated$direct_deaths_1_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_all = (df_aggregated$direct_deaths_2_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_all = (df_aggregated$direct_death_rate_1_all + df_aggregated$direct_death_rate_2_all) / 2
    df_aggregated$all_death_rate_1_all = (df_aggregated$all_deaths_1_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_all = (df_aggregated$all_deaths_2_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_all = (df_aggregated$all_death_rate_1_all + df_aggregated$all_death_rate_2_all) / 2
    df_aggregated$direct_death_rate_1_U5 = (df_aggregated$direct_deaths_1_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_U5 = (df_aggregated$direct_deaths_2_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_U5 = (df_aggregated$direct_death_rate_1_U5 + df_aggregated$direct_death_rate_2_U5) / 2
    df_aggregated$all_death_rate_1_U5 = (df_aggregated$all_deaths_1_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_U5 = (df_aggregated$all_deaths_2_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_U5 = (df_aggregated$all_death_rate_1_U5 + df_aggregated$all_death_rate_2_U5) / 2
    df_aggregated$average_PfPR_all = (df_aggregated$positives_all_sum  / (df_aggregated$pop_all_sum))
    df_aggregated$average_PfPR_U5 = (df_aggregated$positives_U5_sum  / (df_aggregated$pop_U5_sum))
    df_aggregated$annual_num_mLBW = (df_aggregated$mLBW_sum / (end_year - start_year + 1))
    df_aggregated$annual_num_mStill = (df_aggregated$mStill_sum / (end_year - start_year + 1))
    
    
    df_aggregated = df_aggregated[,which(colnames(df_aggregated) %in% c('Run_Number','cases_all_sum','cases_U5_sum', 
                                                                        'mLBW_sum', 'mStill_sum', 'incidence_all', 'incidence_U5', 'average_PfPR_all','average_PfPR_U5',  'annual_num_mLBW', 'annual_num_mStill', 
                                                                        'direct_deaths_1_all_sum', 'direct_deaths_2_all_sum', 'direct_deaths_1_U5_sum', 'direct_deaths_2_U5_sum',
                                                                        'all_deaths_1_all_sum', 'all_deaths_2_all_sum', 'all_deaths_1_U5_sum', 'all_deaths_2_U5_sum', 
                                                                        'direct_death_rate_1_all', 'direct_death_rate_2_all',  'all_death_rate_1_all', 'all_death_rate_2_all', 
                                                                        'direct_death_rate_1_U5', 'direct_death_rate_2_U5', 'all_death_rate_1_U5', 'all_death_rate_2_U5', 
                                                                        'direct_death_rate_mean_all', 'direct_death_rate_mean_U5', 'all_death_rate_mean_all', 'all_death_rate_mean_U5'))]
    write.csv(df_aggregated, output_filename, row.names=FALSE)
  }
  return(df_aggregated)
}



# total over a time period for each state
get_cumulative_burden_by_state = function(sim_output_filepath, experiment_name, start_year, end_year, admin_pop, group_col='State', LLIN2y_flag=FALSE, overwrite_files=FALSE, mean_across_seeds=FALSE){
  #'  @param group_col column in admin_pop used to group LGAs (default 'State'). Set to 'Funder' (or any
  #'        admin-level column present in admin_pop) to aggregate by that grouping instead. Population
  #'        weighting is unaffected: rates are Sum(rate*pop)/Sum(pop) over the LGAs within each group.
  #'  @description get cumulative U5 and all-age burden over specified years in each state (for all malaria metrics, separate values for each seed)
  #'  @return save and return data frame where each row is a seed x state and each column is the total over all included years of different burden metrics:
  #'      sum of:
  #'         - clinical cases (all ages)
  #'         - clinical cases (U5)
  #'         - deaths (all ages) - upper, lower, average parameter estimates
  #'         - deaths (U5) - upper, lower, average parameter estimates
  #'         - mLBW
  #'         - malaria-attributable stillbirths
  #'       average of annual (population weighted):
  #'         - PfPR (all ages)
  #'         - PfPR (U5)
  #'         - incidence (all ages)
  #'         - incidence (U5)
  #'         - death rate (all ages)
  #'         - death rate (U5)
  #'         - mLBW
  #'         - stillbirths
  
  
  output_filename = paste0(sim_output_filepath, '/', experiment_name, '/cumulativeBurden_', start_year, '_', end_year, '_by', group_col, '.csv')
  if(file.exists(output_filename) & !overwrite_files){
    df_aggregated = read.csv(output_filename)
  }else{
    cur_file = fread(file = paste0(sim_output_filepath, '/', experiment_name, '/malariaBurden_withAdjustments.csv'), check.names=TRUE)
    # filter to relevant years
    cur_file = cur_file[cur_file$year <= end_year,]
    cur_file = cur_file[cur_file$year >= start_year,]
    # merge population sizes in each admin
    df = merge(cur_file, admin_pop, by='admin_name')
    
    # all age metrics - rescaled to full population
    df$positives_all_ages = df$PfPR_MiP_adjusted * df$pop_size
    df$cases_all_ages = df$New_Clinical_Cases * (df$pop_size / df$Statistical_Population)
    df$severe_all_ages = df$severe_total * (df$pop_size / df$Statistical_Population)
    df$direct_deaths_1_all_ages = df$direct_mortality_nonMiP_1 * df$pop_size / df$Statistical_Population
    df$direct_deaths_2_all_ages = df$direct_mortality_nonMiP_2 * df$pop_size / df$Statistical_Population
    df$all_deaths_1_all_ages = df$total_mortality_1 * df$pop_size / df$Statistical_Population
    df$all_deaths_2_all_ages = df$total_mortality_2 * df$pop_size / df$Statistical_Population
    df$num_mLBW = df$mLBW_births * df$pop_size / df$Statistical_Population
    df$num_mStillbirths = df$MiP_stillbirths * df$pop_size / df$Statistical_Population
    # U5 metrics - rescaled to full population
    df$pop_size_U5 = df$pop_size * (df$Pop_U5 / df$Statistical_Population)  # assumes fraction of individual U5 in simulation is same as fraction in full population
    df$positives_U5 = df$PfPR_U5 * df$pop_size_U5
    df$cases_U5 = df$New_clinical_cases_U5 * df$pop_size_U5 / df$Pop_U5
    df$severe_U5 = df$Severe_cases_U5 * df$pop_size_U5 / df$Pop_U5
    df$direct_deaths_1_U5 = df$direct_mortality_nonMiP_U5_1 * df$pop_size_U5 / df$Pop_U5
    df$direct_deaths_2_U5 = df$direct_mortality_nonMiP_U5_2 * df$pop_size_U5 / df$Pop_U5
    df$all_deaths_1_U5 = df$total_mortality_U5_1 * df$pop_size_U5 / df$Pop_U5
    df$all_deaths_2_U5 = df$total_mortality_U5_2 * df$pop_size_U5 / df$Pop_U5
    
    df_aggregated = df %>% group_by(Run_Number, across(all_of(group_col))) %>%
      dplyr::summarize(pop_all_sum = sum(pop_size),
                       pop_U5_sum = sum(pop_size_U5),
                       cases_all_sum = sum(cases_all_ages),
                       cases_U5_sum = sum(cases_U5),
                       severe_all_sum = sum(severe_all_ages),
                       severe_U5_sum = sum(severe_U5),
                       positives_all_sum = sum(positives_all_ages),
                       positives_U5_sum = sum(positives_U5),
                       direct_deaths_1_all_sum = sum(direct_deaths_1_all_ages),
                       direct_deaths_2_all_sum = sum(direct_deaths_2_all_ages),
                       all_deaths_1_all_sum = sum(all_deaths_1_all_ages),
                       all_deaths_2_all_sum = sum(all_deaths_2_all_ages),
                       direct_deaths_1_U5_sum = sum(direct_deaths_1_U5),
                       direct_deaths_2_U5_sum = sum(direct_deaths_2_U5),
                       all_deaths_1_U5_sum = sum(all_deaths_1_U5),
                       all_deaths_2_U5_sum = sum(all_deaths_2_U5),
                       mLBW_sum = sum(num_mLBW),
                       mStill_sum = sum(num_mStillbirths),
                       num_values_grouped = n())
    
    
    # clinical incidence (annual): sum of number of cases over all months / (pop size) / (number of years) * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months /  number of months) / (number of months / 12)  * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months / 12)  * 1000
    df_aggregated$incidence_all = (df_aggregated$cases_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$incidence_U5 = (df_aggregated$cases_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$severe_incidence_all = (df_aggregated$severe_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$severe_incidence_U5 = (df_aggregated$severe_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_1_all = (df_aggregated$direct_deaths_1_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_all = (df_aggregated$direct_deaths_2_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_all = (df_aggregated$direct_death_rate_1_all + df_aggregated$direct_death_rate_2_all) / 2
    df_aggregated$all_death_rate_1_all = (df_aggregated$all_deaths_1_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_all = (df_aggregated$all_deaths_2_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_all = (df_aggregated$all_death_rate_1_all + df_aggregated$all_death_rate_2_all) / 2
    df_aggregated$direct_death_rate_1_U5 = (df_aggregated$direct_deaths_1_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_U5 = (df_aggregated$direct_deaths_2_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_U5 = (df_aggregated$direct_death_rate_1_U5 + df_aggregated$direct_death_rate_2_U5) / 2
    df_aggregated$all_death_rate_1_U5 = (df_aggregated$all_deaths_1_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_U5 = (df_aggregated$all_deaths_2_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_U5 = (df_aggregated$all_death_rate_1_U5 + df_aggregated$all_death_rate_2_U5) / 2
    df_aggregated$average_PfPR_all = (df_aggregated$positives_all_sum  / (df_aggregated$pop_all_sum))
    df_aggregated$average_PfPR_U5 = (df_aggregated$positives_U5_sum  / (df_aggregated$pop_U5_sum))
    df_aggregated$annual_num_mLBW = (df_aggregated$mLBW_sum / (end_year - start_year + 1))
    df_aggregated$annual_num_mStill = (df_aggregated$mStill_sum / (end_year - start_year + 1))
    
    
    df_aggregated = df_aggregated[,which(colnames(df_aggregated) %in% c(group_col,'Run_Number','cases_all_sum','cases_U5_sum','severe_all_sum','severe_U5_sum',
                                                                        'mLBW_sum', 'mStill_sum', 'incidence_all', 'incidence_U5', 'severe_incidence_all', 'severe_incidence_U5', 'average_PfPR_all','average_PfPR_U5',  'annual_num_mLBW', 'annual_num_mStill', 
                                                                        'direct_deaths_1_all_sum', 'direct_deaths_2_all_sum', 'direct_deaths_1_U5_sum', 'direct_deaths_2_U5_sum',
                                                                        'all_deaths_1_all_sum', 'all_deaths_2_all_sum', 'all_deaths_1_U5_sum', 'all_deaths_2_U5_sum', 
                                                                        'direct_death_rate_1_all', 'direct_death_rate_2_all',  'all_death_rate_1_all', 'all_death_rate_2_all', 
                                                                        'direct_death_rate_1_U5', 'direct_death_rate_2_U5', 'all_death_rate_1_U5', 'all_death_rate_2_U5', 
                                                                        'direct_death_rate_mean_all', 'direct_death_rate_mean_U5', 'all_death_rate_mean_all', 'all_death_rate_mean_U5'))]
    
    write.csv(df_aggregated, output_filename, row.names=FALSE)
  }
  if(mean_across_seeds){
    df_aggregated = df_aggregated %>% group_by(across(all_of(group_col))) %>%
      summarise_all(mean, na.rm=TRUE)
  }
  return(df_aggregated)
}






# total over a time period for each admin
get_cumulative_burden_by_admin = function(sim_output_filepath, experiment_name, start_year, end_year, admin_pop, LLIN2y_flag=FALSE, overwrite_files=FALSE, mean_across_seeds=FALSE){
  #'  @description get cumulative U5 and all-age burden over specified years in each state (for all malaria metrics, separate values for each seed)
  #'  @return save and return data frame where each row is a seed x state and each column is the total over all included years of different burden metrics:
  #'      sum of:
  #'         - clinical cases (all ages)
  #'         - clinical cases (U5)
  #'         - deaths (all ages) - upper, lower, average parameter estimates
  #'         - deaths (U5) - upper, lower, average parameter estimates
  #'         - mLBW
  #'         - malaria-attributable stillbirths
  #'       average of annual (population weighted):
  #'         - PfPR (all ages)
  #'         - PfPR (U5)
  #'         - incidence (all ages)
  #'         - incidence (U5)
  #'         - death rate (all ages)
  #'         - death rate (U5)
  #'         - mLBW
  #'         - stillbirths
  
  
  output_filename = paste0(sim_output_filepath, '/', experiment_name, '/cumulativeBurden_', start_year, '_', end_year, '_byAdmin.csv')
  if(file.exists(output_filename) & !overwrite_files){
    df_aggregated = read.csv(output_filename)
  }else{
    cur_file = fread(file = paste0(sim_output_filepath, '/', experiment_name, '/malariaBurden_withAdjustments.csv'), check.names=TRUE)
    # filter to relevant years
    cur_file = cur_file[cur_file$year <= end_year,]
    cur_file = cur_file[cur_file$year >= start_year,]
    # merge population sizes in each admin
    df = merge(cur_file, admin_pop, by='admin_name')
    
    # all age metrics - rescaled to full population
    df$positives_all_ages = df$PfPR_MiP_adjusted * df$pop_size
    df$cases_all_ages = df$New_Clinical_Cases * (df$pop_size / df$Statistical_Population)
    df$severe_all_ages = df$severe_total * (df$pop_size / df$Statistical_Population)
    df$direct_deaths_1_all_ages = df$direct_mortality_nonMiP_1 * df$pop_size / df$Statistical_Population
    df$direct_deaths_2_all_ages = df$direct_mortality_nonMiP_2 * df$pop_size / df$Statistical_Population
    df$all_deaths_1_all_ages = df$total_mortality_1 * df$pop_size / df$Statistical_Population
    df$all_deaths_2_all_ages = df$total_mortality_2 * df$pop_size / df$Statistical_Population
    df$num_mLBW = df$mLBW_births * df$pop_size / df$Statistical_Population
    df$num_mStillbirths = df$MiP_stillbirths * df$pop_size / df$Statistical_Population
    # U5 metrics - rescaled to full population
    df$pop_size_U5 = df$pop_size * (df$Pop_U5 / df$Statistical_Population)  # assumes fraction of individual U5 in simulation is same as fraction in full population
    df$positives_U5 = df$PfPR_U5 * df$pop_size_U5
    df$cases_U5 = df$New_clinical_cases_U5 * df$pop_size_U5 / df$Pop_U5
    df$severe_U5 = df$Severe_cases_U5 * df$pop_size_U5 / df$Pop_U5
    df$direct_deaths_1_U5 = df$direct_mortality_nonMiP_U5_1 * df$pop_size_U5 / df$Pop_U5
    df$direct_deaths_2_U5 = df$direct_mortality_nonMiP_U5_2 * df$pop_size_U5 / df$Pop_U5
    df$all_deaths_1_U5 = df$total_mortality_U5_1 * df$pop_size_U5 / df$Pop_U5
    df$all_deaths_2_U5 = df$total_mortality_U5_2 * df$pop_size_U5 / df$Pop_U5
    
    df_aggregated = df %>% group_by(Run_Number, admin_name) %>%
      dplyr::summarize(pop_all_sum = sum(pop_size),
                       pop_U5_sum = sum(pop_size_U5),
                       cases_all_sum = sum(cases_all_ages),
                       cases_U5_sum = sum(cases_U5),
                       severe_all_sum = sum(severe_all_ages),
                       severe_U5_sum = sum(severe_U5),
                       positives_all_sum = sum(positives_all_ages),
                       positives_U5_sum = sum(positives_U5),
                       direct_deaths_1_all_sum = sum(direct_deaths_1_all_ages),
                       direct_deaths_2_all_sum = sum(direct_deaths_2_all_ages),
                       all_deaths_1_all_sum = sum(all_deaths_1_all_ages),
                       all_deaths_2_all_sum = sum(all_deaths_2_all_ages),
                       direct_deaths_1_U5_sum = sum(direct_deaths_1_U5),
                       direct_deaths_2_U5_sum = sum(direct_deaths_2_U5),
                       all_deaths_1_U5_sum = sum(all_deaths_1_U5),
                       all_deaths_2_U5_sum = sum(all_deaths_2_U5),
                       mLBW_sum = sum(num_mLBW),
                       mStill_sum = sum(num_mStillbirths),
                       num_values_grouped = n())
    
    
    # clinical incidence (annual): sum of number of cases over all months / (pop size) / (number of years) * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months /  number of months) / (number of months / 12)  * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months / 12)  * 1000
    df_aggregated$incidence_all = (df_aggregated$cases_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$incidence_U5 = (df_aggregated$cases_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$severe_incidence_all = (df_aggregated$severe_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$severe_incidence_U5 = (df_aggregated$severe_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_1_all = (df_aggregated$direct_deaths_1_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_all = (df_aggregated$direct_deaths_2_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_all = (df_aggregated$direct_death_rate_1_all + df_aggregated$direct_death_rate_2_all) / 2
    df_aggregated$all_death_rate_1_all = (df_aggregated$all_deaths_1_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_all = (df_aggregated$all_deaths_2_all_sum / (df_aggregated$pop_all_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_all = (df_aggregated$all_death_rate_1_all + df_aggregated$all_death_rate_2_all) / 2
    df_aggregated$direct_death_rate_1_U5 = (df_aggregated$direct_deaths_1_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_U5 = (df_aggregated$direct_deaths_2_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_U5 = (df_aggregated$direct_death_rate_1_U5 + df_aggregated$direct_death_rate_2_U5) / 2
    df_aggregated$all_death_rate_1_U5 = (df_aggregated$all_deaths_1_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_U5 = (df_aggregated$all_deaths_2_U5_sum / (df_aggregated$pop_U5_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_U5 = (df_aggregated$all_death_rate_1_U5 + df_aggregated$all_death_rate_2_U5) / 2
    df_aggregated$average_PfPR_all = (df_aggregated$positives_all_sum  / (df_aggregated$pop_all_sum))
    df_aggregated$average_PfPR_U5 = (df_aggregated$positives_U5_sum  / (df_aggregated$pop_U5_sum))
    df_aggregated$annual_num_mLBW = (df_aggregated$mLBW_sum / (end_year - start_year + 1))
    df_aggregated$annual_num_mStill = (df_aggregated$mStill_sum / (end_year - start_year + 1))
    
    
    df_aggregated = df_aggregated[,which(colnames(df_aggregated) %in% c('admin_name','Run_Number','cases_all_sum','cases_U5_sum','severe_all_sum','severe_U5_sum', 
                                                                        'mLBW_sum', 'mStill_sum', 'incidence_all', 'incidence_U5', 'severe_incidence_all', 'severe_incidence_U5', 'average_PfPR_all','average_PfPR_U5',  'annual_num_mLBW', 'annual_num_mStill', 
                                                                        'direct_deaths_1_all_sum', 'direct_deaths_2_all_sum', 'direct_deaths_1_U5_sum', 'direct_deaths_2_U5_sum',
                                                                        'all_deaths_1_all_sum', 'all_deaths_2_all_sum', 'all_deaths_1_U5_sum', 'all_deaths_2_U5_sum', 
                                                                        'direct_death_rate_1_all', 'direct_death_rate_2_all',  'all_death_rate_1_all', 'all_death_rate_2_all', 
                                                                        'direct_death_rate_1_U5', 'direct_death_rate_2_U5', 'all_death_rate_1_U5', 'all_death_rate_2_U5', 
                                                                        'direct_death_rate_mean_all', 'direct_death_rate_mean_U5', 'all_death_rate_mean_all', 'all_death_rate_mean_U5'))]
    
    if(mean_across_seeds){
      df_aggregated = df_aggregated %>% group_by(admin_name) %>%
        summarise_all(mean, na.rm=TRUE)
    }
    write.csv(df_aggregated, output_filename, row.names=FALSE)
  }
  if(mean_across_seeds){
    df_aggregated = df_aggregated %>% group_by(admin_name) %>%
      summarise_all(mean, na.rm=TRUE)
  }
  return(df_aggregated)
}









get_cumulative_U1_burden = function(sim_output_filepath, experiment_name, start_year, end_year, admin_pop, district_subset=district_subset, cur_admins='all', LLIN2y_flag=FALSE, overwrite_files=FALSE){
  #'  @description get cumulative U1 burden over specified years in specified districts (for all malaria metrics, separate values for each seed)
  #'  @return save and return data frame where each row is a seed and each column is the total over all included years of different burden metrics:
  #'      sum of:
  #'         - clinical cases (U1)
  #'         - deaths (U1) - upper, lower, average parameter estimates
  #'       average of annual (population weighted):
  #'         - PfPR (U1)
  #'         - incidence (U1)
  #'         - death rate (U1)

  

  output_filename = paste0(sim_output_filepath, '/', experiment_name, '/cumulativeBurden_IPTi_', start_year, '_', end_year, '_', district_subset, '.csv')
  if(file.exists(output_filename) & !overwrite_files){
    df_aggregated = read.csv(output_filename)
  }else{
    # if we include all admins, get list of names from population size dataframe
    if(cur_admins[1] == 'all'){
      cur_admins = unique(admin_pop$admin_name)
    }
    
    cur_file = fread(file = paste0(sim_output_filepath, '/', experiment_name, '/malariaBurden_withAdjustments.csv'), check.names=TRUE)
    # filter to relevant years
    cur_file = cur_file[cur_file$year <= end_year,]
    cur_file = cur_file[cur_file$year >= start_year,]
    # merge population sizes in each admin
    df = merge(cur_file, admin_pop, by='admin_name')
    # subset to appropriate admins
    df = df[df$admin_name %in% cur_admins,]  
    
    # U1 metrics - rescaled to full population
    df$pop_size_U1 = df$pop_size * (df$Pop_U1 / df$Statistical_Population)  # assumes fraction of individual U1 in simulation is same as fraction in full population
    df$positives_U1 = df$PfPR_U1 * df$pop_size_U1
    df$cases_U1 = df$New_clinical_cases_U1 * df$pop_size_U1 / df$Pop_U1
    df$direct_deaths_1_U1 = df$direct_mortality_nonMiP_U1_1 * df$pop_size_U1 / df$Pop_U1
    df$direct_deaths_2_U1 = df$direct_mortality_nonMiP_U1_2 * df$pop_size_U1 / df$Pop_U1
    df$all_deaths_1_U1 = df$total_mortality_U1_1 * df$pop_size_U1 / df$Pop_U1
    df$all_deaths_2_U1 = df$total_mortality_U1_2 * df$pop_size_U1 / df$Pop_U1
    
    df_aggregated = df %>% dplyr::group_by(Run_Number) %>%
      dplyr::summarize(pop_U1_sum = sum(pop_size_U1),
                       cases_U1_sum = sum(cases_U1),
                       positives_U1_sum = sum(positives_U1),
                       direct_deaths_1_U1_sum = sum(direct_deaths_1_U1),
                       direct_deaths_2_U1_sum = sum(direct_deaths_2_U1),
                       all_deaths_1_U1_sum = sum(all_deaths_1_U1),
                       all_deaths_2_U1_sum = sum(all_deaths_2_U1),
                       num_values_grouped = n())
    
    
    # clinical incidence (annual): sum of number of cases over all months / (pop size) / (number of years) * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months /  number of months) / (number of months / 12)  * 1000
    #     = sum of number of cases over all months / (sum of pop size over all months / 12)  * 1000
    df_aggregated$incidence_U1 = (df_aggregated$cases_U1_sum / (df_aggregated$pop_U1_sum / 12) * 1000)
    df_aggregated$direct_death_rate_1_U1 = (df_aggregated$direct_deaths_1_U1_sum / (df_aggregated$pop_U1_sum / 12) * 1000)
    df_aggregated$direct_death_rate_2_U1 = (df_aggregated$direct_deaths_2_U1_sum / (df_aggregated$pop_U1_sum / 12) * 1000)
    df_aggregated$direct_death_rate_mean_U1 = (df_aggregated$direct_death_rate_1_U1 + df_aggregated$direct_death_rate_2_U1) / 2
    df_aggregated$all_death_rate_1_U1 = (df_aggregated$all_deaths_1_U1_sum / (df_aggregated$pop_U1_sum / 12) * 1000)
    df_aggregated$all_death_rate_2_U1 = (df_aggregated$all_deaths_2_U1_sum / (df_aggregated$pop_U1_sum / 12) * 1000)
    df_aggregated$all_death_rate_mean_U1 = (df_aggregated$all_death_rate_1_U1 + df_aggregated$all_death_rate_2_U1) / 2
    df_aggregated$average_PfPR_U1 = (df_aggregated$positives_U1_sum  / (df_aggregated$pop_U1_sum))
    
    df_aggregated = df_aggregated[,which(colnames(df_aggregated) %in% c('Run_Number','cases_U1_sum', 'incidence_U1', 'average_PfPR_U1', 
                                                                        'direct_deaths_1_U1_sum', 'direct_deaths_2_U1_sum','direct_death_rate_1_U1', 'direct_death_rate_2_U1', 'direct_death_rate_mean_U1',
                                                                        'all_deaths_1_U1_sum', 'all_deaths_2_U1_sum','all_death_rate_1_U1', 'all_death_rate_2_U1', 'all_death_rate_mean_U1'))]
    write.csv(df_aggregated, output_filename, row.names=FALSE)
  }
  return(df_aggregated)
}





####################################################################################
# relative simulation burden between two experiments over specified time interval
####################################################################################

get_relative_burden = function(sim_output_filepath, reference_experiment_name, comparison_experiment_name, comparison_scenario_name, start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all', 
                               LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE ){
  #'  @description get relative change in U5 and all-age burden when comparing between two simulations in specified years and in specified districts (for all malaria metrics, separate values for each seed)
  #'  @return data frame where each row is a seed and each column is the relative change of different burden metrics, calculated as (reference-comparison) / reference:

  reference_df = get_cumulative_burden(sim_output_filepath=sim_output_filepath, experiment_name=reference_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_burden(sim_output_filepath=sim_output_filepath, experiment_name=comparison_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  
  if(align_seeds){  # compare one run seed against the matching run seed in the other experiment
    # align seeds
    reference_df = reference_df[order(reference_df$Run_Number),]
    comparison_df = comparison_df[order(comparison_df$Run_Number),]
  } else{  # compare the averages across all seeds from one experiment against the average from the other experiment
    reference_df = reference_df %>% summarise_all(mean)
    comparison_df = comparison_df %>% summarise_all(mean)
  }
  relative_burden_df = data.frame('Run_Number' = reference_df$Run_Number)
  # iterate through burden indicators, calculating relative burden and adding to dataframe
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) == 'Run_Number')]
  for(bb in 1:length(burden_indicators)){
    # relative_burden_cur = (comparison_df[[burden_indicators[bb]]] - reference_df[[burden_indicators[bb]]]) / reference_df[[burden_indicators[bb]]]
    relative_burden_cur = (reference_df[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]]) / reference_df[[burden_indicators[bb]]]
    relative_burden_df[[burden_indicators[bb]]] = relative_burden_cur
  }
  relative_burden_df$scenario = comparison_scenario_name
  return(relative_burden_df)
}



# get relative burden by state for grid plot
get_relative_burden_by_state = function(sim_output_filepath, reference_experiment_name, comparison_experiment_name, comparison_scenario_name, start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all', group_col='State',
                               LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE ){
  #'  @description get relative change in U5 and all-age burden when comparing between two simulations in specified years and in specified districts (for all malaria metrics, separate values for each seed)
  #'  @param group_col column in admin_pop used to group LGAs (default 'State'; e.g. 'Funder' to group by funder)
  #'  @return data frame where each row is a seed and each column is the relative change of different burden metrics, calculated as (reference-comparison) / reference:

  reference_df = get_cumulative_burden_by_state(sim_output_filepath=sim_output_filepath, experiment_name=reference_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, group_col=group_col, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_burden_by_state(sim_output_filepath=sim_output_filepath, experiment_name=comparison_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, group_col=group_col, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)

  if(align_seeds){  # compare one run seed against the matching run seed in the other experiment
    # align seeds and groups into same order
    reference_df = reference_df[order(reference_df[[group_col]], reference_df$Run_Number),]
    comparison_df = comparison_df[order(comparison_df[[group_col]], comparison_df$Run_Number),]
  } else{  # compare the averages across all seeds from one experiment against the average from the other experiment
    reference_df = reference_df %>% group_by(across(all_of(group_col))) %>% summarise_all(mean)
    comparison_df = comparison_df %>% group_by(across(all_of(group_col))) %>% summarise_all(mean)
    # align groups into same order
    reference_df = reference_df[order(reference_df[[group_col]]),]
    comparison_df = comparison_df[order(comparison_df[[group_col]]),]
  }
  relative_burden_df = data.frame(Run_Number=reference_df$Run_Number)
  relative_burden_df[[group_col]] = reference_df[[group_col]]
  # iterate through burden indicators, calculating relative burden and adding to dataframe
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) %in% c(group_col, 'Run_Number'))]
  for(bb in 1:length(burden_indicators)){
    # relative_burden_cur = (comparison_df[[burden_indicators[bb]]] - reference_df[[burden_indicators[bb]]]) / reference_df[[burden_indicators[bb]]]
    relative_burden_cur = (reference_df[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]]) / reference_df[[burden_indicators[bb]]]
    relative_burden_df[[burden_indicators[bb]]] = relative_burden_cur
  }
  relative_burden_df$scenario = comparison_scenario_name
  return(relative_burden_df)
}




get_relative_U1_burden = function(sim_output_filepath, reference_experiment_name, comparison_experiment_name, comparison_scenario_name, start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all', LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE){
  #'  @description get relative change in U1 burden when comparing between two simulations in specified years and in specified districts (for all malaria metrics, separate values for each seed)
  #'  @return data frame where each row is a seed and each column is the relative change of different burden metrics, calculated as (reference-comparison) / reference:
  
  reference_df = get_cumulative_U1_burden(sim_output_filepath=sim_output_filepath, experiment_name=reference_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_U1_burden(sim_output_filepath=sim_output_filepath, experiment_name=comparison_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  
  if(align_seeds){  # compare one run seed against the matching run seed in the other experiment
    # align seeds
    reference_df = reference_df[order(reference_df$Run_Number),]
    comparison_df = comparison_df[order(comparison_df$Run_Number),]
  } else{  # compare the averages across all seeds from one experiment against the average from the other experiment
    reference_df = reference_df %>% summarise_all(mean)
    comparison_df = comparison_df %>% summarise_all(mean)
  }
  relative_burden_df = data.frame('Run_Number' = reference_df$Run_Number)
  # iterate through burden indicators, calculating relative burden and adding to dataframe
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) == 'Run_Number')]
  for(bb in 1:length(burden_indicators)){
    # relative_burden_cur = (comparison_df[[burden_indicators[bb]]] - reference_df[[burden_indicators[bb]]]) / reference_df[[burden_indicators[bb]]]
    relative_burden_cur = (reference_df[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]]) / reference_df[[burden_indicators[bb]]]
    relative_burden_df[[burden_indicators[bb]]] = relative_burden_cur
  }
  relative_burden_df$scenario = comparison_scenario_name
  return(relative_burden_df)
}





####################################################################################
# absolute difference in simulation burden between two experiments over specified time interval
####################################################################################

get_difference_burden = function(sim_output_filepath, reference_experiment_name, comparison_experiment_name, comparison_scenario_name, start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all', 
                               LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE ){
  #'  @description get relative change in U5 and all-age burden when comparing between two simulations in specified years and in specified districts (for all malaria metrics, separate values for each seed)
  #'  @return data frame where each row is a seed and each column is the relative change of different burden metrics, calculated as (reference-comparison) / reference:
  
  reference_df = get_cumulative_burden(sim_output_filepath=sim_output_filepath, experiment_name=reference_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_burden(sim_output_filepath=sim_output_filepath, experiment_name=comparison_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  
  if(align_seeds){  # compare one run seed against the matching run seed in the other experiment
    # align seeds
    reference_df = reference_df[order(reference_df$Run_Number),]
    comparison_df = comparison_df[order(comparison_df$Run_Number),]
  } else{  # compare the averages across all seeds from one experiment against the average from the other experiment
    reference_df = reference_df %>% summarise_all(mean)
    comparison_df = comparison_df %>% summarise_all(mean)
  }
  difference_burden_df = data.frame('Run_Number' = reference_df$Run_Number)
  # iterate through burden indicators, calculating relative burden and adding to dataframe
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) == 'Run_Number')]
  for(bb in 1:length(burden_indicators)){
    difference_burden_cur = (reference_df[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]])
    difference_burden_df[[burden_indicators[bb]]] = difference_burden_cur
  }
  difference_burden_df$scenario = comparison_scenario_name
  return(difference_burden_df)
}



# get relative burden by state for grid plot
get_difference_burden_by_state = function(sim_output_filepath, reference_experiment_name, comparison_experiment_name, comparison_scenario_name, start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all', 
                                        LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE ){
  #'  @description get relative change in U5 and all-age burden when comparing between two simulations in specified years and in specified districts (for all malaria metrics, separate values for each seed)
  #'  @return data frame where each row is a seed and each column is the relative change of different burden metrics, calculated as (reference-comparison) / reference:
  
  reference_df = get_cumulative_burden_by_state(sim_output_filepath=sim_output_filepath, experiment_name=reference_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_burden_by_state(sim_output_filepath=sim_output_filepath, experiment_name=comparison_experiment_name, start_year=start_year, end_year=end_year, admin_pop=admin_pop, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  
  if(align_seeds){  # compare one run seed against the matching run seed in the other experiment
    # align seeds and states into same order
    reference_df = reference_df[order(reference_df$State, reference_df$Run_Number),]
    comparison_df = comparison_df[order(comparison_df$State, comparison_df$Run_Number),]
  } else{  # compare the averages across all seeds from one experiment against the average from the other experiment
    reference_df = reference_df %>% group_by(State) %>% summarise_all(mean)
    comparison_df = comparison_df %>% group_by(State) %>% summarise_all(mean)
    # align states into same order
    reference_df = reference_df[order(reference_df$State),]
    comparison_df = comparison_df[order(comparison_df$State),]
  }
  difference_burden_df = data.frame('State'=reference_df$State, 'Run_Number'=reference_df$Run_Number)
  # iterate through burden indicators, calculating relative burden and adding to dataframe
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) %in% c('State', 'Run_Number'))]
  for(bb in 1:length(burden_indicators)){
    difference_burden_cur = (reference_df[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]])
    difference_burden_df[[burden_indicators[bb]]] = difference_burden_cur
  }
  difference_burden_df$scenario = comparison_scenario_name
  return(difference_burden_df)
}






####################################################################################
# timeseries of simulation burden and/or interventions
####################################################################################

# Helper: flag years built from fewer than the expected number of months (max observed months/year).
#   set_na=TRUE  -> blank those years (used for SUM-based annual totals, which would be undercounted and
#                   misleadingly read as 'lower burden'); shown as a gap rather than a low value.
#   set_na=FALSE -> warn only (used for MEAN-based annual values, which are not undercounted but may be
#                   seasonally biased if months are missing).
# Expects a data frame with 'year' and 'n_months' columns; returns it with 'n_months' dropped.
.flag_incomplete_annual = function(df, exp_name, burden_colname, set_na){
  exp_mo = suppressWarnings(max(df$n_months, na.rm=TRUE))
  if(is.finite(exp_mo)){
    inc = which(df$n_months < exp_mo)
    if(length(inc)){
      it = unique(df[inc, c('year','n_months')]); it = it[order(it$year),]
      verb = if(set_na) 'set to NA (annual sum would be undercounted)' else 'averaged over partial months (possible seasonal bias)'
      warning(sprintf('get_burden_timeseries_exp [%s | %s]: incomplete year(s) %s: %s; expected %d months/year.',
                      exp_name, burden_colname, verb,
                      paste(sprintf('%s = %d mo', it$year, it$n_months), collapse='; '), exp_mo), call.=FALSE)
      if(set_na) df$burden[inc] = NA
    }
  }
  df$n_months = NULL
  df
}

get_burden_timeseries_exp = function(exp_filepath, exp_name, district_subset, cur_admins, pop_sizes, min_year, max_year, burden_colname, age_plotted, plot_by_month=FALSE, overwrite_files=FALSE){
  #'  @description subset simulation output to appropriate admin and time period, and calculate monthly or annual mean, min, max burdens (for a specified malaria burden metric) across all runs
  #'  @return data frame where each row is a time point and there are columns for the mean, minimum, and maximum burden value across seeds, and also a column for the scenario name
  
  # check whether file already exists, otherwise create new dataframe
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  output_filename = paste0(exp_filepath, '/timeseries_', burden_colname, '_', district_subset, '_', time_string, '_', min_year, '_', max_year, '.csv')
  if(file.exists(output_filename) & !overwrite_files){
    cur_sim_output_agg = read.csv(output_filename)
  } else{
    # read in simulation information, subset to appropriate years
    cur_sim_output = fread(file = paste0(exp_filepath, '/malariaBurden_withAdjustments.csv'))
    cur_sim_output = cur_sim_output[intersect(which(cur_sim_output$year >= min_year), which(cur_sim_output$year <= max_year)),]
    # subset to appropriate admins
    cur_sim_output = cur_sim_output[cur_sim_output$admin_name %in% cur_admins,]  
    
    # merge to get real-world population sizes in each admin
    cur_sim_output = merge(cur_sim_output, pop_sizes, by='admin_name')
    
    # get simulation population denominator and the real-world population size in each admin
    if(age_plotted == 'U1'){
      cur_sim_output$population = cur_sim_output$Pop_U1
      cur_sim_output$true_population = cur_sim_output$pop_size * cur_sim_output$Pop_U1 / cur_sim_output$Statistical_Population
    } else if(age_plotted == 'U5'){
      cur_sim_output$population = cur_sim_output$Pop_U5
      cur_sim_output$true_population = cur_sim_output$pop_size * cur_sim_output$Pop_U5 / cur_sim_output$Statistical_Population
    } else if(age_plotted == 'births'){
      cur_sim_output$population = cur_sim_output$num_all_births
      cur_sim_output$true_population = cur_sim_output$pop_size * cur_sim_output$num_all_births / cur_sim_output$Statistical_Population
    } else{
      cur_sim_output$population = cur_sim_output$Statistical_Population
      cur_sim_output$true_population = cur_sim_output$pop_size
    }

    
    # calculate average over included admins
    if(grepl('PfPR', burden_colname)){
      # get total number of positives (PfPR * true population for all included admins), then divide by total population size of all admins
      cur_sim_output$positives = cur_sim_output[[burden_colname]] * cur_sim_output$true_population
      cur_sim_output_agg_admin = as.data.frame(cur_sim_output) %>% dplyr::select(month, year, date, admin_name, positives, true_population, Run_Number) %>%
        dplyr::group_by(month, year, date, Run_Number) %>%  # take sum of positives and population across admins
        dplyr::summarise(total_positives = sum(positives),
                         total_population = sum(true_population))
      cur_sim_output_agg_admin$burden = cur_sim_output_agg_admin$total_positives / cur_sim_output_agg_admin$total_population
      if(!plot_by_month){
        cur_sim_output_agg_admin = cur_sim_output_agg_admin %>% dplyr::group_by(year, Run_Number) %>% # take average PfPR across months
          dplyr::summarise(burden = mean(burden), n_months = dplyr::n(), .groups='drop')
        cur_sim_output_agg_admin = .flag_incomplete_annual(cur_sim_output_agg_admin, exp_name, burden_colname, set_na=FALSE)
      }
      
    } else if((grepl('mLBW', burden_colname)) | (grepl('stillbirth', burden_colname))){
      # rescale to number present in full admin (with true population instead of simulated population size)
      cur_sim_output$true_burden = cur_sim_output[[burden_colname]] * (cur_sim_output$true_population / cur_sim_output$population)
      # take sum of the number of mLBWs or stillbirths and births across all included admins
      select_col_names = c('true_burden', 'month', 'year', 'date', 'true_population', 'Run_Number')
      cur_sim_output_agg_admin = as.data.frame(cur_sim_output) %>% dplyr::select(match(select_col_names, names(.))) %>%
        dplyr::group_by(month, year, date, Run_Number) %>%
        dplyr::summarise_all(sum)
      cur_sim_output_agg_admin$burden = cur_sim_output_agg_admin$true_burden / cur_sim_output_agg_admin$true_population * 1000
      if(!plot_by_month){
        cur_sim_output_agg_admin = cur_sim_output_agg_admin %>% dplyr::group_by(year, Run_Number) %>% # take average mLBW or stillbirths per birth across months to get annual average rate
          dplyr::summarise(burden = mean(burden), n_months = dplyr::n(), .groups='drop')
        cur_sim_output_agg_admin = .flag_incomplete_annual(cur_sim_output_agg_admin, exp_name, burden_colname, set_na=FALSE)
      }
    }else{
      # rescale case (or death) numbers to number present in full admin (with true population instead of simulated population size)
      cur_sim_output$true_burden = cur_sim_output[[burden_colname]] * (cur_sim_output$true_population / cur_sim_output$population)
      # take sum of the number of cases (or deaths) and population sizes across all included admins
      select_col_names = c('true_burden', 'month', 'year', 'date', 'true_population', 'Run_Number')
      cur_sim_output_agg_admin = as.data.frame(cur_sim_output) %>% dplyr::select(match(select_col_names, names(.))) %>%
        dplyr::group_by(month, year, date, Run_Number) %>%
        dplyr::summarise_all(sum)
      cur_sim_output_agg_admin$burden = cur_sim_output_agg_admin$true_burden / cur_sim_output_agg_admin$true_population * 1000
      if(!plot_by_month){
        cur_sim_output_agg_admin = cur_sim_output_agg_admin %>% dplyr::group_by(year, Run_Number) %>% # take sum across months to get annual burden
          dplyr::summarise(burden = sum(burden), n_months = dplyr::n(), .groups='drop')
        cur_sim_output_agg_admin = .flag_incomplete_annual(cur_sim_output_agg_admin, exp_name, burden_colname, set_na=TRUE)
      }
    }
    
    # take average, max, and min burdens across simulation seeds. Track how many seeds contributed to
    # each time point: groups with NO data become NA (a visible gap), and points summarised over only
    # SOME of the seeds are warned about (the mean is then over the available seeds, not all of them).
    n_runs_total = length(unique(cur_sim_output_agg_admin$Run_Number))
    if(plot_by_month){
      cur_sim_output_agg = as.data.frame(cur_sim_output_agg_admin) %>% dplyr::select(month, year, date, burden) %>%
        dplyr::group_by(month, year, date) %>%
        dplyr::summarise(n_runs_present = sum(!is.na(burden)),
                         mean_burden = if(all(is.na(burden))) NA_real_ else mean(burden, na.rm=TRUE),
                         max_burden  = if(all(is.na(burden))) NA_real_ else max(burden, na.rm=TRUE),
                         min_burden  = if(all(is.na(burden))) NA_real_ else min(burden, na.rm=TRUE), .groups='drop')
    } else{
      cur_sim_output_agg = as.data.frame(cur_sim_output_agg_admin) %>% dplyr::select(year, burden) %>%
        dplyr::group_by(year) %>%
        dplyr::summarise(n_runs_present = sum(!is.na(burden)),
                         mean_burden = if(all(is.na(burden))) NA_real_ else mean(burden, na.rm=TRUE),
                         max_burden  = if(all(is.na(burden))) NA_real_ else max(burden, na.rm=TRUE),
                         min_burden  = if(all(is.na(burden))) NA_real_ else min(burden, na.rm=TRUE), .groups='drop')
    }
    partial = which(cur_sim_output_agg$n_runs_present > 0 & cur_sim_output_agg$n_runs_present < n_runs_total)
    if(length(partial)){
      warning(sprintf('get_burden_timeseries_exp [%s | %s]: %d time point(s) summarised over only some of %d runs (mean over available runs): year(s) %s.',
                      exp_name, burden_colname, length(partial), n_runs_total,
                      paste(sprintf('%s (%d run[s])', cur_sim_output_agg$year[partial], cur_sim_output_agg$n_runs_present[partial]), collapse='; ')), call.=FALSE)
    }
    cur_sim_output_agg$n_runs_present = NULL
    
    cur_sim_output_agg$scenario = exp_name
    write.csv(cur_sim_output_agg, output_filename, row.names=FALSE)
  }
  return(cur_sim_output_agg)
}



get_burden_timeseries_by_lga = function(exp_filepath, exp_name, pop_filepath, overwrite_files=FALSE){
  #'  @description subset simulation output to appropriate admin and time period, and calculate annual mean burden (for a set of malaria burden metric) across all runs
  #'  @return data frame where each row is a time point and there are columns for the mean, minimum, and maximum burden value across seeds, and also a column for the scenario name
  # Assumes malariaBurden_withAdjustments.csv contains the following burden colnames: c('PfPR_U5', 'PfPR_MiP_adjusted', 'New_clinical_cases_U5', 'New_Clinical_Cases', 'direct_mortality_nonMiP_U5_mean', 'direct_mortality_nonMiP_mean', 'total_mortality_U5_mean', 'total_mortality_mean')    
  
  # check whether file already exists, otherwise create new dataframe
  output_filename = paste0(exp_filepath, '/timeseries_burden_annual_by_LGA.csv')
  if(file.exists(output_filename) & !overwrite_files){
    cur_sim_output_agg = read.csv(output_filename)
  } else{
    # read in information about LGAs
    admin_info = read.csv(pop_filepath)
    admin_info = admin_info[,c('admin_name','pop_size','State')]
    
    # read in simulation information, subset to appropriate years
    cur_sim_output = fread(file = paste0(exp_filepath, '/malariaBurden_withAdjustments.csv'))

    # merge to get real-world population sizes in each admin and the State each admin belongs to
    cur_sim_output = merge(cur_sim_output, admin_info, by='admin_name')
    # get simulation population denominator and the real-world population size in each admin
    cur_sim_output$true_population_U5 = cur_sim_output$pop_size * cur_sim_output$Pop_U5 / cur_sim_output$Statistical_Population
    cur_sim_output$true_population_all = cur_sim_output$pop_size
    
    # process simulation output to get total numbers in each state:
    #  - get total numbers in each admin, scaled to appropriate LGA population size from simulation population size
    #  - subset to relevant columns
    #  - get total within each admin-year-Run (across months)
    cur_sim_output_a = cur_sim_output %>% mutate(
      positives_U5 = PfPR_U5 * true_population_U5,
      positives_all = PfPR_MiP_adjusted * true_population_all,
      num_cases_U5 = New_clinical_cases_U5 / Pop_U5 * true_population_U5, 
      num_cases_all = New_Clinical_Cases / Statistical_Population * true_population_all,
      num_direct_mortality_U5 = direct_mortality_nonMiP_U5_mean / Pop_U5 * true_population_U5, 
      num_direct_mortality_all = direct_mortality_nonMiP_mean / Statistical_Population * true_population_all,
      num_total_mortality_U5 = total_mortality_U5_mean / Pop_U5 * true_population_U5, 
      num_total_mortality_all = total_mortality_mean / Statistical_Population * true_population_all
    ) %>% 
      dplyr::select(Run_Number, admin_name, State, year, true_population_U5, true_population_all, positives_U5, positives_all, num_cases_U5, num_cases_all, num_direct_mortality_U5, num_direct_mortality_all, num_total_mortality_U5, num_total_mortality_all) %>%
      group_by(Run_Number, admin_name, State, year) %>%
      summarise_all(sum) %>%
      ungroup()  %>% 
      mutate(  # take average population over 12 months rather than sum of population sizes in each month
        true_population_U5 = true_population_U5 / 12,
        true_population_all = true_population_all / 12
      )
    # process simulation output to get per-capita values and average across runs:
    #  - get per-capita or rate values
    #  - find average across runs
    cur_sim_output_agg = cur_sim_output_a %>% mutate(
      PfPR_U5 = positives_U5 / true_population_U5 / 12,  # divide by number of months in year to get population-weighted average prevalence across all months rather than sum of positives across months
      PfPR_all = positives_all / true_population_all / 12,
      incidence_pp_U5 = num_cases_U5 / true_population_U5, 
      incidence_pp_all = num_cases_all / true_population_all,
      direct_mortality_pp_U5 = num_direct_mortality_U5 / true_population_U5, 
      direct_mortality_pp_all = num_direct_mortality_all / true_population_all,
      total_mortality_pp_U5 = num_total_mortality_U5 / true_population_U5, 
      total_mortality_pp_all = num_total_mortality_all / true_population_all
    ) %>%
      group_by(admin_name, State, year) %>%
      summarise_all(mean) %>%
      ungroup()
    
    # save result
    cur_sim_output_agg$scenario = exp_name
    write.csv(cur_sim_output_agg, output_filename, row.names=FALSE)
  }
  return(cur_sim_output_agg)
}





get_burden_timeseries_by_state = function(exp_filepath, exp_name, pop_filepath, group_col='State', overwrite_files=FALSE){
  #'  @description subset simulation output to appropriate admin and time period, and calculate annual mean burden (for a set of malaria burden metric) across all runs
  #'  @param group_col column in the pop file used to group LGAs (default 'State'; e.g. 'Funder'). Rates are
  #'        population-weighted within each group (Sum of numerators / Sum of population). Cache filename uses
  #'        tolower(group_col) so 'State' keeps the existing 'by_state' name.
  #'  @return data frame where each row is a time point and there are columns for the mean, minimum, and maximum burden value across seeds, and also a column for the scenario name
  # Assumes malariaBurden_withAdjustments.csv contains the following burden colnames: c('PfPR_U5', 'PfPR_MiP_adjusted', 'New_clinical_cases_U5', 'New_Clinical_Cases', 'direct_mortality_nonMiP_U5_mean', 'direct_mortality_nonMiP_mean', 'total_mortality_U5_mean', 'total_mortality_mean')

  # check whether file already exists, otherwise create new dataframe
  output_filename = paste0(exp_filepath, '/timeseries_burden_annual_by_', tolower(group_col), '.csv')
  if(file.exists(output_filename) & !overwrite_files){
    cur_sim_output_agg = read.csv(output_filename)
  } else{
    # read in information about LGAs
    admin_info = read.csv(pop_filepath)
    admin_info = admin_info[,c('admin_name','pop_size',group_col)]

    # read in simulation information, subset to appropriate years
    cur_sim_output = fread(file = paste0(exp_filepath, '/malariaBurden_withAdjustments.csv'))

    # merge to get real-world population sizes in each admin and the State each admin belongs to
    cur_sim_output = merge(cur_sim_output, admin_info, by='admin_name')
    # get simulation population denominator and the real-world population size in each admin
    cur_sim_output$true_population_U5 = cur_sim_output$pop_size * cur_sim_output$Pop_U5 / cur_sim_output$Statistical_Population
    cur_sim_output$true_population_all = cur_sim_output$pop_size

    # process simulation output to get total numbers in each state:
    #  - get total numbers in each admin, scaled to appropriate LGA population size from simulation population size
    #  - subset to relevant columns
    #  - get total within each state-year-Run (across months and LGAs)
    cur_sim_output_a = cur_sim_output %>% mutate(
        positives_U5 = PfPR_U5 * true_population_U5,
        positives_all = PfPR_MiP_adjusted * true_population_all,
        num_cases_U5 = New_clinical_cases_U5 / Pop_U5 * true_population_U5,
        num_cases_all = New_Clinical_Cases / Statistical_Population * true_population_all,
        num_direct_mortality_U5 = direct_mortality_nonMiP_U5_mean / Pop_U5 * true_population_U5,
        num_direct_mortality_all = direct_mortality_nonMiP_mean / Statistical_Population * true_population_all,
        num_total_mortality_U5 = total_mortality_U5_mean / Pop_U5 * true_population_U5,
        num_total_mortality_all = total_mortality_mean / Statistical_Population * true_population_all
      ) %>%
      dplyr::select(Run_Number, all_of(group_col), year, true_population_U5, true_population_all, positives_U5, positives_all, num_cases_U5, num_cases_all, num_direct_mortality_U5, num_direct_mortality_all, num_total_mortality_U5, num_total_mortality_all) %>%
      group_by(Run_Number, across(all_of(group_col)), year) %>%
      summarise_all(sum) %>%
      ungroup() %>%
      mutate(  # take average population over 12 months rather than sum over all months
        true_population_U5 = true_population_U5 / 12,
        true_population_all = true_population_all / 12
      )
    # process simulation output to get per-capita values and average across runs:
    #  - get per-capita or rate values
    #  - find average across runs
    cur_sim_output_agg = cur_sim_output_a %>% mutate(
      PfPR_U5 = positives_U5 / true_population_U5 / 12,  # divide by number of months in year to get population-weighted average prevalence across all months rather than sum of positives across months
      PfPR_all = positives_all / true_population_all / 12,
      incidence_pp_U5 = num_cases_U5 / true_population_U5,
      incidence_pp_all = num_cases_all / true_population_all,
      direct_mortality_pp_U5 = num_direct_mortality_U5 / true_population_U5,
      direct_mortality_pp_all = num_direct_mortality_all / true_population_all,
      total_mortality_pp_U5 = num_total_mortality_U5 / true_population_U5,
      total_mortality_pp_all = num_total_mortality_all / true_population_all
      ) %>%
      group_by(across(all_of(group_col)), year) %>%
      summarise_all(mean) %>%
      ungroup()

    # save result
    cur_sim_output_agg$scenario = exp_name
    write.csv(cur_sim_output_agg, output_filename, row.names=FALSE)
  }
  return(cur_sim_output_agg)
}




# ----- LLIN, IRS, vaccine, PMC intervention coverage by state ----- #
get_intervention_use_timeseries_state = function(exp_filepath, exp_name, pop_sizes, min_year, max_year, indoor_protection_fraction, plot_by_month=FALSE){
  #'  @description subset simulation output to appropriate admin(s) and time period, and calculate monthly net usage, net distribution, and IRS coverage across all runs
  #'  @return data frame where each row is a time point and there are columns for the mean, minimum, and maximum ITN coverage across seeds, as well as columns for mean nets distributed and IRS coverages, and also a column for the scenario name
  
  # read in simulation information, merge to single dataframe, subset to appropriate years
  interv_use_all = fread(file = paste0(exp_filepath, '/MonthlyUsageLLIN.csv'))
  interv_dist_all = fread(file = paste0(exp_filepath, '/monthly_Event_Count.csv'))
  intervention_use_columns_orig = c('Bednet_Using')
  intervention_use_columns_new_name = c('itn_used_per_cap')
  intervention_use_columns = intervention_use_columns_orig[intervention_use_columns_orig %in% colnames(interv_use_all)]
  intervention_distribution_columns_orig = c('Received_IRS', 'Received_Vaccine', 'Received_Campaign_Drugs', 'Received_PMC_VaccDrug', 'Bednet_Got_New_One')
  intervention_distribution_columns_new_name = c('irs_per_cap', 'vacc_per_cap', 'smc_per_cap', 'pmc_per_cap', 'new_net_per_cap')
  intervention_distribution_columns = intervention_distribution_columns_orig[intervention_distribution_columns_orig %in% colnames(interv_dist_all)]
  intervention_columns = c(intervention_use_columns, intervention_distribution_columns)
  interv_dist_all = interv_dist_all %>% dplyr::select(one_of(c('admin_name', 'date', 'Run_Number', intervention_distribution_columns)))
  interv_all = merge(interv_use_all, interv_dist_all, by=c('admin_name', 'date', 'Run_Number'))
  interv_all$date = as.Date(interv_all$date)
  interv_all$year = lubridate::year(interv_all$date)
  interv_all = interv_all[intersect(which(interv_all$year >= min_year), which(interv_all$year <= max_year)),]
  colnames(interv_all) = gsub(' ','.',colnames(interv_all))
  
  # rescale numbers to the true population in each admin, rather than the simulated population size
  interv_all = merge(interv_all, pop_sizes, by='admin_name')
  colnames(interv_all)[colnames(interv_all) == 'pop_size'] = 'true_population'
  interv_pop_scaled = as.data.frame(interv_all)
  interv_pop_scaled[,intervention_columns] = interv_pop_scaled[,intervention_columns]* (interv_pop_scaled$true_population / interv_pop_scaled$Statistical.Population)
  
  # get sum of numbers across all included admins (keeping runs and months separate)
  interv_sums = interv_pop_scaled %>% dplyr::select(one_of('year', 'date', 'Run_Number', 'State', 'true_population', intervention_columns)) %>% dplyr::group_by(date, year, State, Run_Number) %>% 
    dplyr::summarise_all(sum) %>% ungroup()
  interv_per_capita = interv_sums
  interv_per_capita[,intervention_columns] = interv_per_capita[,intervention_columns] / interv_per_capita$true_population
  
  # for bednets, rescale use rates to un-adjust for the fraction of time spent under a net
  if('Bednet_Using' %in% colnames(interv_per_capita)) interv_per_capita$Bednet_Using = interv_per_capita$Bednet_Using / indoor_protection_fraction
  
  # get average coverage across months in a year for the annual report
  if(!plot_by_month){
    interv_use_sums = interv_per_capita %>% dplyr::select(one_of('year', 'Run_Number', 'State', intervention_use_columns)) %>% dplyr::group_by(year, State, Run_Number) %>% # take average across months for annual coverage
      dplyr::summarise_all(mean) %>% ungroup()
    interv_dist_sums = interv_per_capita %>% dplyr::select(one_of('year', 'Run_Number', 'State', intervention_distribution_columns)) %>% dplyr::group_by(year, State, Run_Number) %>% # take sum across months for annual total
      dplyr::summarise_all(sum) %>% ungroup()
    interv_per_capita = merge(interv_use_sums, interv_dist_sums, all=TRUE)
  } 
  
  # take average across simulation seeds
  if(plot_by_month){
    interv_agg = interv_per_capita %>% dplyr::select(one_of('year', 'date', 'State', intervention_columns)) %>% 
      dplyr::group_by(year, date, State) %>%
      dplyr::summarise_all(mean)
  } else{
    interv_agg = interv_per_capita %>% dplyr::select(one_of('year', 'State', intervention_columns)) %>%  # , new_net_per_cap
      dplyr::group_by(year, State) %>%
      dplyr::summarise_all(mean)
  }
  
  # update column names for plotting
  for(cc in intervention_columns){
    if(cc %in% intervention_use_columns){
      colnames(interv_agg)[colnames(interv_agg)==cc] = intervention_use_columns_new_name[which(intervention_use_columns_orig==cc)]
    } else{
      colnames(interv_agg)[colnames(interv_agg)==cc] = intervention_distribution_columns_new_name[which(intervention_distribution_columns_orig==cc)]
    }
  }
  
  interv_agg$scenario = exp_name
  return(interv_agg)
}


# ----- LLIN, IRS, vaccine, PMC intervention coverage nationally ----- #
get_intervention_use_timeseries_exp = function(exp_filepath, exp_name, cur_admins, pop_sizes, min_year, max_year, indoor_protection_fraction, plot_by_month=TRUE){
  #'  @description subset simulation output to appropriate admin(s) and time period, and calculate monthly net usage, net distribution, and IRS coverage across all runs
  #'  @return data frame where each row is a time point and there are columns for the mean, minimum, and maximum ITN coverage across seeds, as well as columns for mean nets distributed and IRS coverages, and also a column for the scenario name
  
  # read in simulation information, merge to single dataframe, subset to appropriate years
  interv_use_all = fread(file = paste0(exp_filepath, '/MonthlyUsageLLIN.csv'))
  interv_dist_all = fread(file = paste0(exp_filepath, '/monthly_Event_Count.csv'))
  intervention_use_columns_orig = c('Bednet_Using')
  intervention_use_columns_new_name = c('coverage')
  intervention_use_columns = intervention_use_columns_orig[intervention_use_columns_orig %in% colnames(interv_use_all)]
  intervention_distribution_columns_orig = c('Received_IRS', 'Received_Vaccine', 'Received_Campaign_Drugs', 'Received_PMC_VaccDrug', 'Bednet_Got_New_One')
  intervention_distribution_columns_new_name = c('irs_per_cap', 'vacc_per_cap', 'smc_per_cap', 'pmc_per_cap', 'new_net_per_cap')
  intervention_distribution_columns = intervention_distribution_columns_orig[intervention_distribution_columns_orig %in% colnames(interv_dist_all)]
  intervention_columns = c(intervention_use_columns, intervention_distribution_columns)
  interv_dist_all = interv_dist_all %>% dplyr::select(one_of(c('admin_name', 'date', 'Run_Number', intervention_distribution_columns)))
  interv_all = merge(interv_use_all, interv_dist_all, by=c('admin_name', 'date', 'Run_Number'))
  interv_all$date = as.Date(interv_all$date)
  interv_all$year = lubridate::year(interv_all$date)
  interv_all = interv_all[intersect(which(interv_all$year >= min_year), which(interv_all$year <= max_year)),]

  # subset to appropriate admins
  interv_all = interv_all[interv_all$admin_name %in% cur_admins,]  
  colnames(interv_all) = gsub(' ','.',colnames(interv_all))
  
  # rescale numbers to the true population in each admin, rather than the simulated population size
  interv_all = merge(interv_all, pop_sizes, by='admin_name')
  colnames(interv_all)[colnames(interv_all) == 'pop_size'] = 'true_population'
  interv_pop_scaled = as.data.frame(interv_all)
  interv_pop_scaled[,intervention_columns] = interv_pop_scaled[,intervention_columns]* (interv_pop_scaled$true_population / interv_pop_scaled$Statistical.Population)
  
  # get sum of numbers across all included admins (keeping runs and months separate)
  interv_sums = interv_pop_scaled %>% dplyr::select(one_of('year', 'date', 'Run_Number', 'true_population', intervention_columns)) %>% dplyr::group_by(date, year, Run_Number) %>% 
    dplyr::summarise_all(sum) %>% ungroup()
  interv_per_capita = interv_sums
  interv_per_capita[,intervention_columns] = interv_per_capita[,intervention_columns] / interv_per_capita$true_population
  
  # for bednets, rescale use rates to un-adjust for the fraction of time spent under a net
  if('Bednet_Using' %in% colnames(interv_per_capita)) interv_per_capita$Bednet_Using = interv_per_capita$Bednet_Using / indoor_protection_fraction
  
  # get average coverage across months in a year for the annual report
  if(!plot_by_month){
    interv_use_sums = interv_per_capita %>% dplyr::select(one_of('year', 'Run_Number', intervention_use_columns)) %>% dplyr::group_by(year, Run_Number) %>% # take average across months for annual coverage
      dplyr::summarise_all(mean) %>% ungroup()
    interv_dist_sums = interv_per_capita %>% dplyr::select(one_of('year', 'Run_Number', intervention_distribution_columns)) %>% dplyr::group_by(year, Run_Number) %>% # take sum across months for annual total
      dplyr::summarise_all(sum) %>% ungroup()
    interv_per_capita = merge(interv_use_sums, interv_dist_sums, all=TRUE)
  } 
  
  # take average across simulation seeds
  if(plot_by_month){
    interv_agg = interv_per_capita %>% dplyr::select(one_of('year', 'date', intervention_columns)) %>% 
      dplyr::group_by(year, date) %>%
      dplyr::summarise_all(mean)
  } else{
    interv_agg = interv_per_capita %>% dplyr::select(one_of('year', intervention_columns)) %>%  # , new_net_per_cap
      dplyr::group_by(year) %>%
      dplyr::summarise_all(mean)
  }
  
  # update column names for plotting
  for(cc in intervention_columns){
    if(cc %in% intervention_use_columns){
      colnames(interv_agg)[colnames(interv_agg)==cc] = intervention_use_columns_new_name[which(intervention_use_columns_orig==cc)]
    } else{
      colnames(interv_agg)[colnames(interv_agg)==cc] = intervention_distribution_columns_new_name[which(intervention_distribution_columns_orig==cc)]
    }
  }
  
  interv_agg$scenario = exp_name
  return(interv_agg)
}


# ----- CM intervention coverage by state ----- #
get_cm_timeseries_by_state = function(cm_filepath, admin_info, end_year, exp_name, min_year, plot_by_month=TRUE){
  
  cm_input = read.csv(cm_filepath)
  if(!('seed' %in% cm_input)) cm_input$seed = 1
  cm_input = merge(cm_input, admin_info, by='admin_name')
  
  # CM is sometimes repeated for several years but only listed once; change to repeate the appropriate number of times
  cm_input$years_repeated = cm_input$duration/365
  if(any(cm_input$years_repeated>1)){
    cur_cm_years = unique(cm_input$year)
    for(yy in cur_cm_years){
      # get first instance of this year
      cur_year = cm_input[cm_input$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          cm_input = rbind(cm_input, temp_year)
        }
      }
    }
  }
  if(any(cm_input$duration==-1) & (max(cm_input$year)<end_year)){
    cm_repeated = cm_input[cm_input$duration == -1,]
    for(rr in 1:(end_year - cm_repeated$year[1])){
      temp_year = cm_repeated
      temp_year$year = cm_repeated$year + rr
      temp_year$simday = cm_repeated$simday + rr*365
      cm_input = rbind(cm_input, temp_year)
    }
  }
  # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
  cm_input <- cm_input %>% group_by(year, admin_name, State, seed) %>%
    summarise_all(mean) %>% ungroup()
  
  # cm_input = cm_input[intersect(which(cm_input$year >= min_year), which(cm_input$year <= end_year)),]
  
  
  # get population-weighted CM coverage across admins
  cm_input$multiplied_U5_cm = cm_input$U5_coverage * cm_input$pop_size
  
  # get sum of population sizes and multiplied CM coverage across included admins
  cm_input_agg_admin <- cm_input %>% dplyr::select(year, State, seed, multiplied_U5_cm, pop_size) %>% group_by(year, State, seed) %>%
    summarise_all(sum) %>% ungroup()
  # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
  cm_input_agg_admin$U5_coverage = cm_input_agg_admin$multiplied_U5_cm / cm_input_agg_admin$pop_size
  
  
  # take average, max, and min burdens across simulation seeds
  if(plot_by_month){
    # subdivide year values and add dates
    # date dataframe
    included_years = unique(cm_input_agg_admin$year)
    all_months = as.Date(paste0(rep(included_years, each=12),'-',c('01','02','03','04','05','06','07','08','09','10','11','12'), '-01' ))
    date_df = data.frame(year=rep(included_years, each=12), date=all_months)
    cm_input_agg_admin_monthly = merge(cm_input_agg_admin, date_df, by='year', all.x=TRUE, all.y=TRUE)
    cm_agg = as.data.frame(cm_input_agg_admin_monthly) %>% dplyr::select(year, State, date, U5_coverage) %>%
      dplyr::group_by(year, State, date) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    cm_agg = as.data.frame(cm_input_agg_admin) %>% dplyr::select(year, State, U5_coverage) %>%
      dplyr::group_by(year, State) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  cm_agg$scenario = exp_name
  return(cm_agg)
}


# ----- CM intervention coverage ----- #
get_cm_timeseries_exp = function(cm_filepath, pop_sizes, end_year, exp_name, cur_admins, min_year, plot_by_month=TRUE){
  
  cm_input = read.csv(cm_filepath)
  # subset to appropriate admins
  cm_input = cm_input[cm_input$admin_name %in% cur_admins,] 
  if(!('seed' %in% cm_input)) cm_input$seed = 1
  
  # CM is sometimes repeated for several years but only listed once; change to repeate the appropriate number of times
  cm_input$years_repeated = cm_input$duration/365
  if(any(cm_input$years_repeated>1)){
    cur_cm_years = unique(cm_input$year)
    for(yy in cur_cm_years){
      # get first instance of this year
      cur_year = cm_input[cm_input$year==yy,]
      if(cur_year$years_repeated[1]>1){
        for(rr in 1:(cur_year$years_repeated[1] - 1)){
          temp_year = cur_year
          temp_year$year = cur_year$year + rr
          temp_year$simday = cur_year$simday + rr*365
          cm_input = rbind(cm_input, temp_year)
        }
      }
    }
  }
  if(any(cm_input$duration==-1) & (max(cm_input$year)<end_year)){
    cm_repeated = cm_input[cm_input$duration == -1,]
    for(rr in 1:(end_year - cm_repeated$year[1])){
      temp_year = cm_repeated
      temp_year$year = cm_repeated$year + rr
      temp_year$simday = cm_repeated$simday + rr*365
      cm_input = rbind(cm_input, temp_year)
    }
  }
  # if there are multiple values in a single year (for a single DS/LGA), take the mean of those values
  cm_input <- cm_input %>% group_by(year, admin_name, seed) %>%
    summarise_all(mean) %>% ungroup()
  
  cm_input = cm_input[intersect(which(cm_input$year >= min_year), which(cm_input$year <= end_year)),]
  
  
  # get population-weighted CM coverage across admins
  cm_input = merge(cm_input, pop_sizes, by='admin_name')
  cm_input$multiplied_U5_cm = cm_input$U5_coverage * cm_input$pop_size
  
  # get sum of population sizes and multiplied CM coverage across included admins
  cm_input_agg_admin <- cm_input %>% dplyr::select(year, seed, multiplied_U5_cm, pop_size) %>% group_by(year, seed) %>%
    summarise_all(sum) %>% ungroup()
  # get population-weighted U5 coverage across all included admin by dividing by sum of population sizes
  cm_input_agg_admin$U5_coverage = cm_input_agg_admin$multiplied_U5_cm / cm_input_agg_admin$pop_size
  
  
  # take average, max, and min burdens across simulation seeds
  if(plot_by_month){
    # subdivide year values and add dates
    # date dataframe
    included_years = unique(cm_input_agg_admin$year)
    all_months = as.Date(paste0(rep(included_years, each=12),'-',c('01','02','03','04','05','06','07','08','09','10','11','12'), '-01' ))
    date_df = data.frame(year=rep(included_years, each=12), date=all_months)
    cm_input_agg_admin_monthly = merge(cm_input_agg_admin, date_df, by='year', all.x=TRUE, all.y=TRUE)
    cm_agg = as.data.frame(cm_input_agg_admin_monthly) %>% dplyr::select(year, date, U5_coverage) %>%
      dplyr::group_by(year, date) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  } else{
    cm_agg = as.data.frame(cm_input_agg_admin) %>% dplyr::select(year, U5_coverage) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(mean_coverage = mean(U5_coverage),
                       max_coverage = max(U5_coverage),
                       min_coverage = min(U5_coverage))
  }
  
  cm_agg$scenario = exp_name
  return(cm_agg)
  
}








# calculate average burden over specified timeperiod for each state
get_ave_burden_by_state = function(sim_future_output_dir, pop_filepath,
                                   barplot_start_year, barplot_end_year, 
                                   pyr, chw_cov,
                                   scenario_names, experiment_names, scenario_palette, group_col='State', LLIN2y_flag=FALSE, overwrite_files=FALSE,
                                   burden_metric_subset=c(), include_to_present=TRUE, file_suffix=''){
  # group_col: admin_pop column used to group LGAs (default 'State'; e.g. 'Funder'). Output filenames use
  # tolower(group_col) so 'State' reproduces the existing '_by_state_' names and 'Funder' -> '_by_funder_'.
  group_lc = tolower(group_col)
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')#, 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')#, 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }
  
  # if first experiment name is to-present, skip it
  if(include_to_present){
    exp_start_index = 2
  } else{
    exp_start_index = 1
  }
  # iterate through scenarios, calculating the burden for all metrics. Combine all scenarios into a dataframe (each scenario set in separate rows)
  burden_all_df = data.frame()
  for(ss in exp_start_index:length(scenario_names)){
    scenario_name_cur = scenario_names[ss]
    burden_df = get_cumulative_burden_by_state(sim_output_filepath=sim_future_output_dir, experiment_name=scenario_name_cur, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, group_col=group_col, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
    burden_df$scenario = scenario_name_cur

    # only save relevant columns for plotting
    burden_df = burden_df[,which(colnames(burden_df) %in% c('scenario', 'Run_Number', group_col, burden_colnames))]
    if(nrow(burden_all_df) == 0){
      burden_all_df = burden_df
    }else{
      burden_all_df = rbind(burden_all_df, burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  burden_all_df$scenario = factor(burden_all_df$scenario, levels=scenario_names[exp_start_index:length(scenario_names)])
  
  # get average across runs
  burden_average_df = burden_all_df %>% group_by(across(all_of(group_col)), scenario) %>%
    summarise_all(mean, na.rm=TRUE)

  write.csv(burden_average_df, paste0(sim_future_output_dir, '/_plots/','average_burden_by_', group_lc, '_', barplot_start_year,'_', barplot_end_year, file_suffix,'.csv'), row.names=FALSE)

  # if the burden metrics include PfPR, create separate table with groups as rows and scenarios as columns
  if('average_PfPR_all' %in% colnames(burden_all_df)){
    pfpr_df = burden_average_df %>% dplyr::select(all_of(group_col), scenario, average_PfPR_all) %>%
      mutate(average_PfPR_all = round(average_PfPR_all,2)) %>%
      pivot_wider(id_cols=all_of(group_col), names_from=scenario, values_from=average_PfPR_all)
    write.csv(pfpr_df, paste0(sim_future_output_dir, '/_plots/','average_PfPR_allAge_by_', group_lc, '_', barplot_start_year,'_', barplot_end_year, file_suffix,'.csv'), row.names=FALSE)
    # percentage-formatted version (e.g., '5.2%' instead of 0.05)
    pfpr_pct_df = burden_average_df %>% dplyr::select(all_of(group_col), scenario, average_PfPR_all) %>%
      mutate(average_PfPR_all = paste0(round(average_PfPR_all * 100, 1), '%')) %>%
      pivot_wider(id_cols=all_of(group_col), names_from=scenario, values_from=average_PfPR_all)
    write.csv(pfpr_pct_df, paste0(sim_future_output_dir, '/_plots/','average_PfPR_allAge_by_', group_lc, '_', barplot_start_year,'_', barplot_end_year, file_suffix,'_percent.csv'), row.names=FALSE)
  }
  if('average_PfPR_U5' %in% colnames(burden_all_df)){
    pfpr_df = burden_average_df %>% dplyr::select(all_of(group_col), scenario, average_PfPR_U5) %>%
      mutate(average_PfPR_U5 = round(average_PfPR_U5,2)) %>%
      pivot_wider(id_cols=all_of(group_col), names_from=scenario, values_from=average_PfPR_U5)
    write.csv(pfpr_df, paste0(sim_future_output_dir, '/_plots/','average_PfPR_U5_by_', group_lc, '_', barplot_start_year,'_', barplot_end_year, file_suffix,'.csv'), row.names=FALSE)
    # percentage-formatted version (e.g., '5.2%' instead of 0.05)
    pfpr_pct_df = burden_average_df %>% dplyr::select(all_of(group_col), scenario, average_PfPR_U5) %>%
      mutate(average_PfPR_U5 = paste0(round(average_PfPR_U5 * 100, 1), '%')) %>%
      pivot_wider(id_cols=all_of(group_col), names_from=scenario, values_from=average_PfPR_U5)
    write.csv(pfpr_pct_df, paste0(sim_future_output_dir, '/_plots/','average_PfPR_U5_by_', group_lc, '_', barplot_start_year,'_', barplot_end_year, file_suffix,'_percent.csv'), row.names=FALSE)
  }

}



####################################################################################
# relative/difference burden for a scenario compared to a reference experiment
#   over a different year range (e.g., future scenario vs. to-present baseline year)
####################################################################################

get_relative_burden_vs_ref_year = function(reference_sim_output_filepath, reference_experiment_name, ref_start_year, ref_end_year,
                                           comparison_sim_output_filepath, comparison_experiment_name, comparison_scenario_name,
                                           start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all',
                                           LLIN2y_flag=FALSE, overwrite_files=FALSE){
  #'  @description get relative change in burden comparing a future scenario to a reference experiment over a reference year range
  #'  @return data frame where each row is a comparison seed, columns are (ref_mean - comp_seed) / ref_mean

  reference_df = get_cumulative_burden(sim_output_filepath=reference_sim_output_filepath, experiment_name=reference_experiment_name,
                                       start_year=ref_start_year, end_year=ref_end_year, admin_pop=admin_pop,
                                       district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_burden(sim_output_filepath=comparison_sim_output_filepath, experiment_name=comparison_experiment_name,
                                        start_year=start_year, end_year=end_year, admin_pop=admin_pop,
                                        district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)

  # average reference across seeds (cross-year comparison makes seed alignment inapplicable)
  reference_mean = reference_df %>% summarise_all(mean)
  comparison_df = comparison_df[order(comparison_df$Run_Number),]

  relative_burden_df = data.frame('Run_Number' = comparison_df$Run_Number)
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) == 'Run_Number')]
  for(bb in 1:length(burden_indicators)){
    relative_burden_df[[burden_indicators[bb]]] = (reference_mean[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]]) / reference_mean[[burden_indicators[bb]]]
  }
  relative_burden_df$scenario = comparison_scenario_name
  return(relative_burden_df)
}


get_difference_burden_vs_ref_year = function(reference_sim_output_filepath, reference_experiment_name, ref_start_year, ref_end_year,
                                             comparison_sim_output_filepath, comparison_experiment_name, comparison_scenario_name,
                                             start_year, end_year, admin_pop, district_subset='allDistricts', cur_admins='all',
                                             LLIN2y_flag=FALSE, overwrite_files=FALSE){
  #'  @description get absolute difference in burden comparing a future scenario to a reference experiment over a reference year range
  #'  @return data frame where each row is a comparison seed, columns are ref_mean - comp_seed

  reference_df = get_cumulative_burden(sim_output_filepath=reference_sim_output_filepath, experiment_name=reference_experiment_name,
                                       start_year=ref_start_year, end_year=ref_end_year, admin_pop=admin_pop,
                                       district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
  comparison_df = get_cumulative_burden(sim_output_filepath=comparison_sim_output_filepath, experiment_name=comparison_experiment_name,
                                        start_year=start_year, end_year=end_year, admin_pop=admin_pop,
                                        district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)

  reference_mean = reference_df %>% summarise_all(mean)
  comparison_df = comparison_df[order(comparison_df$Run_Number),]

  difference_burden_df = data.frame('Run_Number' = comparison_df$Run_Number)
  burden_indicators = colnames(reference_df)[-which(colnames(reference_df) == 'Run_Number')]
  for(bb in 1:length(burden_indicators)){
    difference_burden_df[[burden_indicators[bb]]] = (reference_mean[[burden_indicators[bb]]] - comparison_df[[burden_indicators[bb]]])
  }
  difference_burden_df$scenario = comparison_scenario_name
  return(difference_burden_df)
}


# national-level burden summary for all scenarios saved to CSV
# accepts per-scenario output directories and year ranges so reference and future scenarios can differ
# ref_scenario_name: scenario used as denominator for percent-reduction CSV; defaults to scenario_names[1]
# national_population: if provided, generates a numbers_averted CSV for all-age incidence and mortality (counts = rate_per_1000 * population / 1000)
create_csv_national_burden_all_scenarios = function(scenario_names, experiment_names, sim_output_dirs,
                                                    start_years, end_years,
                                                    pop_filepath, district_subset='allDistricts', cur_admins='all',
                                                    output_dir, LLIN2y_flag=FALSE, overwrite_files=FALSE,
                                                    burden_metric_subset=c(), file_suffix='',
                                                    ref_scenario_name=NULL, national_population=NULL){
  admin_pop = read.csv(pop_filepath)

  burden_metrics  = c('PfPR','PfPR','incidence','incidence','directMortality','directMortality','allMortality','allMortality')
  burden_colnames = c('average_PfPR_U5','average_PfPR_all','incidence_U5','incidence_all',
                      'direct_death_rate_mean_U5','direct_death_rate_mean_all','all_death_rate_mean_U5','all_death_rate_mean_all')
  if(length(burden_metric_subset) >= 1){
    burden_colnames = burden_colnames[burden_metrics %in% burden_metric_subset]
  }

  burden_all_df = data.frame()
  for(ss in 1:length(scenario_names)){
    burden_df = get_cumulative_burden(sim_output_filepath=sim_output_dirs[ss], experiment_name=experiment_names[ss],
                                      start_year=start_years[ss], end_year=end_years[ss],
                                      admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins,
                                      LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
    burden_df$scenario = scenario_names[ss]
    burden_df = burden_df[, which(colnames(burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(burden_all_df) == 0) burden_all_df = burden_df else burden_all_df = rbind(burden_all_df, burden_df)
  }

  burden_all_df$scenario = factor(burden_all_df$scenario, levels=scenario_names)
  burden_mean = burden_all_df %>% dplyr::select(-Run_Number) %>% group_by(scenario) %>% summarise_all(mean, na.rm=TRUE)

  pfpr_cols = intersect(c('average_PfPR_U5', 'average_PfPR_all'), colnames(burden_mean))
  if(length(pfpr_cols) > 0){
    burden_mean = burden_mean %>% mutate(across(all_of(pfpr_cols), ~ round(. * 100, 1)))
  }

  # write.csv(burden_all_df, paste0(output_dir, '/national_burden_all_scenarios_bySeed', file_suffix, '.csv'), row.names=FALSE)
  write.csv(burden_mean, paste0(output_dir, '/national_burden_all_scenarios_summary', file_suffix, '.csv'), row.names=FALSE)

  # percent reduction and absolute difference relative to reference scenario (defaults to first scenario, expected to be the ref year)
  if(is.null(ref_scenario_name)) ref_scenario_name = as.character(scenario_names[1])
  ref_row = burden_mean[as.character(burden_mean$scenario) == ref_scenario_name, ]
  if(nrow(ref_row) == 1){
    numeric_cols = names(burden_mean)[sapply(burden_mean, is.numeric)]
    pct_reduction_df = burden_mean
    burden_averted_df = burden_mean
    for(col in numeric_cols){
      ref_val = as.numeric(ref_row[[col]])
      comp_vals = as.numeric(burden_mean[[col]])
      burden_averted_df[[col]] = round(ref_val - comp_vals, 2)
      if(!is.na(ref_val) && ref_val != 0){
        pct_reduction_df[[col]] = round((ref_val - comp_vals) / ref_val * 100, 1)
      } else {
        pct_reduction_df[[col]] = NA_real_
      }
    }
    write.csv(pct_reduction_df,  paste0(output_dir, '/national_burden_all_scenarios_pct_reduction', file_suffix, '.csv'), row.names=FALSE)
    write.csv(burden_averted_df, paste0(output_dir, '/national_burden_all_scenarios_burden_averted', file_suffix, '.csv'), row.names=FALSE)

    # numbers averted (total counts) for all-age incidence and mortality
    if(!is.null(national_population)){
      count_cols = intersect(c('incidence_all', 'all_death_rate_mean_all'), numeric_cols)
      if(length(count_cols) > 0){
        ref_idx = which(as.character(burden_mean$scenario) == ref_scenario_name)
        numbers_averted_df = data.frame(scenario=as.character(burden_mean$scenario), stringsAsFactors=FALSE)
        for(col in count_cols){
          counts = as.numeric(burden_mean[[col]]) * national_population / 1000
          ref_count = if(length(ref_idx) == 1) counts[ref_idx] else NA_real_
          numbers_averted_df[[col]] = round(ref_count - counts, 0)
        }
        write.csv(numbers_averted_df, paste0(output_dir, '/national_burden_all_scenarios_numbers_averted', file_suffix, '.csv'), row.names=FALSE)
      }
    }
  }

  return(invisible(burden_mean))
}
