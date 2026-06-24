#################################################################################################
# High burden, high impact 
# malariaInPregnancy_functions.R#
# Monique Ambrose
# February 2020
#
# 
# Functions for post-processing simulations to account for the impacts of malaria in pregnancy and IPTp.
#
# This script contains functions needed to prepare the data, perform calculations, and update 
#    simulation output
#
#
#################################################################################################

library(data.table)
library(lubridate)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
###                                     functions
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###

get_pop_size_infection_from_sim = function(first_year=2010, last_year=2019, admin_names, sim_output){
  #' calculate the population size in each DS and month as well as the probability of escaping infection during a month
  #' 
  #'  @param first_year first year of simulation to be included in current calculations
  #'  @param last_year last year of simulation to be included in current calculations
  #'  @param admin_names vector of names for each of the health districts simulated
  #'  @param sim_output data table of simulation results across all health districts and months from a single run
  #'  
  #' Return a list of data tables. Each data table, columns correspond to DS (same order as iptp coverage file) and rows correspond to year/month combinations:
  #'   pop_size - mean population size across all days in a month and ds
  #'   popUnder15_size - mean population size for the under 15age bin across all days in a month and ds
  #'   pop1530_size - mean population size for the 15-30 age bin across all days in a month and ds
  #'   pop3050_size- mean population size for the 30-50 age bin across all days in a month and ds
  #'   popOver50_size - mean population size for the over 50 age bin across all days in a month and ds
  #'   prob_no_infection_1530- the probability of individual in 15-30 age group escaping infection across 
  #'       all days in that month, assuming uniform risk for all individuals in each age group.
  #'   prob_no_infection_3050- the probability of individual in 30-50 age group escaping infection across 
  #'       all days in that month, assuming uniform risk for all individuals in each age group.
  #'   pfpr_under15 - parasite prevalence for each month, directly from simulation output
  #'   pfpr_1530 - parasite prevalence for each month, directly from simulation output
  #'   pfpr_3050 - parasite prevalence for each month, directly from simulation output
  #'   pfpr_over50 - parasite prevalence for each month, directly from simulation output
  
  # replace spaces with periods in column names
  colnames(sim_output) = gsub(' ', '.', colnames(sim_output))
  sim_output = setorder(sim_output, admin_name,year,month)
  
  # each row corresponds to a month in the simulation
  months = rep(1:12, times=length(first_year:last_year))
  years = rep(first_year:last_year, each=12)
  days_in_sim_month = c(rep(30,11), 35)  # each month assumed to have 30 days except December, which has the extra 5

  # initialize data tables - rows are month of simulation, columns are DS
  pop_size = data.table('year'=years,'month'=months)
  popUnder15_size = data.table('year'=years,'month'=months)
  pop1530_size = data.table('year'=years,'month'=months)
  pop3050_size = data.table('year'=years,'month'=months)
  popOver50_size = data.table('year'=years,'month'=months)
  prob_no_infection_1530 = data.table('year'=years,'month'=months)
  prob_no_infection_3050 = data.table('year'=years,'month'=months)
  pfpr_under15 = data.table('year'=years,'month'=months)
  pfpr_1530 = data.table('year'=years,'month'=months)
  pfpr_3050 = data.table('year'=years,'month'=months)
  pfpr_over50 = data.table('year'=years,'month'=months)
  
  # subset to under15, 15-30, 30-50, and over50 age groups
  sim_output_under15 = sim_output[sim_output$AgeGroup=='Under15',]
  sim_output_1530 = sim_output[sim_output$AgeGroup=='15to30',]
  sim_output_3050 = sim_output[sim_output$AgeGroup=='30to50',]
  sim_output_over50 = sim_output[sim_output$AgeGroup=='50plus',]
  
  for(i_ds in 1:length(admin_names)){
    if(i_ds %% 40 == 0){
      print(paste('currently calculating population sizes and probability of avoiding infection for ds#', i_ds, 'out of', length(admin_names)))
    }
    # subset data table to section for this DS
    ds_rows_under15 = which(toupper(sim_output_under15$admin_name) == toupper(admin_names[i_ds]))
    ds_rows_1530 = which(toupper(sim_output_1530$admin_name) == toupper(admin_names[i_ds]))
    ds_rows_3050 = which(toupper(sim_output_3050$admin_name) == toupper(admin_names[i_ds]))
    ds_rows_over50 = which(toupper(sim_output_over50$admin_name) == toupper(admin_names[i_ds]))
    ds_rows_all = which(toupper(sim_output$admin_name) == toupper(admin_names[i_ds]))
    if(length(ds_rows_1530)>0){
      sim_output_under15_ds = sim_output_under15[ds_rows_under15, ]
      sim_output_1530_ds = sim_output_1530[ds_rows_1530, ]
      sim_output_3050_ds = sim_output_3050[ds_rows_3050, ]
      sim_output_over50_ds = sim_output_over50[ds_rows_over50, ]
      sim_output_ds = sim_output[ds_rows_all, ]
      
      # calculate mean population size for all age groups (sum across age bins) within a month and year (order should be the same as months and years vectors from above)
      pop_month_full = aggregate(sim_output_ds$Statistical.Population, by=list(sim_output_ds$month, sim_output_ds$year), sum)
      # check that the order of years/months in the subsets of sim_output are consistent
      if(all(pop_month_full$Group.1 == months) & all(pop_month_full$Group.2 == years)){
        pop_month = pop_month_full$x
      } else{
        warning(paste('month and year ordering do not match what was expected from population size calculation grouping. need to add code to re-order'))
        (paste('month and year ordering do not match what was expected from population size calculation grouping. need to add code to re-order'))
        pop_month = rep(NA, length(months))
      }
      # add the column for population sizes for this DS to the data table
      pop_size[, admin_names[i_ds]:= pop_month]

      # check that the order of years/months in the subsets of sim_output are consistent
      if(all(sim_output_1530_ds$month == months) & all(sim_output_1530_ds$year == years) & all(sim_output_3050_ds$month == months) & all(sim_output_3050_ds$year == years)){
        
        # add the column for population sizes for this DS to the data tables
        popUnder15_size[, admin_names[i_ds]:= sim_output_under15_ds$Statistical.Population]
        pop1530_size[, admin_names[i_ds]:= sim_output_1530_ds$Statistical.Population]
        pop3050_size[, admin_names[i_ds]:= sim_output_3050_ds$Statistical.Population]
        popOver50_size[, admin_names[i_ds]:= sim_output_over50_ds$Statistical.Population]
        pfpr_under15[, admin_names[i_ds]:= sim_output_under15_ds$PfPR]
        pfpr_1530[, admin_names[i_ds]:= sim_output_1530_ds$PfPR]
        pfpr_3050[, admin_names[i_ds]:= sim_output_3050_ds$PfPR]
        pfpr_over50[, admin_names[i_ds]:= sim_output_over50_ds$PfPR]

        # Probability an individual avoids infectious exposure
        # Explanation of approach: 
        # We need to calculate the probability a given individual avoids infectious exposure in a month using the number of individuals who get new infections 
        #    in that month (num_new) and the population size (pop_size). However, if individuals recover quickly from earlier infections, they could 
        #    be counted as having a new infection multiple times within a month. This double counting is not likely to make much difference to the 
        #    overall probability of infection when the number of new infections is much smaller than the population size, but becomes more and more 
        #    important as the number of new infections approaches (or exceeds) the population size. Ignoring the double counting will lead to 
        #    underestimating the probability an individual gets 0 infections in the month if many people are getting multiple infections. In addition, 
        #    since we are interested in infectious exposure rather than any new infection, the denominator is really the number of people who are not 
        #    currently infected, ignoring this will lead to overestimates of the probability an individual avoids infectious exposure. 
        #      --> est_prob_avoid_infection1a = (1-num_new/pop_size)  # could be negative
        #      --> est_prob_avoid_infection1b = (1-num_new/(pop_size*(1-prev)))  # could be negative
        #      --> est_prob_avoid_infection1b <= true_prob_avoid_infection
        #    Another approach would be to look at the probability an individual avoids infection on every day of the month, by approximating 
        #    the number of new infections on each day as the average daily number within the month. However, this won't give the appropriate value 
        #    either, because it isn't the number of infectious exposures (which would give probability of avoiding infection in the month), it's the 
        #    number of infectious exposures among individuals who aren't currently infected. The denominator of what is being observed is 
        #    actually (1-prev)*pop_size rather than pop_size, so estimates will be too large.
        #      --> est_prob_avoid_infection2a = (1-num_new/30/pop_size)^30
        #      --> est_prob_avoid_infection2a >= true_prob_avoid_infection 
        #    However, since we know the prevalence during that month, we can use that to calculate the true denominator of people who could be 
        #    observed to have a new infection. pop_size_adjusted = pop_size*(1-prev). If prevalence and rate of new infections are relatively 
        #    uniform throughout the month, this should give a reasonable approximation of the probability any given individual avoids getting a 
        #    new infectious exposure during the month. 
        #      --> est_prob_avoid_infection2b = (1-num_new/30/(pop_size*(1-prev)))^30
        #      --> est_prob_avoid_infection2b ~= true_prob_avoid_infection 
        #
        # # plot to compare estimates and check logic. should see that blue dots are <= black dots and red dots are >= black dots
        # prob_no_infection_1530_1a = (1 - sim_output_1530_ds$New.Infections / sim_output_1530_ds$Statistical.Population )
        # prob_no_infection_1530_1b = (1 - sim_output_1530_ds$New.Infections / (sim_output_1530_ds$Statistical.Population * (1-sim_output_1530_ds$PfPR)) )
        # prob_no_infection_1530_2a = (1 - sim_output_1530_ds$New.Infections/ days_in_sim_month[sim_output_1530_ds$month] / sim_output_1530_ds$Statistical.Population ) ^ days_in_sim_month[sim_output_1530_ds$month]
        # prob_no_infection_1530_2b = (1 - sim_output_1530_ds$New.Infections / days_in_sim_month[sim_output_1530_ds$month] / (sim_output_1530_ds$Statistical.Population * (1-sim_output_1530_ds$PfPR))) ^ days_in_sim_month[sim_output_1530_ds$month]
        # plot(prob_no_infection_1530_2b, type='p', pch=20, col='black', ylim=c(0,1), bty='L')
        # points(prob_no_infection_1530_1a, pch=20, col='green')
        # points(prob_no_infection_1530_1b, pch=20, col='blue')
        # points(prob_no_infection_1530_2a, pch=20, col='red')
        # points(sim_output_1530_ds$PfPR, pch='-')
        
        # add the columns for this DS to the data tables. maximum probability of being infected in a day is one.
        prob_no_infection_1530[, admin_names[i_ds] := ((1 - sapply(sim_output_1530_ds$New.Infections / days_in_sim_month[sim_output_1530_ds$month] / (sim_output_1530_ds$Statistical.Population * (1-sim_output_1530_ds$PfPR)), min,1)) ^ days_in_sim_month[sim_output_1530_ds$month])]
        prob_no_infection_3050[, admin_names[i_ds] := ((1 - sapply(sim_output_3050_ds$New.Infections / days_in_sim_month[sim_output_3050_ds$month] / (sim_output_3050_ds$Statistical.Population * (1-sim_output_3050_ds$PfPR)), min,1)) ^ days_in_sim_month[sim_output_3050_ds$month])]

      } else{
        warning(paste('month and year ordering do not match what was expected. need to add code to re-order'))
        (paste('month and year ordering do not match what was expected. need to add code to re-order'))
      }
    } else{
      warning(paste('ds name not found in simulation:', admin_names[i_ds]))
      (paste('ds name not found in simulation:', admin_names[i_ds]))
    }
  }
  return(list(pop_size, popUnder15_size, pop1530_size, pop3050_size, popOver50_size, prob_no_infection_1530, prob_no_infection_3050, pfpr_under15, pfpr_1530, pfpr_3050, pfpr_over50))
}






calc_mStill_mLBW_death = function(preg_withoutNMStill_1530, preg_withoutNMStill_3050, prob_infect_in_preg_1530, prob_infect_in_preg_3050, f_iptp_t, pm_LBW_iptp, pm_still_iptp, d0, d1, frac_first_second_birth){
  #' calculate number of malaria-attributed stillbirths, livebirths, mLBWs, and mLBW deaths with the estimated levels of IPTp
  #' 
  #' @param preg_withoutNMStill_1530 number of pregnancies scheduled for delibery in month t in 15-30 age group, excluding preganancies that result in stillbirth from non-malarial causes
  #' @param preg_withoutNMStill_3050 number of pregnancies scheduled for delibery in month t in 30-50 age group, excluding preganancies that result in stillbirth from non-malarial causes
  #' @param prob_infect_in_preg_1530 fraction of individuals aged 15-30 infected in any of the nine months before t
  #' @param prob_infect_in_preg_3050 fraction of individuals aged 30-50 infected in any of the nine months before t
  #' @param f_iptp_t fraction of pregnancies where birth occurs in month t with 0, 1, 2, or 3 IPTp doses
  #' @param pm_LBW_iptp probability of a mLBW infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #' @param pm_still_iptp probability of a stillborn infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #' @param d0  probability of death for a non-LBW infant
  #' @param d1  probability of death for a LBW infant
  #' @param frac_first_second_birth fraction of all births that are of the woman's first or second child
  #' 
  #' Return vector containing:
  #'  - number of malaria-attributable stillbirths
  #'  - number of live first or second births in month t
  #'  - number of mLBW births
  #'  - number of mLBW-attributable deaths
  #'  - number of mLBW-attributable deaths if not IPTp were used
  #'  - fraction of livebirths that occur thanks to IPTp (that otherwise would have been stillbirths)

  # estimate number of stillbirths
  N_nonstill_1530_exposed = preg_withoutNMStill_1530 * prob_infect_in_preg_1530  # number of births in month t with malaria exposure (excluding stillbirths from non-malarial causes) - includes stillbirths from malaria
  N_nonstill_3050_exposed = preg_withoutNMStill_3050 * prob_infect_in_preg_3050  # number of births in month t with malaria exposure (excluding stillbirths from non-malarial causes) - includes stillbirths from malaria
  D_still_1530 = sum(N_nonstill_1530_exposed * f_iptp_t * pm_still_iptp)  # number of malaria-attributable stillbirths in month t with the estimated levels of IPTp
  D_still_3050 = sum(N_nonstill_3050_exposed * f_iptp_t * pm_still_iptp)  # number of malaria-attributable stillbirths in month t with the estimated levels of IPTp
  D_still_all = D_still_1530 + D_still_3050
  # stillbirth numbers if there were no IPTp
  D_still_1530_noIPTp = sum(N_nonstill_1530_exposed * pm_still_iptp[1])  # number of malaria-attributable stillbirths in month t with no IPTp
  D_still_3050_noIPTp = sum(N_nonstill_3050_exposed * pm_still_iptp[1])  # number of malaria-attributable stillbirths in month t with no IPTp
  D_still_all_noIPTp = D_still_1530_noIPTp + D_still_3050_noIPTp
  
  # estimate mLBW numbers
  N_1or2_livebirth = ((preg_withoutNMStill_1530 - D_still_1530) + (preg_withoutNMStill_3050 - D_still_3050)) * frac_first_second_birth  # number of livebirths from first or second pregnancies in the population (excludes stillbirths)
  N_mLBW =  sum(N_1or2_livebirth * prob_infect_in_preg_1530 * f_iptp_t * pm_LBW_iptp)  # the number of mLBW infants during month t with the estimated levels of IPTp
  D_mLBW = N_mLBW * (d1 - d0)  # expected number of deaths due to mLBW in month t with the estimated levels of IPTp
  
  # fraction of livebirths that occur thanks to IPTp (they otherwise would have been stillbirths - we use this to see whether IPTi population fractions need to be adjusted)
  frac_birth_averted_mStill = (D_still_all_noIPTp - D_still_all) / ((preg_withoutNMStill_1530 - D_still_1530) + (preg_withoutNMStill_3050 - D_still_3050))

  return(c(D_still_all, N_1or2_livebirth, N_mLBW, D_mLBW, D_still_all_noIPTp, frac_birth_averted_mStill))
}








calculate_MiP_numbers = function(pop_size, prob_no_infection_1530, prob_no_infection_3050, iptp_coverage_df, admin_names, pm_LBW_iptp, pm_still_iptp, name_for_pm_still_iptp, d0, d1, births_per_person_per_month, 
                                 non_mal_stillbirth_fraction, frac_preg_under_30, frac_first_second_birth, frac_iptp_1, frac_iptp_2, frac_iptp_3, frac_iptp_years, sim_output_base_filepath, cur_run, coverage_string, overwrite_files_flag){
  #' save data tables that contain information necessary for calculating MiP outcomes, including PfPR adjustments from IPTp, Severe case adjustment from maternal severe anemia, stillbirths from malaria, and mLBW and mLBW-deaths
  #'
  #'  @param pop_size mean population size across all days in a month and ds
  #'  @param prob_no_infection_1530 the probability of individual in 15-30 age group escaping infection across 
  #'       all days in that month, assuming uniform risk for all individuals in each age group.
  #'  @param prob_no_infection_3050 the probability of individual in 30-50 age group escaping infection across 
  #'       all days in that month, assuming uniform risk for all individuals in each age group.
  #'  @param iptp_coverage_df fraction of pregnant individuals receiving >=1 IPTp dose in each DS (row) and year (column)
  #'  @param admin_names vector of names for each of the health districts simulated
  #'  @param pm_LBW_iptp probability of a mLBW infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #'  @param pm_still_iptp probability of a stillborn infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #'  @param name_for_pm_still_iptp string tag describing how protective IPTp is against stillbirth compared to protectiveness against mLBW
  #'  @param d0  probability of death for a non-LBW infant
  #'  @param d1  probability of death for a LBW infant
  #'  @param births_per_person_per_month number of births per person per month
  #'  @param non_mal_stillbirth_fraction fraction of births that result in stillbirth from causes other than malaria
  #'  @param frac_preg_under_30 fraction of pregnancies that occur in individuals under the age of 30
  #'  @param frac_first_second_birth fraction of all births that are of the woman's first or second child
  #'  @param frac_iptp_1 fraction of individuals with >=1 IPTp dose who receive 1 dose in each year
  #'  @param frac_iptp_2 fraction of individuals with >=1 IPTp dose who receive 2 doses in each year
  #'  @param frac_iptp_3 fraction of individuals with >=1 IPTp dose who receive 3 doses in each year
  #'  @param frac_iptp_years vector giving the years corresponding to each entry of the frac_iptp_i variables
  #'  @param sim_output_base_filepath filepath where simulation output is saved
  #'  @param cur_run current run or seed of the simulation experiment that is being analyzed
  #'  @param coverage_string description of the assumed IPTp coverage
  #'  @param overwrite_files_flag boolean that describes whether or not analyses should be run and saved even if output files already exist
  #'   
  #' The following ata tables saved as csv files (each data tables  has DS as column, month of simulation as row):
  #'   - probability of infection during course of pregnancy among 15-30 age group <-- needed for maternal severe case calculations, calculating stillbirths, and calculating mLBW
  #'   - probability of infection during course of pregnancy among 30-50 age group <-- needed for calculating stillbirths
  #'   - number of pregnancies in each month (including all types of stillbirths) among 15-30 age group  <-- needed for PfPR adjustments
  #'   - number of pregnancies in each month (including all types of stillbirths) among 30-50 age group  <-- needed for PfPR adjustments
  #'   - number of first or second pregnancies in each month (including all types of stillbirths). Current assumption is that all are in 15-30 age group.  <-- needed for maternal severe case calculations
  #'   - number of pregnancies in each month (*excluding* stillbirths from non-malarial causes) among 15-30 age group  <-- needed for calculating stillbirths
  #'   - number of pregnancies in each month (*excluding* stillbirths from non-malarial causes) among 30-50 age group  <-- needed for calculating stillbirths
  #'   - number of first or second pregnancies in each month (*excluding* stillbirths from non-malarial causes). Current assumption is that all are in 15-30 age group.  <-- needed for maternal mLBW calculations
  #'   - number of stillbirths with estimated IPTp
  #'   - number of mLBW with estimated IPTp
  #'   - number of mLBW-deaths with estimated IPTp
  #'   - number of stillbirths without IPTp
  #'   - fraction of livebirths that occur thanks to IPTp (that otherwise would have been stillbirths)
  #'   - number of IPTp doses used in each DS and each month
  
  if(overwrite_files_flag | !file.exists(paste(sim_output_base_filepath, '/mortality_LBW/prob_infected_in_pregnancy_3050_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep='')) | !file.exists(paste(sim_output_base_filepath, '/mortality_LBW/mLBW_deaths_each_DS_month_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep='')) | !file.exists(paste(sim_output_base_filepath, '/mortality_LBW/num_doses_IPTp_each_DS_month',coverage_string,'_run',cur_run, '.csv', sep=''))){
    
    # set up empty data tables which will store probabilities and/or counts for each DS and month
    empty_values_table = copy(pop_size)
    set(empty_values_table, 1:dim(empty_values_table)[1], 3:dim(empty_values_table)[2], NA)
    prob_infected_in_pregnancy_1530 = copy(empty_values_table)
    prob_infected_in_pregnancy_3050 = copy(empty_values_table)
    num_preg_1530_withAllStill = copy(empty_values_table)
    num_preg_3050_withAllStill = copy(empty_values_table)
    num_1or2_preg_withAllStill = copy(empty_values_table)
    num_preg_1530_withoutNMStill = copy(empty_values_table)
    num_preg_3050_withoutNMStill = copy(empty_values_table)
    num_1or2_preg_withoutNMStill = copy(empty_values_table)
    mStill_each_DS_month = copy(empty_values_table)
    mLBW_each_DS_month = copy(empty_values_table)
    mLBW_deaths_each_DS_month = copy(empty_values_table)
    mStill_noIPTp_each_DS_month = copy(empty_values_table)
    frac_birth_averted_mStill_each_DS_month = copy(empty_values_table)
    num_doses_IPTp_each_DS_month = copy(empty_values_table)
    
    # iterate though DS and month, calculating number of mLBW with no IPTp and with estimated IPTp coverage
    for(i_ds in 1:length(admin_names)){
      if(i_ds %% 40 == 0){
        print(paste('currently calculating number of pregnancies and probability of infection for ds#', i_ds, 'out of', length(admin_names)))
      }
      
      # get the column number for this ds for the simulation output data tables (population size and probability of avoiding infection)
      sim_col_num = which(colnames(pop_size) == admin_names[i_ds])
      
      # get the row number for the IPTp coverage data frame
      iptp_row_num = which(iptp_coverage_df$admin_name == admin_names[i_ds])
      
      # iterate through months in the simulation
      for(i_m in 13:dim(prob_no_infection_1530)[1]){  # start at 2nd year because don't have full simulated infection information on pregnancies for births before the 9th month of the 1st year
        # find the iptp coverage for the current year
        # get the current year
        cur_year = prob_no_infection_1530$year[i_m]
        # get the column corresponding to this year in the iptp coverage dataframe
        iptp_col_num = grep(cur_year, colnames(iptp_coverage_df))
        # get the IPTp coverage for this DS and this year
        iptp_coverage_cur = iptp_coverage_df[iptp_row_num, iptp_col_num]
  
        # determine what fraction of individuals get each number of IPTp doses (0->3)
        # which ratio of doses corresponds to this year
        iptp_dose_year_index = which(frac_iptp_years == cur_year)
        f_iptp_t = c((1-iptp_coverage_cur), iptp_coverage_cur * frac_iptp_1[iptp_dose_year_index], iptp_coverage_cur * frac_iptp_2[iptp_dose_year_index], iptp_coverage_cur * frac_iptp_3[iptp_dose_year_index])
        
        # calculate the fraction of individuals giving birth in this month that were infected with malaria (or would have been without IPTp)
        prob_infect_in_preg_1530 = 1 - prod(prob_no_infection_1530[(i_m-8):i_m, ..sim_col_num])
        prob_infect_in_preg_3050 = 1 - prod(prob_no_infection_3050[(i_m-8):i_m, ..sim_col_num])
        
        # total pregnancies in each age group, total pregnancies excluding stillbirths from non-malarial causes, and total first-and-second pregnancies
        total_preg_withAllStill = as.numeric(pop_size[i_m, ..sim_col_num])*births_per_person_per_month
        preg_withAllStill_1530 = total_preg_withAllStill * frac_preg_under_30
        preg_withAllStill_3050 = total_preg_withAllStill * (1-frac_preg_under_30)
        total_preg_withoutNMStill = total_preg_withAllStill * (1-non_mal_stillbirth_fraction)
        preg_withoutNMStill_1530 = total_preg_withoutNMStill * frac_preg_under_30
        preg_withoutNMStill_3050 = total_preg_withoutNMStill * (1 - frac_preg_under_30)
        
        # calculate the stillbirth, mLBW, and mLBW mortality values for this month and DS
        mStill_mLBW_ds_m = calc_mStill_mLBW_death(preg_withoutNMStill_1530=preg_withoutNMStill_1530, preg_withoutNMStill_3050=preg_withoutNMStill_3050, prob_infect_in_preg_1530=prob_infect_in_preg_1530, prob_infect_in_preg_3050=prob_infect_in_preg_3050, f_iptp_t=f_iptp_t, pm_LBW_iptp=pm_LBW_iptp, pm_still_iptp=pm_still_iptp, d0=d0, d1=d1, frac_first_second_birth=frac_first_second_birth)
        
        # calculate number of IPTp doses used for this DS and month
        num_IPTp = sum((total_preg_withAllStill * f_iptp_t) * c(0, 1, 2, 3))
          
        # store values in data frames for this DS and month
        set(prob_infected_in_pregnancy_1530, i_m, sim_col_num, prob_infect_in_preg_1530)
        set(prob_infected_in_pregnancy_3050, i_m, sim_col_num, prob_infect_in_preg_3050)
        set(num_preg_1530_withAllStill, i_m, sim_col_num, preg_withAllStill_1530)
        set(num_preg_3050_withAllStill, i_m, sim_col_num, preg_withAllStill_3050)
        set(num_1or2_preg_withAllStill, i_m, sim_col_num, total_preg_withAllStill*frac_first_second_birth)  # assumes the fraction of first and second pregnancies is the same as the fraction of first and second births (this may not be quite right if stillbirths occur more in earlier or later pregnancies)
        set(num_preg_1530_withoutNMStill, i_m, sim_col_num, preg_withoutNMStill_1530)
        set(num_preg_3050_withoutNMStill, i_m, sim_col_num, preg_withoutNMStill_3050)
        set(num_1or2_preg_withoutNMStill, i_m, sim_col_num, total_preg_withoutNMStill*frac_first_second_birth)
        set(mStill_each_DS_month, i_m, sim_col_num, mStill_mLBW_ds_m[1])
        set(mLBW_each_DS_month, i_m, sim_col_num, mStill_mLBW_ds_m[3])
        set(mLBW_deaths_each_DS_month, i_m, sim_col_num, mStill_mLBW_ds_m[4])
        set(mStill_noIPTp_each_DS_month, i_m, sim_col_num, mStill_mLBW_ds_m[5])
        set(frac_birth_averted_mStill_each_DS_month, i_m, sim_col_num, mStill_mLBW_ds_m[6])
        set(num_doses_IPTp_each_DS_month, i_m, sim_col_num, num_IPTp)
      }
    }
    # save data tables
    fwrite(prob_infected_in_pregnancy_1530, paste(sim_output_base_filepath, '/mortality_LBW/prob_infected_in_pregnancy_1530_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(prob_infected_in_pregnancy_3050, paste(sim_output_base_filepath, '/mortality_LBW/prob_infected_in_pregnancy_3050_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_preg_1530_withAllStill, paste(sim_output_base_filepath, '/mortality_LBW/num_preg_1530_withAllStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_preg_3050_withAllStill, paste(sim_output_base_filepath, '/mortality_LBW/num_preg_3050_withAllStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_1or2_preg_withAllStill, paste(sim_output_base_filepath, '/mortality_LBW/num_1or2_preg_withAllStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_preg_1530_withoutNMStill, paste(sim_output_base_filepath, '/mortality_LBW/num_preg_1530_withoutNMStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_preg_3050_withoutNMStill, paste(sim_output_base_filepath, '/mortality_LBW/num_preg_3050_withoutNMStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_1or2_preg_withoutNMStill, paste(sim_output_base_filepath, '/mortality_LBW/num_1or2_preg_withoutNMStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(mStill_each_DS_month, paste(sim_output_base_filepath, '/mortality_LBW/mStill_each_DS_month',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(mLBW_each_DS_month, paste(sim_output_base_filepath, '/mortality_LBW/mLBW_each_DS_month_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(mLBW_deaths_each_DS_month, paste(sim_output_base_filepath, '/mortality_LBW/mLBW_deaths_each_DS_month_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(mStill_noIPTp_each_DS_month, paste(sim_output_base_filepath, '/mortality_LBW/mStill_noIPTp_each_DS_month',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(frac_birth_averted_mStill_each_DS_month, paste(sim_output_base_filepath, '/mortality_LBW/frac_birth_averted_mStill_each_DS_month',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    fwrite(num_doses_IPTp_each_DS_month, paste(sim_output_base_filepath, '/mortality_LBW/num_doses_IPTp_each_DS_month',coverage_string,'_run',cur_run, '.csv', sep=''))
  }
}






calculate_pregnancy_values_across_runs = function(sim_output_all_runs, first_year, last_year, admin_names, iptp_coverage_df, pm_LBW_iptp, list_of_pm_still_iptp, names_for_pm_still_iptp, d0, d1, 
                                                  births_per_person_per_month, non_mal_stillbirth_fraction, frac_first_second_birth, frac_preg_under_30, frac_iptp_1, frac_iptp_2, frac_iptp_3, frac_iptp_years, coverage_string,
                                                  sim_output_base_filepath, sim_output_base_filepath_2010_allInter, future_projection_flag, admin_subset_flag=FALSE, admin_subset=NA, overwrite_files_flag=FALSE){
  #' calculate MiP related counts and probabilities for each run in simulation output. Results are saved as csv files.
  #'
  #'  @param sim_output_all_runs data table of simulation results across all health districts and months from all runs
  #'  @param first_year first year of simulation to be included in current calculations
  #'  @param last_year last year of simulation to be included in current calculations
  #'  @param admin_names vector of names for each of the health districts simulated
  #'  @param iptp_coverage_df fraction of pregnant individuals receiving >=1 IPTp dose in each DS (row) and year (column)
  #'  @param pm_LBW_iptp probability of a mLBW infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #'  @param list_of_pm_still_iptp list of probabilities of a stillborn infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #'  @param names_for_pm_still_iptp vector of string tags describing how protective IPTp is against stillbirth compared to protectiveness against mLBW  (correspond to entries in list_of_pm_still_iptp)
  #'  @param d0  probability of death for a non-LBW infant
  #'  @param d1  probability of death for a LBW infant
  #'  @param births_per_person_per_month number of births per person per month
  #'  @param non_mal_stillbirth_fraction fraction of births that result in stillbirth from causes other than malaria
  #'  @param frac_first_second_birth fraction of all births that are of the woman's first or second child
  #'  @param frac_preg_under_30 fraction of pregnancies that occur in individuals under the age of 30
  #'  @param frac_iptp_1 fraction of individuals with >=1 IPTp dose who receive 1 dose in each year
  #'  @param frac_iptp_2 fraction of individuals with >=1 IPTp dose who receive 2 doses in each year
  #'  @param frac_iptp_3 fraction of individuals with >=1 IPTp dose who receive 3 doses in each year
  #'  @param frac_iptp_years vector giving the years corresponding to each entry of the frac_iptp_i variables
  #'  @param coverage_string description of the assumed IPTp coverage
  #'  @param sim_output_base_filepath filepath where simulation output is saved
  #'  @param sim_output_base_filepath_2010_allInter filepath where simulation output from the 2010-2020 simulations is saved
  #'  @param future_projection_flag indicates whether or not the simulation is a projection of future scenarios
  #'  @param admin_subset_flag boolean indicating whether admins should be subsetted to areas that received IRS or IPTi (for simulations run with a subset of admins)
  #'  @param admin_subset vector of admin names included in the irs-only simulation. Must be specified if admin_subset_flag is TRUE
  #'  @param overwrite_files_flag boolean that describes whether or not analyses should be run and saved even if output files already exist


  all_runs = unique(sim_output_all_runs$Run_Number)
  # iterate through runs, saving information about population sizes, number of pregnancies, probability of infection during pregnancy, number of stillbirths, and number of mLBW births and mLBW deaths
  for(rr in 1:length(all_runs)){
    print(paste('currently on run index', rr, 'out of', length(all_runs)))
    
    cur_run = all_runs[rr]
    # subset data table of simulation outputs to current run
    sim_output = sim_output_all_runs[sim_output_all_runs$Run_Number == cur_run,]
    
    # get data tables with monthly simulation population sizes and probabilities of avoiding infection
    # check whether population files already exist
    if(overwrite_files_flag | !file.exists(paste(sim_output_base_filepath, '/mortality_LBW/pop_size_run',cur_run, '.csv', sep='')) | !file.exists(paste(sim_output_base_filepath, '/mortality_LBW/pop3050_size_run',cur_run, '.csv', sep=''))){
      monthly_sim_output = get_pop_size_infection_from_sim(first_year=first_year, last_year=last_year, admin_names=admin_names, sim_output=sim_output)
      pop_size = monthly_sim_output[[1]]
      popUnder15_size = monthly_sim_output[[2]]
      pop1530_size = monthly_sim_output[[3]]
      pop3050_size = monthly_sim_output[[4]]
      popOver50_size = monthly_sim_output[[5]]
      prob_no_infection_1530 = monthly_sim_output[[6]]
      prob_no_infection_3050 = monthly_sim_output[[7]]
      pfpr_under15 = monthly_sim_output[[8]]
      pfpr_1530 = monthly_sim_output[[9]]
      pfpr_3050 = monthly_sim_output[[10]]
      pfpr_over50 = monthly_sim_output[[11]]
      
      # if running future projections AND the monthly probability of avoiding malaria infection were already calculated for the last nine months of the 'past' simulations,
      #    grab the last nine months from the 'past' simulations so that the probability of being infected during teh nine-month pregnancy can be calculated for the 
      #    future simulations in year 2020
      if(future_projection_flag & file.exists(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/prob_no_infection_3050_run',cur_run, '.csv', sep='')) & file.exists(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/prob_no_infection_1530_run',cur_run, '.csv', sep=''))){
        # read in the population sizes and probabilities of not being infected in each month from the 2010-2019 simulations (the 'all interventions' scenario)
        pop_size_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pop_size_run',cur_run, '.csv', sep=''))
        popUnder15_size_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/popUnder15_size_run',cur_run, '.csv', sep=''))
        pop1530_size_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pop1530_size_run',cur_run, '.csv', sep=''))
        pop3050_size_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pop3050_size_run',cur_run, '.csv', sep=''))
        popOver50_size_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/popOver50_size_run',cur_run, '.csv', sep=''))
        prob_no_infection_1530_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/prob_no_infection_1530_run',cur_run, '.csv', sep=''))
        prob_no_infection_3050_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/prob_no_infection_3050_run',cur_run, '.csv', sep=''))
        pfpr_under15_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pfpr_under15_run',cur_run, '.csv', sep=''))
        pfpr_1530_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pfpr_1530_run',cur_run, '.csv', sep=''))
        pfpr_3050_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pfpr_3050_run',cur_run, '.csv', sep=''))
        pfpr_over50_previous = fread(paste(sim_output_base_filepath_2010_allInter, '/mortality_LBW/pfpr_over50_run',cur_run, '.csv', sep=''))
        
        # copy the rows from the year 2019 from the 2010-2019 simulations as the first 12 rows of the 2020-2025 simulations
          rows_previous_year = which(pop_size_previous$year == (first_year-1))
        if(admin_subset_flag){
          cols_previous_admin_subset = which(colnames(pop_size_previous) %in% c('year','month',admin_subset))
          pop_size = rbindlist(list(copy(pop_size_previous[rows_previous_year, ..cols_previous_admin_subset]), pop_size))
          popUnder15_size = rbindlist(list(copy(popUnder15_size_previous[rows_previous_year, ..cols_previous_admin_subset]), popUnder15_size))
          pop1530_size = rbindlist(list(copy(pop1530_size_previous[rows_previous_year, ..cols_previous_admin_subset]), pop1530_size))
          pop3050_size = rbindlist(list(copy(pop3050_size_previous[rows_previous_year, ..cols_previous_admin_subset]), pop3050_size))
          popOver50_size = rbindlist(list(copy(popOver50_size_previous[rows_previous_year, ..cols_previous_admin_subset]), popOver50_size))
          prob_no_infection_1530 = rbindlist(list(copy(prob_no_infection_1530_previous[rows_previous_year, ..cols_previous_admin_subset]), prob_no_infection_1530))
          prob_no_infection_3050 = rbindlist(list(copy(prob_no_infection_3050_previous[rows_previous_year, ..cols_previous_admin_subset]), prob_no_infection_3050))
          pfpr_under15 = rbindlist(list(copy(pfpr_under15_previous[rows_previous_year, ..cols_previous_admin_subset]), pfpr_under15))
          pfpr_1530 = rbindlist(list(copy(pfpr_1530_previous[rows_previous_year, ..cols_previous_admin_subset]), pfpr_1530))
          pfpr_3050 = rbindlist(list(copy(pfpr_3050_previous[rows_previous_year, ..cols_previous_admin_subset]), pfpr_3050))
          pfpr_over50 = rbindlist(list(copy(pfpr_over50_previous[rows_previous_year, ..cols_previous_admin_subset]), pfpr_over50))
        } else{
          pop_size = rbindlist(list(copy(pop_size_previous[rows_previous_year,]), pop_size))
          popUnder15_size = rbindlist(list(copy(popUnder15_size_previous[rows_previous_year,]), popUnder15_size))
          pop1530_size = rbindlist(list(copy(pop1530_size_previous[rows_previous_year,]), pop1530_size))
          pop3050_size = rbindlist(list(copy(pop3050_size_previous[rows_previous_year,]), pop3050_size))
          popOver50_size = rbindlist(list(copy(popOver50_size_previous[rows_previous_year,]), popOver50_size))
          prob_no_infection_1530 = rbindlist(list(copy(prob_no_infection_1530_previous[rows_previous_year,]), prob_no_infection_1530))
          prob_no_infection_3050 = rbindlist(list(copy(prob_no_infection_3050_previous[rows_previous_year,]), prob_no_infection_3050))
          pfpr_under15 = rbindlist(list(copy(pfpr_under15_previous[rows_previous_year,]), pfpr_under15))
          pfpr_1530 = rbindlist(list(copy(pfpr_1530_previous[rows_previous_year,]), pfpr_1530))
          pfpr_3050 = rbindlist(list(copy(pfpr_3050_previous[rows_previous_year,]), pfpr_3050))
          pfpr_over50 = rbindlist(list(copy(pfpr_over50_previous[rows_previous_year,]), pfpr_over50))
        }

        
        first_year_plotted = first_year
      } else{
        first_year_plotted = first_year + 1  # only start plotting a year after the first simulated year because we don't observe infections during pregnancies for the first nine birth months of the year
      }
      
      # save data tables
      fwrite(pop_size, paste(sim_output_base_filepath, '/mortality_LBW/pop_size_run',cur_run, '.csv', sep=''))
      fwrite(popUnder15_size, paste(sim_output_base_filepath, '/mortality_LBW/popUnder15_size_run',cur_run, '.csv', sep=''))
      fwrite(pop1530_size, paste(sim_output_base_filepath, '/mortality_LBW/pop1530_size_run',cur_run, '.csv', sep=''))
      fwrite(pop3050_size, paste(sim_output_base_filepath, '/mortality_LBW/pop3050_size_run',cur_run, '.csv', sep=''))
      fwrite(popOver50_size, paste(sim_output_base_filepath, '/mortality_LBW/popOver50_size_run',cur_run, '.csv', sep=''))
      fwrite(prob_no_infection_1530, paste(sim_output_base_filepath, '/mortality_LBW/prob_no_infection_1530_run',cur_run, '.csv', sep=''))
      fwrite(prob_no_infection_3050, paste(sim_output_base_filepath, '/mortality_LBW/prob_no_infection_3050_run',cur_run, '.csv', sep=''))
      fwrite(pfpr_under15, paste(sim_output_base_filepath, '/mortality_LBW/pfpr_under15_run',cur_run, '.csv', sep=''))
      fwrite(pfpr_1530, paste(sim_output_base_filepath, '/mortality_LBW/pfpr_1530_run',cur_run, '.csv', sep=''))
      fwrite(pfpr_3050, paste(sim_output_base_filepath, '/mortality_LBW/pfpr_3050_run',cur_run, '.csv', sep=''))
      fwrite(pfpr_over50, paste(sim_output_base_filepath, '/mortality_LBW/pfpr_over50_run',cur_run, '.csv', sep=''))
      
    } else{  
      # read in existing data tables
      pop_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/pop_size_run',cur_run, '.csv', sep=''))
      prob_no_infection_1530 = fread(paste(sim_output_base_filepath, '/mortality_LBW/prob_no_infection_1530_run',cur_run, '.csv', sep=''))
      prob_no_infection_3050 = fread(paste(sim_output_base_filepath, '/mortality_LBW/prob_no_infection_3050_run',cur_run, '.csv', sep=''))
    }

    # iterate through IPTp protectiveness against stillbirth list, storing data table csv files with MiP outcomes (number of pregnancies, probability of infection, adverse outcomes, etc.)
    for (i_still in 1:length(list_of_pm_still_iptp)){
      print(paste('Starting calculations for number of pregnancies, probability of infection, and number of MiP outcomes for i_still', i_still, 'out of', length(list_of_pm_still_iptp)))
      # convert monthly simulation information into counts of mLBW, mortality, etc
      calculate_MiP_numbers(pop_size, prob_no_infection_1530, prob_no_infection_3050, iptp_coverage_df, admin_names, pm_LBW_iptp, 
                            pm_still_iptp=list_of_pm_still_iptp[[i_still]], name_for_pm_still_iptp=names_for_pm_still_iptp[i_still], 
                            d0, d1, births_per_person_per_month, non_mal_stillbirth_fraction, frac_preg_under_30, frac_first_second_birth, 
                            frac_iptp_1, frac_iptp_2, frac_iptp_3, frac_iptp_years, sim_output_base_filepath, cur_run, coverage_string, overwrite_files_flag)
    }
  }
}





convert_IPTp_coverage_format = function(iptp_coverage_df, pop_size){
  #' reformat IPTp coverage dataframe into a data table with the same rows and columns as MiP output files.
  #' 
  #'  @param iptp_coverage_df fraction of pregnant individuals receiving >=1 IPTp dose in each DS (row) and year (column)  # set up empty data frame to store IPTp coverages in same format as other MiP output files (rows=month, col=DS)
  #'  @param pop_size mean population size across all days in a month and ds

  iptp_coverage_table = copy(pop_size)
  set(iptp_coverage_table, 1:dim(iptp_coverage_table)[1], 3:dim(iptp_coverage_table)[2], NA)

  admin_names = colnames(iptp_coverage_table)[3:dim(iptp_coverage_table)[2]]
  years = unique(iptp_coverage_table$year)
  
  # iterate through DS
  for(col_index in 3:dim(iptp_coverage_table)[2]){
    cur_admin_name = colnames(iptp_coverage_table)[col_index]
    iptp_df_ds_row_num = which(toupper(iptp_coverage_df$admin_name) == toupper(cur_admin_name))
    # iterate through years, finding appropriate IPTp coverage value for the year and assigning to all months
    for(i_y in years){
      # find column in iptp_coverage_df for this year
      iptp_df_year_col_num = grep(i_y, colnames(iptp_coverage_df))
      
      # get IPTp coverage for this DS and year
      if((length(iptp_df_ds_row_num) == 1) & (length(iptp_df_year_col_num) == 1)){
        iptp_cov_cur = iptp_coverage_df[iptp_df_ds_row_num, iptp_df_year_col_num]
      }else{
        iptp_cov_cur = NA
        # for the 2020-2025 simulations, the 2019 values will be included as empty/blank rows when the 2010-2020 simulations have been used
        #    to get 2019 MiP values. No need to print out an error statement for not finding them in the 2020-2025 simulation output.
        if(i_y != years[1]) warning(paste(cur_admin_name, 'or year', i_y, 'not found in IPTp dataframe'))  # TODO: change years[1] to be the year before the simulation start year (was 2019, but had to generalize for simulations starting in 2021)
      }
      
      # find rows where this IPTp coverage should be inserted in the iptp_coverage_table
      iptp_table_rows_cur = which(iptp_coverage_table$year == i_y)
      # insert IPTp coverage in appropriate DS and for all months of the current year
      set(iptp_coverage_table, iptp_table_rows_cur, col_index, iptp_cov_cur)
    }
  } 
  return(iptp_coverage_table)
}





adjust_sim_output_for_MiP = function(prob_severe_from_MiP=0.57,  # probability of severe anemia attributable to MiP for first&second pregnancies
                                     reduced_prob_severe_IPTp=0.61, # multiplies prob_severe_from_MiP to give probability of severe malaria with IPTp
                                     prob_severe_MiP_treated=0.5,  # probability a woman with severe disease from MiP gets treatment
                                     weeks_protection_each_dose=10,  # number of weeks of protection from new infection from one IPTp dose
                                     frac_iptp_1,
                                     frac_iptp_2,
                                     frac_iptp_3,
                                     iptp_coverage_df,
                                     pm_still_iptp,
                                     name_for_pm_still_iptp,
                                     coverage_string,
                                     sim_output_base_filepath,
                                     sim_output_allAgeMonthly_filename_cur,
                                     sim_output_allAgeMonthly_filename_MiP,
                                     last_year
                                     ){
  #' read in simulation output and calculate updated numbers taking effects of malaria-in-pregnancy (MiP) and IPTp usage into account
  #' 
  #'  @param prob_severe_from_MiP=0.57 probability of severe anemia attributable to MiP for first&second pregnancies
  #'  @param reduced_prob_severe_IPTp=0.61 multiplies prob_severe_from_MiP to give probability of severe malaria with IPTp
  #'  @param prob_severe_MiP_treated=0.5 probability a woman with severe disease from MiP gets treatment
  #'  @param weeks_protection_each_dose=10 number of weeks of protection from new infection from one IPTp dose
  #'  @param frac_iptp_1 fraction of individuals with >=1 IPTp dose who receive 1 dose in each year
  #'  @param frac_iptp_2 fraction of individuals with >=1 IPTp dose who receive 2 doses in each year
  #'  @param frac_iptp_3 fraction of individuals with >=1 IPTp dose who receive 3 doses in each year
  #'  @param iptp_coverage_df fraction of pregnant individuals receiving >=1 IPTp dose in each DS (row) and year (column)
  #'  @param pm_still_iptp probability of a stillborn infant given that the mother was infected with malaria in her first or second pregnancy - 0, 1, 2, or 3 IPTp doses
  #'  @param name_for_pm_still_iptp string tag describing how protective IPTp is against stillbirth compared to protectiveness against mLBW
  #'  @param coverage_string description of the assumed IPTp coverage
  #'  @param sim_output_base_filepath filepath where simulation output is saved
  #'  @param sim_output_allAgeMonthly_filename_cur name for data frame to be modified with MiP and IPTp adjustments
  #'  @param sim_output_allAgeMonthly_filename_MiP name to use for saving modified data frame
  #'  @param last_year final year included in simulation
  
  # read in simulation files that will be adjusted
  sim_output_allAgeMonthly = fread(sim_output_allAgeMonthly_filename_cur)
  # replace spaces with periods in column names
  colnames(sim_output_allAgeMonthly) = gsub(' ', '.', colnames(sim_output_allAgeMonthly))
  
  if('LGA' %in% colnames(sim_output_allAgeMonthly)){
    sim_ds_colname = 'LGA'
  } else{
    sim_ds_colname = 'admin_name'
  }
  
  # remove any values from after the final included year in simulations
  sim_output_allAgeMonthly = sim_output_allAgeMonthly[lubridate::year(sim_output_allAgeMonthly$date) <= last_year,]

  # create new data tables with additional columns
  empty_column_values = as.numeric(rep(NA, dim(sim_output_allAgeMonthly)[1]))
  adjusted_allAgeMonthly = cbind(sim_output_allAgeMonthly, 
                                 data.table(PfPR_MiP_adjusted=empty_column_values,
                                 severe_maternal=empty_column_values,
                                 severe_total=empty_column_values,
                                 severe_maternal_treated=empty_column_values,
                                 severe_total_treated=empty_column_values,
                                 mLBW_births=empty_column_values,
                                 mLBW_deaths=empty_column_values,
                                 MiP_stillbirths=empty_column_values,
                                 MiP_stillbirths_noIPTp=empty_column_values
                                 ))

  # add columns for the year and month
  setDT(adjusted_allAgeMonthly)[, year := format(as.Date(date), "%Y") ]
  setDT(adjusted_allAgeMonthly)[, month := format(as.Date(date), "%m") ]

  
  all_runs = unique(sim_output_allAgeMonthly$Run_Number)
  
  # iterate through runs, adding entries for new columns in the relevant run & month & ds row
  for(rr in 1:length(all_runs)){
    print(paste('currently adding new MiP column entries for run index', rr))
    
    cur_run = all_runs[rr]
  
    # read in pregnancy-related values for this simulation
    pop_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/pop_size_run',cur_run, '.csv', sep=''))
    popUnder15_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/popUnder15_size_run',cur_run, '.csv', sep=''))
    pop1530_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/pop1530_size_run',cur_run, '.csv', sep=''))
    pop3050_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/pop3050_size_run',cur_run, '.csv', sep=''))
    popOver50_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/popOver50_size_run',cur_run, '.csv', sep=''))
    pfpr_under15 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_under15_run',cur_run, '.csv', sep=''))
    pfpr_1530 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_1530_run',cur_run, '.csv', sep=''))
    pfpr_3050 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_3050_run',cur_run, '.csv', sep=''))
    pfpr_over50 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_over50_run',cur_run, '.csv', sep=''))
    prob_infected_in_pregnancy_1530 = fread(paste(sim_output_base_filepath, '/mortality_LBW/prob_infected_in_pregnancy_1530_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    prob_infected_in_pregnancy_3050 = fread(paste(sim_output_base_filepath, '/mortality_LBW/prob_infected_in_pregnancy_3050_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    num_preg_1530_withAllStill = fread(paste(sim_output_base_filepath, '/mortality_LBW/num_preg_1530_withAllStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    num_preg_3050_withAllStill = fread(paste(sim_output_base_filepath, '/mortality_LBW/num_preg_3050_withAllStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    num_1or2_preg_withAllStill = fread(paste(sim_output_base_filepath, '/mortality_LBW/num_1or2_preg_withAllStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    mStill_each_DS_month = fread(paste(sim_output_base_filepath, '/mortality_LBW/mStill_each_DS_month',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    mLBW_each_DS_month = fread(paste(sim_output_base_filepath, '/mortality_LBW/mLBW_each_DS_month_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    mLBW_deaths_each_DS_month = fread(paste(sim_output_base_filepath, '/mortality_LBW/mLBW_deaths_each_DS_month_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    mStill_noIPTp_each_DS_month = fread(paste(sim_output_base_filepath, '/mortality_LBW/mStill_noIPTp_each_DS_month',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    # prob_no_infection_1530 = fread(paste(sim_output_base_filepath, '/mortality_LBW/prob_no_infection_1530_run',cur_run, '.csv', sep=''))
    # prob_no_infection_3050 = fread(paste(sim_output_base_filepath, '/mortality_LBW/prob_no_infection_3050_run',cur_run, '.csv', sep='')) 
    # num_preg_1530_withoutNMStill = fread(paste(sim_output_base_filepath, '/mortality_LBW/num_preg_1530_withoutNMStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    # num_preg_3050_withoutNMStill = fread(paste(sim_output_base_filepath, '/mortality_LBW/num_preg_3050_withoutNMStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))
    # num_1or2_preg_withoutNMStill = fread(paste(sim_output_base_filepath, '/mortality_LBW/num_1or2_preg_withoutNMStill_',coverage_string,'_run',cur_run, '_', name_for_pm_still_iptp, '.csv', sep=''))

    # get the IPTp coverage for each DS and month in the same format as the MiP files (rows=month, col=DS)
    iptp_coverage_table = convert_IPTp_coverage_format(iptp_coverage_df, pop_size)
    
    # sometimes, the data tables of MiP outcomes will include one additional early year that is not included in sim_output_allAgeMonthly.
    #   (this occurs when we append the 2019 simulation results from the 2010-2019 simulations so that we can calculate MiP probabilities 
    #   in 2020 even for pregnancies that occurred in 2019.)
    # check whether this is the case for the current set of simulations, and, if so, remove the first year from the data tables.
    # note: sometimes the final month recorded in adjusted_allAgeMonthly is the first month of the year after the data tables end; this is okay.
    sorted_sim_years = sort(unique(as.integer(adjusted_allAgeMonthly$year)))
    if(!(pop_size$year[1] %in% sorted_sim_years) & all(sort(unique(pop_size$year[13:length(pop_size$year)])) == sorted_sim_years[1:(length(unique(pop_size$year[13:length(pop_size$year)])))])){
        # remove the first year (first 12 rows)
        pop_size = pop_size[-(1:12),]
        popUnder15_size = popUnder15_size[-(1:12),]
        pop1530_size = pop1530_size[-(1:12),]
        pop3050_size = pop3050_size[-(1:12),]
        popOver50_size = popOver50_size[-(1:12),]
        pfpr_under15 = pfpr_under15[-(1:12),]
        pfpr_1530 = pfpr_1530[-(1:12),]
        pfpr_3050 = pfpr_3050[-(1:12),]
        pfpr_over50 = pfpr_over50[-(1:12),]
        prob_infected_in_pregnancy_1530 = prob_infected_in_pregnancy_1530[-(1:12),]
        prob_infected_in_pregnancy_3050 = prob_infected_in_pregnancy_3050[-(1:12),]
        num_preg_1530_withAllStill = num_preg_1530_withAllStill[-(1:12),]
        num_preg_3050_withAllStill = num_preg_3050_withAllStill[-(1:12),]
        num_1or2_preg_withAllStill = num_1or2_preg_withAllStill[-(1:12),]
        mStill_each_DS_month = mStill_each_DS_month[-(1:12),]
        mLBW_each_DS_month = mLBW_each_DS_month[-(1:12),]
        mLBW_deaths_each_DS_month = mLBW_deaths_each_DS_month[-(1:12),]
        mStill_noIPTp_each_DS_month = mStill_noIPTp_each_DS_month[-(1:12),]
        iptp_coverage_table = iptp_coverage_table[-(1:12),]
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # adjustment of PfPR for IPTp effects
    # (occurs in 15to30 and 30to50 age groups)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    # given high number of new infections in simulations, I think a reasonable approach is to reduce PfPR for individuals actively 
    #   protected by IPTp and then assume they are infected again once protection ends rather than calculating an additional delay 
    #   until new infection
    
    # at any given point in time, what fraction of individuals in a particular age group are expected to be protected due to IPTp usage?
    
    # assume 10 weeks complete protection for each dose, in line with SP being found to be protective for 8 weeks postpartun in Mozambique study
    #   if only one dose given, assume 10 weeks protection, if 2 doses (assume within 1 month), assume 10+4 weeks, for >=3 doses assume 10+4+4 weeks
    # default value: weeks_protection_each_dose = 10
    weeks_protection_IPTp1 = weeks_protection_each_dose
    weeks_protection_IPTp2 = weeks_protection_each_dose + 4
    weeks_protection_IPTp3 = weeks_protection_each_dose + 8
    frac_time_unprotected_IPTp1 = 1-weeks_protection_IPTp1/(9*4)
    frac_time_unprotected_IPTp2 = 1-weeks_protection_IPTp2/(9*4)
    frac_time_unprotected_IPTp3 = 1-weeks_protection_IPTp3/(9*4)
    
    # fraction of people in 15-30 age group pregnant
    # assuming a fairly constant birth rate through time (which seems reasonable given how all upstream calculations are done), and assuming 
    #    pregnancies last for nine months, the number of pregnancies during a given month will be approximately the number of pregnancies 
    #    finishing at that month times nine.
    num_ongoing_preg_1530 = num_preg_1530_withAllStill * 9
    frac_preg_1530 = num_ongoing_preg_1530 / pop1530_size
    set(frac_preg_1530, 1:dim(pop1530_size)[1], 1:2, pop1530_size[,1:2])
    # fraction of people in 30-50 age group pregnant
    num_ongoing_preg_3050 = num_preg_3050_withAllStill * 9
    frac_preg_3050 = num_ongoing_preg_3050 / pop3050_size
    set(frac_preg_3050, 1:dim(pop3050_size)[1], 1:2, pop3050_size[,1:2])
    
    
    # number of births (including stillbirths) in each month
    num_all_births = num_preg_1530_withAllStill + num_preg_3050_withAllStill
    set(num_all_births, 1:dim(num_all_births)[1], 1:2, pop1530_size[,1:2])

    # at any given point in time, what fraction of people protected? depends on:
    #  - what fraction of individuals in the age group are pregnant at any given time
    #  - what fraction of them get each number of doses
    #  - what fraction of time during their pregnancy is protected
    # pfpr = fraction_not_pregnant * pfpr_sim + 
    #        fraction_pregnant * (doses_0 * pfpr_sim * 1 + 
    #                             doses_1 * pfpr_sim * frac_time_unprotected_IPTp1/(9*4) +
    #                             doses_2 * pfpr_sim * frac_time_unprotected_IPTp2/(9*4) +
    #                             doses_3 * pfpr_sim * frac_time_unprotected_IPTp3/(9*4))
    
    frac_iptp_1_month = rep(frac_iptp_1, each=12)  # vectors of the fraction of individuals each number of IPTp doses (among people with >=1 IPTp dose) for each month.  
    frac_iptp_2_month = rep(frac_iptp_2, each=12)
    frac_iptp_3_month = rep(frac_iptp_3, each=12)
    pfpr_1530_withIPTp = (1-frac_preg_1530) * pfpr_1530 + 
      frac_preg_1530 * pfpr_1530 * ((1 - iptp_coverage_table) + 
                                      iptp_coverage_table * frac_iptp_1_month * frac_time_unprotected_IPTp1 + 
                                      iptp_coverage_table * frac_iptp_2_month * frac_time_unprotected_IPTp2 + 
                                      iptp_coverage_table * frac_iptp_3_month * frac_time_unprotected_IPTp3)
    set(pfpr_1530_withIPTp, 1:dim(pop1530_size)[1], 1:2, pop1530_size[,1:2])
    
    pfpr_3050_withIPTp = (1-frac_preg_3050) * pfpr_3050 + 
      frac_preg_3050 * pfpr_3050 * ((1 - iptp_coverage_table) + 
                                      iptp_coverage_table * frac_iptp_1_month * frac_time_unprotected_IPTp1 + 
                                      iptp_coverage_table * frac_iptp_2_month * frac_time_unprotected_IPTp2 + 
                                      iptp_coverage_table * frac_iptp_3_month * frac_time_unprotected_IPTp3)
    set(pfpr_3050_withIPTp, 1:dim(pop1530_size)[1], 1:2, pop1530_size[,1:2])
    
    
    # calculate overall PfPR for population with the adjusted PfPRs among these age groups
    pfpr_allAges_withIPTp = (pfpr_under15 * popUnder15_size + 
                               pfpr_1530_withIPTp * pop1530_size + 
                               pfpr_3050_withIPTp * pop3050_size + 
                               pfpr_over50 * popOver50_size) / 
      (popUnder15_size + pop1530_size + pop3050_size + popOver50_size)
    
    
    # calculate PfPR for population, unadjusted for IPTp
    pfpr_allAges_noIPTp = (pfpr_under15 * popUnder15_size + 
                             pfpr_1530 * pop1530_size + 
                             pfpr_3050 * pop3050_size + 
                             pfpr_over50 * popOver50_size) / 
      (popUnder15_size + pop1530_size + pop3050_size + popOver50_size)
    
    
    
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # adjustment of severe cases for malaria in pregnancy and IPTp
    # (occurs in 15to30 and 30to50 age groups)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    # women with malaria in pregnancy are more likely to have severe anemia (a type of severe malaria) than non-pregnant adults.
    # Schantz-Dunn and Nour (2009. 'Malaria and Pregnancy: A Global Health Perspective') suggest that the risk of severe disease is
    #   three times higher during pregnancy. They cite two books I can't access (Infectious Disease in Obstetrics and Gynecology and
    #   WHO, Guidelines for the Treatment of Malaria), but I'm not certain where the number is from or whether it is actually very robust.
    #   Actually, a number of groups cite the same number and sources.
    # Another approach is to look at the values reported in Desai et al. They report that 5-10% of pregnant women develop severe anaemia
    #   and around 26% of these instances are attributable to malaria (averaged across all gravidities, though paucigravidae).
    #   They also say that preventing malaria reduces risk of severe maternal anaemia by 38%
    # I think the best source is Guyatt and Snow 2001 (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.102.468&rep=rep1&type=pdf)
    #    Their lit review found that in Banfora, Burkina Faso, the  15.7% of first-time pregnant women with malaria had severe anemia compared to
    #    10.0% without malaria. As a guess from these numbers, we could say that 5.7% of women in their first pregnancy who are infected with
    #    malaria and who would not otherwise have had severe anemia get severe anemia as a result of malaria infection.
    # Because I don't have numbers on second pregnancies or later, I think I'll take a similar approach as for mLBW and assume that the first
    #    and second pregnancies have the same risk and that there is no risk for later pregnancies. This is probably too much risk for the second
    #    pregnancy and too little risk for the later ones, but lacking numbers to guide estimates, hopefully they will balance one another.
    # default value: prob_severe_from_MiP = 0.057 
    # 
    # however, we expect fewer individuals to get severe disease from MiP if they are using IPTps. Shulman et al. 1999 ('Intermittent
    #   sulphadoxine-pyrimethamine to prevent severe anaemia secondary to malaria in pregnancy: a randomised placebo-controlled trial')
    #   found a protective efficacy of 39% among women with SP (group included 1, 2, or 3 doses) compared to placebo for severe anaemia.
    #   This is among women with and without malaria, so actual protection would likely be higher among individuals with malaria.
    #    However, the study was done in an area of 74% prevalence and an area of 49% prevalence in 1-9 year old children.
    #   there didn't seem to be a difference depending on number of doses, so perhaps I'll use just two values: one for if IPTp taken
    #   and one for no IPTp.
    #   number from Desai et al. 2007: 'prevention of [malaria] infections [during pregnancy] reduces the risk of severe maternal anaemia by 38%'
    # default value: reduced_prob_severe_IPTp = 1-0.39 = 0.61
    prob_severe_from_MiP_withIPTp = prob_severe_from_MiP * reduced_prob_severe_IPTp
    
    # calculate the number of individuals who are giving birth in each month from 1st or 2nd pregnancy who were infected at some point during pregnancy
    num_1or2_preg_withAllStill_infected = num_1or2_preg_withAllStill * prob_infected_in_pregnancy_1530
    # replace year and month columns with original values (they were squared in previous step)
    set(num_1or2_preg_withAllStill_infected, 1:dim(num_1or2_preg_withAllStill)[1], 1:2, num_1or2_preg_withAllStill[,1:2])
    
    # calculate the number of these infectious-exposure-during-1st-or-2nd-pregnancies have severe disease using fraction of individuals with and without IPTp
    #   [# severe cases] = [# infected pregnancies without IPTp] * prob_severe_from_MiP + [# infected pregnancies with IPTp] * prob_severe_from_MiP_withIPTp
    #                    = [# infected pregnancies] * (1-[IPTp coverage]) * prob_severe_from_MiP + [# infected pregnancies] * [IPTp coverage] * prob_severe_from_MiP_withIPTp
    
    num_severe_from_MiP = num_1or2_preg_withAllStill_infected * (((1-iptp_coverage_table) * prob_severe_from_MiP) + (iptp_coverage_table * prob_severe_from_MiP_withIPTp))
    # replace year and month columns with original values (they were multiplied by coverage and probability of severe infection in previous step)
    set(num_severe_from_MiP, 1:dim(num_1or2_preg_withAllStill)[1], 1:2, num_1or2_preg_withAllStill[,1:2])
    
    
    
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # adjustment of clinical cases for malaria in pregnancy - NONE DONE CURRENTLY
    # (occurs in 15to30 and 30to50 age groups)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # several factors might cause MiP and IPTp to affect the number of clinical cases among pregnant women. First, taking IPTp may prevent 
    #    new malaria cases for the duration of protection (presumably contributing to a decrease in new clinical infections). However, in 
    #    high-transmission settings, most new infections in adults likely do not become clinical. Furthermore, by clearing existing 
    #    infections, there is potential for new infections to be clinical (at least according to how our current model would work if we 
    #    were to model it exlicitly instead of adjusting in post-processing). 
    # Second, given that women with malaria in pregnancy are more likely to have severe malaria, it seems reasonable to imagine that they  
    #    may be more likely to have a clinical case, but I haven't found any numbers on this. In high-transmission regions, many articles 
    #    I read indicated that even pregnant women weren't likely to have malaria symptoms, so for now I'll only count an increase in 
    #    severe disease due to anemia rather than also in clinical cases. In theory, this could lead to a larger number of severe cases
    #    in a month than clinical cases, but given the numbers I think it's unlikely.
    
    
    
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # adjustment of mortality from malaria in pregnancy - DONE FOR PLOTTING
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    # since I'm modifying severe malaria to account for anemia, I think I'll leave the death rate given severe disease alone. I'll specify the fraction
    #    of women with the treated death rate and with the untreated death rates. The untreated death rate for all cases would be more in line with 
    #    the estimate that almost 50% of women with severe anemia from malaria in pregnancy die from (Schantz-Dunn and Nour. 2009. Malaria and 
    #    Pregnancy: A Global Health Perspective). This might be an overestimate, since they cite the same sources I couldn't access before so I 
    #    don't know if they are considering treatment.
    # could probably use the CM values for adults with severe disease from corresponding simulations
    # default value: prob_severe_MiP_treated = 0.5

    
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #                     add columns to All_Age_Monthly_Cases.csv with updated MiP and IPTp outcomes
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # rows for the current simulation run
    rows_adjusted_allAge_rr = which(adjusted_allAgeMonthly$Run_Number == cur_run)
    # rows_adjusted_U5PfPR_rr = which(adjusted_U5PfPR$Run_Number == cur_run)

    print('Inserting new/adjusted values into existing simulation output All_Ages table...')
    
    for(i_MiP_col in 3:(dim(pop1530_size)[2])){
      if((i_MiP_col %% 40) == 0){
        print(paste('currently working on DS ', (i_MiP_col-2), 'out of ', (dim(pop1530_size)[2]-2)), sep='')
      }
      rows_adjusted_allAge_rr_ds = rows_adjusted_allAge_rr[which(toupper(adjusted_allAgeMonthly[[sim_ds_colname]][rows_adjusted_allAge_rr]) == toupper(colnames(pop1530_size)[i_MiP_col]))]

      # check whether the months and years are in the same order in both data tables (if so, can slot in new values directly without having to search for exact row)
      if(all(as.integer(adjusted_allAgeMonthly$month[rows_adjusted_allAge_rr_ds][1:length(pop1530_size$month)]) == pop1530_size$month) & all(as.integer(adjusted_allAgeMonthly$year[rows_adjusted_allAge_rr_ds][1:length(pop1530_size$year)]) == pop1530_size$year)){
        new_pfpr = as.vector(as.matrix(pfpr_allAges_withIPTp[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'PfPR_MiP_adjusted' := new_pfpr]
        
        unadjusted_pfpr = as.vector(as.matrix(pfpr_allAges_noIPTp[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'PfPR_unadjusted' := unadjusted_pfpr]
        
        new_severe_m = as.vector(as.matrix(num_severe_from_MiP[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'severe_maternal' := new_severe_m]

        new_severe_tot = as.vector(as.numeric(adjusted_allAgeMonthly$severe_maternal[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)]]) + as.numeric(adjusted_allAgeMonthly$New.Severe.Cases[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)]]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'severe_total' := new_severe_tot]

        new_severe_t = as.vector(as.numeric(adjusted_allAgeMonthly$severe_maternal[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)]]) * prob_severe_MiP_treated)
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'severe_maternal_treated' := new_severe_t]

        new_severe_t_tot = as.vector(as.numeric(adjusted_allAgeMonthly$severe_maternal_treated[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)]]) + as.numeric(adjusted_allAgeMonthly$Received_Severe_Treatment[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)]]))
        if(length(new_severe_t_tot) > 0) {
          adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'severe_total_treated' := new_severe_t_tot]
        } else{
          adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'severe_total_treated' := 0]
        }

        new_mLBW = as.vector(as.matrix(mLBW_each_DS_month[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'mLBW_births' := new_mLBW]

        new_mLBW_d = as.vector(as.matrix(mLBW_deaths_each_DS_month[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'mLBW_deaths' := new_mLBW_d]

        new_still = as.vector(as.matrix(mStill_each_DS_month[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'MiP_stillbirths' := new_still]

        new_still_noIPTp = as.vector(as.matrix(mStill_noIPTp_each_DS_month[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'MiP_stillbirths_noIPTp' := new_still_noIPTp]
        
        new_all_births = as.vector(as.matrix(num_all_births[, ..i_MiP_col]))
        adjusted_allAgeMonthly[rows_adjusted_allAge_rr_ds[1:length(pop1530_size$month)], 'num_all_births' := new_all_births]
        
      } else{
        print('order of months/years within simulation All_Ages output not the same as in existing data table so each row is being matched individually. This will take a while...')
        
        for(i_MiP_row in 13:dim(pop1530_size)[1]){
          # find the row corresponding to this Run_Number, DS, and date
          row_adjusted_allAge_rr_ds_date = rows_adjusted_allAge_rr_ds[intersect(which(as.integer(adjusted_allAgeMonthly$month[rows_adjusted_allAge_rr_ds])==pop1530_size$month[i_MiP_row]), which(adjusted_allAgeMonthly$year[rows_adjusted_allAge_rr_ds] == pop1530_size$year[i_MiP_row]))]

          # add new values in this row for the All_Age_Monthly csv file
          if(length(row_adjusted_allAge_rr_ds_date) == 1){
            adjusted_allAgeMonthly$PfPR_MiP_adjusted[row_adjusted_allAge_rr_ds_date] = pfpr_allAges_withIPTp[i_MiP_row, ..i_MiP_col]
            adjusted_allAgeMonthly$severe_maternal[row_adjusted_allAge_rr_ds_date] = num_severe_from_MiP[i_MiP_row, ..i_MiP_col]
            adjusted_allAgeMonthly$severe_total[row_adjusted_allAge_rr_ds_date] = as.numeric(adjusted_allAgeMonthly$severe_maternal[row_adjusted_allAge_rr_ds_date]) + as.numeric(adjusted_allAgeMonthly$New.Severe.Cases[row_adjusted_allAge_rr_ds_date])
            adjusted_allAgeMonthly$severe_maternal_treated[row_adjusted_allAge_rr_ds_date] = as.numeric(adjusted_allAgeMonthly$severe_maternal[row_adjusted_allAge_rr_ds_date]) * prob_severe_MiP_treated
            adjusted_allAgeMonthly$severe_total_treated[row_adjusted_allAge_rr_ds_date] = as.numeric(adjusted_allAgeMonthly$severe_maternal_treated[row_adjusted_allAge_rr_ds_date]) + as.numeric(adjusted_allAgeMonthly$Received_Severe_Treatment[row_adjusted_allAge_rr_ds_date])
            adjusted_allAgeMonthly$mLBW_births[row_adjusted_allAge_rr_ds_date] = mLBW_each_DS_month[i_MiP_row, ..i_MiP_col]
            adjusted_allAgeMonthly$mLBW_deaths[row_adjusted_allAge_rr_ds_date] = mLBW_deaths_each_DS_month[i_MiP_row, ..i_MiP_col]
            adjusted_allAgeMonthly$MiP_stillbirths[row_adjusted_allAge_rr_ds_date] = mStill_each_DS_month[i_MiP_row, ..i_MiP_col]
            adjusted_allAgeMonthly$MiP_stillbirths_noIPTp[row_adjusted_allAge_rr_ds_date] = mStill_noIPTp_each_DS_month[i_MiP_row, ..i_MiP_col]
            adjusted_allAgeMonthly$num_all_births[row_adjusted_allAge_rr_ds_date] = num_all_births[i_MiP_row, ..i_MiP_col]
          } else{
            warning(paste('match not found in All_Age for ds ', colnames(pop1530_size)[i_MiP_col], '; run num ', cur_run, '; month ', pop1530_size$month[i_MiP_row], '; year ', pop1530_size$year[i_MiP_row], sep=''))
          }
        }
      }
    }

  } # finish iterating through runs in simulation output
  fwrite(adjusted_allAgeMonthly, sim_output_allAgeMonthly_filename_MiP)
  # fwrite(adjusted_U5PfPR, sim_output_U5PfPR_filename_MiP)
}





