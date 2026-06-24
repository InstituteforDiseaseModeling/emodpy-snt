# malariaMortality_functions.R
library(data.table)


# read in simulation output and calculate updated numbers from multiple sources of mortality

adjust_sim_output_mortality = function(cfr_severe_treated = 0.097,  # case fatality rate for severe, treated cases
                                       cfr_severe_untreated_1 = 0.539,  # case fatality rate for severe, untreated cases (upper estimate)
                                       cfr_severe_untreated_2 = 0.177,  # case fatality rate for severe, untreated cases (lower estimate)
                                       starting_prob_1 = 0.037,  # (QD from Ross et al, upper value)
                                       starting_prob_2 = 0.01,  # (QD from Ross et al, lower value)
                                       age_shape_param = 0.117,  # (aF* from Ross et al)
                                       sim_output_U1PfPR_filename_old,
                                       sim_output_U5PfPR_filename_old,
                                       sim_output_allAgeMonthly_filename_old,
                                       sim_output_MiP_mortality_filename, 
                                       pop_size_U1_colname = 'Pop.U1',
                                       pop_size_U5_colname = 'Pop.U5',
                                       cases_U1_colname = 'Cases.U1',
                                       cases_U5_colname = 'Cases.U5',
                                       cases_all_colname = 'New.Clinical.Cases',
                                       severe_U1_colname = 'Severe.cases.U1',
                                       severe_U5_colname = 'Severe.cases.U5',
                                       severe_all_colname = 'New.Severe.Cases',  # does not include maternal severe (due to order of rescaling between IPTp and IPTi)
                                       severe_treated_U1_colname = 'Num_U1_Received_Severe_Treatment',
                                       severe_treated_U5_colname = 'Num_U5_Received_Severe_Treatment',
                                       severe_treated_all_colname = 'severe_total_treated',  # includes treatment from non-maternal and maternal severe malaria
                                       severe_maternal_colname = 'severe_maternal'  # cases of severe disease due to maternal anemia (needs to be added to New Severe Cases)
                                       ){
  #' read in simulation output and calculate updated numbers from multiple sources of mortality
  #' 
  #'  @param cfr_severe_treated  case fatality rate for severe, treated cases
  #'  @param cfr_severe_untreated_1 case fatality rate for severe, treated cases (upper estimate)
  #'  @param cfr_severe_untreated_2 case fatality rate for severe, treated cases (lower estimate)
  #'  @param starting_prob_1 QD from Ross et al, upper value
  #'  @param starting_prob_2 QD from Ross et al, lower value
  #'  @param age_shape_param aF* from Ross et al
  #'  @param sim_output_U1PfPR_filename_old filepath for the U1 simulation output file that will be modified
  #'  @param sim_output_U5PfPR_filename_old filepath for the U5 simulation output file that will be modified
  #'  @param sim_output_allAgeMonthly_filename_old filepath for the all-age simulation output file that will be modified
  #'  @param sim_output_MiP_mortality_filename filepath for the adjusted simulation output file that is created
  #'  @param pop_size_U1_colname column name corresponding to 'Pop U1'
  #'  @param pop_size_U5_colname column name corresponding to 'Pop U5'
  #'  @param cases_U1_colname column name corresponding to 'Cases U1'
  #'  @param cases_U5_colname column name corresponding to 'Cases U5'
  #'  @param cases_all_colname column name corresponding to 'New Clinical Cases'
  #'  @param severe_U1_colname column name corresponding to 'Severe cases U1'
  #'  @param severe_U5_colname column name corresponding to 'Severe cases U5'
  #'  @param severe_all_colname column name corresponding to 'New Severe Cases' - does not include maternal severe (due to order of rescaling between IPTp and IPTi)
  #'  @param severe_treated_U1_colname  column name corresponding to 'Num_U1_Received_Severe_Treatment'
  #'  @param severe_treated_U5_colname  column name corresponding to 'Num_U5_Received_Severe_Treatment'
  #'  @param severe_treated_all_colname column name corresponding to 'severe_total_treated' - includes treatment from non-maternal and maternal severe malaria
  #'  @param severe_maternal_colname column name corresponding to 'severe_maternal' - cases of severe disease due to maternal anemia (needs to be added to New Severe Cases)
  
  if(file.exists(sim_output_U1PfPR_filename_old)){
    include_U1=TRUE
  } else include_U1 = FALSE
  
  # read in simulation files that will be adjusted
  sim_output_allAgeMonthly = fread(sim_output_allAgeMonthly_filename_old, check.names=TRUE)
  sim_output_U5PfPR  = fread(sim_output_U5PfPR_filename_old, check.names=TRUE)
  if(include_U1) sim_output_U1PfPR  = fread(sim_output_U1PfPR_filename_old, check.names=TRUE)
  
  if('admin_name' %in% colnames(sim_output_allAgeMonthly)){
    sim_ds_colname = 'admin_name'
  } else if('LGA' %in% colnames(sim_output_allAgeMonthly)){
    sim_ds_colname = 'LGA'
  } else{
    sim_ds_colname = 'DS_Name'
  }
  # replace any spaces in column names with periods
  colnames(sim_output_allAgeMonthly) = gsub(' ', '.', colnames(sim_output_allAgeMonthly))
  colnames(sim_output_U5PfPR) = gsub(' ', '.', colnames(sim_output_U5PfPR))
  if(include_U1) colnames(sim_output_U1PfPR) = gsub(' ', '.', colnames(sim_output_U1PfPR))
  

  # create new data tables with additional columns
  if(include_U1){
    empty_column_values = as.numeric(rep(NA, dim(sim_output_U1PfPR)[1]))
    adjusted_U1PfPR = cbind(sim_output_U1PfPR, 
                            data.table(New_clinical_cases_U1=empty_column_values,
                                       Severe_cases_U1=empty_column_values,
                                       direct_mortality_nonMiP_U1_1=empty_column_values,
                                       direct_mortality_nonMiP_U1_2=empty_column_values,
                                       indirect_mortality_nonMiP_U1_1=empty_column_values,
                                       indirect_mortality_nonMiP_U1_2=empty_column_values,
                                       total_mortality_U1_1=empty_column_values,
                                       total_mortality_U1_2=empty_column_values
                            ))
  }
  
  empty_column_values = as.numeric(rep(NA, dim(sim_output_U5PfPR)[1]))
  adjusted_U5PfPR = cbind(sim_output_U5PfPR, 
                          data.table(New_clinical_cases_U5=empty_column_values,
                                     Severe_cases_U5=empty_column_values,
                                     direct_mortality_nonMiP_U5_1=empty_column_values,
                                     direct_mortality_nonMiP_U5_2=empty_column_values,
                                     indirect_mortality_nonMiP_U5_1=empty_column_values,
                                     indirect_mortality_nonMiP_U5_2=empty_column_values,
                                     total_mortality_U5_1=empty_column_values,
                                     total_mortality_U5_2=empty_column_values
                          ))
  
  empty_column_values = as.numeric(rep(NA, dim(sim_output_allAgeMonthly)[1]))
  adjusted_allAgeMonthly = cbind(sim_output_allAgeMonthly, 
                                 data.table(direct_mortality_nonMiP_1=empty_column_values,
                                            direct_mortality_nonMiP_2=empty_column_values,
                                            indirect_mortality_nonMiP_1=empty_column_values,
                                            indirect_mortality_nonMiP_2=empty_column_values,
                                            total_mortality_1=empty_column_values,
                                            total_mortality_2=empty_column_values
                                            ))

  
  # add columns for the year and month
  setDT(adjusted_allAgeMonthly)[, year := format(as.Date(date), "%Y") ]
  adjusted_allAgeMonthly$year = as.numeric(adjusted_allAgeMonthly$year)
  setDT(adjusted_allAgeMonthly)[, month := format(as.Date(date), "%m") ]
  adjusted_allAgeMonthly$month = as.numeric(adjusted_allAgeMonthly$month)
  

  # merge data tables to get a single data frame that has outcome information for all age groups
  adjusted_output = merge(adjusted_allAgeMonthly, adjusted_U5PfPR, by=c(sim_ds_colname, 'Run_Number', 'month', 'year'), all=TRUE)
  if(include_U1) adjusted_output = merge(adjusted_output, adjusted_U1PfPR, by=c(sim_ds_colname, 'Run_Number', 'month', 'year'), all=TRUE)
  
  
  # NA values in the column for number of U1 or U5 receiving severe treatment should be zero
  if(include_U1) adjusted_output[[severe_treated_U1_colname]][is.na(adjusted_output[[severe_treated_U1_colname]])] = 0
  adjusted_output[[severe_treated_U5_colname]][is.na(adjusted_output[[severe_treated_U5_colname]])] = 0
  # rename the severe cases column name for clarity
  if(include_U1) colnames(adjusted_output)[which(colnames(adjusted_output)==severe_U1_colname)] = 'severe_case_U1_rate'
  colnames(adjusted_output)[which(colnames(adjusted_output)==severe_U5_colname)] = 'severe_case_U5_rate'
  # calculate the number of clinical cases and severe clinical cases among children U5 (from the per-person, per-month daily rate)
  if(include_U1) adjusted_output$New_clinical_cases_U1 = adjusted_output[[cases_U1_colname]] * adjusted_output[[pop_size_U1_colname]] * 30 / 365
  adjusted_output$New_clinical_cases_U5 = adjusted_output[[cases_U5_colname]] * adjusted_output[[pop_size_U5_colname]] * 30 / 365
  if(include_U1) adjusted_output$Severe_cases_U1 = adjusted_output$severe_case_U1_rate * adjusted_output[[pop_size_U1_colname]] * 30 / 365
  adjusted_output$Severe_cases_U5 = adjusted_output$severe_case_U5_rate * adjusted_output[[pop_size_U5_colname]] * 30 / 365
  
  
  # add direct mortality
  if(include_U1) adjusted_output$direct_mortality_nonMiP_U1_1 = (adjusted_output[[severe_treated_U1_colname]] * cfr_severe_treated) + sapply((adjusted_output$Severe_cases_U1 - adjusted_output[[severe_treated_U1_colname]]), max, 0) * cfr_severe_untreated_1
  if(include_U1) adjusted_output$direct_mortality_nonMiP_U1_2 = (adjusted_output[[severe_treated_U1_colname]] * cfr_severe_treated) + sapply((adjusted_output$Severe_cases_U1 - adjusted_output[[severe_treated_U1_colname]]), max, 0) * cfr_severe_untreated_2
  adjusted_output$direct_mortality_nonMiP_U5_1 = (adjusted_output[[severe_treated_U5_colname]] * cfr_severe_treated) + sapply((adjusted_output$Severe_cases_U5 - adjusted_output[[severe_treated_U5_colname]]), max, 0) * cfr_severe_untreated_1
  adjusted_output$direct_mortality_nonMiP_U5_2 = (adjusted_output[[severe_treated_U5_colname]] * cfr_severe_treated) + sapply((adjusted_output$Severe_cases_U5 - adjusted_output[[severe_treated_U5_colname]]), max, 0) * cfr_severe_untreated_2
  adjusted_output$direct_mortality_nonMiP_1 = (adjusted_output[[severe_treated_all_colname]] * cfr_severe_treated) + sapply((adjusted_output[[severe_all_colname]] + adjusted_output[[severe_maternal_colname]] - adjusted_output[[severe_treated_all_colname]]), max, 0) * cfr_severe_untreated_1
  adjusted_output$direct_mortality_nonMiP_2 = (adjusted_output[[severe_treated_all_colname]] * cfr_severe_treated) + sapply((adjusted_output[[severe_all_colname]] + adjusted_output[[severe_maternal_colname]] - adjusted_output[[severe_treated_all_colname]]), max, 0) * cfr_severe_untreated_2

  
  # add indirect mortality - since ages aren't broken out more finely in the All_Age_Monthly and the U5_PfPR files, and the probability of death asymptotes fairly quickly after age 5, use two parameter values: one for all under 5 and one for all over 5
  # note that if Clinical cases among U5 are adjusted for IPTi, those adjusted numbers should be used here.
  prob_indirect_death_U1_1 = mean(starting_prob_1/(1+seq(0,1,0.1)/age_shape_param)) # ages 0-1
  prob_indirect_death_U1_2 = mean(starting_prob_2/(1+seq(0,1,0.1)/age_shape_param))
  prob_indirect_death_U5_1 = mean(starting_prob_1/(1+seq(0,5,0.1)/age_shape_param)) # ages 0-5
  prob_indirect_death_U5_2 = mean(starting_prob_2/(1+seq(0,5,0.1)/age_shape_param))
  prob_indirect_death_O5_1 = mean(starting_prob_1/(1+seq(5,80,0.1)/age_shape_param)) # ages >=5
  prob_indirect_death_O5_2 = mean(starting_prob_2/(1+seq(5,80,0.1)/age_shape_param))
  
  # will use any clinical case, including severe - so a very small fraction of deaths will be double counted (when an individual is estimated to die both from severe disease and 
  #    from indirect mortality, but given the numbers we're using, this value should be so small as to be negligable). the alternative would be to exclude severe disease (assume no one
  #    with severe disease dies from indirect mortality), which would error in the other direction
  if(include_U1) adjusted_output$indirect_mortality_nonMiP_U1_1 = adjusted_output$New_clinical_cases_U1 * prob_indirect_death_U1_1
  if(include_U1) adjusted_output$indirect_mortality_nonMiP_U1_2 = adjusted_output$New_clinical_cases_U1 * prob_indirect_death_U1_2
  adjusted_output$indirect_mortality_nonMiP_U5_1 = adjusted_output$New_clinical_cases_U5 * prob_indirect_death_U5_1
  adjusted_output$indirect_mortality_nonMiP_U5_2 = adjusted_output$New_clinical_cases_U5 * prob_indirect_death_U5_2
  # adults age >=5 only (exclude U5)
  O5_mortality_nonMiP_1 = (adjusted_output[[cases_all_colname]] - adjusted_output$New_clinical_cases_U5) * prob_indirect_death_O5_1
  O5_mortality_nonMiP_2 = (adjusted_output[[cases_all_colname]] - adjusted_output$New_clinical_cases_U5) * prob_indirect_death_O5_2
  # all ages: add U5 and O5 deaths. Note that we only consider the indirect deaths for the entire age group 0-5 here, rather than calculating 0-1 and 1-5 separately. This may mean there could be inconsistencies with the indirect U1 death counts.
  adjusted_output$indirect_mortality_nonMiP_1 = adjusted_output$indirect_mortality_nonMiP_U5_1 + O5_mortality_nonMiP_1
  adjusted_output$indirect_mortality_nonMiP_2 = adjusted_output$indirect_mortality_nonMiP_U5_2 + O5_mortality_nonMiP_2
  
  
  # total mortality - does not include stillbirths, does include mortality from severe disease, indirect mortality, maternal mortality from severe anemia, and mLBW deaths
  if(include_U1) adjusted_output$total_mortality_U1_1 = adjusted_output$indirect_mortality_nonMiP_U1_1 + adjusted_output$direct_mortality_nonMiP_U1_1 + adjusted_output$mLBW_deaths
  if(include_U1) adjusted_output$total_mortality_U1_2 = adjusted_output$indirect_mortality_nonMiP_U1_2 + adjusted_output$direct_mortality_nonMiP_U1_2 + adjusted_output$mLBW_deaths
  adjusted_output$total_mortality_U5_1 = adjusted_output$indirect_mortality_nonMiP_U5_1 + adjusted_output$direct_mortality_nonMiP_U5_1 + adjusted_output$mLBW_deaths
  adjusted_output$total_mortality_U5_2 = adjusted_output$indirect_mortality_nonMiP_U5_2 + adjusted_output$direct_mortality_nonMiP_U5_2 + adjusted_output$mLBW_deaths
  adjusted_output$total_mortality_1 = adjusted_output$direct_mortality_nonMiP_1 + adjusted_output$indirect_mortality_nonMiP_1 + adjusted_output$mLBW_deaths  # includes severe maternal cases from MiP
  adjusted_output$total_mortality_2 = adjusted_output$direct_mortality_nonMiP_2 + adjusted_output$indirect_mortality_nonMiP_2 + adjusted_output$mLBW_deaths  # includes severe maternal cases from MiP
  
  # replace periods or spaces in column names with underscores
  colnames(adjusted_output) = str_replace_all(colnames(adjusted_output),'\\.','_')
  colnames(adjusted_output) = str_replace_all(colnames(adjusted_output),' ','_')
  
  # final check for multiple rows (either identical duplicates or repeats of what should be a single row with different values)
  check_column_names = c(sim_ds_colname,'year','Run_Number','month')
  if (any(duplicated(adjusted_output[,..check_column_names]))){
    warning('YIKES! Multiple rows detected during final stage of processing for at least one admin-year-month-run combination. This may result in the final burden csv having more than one row per output.')
  } 
  
  # calculate mean of two mortality estimates as well as returning each estimate alone
  adjusted_output$direct_mortality_nonMiP_mean = (adjusted_output$direct_mortality_nonMiP_1 +adjusted_output$direct_mortality_nonMiP_2)/2
  adjusted_output$indirect_mortality_nonMiP_mean = (adjusted_output$indirect_mortality_nonMiP_1 +adjusted_output$indirect_mortality_nonMiP_2)/2
  adjusted_output$total_mortality_mean = (adjusted_output$total_mortality_1 +adjusted_output$total_mortality_2)/2
  adjusted_output$direct_mortality_nonMiP_U5_mean = (adjusted_output$direct_mortality_nonMiP_U5_1 +adjusted_output$direct_mortality_nonMiP_U5_2)/2
  adjusted_output$indirect_mortality_nonMiP_U5_mean = (adjusted_output$indirect_mortality_nonMiP_U5_1 +adjusted_output$indirect_mortality_nonMiP_U5_2)/2
  adjusted_output$total_mortality_U5_mean = (adjusted_output$total_mortality_U5_1 +adjusted_output$total_mortality_U5_2)/2
  adjusted_output$direct_mortality_nonMiP_U1_mean = (adjusted_output$direct_mortality_nonMiP_U1_1 +adjusted_output$direct_mortality_nonMiP_U1_2)/2
  adjusted_output$indirect_mortality_nonMiP_U1_mean = (adjusted_output$indirect_mortality_nonMiP_U1_1 +adjusted_output$indirect_mortality_nonMiP_U1_2)/2
  adjusted_output$total_mortality_U1_mean = (adjusted_output$total_mortality_U1_1 +adjusted_output$total_mortality_U1_2)/2

  # write output to csv file
  fwrite(adjusted_output, sim_output_MiP_mortality_filename)
}



# Description of select columns in output file (stars next to those relevant for plotting output):
#   *Statistical_Population - among all ages (original simulation output)
#   PfHRP2_Prevalence - original simulation output, among all ages, does not include adjustments for IPTp
#   Received_Severe_treatment - original simulation output - number of individuals with severe disease who were treated, does not include treatment for severe disease from MiP
#   *New_Clinical_Cases - adjusted for IPTi (no change from MiP or IPTp)
#   New_Severe_Cases - original simulation output - does not include severe disease from MiP
#   *PfPR_MiP_adjusted - PfPR for all ages, adjusted to account for IPTp (and IPTi if relevant) - PfPR comes from MalariaSummaryReport, so I believe is measured by microscopy
#   severe_maternal - expected number of individuals with a new occurrence of severe anemia due to MiP
#   severe_total - sum of 'New Severe Cases' from original simulation output and severe_maternal. represents number of individuals with severe disease due to malaria for the entire population.
#   severe_maternal_treated - expected number of individuals treated among those with a new occurrence of severe anemia due to MiP
#   mLBW_births - number of infants expected to be born with low birth rate due to MiP
#   mLBW_deaths - number of mLBW infants expected to die within their first year who would have survived in the absence of MiP
#   *MiP_stillbirths - number of stillbirths attributable to MiP with assumed IPTp coverage
#   MiP_stillbirths_noIPTp - number of stillbirths attributable to MiP if no IPTp were used in the population
#   direct_mortality_nonMiP_1 - number of deaths due to severe disease using larger estimate for CFR for untreated severe cases (includes mortality from severe maternal anemia, does not include mLBW or stillbirths)
#   direct_mortality_nonMiP_2 - number of deaths due to severe disease using smaller estimate for CFR for untreated severe cases (includes mortality from severe maternal anemia, does not include mLBW or stillbirths)
#   indirect_mortality_nonMiP_1 - number of deaths that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses larger estimate for CFR
#   indirect_mortality_nonMiP_2 - number of deaths that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses smaller estimate for CFR
#   *total_mortality_1 - number of deaths in population directly or indirectly attributable to malaria (includes severe materal anemia, deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses larger estimates for CFRs.
#   *total_mortality_2 - number of deaths in population directly or indirectly attributable to malaria (includes severe materal anemia, deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses smaller estimates for CFRs.
#   *PfPR_U5 - PfPR among U5, adjusted for IPTi - PfPR comes from MalariaSummaryReport, so I believe is measured by microscopy
#   *Pop_U5 - population size among individuals U5 (original simulation output)
#   Num_U5_Received_Severe_Treatment - number of U5 individuals with severe disease who were treated from simulation (does not include mLBW or stillbirths)
#   *New_clinical_cases_U5 - number of U5 individuals with clinical cases (does not include mLBW or stillbirths), adjusted for IPTi
#   Severe_cases_U5 - number of U5 individuals with severe disease (does not include mLBW or stillbirths)
#   direct_mortality_nonMiP_U5_1 - number of deaths in U5 due to severe disease using larger estimate for CFR for untreated severe cases (does not include mLBW or stillbirths)
#   direct_mortality_nonMiP_U5_2 - number of deaths in U5 due to severe disease using smaller estimate for CFR for untreated severe cases (does not include mLBW or stillbirths)
#   indirect_mortality_nonMiP_U5_1 - number of deaths in U5 that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses larger estimate for CFR
#   indirect_mortality_nonMiP_U5_2 - number of deaths in U5 that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses smaller estimate for CFR
#   *total_mortality_U5_1 - number of deaths in U5 directly or indirectly attributable to malaria (includes deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses larger estimates for CFRs.
#   *total_mortality_U5_2 - number of deaths in U5 directly or indirectly attributable to malaria (includes deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses smaller estimates for CFRs.
#    direct_mortality_nonMiP_mean - number of deaths due to severe disease using mean of CFR estimates for untreated severe cases (includes mortality from severe maternal anemia, does not include mLBW or stillbirths)
#    indirect_mortality_nonMiP_mean - number of deaths that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses mean of CFR estimates
#    total_mortality_mean - number of deaths in population directly or indirectly attributable to malaria (includes severe materal anemia, deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses mean of CFR estimates
#    direct_mortality_nonMiP_U5_mean - number of deaths in U5 due to severe disease using mean of CFR estimates for untreated severe cases (does not include mLBW or stillbirths)
#    indirect_mortality_nonMiP_U5_mean - number of deaths in U5 that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses mean of CFR estimates
#    total_mortality_U5_mean - number of deaths in U5 directly or indirectly attributable to malaria (includes deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses mean of CFR estimates
#    direct_mortality_nonMiP_U1_mean - number of deaths in U1 due to severe disease using mean of CFR estimates for untreated severe cases (does not include mLBW or stillbirths)
#    indirect_mortality_nonMiP_U1_mean - number of deaths in U1 that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses mean of CFR estimates
#    total_mortality_U1_mean - number of deaths in U1 directly or indirectly attributable to malaria (includes deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses mean of CFR estimates
#    num_all_births - number of births (both still and live births) occurring in the month









