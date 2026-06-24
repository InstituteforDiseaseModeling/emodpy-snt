# rescale_IPTp_based_on_dose_trajectory.R

# We have IPTp coverage from DHS, but only up until 2016. We also have data from the PNILP on number of people receiving IPTp in each month and year from 2017 onwards.
#   To project coverage increases after 2016, we look at how number of IPTp treatments increases and rescale accordingly (assuming the 2017 number data represents the 2016 coverage)

library(readxl)
library(dplyr)
library(ggplot2)


# helper function to extract the total number of IPTp3 doses given in each DS for a given year (summed across health facilities)
get_admin_sums = function(cur_sheet, colname_admin, colname_total){
  # remove rows without DS and unneeded columns
  cur_sheet = cur_sheet[,(grepl(colname_admin, colnames(cur_sheet))) | grepl(colname_total, colnames(cur_sheet))]
  if(any(is.na(cur_sheet[[colname_admin]]))) cur_sheet = cur_sheet[-(which(is.na(cur_sheet[[colname_admin]]))),]
  
  # get sum of number of doses (across all months and all health facilities in a DS)
  df_sums0 = cur_sheet %>% 
    group_by_at(all_of(colname_admin)) %>%
    summarise_all(sum, na.rm=TRUE)
  df_sums = data.frame(admin_name=df_sums0[[colname_admin]], count=rowSums(df_sums0[,-which(grepl(colname_admin, colnames(df_sums0)))]))
  
  # remove the 'DS ' from the name and standardize
  df_sums$admin_name = gsub('DS ', '', df_sums$admin_name)
  # standardize shapefile names
  df_sums$admin_name = standardize_admin_names_in_vector(target_names=archetype_info$admin_name, origin_names=df_sums$admin_name)
  
  return(df_sums[,c('admin_name', 'count')])
}



rescale_IPTp_later_years = function(hbhi_dir, script_dir, ds_pop_df_filename, rescale_iptp_dose_filename, maximum_coverage_iptp){
  source(paste0(script_dir,'/standardize_admin_names.R'))
  
  # filepath for original (and updated) IPTp csv file
  iptp_input_filename = paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/iptp_2010_toPresent.csv')
  
  # get the names of each sheet (the year).
  sheet_names = excel_sheets(path=rescale_iptp_dose_filename)
  years_included = sort(as.integer(sheet_names))
  # specify which columns are used
  colname_admin = 'District'
  # value_column_substring = 'CPN_Femmes'  # will sum across all columns containing this substring
  colname_total = 'Total'
  
  # get dataframe with standardized admin names, states, archetype, and population information
  archetype_info = read.csv(ds_pop_df_filename)
  
  
  # get results for each year - for IPTp3
  iptp_counts = data.frame()
  for (yy in years_included){
    cur_sheet = read_excel(path = rescale_iptp_dose_filename, sheet = sheet_names[grepl(yy, sheet_names)])
    cur_counts = get_admin_sums(cur_sheet=cur_sheet, colname_admin=colname_admin, colname_total=colname_total)
    cur_counts$year = yy
    if(nrow(iptp_counts)>0){
      iptp_counts = merge(iptp_counts, cur_counts, all=TRUE)
    } else{
      iptp_counts = cur_counts
    }
  }
  
  
  # get the trajectory from >=3 doses in each DS
  iptp_counts$dose_num='3 or more'
  # divide by 2017 value to get comparable magnitudes (and since we'll be rescaling based on 2017 value later)
  iptp_counts_2017 = iptp_counts[iptp_counts$year == 2017,]
  colnames(iptp_counts_2017)[colnames(iptp_counts_2017)=='count'] = 'count_2017'
  iptp_counts = merge(iptp_counts, iptp_counts_2017[,c('admin_name', 'count_2017')], all=TRUE)
  iptp_counts$count_rel = iptp_counts$count / iptp_counts$count_2017
  iptp_counts$coverage_multiplier = iptp_counts$count_rel
  
  
  # read in old IPTp coverage file and rescale appropriately based on changes in dose numbers
  # If X1 is the number of doses needed for P1 coverage, then what coverage is achieved with X2 doses?
  #     P1=X1/N  (N=target population size, assume this is relatively constant between years)
  #     P2=X2/N 
  #     --> N=X1/P1 --> P2=P1*X2/X1
  iptp_sampled_coverage = read.csv(iptp_input_filename)
  # check that the sampled coverage is the same after min(iptp_counts$year)
  if(length(unique(iptp_sampled_coverage$coverage[iptp_sampled_coverage$admin_name==iptp_sampled_coverage$admin_name[1] & iptp_sampled_coverage$year>=min(iptp_counts$year) & iptp_sampled_coverage$seed==1]))==1){
    write.csv(iptp_sampled_coverage, gsub('.csv','_original.csv',iptp_input_filename), row.names=FALSE)
    
    # rescale future values
    iptp_sampled_coverage = merge(iptp_sampled_coverage, iptp_counts[,c('admin_name', 'year', 'coverage_multiplier')], all=TRUE)
    iptp_sampled_coverage$coverage_multiplier[iptp_sampled_coverage$year<min(iptp_counts$year)] = 1
    iptp_sampled_coverage_updated = iptp_sampled_coverage
    iptp_sampled_coverage_updated$coverage = iptp_sampled_coverage$coverage * iptp_sampled_coverage$coverage_multiplier
    
    # set maximum realistic coverage
    iptp_sampled_coverage_updated$coverage = sapply(iptp_sampled_coverage_updated$coverage, min, maximum_coverage_iptp)
    
    # create plot of results
    if(!(dir.exists(file.path(hbhi_dir, 'simulation_inputs','IPTp')))) dir.creage(file.path(hbhi_dir, 'simulation_inputs','IPTp'))
    gg=ggplot(iptp_sampled_coverage_updated[iptp_sampled_coverage_updated$seed==1,], aes(x=year, y=coverage, col=admin_name))+
      geom_line()+
      theme(legend.position = 'none')
    ggsave(file.path(hbhi_dir, 'simulation_inputs','IPTp','rescaled_coverage_timeseries.png'),gg,width=4, height=3)
    
    # add any extra years as repeat of the final year - NOT CURRENTLY INCLUDED - add error check instead
    if(max(iptp_counts$year)<max(iptp_sampled_coverage_updated$year)){
      warning('PROBLEM: the maximum year in the IPTp coverage dataframe is greater than that in the dose number data. Need to add support for this in the code (e.g., continuing the coverage from the latest year)')
    }
    write.csv(iptp_sampled_coverage_updated, iptp_input_filename, row.names=FALSE)
  }
}



# ##############################
# ##### Previous work comparing the difference of using >=3 doses, 1-2 doses, or the mean of the two for hte trajectories
# ##############################
# # more recent data on number of nets distributed through each channel in each year (all year-channel combinations saved as different sheets in excel file)
# new_iptp_1or2_filename = paste0(dta_dir, '/burundi/WHO/snt_2022/Donn�es sur TPIg 1et2.xlsx')
# # get the names of each sheet (the year).
# sheet_names_a = excel_sheets(path=new_iptp_1or2_filename)
# 
# # get results for each year - for IPTp1or2
# iptp_counts_a = data.frame()
# for (yy in years_included){
#   cur_sheet = read_excel(path = new_iptp_1or2_filename, sheet = sheet_names_a[grepl(yy, sheet_names_a)])
#   cur_counts = get_admin_sums(cur_sheet=cur_sheet, colname_admin=colname_admin, colname_total=colname_total)
#   cur_counts$year = yy
#   if(nrow(iptp_counts_a)>0){
#     iptp_counts_a = merge(iptp_counts_a, cur_counts, all=TRUE)
#   } else{
#     iptp_counts_a = cur_counts
#   }
# }
# 
# # combine the trajectories from 1-2 doses and >=3 doses
# iptp_counts_a$dose_num='1 or 2'
# # divide by 2017 value to get comparable magnitudes (and since we'll be rescaling based on 2017 value later)
# iptp_counts_a_2017 = iptp_counts_a[iptp_counts_a$year == 2017,]
# colnames(iptp_counts_a_2017)[colnames(iptp_counts_a_2017)=='count'] = 'count_2017'
# iptp_counts_a = merge(iptp_counts_a, iptp_counts_a_2017[,c('admin_name', 'count_2017')], all=TRUE)
# iptp_counts_a$count_rel = iptp_counts_a$count / iptp_counts_a$count_2017
# # merge the 1-2 and >3 results
# iptp_counts_all = merge(iptp_counts_a, iptp_counts, all=TRUE)
# iptp_counts_comb = iptp_counts_all %>% group_by(admin_name, year) %>%
#   summarise(coverage_multiplier=mean(count_rel))
# 
# # plot separately
# ggplot(iptp_counts_all, aes(x=year, y=count_rel, col=admin_name, linetype=dose_num))+
#   geom_line()+
#   theme(legend.position = 'none')
# # plot together
# ggplot(iptp_counts_comb, aes(x=year, y=coverage_multiplier, col=admin_name))+
#   geom_line()+
#   theme(legend.position = 'none')
# 
# # plan to use just the IPTp3 dose rates for rescaling, which will likely slightly overestimate the rate of 1 and 2 doses (since we're assuming the 
# #     same dose distribution), but since we are skipping the growth between 2016 and 2017, that's probably a minimal impact
# iptp_counts %>% group_by(year) %>% summarise(median_multiplier = median(count_rel), mean_multiplier = mean(count_rel))
# iptp_counts_a %>% group_by(year) %>% summarise(median_multiplier = median(count_rel), mean_multiplier = mean(count_rel))
# iptp_counts_comb %>% group_by(year) %>% summarise(median_multiplier = median(coverage_multiplier), mean_multiplier = mean(coverage_multiplier))
