# add_permethrin_mortality_to_itn_coverage_files.R
# Contact: mambrose
# Created: 2021
# Contents: functions that add permethrin bioassay mortality values to all ITN coverage in a project

library(reshape2)


add_mortality_to_file = function(itn_coverage_filename, mortality_df_name, itn_coverage_dir, itn_coverage_mort_dir){
  #' Given an ITN coverage file with 'admin_name' and 'year' columns, find the matching permethrin resistance value for that location and time and add it as a column to the ITN file
  #' @param itn_coverage_filename Full filepath of ITN coverage file. 
  #'                              Assumes the year of the distribution is in a column called 'dist_year' or 'year' and the location (either the name of the admin or the archetype) 
  #'                              is in a column named 'admin_name' or 'seasonality_archetype', respectively
  #' @param mortality_df_name Full filepath to file giving bioassay mortality values. 
  #'                          Assumes admin name is in a column called 'DS' and the year is in a column called 'year'
  #' @param itn_coverage_dir Filepath to intermediate directory for temporary files that contain ITN coverage information (but yet not mortality/killing/blocking)
  #' @param itn_coverage_mort_dir Filepath to intermediate directory for temporary files that contain ITN coverage and mortality information (but yet not killing/blocking)
  
  # read in ITN coverage csv
  itn_cov_df = read.csv(paste0(itn_coverage_dir, '/', itn_coverage_filename))
  colnames(itn_cov_df)[colnames(itn_cov_df) == 'dist_year'] = 'year'
  
  # read-in and reformat mortality dataframe
  mortality_df_0 = read.csv(mortality_df_name)
  colnames(mortality_df_0)[which(colnames(mortality_df_0) %in% c('DS', 'LGA'))] = 'admin_name'
  colnames(mortality_df_0) = gsub('X','',colnames(mortality_df_0))  # get rid of the X in front of years in column names
  # get rid of extra indexing column if there is one
  if(all(mortality_df_0[,1] == 1:nrow(mortality_df_0))){
    mortality_df_0 = mortality_df_0[,-1]
  }
  mortality_df = reshape2::melt(mortality_df_0, id.vars = 'admin_name')
  colnames(mortality_df)[which(colnames(mortality_df)=='variable')] = 'year'
  colnames(mortality_df)[which(colnames(mortality_df)=='value')] = 'bio_mortality'
  mortality_df$year = as.integer(as.character(mortality_df$year))

  # if there are years in the itn coverage df that are not in the mortality df, find the closest year from the mortality df and use those values
  mort_years = unique(mortality_df$year)
  for(yy in unique(itn_cov_df$year)){
    if(!(yy %in% mort_years)){
      closest_year = mort_years[which.min(abs(mort_years - yy))]
      # add new rows to mortality df for the best match to this year
      mort_yy = mortality_df[mortality_df$year == closest_year,]
      mort_yy$year = yy
      mortality_df = rbind(mortality_df, mort_yy)
    }
  }

  # merge the two dataframes
  if('admin_name' %in% colnames(itn_cov_df)){
    itn_cov_mort = merge(itn_cov_df, mortality_df, by=c('year', 'admin_name'))
  } else{
    itn_cov_mort = merge(itn_cov_df, mortality_df, by.x=c('year', 'seasonality_archetype'), by.y=c('year', 'admin_name'))
  }
  
  # save csv
  itn_coverage_mort_filename = gsub('_coverage_dates_LGA','', itn_coverage_filename)
  itn_coverage_mort_filename = gsub('_coverages','', itn_coverage_mort_filename)
  itn_coverage_mort_filename = gsub('_coverage','', itn_coverage_mort_filename)
  write.csv(itn_cov_mort, paste0(itn_coverage_mort_dir, '/', itn_coverage_mort_filename), row.names=FALSE)
}


add_bioassay_mortality_to_itn_file = function(hbhi_dir, mortality_df_name){
  #' Given a directory where there all ITN-coverage files are stored, iterate through all rows of all files, adding the appropriate bioassay mortality for the location and year.
  #' @param hbhi_dir Filepath to the base hbhi project directory
  #' @param mortality_df_name Full filepath to file giving bioassay mortality values. 
  
  itn_coverage_dir = paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')  # intermediate files that only have coverage information
  itn_coverage_mort_dir = paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage_mort')  # intermediate files that have coverage and mortality (but not killing and blocking)
  ifelse(!dir.exists(itn_coverage_mort_dir), dir.create(itn_coverage_mort_dir), FALSE)
  itn_files = list.files(itn_coverage_dir, pattern = '.csv')
  # remove state-level file created for records
  itn_files = itn_files[!grepl('state',itn_files)]
  
  for(ii in 1:length(itn_files)){
    itn_coverage_filename = itn_files[ii]
    add_mortality_to_file(itn_coverage_filename, mortality_df_name, itn_coverage_dir, itn_coverage_mort_dir)
  }
}
