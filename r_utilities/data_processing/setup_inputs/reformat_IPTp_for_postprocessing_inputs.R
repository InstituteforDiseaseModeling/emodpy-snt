# reformat_IPTp_for_postprocessing_inputs.R
library(reshape2)

reformat_IPTp_input = function(hbhi_dir, ds_pop_df_filename){
  ds_pop_df = read.csv(ds_pop_df_filename)
  
  # read in seed 1 (best-guess estiamte) from the sampled coverage from DHS observations for number of individuals with >=1 IPTp doses
  iptp_sampled_coverage = read.csv(paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/iptp_2010_toPresent.csv'))
  if('seed' %in% colnames(iptp_sampled_coverage)){
    iptp_coverage_est = iptp_sampled_coverage[iptp_sampled_coverage$seed == 1,]
  } else iptp_coverage_est = iptp_sampled_coverage
  iptp_coverage_est = iptp_coverage_est[,which(colnames(iptp_coverage_est) %in% c('admin_name', 'year','coverage'))]

   # reformat so that admins are rows and years are columns
  iptp_coverage_est_wide = reshape2::dcast(iptp_coverage_est, admin_name ~ year, value.var="coverage")
  
  # if some admin don't have IPTp, add rows with zero coverage
  if(!all(ds_pop_df$admin_name %in% iptp_coverage_est_wide$admin_name)){
    missing_admins = ds_pop_df$admin_name[which(!(ds_pop_df$admin_name %in% iptp_coverage_est_wide$admin_name))]
    for(ii in 1:length(missing_admins)){
      new_row = iptp_coverage_est_wide[1,]
      new_row[,2:ncol(new_row)] = 0
      new_row$admin_name = missing_admins[ii]
      iptp_coverage_est_wide = rbind(iptp_coverage_est_wide, new_row)
    }
  }
  if(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/IPTp'))) dir.create(paste0(hbhi_dir, '/simulation_inputs/IPTp'))
  write.csv(iptp_coverage_est_wide, paste0(hbhi_dir,'/simulation_inputs/IPTp/estimated_past_IPTp_each_DS.csv'))
}






