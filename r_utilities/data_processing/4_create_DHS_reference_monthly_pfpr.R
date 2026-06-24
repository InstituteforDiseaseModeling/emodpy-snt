#  4_create_DHS_reference_monthly_pfpr.R


# create data frame with a row entry for each admin-month-year present in the DHS dataset, giving the total number tested and the number microscopy positive
create_DHS_reference_monthly_pfpr = function(hbhi_dir, dta_dir, admin_shape, ds_pop_df_filename, pfpr_dhs_ref_years=c(2012, 2016), min_num_total=30, pfpr_measure='mic'){
  if(pfpr_measure=='mic'){
    pfpr_measure_name = 'microscopy'
  }else if(pfpr_measure=='rdt'){
    pfpr_measure_name = 'RDT'
  }else{
    warning('Name of PfPR metric not recognized. Assuming microscopy.')
    pfpr_measure='mic'
    pfpr_measure_name = 'microscopy'
  }
  
  for (yy in 1:length(pfpr_dhs_ref_years)){
    year = pfpr_dhs_ref_years[yy]
    DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',year,'_files_recodes_for_sims.csv'))
    location_index = which(DHS_file_recode_df$variable == 'locations')
    locations_shp = shapefile(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[location_index], '/', DHS_file_recode_df$filename[location_index]))
    locations = data.frame(clusterid = locations_shp$DHSCLUST, latitude=locations_shp$LATNUM, longitude=locations_shp$LONGNUM)
    MIS_outputs = locations
    
    var_index = which(DHS_file_recode_df$variable == pfpr_measure)
    if(!is.na(DHS_file_recode_df$filename[var_index])){
      cur_dta = read.dta(paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index]))
      
      cur_dta$pos = NA
      cur_dta$pos[which(cur_dta[,grep(colnames(cur_dta), pattern = DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$pos_pattern[var_index])] = 1
      cur_dta$pos[which(cur_dta[,grep(colnames(cur_dta), pattern = DHS_file_recode_df$code[var_index])] == DHS_file_recode_df$neg_pattern[var_index])] = 0      
      cur_dta$month = cur_dta[,grep(colnames(cur_dta), pattern = DHS_file_recode_df$month_code[var_index])]
      cur_dta$year = cur_dta[,grep(colnames(cur_dta), pattern = DHS_file_recode_df$year_code[var_index])]
      
      dta_cluster_0 = cur_dta  %>%
        filter(!is.na(cur_dta$pos)) %>%
        group_by_at(c(DHS_file_recode_df$cluster_id_code[var_index], 'month', 'year')) %>%
        summarize(rate = mean(pos, na.rm = TRUE),
                  num_pos = sum(pos),
                  num_tested = n()) 
      MIS_outputs = merge(MIS_outputs, dta_cluster_0, by.y=DHS_file_recode_df$cluster_id_code[var_index], by.x='clusterid', all=TRUE)
      
      
      ####=========================================================================================================####
      # determine which clusters are in which admins
      ####=========================================================================================================####
      # turn MIS output data frame into spatial points data frame
      points_crs = crs(admin_shape)
      MIS_shape = SpatialPointsDataFrame(MIS_outputs[,c('longitude', 'latitude')],
                                         MIS_outputs,
                                         proj4string = points_crs)
      # find which admins each cluster belongs to
      MIS_locations = over(MIS_shape, admin_shape)
      if(nrow(MIS_locations) == nrow(MIS_shape)){
        MIS_shape$admin_name = MIS_locations$NOMDEP
        MIS_shape$NOMREGION = MIS_locations$NOMREGION
        MIS_shape$NAME_1 = MIS_locations$NAME_1
      }
      
      # combine this year's dataframe with previous years
      if(yy ==1){
        ds_pfpr_all_years = as.data.frame(MIS_shape)
      } else{
        ds_pfpr_all_years = rbind(ds_pfpr_all_years, as.data.frame(MIS_shape))
      }
    }else warning(paste0('filename not specified for ', pfpr_measure_name, ' year ', pfpr_dhs_ref_years[yy]))
  }
  ds_pfpr_all_years = ds_pfpr_all_years[!is.na(ds_pfpr_all_years$admin_name),]
  
  # combine all entries from clusters in the same month-year-admin
  ds_pfpr_all_years_agg = ds_pfpr_all_years %>% 
    group_by_at(c('admin_name', 'NOMREGION', 'month', 'year')) %>%
    summarize(num_pos = sum(num_pos),
              num_tested = sum(num_tested))
  
    
  
  
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # add any admins that did not have any DHS clusters
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  archetype_info = read.csv(ds_pop_df_filename)
  colnames(archetype_info)[colnames(archetype_info)=='State'] = 'NOMREGION'

  # standardize admin names
  ds_pfpr_all_years_agg = standardize_admin_names_in_df(target_names_df=archetype_info, origin_names_df=ds_pfpr_all_years_agg, target_names_col='admin_name', origin_names_col='admin_name')
  
  # merge to add in missing admins
  archetype_info = archetype_info[,c('admin_name', 'NOMREGION', 'Archetype')]
  ds_pfpr_all_years_agg_expanded = merge(ds_pfpr_all_years_agg, archetype_info, all=TRUE)
  # check that all names have been matched successfully
  if(length(ds_pfpr_all_years_agg_expanded$admin_name[which(!(ds_pfpr_all_years_agg_expanded$admin_name %in% archetype_info$admin_name))])>0) warning('Not all LGA names from the shapefile were matched with archetype file')
  if(length(ds_pfpr_all_years_agg_expanded$NOMREGION[which(!(ds_pfpr_all_years_agg_expanded$NOMREGION %in% archetype_info$NOMREGION))])>0) warning('Not all state names from the shapefile were matched with archetype file')

    # change to zero sample size for admins that were not included in DHS
  ds_pfpr_all_years_agg_expanded$num_pos[is.na(ds_pfpr_all_years_agg_expanded$num_pos)] = 0
  ds_pfpr_all_years_agg_expanded$num_tested[is.na(ds_pfpr_all_years_agg_expanded$num_tested)] = 0
  
  admin_pfpr_sums = ds_pfpr_all_years_agg_expanded
  admin_pfpr_sums$data_spatial_level = 'admin'
  write.csv(as.data.frame(admin_pfpr_sums), paste0(hbhi_dir, '/estimates_from_DHS/DHS_monthly_', pfpr_measure_name, '_adminLevelDataOnly.csv'), row.names=FALSE)
  
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # when the number surveyed in a admin is lower than the threshold, use the region value instead
  #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # region level values: combine all entries from clusters in the same month-year-region
  region_pfpr_sums_agg = admin_pfpr_sums %>% 
    group_by_at(c('NOMREGION', 'month', 'year', 'Archetype')) %>%
    summarize(num_pos = sum(num_pos),
              num_tested = sum(num_tested))
  
  # determine the total number in each admin across all years/months of sampling
  agg_to_admin = admin_pfpr_sums %>% group_by(admin_name, NOMREGION) %>%
    summarise(total_positive = sum(num_pos),
              total_tested = sum(num_tested))
  # hist(agg_to_admin$total_tested, breaks=60)
  
  # what are the total number sampled in each state/region when aggregated across years/months?
  agg_to_region = region_pfpr_sums_agg %>% group_by(NOMREGION) %>%
    summarise(total_positive = sum(num_pos),
              total_tested = sum(num_tested))
  # hist(agg_to_region$total_tested, breaks=60)
  if(any(agg_to_region$total_tested < min_num_total)) warning('Some regions do not have the minimum number of required U5 tested with ', pfpr_measure_name,', need to modify code to aggregate further.')
  
  admin_pfpr_sums_adjMin = admin_pfpr_sums
  # if the total number of U5s surveyed across all months/years/clusters is less than the threshold min_num_total, replace the admin-level data with the aggregated state-level data instead
  for(i_admin in 1:nrow(agg_to_admin)){
    cur_admin = agg_to_admin$admin_name[i_admin]
    cur_region = agg_to_admin$NOMREGION[i_admin]
    if(agg_to_admin$total_tested[i_admin] < min_num_total){
      # remove all of the entries for this admin and replace them with the corresponding entries from the region/state
      admin_rows_to_remove = which(admin_pfpr_sums_adjMin$admin_name == cur_admin)
      region_rows_to_add = which(region_pfpr_sums_agg$NOMREGION == cur_region)
      added_df = region_pfpr_sums_agg[region_rows_to_add,]
      added_df$data_spatial_level = 'region'
      added_df$admin_name = cur_admin
      added_df = added_df[added_df$num_tested>0,]

      # remove the old rows from admin_pfpr_sums_adjMin and add in the new rows
      admin_pfpr_sums_adjMin = admin_pfpr_sums_adjMin[-admin_rows_to_remove,]
      admin_pfpr_sums_adjMin = merge(admin_pfpr_sums_adjMin, added_df, all=TRUE)
    }
  }   
  admin_pfpr_sums_adjMin = admin_pfpr_sums_adjMin[(admin_pfpr_sums_adjMin$num_tested>0) & !is.na(admin_pfpr_sums_adjMin$year),]
  write.csv(as.data.frame(admin_pfpr_sums_adjMin), paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_monthly_',pfpr_measure_name,'.csv'), row.names=FALSE)
}

# # plot distribution of aggregated prevalences across LGAs
# hist(agg_to_region$total_positive/agg_to_region$total_tested, breaks=40)
# hist(admin_pfpr_sums$num_pos[admin_pfpr_sums$num_tested>20]/admin_pfpr_sums$num_tested[admin_pfpr_sums$num_tested>20])

