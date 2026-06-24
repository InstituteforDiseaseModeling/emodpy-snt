# plot_sim_input_functions.R

# library(rgdal)
library(raster)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(sf)
library(reshape2)
library(data.table)
library(dplyr)
library(geofacet)
library(ggpattern)


separate_plot_text_size=12
text_size = 15
save_plots = TRUE


# Returns TRUE if the cached state-level timeseries CSV at `cache_path` exists,
# overwrite is not forced, AND its set of scenarios exactly matches the
# requested `expected_scenarios`. The plot_state_grid_* functions use this
# instead of a bare file.exists() check so that adding a new scenario to the
# coordinator (or renaming one) automatically invalidates the cached
# per-state aggregate -- otherwise the new scenario silently fails to appear
# in the geofacet plot.
cache_matches_scenarios = function(cache_path, expected_scenarios, overwrite_files){
  if(overwrite_files) return(FALSE)
  if(!file.exists(cache_path)) return(FALSE)
  cached = tryCatch(read.csv(cache_path, stringsAsFactors = FALSE),
                    error = function(e) NULL)
  if(is.null(cached) || !('scenario' %in% colnames(cached))) return(FALSE)
  cached_scenarios = unique(cached$scenario)
  # 'to-present' is appended at write time from per-scenario rows; allow it
  # to be present in the cache even when not in the caller's scenario list.
  cached_extras = setdiff(cached_scenarios, c(expected_scenarios, 'to-present'))
  missing_from_cache = setdiff(expected_scenarios, cached_scenarios)
  if(length(missing_from_cache) > 0 || length(cached_extras) > 0){
    message(sprintf(
      "  Cache %s out of date (missing: %s; extra: %s) -- regenerating.",
      basename(cache_path),
      if(length(missing_from_cache) > 0) paste(missing_from_cache, collapse = ', ') else '-',
      if(length(cached_extras) > 0) paste(cached_extras, collapse = ', ') else '-'))
    return(FALSE)
  }
  TRUE
}




#####################################################################################################
# ================================================================================================= #
# create state-grid-arranged plot panel with intervention timeseries for included scenarios
# ================================================================================================= #
#####################################################################################################

##############
# CM
##############
plot_state_grid_cm = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                              min_year, max_year, sim_end_years, 
                              scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                              separate_admin_lines_flag = FALSE,  act_adherence_effective_multiplier=1, overwrite_files=FALSE){

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]

  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether CM output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_cm_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_cm_timeseries_by_state(cm_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_cm_timeseries_by_state(cm_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                        min_year=min_year, get_state_level=TRUE)
      }

      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  # adjust for ACT adherence to get coverage of individuals taking ACTs (not reduced effective coverage due to not completing the regimen)
  if(act_adherence_effective_multiplier<1){
    timeseries_df$mean_coverage = timeseries_df$mean_coverage / act_adherence_effective_multiplier
    timeseries_df$min_coverage = timeseries_df$min_coverage / act_adherence_effective_multiplier
    timeseries_df$max_coverage = timeseries_df$max_coverage / act_adherence_effective_multiplier
  }
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Effective treatment rate (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 

  } else{
      gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
        scale_fill_manual(values = scenario_palette) + 
        geom_line(linewidth=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Effective treatment rate (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_CM_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}



##############
# ITN ANC
##############
plot_state_grid_itn_anc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                    min_year, max_year, sim_end_years, 
                                    scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                    separate_admin_lines_flag = FALSE, overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether CM output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_itn_anc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input ITN ANC coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$ANC_ITN_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_itn_anc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                             min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_itn_anc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                             min_year=min_year, get_state_level=TRUE)
      }

      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }

  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.5) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ANC ITN coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ANC ITN coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_ITN_ANC_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}




##############
# ITN EPI
##############
plot_state_grid_itn_epi = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                   min_year, max_year, sim_end_years, 
                                   scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                   separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether CM output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_itn_epi_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input ITN epi coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$EPI_ITN_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_itn_epi_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                             min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_itn_epi_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                             min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('EPI ITN coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('EPI ITN coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_ITN_EPI_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}






##############
# SMC
##############

plot_state_grid_smc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                               min_year, max_year, sim_end_years, 
                               scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                               separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_smc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$SMC_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_smc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_smc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC coverage (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
        # scale_fill_manual(values = scenario_palette) + 
        geom_line(linewidth=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC coverage (U5)')) + 
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+ 
        # theme(legend.position = "none", text = element_text(size = text_size))+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_SMC_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}




##############
# SPAQ (combined SMC + IPTsc) -- one state-grid plot per age band
##############
# Mirrors plot_state_grid_smc but reads SPAQ_filename and plots a single age band
# (band_min/band_max), writing band-specific output files. Effective coverage per band
# is computed in get_spaq_timeseries_by_state.
plot_state_grid_spaq = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations,
                                min_year, max_year, sim_end_years,
                                scenario_names, scenario_input_references, experiment_names, scenario_palette,
                                band_min, band_max, band_label = NULL,
                                separate_admin_lines_flag = FALSE, overwrite_files=FALSE){

  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  if(is.null(band_label)) band_label = paste0(band_min, '-', band_max)
  band_file_tag = gsub('[^A-Za-z0-9]+', '', gsub('\\.', 'p', band_label))  # filename-safe, e.g. 0p25to5

  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }

  # check whether SPAQ output already exists for this comparison + band
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_spaq_age', band_file_tag,
                               '_state_', time_string, 'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input SPAQ coverages for this band
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$SPAQ_filename[cur_int_row], '.csv')

      cur_timeseries_agg = get_spaq_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info,
                                                        end_year=end_year, exp_name=scenario_names[ee],
                                                        min_year=min_year, band_min=band_min, band_max=band_max,
                                                        get_state_level=!separate_admin_lines_flag)

      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }

    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }


  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  ylab_string = paste0('SPAQ coverage (ages ', band_label, ')')

  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) +
        scale_color_manual(values = scenario_palette) +
        xlab('year') +
        ylab(ylab_string) +
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free')

  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(linewidth=1) +
        scale_color_manual(values = scenario_palette) +
        xlab('year') +
        ylab(ylab_string) +
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free')
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/', time_string, 'Timeseries_SPAQ_age', band_file_tag,
                '_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}





##############
# ITN mass campaign
##############

plot_state_grid_itn = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                               min_year, max_year, sim_end_years, 
                               scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                               separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_itn_mass_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$ITN_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_itn_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_itn_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ITN initial use rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      # scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('ITN initial use rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_ITN_mass_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}








##############
# Vaccines
##############

plot_state_grid_vacc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                min_year, max_year, sim_end_years, 
                                scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_vacc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$vacc_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_vacc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_vacc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                         min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('EPI malaria vaccine coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      # scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('EPI malaria vaccine coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_vacc_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}







##############
# PMC
##############

plot_state_grid_pmc = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                                min_year, max_year, sim_end_years, 
                                scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }
  
  # check whether SMC output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_pmc_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$PMC_filename[cur_int_row], '.csv')
      
      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_pmc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                          min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_pmc_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                                          min_year=min_year, get_state_level=TRUE)
      }
      
      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }
    
    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      # join past and future simulation trajectories
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  timeseries_df = timeseries_df[timeseries_df$scenario %in% scenario_names,]
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State
  
  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('PMC coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
    
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      # scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('PMC coverage')) + 
      coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_PMC_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}









#####################################################################
# plot map of admin subsets
#####################################################################
# NOTE: kept in sync with the copy in plot_sim_output_functions.R (which is the one the GC8 results
# script sources). Thin LGA borders + bold state borders (state borders dissolved from the LGA polygons
# by the 'State' column of pop_filepath).
plot_included_admin_map = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, admin_shapefile_filepath, shapefile_admin_colname,
                                   save_width=4.8, save_height=4.8, lga_linewidth=0.12, state_linewidth=0.5,
                                   lga_border_color='grey50', state_border_color='black'){
  admin_pop = read.csv(pop_filepath)
  admin_shapefile = st_read(admin_shapefile_filepath)
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=admin_pop$admin_name, origin_names=admin_shapefile[[shapefile_admin_colname]])

  admin_in_map = data.frame(admin_name = admin_pop$admin_name, admin_included='no')
  admin_in_map$admin_included[admin_in_map$admin_name %in% cur_admins] = 'yes'
  included_colors = c('yes'='#006692', 'no'='grey96')

  # shapefile already carries its own 'State' column, so join only inclusion status (merging State
  # from pop would collide -> State.x/State.y and break the dissolve below).
  admin_cur = admin_shapefile %>%
    dplyr::left_join(admin_in_map, by=c('NOMDEP' = 'admin_name'))

  gg_map = ggplot(admin_cur) +
    geom_sf(aes(fill=admin_included), linewidth=lga_linewidth, color=lga_border_color) +   # thin LGA borders
    scale_fill_manual(values=included_colors, drop=FALSE, na.value='grey96') +
    theme_void() +
    theme(legend.position = 'none')

  # st_make_valid() first: the shapefile has minor invalid geometries that otherwise break st_union.
  if('State' %in% names(admin_cur)){
    state_sf = sf::st_make_valid(admin_cur[!is.na(admin_cur$State), ]) %>%
      dplyr::group_by(State) %>%
      dplyr::summarise(.groups = 'drop')
    gg_map = gg_map + geom_sf(data = state_sf, fill = NA, color = state_border_color, linewidth = state_linewidth)
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/map_admins_included_', district_subset, '.png'), gg_map, dpi=600, width=save_width, height=save_height, units='in')
}




##############
# IRS
##############

plot_state_grid_irs = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations,
                               min_year, max_year, sim_end_years,
                               scenario_names, scenario_input_references, experiment_names, scenario_palette,
                               separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]

  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }

  # check whether IRS output already exists for this comparison
  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_irs_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$IRS_filename[cur_int_row], '.csv')

      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_irs_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee],
                                                         min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_irs_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee],
                                                         min_year=min_year, get_state_level=TRUE)
      }

      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }

    if(any(grepl('to-present', timeseries_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }


  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State

  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) +
        scale_color_manual(values = scenario_palette) +
        xlab('year') +
        ylab(paste0('IRS effective coverage')) +
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free')

  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(linewidth=1) +
        scale_color_manual(values = scenario_palette) +
        xlab('year') +
        ylab(paste0('IRS effective coverage')) +
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free')
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_IRS_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}


##############
# IPTp
##############

plot_state_grid_iptp = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations,
                                 min_year, max_year, sim_end_years,
                                 scenario_names, scenario_input_references, experiment_names, scenario_palette,
                                 separate_admin_lines_flag = FALSE,  overwrite_files=FALSE){

  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]

  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  time_string = 'annual'
  if(separate_admin_lines_flag){
    separate_admin_string = '_separated_admins'
  } else{
    separate_admin_string = ''
  }

  timeseries_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_iptp_state_',time_string,'Timeseries', separate_admin_string, '.csv')
  if(cache_matches_scenarios(timeseries_filepath, scenario_names, overwrite_files)){
    timeseries_df = read.csv(timeseries_filepath)
  } else{
    timeseries_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      input_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$IPTp_filename[cur_int_row], '.csv')

      if(separate_admin_lines_flag){
        cur_timeseries_agg = get_iptp_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee],
                                                          min_year=min_year, get_state_level=FALSE)
      } else{
        cur_timeseries_agg = get_iptp_timeseries_by_state(input_filepath=input_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee],
                                                          min_year=min_year, get_state_level=TRUE)
      }

      if(nrow(timeseries_df)==0){
        timeseries_df = cur_timeseries_agg
      } else{
        timeseries_df = rbind(timeseries_df, cur_timeseries_agg)
      }
    }

    if(any(grepl('to-present', timeseries_df$scenario))){
      to_present_df = timeseries_df[timeseries_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        timeseries_df = rbind(timeseries_df, final_to_present_row)
      }
    }
    write.csv(timeseries_df, timeseries_filepath, row.names=FALSE)
  }

  timeseries_df$scenario = factor(timeseries_df$scenario, levels=rev(scenario_names))
  timeseries_df$code = timeseries_df$State

  if(separate_admin_lines_flag){
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(aes(group=interaction(admin_name, scenario), color=scenario), linewidth=0.8) +
        scale_color_manual(values = scenario_palette) +
        xlab('year') +
        ylab(paste0('IPTp coverage (>=1 dose)')) +
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free')
  } else{
    gg = ggplot(timeseries_df, aes(x=year, y=mean_coverage, color=scenario)) +
        geom_line(linewidth=1) +
        scale_color_manual(values = scenario_palette) +
        xlab('year') +
        ylab(paste0('IPTp coverage (>=1 dose)')) +
        coord_cartesian(xlim=c(min_year, max_year), ylim=c(0,1))+
        scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
        theme_bw()+
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free')
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_IPTp_by_state', separate_admin_string, '.png'), gg, dpi=600, width=12, height=10, units='in')
}
