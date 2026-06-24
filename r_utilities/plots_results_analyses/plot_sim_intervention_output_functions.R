# plot_sim_output_functions.R

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




######################################################################
# create plot panel with intervention distribution info
######################################################################

plot_simulation_intervention_output_timeseries_by_state = function(sim_output_dir, pop_filepath, grid_layout_state_locations,
                                               plot_by_month, min_year, max_year, sim_end_years, 
                                               scenario_filepaths, scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                               indoor_protection_fraction=0.75, remove_exp_name_substring='', overwrite_files=FALSE){
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  pop_sizes = read.csv(pop_filepath)
  pop_sizes = pop_sizes[,c('admin_name', 'State', 'pop_size')]
  
  # create output directories
  if(!dir.exists(paste0(sim_output_dir, '/_plots'))) dir.create(paste0(sim_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  

  # ----- LLIN, vaccine, and IRS intervention coverage ----- #
  
  # check whether LLIN/IRS output already exists for this comparison
  inter_df_filepath = paste0(sim_output_dir, '/_plots/timeseries_dfs/df_intervention_use_',time_string,'Timeseries.csv')
  if(file.exists(inter_df_filepath) & !overwrite_files){
    inter_use_df = read.csv(inter_df_filepath)
  } else{
    # iterate through scenarios, storing relevant output
    inter_use_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_inter_use = get_intervention_use_timeseries_state(exp_filepath = scenario_filepaths[ee],
                                                        exp_name = scenario_names[ee], 
                                                        pop_sizes=pop_sizes, min_year=min_year, max_year=max_year, indoor_protection_fraction=indoor_protection_fraction, plot_by_month=plot_by_month)
      if(nrow(inter_use_df)==0){
        inter_use_df = cur_inter_use
      } else{
        inter_use_df = merge(inter_use_df, cur_inter_use, all=TRUE)
      }
    }
    write.csv(inter_use_df, inter_df_filepath, row.names=FALSE)
  }
  
  
  # ----- Case management ----- #
  
  # check whether CM output already exists for this comparison
  cm_df_filepath = paste0(sim_output_dir, '/_plots/timeseries_dfs/df_cm_',time_string,'Timeseries.csv')
  if(file.exists(cm_df_filepath)){
    cm_df = read.csv(cm_df_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    cm_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath) %>%
        mutate(ScenarioName = gsub(remove_exp_name_substring, '', ScenarioName))
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      cm_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      cur_cm_agg = get_cm_timeseries_by_state(cm_filepath=cm_filepath, admin_info=pop_sizes, end_year=end_year, exp_name= scenario_names[ee],
                                 min_year=min_year, plot_by_month=plot_by_month)
        
      if(nrow(cm_df)==0){
        cm_df = cur_cm_agg
      } else{
        cm_df = rbind(cm_df, cur_cm_agg)
      }
    }
    write.csv(cm_df, cm_df_filepath, row.names=FALSE)
  }
  
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots for each intervention, faceted by state
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  inter_use_df$code = inter_use_df$State
  # ----- LLIN use and distribution ----- #
  # plot net use through time
  if(plot_by_month){
    g_net_use = ggplot(inter_use_df, aes(x=as.Date(date), y=itn_used_per_cap, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('LLIN use (all ages)')) + 
      coord_cartesian(ylim=c(0,NA))+
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') +
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  } else{
    g_net_use = ggplot(inter_use_df, aes(x=year, y=itn_used_per_cap, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      # geom_hline(yintercept=0.22, alpha=0.1)+
      # geom_hline(yintercept=0.39, alpha=0.1)+
      xlab('year') + 
      ylab(paste0('LLIN use (all ages)')) + 
      coord_cartesian(ylim=c(0,NA))+
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  }
  ggsave(paste0(sim_output_dir, '/_plots/',time_string,'Timeseries_ITN_by_state.png'), g_net_use, dpi=600, width=12, height=10, units='in')
  
  
  
  # ----- vaccine ----- #
  if('vacc_per_cap' %in% colnames(inter_use_df)){
    if(plot_by_month){
      g_vacc = ggplot(inter_use_df, aes(x=as.Date(date), y=vacc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('Vaccines (primary series + booster) per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_vacc = ggplot(inter_use_df, aes(x=year, y=vacc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Vaccines (primary series + booster) per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    }
    ggsave(paste0(sim_output_dir, '/_plots/',time_string,'Timeseries_Vacc_by_state.png'), g_vacc, dpi=600, width=12, height=10, units='in')
  }
  
  
  # ----- PMC ----- #
  if('pmc_per_cap' %in% colnames(inter_use_df)){
    if(plot_by_month){
      g_pmc = ggplot(inter_use_df, aes(x=as.Date(date), y=pmc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('PMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_pmc = ggplot(inter_use_df, aes(x=year, y=pmc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('PMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    }
    ggsave(paste0(sim_output_dir, '/_plots/',time_string,'Timeseries_PMC_by_state.png'), g_pmc, dpi=600, width=12, height=10, units='in')
}
  
  
  # ----- SMC ----- #
  if('smc_per_cap' %in% colnames(inter_use_df)){
    if(plot_by_month){
      g_smc = ggplot(inter_use_df, aes(x=as.Date(date), y=smc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('SMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_smc = ggplot(inter_use_df, aes(x=year, y=smc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    }
    ggsave(paste0(sim_output_dir, '/_plots/',time_string,'Timeseries_SMC_by_state.png'), g_smc, dpi=600, width=12, height=10, units='in')
  }
  
  
  # ----- IRS ----- #
  if('irs_per_cap' %in% colnames(inter_use_df)){
    if(plot_by_month){
      g_irs = ggplot(inter_use_df, aes(x=as.Date(date), y=irs_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('IRS rounds per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_irs = ggplot(inter_use_df, aes(x=year, y=irs_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('IRS per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    }
    ggsave(paste0(sim_output_dir, '/_plots/',time_string,'Timeseries_IRS_by_state.png'), g_irs, dpi=600, width=12, height=10, units='in')
  }
  
  
  # ----- Case management ----- #
  cm_df$code = cm_df$State
  if(plot_by_month){
    g_cm = ggplot(cm_df, aes(x=as.Date(date), y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(ylim=c(0,NA))+
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  } else{
    g_cm = ggplot(cm_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(ylim=c(0,NA))+
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='fixed') +
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  }
  ggsave(paste0(sim_output_dir, '/_plots/',time_string,'Timeseries_CM_by_state.png'), g_cm, dpi=600, width=12, height=10, units='in')
  
}







