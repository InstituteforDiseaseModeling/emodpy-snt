# functions_counterfactual_timeseries.R



######################################################################
# create plot panel with all burden metrics, no intervention info
######################################################################

plot_burden_reduction_attribution = function(sim_output_dir, pop_filepath, district_subset, cur_admins, 
                                             plot_by_month, min_year, max_year, 
                                             scenario_filepaths, scenario_names, experiment_names, inter_ordered, scenario_palette, scenario_palette_stacked, 
                                             burden_metrics='PfPR', burden_metric_names='PfPR (U5)', burden_colnames='PfPR_U5', 
                                             LLIN2y_flag=FALSE, overwrite_files=FALSE, noInterName='noInter', allInterName='allInter',
                                             firstInPrefix='only',lastInPrefix='allInterExcept', create_plots_flag=TRUE){
  
  
  
  pop_sizes = read.csv(pop_filepath)
  pop_sizes = pop_sizes[,c('admin_name','pop_size')]
  # if we include all admins, get list of names from population size dataframe
  if(cur_admins[1] == 'all'){
    cur_admins = unique(pop_sizes$admin_name)
  }
  
  # create output directories
  if(!dir.exists(paste0(sim_output_dir, '/_plots'))) dir.create(paste0(sim_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  
  # ----- malaria burden ----- #
  
  for(bb in 1:length(burden_colnames)){
    burden_metric_name = burden_metric_names[bb]
    burden_colname = burden_colnames[bb]
    burden_metric = burden_metrics[bb]
    
    
    if(grepl('U5', burden_metric_name)){
      age_plotted = 'U5'
    } else age_plotted = 'all'
    
    # check whether burden output already exists for this comparison
    if(LLIN2y_flag){
      llin2y_string = '_2yLLIN'
    } else{
      llin2y_string = ''
    }
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # create data frame with annual burden estimate for all experiments
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    
    burden_df_filepath = paste0(sim_output_dir, '/_plots/timeseries_dfs/df_burden_',time_string,'Timeseries_', burden_metric, '_', age_plotted, '_',district_subset, llin2y_string,'.csv')
    if(file.exists(burden_df_filepath) & !overwrite_files){
      burden_df = read.csv(burden_df_filepath)
    } else{
      # iterate through scenarios, storing relevant output
      burden_df = data.frame()
      for(ee in 1:length(scenario_filepaths)){
        cur_sim_output_agg = get_burden_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                       exp_name = scenario_names[ee], district_subset=district_subset,
                                                       cur_admins=cur_admins, pop_sizes=pop_sizes, min_year=min_year, max_year=max_year, burden_colname=burden_colname, age_plotted=age_plotted, plot_by_month=plot_by_month)
        if(nrow(burden_df)==0){
          burden_df = cur_sim_output_agg
        } else{
          burden_df = rbind(burden_df, cur_sim_output_agg)
        }
      }
      write.csv(burden_df, burden_df_filepath, row.names=FALSE)
    }
    
    
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # calculate burden averted by intervention when first-in and last-out
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    burden_averted_filepath = paste0(sim_output_dir, '/_plots/timeseries_dfs/df_burden_averted_',time_string,'Timeseries_', burden_metric, '_', age_plotted, '_',district_subset, llin2y_string,'.csv')
    if(file.exists(burden_averted_filepath) & !overwrite_files){
      averted_df = read.csv(burden_averted_filepath)
    } else{
      ### look at impact attributed to each intervention when it is the first (only) intervention added on top of a simulation without interventions
      first_in_df = burden_df[grepl(firstInPrefix, burden_df$scenario),]
      first_in_df$intervention = sapply(first_in_df$scenario, gsub, pattern=firstInPrefix, replacement='')
      # format the no-intervention reference data frame
      noInter_df = burden_df[burden_df$scenario == noInterName,]
      noInter_df = noInter_df %>% rename(ref_mean_burden=mean_burden) %>%  # , ref_max_burden=max_burden, ref_min_burden=min_burden
        dplyr::select(year, ref_mean_burden)  # , ref_max_burden, ref_min_burden
      # get difference between intervention and no-intervention simulation
      first_in_df = merge(first_in_df, noInter_df, all=TRUE)
      first_in_df$burden_averted_first_in = first_in_df$ref_mean_burden - first_in_df$mean_burden
      
      ### look at impact attributed to each intervention when it is the final intervention added (i.e. the first subtracted from a simulation with all interventions)
      last_in_df = burden_df[grepl(lastInPrefix, burden_df$scenario),]
      last_in_df$intervention = sapply(last_in_df$scenario, gsub, pattern=lastInPrefix, replacement='')
      # format the no-intervention reference data frame
      allInter_df = burden_df[burden_df$scenario == allInterName,]
      allInter_df = allInter_df %>% rename(ref_mean_burden=mean_burden) %>%  # , ref_max_burden=max_burden, ref_min_burden=min_burden
        dplyr::select(year, ref_mean_burden)  # , ref_max_burden, ref_min_burden
      # get difference between intervention and no-intervention simulation
      last_in_df = merge(last_in_df, allInter_df, all=TRUE)
      last_in_df$burden_averted_last_in = last_in_df$mean_burden - last_in_df$ref_mean_burden
      
      averted_df = merge(first_in_df[,c('intervention','year','burden_averted_first_in')], last_in_df[,c('intervention','year','burden_averted_last_in')], all=TRUE)
      write.csv(averted_df, burden_averted_filepath, row.names=FALSE)
    }
    
    if(create_plots_flag){
      gg = ggplot(averted_df, aes(x=burden_averted_first_in, y=burden_averted_last_in, color=intervention, size=as.factor(year)))+
        geom_point()+
        geom_abline()+
        ggtitle(paste0(burden_metric_name, ' averted by first-in versus last-out interventions'))+
        ylab('burden averted when intervention is last-in')+
        xlab('burden averted when intervention is first-in')+
        scale_color_manual(values = scenario_palette) +
        theme_classic()
      ggsave(paste0(sim_output_dir, '/_plots/firstIn_versus_lastOut_burden_attribution_', burden_metric, '_', age_plotted, '_',district_subset,'.png'), gg, width=8, height=7, units='in')
    }
    
    # get average between the two burden-averted calculations
    averted_df$burden_averted_mean = (averted_df$burden_averted_first_in + averted_df$burden_averted_last_in)/2
    
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # create plots of timeseries (attributed burden averted)
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    
    # create dataframe with stacked/cumulative burden averted for plotting
    stacked_df = burden_df[burden_df$scenario %in% c(noInterName, allInterName),c('year','mean_burden','scenario')]
    stacked_df$interventions_included = 'none'
    stacked_df$interventions_included[stacked_df$scenario==allInterName] = 'all'
    stacked_df = stacked_df[,c('year','mean_burden','interventions_included')]
    stacked_df$ribbon_min = stacked_df$mean_burden
    cur_inter_combo = ''
    cur_stacked_values = stacked_df[stacked_df$interventions_included == 'none', c('year','mean_burden', 'ribbon_min')]
    for(i_inter in 1:length(inter_ordered)){
      cur_inter_combo = paste0(cur_inter_combo, inter_ordered[i_inter], '_')
      # get the burden averted by current intervention
      averted_cur = averted_df[averted_df$intervention==inter_ordered[i_inter],]
      # calculate new stacked burden value and add rows to stacked_df
      new_stacked_values = merge(averted_cur, cur_stacked_values, all=TRUE)
      new_stacked_values$ribbon_min = new_stacked_values$mean_burden
      new_stacked_values$mean_burden = new_stacked_values$mean_burden - new_stacked_values$burden_averted_mean
      new_stacked_values$interventions_included = cur_inter_combo
      stacked_df = merge(stacked_df, new_stacked_values[,c('year','mean_burden', 'ribbon_min','interventions_included')], all=TRUE)
      # update current burden df
      cur_stacked_values = new_stacked_values[, c('year','mean_burden', 'ribbon_min')]
    }
    # make sure ordering is correct
    stacked_df$interventions_included = factor(stacked_df$interventions_included, levels=names(scenario_palette_stacked))
    
    # plot timeseries with burden attributions
    gg = ggplot(stacked_df, aes(x=year, y=mean_burden, color=interventions_included)) +
      geom_ribbon(aes(ymin=ribbon_min, ymax=mean_burden, fill=interventions_included), alpha=0.8, color=NA)+
      scale_fill_manual(values = scenario_palette_stacked) +
      geom_line(data=stacked_df[stacked_df$interventions_included == 'none',], aes(x=year, y=mean_burden, group=interventions_included), size=1, color='darkgrey') + 
      geom_line(data=stacked_df[stacked_df$interventions_included == cur_inter_combo,], aes(x=year, y=mean_burden, group=interventions_included), size=1, color='black') + 
      scale_color_manual(values = scenario_palette_stacked) +
      xlab('year') + 
      ylab(burden_metric_name) + 
      ylim(min(0, stacked_df$mean_burden),NA) +
      coord_cartesian(xlim=c(min_year, max_year))+
      theme_classic()+ 
      theme(legend.position = "right", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = 10))
    if(create_plots_flag){
      ggsave(paste0(sim_output_dir, '/_plots/timeseries_burden_attribution_', burden_metric, '_', age_plotted, '_',district_subset,'.png'), gg, width=10*0.8, height=4.5*0.8, units='in')
    }
    
    if(burden_colname == 'PfPR_U5'){
      returned_gg = gg
    }
  }
  return(returned_gg)
}


###################################################################################################################
# create plot showing intervention coverage from csv inputs (CM, SMC, IPTp) and/or simulation output (ITN use)
###################################################################################################################

# helper function to get weighted average coverage based on population size in included admins
get_weighted_coverage = function(coverage_df, admin_info){
  coverage_df = merge(coverage_df, admin_info[,c('admin_name','pop_size')], by='admin_name', all.x=TRUE)
  coverage_df$cov_x_pop = coverage_df$coverage * coverage_df$pop_size
  wmean_coverage = coverage_df %>% group_by(simday, duration) %>%
    summarise(sum_cov_x_pop = sum(cov_x_pop))
  total_pop = sum(admin_info$pop_size)
  wmean_coverage$coverage = wmean_coverage$sum_cov_x_pop/total_pop
  return(wmean_coverage)
}

# create main plot of all intervention inputs
plot_input_intervention_coverages = function(base_sim_input_dir, sim_output_dir, pop_filepath, intervention_coordinator_filepath, cur_admins, sim_start_year, min_year, max_year, indoor_protection_fraction, all_inter_exp_name='NGA_toPresent_allInter'){
  admin_info = read.csv(pop_filepath)
  if(cur_admins[1] == 'all'){
    cur_admins = unique(admin_info$admin_name)
    admin_info_cur = admin_info
  } else{
    admin_info_cur = admin_info[admin_info$admin_name %in% cur_admins,]
  }
  intervention_coordinator = read.csv(intervention_coordinator_filepath)
  
  # CM
  cm_input_df = read.csv(paste0(base_sim_input_dir, '/', intervention_coordinator$CM_filename[1],'.csv'))
  cm_input_df = cm_input_df[cm_input_df$admin_name %in% cur_admins,]
  cm_input_df$coverage = cm_input_df$U5_coverage
  mean_cm_coverage = get_weighted_coverage(coverage_df=cm_input_df, admin_info=admin_info_cur)
  mean_cm_coverage$year = mean_cm_coverage$simday/365 + sim_start_year
  mean_cm_coverage = mean_cm_coverage[(mean_cm_coverage$year>=min_year) & (mean_cm_coverage$year<=max_year),]
  sim_year_end_to_present = mean_cm_coverage[mean_cm_coverage$duration==-1,]
  sim_year_end_to_present$duration = round(max_year - sim_year_end_to_present$year)
  
  # SMC
  smc_input_df = read.csv(paste0(base_sim_input_dir, '/', intervention_coordinator$SMC_filename[1],'.csv'))
  smc_input_df = smc_input_df[smc_input_df$admin_name %in% cur_admins,]
  if(nrow(smc_input_df)>0){
    smc_input_df$coverage = smc_input_df$coverage_high_access_U5 * smc_input_df$high_access_U5 + smc_input_df$coverage_low_access_U5 * (1-smc_input_df$high_access_U5)
    # plot average across all four rounds
    smc_input_df = smc_input_df %>% group_by(admin_name, year) %>%
      summarise(coverage = mean(coverage),
                simday = mean(simday))
    smc_input_df$duration = 0
    # set to same simday (otherwise will have daily coverage instead of yearly when distributions occur on different days)
    smc_input_df$simday = (smc_input_df$year-sim_start_year)*365
    mean_smc_coverage = get_weighted_coverage(coverage_df=smc_input_df, admin_info=admin_info_cur)
    mean_smc_coverage$year = mean_smc_coverage$simday/365 + sim_start_year
    mean_smc_coverage = mean_smc_coverage[(mean_smc_coverage$year>=min_year) & (mean_smc_coverage$year<=max_year),]
  } else mean_smc_coverage = data.frame(year=c(min_year, max_year), coverage=c(0,0))
  
  # LLIN use
  cur_net_agg = get_intervention_use_timeseries_exp(exp_filepath = paste0(sim_output_dir, '/', all_inter_exp_name),
                                                    exp_name = all_inter_exp_name, 
                                                    cur_admins=cur_admins, pop_sizes=admin_info, min_year=min_year, max_year=max_year, indoor_protection_fraction=indoor_protection_fraction, plot_by_month=FALSE)
  
  
  # create plot
  # LLIN and CM
  g_all_inter = ggplot() +
    geom_line(data=cur_net_agg, aes(x=year, y=coverage, color="ITN use"), size=1) +
    geom_line(data=mean_cm_coverage, aes(x=year+0.5, y=coverage, color="Effective treatment with ACT"), size=1) +
    # geom_segment(data=mean_cm_coverage, aes(x=year, xend=year+duration/365, y=coverage), color=rgb(0.1,0.9,0.4), size=1)+
    geom_segment(data=sim_year_end_to_present, aes(x=year+0.5, xend=year+duration, y=coverage, color="Effective treatment with ACT"), size=1)+
    scale_color_manual("",
                       breaks=c('Effective treatment with ACT', 'ITN use','SMC (U5)'),
                       values=c('#FFC145', '#A5351F', '#0090A4'))+
    xlab('year') + 
    ylab(paste0('intervention coverage')) + 
    coord_cartesian(ylim=c(0,0.9))+
    theme_bw()
  # SMC
  if(nrow(mean_smc_coverage)>0){
    g_all_inter = g_all_inter +
      geom_line(data=mean_smc_coverage, aes(x=year, y=coverage, color="SMC (U5)"), size=1)
  }
  return(g_all_inter)
}


################################################
# add DHS/MIS U5 prevalence estimate to plot
################################################
add_state_prev_points_from_dhs = function(gg, hbhi_dir, dhs_years, state_name){
  for(dd in 1:length(dhs_years)){
    dhs_prev = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/state_prev_results_',dhs_years[dd], '.csv'))
    dhs_prev$state = gsub('-','',gsub(' ', '', toupper(dhs_prev$region)))
    
    if(toupper(state_name) %in% dhs_prev$state){
      gg = gg + 
        geom_point(data = dhs_prev[dhs_prev$state==toupper(state_name),], x=dhs_years[dd], aes(y=mic_rate), color='black', size=3)
    } else{
      warning(paste0('State name (',state_name,') did not match any name in ',dhs_years[dd],' DHS prevalence file.'))
    }
  }
  return(gg)
}
