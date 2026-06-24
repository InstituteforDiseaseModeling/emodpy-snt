# get_routine_itn_coverages.R

# script contains several functions that have different approaches for creating routine ITN input files depending on the specific data files provided by the country


library(ggplot2)
library(dplyr)


add_routine_itn_decay = function(hbhi_dir, routine_itn_coverage_filename, sim_start_year, last_sim_year, season_calib_start_year, include_trans_calib=FALSE){
  # read in admin and archetype info
  ds_pop_df_filename = paste0(hbhi_dir, '/admin_pop_archetype.csv')
  admin_pop_dataframe = read.csv(ds_pop_df_filename)
  
  # read in coverage and date information
  sim_input = read.csv(routine_itn_coverage_filename)
 
  # read in sampled retention lognormal mus and sigmas
  net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
  net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
  num_seeds = max(1,max(net_discard_decay$seed, na.rm=TRUE), na.rm=TRUE)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #####  create main simulation input file  ####
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  if(num_seeds>1){
    # combine coverage and net parameters (need separate row for each seed)
    sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=num_seeds))
    sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
    sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
  } else{
    net_discard_decay_no_seed = net_discard_decay %>% dplyr::select(-seed)
    if(nrow(net_discard_decay_no_seed) != 1){
      stop('Expected exactly one row of ITN discard decay parameters when num_seeds == 1 in transmission calibration input creation.')
    }
    sim_input_with_net_params = sim_input
    for(col_name in colnames(net_discard_decay_no_seed)){
      sim_input_with_net_params[[col_name]] = net_discard_decay_no_seed[[col_name]][1]
    }
  }
  
  # check that there are separate rows for each year (so that permethrin mortality can change in response to insecticide resistance, added in a later step)
  if(!(all(min(sim_input_with_net_params$year, na.rm=TRUE):last_sim_year %in% sim_input_with_net_params$year))) warning('The routine ITN file does not list all years separately, so if insecticide resistance changes by year, not all changes will be captured by ITN distributions. Please update input-creation scripts!')
  write.csv(sim_input_with_net_params, gsub('.csv', '_withDecay.csv', routine_itn_coverage_filename), row.names=FALSE)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  ######   calibration input files   ##########
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # difference in simday from main simulation and calibration (simday_calib = simday_main - this_value)
  if(season_calib_start_year > sim_start_year){
    simday_difference = (season_calib_start_year - sim_start_year)*365
  } else{
    warning('For setting simday in calibration input, we assume that calibration start year is later than to-present start year. However, this is not the case here. Please update code accordingly.')
  }
  if('birthday_age' %in% colnames(sim_input)){
    # get mean estimate within archetype and assign value to representative admin
    archetype_rep_admins = unique(admin_pop_dataframe$seasonality_archetype)
    arch_sim_input = sim_input %>% merge(admin_pop_dataframe[,c('admin_name','seasonality_archetype', 'pop_size')]) %>%
      mutate(pop_mult_coverage = coverage * pop_size) %>%
      group_by(seasonality_archetype, year, birthday_age) %>%
      summarise(pop_mult_coverage_sum = sum(pop_mult_coverage),
                # mean_coverage = mean(coverage),
                # median_coverage = median(coverage),
                pop_size = sum(pop_size),
                simday = mean(simday), 
                duration = mean(duration)) %>%
      ungroup() %>%
      mutate(coverage = pop_mult_coverage_sum / pop_size) %>%
      mutate(admin_name = seasonality_archetype) %>%
      mutate(simday = simday - simday_difference) %>%
      dplyr::filter(simday>=0)
  } else{
    # get mean estimate within archetype and assign value to representative admin
    archetype_rep_admins = unique(admin_pop_dataframe$seasonality_archetype)
    arch_sim_input = sim_input %>% merge(admin_pop_dataframe[,c('admin_name','seasonality_archetype', 'pop_size')]) %>%
      mutate(pop_mult_coverage = coverage * pop_size) %>%
      group_by(seasonality_archetype, year) %>%
      summarise(pop_mult_coverage_sum = sum(pop_mult_coverage),
                # mean_coverage = mean(coverage),
                # median_coverage = median(coverage),
                pop_size = sum(pop_size),
                simday = mean(simday), 
                duration = mean(duration)) %>%
      ungroup() %>%
      mutate(coverage = pop_mult_coverage_sum / pop_size) %>%
      mutate(admin_name = seasonality_archetype) %>%
      mutate(simday = simday - simday_difference) %>%
      dplyr::filter(simday>=0)
    }
                                       
  # add the mean net retention parameters
  arch_sim_input$net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[1]
  arch_sim_input$net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[1] 
  
  # check that there are separate rows for each year (so that permethrin mortality can change in response to insecticide resistance, added in a later step)
  if(!(all(min(arch_sim_input$year, na.rm=TRUE):last_sim_year %in% arch_sim_input$year))) warning('The routine ITN file does not list all years separately, so if insecticide resistance changes by year, not all changes will be captured by ITN distributions. Please update input-creation scripts!')
  write.csv(arch_sim_input, gsub('.csv', '_seasonCalib_withDecay.csv', routine_itn_coverage_filename), row.names=FALSE)
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #####  create transmission calibration input file  ####
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  if(include_trans_calib){
    # read in sampled retention lognormal mus and sigmas
    net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'))
    net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
    num_seeds = max(1,max(net_discard_decay$seed, na.rm=TRUE), na.rm=TRUE)
    
    if(num_seeds>1){
      # combine coverage and net parameters (need separate row for each seed)
      sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=num_seeds))
      sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
      sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
    } else{
      net_discard_decay_no_seed = net_discard_decay %>% dplyr::select(-seed)
      if(nrow(net_discard_decay_no_seed) != 1){
        stop('Expected exactly one row of ITN discard decay parameters when num_seeds == 1 in transmission calibration input creation.')
      }
      sim_input_with_net_params = sim_input
      for(col_name in colnames(net_discard_decay_no_seed)){
        sim_input_with_net_params[[col_name]] = net_discard_decay_no_seed[[col_name]][1]
      }
    }
    
    # check that there are separate rows for each year (so that permethrin mortality can change in response to insecticide resistance, added in a later step)
    if(!(all(min(sim_input_with_net_params$year, na.rm=TRUE):last_sim_year %in% sim_input_with_net_params$year))) warning('The routine ITN file does not list all years separately, so if insecticide resistance changes by year, not all changes will be captured by ITN distributions. Please update input-creation scripts!')
    write.csv(sim_input_with_net_params, gsub('.csv', '_transCalib_withDecay.csv', routine_itn_coverage_filename), row.names=FALSE)
  }
}





# #######################################################################################
# # from inputs on EPI ITN coverage (same value used for all districts, value from PNI report) and district-specific ANC coverage from WHO file, create two routine ITN simulation input files: 
# #   1) ITN-coverage-at-birth by admin in file format needed for simulations (birth-triggered distribution of ITNs to infants and GaveBirth triggered distribution to mothers) - based on ANC coverage
# #   2) ITN-coverage-at-18-months from EPI ("Expanded Programme on Immunizations") ITN distributions. Note that we often may inflate the coverage and probability of use given coverage for 18-month-old children to compensate for the fact that in practice, these nets would likely also protect some f the siblings of the children receiving the nets.
# #          We assume the same EPI coverage across all years and all districts
# 
# 
# # provides nets on birthdays, assumes the same coverage across all years and all admins
# create_epi_itn_inputs = function(hbhi_dir, birthday_ages=c(1,3), coverage_each_birthday=c(0.915, 0.732), epi_itn_use_given_access=0.9, 
#                                  last_sim_year=2021, sim_start_year=2010, season_calib_start_year=2011, trans_calib_start_year=2010, include_trans_calib=FALSE){
#   
#   ds_pop_df_filename = paste0(hbhi_dir, '/admin_pop_archetype.csv')
#   admin_pop_dataframe = read.csv(ds_pop_df_filename)
#   admin_names = admin_pop_dataframe$admin_name
#   
#   # iterate through birthday ages, adding coverages for each
#   for(bb in 1:length(birthday_ages)){
#     sim_input_cur = data.frame(admin_name = admin_names)
#     sim_input_cur$coverage = coverage_each_birthday[bb] * epi_itn_use_given_access
#     sim_input_cur$birthday_age = birthday_ages[bb]
#     if(bb==1){
#       sim_input = sim_input_cur
#     } else{
#       sim_input = rbind(sim_input, sim_input_cur)
#     }
#   }  
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   #####  create main simulation input file  ####
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   
#   # read in sampled retention lognormal mus and sigmas
#   net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
#   net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#   # TODO: make sure that simulations are using seeds, otherwise if num seeds<=1, don't include seed column
#   # combine coverage and net parameters (need separate row for each seed)
#   sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=max(net_discard_decay$seed)))
#   sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
#   sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
#   
#   # create separate rows for each year (so that permethrin mortality can change in response to insecticide resistance, added in a later step)
#   sim_input_with_year = sim_input_with_net_params %>% slice(rep(1:n(), each=length(sim_start_year:last_sim_year)))
#   sim_input_with_year$year = rep(sim_start_year:last_sim_year, times=nrow(sim_input_with_net_params))
#   sim_input_with_year$simday = (sim_input_with_year$year - sim_start_year) * 365 + 1
#   sim_input_with_year$duration = 365
#   
#   # write to csv
#   ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
#   ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
#   write.csv(sim_input_with_year, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/epi_itn_use_coverages_2010_2020.csv'), row.names=FALSE)
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   ######   calibration input files   ##########
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   # get archetype representative admin names
#   archetype_rep_admins = unique(admin_pop_dataframe$seasonality_archetype)
#   arch_sim_input = sim_input[sim_input$admin_name %in% archetype_rep_admins,] %>%
#     mutate(admin_name = seasonality_archetype)
#   
#   # add the mean net retention parameters
#   arch_sim_input$net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[1]
#   arch_sim_input$net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[1] 
#   
#   # create separate rows for each year (so that permethrin mortality can change in response to insecticide resistance, added in a later step)
#   arch_sim_input_with_year = arch_sim_input %>% slice(rep(1:n(), each=length(season_calib_start_year:last_sim_year)))
#   arch_sim_input_with_year$year = rep(season_calib_start_year:last_sim_year, times=nrow(arch_sim_input))
#   arch_sim_input_with_year$simday = (arch_sim_input_with_year$year - season_calib_start_year) * 365 + 1
#   arch_sim_input_with_year$duration = 365
#   
#   write.csv(arch_sim_input_with_year, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/epi_itn_use_coverages_season_calib.csv'), row.names=FALSE)
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   #####  create transmission calibration input file  ####
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   if(include_trans_calib){
#     # read in sampled retention lognormal mus and sigmas
#     net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'))
#     net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     # TODO: make sure that simulations are using seeds, otherwise if num seeds<=1, don't include seed column
#     # combine coverage and net parameters (need separate row for each seed)
#     sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=max(net_discard_decay$seed)))
#     sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
#     sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
#     
#     # create separate rows for each year (so that permethrin mortality can change in response to insecticide resistance, added in a later step)
#     sim_input_with_year = sim_input_with_net_params %>% slice(rep(1:n(), each=length(trans_calib_start_year:last_sim_year)))
#     sim_input_with_year$year = rep(trans_calib_start_year:last_sim_year, times=nrow(sim_input_with_net_params))
#     sim_input_with_year$simday = (sim_input_with_year$year - trans_calib_start_year) * 365 + 1
#     sim_input_with_year$duration = 365
#     
#     # write to csv
#     write.csv(sim_input_with_year, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/epi_itn_use_coverages_trans_calib.csv'), row.names=FALSE)
#   }
# }
# # # CHECKING: does the access rate for the final year make sense given the number of births in that year and the number of nets 
# # # do these numbers align well with number of nets allocated to ANC compared with the number of births in Burundi?
# # birth_rate = 0.039  # per person per year
# # national_pop_size = 11500000
# # annual_births = national_pop_size * birth_rate
# # anc_nets_reported_2017 = 447000  # from Table 3 in PMI MOP 2018
# # anc_access = anc_nets_distributed_2017 / annual_births
# 
# 
# 
# 
# 
# 
# create_anc_itn_inputs_from_WHO = function(anc_itn_access_filename, hbhi_dir, ds_pop_df_filename, max_anc_itn_access=0.97, anc_itn_use_given_access=0.9, sim_start_year=2010, calib_anc_itn_year=2016, num_samples=0){
#   
#   data_file = read.csv(anc_itn_access_filename)
#   # head(data_file)
#   data_file = data_file[!is.na(data_file$cov_llins_anc),]  # remove NA rows
#   
#   admin_pop_dataframe = read.csv(ds_pop_df_filename)
#   
#   # what is the first and final year in the ANC LLIN file
#   first_year = min(data_file$year)
#   last_year = max(data_file$year)
#   # assume that the final distribution is the rate at which ANCs continue to be delivered even if there are no data on following years?
#   continue_coverage = TRUE   # ( this will set the duration of the final year to be -1)
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   #####  create main simulation input file  ####
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   if(all(data_file$adm2 %in% admin_pop_dataframe$admin_name)){  # if all the names match up perfectly, don't even need to merge / correct
#     sim_input = data_file[,which(colnames(data_file) %in% c('adm2', 'year', 'cov_llins_anc'))]
#     colnames(sim_input)[colnames(sim_input) == 'adm2'] = 'admin_name'
#     if(mean(sim_input$cov_llins_anc, na.rm=TRUE)>1){ # reported in terms of percentages
#       sim_input$access = sim_input$cov_llins_anc / 100
#     } else{
#       sim_input$access = sim_input$cov_llins_anc
#     }
#     sim_input$simday = (sim_input$year - sim_start_year) * 365
#     sim_input$duration = 365
#     
#     # check whether all of the admins are included in the data file for the final year (otherwise might need to use the second-to-last year for some)
#     if (all(sort(sim_input$admin_name[sim_input$year == last_year]) == sort(admin_pop_dataframe$admin_name))){
#       sim_input$duration[sim_input$year == last_year] = -1  # these coverages will continue if the same simulation is run longer than the last year in the ANC file (but not across serialization)
#     } else warning('The admins in the final year of the ANC LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     
#     
#     # according to the PMI MOP, ANC LLINs 'continued' in 2016 (implying they were also distributed before 2016) and 'during this period' (after 2009), 'USAID support continued to focus on the procurement and distribution of ITNs for pregnant women at ANC clinics and children during routine immunization services.'
#     #   lacking any other information, extend ANC access backwards assuming linear decline to 0 coverage in 2009
#     # check whether all of the admins are included in the data file for the first year (otherwise might need to use the second year for some)
#     if(2010<first_year){
#       if (all(sort(sim_input$admin_name[sim_input$year == first_year]) == sort(admin_pop_dataframe$admin_name))){
#         years_projected = 2009:(first_year-1)
#         for(yy in max(2010, sim_start_year):(first_year-1)){
#           yy_df = sim_input[sim_input$year == first_year, ]  # set equal to the first year ANC recorded (access will then be multiplied to get linear decrease from this starting point back in time)
#           yy_df$access = yy_df$access * ( 1 - (first_year - yy) / length(2009:(first_year-1)) )
#           yy_df$year = yy
#           yy_df$duration = 365
#           yy_df$simday = (yy - sim_start_year) * 365
#           sim_input = rbind(sim_input, yy_df)
#         }
#       } else warning('The admins in the first year of the ANC LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     }
#     
#     # set a maximum reasonable coverage
#     sim_input$access = sapply(sim_input$access, min, max_anc_itn_access)
#     
#     # convert all values from net access to net use (at peak use season)
#     sim_input$coverage = sim_input$access * anc_itn_use_given_access
#     sim_input = sim_input[,-which(colnames(sim_input) == 'cov_llins_anc')]
#     
#     # read in sampled retention lognormal mus and sigmas
#     net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
#     net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     # TODO: make sure that simulations are using seeds, otherwise if num seeds<=1, don't include seed column
#     # combine coverage and net parameters (need separate row for each seed)
#     sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=max(net_discard_decay$seed)))
#     sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
#     sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
#     # remove seed if there is only one
#     if(length(unique(sim_input_with_net_params$seed))==1){
#       sim_input_with_net_params = sim_input_with_net_params[,-which(colnames(sim_input_with_net_params)=='seed')]
#     }
#     
#     # write to csv
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
#     write.csv(sim_input_with_net_params, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/anc_itn_use_coverages_2010_toPresent.csv'), row.names=FALSE)
#     
#     
#     # create plots showing ANC use coverages through time
#     if(!dir.exists(paste0(hbhi_dir,'/simulation_inputs/plots'))) dir.create(paste0(hbhi_dir,'/simulation_inputs/plots'))
#     gg = ggplot(sim_input, aes(x=year, y=coverage, color=admin_name)) + 
#       geom_line() + 
#       ylab('ANC LLIN use coverage') +
#       theme_classic()+
#       theme(legend.position='none')
#     ggsave(filename=paste0(hbhi_dir,'/simulation_inputs/plots/itn_anc_coverage_through_time.png'), gg, width=4, height=3, units='in')
#     
#   } else warning('Not all admin names from the ANC LLIN file match perfectly with the reference population admin file. If this warning appears, it will be necessary to write code that can handle this.')
#   
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   ######   calibration input files   ##########
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   # in each archetype, use a single coverage across all years calculated by 
#   #   1) averaging coverage within each admin across years reported in the routine dataset, 
#   #   2) taking the population-weighted average across all admins in each archetype
#   
#   # get average across years included in the routine dataset
#   admin_averages = sim_input[sim_input$year >= first_year,] %>% group_by(admin_name) %>% summarise(mean_coverage = mean(coverage))
#   
#   # get population-size weighted coverages for archetypes used in seasonality calibrations
#   admin_pop_dataframe_coverage = merge(admin_pop_dataframe, admin_averages, by='admin_name')
#   admin_pop_dataframe_coverage$product_cov_pop = admin_pop_dataframe_coverage$mean_coverage * admin_pop_dataframe_coverage$pop_size
#   arch_sums = admin_pop_dataframe_coverage %>%  
#     group_by(seasonality_archetype) %>%
#     summarise(sum_pop = sum(pop_size),
#               sum_product = sum(product_cov_pop))
#   arch_sums$coverage = arch_sums$sum_product / arch_sums$sum_pop
#   arch_anc_itn = arch_sums[,which(colnames(arch_sums) %in% c('seasonality_archetype','coverage'))]
#   colnames(arch_anc_itn)[colnames(arch_anc_itn) == 'seasonality_archetype'] = 'admin_name'
#   arch_anc_itn$duration = -1
#   arch_anc_itn$simday = 0
#   arch_anc_itn$year = calib_anc_itn_year
#   # add the mean net retention parameters
#   arch_anc_itn$net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[1]
#   arch_anc_itn$net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[1] 
#   
#   
#   write.csv(arch_anc_itn, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/anc_itn_use_coverages_season_calib.csv'), row.names=FALSE)
#   
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   #####  create transmission calibration input file  ####
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   if(all(data_file$adm2 %in% admin_pop_dataframe$admin_name)){  # if all the names match up perfectly, don't even need to merge / correct
#     sim_input = data_file[,which(colnames(data_file) %in% c('adm2', 'year', 'cov_llins_anc'))]
#     colnames(sim_input)[colnames(sim_input) == 'adm2'] = 'admin_name'
#     sim_input$access = sim_input$cov_llins_anc / 100
#     sim_input$simday = (sim_input$year - sim_start_year) * 365
#     sim_input$duration = 365
#     
#     # check whether all of the admins are included in the data file for the final year (otherwise might need to use the second-to-last year for some)
#     if (all(sort(sim_input$admin_name[sim_input$year == last_year]) == sort(admin_pop_dataframe$admin_name))){
#       sim_input$duration[sim_input$year == last_year] = -1  # these coverages will continue if the same simulation is run longer than the last year in the ANC file (but not across serialization)
#     } else warning('The admins in the final year of the ANC LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     
#     
#     # according to the PMI MOP, ANC LLINs 'continued' in 2016 (implying they were also distributed before 2016) and 'during this period' (after 2009), 'USAID support continued to focus on the procurement and distribution of ITNs for pregnant women at ANC clinics and children during routine immunization services.'
#     #   lacking any other information, extend ANC access backwards assuming linear decline to 0 coverage in 2009
#     # check whether all of the admins are included in the data file for the first year (otherwise might need to use the second year for some)
#     if(2010<first_year){
#       if (all(sort(sim_input$admin_name[sim_input$year == first_year]) == sort(admin_pop_dataframe$admin_name))){
#         years_projected = 2009:(first_year-1)
#         for(yy in max(2010, sim_start_year):(first_year-1)){
#           yy_df = sim_input[sim_input$year == first_year, ]  # set equal to the first year ANC recorded (access will then be multiplied to get linear decrease from this starting point back in time)
#           yy_df$access = yy_df$access * ( 1 - (first_year - yy) / length(2009:(first_year-1)) )
#           yy_df$year = yy
#           yy_df$duration = 365
#           yy_df$simday = (yy - sim_start_year) * 365
#           sim_input = rbind(sim_input, yy_df)
#         }
#       } else warning('The admins in the first year of the ANC LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     }
#     
#     # set a maximum reasonable coverage
#     sim_input$access = sapply(sim_input$access, min, max_anc_itn_access)
#     
#     # convert all values from net access to net use (at peak use season)
#     sim_input$coverage = sim_input$access * anc_itn_use_given_access
#     sim_input = sim_input[,-which(colnames(sim_input) == 'cov_llins_anc')]
#     
#     if(num_samples>0){ # if there are uncertainty samples, use quantiles, otherwise can use main file which contains the mean value
#       # read in sampled retention lognormal mus and sigmas
#       net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'))
#       net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     } else{
#       # read in sampled retention lognormal mus and sigmas
#       net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
#       net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     }
#     # TODO: make sure that simulations are using seeds, otherwise if num seeds<=1, don't include seed column
#     # combine coverage and net parameters (need separate row for each seed)
#     sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=max(net_discard_decay$seed)))
#     sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
#     sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
#     # remove seed if there is only one
#     if(length(unique(sim_input_with_net_params$seed))==1){
#       sim_input_with_net_params = sim_input_with_net_params[,-which(colnames(sim_input_with_net_params)=='seed')]
#     }
#     
#     # write to csv
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
#     write.csv(sim_input_with_net_params, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/anc_itn_use_coverages_trans_calib.csv'), row.names=FALSE)
#     
#     
#     # # create plots showing ANC use coverages through time
#     # ggplot(sim_input, aes(x=year, y=coverage, color=admin_name)) + 
#     #   geom_line() + 
#     #   ylab('ANC LLIN use coverage') +
#     #   theme_classic()
#     
#   } else warning('Not all admin names from the ANC LLIN file match perfectly with the reference population admin file. If this warning appears, it will be necessary to write code that can handle this.')
# }
# 
# 
# 
# 
# 
# # add EPI nets for both the child receiving the EPI ITN and sometimes a sibling
# create_epi_itn_inputs_from_WHO = function(epi_itn_access_filename, hbhi_dir, ds_pop_df_filename, birthday_ages, coverage_age_rel_age1, max_epi_itn_access, epi_itn_use_given_access, 
#                                           last_sim_year, sim_start_year=2010, calib_epi_itn_year=2016, num_samples=0){
#   
#   data_file = read.csv(epi_itn_access_filename)
#   data_file = data_file[!is.na(data_file$coverage),]  # remove NA rows
#   
#   admin_pop_dataframe = read.csv(ds_pop_df_filename)
#   
#   # what is the first and final year in the EPI LLIN file
#   first_year = min(data_file$year)
#   last_year = max(data_file$year)
#   # assume that the final distribution is the rate at which EPIs continue to be delivered even if there are no data on following years?
#   continue_coverage = TRUE   # ( this will set the duration of the final year to be -1)
#   
#   # update column names and units (from percentage to decimal, if necessary)
#   colnames(data_file)[colnames(data_file) == 'adm2'] = 'admin_name'
#   if(mean(data_file$coverage, na.rm=TRUE)>5){ # proxy for determining if reported in terms of percentages
#     data_file$coverage = data_file$coverage / 100
#   } 
#   
#   # adjust access coverage to have maximum rate of individuals receiving EPI ITNs
#   data_file$coverage = sapply(data_file$coverage, min, max_epi_itn_access)
#   # convert all values from net access to net use (at peak use season)
#   data_file$coverage = data_file$coverage * epi_itn_use_given_access
#   # iterate through birthday ages, adding coverages for each
#   data_file$birthday_age = birthday_ages[1]
#   if(length(birthday_ages)>1){
#     data_file_original = data_file
#     for(bb in 2:length(birthday_ages)){
#       data_file_cur_age = data_file_original
#       data_file_cur_age$coverage = data_file_original$coverage * coverage_age_rel_age1[bb]
#       data_file_cur_age$birthday_age = birthday_ages[bb]
#       data_file = rbind(data_file, data_file_cur_age)
#     }  
#   }
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   #####  create main simulation input file  ####
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   if(all(data_file$admin_name %in% admin_pop_dataframe$admin_name)){  # if all the names match up perfectly, don't even need to merge / correct
#     sim_input = data_file[,which(colnames(data_file) %in% c('admin_name', 'year', 'coverage', 'birthday_age'))]
#     sim_input$simday = (sim_input$year - sim_start_year) * 365
#     sim_input$duration = 365
#     
#     # check whether all of the admins are included in the data file for the final year (otherwise might need to use the second-to-last year for some)
#     if (all(sort(sim_input$admin_name[(sim_input$year == last_year) & (sim_input$birthday_age == birthday_ages[1])]) == sort(admin_pop_dataframe$admin_name))){
#       sim_input$duration[sim_input$year == last_year] = -1  # these coverages will continue if the same simulation is run longer than the last year in the EPI file (but not across serialization)
#     } else warning('The admins in the final year of the EPI LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     
#     # extend EPI ITN access backwards assuming linear decline to 0 coverage in 2009
#     # check whether all of the admins are included in the data file for the first year (otherwise might need to use the second year for some)
#     if(2010<first_year){
#       if (all(sort(sim_input$admin_name[(sim_input$year == first_year) & (sim_input$birthday_age == birthday_ages[1])]) == sort(admin_pop_dataframe$admin_name))){
#         years_projected = 2009:(first_year-1)
#         for(yy in max(2010, sim_start_year):(first_year-1)){
#           yy_df = sim_input[sim_input$year == first_year, ]  # set equal to the first year EPI recorded (access will then be multiplied to get linear decrease from this starting point back in time)
#           yy_df$coverage = yy_df$coverage * ( 1 - (first_year - yy) / length(2009:(first_year-1)) )
#           yy_df$year = yy
#           yy_df$duration = 365
#           yy_df$simday = (yy - sim_start_year) * 365
#           sim_input = rbind(sim_input, yy_df)
#         }
#       } else warning('The admins in the first year of the EPI LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     }
#         
#     # read in sampled retention lognormal mus and sigmas
#     net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
#     net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     # TODO: make sure that simulations are using seeds, otherwise if num seeds<=1, don't include seed column
#     # combine coverage and net parameters (need separate row for each seed)
#     sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=max(net_discard_decay$seed)))
#     sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
#     sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
#     # remove seed if there is only one
#     if(length(unique(sim_input_with_net_params$seed))==1){
#       sim_input_with_net_params = sim_input_with_net_params[,-which(colnames(sim_input_with_net_params)=='seed')]
#     }
#     
#     # write to csv
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
#     write.csv(sim_input_with_net_params, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/epi_itn_use_coverages_2010_toPresent.csv'), row.names=FALSE)
#     
#     
#     # create plots showing EPI use coverages through time
#     if(!dir.exists(paste0(hbhi_dir,'/simulation_inputs/plots'))) dir.create(paste0(hbhi_dir,'/simulation_inputs/plots'))
#     gg = ggplot(sim_input, aes(x=year, y=coverage, color=admin_name, linetype=factor(birthday_age))) + 
#       geom_line() + 
#       ylab('EPI LLIN use coverage') +
#       theme_classic()+
#       theme(legend.position='none')
#     ggsave(filename=paste0(hbhi_dir,'/simulation_inputs/plots/itn_epi_coverage_through_time.png'), gg, width=4, height=3, units='in')
#     
#   } else warning('Not all admin names from the EPI LLIN file match perfectly with the reference population admin file. If this warning appears, it will be necessary to write code that can handle this.')
#   
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   ######   calibration input files   ##########
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   # in each archetype, use a single coverage across all years calculated by 
#   #   1) averaging coverage within each admin across years reported in the routine dataset, 
#   #   2) taking the population-weighted average across all admins in each archetype
#   
#   # get average across years included in the routine dataset
#   admin_averages = sim_input[sim_input$year >= first_year,] %>% group_by(admin_name, birthday_age) %>% summarise(mean_coverage = mean(coverage))
#   
#   # get population-size weighted coverages for archetypes used in seasonality calibrations
#   admin_pop_dataframe_coverage = merge(admin_pop_dataframe, admin_averages, by='admin_name')
#   admin_pop_dataframe_coverage$product_cov_pop = admin_pop_dataframe_coverage$mean_coverage * admin_pop_dataframe_coverage$pop_size
#   arch_sums = admin_pop_dataframe_coverage %>%  
#     group_by(seasonality_archetype, birthday_age) %>%
#     summarise(sum_pop = sum(pop_size),
#               sum_product = sum(product_cov_pop))
#   arch_sums$coverage = arch_sums$sum_product / arch_sums$sum_pop
#   arch_epi_itn = arch_sums[,which(colnames(arch_sums) %in% c('seasonality_archetype', 'coverage', 'birthday_age'))]
#   colnames(arch_epi_itn)[colnames(arch_epi_itn) == 'seasonality_archetype'] = 'admin_name'
#   arch_epi_itn$duration = -1
#   arch_epi_itn$simday = 0
#   arch_epi_itn$year = calib_epi_itn_year
#   # add the mean net retention parameters
#   arch_epi_itn$net_life_lognormal_mu = net_discard_decay$net_life_lognormal_mu[1]
#   arch_epi_itn$net_life_lognormal_sigma = net_discard_decay$net_life_lognormal_sigma[1] 
#   
#   
#   write.csv(arch_epi_itn, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/epi_itn_use_coverages_season_calib.csv'), row.names=FALSE)
#   
#   
#   
#   
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   #####  create transmission calibration input file  ####
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   if(all(data_file$admin_name %in% admin_pop_dataframe$admin_name)){  # if all the names match up perfectly, don't even need to merge / correct
#     sim_input = data_file[,which(colnames(data_file) %in% c('admin_name', 'year', 'coverage', 'birthday_age'))]
#     sim_input$simday = (sim_input$year - sim_start_year) * 365
#     sim_input$duration = 365
#     
#     # check whether all of the admins are included in the data file for the final year (otherwise might need to use the second-to-last year for some)
#     if (all(sort(sim_input$admin_name[(sim_input$year == last_year) & (sim_input$birthday_age == birthday_ages[1])]) == sort(admin_pop_dataframe$admin_name))){
#       sim_input$duration[sim_input$year == last_year] = -1  # these coverages will continue if the same simulation is run longer than the last year in the EPI file (but not across serialization)
#     } else warning('The admins in the final year of the EPI LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     
#     
#     # extend EPI access backwards assuming linear decline to 0 coverage in 2009
#     # check whether all of the admins are included in the data file for the first year (otherwise might need to use the second year for some)
#     if(2010<first_year){
#       if (all(sort(sim_input$admin_name[(sim_input$year == first_year) & (sim_input$birthday_age == birthday_ages[1])]) == sort(admin_pop_dataframe$admin_name))){
#         years_projected = 2009:(first_year-1)
#         for(yy in max(2010, sim_start_year):(first_year-1)){
#           yy_df = sim_input[sim_input$year == first_year, ]  # set equal to the first year EPI recorded (access will then be multiplied to get linear decrease from this starting point back in time)
#           yy_df$coverage = yy_df$coverage * ( 1 - (first_year - yy) / length(2009:(first_year-1)) )
#           yy_df$year = yy
#           yy_df$duration = 365
#           yy_df$simday = (yy - sim_start_year) * 365
#           sim_input = rbind(sim_input, yy_df)
#         }
#       } else warning('The admins in the first year of the EPI LLIN file are not a perfect match to all admins listed in the population file. If this warning appears, it will be necessarry to write code that can handle this.')
#     }
#     
#     
#     if(num_samples>0){ # if there are uncertainty samples, use quantiles, otherwise can use main file which contains the mean value
#       # read in sampled retention lognormal mus and sigmas
#       net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_quantile_discard_decay_params.csv'))
#       net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     } else{
#       # read in sampled retention lognormal mus and sigmas
#       net_discard_decay = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))
#       net_discard_decay = net_discard_decay[,which(colnames(net_discard_decay) %in% c('seed', 'net_life_lognormal_mu', 'net_life_lognormal_sigma'))]
#     }
#     # TODO: make sure that simulations are using seeds, otherwise if num seeds<=1, don't include seed column
#     # combine coverage and net parameters (need separate row for each seed)
#     sim_input_with_seed = sim_input %>% slice(rep(1:n(), each=max(net_discard_decay$seed)))
#     sim_input_with_seed$seed = rep(1:max(net_discard_decay$seed), times=nrow(sim_input))
#     sim_input_with_net_params = merge(sim_input_with_seed, net_discard_decay, by=c('seed'), all.x=TRUE)
#     # remove seed if there is only one
#     if(length(unique(sim_input_with_net_params$seed))==1){
#       sim_input_with_net_params = sim_input_with_net_params[,-which(colnames(sim_input_with_net_params)=='seed')]
#     }
#     
#     # write to csv
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files')), FALSE)
#     ifelse(!dir.exists(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), dir.create(paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage')), FALSE)
#     write.csv(sim_input_with_net_params, paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/epi_itn_use_coverages_trans_calib.csv'), row.names=FALSE)
#     
#     
#     # # create plots showing EPI use coverages through time
#     # ggplot(sim_input, aes(x=year, y=coverage, color=admin_name)) + 
#     #   geom_line() + 
#     #   ylab('EPI LLIN use coverage') +
#     #   theme_classic()
#     
#   } else warning('Not all admin names from the EPI LLIN file match perfectly with the reference population admin file. If this warning appears, it will be necessary to write code that can handle this.')
# }






