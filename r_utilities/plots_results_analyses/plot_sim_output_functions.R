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



# for a plot with horizontal gridlines but not the box around the plot
theme_gridlines_no_box <- function(...) {
  theme_bw(...) +
    theme(
      panel.border        = element_blank(),
      axis.line           = element_line(colour = "black"),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor.x  = element_blank()
    )
}

####################################################################################
# barplots for burden relative to BAU: percent reduction
####################################################################################

plot_relative_burden_barplots = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                         barplot_start_year, barplot_end_year, 
                                         pyr, chw_cov,
                                         scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, separate_plots_flag=FALSE, standard_max_y = 0.1, show_error_bar=TRUE, align_seeds=TRUE,
                                         include_to_present=TRUE, burden_metric_subset=c(), scenario_barfill=NULL){
  # scenario_barfill (optional named vector scenario -> ggpattern pattern, e.g. 'stripe'/'none') hatches
  # selected bars; NULL = plain solid bars (default, so existing callers are unaffected).
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')#, 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')#, 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }

  # first comparison name is to-present (skip it), second is BAU (use as reference), comparison scenarios start at the third index
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2  
  }
  # iterate through comparison scenarios, calculating the burden reduction of all metrics relative to BAU (seedwise comparisons, so one output for each run). Combine all scenario reductions into a dataframe (each scenario set in separate rows)
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, 
                                             start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    # only save relevant columns for plotting
    relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    }else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_names[comparison_start_index:length(scenario_names)])
  
  # guard against non-finite % reductions: a near-zero reference makes (ref-comp)/ref blow up to
  # +/-Inf, which (left unhandled) makes the y-limits non-finite and crashes axis-break computation
  # (Error in seq.default: invalid '(to - from)/by'). Report which scenario x metric were affected
  # (so the offending experiment / cumulativeBurden cache can be checked/deleted), then set them to NA.
  for(cc in burden_colnames){
    bad = !is.finite(relative_burden_all_df[[cc]])
    if(any(bad)){
      for(sc in unique(as.character(relative_burden_all_df$scenario[bad]))){
        n = sum(bad & as.character(relative_burden_all_df$scenario) == sc)
        message(sprintf('plot_relative_burden_barplots: %d non-finite %% reduction value(s) set to NA  [metric: %s | scenario: %s]', n, cc, sc))
      }
      relative_burden_all_df[[cc]][bad] = NA
    }
  }
  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_y = 0
  finite_vals = unlist(relative_burden_all_df[, burden_colnames, drop=FALSE], use.names=FALSE)
  finite_vals = finite_vals[is.finite(finite_vals)]
  cur_min = if(length(finite_vals)) min(finite_vals) else 0
  cur_max = if(length(finite_vals)) max(finite_vals) else standard_max_y
  if(cur_min < standard_min_y) standard_min_y = cur_min
  if(cur_max > standard_max_y) standard_max_y = cur_max
  
  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name = burden_metric_names[bb]
    select_col_names = c(current_burden_name, 'scenario')
    # get mean, min, and max among all runs for this burden metric
    rel_burden_agg = as.data.frame(relative_burden_all_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name), na.rm=TRUE),
                       max_rel = max(get(current_burden_name), na.rm=TRUE),
                       min_rel = min(get(current_burden_name), na.rm=TRUE))

    bar_layer = if(is.null(scenario_barfill)){
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity')
    } else {
      ggpattern::geom_bar_pattern(aes(x=scenario, y=mean_rel, fill=scenario, pattern=scenario), stat='identity',
                                  pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                  pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5)
    }
    gg_list[[bb]] = ggplot(rel_burden_agg) +
      bar_layer +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
      ylab('Percent reduction') +
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values = scenario_palette) +
      (if(!is.null(scenario_barfill)) scale_pattern_manual(values=scenario_barfill)) +
      (if(!is.null(scenario_barfill)) guides(pattern=guide_legend(override.aes=list(fill='white')), fill=guide_legend(override.aes=list(pattern='none')))) +
      # theme_classic()+
      theme_gridlines_no_box()+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size),
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg_list[[bb]] = gg_list[[bb]] +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour="black", alpha=0.9, size=1) 
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] + 
        ylab('Percent reduction in burden \n ((Current - Plan) / Current) * 100') + 
        theme(legend.position='none', plot.title = element_blank(), text=element_text(size =separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/','barplot_percent_reduction_', burden_metric_name,'_',district_subset,'.png'), separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }
  # for each burden type, 
  # get mean, min, and max among all runs for each burden metric, each saved as a separate column
  # create barplot for each burden type (using columns of dataframe, separate bar for each scenario)
  
  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  # remove legend from main plots
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position = "none")  + theme(text = element_text(size = text_size))   
  }
  
  if(save_plots){
    # arrangeGrob (not grid.arrange): build the grob WITHOUT drawing to the active device. grid.arrange
    # draws immediately to the current device, and on some interactive devices ggpattern's crosshatch
    # renderer crashes (Error in seq.default: invalid '(to - from)/by'); ggsave then draws to PNG safely.
    gg_saved = arrangeGrob(grobs = gg_list[-1], layout_matrix = matrix(c(1:(length(burden_colnames))), nrow=2, byrow=FALSE))
    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_percent_reduction_burden_', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg_saved, dpi=600, width=14, height=7, units='in')
  }

  # ----- combine all burden plots ----- #
  # gg = grid.arrange(grobs = gg_list, layout_matrix = matrix(c(1,1,2:(length(burden_colnames)+1)), ncol=2, byrow=TRUE))
  gg = arrangeGrob(grobs = gg_list, layout_matrix = rbind(matrix(rep(1, length(burden_colnames)/2), nrow=1), matrix(2:(length(burden_colnames)+1), nrow=2, byrow=FALSE)))  # arrangeGrob, not grid.arrange (see note above) -- caller ggsave's the returned grob

  return(gg)
}


####################################################################################
# Grouped barplots: scenarios as grouped bars within each burden metric, U5 over all-ages
####################################################################################
# Companion to plot_relative_burden_barplots(). When only a few scenarios are being
# compared, the per-metric panel layout looks sparse; this version instead returns just
# two panels (U5 on top, all-age on bottom), where the x-axis groups are burden metrics
# and bars within each group are colored by scenario (dodged).
# Signature matches plot_relative_burden_barplots() for easy swap at call sites.
plot_relative_burden_grouped_barplots = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                                                 barplot_start_year, barplot_end_year,
                                                 pyr, chw_cov,
                                                 scenario_names, experiment_names, scenario_palette,
                                                 LLIN2y_flag=FALSE, overwrite_files=FALSE, separate_plots_flag=FALSE,
                                                 standard_max_y=0.1, show_error_bar=TRUE, align_seeds=TRUE,
                                                 include_to_present=TRUE, burden_metric_subset=c(),
                                                 font_scale=1, legend_scale=0.6, scenario_barfill=NULL){
  # scenario_barfill (optional named vector scenario -> ggpattern pattern, e.g. 'stripe'/'none')
  # hashes selected bars (e.g. scaled-up routine-ITN variants); NULL = plain bars (default).
  # font_scale multiplies all in-plot text; legend_scale shrinks the legend relative to body
  # text (legends were getting cut off). Defaults preserve sizing for existing callers except
  # the now-smaller legend. (ggsave dimensions are set at the call site in the main script.)
  admin_pop = read.csv(pop_filepath)

  # burden metrics (same set as plot_relative_burden_barplots)
  burden_metrics      = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')
  burden_colnames     = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  if(length(burden_metric_subset) >= 1){
    keep_idx = which(burden_metrics %in% burden_metric_subset)
    burden_metrics      = burden_metrics[keep_idx]
    burden_colnames     = burden_colnames[keep_idx]
    burden_metric_names = burden_metric_names[keep_idx]
  }

  # determine reference scenario (same logic as plot_relative_burden_barplots)
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2
  }

  # compute % reductions for each comparison scenario relative to reference
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir,
                                             reference_experiment_name=reference_experiment_name,
                                             comparison_experiment_name=comparison_experiment_name,
                                             comparison_scenario_name=comparison_scenario_name,
                                             start_year=barplot_start_year, end_year=barplot_end_year,
                                             admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins,
                                             LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    relative_burden_df = relative_burden_df[, which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    } else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }

  # keep scenario factor order matching scenario_names
  scenario_order = scenario_names[comparison_start_index:length(scenario_names)]
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_order)

  # reshape to long form: scenario, burden_col, rel_reduction
  long_df = relative_burden_all_df %>%
    tidyr::pivot_longer(cols = all_of(burden_colnames), names_to = 'burden_col', values_to = 'rel_reduction')

  # attach metric family name (without age suffix) and age group
  col_to_metric_name = setNames(burden_metric_names, burden_colnames)
  long_df$burden_metric_name = col_to_metric_name[long_df$burden_col]
  long_df$age_group     = ifelse(grepl('\\(U5\\)', long_df$burden_metric_name), 'U5', 'all ages')
  long_df$metric_family = gsub(' \\(U5\\)| \\(all ages\\)', '', long_df$burden_metric_name)
  # preserve order of metric families as they first appear in burden_metric_names
  metric_family_levels = unique(gsub(' \\(U5\\)| \\(all ages\\)', '', burden_metric_names))
  long_df$metric_family = factor(long_df$metric_family, levels=metric_family_levels)

  # aggregate across runs: mean / min / max per scenario x metric x age group
  agg_df = long_df %>%
    dplyr::group_by(scenario, age_group, metric_family) %>%
    dplyr::summarise(mean_rel = mean(rel_reduction),
                     min_rel  = min(rel_reduction),
                     max_rel  = max(rel_reduction),
                     .groups = 'drop')

  # y-axis range: widen if data exceeds the standard window
  standard_min_y_use = min(0, min(agg_df$mean_rel))
  standard_max_y_use = max(standard_max_y, max(agg_df$mean_rel))

  # one grouped barplot per age group
  age_groups = c('U5', 'all ages')
  age_titles = c('U5' = 'Under-5', 'all ages' = 'All ages')
  gg_list = list()
  for(aa in seq_along(age_groups)){
    cur_age = age_groups[aa]
    cur_df  = agg_df[agg_df$age_group == cur_age, ]
    # bar layer: plain bars, or hashed bars when scenario_barfill (scenario -> pattern) is supplied
    bar_layer = if(is.null(scenario_barfill)){
      geom_bar(stat='identity', position=position_dodge(width=0.8), width=0.7)
    } else {
      ggpattern::geom_bar_pattern(aes(pattern=scenario), stat='identity', position=position_dodge(width=0.8), width=0.7,
                                  pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                  pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5)
    }
    gg = ggplot(cur_df, aes(x=metric_family, y=mean_rel, fill=scenario)) +
      bar_layer +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_y_use, standard_max_y_use)) +
      scale_fill_manual(values = scenario_palette) +
      (if(!is.null(scenario_barfill)) scale_pattern_manual(values=scenario_barfill)) +
      (if(!is.null(scenario_barfill)) guides(pattern=guide_legend(override.aes=list(fill='white')), fill=guide_legend(override.aes=list(pattern='none')))) +
      geom_hline(yintercept=0, color='black') +
      ylab('Percent reduction') +
      ggtitle(age_titles[cur_age]) +
      theme_gridlines_no_box() +
      theme(legend.position = 'top', legend.box='horizontal', legend.title=element_blank(),
            text=element_text(size=text_size*font_scale),
            legend.text=element_text(size=text_size*legend_scale*font_scale),
            legend.key.size=unit(legend_scale*font_scale, 'lines'),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle=20, vjust=1, hjust=0.8),
            plot.title = element_text(size=rel(1)),
            plot.margin = unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg = gg + geom_errorbar(aes(x=metric_family, ymin=min_rel, ymax=max_rel, group=scenario),
                              position=position_dodge(width=0.8), width=0.3, colour='black', alpha=0.9, size=0.6)
    }
    gg_list[[aa]] = gg
  }

  # extract a single legend, strip from individual panels
  legend_grob = ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))
  for(aa in seq_along(gg_list)){
    gg_list[[aa]] = gg_list[[aa]] + theme(legend.position='none')
  }

  # layout: legend row (spans full width), U5 panel left | all-age panel right
  gg_combined = grid.arrange(grobs = c(list(legend_grob), gg_list),
                             layout_matrix = rbind(c(1, 1), c(2, 3)),
                             heights = c(0.6, 3),
                             widths = c(1, 1))
  return(gg_combined)
}




# ---------------------------------------------------------------------------
# Companion to plot_relative_burden_grouped_barplots(): same data pipeline,
# layout, scenario order and colors, but designed for the case where some
# comparison scenarios have HIGHER burden than the reference (negative percent
# reduction). Adds directional zone shading + labels so negative bars read
# clearly WITHOUT renaming the familiar "Percent reduction" axis:
#   y > 0  (blue zone): scenario has LOWER burden than the modeled
#                       continue-current-implementation projection
#   y < 0  (red zone):  scenario has HIGHER burden (e.g. reduced-investment
#                       scenarios) than that projection
# The reference is itself a modeled scenario (NOT observed / "today's" burden);
# the caption states this explicitly. Blue/red is a colorblind-safe diverging
# pair (RdBu), unlike red/green.
# Drop-in: identical signature to plot_relative_burden_grouped_barplots().
# ---------------------------------------------------------------------------
plot_relative_burden_grouped_barplots_diverging = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                                                           barplot_start_year, barplot_end_year,
                                                           pyr, chw_cov,
                                                           scenario_names, experiment_names, scenario_palette,
                                                           LLIN2y_flag=FALSE, overwrite_files=FALSE, separate_plots_flag=FALSE,
                                                           standard_min_y=0, standard_max_y=0.1, show_error_bar=TRUE, align_seeds=TRUE,
                                                           include_to_present=TRUE, burden_metric_subset=c(),
                                                           zone_pos_color='#2166AC', zone_neg_color='#B2182B', zone_alpha=0.09,
                                                           zone_pos_label='lower malaria burden\nthan current implementation',
                                                           zone_neg_label='higher malaria burden\nthan current implementation',
                                                           zone_label_size=3, zone_label_pad_frac=0.08, zone_label_lineheight=0.8,
                                                           reference_label='continue current implementation',
                                                           font_scale=1, legend_scale=0.6, scenario_barfill=NULL){
  # scenario_barfill (optional named vector scenario -> ggpattern pattern) hashes selected bars
  # (e.g. scaled-up routine-ITN variants); NULL = plain bars (default).
  # font_scale multiplies all in-plot text; legend_scale shrinks the legend relative to body text.
  # Defaults preserve sizing for existing callers except the now-smaller legend. (ggsave dims set
  # at the call site in the main script.)
  admin_pop = read.csv(pop_filepath)

  # burden metrics (same set as plot_relative_burden_grouped_barplots)
  burden_metrics      = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')
  burden_colnames     = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  if(length(burden_metric_subset) >= 1){
    keep_idx = which(burden_metrics %in% burden_metric_subset)
    burden_metrics      = burden_metrics[keep_idx]
    burden_colnames     = burden_colnames[keep_idx]
    burden_metric_names = burden_metric_names[keep_idx]
  }

  # determine reference scenario (same logic as plot_relative_burden_grouped_barplots)
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2
  }

  # compute % reductions for each comparison scenario relative to reference
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir,
                                             reference_experiment_name=reference_experiment_name,
                                             comparison_experiment_name=comparison_experiment_name,
                                             comparison_scenario_name=comparison_scenario_name,
                                             start_year=barplot_start_year, end_year=barplot_end_year,
                                             admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins,
                                             LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    relative_burden_df = relative_burden_df[, which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    } else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }

  # keep scenario factor order matching scenario_names (related scenarios stay grouped)
  scenario_order = scenario_names[comparison_start_index:length(scenario_names)]
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_order)

  # reshape to long form: scenario, burden_col, rel_reduction
  long_df = relative_burden_all_df %>%
    tidyr::pivot_longer(cols = all_of(burden_colnames), names_to = 'burden_col', values_to = 'rel_reduction')

  # attach metric family name (without age suffix) and age group
  col_to_metric_name = setNames(burden_metric_names, burden_colnames)
  long_df$burden_metric_name = col_to_metric_name[long_df$burden_col]
  long_df$age_group     = ifelse(grepl('\\(U5\\)', long_df$burden_metric_name), 'U5', 'all ages')
  long_df$metric_family = gsub(' \\(U5\\)| \\(all ages\\)', '', long_df$burden_metric_name)
  metric_family_levels  = unique(gsub(' \\(U5\\)| \\(all ages\\)', '', burden_metric_names))
  long_df$metric_family = factor(long_df$metric_family, levels=metric_family_levels)

  # aggregate across runs: mean / min / max per scenario x metric x age group
  agg_df = long_df %>%
    dplyr::group_by(scenario, age_group, metric_family) %>%
    dplyr::summarise(mean_rel = mean(rel_reduction),
                     min_rel  = min(rel_reduction),
                     max_rel  = max(rel_reduction),
                     .groups = 'drop')

  # y-range: always include 0; widen to most negative/positive shown (whiskers if drawn)
  y_lo_source = if(show_error_bar) min(agg_df$min_rel) else min(agg_df$mean_rel)
  y_hi_source = if(show_error_bar) max(agg_df$max_rel) else max(agg_df$mean_rel)
  standard_min_y_use = min(standard_min_y, y_lo_source)
  standard_max_y_use = max(standard_max_y, y_hi_source)
  has_neg = standard_min_y_use < 0
  # reserve a gap below the lowest bar so the 'higher burden' label sits in clear space
  # (scaled to the plotted y-range; tune via zone_label_pad_frac)
  label_pad = zone_label_pad_frac * (standard_max_y_use - standard_min_y_use)
  y_floor   = if(has_neg) standard_min_y_use - label_pad else standard_min_y_use

  # one grouped barplot per age group
  age_groups = c('U5', 'all ages')
  age_titles = c('U5' = 'Under-5', 'all ages' = 'All ages')
  gg_list = list()
  for(aa in seq_along(age_groups)){
    cur_age  = age_groups[aa]
    cur_df   = agg_df[agg_df$age_group == cur_age, ]
    n_family = length(levels(cur_df$metric_family))
    x_right  = n_family + 0.45   # near the right edge of the panel, for corner-tucked zone labels
    # bar layer: plain bars, or hashed bars when scenario_barfill (scenario -> pattern) is supplied
    bar_layer = if(is.null(scenario_barfill)){
      geom_bar(stat='identity', position=position_dodge(width=0.8), width=0.7)
    } else {
      ggpattern::geom_bar_pattern(aes(pattern=scenario), stat='identity', position=position_dodge(width=0.8), width=0.7,
                                  pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                  pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5)
    }
    gg = ggplot(cur_df, aes(x=metric_family, y=mean_rel, fill=scenario)) +
      # directional zone shading -- drawn first, sits behind the bars
      annotate('rect', xmin=-Inf, xmax=Inf, ymin=0,    ymax=Inf, fill=zone_pos_color, alpha=zone_alpha) +
      annotate('rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0,   fill=zone_neg_color, alpha=zone_alpha) +
      bar_layer +
      scale_y_continuous(labels=percent_format(), limits=c(y_floor, standard_max_y_use)) +
      scale_fill_manual(values = scenario_palette) +
      (if(!is.null(scenario_barfill)) scale_pattern_manual(values=scenario_barfill)) +
      (if(!is.null(scenario_barfill)) guides(pattern=guide_legend(override.aes=list(fill='white')), fill=guide_legend(override.aes=list(pattern='none')))) +
      geom_hline(yintercept=0, color='black') +
      # zone guide label (top-right corner of blue / lower-burden zone)
      annotate('text', x=x_right, y=standard_max_y_use, label=zone_pos_label,
               hjust=1, vjust=1.1, lineheight=zone_label_lineheight, color=zone_pos_color, fontface='plain', size=zone_label_size, alpha=0.9) +
      ylab('Percent reduction') +
      ggtitle(age_titles[cur_age]) +
      theme_gridlines_no_box() +
      theme(legend.position = 'top', legend.box='horizontal', legend.title=element_blank(),
            text=element_text(size=text_size*font_scale),
            legend.text=element_text(size=text_size*legend_scale*font_scale),
            legend.key.size=unit(legend_scale*font_scale, 'lines'),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle=20, vjust=1, hjust=0.8),
            plot.title = element_text(size=rel(1)),
            plot.margin = unit(c(0,1,1,0), 'cm'))
    if(has_neg){
      # zone guide label (bottom-right corner of red / higher-burden zone) -- only when negatives
      # are present; tucked into the reserved gap below the lowest bar, near the bottom edge
      gg = gg + annotate('text', x=x_right, y=y_floor, label=zone_neg_label,
                         hjust=1, vjust=0, lineheight=zone_label_lineheight, color=zone_neg_color, fontface='plain', size=zone_label_size, alpha=0.9)
    }
    if(show_error_bar){
      gg = gg + geom_errorbar(aes(x=metric_family, ymin=min_rel, ymax=max_rel, group=scenario),
                              position=position_dodge(width=0.8), width=0.3, colour='black', alpha=0.9, size=0.6)
    }
    gg_list[[aa]] = gg
  }

  # extract a single legend, strip from individual panels
  legend_grob = ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))
  for(aa in seq_along(gg_list)){
    gg_list[[aa]] = gg_list[[aa]] + theme(legend.position='none')
  }

  # caption: make clear the reference is a modeled scenario (not observed burden), and state the convention
  caption_grob = grid::textGrob(
    paste0('Bars compare each modeled scenario to the modeled "', reference_label,
           '" projection (a simulated scenario, not observed/current burden). ',
           'Above 0 = lower malaria burden than that scenario; below 0 = higher.'),
    gp = grid::gpar(fontsize = text_size*0.7, fontface='italic'), x=0.01, hjust=0)

  # layout: legend row (spans full width), U5 panel left | all-age panel right, caption underneath
  gg_combined = grid.arrange(grobs = c(list(legend_grob), gg_list),
                             layout_matrix = rbind(c(1, 1), c(2, 3)),
                             heights = c(0.6, 3),
                             widths  = c(1, 1),
                             bottom  = caption_grob)
  return(gg_combined)
}




####################################################################################
# barplots for burden relative to BAU at the state level, displayed in a state grid
####################################################################################

plot_relative_burden_barplots_by_state = function(sim_future_output_dir, pop_filepath,grid_layout_state_locations,
                                                 barplot_start_year, barplot_end_year, 
                                                 pyr, chw_cov,
                                                 scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                                 burden_metric_subset=c(), include_to_present=TRUE, file_suffix='', scenario_barfill=NULL){
  # scenario_barfill (optional named vector scenario -> ggpattern pattern, e.g. 'stripe'/'none') hatches
  # selected bars; NULL = plain solid bars (default, so existing callers are unaffected).
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')#, 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')#, 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }
  
  # first comparison name is to-present (skip it), second is BAU (use as reference), comparison scenarios start at the third index
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2  
  }
  # iterate through comparison scenarios, calculating the burden reduction of all metrics relative to BAU (seedwise comparisons, so one output for each run). Combine all scenario reductions into a dataframe (each scenario set in separate rows)
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden_by_state(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, 
                                                   start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    # only save relevant columns for plotting
    relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', 'State', burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    }else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_names[comparison_start_index:length(scenario_names)])

  # guard against non-finite % reductions (near-zero reference -> (ref-comp)/ref = +/-Inf), which would
  # otherwise crash axis-break computation (Error in seq.default). Report affected scenario x metric x
  # state(s) (so the offending experiment / cumulativeBurden cache can be checked/deleted), then set to NA.
  for(cc in burden_colnames){
    bad = !is.finite(relative_burden_all_df[[cc]])
    if(any(bad)){
      for(sc in unique(as.character(relative_burden_all_df$scenario[bad]))){
        sel = bad & as.character(relative_burden_all_df$scenario) == sc
        sts = sort(unique(as.character(relative_burden_all_df$State[sel])))
        message(sprintf('plot_relative_burden_barplots_by_state: %d non-finite %% reduction value(s) set to NA  [metric: %s | scenario: %s | state(s): %s]', sum(sel), cc, sc, paste(sts, collapse=', ')))
      }
      relative_burden_all_df[[cc]][bad] = NA
    }
  }

  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name = burden_metric_names[bb]
    select_col_names = c(current_burden_name, 'scenario', 'State')
    # get mean, min, and max among all runs for this burden metric
    rel_burden_agg = as.data.frame(relative_burden_all_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
      dplyr::group_by(scenario, State) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name), na.rm=TRUE),
                       max_rel = max(get(current_burden_name), na.rm=TRUE),
                       min_rel = min(get(current_burden_name), na.rm=TRUE))
    
    rel_burden_agg$code = rel_burden_agg$State
    bar_layer = if(is.null(scenario_barfill)){
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity')
    } else {
      ggpattern::geom_bar_pattern(aes(x=scenario, y=mean_rel, fill=scenario, pattern=scenario), stat='identity',
                                  pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                  pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5)
    }
    gg = ggplot(rel_burden_agg) +
      bar_layer +
      scale_y_continuous(labels=percent_format(), n.breaks=4) + #,limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
      # ylab('Percent reduction in burden \n ((Baseline - Plan) / Baseline) * 100') +
      ylab(paste0('Percent reduction in ', burden_metric_name)) +
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values = scenario_palette) +
      (if(!is.null(scenario_barfill)) scale_pattern_manual(values=scenario_barfill)) +
      (if(!is.null(scenario_barfill)) guides(pattern=guide_legend(override.aes=list(fill='white')), fill=guide_legend(override.aes=list(pattern='none')))) +
      theme_classic()+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size),
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm')) +
      facet_geo(~code, grid = grid_layout_state_locations, label="name") #, scales='free')
    
    if(show_error_bar){
      gg = gg +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour="black", alpha=0.9, size=1) 
    }
    ggsave(paste0(sim_future_output_dir, '/_plots/','barplot_percent_reduction_', burden_metric_name,'_stateGrid',file_suffix,'.png'), gg, dpi=600, width=18*.6, height=12*.6, units='in')  # , width=18, height=12, units='in'
  }
}


####################################################################################
# Grouped state-grid barplots: bar groups = burden metric, fill = scenario, faceted by state
####################################################################################
# Companion to plot_relative_burden_barplots_by_state(). When only a few scenarios are
# compared, this avoids producing 6-8 near-empty per-metric files. Instead saves TWO PNGs:
#   - one for U5 metrics (all burden metrics plotted as grouped bars within each state)
#   - one for all-age metrics (same layout)
# Within each state facet: x-axis is burden metric, fill is scenario (dodged bars).
# X-axis text is hidden in every facet; metric order is communicated via the plot subtitle.
# Signature matches plot_relative_burden_barplots_by_state() for easy swap at call sites.
plot_relative_burden_grouped_barplots_by_state = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations,
                                                          barplot_start_year, barplot_end_year,
                                                          pyr, chw_cov,
                                                          scenario_names, experiment_names, scenario_palette, group_col='State', group_levels=NULL,
                                                          LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                                          burden_metric_subset=c(), include_to_present=TRUE, file_suffix='', scenario_barfill=NULL,
                                                          diverging=FALSE, zone_pos_color='#2166AC', zone_neg_color='#B2182B', zone_alpha=0.09,
                                                          zone_pos_label='lower malaria burden than current implementation',
                                                          zone_neg_label='higher malaria burden than current implementation',
                                                          font_scale=1, legend_scale=0.6, save_width=18*.66, save_height=12*.77){
  # scenario_barfill (optional named vector scenario -> ggpattern pattern, e.g. 'stripe'/'none') hatches
  # selected bars; NULL = plain solid bars (default, so existing callers are unaffected).
  # diverging=TRUE shades the panel background blue above 0 (lower burden than current implementation)
  # and red below 0 (higher burden), and adds a caption stating the convention -- making negative
  # 'burden averted' (i.e. worse than current implementation) immediately visible. Filenames get a
  # '_diverging' tag so they don't overwrite the standard version.
  # font_scale multiplies all in-plot text; legend_scale shrinks the legend; save_width/save_height
  # set the saved PNG dimensions. Defaults preserve prior behaviour (except the smaller legend).
  # group_col: column used to group/facet LGAs (default 'State' -> geofacet; any other value, e.g. 'Funder',
  #   uses facet_wrap with fixed shared scales). PNG names use 'stateGrid' for State, '<group_col>Grid' otherwise.
  grid_tag = if(group_col == 'State') 'stateGrid' else paste0(group_col, 'Grid')
  admin_pop = read.csv(pop_filepath)

  # burden metrics (same set as plot_relative_burden_barplots_by_state)
  burden_metrics      = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')
  burden_colnames     = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  if(length(burden_metric_subset) >= 1){
    keep_idx = which(burden_metrics %in% burden_metric_subset)
    burden_metrics      = burden_metrics[keep_idx]
    burden_colnames     = burden_colnames[keep_idx]
    burden_metric_names = burden_metric_names[keep_idx]
  }

  # determine reference scenario (same logic as plot_relative_burden_barplots_by_state)
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2
  }

  # compute % reductions for each comparison scenario relative to reference
  relative_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    relative_burden_df = get_relative_burden_by_state(sim_output_filepath=sim_future_output_dir,
                                                      reference_experiment_name=reference_experiment_name,
                                                      comparison_experiment_name=comparison_experiment_name,
                                                      comparison_scenario_name=comparison_scenario_name,
                                                      start_year=barplot_start_year, end_year=barplot_end_year,
                                                      admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, group_col=group_col,
                                                      LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    relative_burden_df = relative_burden_df[, which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', group_col, burden_colnames))]
    if(nrow(relative_burden_all_df) == 0){
      relative_burden_all_df = relative_burden_df
    } else{
      relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
    }
  }

  # keep scenario factor order matching scenario_names
  scenario_order = scenario_names[comparison_start_index:length(scenario_names)]
  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_order)

  # reshape to long form: scenario, State, burden_col, rel_reduction
  long_df = relative_burden_all_df %>%
    tidyr::pivot_longer(cols = all_of(burden_colnames), names_to = 'burden_col', values_to = 'rel_reduction')

  # attach metric family name (no age suffix) and age group
  col_to_metric_name = setNames(burden_metric_names, burden_colnames)
  long_df$burden_metric_name = col_to_metric_name[long_df$burden_col]
  long_df$age_group     = ifelse(grepl('\\(U5\\)', long_df$burden_metric_name), 'U5', 'all ages')
  long_df$metric_family = gsub(' \\(U5\\)| \\(all ages\\)', '', long_df$burden_metric_name)
  metric_family_levels  = unique(gsub(' \\(U5\\)| \\(all ages\\)', '', burden_metric_names))
  long_df$metric_family = factor(long_df$metric_family, levels=metric_family_levels)

  # aggregate across runs: mean / min / max per scenario x State x metric x age group
  agg_df = long_df %>%
    dplyr::group_by(scenario, across(all_of(group_col)), age_group, metric_family) %>%
    dplyr::summarise(mean_rel = mean(rel_reduction),
                     min_rel  = min(rel_reduction),
                     max_rel  = max(rel_reduction),
                     .groups = 'drop')
  agg_df$code = agg_df[[group_col]]
  if(!is.null(group_levels)) agg_df$code = factor(agg_df$code, levels=group_levels)  # custom facet order

  # subtitle text communicates which metric maps to which x-position
  subtitle_text = paste0('Within each ', tolower(group_col), ' (left -> right): ', paste(metric_family_levels, collapse = '  |  '))

  # one PNG per age group
  age_groups   = c('U5', 'all ages')
  age_titles   = c('U5' = 'Under-5', 'all ages' = 'All ages')
  age_filetags = c('U5' = 'U5', 'all ages' = 'allAges')
  for(aa in seq_along(age_groups)){
    cur_age = age_groups[aa]
    cur_df  = agg_df[agg_df$age_group == cur_age, ]

    gg = ggplot(cur_df, aes(x=metric_family, y=mean_rel, fill=scenario))
    if(diverging){
      # directional zone shading (behind bars): blue above 0 = lower burden, red below 0 = higher burden
      gg = gg +
        annotate('rect', xmin=-Inf, xmax=Inf, ymin=0,    ymax=Inf, fill=zone_pos_color, alpha=zone_alpha) +
        annotate('rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0,   fill=zone_neg_color, alpha=zone_alpha)
    }
    bar_layer = if(is.null(scenario_barfill)){
      geom_bar(stat='identity', position=position_dodge(width=0.55), width=0.5)
    } else {
      ggpattern::geom_bar_pattern(aes(pattern=scenario), stat='identity', position=position_dodge(width=0.55), width=0.5,
                                  pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                  pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5)
    }
    gg = gg +
      bar_layer +
      scale_y_continuous(labels=percent_format(), n.breaks=4) +
      scale_fill_manual(values = scenario_palette) +
      (if(!is.null(scenario_barfill)) scale_pattern_manual(values=scenario_barfill)) +
      (if(!is.null(scenario_barfill)) guides(pattern=guide_legend(override.aes=list(fill='white')), fill=guide_legend(override.aes=list(pattern='none')))) +
      geom_hline(yintercept=0, color='black') +
      ylab('Percent reduction') +
      ggtitle(paste0('Percent reduction in burden (', age_titles[cur_age], ')'),
              subtitle = subtitle_text) +
      theme_classic() +
      theme(legend.position = 'top', legend.box='horizontal', legend.title=element_blank(),
            text=element_text(size=text_size*font_scale),
            legend.text=element_text(size=text_size*legend_scale*font_scale),
            legend.key.size=unit(legend_scale*font_scale, 'lines'),
            plot.subtitle = element_text(size=text_size*font_scale, hjust=0.5),
            axis.title.x=element_blank(),
            axis.text.x=element_text(angle=20, vjust=1, hjust=0.8),
            plot.margin=unit(c(0,1,1,0), 'cm')) +
      (if(group_col == 'State') facet_geo(~code, grid = grid_layout_state_locations, label='name') else facet_wrap(~code))

    if(diverging){
      # caption states the higher-vs-lower-burden convention (zone labels in every facet would clutter)
      gg = gg +
        labs(caption = paste0('Blue (above 0): ', zone_pos_label, '.    Red (below 0): ', zone_neg_label, '.')) +
        theme(plot.caption = element_text(hjust=0, size=text_size*0.75*font_scale))
    }

    if(show_error_bar){
      gg = gg + geom_errorbar(aes(x=metric_family, ymin=min_rel, ymax=max_rel, group=scenario),
                              position=position_dodge(width=0.55), width=0.25, colour='black', alpha=0.9, size=0.6)
    }

    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_percent_reduction_grouped_',
                  age_filetags[cur_age], '_', grid_tag, if(diverging) '_diverging' else '', file_suffix, '.png'),
           gg, dpi=600, width=save_width, height=save_height, units='in')
  }
  invisible(NULL)
}




####################################################################################
# Across-funder comparison barplots WITHIN a scenario
####################################################################################
# For a single scenario (e.g. Scenario 1), compare % reduction in burden relative to current
# implementation ACROSS groups: National (all LGAs) + each funder. Bars are COLOURED BY GROUP
# (funder); the figure is a facet grid with burden metric across columns and age group down rows
# (U5 top, all-ages bottom). Current vs PF target coverage are dodged bars, with current hatched
# (coverage_barfill). One PNG per scenario.
#   reference_experiment_name : current-implementation experiment (the comparison baseline)
#   coverage_experiments      : named char vector, e.g. c(current='s1_curCov', target='s1_targetCov')
#                               (drop an element to show only one coverage)
#   funder_palette            : named colours for 'National' + each funder value
#   group_levels              : x-axis order of groups (default: National first, then funders)
plot_relative_burden_across_funders = function(sim_future_output_dir, pop_filepath,
        reference_experiment_name, coverage_experiments, scenario_label,
        barplot_start_year, barplot_end_year,
        funder_palette, funder_col='Funder', group_levels=NULL,
        burden_metric_subset=c('PfPR','incidence','allMortality'),
        overwrite_files=FALSE, align_seeds=TRUE,
        coverage_barfill=c(current='stripe', target='none'),
        bar_width=0.6, national_width_factor=2,   # National bars drawn national_width_factor x wider (it spans all funders)
        diverging=TRUE, zone_pos_color='#2166AC', zone_neg_color='#B2182B', zone_alpha=0.09,
        zone_pos_label='lower malaria burden than current implementation',
        zone_neg_label='higher malaria burden than current implementation',
        file_suffix='', font_scale=1, legend_scale=0.6, save_width=10, save_height=7){

  admin_pop = read.csv(pop_filepath)
  admin_pop$National_group = 'National'   # constant column -> National aggregated like a funder group

  # burden metric column -> family/age mapping (subset to requested metrics)
  metric_map = data.frame(
    col    = c('average_PfPR_U5','average_PfPR_all','incidence_U5','incidence_all','direct_death_rate_mean_U5','direct_death_rate_mean_all','all_death_rate_mean_U5','all_death_rate_mean_all'),
    metric = c('PfPR','PfPR','incidence','incidence','directMortality','directMortality','allMortality','allMortality'),
    family = c('PfPR','PfPR','incidence','incidence','mortality','mortality','mortality','mortality'),
    age    = c('U5','all ages','U5','all ages','U5','all ages','U5','all ages'),
    stringsAsFactors=FALSE)
  metric_map = metric_map[metric_map$metric %in% burden_metric_subset, ]
  metric_cols = metric_map$col

  # assemble mean % reduction per group x coverage x metric (National + each funder, pop-weighted)
  rows = data.frame()
  for(cov in names(coverage_experiments)){
    comp = coverage_experiments[[cov]]
    for(gcol in c(funder_col, 'National_group')){
      rb = get_relative_burden_by_state(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name,
              comparison_experiment_name=comp, comparison_scenario_name=scenario_label,
              start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop,
              group_col=gcol, overwrite_files=overwrite_files, align_seeds=align_seeds)
      m = aggregate(rb[, metric_cols, drop=FALSE], by=list(group=rb[[gcol]]), FUN=function(x) mean(x, na.rm=TRUE))
      m$coverage = cov
      rows = rbind(rows, m)
    }
  }

  # long form + metric family/age labels
  long = tidyr::pivot_longer(rows, cols=all_of(metric_cols), names_to='col', values_to='rel_reduction')
  long = merge(long, metric_map[, c('col','family','age')], by='col')
  long$family   = factor(long$family, levels=intersect(c('PfPR','incidence','mortality'), unique(metric_map$family)))
  long$age      = factor(long$age, levels=c('U5','all ages'), labels=c('Under-5','All ages'))   # U5 top row, all-ages bottom
  if(is.null(group_levels)) group_levels = c('National', setdiff(unique(rows$group), 'National'))
  long$group    = factor(long$group, levels=group_levels)
  long$coverage = factor(long$coverage, levels=intersect(c('current','target'), names(coverage_experiments)))

  # numeric x-positions: each group center placed so edge-to-edge gap is equal for all adjacent pairs
  nat_width  = bar_width * national_width_factor
  grp_widths = ifelse(group_levels == 'National', nat_width, bar_width)
  names(grp_widths) = group_levels
  gap        = 1 - bar_width     # matches default funder-funder spacing from scale_x_discrete
  x_centers  = numeric(length(group_levels))
  names(x_centers) = group_levels
  for(i in seq_along(group_levels)[-1])
    x_centers[i] = x_centers[i-1] + grp_widths[i-1]/2 + gap + grp_widths[i]/2
  long$x_pos = x_centers[as.character(long$group)]

  has_neg = any(long$rel_reduction < 0, na.rm=TRUE)
  gg = ggplot(long, aes(x=x_pos, y=rel_reduction, fill=group, pattern=coverage, group=coverage))
  if(diverging && has_neg){
    # no blue zone above 0; shade red below 0 ONLY when some bar is negative (worse than current implementation)
    gg = gg +
      annotate('rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, fill=zone_neg_color, alpha=zone_alpha)
  }
  gg = gg +
    # funder bars (standard width); National bars (wider) drawn as a second layer
    ggpattern::geom_bar_pattern(data = function(d) d[d$group != 'National', ],
                                stat='identity', position=position_dodge(width=bar_width), width=bar_width,
                                pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5) +
    ggpattern::geom_bar_pattern(data = function(d) d[d$group == 'National', ],
                                stat='identity', position=position_dodge(width=nat_width), width=nat_width, show.legend=FALSE,
                                pattern_fill='white', pattern_colour=NA, pattern_angle=45,
                                pattern_density=0.3, pattern_spacing=0.03, pattern_key_scale_factor=0.5) +
    scale_fill_manual(values=funder_palette) +
    scale_pattern_manual(values=coverage_barfill,
                         labels=c(current='current coverage', target='PF target coverage')) +
    scale_x_continuous(breaks=x_centers, labels=names(x_centers)) +
    guides(fill='none', pattern=guide_legend(title=NULL, override.aes=list(fill='grey70'))) +   # x-axis already labels the groups
    scale_y_continuous(labels=percent_format(), n.breaks=4) +
    geom_hline(yintercept=0, color='black') +
    facet_grid(age ~ family) +                                                                       # rows = age (U5 top), cols = metric
    ylab('Percent reduction') +
    ggtitle(paste0('Burden averted by ', scenario_label, ' relative to continued current implementation')) +
    theme_classic() +
    theme(legend.position='top', legend.justification='right', legend.box='vertical',
          text=element_text(size=text_size*font_scale),
          legend.text=element_text(size=text_size*font_scale),
          legend.key.size=unit(legend_scale*font_scale,'lines'),
          axis.title.x=element_blank(), axis.text.x=element_text(angle=30, vjust=1, hjust=1),
          plot.margin=unit(c(0,1,1,0),'cm'))
  if(diverging && has_neg){
    # caption only when the red zone is shown
    gg = gg + labs(caption=paste0('Red (below 0): ', zone_neg_label, '.')) +
      theme(plot.caption=element_text(hjust=0, size=text_size*0.75*font_scale))
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/barplot_percent_reduction_acrossFunders_', gsub('[^A-Za-z0-9]+','',scenario_label), file_suffix, '.png'),
         gg, dpi=600, width=save_width, height=save_height, units='in')
  invisible(NULL)
}




####################################################################################
# barplots for burden reduction relative to BAU: difference
####################################################################################

plot_difference_burden_barplots = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                         barplot_start_year, barplot_end_year, 
                                         pyr, chw_cov,
                                         scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, separate_plots_flag=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                         include_to_present=TRUE, burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)
  
  # burden metrics
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')#, 'mLBW_deaths', 'MiP_stillbirths')
  burden_colnames = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')#, 'annual_num_mLBW', 'annual_num_mStill')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
  }
  
  # first comparison name is to-present (skip it), second is BAU (use as reference), comparison scenarios start at the third index
  if(include_to_present){
    reference_experiment_name = experiment_names[2]
    comparison_start_index = 3
  } else{
    reference_experiment_name = experiment_names[1]
    comparison_start_index = 2  
  }
  # iterate through comparison scenarios, calculating the burden reduction of all metrics relative to BAU (seedwise comparisons, so one output for each run). Combine all scenario reductions into a dataframe (each scenario set in separate rows)
  difference_burden_all_df = data.frame()
  for(ss in comparison_start_index:length(scenario_names)){
    comparison_experiment_name = experiment_names[ss]
    comparison_scenario_name = scenario_names[ss]
    difference_burden_df = get_difference_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, 
                                             start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
    # only save relevant columns for plotting
    difference_burden_df = difference_burden_df[,which(colnames(difference_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
    if(nrow(difference_burden_all_df) == 0){
      difference_burden_all_df = difference_burden_df
    }else{
      difference_burden_all_df = rbind(difference_burden_all_df, difference_burden_df)
    }
  }
  
  # get factors in the correct order (rather than alphabetical)
  difference_burden_all_df$scenario = factor(difference_burden_all_df$scenario, levels=scenario_names[comparison_start_index:length(scenario_names)])
  
  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_y = 0
  standard_max_y = 0.1
  cur_min = min(difference_burden_all_df[,2:(1+length(burden_colnames))])
  cur_max = max(difference_burden_all_df[,2:(1+length(burden_colnames))])
  if(cur_min < standard_min_y) standard_min_y = cur_min
  if(cur_max > standard_max_y) standard_max_y = cur_max
  
  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name = burden_metric_names[bb]
    select_col_names = c(current_burden_name, 'scenario')
    # get mean, min, and max among all runs for this burden metric
    rel_burden_agg = as.data.frame(difference_burden_all_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                       max_rel = max(get(current_burden_name)),
                       min_rel = min(get(current_burden_name)))
    
    gg_list[[bb]] = ggplot(rel_burden_agg) + 
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity') +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
      ylab('Burden averted') + 
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values = scenario_palette) + 
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size), 
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg_list[[bb]] = gg_list[[bb]] +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour="black", alpha=0.9, size=1) 
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] + 
        ylab('Burden averted \n ((Baseline - Plan)') + 
        theme(legend.position='none', plot.title = element_blank(), text=element_text(size =separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/','barplot_burden_averted_', burden_metric_name,'_',district_subset,'.png'), separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }
  # for each burden type, 
  # get mean, min, and max among all runs for each burden metric, each saved as a separate column
  # create barplot for each burden type (using columns of dataframe, separate bar for each scenario)
  
  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  # remove legend from main plots
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position = "none")  + theme(text = element_text(size = text_size))   
  }
  
  if(save_plots){
    gg_saved = grid.arrange(grobs = gg_list[-1], layout_matrix = matrix(c(1:(length(burden_colnames))), nrow=2, byrow=FALSE))
    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_burden_averted_', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg_saved, dpi=600, width=14, height=7, units='in')
  }
  
  # ----- combine all burden plots ----- #
  # gg = grid.arrange(grobs = gg_list, layout_matrix = matrix(c(1,1,2:(length(burden_colnames)+1)), ncol=2, byrow=TRUE))
  gg = grid.arrange(grobs = gg_list, layout_matrix = rbind(matrix(rep(1, length(burden_colnames)/2), nrow=1), matrix(2:(length(burden_colnames)+1), nrow=2, byrow=FALSE)))
  
  return(gg)
}






####################################################################################################################################
# barplot of the impact a specific intervention has in relevant admins 
#  (percent reduction when intervention is included versus matched simulation without the intervention)
####################################################################################################################################


plot_barplot_impact_specific_intervention = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                              barplot_start_year, barplot_end_year, 
                                              pyr, chw_cov,
                                              experiment_names_without, experiment_names_with, scenario_palette, scenario_barfill=NA, intervention_name='PMC', age_group = 'U1', LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                              burden_metric_subset=c(), default_ylim_max=0.03,
                                              font_scale=1, legend_scale=0.6, save_width=4.8, save_height=4.8){
  # font_scale multiplies all in-plot text; legend_scale shrinks the legend; save_width/save_height
  # apply to the optional single-plot save. (Combo panels are assembled + saved by the caller.)
  admin_pop = read.csv(pop_filepath)
  comparison_scenario_name = intervention_name
  
  # iterate through the matched pairs of experiments without / with the intervention
  rel_burden_agg = data.frame()
  if(length(experiment_names_without) == length(experiment_names_with)){
    for(ii in 1:length(experiment_names_without)){
      # first experiment is without interventions, second experiment is with intervention
      reference_experiment_name = experiment_names_without[ii]
      # calculating the burden reduction of all metrics relative to no intervention (seedwise comparisons, so one output for each run). 
      comparison_experiment_name = experiment_names_with[ii]  
      
      # set which burden metrics are relevant and get relative burden between simulations
      burden_metrics = c('PfPR', 'incidence', 'directMortality', 'allMortality')
      if(age_group=='U1'){
        burden_colnames = c('average_PfPR_U1', 'incidence_U1', 'direct_death_rate_mean_U1', 'all_death_rate_mean_U1')
        burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'direct mortality (U1)', 'mortality (U1)')
        relative_burden_df = get_relative_U1_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      } else if (age_group=='U5'){
        burden_colnames = c('average_PfPR_U5', 'incidence_U5', 'direct_death_rate_mean_U5', 'all_death_rate_mean_U5')
        burden_metric_names = c('PfPR (U5)', 'incidence (U5)', 'direct mortality (U5)', 'mortality (U5)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }else{
        burden_metrics = c(burden_metrics)#, 'mLBW_deaths', 'MiP_stillbirths')
        burden_colnames = c('average_PfPR_all', 'incidence_all', 'direct_death_rate_mean_all', 'all_death_rate_mean_all')#, 'annual_num_mLBW', 'annual_num_mStill')
        burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'direct mortality (all ages)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }
      
  
      # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
      if((length(burden_metric_subset)>=1)){
        burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
        burden_colnames = burden_colnames[burden_metrics_subset_indices]
        burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
      }
      
      # only save relevant columns for plotting
      relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
      
      for(bb in 1:length(burden_colnames)){
        current_burden_name = burden_colnames[bb]
        burden_metric_name = burden_metric_names[bb]
        select_col_names = c(current_burden_name, 'scenario')
        # get mean, min, and max among all runs for this burden metric
        rel_burden_agg_bb = as.data.frame(relative_burden_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
          dplyr::group_by(scenario) %>%
          dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                           max_rel = max(get(current_burden_name)),
                           min_rel = min(get(current_burden_name)))
        rel_burden_agg_bb$burden_metric = burden_metric_name
        rel_burden_agg_bb$scenario_name = experiment_names_with[ii]
        if(nrow(rel_burden_agg)<1){
          rel_burden_agg = rel_burden_agg_bb
        } else{
          rel_burden_agg = merge(rel_burden_agg, rel_burden_agg_bb, all=TRUE)
        }
      }
    }
  }

  rel_burden_agg$burden_metric = gsub('\\(births\\)', '', rel_burden_agg$burden_metric)
  rel_burden_agg$burden_metric = factor(rel_burden_agg$burden_metric, levels=gsub('\\(births\\)', '', burden_metric_names))
  rel_burden_agg$scenario_name = factor(rel_burden_agg$scenario_name, levels=experiment_names_with)

  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_y = 0
  standard_max_y = default_ylim_max
  cur_min = min(rel_burden_agg[,2:4])
  cur_max = max(rel_burden_agg[,2:4])
  if(cur_min < standard_min_y) standard_min_y = cur_min
  if(cur_max > standard_max_y) standard_max_y = cur_max
  if(any(is.na(scenario_barfill))){
    scenario_barfill = rep('none', length(unique(rel_burden_agg$scenario_name)))
    names(scenario_barfill) = unique(rel_burden_agg$scenario_name)
  }
  # original without shading:
  # gg = ggplot(rel_burden_agg) +
  #   geom_bar(aes(x=burden_metric, y=mean_rel, fill=scenario_name), stat='identity', position="dodge") +
  #   scale_y_continuous(labels=percent_format(), limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
  #   ylab(paste0('Percent reduction in burden \n ((without ', intervention_name, ' - with ', intervention_name, ') / without ', intervention_name, ') * 100')) +
  #   geom_hline(yintercept=0, color='black') +
  #   ggtitle(paste0('Comparison of burden in proposed ', intervention_name, ' districts')) +
  #   scale_fill_manual(values = scenario_palette) +
  #   theme_classic()+
  #   theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size),
  #         axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
  #         plot.margin=unit(c(0,1,1,0), 'cm'))
  
  gg = ggplot(rel_burden_agg) + 
    # geom_bar(aes(x=burden_metric, y=mean_rel, fill=scenario_name, pattern=scenario_name), stat='identity', position="dodge") +
    scale_y_continuous(labels=percent_format(), limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
    ylab(paste0('Percent reduction in burden \n ((without ', intervention_name, ' - with ', intervention_name, ') / without ', intervention_name, ') * 100')) + 
    geom_hline(yintercept=0, color='black') +
    ggtitle(paste0('Comparison of burden in proposed ', intervention_name, ' districts')) + 
    scale_fill_manual(values = scenario_palette) + 
    geom_bar_pattern(aes(x=burden_metric, y=mean_rel, fill=scenario_name, pattern=scenario_name), stat='identity', position="dodge", #position = position_dodge(preserve = "single"),
                     # color = "white",
                     pattern_fill = "white",
                     pattern_linetype=0,
                     pattern_angle = 45,
                     pattern_density = 0.35,
                     pattern_spacing = 0.06,# 0.025,
                     pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = scenario_barfill) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    # theme_classic()+ 
    theme_gridlines_no_box()+ 
    theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size*font_scale), legend.text=element_text(size = text_size*legend_scale*font_scale), legend.key.size=unit(legend_scale*font_scale, 'lines'),
          axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
          plot.margin=unit(c(0,1,1,0), 'cm'))
  
  if(show_error_bar){
    gg = gg +
      geom_errorbar(aes(x=burden_metric, ymin=min_rel, ymax=max_rel, group=scenario_name), position='dodge',  colour="black", alpha=0.9, size=1) # width=0.4,
  }

  if(save_plots){
    ggsave(paste0(sim_future_output_dir, '/_plots/barplot_', intervention_name, '_percent_reduction_burden_', age_group, '', barplot_start_year, '_', barplot_end_year, '.png'), gg, dpi=600, width=save_width, height=save_height, units='in')
  }
  
  return(gg)
}






####################################################################################################################################
# barplot of the impact a specific intervention has in relevant admins 
#  (burden reduction when intervention is included versus matched simulation without the intervention)
####################################################################################################################################


table_difference_impact_specific_intervention = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                                     barplot_start_year, barplot_end_year, 
                                                     pyr, chw_cov,
                                                     experiment_names_without, experiment_names_with, intervention_name='PMC', age_group = 'U1', LLIN2y_flag=FALSE, overwrite_files=FALSE, align_seeds=TRUE,
                                                     burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)
  comparison_scenario_name = intervention_name
  
  # iterate through the matched pairs of experiments without / with the intervention
  rel_burden_agg = data.frame()
  if(length(experiment_names_without) == length(experiment_names_with)){
    for(ii in 1:length(experiment_names_without)){
      # first experiment is without interventions, second experiment is with intervention
      reference_experiment_name = experiment_names_without[ii]
      # calculating the burden reduction of all metrics relative to no intervention (seedwise comparisons, so one output for each run). 
      comparison_experiment_name = experiment_names_with[ii]  
      
      # set which burden metrics are relevant and get relative burden between simulations
      burden_metrics = c('PfPR', 'incidence', 'directMortality', 'allMortality')
      if(age_group=='U1'){
        warning('Have not yet added support for burden reduction for U1... need to add relevant functions')
        # burden_colnames = c('average_PfPR_U1', 'incidence_U1', 'direct_death_rate_mean_U1', 'all_death_rate_mean_U1')
        # burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'direct mortality (U1)', 'mortality (U1)')
        # relative_burden_df = get_relative_U1_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      } else if (age_group=='U5'){
        burden_colnames = c('average_PfPR_U5', 'incidence_U5', 'direct_death_rate_mean_U5', 'all_death_rate_mean_U5')
        burden_metric_names = c('PfPR (U5)', 'incidence (U5)', 'direct mortality (U5)', 'mortality (U5)')
        relative_burden_df = get_difference_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }else{
        burden_metrics = c(burden_metrics)#, 'mLBW_deaths', 'MiP_stillbirths')
        burden_colnames = c('average_PfPR_all', 'incidence_all', 'direct_death_rate_mean_all', 'all_death_rate_mean_all')#, 'annual_num_mLBW', 'annual_num_mStill')
        burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'direct mortality (all ages)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
        relative_burden_df = get_difference_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }
      
      
      # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
      if((length(burden_metric_subset)>=1)){
        burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
        burden_colnames = burden_colnames[burden_metrics_subset_indices]
        burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
      }
      
      # only save relevant columns for plotting
      relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]
      
      for(bb in 1:length(burden_colnames)){
        current_burden_name = burden_colnames[bb]
        burden_metric_name = burden_metric_names[bb]
        select_col_names = c(current_burden_name, 'scenario')
        # get mean, min, and max among all runs for this burden metric
        rel_burden_agg_bb = as.data.frame(relative_burden_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
          dplyr::group_by(scenario) %>%
          dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                           max_rel = max(get(current_burden_name)),
                           min_rel = min(get(current_burden_name)))
        rel_burden_agg_bb$burden_metric = burden_metric_name
        rel_burden_agg_bb$scenario_name = experiment_names_with[ii]
        if(nrow(rel_burden_agg)<1){
          rel_burden_agg = rel_burden_agg_bb
        } else{
          rel_burden_agg = merge(rel_burden_agg, rel_burden_agg_bb, all=TRUE)
        }
      }
    }
  }
  
  rel_burden_agg$burden_metric = gsub('\\(births\\)', '', rel_burden_agg$burden_metric)
  rel_burden_agg$burden_metric = factor(rel_burden_agg$burden_metric, levels=gsub('\\(births\\)', '', burden_metric_names))
  rel_burden_agg$scenario_name = factor(rel_burden_agg$scenario_name, levels=experiment_names_with)
  
  return(rel_burden_agg)
}







plot_barplot_impact_two_specific_interventions = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                                                     barplot_start_year, barplot_end_year,
                                                     pyr, chw_cov,
                                                     experiment_names_without, experiment_names_with, scenario_palette, intervention_name='PMC', age_group = 'U1', LLIN2y_flag=FALSE, overwrite_files=FALSE, show_error_bar=TRUE, align_seeds=TRUE,
                                                     intervention_strings = c('Vacc')){
  admin_pop = read.csv(pop_filepath)
  comparison_scenario_name = intervention_name

  # iterate through the matched pairs of experiments without / with the intervention
  rel_burden_agg = data.frame()
  if(length(experiment_names_without) == length(experiment_names_with)){
    for(ii in 1:length(experiment_names_without)){
      # first experiment is without interventions, second experiment is with intervention
      reference_experiment_name = experiment_names_without[ii]
      # calculating the burden reduction of all metrics relative to no intervention (seedwise comparisons, so one output for each run).
      comparison_experiment_name = experiment_names_with[ii]

      # set which burden metrics are relevant and get relative burden between simulations
      if(age_group=='U1'){
        burden_colnames = c('average_PfPR_U1', 'incidence_U1', 'direct_death_rate_mean_U1', 'all_death_rate_mean_U1')
        burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'direct mortality (U1)', 'mortality (U1)')
        relative_burden_df = get_relative_U1_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      } else if (age_group=='U5'){
        burden_colnames = c('average_PfPR_U5', 'incidence_U5', 'direct_death_rate_mean_U5', 'all_death_rate_mean_U5')
        burden_metric_names = c('PfPR (U5)', 'incidence (U5)', 'direct mortality (U5)', 'mortality (U5)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }else{
        burden_colnames = c('average_PfPR_all', 'incidence_all', 'direct_death_rate_mean_all', 'all_death_rate_mean_all')
        burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'direct mortality (all ages)', 'mortality (all ages)')
        relative_burden_df = get_relative_burden(sim_output_filepath=sim_future_output_dir, reference_experiment_name=reference_experiment_name, comparison_experiment_name=comparison_experiment_name, comparison_scenario_name=comparison_scenario_name, start_year=barplot_start_year, end_year=barplot_end_year, admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins, LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files, align_seeds=align_seeds)
      }

      # only save relevant columns for plotting
      relative_burden_df = relative_burden_df[,which(colnames(relative_burden_df) %in% c('scenario', 'Run_Number', burden_colnames))]

      for(bb in 1:length(burden_colnames)){
        current_burden_name = burden_colnames[bb]
        burden_metric_name = burden_metric_names[bb]
        select_col_names = c(current_burden_name, 'scenario')
        # get mean, min, and max among all runs for this burden metric
        rel_burden_agg_bb = as.data.frame(relative_burden_df) %>% dplyr::select(match(select_col_names, names(.))) %>%
          dplyr::group_by(scenario) %>%
          dplyr::summarise(mean_rel = mean(get(current_burden_name)),
                           max_rel = max(get(current_burden_name)),
                           min_rel = min(get(current_burden_name)))
        rel_burden_agg_bb$burden_metric = burden_metric_name
        rel_burden_agg_bb$scenario_name = experiment_names_with[ii]
        if(length(intervention_strings)>1){
          rel_burden_agg_bb$intervention_info = NA
          for(jj in 1:length(intervention_strings)){
            if(grepl(intervention_strings[jj], experiment_names_without[ii])){
              rel_burden_agg_bb$intervention_info = intervention_strings[jj]
            }
          }
        }
        if(nrow(rel_burden_agg)<1){
          rel_burden_agg = rel_burden_agg_bb
        } else{
          rel_burden_agg = merge(rel_burden_agg, rel_burden_agg_bb, all=TRUE)
        }
      }
    }
  }

  rel_burden_agg$burden_metric = factor(rel_burden_agg$burden_metric, levels=burden_metric_names)

  # get minimum and maximum reductions - these will be used if they are smaller / greater than the current min/max
  standard_min_y = 0
  standard_max_y = 0.2
  cur_min = min(rel_burden_agg[,2:4])
  cur_max = max(rel_burden_agg[,2:4])
  if(cur_min < standard_min_y) standard_min_y = cur_min
  if(cur_max > standard_max_y) standard_max_y = cur_max

  
  # create list where each element is a barplot corresponding to one of the interventions in intervention_strings
  gg_list = list()
  for(jj in 1:length(intervention_strings)){
    gg = ggplot(rel_burden_agg[rel_burden_agg$intervention_info == intervention_strings[jj],]) +
      geom_bar(aes(x=burden_metric, y=mean_rel, fill=scenario_name), stat='identity', position="dodge") +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_y, standard_max_y)) +   # turn into percent reduction
      ylab(paste0('Percent reduction in burden \n ((without ', intervention_strings[jj], ' - with ', intervention_name, ') / without ', intervention_strings[jj], ') * 100')) +
      geom_hline(yintercept=0, color='black') +
      ggtitle(paste0('Comparison of burden in proposed ', intervention_name, ' districts')) +
      scale_fill_manual(values = scenario_palette) +
      theme_classic()+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size),
            axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    
    if(show_error_bar){
      gg = gg +
        geom_errorbar(aes(x=burden_metric, ymin=min_rel, ymax=max_rel, group=scenario_name), position='dodge',  colour="black", alpha=0.9, size=1) # width=0.4,
    }
    gg_list[[jj]] = gg
  }
  return(gg_list)
}







######################################################################
# create plot panel with all burden metrics, no intervention info (either showing burden or burden relative to burden in specified year)
######################################################################

plot_simulation_output_burden_all = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                             plot_by_month, min_year, max_year, sim_end_years, relative_year=NA,
                                             pyr='', chw_cov='',
                                             scenario_filepaths, scenario_names, experiment_names, scenario_palette, LLIN2y_flag=FALSE, overwrite_files=FALSE, 
                                             separate_plots_flag=FALSE, extend_past_timeseries_year=NA, scenario_linetypes=NA, plot_CI=TRUE, include_U1=FALSE,
                                             burden_metric_subset=c(), ymax_each_burden=NA,
                                             font_scale=1, legend_scale=0.6, display_xlim=NULL){
  # font_scale multiplies all in-plot text; legend_scale shrinks the legend. Defaults preserve
  # prior sizing except the now-smaller legend. (ggsave dims are set at the call site.)
  # display_xlim: optional c(lo, hi) to CROP the x-axis view (annual plots only) without dropping data
  #   -- e.g. display_xlim=c(2026, 2028.2) zooms to mid-2028 while keeping the 2028->2029 line segment.
  #   Keep max_year at the data's last year (e.g. 2029) so that segment exists to be cropped.
  if (!is.na(relative_year)){ if(relative_year<min_year){
    warning('specified minimum year must be <= relative year. Setting min_year to relative_year.')
    min_year = relative_year
  }}
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  pop_sizes = read.csv(pop_filepath)
  pop_sizes = pop_sizes[,c('admin_name','pop_size')]
  # if we include all admins, get list of names from population size dataframe
  if(cur_admins[1] == 'all'){
    cur_admins = unique(pop_sizes$admin_name)
  }
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  
  # ----- malaria burden ----- #
  burden_metrics = c('PfPR', 'PfPR', 'incidence', 'incidence', 'directMortality', 'directMortality', 'allMortality', 'allMortality')#, 'mLBW_deaths', 'MiP_stillbirths')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
  burden_colnames = c('PfPR_U5', 'PfPR_MiP_adjusted', 'New_clinical_cases_U5', 'New_Clinical_Cases', 'direct_mortality_nonMiP_U5_mean', 'direct_mortality_nonMiP_mean', 'total_mortality_U5_mean', 'total_mortality_mean')#, 'mLBW_deaths', 'MiP_stillbirths')    
  if(include_U1){
    burden_metrics = c('PfPR', 'PfPR', 'PfPR', 'incidence','incidence', 'incidence', 'directMortality', 'directMortality', 'directMortality', 'allMortality', 'allMortality', 'allMortality')#, 'mLBW_deaths', 'MiP_stillbirths')
    burden_metric_names = c('PfPR (U1)', 'PfPR (U5)', 'PfPR (all ages)', 'incidence (U1)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U1)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U1)', 'mortality (U5)', 'mortality (all ages)')#, 'mLBW mortality (births)', 'stillbirths (births)')
    burden_colnames = c('PfPR_U1', 'PfPR_U5', 'PfPR_MiP_adjusted', 'New_clinical_cases_U1', 'New_clinical_cases_U5', 'New_Clinical_Cases', 'direct_mortality_nonMiP_U1_mean', 'direct_mortality_nonMiP_U5_mean', 'direct_mortality_nonMiP_mean', 'total_mortality_U1_mean', 'total_mortality_U5_mean', 'total_mortality_mean')#, 'mLBW_deaths', 'MiP_stillbirths')
  }
  # allow subsetting of which burden metrics plotted (based on burden_metric_subset argument)
  if((length(burden_metric_subset)>=1)){
    burden_metrics_subset_indices = which(burden_metrics %in% burden_metric_subset)
    burden_colnames = burden_colnames[burden_metrics_subset_indices]
    burden_metric_names = burden_metric_names[burden_metrics_subset_indices]
    burden_metrics = burden_metrics[burden_metrics_subset_indices]
  }
    
  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    burden_metric_name = burden_metric_names[bb]
    burden_colname = burden_colnames[bb]
    burden_metric = burden_metrics[bb]
    
    if(grepl('U1', burden_metric_name)){
      age_plotted = 'U1'
    } else if(grepl('U5', burden_metric_name)){
      age_plotted = 'U5'
    } else if(grepl('births', burden_metric_name)){
      age_plotted = 'births'
    } else age_plotted = 'all'
    
    
    # burden timeseries are built in memory below from the per-experiment caches (no combined df saved)
    if(LLIN2y_flag){
      llin2y_string = '_2yLLIN'
    } else{
      llin2y_string = ''
    }
    
    # iterate through scenarios, storing relevant output
    burden_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_sim_output_agg = get_burden_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                     exp_name = scenario_names[ee], district_subset=district_subset,
                                                     cur_admins=cur_admins, pop_sizes = pop_sizes, min_year=min_year, max_year=max_year, burden_colname=burden_colname, age_plotted=age_plotted, plot_by_month=plot_by_month)
      if(nrow(burden_df)==0){
        burden_df = cur_sim_output_agg
      } else{
        burden_df = rbind(burden_df, cur_sim_output_agg)
      }
    }
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # Connect the 'to-present' and 'future-projection' simulations in the plot
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # Two alternatives for how this is done, controlled by extend_past_timeseries:
    #   - (FALSE) extend the 'future-projection' lines all back to the end of the 'to-present' simulations, which is desirable if the future projection scenarios separate right away
    #   - (TRUE) extend the end of the 'to-present' line up to the specified point in the 'future-projections' timeseries. This is only desirable if all 'future-projections' are 
    #            identical up to that point (e.g., 'to-present' simulations only run to 2020 and we are currently in 2023, so 2021-2022 are identical in all 'future projection' scenarios)
    if('to-present' %in% burden_df$scenario){
      connect_future_with_past = TRUE
      similarity_threshold = 0.15
      if(!is.na(extend_past_timeseries_year) & (extend_past_timeseries_year %in% burden_df$year[burden_df$scenario != 'to-present']) & time_string=='annual'){
        # check whether the future projections are all nearly identical (minus stochasticity) for the initial year (otherwise, use version that extends future-projection lines back to past)
        future_df = burden_df[burden_df$scenario != 'to-present',]
        earliest_future_year = min(future_df$year)
        compare_burdens = future_df$mean_burden[future_df$year == earliest_future_year]
        if(all(compare_burdens<(compare_burdens[1]*(1+similarity_threshold))) & all(compare_burdens>(compare_burdens[1]*(1-similarity_threshold)))){
          connect_future_with_past = FALSE
          merge_years = earliest_future_year
          if(extend_past_timeseries_year > earliest_future_year){
            # check which years (up to a maximum of extend_past_timeseries_year) should be included in the to-present line
            yy = earliest_future_year + 1
            while(yy <= extend_past_timeseries_year){
              compare_burdens = future_df$mean_burden[future_df$year == (yy)]
              if(all(compare_burdens<(compare_burdens[1]*1.05)) & all(compare_burdens>(compare_burdens[1]*0.95))){
                merge_years = c(merge_years, yy)
                yy = yy+1
              } else{  # as soon as they don't match for a year, stop trying to match any future years
                yy=99999999
              }
            }
          }
          # get the mean value from the 'future-projection' rows so that it can be added to the 'to-present' scenario
          past_from_future_df = future_df[future_df$year %in% merge_years,]
          past_from_future_df_means = past_from_future_df %>% dplyr::select(-scenario) %>% group_by(year) %>%
            summarise_all(mean) %>% ungroup()
          past_from_future_df_means$scenario = 'to-present'
          # delete the old 'future-projection' rows for all but the final of these years
          delete_future_years = merge_years[merge_years != max(merge_years)]
          if(length(delete_future_years)>0) burden_df = burden_df[-which(burden_df$year %in% merge_years),]
          # add the rows to the 'to-present' scenario in the data frame
          burden_df = merge(burden_df, past_from_future_df_means, all=TRUE)
        }else{
          connect_future_with_past = TRUE
        }
      } 
      if(connect_future_with_past){
        # add the final 'to-present' row to all future simulations for a continuous plot
        to_present_df = burden_df[burden_df$scenario == 'to-present',]
        if(plot_by_month){
          final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
          for(ss in 2:length(scenario_names)){
            final_to_present_row$scenario = scenario_names[ss]
            burden_df = rbind(burden_df, final_to_present_row)
          }
        } else{
          final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
          for(ss in 2:length(scenario_names)){
            final_to_present_row$scenario = scenario_names[ss]
            burden_df = rbind(burden_df, final_to_present_row)
          }
        }
      }
    }
      
    
    
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # if plotting burden relative to specified year, calculate relative values
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    if (!is.na(relative_year) & time_string=='annual'){
      # if the reference year is in the to-present simulation, use the same reference for all scenarios
      # if the reference year is not in the to-present simulation, use average value across scenarios and check that all scenarios have similar values for the reference year (if they do not, send a warning)
      if(('to-present' %in% burden_df$scenario) & (relative_year %in% unique(burden_df$year[burden_df$scenario=='to-present']))){
        # get the burden in the reference year
        reference_burden_cur = burden_df$mean_burden[burden_df$scenario=='to-present' & burden_df$year == relative_year]
      } else{
        similarity_threshold = 0.1
        all_ref_year_burdens = burden_df$mean_burden[burden_df$year == relative_year]
        reference_burden_cur = mean(all_ref_year_burdens)
        if(any(all_ref_year_burdens>(reference_burden_cur*(1+similarity_threshold))) | any(all_ref_year_burdens<(reference_burden_cur*(1-similarity_threshold)))){
          warning(paste0('in the reference year, some scenarios have different burdens for ',burden_metric_name))
        }
      }
      # calculate all relative burden values as mean_burden / reference_burden_cur: this will be referred to as 'burden relative to burden in relative_year'
      burden_df$mean_burden = burden_df$mean_burden / reference_burden_cur
      burden_df$max_burden = NA
      burden_df$min_burden = NA
        
      ylab_add_component = paste0('\n relative to ', relative_year)
      relative_string = paste0('_relativeTo', relative_year)
    } else{
      ylab_add_component = ''
      relative_string = ''
    }
      
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # create scenario-comparison plots
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    
    # subset to relevant scenarios currently being compared
    burden_df = burden_df[burden_df$scenario %in% scenario_names,]
    # get factors in the correct order (rather than alphabetical)
    burden_df$scenario = factor(burden_df$scenario, levels=rev(scenario_names))
    
    if(is.na(scenario_linetypes[1])){
      scenario_linetypes = rep(1, length(unique(burden_df$scenario)))
      names(scenario_linetypes) = unique(burden_df$scenario)
    }
    # ----- plot malaria burden ----- #
    if(length(ymax_each_burden)<bb){
      ylim_max = NA
    } else{
      ylim_max = ymax_each_burden[bb]
    }
    if(plot_by_month){
      gg_list[[bb]] = ggplot(burden_df, aes(x=as.Date(date), y=mean_burden, color=scenario)) +
        geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
        scale_fill_manual(values = rev(scenario_palette)) + 
        geom_line(linewidth=1) + 
        scale_color_manual(values = rev(scenario_palette)) + 
        xlab('date') + 
        ylab(paste0(gsub('\\(births\\)', '', burden_metric_name),ylab_add_component)) + 
        xlim(as.Date(paste0(min_year, '-01-01')), as.Date(paste0(max_year, '-01-01'))) +
        coord_cartesian(ylim=c(0, ylim_max)) +
        theme_classic()+ 
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text=element_text(size = text_size*font_scale), legend.text=element_text(size = text_size*legend_scale*font_scale), legend.key.size=unit(legend_scale*font_scale, 'lines'))
    } else{
      gg_list[[bb]] = ggplot(burden_df, aes(x=year, y=mean_burden, color=scenario, linetype=scenario))  +
        geom_line(linewidth=1) + 
        scale_linetype_manual(values=rev(scenario_linetypes)) +
        scale_color_manual(values = rev(scenario_palette)) + 
        xlab('year') +
        ylab(paste0(gsub('\\(births\\)', '', burden_metric_name), ylab_add_component)) +
        scale_x_continuous(breaks = function(b){ br <- scales::breaks_extended()(b); ib <- br[br == floor(br)]; if(length(ib)) ib else seq(floor(b[1]), ceiling(b[2])) }) +   # whole years only, but spaced adaptively (every 1/2/5/... yr) so long series don't crowd; avoids pretty_breaks()'s 2026.5 half-years
        coord_cartesian(xlim=if(is.null(display_xlim)) c(min_year, max_year) else display_xlim, ylim=c(0, ylim_max)) +   # x-range via coord (display_xlim crops the view without dropping data); using xlim() + scale_x_continuous() together triggered 'Scale for x already present'

        theme_classic()+ 
        theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text=element_text(size = text_size*font_scale), legend.text=element_text(size = text_size*legend_scale*font_scale), legend.key.size=unit(legend_scale*font_scale, 'lines'))
    }
    if(burden_metric == 'PfPR'){
      gg_list[[bb]] = gg_list[[bb]] + scale_y_continuous(labels = percent_format(accuracy = 1))
    }
    if(plot_CI){
      gg_list[[bb]] =  gg_list[[bb]] +
        geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
        scale_fill_manual(values = scenario_palette)
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] + theme(legend.position='none', text=element_text(size =separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_', relative_string, burden_metric_name,'_',district_subset,'.png'), separate_plot, dpi=600, width=4, height=3, units='in')  # width=4*0.9, height=3*0.7, # width=4*0.7, height=3*1.2,
    }
  }
  # gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]])), (ggplot() + theme_void())), gg_list)
  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  # remove legend from main plots
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position = "none")  + theme(text = element_text(size = text_size*font_scale))
  }
  # ----- combine all burden plots ----- #
  # gg = grid.arrange(grobs = gg_list, layout_matrix = matrix(c(1,1,2:(length(burden_colnames)+1)), ncol=2, byrow=TRUE))  # other orientation
  nrow_plot = 2
  if(include_U1) nrow_plot = 3
  num_in_matrix = ceiling(length(burden_colnames)/nrow_plot)*nrow_plot
  gg = grid.arrange(grobs = gg_list, layout_matrix = rbind(matrix(rep(1, ceiling(length(burden_colnames)/nrow_plot)), nrow=1), matrix(2:(num_in_matrix+1), nrow=nrow_plot, byrow=FALSE)))
  
  if(save_plots){
    ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_burden', relative_string,'_pyr', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg, dpi=600, width=9, height=3*nrow_plot, units='in')
  }
  
  return(gg)
}





######################################################################
# create state grid plot with timeseries of burden, no intervention info
######################################################################

plot_simulation_output_burden_by_state = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations,
                                             min_year, max_year, sim_end_years, relative_year=NA,
                                             scenario_filepaths, scenario_names, experiment_names, scenario_palette, group_col='State', group_levels=NULL, LLIN2y_flag=FALSE, overwrite_files=FALSE,
                                             extend_past_timeseries_year=NA, scenario_linetypes=NA,filename_suffix='',
                                             font_scale=1, legend_scale=0.6, save_width=12*0.9, save_height=8.5*0.9){
  # font_scale multiplies all in-plot text; legend_scale shrinks the legend; save_width/save_height
  # set the saved PNG size. Defaults preserve prior behaviour (except the smaller legend).
  # group_col: column used to group/facet LGAs (default 'State' -> geofacet via grid_layout_state_locations;
  #   any other value, e.g. 'Funder', uses facet_wrap with fixed shared scales). Filenames use tolower(group_col).
  if (!is.na(relative_year)){
    if(relative_year<min_year){
      warning('specified minimum year must be <= relative year. Setting min_year to relative_year.')
      min_year = relative_year
  }}
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))


  
  # iterate through scenarios, storing relevant output
  burden_df_all = data.frame()
  for(ee in 1:length(scenario_filepaths)){
    exp_filepath = scenario_filepaths[ee]
    exp_name = scenario_names[ee]
    cur_sim_output_agg = get_burden_timeseries_by_state(exp_filepath=exp_filepath, exp_name=exp_name, pop_filepath=pop_filepath, group_col=group_col, overwrite_files=overwrite_files)
    if(nrow(burden_df_all)==0){
      burden_df_all = cur_sim_output_agg
    } else{
      burden_df_all = rbind(burden_df_all, cur_sim_output_agg)
    }
  }

  # subset to relevant scenarios currently being compared
  burden_df_all = burden_df_all[burden_df_all$scenario %in% scenario_names,]

  # ----- malaria burden ----- #
  burden_metrics = c( 'PfPR_U5', 'PfPR_all', 'incidence_pp_U5', 'incidence_pp_all',  'total_mortality_pp_U5', 'total_mortality_pp_all')  # ,  'direct_mortality_pp_U5', 'direct_mortality_pp_all'
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence per person (U5)', 'incidence per person (all ages)', 'total mortality per person (U5)','total mortality per person (all ages)')  # ,  'direct mortality per person (U5)','direct mortality per person (all ages)'
  
  for(bb in 1:length(burden_metrics)){
    burden_metric = burden_metrics[bb]
    burden_metric_name = burden_metric_names[bb]
    burden_df = burden_df_all[,c(group_col,'year', 'scenario', burden_metric)]
    burden_df$mean_burden = burden_df[[burden_metric]]

    # connect the 'to-present' and 'future-projection' simulations in the plot. Two alternatives for how this is done, controlled by extend_past_timeseries:
    #   - (FALSE) extend the 'future-projection' lines all back to the end of the 'to-present' simulations, which is desirable if the future projection scenarios separate right away
     if('to-present' %in% burden_df$scenario){
        # add the final 'to-present' row to all future simulations for a continuous plot
        to_present_df = burden_df[burden_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          burden_df = rbind(burden_df, final_to_present_row)
        }
    }
    
    
    
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # if plotting burden relative to specified year, calculate relative values
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    if (!is.na(relative_year)){
      # if the reference year is in the to-present simulation, use the same reference for all scenarios
      # if the reference year is not in the to-present simulation, use average value across scenarios and check that all scenarios have similar values for the reference year (if they do not, send a warning)
      if(('to-present' %in% burden_df$scenario) & (relative_year %in% unique(burden_df$year[burden_df$scenario=='to-present']))){
        # get the burden in the reference year
        reference_burden_cur = burden_df[burden_df$scenario=='to-present' & burden_df$year == relative_year, c(group_col,'mean_burden')]
      } else{
        similarity_threshold = 0.1
        all_ref_year_burdens = burden_df[burden_df$year == relative_year, c(group_col,'mean_burden')]
        reference_burden_cur = all_ref_year_burdens %>% group_by(across(all_of(group_col))) %>%
          summarise(mean_burden = mean(mean_burden)) %>% ungroup()
        # if(any(all_ref_year_burdens>(reference_burden_cur*(1+similarity_threshold))) | any(all_ref_year_burdens<(reference_burden_cur*(1-similarity_threshold)))){
        #   warning(paste0('in the reference year, some scenarios have different burdens for ',burden_metric_name))
        # }
      }
      # calculate all relative burden values as mean_burden / reference_burden_cur: this will be referred to as 'burden relative to burden in relative_year'
      colnames(reference_burden_cur)[colnames(reference_burden_cur)=='mean_burden'] = 'ref_mean_burden'
      burden_df = merge(burden_df, reference_burden_cur, all=TRUE)
      burden_df$mean_burden = burden_df$mean_burden / burden_df$ref_mean_burden
      
      burden_ylab = paste0(gsub(' per person', '', burden_metric_name), ' relative to ', relative_year)
      relative_string = paste0('_relativeTo', relative_year)
    } else{
      burden_ylab = burden_metric_name
      relative_string = ''
    }
    
      
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    # create scenario-comparison plots
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
    
    # get factors in the correct order (rather than alphabetical)
    burden_df$scenario = factor(burden_df$scenario, levels=rev(scenario_names))
    
    if(is.na(scenario_linetypes[1])){
      scenario_linetypes = rep(1, length(unique(burden_df$scenario)))
      names(scenario_linetypes) = unique(burden_df$scenario)
    }

    burden_df$code = burden_df[[group_col]]
    if(!is.null(group_levels)) burden_df$code = factor(burden_df$code, levels=group_levels)  # custom facet order
    gg = ggplot(burden_df, aes(x=year, y=mean_burden, color=scenario, linetype=scenario))+
      geom_line(linewidth=1) + 
      scale_linetype_manual(values=scenario_linetypes) +
      scale_color_manual(values = scenario_palette) +
      xlab('year') + 
      ylab(burden_ylab) + 
      coord_cartesian(ylim=c(0, ifelse(!is.na(relative_year),2,NA)), xlim=c(min_year, max_year)) + 
      # scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      # scale_y_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      # skip first and last year so adjacent facets' labels don't run together
      scale_x_continuous(breaks = unique(c(min_year + 1, max_year - 1))) +
      scale_y_continuous(n.breaks= 3, labels = if(grepl('PfPR', burden_metric)) percent_format(accuracy = 1) else waiver()) +
      theme_bw()+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text=element_text(size = text_size*font_scale), legend.text=element_text(size = text_size*legend_scale*font_scale), legend.key.size=unit(legend_scale*font_scale, 'lines')) +
      (if(group_col == 'State') facet_geo(~code, grid = grid_layout_state_locations, label="name") else facet_wrap(~code))#, scales='free')

      ggsave(paste0(sim_future_output_dir, '/_plots/Timeseries_burden',relative_string,'_',tolower(group_col),'_grid_',burden_metric,filename_suffix,'.png'), gg, dpi=600, width=save_width, height=save_height, units='in')
  }
  return(gg)
}





######################################################################
# create csv of outputs for partners with timeseries of burden across metrics and years and scenarios
######################################################################
# create csv with annual average burden across specified timeseries for each state and scenario (requested by Maikore)
# columns are scenario, State, year, PfPR_all, PfPR_U5, incidence_all, incidence_U5
create_csv_timeseries_state_burden_each_scenario = function(sim_future_output_dir, pop_filepath,
                                                  min_year, max_year, 
                                                  scenario_filepaths, scenario_names, experiment_names, LLIN2y_flag=FALSE, overwrite_files=FALSE, 
                                                  filename_suffix=''){
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  
  # iterate through scenarios, storing relevant output
  burden_df_all = data.frame()
  for(ee in 1:length(scenario_filepaths)){
    exp_filepath = scenario_filepaths[ee]
    exp_name = scenario_names[ee]
    cur_sim_output_agg = get_burden_timeseries_by_state(exp_filepath=exp_filepath, exp_name=exp_name, pop_filepath=pop_filepath, overwrite_files=overwrite_files)
    if(nrow(burden_df_all)==0){
      burden_df_all = cur_sim_output_agg
    } else{
      burden_df_all = rbind(burden_df_all, cur_sim_output_agg)
    }
  }
  
  # subset to relevant scenarios currently being compared
  burden_df_all = burden_df_all %>% filter(scenario %in% scenario_names, year>=min_year, year<=max_year)
  
  # ----- malaria burden metrics to include ----- #
  burden_metrics = c( 'PfPR_U5', 'PfPR_all', 'incidence_pp_U5', 'incidence_pp_all')#,  'total_mortality_pp_U5', 'total_mortality_pp_all')  # ,  'direct_mortality_pp_U5', 'direct_mortality_pp_all'
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence per person (U5)', 'incidence per person (all ages)', 'total mortality per person (U5)','total mortality per person (all ages)')  # ,  'direct mortality per person (U5)','direct mortality per person (all ages)'
  colnames_include = c('scenario','State','year', burden_metrics)
  burden_df_all = burden_df_all %>% dplyr::select(all_of(colnames_include)) %>%
    mutate(incidence_U5 = incidence_pp_U5 * 1000,
           incidence_all = incidence_pp_all * 1000) %>%
    dplyr::select(-c(incidence_pp_U5, incidence_pp_all))
  
  write.csv(burden_df_all, paste0(sim_future_output_dir, '/_plots/timeseries_dfs/Timeseries_state_burden_each_scenario_', min_year, '_', max_year, filename_suffix,'.csv'), row.names=FALSE)
}






######################################################################
# create plot panel with selected burden metric and intervention info
######################################################################
# note: plot of ITN use rates through time is for the entire population (always shows all-age, even when burden plot shows U5)

plot_simulation_intervention_output = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
                                               plot_by_month, min_year, max_year, sim_end_years, 
                                               burden_metric, age_plotted, 
                                               pyr, chw_cov,
                                               scenario_filepaths, scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                                               indoor_protection_fraction=0.75, remove_exp_name_substring='', LLIN2y_flag=FALSE, overwrite_files=FALSE){
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  pop_sizes = read.csv(pop_filepath)
  pop_sizes = pop_sizes[,c('admin_name','pop_size')]
  # if we include all admins, get list of names from population size dataframe
  if(cur_admins[1] == 'all'){
    cur_admins = unique(pop_sizes$admin_name)
  }
  
  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  
  # ----- malaria burden ----- #
  
  # Get output column name for specified burden metric
  # Note: need to divide by pop size and multiply by 1000 if not PfPR
  burden_colname = NA
  if(burden_metric == 'PfPR'){
    if(age_plotted == 'U5'){
      burden_colname = 'PfPR_U5'
    } else if(age_plotted == 'all'){
      burden_colname = 'PfPR_MiP_adjusted'
    }
  } else if(burden_metric == 'incidence'){
    if(age_plotted == 'U5'){
      burden_colname = 'New_clinical_cases_U5'
    } else if(age_plotted == 'all'){
      burden_colname = 'New_Clinical_Cases'
    }
  } else if(burden_metric == 'directMortality'){
    if(age_plotted == 'U5'){
      burden_colname = 'direct_mortality_nonMiP_U5_mean'
    } else if(age_plotted == 'all'){
      burden_colname = 'direct_mortality_nonMiP_mean'
    }
  } else if(burden_metric == 'allMortality'){
    if(age_plotted == 'U5'){
      burden_colname = 'total_mortality_U5_mean'
    } else if(age_plotted == 'all'){
      burden_colname = 'total_mortality_mean'
    }
  } 
  if(is.na(burden_colname)){
    warning('PROBLEM DETECTED: name of burden metric or age group not currently supported')
  }
  
  # (the combined burden / intervention-coverage timeseries below are built in memory each call,
  #  not saved to or reloaded from timeseries_dfs -- see the if(FALSE) notes below)
  if(LLIN2y_flag){
    llin2y_string = '_2yLLIN'
  } else{
    llin2y_string = ''
  }
  burden_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_burden_',time_string,'Timeseries_', burden_metric, '_', age_plotted, '_pyr', pyr, '_', chw_cov, 'CHW_',district_subset, llin2y_string,'.csv')
  if(FALSE){  # always regenerate in memory; never reload a shared scenario-combined df from
              # timeseries_dfs (those are keyed by district_subset, NOT the scenario set, so reloading
              # would silently mix scenarios across runs/sections, e.g. combo base vs _vacc2).
    burden_df = read.csv(burden_df_filepath)
  } else{
    # iterate through scenarios, storing relevant output
    burden_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_sim_output_agg = get_burden_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                     exp_name = scenario_names[ee],  district_subset=district_subset,
                                                     cur_admins=cur_admins, pop_sizes = pop_sizes, min_year=min_year, max_year=max_year, burden_colname=burden_colname, age_plotted=age_plotted, plot_by_month=plot_by_month)
      if(nrow(burden_df)==0){
        burden_df = cur_sim_output_agg
      } else{
        burden_df = merge(burden_df, cur_sim_output_agg, all=TRUE)
      }
    }
    
    # add the final 'to-present' row to all future simulations for a continuous plot
    to_present_df = burden_df[burden_df$scenario == 'to-present',]
    if(nrow(to_present_df)>0){
      if(plot_by_month){
        final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          burden_df = rbind(burden_df, final_to_present_row)
        }
      } else{
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          burden_df = rbind(burden_df, final_to_present_row)
        }
      }
    }
    # (not saved: regenerated each call from the per-experiment caches, to avoid the
    #  district_subset-keyed shared-cache hazard)
  }




  # ----- LLIN, vaccine, and IRS intervention coverage ----- #
  
  # build the LLIN/IRS coverage timeseries in memory each call (not saved/reloaded -- see if(FALSE))
  llin_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_llin_irs_',time_string,'Timeseries', '_pyr', pyr, '_', chw_cov, 'CHW_',district_subset, llin2y_string,'.csv')
  if(FALSE){  # always regenerate in memory (see note above); never reload the shared df
    net_use_df = read.csv(llin_df_filepath)
  } else{
    # iterate through scenarios, storing relevant output
    net_use_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      cur_net_agg = get_intervention_use_timeseries_exp(exp_filepath = scenario_filepaths[ee],
                                                        exp_name = scenario_names[ee], 
                                                        cur_admins=cur_admins, pop_sizes=pop_sizes, min_year=min_year, max_year=max_year, indoor_protection_fraction=indoor_protection_fraction, plot_by_month=plot_by_month)
      if(nrow(net_use_df)==0){
        net_use_df = cur_net_agg
      } else{
        net_use_df = merge(net_use_df, cur_net_agg, all=TRUE)
      }
    }
    
    # first, remove the final 'to-present' month or year - it should was overwritten in the pick-up from burn-in
    # then, add the final 'to-present' row to all future simulations for a continuous plot
    if(plot_by_month){
      if(any(net_use_df$scenario == 'to-present') & length(scenario_names)>1){
        # remove excess month from to-present simulation
        max_to_present_date = max(net_use_df$date[net_use_df$scenario == 'to-present'])
        row_to_remove = intersect(which(net_use_df$scenario == 'to-present'), which(net_use_df$date == max_to_present_date))
        net_use_df = net_use_df[-row_to_remove,]
        
        # join past and future simulation trajectories
        to_present_df = net_use_df[net_use_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
        
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          net_use_df = rbind(net_use_df, final_to_present_row)
        }
      }
    } else{
      # remove excess year from to-present simulation
      if(any(net_use_df$scenario == 'to-present') & length(scenario_names)>1){
        max_to_present_date = max(net_use_df$year[net_use_df$scenario == 'to-present'])
        min_projection_date = min(net_use_df$year[net_use_df$scenario != 'to-present'])
        if(max_to_present_date >= min_projection_date){
          row_to_remove = intersect(which(net_use_df$scenario == 'to-present'), which(net_use_df$year == max_to_present_date))
          net_use_df = net_use_df[-row_to_remove,]
        }

        # join past and future simulation trajectories
        to_present_df = net_use_df[net_use_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          net_use_df = rbind(net_use_df, final_to_present_row)
        }
      }
    }
    # (not saved: regenerated each call -- shared df is keyed by district_subset, not scenario set)
  }
  
  
  
  
  # ----- Case management ----- #
  
  # build the CM coverage timeseries in memory each call (not saved/reloaded -- see if(FALSE))
  cm_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_cm_',time_string,'Timeseries', '_pyr', pyr, '_', chw_cov, 'CHW_',district_subset, llin2y_string,'.csv')
  if(FALSE){  # always regenerate in memory; never reload the shared (district_subset-keyed) df
    cm_df = read.csv(cm_df_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    cm_df = data.frame()
    for(ee in 1:length(scenario_filepaths)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      if(experiment_intervention_name %in% intervention_file_info$ScenarioName){
        cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      } else{
        cur_int_row = which(gsub(remove_exp_name_substring, '', intervention_file_info$ScenarioName) == experiment_intervention_name)
      }
      
      # read in intervention files
      cm_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      
      cur_cm_agg = get_cm_timeseries_exp(cm_filepath=cm_filepath, pop_sizes=pop_sizes, end_year=end_year, exp_name = scenario_names[ee], 
                                         cur_admins=cur_admins, min_year=min_year, plot_by_month=plot_by_month)
      
      if(nrow(cm_df)==0){
        cm_df = cur_cm_agg
      } else{
        cm_df = rbind(cm_df, cur_cm_agg)
      }
    }
    
    # add the final 'to-present' row to all future simulations for a continuous plot
    if(plot_by_month){
      # join past and future simulation trajectories
      to_present_df = cm_df[cm_df$scenario == 'to-present',]
      final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
      for(ss in 2:length(scenario_names)){
        final_to_present_row$scenario = scenario_names[ss]
        cm_df = rbind(cm_df, final_to_present_row)
      }
    } else{
      # join past and future simulation trajectories
      to_present_df = cm_df[cm_df$scenario == 'to-present',]
      if(nrow(to_present_df)>0){
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          cm_df = rbind(cm_df, final_to_present_row)
        }
      }
    }
    # (not saved: regenerated each call -- shared df is keyed by district_subset, not scenario set)
  }
  
  
  
  
  
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  burden_df$scenario = factor(burden_df$scenario, levels=scenario_names)
  
  # ----- malaria burden ----- #
  
  if(plot_by_month){
    g_burden = ggplot(burden_df, aes(x=as.Date(date), y=mean_burden, color=scenario)) +
      geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0(burden_metric, ' - ', age_plotted)) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size))
  } else{
    g_burden = ggplot(burden_df, aes(x=year, y=mean_burden, color=scenario)) +
      geom_ribbon(aes(ymin=min_burden, ymax=max_burden, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0(burden_metric, ' - ', age_plotted)) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), text = element_text(size = text_size), legend.text=element_text(size = text_size))
  }
  
  inter_plot_list = list()
  # ----- LLIN use and distribution ----- #
  # plot net use through time
  if(plot_by_month){
    g_net_use = ggplot(net_use_df, aes(x=as.Date(date), y=coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('LLIN use (all ages)')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  } else{
    g_net_use = ggplot(net_use_df, aes(x=year, y=coverage, color=scenario)) +
      # geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      # geom_hline(yintercept=0.22, alpha=0.1)+
      # geom_hline(yintercept=0.39, alpha=0.1)+
      xlab('year') + 
      ylab(paste0('LLIN use (all ages)')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
    
    g_all_inter = ggplot() +
      geom_line(data=net_use_df, aes(x=year, y=coverage), color=rgb(1,0.6,1), linewidth=1) + 
      xlab('year') + 
      ylab(paste0('coverage metric')) + 
      coord_cartesian(ylim=c(0,NA))+
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
  }
  inter_plot_list = append(inter_plot_list, list(g_net_use))
  # # plot net distribution numbers through time (how many nets distributed in each month or year per person?)
  # if(plot_by_month){
  #   g_net_dist = ggplot(net_use_df, aes(x=as.Date(date), y=new_net_per_cap, color=scenario)) +
  #     geom_point(size=1) + 
  #     scale_color_manual(values = scenario_palette) + 
  #     xlab('date') + 
  #     ylab(paste0('LLINs distributed per person')) + 
  #     theme_classic()+ 
  #     theme(legend.position = "none", text = element_text(size = text_size))
  # } else{
  #   g_net_dist = ggplot(net_use_df, aes(x=year, y=new_net_per_cap, color=scenario)) +
  #     geom_point(size=2) + 
  #     geom_line(alpha=0.2, linewidth=2) +
  #     scale_color_manual(values = scenario_palette) + 
  #     xlab('year') + 
  #     ylab(paste0('LLINs distributed per person')) + 
  #     theme_classic()+ 
  #     theme(legend.position = "none", text = element_text(size = text_size))
  # }
  
  
  # ----- vaccine ----- #
  if('vacc_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_vacc = ggplot(net_use_df, aes(x=as.Date(date), y=vacc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('Vaccines (primary series + booster) per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_vacc = ggplot(net_use_df, aes(x=year, y=vacc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('Vaccines (primary series + booster) per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=vacc_per_cap), color=rgb(0,0.3,0), linewidth=1)
        
    }
    inter_plot_list = append(inter_plot_list, list(g_vacc))
  }
  
  
  # ----- PMC ----- #
  if('pmc_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_pmc = ggplot(net_use_df, aes(x=as.Date(date), y=pmc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('PMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_pmc = ggplot(net_use_df, aes(x=year, y=pmc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('PMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=pmc_per_cap), color=rgb(0.0,0.4,1), linewidth=1)
    }
    inter_plot_list = append(inter_plot_list, list(g_pmc))
  }
  
  
  # ----- SMC ----- #
  if('smc_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_smc = ggplot(net_use_df, aes(x=as.Date(date), y=smc_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('SMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_smc = ggplot(net_use_df, aes(x=year, y=smc_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('SMC doses per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=smc_per_cap), color=rgb(0.0,0.4,1), linewidth=1)
    }
    inter_plot_list = append(inter_plot_list, list(g_smc))
  }
  
  
  # ----- IRS ----- #
  if('irs_per_cap' %in% colnames(net_use_df)){
    if(plot_by_month){
      g_irs = ggplot(net_use_df, aes(x=as.Date(date), y=irs_per_cap, color=scenario)) +
        geom_point(size=1) + 
        scale_color_manual(values = scenario_palette) + 
        xlab('date') + 
        ylab(paste0('IRS rounds per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
    } else{
      g_irs = ggplot(net_use_df, aes(x=year, y=irs_per_cap, color=scenario)) +
        geom_point(size=2) + 
        geom_line(alpha=0.2, linewidth=2) +
        scale_color_manual(values = scenario_palette) + 
        xlab('year') + 
        ylab(paste0('IRS per person')) + 
        coord_cartesian(ylim=c(0,NA))+
        theme_classic()+ 
        theme(legend.position = "none", text = element_text(size = text_size))
      
      g_all_inter = g_all_inter +
        geom_line(data=net_use_df, aes(x=year, y=irs_per_cap), color=rgb(1,0,1), linewidth=1)
    }
    inter_plot_list = append(inter_plot_list, list(g_irs))
  }
  
  
  # ----- Case management ----- #
  if(plot_by_month){
    g_cm = ggplot(cm_df, aes(x=as.Date(date), y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(ylim=c(0,NA))+
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
      theme_classic()+ 
      theme(legend.position = "none", text = element_text(size = text_size))
    
    g_all_inter = g_all_inter +
      geom_line(data=cm_df, aes(x=year, y=mean_coverage), color=rgb(0.1,0.9,0.4), linewidth=1)
  }
  inter_plot_list = append(inter_plot_list, list(g_cm))
  
  
  
  # ----- combine burden and intervention plots ----- #
  gg_leg = ggpubr::as_ggplot(ggpubr::get_legend(g_burden))
  g_burden = g_burden + theme(legend.position = "none")
  gg = plot_grid(plotlist=append(list(gg_leg, g_burden), inter_plot_list), ncol=1, nrow=(2+length(inter_plot_list)), align='vh', axis='lrtb')  # (gg_leg, g_burden, plot_list)
  
  if(save_plots){
    ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_', burden_metric, '_', age_plotted, '_versusInterventions_pyr', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), gg, dpi=600, width=7, height=4*(2+length(inter_plot_list)), units='in')
    
    if(!plot_by_month){
      ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_interventions_pyr', pyr, '_', chw_cov, 'CHW_',district_subset,'.png'), g_all_inter, dpi=600, width=7, height=5, units='in')
    }
  }
  return(gg)
}





######################################################################
# create plot panel with CM timeseries for included scenarios in each state
######################################################################

plot_state_grid_cm = function(sim_future_output_dir, pop_filepath, grid_layout_state_locations, 
                               plot_by_month, min_year, max_year, sim_end_years, 
                               scenario_names, scenario_input_references, experiment_names, scenario_palette, 
                               overwrite_files=FALSE){

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # combine simulation output from multiple scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  admin_info = read.csv(pop_filepath)
  admin_info = admin_info[,c('admin_name','pop_size','State')]

  # create output directories
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots'))) dir.create(paste0(sim_future_output_dir, '/_plots'))
  if(!dir.exists(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))) dir.create(paste0(sim_future_output_dir, '/_plots/timeseries_dfs'))
  if(plot_by_month){
    time_string = 'monthly'
  } else time_string = 'annual'
  
  # build the CM-by-state timeseries in memory each call (not saved/reloaded -- see if(FALSE))
  cm_df_filepath = paste0(sim_future_output_dir, '/_plots/timeseries_dfs/df_cm_state_',time_string,'Timeseries.csv')
  if(FALSE){  # always regenerate in memory; never reload the shared (district_subset-keyed) df
    cm_df = read.csv(cm_df_filepath)
  } else{
    # iterate through scenarios, storing input CM coverages
    cm_df = data.frame()
    for(ee in 1:length(experiment_names)){
      intervention_csv_filepath = scenario_input_references[ee]
      intervention_file_info = read.csv(intervention_csv_filepath)
      experiment_intervention_name = experiment_names[ee]
      end_year = sim_end_years[ee]
      cur_int_row = which(intervention_file_info$ScenarioName == experiment_intervention_name)
      # read in intervention files
      cm_filepath = paste0(hbhi_dir, '/simulation_inputs/', intervention_file_info$CM_filename[cur_int_row], '.csv')
      
      cur_cm_agg = get_cm_timeseries_by_state(cm_filepath=cm_filepath, admin_info=admin_info, end_year=end_year, exp_name = scenario_names[ee], 
                                              min_year=min_year, plot_by_month=plot_by_month)
      
      if(nrow(cm_df)==0){
        cm_df = cur_cm_agg
      } else{
        cm_df = rbind(cm_df, cur_cm_agg)
      }
    }
    
    if(any(grepl('to-present', cm_df$scenario))){
      # add the final 'to-present' row to all future simulations for a continuous plot
      if(plot_by_month){
        # join past and future simulation trajectories
        to_present_df = cm_df[cm_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[as.Date(to_present_df$date) == max(as.Date(to_present_df$date)),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          cm_df = rbind(cm_df, final_to_present_row)
        }
      } else{
        # join past and future simulation trajectories
        to_present_df = cm_df[cm_df$scenario == 'to-present',]
        final_to_present_row = to_present_df[to_present_df$year == max(to_present_df$year),]
        for(ss in 2:length(scenario_names)){
          final_to_present_row$scenario = scenario_names[ss]
          cm_df = rbind(cm_df, final_to_present_row)
        }
      }
    }
    # (not saved: regenerated each call -- shared df is keyed by district_subset, not scenario set)
  }
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # create scenario-comparison plots
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  # get factors in the correct order (rather than alphabetical)
  cm_df$scenario = factor(cm_df$scenario, levels=rev(scenario_names))
  cm_df$code = cm_df$State
  
  if(plot_by_month){
    g_cm = ggplot(cm_df, aes(x=as.Date(date), y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('date') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year))+
      theme_bw()+ 
      theme(legend.position = "none", text = element_text(size = text_size))+
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  } else{
    g_cm = ggplot(cm_df, aes(x=year, y=mean_coverage, color=scenario)) +
      geom_ribbon(aes(ymin=min_coverage, ymax=max_coverage, fill=scenario), alpha=0.1, color=NA)+
      scale_fill_manual(values = scenario_palette) + 
      geom_line(linewidth=1) + 
      scale_color_manual(values = scenario_palette) + 
      xlab('year') + 
      ylab(paste0('Effective treatment rate (U5)')) + 
      coord_cartesian(xlim=c(min_year, max_year))+
      scale_x_continuous(breaks= pretty_breaks(), guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()+ 
      # theme(legend.position = "none", text = element_text(size = text_size))+
      theme(legend.position = "top", legend.box='horizontal', legend.title = element_blank(), legend.text=element_text(size = text_size)) +  # legend.position = "none"
      facet_geo(~code, grid = grid_layout_state_locations, label="name", scales='free') 
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/',time_string,'Timeseries_CM_by_state.png'), g_cm, dpi=600, width=12, height=10, units='in')
}














#####################################################################
# plot map of admin subsets
#####################################################################
plot_included_admin_map = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, admin_shapefile_filepath, shapefile_admin_colname,
                                   save_width=4.8, save_height=4.8, lga_linewidth=0.12, state_linewidth=0.5,
                                   lga_border_color='grey50', state_border_color='black'){
  # save_width/save_height set the saved PNG size (map is theme_void(), so no font_scale needed).
  # LGA boundaries are drawn thin (lga_linewidth, lga_border_color); STATE boundaries are overlaid
  # thicker/darker (state_linewidth, state_border_color) by dissolving the LGA polygons on the
  # shapefile's own 'State' column -- so the LGA-level shapefile is sufficient (no separate state file).
  admin_pop = read.csv(pop_filepath)
  admin_shapefile = st_read(admin_shapefile_filepath)
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=admin_pop$admin_name, origin_names=admin_shapefile[[shapefile_admin_colname]])

  admin_in_map = data.frame(admin_name = admin_pop$admin_name, admin_included='no')
  admin_in_map$admin_included[admin_in_map$admin_name %in% cur_admins] = 'yes'
  included_colors = c('yes'='#006692', 'no'='grey96')

  # join inclusion status onto the LGA polygons. The shapefile already carries its own 'State' column,
  # so do NOT merge State from pop here -- that collides (-> State.x/State.y) and breaks the dissolve.
  admin_cur = admin_shapefile %>%
    dplyr::left_join(admin_in_map, by=c('NOMDEP' = 'admin_name'))

  gg_map = ggplot(admin_cur) +
    geom_sf(aes(fill=admin_included), linewidth=lga_linewidth, color=lga_border_color) +   # thin LGA borders
    scale_fill_manual(values=included_colors, drop=FALSE, na.value='grey96') +
    theme_void() +
    theme(legend.position = 'none')

  # overlay bolder state boundaries by dissolving the LGA polygons on the shapefile's 'State' column.
  # st_make_valid() first: the shapefile has minor invalid geometries that otherwise break st_union.
  if('State' %in% names(admin_cur)){
    state_sf = sf::st_make_valid(admin_cur[!is.na(admin_cur$State), ]) %>%
      dplyr::group_by(State) %>%
      dplyr::summarise(.groups = 'drop')
    gg_map = gg_map + geom_sf(data = state_sf, fill = NA, color = state_border_color, linewidth = state_linewidth)
  }
  ggsave(paste0(sim_future_output_dir, '/_plots/map_admins_included_', district_subset, '.png'), gg_map, dpi=600, width=save_width, height=save_height, units='in')
}



#####################################################################
# plot maps of burden with and without the intervention
#####################################################################
plot_burden_maps = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                            barplot_start_year, barplot_end_year,
                            pyr, chw_cov,
                            scenario_names, experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', LLIN2y_flag=FALSE,
                            filename_sufffix='', overwrite_files=FALSE){


  admin_pop = read.csv(pop_filepath)
  if(!(cur_admins[1] == 'all')){
    admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),]
  }
  admin_shapefile = shapefile(admin_shapefile_filepath)
  # standardize shapefile names
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=archetype_info$LGA, origin_names=admin_shapefile$NOMDEP)
  
  years_included = barplot_end_year - barplot_start_year + 1

  # burden metrics
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'death_rate_mean_U5', 'death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  # burden_colnames_for_map = c('pfpr_all', 'incidence_all', 'mortality_rate_all')
  # burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'mortality (all)')
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  burden_colnames_for_map = c('pfpr_u5', 'pfpr_all', 'incidence_u5', 'incidence_all', 'mortality_rate_u5',  'mortality_rate_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  

  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #   iterate through scenarios, creating dataframe including all burden metrics
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  num_scenarios = length(experiment_names)
  burden_df_all = data.frame()
  for(ee in 1:num_scenarios){
    experiment_name = experiment_names[ee]
    cur_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop, comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
    cur_burden_df$scenario_name = scenario_names[ee]
    if(nrow(burden_df_all) == 0){
      burden_df_all = cur_burden_df
    } else{
      burden_df_all = rbind(burden_df_all, cur_burden_df)
    }
  }



  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #      create maps showing each burden metric for all scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  if(LLIN2y_flag){
    llin2y_string = '_2yLLIN'
  } else{
    llin2y_string = ''
  }
  num_colors = 40
  # colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
  colorscale = colorRampPalette(brewer.pal(9, 'YlOrRd'))(num_colors)
  


  # iterate through burden metrics, creating plots for each
  for(cc in 1:length(burden_colnames_for_map)){

    if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_', burden_colnames_for_map[cc], '_', pyr, '_', chw_cov, 'CHW_', district_subset, llin2y_string, filename_sufffix, '.png'), res=600, width=(num_scenarios*3+2)*3/4, height=3, units='in')
    par(mar=c(0,1,2,0))
    # set layout for panel of maps
    layout_matrix = matrix(rep(c(rep(1:num_scenarios, each=3),rep((num_scenarios+1),2)),2), nrow=2, byrow=TRUE)
    layout(mat = layout_matrix)

    cur_colname = burden_colnames_for_map[cc]
    min_value = min(min(burden_df_all[[cur_colname]], na.rm=TRUE), 0)
    max_value = max(max(burden_df_all[[cur_colname]], na.rm=TRUE), 0.65)

    # iterate through scenarios
    for(ee in 1:num_scenarios){
      cur_burden_df = burden_df_all[burden_df_all$scenario_name == scenario_names[ee],]
      vals_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
      for (i_ds in 1:length(vals_ordered$ds_ordered)){
        cur_ds = vals_ordered$ds_ordered[i_ds]
        if(toupper(cur_ds) %in% toupper(cur_burden_df$admin_name)){
          vals_ordered$value[i_ds] = cur_burden_df[which(toupper(cur_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
        }
      }

      col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
      col_cur[is.na(col_cur)] = 'grey'
      plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=scenario_names[ee])
    }
    # legend
    legend_label_vals = seq(min_value, max_value, length.out=5)
    legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
    text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
    rasterImage(legend_image, 0, 0, 1,1)
    # fourth blank plot
    # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
    par(mfrow=c(1,1), mar=c(5,4,4,2))
    if(save_plots) dev.off()
  }

}





# plot map with the reduction in burden relative to the first scenario
plot_burden_relative_reduction_maps = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins,
                            barplot_start_year, barplot_end_year,
                            pyr, chw_cov,
                            scenario_names, experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', LLIN2y_flag=FALSE,
                            overwrite_files=FALSE){
  
  
  admin_pop = read.csv(pop_filepath)
  if(!(cur_admins[1] == 'all')){
    admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),]
  }
  admin_shapefile = shapefile(admin_shapefile_filepath)
  # standardize shapefile names
  admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=archetype_info$LGA, origin_names=admin_shapefile$NOMDEP)
  
  years_included = barplot_end_year - barplot_start_year + 1
  
  # burden metrics
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'death_rate_mean_U5', 'death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  # burden_colnames_for_map = c('pfpr_all', 'incidence_all', 'mortality_rate_all')
  # burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'mortality (all)')
  # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'direct_death_rate_mean_U5', 'direct_death_rate_mean_all', 'all_death_rate_mean_U5', 'all_death_rate_mean_all')
  # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'direct mortality (U5)', 'direct mortality (all ages)', 'mortality (U5)', 'mortality (all ages)')
  burden_colnames_for_map = c('pfpr_u5', 'pfpr_all', 'incidence_u5', 'incidence_all', 'mortality_rate_u5',  'mortality_rate_all')
  burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #   iterate through scenarios, creating dataframe including all burden metrics
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  num_scenarios = length(experiment_names)
  burden_df_all = data.frame()
  for(ee in 1:num_scenarios){
    experiment_name = experiment_names[ee]
    cur_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop, comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
    cur_burden_df$scenario_name = scenario_names[ee]
    if(nrow(burden_df_all) == 0){
      burden_df_all = cur_burden_df
    } else{
      burden_df_all = rbind(burden_df_all, cur_burden_df)
    }
  }
  
  
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #      create maps showing each burden metric for all scenarios
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  if(LLIN2y_flag){
    llin2y_string = '_2yLLIN'
  } else{
    llin2y_string = ''
  }
  num_colors = 40
  colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
  
  
  
  # iterate through burden metrics, creating plots for each
  for(cc in 1:length(burden_colnames_for_map)){
    
    if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_rel_reduction_', burden_colnames_for_map[cc], '_', pyr, '_', chw_cov, 'CHW_', district_subset, llin2y_string, '.png'), res=600, width=(num_scenarios*3+2)*3/4, height=3, units='in')
    par(mar=c(0,1,2,0))
    # set layout for panel of maps
    layout_matrix = matrix(rep(c(rep(1:num_scenarios, each=3),rep((num_scenarios+1),2)),2), nrow=2, byrow=TRUE)
    layout(mat = layout_matrix)
    
    cur_colname = burden_colnames_for_map[cc]
    
    # get reference column for this burden metric
    burden_df_ref = burden_df_all[burden_df_all$scenario_name==scenario_names[1],c('admin_name',cur_colname)]
    colnames(burden_df_ref)[which(colnames(burden_df_ref)==cur_colname)]='reference_value'
    burden_df_relative = merge(burden_df_all, burden_df_ref, all=TRUE)
    burden_df_relative$rel_reduction = (burden_df_relative$reference_value - burden_df_relative[[cur_colname]]) / burden_df_relative$reference_value
    
    # set minimum and maximum plotted in legend
    abs_max = max(abs( burden_df_relative$rel_reduction), na.rm=TRUE)
    min_value = -1 * abs_max
    max_value = abs_max
    
    # iterate through scenarios
    for(ee in 2:num_scenarios){
      cur_burden_df = burden_df_relative[burden_df_relative$scenario_name == scenario_names[ee],]
      vals_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
      for (i_ds in 1:length(vals_ordered$ds_ordered)){
        cur_ds = vals_ordered$ds_ordered[i_ds]
        if(toupper(cur_ds) %in% toupper(cur_burden_df$admin_name)){
          vals_ordered$value[i_ds] = cur_burden_df[which(toupper(cur_burden_df$admin_name) == toupper(cur_ds)), 'rel_reduction']
        }
      }
      
      col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
      col_cur[is.na(col_cur)] = 'grey'
      plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=scenario_names[ee])
    }
    # legend
    legend_label_vals = seq(min_value, max_value, length.out=5)
    legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
    text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
    rasterImage(legend_image, 0, 0, 1,1)
    # fourth blank plot
    # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
    par(mfrow=c(1,1), mar=c(5,4,4,2))
    if(save_plots) dev.off()
  }
  
}
# 
# 
# 
# 
# 
# #####################################################################
# # plot maps of burden with and without IPTi
# #####################################################################
# 
# line2user <- function(line, side) {
#   lh <- par('cin')[2] * par('cex') * par('lheight')
#   x_off <- diff(grconvertX(0:1, 'inches', 'user'))
#   y_off <- diff(grconvertY(0:1, 'inches', 'user'))
#   switch(side,
#          `1` = par('usr')[3] - line * y_off * lh,
#          `2` = par('usr')[1] - line * x_off * lh,
#          `3` = par('usr')[4] + line * y_off * lh,
#          `4` = par('usr')[2] + line * x_off * lh,
#          stop("side must be 1, 2, 3, or 4", call.=FALSE))
# }
# 
# 
# 
# plot_IPTi_burden_maps = function(sim_future_output_dir, pop_filepath, district_subset, cur_admins, 
#                                  barplot_start_year, barplot_end_year, 
#                                  experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', overwrite_files=FALSE){
#   
#   admin_pop = read.csv(pop_filepath)
#   admin_shapefile = shapefile(admin_shapefile_filepath)
#   
#   years_included = barplot_end_year - barplot_start_year + 1 
#   
#   # burden metrics
#   burden_metric_names = c('PfPR (U1)', 'incidence (U1)', 'mortality (U1)')
#   burden_colnames_for_map = c('pfpr_u1', 'incidence_u1', 'mortality_rate_u1')
#   
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      read in and format malaria burden simulation output in IPTi admins
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   experiment_names_descriptions = c('noIPTi', 'IPTi')
#   pfpr_u1_df = data.frame(admin_name=cur_admins)
#   deaths_u1_df = data.frame(admin_name=cur_admins)
#   clinical_cases_u1_df = data.frame(admin_name=cur_admins)
#   pop_u1_df = data.frame(admin_name=admin_pop$admin_name)
#   
#   # no-IPTi burden df
#   experiment_name = experiment_names[1]
#   option_name = experiment_names_descriptions[1]
#   noIPTi_burden_df = get_total_U1_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   # IPTi burden df
#   experiment_name = experiment_names[2]
#   option_name = experiment_names_descriptions[2]
#   IPTi_burden_df = get_total_U1_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   
#   
#   
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      create panel of maps showing all burden metrics
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_IPTi_burden_', pyr, '_', chw_cov, 'CHW.png'), res=600, width=6, height=3*length(burden_metric_names), units='in')
#   par(mar=c(0,1,2,0))
#   
#   num_colors = 40
#   colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
#   
#   # set layout for panel of maps
#   base_matrix = matrix(c(1,1,1,2,2,2,3,3, 1,1,1,2,2,2,3,3), nrow=2, byrow=TRUE)
#   # add rows for each burden metric
#   layout_matrix = base_matrix
#   for(cc in 2:length(burden_colnames_for_map)){
#     layout_matrix = rbind(layout_matrix, base_matrix + 3*(cc-1))
#   }
#   # add row for title
#   layout_matrix = layout_matrix + 1
#   layout_matrix = rbind(rep(1, ncol(layout_matrix)), layout_matrix)
#   layout(mat = layout_matrix)
#   
#   # title
#   plot.new()
#   text(0.5,0.5,"Malaria burden in each health district",cex=2.5,font=1)
#   # text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
#   #      line2user(line=4, side=3), "Malaria burden in each health district", xpd=NA, cex=2, font=2)
#   
#   # iterate through burden metrics, creating plots of each
#   for(cc in 1:length(burden_colnames_for_map)){
#     cur_colname = burden_colnames_for_map[cc]
#     vals_ordered_noipti = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     vals_ordered_ipti = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     for (i_ds in 1:length(vals_ordered_noipti$ds_ordered)){
#       cur_ds = vals_ordered_noipti$ds_ordered[i_ds]
#       if(toupper(cur_ds) %in% toupper(noIPTi_burden_df$admin_name)){
#         vals_ordered_noipti$value[i_ds] = noIPTi_burden_df[which(toupper(noIPTi_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#         vals_ordered_ipti$value[i_ds] = IPTi_burden_df[which(toupper(IPTi_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#       }
#     }
#     min_value = min(c(vals_ordered_noipti$value, vals_ordered_ipti$value), na.rm=TRUE)
#     max_value = max(c(vals_ordered_noipti$value, vals_ordered_ipti$value), na.rm=TRUE)
#     # without IPTi
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_noipti$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - without IPTi'))
#     # with IPTi
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_ipti$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - with IPTi'))
#     # legend
#     legend_label_vals = seq(min_value, max_value, length.out=5)
#     legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
#     plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
#     text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
#     rasterImage(legend_image, 0, 0, 1,1)
#     # fourth blank plot
#     # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
#     
#   }
#   par(mfrow=c(1,1), mar=c(5,4,4,2))
#   if(save_plots) dev.off()
#   
#   
# }
# 
# 
# 
# #####################################################################
# # plot maps of burden with and without the intervention
# #####################################################################
# 
# plot_burden_maps_with_without_inter = function(sim_future_output_dir, pop_filepath, cur_admins, district_subset,
#                                                barplot_start_year, barplot_end_year,
#                                                experiment_names, admin_shapefile_filepath, shapefile_admin_colname='NOMDEP', inter_name='PBO',
#                                                overwrite_files=FALSE){
#   
#   admin_pop = read.csv(pop_filepath)
#   admin_shapefile = shapefile(admin_shapefile_filepath)
#   
#   years_included = barplot_end_year - barplot_start_year + 1
#   
#   # burden metrics
#   # burden_colnames_for_map = c('average_PfPR_U5', 'average_PfPR_all', 'incidence_U5', 'incidence_all', 'death_rate_mean_U5', 'death_rate_mean_all')
#   # burden_metric_names = c('PfPR (U5)', 'PfPR (all ages)', 'incidence (U5)', 'incidence (all ages)', 'mortality (U5)', 'mortality (all ages)')
#   burden_colnames_for_map = c('pfpr_all', 'incidence_all', 'mortality_rate_all')
#   burden_metric_names = c('PfPR (all ages)', 'incidence (all ages)', 'mortality (all)')
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      read in and format malaria burden simulation output in specified set of admins
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   experiment_names_descriptions = paste0(c('no', ''), inter_name)
#   
#   # no-intervention burden df
#   experiment_name = experiment_names[1]
#   option_name = experiment_names_descriptions[1]
#   noPBO_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   # with intervention burden df
#   experiment_name = experiment_names[2]
#   option_name = experiment_names_descriptions[2]
#   PBO_burden_df = get_total_burden(sim_output_filepath=sim_future_output_dir, experiment_name=experiment_name, admin_pop=admin_pop[which(admin_pop$admin_name %in% cur_admins),], comparison_start_year=barplot_start_year, comparison_end_year=barplot_end_year, district_subset=district_subset, cur_admins=cur_admins, overwrite_files=overwrite_files)
#   
#   # # increase in burden in each admin: (without intervention - with intervention)  /  with intervention
#   # rel_burd_increase = (noPBO_burden_df$incidence_all - PBO_burden_df$incidence_all ) / PBO_burden_df$incidence_all
#   # min(rel_burd_increase)
#   # max(rel_burd_increase)
#   
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   #      create panel of maps showing all burden metrics
#   ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#   if(save_plots) png(paste0(sim_future_output_dir, '/_plots/map_withWithout', inter_name,'_burden_', pyr, '_', chw_cov, 'CHW.png'), res=600, width=6, height=3*length(burden_metric_names), units='in')
#   par(mar=c(0,1,2,0))
#   
#   num_colors = 40
#   colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
#   
#   # set layout for panel of maps
#   base_matrix = matrix(c(1,1,1,2,2,2,3,3, 1,1,1,2,2,2,3,3), nrow=2, byrow=TRUE)
#   # add rows for each burden metric
#   layout_matrix = base_matrix
#   for(cc in 2:length(burden_colnames_for_map)){
#     layout_matrix = rbind(layout_matrix, base_matrix + 3*(cc-1))
#   }
#   # add row for title
#   layout_matrix = layout_matrix + 1
#   layout_matrix = rbind(rep(1, ncol(layout_matrix)), layout_matrix)
#   layout(mat = layout_matrix)
#   
#   # title
#   plot.new()
#   text(0.5,0.5,"Malaria burden in each health district",cex=2.5,font=1)
#   # text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
#   #      line2user(line=4, side=3), "Malaria burden in each health district", xpd=NA, cex=2, font=2)
#   
#   # iterate through burden metrics, creating plots of each
#   for(cc in 1:length(burden_colnames_for_map)){
#     cur_colname = burden_colnames_for_map[cc]
#     vals_ordered_noPBO = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     vals_ordered_PBO = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
#     for (i_ds in 1:length(vals_ordered_noPBO$ds_ordered)){
#       cur_ds = vals_ordered_noPBO$ds_ordered[i_ds]
#       if(toupper(cur_ds) %in% toupper(noPBO_burden_df$admin_name)){
#         vals_ordered_noPBO$value[i_ds] = noPBO_burden_df[which(toupper(noPBO_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#         vals_ordered_PBO$value[i_ds] = PBO_burden_df[which(toupper(PBO_burden_df$admin_name) == toupper(cur_ds)), cur_colname]
#       }
#     }
#     min_value = min(c(vals_ordered_noPBO$value, vals_ordered_PBO$value), na.rm=TRUE)
#     max_value = max(c(vals_ordered_noPBO$value, vals_ordered_PBO$value), na.rm=TRUE)
#     # without intervention
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_noPBO$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - without ', inter_name))
#     # with intervention
#     col_cur = colorscale[sapply(floor((num_colors)*(vals_ordered_PBO$value - min_value) / (max_value - min_value))+1, min, num_colors)]
#     col_cur[is.na(col_cur)] = 'grey'
#     plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=paste0(burden_metric_names[cc], ' - with ', inter_name))
#     # legend
#     legend_label_vals = seq(min_value, max_value, length.out=5)
#     legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
#     plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = burden_metric_names[cc])
#     text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
#     rasterImage(legend_image, 0, 0, 1,1)
#     # fourth blank plot
#     # plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
#     
#   }
#   par(mfrow=c(1,1), mar=c(5,4,4,2))
#   if(save_plots) dev.off()
#
#
# }




####################################################################################
# barplots: percent reduction in burden relative to a reference experiment at a
#   specified reference year (e.g., to-present simulation at 2025)
####################################################################################

plot_relative_burden_barplots_vs_ref_year = function(sim_future_output_dir, reference_sim_output_dir, pop_filepath,
                                                     district_subset, cur_admins,
                                                     barplot_start_year, barplot_end_year,
                                                     reference_experiment_name, ref_year,
                                                     pyr='', chw_cov='',
                                                     scenario_names, experiment_names, scenario_palette,
                                                     LLIN2y_flag=FALSE, overwrite_files=FALSE,
                                                     separate_plots_flag=FALSE, standard_max_y=0.95,
                                                     show_error_bar=FALSE, burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)

  burden_metrics      = c('PfPR','PfPR','incidence','incidence','directMortality','directMortality','allMortality','allMortality')
  burden_colnames     = c('average_PfPR_U5','average_PfPR_all','incidence_U5','incidence_all',
                          'direct_death_rate_mean_U5','direct_death_rate_mean_all','all_death_rate_mean_U5','all_death_rate_mean_all')
  burden_metric_names = c('PfPR (U5)','PfPR (all ages)','incidence (U5)','incidence (all ages)',
                          'direct mortality (U5)','direct mortality (all ages)','mortality (U5)','mortality (all ages)')
  if(length(burden_metric_subset) >= 1){
    idx = which(burden_metrics %in% burden_metric_subset)
    burden_colnames     = burden_colnames[idx]
    burden_metric_names = burden_metric_names[idx]
  }

  relative_burden_all_df = data.frame()
  for(ss in 1:length(scenario_names)){
    relative_burden_df = get_relative_burden_vs_ref_year(
      reference_sim_output_filepath=reference_sim_output_dir,
      reference_experiment_name=reference_experiment_name,
      ref_start_year=ref_year, ref_end_year=ref_year,
      comparison_sim_output_filepath=sim_future_output_dir,
      comparison_experiment_name=experiment_names[ss],
      comparison_scenario_name=scenario_names[ss],
      start_year=barplot_start_year, end_year=barplot_end_year,
      admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins,
      LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
    relative_burden_df = relative_burden_df[, which(colnames(relative_burden_df) %in% c('scenario','Run_Number',burden_colnames))]
    if(nrow(relative_burden_all_df) == 0) relative_burden_all_df = relative_burden_df else relative_burden_all_df = rbind(relative_burden_all_df, relative_burden_df)
  }

  relative_burden_all_df$scenario = factor(relative_burden_all_df$scenario, levels=scenario_names)

  standard_min_y = 0
  cur_min = min(relative_burden_all_df[, which(colnames(relative_burden_all_df) %in% burden_colnames)], na.rm=TRUE)
  cur_max = max(relative_burden_all_df[, which(colnames(relative_burden_all_df) %in% burden_colnames)], na.rm=TRUE)
  if(cur_min < standard_min_y) standard_min_y = cur_min
  if(cur_max > standard_max_y) standard_max_y = cur_max

  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name  = burden_metric_names[bb]
    rel_burden_agg = as.data.frame(relative_burden_all_df) %>%
      dplyr::select(match(c(current_burden_name,'scenario'), names(.))) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(mean_rel = mean(get(current_burden_name), na.rm=TRUE),
                       max_rel  = max(get(current_burden_name),  na.rm=TRUE),
                       min_rel  = min(get(current_burden_name),  na.rm=TRUE))
    gg_list[[bb]] = ggplot(rel_burden_agg) +
      geom_bar(aes(x=scenario, y=mean_rel, fill=scenario), stat='identity') +
      scale_y_continuous(labels=percent_format(), limits=c(standard_min_y, standard_max_y)) +
      ylab(paste0('% change vs ', ref_year)) +
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values=scenario_palette) +
      theme_gridlines_no_box() +
      theme(legend.position='top', legend.box='horizontal', legend.title=element_blank(),
            text=element_text(size=text_size), legend.text=element_text(size=text_size),
            axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg_list[[bb]] = gg_list[[bb]] +
        geom_errorbar(aes(x=scenario, ymin=min_rel, ymax=max_rel), width=0.4, colour='black', alpha=0.9, size=1)
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] +
        theme(legend.position='none', plot.title=element_blank(), text=element_text(size=separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/barplot_percent_reduction_vs', ref_year, '_', burden_metric_name, '_', district_subset, '.png'),
             separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }

  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position='none') + theme(text=element_text(size=text_size))
  }
  gg = grid.arrange(grobs=gg_list, layout_matrix=rbind(matrix(rep(1, ceiling(length(burden_colnames)/2)), nrow=1),
                                                       matrix(2:(length(burden_colnames)+1), nrow=2, byrow=FALSE)))
  return(gg)
}


####################################################################################
# barplots: absolute burden averted relative to a reference experiment at a
#   specified reference year (e.g., to-present simulation at 2025)
#   y-axis is metric-specific (no shared scale) since units differ across metrics
####################################################################################

plot_difference_burden_barplots_vs_ref_year = function(sim_future_output_dir, reference_sim_output_dir, pop_filepath,
                                                       district_subset, cur_admins,
                                                       barplot_start_year, barplot_end_year,
                                                       reference_experiment_name, ref_year,
                                                       pyr='', chw_cov='',
                                                       scenario_names, experiment_names, scenario_palette,
                                                       LLIN2y_flag=FALSE, overwrite_files=FALSE,
                                                       separate_plots_flag=FALSE,
                                                       show_error_bar=FALSE, burden_metric_subset=c()){
  admin_pop = read.csv(pop_filepath)

  burden_metrics      = c('PfPR','PfPR','incidence','incidence','directMortality','directMortality','allMortality','allMortality')
  burden_colnames     = c('average_PfPR_U5','average_PfPR_all','incidence_U5','incidence_all',
                          'direct_death_rate_mean_U5','direct_death_rate_mean_all','all_death_rate_mean_U5','all_death_rate_mean_all')
  burden_metric_names = c('PfPR (U5)','PfPR (all ages)','incidence (U5)','incidence (all ages)',
                          'direct mortality (U5)','direct mortality (all ages)','mortality (U5)','mortality (all ages)')
  burden_metric_units = c('prevalence (%)','prevalence (%)','cases per 1000/year','cases per 1000/year',
                          'deaths per 1000/year','deaths per 1000/year','deaths per 1000/year','deaths per 1000/year')
  if(length(burden_metric_subset) >= 1){
    idx = which(burden_metrics %in% burden_metric_subset)
    burden_colnames     = burden_colnames[idx]
    burden_metric_names = burden_metric_names[idx]
    burden_metric_units = burden_metric_units[idx]
  }

  difference_burden_all_df = data.frame()
  for(ss in 1:length(scenario_names)){
    difference_burden_df = get_difference_burden_vs_ref_year(
      reference_sim_output_filepath=reference_sim_output_dir,
      reference_experiment_name=reference_experiment_name,
      ref_start_year=ref_year, ref_end_year=ref_year,
      comparison_sim_output_filepath=sim_future_output_dir,
      comparison_experiment_name=experiment_names[ss],
      comparison_scenario_name=scenario_names[ss],
      start_year=barplot_start_year, end_year=barplot_end_year,
      admin_pop=admin_pop, district_subset=district_subset, cur_admins=cur_admins,
      LLIN2y_flag=LLIN2y_flag, overwrite_files=overwrite_files)
    difference_burden_df = difference_burden_df[, which(colnames(difference_burden_df) %in% c('scenario','Run_Number',burden_colnames))]
    if(nrow(difference_burden_all_df) == 0) difference_burden_all_df = difference_burden_df else difference_burden_all_df = rbind(difference_burden_all_df, difference_burden_df)
  }

  difference_burden_all_df$scenario = factor(difference_burden_all_df$scenario, levels=scenario_names)

  gg_list = list()
  for(bb in 1:length(burden_colnames)){
    current_burden_name = burden_colnames[bb]
    burden_metric_name  = burden_metric_names[bb]
    burden_metric_unit  = burden_metric_units[bb]
    diff_burden_agg = as.data.frame(difference_burden_all_df) %>%
      dplyr::select(match(c(current_burden_name,'scenario'), names(.))) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(mean_diff = mean(get(current_burden_name), na.rm=TRUE),
                       max_diff  = max(get(current_burden_name),  na.rm=TRUE),
                       min_diff  = min(get(current_burden_name),  na.rm=TRUE))
    gg_list[[bb]] = ggplot(diff_burden_agg) +
      geom_bar(aes(x=scenario, y=mean_diff, fill=scenario), stat='identity') +
      scale_y_continuous(labels = if(grepl('PfPR', current_burden_name)) percent_format(accuracy = 0.1) else comma_format()) +
      ylab(paste0('Averted vs ', ref_year, '\n(', burden_metric_unit, ')')) +
      geom_hline(yintercept=0, color='black') +
      ggtitle(gsub('\\(births\\)', '', burden_metric_name)) +
      scale_fill_manual(values=scenario_palette) +
      theme_gridlines_no_box() +
      theme(legend.position='top', legend.box='horizontal', legend.title=element_blank(),
            text=element_text(size=text_size), legend.text=element_text(size=text_size),
            axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), axis.line.x=element_blank(),
            plot.margin=unit(c(0,1,1,0), 'cm'))
    if(show_error_bar){
      gg_list[[bb]] = gg_list[[bb]] +
        geom_errorbar(aes(x=scenario, ymin=min_diff, ymax=max_diff), width=0.4, colour='black', alpha=0.9, size=1)
    }
    if(separate_plots_flag){
      separate_plot = gg_list[[bb]] +
        theme(legend.position='none', plot.title=element_blank(), text=element_text(size=separate_plot_text_size))
      ggsave(paste0(sim_future_output_dir, '/_plots/barplot_burden_averted_vs', ref_year, '_', burden_metric_name, '_', district_subset, '.png'),
             separate_plot, dpi=600, width=4, height=3, units='in')
    }
  }

  gg_list = append(list(ggpubr::as_ggplot(ggpubr::get_legend(gg_list[[1]]))), gg_list)
  for(bb in 2:(length(burden_colnames)+1)){
    gg_list[[bb]] = gg_list[[bb]] + theme(legend.position='none') + theme(text=element_text(size=text_size))
  }
  gg = grid.arrange(grobs=gg_list, layout_matrix=rbind(matrix(rep(1, ceiling(length(burden_colnames)/2)), nrow=1),
                                                       matrix(2:(length(burden_colnames)+1), nrow=2, byrow=FALSE)))
  return(gg)
}




