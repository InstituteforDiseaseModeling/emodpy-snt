# plot_map_intervention_coverage_input.R

# library(rgdal)
library(raster)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(tidyverse)
library(sf)
library(reshape2)
library(pals)
library(prettyGraphs)


##########################################################
# setup
##########################################################
user = Sys.getenv("USERNAME")
user_path = file.path("C:/Users",user)

country = 'NGA'  #'SLE'  # 'BDI'
dta_dir = 'C:/Users/moniqueam/IDM Dropbox/Malaria Team Folder/data'
script_dir = 'C:/Users/moniqueam/Documents/emodpy-snt/r_utilities'


if(country =='BDI'){
  hbhi_dir = 'C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder/projects/burundi_hbhi/snt_2023'
  admin_shapefile_filepath = (paste0(hbhi_dir, '/SpatialClustering/reference_rasters_shapefiles/bdi_adm2.shp'))
  data_dir = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/data/Burundi'
} else if(country=='NGA'){
  hbhi_dir = 'C:/Users/moniqueam/IDM Dropbox/Malaria Team Folder/projects/snt/Nigeria/snt_2023'
  admin_shapefile_filepath = paste0(hbhi_dir, '/SpatialClustering/reference_rasters_shapefiles/NGA_DS_clusteringProjection.shp')
  
  # base_filepath = paste0(user_path, '/IDM Dropbox/NU_collaboration')
  # hbhi_dir = paste0(user_path, '/IDM Dropbox/Monique Ambrose/NU_collaboration/hbhi_nigeria/snt_2022')
  # admin_shapefile_filepath = (paste0(base_filepath, '/hbhi_nigeria/SpatialClustering/reference_rasters_shapefiles/NGA_DS_clusteringProjection.shp'))
}


sim_output_dir = paste0(hbhi_dir, '/simulation_output')
pop_filepath = paste0(hbhi_dir, '/admin_pop_archetype.csv')
admin_shapefile = st_read(admin_shapefile_filepath)

base_sim_input_dir = paste0(hbhi_dir, '/simulation_inputs')
intervention_coordinator = read.csv(paste0(base_sim_input_dir, '/_intervention_file_references/Interventions_to_present.csv'))
# intervention_coordinator = read.csv(paste0(base_sim_input_dir, '/_intervention_file_references/Interventions_for_projections.csv'))
scenario_row = 1  #5
if(scenario_row>nrow(intervention_coordinator)) scenario_row = nrow(intervention_coordinator)

source(paste0(script_dir,'/standardize_admin_names.R'))

# get dataframe with standardized admin names, states, archetype, and population information
admin_pop = read.csv(pop_filepath)
# standardize shapefile names
admin_shapefile$NOMDEP = standardize_admin_names_in_vector(target_names=admin_pop$admin_name, origin_names=admin_shapefile$NOMDEP)

if(!dir.exists(file.path(base_sim_input_dir, 'plots'))) dir.create(file.path(base_sim_input_dir, 'plots'))
if(!dir.exists(file.path(base_sim_input_dir, 'plots', 'interventions_2010_toPresent'))) dir.create(file.path(base_sim_input_dir, 'plots', 'interventions_2010_toPresent'))
if(!dir.exists(file.path(base_sim_input_dir, 'plots', 'interventions_projections'))) dir.create(file.path(base_sim_input_dir, 'plots', 'interventions_projections'))



##########################################################
# functions for ggplot version
##########################################################

# function to combine multiple plots for same intervention that share a legend
grid_arrange_shared_legend_plotlist =function(..., 
                                              plotlist=NULL,
                                              ncol = length(list(...)),
                                              nrow = NULL,
                                              position = c("bottom", "right")) {
  
  plots <- c(list(...), plotlist)
  
  if (is.null(nrow)) nrow = ceiling(length(plots)/ncol)
  
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position) + guides(fill = guide_legend(reverse=T)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}

change_legend_size <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

# function to create and save maps of intervention coverages
create_coverage_input_maps = function(inter_input, inter_years, output_filename, colorscale, min_value, max_value, num_colors){
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed==1,]
  plot_list = list()
  for(yy in 1:length(inter_years)){
    inter_input_cur = inter_input[inter_input$year == inter_years[yy],]
    inter_input_cur$output_value = inter_input_cur[[coverage_colname]]
    admin_cur = admin_shapefile %>%
      left_join(inter_input_cur, by=c('NOMDEP' = 'admin_name')) %>%
      mutate(binned_values=cut(output_value,
                               breaks=round(seq((min_value), (max_value), length.out = (num_colors+1)),2)))
    plot_list[[yy]] = ggplot(admin_cur) +
      geom_sf(aes(fill=binned_values), size=0.1) +
      scale_fill_manual(values=setNames(colorscale, levels(admin_cur$binned_values)), drop=FALSE, name='coverage', na.value='grey96') + 
      ggtitle(inter_years[yy]) +
      guides(fill = guide_legend(reverse=T)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5)) 
    plot_list[[yy]] = change_legend_size(plot_list[[yy]], pointSize=10, textSize=10, spaceLegend=1)
  }
  
  gg = grid_arrange_shared_legend_plotlist(plotlist=plot_list, ncol=length(plot_list), position='right')
  ggsave(output_filename, gg, width = (length(inter_years)+1)*1.8, height=2.3, units='in', dpi=800)
}

##########################################################
# create ggplot2 maps
##########################################################
# num_colors=9
# min_value = 0#min(inter_input[[coverage_colname]], na.rm=TRUE)
# max_value = 0.9#max(inter_input[[coverage_colname]], na.rm=TRUE)
# colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
# plot(1:num_colors, col=colorscale, pch=20, cex=3)
# ### updated colors for consistency with plotting DHS outputs
# # colors_range_0_to_1 = add.alpha(pals::parula(101), alpha=0.5)  # specify colorscale. To go from a value vv to a color, take colors_range_0_to_1[1+round(vv*100)]
# # plot(1:101, col=add.alpha(pals::parula(101), alpha=0.5), pch=20 )
# # colorscale = add.alpha(pals::parula(num_colors), alpha=0.5)  # updated alpha from DHS plot to go from more to less transparent when plotting simulated coverage (easier to interpret empty spaces)
num_colors=10
# colorscale = unname(rev(add.alpha(pals::parula(num_colors), alpha=0.3)))  # updated alpha from DHS plot to go from more to less transparent when plotting simulated coverage (easier to interpret empty spaces)
max_value = 1
min_value=0
# colors_range_0_to_1 = rev(add.alpha(pals::parula(101), alpha=0.5))  # specify colorscale. To go from a value vv to a color, take colors_range_0_to_1[1+round(vv*100)]
# num_colors=10
# min_value = 0
# max_value = 1
# color_vals = rev(pals::parula(num_colors))
alpha_vals=seq(0.3,1,length.out=num_colors)
color_vals = unname(rev(pals::parula(num_colors)))  # updated alpha from DHS plot to go from more to less transparent when plotting simulated coverage (easier to interpret empty spaces)
colorscale = color_vals
for(ii in 1:num_colors){
  colorscale[ii] = add.alpha(color_vals[ii], alpha=alpha_vals[ii])
}
# plot(1:num_colors, col=colorscale, pch=20, cex=3)


###### plot itn mass distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$ITN_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$ITN_filename[scenario_row], '.png')
coverage_colname = 'itn_u5'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = sort(unique(inter_input$year))
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}

# plot for the pyrethroid-only DS
intervention_filename = paste0(base_sim_input_dir, '/interventions_projections/llin_mass_pri1.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/interventions_projections/llin_mass_pri1_adminSubset.png')
coverage_colname = 'itn_u5'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  inter_input = inter_input[inter_input$llin_type=='standard',]
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = 2022
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}

###### plot itn ANC distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$ANC_ITN_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$ANC_ITN_filename[scenario_row], '.png')
coverage_colname = 'coverage'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = sort(unique(inter_input$year))
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}


###### plot itn EPI distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$EPI_ITN_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$EPI_ITN_filename[scenario_row], '.png')
coverage_colname = 'coverage'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  inter_input = inter_input[inter_input$birthday_age == 1,]
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = sort(unique(inter_input$year))
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}


###### plot itn community distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$CHW_ITN_annual_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$CHW_ITN_annual_filename[scenario_row], '.png')
coverage_colname = 'itn_u5'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = sort(unique(inter_input$year))
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}



###### plot IRS distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$IRS_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$IRS_filename[scenario_row], '.png')
coverage_colname = 'effective_coverage'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = unique(inter_input$year)
  inter_years = sort(unique(inter_input$year))
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}


###### plot CM distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$CM_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$CM_filename[scenario_row], '.png')
coverage_colname = 'U5_coverage'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = sort(unique(inter_input$year))
  inter_years = inter_years[seq(2,length(inter_years), by=2)]
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}

# ###### plot CM in GF states only #######
# # scenario file describing which LGAs are in each plan and standardize names
# scenario_nmcp_df = read.csv(paste0(user_path, '/Dropbox (IDM)/NU_collaboration/nigeria_who/NGA_2022_SNT/future_projection_scenarios_20230225.csv'))
# scenario_nmcp_df$admin_name = standardize_admin_names_in_vector(target_names=admin_pop$LGA, origin_names=scenario_nmcp_df$adm2)
# # get subset of admin to plot
# cur_admins_gf = unique(scenario_nmcp_df$admin_name[grepl('Global Fund', scenario_nmcp_df$funder)])
# intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$CM_filename[scenario_row],'.csv')
# output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$CM_filename[scenario_row], '_subsetAdminsGF.png')
# coverage_colname = 'U5_coverage'
# if(file.exists(intervention_filename)){
#   inter_input = read.csv(intervention_filename)
#   if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
#   inter_input = inter_input[inter_input$year == max(inter_input$year),]
#   inter_input = inter_input[inter_input$admin_name %in% cur_admins_gf,]
#   inter_years = sort(unique(inter_input$year))
#   create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
# }





###### plot SMC distribution coverage #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$SMC_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$SMC_filename[scenario_row], '.png')
coverage_colname = 'coverage_U5_total'
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  if('round' %in% colnames(inter_input)) inter_input = inter_input[inter_input$round == 1,]
  inter_years = sort(unique(inter_input$year))
  inter_input$coverage_U5_total = inter_input$high_access_U5 * inter_input$coverage_high_access_U5 + (1-inter_input$high_access_U5) * inter_input$coverage_low_access_U5
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}


###### plot RTSS vaccine distribution coverage #######
if('vacc_filename' %in% colnames(intervention_coordinator)){
  intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$vacc_filename[scenario_row],'.csv')
  output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$vacc_filename[scenario_row], '.png')
  coverage_colname = 'coverage'
  if(file.exists(intervention_filename)){
    inter_input = read.csv(intervention_filename)
    inter_input = inter_input[inter_input$vaccine=='primary',]
    if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
    inter_input$year = 2024
    inter_years = sort(unique(inter_input$year))
    create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
  }
}



###### plot IPTp #######
intervention_filename = paste0(hbhi_dir,'/simulation_inputs/IPTp/estimated_past_IPTp_each_DS.csv')
if(file.exists(intervention_filename)){
  inter_input_raw = read.csv(intervention_filename)[,-1]
  inter_input = reshape2::melt(inter_input_raw, id.vars=c("admin_name"))
  inter_input$year = as.numeric(gsub('X','',inter_input$variable))
  coverage_colname = 'value'
  output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/sim_input_map_IPTp.png')
  inter_years = sort(unique(inter_input$year))
  create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
}


# ###### plot IPTp in GF states only #######
# # scenario file describing which LGAs are in each plan and standardize names
# scenario_nmcp_df = read.csv(paste0(user_path, '/Dropbox (IDM)/NU_collaboration/nigeria_who/NGA_2022_SNT/future_projection_scenarios_20230225.csv'))
# scenario_nmcp_df$admin_name = standardize_admin_names_in_vector(target_names=admin_pop$LGA, origin_names=scenario_nmcp_df$adm2)
# # get subset of admin to plot
# cur_admins_gf = unique(scenario_nmcp_df$admin_name[grepl('Global Fund', scenario_nmcp_df$funder)])
# intervention_filename = paste0(hbhi_dir,'/simulation_inputs/IPTp/estimated_past_IPTp_each_DS.csv')
# if(file.exists(intervention_filename)){
#   inter_input_raw = read.csv(intervention_filename)[,-1]
#   inter_input = reshape2::melt(inter_input_raw, id.vars=c("admin_name"))
#   inter_input$year = as.numeric(gsub('X','',inter_input$variable))
#   inter_input = inter_input[inter_input$year == max(inter_input$year),]
#   inter_input = inter_input[inter_input$admin_name %in% cur_admins_gf,]
#   coverage_colname = 'value'
#   output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/sim_input_map_IPTp_subsetAdminsGF.png')
#   inter_years = sort(unique(inter_input$year))
#   create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
# }



###### plot IPTi #######
if('PMC_filename' %in% colnames(intervention_coordinator)){
  intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$PMC_filename[scenario_row],'.csv')
  output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$PMC_filename[scenario_row], '.png')
  if(file.exists(intervention_filename)){
    inter_input = read.csv(intervention_filename)
    # subset to first touchpoint, based on vacc_dpt3
    inter_input = inter_input[inter_input$vaccine == 'vacc_dpt3',]
    inter_input$year = 2024
    coverage_colname = 'coverage'
    inter_years = unique(inter_input$year)
    create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
  }
}




####################################################
# categorical inputs
####################################################

###### plot net types for itn mass distribution #######
intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$ITN_filename[scenario_row],'.csv')
output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$ITN_filename[scenario_row], '_netType.png')
coverage_colname = 'lin_type'
start_day = (2024-2022)*365
if(file.exists(intervention_filename)){
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_input = inter_input[inter_input$simday > start_day,]
  inter_input = inter_input[,c('admin_name', 'llin_type')]
  admin_cur = admin_shapefile %>%
    left_join(inter_input, by=c('NOMDEP' = 'admin_name'))
  net_palette = c('#D1AB8F', '#9CE4B2','#C78FFF')
  names(net_palette) = c('standard', 'PBO', 'IG2')
  gg = ggplot(admin_cur) +
    geom_sf(aes(fill=llin_type)) +
    scale_fill_manual(values=net_palette, name='Net type', drop=FALSE, na.value='white') +
    theme_void() 
  ggsave(filename=output_filename, plot=gg, width=2*1.8, height=2.3, units='in', dpi=800)
}


# ###### plot net types for itn mass distribution - subset to LGAs with PBO in prioritized plan #######
# # scenario file describing which LGAs are in each plan and standardize names
# scenario_nmcp_df = read.csv(paste0(user_path, '/Dropbox (IDM)/NU_collaboration/nigeria_who/NGA_2022_SNT/future_projection_scenarios_20230225.csv'))
# scenario_nmcp_df$admin_name = standardize_admin_names_in_vector(target_names=admin_pop$LGA, origin_names=scenario_nmcp_df$adm2)
# # get subset of admin to plot
# cur_admins_pbo = unique(scenario_nmcp_df$admin_name[grepl('PBO', scenario_nmcp_df$mass_llins_fund)])
# # read in and plot file for this scenario_row
# intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$ITN_filename[scenario_row],'.csv')
# output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$ITN_filename[scenario_row], '_netType_subsetAdmins.png')
# coverage_colname = 'lin_type'
# start_day = (2024-2022)*365
# if(file.exists(intervention_filename)){
#   inter_input = read.csv(intervention_filename)
#   if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
#   inter_input = inter_input[inter_input$simday > start_day,]
#   inter_input = inter_input[inter_input$admin_name %in% cur_admins_pbo,]
#   inter_input = inter_input[,c('admin_name', 'State', 'llin_type')]
#   admin_cur = admin_shapefile %>%
#     left_join(inter_input, by=c('NOMDEP' = 'admin_name'))
#   net_palette = c('#D1AB8F', '#9CE4B2','#C78FFF')
#   names(net_palette) = c('standard', 'PBO', 'IG2')
#   gg = ggplot(admin_cur) +
#     geom_sf(aes(fill=llin_type)) +
#     scale_fill_manual(values=net_palette, name='Net type', drop=FALSE, na.value='white') +
#     theme_void() 
#   ggsave(filename=output_filename, plot=gg, width=2*1.8, height=2.3, units='in', dpi=800)
# }





if(FALSE){
  
  # permethrin resistance colorscale
  num_colors=10
  max_value = 1
  min_value=0
  color_vals_IR = unname((pals::coolwarm(num_colors)))  # updated alpha from DHS plot to go from more to less transparent when plotting simulated coverage (easier to interpret empty spaces)
  colorscale_IR = color_vals_IR

  ###### plot permethrin resistance #######
  intervention_filename = paste0(hbhi_dir,'/simulation_inputs/intermediate_files/insecticide_resistance/permethrin_mortality_admin_estimates.csv')
  intervention_filename = paste0(hbhi_dir,'/ento/insecticide_resistance_DS_means/insecticide_Permethrin_mortality_DS_means.csv')
  if(file.exists(intervention_filename)){
    inter_input_raw = read.csv(intervention_filename)[,-1]
    if('DS' %in% colnames(inter_input_raw)) colnames(inter_input_raw)[which(colnames(inter_input_raw)=='DS')] = 'admin_name'
    inter_input = melt(inter_input_raw, id.vars=c("admin_name"))
    inter_input$year = as.numeric(gsub('X','',inter_input$variable))
    inter_input$resistance = 1-inter_input$value
    coverage_colname = 'resistance'
    # coverage_colname = 'value'
    output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/sim_input_map_perm_resistance.png') # mortality
    inter_years = seq(2011,2017)
    create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale_IR, min_value=min_value, max_value=max_value, num_colors=num_colors)
  }
  
  # alternative assumption about resistance
  alt_mort=0.5
  intervention_filename = paste0(hbhi_dir,'/ento/insecticide_resistance_DS_means/insecticide_Permethrin_mortality_DS_means.csv')
  if(file.exists(intervention_filename)){
    inter_input_raw = read.csv(intervention_filename)[,-1]
    if('DS' %in% colnames(inter_input_raw)) colnames(inter_input_raw)[which(colnames(inter_input_raw)=='DS')] = 'admin_name'
    inter_input = melt(inter_input_raw, id.vars=c("admin_name"))
    inter_input$year = 2023
    inter_input$resistance = 1-alt_mort
    coverage_colname = 'resistance'
    # coverage_colname = 'value'
    output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/sim_input_map_perm_resistance_mort',round(100*alt_mort),'.png') # mortality
    inter_years = 2023
    create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale_IR, min_value=min_value, max_value=max_value, num_colors=num_colors)
  }
  
  
  ###### plot permethrin resistance for GF LGAs #######
  # scenario file describing which LGAs are in each plan and standardize names
  scenario_nmcp_df = read.csv(paste0(user_path, '/Dropbox (IDM)/NU_collaboration/nigeria_who/NGA_2022_SNT/future_projection_scenarios_20230225.csv'))
  scenario_nmcp_df$admin_name = standardize_admin_names_in_vector(target_names=admin_pop$LGA, origin_names=scenario_nmcp_df$adm2)
  # get subset of admin to plot
  cur_admins_gf = unique(scenario_nmcp_df$admin_name[grepl('Global Fund', scenario_nmcp_df$funder)])
  intervention_filename = paste0(hbhi_dir,'/simulation_inputs/intermediate_files/insecticide_resistance/permethrin_mortality_admin_estimates.csv')
  if(file.exists(intervention_filename)){
    inter_input_raw = read.csv(intervention_filename)[,-1]
    if('DS' %in% colnames(inter_input_raw)) colnames(inter_input_raw)[which(colnames(inter_input_raw)=='DS')] = 'admin_name'
    inter_input = melt(inter_input_raw, id.vars=c("admin_name"))
    inter_input$year = as.numeric(gsub('X','',inter_input$variable))
    inter_input = inter_input[inter_input$admin_name %in% cur_admins_gf,]
    inter_input$survival = 1-inter_input$value
    coverage_colname = 'survival'
    # coverage_colname = 'value'
    output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/sim_input_map_perm_bioassay_survival_subsetAdminGF.png')
    inter_years = seq(2014,2018, by=2)
    create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale_IR, min_value=min_value, max_value=max_value, num_colors=num_colors)
  }
  
  
  ######  plot EMOD LLIN kill parameter for each mass distribution  ###### 
  intervention_filename = paste0(base_sim_input_dir, '/', intervention_coordinator$ITN_filename[scenario_row],'.csv')
  output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/',intervention_coordinator$ITN_filename[scenario_row], '_KillInitial.png')
  coverage_colname = 'kill_initial'
  if(file.exists(intervention_filename)){
    inter_input = read.csv(intervention_filename)
    if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
    inter_years = sort(unique(inter_input$year))
    create_coverage_input_maps(inter_input=inter_input, inter_years=inter_years, output_filename=output_filename, colorscale=colorscale, min_value=min_value, max_value=max_value, num_colors=num_colors)
  }
  
  
  
  ###### plot categorical packages #######
  GF_package_filename = paste0(data_dir,'/WHO/scenario_plans/BDI_GF.csv')
  if(file.exists(GF_package_filename)){
    packages_df = read.csv(GF_package_filename)
    coverage_colname = 'plan_pri'
    output_filename = paste0(hbhi_dir, '/simulation_inputs/plots/sim_input_map_GF_packages.png')
    # update names without iCCM or PECADOM
    packages_df[[coverage_colname]][packages_df[[coverage_colname]] == 'Mass&Rout-Pyr'] = 'contCM+Mass&Rout-Pyr'
    
    
    # PBO version
    packages_df_PBO = packages_df
    packages_df_PBO[[coverage_colname]] = gsub('Pyr', 'PBO', packages_df_PBO[[coverage_colname]])
    
    # BAU scenario
    irs_districts = unique(read.csv(paste0(hbhi_dir, '/simulation_inputs/interventions_projections/irs_bau.csv'))$admin_name)
    packages_df_BAU = packages_df
    packages_df_BAU[[coverage_colname]] = 'contCM+Mass&Rout-Pyr'
    packages_df_BAU[[coverage_colname]][packages_df_BAU$adm2 %in% irs_districts] = 'contCM+Mass&Rout-Pyr+IRS'
    
    
    package_names_temp = sort(unique(as.character(c(packages_df[[coverage_colname]], packages_df_PBO[[coverage_colname]], packages_df_BAU[[coverage_colname]]))))
    package_names = package_names_temp
    package_names[1] = package_names_temp[2]
    package_names[2] = package_names_temp[1]
    packages_df[[coverage_colname]] = factor(packages_df[[coverage_colname]], levels=package_names)
    packages_df_PBO[[coverage_colname]] = factor(packages_df_PBO[[coverage_colname]], levels=package_names)
    packages_df_BAU[[coverage_colname]] = factor(packages_df_BAU[[coverage_colname]], levels=package_names)
    
    
    # attach packages to map
    packages_df$output_value = packages_df[[coverage_colname]]
    admin_cur = admin_shapefile %>%
      left_join(packages_df, by=c('NOMDEP' = 'adm2'))
    packages_df_PBO$output_value = packages_df_PBO[[coverage_colname]]
    admin_cur_PBO = admin_shapefile %>%
      left_join(packages_df_PBO, by=c('NOMDEP' = 'adm2'))
    packages_df_BAU$output_value = packages_df_BAU[[coverage_colname]]
    admin_cur_BAU = admin_shapefile %>% 
      left_join(packages_df_BAU, by=c('NOMDEP' = 'adm2'))
    
    plot_list = list()
    plot_list[[1]] = ggplot(admin_cur_BAU) +
      geom_sf(aes(fill=output_value)) +
      scale_fill_brewer(palette = 'Set3', name='Interventions', drop=FALSE) +
      theme_void() 
    plot_list[[2]] = ggplot(admin_cur) +
      geom_sf(aes(fill=output_value)) +
      scale_fill_brewer(palette = 'Set3', name='Interventions', drop=FALSE) +
      theme_void() 
    plot_list[[3]] = ggplot(admin_cur_PBO) +
      geom_sf(aes(fill=output_value)) +
      scale_fill_brewer(palette = 'Set3', name='Interventions', drop=FALSE) +
      theme_void() 
    
    
    gg = grid_arrange_shared_legend_plotlist(plotlist=plot_list, ncol=length(plot_list), position='right')
    
    ggsave(output_filename, gg, width=(length(plot_list)+1)*2.3, height=3, units='in', dpi=800)
    
  }
  
  
  
  
  ##########################################################
  ###### plot population sizes in each admin #######
  pop_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
  for (i_ds in 1:length(pop_ordered$ds_ordered)){
    cur_ds = pop_ordered$ds_ordered[i_ds]
    pop_ordered$value[i_ds] = admin_pop$pop_size[which(toupper(admin_pop$admin_name) == toupper(cur_ds))]
  }
  num_colors=8
  min_value = min(c(pop_ordered$value), na.rm=TRUE)
  max_value = max(c(pop_ordered$value), na.rm=TRUE)
  colorscale = colorRampPalette(brewer.pal(9, 'BuGn'))(num_colors)
  
  
  admin_cur = admin_shapefile %>%
    left_join(admin_pop, by=c('NOMDEP' = 'admin_name')) %>%
    mutate(binned_values=cut(pop_size,
                             breaks=round(seq(floor(min_value), ceiling(max_value), length.out = num_colors)),
                             include.lowest=TRUE,
                             dig.lab=50))
  ggplot(admin_cur) +
    geom_sf(aes(fill=binned_values)) +
    scale_fill_manual(values=setNames(colorscale, levels(admin_cur$binned_values))) + 
    theme_map()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ################################################################################
  # plots from base R instead of ggplot2
  ################################################################################
  
  hbhi_dir = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/burundi_hbhi'
  sim_output_dir = paste0(hbhi_dir, '/simulation_output')
  pop_filepath = paste0(hbhi_dir, '/admin_pop_archetype.csv')
  admin_shapefile_filepath = paste0(hbhi_dir, '/SpatialClustering/reference_rasters_shapefiles/bdi_adm2.shp')
  shapefile_admin_colname = 'NOMDEP'
  admin_shapefile = shapefile(admin_shapefile_filepath)
  admin_shapefile = st_read(admin_shapefile_filepath)
  
  
  
  ###### plot population sizes in each admin #######
  admin_pop = read.csv(pop_filepath)
  pop_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
  for (i_ds in 1:length(pop_ordered$ds_ordered)){
    cur_ds = pop_ordered$ds_ordered[i_ds]
    pop_ordered$value[i_ds] = admin_pop$pop_size[which(toupper(admin_pop$admin_name) == toupper(cur_ds))]
  }
  num_colors=8
  min_value = min(c(pop_ordered$value), na.rm=TRUE)
  max_value = max(c(pop_ordered$value), na.rm=TRUE)
  colorscale = colorRampPalette(brewer.pal(9, 'BuGn'))(num_colors)
  col_cur = colorscale[sapply(floor((num_colors)*(pop_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
  col_cur[is.na(col_cur)] = 'grey'
  plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3))
  
  # legend
  legend_label_vals = seq(min_value/1000, max_value/1000, length.out=5)
  legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = c('population size', '(thousands)'))
  text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals))
  rasterImage(legend_image, 0, 0, 1,1)
  
  admin_cur = admin_shapefile %>%
    left_join(admin_pop, by=c('NOMDEP' = 'admin_name')) %>%
    mutate(binned_values=cut(pop_size,
                             breaks=round(seq(floor(min_value), ceiling(max_value), length.out = num_colors)),
                             include.lowest=TRUE,
                             dig.lab=50))
  ggplot(admin_cur) +
    geom_sf(aes(fill=binned_values)) +
    scale_fill_manual(values=setNames(colorscale, levels(admin_cur$binned_values))) + 
    theme_map()
  
  
  
  
  
  num_colors=8
  min_value = 0.2#min(inter_input[[coverage_colname]], na.rm=TRUE)
  max_value = 1#max(inter_input[[coverage_colname]], na.rm=TRUE)
  colorscale = colorRampPalette(brewer.pal(9, 'YlGnBu'))(num_colors)
  
  
  # legend
  legend_label_vals = seq(min_value, max_value, length.out=5)
  legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = c('coverage'))
  text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
  par(mfrow=c(1,1), mar=c(5,4,4,2))
  
  
  
  ###### plot itn mass distribution coverage #######
  intervention_filename = paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/itn_mass_coverages_mort_2010_toPresent.csv')
  coverage_colname = 'itn_u5'
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = unique(inter_input$year)
  
  
  png(paste0(hbhi_dir, '/simulation_inputs/plots/itn_mass_sim_input.png'), width = (length(inter_years)+1)*1.8, height=2.3, units='in', res=800)
  par(mfrow=c(1,(length(inter_years)+1)), mar=c(0,1,2,0))
  
  for(yy in 1:length(inter_years)){
    inter_input_cur = inter_input[inter_input$year == inter_years[yy],]
    
    inter_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
    for (i_ds in 1:length(inter_ordered$ds_ordered)){
      cur_ds = inter_ordered$ds_ordered[i_ds]
      inter_ordered$value[i_ds] = inter_input_cur[[coverage_colname]][which(toupper(inter_input_cur$admin_name) == toupper(cur_ds))]
    }
    col_cur = colorscale[sapply(floor((num_colors)*(inter_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
    col_cur[is.na(col_cur)] = 'grey'
    plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=inter_years[yy])
  }
  
  # legend
  legend_label_vals = seq(min_value, max_value, length.out=5)
  legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = c('coverage'))
  text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
  par(mfrow=c(1,1), mar=c(5,4,4,2))
  
  
  
  
  
  ###### plot itn ANC distribution coverage #######
  intervention_filename = paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/anc_itn_use_coverages_mort_2010_toPresent.csv')
  coverage_colname = 'coverage'
  inter_input = read.csv(intervention_filename)
  if('seed' %in% colnames(inter_input)) inter_input = inter_input[inter_input$seed == 1,]
  inter_years = c(2016, 2019)  #unique(inter_input$year)
  
  png(paste0(hbhi_dir, '/simulation_inputs/plots/itn_anc_sim_input.png'), width = (length(inter_years)+1)*1.8, height=2.3, units='in', res=800)
  par(mfrow=c(1,(length(inter_years)+1)), mar=c(0,1,2,0))
  for(yy in 1:length(inter_years)){
    inter_input_cur = inter_input[inter_input$year == inter_years[yy],]
    
    inter_ordered = data.frame('ds_ordered'=admin_shapefile[[shapefile_admin_colname]], 'value'=rep(NA, length(admin_shapefile[[shapefile_admin_colname]])))
    for (i_ds in 1:length(inter_ordered$ds_ordered)){
      cur_ds = inter_ordered$ds_ordered[i_ds]
      inter_ordered$value[i_ds] = inter_input_cur[[coverage_colname]][which(toupper(inter_input_cur$admin_name) == toupper(cur_ds))]
    }
    col_cur = colorscale[sapply(floor((num_colors)*(inter_ordered$value - min_value) / (max_value - min_value))+1, min, num_colors)]
    col_cur[is.na(col_cur)] = 'grey'
    plot(admin_shapefile, col=col_cur, border=rgb(0.3,0.3,0.3), main=inter_years[yy])
  }
  
  # legend
  legend_label_vals = seq(min_value, max_value, length.out=5)
  legend_image = as.raster(matrix(rev(colorscale[sapply(floor((num_colors)*(legend_label_vals - min_value) / (max_value - min_value))+1, min, num_colors)]), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = c('coverage'))
  text(x=1.5, y = seq(0,1,length.out=5), labels = round(legend_label_vals,2))
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
  par(mfrow=c(1,1), mar=c(5,4,4,2))
  
  
  
  
  
  
}






