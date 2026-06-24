# get_MiP_protection_coverage.R

# functions to calculate what coverage should be assumed for IPTp and LLINs distributed during pregnancy.
#   This information is used during simulation post-processing to adjust simulation output for the impact of MiP, IPTp, LLINp, and mortality



#' New-style IPTp input loader
#'
#' Reads the two long-format CSVs produced by future_projection_intervention_
#' inputs_NMSP_coverage.R (future scenarios) and 0_main_DHS_data_to_sim_inputs.R
#' (to-present scenarios). The coordinator CSV stores a single basename
#' ('IPTp_filename'); the doses CSV is found by appending '_doses' to that
#' basename, matching the pairing convention used everywhere else in the
#' pipeline. Returns the data in the shape the downstream postprocessing
#' expects, so the call site can be a near drop-in replacement for
#' get_IPTp_coverages().
#'
#' @param iptp_filepath full path to the coverage CSV ('admin_name', 'year',
#'                      'coverage' columns; long format)
#' @param iptp_doses_filepath full path to the dose-distribution CSV
#'                      ('year', 'doses_1', 'doses_2', 'doses_3' columns)
#' @param first_year first simulation year
#' @param last_year last simulation year
#' @param admin_subset optional vector of admin_names to restrict the coverage
#'                      frame to (used when sims were run on a subset of LGAs)
#'
#' @return list with:
#'   coverage  -- data frame: admin_name, year, coverage  (long)
#'   doses     -- data frame: year, doses_1, doses_2, doses_3
load_iptp_inputs = function(iptp_filepath, iptp_doses_filepath,
                             first_year, last_year, admin_subset = NULL) {
  if (!file.exists(iptp_filepath))
    stop("IPTp coverage file not found: ", iptp_filepath)
  if (!file.exists(iptp_doses_filepath))
    stop("IPTp doses file not found: ", iptp_doses_filepath)

  cov   = read.csv(iptp_filepath,        as.is = TRUE)
  doses = read.csv(iptp_doses_filepath,  as.is = TRUE)

  # Normalise admin names (replace '/' with '-', matching what the legacy
  # get_IPTp_coverages did).
  if ('admin_name' %in% colnames(cov)) {
    cov$admin_name = gsub('/', '-', cov$admin_name)
  }
  if (!is.null(admin_subset)) {
    admin_subset = gsub('/', '-', admin_subset)
    cov = cov[cov$admin_name %in% admin_subset, , drop = FALSE]
  }

  # Filter to the simulation window.
  cov   = cov[cov$year   >= first_year & cov$year   <= last_year, , drop = FALSE]
  doses = doses[doses$year >= first_year & doses$year <= last_year, , drop = FALSE]

  # Basic sanity checks.
  required_cov = c('admin_name', 'year', 'coverage')
  missing_cov  = setdiff(required_cov, colnames(cov))
  if (length(missing_cov) > 0)
    stop("IPTp coverage file ", basename(iptp_filepath),
         " is missing required column(s): ",
         paste(missing_cov, collapse = ', '))

  required_doses = c('year', 'doses_1', 'doses_2', 'doses_3')
  missing_doses  = setdiff(required_doses, colnames(doses))
  if (length(missing_doses) > 0)
    stop("IPTp doses file ", basename(iptp_doses_filepath),
         " is missing required column(s): ",
         paste(missing_doses, collapse = ', '))

  dose_sums = doses$doses_1 + doses$doses_2 + doses$doses_3
  if (any(abs(dose_sums - 1) > 1e-3, na.rm = TRUE))
    warning("IPTp doses file ", basename(iptp_doses_filepath),
            ": doses_1 + doses_2 + doses_3 != 1 for some year(s). Max abs diff = ",
            signif(max(abs(dose_sums - 1), na.rm = TRUE), 3))

  list(coverage = cov, doses = doses)
}

#' Derive the doses-file path from the coverage-file path by appending '_doses'
#' before the '.csv' extension. Mirrors the convention enforced by the
#' generator scripts: the coordinator stores only the coverage basename.
iptp_doses_path_from_coverage_path = function(iptp_filepath) {
  sub("\\.csv$", "_doses.csv", iptp_filepath)
}




get_IPTp_coverages = function(iptp_estimates_filename, iptp_dose_number_filename, future_projection_flag, first_year, last_year, coverage_string, iptp_estimates_ci_l_filename=NA, iptp_estimates_ci_u_filename=NA, iptp_project_trajectory_filename=NA,
                              admin_subset_flag=FALSE, admin_subset=NA){
  #'  set IPTp coverage (and number of doses) through time for a particular scenario. Uses information on whether simulation 
  #'     is a past scenario or future projection as well as the character string describing the IPTp scenario.
  #'  
  #' @param iptp_estimates_filename filepath for estimates of IPTp coverage (>=1 dose) in each DS in past years
  #' @param iptp_dose_number_filename filepath for estimates of number of doses taken by people who receive IPTp
  #' @param future_projection_flag boolean indicating whether the simulation is a future projection
  #' @param first_year first year of simulation
  #' @param last_year final year of simulation
  #' @param coverage_string character string indicating which of several possible IPTp coverage scenarios should be applied
  #' @param iptp_estimates_ci_l_filename filename for lower bound estimates of IPTp coverage (>=1 dose) in each DS in past years (only used when coverage_string == 'estimateWorseCoverage')
  #' @param iptp_estimates_ci_u_filename filename for upper bound estimates of IPTp coverage (>=1 dose) in each DS in past years (only used when coverage_string == 'estimateBetterCoverage')
  #' @param admin_subset_flag boolean indicating whether admins should be subsetted to areas that received IRS or IPTi (for simulations run with a subset of admins)
  #' @param admin_subset vector of admin names included in the irs- or IPTi-only simulation. Must be specified if admin_subset_flag is TRUE
  #' 
  #' return list with:
  #'   1) IPTp coverage (>=1 doses) through time in each DS and 
  #'   2) fraction of individuals with 1,2, or 3 doses among individuals who receive IPTp in each year
  #'   3) names of DS/LGAs

  # read in estimated IPTp coverage
  iptp_coverage_df = read.csv(iptp_estimates_filename, as.is=TRUE)
  # get rid of first column if it is just row numbers
  if(all(iptp_coverage_df[1:10,1] == 1:10)){
    iptp_coverage_df = iptp_coverage_df[,-1]
  }
  admin_names = iptp_coverage_df$admin_name
  # replace slashes with dashes in admin names from admin_names and from iptp_coverage_df to match simulation output
  admin_names = gsub('/', '-', admin_names)
  iptp_coverage_df$admin_name = gsub('/', '-', iptp_coverage_df$admin_name)
  # check that all of the admin names from the IPTp file match with the simulation output
  # for(i_ds in 1:length(admin_names)){
  #   if(length(grep(toupper(admin_names[i_ds]), toupper(sim_output_all_runs$admin_name))) == 0){
  #     warning(print(paste('DS name not matched:', admin_names[i_ds])))
  #   }
  # }
  # if only analyzing admins where IRS was administered, subset appropriately
  if(admin_subset_flag){
    iptp_coverage_df = iptp_coverage_df[iptp_coverage_df$admin_name %in% admin_subset, ]
    admin_names = iptp_coverage_df$admin_name
  }
  
  # read in the number of IPTp doses taken by people who receive IPTp
  iptp_dose_number = read.csv(iptp_dose_number_filename, row.names=1)

  ## - - - - - - - - - - - - - - - - - - - - ##
  # simulations of historical transmission
  ## - - - - - - - - - - - - - - - - - - - - ##
  
  if (!future_projection_flag){
    if(coverage_string %in% c('noCoverage', 'none')){
      project_coverage = rep(0, length(iptp_coverage_df[,dim(iptp_coverage_df)[2]]))
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp coverage (>=1 dose)
        new_coverage_df = data.frame(project_coverage)
        colnames(new_coverage_df) = yy
        iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
      
    } else if(coverage_string == 'estimateWorseCoverage'){
      # read in estimated IPTp coverage
      iptp_coverage_df = read.csv(iptp_estimates_ci_l_filename, as.is=TRUE)[,-1]
      iptp_coverage_df$admin_name = gsub('/', '-', iptp_coverage_df$admin_name)
      iptp_relative_risk = iptp_relative_risk_estimateWorseProtection
      
    } else if(coverage_string == 'estimateBetterCoverage'){
      # read in estimated IPTp coverage
      iptp_coverage_df = read.csv(iptp_estimates_ci_u_filename, as.is=TRUE)[,-1]
      iptp_coverage_df$admin_name = gsub('/', '-', iptp_coverage_df$admin_name)
      iptp_relative_risk = iptp_relative_risk_estimateBetterProtection
      
    } else if(coverage_string != 'curCoverage') warning('coverage string not recognized for historical simulations... assuming estimated true coverage')
    
    # make sure that values are included for all years
    if (first_year < min(as.numeric(gsub('X','', colnames(iptp_coverage_df)[-1])))){
      warning('The first year of the IPTp coverage dataframe is after the first year of the simulation. Need to add more IPTp data/assumptions.')
    }
    if(last_year > max(as.numeric(gsub('X','', colnames(iptp_coverage_df)[-1])))){
      warning('The simulation continues for more years than supplied in the IPTp coverage dataframe. Assuming continuation of lastest coverage. To use a different assumption, need to update IPTp input file.')
      # use coverage from latest year for projections
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      
      for(yy in (max(as.numeric(gsub('X','', colnames(iptp_coverage_df)[-1])))+1):last_year){
        # update iptp coverage (>=1 dose)
        new_coverage_df = data.frame(project_coverage)
        colnames(new_coverage_df) = yy
        iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
    }
    
  ## - - - - - - - - - - - - - - - - - - - - ##
  # simulations of future transmission
  ## - - - - - - - - - - - - - - - - - - - - ##
  } else{
    constant_future_values = TRUE
    # replace entries with projection scenario
    if(coverage_string %in% c('noCoverage', 'none')){
      project_coverage = rep(0, length(iptp_coverage_df[,dim(iptp_coverage_df)[2]]))
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp coverage (>=1 dose)
        new_coverage_df = data.frame(project_coverage)
        colnames(new_coverage_df) = yy
        iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
    } else if(coverage_string =='curCoverage'){
      # use coverage from latest year for projections
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
    } else if(coverage_string =='increase_by10'){
      # increase probability of getting at least one IPTp does by 10%, up to a maximum 80%
      project_coverage = sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 1.1, min, 0.8)
      project_coverage[iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8] = iptp_coverage_df[which(iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8), dim(iptp_coverage_df)[2]]
      # also increase probability of getting 3 IPTp doses (given >=1 dose) by 10% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] * 1.1)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    }  else if(coverage_string =='increase_by20'){
      # increase probability of getting at least one IPTp does by 20%, up to a maximum 80%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 1.2, min, 0.8)
      project_coverage[iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8] = iptp_coverage_df[which(iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8), dim(iptp_coverage_df)[2]]
      # also increase probability of getting 3 IPTp doses (given >=1 dose) by 20% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] * 1.2)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    }  else if(coverage_string =='increase_by30'){
      # increase probability of getting at least one IPTp does by 30%, up to a maximum 80%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 1.3, min, 0.8)
      project_coverage[iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8] = iptp_coverage_df[which(iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8), dim(iptp_coverage_df)[2]]
      # also increase probability of getting 3 IPTp doses (given >=1 dose) by 30% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] * 1.3)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    } else if(coverage_string =='increase_10'){
      # increase probability of getting at least one IPTp 10%, up to a maximum 80%
      project_coverage = sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]] + 0.1, min, 0.8)
      project_coverage[iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8] = iptp_coverage_df[which(iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8), dim(iptp_coverage_df)[2]]
      # also increase probability of getting 3 IPTp doses (given >=1 dose)  10% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.1)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    }  else if(coverage_string =='increase_20'){
      # increase probability of getting at least one IPTp 20%, up to a maximum 80%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]] + 0.2, min, 0.8)
      project_coverage[iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8] = iptp_coverage_df[which(iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8), dim(iptp_coverage_df)[2]]
      # also increase probability of getting 3 IPTp doses (given >=1 dose) 20% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.2)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    }  else if(coverage_string =='increase_30'){
      # increase probability of getting at least one IPTp 30%, up to a maximum 80%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]]+ 0.3, min, 0.8)
      project_coverage[iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8] = iptp_coverage_df[which(iptp_coverage_df[,dim(iptp_coverage_df)[2]]>0.8), dim(iptp_coverage_df)[2]]
      # also increase probability of getting 3 IPTp doses (given >=1 dose) 30% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.3)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    } else if(coverage_string =='increase_to_50'){
      # increase probability of getting at least one IPTp up to 50%, unless already over 50%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]], max, 0.5)
      # also increase probability of getting 3 IPTp doses (given >=1 dose) 20% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.2)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    } else if(coverage_string =='increase_to_60'){
      # increase probability of getting at least one IPTp up to 60%, unless already over 60%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]], max, 0.6)
      # also increase probability of getting 3 IPTp doses (given >=1 dose) 20% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.2)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    } else if(coverage_string =='increase_to_70'){
      # increase probability of getting at least one IPTp up to 70%, unless already over 70%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]], max, 0.7)
      # also increase probability of getting 3 IPTp doses (given >=1 dose) 20% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.2)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    } else if(coverage_string =='increase_to_80'){
      # increase probability of getting at least one IPTp up to 80%, unless already over 80%
      project_coverage =  sapply(iptp_coverage_df[,dim(iptp_coverage_df)[2]], max, 0.8)
      # also increase probability of getting 3 IPTp doses (given >=1 dose) 30% and split remaining probability evenly between 1 and 2 doses
      new_three_dose_fraction = min(1,iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.3)
      project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
    } else if(coverage_string =='betterCoverage'){
      # specify coverage to use for projections
      project_coverage = rep(0.95, length(admin_names))
      project_dose_number = c(0.1, 0.2, 0.7)
    } else if(coverage_string =='betterCoverage2'){
      # specify coverage to use for projections
      if(country_name == 'Burkina'){
        # since coverage with >=1 IPTp is already higher than 80% in all Burkina admins (the average fraction of people with >=1 IPTp across admins is 0.9), 
        #   for this scenario assume everyone gets >=1 IPTp, and increase the fraction of people with higher number of doses from (0.12, 0.27, 0.6) 
        #   to (0.08, 0.18, 0.74). Overall, the fraction of people with >=3 IPTp increases from 0.9*0.6 = 0.54 to 1*0.74, so by 20%, in line with 
        #   that scenario.
        project_coverage = rep(1, length(admin_names))
        project_dose_number = c(0.08, 0.18, 0.74)
      } else if(country_name == 'Nigeria'){
        # coverage with >=1 IPTp is not as high in Nigeria as was reported for Burkina (the average fraction of people with >=1 IPTp across admin is 0.64).
        #   If we increase the fraction of people who get IPTp by 20%, and increase the fraction of IPTp-treated people who get three doses to 
        #   (0.64*d3+0.2)/(0.64+0.2)=0.44 , then the fraction of individuals with >=3 doses increases 20% (d3 represents the original fraction of 
        #   individuals with >=3 doses).This won't work exactly perfectly, because the coverage isn't equal in all admins, so some admins may have a larger 
        #   or smaller increase in the fraction of individuals with >=3 doses, but overall, I think it should be close to the requested scenario. 
        #   Also, when adding 20% to the >=1 IPTp coverage in each admin, some admins will truncate at 100% coverage, and that will bring the intervention 
        #   mean down a bit. 
        project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]] + 0.2
        project_coverage[project_coverage>1] = 1
        new_three_dose_fraction = (mean(iptp_coverage_df[,dim(iptp_coverage_df)[2]])*iptp_dose_number[3,dim(iptp_dose_number)[2]] + 0.2)/(0.2+mean(iptp_coverage_df[,dim(iptp_coverage_df)[2]]))
        project_dose_number = c((1-new_three_dose_fraction)/2, (1-new_three_dose_fraction)/2, new_three_dose_fraction)
        # # old fraction of individuals getting >=3 IPTp doses (mean across admins):  # 0.1716
        # mean(iptp_coverage_df[,dim(iptp_coverage_df)[2]]*iptp_dose_number[3,dim(iptp_dose_number)[2]])
        # # new fraction of individuals getting >=3 IPTp doses (mean across admins):  # 0.3667  <- around 0.2 greater than the original coverage, so looks good!
        # mean(project_coverage * project_dose_number[3])
      }
    } else if(coverage_string =='betterCoverage3'){
      # specify coverage to use for projections
      project_coverage = rep(1, length(admin_names))
      project_dose_number = c(0.05, 0.15, 0.8)
    # } else if(coverage_string == 'chwCoverage'){
    #   # for admins without CHW communication, use coverage from latest year for projections
    #   project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
    #   project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
    #   # for admins with communication, move frac_chw_reached from the 0 IPTp dose group into the >=1 IPTp group
    #   frac_chw_reached = 0.5
    #   # read in the csv file that tells which DS will be targeted
    #   CHW_communication_df = read.csv(projection_CHW_communication_filename, as.is=TRUE)
    #   # several DS have slightly different spellings... alter these now
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='BITOU'] = 'BITTOU'
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='GOROM'] = 'GOROM-GOROM'
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='NDOROLA'] = "N'DOROLA"
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='NONGR-MASSOM'] = 'NONGR-MASSOUM'
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='SIG-NOGHIN'] = 'SIG-NONGHIN'
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='MANNI'] = 'MANI'
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='FADA'] = "FADA N'GOURMA"
    #   CHW_communication_df$adm2[CHW_communication_df$adm2=='KARANGASSO VIGUE'] = 'KARANGASSO - VIGUE'
    #   # iterate through the DS, matching the admin_names value with the row from CHW_communication_df and increasing coverage if appropriate
    #   for(i_ds in 1:length(admin_names)){
    #     # find the row in CHW_communication_df
    #     cur_row = which(toupper(CHW_communication_df$adm2) == toupper(admin_names[i_ds]))
    #     if(length(cur_row) == 1){
    #       # check whether this DS gets enhanced CHW communication
    #       if(iconv(CHW_communication_df$com1[cur_row],to="ASCII//TRANSLIT") == 'Mis en A\"uvre'){
    #         project_coverage[i_ds] = project_coverage[i_ds] + (1-project_coverage[i_ds])*frac_chw_reached
    #       }
    #     } else warning(print(paste('ds name not macthed in chw communication df:', admin_names[i_ds])))
    #   }
    } else if(coverage_string =='coverage75'){
      # use coverage from latest year for projections, with fraction of individuals getting >=1 dose decreased to 75% of original
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 0.75
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
    } else if(coverage_string =='coverage50'){
      # use coverage from latest year for projections, with fraction of individuals getting >=1 dose decreased to 50% of original
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 0.5
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
    } else if(coverage_string =='coverage25'){
      # use coverage from latest year for projections, with fraction of individuals getting >=1 dose decreased to 25% of original
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 0.25
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
    } else if(coverage_string =='noCoverage'){
      # use coverage from latest year for projections, with fraction of individuals getting >=1 dose decreased to 0% of original
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]] * 0
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
    } else if (coverage_string == 'increase_to80'){
      # increase from current values to 80% coverage over 3 years
      constant_future_values = FALSE
      
      # values from past year
      project_coverage_past = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
      project_dose_number_past = iptp_dose_number[,dim(iptp_dose_number)[2]]
      
      # increase so that on thrid year, there is 80% coverage with >=1 dose and 80% of people  with >=1 dose get 3 doses
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(ii in 1:min(3,(last_year-first_year+1))){
        yy = (first_year:last_year)[ii]
        # update iptp coverage (>=1 dose)
        new_coverage_df = data.frame(project_coverage_past + ii/3 * sapply((0.8-project_coverage_past), max, 0))
        colnames(new_coverage_df) = yy
        iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        
        # update iptp doses (of covered individuals, how many doses received?)
        three_dose_values = project_dose_number_past[3] + ii/3 * max((0.8 - project_dose_number_past[3]),0)
        new_dose_df = data.frame(c((1-three_dose_values)/2, (1-three_dose_values)/2, three_dose_values))
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      if((last_year - first_year)>3){
        for(yy in (first_year+3):last_year){
          # remaining years are at the maximum coverage
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(sapply(project_coverage_past, max, 0.8))
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
          
          # update iptp doses (of covered individuals, how many doses received?)
          three_dose_values = max(0.8, project_dose_number_past[3])
          new_dose_df = data.frame(c((1-three_dose_values)/2, (1-three_dose_values)/2, three_dose_values))
          colnames(new_dose_df) = yy
          iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
        }
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
    } else if(coverage_string == 'trajectoryCoverage'){
      constant_future_values = FALSE
      
      iptp_projection_coverage_df = read.csv(iptp_project_trajectory_filename, as.is=TRUE)#[,-1]
      iptp_projection_coverage_df$admin_name = gsub('/', '-', iptp_projection_coverage_df$admin_name)
      
      included_columns = sapply(c('admin_name', first_year:last_year), grep, colnames(iptp_projection_coverage_df))
      iptp_coverage_df = iptp_projection_coverage_df[, included_columns]
      colnames(iptp_coverage_df) = gsub('X', '', colnames(iptp_coverage_df))
      
      # keep same dose number distribution as final year
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
      
    } else if(coverage_string == 'trajectoryCoverage_down25'){
      constant_future_values = FALSE
      
      iptp_projection_coverage_df = read.csv(iptp_project_trajectory_filename, as.is=TRUE)[,-1]
      iptp_projection_coverage_df$admin_name = gsub('/', '-', iptp_projection_coverage_df$admin_name)
      
      included_columns = sapply(c('admin_name', first_year:last_year), grep, colnames(iptp_projection_coverage_df))
      iptp_coverage_df = iptp_projection_coverage_df[, included_columns]
      iptp_coverage_df[2:dim(iptp_coverage_df)[2]] = iptp_coverage_df[2:dim(iptp_coverage_df)[2]] * 0.75
      colnames(iptp_coverage_df) = gsub('X', '', colnames(iptp_coverage_df))
      
      # keep same dose number distribution as final year
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
      
    } else if(coverage_string == 'trajectoryCoverage_down50'){
      constant_future_values = FALSE
      
      iptp_projection_coverage_df = read.csv(iptp_project_trajectory_filename, as.is=TRUE)[,-1]
      iptp_projection_coverage_df$admin_name = gsub('/', '-', iptp_projection_coverage_df$admin_name)
      
      included_columns = sapply(c('admin_name', first_year:last_year), grep, colnames(iptp_projection_coverage_df))
      iptp_coverage_df = iptp_projection_coverage_df[, included_columns]
      iptp_coverage_df[2:dim(iptp_coverage_df)[2]] = iptp_coverage_df[2:dim(iptp_coverage_df)[2]] * 0.5
      colnames(iptp_coverage_df) = gsub('X', '', colnames(iptp_coverage_df))
      
      # keep same dose number distribution as final year
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
      
    } else if(coverage_string == 'trajectoryCoverage_down75'){
      constant_future_values = FALSE
      
      iptp_projection_coverage_df = read.csv(iptp_project_trajectory_filename, as.is=TRUE)[,-1]
      iptp_projection_coverage_df$admin_name = gsub('/', '-', iptp_projection_coverage_df$admin_name)
      
      included_columns = sapply(c('admin_name', first_year:last_year), grep, colnames(iptp_projection_coverage_df))
      iptp_coverage_df = iptp_projection_coverage_df[, included_columns]
      iptp_coverage_df[2:dim(iptp_coverage_df)[2]] = iptp_coverage_df[2:dim(iptp_coverage_df)[2]] * 0.25
      colnames(iptp_coverage_df) = gsub('X', '', colnames(iptp_coverage_df))
      
      # keep same dose number distribution as final year
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
      
    } else if (coverage_string == 'stopAfterOneYear'){
      # keep current values for 1 year, then discontinue
      constant_future_values = FALSE
      
      # use coverage from latest year for first year of projections
      num_years_cur_coverage = 1
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
      project_0_coverage = rep(0, length(iptp_coverage_df[,dim(iptp_coverage_df)[2]]))
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        if(yy < (first_year+num_years_cur_coverage)){
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(project_coverage)
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        } else{
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(project_0_coverage)
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)        
        }
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
    } else if (coverage_string == 'stopAfterTwoYears'){
      # keep current values for 2 years, then discontinue
      constant_future_values = FALSE
      
      # use coverage from latest year for first two years of projections
      num_years_cur_coverage = 2
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
      project_0_coverage = rep(0, length(iptp_coverage_df[,dim(iptp_coverage_df)[2]]))
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        if(yy < (first_year+num_years_cur_coverage)){
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(project_coverage)
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        } else{
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(project_0_coverage)
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)        
        }
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
    } else if (coverage_string == 'stopAfterThreeYears'){
      # keep current values for 2 years, then discontinue
      constant_future_values = FALSE
      
      # use coverage from latest year for first two years of projections
      num_years_cur_coverage = 3
      project_coverage = iptp_coverage_df[,dim(iptp_coverage_df)[2]]
      project_0_coverage = rep(0, length(iptp_coverage_df[,dim(iptp_coverage_df)[2]]))
      project_dose_number = iptp_dose_number[,dim(iptp_dose_number)[2]]
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        if(yy < (first_year+num_years_cur_coverage)){
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(project_coverage)
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        } else{
          # update iptp coverage (>=1 dose)
          new_coverage_df = data.frame(project_0_coverage)
          colnames(new_coverage_df) = yy
          iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)        
        }
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
    } else warning(print('coverage string not recognized'))
    
    if (constant_future_values){
      iptp_coverage_df = data.frame('admin_name' = admin_names)
      iptp_dose_number = data.frame('doses' = c(1,2,3))
      for(yy in first_year:last_year){
        # update iptp coverage (>=1 dose)
        new_coverage_df = data.frame(project_coverage)
        colnames(new_coverage_df) = yy
        iptp_coverage_df = cbind(iptp_coverage_df, new_coverage_df)
        
        # update iptp doses (of covered individuals, how many doses received?)
        new_dose_df = data.frame(project_dose_number)
        colnames(new_dose_df) = yy
        iptp_dose_number = cbind(iptp_dose_number, new_dose_df)
      }
      # remove first column
      iptp_dose_number = iptp_dose_number[,-1]
    }
  }
  write.csv(iptp_coverage_df, gsub('estimated_past_IPTp_each_DS', paste0('sim_assumptions/IPTp_', coverage_string, '_',first_year, '_', last_year), iptp_estimates_filename))
  write.csv(iptp_dose_number, gsub('estimated_past_IPTp_each_DS', paste0('sim_assumptions/IPTp_dose_', coverage_string, '_',first_year, '_', last_year), iptp_estimates_filename))
  
  return(list(iptp_coverage_df, iptp_dose_number, admin_names))
}



get_LLINp_coverages = function(LLINp_filename, admin_names){
  #TODO: finish! this is a placeholder for a function that will return the probability a pregnant individual receives an LLIN through LLIN-in-pregnancy distribution, ordered the same as admin_names
  LLINp_each_ds = read.csv(LLINp_filename, as.is=TRUE)

  return(LLINp_coverage)
}





