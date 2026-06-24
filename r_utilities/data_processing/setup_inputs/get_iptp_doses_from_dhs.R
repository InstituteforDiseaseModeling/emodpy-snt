library(haven)
library(readr)
library(dplyr)

# ---------------------------------------------------------------------------
# Internal helpers: CSPro fixed-width file parsing
# Mirrors logic in 1_DHS_data_extraction.R so this script is self-contained.
# ---------------------------------------------------------------------------

.parse_dcf_layout_iptp <- function(dcf_path) {
  lines <- readLines(dcf_path, warn = FALSE)
  item_starts <- grep("^\\[Item\\]", lines)
  if (length(item_starts) == 0) stop("No [Item] sections found in DCF: ", dcf_path)
  all_section_starts <- sort(c(item_starts,
    grep("^\\[ValueSet\\]|^\\[Level\\]|^\\[Record\\]|^\\[Dictionary\\]|^\\[IdItems\\]", lines)))
  rows <- lapply(seq_along(item_starts), function(i) {
    s <- item_starts[i]
    nexts <- all_section_starts[all_section_starts > s]
    e     <- if (length(nexts) > 0) nexts[1] - 1L else min(s + 40L, length(lines))
    block <- lines[s:e]
    nm <- sub("^Name=",  "", block[grep("^Name=",  block)[1]])
    st <- suppressWarnings(as.integer(sub("^Start=", "", block[grep("^Start=", block)[1]])))
    ln <- suppressWarnings(as.integer(sub("^Len=",   "", block[grep("^Len=",   block)[1]])))
    if (!is.na(nm) && !is.na(st) && !is.na(ln))
      data.frame(Name = nm, Start = st, Len = ln, stringsAsFactors = FALSE)
    else NULL
  })
  layout <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(layout) || nrow(layout) == 0)
    stop("No valid Name/Start/Len entries parsed from DCF: ", dcf_path,
         "\n  Check that the file uses standard CSPro key names (Name=, Start=, Len=).")
  layout
}

.parse_sps_layout_iptp <- function(sps_path) {
  lines <- readLines(sps_path, warn = FALSE)
  dl_line <- grep("(?i)^\\s*DATA\\s+LIST", lines, perl = TRUE)[1]
  if (is.na(dl_line)) stop("No DATA LIST found in SPS file: ", sps_path)
  block <- character(0)
  for (i in dl_line:length(lines)) {
    ln <- lines[i]
    block <- c(block, ln)
    if (grepl("^\\s*\\.\\s*$", ln)) break
    if (i > dl_line && !grepl("^\\s", ln) && grepl("^[A-Za-z]", trimws(ln))) {
      block <- block[-length(block)]; break
    }
  }
  kw <- toupper(c("DATA","LIST","FILE","FIXED","FREE","RECORDS","FORMAT",
                  "TABLE","NOTABLE","END","MISSING","VALUES","VARIABLE",
                  "LABELS","VALUE","FORMATS","HANDLE","NAME"))
  rows <- list()
  for (ln in block) {
    ln <- gsub("/\\d+", " ", ln)
    tokens <- strsplit(trimws(ln), "\\s+")[[1]]
    tokens <- tokens[nchar(tokens) > 0]
    i <- 1
    while (i <= length(tokens)) {
      tok <- tokens[i]
      if (toupper(tok) %in% kw || grepl("^\\d", tok) ||
          grepl("^[(/]", tok) || !grepl("^[A-Za-z_]", tok)) { i <- i + 1; next }
      if (i + 1 <= length(tokens)) {
        col_tok <- tokens[i + 1]
        if (grepl("^\\d+-\\d+$", col_tok)) {
          p     <- strsplit(col_tok, "-")[[1]]
          start <- as.integer(p[1]); end <- as.integer(p[2])
          implicit_dec <- 0L
          if (i + 2 <= length(tokens) && grepl("^\\(\\d+\\)$", tokens[i + 2])) {
            implicit_dec <- as.integer(gsub("[()]", "", tokens[i + 2])); i <- i + 1
          }
          rows[[length(rows) + 1]] <- data.frame(
            Name = tok, Start = start, Len = end - start + 1L,
            ImplicitDec = implicit_dec, stringsAsFactors = FALSE)
          i <- i + 2; next
        } else if (grepl("^\\d+$", col_tok)) {
          start <- as.integer(col_tok)
          rows[[length(rows) + 1]] <- data.frame(
            Name = tok, Start = start, Len = 1L,
            ImplicitDec = 0L, stringsAsFactors = FALSE)
          i <- i + 2; next
        }
      }
      i <- i + 1
    }
  }
  if (length(rows) == 0) stop("No variable definitions parsed from SPS file: ", sps_path)
  do.call(rbind, rows)
}

.read_cspr_iptp <- function(dat_path, dcf_path = NULL) {
  sps_path <- sub("(?i)\\.dat$", ".SPS", dat_path, perl = TRUE)
  if (!file.exists(sps_path))
    sps_path <- sub("(?i)\\.dat$", ".sps", dat_path, perl = TRUE)
  if (file.exists(sps_path)) {
    layout  <- .parse_sps_layout_iptp(sps_path)
    col_pos <- readr::fwf_positions(
      start     = layout$Start,
      end       = layout$Start + layout$Len - 1L,
      col_names = layout$Name)
  } else {
    if (is.null(dcf_path)) {
      dcf_path <- sub("(?i)\\.dat$", ".dcf", dat_path, perl = TRUE)
      if (!file.exists(dcf_path))
        dcf_path <- sub("(?i)\\.dat$", ".DCF", dat_path, perl = TRUE)
    }
    if (!file.exists(dcf_path))
      stop("No .sps or .dcf layout file found for: ", dat_path)
    layout  <- .parse_dcf_layout_iptp(dcf_path)
    layout  <- layout[order(layout$Start), , drop = FALSE]
    col_pos <- readr::fwf_positions(
      start     = layout$Start,
      end       = layout$Start + layout$Len - 1L,
      col_names = layout$Name)
  }
  df <- readr::read_fwf(
    dat_path,
    col_positions  = col_pos,
    col_types      = readr::cols(.default = readr::col_character()),
    progress       = FALSE,
    show_col_types = FALSE)
  as.data.frame(dplyr::mutate(df, dplyr::across(where(is.character), trimws)))
}

.load_survey_file_iptp <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "dta") {
    as.data.frame(haven::read_dta(path))
  } else if (ext == "sav") {
    as.data.frame(haven::read_sav(path))
  } else if (ext == "dat") {
    .read_cspr_iptp(path)
  } else {
    stop("Unsupported file extension: .", ext, "  (supported: .dta, .sav, .dat)")
  }
}

# ---------------------------------------------------------------------------

get_iptp_doses_from_dhs = function(hbhi_dir, dta_dir, years, sim_start_year, last_sim_year, country='NGA'){
  if(!dir.exists(paste0(hbhi_dir,'/simulation_inputs/IPTp'))) dir.create(paste0(hbhi_dir,'/simulation_inputs/IPTp'))

  # calculate fraction of all individuals with each number of IPTp doses
  store_iptp_fractions = matrix(NA, ncol = length(years), nrow=3)
  for(yy in 1:length(years)){
    DHS_file_recode_df = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_',years[yy],'_files_recodes_for_sims.csv'))
    var_index = which(DHS_file_recode_df$variable == 'iptp_doses')
    file_path = paste0(dta_dir, '/', DHS_file_recode_df$folder_dir[var_index], '/', DHS_file_recode_df$filename[var_index])
    cur_dta = .load_survey_file_iptp(file_path)
    code_name = DHS_file_recode_df$code[var_index]
    cur_dta[[code_name]] = suppressWarnings(as.numeric(cur_dta[[code_name]]))

    # get rid of values that are clearly wrong
    cur_dta[[code_name]][cur_dta[[code_name]]> 10] = NA
    denom = length(which(cur_dta[[code_name]] >=1))
    # only include results if there are 30 or more datapoints
    if(denom>30){
      store_iptp_fractions[, yy] = c(  length(which(cur_dta[[code_name]] ==1))/denom,
                                       length(which(cur_dta[[code_name]] ==2))/denom,
                                       length(which(cur_dta[[code_name]] >=3))/denom)
    } else{
      print(paste0('Fewer than 30 individuals included in sample for year ', years[yy], '. This year excluded from fit.'))
    }
  }
  # assume same dose distribution in last_sim_year as in last observation year
  if(last_sim_year > max(years)){
    store_iptp_fractions = cbind(store_iptp_fractions, store_iptp_fractions[,dim(store_iptp_fractions)[2]])
    colnames_iptp =  c(years, last_sim_year)
  } else{
    colnames_iptp = years
  }
  # label rows and columns of matrix
  colnames(store_iptp_fractions) = colnames_iptp
  rownames(store_iptp_fractions) = c('1 dose', '2 doses', '3 doses')
  
  # remove NA columns (survey years where fewer than 30 individuals were sampled)
  na_cols = which(is.na(store_iptp_fractions[1,]))
  if(length(na_cols) >= 1) {
    store_iptp_fractions = store_iptp_fractions[,-na_cols]
    years_included = colnames_iptp[-na_cols]
  } else{
    years_included = colnames_iptp
  }
  
  if(country !='NGA'){
    # smoothed function through these points to get the fraction given each dosage for each year sim_start_year-last_sim_year
    # xx=years_included
    iptp_doses_all_years = matrix(NA, nrow=3, ncol=length(sim_start_year:last_sim_year))
    for(i_dose in 1:dim(store_iptp_fractions)[1]){
      iptp_doses_all_years[i_dose, ] = splinefun(years_included, store_iptp_fractions[i_dose,], method="hyman")(sim_start_year:last_sim_year)
    }
    
    # check result
    if(!all(apply(iptp_doses_all_years, 2, sum)==1)){  # should all be equal to 1, if not, need to divide each column by column sum
      warning('PROBLEM DETECTED: in fraction of individuals receiving 1, 2, or 3 doses. Sum is not 1.')
    }
    # label rows and columns of matrix
    colnames(iptp_doses_all_years) = sim_start_year:last_sim_year
    rownames(iptp_doses_all_years) = c('1 dose', '2 doses', '3 doses')
    
    
    # write results to csv
    if(!is.na(iptp_doses_all_years[1,1])){
      write.csv(iptp_doses_all_years, paste(hbhi_dir,'/simulation_inputs/IPTp/estimated_past_num_doses.csv', sep=''))
      
    } else{
      warning('unknown reason for NAs in smoothed IPTp dose matrix - not writing csv output file.')
    }
  } else{
    # For Nigeria: use a two-epoch, piecewise-smoothed profile.
    # Year-to-year fluctuations in dose distribution across DHS/MIS surveys are unlikely
    # to reflect true trends (DHS and MIS surveys show different apparent coverage levels).
    # Rather than a single grand mean, split surveys into two epochs:
    #   "recent" = final survey year and any within 5 years of it
    #   "early"  = all remaining surveys
    # Within each epoch, combine the survey-specific dose fractions by summing the
    # per-survey fractions (rowSums across survey columns) and renormalising, i.e.
    # an unweighted aggregation of survey-level fractions rather than summing counts.
    # Simulation years are assigned the early profile for years <= the last early
    # survey year, the recent profile for years >= the first recent survey year
    # (including all future years), and a linear interpolation in between. This
    # captures big-picture directional change while suppressing noise.

    # Only include actual survey years (exclude any synthetic last_sim_year column)
    survey_yr_cols  = as.numeric(colnames(store_iptp_fractions)) %in% as.numeric(years)
    fracs_surveys   = store_iptp_fractions[, survey_yr_cols, drop = FALSE]
    survey_yrs_used = as.numeric(colnames(fracs_surveys))

    max_survey_yr = max(survey_yrs_used)
    recent_mask   = survey_yrs_used >= (max_survey_yr - 5)
    early_mask    = !recent_mask

    if (sum(early_mask) == 0) {
      warning('All IPTp survey years are within 5 years of the final survey; using single mean for NGA.')
      early_mask = recent_mask
    }

    recent_mean = rowSums(fracs_surveys[, recent_mask, drop = FALSE])
    recent_mean = recent_mean / sum(recent_mean)
    early_mean  = rowSums(fracs_surveys[, early_mask,  drop = FALSE])
    early_mean  = early_mean  / sum(early_mean)

    year_early  = max(survey_yrs_used[early_mask])   # early mean holds through the last early survey year
    year_recent = min(survey_yrs_used[recent_mask])  # recent mean applies from the first recent survey year

    sim_years = sim_start_year:last_sim_year
    iptp_doses_all_years = matrix(NA, nrow=3, ncol=length(sim_years))
    for(i_dose in 1:3){
      for(yy in seq_along(sim_years)){
        yr = sim_years[yy]
        if(yr <= year_early){
          iptp_doses_all_years[i_dose, yy] = early_mean[i_dose]
        } else if(yr >= year_recent){
          iptp_doses_all_years[i_dose, yy] = recent_mean[i_dose]
        } else{
          t = (yr - year_early) / (year_recent - year_early)
          iptp_doses_all_years[i_dose, yy] = (1 - t) * early_mean[i_dose] + t * recent_mean[i_dose]
        }
      }
    }
    # label rows and columns of matrix
    colnames(iptp_doses_all_years) = sim_start_year:last_sim_year
    rownames(iptp_doses_all_years) = c('1 dose', '2 doses', '3 doses')
    write.csv(iptp_doses_all_years, paste(hbhi_dir,'/simulation_inputs/IPTp/estimated_past_num_doses.csv', sep=''))
  }

  # ── Long-format doses companion file ─────────────────────────────────────
  # Matches the convention used by the new IPTp pipeline: one row per year
  # with columns doses_1 / doses_2 / doses_3. Filename pairs with the long-
  # format coverage file written by create_iptp_input_from_DHS, so the
  # coordinator can store a single 'iptp_2010_toPresent' basename and the
  # consumer constructs the doses path by appending '_doses'.
  if(exists('iptp_doses_all_years') && !is.null(iptp_doses_all_years)){
    doses_long = data.frame(
      year    = as.integer(colnames(iptp_doses_all_years)),
      doses_1 = as.numeric(iptp_doses_all_years[1, ]),
      doses_2 = as.numeric(iptp_doses_all_years[2, ]),
      doses_3 = as.numeric(iptp_doses_all_years[3, ])
    )
    out_dir = paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent')
    if(!dir.exists(out_dir)) dir.create(out_dir)
    write.csv(doses_long,
              paste0(out_dir, '/iptp_2010_toPresent_doses.csv'),
              row.names = FALSE)
  }
 
  # # plot the number of IPTp doses reported in different surveys
  # par(mar=c(5,5,4,3))
  # layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
  # barplot(store_iptp_fractions, beside=FALSE, xlab='year', ylab='fraction of individuals receiving dose', col= c(rgb(145/255, 48/255, 88/255, 0.5), rgb(0/255, 160/255, 138/255, 0.5), rgb(83/255, 147/255, 195/255)), cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
  # plot(NA, xlim=c(2008,2018), ylim=c(0,0.6), xlab='', ylab='', axes=FALSE)
  # legend('center', legend=rev(c('1', '2', '>=3')), fill=rev(c(rgb(145/255, 48/255, 88/255, 0.5), rgb(0/255, 160/255, 138/255, 0.5), rgb(83/255, 147/255, 195/255))), title='IPTp doses taken', bty='n', cex=1.8)
  # par(mfrow=c(1,1))
  
  png(paste0(hbhi_dir, '/simulation_inputs/plots/IPTp_barchart_fraction_num_doses_DHSMIS.png'), width=7, height=4, units='in', res=900)
  layout(matrix(c(1,1,2), nrow=1))
  if(last_sim_year > max(years)){
    bp = barplot(store_iptp_fractions[,1:(dim(store_iptp_fractions)[2]-1)], col=(c(rgb(0,0,0,0.75), rgb(0.5,0.5,0.5,0.5), rgb(0.7,0.9, 1, 1))), xlab='year', ylab='IPTp doses', names.arg=years_included[1:(dim(store_iptp_fractions)[2]-1)], cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
  } else{
    bp = barplot(store_iptp_fractions[,1:(dim(store_iptp_fractions)[2])], col=(c(rgb(0,0,0,0.75), rgb(0.5,0.5,0.5,0.5), rgb(0.7,0.9, 1, 1))), xlab='year', ylab='IPTp doses', names.arg=years_included[1:(dim(store_iptp_fractions)[2])], cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
  }
  if(country =='NGA'){
    # draw a segment over each bar showing the assumed dose boundaries used in simulations
    for(ii in 1:length(bp)){
      yr = as.numeric(years_included[ii])
      if(yr <= year_early){
        assumed = early_mean
      } else if(yr >= year_recent){
        assumed = recent_mean
      } else{
        tt = (yr - year_early) / (year_recent - year_early)
        assumed = (1 - tt) * early_mean + tt * recent_mean
      }
      cum_vals = cumsum(assumed)[1:2]  # boundaries between 1/2 doses and 2/3+ doses
      segments(x0=bp[ii]-0.5, x1=bp[ii]+0.5, y0=cum_vals, y1=cum_vals, col='red', lwd=2)
    }
  }
  plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE, ylab='', xlab='')
  if(country =='NGA'){
    legend(-0.01,0.95, c('>=3 doses', '2 doses', '1 dose', 'sim assumed'), col=c(rev(c(rgb(0,0,0,0.75), rgb(0.5,0.5,0.5,0.5), rgb(0.7,0.9, 1, 1))), 'red'), lwd=c(5,5,5,2), bty='n', cex=1.5)
  } else{
    legend(-0.01,0.95, c('>=3 doses', '2 doses', '1 dose'), col=rev(c(rgb(0,0,0,0.75), rgb(0.5,0.5,0.5,0.5), rgb(0.7,0.9, 1, 1))), lwd=5, bty='n', cex=1.5)
  }
  dev.off()
}






