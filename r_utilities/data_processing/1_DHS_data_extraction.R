# =============================================================================
# 1_DHS_data_extraction.R
#
# Extracts, aggregates, and exports survey data from DHS, MIS, and MICS
# household surveys for use as malaria simulation inputs.  Handles multiple
# file formats, applies hierarchical sample-size fallbacks, and produces
# standardised CSV outputs and diagnostic plots.
#
# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------
#   extract_DHS_data                    — main extraction: all variables, all years
#   extract_archetype_level_DHS_data    — aggregate to archetype level
#   extract_country_level_DHS_data      — aggregate to national level
#   extract_state_level_DHS_data        — aggregate to state/region level
#   extract_vaccine_DHS_data            — vaccination coverage from DHS
#   extract_vaccine_MICS_data           — vaccination coverage from MICS
#   extract_frac_itn_from_campaign_by_state — fraction of ITNs from campaigns
#   plot_extracted_DHS_data             — choropleth maps of extracted estimates
#   plot_state_level_DHS_data           — state-level choropleth maps
#   create_DHS_reference_monthly_pfpr   — monthly PfPR reference curves by admin
#
# -----------------------------------------------------------------------------
# SUPPORTED FILE FORMATS (auto-detected from extension)
# -----------------------------------------------------------------------------
#   .dta         — Stata       (classic DHS / MIS)
#   .sav         — SPSS        (MICS)
#   .dat + .dcf  — CSPro fixed-width  (MIS 2025+)
#
# -----------------------------------------------------------------------------
# RECODE CSV
# -----------------------------------------------------------------------------
# Each survey year requires a recode CSV at:
#   {hbhi_dir}/estimates_from_DHS/DHS_{year}_files_recodes_for_sims.csv
#
# Required columns:
#   variable            — short name used throughout the pipeline
#   folder_dir          — path to data file, relative to dta_dir
#   filename            — data file name (e.g. PERSONS.DAT, hh.dta)
#   code                — column name in the data file for the outcome
#   pos_pattern         — value representing a positive result
#   neg_pattern         — value representing a negative result
#   cluster_id_code     — column name for the cluster/PSU identifier
#   month_code          — column name for interview month
#   year_code           — column name for interview year
#   age_code            — column name for individual age (where applicable)
#   state_code          — column name for state/region (where applicable)
#   min_age_months_to_include — minimum age (months) for inclusion
#
# Columns used by the "locations" row (GPS coordinates):
#   latitude_code       — column name for GPS latitude
#   longitude_code      — column name for GPS longitude
#   household_id_code   — column name for household identifier
#
# Optional columns:
#   admin_code          — column name in the survey data file containing the
#                         MIS-recorded LGA/admin name.  When provided in the
#                         "locations" row (CSPro format), this name is used
#                         inside extract_gps_from_nets() as a cross-check
#                         against GPS-derived admin assignments: clusters whose
#                         GPS-assigned state does not match the MIS-recorded
#                         state are flagged and filtered before the cluster
#                         centroid is computed.  It does NOT replace the
#                         spatial join — admin unit assignment always uses the
#                         shapefile spatial join when GPS coordinates are
#                         available.
#   region_code         — column name in the survey data file containing the
#                         MIS-recorded state/region name (companion to
#                         admin_code, used in the same cross-check).
#
# NOTE on CSPro format:
#   pos_pattern / neg_pattern must be the numeric codes stored in the DAT file
#   (e.g. "1", "2"), not string labels (e.g. "positive", "yes").
#   For .dta / .sav files, string labels work as before.
#
# NOTE on admin_shape:
#   All externally-called functions require an admin_shape
#   SpatialPolygonsDataFrame with NOMDEP (admin name) and NOMREGION
#   (state/region name) columns.  It is used for spatial join (cluster-to-
#   admin assignment), state-consistency cross-checking in GPS processing,
#   GPS diagnostic plots, and all map outputs.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)      # read_dta, read_sav
  library(dplyr)
  library(stringr)
  library(readr)      # read_fwf — CSPro DAT files
  library(pals)
  library(prettyGraphs)
  library(raster)     # shapefile, crs, over — kept for back-compat with sp workflow
  library(sp)
})

# Internal null-coalescing operator
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x[1])) x else y


# =============================================================================
# SECTION 1: FILE LOADING & PARSING
# Read survey data files (.dta, .sav, .dat) and their layout dictionaries.
# Provides helpers to build file paths and normalise date columns for CSV output.
# =============================================================================
# Parse a CSPro .dcf data-dictionary file into a variable layout table.
# Inputs: dcf_path -- path to the .dcf file
# Returns: data frame with columns Name, Start, Len (one row per [Item] entry).
# Only [Item] entries are extracted; sub-items and occurrences are ignored.
parse_dcf_layout <- function(dcf_path) {
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
  do.call(rbind, rows)
}
# Parse an SPSS syntax (.sps) DATA LIST block to extract variable positions.
# Inputs: sps_path -- path to the .sps file
# Returns: data frame with columns Name, Start, Len, ImplicitDec.
# ImplicitDec > 0: stored integer must be divided by 10^ImplicitDec.
# implicit decimals → -12.345678).
parse_sps_layout <- function(sps_path) {
  lines <- readLines(sps_path, warn = FALSE)
  
  dl_line <- grep("(?i)^\\s*DATA\\s+LIST", lines, perl = TRUE)[1]
  if (is.na(dl_line)) stop("No DATA LIST found in SPS file: ", sps_path)
  
  # Collect all lines in the DATA LIST block until end-of-command (bare period)
  # or a new top-level command (non-whitespace start, after dl_line).
  block <- character(0)
  for (i in dl_line:length(lines)) {
    ln <- lines[i]
    block <- c(block, ln)
    if (grepl("^\\s*\\.\\s*$", ln)) break
    if (i > dl_line && !grepl("^\\s", ln) && grepl("^[A-Za-z]", trimws(ln))) {
      block <- block[-length(block)]; break
    }
  }
  
  # Keywords that can appear on DATA LIST header line — not variable names
  kw <- toupper(c("DATA","LIST","FILE","FIXED","FREE","RECORDS","FORMAT",
                  "TABLE","NOTABLE","END","MISSING","VALUES","VARIABLE",
                  "LABELS","VALUE","FORMATS","HANDLE","NAME"))
  
  rows <- list()
  for (ln in block) {
    ln <- gsub("/\\d+", " ", ln)          # strip record indicators (/1, /2 …)
    tokens <- strsplit(trimws(ln), "\\s+")[[1]]
    tokens <- tokens[nchar(tokens) > 0]
    i <- 1
    while (i <= length(tokens)) {
      tok <- tokens[i]
      # Skip SPSS keywords, column-spec tokens, format tokens, or non-name tokens
      if (toupper(tok) %in% kw || grepl("^\\d", tok) ||
          grepl("^[(/]", tok) || !grepl("^[A-Za-z_]", tok)) {
        i <- i + 1; next
      }
      # Next token should be column spec  start-end  or  start
      if (i + 1 <= length(tokens)) {
        col_tok <- tokens[i + 1]
        if (grepl("^\\d+-\\d+$", col_tok)) {
          p     <- strsplit(col_tok, "-")[[1]]
          start <- as.integer(p[1]); end <- as.integer(p[2])
          # Optional format token: (N) or (A) immediately following
          implicit_dec <- 0L
          if (i + 2 <= length(tokens) && grepl("^\\(\\d+\\)$", tokens[i + 2])) {
            implicit_dec <- as.integer(gsub("[()]", "", tokens[i + 2]))
            i <- i + 1
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

# Parse VALUE LABELS blocks from an SPS syntax file.
# Inputs: sps_path -- path to the .sps file
# Returns: named list VARNAME -> named character vector of code -> label mappings.
# Used to decode numeric codes into human-readable labels (e.g. state, LGA names).
parse_sps_value_labels <- function(sps_path) {
  lines    <- readLines(sps_path, warn = FALSE, encoding = "latin1")
  result   <- list()
  in_block <- FALSE
  cur_var  <- NULL
  
  new_cmd <- paste0("^(",
                    "VARIABLE +LABELS?|MISSING +VALUES?|EXECUTE|DATA +LIST|",
                    "FREQUENCIES|FORMATS|RECODE|COMPUTE|IF\\b|FILTER|WEIGHT|",
                    "DESCRIPTIVES|CROSSTABS|SELECT|SAVE|GET|SET|SHOW|TITLE|",
                    "\\*|\\$)",
                    collapse = "")
  
  for (line in lines) {
    lt <- trimws(line)
    lu <- toupper(lt)
    
    # Detect start of a VALUE LABELS block (also ADD VALUE LABELS)
    if (grepl("^(ADD +)?VALUE +LABELS?\\b", lu, perl = TRUE)) {
      in_block <- TRUE
      next
    }
    # Detect end: a new SPSS command at the start of a line
    if (in_block && nchar(lt) > 0 && grepl(new_cmd, lu, perl = TRUE)) {
      in_block <- FALSE; cur_var <- NULL
    }
    if (!in_block || nchar(lt) == 0) next
    
    if (startsWith(lt, "/")) {
      # New variable section: /VARNAME  (may have multiple names or a description after)
      var_part <- trimws(sub("^/\\s*", "", lt))
      cur_var  <- strsplit(var_part, "[\\s\"']")[[1]][1]   # take only the name token
      if (nchar(cur_var) > 0 && !cur_var %in% names(result))
        result[[cur_var]] <- character(0)
    } else if (grepl("^[A-Za-z_][A-Za-z0-9_$@.]*$", lt)) {
      # Variable name on its own line (no leading /) — e.g. the indented format:
      #   VALUE LABELS
      #     QHREGION
      #        1 "SOKOTO"
      cur_var <- lt
      if (!cur_var %in% names(result)) result[[cur_var]] <- character(0)
    } else if (!is.null(cur_var)) {
      # Value-label line:  CODE 'Label'  or  CODE "Label"
      m <- regexpr("^(\\S+)\\s+['\"](.+?)['\"]", lt, perl = TRUE)
      if (m > 0) {
        cs <- attr(m, "capture.start"); cl <- attr(m, "capture.length")
        code  <- substr(lt, cs[1], cs[1] + cl[1] - 1)
        label <- substr(lt, cs[2], cs[2] + cl[2] - 1)
        result[[cur_var]] <- c(result[[cur_var]], setNames(label, code))
      }
    }
  }
  result
}

# Read a CSPro fixed-width .dat file using its companion .sps (preferred) or .dcf layout.
# Inputs: dat_path -- path to .dat file; dcf_path -- optional fallback if no .sps found
# Returns: plain data frame with character columns (whitespace trimmed).
# Note: SPS positions are absolute column offsets; DCF positions are record-relative.
# Returns a plain data frame with character columns (whitespace trimmed).
read_cspr <- function(dat_path, dcf_path = NULL) {
  sps_path <- sub("(?i)\\.dat$", ".SPS", dat_path, perl = TRUE)
  
  if (file.exists(sps_path)) {
    layout  <- parse_sps_layout(sps_path)   # returns Name, Start, Len, ImplicitDec
    col_pos <- readr::fwf_positions(
      start     = layout$Start,
      end       = layout$Start + layout$Len - 1L,
      col_names = layout$Name
    )
  } else {
    if (is.null(dcf_path))
      dcf_path <- sub("(?i)\\.dat$", ".dcf", dat_path, perl = TRUE)
    if (!file.exists(dcf_path))
      stop("No .sps or .dcf layout file found for: ", dat_path)
    layout  <- parse_dcf_layout(dcf_path)   # returns Name, Start, Len
    col_pos <- readr::fwf_widths(layout$Len, layout$Name)
  }
  
  df <- readr::read_fwf(
    dat_path,
    col_positions  = col_pos,
    col_types      = readr::cols(.default = readr::col_character()),
    progress       = FALSE,
    show_col_types = FALSE
  )
  mutate(df, across(where(is.character), trimws))
}
# Load a survey data file, dispatching on extension (.dta, .sav, .dat).
# Inputs: path -- full path to the survey data file
# Returns: plain data frame; haven_labelled columns are preserved for code_positivity().
load_survey_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "dta") {
    df <- haven::read_dta(path)
  } else if (ext == "sav") {
    df <- haven::read_sav(path)
  } else if (ext == "dat") {
    return(as.data.frame(read_cspr(path)))
  } else {
    stop("Unsupported file extension: .", ext, "  (supported: .dta, .sav, .dat+.dcf)")
  }
  # Keep haven_labelled columns as-is. code_positivity() resolves patterns
  # against the label map, so it handles both numeric patterns (e.g. "1"/"0"
  # used in older recode CSVs) and text-label patterns (e.g. "positive"/"negative"
  # used in newer recode CSVs) without any changes to the recode CSVs.
  as.data.frame(df)
}

# Construct a full file path from dta_dir, folder_dir, and filename.
# Returns: file.path string; omits folder_dir if it is NA or empty.
build_file_path <- function(dta_dir, folder_dir, filename) {
  if (!is.na(folder_dir) && nchar(trimws(folder_dir)) > 0)
    file.path(dta_dir, folder_dir, filename)
  else
    file.path(dta_dir, filename)
}

# Format date columns as mm/dd/yyyy strings before write.csv, for compatibility
# with downstream consumers (e.g. 3_sample_itn_mass_distribution_coverages.R).
# Handles Date/POSIXt objects and numeric values (dplyr mean() on Date columns
# returns numeric days-since-epoch, which also needs converting).
# Inputs: df -- data frame whose column names contain 'date'
# Returns: df with date columns as character "mm/dd/yyyy" strings.
format_dates_for_csv <- function(df) {
  date_cols <- grep("date", colnames(df), value = TRUE, ignore.case = TRUE)
  for (col in date_cols) {
    x <- df[[col]]
    if (inherits(x, "Date") || inherits(x, "POSIXt")) {
      df[[col]] <- format(as.Date(x), "%m/%d/%Y")
    } else if (is.numeric(x)) {
      # dplyr mean() on a Date column yields numeric days since 1970-01-01
      df[[col]] <- format(as.Date(x, origin = "1970-01-01"), "%m/%d/%Y")
    }
  }
  df
}

# Parse a mean_date column read back from CSV.
# Handles both the current mm/dd/yyyy string format and legacy numeric
# days-since-1970-01-01 values, so scripts work with both old and new CSVs.
.parse_csv_date <- function(x) {
  if (is.numeric(x)) {
    as.Date(x, origin = "1970-01-01")
  } else {
    as.Date(as.character(x), tryFormats = c("%m/%d/%Y", "%m-%d-%Y", "%Y-%m-%d"))
  }
}


# =============================================================================
# SECTION 2: ADMIN NAME & ARCHETYPE HELPERS
# Normalise admin/state name strings for fuzzy matching, rename archetype CSV
# columns to canonical names, and expand admin_sums to include all admins.
# =============================================================================

# Standardise an admin (LGA) name: uppercase, spaces/slashes replaced with hyphens.
# Inputs: name -- raw admin name string
# Returns: normalised character string for use with match/merge operations.
normalize_admin_name <- function(name) {
  name <- toupper(name)
  name <- gsub(" ", "-", name, fixed = TRUE)
  name <- gsub("/", "-", name, fixed = TRUE)
  name
}

# Rename common column-name variants in an archetype/pop CSV to canonical names.
# Canonical targets: NOMDEP (LGA), NOMREGION (state), ZONE (geopolitical zone).
# Returns: df with standardised column names.
normalize_archetype_info <- function(df) {
  rename_once <- function(df, candidates, target) {
    for (c in candidates) {
      if (c %in% colnames(df) && !target %in% colnames(df)) {
        colnames(df)[colnames(df) == c] <- target; break
      }
    }
    df
  }
  df <- rename_once(df, c("LGA","DS","admin_name"),          "NOMDEP")
  df <- rename_once(df, c("State"),                           "NOMREGION")
  df <- rename_once(df, c("Geopolitical.zone","Geopolitical.Zone","Zone","ZONE"), "ZONE")
  df
}

# Build a deduplicated NOMREGION → grouping_col lookup table (one row per state).
# If archetype_info has multiple grouping values for one state (e.g. multiple
# archetypes per state), merging directly inflates cluster_df row counts, causing
# inflated arch/national sums.  This helper keeps only the first occurrence per
# state and emits a warning when inconsistencies are found.
.safe_zone_map <- function(archetype_info, grouping_col) {
  zm <- dplyr::distinct(archetype_info[, c("NOMREGION", grouping_col)])
  dups <- duplicated(zm[["NOMREGION"]])
  if (any(dups)) {
    warning("[grouping] Multiple '", grouping_col, "' values per state in archetype_info — ",
            "keeping first per state. States affected: ",
            paste(unique(zm$NOMREGION[dups]), collapse = ", "))
    zm <- zm[!dups, , drop = FALSE]
  }
  zm
}

# Extract and scale the survey weight column (hv005, v005, or QWEIGHT) to fractions.
# Returns: numeric vector of weights divided by 1e6, or NULL if no weight column found.
get_survey_weight <- function(df) {
  # DHS/MIS Stata files: hv005 / v005 store weight as integer × 1e6
  dhs_col <- intersect(c("hv005", "v005"), colnames(df))[1]
  if (!is.na(dhs_col)) return(as.numeric(df[[dhs_col]]) / 1e6)
  # CSPro MIS files: QWEIGHT stored as integer × 1e6 (6 implicit decimal places)
  cspr_col <- intersect(c("QWEIGHT", "qweight"), colnames(df))[1]
  if (!is.na(cspr_col)) return(as.numeric(df[[cspr_col]]) / 1e6)
  NULL
}

# Derive a per-person ITN-use column from occupant line-number suffix columns.
# Used when {code} is absent but {code}A / {code}B / ... columns are present.
# Inputs: df; code -- target column name; line_col -- respondent line number;
#   pos_val/neg_val -- values for match/no-match.
# Returns: df with {code} column added (or unchanged if already present).
derive_itn_col <- function(df, code, line_col, pos_val, neg_val) {
  if (code %in% colnames(df)) return(df)
  
  # Find all present suffix columns: A, B, C, D, ... (up to Z)
  suffix_cols <- paste0(code, LETTERS)
  suffix_cols <- suffix_cols[suffix_cols %in% colnames(df)]
  
  if (length(suffix_cols) == 0) {
    warning("ITN column '", code, "' not found and no '", code,
            "A/B/C/D' suffix columns are present — ITN positivity will be NA.")
    df[[code]] <- NA_character_
    return(df)
  }
  if (!line_col %in% colnames(df))
    stop("household_line_number_code column '", line_col, "' not found in ITN data.")
  
  message("  ITN: '", code, "' absent — deriving from ",
          paste(suffix_cols, collapse = "/"), " matched against '", line_col, "'.")
  
  # Integer comparison is safer than character for haven_labelled columns
  line_vals <- suppressWarnings(as.integer(df[[line_col]]))
  occ_match <- Reduce(`|`, lapply(suffix_cols, function(sc) {
    occ <- suppressWarnings(as.integer(df[[sc]]))
    !is.na(line_vals) & !is.na(occ) & occ == line_vals
  }))
  df[[code]] <- ifelse(occ_match, as.character(pos_val), as.character(neg_val))
  df
}


# =============================================================================
# SECTION 3: SURVEY WEIGHT & ITN HELPERS
# Multi-file ITN calculation and per-age-group aggregation for CSPro MIS surveys.
# =============================================================================

# Compute ITN use rates from separate PERSONS.DAT and NETS.DAT files (CSPro MIS).
# Joins unique persons to net-occupant records; aggregates by cluster and age group.
# Inputs: persons_df/nets_df -- loaded data frames; *_col -- column name arguments;
#   age_groups -- named list of c(lo, hi) age bounds in years.
# Returns: named list (one entry per age group) of cluster-level rate data frames.
compute_itn_from_persons_and_nets <- function(
    persons_df, nets_df,
    person_cluster_col, person_hh_col, person_line_col, person_age_col,
    net_cluster_col,    net_hh_col,
    age_groups
) {
  # ---- 1. Unique person roster from PERSONS.DAT ----
  # PERSONS.DAT has one row per (person × net) because NETS is nested inside
  # PERSONS in the CSPro questionnaire structure. Deduplicate by (cluster, hh,
  # person line number) so each individual is counted exactly once.
  wt_col <- intersect(c("QHWEIGHT", "itn_weights", "hv005", "v005"), colnames(persons_df))[1]
  select_cols <- c(person_cluster_col, person_hh_col, person_line_col, person_age_col)
  if (!is.na(wt_col)) select_cols <- c(select_cols, wt_col)
  
  unique_persons <- persons_df %>%
    dplyr::select(all_of(select_cols)) %>%
    dplyr::mutate(
      .cluster = suppressWarnings(as.integer(.data[[person_cluster_col]])),
      .hh      = suppressWarnings(as.integer(.data[[person_hh_col]])),
      .line    = suppressWarnings(as.integer(.data[[person_line_col]]))
    ) %>%
    dplyr::distinct(.cluster, .hh, .line, .keep_all = TRUE)
  
  message("  ITN multifile: ", nrow(unique_persons), " unique persons from PERSONS.DAT")
  
  # ---- 2. Net occupants from NETS.DAT ----
  # Find all present QH129A/B/C/D columns in nets_df
  occ_cols <- paste0("QH129", LETTERS)
  occ_cols <- occ_cols[occ_cols %in% colnames(nets_df)]
  if (length(occ_cols) == 0)
    stop("compute_itn_from_persons_and_nets: no QH129A/B/C/D columns found in NETS file.")
  
  message("  ITN multifile: occupant columns found in NETS.DAT: ",
          paste(occ_cols, collapse = ", "))
  
  # Pivot occupant columns to long: one row per (cluster, hh, person_line)
  # who slept under any net. Do this without tidyr for minimal dependencies.
  net_sleepers <- do.call(rbind, lapply(occ_cols, function(col) {
    tmp <- data.frame(
      cluster = suppressWarnings(as.integer(nets_df[[net_cluster_col]])),
      hh      = suppressWarnings(as.integer(nets_df[[net_hh_col]])),
      line    = suppressWarnings(as.integer(nets_df[[col]])),
      stringsAsFactors = FALSE
    )
    tmp
  })) %>%
    dplyr::filter(!is.na(line), line > 0) %>%
    dplyr::distinct(cluster, hh, line) %>%
    dplyr::mutate(slept_under_net = 1L)
  
  message("  ITN multifile: ", nrow(net_sleepers),
          " unique (cluster, hh, person) records in NETS.DAT occupant columns")
  
  # ---- 3. Join persons to net-sleepers ----
  persons_itn <- unique_persons %>%
    dplyr::left_join(net_sleepers,
                     by = c(".cluster" = "cluster", ".hh" = "hh", ".line" = "line")) %>%
    dplyr::mutate(
      pos     = ifelse(is.na(slept_under_net), 0L, 1L),
      age_num = parse_dhs_age(.data[[person_age_col]], person_age_col)
    )
  
  message("  ITN multifile: ", sum(persons_itn$pos, na.rm = TRUE), " / ",
          nrow(persons_itn), " persons slept under a net (all ages, all clusters)")
  
  # Survey weight: QHWEIGHT in CSPro is stored as integer × 1e6
  has_weight <- !is.na(wt_col) && wt_col %in% colnames(persons_itn)
  if (has_weight) {
    raw_wt <- suppressWarnings(as.numeric(persons_itn[[wt_col]]))
    persons_itn$itn_weights <- ifelse(raw_wt > 1000, raw_wt / 1e6, raw_wt)
  }
  
  # ---- 4. Aggregate per cluster per age group ----
  results <- list()
  for (grp in names(age_groups)) {
    lo <- age_groups[[grp]][1];  hi <- age_groups[[grp]][2]
    if (is.infinite(lo) && is.infinite(hi)) {
      sub <- persons_itn
    } else {
      sub <- persons_itn[!is.na(persons_itn$age_num) &
                           persons_itn$age_num >  lo &
                           persons_itn$age_num <= hi, ]
    }
    if (has_weight && grp == "itn_all") {
      agg <- sub %>%
        dplyr::filter(!is.na(pos)) %>%
        dplyr::group_by(across(all_of(person_cluster_col))) %>%
        dplyr::summarise(
          rate       = mean(pos, na.rm = TRUE),
          num_pos    = sum(pos, na.rm = TRUE),
          num_tested = dplyr::n(),
          itn_weights = mean(itn_weights, na.rm = TRUE),
          .groups    = "drop"
        )
    } else {
      agg <- sub %>%
        dplyr::filter(!is.na(pos)) %>%
        dplyr::group_by(across(all_of(person_cluster_col))) %>%
        dplyr::summarise(
          rate       = mean(pos, na.rm = TRUE),
          num_pos    = sum(pos, na.rm = TRUE),
          num_tested = dplyr::n(),
          .groups    = "drop"
        )
    }
    results[[grp]] <- agg
  }
  results
}


# =============================================================================
# SECTION 4: POSITIVITY CODING & AGE PARSING
# =============================================================================

# Convert a DHS household-member age column to numeric, applying DHS special codes.
# Inputs: age_char_vec -- raw age column; age_col_name -- for warning messages
# Returns: numeric vector; codes 96/97/99+ -> NA with message; 98 (DK) -> NA silently.
parse_dhs_age <- function(age_char_vec, age_col_name) {
  if (haven::is.labelled(age_char_vec)) {
    # Numeric codes stored directly; no gsub needed. Using as.numeric() preserves
    # the raw DHS codes (96/97/98/99+) so the validity checks below fire correctly.
    raw <- suppressWarnings(as.numeric(age_char_vec))
  } else {
    raw <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", age_char_vec)))
  }
  raw[!is.na(raw) & raw == 98] <- NA   # don't know — silent
  invalid <- !is.na(raw) & (raw %in% c(96, 97) | raw >= 99)
  if (any(invalid)) {
    bad_codes <- sort(unique(raw[invalid]))
    message(sprintf(
      "  [WARNING] Age column '%s': %d record(s) with invalid age code(s) [%s] treated as NA.",
      age_col_name, sum(invalid), paste(bad_codes, collapse = ", ")
    ))
    raw[invalid] <- NA
  }
  raw
}

# Code a binary positivity indicator 'pos' onto a data frame.
# pos = 1 for pos_pattern + alternate_pos_patterns; 0 for neg_pattern; NA otherwise.
# Handles haven_labelled (label lookup), numeric, and character column types.
# Inputs: df; var_code -- column name; pos_pattern/neg_pattern; alternate_pos_patterns.
# Returns: df with pos column added.
code_positivity <- function(df, var_code, pos_pattern, neg_pattern,
                            alternate_pos_patterns = character(0)) {
  if (!var_code %in% colnames(df)) {
    warning("Variable '", var_code, "' not found in data — pos set to NA.")
    df$pos <- NA_real_
    return(df)
  }
  col <- df[[var_code]]
  
  if (haven::is.labelled(col)) {
    # Resolve a pattern string to the numeric code stored in the column.
    # Tries label-name lookup first (e.g. "positive" → 1), then numeric coercion.
    labs <- attr(col, "labels")   # named numeric vector: name = label text, value = code
    resolve <- function(p) {
      m <- match(as.character(p), names(labs))
      if (!is.na(m)) return(labs[[m]])
      suppressWarnings(as.numeric(p))   # numeric pattern (e.g. "1") used directly
    }
    cmp <- as.numeric(col)
    match1 <- function(p) { rv <- resolve(p); !is.na(cmp) & !is.na(rv) & cmp == rv }
    df$pos <- NA_real_
    df$pos[match1(pos_pattern)] <- 1
    df$pos[match1(neg_pattern)] <- 0
    for (alt in alternate_pos_patterns) df$pos[match1(alt)] <- 1
  } else {
    coerce <- function(p) if (is.numeric(col)) suppressWarnings(as.numeric(p)) else as.character(p)
    df$pos <- NA_real_
    df$pos[col == coerce(pos_pattern)] <- 1
    df$pos[col == coerce(neg_pattern)] <- 0
    for (alt in alternate_pos_patterns) df$pos[col == coerce(alt)] <- 1
  }
  df
}


# =============================================================================
# SECTION 5: CLUSTER-LEVEL AGGREGATION
# Aggregate individual-level positivity to cluster-level rates and counts.
# merge_cluster_result joins results back into the running cluster data frame.
# =============================================================================

# Aggregate coded 'pos' column to cluster level (unweighted or ITN-weighted).
# Returns: cluster_id_col | rate | num_pos | num_tested [| itn_weights].
aggregate_to_clusters <- function(df, cluster_id_col, include_weight = FALSE) {
  if (include_weight) {
    wt <- get_survey_weight(df)
    df$itn_weights <- if (!is.null(wt)) wt else rep(1, nrow(df))
    df %>%
      filter(!is.na(pos)) %>%
      group_by(across(all_of(cluster_id_col))) %>%
      summarise(rate = mean(pos, na.rm = TRUE), num_pos = sum(pos, na.rm = TRUE),
                num_tested = n(), itn_weights = mean(itn_weights, na.rm = TRUE),
                .groups = "drop")
  } else {
    df %>%
      filter(!is.na(pos)) %>%
      group_by(across(all_of(cluster_id_col))) %>%
      summarise(rate = mean(pos, na.rm = TRUE), num_pos = sum(pos, na.rm = TRUE),
                num_tested = n(), .groups = "drop")
  }
}

# Aggregate 'pos' to cluster x month x year level for monthly PfPR reference.
# Returns: cluster_id_col | month_col | year_col | rate | num_pos | num_tested.
aggregate_to_clusters_monthly <- function(df, cluster_id_col, month_col, year_col) {
  df %>%
    filter(!is.na(pos)) %>%
    group_by(across(all_of(c(cluster_id_col, month_col, year_col)))) %>%
    summarise(rate = mean(pos, na.rm = TRUE), num_pos = sum(pos, na.rm = TRUE),
              num_tested = n(), .groups = "drop")
}

# Vaccine aggregation for DHS: compute child age from survey/birth dates,
# exclude children below min_age_months, then aggregate to cluster level.
aggregate_to_clusters_vacc_dhs <- function(df, cluster_id_col,
                                           survey_month_col, survey_year_col,
                                           birth_month_col,  birth_year_col,
                                           min_age_months) {
  df$survey_date <- as.Date(paste0(as.integer(df[[survey_year_col]]), "-",
                                   as.integer(df[[survey_month_col]]), "-28"))
  df$birth_date  <- as.Date(paste0(as.integer(df[[birth_year_col]]),  "-",
                                   as.integer(df[[birth_month_col]]),  "-01"))
  df$age_days    <- as.numeric(df$survey_date - df$birth_date)
  df$pos[df$age_days <= min_age_months * 30.4] <- NA
  df %>%
    filter(!is.na(pos)) %>%
    group_by(across(all_of(cluster_id_col))) %>%
    summarise(rate = mean(pos, na.rm = TRUE), num_pos = sum(pos, na.rm = TRUE),
              num_tested = n(), .groups = "drop")
}

# Vaccine aggregation for MICS: age column already available; filter then aggregate.
# Returns: cluster_id_col | rate | num_pos | num_tested.
aggregate_to_clusters_vacc_mics <- function(df, cluster_id_col, age_col, min_age_months) {
  df$pos[df[[age_col]] <= min_age_months * 30.4] <- NA
  df %>%
    filter(!is.na(pos)) %>%
    group_by(across(all_of(cluster_id_col))) %>%
    summarise(rate = mean(pos, na.rm = TRUE), num_pos = sum(pos, na.rm = TRUE),
              num_tested = n(), .groups = "drop")
}

# Merge a single-variable cluster result into the running cluster data frame.
# Renames rate/num_pos/num_tested to {var_name}_rate/_num_true/_num_total.
merge_cluster_result <- function(MIS_outputs, cluster_result, cluster_id_col, var_name) {
  merged <- merge(MIS_outputs, cluster_result,
                  by.x = "clusterid", by.y = cluster_id_col, all = TRUE)
  merged <- dplyr::rename(merged,
                          !!paste0(var_name, "_rate")      := rate,
                          !!paste0(var_name, "_num_true")  := num_pos,
                          !!paste0(var_name, "_num_total") := num_tested
  )
  merged
}

# Resolve a column name that may use CSPro occurrence suffixes ({col}$1, {col}$2, ...).
# Returns: the resolved column name, or NULL if neither bare name nor {col}$1 is found.
# Returns NULL if neither the bare name nor $1 is found.
resolve_occurrence_col <- function(df, col) {
  if (is.na(col) || !nzchar(trimws(col))) return(NULL)
  if (col %in% colnames(df)) return(col)
  occ1 <- paste0(col, "$1")
  if (occ1 %in% colnames(df)) return(occ1)
  NULL
}


# Compute mean interview date per cluster from month/year columns.
# Handles haven_labelled months and month-name strings (e.g. 'January').
# Returns: data frame with cluster_id_col | mean_date, or NULL if date codes are NA.
extract_cluster_dates <- function(df, cluster_id_col, month_col, year_col) {
  if (is.na(month_col) || !nzchar(trimws(month_col)) ||
      is.na(year_col)  || !nzchar(trimws(year_col)))  return(NULL)
  month_col <- resolve_occurrence_col(df, month_col)
  year_col  <- resolve_occurrence_col(df, year_col)
  if (is.null(month_col) || is.null(year_col)) return(NULL)
  month_vals <- df[[month_col]]
  # Handle haven_labelled month columns (e.g. labelled 1–12 or with month names).
  # as.integer() strips the label map and returns the underlying numeric code.
  if (haven::is.labelled(month_vals)) {
    month_vals <- as.integer(month_vals)
  } else if (is.character(month_vals) && any(month_vals %in% month.name, na.rm = TRUE)) {
    month_vals <- match(month_vals, month.name)
  }
  year_vals <- if (haven::is.labelled(df[[year_col]])) as.integer(df[[year_col]]) else df[[year_col]]
  df$date <- as.Date(paste0(year_vals, "-", month_vals, "-01"),
                     tryFormats = c("%Y-%m-%d", "%y-%m-%d"))
  df %>%
    filter(!is.na(date)) %>%
    group_by(across(all_of(cluster_id_col))) %>%
    summarise(mean_date = mean(date, na.rm = TRUE), .groups = "drop")
}


# =============================================================================
# SECTION 6: TREATMENT & ACT CODING
# Build per-person treatment-seeking and ACT indicator columns from raw survey data.
# Supports both classic per-column DHS format and CSPro concatenated-string format.
# =============================================================================

# Positive-response values used for any yes/no treatment column.
# Internal yes/no value sets and row-level matcher used by build_treatment_columns.
.YES_VALS <- c("yes","Yes","YES","Y","1",1L,"T","TRUE","True")
.ANY_VALS <- c(.YES_VALS, "no","No","NO","N","0",0L,"F","FALSE","False")
.any_match <- function(row, vals) any(row %in% vals, na.rm = TRUE)

# Add treatment-seeking columns to an individual-level data frame.
# Produces: sought_treatment, public_treatment, private_formal_treatment,
#   private_informal_treatment.
# Two modes: classic DHS (h32x columns) and CSPro MIS (Q408 concatenated string).
# Inputs: df; year -- survey year (affects private-informal column mapping).
build_treatment_columns <- function(df, year,
                                    codes_public           = paste0("h32", c("a","b","c","d","e","f","g","h","i")),
                                    codes_private_formal   = paste0("h32", c("j","k","l","m","n","o","p","na","nb")),
                                    codes_private_informal = NULL,
                                    treat_concat_col         = "Q408",
                                    treat_sought_col         = "Q407",   # yes/no sought-treatment indicator (asked of febrile children)
                                    public_letters           = c("A","B","C","D","E","F","U"),
                                    private_formal_letters   = c("G","H","I","J","K","L","M","N","O"),
                                    private_informal_letters = c("P","Q","S","T","X")) {
  
  # --- CSPro concatenated-string mode ---
  # Triggered when no h32x columns are present but treat_concat_col exists.
  #
  # Denominator: febrile children — those for whom treat_sought_col (Q407) is
  # non-NA (the question is only asked when the child had a fever).
  # sought_treatment numerator: Q407 == 1 (yes, sought advice/treatment).
  # Sector columns (public/private): Q408 is only filled for children who sought
  # treatment (Q407=1); for all other febrile children it is blank.  Sector flags
  # are therefore 0 for febrile children who did not seek treatment and 1 for
  # those who sought treatment at that type of facility.
  has_h32 <- any(intersect(c(codes_public, codes_private_formal), colnames(df)) != character(0))
  if (!has_h32 && !is.null(treat_concat_col) && treat_concat_col %in% colnames(df)) {
    s <- trimws(df[[treat_concat_col]])
    
    # Febrile children: Q407 is 1 (yes) or 2 (no); non-febrile children have NA
    if (!is.null(treat_sought_col) && treat_sought_col %in% colnames(df)) {
      q407      <- suppressWarnings(as.integer(df[[treat_sought_col]]))
      has_fever <- !is.na(q407)          # TRUE = febrile child (in denominator)
      sought    <- has_fever & q407 == 1 # TRUE = sought treatment
    } else {
      # Fallback if Q407 absent: use Q408 non-empty as the sought flag (original behaviour)
      message("  build_treatment_columns: '", treat_sought_col,
              "' not found — falling back to Q408 non-empty for sought_treatment denominator.")
      has_fever <- nchar(s) > 0 & !is.na(s)
      sought    <- has_fever
    }
    
    has_pub  <- grepl(paste0("[", paste(public_letters,           collapse=""), "]"), s)
    has_prf  <- grepl(paste0("[", paste(private_formal_letters,   collapse=""), "]"), s)
    has_pri  <- grepl(paste0("[", paste(private_informal_letters, collapse=""), "]"), s)
    
    # All four columns share the same denominator: febrile children (has_fever).
    # Sector columns are 1 if the child sought treatment at that sector, 0 otherwise.
    df$sought_treatment           <- ifelse(has_fever, as.numeric(sought),              NA_real_)
    df$public_treatment           <- ifelse(has_fever, as.numeric(sought & has_pub),    NA_real_)
    df$private_formal_treatment   <- ifelse(has_fever, as.numeric(sought & has_prf),    NA_real_)
    df$private_informal_treatment <- ifelse(has_fever, as.numeric(sought & has_pri),    NA_real_)
    return(df)
  }
  
  # --- Classic per-column mode ---
  if (is.null(codes_private_informal)) {
    # 2013: h32v = traditional (excluded), h32t = private informal (included)
    # non-2013: h32t = traditional (excluded), h32v = private informal (included)
    codes_private_informal <- if (year == 2013)
      paste0("h32", c("q","r","s","u","t","w","x","nc","nd","ne"))
    else
      paste0("h32", c("q","r","s","u","v","w","x","nc","nd","ne"))
  }
  codes_all <- intersect(
    c(codes_public, codes_private_formal, codes_private_informal), colnames(df))
  codes_public           <- intersect(codes_public,           colnames(df))
  codes_private_formal   <- intersect(codes_private_formal,   colnames(df))
  codes_private_informal <- intersect(codes_private_informal, colnames(df))
  
  # Convert any haven_labelled h32* columns to character label text before .YES_VALS
  # comparison.  as_factor() converts via the label map first.
  for (cc in codes_all) {
    if (haven::is.labelled(df[[cc]]))
      df[[cc]] <- as.character(haven::as_factor(df[[cc]]))
  }
  
  if (length(codes_all) > 1) {
    df$sought_treatment           <- as.numeric(apply(df[, codes_all,              drop=FALSE], 1, .any_match, vals=.YES_VALS))
    df$responded_sought_treatment <- as.numeric(apply(df[, codes_all,              drop=FALSE], 1, .any_match, vals=.ANY_VALS))
    df$public_treatment           <- as.numeric(apply(df[, codes_public,           drop=FALSE], 1, .any_match, vals=.YES_VALS))
    df$private_formal_treatment   <- as.numeric(apply(df[, codes_private_formal,   drop=FALSE], 1, .any_match, vals=.YES_VALS))
    df$private_informal_treatment <- as.numeric(apply(df[, codes_private_informal, drop=FALSE], 1, .any_match, vals=.YES_VALS))
    # Mask non-responders
    no_resp <- df$responded_sought_treatment != 1
    df$sought_treatment[no_resp]           <- NA
    df$public_treatment[no_resp]           <- NA
    df$private_formal_treatment[no_resp]   <- NA
    df$private_informal_treatment[no_resp] <- NA
  } else {
    df$sought_treatment <- df$public_treatment <-
      df$private_formal_treatment <- df$private_informal_treatment <- NA_real_
  }
  df
}

# Add ACT indicator columns; requires build_treatment_columns() to have been called.
# Produces: art_given_antimalarial, act_given_antimal_{sector}, act_given_treat_{sector}.
# Two modes: classic per-drug columns (art_codes = character vector) or CSPro
#   concatenated string (art_codes = list with $drug_col and $act_letters).
build_act_columns <- function(df,
                              art_codes     = c("ml13e"),
                              non_art_codes = c("ml13a","ml13b","ml13c","ml13d","ml13aa","ml13da","ml13h")) {
  
  # --- CSPro concatenated-string mode ---
  # Detected when art_codes is a list containing a $drug_col element.
  if (is.list(art_codes) && "drug_col" %in% names(art_codes)) {
    drug_col    <- art_codes$drug_col
    act_letters <- art_codes$act_letters %||% c("A")
    non_art_letters <- if (is.list(non_art_codes) && "non_art_letters" %in% names(non_art_codes))
      non_art_codes$non_art_letters
    else c("B","C","D","E","F","G","H","I")
    if (!drug_col %in% colnames(df)) {
      warning("Drug concatenation column '", drug_col, "' not found; ACT columns set to NA.")
      df$art_given_antimalarial <- NA_real_
      return(df)
    }
    s <- trimws(df[[drug_col]])
    all_letters <- c(act_letters, non_art_letters)
    has_act     <- grepl(paste0("[", paste(act_letters,  collapse = ""), "]"), s)
    has_antimal <- grepl(paste0("[", paste(all_letters,  collapse = ""), "]"), s)
    df$any_antimalarial       <- ifelse(has_antimal, 1, NA_real_)
    df$art_antimalarial       <- as.numeric(has_act)
    df$art_given_antimalarial <- df$art_antimalarial / df$any_antimalarial
    df <- df %>% mutate(
      act_given_antimal_public           = if_else(public_treatment          == 1, art_given_antimalarial, NA_real_),
      act_given_antimal_private_formal   = if_else(private_formal_treatment  == 1, art_given_antimalarial, NA_real_),
      act_given_antimal_private_informal = if_else(private_informal_treatment == 1, art_given_antimalarial, NA_real_),
      act_given_treat_public             = if_else(public_treatment          == 1, art_antimalarial, NA_real_),
      act_given_treat_private_formal     = if_else(private_formal_treatment  == 1, art_antimalarial, NA_real_),
      act_given_treat_private_informal   = if_else(private_informal_treatment == 1, art_antimalarial, NA_real_)
    )
    return(df)
  }
  
  # --- Classic per-drug column mode ---
  art_codes     <- intersect(art_codes,     colnames(df))
  non_art_codes <- intersect(non_art_codes, colnames(df))
  all_antimal   <- c(art_codes, non_art_codes)
  
  # Same haven_labelled → label text conversion as in build_treatment_columns.
  for (cc in all_antimal) {
    if (haven::is.labelled(df[[cc]]))
      df[[cc]] <- as.character(haven::as_factor(df[[cc]]))
  }
  
  if (length(all_antimal) == 0) {
    warning("No antimalarial columns found; ACT columns set to NA.")
    df$art_given_antimalarial <- NA_real_
    return(df)
  }
  
  df$any_antimalarial <- as.numeric(apply(df[, all_antimal, drop=FALSE], 1, .any_match, vals=.YES_VALS))
  df$any_antimalarial[df$any_antimalarial != 1] <- NA  # denominator = those receiving any antimalarial
  
  df$art_antimalarial <- if (length(art_codes) > 1)
    as.numeric(apply(df[, art_codes, drop=FALSE], 1, .any_match, vals=.YES_VALS))
  else if (length(art_codes) == 1)
    as.numeric(df[[art_codes]] %in% .YES_VALS)
  else NA_real_
  
  df$art_given_antimalarial <- df$art_antimalarial / df$any_antimalarial
  
  df <- df %>% mutate(
    act_given_antimal_public           = if_else(public_treatment          == 1, art_given_antimalarial, NA_real_),
    act_given_antimal_private_formal   = if_else(private_formal_treatment  == 1, art_given_antimalarial, NA_real_),
    act_given_antimal_private_informal = if_else(private_informal_treatment == 1, art_given_antimalarial, NA_real_),
    act_given_treat_public             = if_else(public_treatment          == 1, art_antimalarial, NA_real_),
    act_given_treat_private_formal     = if_else(private_formal_treatment  == 1, art_antimalarial, NA_real_),
    act_given_treat_private_informal   = if_else(private_informal_treatment == 1, art_antimalarial, NA_real_)
  )
  df
}


# =============================================================================
# SECTION 7: GPS & SPATIAL JOIN
# Assign survey clusters to administrative units via point-in-polygon spatial join.
# extract_gps_from_nets builds cluster GPS from household-level CSPro .dat files,
# with outlier removal, state-consistency filtering, and outside-boundary snapping.
# =============================================================================

# Assign cluster GPS points to admin units via point-in-polygon join on admin_shape.
# Inputs: cluster_df with latitude/longitude; admin_shape -- SpatialPolygonsDataFrame
# Returns: cluster_df with NOMREGION and NOMDEP columns added.
# Uses sf::st_join to avoid sp/raster/terra CRS chain which can fail when GDAL state
# has not been primed by a prior terra::vect() or sf::st_as_sf() call in the session.
assign_clusters_to_admins <- function(cluster_df, admin_shape,
                                      region_col = "NOMREGION",
                                      admin_col  = "NOMDEP") {
  has_coords <- !is.na(cluster_df$latitude) & !is.na(cluster_df$longitude)
  n_drop <- sum(!has_coords)
  if (n_drop > 0) {
    warning("assign_clusters_to_admins: dropping ", n_drop,
            " cluster(s) with NA coordinates before spatial join.")
    cluster_df <- cluster_df[has_coords, , drop = FALSE]
  }
  admin_sf <- sf::st_as_sf(admin_shape)
  if (is.na(sf::st_crs(admin_sf)) ||
      !identical(sf::st_crs(admin_sf), sf::st_crs(4326)))
    admin_sf <- sf::st_transform(admin_sf, 4326)
  # Drop any pre-existing admin-name columns to prevent sf::st_join from adding
  # .x/.y suffixes that would break the column lookup below.
  for (.col in unique(c(region_col, admin_col, "NOMREGION", "NOMDEP"))) {
    if (.col %in% colnames(cluster_df)) cluster_df[[.col]] <- NULL
  }
  pts_sf <- sf::st_as_sf(cluster_df,
                         coords = c("longitude", "latitude"),
                         crs    = 4326,
                         remove = FALSE)
  # Tag each input row so we can detect polygon-overlap duplicates after the join.
  pts_sf$.row_id. <- seq_len(nrow(pts_sf))
  joined_full <- suppressMessages(
    sf::st_join(pts_sf, admin_sf[, c(region_col, admin_col)])
  )
  joined_df <- as.data.frame(joined_full)
  joined_df$geometry <- NULL

  # ---- Detect ambiguous clusters (matched >1 polygon) ----
  row_counts    <- table(joined_df$.row_id.)
  ambiguous_ids <- as.integer(names(row_counts[row_counts > 1]))
  ambiguous_log <- NULL
  if (length(ambiguous_ids) > 0) {
    message("  [admin_assign] ", length(ambiguous_ids),
            " cluster(s) matched more than one admin polygon — keeping first match.")
    clust_col <- if ("clusterid" %in% colnames(joined_df)) "clusterid" else NULL
    amb_rows  <- joined_df[joined_df$.row_id. %in% ambiguous_ids, ]
    ambiguous_log <- do.call(rbind, lapply(ambiguous_ids, function(rid) {
      rows <- amb_rows[amb_rows$.row_id. == rid, ]
      data.frame(
        latitude            = rows$latitude[1],
        longitude           = rows$longitude[1],
        clusterid           = if (!is.null(clust_col)) rows[[clust_col]][1] else NA_integer_,
        n_matches           = nrow(rows),
        all_matched_admins  = paste(rows[[admin_col]],  collapse = " / "),
        all_matched_regions = paste(rows[[region_col]], collapse = " / "),
        selected_admin      = rows[[admin_col]][1],
        selected_region     = rows[[region_col]][1],
        stringsAsFactors    = FALSE
      )
    }))
  }

  # Keep only the first polygon match per input cluster point.
  result          <- joined_df[!duplicated(joined_df$.row_id.), , drop = FALSE]
  result$.row_id. <- NULL
  result$NOMREGION <- result[[region_col]]
  result$NOMDEP    <- result[[admin_col]]
  if (region_col != "NOMREGION") result[[region_col]] <- NULL
  if (admin_col  != "NOMDEP")    result[[admin_col]]  <- NULL

  list(cluster_df = result, ambiguous_log = ambiguous_log)
}

# Map page showing clusters that matched >1 admin polygon during the spatial join.
# Each point is labelled with the selected admin and the alternative(s); a caption
# gives the total count.  Saved as a standalone PDF in the same plots directory as
# the GPS diagnostics PDF.
.plot_ambiguous_admin_assignments <- function(ambiguous_log, admin_shape, plot_file) {
  if (is.null(ambiguous_log) || nrow(ambiguous_log) == 0) return(invisible(NULL))
  for (pkg in c("ggplot2", "sf")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning("Package '", pkg, "' needed for ambiguous-polygon map — skipping.")
      return(invisible(NULL))
    }
  }
  admin_sf <- sf::st_as_sf(admin_shape)
  if (is.na(sf::st_crs(admin_sf)) || !identical(sf::st_crs(admin_sf), sf::st_crs(4326)))
    admin_sf <- sf::st_transform(admin_sf, 4326)

  # Build a short label: "selected / alt" — truncate long multi-match strings
  ambiguous_log$label <- paste0(
    ambiguous_log$selected_admin, "\n[alt: ",
    gsub(paste0("^", ambiguous_log$selected_admin, " / "), "", ambiguous_log$all_matched_admins),
    "]"
  )
  # For cases where selected is not at the start of all_matched_admins, fall back
  ambiguous_log$alt <- mapply(function(all, sel) {
    parts <- strsplit(all, " / ")[[1]]
    alts  <- parts[parts != sel]
    if (length(alts) == 0) "—" else paste(alts, collapse = " / ")
  }, ambiguous_log$all_matched_admins, ambiguous_log$selected_admin, USE.NAMES = FALSE)
  ambiguous_log$label <- paste0(ambiguous_log$selected_admin, "\n[alt: ", ambiguous_log$alt, "]")

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = admin_sf, fill = "grey97", colour = "grey70", linewidth = 0.25) +
    ggplot2::geom_point(
      data = ambiguous_log,
      ggplot2::aes(longitude, latitude, colour = selected_admin),
      size = 3, alpha = 0.85, show.legend = FALSE) +
    ggplot2::geom_text(
      data = ambiguous_log,
      ggplot2::aes(longitude, latitude, label = label),
      size = 2.2, vjust = -0.6, lineheight = 0.85, check_overlap = FALSE) +
    ggplot2::labs(
      title    = "Admin assignment: clusters matched to >1 polygon",
      subtitle = paste0(
        nrow(ambiguous_log), " cluster(s) fell on or near an admin boundary and matched ",
        "multiple polygons.\nThe FIRST matched polygon was used (shown). ",
        "Verify that the selected LGA is correct."),
      x = NULL, y = NULL,
      caption = paste0("Selected admin shown as point colour and label; alternative(s) in brackets.\n",
                       "If incorrect, review the admin shapefile for boundary overlaps.")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "grey50"))

  grDevices::pdf(plot_file, width = 11, height = 8.5)
  print(p)
  grDevices::dev.off()
  message("  [admin_assign] Ambiguous-polygon map saved to ", basename(plot_file))
}

# Great-circle distance in km between two lat/lon pairs (vectorised).
# Internal helper used by .max_spread_km and the outside-Nigeria snap step.
.haversine_km <- function(lat1, lon1, lat2, lon2) {
  R    <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a    <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# Maximum distance (km) from the group centroid to any point in the set.
# Returns 0 when fewer than 2 points are supplied.
.max_spread_km <- function(lats, lons) {
  if (length(lats) < 2) return(0)
  max(.haversine_km(lats, lons, mean(lats, na.rm = TRUE), mean(lons, na.rm = TRUE)),
      na.rm = TRUE)
}

# Produce a one-page arrow map documenting GPS snap corrections.
# Shows clusters moved from outside-Nigeria GPS to nearest LGA centroid.
# Inputs: snap_df -- data frame of snapped clusters; admin_sf -- sf polygon object.
# Outputs: PDF saved to plot_file.
.plot_snap_corrections <- function(snap_df, admin_sf, plot_file) {
  for (pkg in c("ggplot2", "sf")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning("Package '", pkg, "' needed for snap correction plot — skipping.")
      return(invisible(NULL))
    }
  }
  
  pad <- 1.0  # degrees padding around old + new points
  all_lons <- c(snap_df$old_lon, snap_df$new_lon)
  all_lats <- c(snap_df$old_lat, snap_df$new_lat)
  xlim_s <- range(all_lons, na.rm = TRUE) + c(-pad, pad)
  ylim_s <- range(all_lats, na.rm = TRUE) + c(-pad, pad)
  
  lga_ctrd_df <- tryCatch({
    ctrd   <- suppressWarnings(sf::st_centroid(admin_sf))
    coords <- as.data.frame(sf::st_coordinates(ctrd))
    data.frame(NOMDEP = admin_sf$NOMDEP, X = coords$X, Y = coords$Y,
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  
  snap_df$label <- paste0("Cluster ", snap_df$cluster,
                          "\n-> ", snap_df$snap_lga,
                          " (", round(snap_df$border_dist_miles, 1), " mi)")
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = admin_sf, fill = "grey95", colour = "grey65", linewidth = 0.3) +
    {if (!is.null(lga_ctrd_df))
      ggplot2::geom_text(data = lga_ctrd_df,
                         ggplot2::aes(X, Y, label = NOMDEP),
                         size = 2.3, colour = "grey45", check_overlap = TRUE)} +
    ggplot2::geom_segment(
      data = snap_df,
      ggplot2::aes(x = old_lon, y = old_lat, xend = new_lon, yend = new_lat),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed"),
      colour    = "navy", linewidth = 0.9) +
    ggplot2::geom_point(
      data = snap_df,
      ggplot2::aes(old_lon, old_lat),
      colour = "red3", size = 3.5, shape = 4, stroke = 1.8) +
    ggplot2::geom_point(
      data = snap_df,
      ggplot2::aes(new_lon, new_lat),
      colour = "darkgreen", size = 3, shape = 16) +
    ggplot2::geom_text(
      data = snap_df,
      ggplot2::aes(old_lon, old_lat, label = label),
      size = 2.8, hjust = -0.08, vjust = 0.5, colour = "black") +
    ggplot2::coord_sf(xlim = xlim_s, ylim = ylim_s, expand = FALSE) +
    ggplot2::labs(
      title    = "GPS snap corrections: clusters outside Nigeria moved to nearest LGA centroid",
      subtitle = "Red X = original GPS (outside Nigeria); green dot = LGA centroid (snapped to); arrow = correction applied",
      x = NULL, y = NULL) +
    ggplot2::theme_minimal()
  
  grDevices::pdf(plot_file, width = 11, height = 8.5)
  print(p)
  grDevices::dev.off()
  message("  GPS: snap correction map saved to ", basename(plot_file))
}


# Produce a multi-page GPS diagnostic PDF (before/after maps, spread scatter,
# LGA/state inconsistency map, optional snap-correction page).
# Inputs: hh_means_all, clust_original, clust_after -- household/cluster data frames;
#   admin_shape; snap_log -- NULL or snap data frame from the snap step.
# Outputs: PDF saved to plot_file.
.plot_gps_diagnostics <- function(
    hh_means_all,        # HH data: post-snap/pre-null snapshot with all flag columns
    clust_original,      # cluster centroids computed from pre-null hh_means (post-snap)
    clust_after,         # final cluster centroids after GPS nulling + Step 2
    admin_shape,
    max_clust_spread_km,
    plot_file,
    snap_log = NULL      # data.frame from Step 1.7 snap corrections, or NULL
) {
  for (pkg in c("ggplot2", "sf")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning("Package '", pkg, "' needed for GPS diagnostic plots — skipping.")
      return(invisible(NULL))
    }
  }
  
  admin_sf <- sf::st_as_sf(admin_shape)
  if (is.na(sf::st_crs(admin_sf)) || !identical(sf::st_crs(admin_sf), sf::st_crs(4326)))
    admin_sf <- sf::st_transform(admin_sf, 4326)
  
  # Shared helper: LGA centroid labels (used on several pages)
  .lga_labels <- function() {
    tryCatch({
      ctrd   <- suppressWarnings(sf::st_centroid(admin_sf))
      coords <- as.data.frame(sf::st_coordinates(ctrd))
      data.frame(NOMDEP = admin_sf$NOMDEP, X = coords$X, Y = coords$Y,
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
  }
  
  # Restore original outside-Nigeria lat/lon for BEFORE map.
  # hh_means_all and clust_original were captured post-snap, so snapped clusters
  # show LGA centroid positions; replace these with old_lat/old_lon from snap_log.
  hh_before    <- hh_means_all
  clust_before <- clust_original
  if (!is.null(snap_log) && nrow(snap_log) > 0) {
    for (i in seq_len(nrow(snap_log))) {
      cl <- snap_log$cluster[i]
      hh_idx <- hh_before$cluster == cl
      hh_before$lat[hh_idx] <- snap_log$old_lat[i]
      hh_before$lon[hh_idx] <- snap_log$old_lon[i]
      cl_idx <- clust_before$cluster == cl
      clust_before$latitude[cl_idx]  <- snap_log$old_lat[i]
      clust_before$longitude[cl_idx] <- snap_log$old_lon[i]
    }
  }
  
  # Ensure needed flag columns exist (defaults to FALSE if absent)
  for (col in c("spread_flagged", "state_filter_removed", "outlier_removed",
                "cross_lga_in_cluster", "was_outside_nigeria")) {
    if (!col %in% colnames(hh_before))
      hh_before[[col]] <- FALSE
  }
  
  # ---- Cluster-level issue summary ----
  clust_issues <- hh_before %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      has_outside      = any(was_outside_nigeria %in% TRUE, na.rm = TRUE),
      has_state_filter = any(state_filter_removed %in% TRUE, na.rm = TRUE),
      has_cross_lga    = any(cross_lga_in_cluster %in% TRUE, na.rm = TRUE),
      has_spread       = any(outlier_removed %in% TRUE, na.rm = TRUE),
      .groups          = "drop"
    )
  clust_issues$issue_type <- with(clust_issues, ifelse(
    has_outside,      "Outside Nigeria",
    ifelse(has_state_filter, "Cross-state (GPS removed)",
           ifelse(has_cross_lga,    "Cross-LGA (flagged)",
                  ifelse(has_spread,       "Spread outlier (GPS removed)",
                         "No issues")))))
  clust_issues$issue_type <- factor(clust_issues$issue_type,
                                    levels = c("Outside Nigeria", "Cross-state (GPS removed)",
                                               "Cross-LGA (flagged)", "Spread outlier (GPS removed)", "No issues"))
  
  clust_before_issues <- merge(clust_before, clust_issues, by = "cluster", all.x = TRUE)
  clust_before_issues$issue_type[is.na(clust_before_issues$issue_type)] <- "No issues"
  
  issue_colors <- c(
    "Outside Nigeria"              = "red3",
    "Cross-state (GPS removed)"    = "darkorange",
    "Cross-LGA (flagged)"          = "goldenrod3",
    "Spread outlier (GPS removed)" = "steelblue3",
    "No issues"                    = "grey70"
  )
  problem_ids <- clust_issues$cluster[clust_issues$issue_type != "No issues"]
  hh_problems <- hh_before[hh_before$cluster %in% problem_ids, ]
  
  n_outside <- sum(clust_issues$has_outside,      na.rm = TRUE)
  n_state   <- sum(clust_issues$has_state_filter, na.rm = TRUE)
  n_lga     <- sum(clust_issues$has_cross_lga,    na.rm = TRUE)
  n_spread  <- sum(clust_issues$has_spread,        na.rm = TRUE)
  
  # ======================================================================== #
  # PAGE 1: BEFORE — all issues visible                                       #
  # ======================================================================== #
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = admin_sf, fill = "grey97", colour = "grey70", linewidth = 0.2) +
    ggplot2::geom_point(
      data = hh_problems,
      ggplot2::aes(lon, lat),
      colour = "grey55", size = 0.8, alpha = 0.45) +
    ggplot2::geom_point(
      data = clust_before_issues[clust_before_issues$issue_type == "No issues", ],
      ggplot2::aes(longitude, latitude),
      colour = "grey70", size = 0.9, alpha = 0.6) +
    ggplot2::geom_point(
      data = clust_before_issues[clust_before_issues$issue_type != "No issues", ],
      ggplot2::aes(longitude, latitude, colour = issue_type),
      size = 2.5, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = issue_colors, name = "Issue type", drop = FALSE) +
    ggplot2::labs(
      title    = "GPS - BEFORE cleaning: all identified issues",
      subtitle = paste0(
        "Outside Nigeria: ", n_outside, "  |  Cross-state: ", n_state,
        "  |  Cross-LGA: ", n_lga, "  |  Spread outliers: ", n_spread,
        "  |  Grey dots = individual HHs in problem clusters"),
      x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  
  # ======================================================================== #
  # PAGE 2: Spread outlier resolution                                         #
  # ======================================================================== #
  flagged_cl <- unique(hh_before$cluster[!is.na(hh_before$dist_from_median_km)])
  hh_flagged <- hh_before[hh_before$cluster %in% flagged_cl, ]
  p2 <- NULL
  if (nrow(hh_flagged) > 0) {
    clust_order <- hh_flagged %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(max_d = max(dist_from_median_km, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(max_d)) %>%
      dplyr::pull(cluster)
    hh_flagged$cluster_f <- factor(hh_flagged$cluster, levels = clust_order)
    
    p2 <- ggplot2::ggplot(hh_flagged,
                          ggplot2::aes(cluster_f, dist_from_median_km, colour = outlier_removed)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::geom_hline(yintercept = max_clust_spread_km,
                          linetype = "dashed", colour = "red3", linewidth = 0.8) +
      ggplot2::scale_colour_manual(
        values = c("FALSE" = "steelblue3", "TRUE" = "red3"),
        labels = c("Kept", "GPS nulled"), name = "Household") +
      ggplot2::labs(
        title    = "Spread outlier resolution: household distance from cluster median",
        subtitle = paste0(
          "Check: households with within-cluster spread > ", max_clust_spread_km, " km flagged.\n",
          "Resolution: households above dashed threshold have GPS coordinates set to NA",
          " (non-GPS data preserved)."),
        x = "Cluster (sorted by max HH distance)", y = "Distance from cluster median (km)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7))
  }
  
  # ======================================================================== #
  # PAGE 3: LGA / state inconsistency resolution                              #
  # ======================================================================== #
  state_cl   <- unique(hh_before$cluster[hh_before$state_filter_removed %in% TRUE])
  cross_lga_cl <- unique(hh_before$cluster[hh_before$cross_lga_in_cluster %in% TRUE])
  affected_cl  <- unique(c(state_cl, cross_lga_cl))
  p3 <- NULL
  if (length(affected_cl) > 0) {
    hh_lga <- hh_before[hh_before$cluster %in% affected_cl, ]
    hh_lga$status <- ifelse(
      hh_lga$state_filter_removed %in% TRUE, "Cross-state: GPS nulled",
      ifelse(hh_lga$outlier_removed %in% TRUE, "Spread outlier: GPS nulled", "Kept"))
    hh_lga$status <- factor(hh_lga$status,
                            levels = c("Kept", "Cross-state: GPS nulled", "Spread outlier: GPS nulled"))
    
    lga_lbl <- .lga_labels()
    
    p3 <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = admin_sf, fill = "grey97", colour = "grey70", linewidth = 0.2) +
      {if (!is.null(lga_lbl))
        ggplot2::geom_text(data = lga_lbl,
                           ggplot2::aes(X, Y, label = NOMDEP),
                           size = 1.9, colour = "grey50", check_overlap = TRUE)} +
      ggplot2::geom_point(data = hh_lga,
                          ggplot2::aes(lon, lat, colour = status, shape = status),
                          size = 2.5, alpha = 0.85) +
      ggplot2::scale_colour_manual(
        values = c("Kept" = "steelblue3",
                   "Cross-state: GPS nulled"    = "darkorange",
                   "Spread outlier: GPS nulled" = "red3"),
        name = "Household status") +
      ggplot2::scale_shape_manual(
        values = c("Kept" = 16L,
                   "Cross-state: GPS nulled"    = 4L,
                   "Spread outlier: GPS nulled" = 4L),
        name = "Household status") +
      ggplot2::labs(
        title    = "LGA / state inconsistency resolution",
        subtitle = paste0(
          "Check: households within the same cluster assigned to different states or LGAs.\n",
          "Resolution: where GPS-assigned state differs from MIS-recorded state, households",
          " not matching MIS state have GPS nulled (", length(state_cl), " cluster(s)).\n",
          "Cross-LGA within same state is flagged only (", length(cross_lga_cl),
          " cluster(s)); cluster median used for LGA assignment."),
        x = NULL, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "right")
  }
  
  # ======================================================================== #
  # PAGE 4: Outside-Nigeria snap resolution (omitted if no snaps)             #
  # ======================================================================== #
  p4 <- NULL
  if (!is.null(snap_log) && nrow(snap_log) > 0) {
    pad   <- 1.0
    xlim4 <- range(c(snap_log$old_lon, snap_log$new_lon), na.rm = TRUE) + c(-pad, pad)
    ylim4 <- range(c(snap_log$old_lat, snap_log$new_lat), na.rm = TRUE) + c(-pad, pad)
    lga_lbl4 <- .lga_labels()
    snap_log$label <- paste0("Cluster ", snap_log$cluster,
                             "\n-> ", snap_log$snap_lga, ", ", snap_log$snap_state,
                             " (", round(snap_log$border_dist_miles, 1), " mi)")
    
    p4 <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = admin_sf, fill = "grey95", colour = "grey65", linewidth = 0.3) +
      {if (!is.null(lga_lbl4))
        ggplot2::geom_text(data = lga_lbl4,
                           ggplot2::aes(X, Y, label = NOMDEP),
                           size = 2.3, colour = "grey45", check_overlap = TRUE)} +
      ggplot2::geom_segment(
        data = snap_log,
        ggplot2::aes(x = old_lon, y = old_lat, xend = new_lon, yend = new_lat),
        arrow  = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed"),
        colour = "navy", linewidth = 0.9) +
      ggplot2::geom_point(data = snap_log,
                          ggplot2::aes(old_lon, old_lat),
                          colour = "red3", size = 3.5, shape = 4, stroke = 1.8) +
      ggplot2::geom_point(data = snap_log,
                          ggplot2::aes(new_lon, new_lat),
                          colour = "darkgreen", size = 3, shape = 16) +
      ggplot2::geom_text(data = snap_log,
                         ggplot2::aes(old_lon, old_lat, label = label),
                         size = 2.8, hjust = -0.1, colour = "black") +
      ggplot2::coord_sf(xlim = xlim4, ylim = ylim4, expand = FALSE) +
      ggplot2::labs(
        title    = "Outside-Nigeria cluster resolution",
        subtitle = paste0(
          "Check: cluster(s) where ALL households fell outside Nigeria's shapefile boundaries.\n",
          "Resolution: nearest LGA by border distance identified; if border within ",
          "max_outside_snap_miles threshold, cluster GPS replaced with LGA centroid (green dot).\n",
          "Red X = original GPS; arrow = correction applied."),
        x = NULL, y = NULL) +
      ggplot2::theme_minimal()
  }
  
  # ======================================================================== #
  # PAGE 5: AFTER — final cluster positions                                   #
  # ======================================================================== #
  clust_after$flag_label <- ifelse(
    clust_after$clust_spread_flag %in% TRUE, "Remaining spread flag", "Clean")
  
  p5 <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = admin_sf, fill = "grey97", colour = "grey70", linewidth = 0.2) +
    ggplot2::geom_point(
      data = clust_after,
      ggplot2::aes(longitude, latitude, colour = flag_label),
      size = 1.3, alpha = 0.8) +
    ggplot2::scale_colour_manual(
      values = c("Clean" = "steelblue3", "Remaining spread flag" = "darkorange"),
      name = NULL) +
    ggplot2::labs(
      title    = "GPS - AFTER cleaning: final cluster positions",
      subtitle = paste0(
        nrow(clust_after), " clusters total; ",
        sum(clust_after$clust_spread_flag %in% TRUE, na.rm = TRUE),
        " still spread-flagged after cleaning (orange)."),
      x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  
  grDevices::pdf(plot_file, width = 11, height = 8.5)
  print(p1)
  if (!is.null(p2)) print(p2)
  if (!is.null(p3)) print(p3)
  if (!is.null(p4)) print(p4)
  print(p5)
  grDevices::dev.off()
  message("  GPS: diagnostic plots (", 2L + !is.null(p2) + !is.null(p3) + !is.null(p4),
          " pages) saved to ", basename(plot_file))
}

# Build a cluster-level GPS data frame from a CSPro household-level .dat GPS file.
# Steps: read layout, remove invalid coords, compute HH centroids, cluster centroids,
#   flag/remove outlier HHs, apply state-consistency filter, snap outside-boundary
#   clusters to nearest LGA centroid.
# Inputs: dat_path -- path to GPS .dat file; column name overrides; admin_shape.
# Returns: data frame with clusterid, latitude, longitude, spread flags.
extract_gps_from_nets <- function(
    dat_path,
    sps_path            = NULL,
    cluster_col         = "QHCLUST",
    hh_col              = "QHNUMBER",
    lat_col             = "GHLATITUDE",
    lon_col             = "GHLONGITUDE",
    max_hh_spread_km    = 1,
    max_clust_spread_km    = 5,
    max_outside_snap_miles = 20,  # snap entirely-outside-Nigeria clusters to nearest LGA centroid if border is within this many miles
    gps_diag_file          = NULL,  # if provided, write per-household detail for spread-flagged clusters
    mis_state_col       = NULL,  # column name in NETS.DAT for MIS-recorded state  (from recode_df state_code)
    mis_lga_col         = NULL,  # column name in NETS.DAT for MIS-recorded LGA    (from recode_df admin_code)
    admin_shape         = NULL,  # SpatialPolygonsDataFrame; used for LGA assignment in CSV + plots
    create_plots        = FALSE, # produce before/after maps, distance plot, and LGA consistency plot
    plot_file           = NULL   # path for PDF output
) {
  # --- find layout file (.sps preferred, .dcf as fallback) ---
  if (is.null(sps_path)) {
    sps_path <- sub("(?i)\\.dat$", ".sps", dat_path, perl = TRUE)
    if (!file.exists(sps_path)) sps_path <- NULL
  }
  
  if (!is.null(sps_path)) {
    layout       <- parse_sps_layout(sps_path)
    value_labels <- parse_sps_value_labels(sps_path)
  } else {
    dcf_path <- sub("(?i)\\.dat$", ".dcf", dat_path, perl = TRUE)
    if (!file.exists(dcf_path))
      stop("No .sps or .dcf layout file found for: ", dat_path)
    layout       <- parse_dcf_layout(dcf_path)
    layout$ImplicitDec <- 0L
    value_labels <- list()
  }
  
  # --- validate needed columns exist in layout, then read only those columns ---
  required_cols     <- c(cluster_col, hh_col, lat_col, lon_col)
  missing_in_layout <- required_cols[!required_cols %in% layout$Name]
  if (length(missing_in_layout) > 0)
    stop("Column(s) not found in GPS file layout: ",
         paste(missing_in_layout, collapse = ", "),
         "\n  Available: ", paste(layout$Name, collapse = ", "))
  
  # Optional MIS state/LGA columns — warn and skip if absent
  extra_cols    <- na.omit(unique(c(mis_state_col, mis_lga_col)))
  extra_present <- extra_cols[extra_cols %in% layout$Name]
  extra_missing <- setdiff(extra_cols, extra_present)
  if (length(extra_missing) > 0)
    warning("MIS state/LGA column(s) not found in NETS.DAT layout, skipping: ",
            paste(extra_missing, collapse = ", "))
  needed_cols <- c(required_cols, extra_present)
  
  layout_sub <- layout[layout$Name %in% needed_cols, , drop = FALSE]
  df <- readr::read_fwf(
    dat_path,
    col_positions  = readr::fwf_positions(
      start     = layout_sub$Start,
      end       = layout_sub$Start + layout_sub$Len - 1L,
      col_names = layout_sub$Name
    ),
    col_types      = readr::cols(.default = readr::col_character()),
    progress       = FALSE,
    show_col_types = FALSE
  )
  df <- mutate(df, across(where(is.character), trimws))
  
  # Apply implicit decimal scaling for non-GPS columns that need it.
  # GPS lat/lon columns are always stored as decimal degrees — the SPS (6) token
  # describes storage precision, not integer scaling, so we skip them here.
  gps_cols <- c(lat_col, lon_col)
  for (i in seq_len(nrow(layout_sub))) {
    dec <- layout_sub$ImplicitDec[i]
    if (!is.na(dec) && dec > 0 && !layout_sub$Name[i] %in% gps_cols)
      df[[layout_sub$Name[i]]] <- suppressWarnings(as.numeric(df[[layout_sub$Name[i]]])) / (10^dec)
  }
  
  # --- convert to appropriate types ---
  df[[lat_col]]     <- suppressWarnings(as.numeric(df[[lat_col]]))
  df[[lon_col]]     <- suppressWarnings(as.numeric(df[[lon_col]]))
  df[[cluster_col]] <- trimws(as.character(df[[cluster_col]]))
  df[[hh_col]]      <- trimws(as.character(df[[hh_col]]))
  
  # --- remove invalid GPS rows ---
  invalid <- is.na(df[[lat_col]]) | is.na(df[[lon_col]]) |
    abs(df[[lat_col]]) > 90 | abs(df[[lon_col]]) > 180 |
    (df[[lat_col]] == 0 & df[[lon_col]] == 0)
  n_invalid <- sum(invalid)
  if (n_invalid > 0)
    message("  GPS: removed ", n_invalid, " row(s) with invalid coordinates.")
  df <- df[!invalid, , drop = FALSE]
  
  if (nrow(df) == 0) stop("No valid GPS coordinates remain after filtering in: ", dat_path)
  
  # Rename to standard names before grouping
  gps <- data.frame(
    cluster = df[[cluster_col]],
    hh      = df[[hh_col]],
    lat     = df[[lat_col]],
    lon     = df[[lon_col]],
    stringsAsFactors = FALSE
  )
  .decode_col <- function(vals, col_name) {
    raw <- trimws(as.character(vals))
    if (col_name %in% names(value_labels)) {
      lbl <- value_labels[[col_name]]
      decoded <- lbl[raw]               # named-vector lookup by code
      ifelse(is.na(decoded), raw, decoded)   # fall back to raw code if no label found
    } else raw
  }
  if (!is.null(mis_state_col) && mis_state_col %in% colnames(df))
    gps$MIS_state <- .decode_col(df[[mis_state_col]], mis_state_col)
  if (!is.null(mis_lga_col)   && mis_lga_col   %in% colnames(df))
    gps$MIS_LGA   <- .decode_col(df[[mis_lga_col]],   mis_lga_col)
  
  # --- Step 1: household-level means ---
  hh_means <- gps %>%
    group_by(cluster, hh) %>%
    summarise(
      lat          = mean(lat, na.rm = TRUE),
      lon          = mean(lon, na.rm = TRUE),
      hh_spread_km = .max_spread_km(lat, lon),
      .groups      = "drop"
    )
  
  # Attach MIS state/LGA — constant within a household so take the first value
  for (mis_col in intersect(c("MIS_state", "MIS_LGA"), colnames(gps))) {
    first_vals <- gps %>%
      group_by(cluster, hh) %>%
      summarise(!!mis_col := first(.data[[mis_col]]), .groups = "drop")
    hh_means <- merge(hh_means, first_vals, by = c("cluster", "hh"), all.x = TRUE)
  }
  
  n_hh_flagged <- sum(hh_means$hh_spread_km > max_hh_spread_km, na.rm = TRUE)
  if (n_hh_flagged > 0)
    message("  GPS: ", n_hh_flagged, " household(s) with within-household spread > ",
            max_hh_spread_km, " km — averaged anyway.")
  
  hh_means$hh_spread_flag <- hh_means$hh_spread_km > max_hh_spread_km
  
  # --- Step 1.5: remove outlier households from spread-flagged clusters ---
  # For any cluster whose household spread would exceed max_clust_spread_km,
  # compute each household's distance from the cluster median centroid and
  # drop households farther than max_clust_spread_km from that median.
  # This targets individual GPS errors rather than discarding entire clusters.
  init_spread <- hh_means %>%
    group_by(cluster) %>%
    summarise(
      init_spread_km = .max_spread_km(lat, lon),
      med_lat        = median(lat, na.rm = TRUE),
      med_lon        = median(lon, na.rm = TRUE),
      .groups        = "drop"
    )
  
  flagged_clusters <- init_spread$cluster[init_spread$init_spread_km > max_clust_spread_km]
  
  hh_means <- merge(hh_means,
                    init_spread[, c("cluster", "med_lat", "med_lon", "init_spread_km")],
                    by = "cluster", all.x = TRUE)
  
  # Distance from cluster median — only computed for households in flagged clusters
  hh_means$dist_from_median_km <- ifelse(
    hh_means$cluster %in% flagged_clusters,
    .haversine_km(hh_means$lat, hh_means$lon, hh_means$med_lat, hh_means$med_lon),
    NA_real_
  )
  hh_means$outlier_removed <- !is.na(hh_means$dist_from_median_km) &
    hh_means$dist_from_median_km > max_clust_spread_km
  
  hh_means$state_filter_removed <- FALSE
  
  # --- Step 1.6: spatial join + state-consistency filter ---
  # Join each household centroid to its shapefile LGA/state (admin_name / region).
  # Then, for clusters where the GPS-assigned state varies across households, keep
  # only households whose region fuzzy-matches the MIS-recorded state (MIS_state).
  # Non-matching households have their GPS nulled (coordinates → NA) rather than
  # being dropped — their non-GPS data remains valid.
  if (!is.null(admin_shape) && requireNamespace("sf", quietly = TRUE)) {
    admin_sf_join <- sf::st_as_sf(admin_shape)
    if (is.na(sf::st_crs(admin_sf_join)) ||
        !identical(sf::st_crs(admin_sf_join), sf::st_crs(4326)))
      admin_sf_join <- sf::st_transform(admin_sf_join, 4326)
    valid_hh <- !is.na(hh_means$lat) & !is.na(hh_means$lon)
    hh_sf_join <- sf::st_as_sf(hh_means[valid_hh, ],
                               coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    joined_16   <- suppressMessages(
      sf::st_join(hh_sf_join, admin_sf_join[, c("NOMDEP", "NOMREGION")])
    )
    jdf_16 <- as.data.frame(joined_16)[, c("cluster", "hh", "NOMDEP", "NOMREGION")]
    names(jdf_16)[names(jdf_16) == "NOMDEP"]    <- "admin_name"
    names(jdf_16)[names(jdf_16) == "NOMREGION"] <- "region"
    hh_means <- merge(hh_means, jdf_16, by = c("cluster", "hh"), all.x = TRUE)
    
    if ("MIS_state" %in% colnames(hh_means)) {
      .norm_str <- function(x) tolower(trimws(gsub("[^a-zA-Z0-9 ]", "", as.character(x))))
      # clusters where GPS-assigned state is not uniform
      clusters_multi_state <- unique(hh_means$cluster[
        !is.na(hh_means$region) &
          as.integer(ave(hh_means$region, hh_means$cluster,
                         FUN = function(r) length(unique(na.omit(r))))) > 1L
      ])
      for (cl in clusters_multi_state) {
        cl_idx  <- which(hh_means$cluster == cl)
        mis_val <- na.omit(hh_means$MIS_state[cl_idx])[1]
        if (length(mis_val) == 0 || !nzchar(trimws(mis_val))) next
        mis_norm <- .norm_str(mis_val)
        matches  <- vapply(hh_means$region[cl_idx], function(r) {
          if (is.na(r)) return(FALSE)
          isTRUE(agrepl(mis_norm, .norm_str(r), max.distance = 0.2, ignore.case = TRUE))
        }, logical(1))
        # Only filter if at least one household matches; if none match keep all
        if (any(matches, na.rm = TRUE))
          hh_means$state_filter_removed[cl_idx[!matches]] <- TRUE
      }
      n_sf_rem <- sum(hh_means$state_filter_removed, na.rm = TRUE)
      if (n_sf_rem > 0)
        message("  GPS: ", n_sf_rem, " household(s) in ",
                length(clusters_multi_state), " cluster(s) removed by state-consistency filter.")
    }
    # Flag households originally outside Nigeria BEFORE snapping modifies admin_name
    hh_means$was_outside_nigeria <- is.na(hh_means$admin_name)
  }
  
  # --- Step 1.7: snap entirely-outside-Nigeria clusters to nearest LGA ---
  # For clusters where ALL households fall outside the shapefile, find the LGA whose
  # border is closest to the cluster centroid.  If that border is within
  # max_outside_snap_miles, replace every household's GPS with the centroid of that LGA
  # (so Step 2 will place the cluster there) and update admin_name / region accordingly.
  # Writes a snap log CSV and arrow-map PDF alongside gps_diag_file for documentation.
  snap_log <- NULL   # will hold data.frame of snapped clusters for plots/CSV
  if (!is.null(admin_shape) && requireNamespace("sf", quietly = TRUE) &&
      "admin_name" %in% colnames(hh_means)) {
    # Clusters where every household is outside Nigeria
    all_outside_clusters <- hh_means %>%
      group_by(cluster) %>%
      summarise(
        all_outside = all(is.na(admin_name)),
        mean_lat    = mean(lat, na.rm = TRUE),
        mean_lon    = mean(lon, na.rm = TRUE),
        .groups     = "drop"
      ) %>%
      dplyr::filter(all_outside, !is.na(mean_lat), !is.na(mean_lon))
    
    if (nrow(all_outside_clusters) > 0) {
      admin_sf_snap <- sf::st_as_sf(admin_shape)
      if (!identical(sf::st_crs(admin_sf_snap), sf::st_crs(4326)))
        admin_sf_snap <- sf::st_transform(admin_sf_snap, 4326)
      
      snap_rows <- list()
      for (i in seq_len(nrow(all_outside_clusters))) {
        cl      <- all_outside_clusters$cluster[i]
        pt      <- sf::st_sfc(sf::st_point(c(all_outside_clusters$mean_lon[i],
                                             all_outside_clusters$mean_lat[i])),
                              crs = 4326)
        # st_distance to polygon returns distance to nearest border when point is outside
        dists_m     <- as.numeric(sf::st_distance(pt, admin_sf_snap))
        nearest_idx <- which.min(dists_m)
        dist_km     <- dists_m[nearest_idx] / 1000
        dist_miles  <- dist_km / 1.60934
        nearest_lga   <- as.character(admin_sf_snap$NOMDEP[nearest_idx])
        nearest_state <- as.character(admin_sf_snap$NOMREGION[nearest_idx])
        
        if (dist_miles <= max_outside_snap_miles) {
          # Snap target: centroid of the nearest LGA polygon
          ctrd        <- suppressWarnings(sf::st_centroid(admin_sf_snap[nearest_idx, ]))
          ctrd_coords <- sf::st_coordinates(ctrd)
          new_lon     <- ctrd_coords[1, "X"]
          new_lat     <- ctrd_coords[1, "Y"]
          
          cl_idx <- which(hh_means$cluster == cl)
          snap_rows[[length(snap_rows) + 1]] <- data.frame(
            cluster    = cl,
            old_lat    = all_outside_clusters$mean_lat[i],
            old_lon    = all_outside_clusters$mean_lon[i],
            new_lat    = new_lat,
            new_lon    = new_lon,
            snap_lga   = nearest_lga,
            snap_state = nearest_state,
            border_dist_km    = round(dist_km,    2),
            border_dist_miles = round(dist_miles, 2),
            stringsAsFactors  = FALSE
          )
          hh_means$lat[cl_idx]        <- new_lat
          hh_means$lon[cl_idx]        <- new_lon
          hh_means$admin_name[cl_idx] <- nearest_lga
          hh_means$region[cl_idx]     <- nearest_state
          
          message(sprintf(
            "  GPS: cluster %s (all HHs outside Nigeria) snapped to centroid of %s, %s (nearest border: %.1f km / %.1f miles)",
            cl, nearest_lga, nearest_state, dist_km, dist_miles))
        } else {
          message(sprintf(
            "  GPS: cluster %s is %.1f miles from nearest LGA (%s) - exceeds max_outside_snap_miles=%.0f; not snapped.",
            cl, dist_miles, nearest_lga, max_outside_snap_miles))
        }
      }
      
      if (length(snap_rows) > 0) {
        snap_log <- do.call(rbind, snap_rows)
        # Write snap log CSV
        if (!is.null(gps_diag_file)) {
          snap_csv <- file.path(dirname(gps_diag_file), "gps_snap_log.csv")
          write.csv(snap_log, snap_csv, row.names = FALSE)
          message("  GPS: snap log written to ", basename(snap_csv))
        }
        # Write snap arrow-map PDF
        snap_pdf_path <- if (!is.null(gps_diag_file))
          file.path(dirname(gps_diag_file), "gps_snap_corrections.pdf")
        else if (!is.null(plot_file))
          file.path(dirname(plot_file), "gps_snap_corrections.pdf")
        else NULL
        if (!is.null(snap_pdf_path) && requireNamespace("ggplot2", quietly = TRUE)) {
          .plot_snap_corrections(snap_log, admin_sf_snap, snap_pdf_path)
        }
      }
    }
  }
  
  # Write per-household diagnostics before removing outliers.
  # Includes all households in any cluster that fails at least one GPS check:
  #   spread_flagged        — cluster household spread exceeds max_clust_spread_km
  #   cross_lga_in_cluster  — households in this cluster map to >1 LGA
  #   outside_nigeria       — this specific household falls outside Nigeria's boundaries
  if (!is.null(gps_diag_file)) {
    # admin_name / region were attached in Step 1.6 when admin_shape is available.
    # Here we add cluster-level consistency flags and print console diagnostics.
    if ("admin_name" %in% colnames(hh_means)) {
      # Cluster-level LGA consistency flags
      lga_consist <- hh_means %>%
        group_by(cluster) %>%
        summarise(
          cross_lga_in_cluster = n_distinct(admin_name, na.rm = TRUE) > 1,
          any_outside_nigeria  = any(is.na(admin_name)),
          .groups = "drop"
        )
      hh_means <- merge(hh_means, lga_consist, by = "cluster", all.x = TRUE)
      hh_means$outside_nigeria <- is.na(hh_means$admin_name)
      hh_means$spread_flagged  <- hh_means$cluster %in% flagged_clusters
      
      # Console report: GPS-removed households in a different LGA from the kept majority
      any_removed <- hh_means$outlier_removed | hh_means$state_filter_removed
      kept_lgas <- hh_means[!any_removed, c("cluster", "admin_name")]
      modal_lga <- do.call(rbind, lapply(split(kept_lgas, kept_lgas$cluster), function(x) {
        lgas <- na.omit(x$admin_name)
        data.frame(cluster   = x$cluster[1],
                   modal_lga = if (length(lgas)) names(sort(table(lgas), decreasing = TRUE))[1]
                   else NA_character_,
                   stringsAsFactors = FALSE)
      }))
      hh_means <- merge(hh_means, modal_lga, by = "cluster", all.x = TRUE)
      diff_lga <- hh_means[any_removed &
                             !is.na(hh_means$admin_name) & !is.na(hh_means$modal_lga) &
                             hh_means$admin_name != hh_means$modal_lga, ]
      if (nrow(diff_lga) > 0) {
        message("  GPS: ", nrow(diff_lga),
                " removed household(s) map to a different LGA than the rest of their cluster:")
        for (i in seq_len(nrow(diff_lga)))
          message(sprintf("    Cluster %s HH %s: removed LGA = '%s',  cluster LGA = '%s'",
                          diff_lga$cluster[i], diff_lga$hh[i],
                          diff_lga$admin_name[i], diff_lga$modal_lga[i]))
      }
      hh_means$modal_lga <- NULL
    } else {
      hh_means$admin_name           <- NA_character_
      hh_means$region               <- NA_character_
      hh_means$cross_lga_in_cluster <- NA
      hh_means$any_outside_nigeria  <- NA
      hh_means$outside_nigeria      <- NA
    }
    # MIS_LGA consistency: TRUE = all households in cluster share the same code
    if ("MIS_LGA" %in% colnames(hh_means)) {
      lga_code_consist <- hh_means %>%
        group_by(cluster) %>%
        summarise(
          mis_lga_consistent = n_distinct(MIS_LGA[!is.na(MIS_LGA)], na.rm = FALSE) <= 1,
          .groups = "drop"
        )
      hh_means <- merge(hh_means, lga_code_consist, by = "cluster", all.x = TRUE)
    } else {
      hh_means$mis_lga_consistent <- NA
    }
    
    # Collect all clusters failing at least one check
    problem_clusters <- unique(hh_means$cluster[
      hh_means$spread_flagged |
        hh_means$state_filter_removed |
        (is.logical(hh_means$cross_lga_in_cluster) & !is.na(hh_means$cross_lga_in_cluster) & hh_means$cross_lga_in_cluster) |
        (is.logical(hh_means$any_outside_nigeria)  & !is.na(hh_means$any_outside_nigeria)  & hh_means$any_outside_nigeria) |
        (is.logical(hh_means$mis_lga_consistent)   & !is.na(hh_means$mis_lga_consistent)   & !hh_means$mis_lga_consistent)
    ])
    
    if (length(problem_clusters) > 0) {
      diag_cols <- c("cluster", "hh", "lat", "lon",
                     "init_spread_km", "dist_from_median_km",
                     "outlier_removed", "state_filter_removed",
                     "spread_flagged", "cross_lga_in_cluster", "outside_nigeria",
                     "mis_lga_consistent",
                     intersect(c("MIS_state", "MIS_LGA", "admin_name", "region"),
                               colnames(hh_means)))
      diag_hh <- hh_means[hh_means$cluster %in% problem_clusters,
                          diag_cols[diag_cols %in% colnames(hh_means)], drop = FALSE]
      write.csv(diag_hh, gps_diag_file, row.names = FALSE)
      n_cross_lga <- if (exists("lga_consist"))
        sum(lga_consist$cross_lga_in_cluster | lga_consist$any_outside_nigeria, na.rm = TRUE)
      else NA_integer_
      message("  GPS: household-level diagnostics for ", length(problem_clusters),
              " cluster(s) written to ", basename(gps_diag_file),
              " (spread: ", sum(unique(hh_means$cluster[hh_means$spread_flagged]) %in% problem_clusters),
              ", state-filter: ", sum(hh_means$state_filter_removed, na.rm = TRUE),
              ", cross-LGA: ", n_cross_lga, ")")
    } else {
      message("  GPS: no clusters failed any GPS checks — diagnostic CSV not written.")
    }
  }
  
  n_outlier_hh      <- sum(hh_means$outlier_removed, na.rm = TRUE)
  n_outlier_clusters <- length(unique(hh_means$cluster[hh_means$outlier_removed]))
  if (n_outlier_hh > 0)
    message("  GPS: nulling GPS for ", n_outlier_hh, " spread-outlier household(s) across ",
            n_outlier_clusters, " cluster(s).")
  
  # Track removal counts per cluster
  removal_counts <- hh_means %>%
    group_by(cluster) %>%
    summarise(
      n_outlier_hh_removed      = sum(outlier_removed,      na.rm = TRUE),
      n_state_filter_hh_removed = sum(state_filter_removed, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Snapshot for plots: taken BEFORE Step 1.7 snapping so original outside-Nigeria
  # coordinates are preserved.  hh_means_for_plots has all flag columns from Steps
  # 1.5, 1.6 and the diagnostic block, but lat/lon still reflect raw survey GPS.
  if (create_plots) {
    hh_means_for_plots <- hh_means
    clust_original <- hh_means %>%
      group_by(cluster) %>%
      summarise(
        latitude        = mean(lat, na.rm = TRUE),
        longitude       = mean(lon, na.rm = TRUE),
        clust_spread_km = .max_spread_km(lat, lon),
        .groups         = "drop"
      ) %>%
      as.data.frame()
  }
  
  # Null out GPS coordinates for removed households rather than dropping rows.
  # The household still exists (non-GPS data is valid); only its GPS is untrusted.
  gps_removed <- hh_means$outlier_removed | hh_means$state_filter_removed
  hh_means$lat[gps_removed] <- NA_real_
  hh_means$lon[gps_removed] <- NA_real_
  
  # --- Step 2: cluster-level means (using cleaned household set) ---
  clust_means <- hh_means %>%
    group_by(cluster) %>%
    summarise(
      latitude         = mean(lat, na.rm = TRUE),
      longitude        = mean(lon, na.rm = TRUE),
      clust_spread_km  = .max_spread_km(lat, lon),
      n_households     = n(),
      n_hh_spread_flag = sum(hh_spread_flag, na.rm = TRUE),
      .groups          = "drop"
    )
  
  clust_means <- merge(clust_means, removal_counts, by = "cluster", all.x = TRUE)
  clust_means$n_outlier_hh_removed[is.na(clust_means$n_outlier_hh_removed)]           <- 0L
  clust_means$n_state_filter_hh_removed[is.na(clust_means$n_state_filter_hh_removed)] <- 0L
  
  n_clust_flagged <- sum(clust_means$clust_spread_km > max_clust_spread_km, na.rm = TRUE)
  if (n_clust_flagged > 0)
    message("  GPS: ", n_clust_flagged, " cluster(s) still have spread > ",
            max_clust_spread_km, " km after outlier removal.")
  
  clust_means$clust_spread_flag <- clust_means$clust_spread_km > max_clust_spread_km
  clust_means_df <- as.data.frame(rename(clust_means, clusterid = cluster))
  
  if (create_plots) {
    if (is.null(admin_shape))
      warning("create_plots = TRUE but admin_shape not provided — skipping GPS plots.")
    else {
      pf <- if (!is.null(plot_file)) plot_file else
        sub("\\.csv$", "_plots.pdf",
            if (!is.null(gps_diag_file)) gps_diag_file
            else file.path(dirname(dat_path), "DHS_cluster_gps_diagnostics.pdf"))
      .plot_gps_diagnostics(hh_means_for_plots, clust_original, clust_means_df,
                            admin_shape, max_clust_spread_km, pf,
                            snap_log = snap_log)
    }
  }
  
  clust_means_df
}


# Initialise cluster_df from survey data when no GPS shapefile is available.
# Reads admin_code (LGA) and optionally region_code (state) from the recode CSV.
# Returns: data frame with clusterid, NOMDEP, and optionally NOMREGION.
init_cluster_df_from_lga <- function(recode_df, dta_dir) {
  if (!"admin_code" %in% colnames(recode_df))
    stop("No GPS 'locations' row in recode CSV and no 'admin_code' column found. ",
         "Add admin_code (and optionally region_code) to the recode CSV for this year.")
  
  use_idx <- which(!is.na(recode_df$admin_code) &
                     nchar(trimws(as.character(recode_df$admin_code))) > 0)
  if (length(use_idx) == 0)
    stop("No GPS 'locations' row and no non-empty admin_code values in recode CSV. ",
         "Cannot initialise cluster_df.")
  
  idx        <- use_idx[1]
  path       <- build_file_path(dta_dir, recode_df$folder_dir[idx], recode_df$filename[idx])
  df         <- load_survey_file(path)
  cid_col    <- recode_df$cluster_id_code[idx]
  admin_col  <- as.character(recode_df$admin_code[idx])
  
  if (!cid_col   %in% colnames(df)) stop("cluster_id_code '", cid_col,   "' not in ", basename(path))
  if (!admin_col %in% colnames(df)) stop("admin_code '",      admin_col, "' not in ", basename(path))
  
  # One row per cluster — take the first occurrence (LGA should be constant within a cluster)
  cluster_df <- df[!duplicated(df[[cid_col]]), c(cid_col, admin_col), drop = FALSE]
  colnames(cluster_df)[colnames(cluster_df) == cid_col]   <- "clusterid"
  colnames(cluster_df)[colnames(cluster_df) == admin_col] <- "NOMDEP"
  
  # Optionally attach region/state
  has_region <- "region_code" %in% colnames(recode_df) &&
    !is.na(recode_df$region_code[idx]) &&
    nchar(trimws(as.character(recode_df$region_code[idx]))) > 0
  if (has_region) {
    region_col <- as.character(recode_df$region_code[idx])
    if (!region_col %in% colnames(df)) {
      warning("region_code '", region_col, "' not found in ", basename(path),
              "; NOMREGION will be NA.")
    } else {
      region_lookup <- df[!duplicated(df[[cid_col]]), c(cid_col, region_col), drop = FALSE]
      colnames(region_lookup) <- c("clusterid", "NOMREGION")
      cluster_df <- merge(cluster_df, region_lookup, by = "clusterid", all.x = TRUE)
    }
  }
  
  rownames(cluster_df) <- NULL
  cluster_df
}
# =============================================================================
# SECTION 8: ADMIN/REGION/ARCHETYPE AGGREGATION
# Sum cluster-level counts to admin, state, zone, and national levels.
# apply_min_sample_fallback cascades through fallback levels for small samples.
# =============================================================================

# Aggregate cluster-level _num_true/_num_total counts to a spatial grouping level.
# Inputs: cluster_df; group_cols -- character vector (empty = national); variables.
# Returns: data frame with summed counts and recomputed rates.
compute_level_sums <- function(cluster_df, group_cols, variables) {
  num_cols <- c(grep("num_total$", names(cluster_df), value = TRUE),
                grep("num_true$",  names(cluster_df), value = TRUE))
  cols     <- intersect(c(group_cols, num_cols), colnames(cluster_df))
  
  if (length(group_cols) > 0) {
    sums <- cluster_df[, cols] %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  } else {
    sums <- cluster_df[, intersect(num_cols, colnames(cluster_df))] %>%
      summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
  }
  for (v in variables) {
    tc <- paste0(v, "_num_true"); oc <- paste0(v, "_num_total"); rc <- paste0(v, "_rate")
    if (tc %in% colnames(sums) && oc %in% colnames(sums))
      sums[[rc]] <- sums[[tc]] / sums[[oc]]
  }
  sums
}

# Apply hierarchical minimum-sample-size fallback to admin_sums.
# For each variable, admins below min_num_total cascade through fallback levels
#   (state -> zone/archetype -> national) until sufficient data is found.
# Inputs: admin_sums; fallback_list -- ordered list of {join_col, data} entries.
# Returns: admin_sums with small-sample cells replaced by fallback estimates.
apply_min_sample_fallback <- function(admin_sums, fallback_list, variables,
                                      min_num_total, include_date = TRUE) {
  for (v in variables) {
    oc <- paste0(v, "_num_total"); tc <- paste0(v, "_num_true"); rc <- paste0(v, "_rate")
    if (!oc %in% colnames(admin_sums)) next
    
    # Track which rows still need filling; process fallback levels in priority order.
    needs_fill <- is.na(admin_sums[[oc]]) | admin_sums[[oc]] < min_num_total
    
    for (fb in fallback_list) {
      if (!any(needs_fill)) break
      fb_df <- fb$data
      jc    <- fb$join_col   # NULL = national (single-row) fallback
      if (!oc %in% colnames(fb_df)) next
      
      if (is.null(jc)) {
        # National fallback applies to all remaining unfilled rows in one pass
        admin_sums[[oc]][needs_fill] <- fb_df[[oc]][1]
        admin_sums[[tc]][needs_fill] <- fb_df[[tc]][1]
        admin_sums[[rc]][needs_fill] <- fb_df[[rc]][1]
        if (include_date && grepl("itn", v) &&
            "mean_date" %in% colnames(admin_sums) && "mean_date" %in% colnames(fb_df))
          admin_sums$mean_date[needs_fill] <- fb_df$mean_date[1]
        needs_fill[] <- FALSE
      } else {
        if (!jc %in% colnames(admin_sums) || !jc %in% colnames(fb_df)) next
        # Align fallback rows to admin_sums rows via match, then fill in batch
        fb_aligned <- fb_df[match(admin_sums[[jc]], fb_df[[jc]]), , drop = FALSE]
        fill_rows  <- needs_fill & !is.na(fb_aligned[[oc]])
        if (!any(fill_rows)) next
        admin_sums[[oc]][fill_rows] <- fb_aligned[[oc]][fill_rows]
        admin_sums[[tc]][fill_rows] <- fb_aligned[[tc]][fill_rows]
        admin_sums[[rc]][fill_rows] <- fb_aligned[[rc]][fill_rows]
        if (include_date && grepl("itn", v) &&
            "mean_date" %in% colnames(admin_sums) && "mean_date" %in% colnames(fb_aligned))
          admin_sums$mean_date[fill_rows] <- fb_aligned$mean_date[fill_rows]
        # Only mark as done when the filled value itself meets the threshold.
        # If the fallback level also has < min_num_total, leave needs_fill TRUE
        # so the next (higher) level is tried — matching the old code behaviour of
        # cascading state → zone → national until sufficient data is found.
        newly_sufficient <- fill_rows & !is.na(admin_sums[[oc]]) &
          admin_sums[[oc]] >= min_num_total
        needs_fill[newly_sufficient] <- FALSE
      }
    }
  }
  admin_sums
}

# Merge archetype/pop file into admin_sums, expand to all admins, zero-fill missing.
# Emits warnings for LGA or state names not matching the archetype file.
# Returns: expanded data frame aligned to archetype_info.
expand_and_merge_archetype <- function(admin_sums, archetype_info,
                                       include_zone = TRUE) {
  archetype_info$name_match <- vapply(archetype_info$NOMDEP, normalize_admin_name, character(1))
  admin_sums$name_match     <- vapply(admin_sums$NOMDEP,     normalize_admin_name, character(1))
  
  unmatched_adm <- admin_sums$name_match[!admin_sums$name_match %in% archetype_info$name_match]
  unmatched_reg <- admin_sums$NOMREGION[ !admin_sums$NOMREGION  %in% archetype_info$NOMREGION]
  if (length(unmatched_adm) > 0) warning("Unmatched LGA names (", length(unmatched_adm), "): ",
                                         paste(unique(unmatched_adm), collapse=", "))
  if (length(unmatched_reg) > 0) warning("Unmatched state names (", length(unmatched_reg), "): ",
                                         paste(unique(unmatched_reg), collapse=", "))
  
  colnames(admin_sums)[colnames(admin_sums) == "NOMDEP"] <- "NOMDEP_dhs_orig"
  
  arch_cols <- c("name_match","NOMDEP","NOMREGION","Archetype")
  if (include_zone && "ZONE" %in% colnames(archetype_info)) arch_cols <- c(arch_cols, "ZONE")
  
  expanded <- merge(admin_sums, archetype_info[, arch_cols], all = TRUE)
  expanded  <- expanded[, !colnames(expanded) %in% c("NOMDEP_dhs_orig","name_match")]
  
  # Zero-fill total columns for admins with no survey presence
  oc_cols <- grep("num_total", colnames(expanded))
  expanded[, oc_cols][is.na(expanded[, oc_cols])] <- 0
  expanded
}


# Add act_middle_rate_* columns as a weighted blend of low and high ACT estimates.
# Inputs: df; w -- weight assigned to the high estimate (0-1).
# Returns: df with act_middle_rate_{public,private_formal,private_informal} added.
.add_act_middle_rates <- function(df, w) {
  mutate(df,
         act_middle_rate_public           = act_low_public_rate           * (1 - w) + act_high_public_rate           * w,
         act_middle_rate_private_formal   = act_low_private_formal_rate   * (1 - w) + act_high_private_formal_rate   * w,
         act_middle_rate_private_informal = act_low_private_informal_rate * (1 - w) + act_high_private_informal_rate * w
  )
}

# Rename art_(high|low)* columns to act_(high|low)* for consistent output naming.
# Returns: df with column names updated in-place.
.rename_art_to_act <- function(df) {
  names(df) <- sub("^art_(high|low)", "act_\\1", names(df))
  df
}


# Aggregate cluster-level ACT rates to a given spatial grouping (or national).
# Uses survey weights when present; falls back to raw count ratios.
# Inputs: cluster_df; group_col -- NULL for national; national_fallback -- optional
#   prior national result used to replace zone estimates with small sample sizes.
aggregate_act_rates <- function(cluster_df, group_col = NULL,
                                min_num_total = 30, national_fallback = NULL) {
  rate_cols  <- grep("^art_(high|low).*_rate$",      names(cluster_df), value = TRUE)
  count_cols <- grep("^art_(high|low).*_num_total$",  names(cluster_df), value = TRUE)
  true_cols  <- grep("^art_(high|low).*_num_true$",   names(cluster_df), value = TRUE)
  if (length(rate_cols) == 0) return(NULL)
  
  use_wt <- "itn_weights" %in% colnames(cluster_df)
  
  summarise_act <- function(df) {
    if (use_wt) {
      df %>% summarise(
        across(all_of(rate_cols),  ~ stats::weighted.mean(.x, itn_weights, na.rm=TRUE)),
        across(all_of(count_cols), ~ sum(.x, na.rm=TRUE)),
        across(all_of(true_cols),  ~ sum(.x, na.rm=TRUE)),
        .groups = "drop")
    } else {
      df %>% summarise(
        across(all_of(c(count_cols, true_cols)), ~ sum(.x, na.rm=TRUE)),
        .groups = "drop") %>%
        mutate(across(all_of(rate_cols), ~ {
          tc <- sub("_rate$", "_num_true",  cur_column())
          oc <- sub("_rate$", "_num_total", cur_column())
          if (tc %in% colnames(.) && oc %in% colnames(.)) .[[tc]] / .[[oc]] else NA_real_
        }))
    }
  }
  
  if (is.null(group_col)) {
    result <- summarise_act(cluster_df)
  } else {
    result <- cluster_df %>%
      group_by(across(all_of(group_col))) %>%
      group_modify(~ summarise_act(.x)) %>%
      ungroup()
    
    # Replace zone estimates with national where sample too small
    if (!is.null(national_fallback)) {
      for (oc in count_cols) {
        rc <- sub("_num_total$", "_rate", oc)
        if (!rc %in% colnames(result)) next
        small <- !is.na(result[[oc]]) & result[[oc]] < min_num_total
        if (any(small) && rc %in% colnames(national_fallback))
          result[[rc]][small] <- national_fallback[[rc]][1]
      }
    }
  }
  result
}


# =============================================================================
# SECTION 9: CORE EXTRACTION -- extract_DHS_data
# Main entry point: loops over years, builds cluster-level outputs, applies spatial
# join, aggregates to admin/archetype/national level, writes standardised CSVs.
# =============================================================================

extract_DHS_data <- function(
    hbhi_dir, dta_dir, years, ds_pop_df_filename, admin_shape = NULL,
    min_num_total          = 30,
    variables              = c("mic","rdt","itn_all","itn_u5","itn_5_10","itn_10_15",
                               "itn_15_20","itn_o20","iptp","cm","blood_test",
                               "public_treatment","private_formal_treatment",
                               "private_informal_treatment","art_high","art_high_public",
                               "art_high_private_formal","art_high_private_informal",
                               "art_low","art_low_public","art_low_private_formal",
                               "art_low_private_informal"),
    cm_high_estimate_weight = 2/3,
    overwrite_files_flag   = FALSE,
    # If TRUE, ACT rates (by sector) are estimated at the national level and applied
    # uniformly to all LGAs. If FALSE, zone/archetype-level rates are used instead
    # (noisier but more spatially granular). National is the default based on
    # consistency testing showing too much noise at zone level.
    use_national_act = TRUE,
    # ACT codes may vary by year. Supply a named list keyed by year (as character)
    # with a "default" fallback, e.g.:
    #   art_codes = list("2024" = c("h37e"), default = c("ml13e"))
    art_codes     = c("ml13e"),
    non_art_codes = c("ml13a","ml13b","ml13c","ml13d","ml13aa","ml13da","ml13h"),
    create_gps_plots = FALSE  # produce GPS diagnostic PDF (requires ggplot2 + sf)
) {
  out_dir <- file.path(hbhi_dir, "estimates_from_DHS")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Helper: load a variable from its recode CSV row, code positivity, aggregate
  # to clusters, and merge into cluster_df. Defined once here — not per year.
  add_var <- function(cluster_df, var_name, recode_df,
                      alt_pos = character(0), age_filter = NULL,
                      include_weight = FALSE, verbose = TRUE) {
    idx <- which(recode_df$variable == var_name)
    if (length(idx) == 0) return(cluster_df)
    fn <- recode_df$filename[idx]
    if (is.na(fn) || !nzchar(trimws(fn))) return(cluster_df)
    path <- build_file_path(dta_dir, recode_df$folder_dir[idx], recode_df$filename[idx])
    df   <- load_survey_file(path)
    if (!is.null(age_filter)) {
      age_col  <- recode_df$age_code[idx]
      age_vals <- parse_dhs_age(df[[age_col]], age_col)
      df <- df[!is.na(age_vals) & age_vals > age_filter[1] & age_vals <= age_filter[2], ]
    }
    
    code_col   <- recode_df$code[idx]
    pos_pat    <- recode_df$pos_pattern[idx]
    neg_pat    <- recode_df$neg_pattern[idx]
    cid_col    <- recode_df$cluster_id_code[idx]
    
    if (is.na(cid_col) || !nzchar(trimws(cid_col))) {
      warning("add_var: cluster_id_code is missing for variable '", var_name,
              "' — please add the cluster ID column name to the recode CSV. Skipping.")
      return(cluster_df)
    }
    if (is.na(code_col) || !nzchar(trimws(code_col))) {
      message("  add_var: code column is empty for '", var_name,
              "' in the recode CSV — skipping (cannot compute positivity without a column name).")
      return(cluster_df)
    }
    
    if (verbose) {
      col_present <- code_col %in% colnames(df)
      col_class   <- if (col_present) class(df[[code_col]]) else "ABSENT"
      col_vals    <- if (col_present) paste(head(unique(df[[code_col]]), 8), collapse = ", ") else "N/A"
      message(sprintf("  [diag] %s: code='%s' class=%s | pos='%s' neg='%s' | sample vals: %s",
                      var_name, code_col, col_class, pos_pat, neg_pat, col_vals))
    }
    
    df  <- code_positivity(df, code_col, pos_pat, neg_pat, alt_pos)
    
    if (verbose) {
      message(sprintf("  [diag] %s: pos=1:%d  pos=0:%d  pos=NA:%d",
                      var_name,
                      sum(df$pos == 1, na.rm = TRUE),
                      sum(df$pos == 0, na.rm = TRUE),
                      sum(is.na(df$pos))))
    }
    
    res <- aggregate_to_clusters(df, cid_col, include_weight = include_weight)
    
    if (verbose) {
      message(sprintf("  [diag] %s: aggregated to %d clusters | cluster_df has %d clusters",
                      var_name, nrow(res), nrow(cluster_df)))
      if (nrow(res) > 0 && nrow(cluster_df) > 0) {
        n_match <- sum(res[[cid_col]] %in% cluster_df$clusterid)
        message(sprintf("  [diag] %s: %d / %d cluster IDs match cluster_df (type survey=%s, gps=%s)",
                        var_name, n_match, nrow(res),
                        class(res[[cid_col]]), class(cluster_df$clusterid)))
      }
    }
    
    cluster_df <- merge_cluster_result(cluster_df, res, cid_col, var_name)
    
    if (verbose) {
      rate_col <- paste0(var_name, "_rate")
      n_na <- if (rate_col %in% colnames(cluster_df)) sum(is.na(cluster_df[[rate_col]])) else NA
      message(sprintf("  [diag] %s: after merge — %s NA in %d rows",
                      var_name, n_na, nrow(cluster_df)))
    }
    
    if (!"mean_date" %in% colnames(cluster_df)) {
      dr <- extract_cluster_dates(df, cid_col,
                                  recode_df$month_code[idx], recode_df$year_code[idx])
      if (!is.null(dr))
        cluster_df <- merge(cluster_df, dr,
                            by.x = "clusterid", by.y = cid_col,
                            all = TRUE)
    }
    cluster_df
  }
  
  # Resolve year-specific ACT code lists. Supply a named list keyed by year
  # (as character) with a "default" entry, or a plain vector for all years.
  resolve_codes <- function(codes_arg, year) {
    if (is.list(codes_arg))
      codes_arg[[as.character(year)]] %||% codes_arg[["default"]] %||% codes_arg[[1]]
    else codes_arg
  }
  
  for (year in years) {
    message("\n=== Processing year: ", year, " ===")
    cluster_file <- file.path(out_dir, paste0("DHS_cluster_outputs_", year, ".csv"))
    
    # ------------------------------------------------------------------ #
    # STEP 1 — Build cluster-level output data frame                      #
    # ------------------------------------------------------------------ #
    if (file.exists(cluster_file) && !overwrite_files_flag) {
      message("  Loading existing cluster file.")
      cluster_df <- read.csv(cluster_file)[, -1]
    } else {
      recode_df <- read.csv(file.path(out_dir, paste0("DHS_", year, "_files_recodes_for_sims.csv")))
      plots_dir <- file.path(out_dir, "plots")
      if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

      # Cluster initialization:
      #   If recode CSV has a "locations" row pointing to a GPS shapefile, use that
      #   (classic DHS/MIS approach). Otherwise extract cluster IDs and LGA/state
      #   names directly from the data (new CSPro MIS format with admin_code column).
      loc_idx <- which(recode_df$variable == "locations")
      has_gps <- length(loc_idx) > 0 && !is.na(recode_df$filename[loc_idx])
      if (has_gps) {
        loc_path <- build_file_path(dta_dir, recode_df$folder_dir[loc_idx],
                                    recode_df$filename[loc_idx])
        loc_ext  <- tolower(tools::file_ext(loc_path))
        if (loc_ext == "dat") {
          # New format: household GPS coordinates in a .dat file (e.g. NETS.DAT + NETS.SPS).
          # Averages within household then within cluster; flags widely-spread clusters.
          message("  Building cluster GPS from household coordinates in ", basename(loc_path))
          # Column names in the GPS .dat file — read from recode CSV if provided,
          # otherwise fall back to the standard MIS CSPro variable names.
          # Recode CSV columns: cluster_id_code (existing), household_id_code,
          #                     latitude_code, longitude_code
          cid_col <- recode_df$cluster_id_code[loc_idx]
          hh_col  <- if ("household_id_code"   %in% colnames(recode_df) &&
                         !is.na(recode_df$household_id_code[loc_idx]))
            recode_df$household_id_code[loc_idx]   else "QHNUMBER"
          lat_col <- if ("latitude_code"  %in% colnames(recode_df) &&
                         !is.na(recode_df$latitude_code[loc_idx]))
            recode_df$latitude_code[loc_idx]  else "GHLATITUDE"
          lon_col <- if ("longitude_code" %in% colnames(recode_df) &&
                         !is.na(recode_df$longitude_code[loc_idx]))
            recode_df$longitude_code[loc_idx] else "GHLONGITUDE"
          # Read optional MIS state/LGA column names from the locations row.
          # Handles duplicate column names: R renames repeated headers to col.1, col.2, …
          # so we scan all matching columns and return the first non-NA value.
          .get_recode_col <- function(col) {
            pattern  <- paste0("^", col, "(\\.\\d+)?$")
            cols_hit <- grep(pattern, colnames(recode_df), value = TRUE, perl = TRUE)
            for (cn in cols_hit) {
              v <- recode_df[[cn]][loc_idx]
              if (!is.na(v) && nzchar(trimws(as.character(v)))) return(trimws(as.character(v)))
            }
            NULL
          }
          cluster_df <- extract_gps_from_nets(loc_path, cluster_col = cid_col,
                                              hh_col = hh_col,
                                              lat_col = lat_col, lon_col = lon_col,
                                              mis_state_col = .get_recode_col("state_code"),
                                              mis_lga_col   = .get_recode_col("admin_code"),
                                              gps_diag_file = file.path(plots_dir,
                                                                        paste0("DHS_cluster_gps_spread_", year, ".csv")),
                                              admin_shape   = admin_shape,
                                              create_plots  = create_gps_plots,
                                              plot_file     = file.path(plots_dir,
                                                                        paste0("DHS_cluster_gps_diagnostics_", year, ".pdf")))
        } else {
          # Classic format: GPS shapefile (DHSCLUST / LATNUM / LONGNUM).
          loc_shp    <- raster::shapefile(loc_path)
          cluster_df <- data.frame(clusterid = loc_shp$DHSCLUST,
                                   latitude  = loc_shp$LATNUM,
                                   longitude = loc_shp$LONGNUM)
        }
      } else {
        message("  No GPS locations file found; extracting LGA from survey data.")
        cluster_df <- init_cluster_df_from_lga(recode_df, dta_dir)
      }
      
      # PfPR
      cluster_df <- add_var(cluster_df, "mic",  recode_df)
      cluster_df <- add_var(cluster_df, "rdt",  recode_df)
      
      # ITNs — load once, filter by age group
      itn_idx <- which(recode_df$variable == "itn")
      if (length(itn_idx) > 0 && !is.na(recode_df$filename[itn_idx])) {
        itn_path <- build_file_path(dta_dir, recode_df$folder_dir[itn_idx], recode_df$filename[itn_idx])
        itn_df   <- load_survey_file(itn_path)
        
        itn_line_col <- if ("household_line_number_code" %in% colnames(recode_df) &&
                            !is.na(recode_df$household_line_number_code[itn_idx]))
          recode_df$household_line_number_code[itn_idx] else NULL
        cid      <- recode_df$cluster_id_code[itn_idx]
        age_col  <- recode_df$age_code[itn_idx]
        
        age_groups <- list(itn_all=c(-Inf,Inf), itn_u5=c(-Inf,5), itn_5_10=c(5,10),
                           itn_10_15=c(10,15), itn_15_20=c(15,20), itn_o20=c(20,Inf),
                           itn_o5=c(5,Inf))
        
        # --- CSPro multi-file path ---
        # When itn_line_col is set AND the recode CSV has a "locations" row pointing
        # to a .dat NETS file, use compute_itn_from_persons_and_nets():
        #   - PERSONS.DAT provides the unique person roster and ages
        #   - NETS.DAT provides per-net occupant line numbers (QH129A/B/C/D)
        # This is the authoritative approach for CSPro MIS surveys where NETS is
        # a sibling repeating group of PERSONS under the household.
        loc_idx <- which(recode_df$variable == "locations")
        nets_dat_path <- if (length(loc_idx) > 0 &&
                             !is.na(recode_df$filename[loc_idx]) &&
                             nzchar(trimws(recode_df$filename[loc_idx])) &&
                             tolower(tools::file_ext(recode_df$filename[loc_idx])) == "dat") {
          build_file_path(dta_dir, recode_df$folder_dir[loc_idx], recode_df$filename[loc_idx])
        } else NULL
        
        use_multifile_itn <- !is.null(itn_line_col) && !is.null(nets_dat_path)
        
        if (use_multifile_itn) {
          message("  ITN: using multi-file approach (PERSONS.DAT + NETS.DAT)")
          nets_df   <- load_survey_file(nets_dat_path)
          net_cid   <- recode_df$cluster_id_code[loc_idx]
          net_hh_col <- intersect(c("QHNUMBER", "hv002"), colnames(nets_df))[1]
          per_hh_col <- intersect(c("QHNUMBER", "hv002"), colnames(itn_df))[1]
          
          itn_results <- compute_itn_from_persons_and_nets(
            persons_df         = itn_df,
            nets_df            = nets_df,
            person_cluster_col = cid,
            person_hh_col      = per_hh_col,
            person_line_col    = itn_line_col,
            person_age_col     = age_col,
            net_cluster_col    = net_cid,
            net_hh_col         = net_hh_col,
            age_groups         = age_groups
          )
          for (grp in names(age_groups)) {
            if (!is.null(itn_results[[grp]]))
              cluster_df <- merge_cluster_result(cluster_df, itn_results[[grp]], cid, grp)
          }
          
        } else {
          # --- Classic path (pre-computed per-person column, e.g. old DHS .dta) ---
          # derive_itn_col handles the rare case where the column is absent and
          # must be inferred from suffix columns within the same file.
          if (!is.null(itn_line_col))
            itn_df <- derive_itn_col(itn_df, recode_df$code[itn_idx], itn_line_col,
                                     recode_df$pos_pattern[itn_idx],
                                     recode_df$neg_pattern[itn_idx])
          
          itn_df <- code_positivity(itn_df, recode_df$code[itn_idx],
                                    recode_df$pos_pattern[itn_idx],
                                    recode_df$neg_pattern[itn_idx])
          
          itn_age_vals <- parse_dhs_age(itn_df[[age_col]], age_col)
          for (grp in names(age_groups)) {
            lo <- age_groups[[grp]][1]; hi <- age_groups[[grp]][2]
            if (is.infinite(lo) && is.infinite(hi)) {
              sub <- itn_df
            } else {
              sub <- itn_df[!is.na(itn_age_vals) & itn_age_vals > lo & itn_age_vals <= hi, ]
            }
            res <- aggregate_to_clusters(sub, cid, include_weight = (grp == "itn_all"))
            cluster_df <- merge_cluster_result(cluster_df, res, cid, grp)
          }
        }
        
        if (!"mean_date" %in% colnames(cluster_df)) {
          dr <- extract_cluster_dates(itn_df, cid, recode_df$month_code[itn_idx],
                                      recode_df$year_code[itn_idx])
          if (!is.null(dr)) cluster_df <- merge(cluster_df, dr,
                                                by.x="clusterid", by.y=cid, all=TRUE)
        }
      }
      
      # Other binary indicators
      cluster_df <- add_var(cluster_df, "iptp",              recode_df)
      cluster_df <- add_var(cluster_df, "cm",                recode_df)
      cluster_df <- add_var(cluster_df, "received_treatment", recode_df)
      cluster_df <- add_var(cluster_df, "blood_test",         recode_df)
      
      # Treatment-seeking and ACT (loaded once, derived columns built)
      treat_idx <- which(recode_df$variable == "sought_treatment")
      art_idx   <- which(recode_df$variable == "art_given_antimal")
      
      cur_art     <- resolve_codes(art_codes,     year)
      cur_non_art <- resolve_codes(non_art_codes, year)
      
      if (length(treat_idx) > 0 && !is.na(recode_df$filename[treat_idx])) {
        treat_path <- build_file_path(dta_dir, recode_df$folder_dir[treat_idx],
                                      recode_df$filename[treat_idx])
        treat_df <- load_survey_file(treat_path)
        treat_df <- build_treatment_columns(treat_df, year)
        withCallingHandlers(
          treat_df <- build_act_columns(treat_df, cur_art, cur_non_art),
          warning = function(w) {
            warning("Year ", year, ": ", conditionMessage(w), call. = FALSE)
            invokeRestart("muffleWarning")
          }
        )
        cid      <- recode_df$cluster_id_code[treat_idx]
        
        # Single-sector subset for ACT calculations
        act_df <- treat_df %>%
          filter(!is.na(public_treatment), !is.na(private_formal_treatment),
                 !is.na(private_informal_treatment)) %>%
          mutate(n_sectors = public_treatment + private_formal_treatment +
                   private_informal_treatment) %>%
          filter(n_sectors == 1)
        
        # Treatment-seeking rates by sector
        for (tv in c("sought_treatment","public_treatment",
                     "private_formal_treatment","private_informal_treatment")) {
          treat_df$pos <- treat_df[[tv]]
          res <- aggregate_to_clusters(treat_df, cid)
          cluster_df <- merge_cluster_result(cluster_df, res, cid, tv)
        }
        
        # ACT rate columns → output names
        # art_given_antimalarial (overall ACT rate) is aggregated from treat_df —
        # it does not require sector assignment and must not be filtered by act_df
        # (which would be empty when sector columns are unavailable, e.g. MIS 2025).
        # Sector-disaggregated columns require act_df (single-sector rows only).
        overall_act_map <- c(art_given_antimalarial = "art_high")
        sector_act_map  <- c(
          act_given_antimal_public           = "art_high_public",
          act_given_antimal_private_formal   = "art_high_private_formal",
          act_given_antimal_private_informal = "art_high_private_informal",
          act_given_treat_public             = "art_low_public",
          act_given_treat_private_formal     = "art_low_private_formal",
          act_given_treat_private_informal   = "art_low_private_informal"
        )
        for (src_col in names(overall_act_map)) {
          if (!src_col %in% colnames(treat_df)) next
          treat_df$pos <- treat_df[[src_col]]
          res <- aggregate_to_clusters(treat_df, cid)
          cluster_df <- merge_cluster_result(cluster_df, res, cid, overall_act_map[[src_col]])
        }
        for (src_col in names(sector_act_map)) {
          if (!src_col %in% colnames(act_df)) next
          act_df$pos <- act_df[[src_col]]
          res <- aggregate_to_clusters(act_df, cid)
          cluster_df <- merge_cluster_result(cluster_df, res, cid, sector_act_map[[src_col]])
        }
        
      } else if (length(art_idx) > 0 && !is.na(recode_df$filename[art_idx])) {
        # Fallback: no sector breakdown available — use aggregate ACT rate
        cluster_df <- add_var(cluster_df, "art_given_antimal", recode_df)
      }
      
      # Assign clusters to admin units via spatial join.
      # When GPS coordinates are present (all paths except init_cluster_df_from_lga),
      # the spatial join is always performed using admin_shape — admin_code in the
      # recode CSV is used only as a cross-check inside extract_gps_from_nets, not
      # to bypass this step.
      # The spatial join is skipped only when there are no GPS coordinates (i.e.,
      # init_cluster_df_from_lga was used because the recode CSV had no locations
      # row), in which case NOMDEP/NOMREGION were read directly from the survey data
      # and a spatial join is not possible.
      if (all(c("latitude", "longitude") %in% colnames(cluster_df))) {
        if (is.null(admin_shape))
          stop("admin_shape is required to assign clusters to admin units via spatial join.")
        .join <- assign_clusters_to_admins(cluster_df, admin_shape)
        cluster_df <- .join$cluster_df
        if (!is.null(.join$ambiguous_log)) {
          amb_csv <- file.path(plots_dir, paste0("admin_assignment_ambiguity_", year, ".csv"))
          amb_pdf <- file.path(plots_dir, paste0("admin_assignment_ambiguity_", year, ".pdf"))
          write.csv(.join$ambiguous_log, amb_csv, row.names = FALSE)
          message("  [admin_assign] Ambiguity log (", nrow(.join$ambiguous_log),
                  " cluster(s)) saved to ", basename(amb_csv))
          .plot_ambiguous_admin_assignments(.join$ambiguous_log, admin_shape, amb_pdf)
        } else {
          message("  [admin_assign] All clusters matched exactly one polygon.")
        }
      } else if (!all(c("NOMDEP", "NOMREGION") %in% colnames(cluster_df))) {
        stop("No GPS coordinates and no admin names available — cannot assign clusters to admin units.")
      }
      write.csv(format_dates_for_csv(cluster_df), cluster_file)
    }
    
    # ------------------------------------------------------------------ #
    # STEP 2 — Admin-level aggregation                                    #
    # ------------------------------------------------------------------ #
    cluster_df <- cluster_df[!is.na(cluster_df$NOMREGION), ]
    
    present_vars <- variables[paste0(variables, "_num_total") %in% colnames(cluster_df)]
    admin_sums  <- compute_level_sums(cluster_df, c("NOMREGION","NOMDEP"), present_vars)
    region_sums <- compute_level_sums(cluster_df, "NOMREGION", present_vars)
    
    if ("mean_date" %in% colnames(cluster_df)) {
      cluster_df$mean_date <- .parse_csv_date(cluster_df$mean_date)
      admin_sums  <- merge(admin_sums,
                           cluster_df %>% group_by(NOMREGION, NOMDEP) %>%
                             summarise(mean_date = mean(mean_date, na.rm=TRUE), .groups="drop"),
                           all.x = TRUE)
      region_sums <- merge(region_sums,
                           cluster_df %>% group_by(NOMREGION) %>%
                             summarise(mean_date = mean(mean_date, na.rm=TRUE), .groups="drop"),
                           all.x = TRUE)
    }
    write.csv(format_dates_for_csv(as.data.frame(dplyr::select(admin_sums, -matches("^art_(high|low)")))),
              file.path(out_dir, paste0("DHS_admin_", year, "_0.csv")))
    
    # Expand to include all admins; merge archetype info
    archetype_info <- normalize_archetype_info(read.csv(ds_pop_df_filename))
    admin_sums     <- expand_and_merge_archetype(admin_sums, archetype_info)
    write.csv(format_dates_for_csv(as.data.frame(dplyr::select(admin_sums, -matches("^art_(high|low)")))),
              file.path(out_dir, paste0("DHS_admin_", year, ".csv")))
    
    grouping_col <- if ("ZONE" %in% colnames(archetype_info)) "ZONE" else "Archetype"

    # Bring grouping column into cluster_df for ACT aggregation (NOMREGION → ZONE/Archetype).
    # .safe_zone_map deduplicates to one row per state — merging a non-deduplicated table
    # would create one cluster row per archetype, inflating all arch/national sums.
    zone_map   <- .safe_zone_map(archetype_info, grouping_col)
    cluster_df <- merge(cluster_df, zone_map, by = "NOMREGION", all.x = TRUE)
    
    # Bring seasonality_archetype into cluster_df for archetype-level aggregation.
    # This is an LGA-level property so must be merged by NOMDEP, not NOMREGION.
    if ("seasonality_archetype" %in% colnames(archetype_info) &&
        "NOMDEP" %in% colnames(archetype_info) &&
        "NOMDEP" %in% colnames(cluster_df)) {
      sa_map_lga <- dplyr::distinct(archetype_info[, c("NOMDEP", "seasonality_archetype")])
      cluster_df <- merge(cluster_df, sa_map_lga, by = "NOMDEP", all.x = TRUE)
    }
    
    # ------------------------------------------------------------------ #
    # STEP 3 (part A) — ACT rate aggregation                             #
    # Must run before art columns are dropped from cluster_df below.     #
    # ------------------------------------------------------------------ #
    has_sector_act <- any(grepl("^art_high_public_num_total", colnames(cluster_df)))
    
    if (has_sector_act) {
      national_act_raw <- aggregate_act_rates(cluster_df)  # keep art_ prefix for use as fallback
      national_act <- .rename_art_to_act(national_act_raw)
      if ("act_high_rate" %in% colnames(national_act))
        national_act <- dplyr::rename(national_act, act_high_rate_all_sectors = act_high_rate)
      national_act <- .add_act_middle_rates(national_act, cm_high_estimate_weight)
      
      # Always write zone/archetype file with a National row appended for comparison.
      if (grouping_col %in% colnames(cluster_df)) {
        zone_act <- .rename_art_to_act(aggregate_act_rates(cluster_df, grouping_col, min_num_total, national_act_raw))
        if ("act_high_rate" %in% colnames(zone_act))
          zone_act <- dplyr::rename(zone_act, act_high_rate_all_sectors = act_high_rate)
        zone_act <- .add_act_middle_rates(zone_act, cm_high_estimate_weight)
        national_row <- national_act
        national_row[[grouping_col]] <- "National"
        write.csv(dplyr::bind_rows(zone_act, national_row),
                  file.path(out_dir, paste0("DHS_ACT_prob_estimates_by_sector_and_zone_", year, ".csv")))
      } else {
        write.csv(national_act,
                  file.path(out_dir, paste0("DHS_ACT_prob_estimates_by_sector_and_zone_", year, ".csv")))
      }
      
      # Choose which level to use for downstream effective treatment calculation.
      act_rates <- if (!use_national_act && grouping_col %in% colnames(cluster_df)) zone_act else national_act
    }
    
    cluster_df <- cluster_df %>% dplyr::select(-matches("^art"))  # exclude raw ACT cluster cols
    
    # Archetype-level aggregation uses seasonality_archetype if available (LGA-level
    # grouping merged in above), otherwise falls back to the zone/archetype grouping.
    arch_col  <- if ("seasonality_archetype" %in% colnames(cluster_df)) "seasonality_archetype"
    else grouping_col
    arch_sums     <- compute_level_sums(cluster_df, arch_col, present_vars)
    national_sums <- compute_level_sums(cluster_df, character(0), present_vars)
    # Zone-level sums for the fallback mechanism: always indexed by grouping_col (ZONE/Archetype).
    # arch_sums may use seasonality_archetype (LGA-level), which has no ZONE column and would
    # cause apply_min_sample_fallback to skip the zone level and fall through to national.
    zone_fallback_sums <- if (arch_col != grouping_col) {
      compute_level_sums(cluster_df, grouping_col, present_vars)
    } else {
      arch_sums
    }
    
    if ("mean_date" %in% colnames(cluster_df)) {
      arch_sums <- merge(arch_sums,
                         cluster_df %>% group_by(across(all_of(arch_col))) %>%
                           summarise(mean_date = mean(mean_date, na.rm=TRUE), .groups="drop"),
                         all.x = TRUE)
      national_sums$mean_date <- mean(cluster_df$mean_date, na.rm=TRUE)
    }
    
    # Write archetype- and national-level rate files from in-memory aggregates.
    # These use the same counts as the fallback mechanism (raw cluster sums before
    # any minN substitution), so they are unbiased totals.
    arch_out <- arch_sums
    arch_out$year <- year
    existing_arch <- file.path(out_dir, "DHS_archetype_rates.csv")
    if (file.exists(existing_arch)) {
      prev_arch <- read.csv(existing_arch)
      prev_arch <- prev_arch[prev_arch$year != year, ]
      if ("mean_date" %in% colnames(prev_arch))
        prev_arch$mean_date <- .parse_csv_date(prev_arch$mean_date)
      arch_out <- dplyr::bind_rows(prev_arch, arch_out)
    }
    write.csv(format_dates_for_csv(as.data.frame(arch_out)),
              file.path(out_dir, "DHS_archetype_rates.csv"), row.names = FALSE)
    
    national_out <- national_sums
    national_out$year <- year
    existing_nat <- file.path(out_dir, "DHS_national_rates.csv")
    if (file.exists(existing_nat)) {
      prev <- read.csv(existing_nat)
      prev <- prev[prev$year != year, ]
      if ("mean_date" %in% colnames(prev))
        prev$mean_date <- .parse_csv_date(prev$mean_date)
      national_out <- dplyr::bind_rows(prev, national_out)
    }
    write.csv(format_dates_for_csv(as.data.frame(national_out)),
              file.path(out_dir, "DHS_national_rates.csv"), row.names = FALSE)
    
    # ZONE is already in admin_sums from expand_and_merge_archetype — no merge needed.
    fallback_list <- list(
      list(join_col = "NOMREGION",  data = region_sums),
      list(join_col = grouping_col, data = zone_fallback_sums),
      list(join_col = NULL,         data = national_sums)
    )
    admin_sums <- apply_min_sample_fallback(admin_sums, fallback_list,
                                            present_vars, min_num_total)
    
    # ------------------------------------------------------------------ #
    # STEP 3 (part B) — Apply ACT rates to admin_sums                   #
    # ------------------------------------------------------------------ #
    if (has_sector_act) {
      admin_sums <- admin_sums %>%
        dplyr::select(-matches("^art_(high|low)")) %>%
        merge(act_rates, all.x = TRUE) %>%
        mutate(
          effective_treatment_low_rate     = public_treatment_rate * act_low_public_rate +
            private_formal_treatment_rate * act_low_private_formal_rate +
            private_informal_treatment_rate * act_low_private_informal_rate,
          effective_treatment_middle_rate  = public_treatment_rate * act_middle_rate_public +
            private_formal_treatment_rate * act_middle_rate_private_formal +
            private_informal_treatment_rate * act_middle_rate_private_informal,
          effective_treatment_high_rate    = public_treatment_rate * act_high_public_rate +
            private_formal_treatment_rate * act_high_private_formal_rate +
            private_informal_treatment_rate * act_high_private_informal_rate,
          effective_treatment_low_num_total    = public_treatment_num_total +
            private_formal_treatment_num_total + private_informal_treatment_num_total,
          effective_treatment_middle_num_total = effective_treatment_low_num_total,
          effective_treatment_high_num_total   = effective_treatment_low_num_total,
          effective_treatment_low_num_true     = effective_treatment_low_num_total    * effective_treatment_low_rate,
          effective_treatment_middle_num_true  = effective_treatment_middle_num_total * effective_treatment_middle_rate,
          effective_treatment_high_num_true    = effective_treatment_high_num_total   * effective_treatment_high_rate
        )
      
    } else if (all(c("cm_rate","blood_test_rate") %in% colnames(admin_sums))) {
      admin_sums <- admin_sums %>% mutate(
        effective_treatment_rate      = (cm_rate + blood_test_rate) / 2,
        effective_treatment_num_total = (cm_num_total + blood_test_num_total) / 2,
        effective_treatment_num_true  = (cm_num_true  + blood_test_num_true)  / 2)
    }
    
    write.csv(format_dates_for_csv(as.data.frame(admin_sums)),
              file.path(out_dir, paste0("DHS_admin_minN", min_num_total, "_", year, ".csv")))
    message("  Done — outputs written for year ", year)
  }
}


# =============================================================================
# SECTION 10: COUNTRY, STATE, AND ARCHETYPE LEVEL EXTRACTION
# Read cluster-level CSVs produced by extract_DHS_data and aggregate to national,
# state, archetype, and admin levels with optional survey weighting.
# =============================================================================

# Compute national-level rates from DHS_cluster_outputs_{year}.csv.
# Outputs: DHS_country_weighted_rates.csv and DHS_country_unweighted_rates.csv.
extract_country_level_DHS_data <- function(hbhi_dir, dta_dir, years) {
  # Reads DHS_cluster_outputs_{year}.csv (written by extract_DHS_data) and computes
  # national-level rates by summing _num_true / _num_total across all clusters.
  # Weighted rates use cluster-level survey weights (itn_weights column, if present).
  # This avoids re-reading raw survey files and correctly handles all derived variables
  # (treatment sectors, ACT rates) that are already computed in the cluster-level output.
  out_dir <- file.path(hbhi_dir, "estimates_from_DHS")
  
  out_wt  <- file.path(out_dir, "DHS_country_weighted_rates.csv")
  out_unw <- file.path(out_dir, "DHS_country_unweighted_rates.csv")
  
  # Seed from existing files, dropping any rows whose year is being recomputed
  read_existing <- function(path) {
    if (!file.exists(path)) return(data.frame())
    prev <- read.csv(path) %>% dplyr::select(-any_of("X"))
    prev[!prev$year %in% years, , drop = FALSE]
  }
  weighted_all   <- read_existing(out_wt)
  unweighted_all <- read_existing(out_unw)
  
  for (year in years) {
    cluster_file <- file.path(out_dir, paste0("DHS_cluster_outputs_", year, ".csv"))
    if (!file.exists(cluster_file)) {
      warning("extract_country_level_DHS_data: cluster file not found for year ", year,
              " — skipping. Run extract_DHS_data first. (", cluster_file, ")")
      next
    }
    cl <- read.csv(cluster_file) %>% dplyr::select(-any_of("X"))
    
    # Identify variable names from _num_total columns; skip raw art_ intermediates
    tot_cols <- grep("_num_total$", colnames(cl), value = TRUE)
    tot_cols <- tot_cols[!grepl("^art_", tot_cols)]
    vars <- sub("_num_total$", "", tot_cols)
    
    has_weight <- "itn_weights" %in% colnames(cl)
    
    wt_row  <- data.frame(admin = "country", year = year)
    unw_row <- data.frame(admin = "country", year = year)
    
    for (v in vars) {
      tc <- paste0(v, "_num_true");  oc <- paste0(v, "_num_total")
      rc <- paste0(v, "_rate")
      if (!tc %in% colnames(cl)) next
      
      # Unweighted: simple ratio of summed counts across clusters
      n_true  <- sum(cl[[tc]], na.rm = TRUE)
      n_total <- sum(cl[[oc]], na.rm = TRUE)
      unw_row[[rc]] <- if (n_total > 0) n_true / n_total else NA_real_
      
      # Weighted: weight each cluster's contribution by itn_weights * num_total
      if (has_weight) {
        sub <- cl[!is.na(cl[[oc]]) & !is.na(cl$itn_weights) & cl[[oc]] > 0, ]
        eff_wt <- sub$itn_weights * sub[[oc]]
        wt_row[[rc]] <- if (sum(eff_wt) > 0)
          sum(sub[[tc]] * sub$itn_weights) / sum(sub[[oc]] * sub$itn_weights)
        else NA_real_
      } else {
        wt_row[[rc]] <- unw_row[[rc]]  # no weights available — fall back to unweighted
      }
    }
    
    weighted_all   <- dplyr::bind_rows(weighted_all,   wt_row)
    unweighted_all <- dplyr::bind_rows(unweighted_all, unw_row)
  }
  
  weighted_all   <- weighted_all[order(weighted_all$year),   , drop = FALSE]
  unweighted_all <- unweighted_all[order(unweighted_all$year), , drop = FALSE]
  
  write.csv(as.data.frame(weighted_all),   out_wt,  row.names = FALSE)
  write.csv(as.data.frame(unweighted_all), out_unw, row.names = FALSE)
}


# Compute state-level rates from DHS_cluster_outputs_{year}.csv files.
# Outputs: DHS_state_{year}.csv (per year) and DHS_state_pfpr.csv (combined, wide).
# Weighting: cluster weight x sample size; falls back to unweighted if absent.

extract_state_level_DHS_data <- function(hbhi_dir, dta_dir, years) {
  out_dir <- file.path(hbhi_dir, "estimates_from_DHS")
  
  state_col_candidates <- c("NOMREGION", "state", "region", "State", "Region")
  
  # Helper: compute weighted rate for one variable across clusters in a group
  .state_rate <- function(sub, tc, oc, has_weight) {
    sub <- sub[!is.na(sub[[oc]]) & sub[[oc]] > 0, ]
    if (nrow(sub) == 0) return(NA_real_)
    if (has_weight && any(!is.na(sub$itn_weights))) {
      sub <- sub[!is.na(sub$itn_weights), ]
      if (nrow(sub) == 0) return(NA_real_)
      sum(sub[[tc]] * sub$itn_weights) / sum(sub[[oc]] * sub$itn_weights)
    } else {
      sum(sub[[tc]], na.rm = TRUE) / sum(sub[[oc]], na.rm = TRUE)
    }
  }
  
  pfpr_list <- list()  # accumulates per-year pfpr for combined output
  
  any_written <- FALSE
  
  for (year in years) {
    cluster_file <- file.path(out_dir, paste0("DHS_cluster_outputs_", year, ".csv"))
    if (!file.exists(cluster_file)) {
      warning("extract_state_level_DHS_data: cluster file not found for year ", year,
              " — skipping. Run extract_DHS_data first. (", cluster_file, ")")
      next
    }
    cl <- read.csv(cluster_file) %>% dplyr::select(-any_of("X"))
    
    # Identify state column
    state_col <- state_col_candidates[state_col_candidates %in% colnames(cl)][1]
    if (is.na(state_col)) {
      warning("extract_state_level_DHS_data: no state column found in ", cluster_file,
              " (tried: ", paste(state_col_candidates, collapse = ", "), ") — skipping year ", year)
      next
    }
    
    has_weight <- "itn_weights" %in% colnames(cl)
    
    # All variables from _num_total columns; skip raw art_ intermediates
    tot_cols <- grep("_num_total$", colnames(cl), value = TRUE)
    tot_cols <- tot_cols[!grepl("^art_", tot_cols)]
    vars <- sub("_num_total$", "", tot_cols)
    
    states <- sort(unique(na.omit(cl[[state_col]])))
    
    # ---- per-year state file ----
    yr_rows <- lapply(states, function(st) {
      sub <- cl[!is.na(cl[[state_col]]) & cl[[state_col]] == st, , drop = FALSE]
      row <- data.frame(state = st, stringsAsFactors = FALSE)
      for (v in vars) {
        tc <- paste0(v, "_num_true");  oc <- paste0(v, "_num_total")
        if (!tc %in% colnames(cl)) next
        row[[paste0(v, "_num_true")]]  <- sum(sub[[tc]], na.rm = TRUE)
        row[[paste0(v, "_num_total")]] <- sum(sub[[oc]], na.rm = TRUE)
        row[[paste0(v, "_rate")]]      <- .state_rate(sub, tc, oc, has_weight)
      }
      row
    })
    yr_df <- dplyr::bind_rows(yr_rows)
    
    out_yr <- file.path(out_dir, paste0("DHS_state_", year, ".csv"))
    write.csv(yr_df, out_yr, row.names = FALSE)
    message("extract_state_level_DHS_data: written ", out_yr)
    any_written <- TRUE
    
    # ---- accumulate pfpr columns for combined file ----
    for (pfpr_var in c("mic", "rdt")) {
      rc <- paste0(pfpr_var, "_rate")
      if (!rc %in% colnames(yr_df)) next
      col_name <- paste0(year, "_", pfpr_var, "_rate")
      tmp <- yr_df[, c("state", rc), drop = FALSE]
      colnames(tmp)[2] <- col_name
      if (is.null(pfpr_list[[pfpr_var]])) {
        pfpr_list[[pfpr_var]] <- tmp
      } else {
        pfpr_list[[pfpr_var]] <- merge(pfpr_list[[pfpr_var]], tmp, by = "state", all = TRUE)
      }
    }
  }
  
  if (!any_written) {
    warning("extract_state_level_DHS_data: no data found — no output written.")
    return(invisible(NULL))
  }
  
  # ---- combined pfpr file (wide across years) ----
  combined <- pfpr_list[["mic"]]
  if (!is.null(pfpr_list[["rdt"]])) {
    combined <- if (is.null(combined)) pfpr_list[["rdt"]] else
      merge(combined, pfpr_list[["rdt"]], by = "state", all = TRUE)
  }
  if (!is.null(combined)) {
    # Reorder columns: state | year1_mic | year1_rdt | year2_mic | year2_rdt | ...
    rate_cols <- colnames(combined)[colnames(combined) != "state"]
    ordered_cols <- c()
    for (yr in sort(years)) {
      for (pfpr_var in c("mic", "rdt")) {
        col <- paste0(yr, "_", pfpr_var, "_rate")
        if (col %in% rate_cols) ordered_cols <- c(ordered_cols, col)
      }
    }
    combined <- combined[, c("state", ordered_cols), drop = FALSE]
    combined <- combined[order(combined$state), ]
    
    out_pfpr <- file.path(out_dir, "DHS_state_pfpr.csv")
    write.csv(combined, out_pfpr, row.names = FALSE)
    message("extract_state_level_DHS_data: written ", out_pfpr)
    invisible(combined)
  }
}


# Aggregate admin-level counts from DHS_admin_{year}.csv to archetype level.
# Uses raw (pre-fallback) counts to avoid inflated archetype totals.
# Outputs: DHS_archetype_rates.csv (combined across all years).
extract_archetype_level_DHS_data <- function(hbhi_dir, dta_dir,
                                             ds_pop_df_filename, years,
                                             min_num_total = 30) {
  out_dir        <- file.path(hbhi_dir, "estimates_from_DHS")
  archetype_info <- read.csv(ds_pop_df_filename) %>%
    normalize_archetype_info() %>%
    dplyr::select(admin_name = NOMDEP, seasonality_archetype, Archetype)
  
  result_all <- data.frame()
  for (year in years) {
    # Read the plain admin file (before apply_min_sample_fallback).
    # The minN file inflates _num_true/_num_total for small admins (replacing them with
    # region/archetype/national counts), so summing across admins in an archetype produces
    # falsely large totals.  The plain admin file has the actual raw counts, matching
    # what the original script used.
    admin_file <- file.path(out_dir, paste0("DHS_admin_", year, ".csv"))
    if (!file.exists(admin_file)) {
      warning("extract_archetype_level_DHS_data: file not found for year ", year,
              " — skipping. (", admin_file, ")")
      next
    }
    counts <- read.csv(admin_file) %>%
      # Drop the unnamed row-index column written by write.csv (read as "X")
      dplyr::select(-any_of(c("X"))) %>%
      dplyr::select(-matches("_rate$")) %>%
      merge(archetype_info %>% rename(NOMDEP = admin_name)) %>%
      dplyr::select(-NOMDEP) %>%
      group_by(seasonality_archetype, Archetype) %>%
      summarise(across(where(is.numeric), ~ sum(.x, na.rm=TRUE)), .groups="drop") %>%
      mutate(year = year)
    
    # Recompute rates from aggregated counts
    vars <- sub("_num_total$", "", grep("_num_total$", colnames(counts), value=TRUE))
    for (v in vars) {
      tc <- paste0(v,"_num_true"); oc <- paste0(v,"_num_total"); rc <- paste0(v,"_rate")
      if (tc %in% colnames(counts) && oc %in% colnames(counts))
        counts[[rc]] <- counts[[tc]] / counts[[oc]]
    }
    result_all <- dplyr::bind_rows(result_all, counts)
  }
  write.csv(as.data.frame(result_all),
            file.path(out_dir, "DHS_archetype_rates.csv"), row.names=FALSE)
}


# =============================================================================
# SECTION 11: PLOTTING
# Choropleth and cluster-point maps of extracted DHS/MIS estimates, including
# state-level and admin-level choropleth maps for all requested variables.
# =============================================================================

# Plot choropleth maps of state-level rates for requested variables and years.
# Inputs: hbhi_dir, year, variables, admin_shape (LGA-level, dissolved to state).
# Outputs: one PNG per variable/weight combination in {hbhi_dir}/estimates_from_DHS/plots/.

plot_state_level_DHS_data <- function(hbhi_dir, year, variables, admin_shape,
                                      weighted = c(TRUE, FALSE)) {
  for (pkg in c("ggplot2", "sf", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("plot_state_level_DHS_data requires the '", pkg, "' package.")
  }
  
  out_dir   <- file.path(hbhi_dir, "estimates_from_DHS")
  plots_dir <- file.path(out_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  state_file <- file.path(out_dir, paste0("DHS_state_", year, ".csv"))
  if (!file.exists(state_file))
    stop("plot_state_level_DHS_data: file not found — run extract_state_level_DHS_data first.\n  ",
         state_file)
  rates_df <- read.csv(state_file, stringsAsFactors = FALSE)
  
  # Dissolve LGA shapefile to state level using NOMREGION
  state_sf <- sf::st_as_sf(admin_shape) %>%
    dplyr::group_by(NOMREGION) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    dplyr::rename(state = NOMREGION)
  
  # Centroids computed once (shared across all variable/weight combinations)
  centroids <- suppressWarnings(sf::st_centroid(state_sf))
  ctrd_coords <- sf::st_coordinates(centroids)
  
  # Fuzzy-join helper: match state names case-insensitively after stripping punctuation
  norm_str <- function(x) tolower(trimws(gsub("[^a-zA-Z0-9 ]", "", as.character(x))))
  sf_states   <- norm_str(state_sf$state)
  rate_states <- norm_str(rates_df$state)
  match_idx   <- match(sf_states, rate_states)
  unmatched   <- state_sf$state[is.na(match_idx)]
  if (length(unmatched) > 0)
    warning("plot_state_level_DHS_data: no rate data for state(s): ",
            paste(unmatched, collapse = ", "))
  
  for (var in variables) {
    rc_weighted <- paste0(var, "_rate")
    tc <- paste0(var, "_num_true");  oc <- paste0(var, "_num_total")
    
    if (!rc_weighted %in% colnames(rates_df)) {
      warning("plot_state_level_DHS_data: column '", rc_weighted, "' not found in ",
              state_file, " — skipping variable '", var, "'.")
      next
    }
    
    # Compute unweighted rate from counts if columns are present
    has_unweighted <- tc %in% colnames(rates_df) && oc %in% colnames(rates_df)
    unweighted_vals <- if (has_unweighted)
      rates_df[[tc]] / rates_df[[oc]]
    else
      rep(NA_real_, nrow(rates_df))
    
    for (wt in weighted) {
      if (!wt && !has_unweighted) {
        message("plot_state_level_DHS_data: num_true/num_total columns absent for '", var,
                "' — cannot produce unweighted map, skipping.")
        next
      }
      
      rate_vals <- if (wt) rates_df[[rc_weighted]] else unweighted_vals
      wt_label  <- if (wt) "cluster-weighted" else "unweighted"
      wt_suffix <- if (wt) "_weighted"        else "_unweighted"
      
      plot_sf <- state_sf
      plot_sf$rate_col <- rate_vals[match_idx]
      
      lbl_df <- data.frame(
        state    = plot_sf$state,
        rate_val = plot_sf$rate_col,
        x        = ctrd_coords[, 1],
        y        = ctrd_coords[, 2],
        stringsAsFactors = FALSE
      )
      lbl_df$label <- paste0(
        lbl_df$state, "\n",
        ifelse(is.na(lbl_df$rate_val), "NA", sprintf("%.3f", lbl_df$rate_val))
      )
      lbl_df$text_col <- ifelse(
        !is.na(lbl_df$rate_val) & lbl_df$rate_val > 0.55, "white", "grey10"
      )
      
      p <- ggplot2::ggplot(plot_sf) +
        ggplot2::geom_sf(ggplot2::aes(fill = rate_col), colour = "white", linewidth = 0.3) +
        ggplot2::scale_fill_gradient(
          low      = "#cce5ff",
          high     = "#003580",
          limits   = c(0, 1),
          na.value = "grey85",
          name     = "Rate"
        ) +
        ggplot2::geom_text(
          data    = lbl_df,
          mapping = ggplot2::aes(x = x, y = y, label = label, colour = text_col),
          size    = 2.2, lineheight = 0.9, fontface = "plain", show.legend = FALSE
        ) +
        ggplot2::scale_colour_identity() +
        ggplot2::labs(
          title    = paste0(var, " — ", year),
          subtitle = paste0("State-level ", wt_label, " rate")
        ) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
          legend.position = "right"
        )
      
      out_file <- file.path(plots_dir,
                            paste0("DHS_state_", var, "_", year, wt_suffix, ".png"))
      ggplot2::ggsave(out_file, p, width = 8, height = 6, dpi = 200)
      message("plot_state_level_DHS_data: written ", out_file)
    }
  }
}


# Plot cluster-level and admin-level choropleth maps for all requested variables.
# Inputs: hbhi_dir, years, admin_shape; optional layout/colour/filter arguments.
# Outputs: PNGs (and optionally PDFs) in {hbhi_dir}/estimates_from_DHS/plots/.
plot_extracted_DHS_data <- function(
    hbhi_dir, years, admin_shape,
    min_num_total  = 30,
    variables      = c("mic","rdt","itn_all","itn_u5","itn_5_10","itn_10_15",
                       "itn_15_20","itn_o20","iptp","cm","received_treatment",
                       "sought_treatment","blood_test"),
    colors_range_0_to_1    = NA,
    all_years_int_plot_panel  = TRUE,
    separate_plots_for_each_var = TRUE,
    plot_separate_pdfs    = FALSE,
    plot_vaccine          = FALSE,
    plot_suffix           = ""
) {
  out_dir   <- file.path(hbhi_dir, "estimates_from_DHS")
  plots_dir <- file.path(out_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir)
  
  vacc_str <- if (plot_vaccine) "vaccine_" else ""
  if (any(is.na(colors_range_0_to_1)))
    colors_range_0_to_1 <- add.alpha(pals::parula(101), alpha=0.5)
  
  vars_panel <- if (!plot_vaccine)
    variables[variables %in% c("mic","rdt","itn_all","itn_u5","iptp","cm",
                               "received_treatment","blood_test",  # ,"sought_treatment"
                               "art", "art_middle_public","effective_treatment_middle")]
  else variables
  
  nyears <- length(years)
  
  # ---- Load all CSVs once, cached by year ----------------------------------------
  # (avoids re-reading the same file once per variable × year, which was the main
  #  speed bottleneck when plotting many variables across many years)
  clust_cache <- setNames(lapply(years, function(y) {
    f <- file.path(out_dir, paste0("DHS_", vacc_str, "cluster_outputs_", y, ".csv"))
    if (!file.exists(f)) { warning("Cluster file not found: ", f); return(NULL) }
    d <- read.csv(f)[, -1]
    d$latitude[d$latitude   == 0] <- NA
    d$longitude[d$longitude == 0] <- NA
    d
  }), as.character(years))
  
  admin_cache <- setNames(lapply(years, function(y) {
    f <- file.path(out_dir, paste0("DHS_", vacc_str, "admin_minN",
                                   min_num_total, "_", y, ".csv"))
    if (!file.exists(f)) { warning("Admin file not found: ", f); return(NULL) }
    read.csv(f)[, -1]
  }), as.character(years))
  
  # Precompute admin polygon reorder index once per year (reused across all variables)
  admin_reorder <- setNames(lapply(as.character(years), function(y) {
    d <- admin_cache[[y]]
    if (is.null(d)) return(NULL)
    match(sapply(admin_shape$NOMDEP, normalize_admin_name),
          sapply(d$NOMDEP,           normalize_admin_name))
  }), as.character(years))
  
  # Cluster plots can only show variables whose _num_total column exists in at
  # least one year's cluster CSV (e.g. effective_treatment_* are admin-only,
  # computed as sector-rate × sector-ACT products during admin aggregation).
  vars_cluster <- vars_panel[sapply(vars_panel, function(v) {
    oc <- paste0(v, "_num_total")
    any(sapply(clust_cache, function(d) !is.null(d) && oc %in% colnames(d)))
  })]
  vars_admin_only <- setdiff(vars_panel, vars_cluster)
  if (length(vars_admin_only) > 0)
    message("plot_extracted_DHS_data: the following variables are admin-level only ",
            "(not in cluster CSV) — shown in admin maps only: ",
            paste(vars_admin_only, collapse = ", "))
  
  # Helper: polygon fill colours for admin map (returns NULL if data missing)
  admin_fill_cols <- function(year, rc) {
    d  <- admin_cache[[as.character(year)]]
    ro <- admin_reorder[[as.character(year)]]
    if (is.null(d) || is.null(ro) || !rc %in% colnames(d)) return(NULL)
    vals <- d[[rc]][ro]
    cols <- rep("grey80", length(vals))
    ok   <- !is.na(vals) & vals >= 0 & vals <= 1
    cols[ok] <- colors_range_0_to_1[1 + round(vals[ok] * 100)]
    cols
  }
  
  # Layout: each map panel 3 columns wide, two legend slots on the right
  base_layout <- matrix(c(rep(1:nyears, each=3), nyears+1,
                          rep(1:nyears, each=3), nyears+2), nrow=2, byrow=TRUE)
  
  # Helpers: legends
  draw_legends <- function(var_label, max_n) {   # cluster: colour + size
    legend_img <- as.raster(matrix(
      rev(colors_range_0_to_1[1 + round(seq(0,1,length.out=20)*100)]), ncol=1))
    plot(c(0,2), c(0,1), type="n", axes=FALSE, xlab="", ylab="", main=var_label)
    text(x=1.5, y=seq(0,1,length.out=5), labels=seq(0,1,length.out=5))
    rasterImage(legend_img, 0, 0, 1, 1)
    plot(rep(0,5), seq(1, max_n, length.out=5),
         cex=seq(1, max_n, length.out=5) / round(max_n/5),
         pch=20, axes=FALSE, xlab="", ylab="sample size"); axis(2)
  }
  draw_color_legend <- function(var_label) {     # admin: colour only
    legend_img <- as.raster(matrix(
      rev(colors_range_0_to_1[1 + round(seq(0,1,length.out=20)*100)]), ncol=1))
    plot(c(0,2), c(0,1), type="n", axes=FALSE, xlab="", ylab="", main=var_label)
    text(x=1.5, y=seq(0,1,length.out=5), labels=seq(0,1,length.out=5))
    rasterImage(legend_img, 0, 0, 1, 1)
    plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
  }
  
  # ---- Separate PDFs per year (cluster points) -----------------------------------
  if (plot_separate_pdfs) {
    for (year in years) {
      d <- clust_cache[[as.character(year)]]
      if (is.null(d)) next
      pdf(file.path(plots_dir, paste0("DHS_", vacc_str, "cluster_observations_",
                                      year, ".pdf")),
          width=7, height=5, useDingbats=FALSE)
      par(mfrow=c(1,1))
      for (var in variables) {
        oc <- paste0(var, "_num_total"); rc <- paste0(var, "_rate")
        if (!oc %in% colnames(d)) next
        max_n <- max(d[[oc]], na.rm=TRUE)
        layout(matrix(c(1,1,1,2, 1,1,1,3), nrow=2, byrow=TRUE))
        plot(admin_shape, main=var)
        points(d$longitude, d$latitude,
               col=colors_range_0_to_1[1+round(d[[rc]]*100)], pch=20,
               cex=d[[oc]] / round(max_n/5))
        draw_legends(var, max_n)
      }
      dev.off()
    }
  }
  
  # ---- Multi-year panel plots ----------------------------------------------------
  if (!all_years_int_plot_panel) return(invisible(NULL))
  
  # Stacked layout matrix for the all-vars-in-one case
  stacked_layout <- base_layout
  for (vv in seq_along(vars_panel)[-1])
    stacked_layout <- rbind(stacked_layout,
                            base_layout + (vv - 1L) * max(base_layout))
  
  panel_w <- 0.5 * 7 * nyears
  panel_h <- 0.5 * 6
  
  if (separate_plots_for_each_var) {
    
    # -- Cluster: one PNG per variable (cluster-available variables only) --
    for (var in vars_cluster) {
      oc <- paste0(var, "_num_total"); rc <- paste0(var, "_rate")
      max_n <- max(sapply(as.character(years), function(y) {
        d <- clust_cache[[y]]
        if (!is.null(d) && oc %in% colnames(d)) max(d[[oc]], na.rm=TRUE) else 0
      }))
      if (max_n == 0) next
      png(file.path(plots_dir, paste0("DHS_", var, "_", vacc_str,
                                      "cluster_observations_all_years",
                                      plot_suffix, ".png")),
          width=panel_w, height=panel_h, units="in", res=300)
      layout(base_layout); par(mar=c(0,0,1,0))
      for (year in years) {
        d <- clust_cache[[as.character(year)]]
        if (!is.null(d) && oc %in% colnames(d)) {
          plot(admin_shape, main=paste0(var, " - ", year), border=rgb(.5,.5,.5,.5))
          points(d$longitude, d$latitude, pch=20,
                 col=colors_range_0_to_1[1+round(d[[rc]]*100)],
                 cex=d[[oc]] / round(max_n/5))
        } else plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
      }
      draw_legends(var, max_n); dev.off()
    }
    
    # -- Admin: one PNG per variable --
    for (var in vars_panel) {
      rc <- paste0(var, "_rate")
      png(file.path(plots_dir, paste0("DHS_", var, "_", vacc_str,
                                      "admin_minN", min_num_total,
                                      "_observations_acrossYears",
                                      plot_suffix, ".png")),
          width=panel_w, height=panel_h, units="in", res=300)
      layout(base_layout); par(mar=c(0,0,1,0))
      for (year in years) {
        fills <- admin_fill_cols(year, rc)
        if (!is.null(fills)) {
          plot(admin_shape, main=paste0(var, " - ", year),
               border=rgb(.5,.5,.5,.5), col=fills)
        } else {
          plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
        }
      }
      draw_color_legend(var); dev.off()
    }
    
  } else {
    
    # -- Cluster: all variables in one PNG (cluster-available variables only) --
    stacked_layout_clust <- base_layout
    for (vv in seq_along(vars_cluster)[-1])
      stacked_layout_clust <- rbind(stacked_layout_clust,
                                    base_layout + (vv - 1L) * max(base_layout))
    png(file.path(plots_dir, paste0("DHS_", vacc_str,
                                    "cluster_observations_all_years",
                                    plot_suffix, ".png")),
        width=panel_w, height=panel_h * max(length(vars_cluster), 1L), units="in", res=300)
    layout(stacked_layout_clust); par(mar=c(0,0,1,0))
    for (var in vars_cluster) {
      oc <- paste0(var, "_num_total"); rc <- paste0(var, "_rate")
      max_n <- max(sapply(as.character(years), function(y) {
        d <- clust_cache[[y]]
        if (!is.null(d) && oc %in% colnames(d)) max(d[[oc]], na.rm=TRUE) else 0
      }))
      for (year in years) {
        d <- clust_cache[[as.character(year)]]
        if (!is.null(d) && oc %in% colnames(d)) {
          plot(admin_shape, main=paste0(var, " - ", year), border=rgb(.5,.5,.5,.5))
          points(d$longitude, d$latitude, pch=20,
                 col=colors_range_0_to_1[1+round(d[[rc]]*100)],
                 cex=if (max_n > 0) d[[oc]] / round(max_n/5) else 1)
        } else plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
      }
      draw_legends(var, max(max_n, 1))
    }
    dev.off()
    
    # -- Admin: all variables in one PNG --
    png(file.path(plots_dir, paste0("DHS_", vacc_str,
                                    "admin_minN", min_num_total,
                                    "_observations_all_years",
                                    plot_suffix, ".png")),
        width=panel_w, height=panel_h * length(vars_panel), units="in", res=300)
    layout(stacked_layout); par(mar=c(0,0,1,0))
    for (var in vars_panel) {
      rc <- paste0(var, "_rate")
      for (year in years) {
        fills <- admin_fill_cols(year, rc)
        if (!is.null(fills)) {
          plot(admin_shape, main=paste0(var, " - ", year),
               border=rgb(.5,.5,.5,.5), col=fills)
        } else {
          plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA)
        }
      }
      draw_color_legend(var)
    }
    dev.off()
  }
}


# =============================================================================
# SECTION 12: VACCINE EXTRACTION
# Extract vaccination coverage from DHS and MICS surveys.
# =============================================================================

# ---- Vaccine: DHS -----------------------------------------------------------

# Extract vaccination coverage from a DHS survey year.
# Inputs: hbhi_dir, dta_dir, year, admin_shape, ds_pop_df_filename.
# Outputs: DHS_vaccine_cluster_outputs_{year}.csv, DHS_vaccine_admin_{year}.csv,
#   DHS_vaccine_admin_minN{min_num_total}_{year}.csv.
extract_vaccine_DHS_data <- function(
    hbhi_dir, dta_dir, year, admin_shape, ds_pop_df_filename,
    min_num_total = 30,
    vaccine_variables = c("vacc_dpt1","vacc_dpt2","vacc_dpt3"),
    vaccine_alternate_positive_patterns = c("reported by mother","vaccination marked on card"),
    survey_month_col = "v006", survey_year_col = "v007",
    birth_month_col  = "b1",   birth_year_col  = "b2"
) {
  out_dir   <- file.path(hbhi_dir, "estimates_from_DHS")
  recode_df <- read.csv(file.path(out_dir, paste0("DHS_", year, "_files_recodes_for_sims.csv")))
  
  loc_idx <- which(recode_df$variable == "locations")
  loc_shp <- raster::shapefile(build_file_path(
    dta_dir, recode_df$folder_dir[loc_idx], recode_df$filename[loc_idx]))
  cluster_df <- data.frame(clusterid = loc_shp$DHSCLUST,
                           latitude  = loc_shp$LATNUM, longitude = loc_shp$LONGNUM)
  
  for (vv in vaccine_variables) {
    idx <- which(recode_df$variable == vv)
    if (length(idx) == 0 || is.na(recode_df$filename[idx])) next
    path <- build_file_path(dta_dir, recode_df$folder_dir[idx], recode_df$filename[idx])
    df   <- load_survey_file(path)
    df   <- code_positivity(df, recode_df$code[idx],
                            recode_df$pos_pattern[idx], recode_df$neg_pattern[idx],
                            vaccine_alternate_positive_patterns)
    res  <- aggregate_to_clusters_vacc_dhs(
      df, recode_df$cluster_id_code[idx],
      survey_month_col, survey_year_col, birth_month_col, birth_year_col,
      min_age_months = recode_df$min_age_months_to_include[idx])
    cluster_df <- merge_cluster_result(cluster_df, res, recode_df$cluster_id_code[idx], vv)
  }
  
  cluster_df <- assign_clusters_to_admins(cluster_df, admin_shape)$cluster_df
  write.csv(as.data.frame(cluster_df),
            file.path(out_dir, paste0("DHS_vaccine_cluster_outputs_", year, ".csv")))
  
  cluster_df <- cluster_df[!is.na(cluster_df$NOMREGION), ]
  admin_sums  <- compute_level_sums(cluster_df, c("NOMREGION","NOMDEP"), vaccine_variables)
  
  archetype_info <- normalize_archetype_info(read.csv(ds_pop_df_filename))
  admin_sums     <- expand_and_merge_archetype(admin_sums, archetype_info)
  
  grouping_col  <- if ("ZONE" %in% colnames(archetype_info)) "ZONE" else "Archetype"
  zone_map      <- .safe_zone_map(archetype_info, grouping_col)

  # Merge grouping column into cluster_df so compute_level_sums can group by it.
  cluster_df <- merge(cluster_df, zone_map, by = "NOMREGION", all.x = TRUE)

  region_sums   <- compute_level_sums(cluster_df, "NOMREGION", vaccine_variables)
  arch_sums     <- compute_level_sums(cluster_df, grouping_col, vaccine_variables)
  national_sums <- compute_level_sums(cluster_df, character(0), vaccine_variables)

  # ZONE is already in admin_sums from expand_and_merge_archetype — no merge needed.
  fallback_list <- list(
    list(join_col = "NOMREGION",  data = region_sums),
    list(join_col = grouping_col, data = arch_sums),
    list(join_col = NULL,         data = national_sums)
  )
  write.csv(as.data.frame(admin_sums),
            file.path(out_dir, paste0("DHS_vaccine_admin_", year, ".csv")))
  
  admin_sums <- apply_min_sample_fallback(admin_sums, fallback_list,
                                          vaccine_variables, min_num_total,
                                          include_date = FALSE)
  write.csv(as.data.frame(admin_sums),
            file.path(out_dir, paste0("DHS_vaccine_admin_minN", min_num_total, "_", year, ".csv")))
}


# ---- Vaccine: MICS ----------------------------------------------------------

# Extract vaccination coverage from a MICS survey year.
# Inputs: hbhi_dir, dta_dir, year, admin_shape, ds_pop_df_filename;
#   use_admin_level -- if TRUE also write LGA-level CSV with fallback applied.
# Outputs: MICS_vaccine_cluster/state/admin CSVs.
extract_vaccine_MICS_data <- function(
    hbhi_dir, dta_dir, year, admin_shape, ds_pop_df_filename,
    min_num_total = 30,
    vaccine_variables = c("vacc_dpt1","vacc_dpt2","vacc_dpt3"),
    vaccine_alternate_positive_patterns = c("reported by mother","vaccination marked on card"),
    use_admin_level = FALSE
) {
  out_dir   <- file.path(hbhi_dir, "estimates_from_DHS")
  recode_df <- read.csv(file.path(out_dir, paste0("MICS_", year, "_files_recodes_for_sims.csv")))
  
  loc_idx <- which(recode_df$variable == "locations")
  loc_shp <- raster::shapefile(build_file_path(
    dta_dir, recode_df$folder_dir[loc_idx], recode_df$filename[loc_idx]))
  cluster_df <- data.frame(clusterid = loc_shp$HH1,
                           latitude  = loc_shp$LATITUDE, longitude = loc_shp$LONGITUDE)
  
  for (vv in vaccine_variables) {
    idx <- which(recode_df$variable == vv)
    if (length(idx) == 0 || is.na(recode_df$filename[idx])) next
    path <- build_file_path(dta_dir, recode_df$folder_dir[idx], recode_df$filename[idx])
    df   <- load_survey_file(path)
    df   <- code_positivity(df, recode_df$code[idx],
                            recode_df$pos_pattern[idx], recode_df$neg_pattern[idx],
                            vaccine_alternate_positive_patterns)
    res  <- aggregate_to_clusters_vacc_mics(df, recode_df$cluster_id_code[idx],
                                            recode_df$age_code[idx],
                                            min_age_months = recode_df$min_age_months_to_include[idx])
    cluster_df <- merge_cluster_result(cluster_df, res, recode_df$cluster_id_code[idx], vv)
  }
  
  cluster_df <- assign_clusters_to_admins(cluster_df, admin_shape,
                                          region_col = "State", admin_col = "GEONAMET")$cluster_df
  cluster_df <- cluster_df %>%
    rename(NOMREGION = NOMREGION, NOMDEP = NOMDEP)  # canonical names already set
  
  write.csv(cluster_df, file.path(out_dir, paste0("MICS_vaccine_cluster_outputs_", year, ".csv")))
  
  archetype_info <- normalize_archetype_info(read.csv(ds_pop_df_filename))
  
  # State-level output (always produced)
  cluster_df <- cluster_df[!is.na(cluster_df$NOMREGION), ]
  region_sums <- compute_level_sums(cluster_df, "NOMREGION", vaccine_variables)
  write.csv(as.data.frame(region_sums),
            file.path(out_dir, paste0("MICS_vaccine_state_", year, ".csv")))
  
  if (use_admin_level) {
    admin_sums <- compute_level_sums(cluster_df, c("NOMREGION","NOMDEP"), vaccine_variables)
    admin_sums <- expand_and_merge_archetype(admin_sums, archetype_info)
    
    grouping_col  <- if ("ZONE" %in% colnames(archetype_info)) "ZONE" else "Archetype"
    zone_map      <- .safe_zone_map(archetype_info, grouping_col)

    # Merge grouping column into cluster_df so compute_level_sums can group by it.
    cluster_df <- merge(cluster_df, zone_map, by = "NOMREGION", all.x = TRUE)

    arch_sums     <- compute_level_sums(cluster_df, grouping_col, vaccine_variables)
    national_sums <- compute_level_sums(cluster_df, character(0), vaccine_variables)

    # ZONE is already in admin_sums from expand_and_merge_archetype — no merge needed.
    write.csv(as.data.frame(admin_sums),
              file.path(out_dir, paste0("MICS_vaccine_admin_", year, ".csv")))
    
    fallback_list <- list(
      list(join_col = "NOMREGION",  data = region_sums),
      list(join_col = grouping_col, data = arch_sums),
      list(join_col = NULL,         data = national_sums)
    )
    admin_sums <- apply_min_sample_fallback(admin_sums, fallback_list,
                                            vaccine_variables, min_num_total,
                                            include_date = FALSE)
    write.csv(as.data.frame(admin_sums),
              file.path(out_dir, paste0("MICS_vaccine_admin_minN", min_num_total, "_", year, ".csv")))
  }
}


# =============================================================================
# SECTION 13: ITN CAMPAIGN FRACTION
# Estimate the fraction of ITNs obtained from mass campaigns vs. other sources.
# =============================================================================

# ---- Campaign ITN fraction --------------------------------------------------
# Compute fraction of ITNs from mass campaigns at state level across survey years.
# Inputs: hbhi_dir, dta_dir, years, archetype_info.
# Outputs: itn_source_state_level_all_years.csv; bar chart plot (returned invisibly).
extract_frac_itn_from_campaign_by_state <- function(hbhi_dir, dta_dir,
                                                    years, archetype_info) {
  out_dir     <- file.path(hbhi_dir, "estimates_from_DHS")
  result_all  <- data.frame()
  
  for (year in years) {
    recode_df <- read.csv(file.path(out_dir, paste0("DHS_", year, "_files_recodes_for_sims.csv")))
    if (!"itn_source" %in% recode_df$variable) next
    
    idx  <- which(recode_df$variable == "itn_source")
    path <- build_file_path(dta_dir, recode_df$folder_dir[idx], recode_df$filename[idx])
    df   <- load_survey_file(path)
    
    # For CSPro .dat files, load_survey_file returns plain numeric columns —
    # haven::as_factor() is a no-op on these, so numeric codes like "1","2","3","4"
    # pass through unchanged and the downstream str_detect filter finds no matches.
    # Use parse_sps_value_labels() to decode numeric codes to their label text
    # before any string operations.  For .dta/.sav files the SPS won't exist and
    # this block is skipped.
    src_col   <- recode_df$code[idx]
    state_col <- recode_df$state_code[idx]
    sps_path  <- sub("\\.[Dd][Aa][Tt]$", ".SPS", path)
    if (!file.exists(sps_path))
      sps_path <- sub("\\.[Dd][Aa][Tt]$", ".sps", path)
    if (file.exists(sps_path)) {
      val_labels <- parse_sps_value_labels(sps_path)
      message(sprintf("  [itn_source] SPS labels found for %d variables; %s in labels: %s; %s in labels: %s",
                      length(val_labels),
                      src_col,   if (src_col   %in% names(val_labels)) "YES" else "NO",
                      state_col, if (!is.na(state_col) && state_col %in% names(val_labels)) "YES" else "NO"))
      for (col in c(src_col, state_col)) {
        if (!is.na(col) && nzchar(trimws(col)) && col %in% names(val_labels)) {
          lbl       <- val_labels[[col]]          # named char: c("1"="label A", ...)
          raw_codes <- as.character(df[[col]])
          decoded   <- lbl[raw_codes]             # NA for unrecognised codes
          df[[col]] <- ifelse(is.na(decoded), raw_codes, as.character(decoded))
          message(sprintf("  [itn_source] decoded %s: %d unique values, e.g. %s",
                          col, length(unique(df[[col]])),
                          paste(head(unique(df[[col]]), 3), collapse=", ")))
        }
      }
    } else {
      message("  [itn_source] no SPS file found at: ", sps_path)
    }
    
    itn_source <- df %>%
      dplyr::select(source = all_of(src_col),
                    State  = all_of(recode_df$state_code[idx])) %>%
      # For haven-labelled (.dta/.sav), as_factor() converts numeric codes to label
      # text.  For CSPro (already decoded above), it is a no-op.
      mutate(source = as.character(haven::as_factor(source)),
             State  = as.character(haven::as_factor(State))) %>%
      mutate(State = str_trim(State, "right"),
             State = str_replace(State, "^[A-Za-z]{2} ", ""),
             State = str_remove(State, " (rural|urban)$")) %>%
      filter(!is.na(source),
             str_detect(source, regex("no|antenatal|ANC|immunization|campaign", ignore_case = TRUE))) %T>%
      { message(sprintf("  [itn_source] rows passing source filter: %d  |  unique source values: %s",
                        nrow(.), paste(head(unique(.$source), 5), collapse=", "))) } %>%
      group_by(State) %>%
      summarise(frac_campaign = mean(str_detect(source, regex("campaign", ignore_case = TRUE))),
                n_total       = n(),
                n_campaign    = sum(str_detect(source, regex("campaign", ignore_case = TRUE))),
                .groups = "drop") %>%
      mutate(year = year)
    
    # Rebuild as a plain data.frame with fully stripped types.
    itn_source_clean <- data.frame(
      State         = as.character(itn_source$State),
      frac_campaign = as.numeric(itn_source$frac_campaign),
      n_total       = as.integer(itn_source$n_total),
      n_campaign    = as.integer(itn_source$n_campaign),
      year          = as.integer(itn_source$year),
      stringsAsFactors = FALSE
    )
    
    # standardize_state_names_in_df uses sapply() internally, which returns a list
    # (triggering "sort.list: x must be atomic") when either the origin or target
    # State vector retains haven_labelled or other non-atomic S3 attributes.
    # Replicate its logic inline with vapply() to enforce character(1) output.
    # iconv to UTF-8 first: CSPro SPS value labels may be Latin-1 encoded, and
    # create_reference_name_match uses perl=TRUE regex which requires valid UTF-8.
    safe_utf8 <- function(x) iconv(as.character(x), from = "latin1", to = "UTF-8", sub = "")
    arch_state_col <- safe_utf8(archetype_info[["State"]])
    arch_state_lookup <- unique(data.frame(
      State_norm = vapply(arch_state_col, create_reference_name_match, character(1)),
      State      = arch_state_col,
      stringsAsFactors = FALSE
    ))
    itn_source_clean$State <- safe_utf8(itn_source_clean$State)
    itn_source_clean$State_norm <- vapply(
      itn_source_clean$State, create_reference_name_match, character(1))
    
    unmatched <- itn_source_clean$State_norm[
      !itn_source_clean$State_norm %in% arch_state_lookup$State_norm]
    if (length(unmatched) > 0)
      warning("extract_frac_itn_from_campaign_by_state: unmatched state names: ",
              paste(unique(unmatched), collapse=", "))
    
    itn_source <- itn_source_clean %>%
      dplyr::select(-State) %>%
      merge(arch_state_lookup, by = "State_norm", all.x = TRUE) %>%
      dplyr::select(-State_norm)
    
    result_all <- dplyr::bind_rows(result_all, itn_source)
  }
  
  ggplot2::ggplot(result_all, ggplot2::aes(x=State, y=frac_campaign,
                                           fill=as.factor(year))) +
    ggplot2::geom_col(position="dodge") +
    ggplot2::labs(x="State", y="Fraction of ITNs from mass campaign",
                  title="Fraction of ITNs from mass campaigns in DHS/MIS") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  
  write.csv(result_all,
            file.path(out_dir, "itn_source_state_level_all_years.csv"), row.names=FALSE)
}


# =============================================================================
# SECTION 14: ITN USE BY AGE
# Compute ITN use by age group relative to the all-age (or U5) rate.
# =============================================================================

# Compute per-age-group ITN use rates relative to all-age and U5 rates.
# Inputs: hbhi_dir, years -- reads DHS_cluster_outputs_{year}.csv.
# Outputs: DHS_weighted/unweighted_relative_ITN_use_by_age.csv,
#   DHS_weighted_ITN_use_by_age_relative_to_u5.csv,
#   plots/DHS_relative_ITN_use_by_age.pdf.

get_relative_itn_use_by_age <- function(hbhi_dir, years) {
  out_dir   <- file.path(hbhi_dir, "estimates_from_DHS")
  plots_dir <- file.path(out_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  age_grps <- c("u5", "5_10", "10_15", "15_20", "o20")
  
  save_weighted   <- matrix(NA, nrow = 0, ncol = length(age_grps))
  save_unweighted <- matrix(NA, nrow = 0, ncol = length(age_grps))
  years_done <- c()
  
  # Collect per-year data for plotting after all years are processed
  plot_data <- list()
  
  for (year in years) {
    f <- file.path(out_dir, paste0("DHS_cluster_outputs_", year, ".csv"))
    if (!file.exists(f)) {
      warning("get_relative_itn_use_by_age: cluster file not found for year ",
              year, " — skipping.")
      next
    }
    cl <- read.csv(f) %>% dplyr::select(-any_of("X"))
    
    # Per-cluster relative rates (age rate / all-age rate)
    for (grp in age_grps)
      cl[[paste0("rel_", grp)]] <- cl[[paste0("itn_", grp, "_rate")]] /
      cl$itn_all_rate
    
    # Drop clusters where any relative rate is Inf (itn_all_rate == 0)
    cl2 <- cl[is.finite(cl$rel_5_10) & is.finite(cl$rel_10_15) &
                is.finite(cl$rel_15_20) & is.finite(cl$rel_o20), ]
    
    rel_cols <- paste0("rel_", age_grps)
    wt       <- cl2$itn_weights
    
    # DHS-weighted mean of cluster-level relative rates
    weighted_means <- sapply(rel_cols, function(rc)
      weighted.mean(cl2[[rc]], wt, na.rm = TRUE))
    
    # Unweighted cluster means (simple mean of per-cluster ratios)
    unweighted_cluster_means <- sapply(rel_cols, function(rc)
      mean(cl2[[rc]], na.rm = TRUE))
    
    # Equal weight to all surveyed individuals (pooled counts across clusters)
    all_rate <- sum(cl$itn_all_num_true,  na.rm = TRUE) /
      sum(cl$itn_all_num_total, na.rm = TRUE)
    unweighted_means <- sapply(age_grps, function(grp) {
      ag_rate <- sum(cl[[paste0("itn_", grp, "_num_true")]],  na.rm = TRUE) /
        sum(cl[[paste0("itn_", grp, "_num_total")]], na.rm = TRUE)
      ag_rate / all_rate
    })
    
    save_weighted   <- rbind(save_weighted,   weighted_means)
    save_unweighted <- rbind(save_unweighted, unweighted_means)
    years_done      <- c(years_done, year)
    
    plot_data[[as.character(year)]] <- list(
      weighted            = weighted_means,
      unweighted_clusters = unweighted_cluster_means,
      unweighted_pooled   = unweighted_means
    )
  }
  
  # ---- Compute downstream relative-to-U5 values (needed for plot below) -------
  last_year    <- years_done[length(years_done)]
  recent_years <- years_done[years_done >= last_year - 5]
  message(sprintf("  DHS_weighted_ITN_use_by_age_relative_to_u5: averaging over years: %s",
                  paste(recent_years, collapse=", ")))
  
  yearly_rel_u5 <- lapply(recent_years, function(yr) {
    cl <- read.csv(file.path(out_dir,
                             paste0("DHS_cluster_outputs_", yr, ".csv"))) %>%
      dplyr::select(-any_of("X"))
    
    for (grp in age_grps)
      cl[[paste0("rel_", grp)]] <- cl[[paste0("itn_", grp, "_rate")]] /
        cl$itn_u5_rate
    
    cl2 <- cl[is.finite(cl$rel_5_10) & is.finite(cl$rel_10_15) &
                is.finite(cl$rel_15_20) & is.finite(cl$rel_o20), ]
    wt <- cl2$itn_weights
    sapply(paste0("rel_", age_grps), function(rc)
      weighted.mean(cl2[[rc]], wt, na.rm = TRUE))
  })
  
  # Average the per-year weighted means across included years
  rel_u5_means <- rowMeans(do.call(cbind, yearly_rel_u5), na.rm = TRUE)
  rel_u5_means <- pmin(rel_u5_means, 1)   # cap at 1 (can exceed due to sampling noise)
  
  # ---- Plot: one panel per year + downstream panel, all in a single PDF ------
  n_years <- length(years_done)
  if (n_years > 0) {
    # Total data panels = per-year panels + 1 downstream panel
    n_panels <- n_years + 1L
    ncols    <- min(3L, n_panels)
    nrows    <- ceiling(n_panels / ncols)
    pw       <- ncols * 3.5
    ph       <- nrows * 3.2 + 1.2   # extra inch at bottom for legend
    age_labels <- c("U5", "5-10", "10-15", "15-20", ">20")
    bar_cols   <- c(rgb(0.3, 0.3, 1),   # unweighted cluster means
                    rgb(0.55, 0.55, 1),  # DHS-weighted cluster means
                    rgb(0.8, 0.8, 1))    # equal-weight pooled
    
    pdf(file.path(plots_dir, "DHS_relative_ITN_use_by_age.pdf"),
        width = pw, height = ph, useDingbats = FALSE)
    layout(
      mat    = rbind(matrix(seq_len(nrows * ncols), nrow = nrows, byrow = TRUE),
                     rep(nrows * ncols + 1L, ncols)),
      heights = c(rep(1, nrows), 0.35)
    )
    par(mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0), cex = 0.85)
    
    # Per-year panels (relative to all-age)
    for (yr in years_done) {
      pd <- plot_data[[as.character(yr)]]
      ratios <- matrix(c(pd$unweighted_clusters, pd$weighted, pd$unweighted_pooled),
                       byrow = TRUE, nrow = 3)
      ymax <- max(ratios, na.rm = TRUE) * 1.15
      barplot(ratios, beside = TRUE,
              col       = bar_cols,
              names.arg = age_labels,
              ylim      = c(0, ymax),
              xlab      = "Age group",
              ylab      = "Rate relative to all-age",
              main      = as.character(yr))
      abline(h = 1, lty = 2, col = "grey50")
    }
    
    # Downstream panel: weighted mean relative to U5, averaged over recent years
    ds_ratios <- matrix(rel_u5_means, nrow = 1)
    ymax_ds   <- max(ds_ratios, na.rm = TRUE) * 1.15
    ds_title  <- if (length(recent_years) == 1) {
      paste0("Used downstream\n(", recent_years, ", rel. to U5)")
    } else {
      paste0("Used downstream\n(avg ", min(recent_years), "\u2013", max(recent_years),
             ", rel. to U5)")
    }
    barplot(ds_ratios, beside = TRUE,
            col       = rgb(0.2, 0.6, 0.2),
            names.arg = age_labels,
            ylim      = c(0, max(ymax_ds, 1.15)),
            xlab      = "Age group",
            ylab      = "Rate relative to U5",
            main      = ds_title)
    abline(h = 1, lty = 2, col = "grey50")
    
    # Fill any empty cells in the grid with blank panels
    n_empty <- nrows * ncols - n_panels
    for (i in seq_len(n_empty)) plot.new()
    
    # Legend panel spanning full bottom row
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
           legend = c("Unweighted mean of clusters",
                      "DHS-weighted mean of clusters",
                      "Equal weight per individual (pooled)",
                      "Downstream value (wt. mean rel. to U5)"),
           fill   = c(bar_cols, rgb(0.2, 0.6, 0.2)),
           border = NA,
           horiz  = TRUE,
           bty    = "n",
           cex    = 0.8)
    
    par(mfrow = c(1, 1))
    dev.off()
  }
  
  # ---- Write multi-year relative-to-all-age CSVs --------------------------------
  col_names <- paste0("relative_use_", age_grps)
  
  colnames(save_weighted)   <- col_names
  wt_df  <- as.data.frame(save_weighted);  wt_df$year  <- years_done
  wt_df  <- wt_df[, c("year", col_names)]
  
  colnames(save_unweighted) <- col_names
  unw_df <- as.data.frame(save_unweighted); unw_df$year <- years_done
  unw_df <- unw_df[, c("year", col_names)]
  
  write.csv(wt_df,  file.path(out_dir, "DHS_weighted_relative_ITN_use_by_age.csv"),
            row.names = FALSE)
  write.csv(unw_df, file.path(out_dir, "DHS_unweighted_relative_ITN_use_by_age.csv"),
            row.names = FALSE)
  
  # ---- Write single-row relative-to-U5 CSV -----------------------------------
  # (rel_u5_means already computed above for the plot)
  
  rel_u5_df <- as.data.frame(matrix(rel_u5_means, nrow = 1))
  colnames(rel_u5_df) <- paste0(age_grps, "_use_relative_to_u5")
  write.csv(rel_u5_df,
            file.path(out_dir, "DHS_weighted_ITN_use_by_age_relative_to_u5.csv"),
            row.names = FALSE)
  
  message("get_relative_itn_use_by_age: outputs written to ", out_dir)
}


# =============================================================================
# SECTION 15: MONTHLY PfPR REFERENCE
# Build monthly PfPR reference curves by admin unit, with region-level fallback.
# =============================================================================

# Build monthly PfPR reference curves by admin unit from multiple DHS/MIS years.
# Inputs: hbhi_dir, dta_dir, admin_shape, ds_pop_df_filename, pfpr_dhs_ref_years.
# Outputs: DHS_monthly_{measure}_adminLevelDataOnly.csv and
#   DHS_admin_monthly_{measure}.csv (with region-level fallback for small admins).
create_DHS_reference_monthly_pfpr <- function(
    hbhi_dir, dta_dir, admin_shape, ds_pop_df_filename,
    pfpr_dhs_ref_years = c(2012, 2016),
    min_num_total = 30,
    pfpr_measure  = "mic"
) {
  out_dir        <- file.path(hbhi_dir, "estimates_from_DHS")
  pfpr_name      <- switch(pfpr_measure, mic="microscopy", rdt="RDT", {
    warning("Unrecognised PfPR measure; defaulting to microscopy.")
    pfpr_measure <<- "mic"; "microscopy"
  })
  
  all_years_df <- data.frame()
  
  for (year in pfpr_dhs_ref_years) {
    recode_df <- read.csv(file.path(out_dir, paste0("DHS_", year, "_files_recodes_for_sims.csv")))
    
    loc_idx  <- which(recode_df$variable == "locations")
    loc_path <- build_file_path(dta_dir, recode_df$folder_dir[loc_idx],
                                recode_df$filename[loc_idx])
    loc_ext  <- tolower(tools::file_ext(loc_path))
    
    if (loc_ext == "shp") {
      # Classic DHS format: GPS shapefile with DHSCLUST / LATNUM / LONGNUM columns
      loc_shp    <- raster::shapefile(loc_path)
      cluster_df <- data.frame(clusterid = loc_shp$DHSCLUST,
                               latitude  = loc_shp$LATNUM,
                               longitude = loc_shp$LONGNUM)
    } else if (loc_ext == "dat") {
      # CSPro format: read GPS from .dat file using recode CSV column names
      cid_col <- as.character(recode_df$cluster_id_code[loc_idx])
      lat_col <- as.character(recode_df$latitude_code[loc_idx])
      lon_col <- as.character(recode_df$longitude_code[loc_idx])
      loc_df  <- as.data.frame(read_cspr(loc_path))
      # Apply implicit-decimal scaling (skip lat/lon columns — already handled by read_cspr SPS)
      cluster_df <- loc_df %>%
        dplyr::select(clusterid = all_of(cid_col),
                      latitude  = all_of(lat_col),
                      longitude = all_of(lon_col)) %>%
        mutate(across(c(latitude, longitude), as.numeric)) %>%
        group_by(clusterid) %>%
        summarise(latitude  = mean(latitude,  na.rm = TRUE),
                  longitude = mean(longitude, na.rm = TRUE),
                  .groups = "drop") %>%
        as.data.frame()
    } else {
      stop("create_DHS_reference_monthly_pfpr: unsupported locations file type '.",
           loc_ext, "' for year ", year)
    }
    
    idx <- which(recode_df$variable == pfpr_measure)
    if (length(idx) == 0 || is.na(recode_df$filename[idx])) {
      warning("No filename for ", pfpr_name, " in year ", year); next
    }
    path <- build_file_path(dta_dir, recode_df$folder_dir[idx], recode_df$filename[idx])
    df   <- load_survey_file(path)
    df   <- code_positivity(df, recode_df$code[idx],
                            recode_df$pos_pattern[idx], recode_df$neg_pattern[idx])
    
    month_col <- recode_df$month_code[idx]
    year_col  <- recode_df$year_code[idx]
    valid_col <- function(x) !is.na(x) && nzchar(trimws(x))
    
    if (!valid_col(month_col))
      stop("create_DHS_reference_monthly_pfpr: month_code is missing for variable '",
           pfpr_measure, "' in year ", year, ".\n",
           "  ACTION REQUIRED: add the interview month column name to the 'month_code' field\n",
           "  in DHS_", year, "_files_recodes_for_sims.csv for the '", pfpr_measure, "' row.")
    if (!valid_col(year_col))
      stop("create_DHS_reference_monthly_pfpr: year_code is missing for variable '",
           pfpr_measure, "' in year ", year, ".\n",
           "  ACTION REQUIRED: add the interview year column name to the 'year_code' field\n",
           "  in DHS_", year, "_files_recodes_for_sims.csv for the '", pfpr_measure, "' row.")
    # Resolve CSPro occurrence-indexed columns (e.g. QHVMONTH → QHVMONTH$1).
    resolved_month <- resolve_occurrence_col(df, month_col)
    resolved_year  <- resolve_occurrence_col(df, year_col)
    if (is.null(resolved_month))
      stop("create_DHS_reference_monthly_pfpr: month_code column '", month_col,
           "' (or '", paste0(month_col, "$1"), "') not found in the survey file",
           " for variable '", pfpr_measure, "' year ", year, ".\n",
           "  ACTION REQUIRED: check the column name in DHS_", year, "_files_recodes_for_sims.csv.")
    if (is.null(resolved_year))
      stop("create_DHS_reference_monthly_pfpr: year_code column '", year_col,
           "' (or '", paste0(year_col, "$1"), "') not found in the survey file",
           " for variable '", pfpr_measure, "' year ", year, ".\n",
           "  ACTION REQUIRED: check the column name in DHS_", year, "_files_recodes_for_sims.csv.")
    month_col <- resolved_month
    year_col  <- resolved_year
    
    month_vals <- df[[month_col]]
    # Coerce haven_labelled to integer first (strips label map, gives raw numeric code).
    # Then handle character month-name labels (e.g. "January" → 1).
    if (haven::is.labelled(month_vals)) month_vals <- as.integer(month_vals)
    if (is.character(month_vals) && any(month_vals %in% month.name, na.rm=TRUE))
      month_vals <- match(month_vals, month.name)
    df$month <- as.integer(month_vals)
    df$year  <- as.integer(df[[year_col]])
    
    res <- aggregate_to_clusters_monthly(df, recode_df$cluster_id_code[idx],
                                         "month", "year")
    cluster_df <- merge(cluster_df, res,
                        by.x = "clusterid", by.y = recode_df$cluster_id_code[idx],
                        all = TRUE)
    cluster_df <- assign_clusters_to_admins(cluster_df, admin_shape)$cluster_df
    cluster_df <- cluster_df %>% rename(admin_name = NOMDEP)
    
    all_years_df <- if (nrow(all_years_df) == 0) as.data.frame(cluster_df)
    else rbind(all_years_df, as.data.frame(cluster_df))
  }
  
  all_years_df <- all_years_df[!is.na(all_years_df$admin_name), ]
  
  # Aggregate to admin × month × year
  admin_monthly <- all_years_df %>%
    group_by(admin_name, NOMREGION, month, year) %>%
    summarise(num_pos = sum(num_pos, na.rm=TRUE),
              num_tested = sum(num_tested, na.rm=TRUE), .groups="drop")
  
  archetype_info <- normalize_archetype_info(read.csv(ds_pop_df_filename)) %>%
    dplyr::select(admin_name = NOMDEP, NOMREGION, Archetype)
  
  admin_monthly <- standardize_admin_names_in_df(
    target_names_df = archetype_info, origin_names_df = admin_monthly,
    target_names_col = "admin_name", origin_names_col = "admin_name")
  
  admin_monthly <- merge(admin_monthly, archetype_info, all = TRUE) %>%
    mutate(data_spatial_level = "admin",
           num_pos    = replace(num_pos,    is.na(num_pos),    0),
           num_tested = replace(num_tested, is.na(num_tested), 0))
  
  write.csv(admin_monthly,
            file.path(out_dir, paste0("DHS_monthly_", pfpr_name, "_adminLevelDataOnly.csv")),
            row.names = FALSE)
  
  # Region-level fallback for admins below min_num_total
  region_monthly <- admin_monthly %>%
    group_by(NOMREGION, Archetype, month, year) %>%
    summarise(num_pos = sum(num_pos), num_tested = sum(num_tested), .groups="drop")
  
  admin_totals <- admin_monthly %>%
    group_by(admin_name, NOMREGION) %>%
    summarise(total_tested = sum(num_tested), .groups="drop")
  
  if (any(region_monthly %>% group_by(NOMREGION) %>%
          summarise(t=sum(num_tested)) %>% pull(t) < min_num_total)) {
    warning("Some regions have fewer than min_num_total U5s tested with ",
            pfpr_name, ". Further aggregation may be needed.")
  }
  
  adjusted <- admin_monthly
  for (i in seq_len(nrow(admin_totals))) {
    adm <- admin_totals$admin_name[i]; reg <- admin_totals$NOMREGION[i]
    if (admin_totals$total_tested[i] < min_num_total) {
      region_rows <- region_monthly[region_monthly$NOMREGION == reg &
                                      region_monthly$num_tested > 0, ]
      region_rows$admin_name        <- adm
      region_rows$data_spatial_level <- "region"
      adjusted <- adjusted[adjusted$admin_name != adm, ]
      adjusted <- merge(adjusted, region_rows, all = TRUE)
    }
  }
  adjusted <- adjusted[adjusted$num_tested > 0 & !is.na(adjusted$year), ]
  write.csv(adjusted,
            file.path(out_dir, paste0("DHS_admin_monthly_", pfpr_name, ".csv")),
            row.names = FALSE)
}
