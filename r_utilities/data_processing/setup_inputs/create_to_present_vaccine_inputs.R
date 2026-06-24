# create_to_present_vaccine_inputs.R

library(dplyr)
library(tidyr)


#=======================================#
# setup vaccine inputs
#=======================================#
# Builds the to-present vaccine input CSVs. LGAs included in the output are
# restricted to the LGAs that appeared in the original pilot file (i.e. those
# observed to receive vaccine in the early rollout). Four coverage assumptions
# are produced in a single call, one CSV per variant:
#
#   1) EPI-capped  — per-LGA DHS EPI coverage (min of DPT3 and measles),
#                    floored at 0.01 and capped at max_r21_epi_cov.
#   2) Kebbi pre-scale-up uniform  — single coverage value applied to all LGAs,
#                                     drawn from the state-level summary CSV at
#                                     unit="Kebbi", period="Pre-Aug 2025",
#                                     Method C (stable r1 x cumulative p31).
#   3) Kebbi full-period uniform   — same as (2) but period="Full period".
#   4) Per-state Method C          — each LGA receives its own state's Method C
#                                     "Full period" coverage (Bayelsa LGAs get the
#                                     Bayelsa value, Kebbi LGAs get the Kebbi value).
#
# Inputs:
#   original_lga_coverages: data.frame from the original LGA-level pilot
#       coverage projections (vacc_dose3_coverages_projected.csv). Only used
#       for the list of LGAs to include in the output (column LGA, with State).
#   state_summary: data.frame from the state-level summary
#       (vacc_p3_summary.csv) with columns unit, period, method, coverage.
#   vacc_dhs_filepath: path to the DHS EPI vaccine coverage CSV.
#
# Output CSVs are written to hbhi_dir/simulation_inputs/interventions_2010_toPresent/:
#   vacc_2010_toPresent_epi{cap_pct}.csv   (cap_pct = round(max_r21_epi_cov*100))
#   vacc_2010_toPresent_kebbiPre.csv
#   vacc_2010_toPresent_kebbiFull.csv
#   vacc_2010_toPresent_byStateC.csv
#
# Each output CSV has two rows per included LGA (primary + booster) with columns:
#   admin_name          -- LGA name (from the original pilot file's inclusion list)
#   coverage            -- fraction receiving the dose. Primary row: variant-specific
#                          (DHS-derived per LGA for variant 1; flat Kebbi value for
#                          variants 2 and 3). Booster row: booster_relative_coverage,
#                          interpreted as conditional on receiving primary.
#   deploy_type         -- 'EPI'
#   RTSS_day            -- simulation day rollout begins (shifted earlier by
#                          primary_series_day to account for birth-triggered timing)
#   vaccine             -- 'primary' or 'booster'
#   rtss_touchpoints    -- child age (days) at this dose: primary_series_day for
#                          primary, booster_day for booster
#   distribution_name   -- 'CONSTANT_DISTRIBUTION'
#   distribution_std    -- 1 (placeholder for CONSTANT_DISTRIBUTION)
#   initial_effect      -- vaccine efficacy at deployment: primary_initial_effect for
#                          primary, booster_initial_effect for booster
#   decay_time_constant -- exponential efficacy decay timescale (days)
#   decay_class         -- 'WaningEffectExponential'
#   duration            -- -1 (persists through end of simulation)

create_vaccine_input_files = function(
    hbhi_dir,
    original_lga_coverages,
    state_summary,
    vacc_dhs_filepath,
    vacc_rollout_day,
    primary_initial_effect,
    booster_initial_effect,
    decay_time_constant,
    primary_series_day,
    booster_day,
    max_r21_epi_cov,
    booster_relative_coverage = 0.8
) {

  # ---- Step 1: LGA inclusion list from the original pilot file ----
  pilot_lgas = original_lga_coverages %>%
    dplyr::distinct(State, LGA) %>%
    dplyr::rename(admin_name = LGA)

  # ---- Step 2: Pull Method C coverage for Kebbi from the state-level summary ----
  # Method label in the summary file is the multi-line factor label, e.g.
  # "C: Stable r1 × cum p31\n(r1_stable × 12 × p31_cum)". Match on the "C:" prefix.
  method_c = state_summary %>%
    dplyr::filter(grepl("^C:", method)) %>%
    dplyr::select(unit, period, coverage)

  get_method_c = function(unit_label, period_label) {
    val = method_c %>%
      dplyr::filter(unit == unit_label, period == period_label) %>%
      dplyr::pull(coverage)
    if (length(val) != 1 || is.na(val))
      stop("Could not find a unique Method C coverage for unit = '", unit_label,
           "', period = '", period_label, "' in state_summary.")
    val
  }
  kebbi_pre    = get_method_c("Kebbi",   "Pre-Aug 2025")
  kebbi_full   = get_method_c("Kebbi",   "Full period")
  bayelsa_full = get_method_c("Bayelsa", "Full period")
  message(sprintf("  [vaccine] Method C coverage: Kebbi pre-Aug 2025 = %.3f, Kebbi full = %.3f, Bayelsa full = %.3f",
                  kebbi_pre, kebbi_full, bayelsa_full))

  # ---- Step 3: Variant 1 coverage — per-LGA DHS EPI, capped at max_r21_epi_cov ----
  var_names = c('vacc_dpt3', 'vacc_measles')
  vacc_dhs = read.csv(vacc_dhs_filepath)
  vacc_dhs$admin_name = vacc_dhs$NOMDEP
  vacc_dhs$coverage = apply(
    vacc_dhs[, which(colnames(vacc_dhs) %in% paste0(var_names, '_rate'))], 1, min)
  vacc_dhs = vacc_dhs %>%
    dplyr::mutate(
      coverage = ifelse(coverage < 0.01, 0.01, coverage),
      coverage = ifelse(coverage > max_r21_epi_cov, max_r21_epi_cov, coverage)
    ) %>%
    dplyr::select(admin_name, coverage)

  variant_epi = pilot_lgas %>%
    dplyr::left_join(vacc_dhs, by = "admin_name") %>%
    # any pilot LGA without a DHS match: floor at 0.01 (same behaviour as the
    # original script's coverage floor)
    dplyr::mutate(coverage = ifelse(is.na(coverage), 0.01, coverage)) %>%
    dplyr::select(admin_name, coverage)

  # ---- Step 4: Variants 2 and 3 — uniform Kebbi-derived coverage ----
  variant_kebbi_pre  = pilot_lgas %>%
    dplyr::mutate(coverage = kebbi_pre)  %>% dplyr::select(admin_name, coverage)
  variant_kebbi_full = pilot_lgas %>%
    dplyr::mutate(coverage = kebbi_full) %>% dplyr::select(admin_name, coverage)

  # ---- Step 4b: Variant 4 — per-state Method C "Full period" ----
  # LGAs in Bayelsa receive the Bayelsa Method C estimate; LGAs in Kebbi receive
  # the Kebbi Method C estimate. Stops if any pilot LGA's state is unmapped.
  state_cov_lookup = tibble::tibble(
    State    = c("Bayelsa",      "Kebbi"),
    coverage = c(bayelsa_full,   kebbi_full)
  )
  variant_by_state = pilot_lgas %>%
    dplyr::left_join(state_cov_lookup, by = "State")
  if (any(is.na(variant_by_state$coverage))) {
    missing = variant_by_state %>%
      dplyr::filter(is.na(coverage)) %>%
      dplyr::pull(State) %>% unique()
    stop("Variant 4 (per-state Method C): no Method C estimate available for state(s): ",
         paste(missing, collapse = ", "))
  }
  variant_by_state = variant_by_state %>% dplyr::select(admin_name, coverage)

  # ---- Step 5: Build primary + booster rows and save each variant ----
  out_dir = file.path(hbhi_dir, 'simulation_inputs', 'interventions_2010_toPresent')
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  build_variant_df = function(cov_df) {
    primary = data.frame(
      admin_name          = cov_df$admin_name,
      coverage            = cov_df$coverage,
      deploy_type         = 'EPI',
      RTSS_day            = vacc_rollout_day,
      vaccine             = 'primary',
      rtss_touchpoints    = primary_series_day,
      distribution_name   = 'CONSTANT_DISTRIBUTION',
      distribution_std    = 1,
      initial_effect      = primary_initial_effect,
      decay_time_constant = decay_time_constant,
      decay_class         = 'WaningEffectExponential',
      duration            = -1
    )
    boost = primary %>%
      dplyr::mutate(
        coverage         = booster_relative_coverage,
        vaccine          = 'booster',
        rtss_touchpoints = booster_day,
        initial_effect   = booster_initial_effect
      )
    rbind(primary, boost)
  }

  write_variant = function(cov_df, suffix) {
    df = build_variant_df(cov_df)
    out_file = file.path(out_dir, paste0('vacc_2010_toPresent_', suffix, '.csv'))
    write.csv(df, out_file, row.names = FALSE)
    message('  [vaccine] Written ', nrow(cov_df), ' LGAs to ', out_file)
  }

  epi_suffix = paste0('epi', sprintf('%g', round(max_r21_epi_cov * 100)))
  write_variant(variant_epi,        epi_suffix)
  write_variant(variant_kebbi_pre,  'kebbiPre')
  write_variant(variant_kebbi_full, 'kebbiFull')
  write_variant(variant_by_state,   'byStateC')

  invisible(list(
    epi        = variant_epi,
    kebbi_pre  = variant_kebbi_pre,
    kebbi_full = variant_kebbi_full,
    by_state   = variant_by_state
  ))
}
