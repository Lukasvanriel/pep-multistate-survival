# ==============================================================================
# covariate_spec.R   (B1: covariate and formula specification)
#
# Single place that says WHICH covariates exist, how they are scaled, and which
# terms each kept transition gets. Plays the role of covariate_config.R in the
# simulation folder (COV_NAMES / rhs_formula_from_covs / beta_mat), but the term
# set differs per transition (fire / pest / harvest terms pass the >= 15 rule
# only for some transitions).
#
# COV_CONT is FINAL: gdd_3, cmi_gs (CMI over the NRCan growing season = window of gdd_3, from 03c step 4b since 2026-10-24; replaced May-Oct and annual cmi_sum), soil_pc1.
#
# Model terms (per kept transition, see term_mask()):
#   main effects : z_gdd_3, z_cmi_gs, z_soil_pc1, fire_f, pest_f, harv_light + harv_heavy (harv_term "3level" = light and heavy
#                  each vs none; there is no "any" level since A2, 2026-11-04).
#                  z_gdd_3 is dropped where a transition has no undisturbed event (A11, n_undist == 0 and >= 1 disturbance term).
#   GDD x disturbance (":" products; GDD only, not CMI). The GDD main effect is present except in A11 transitions (z_gdd_3 dropped):
#                  there the product is the GDD slope WITHIN the disturbance level (unharvested / unburned slope assumed 0, untestable):
#                  gdd_x_fire, gdd_x_pest, gdd_x_harv_light, gdd_x_harv_heavy = z_gdd_3 * flag.
#                  Wherever the disturbance main term exists; for a harvest LEVEL with 0 events the main term only (A6).
#   TSD (one per disturbance type; harvest per LEVEL, A4): tsd_fire, tsd_pest, tsd_harv_light, tsd_harv_heavy.
#                  Model column = flag * f(TSD in years) with f from TSD_FORM (config.R): "log1p" = log(x + 1) (main),
#                  "linear" = x / 10 (decades, S24). Raw TSD in years (tsd_*_y, 0 where the flag is 0) comes from 04; not capped (A8).
#                  Wherever the disturbance main term exists the TSD term exists too (harvest level: only with >= 1 event, A6).
#   Weibull variant 0 (WEIBULL_VARIANT, config.R): S(t) = exp(-e^eta t^alpha); beta is a log hazard ratio. (The fit functions
#   follow in WP3; until then 05a_fit_funs.R still uses its earlier variant.)
#
# Disturbance flags are taken at the END visit of the interval (decision default;
# the *_start versions are still in the interval table).
# ==============================================================================

COV_CONT <- c("gdd_3", "cmi_gs", "soil_pc1")   # final (cmi_gs replaced cmi_sum 2026-10-22; same position)

# Order matters only for P1: P1 uses the FIRST continuous covariate's value of the
# plot's first interval (the simulation's P1 supports one covariate row).

# ------------------------------------------------------------------------------
# make_model_covariates()
# Adds the model columns to an interval table (or a row table that carries the
# same raw columns): scaled continuous (z-scores, centre/scale taken ONCE from the
# interval table so rows and intervals agree), fire, pest, harvest dummies.
# Returns list(data = ..., scaling = data.table(var, mean, sd)).
# ------------------------------------------------------------------------------
make_scaling <- function(iv, cont = COV_CONT) {
  data.table::data.table(var = cont,
                         mean = vapply(cont, function(v) mean(iv[[v]], na.rm = TRUE), numeric(1)),
                         sd   = vapply(cont, function(v) sd(iv[[v]],   na.rm = TRUE), numeric(1)))
}

add_model_covariates <- function(d, scaling, tsd_form = TSD_FORM) {   # TSD_FORM from R/config.R (source it first)
  stopifnot("tsd_form must be log1p or linear" = tsd_form %in% c("log1p", "linear"))
  f_tsd <- if (tsd_form == "log1p") function(x) log(x + 1) else function(x) x / 10
  d <- data.table::copy(d)
  for (i in seq_len(nrow(scaling))) {
    v <- scaling$var[i]
    data.table::set(d, j = paste0("z_", v), value = (d[[v]] - scaling$mean[i]) / scaling$sd[i])
  }
  d[, `:=`(fire_f = as.integer(fire == 1L), pest_f = as.integer(pest == 1L),
           harv_light = as.integer(harv == 1L), harv_heavy = as.integer(harv == 2L))]
  # GDD x disturbance (":" products; main effects are separate terms)
  stopifnot("z_gdd_3" %in% names(d))
  d[, `:=`(gdd_x_fire = z_gdd_3 * fire_f, gdd_x_pest = z_gdd_3 * pest_f,
           gdd_x_harv_light = z_gdd_3 * harv_light, gdd_x_harv_heavy = z_gdd_3 * harv_heavy)]
  # TSD in years (structural zero where the flag is 0) comes from 04 (A5); check it, then build the model columns (A4, A5).
  for (v in c("tsd_fire_y", "tsd_pest_y", "tsd_harv_y")) stopifnot(v %in% names(d), !anyNA(d[[v]]), all(d[[v]] >= 0))
  stopifnot(all(d$tsd_fire_y[d$fire_f == 0L] == 0), all(d$tsd_pest_y[d$pest_f == 0L] == 0), all(d$tsd_harv_y[d$harv == 0L] == 0))
  d[, `:=`(tsd_fire = fire_f * f_tsd(tsd_fire_y), tsd_pest = pest_f * f_tsd(tsd_pest_y),
           tsd_harv_light = harv_light * f_tsd(tsd_harv_y), tsd_harv_heavy = harv_heavy * f_tsd(tsd_harv_y))]
  stopifnot(!anyNA(d[, .(tsd_fire, tsd_pest, tsd_harv_light, tsd_harv_heavy)]))
  d
}

# All possible model terms (columns of the beta matrix).
all_terms <- function(cont = COV_CONT)
  c(paste0("z_", cont), "fire_f", "pest_f", "harv_light", "harv_heavy",
    "gdd_x_fire", "gdd_x_pest", "gdd_x_harv_light", "gdd_x_harv_heavy",
    "tsd_fire", "tsd_pest", "tsd_harv_light", "tsd_harv_heavy")

# ------------------------------------------------------------------------------
# term_mask()
# Logical matrix, one row per KEPT transition, one column per possible term:
# TRUE = the term is in that transition's formula. Zeros in the beta matrix for
# FALSE entries (same role as zeros in the simulation's beta_mat).
# ------------------------------------------------------------------------------
term_mask <- function(tr, cont = COV_CONT) {
  stopifnot("transition_set lacks n_harv1 / n_harv2 / n_undist (rerun 04 of the final pipeline)" = all(c("n_harv1", "n_harv2", "n_undist") %in% names(tr)))
  k <- tr[kept == TRUE]
  terms <- all_terms(cont)
  m <- matrix(FALSE, nrow(k), length(terms), dimnames = list(paste(k$from, k$to, sep = "->"), terms))
  m[, paste0("z_", cont)] <- TRUE                       # continuous terms: every kept transition
  m[, "fire_f"]     <- k$fire_term
  m[, "pest_f"]     <- k$pest_term
  m[, "harv_light"] <- k$harv_term == "3level"
  m[, "harv_heavy"] <- k$harv_term == "3level"
  # A11: no GDD main term where a transition has no undisturbed event and at least one disturbance term
  has_dist <- k$fire_term | k$pest_term | k$harv_term == "3level"
  m[, "z_gdd_3"] <- !(k$n_undist == 0L & has_dist)
  # GDD x disturbance: wherever the disturbance main term exists (no extra threshold)
  m[, "gdd_x_fire"]       <- k$fire_term
  m[, "gdd_x_pest"]       <- k$pest_term
  # A6: per harvest level the GDD x and TSD terms need >= 1 event of that level (a 0-event level keeps its main term only)
  m[, "gdd_x_harv_light"] <- k$harv_term == "3level" & k$n_harv1 > 0L
  m[, "gdd_x_harv_heavy"] <- k$harv_term == "3level" & k$n_harv2 > 0L
  # TSD: one per disturbance type (harvest per level), wherever that type's main term exists (A4)
  m[, "tsd_fire"] <- k$fire_term
  m[, "tsd_pest"] <- k$pest_term
  m[, "tsd_harv_light"] <- k$harv_term == "3level" & k$n_harv1 > 0L
  m[, "tsd_harv_heavy"] <- k$harv_term == "3level" & k$n_harv2 > 0L
  m
}

# RHS string per kept transition (named by "from->to"), for fit_inla_models().
rhs_formulas <- function(mask) {
  out <- apply(mask, 1, function(r) if (any(r)) paste(colnames(mask)[r], collapse = " + ") else "1")
  out
}
