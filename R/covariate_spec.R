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
#   main effects : z_gdd_3, z_cmi_gs, z_soil_pc1, fire_f, pest_f, harv_light + harv_heavy ("3level")
#                  or harv_any ("any")
#   GDD x disturbance (":" products, main effects both present; GDD only, not CMI):
#                  gdd_x_fire, gdd_x_pest, gdd_x_harv_light, gdd_x_harv_heavy, gdd_x_harv_any
#                  = z_gdd_3 * flag.  Allowed wherever the disturbance main term exists
#                  (no extra event threshold).
#   TSD (one per disturbance TYPE): tsd_fire, tsd_pest, tsd_harv (any harvest).
#                  RAW time since the disturbance in DECADES, not centred, not capped.
#                  0 = immediately after the disturbance. Structural zero where the
#                  flag is 0 (built in 04_prepare_inla.R). Wherever the disturbance main
#                  term exists the TSD term exists too.
#   INLA weibullsurv variant 1: S(t) = exp(-(lambda t)^alpha), lambda = exp(eta); the
#   hazard ratio of a coefficient is exp(alpha * beta), NOT exp(beta).
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

add_model_covariates <- function(d, scaling) {
  d <- data.table::copy(d)
  for (i in seq_len(nrow(scaling))) {
    v <- scaling$var[i]
    data.table::set(d, j = paste0("z_", v), value = (d[[v]] - scaling$mean[i]) / scaling$sd[i])
  }
  d[, `:=`(fire_f = as.integer(fire == 1L), pest_f = as.integer(pest == 1L),
           harv_light = as.integer(harv == 1L), harv_heavy = as.integer(harv == 2L),
           harv_any = as.integer(harv >= 1L))]
  # GDD x disturbance (":" products; main effects are separate terms)
  stopifnot("z_gdd_3" %in% names(d))
  d[, `:=`(gdd_x_fire = z_gdd_3 * fire_f, gdd_x_pest = z_gdd_3 * pest_f,
           gdd_x_harv_light = z_gdd_3 * harv_light, gdd_x_harv_heavy = z_gdd_3 * harv_heavy,
           gdd_x_harv_any = z_gdd_3 * harv_any)]
  # TSD columns (decades, structural zero where the flag is 0) come from 04; check them.
  for (v in c("tsd_fire", "tsd_pest", "tsd_harv")) stopifnot(v %in% names(d), !anyNA(d[[v]]), all(d[[v]] >= 0))
  stopifnot(all(d$tsd_fire[d$fire_f == 0L] == 0), all(d$tsd_pest[d$pest_f == 0L] == 0),
            all(d$tsd_harv[d$harv_any == 0L] == 0))
  d
}

# All possible model terms (columns of the beta matrix).
all_terms <- function(cont = COV_CONT)
  c(paste0("z_", cont), "fire_f", "pest_f", "harv_light", "harv_heavy", "harv_any",
    "gdd_x_fire", "gdd_x_pest", "gdd_x_harv_light", "gdd_x_harv_heavy", "gdd_x_harv_any",
    "tsd_fire", "tsd_pest", "tsd_harv")

# ------------------------------------------------------------------------------
# term_mask()
# Logical matrix, one row per KEPT transition, one column per possible term:
# TRUE = the term is in that transition's formula. Zeros in the beta matrix for
# FALSE entries (same role as zeros in the simulation's beta_mat).
# ------------------------------------------------------------------------------
term_mask <- function(tr, cont = COV_CONT) {
  k <- tr[kept == TRUE]
  terms <- all_terms(cont)
  m <- matrix(FALSE, nrow(k), length(terms), dimnames = list(paste(k$from, k$to, sep = "->"), terms))
  m[, paste0("z_", cont)] <- TRUE                       # continuous terms: every kept transition
  m[, "fire_f"]     <- k$fire_term
  m[, "pest_f"]     <- k$pest_term
  m[, "harv_light"] <- k$harv_term == "3level"
  m[, "harv_heavy"] <- k$harv_term == "3level"
  m[, "harv_any"]   <- k$harv_term == "any"
  # GDD x disturbance: wherever the disturbance main term exists (no extra threshold)
  m[, "gdd_x_fire"]       <- k$fire_term
  m[, "gdd_x_pest"]       <- k$pest_term
  m[, "gdd_x_harv_light"] <- k$harv_term == "3level"
  m[, "gdd_x_harv_heavy"] <- k$harv_term == "3level"
  m[, "gdd_x_harv_any"]   <- k$harv_term == "any"
  # TSD: one per disturbance type, wherever that type's main term exists
  m[, "tsd_fire"] <- k$fire_term
  m[, "tsd_pest"] <- k$pest_term
  m[, "tsd_harv"] <- k$harv_term != "none"
  m
}

# RHS string per kept transition (named by "from->to"), for fit_inla_models().
rhs_formulas <- function(mask) {
  out <- apply(mask, 1, function(r) if (any(r)) paste(colnames(mask)[r], collapse = " + ") else "1")
  out
}
