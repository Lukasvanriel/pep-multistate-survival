# ==============================================================================
# 11_within_between_climate.R   RUN ON THE MACHINE WITH INLA (from the project root: Rscript R/sensitivity/11_within_between_climate.R)
#
# QUESTION (concern 2, spatial confounding of climate): is the GDD / CMI effect a climate effect, or a
# north-south gradient of other things (seed sources, land use, ...)? 95% of GDD variance is between plots.
#
# DESIGN (within/between, "Mundlak" decomposition). Each continuous climate covariate is split into
#   *_btw = (plot mean over the plot's intervals - overall mean) / overall SD      -> warm vs cold PLACES (spatial)
#   *_wth = (interval value - plot mean)                       / overall SD      -> warmer vs cooler PERIODS in the same plot (temporal)
#   *_btw + *_wth = z_* exactly, and both are in the SAME unit (1 SD of the original covariate), so their coefficients are directly comparable.
# The within-plot term cannot be confounded by anything that is constant per plot (seed sources, soils, region, land-use regime).
# If beta_wth ~ beta_btw: supports a genuine climate effect. If beta_wth ~ 0 while beta_btw is large: the gradient effect is not reproduced over time.
#
# Caveats to report with the result:
#   * within-plot variation is small (GDD 4.6%, CMI 2.6% of the variance), so beta_wth has wide intervals; and it reflects decadal fluctuation
#     (incl. the warming trend), not long-term climate, so a smaller beta_wth does NOT by itself prove confounding (slow responses / lags).
#   * GDD x disturbance interactions are kept as in the main model (z_gdd_3 x flag, total GDD); soil is unchanged (no within-plot variance).
#   * (WP7, D7) transitions without a GDD main term (A11, z_gdd_3 not in the model) get neither gdd_btw nor gdd_wth (nor d_gdd); the CMI split stays.
#   * the imputations (P1 entries, P2 exits) are REUSED from the main run (extras_<tag>_mXX.rds): they were generated under the main model.
#     Events and their intervals are observed, only their timing inside the bracket is imputed, so this mainly affects timing, not which
#     plots changed. State this as a limitation (sensitivity analysis, not a full re-imputation).
#
# WHAT IT DOES
#   0. Rebuilds the rows of imputation CHECK_M from the stored P1/P2 and refits CHECK_LABELS with the ORIGINAL formula: estimates must
#      reproduce mi_raw_<tag>.csv (stops otherwise). This verifies both the row rebuild and the INLA call.
#   1. Pass 1: refits every transition on its 30 USED imputations (fit_validity_<tag>.csv, used == TRUE) with the within/between formula.
#      Same INLA call as fit_inla_models_pc in 05a (weibullsurv variant WEIBULL_VARIANT = 0, control.fixed N(0, COEF_PRIOR_SD^2) with a flat intercept, adaptive, grid, VB off, config TRUE), plus linear combinations
#      d_gdd = gdd_wth - gdd_btw (where GDD is in the model) and d_cmi = cmi_wth - cmi_btw (derived only; they do not change the fit). Estimates are log hazard ratios: HR = exp(beta); the intercept stays on INLA's decade scale.
#   2. Screen (same rules as 08_final_pool.R, user decision 2026-11-05): no error/NA; alpha-SD ratio >= 0.1; |z_map_alpha| <= 3 and |z_map_intercept| <= 3
#      (INLA mode vs. direct maximisation of the same log-posterior, cap_diagnostics in 05a). No reference table.
#   3. Pass 2: a transition with < 30 valid fits gets replacement fits from the other stored imputations, in increasing m.
#   4. Rubin pooling (pool_all from 05a_fit_funs.R) and the comparison table.
# OUTPUTS (SENS_DIR/wb = data-output/<RUN_TAG>/sensitivity/wb/; inttest version in data-output/inttest/sensitivity/wb/): raw_wb_p1_mXXX.csv / raw_wb_p2_mXXX.csv (per imputation, resumable), validity_wb.csv,
#   pooled_wb.csv, compare_within_between_wb.csv (one row per transition), summary printed.
# PILOT: set PILOT_M <- c(1, 2, 3) to run step 0 and pass 1 for those imputations only and print a time estimate for the full run.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
options(width = 200)

# ---- config ------------------------------------------------------------------
source(here::here("R", "config.R"))   # RUN_TAG (run whose imputations are reused), TABLE_DIR, FIT_DIR, RESULTS_DIR, SENS_DIR
IN_DIR    <- FIT_DIR
OUT_DIR   <- file.path(SENS_DIR, "wb"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
N_WORKERS <- 6L
TARGET    <- 30L; RATIO_MIN <- 0.1; Z_MAP_MAX <- 3
PILOT_M   <- NULL                                    # e.g. c(1, 2, 3) for a timing pilot
CHECK_M   <- 1L
CHECK_LABELS <- c("black_spruce->fir", "white_birch->fir", "fir->white_birch")   # must be VALID fits in CHECK_M (a collapsed original cannot be reproduced)
CHECK_TOL <- 0.1                                     # max |difference| in posterior-SD units for the reproduction check
SEED <- 2026L

SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L   # set before sourcing 05a (its defaults are only used if unset)
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)

# ---- inputs and the within/between covariates --------------------------------------------------------------------------------------------------------
inp <- prepare_fit_inputs(TABLE_DIR)          # same tables, same scaling, same ordering as 05_fit_inla.R
iv  <- inp$iv; sc <- inp$scaling
decomp <- function(v, short) {
  m <- sc[var == v, mean]; s <- sc[var == v, sd]
  iv[, pm__ := mean(get(v)), by = id_pe]
  iv[, (paste0(short, "_btw")) := (pm__ - m) / s]
  iv[, (paste0(short, "_wth")) := (get(v) - pm__) / s]
  iv[, pm__ := NULL]
  stopifnot(isTRUE(all.equal(iv[[paste0(short, "_btw")]] + iv[[paste0(short, "_wth")]], iv[[paste0("z_", v)]])))
}
decomp("gdd_3", "gdd"); decomp("cmi_gs", "cmi")
cat("Variance (in SD units of the original covariate):\n")
print(iv[, .(gdd_btw_sd = sd(gdd_btw), gdd_wth_sd = sd(gdd_wth), cmi_btw_sd = sd(cmi_btw), cmi_wth_sd = sd(cmi_wth))])

m0 <- inp$mask
new_cont <- c("gdd_btw", "gdd_wth", "cmi_btw", "cmi_wth")
has_g <- m0[, "z_gdd_3"]                                       # A11: FALSE where z_gdd_3 is not in the model -> no gdd_btw / gdd_wth there
mask_wb <- cbind(gdd_btw = has_g, gdd_wth = has_g, cmi_btw = TRUE, cmi_wth = TRUE,
                 m0[, setdiff(colnames(m0), c("z_gdd_3", "z_cmi_gs")), drop = FALSE])
inp_wb <- inp
inp_wb$iv <- iv
inp_wb$mask <- mask_wb
inp_wb$terms <- colnames(mask_wb)
inp_wb$rhs <- rhs_formulas(mask_wb)
inp_wb$keep_cols <- unique(c("fire", "harv", "pest", "year_end", inp$terms, inp_wb$terms))
terms_of <- function(mask, lb) colnames(mask)[mask[lb, ]]

# ---- rows of a stored imputation (P1 entries and P2 exits from extras_int_mXX.rds) --------------------------------------------------------------------
extras_file <- function(m) file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))
rebuild_rows <- function(m, x = inp_wb) {
  ex <- readRDS(extras_file(m)); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  fe <- setNames(p1$entry, as.character(p1$id_pe))
  key_iv <- paste(x$iv$id_pe, x$iv$interval)
  ex_row <- rep(NA_real_, nrow(x$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv)
  stopifnot(!anyNA(ii), nrow(p2) == sum(x$iv$is_change))
  ex_row[ii] <- p2$exit_abs
  b_abs <- x$iv$t_end
  p2_fun <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }   # stored exit, clamped to rounding
  rows <- build_rows(x$iv, fe, x$states, x$tr, p2_fun = p2_fun, keep_cols = x$keep_cols)$rows
  rows <- rows[kept == TRUE]; rows[, label := paste(from, to, sep = "->")]
  rows
}

# ---- INLA call: identical settings to the main fits (fit_inla_models_pc in 05a, exact = TRUE: variant WEIBULL_VARIANT, alpha prior pc.alphaw(ALPHA_PRIOR),
# control.fixed N(0, COEF_PRIOR_SD^2) with a flat intercept, one thread), optional derived linear combinations ---------------------------
lc_names <- function(tms, with_lc) if (!with_lc) character(0) else c(if (all(c("gdd_btw", "gdd_wth") %in% tms)) "d_gdd", "d_cmi")
fit_one <- function(d, rhs, with_lc = FALSE, tms = character(0)) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs))
  lc <- NULL; ln <- lc_names(tms, with_lc)
  if (length(ln)) { l <- list(); if ("d_gdd" %in% ln) { l1 <- inla.make.lincomb(gdd_wth = 1, gdd_btw = -1); names(l1) <- "d_gdd"; l <- c(l, l1) }
                    l2 <- inla.make.lincomb(cmi_wth = 1, cmi_btw = -1); names(l2) <- "d_cmi"; lc <- c(l, l2) }
  tryCatch(inla(f, family = "weibullsurv", data = d, lincomb = lc,
                control.family = list(variant = WEIBULL_VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                control.fixed = list(mean = 0, prec = 1 / COEF_PRIOR_SD^2, mean.intercept = 0, prec.intercept = 0),
                control.compute = list(config = TRUE),
                control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                num.threads = "1:1"),
           error = function(e) list(error = conditionMessage(e)))
}
tidy_one <- function(mod, tms, lb, m, n_ev, with_lc) {
  pars <- c("alpha", "(Intercept)", tms, lc_names(tms, with_lc))
  if (!is.null(mod$error)) return(data.table(label = lb, m = m, parameter = pars, est = NA_real_, sd = NA_real_, fit_error = TRUE, n_events = n_ev, error_msg = mod$error[1]))
  fx <- mod$summary.fixed; hy <- mod$summary.hyperpar; lcs <- mod$summary.lincomb.derived
  g <- function(tab, nm, col) if (!is.null(tab) && nm %in% rownames(tab)) tab[nm, col] else NA_real_
  est <- c(g(hy, ALPHA_NAME, "mean"), g(fx, "(Intercept)", "mean"), vapply(tms, function(t) g(fx, t, "mean"), 0),   # intercept on INLA's decade scale (variant 0, B3)
           vapply(lc_names(tms, with_lc), function(t) g(lcs, t, "mean"), 0))
  sdv <- c(g(hy, ALPHA_NAME, "sd"), g(fx, "(Intercept)", "sd"), vapply(tms, function(t) g(fx, t, "sd"), 0),
           vapply(lc_names(tms, with_lc), function(t) g(lcs, t, "sd"), 0))
  data.table(label = lb, m = m, parameter = pars, est = est, sd = sdv, fit_error = FALSE, n_events = n_ev, error_msg = NA_character_)
}
fit_labels <- function(m, labels, mask, rhs, with_lc) {
  rows <- rebuild_rows(m); sp <- split(rows, by = "label", keep.by = TRUE); rm(rows)
  rbindlist(lapply(labels, function(lb) {
    tms <- terms_of(mask, lb); d <- as.data.frame(sp[[lb]][, c("status", "Tstop", "truncation", tms), with = FALSE])
    tb <- tidy_one(fit_one(d, rhs[[lb]], with_lc, tms), tms, lb, m, sum(d$status == 1L), with_lc)
    cap <- cap_diagnostics(d, tms, tb[parameter %in% c("alpha", "(Intercept)", tms)])   # B6 as in 05/07 (05a): MAP vs INLA mode; z_map_alpha / z_map_intercept feed the screen
    tb[, (names(cap)) := as.list(cap)][]
  }))
}

# ---- 0. reproduction check: original formula on stored imputation CHECK_M must reproduce mi_raw ----------------------------------------------------
mi_raw <- fread(file.path(IN_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
fv0 <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
CHECK_LABELS <- intersect(CHECK_LABELS, fv0[m == CHECK_M & valid == TRUE, label]); stopifnot(length(CHECK_LABELS) > 0)
chk <- fit_labels(CHECK_M, CHECK_LABELS, inp$mask, inp$rhs, FALSE)
cmp <- merge(chk[, .(label, parameter, est_new = est, n_new = n_events)], mi_raw[m == CHECK_M, .(label, parameter, est, sd, n_events)], by = c("label", "parameter"))
cmp[, dz := abs(est_new - est) / sd]
cat("\nStep 0, reproduction of mi_raw (imputation", CHECK_M, "):\n"); print(cmp[, .(max_dz = round(max(dz), 4), events_equal = all(n_new == n_events)), by = label])
if (nrow(cmp) == 0 || anyNA(cmp$dz) || any(cmp$dz > CHECK_TOL) || any(cmp$n_new != cmp$n_events)) stop("Row rebuild / INLA call does not reproduce the main run; do not continue.")

# ---- 1. pass 1: the 30 used imputations per transition -------------------------------------------------------------------------------------------------
fv <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
avail_m <- sort(as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG)))))
p1_pairs <- fv[used == TRUE, .(label, m)]
if (!is.null(PILOT_M)) p1_pairs <- p1_pairs[m %in% PILOT_M]
run_pass <- function(pairs, pass) {
  todo <- sort(unique(pairs$m)); done <- todo[file.exists(file.path(OUT_DIR, sprintf("raw_wb_%s_m%03d.csv", pass, todo)))]
  todo <- setdiff(todo, done); cat(sprintf("\n%s: %d imputations to fit (%d already done), %d fits\n", pass, length(todo), length(done), nrow(pairs)))
  t0 <- proc.time()["elapsed"]; labs_by_m <- split(pairs$label, pairs$m)
  for (batch in split(todo, ceiling(seq_along(todo) / N_WORKERS))) {
    out <- mclapply(batch, function(mi) tryCatch(fit_labels(mi, labs_by_m[[as.character(mi)]], inp_wb$mask, inp_wb$rhs, TRUE),
                                                error = function(e) data.table(m = mi, run_error = conditionMessage(e))), mc.cores = N_WORKERS)
    for (r in out) { if ("run_error" %in% names(r)) { cat("  m", r$m, "FAILED:", r$run_error, "\n"); next }
                     fwrite(r, file.path(OUT_DIR, sprintf("raw_wb_%s_m%03d.csv", pass, r$m[1]))) }
    cat(sprintf("  batch up to m=%d done (%.0f s)\n", max(batch), proc.time()["elapsed"] - t0))
  }
  rbindlist(lapply(list.files(OUT_DIR, pattern = sprintf("^raw_wb_%s_m\\d+\\.csv$", pass), full.names = TRUE), fread), fill = TRUE)
}
t_p1 <- proc.time()["elapsed"]
raw1 <- run_pass(p1_pairs, "p1")
t_p1 <- proc.time()["elapsed"] - t_p1
if (!is.null(PILOT_M)) {
  n_full <- uniqueN(fv[used == TRUE, m]); fits_full <- fv[used == TRUE, .N]
  cat(sprintf("\nPILOT: %d fits in %.0f s with %d workers -> full pass 1 (%d fits over %d imputations) roughly %.1f h, plus pass 2 replacements.\n",
              nrow(unique(raw1[, .(label, m)])), t_p1, N_WORKERS, fits_full, n_full, t_p1 / nrow(unique(raw1[, .(label, m)])) * fits_full / 3600))
  quit(save = "no")
}

# ---- 2. screen (same rules as 08_final_pool.R, user decision 2026-11-05) ---------------------------------------------------------------------------------------------------------------------------------------
screen <- function(raw) {
  stopifnot("screen needs z_map_alpha and z_map_intercept in the fit table (cap_diagnostics, 05a)" = all(c("z_map_alpha", "z_map_intercept") %in% names(raw)))
  a <- raw[parameter == "alpha", .(label, m, alpha = est, a_sd = sd, n_events, fit_error, z_map_alpha, z_map_intercept)]
  i <- raw[parameter == "(Intercept)", .(label, m, icpt = est, i_sd = sd)]
  x <- merge(a, i, by = c("label", "m")); x[, ratio := a_sd / (0.78 * alpha / sqrt(n_events))]
  x[, reason := fifelse(fit_error %in% TRUE | is.na(alpha) | is.na(a_sd) | is.na(icpt), "error/NA",
                fifelse(ratio < RATIO_MIN, "alpha_sd_ratio",
                fifelse(is.na(z_map_alpha) | is.na(z_map_intercept) | abs(z_map_alpha) > Z_MAP_MAX | abs(z_map_intercept) > Z_MAP_MAX, "z_map", "ok")))]
  x[, valid := reason == "ok"]; list(x = x)
}
s1 <- screen(raw1)
cat("\nPass 1 screen:\n"); print(s1$x[, .N, by = reason])

# ---- 3. pass 2: replacements from other stored imputations, increasing m ---------------------------------------------------------------------------------
raw_all <- copy(raw1); val <- s1$x[, .(label, m, valid, reason, pass = "p1")]
for (round in 1:20) {
  need <- val[, .(n_valid = sum(valid)), by = label][n_valid < TARGET]
  if (!nrow(need)) break
  pairs <- rbindlist(lapply(seq_len(nrow(need)), function(j) { lb <- need$label[j]; cand <- setdiff(avail_m, val[label == lb, m])
    if (!length(cand)) return(NULL); data.table(label = lb, m = head(cand, TARGET - need$n_valid[j])) }))
  if (is.null(pairs) || !nrow(pairs)) { cat("No candidates left for:", paste(need$label, collapse = ", "), "\n"); break }
  tag <- sprintf("p2r%02d", round); r2 <- run_pass(pairs, tag); r2 <- r2[pairs, on = .(label, m), nomatch = 0]
  s2 <- screen(r2); raw_all <- rbind(raw_all, r2, fill = TRUE); val <- rbind(val, s2$x[, .(label, m, valid, reason, pass = tag)])
}
setorder(val, label, pass, m)                                   # pass-1 fits first, then replacements in increasing m
val[, rank_valid := ifelse(valid, cumsum(valid), NA_integer_), by = label]; val[, used := valid & rank_valid <= TARGET]
fwrite(val, file.path(OUT_DIR, "validity_wb.csv"))
short <- val[, .(used = sum(used)), by = label][used < TARGET]
if (nrow(short)) { cat("\n!! transitions with fewer than", TARGET, "valid fits (pooled over what exists):\n"); print(short) }
cat("\nInvalid fits by reason (all passes):\n"); print(val[valid == FALSE, .N, by = .(label, reason)][order(-N)])

# ---- 4. pooling and comparison -------------------------------------------------------------------------------------------------------------------------------
fin <- raw_all[val[used == TRUE, .(label, m)], on = .(label, m)]
fin[, k := match(label, inp$td$label)]
pooled <- pool_all(fin); fwrite(pooled, file.path(OUT_DIR, "pooled_wb.csv"))

w <- dcast(pooled[parameter %in% c("alpha", new_cont, "d_gdd", "d_cmi")], label ~ parameter, value.var = c("pooled_est", "se", "ci_lo", "ci_hi"))
main <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))[parameter %in% c("z_gdd_3", "z_cmi_gs"), .(label, parameter, est_main = pooled_est)]
main <- dcast(main, label ~ parameter, value.var = "est_main")
w <- merge(w, main, by = "label")
for (v in c("gdd", "cmi")) {
  b <- paste0(v, "_btw"); wi <- paste0(v, "_wth"); d <- paste0("d_", v)
  w[, (paste0("hr_", b)) := exp(get(paste0("pooled_est_", b)))]                                  # variant 0: HR = exp(beta)
  w[, (paste0("hr_", wi)) := exp(get(paste0("pooled_est_", wi)))]
  w[, (paste0("hr_", wi, "_lo")) := exp(get(paste0("ci_lo_", wi)))]; w[, (paste0("hr_", wi, "_hi")) := exp(get(paste0("ci_hi_", wi)))]
  w[, (paste0(v, "_btw_sig")) := get(paste0("ci_lo_", b)) > 0 | get(paste0("ci_hi_", b)) < 0]
  w[, (paste0(v, "_wth_sig")) := get(paste0("ci_lo_", wi)) > 0 | get(paste0("ci_hi_", wi)) < 0]
  w[, (paste0(v, "_same_sign")) := sign(get(paste0("pooled_est_", b))) == sign(get(paste0("pooled_est_", wi)))]
  w[, (paste0(v, "_diff_sig")) := get(paste0("ci_lo_", d)) > 0 | get(paste0("ci_hi_", d)) < 0]
}
w <- merge(w, inp$td[, .(label, n_events)], by = "label"); setorder(w, -n_events)
fwrite(w, file.path(OUT_DIR, "compare_within_between_wb.csv"))

cat("\n=== Within vs between, per transition (per 1 SD of the ORIGINAL covariate; HR = exp(beta); NA GDD rows = transitions without a GDD main term) ===\n")
for (v in c("gdd", "cmi")) {
  cat(sprintf("\n-- %s --\n", toupper(v)))
  print(w[, .(label, n_events, hr_btw = round(get(paste0("hr_", v, "_btw")), 2), hr_wth = round(get(paste0("hr_", v, "_wth")), 2),
              wth_lo = round(get(paste0("hr_", v, "_wth_lo")), 2), wth_hi = round(get(paste0("hr_", v, "_wth_hi")), 2),
              btw_sig = get(paste0(v, "_btw_sig")), wth_sig = get(paste0(v, "_wth_sig")), same_sign = get(paste0(v, "_same_sign")),
              diff_sig = get(paste0(v, "_diff_sig")))])
  s <- w[get(paste0(v, "_btw_sig")) == TRUE]
  cat(sprintf("Transitions with a clear between-plot effect: %d | within-plot same sign: %d | within-plot CI excludes 0: %d | within differs from between (d CI excludes 0): %d\n",
              nrow(s), sum(s[[paste0(v, "_same_sign")]]), sum(s[[paste0(v, "_wth_sig")]]), sum(s[[paste0(v, "_diff_sig")]])))
  fit <- lm(as.formula(sprintf("pooled_est_%s_wth ~ pooled_est_%s_btw", v, v)), data = w, weights = 1 / w[[paste0("se_", v, "_wth")]]^2)
  cs <- summary(fit)$coefficients
  cat(sprintf("Across all transitions, weighted slope of within on between = %.2f (SE %.2f); 1 = same effect over time as across space, 0 = no temporal effect.\n",
              cs[2, 1], cs[2, 2]))
  cat(sprintf("Check: between-plot estimate vs main-run total estimate, median |difference| = %.3f\n",
              median(abs(w[[paste0("pooled_est_", v, "_btw")]] - w[[if (v == "gdd") "z_gdd_3" else "z_cmi_gs"]]), na.rm = TRUE)))
}
cat("\nOutputs in", OUT_DIR, "\n")
