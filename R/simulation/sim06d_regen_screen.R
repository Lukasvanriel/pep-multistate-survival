# ==============================================================================
# R/simulation/sim06d_regen_screen.R   (2026-11-11; header below describes sim06, which this extends)   RUN ON THE MACHINE WITH INLA, from the project root:
#   Rscript R/simulation/sim06_tsd_tau.R 2>&1 | tee data-output/simulation/sim06_run.log      (re-running skips finished fits)
#
# QUESTION: which time-since-disturbance (TSD) form fits best? Main model: effect(TSD) = main + b * log(TSD + 1) (unbounded: a fading
#   effect keeps going and changes sign). Alternatives, all with the same number of coefficients per term:
#     linear       main + b * TSD / 10                                  (S24; also unbounded)
#     exp(tau)     main + b * exp(-TSD / tau)                           (bounded: main = long-run effect, b = extra effect right after the
#                                                                         disturbance that fades on time scale tau; a sign change only if the
#                                                                         long-run effect has the opposite sign)
#   tau on a grid, ONE tau per disturbance type (fire, pest, light harvest, heavy harvest) shared by all transitions (user 2026-11-11),
#   chosen by WAIC summed over transitions; reported also by transition class (into regen / out of regen / between states).
# DESIGN
#   Pass 1, per transition x imputation (all 72 transitions with a TSD term; 5 used imputations each, fit_validity used == TRUE, lowest m):
#     log_all (= main model), linear_all, and per TSD term T: lin_T (only T linear) and exp_T_tau for every tau (only T exponential);
#     all other TSD terms stay log. Paired WAIC differences against log_all on identical rows, SE from the row-wise WAIC contributions.
#   Pass 2 (after pass 1 is complete): per transition x imputation, best_all = every TSD term at its type's best tau (pass-1 profile).
#   Rows are rebuilt from the stored imputations exactly as in R/sensitivity/11 (reproduction check against mi_raw first). The TSD columns
#   are re-transformed row by row: years = expm1(model column) (main form log1p), new column = flag * f(years); flag 0 rows stay 0.
#   INLA call as in 11 / 05a (variant WEIBULL_VARIANT, N(0, COEF_PRIOR_SD^2), flat intercept, pc.alphaw(ALPHA_PRIOR), adaptive, grid, VB off)
#   plus WAIC; no screen (errors are reported).
# CAVEATS: imputations are those of the main (log) model; tau is chosen from a grid (charge ~ +2 WAIC per tuned tau when comparing forms).
# OUTPUT data-output/simulation/sim06_tsd_tau/: fits/ (one .rds per transition x imputation x pass), sim06_waic.csv, sim06_profile.csv,
#   sim06_profile_by_class.csv, sim06_best_tau.csv, sim06_compare.csv, run_log.txt, progress.log.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here("R", "config.R")); options(width = 200)
N_WORKERS <- as.integer(Sys.getenv("SIM06_WORKERS", "6")); N_IMP <- as.integer(Sys.getenv("SIM06_NIMP", "10"))
ONLY      <- Sys.getenv("SIM06_ONLY", "")                     # pilot: comma-separated labels
TAU       <- c(2, 3.5, 6, 10, 15, 25, 40, 60, 100, 160)       # years (user 2026-11-11)
TSD_TERMS <- c(tsd_fire = "fire_f", tsd_pest = "pest_f", tsd_harv_light = "harv_light", tsd_harv_heavy = "harv_heavy")   # TSD column -> its flag
CHECK_M <- 1L; CHECK_LABELS <- c("black_spruce->fir", "white_birch->fir", "fir->white_birch"); CHECK_TOL <- 0.1
stopifnot(TSD_FORM == "log1p")
OUT <- here("data-output", "simulation", "sim06d_regen_screen"); dir.create(file.path(OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
suppressPackageStartupMessages(library(INLA)); RNGkind("L'Ecuyer-CMRG"); set.seed(2026L)
logf <- file.path(OUT, "run_log.txt"); progf <- file.path(OUT, "progress.log")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
prog <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim06d_regen_screen.R")
logl("==== start | R %s | INLA %s | run %s | imputations per transition %d | workers %d | tau %s", getRversion(), packageVersion("INLA"), RUN_TAG, N_IMP, N_WORKERS, paste(TAU, collapse = " "))
logl("script %s md5 %s", self, md5sum(self))

# ---- inputs, transitions, imputations ----
inp <- prepare_fit_inputs(TABLE_DIR); mask <- inp$mask
terms_of <- function(lb) colnames(mask)[mask[lb, ]]
labs <- rownames(mask)[sapply(rownames(mask), function(lb) any(names(TSD_TERMS) %in% terms_of(lb)))]
labs <- intersect(labs, fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[class == "to_regen", label])   # sim06d: into-regen only
if (nzchar(ONLY)) labs <- intersect(labs, strsplit(ONLY, ",")[[1]])
fv <- fread(file.path(RESULTS_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
use <- fv[used == TRUE & label %in% labs][order(label, m)][, head(.SD, N_IMP), by = label][, .(label, m)]
stopifnot(all(use[, .N, by = label]$N == N_IMP))
ar <- fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[, .(label, class)]
logl("transitions %d | TSD terms %d | imputations used %s", length(labs), sum(sapply(labs, function(lb) sum(names(TSD_TERMS) %in% terms_of(lb)))), paste(sort(unique(use$m)), collapse = " "))

extras_file <- function(m) file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))
rebuild_rows <- function(m, x = inp) {                       # as R/sensitivity/11
  ex <- readRDS(extras_file(m)); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  fe <- setNames(p1$entry, as.character(p1$id_pe)); key_iv <- paste(x$iv$id_pe, x$iv$interval)
  ex_row <- rep(NA_real_, nrow(x$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii), nrow(p2) == sum(x$iv$is_change))
  ex_row[ii] <- p2$exit_abs; b_abs <- x$iv$t_end
  p2_fun <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  rows <- build_rows(x$iv, fe, x$states, x$tr, p2_fun = p2_fun, keep_cols = x$keep_cols)$rows
  rows <- rows[kept == TRUE]; rows[, label := paste(from, to, sep = "->")]; rows
}
fit_one <- function(d, tms, waic = TRUE) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(c("1", tms), collapse = " + ")))
  m <- tryCatch(inla(f, family = "weibullsurv", data = d,
                control.family = list(variant = WEIBULL_VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                control.fixed = list(mean = 0, prec = 1 / COEF_PRIOR_SD^2, mean.intercept = 0, prec.intercept = 0),
                control.compute = list(waic = waic),
                control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(m$error)) return(list(error = substr(m$error, 1, 150)))
  fx <- m$summary.fixed
  list(error = NA_character_, waic = if (waic) m$waic$waic else NA_real_, p_eff = if (waic) m$waic$p.eff else NA_real_, local = if (waic) m$waic$local.waic else NULL,
       coef = data.table(parameter = c("alpha", rownames(fx)), mean = c(m$summary.hyperpar[ALPHA_NAME, "mean"], fx$mean), sd = c(m$summary.hyperpar[ALPHA_NAME, "sd"], fx$sd)))
}
# TSD column transform: years from the main (log1p) column; flag 0 rows stay 0
retrans <- function(d, term, form, tau = NA_real_) {
  fl <- d[[TSD_TERMS[[term]]]]; stopifnot(all(d[[term]][fl == 0] == 0)); yrs <- expm1(d[[term]])
  v <- switch(form, log = d[[term]], linear = yrs / 10, exp = exp(-yrs / tau)); v[fl == 0] <- 0; set(d, j = term, value = v); invisible(d)
}


# ==============================================================================
# sim06d: the into-regen transitions of sim06 refitted with the MAIN RUN'S SCREEN (08_final_pool.R / 07_redraw.R): a fit is valid if
#   (1) no INLA error, (2) alpha-SD ratio sd / (0.78 * alpha / sqrt(events)) >= 0.1, (3) |z_map_alpha| <= 3 and |z_map_intercept| <= 3
#   (INLA mode vs direct maximisation of the same log-posterior, cap_diagnostics in 05a). 10 used imputations per transition.
#   Forms as sim06 pass 1 (log_all, lin_T, exp_T_tau; other TSD terms log). A WAIC difference (form minus log_all) is used only where BOTH
#   fits are valid; per transition the mean over valid imputations; profile = sum over transitions; tau chosen per disturbance type for
#   this class. Output data-output/simulation/sim06d_regen_screen/: fits/, sim06d_validity.csv, sim06d_waic.csv, sim06d_profile.csv.
# ==============================================================================
RATIO_MIN <- 0.1; Z_MAP_MAX <- 3
screen_fit <- function(r, d, tms) {                       # returns list(valid, reason, ratio, z_alpha, z_icpt)
  if (!is.na(r$error)) return(list(valid = FALSE, reason = "error", ratio = NA_real_, z_alpha = NA_real_, z_icpt = NA_real_))
  tb <- r$coef[, .(parameter = sub("^\\(Intercept\\)$", "(Intercept)", parameter), est = mean, sd = sd, fit_error = FALSE)]
  a <- tb[parameter == "alpha"]; ratio <- a$sd / (0.78 * a$est / sqrt(sum(d$status)))
  cap <- cap_diagnostics(as.data.frame(d), tms, tb)
  za <- unname(cap["z_map_alpha"]); zi <- unname(cap["z_map_intercept"])
  reason <- if (anyNA(c(ratio, za, zi))) "no_map" else if (ratio < RATIO_MIN) "alpha_sd_ratio" else if (za > Z_MAP_MAX) "z_map_alpha" else if (zi > Z_MAP_MAX) "z_map_intercept" else "ok"
  list(valid = reason == "ok", reason = reason, ratio = ratio, z_alpha = za, z_icpt = zi)
}
fit_forms_s <- function(d, tms) {
  tt <- intersect(names(TSD_TERMS), tms); forms <- list(log_all = list())
  for (T in tt) { forms[[paste0("lin_", T)]] <- setNames("linear", T)
    for (k in seq_along(TAU)) forms[[sprintf("exp_%s_%03d", T, k)]] <- setNames(sprintf("exp:%g", TAU[k]), T) }
  lapply(forms, function(spec) { dd <- copy(d)
    for (T in names(spec)) { s <- spec[[T]]; if (s == "linear") retrans(dd, T, "linear") else retrans(dd, T, "exp", as.numeric(sub("exp:", "", s))) }
    r <- fit_one(as.data.frame(dd), tms); r$spec <- spec; r$screen <- screen_fit(r, dd, tms); r$local_keep <- r$local; r })
}
todo <- use[!file.exists(file.path(OUT, "fits", sprintf("m%02d_%s.rds", m, gsub("->", "_to_", label))))]
prog("[sim06d] %d of %d transition x imputation jobs to run (%d transitions)", nrow(todo), nrow(use), length(labs)); t0 <- Sys.time(); done <- 0L
for (mm in sort(unique(todo$m))) { lb_m <- todo[m == mm, label]; rows <- rebuild_rows(mm); sp <- split(rows[label %in% lb_m], by = "label", keep.by = TRUE); rm(rows); invisible(gc())
  res <- mclapply(lb_m, function(lb) { tms <- terms_of(lb); d <- sp[[lb]][, c("status", "Tstop", "truncation", tms), with = FALSE]
    fits <- fit_forms_s(d, tms); saveRDS(list(label = lb, m = mm, terms = tms, n_events = sum(d$status), fits = fits), file.path(OUT, "fits", sprintf("m%02d_%s.rds", mm, gsub("->", "_to_", lb))))
    sum(!vapply(fits, function(f) f$screen$valid, logical(1))) }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("worker error: ", paste(res[bad], collapse = " | "))
  done <- done + length(lb_m); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  prog("[sim06d] m%02d done (%d jobs, fits failing the screen %d) | %d/%d | elapsed %.0f min, ~%.0f min left", mm, length(lb_m), sum(unlist(res)), done, nrow(todo), el, el / done * (nrow(todo) - done)) }

# ---- summaries ----
FL <- file.path(OUT, "fits", sprintf("m%02d_%s.rds", use$m, gsub("->", "_to_", use$label)))
V <- rbindlist(lapply(FL, function(f) { r <- readRDS(f); rbindlist(lapply(names(r$fits), function(n) { s <- r$fits[[n]]$screen
  data.table(label = r$label, m = r$m, form = n, valid = s$valid, reason = s$reason, ratio = s$ratio, z_alpha = s$z_alpha, z_icpt = s$z_icpt, waic = if (is.null(r$fits[[n]]$waic)) NA_real_ else r$fits[[n]]$waic) })) }))   # failed fits have no WAIC
fwrite(V, file.path(OUT, "sim06d_validity.csv"))
W <- rbindlist(lapply(FL, function(f) { r <- readRDS(f); b <- r$fits$log_all
  rbindlist(lapply(setdiff(names(r$fits), "log_all"), function(n) { x <- r$fits[[n]]; ok <- b$screen$valid && x$screen$valid
    dl <- if (ok) x$local - b$local else NA_real_
    data.table(label = r$label, m = r$m, form = n, term = sub("^(exp|lin)_(tsd_[a-z_]+?)(_[0-9]{3})?$", "\\2", n),
               tau = if (grepl("^exp_", n)) TAU[as.integer(sub(".*_", "", n))] else NA_real_, both_valid = ok,
               dwaic = if (ok) x$waic - b$waic else NA_real_, var_d = if (ok) length(dl) * var(dl) else NA_real_) })) }))
fwrite(W, file.path(OUT, "sim06d_waic.csv"))
PT <- W[, .(n_valid = sum(both_valid), dwaic = mean(dwaic, na.rm = TRUE), var_d = mean(var_d, na.rm = TRUE)), by = .(label, form, term, tau)]
P <- PT[, .(n_transitions = sum(n_valid > 0), n_pairs_valid = sum(n_valid), dwaic_sum = if (any(n_valid > 0)) sum(dwaic, na.rm = TRUE) else NA_real_,
             se = if (any(n_valid > 0)) sqrt(sum(var_d, na.rm = TRUE)) else NA_real_), by = .(term, form, tau)]   # NA, not 0, when nothing is valid
P[, shape := fifelse(grepl("^lin_", form), "linear", "exp")]; fwrite(P, file.path(OUT, "sim06d_profile.csv"))
cat("\nScreen: share of fits valid, by form type:\n"); print(V[, .(fits = .N, valid = round(mean(valid), 3)), by = .(kind = sub("_tsd.*", "", form))])
cat("\nReasons for failing:\n"); print(V[valid == FALSE, .N, by = reason])
cat("\nProfile into regen (valid pairs only; WAIC vs log summed over transitions; negative = better than log):\n")
print(dcast(P[shape == "exp"], term ~ tau, value.var = "dwaic_sum", fun.aggregate = function(x) round(x, 1)))
cat("\nLinear vs log:\n"); print(P[shape == "linear", .(term, n_transitions, n_pairs_valid, dwaic_sum = round(dwaic_sum, 1), se = round(se, 1))])
cat("\nValid pairs per tau (out of transitions x 10):\n"); print(dcast(P[shape == "exp"], term ~ tau, value.var = "n_pairs_valid"))
logl("==== done")
