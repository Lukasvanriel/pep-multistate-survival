# ==============================================================================
# R/simulation/sim06_tsd_tau.R   (2026-11-11)   RUN ON THE MACHINE WITH INLA, from the project root:
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
N_WORKERS <- as.integer(Sys.getenv("SIM06_WORKERS", "6")); N_IMP <- as.integer(Sys.getenv("SIM06_NIMP", "5"))
ONLY      <- Sys.getenv("SIM06_ONLY", "")                     # pilot: comma-separated labels
TAU       <- c(2, 3.5, 6, 10, 15, 25, 40, 60, 100, 160)       # years (user 2026-11-11)
TSD_TERMS <- c(tsd_fire = "fire_f", tsd_pest = "pest_f", tsd_harv_light = "harv_light", tsd_harv_heavy = "harv_heavy")   # TSD column -> its flag
CHECK_M <- 1L; CHECK_LABELS <- c("black_spruce->fir", "white_birch->fir", "fir->white_birch"); CHECK_TOL <- 0.1
stopifnot(TSD_FORM == "log1p")
OUT <- here("data-output", "simulation", "sim06_tsd_tau"); dir.create(file.path(OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
suppressPackageStartupMessages(library(INLA)); RNGkind("L'Ecuyer-CMRG"); set.seed(2026L)
logf <- file.path(OUT, "run_log.txt"); progf <- file.path(OUT, "progress.log")
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time()), s, "\n", file = logf, append = TRUE) }
prog <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, "\n", file = progf, append = TRUE) }
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (!length(self)) self <- here("R", "simulation", "sim06_tsd_tau.R")
logl("==== start | R %s | INLA %s | run %s | imputations per transition %d | workers %d | tau %s", getRversion(), packageVersion("INLA"), RUN_TAG, N_IMP, N_WORKERS, paste(TAU, collapse = " "))
logl("script %s md5 %s", self, md5sum(self))

# ---- inputs, transitions, imputations ----
inp <- prepare_fit_inputs(TABLE_DIR); mask <- inp$mask
terms_of <- function(lb) colnames(mask)[mask[lb, ]]
labs <- rownames(mask)[sapply(rownames(mask), function(lb) any(names(TSD_TERMS) %in% terms_of(lb)))]
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

# ---- 0. reproduction check (log_all = main model) ----
chk_f <- file.path(OUT, "check_ok.txt")
if (!file.exists(chk_f)) {
  mi_raw <- fread(file.path(FIT_DIR, sprintf("mi_raw_%s.csv", RUN_TAG))); lbs <- intersect(CHECK_LABELS, fv[m == CHECK_M & valid == TRUE, label])
  rows <- rebuild_rows(CHECK_M)
  dz <- rbindlist(lapply(lbs, function(lb) { tms <- terms_of(lb); r <- fit_one(as.data.frame(rows[label == lb, c("status", "Tstop", "truncation", tms), with = FALSE]), tms, waic = FALSE)
    x <- merge(r$coef[parameter %in% tms], mi_raw[label == lb & m == CHECK_M, .(parameter, est)], by = "parameter"); x[, .(label = lb, dz = max(abs(mean - est) / sd))] }))
  print(dz); if (any(dz$dz > CHECK_TOL)) stop("rebuild / INLA call does not reproduce the main run")
  writeLines(format(Sys.time()), chk_f); rm(rows); invisible(gc())
}

# ---- pass 1 / pass 2 ----
fit_forms <- function(d, tms, pass, best = NULL) {
  tt <- intersect(names(TSD_TERMS), tms); forms <- list()
  if (pass == 1L) {
    forms$log_all <- list(); forms$linear_all <- setNames(rep("linear", length(tt)), tt)
    for (T in tt) { forms[[paste0("lin_", T)]] <- setNames("linear", T)
      for (k in seq_along(TAU)) forms[[sprintf("exp_%s_%03d", T, k)]] <- setNames(sprintf("exp:%g", TAU[k]), T) }
  } else forms$best_all <- setNames(sprintf("exp:%g", best[tt]), tt)
  lapply(forms, function(spec) { dd <- copy(d)
    for (T in names(spec)) { s <- spec[[T]]; if (s == "linear") retrans(dd, T, "linear") else retrans(dd, T, "exp", as.numeric(sub("exp:", "", s))) }
    r <- fit_one(as.data.frame(dd), tms); r$spec <- spec; r })
}
run_pass <- function(pass, best = NULL) {
  tag <- if (pass == 1L) "p1" else paste0("p2_", paste(sprintf("%s%g", sub("tsd_", "", names(best)), best), collapse = "_"))   # p2 files carry the tau set they used
  todo <- use[!file.exists(file.path(OUT, "fits", sprintf("%s_m%02d_%s.rds", tag, m, gsub("->", "_to_", label))))]
  prog("[sim06] pass %d: %d of %d transition x imputation jobs to run", pass, nrow(todo), nrow(use)); t0 <- Sys.time(); done <- 0L
  for (mm in sort(unique(todo$m))) { lb_m <- todo[m == mm, label]; rows <- rebuild_rows(mm); sp <- split(rows[label %in% lb_m], by = "label", keep.by = TRUE); rm(rows); invisible(gc())
    res <- mclapply(lb_m, function(lb) { tms <- terms_of(lb); d <- sp[[lb]][, c("status", "Tstop", "truncation", tms), with = FALSE]
      fits <- fit_forms(d, tms, pass, best); saveRDS(list(label = lb, m = mm, terms = tms, n_events = sum(d$status), n_rows = nrow(d), fits = fits),
                                                       file.path(OUT, "fits", sprintf("%s_m%02d_%s.rds", tag, mm, gsub("->", "_to_", lb))))
      sum(!is.na(vapply(fits, `[[`, "", "error"))) }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
    bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("worker error: ", paste(res[bad], collapse = " | "))
    done <- done + length(lb_m); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    prog("[sim06] pass %d: imputation m%02d done (%d jobs, INLA errors %d) | %d/%d | elapsed %.0f min, ~%.0f min left", pass, mm, length(lb_m), sum(unlist(res)), done, nrow(todo), el, el / done * (nrow(todo) - done)) }
}
run_pass(1L)

# ---- pass-1 summaries: paired WAIC differences vs log_all, profile per type ----
F1 <- file.path(OUT, "fits", sprintf("p1_m%02d_%s.rds", use$m, gsub("->", "_to_", use$label)))
W <- rbindlist(lapply(F1, function(f) { r <- readRDS(f); b <- r$fits$log_all
  rbindlist(lapply(setdiff(names(r$fits), "log_all"), function(fm) { x <- r$fits[[fm]]; ok <- is.na(b$error) && is.na(x$error)
    dl <- if (ok) x$local - b$local else NA_real_
    data.table(label = r$label, m = r$m, form = fm, term = if (grepl("^(exp|lin)_", fm)) sub("^(exp|lin)_(tsd_[a-z_]+?)(_[0-9]{3})?$", "\\2", fm) else NA_character_,
               tau = if (grepl("^exp_", fm)) TAU[as.integer(sub(".*_", "", fm))] else NA_real_, dwaic = if (ok) x$waic - b$waic else NA_real_,
               var_d = if (ok) length(dl) * var(dl) else NA_real_, error = if (!ok) paste(b$error, x$error) else NA_character_) })) }))
W <- merge(W, ar, by = "label"); fwrite(W, file.path(OUT, "sim06_waic.csv"))
prof <- function(by) W[grepl("^(exp|lin)_", form), .(n_transitions = uniqueN(label), dwaic_sum = sum(dwaic, na.rm = TRUE) / N_IMP,
                                                      se = sqrt(sum(var_d, na.rm = TRUE) / N_IMP) , errors = sum(!is.na(error))), by = c("term", "form", "tau", by)]
P <- prof(NULL); P[, form_short := fifelse(grepl("^lin_", form), "linear", "exp")]; fwrite(P, file.path(OUT, "sim06_profile.csv"))
PC <- prof("class"); fwrite(PC, file.path(OUT, "sim06_profile_by_class.csv"))
BT <- P[form_short == "exp", .SD[which.min(dwaic_sum)], by = term][, .(term, best_tau = tau, dwaic_vs_log = dwaic_sum, se)]
BT <- merge(BT, P[form_short == "linear", .(term, dwaic_linear_vs_log = dwaic_sum)], by = "term"); fwrite(BT, file.path(OUT, "sim06_best_tau.csv"))
cat("\nProfile per TSD type (WAIC difference vs log, summed over transitions, mean over imputations; negative = better than log):\n")
print(dcast(P[form_short == "exp"], term ~ tau, value.var = "dwaic_sum", fun.aggregate = function(x) round(x, 1)))
cat("\nBest tau per type:\n"); print(BT)

# ---- pass 2: every TSD term at its type's best tau ----
best <- setNames(BT$best_tau, BT$term); run_pass(2L, best)
p2tag <- paste0("p2_", paste(sprintf("%s%g", sub("tsd_", "", names(best)), best), collapse = "_"))
F2 <- file.path(OUT, "fits", sprintf("%s_m%02d_%s.rds", p2tag, use$m, gsub("->", "_to_", use$label)))
CMP <- rbindlist(Map(function(f1, f2) { a <- readRDS(f1); b <- readRDS(f2); lg <- a$fits$log_all
  rbindlist(lapply(list(linear_all = a$fits$linear_all, best_all = b$fits$best_all), function(x) { ok <- is.na(lg$error) && is.na(x$error); dl <- if (ok) x$local - lg$local else NA_real_
    data.table(dwaic = if (ok) x$waic - lg$waic else NA_real_, se = if (ok) sqrt(length(dl)) * sd(dl) else NA_real_) }), idcol = "form")[, `:=`(label = a$label, m = a$m, n_tau = sum(names(TSD_TERMS) %in% a$terms))] }, F1, F2))
CS <- CMP[, .(dwaic_median = median(dwaic, na.rm = TRUE), dwaic_min = min(dwaic, na.rm = TRUE), dwaic_max = max(dwaic, na.rm = TRUE), se_median = median(se, na.rm = TRUE), n_tau = n_tau[1]), by = .(label, form)]
CS[, verdict := fifelse(dwaic_max < -2 * se_median & dwaic_max < -4, "better than log (all imputations)",
                fifelse(dwaic_min > 2 * se_median & dwaic_min > 4, "worse than log (all imputations)", "no clear difference"))]
CS <- merge(CS, ar, by = "label"); fwrite(CS, file.path(OUT, "sim06_compare.csv"))
cat("\nPer transition vs log (best_all charged +2 per tau when reading):\n"); print(CS[, .N, by = .(form, verdict)])
logl("==== done")
