# ==============================================================================
# 19_alpha_prior_test.R   RUN ON THE MACHINE WITH INLA, FROM A TERMINAL (~10-15 min with 6 workers)   -- DIAGNOSTIC ONLY
#   cd ~/Desktop/RProjects/pep-multistate-survival && Rscript R/tests/19_alpha_prior_test.R
#
# QUESTION: how much does INLA's default prior on the Weibull shape (pc.alphaw, strength 5, pulls alpha toward 1) drive the results,
#   and does a weaker prior remove the other_temperate->white_birch crashes (flat alpha posterior -> integration grid reaches alpha ~0.13 -> NaN)?
#   Background: in the simulation (p1p2_v2_test/stepB, Steps 8-9) weakening the PC prior from 5 to 0.01 roughly halved the alpha bias of
#   high-alpha transitions, on exact and on imputed data. The earlier "wide prior" test there (loggamma(0.1,0.1)) was invalid because alpha is
#   estimated on an internal scale (alpha = exp(0.1 * theta)). Weakening the prior was never tried on the final fits or on real data.
# PRIORS: default (no hyper given) | pc5 (explicit, same as default: sanity check) | pc1 | pc0.1.
#   control.family = list(variant = 1, hyper = list(alpha = list(prior = "pc.alphaw", param = <strength>))). Lower = weaker pull toward 1.
# DATASETS (stored "int" imputations, rows rebuilt from extras = identical to the fitted data):
#   A. other_temperate->white_birch: the 6 datasets that crash reproducibly (script 18) + 5 healthy ones
#   B. up to N_SPARSE sparse transitions (<= 40 events) whose int estimate of alpha is close to 1 (0.95-1.10): 3 healthy imputations each
#   C. two large transitions (fir->white_birch, black_spruce->regen): 3 healthy imputations each
#   Every dataset is also fitted by maximum likelihood (no prior) as a reference.
# FAILED: INLA error/NA, or alpha-SD ratio < 0.1. SAVES AS IT GOES (one file per dataset); re-running skips finished datasets.
# MUST be run BEFORE 04_prepare_inla.R is re-run (needs the int tables; checked).
# OUTPUT: data/processed/fit_int/prior_test/ (per dataset) + data/processed/fit_int/alpha_prior_test.csv + printed summary.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 220)
RUN_TAG <- "int"; IN_DIR <- here("data", "processed", "fit_int"); N_WORKERS <- 6L; INLA::inla.setOption(num.threads = "1:1")
PRIORS <- list(default = NULL, pc5 = 5, pc1 = 1, pc0.1 = 0.1)
N_SPARSE <- 8L; LARGE <- c("fir->white_birch", "black_spruce->regen"); OTWB <- "other_temperate->white_birch"
OTWB_CRASH <- c(41L, 64L, 78L, 87L, 88L, 111L); OTWB_OK <- 1:5; SANITY_N <- 3L        # pc5 is fitted only on the first SANITY_N datasets
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (!"cmi_sum" %in% names(readRDS(here("data", "processed", "visit_intervals.rds")))) stop("visit_intervals.rds is already the gs version: run this before re-running 04.")
inp <- prepare_fit_inputs(here("data", "processed"))
OUT <- file.path(IN_DIR, "prior_test"); dir.create(OUT, showWarnings = FALSE)
cat("INLA version:", as.character(packageVersion("INLA")), "\n")

key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild_rows <- function(m, lb) {
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  r[paste(from, to, sep = "->") == lb] }

fit_inla <- function(dd, rhs, strength) {
  cf <- list(variant = 1); if (!is.null(strength)) cf$hyper <- list(alpha = list(prior = "pc.alphaw", param = strength))
  tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = dd,
    control.family = cf, control.compute = list(config = TRUE),
    control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) list(error = conditionMessage(e))) }
fit_mle <- function(d, tms) {                                   # S(t) = exp(-(lambda t)^alpha), lambda = exp(eta); theta = (log alpha, intercept, betas)
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(diag(V))
  data.table(parameter = c("alpha", "(Intercept)", tms), est = c(exp(o$par[1]), o$par[-1]), sd = c(exp(o$par[1]) * se[1], se[-1])) }

# ---- datasets ----------------------------------------------------------------------------------------------------------------------------------------------
fv <- fread(file.path(IN_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))
ev <- fread(file.path(IN_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))[parameter == "alpha", .(n_events = median(n_events)), by = label]
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
first_ok <- function(lb, n = 3L) head(fv[label == lb & valid == TRUE & m %in% avail, sort(m)], n)
sp <- merge(fv[valid == TRUE, .(alpha_int = median(alpha)), by = label], ev, by = "label")[n_events <= 40 & alpha_int >= 0.95 & alpha_int <= 1.10 & label != OTWB][order(n_events)]
SPARSE <- head(sp$label, N_SPARSE)
ds <- rbind(data.table(group = "A_otwb", label = OTWB, m = c(OTWB_CRASH, OTWB_OK), crash18 = c(rep(TRUE, length(OTWB_CRASH)), rep(FALSE, length(OTWB_OK)))),
            rbindlist(lapply(SPARSE, function(lb) data.table(group = "B_sparse", label = lb, m = first_ok(lb), crash18 = FALSE))),
            rbindlist(lapply(LARGE, function(lb) data.table(group = "C_large", label = lb, m = first_ok(lb), crash18 = FALSE))))
ds <- ds[m %in% avail]; ds[, id := .I]; ds[, file := file.path(OUT, sprintf("ds_%s_m%03d.rds", gsub("->", "_to_", label), m))]
cat("Sparse transitions selected:\n"); print(sp[label %in% SPARSE]); cat("\nDatasets:\n"); print(ds[, .N, by = .(group, label)])
todo <- ds[!file.exists(file)]; cat(sprintf("Datasets: %d, already done: %d\n", nrow(ds), nrow(ds) - nrow(todo)))

t0 <- proc.time()["elapsed"]
invisible(mclapply(seq_len(nrow(todo)), function(i) { j <- todo[i]
  res <- tryCatch({
    d <- rebuild_rows(j$m, j$label); k <- match(j$label, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
    dd <- as.data.frame(d[, c("status", "Tstop", "truncation", tms), with = FALSE]); nev <- sum(dd$status == 1L)
    pri <- if (j$id <= SANITY_N) names(PRIORS) else setdiff(names(PRIORS), "pc5")
    fits <- rbindlist(lapply(pri, function(p) { mo <- fit_inla(dd, inp$rhs[[k]], PRIORS[[p]])
      if (!is.null(mo$error)) return(data.table(prior = p, parameter = "alpha", est = NA_real_, sd = NA_real_, inla_error = TRUE))
      h <- mo$summary.hyperpar; f <- mo$summary.fixed
      data.table(prior = p, parameter = c("alpha", rownames(f)), est = c(h[ALPHA_NAME, "mean"], f$mean), sd = c(h[ALPHA_NAME, "sd"], f$sd), inla_error = FALSE) }))
    ml <- fit_mle(dd, tms)[, `:=`(prior = "MLE", inla_error = FALSE)]
    cbind(j[, .(group, label, m, crash18)], events = nev, rbind(fits, ml, fill = TRUE)) },
    error = function(e) cbind(j[, .(group, label, m, crash18)], data.table(prior = "run_error", parameter = substr(conditionMessage(e), 1, 150))))
  saveRDS(res, j$file); NULL }, mc.cores = N_WORKERS, mc.preschedule = FALSE))
cat(sprintf("Run time: %.1f min\n", (proc.time()["elapsed"] - t0) / 60))

# ---- summary --------------------------------------------------------------------------------------------------------------------------------------------
res <- rbindlist(lapply(ds[file.exists(file), file], readRDS), fill = TRUE)
fwrite(res, file.path(IN_DIR, "alpha_prior_test.csv"))
a <- res[parameter == "alpha"]
a[, failed := prior != "MLE" & (inla_error | is.na(est) | sd / (0.78 * est / sqrt(events)) < 0.1)]
w <- dcast(a, group + label + m + crash18 + events ~ prior, value.var = c("est", "sd", "failed"))

cat("\n== 1. Sanity check: explicit pc.alphaw(5) must equal the default ==\n")
print(w[!is.na(est_pc5), .(label, m, alpha_default = est_default, alpha_pc5 = est_pc5, diff = signif(est_pc5 - est_default, 3))])

cat("\n== 2. other_temperate->white_birch: failed fits per prior (crash datasets of script 18 | healthy) ==\n")
print(a[group == "A_otwb" & prior != "MLE", .(crash_datasets = sprintf("%d/%d", sum(failed[crash18]), sum(crash18)),
                                               healthy_datasets = sprintf("%d/%d", sum(failed[!crash18]), sum(!crash18))), by = prior])

cat("\n== 3. alpha per transition and prior (median over datasets; MLE = no prior; healthy fits only) ==\n")
s3 <- a[failed == FALSE, .(alpha = round(median(est), 3), sd = round(median(sd), 3)), by = .(group, label, events, prior)]
print(dcast(s3, group + label + events ~ prior, value.var = "alpha")[order(group, events)])
cat("\n   posterior SD of alpha (MLE column = MLE standard error):\n")
print(dcast(s3, group + label + events ~ prior, value.var = "sd")[order(group, events)])

cat("\n== 4. Covariate effects: change in log hazard ratio (alpha * beta) vs the default prior, per term ==\n")
b <- res[!parameter %in% c("alpha", "(Intercept)") & prior %in% names(PRIORS) & !is.na(est)]
b <- merge(b, a[, .(label, m, prior, alpha = est, failed)], by = c("label", "m", "prior"))[failed == FALSE]
b[, `:=`(logHR = alpha * est, logHR_sd = alpha * sd)]
bd <- merge(b[prior != "default"], b[prior == "default", .(label, m, parameter, logHR0 = logHR, sd0 = logHR_sd)], by = c("label", "m", "parameter"))
bd[, shift_sd := (logHR - logHR0) / sd0]
print(bd[, .(terms = .N, median_abs_shift_sd = round(median(abs(shift_sd)), 2), max_abs_shift_sd = round(max(abs(shift_sd)), 2)), by = .(group, label, prior)][order(group, label, prior)])
cat("\n   (shift in units of the default fit's posterior SD of the log hazard ratio, approximated as alpha * SD(beta))\n")
