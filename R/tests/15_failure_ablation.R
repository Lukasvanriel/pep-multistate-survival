# ==============================================================================
# 15_failure_ablation.R   RUN ON THE MACHINE WITH INLA (~30-60 min with 6 workers)   -- DIAGNOSTIC ONLY, changes nothing in the pipeline
#
# QUESTION: WHAT in the data/model of the few failing transitions triggers INLA's collapse? Everything else (88 transitions) never fails.
#   Known so far: (1) the datasets where INLA fails are statistically identical to those where it works (independent MLE: same alpha, SE,
#   intercept, truncation); (2) regen->fir fails (14%) while regen->black_spruce, fitted on the SAME exposure rows, never fails; (3) shortening
#   left truncation cured the regen exits; (4) the failing set flips under small data changes.
# DESIGN: one-factor-at-a-time ABLATION on the stored int imputations (rows rebuilt from extras, identical to the fitted data).
#   For each LABEL: up to N_FAIL imputations that failed in the int run + N_OK that worked. Every variant changes ONE thing and is fitted
#   with the pipeline's free INLA call; a variant that removes the failures in the failing imputations (and keeps the working ones working)
#   points to the trigger. Variants that change the likelihood are diagnostics only.
#     as_is            the data as fitted (control: must reproduce the int failures most of the time; INLA is not fully deterministic)
#     as_is_repeat     the identical fit again (determinism check; with INLA_THREADS = "1:1" it should be identical to as_is)
#   INLA_THREADS (default "1:1"): the int run used INLA's default threading inside 6 parallel workers. Running this once with "1:1" and once with
#     NULL shows whether threading changes which fits fail.
#     drop_<term>      one model term removed (all terms in turn)
#     main_only        interactions and TSD terms removed (the main-effects model of the earlier "diag" run)
#     intercept_only   no covariates
#     no_trunc         truncation set to 0                      (known cure in 06q/06n; likelihood changes)
#     trunc_cap4       truncation capped at 4 decades           (likelihood changes)
#     time_years       time in years instead of decades         (EXACT reparametrisation: alpha unchanged, intercept shifts by log 10)
#     time_century     time in centuries                        (exact)
#     drop_hh_rows     rows flagged heavy harvest removed (+ their terms)   (stand-replacing harvest on old stands)
#     drop_fire_rows   rows flagged fire removed (+ their terms)
#     center_<term>    TSD / continuous terms centred at their mean (exact for the betas; only the intercept changes) -> tests scaling of eta
#   FAILED = INLA error/NA, alpha-SD ratio < 0.1, or |alpha_INLA - alpha_MLE| > 3 MLE-SE (MLE of the SAME variant's data: a reference only).
# MUST be run BEFORE 04_prepare_inla.R is re-run (needs the int tables; checked).
# OUTPUT: data-output/inttest/fits/ablation_fits.csv + summary (failure counts per label x variant, split by int-failed / int-ok imputations).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 220)
RUN_TAG <- "int"; IN_DIR <- here("data-output", "inttest", "fits")
N_FAIL <- 15L; N_OK <- 10L; N_WORKERS <- 6L
INLA_THREADS <- "1:1"   # one thread per INLA call: deterministic numerics and no oversubscription (6 workers x all cores). Set to NULL to use INLA's default
if (!is.null(INLA_THREADS)) INLA::inla.setOption(num.threads = INLA_THREADS)
LABELS <- c("black_spruce->regen", "regen->fir", "fir->sugar_maple")
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
iv0 <- readRDS(here("data", "tables", "inttest", "visit_intervals.rds"))
if (!"cmi_sum" %in% names(iv0)) stop("visit_intervals.rds is already the gs version: run this test before re-running 04 (or restore the int tables).")
inp <- prepare_fit_inputs(here("data", "tables", "inttest"))
stopifnot(all(LABELS %in% inp$td$label))

# ---- rows of a stored imputation -----------------------------------------------------------------------------------------------------------------------
key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild_rows <- function(m) {
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  r[, label := paste(from, to, sep = "->")][label %in% LABELS] }

# ---- INLA call (identical settings to model_fitting.R, exact = TRUE) and the independent MLE ----------------------------------------------------------
fit_inla <- function(d, rhs) tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
  control.family = list(variant = 1), control.compute = list(config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) list(error = conditionMessage(e)))
fit_mle <- function(d, tms) {                                   # S(t) = exp(-(lambda t)^alpha), lambda = exp(eta); theta = (log alpha, intercept, betas)
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); a <- exp(o$par[1])
  list(alpha = a, alpha_se = if (is.null(V)) NA_real_ else a * sqrt(V[1, 1]), icpt = o$par[2], conv = o$convergence) }

# ---- variants: function(d, tms) -> list(d = data, tms = terms) ---------------------------------------------------------------------------------------------
tsd_terms <- c("tsd_fire", "tsd_pest", "tsd_harv"); int_terms <- grep("^gdd_x_", inp$terms, value = TRUE)
make_variants <- function(tms) {
  V <- list(as_is = function(d, t) list(d = d, tms = t),
            as_is_repeat = function(d, t) list(d = d, tms = t),   # same fit again: tests whether INLA is deterministic in this setting
            main_only = function(d, t) list(d = d, tms = setdiff(t, c(tsd_terms, int_terms))),
            intercept_only = function(d, t) list(d = d, tms = character(0)),
            no_trunc = function(d, t) { d$truncation <- 0; list(d = d, tms = t) },
            trunc_cap4 = function(d, t) { d$truncation <- pmin(d$truncation, 4); list(d = d, tms = t) },
            time_years = function(d, t) { d$Tstop <- d$Tstop * 10; d$truncation <- d$truncation * 10; list(d = d, tms = t) },
            time_century = function(d, t) { d$Tstop <- d$Tstop / 10; d$truncation <- d$truncation / 10; list(d = d, tms = t) })
  for (tm in tms) V[[paste0("drop_", tm)]] <- local({ tm_ <- tm; function(d, t) list(d = d, tms = setdiff(t, tm_)) })
  for (tm in intersect(tms, c(tsd_terms, "z_gdd_3", "z_cmi_sum", "z_soil_pc1")))
    V[[paste0("center_", tm)]] <- local({ tm_ <- tm; function(d, t) { d[[tm_]] <- d[[tm_]] - mean(d[[tm_]]); list(d = d, tms = t) } })
  if ("harv_heavy" %in% tms) V$drop_hh_rows <- function(d, t) list(d = d[d$harv_heavy == 0, , drop = FALSE], tms = setdiff(t, c("harv_heavy", "gdd_x_harv_heavy")))
  if ("fire_f" %in% tms) V$drop_fire_rows <- function(d, t) list(d = d[d$fire_f == 0, , drop = FALSE], tms = setdiff(t, c("fire_f", "gdd_x_fire", "tsd_fire")))
  V }
fit_variant <- function(d0, lb, m, vname, vfun, tms0) {
  v <- vfun(d0, tms0); d <- v$d; tms <- v$tms
  # drop terms that became constant (e.g. all-zero after dropping flagged rows)
  tms <- tms[vapply(tms, function(tm) length(unique(d[[tm]])) > 1, logical(1))]
  rhs <- if (length(tms)) paste(tms, collapse = " + ") else "1"
  dd <- as.data.frame(d[, c("status", "Tstop", "truncation", tms), drop = FALSE])
  mo <- fit_inla(dd, rhs); ml <- tryCatch(fit_mle(dd, tms), error = function(e) list(alpha = NA_real_, alpha_se = NA_real_, icpt = NA_real_, conv = NA_integer_))
  a <- if (is.null(mo$error)) mo$summary.hyperpar[ALPHA_NAME, ] else NULL
  data.table(m = m, label = lb, variant = vname, rows = nrow(dd), events = sum(dd$status == 1L), n_terms = length(tms), inla_error = !is.null(mo$error),
             alpha = if (is.null(a)) NA_real_ else a$mean, alpha_sd = if (is.null(a)) NA_real_ else a$sd,
             mle_alpha = ml$alpha, mle_alpha_se = ml$alpha_se, mle_conv = ml$conv) }

# ---- imputations ------------------------------------------------------------------------------------------------------------------------------------------
fv <- fread(file.path(IN_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[label %in% LABELS]
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
jobs <- rbindlist(lapply(LABELS, function(lb) { f <- fv[label == lb & m %in% avail]
  rbind(data.table(label = lb, m = head(f[valid == FALSE, sort(m)], N_FAIL), int_failed = TRUE), data.table(label = lb, m = head(f[valid == TRUE, sort(m)], N_OK), int_failed = FALSE)) }))
cat("Jobs (label x imputation):\n"); print(jobs[, .N, by = .(label, int_failed)])

res <- rbindlist(mclapply(sort(unique(jobs$m)), function(mi) tryCatch({
  r <- rebuild_rows(mi)
  rbindlist(lapply(jobs[m == mi, label], function(lb) {
    k <- match(lb, inp$td$label); tms0 <- colnames(inp$mask)[inp$mask[k, ]]; d0 <- as.data.frame(r[label == lb])
    V <- make_variants(tms0); rbindlist(lapply(names(V), function(vn) fit_variant(d0, lb, mi, vn, V[[vn]], tms0))) })) },
  error = function(e) data.table(m = mi, label = NA_character_, variant = paste("run_error:", conditionMessage(e)))), mc.cores = N_WORKERS), fill = TRUE)
res <- merge(res, jobs, by = c("label", "m"), all.x = TRUE)
res[, `:=`(ratio = alpha_sd / (0.78 * alpha / sqrt(events)), z_mle = (alpha - mle_alpha) / mle_alpha_se)]
res[, failed := inla_error | is.na(alpha) | ratio < 0.1 | abs(z_mle) > 3]
fwrite(res, file.path(IN_DIR, sprintf("ablation_fits_threads_%s.csv", if (is.null(INLA_THREADS)) "default" else gsub(":", "-", INLA_THREADS))))

cat("\n== Failures per variant (int-failed imputations | int-ok imputations) ==\n")
s <- res[!is.na(label), .(fail_in_intfailed = sprintf("%d/%d", sum(failed[int_failed]), sum(int_failed)), fail_in_intok = sprintf("%d/%d", sum(failed[!int_failed]), sum(!int_failed)),
                          share_failed = round(mean(failed), 2), mle_alpha = round(median(mle_alpha), 3), events = round(mean(events))), by = .(label, variant)]
for (lb in LABELS) { cat("\n--", lb, "--\n"); print(s[label == lb][order(share_failed)][, -"label"]) }
