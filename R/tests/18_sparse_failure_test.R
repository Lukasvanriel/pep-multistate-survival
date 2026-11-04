# ==============================================================================
# 18_sparse_failure_test.R   RUN ON THE MACHINE WITH INLA, FROM A TERMINAL (~5-10 min with 6 workers)   -- DIAGNOSTIC ONLY
#   cd ~/Desktop/RProjects/pep-multistate-survival && Rscript R/tests/18_sparse_failure_test.R
#
# QUESTION: why do fir->sugar_maple (alpha jumps to ~3.2) and other_temperate->white_birch (INLA program crashes) fail in some imputations?
#   The data show nothing (unimodal likelihood, no extreme rows, failing and healthy datasets indistinguishable), and refits of a failed
#   dataset usually work -> suspected run-to-run numerical noise from INLA's multi-threading.
# DESIGN: every imputation where the transition failed in the "int" run + N_OK healthy ones; rows rebuilt from extras (identical to the fitted data).
#   Each dataset is fitted three times with the pipeline's INLA call:
#     t1_a, t1_b   one thread per fit ("1:1"), twice  -> is the result deterministic, and does the failure persist?
#     default      INLA's default threading (as in the int run), 6 fits in parallel
#   INLA keeps its working files (incl. its own log) for every fit that FAILS, in data-output/inttest/fits/sparse_test/inla_<label>_m<m>_<setting>/,
#   so a crash can be read there. Files of healthy fits are deleted.
# FAILED: INLA error/NA, or alpha-SD ratio < 0.1 (pipeline screen). The MLE is not used here.
# SAVES AS IT GOES: one file per fit; re-running skips fits already saved.
# MUST be run BEFORE 04_prepare_inla.R is re-run (needs the int tables; checked).
# OUTPUT: data-output/inttest/fits/sparse_failure_test.csv + printed summary (+ the tail of INLA's log for every failed fit).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 200)
RUN_TAG <- "int"; IN_DIR <- here("data-output", "inttest", "fits"); N_WORKERS <- 6L; N_OK <- 5L
LABELS <- c("fir->sugar_maple", "other_temperate->white_birch")
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (!"cmi_sum" %in% names(readRDS(here("data", "tables", "inttest", "visit_intervals.rds")))) stop("visit_intervals.rds is already the gs version: run this before re-running 04.")
inp <- prepare_fit_inputs(here("data", "tables", "inttest")); stopifnot(all(LABELS %in% inp$td$label))
OUT <- file.path(IN_DIR, "sparse_test"); dir.create(OUT, showWarnings = FALSE)
cat("INLA version:", as.character(packageVersion("INLA")), "\n")

key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild_rows <- function(m) {
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  r[, label := paste(from, to, sep = "->")][label %in% LABELS] }

fit_one <- function(dd, rhs, threads, wdir) {
  args <- list(formula = as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = dd,
               control.family = list(variant = 1), control.compute = list(config = TRUE),
               control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
               working.directory = wdir, keep = TRUE)
  if (!is.null(threads)) args$num.threads <- threads
  tryCatch(do.call(inla, args), error = function(e) list(error = conditionMessage(e))) }

# ---- jobs ------------------------------------------------------------------------------------------------------------------------------------------------
fv <- fread(file.path(IN_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[label %in% LABELS]
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
ds <- rbindlist(lapply(LABELS, function(lb) { f <- fv[label == lb & m %in% avail]
  rbind(data.table(label = lb, m = f[valid == FALSE, sort(m)], int_failed = TRUE), data.table(label = lb, m = head(f[valid == TRUE, sort(m)], N_OK), int_failed = FALSE)) }))
jobs <- ds[, .(setting = c("t1_a", "t1_b", "default")), by = .(label, m, int_failed)]
jobs[, tag := sprintf("%s_m%03d_%s", gsub("->", "_to_", label), m, setting)]
jobs[, file := file.path(OUT, paste0("fit_", tag, ".csv"))]
cat("Datasets:\n"); print(ds[, .N, by = .(label, int_failed)])
cat(sprintf("Fits: %d, already saved: %d\n", nrow(jobs), sum(file.exists(jobs$file))))

run_jobs <- function(J) invisible(mclapply(seq_len(nrow(J)), function(i) { j <- J[i]; ts <- proc.time()["elapsed"]
  res <- tryCatch({
    d <- rebuild_rows(j$m)[label == j$label]; k <- match(j$label, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]
    dd <- as.data.frame(d[, c("status", "Tstop", "truncation", tms), with = FALSE])
    wdir <- file.path(OUT, paste0("inla_", j$tag)); unlink(wdir, recursive = TRUE)
    mo <- fit_one(dd, inp$rhs[[k]], if (j$setting == "default") NULL else "1:1", wdir)
    a <- if (is.null(mo$error)) mo$summary.hyperpar[ALPHA_NAME, ] else NULL
    out <- data.table(label = j$label, m = j$m, int_failed = j$int_failed, setting = j$setting, events = sum(dd$status == 1L),
                      inla_error = !is.null(mo$error), alpha = if (is.null(a)) NA_real_ else a$mean, alpha_sd = if (is.null(a)) NA_real_ else a$sd,
                      icpt = if (is.null(mo$error)) mo$summary.fixed["(Intercept)", "mean"] else NA_real_,
                      mode_status = if (is.null(mo$error)) as.integer(mo$mode$mode.status) else NA_integer_)
    out[, ratio := alpha_sd / (0.78 * alpha / sqrt(events))][, failed := inla_error | is.na(alpha) | ratio < 0.1]
    if (!out$failed) unlink(wdir, recursive = TRUE)        # keep INLA's files (incl. its log) only for failed fits
    out }, error = function(e) data.table(label = j$label, m = j$m, int_failed = j$int_failed, setting = j$setting, inla_error = TRUE, failed = TRUE,
                                          run_error = substr(conditionMessage(e), 1, 200)))
  res[, seconds := round(proc.time()["elapsed"] - ts)]; fwrite(res, j$file); NULL }, mc.cores = N_WORKERS, mc.preschedule = FALSE))

t0 <- proc.time()["elapsed"]
run_jobs(jobs[setting != "default" & !file.exists(file)])    # one thread per fit
run_jobs(jobs[setting == "default" & !file.exists(file)])    # INLA default threading, 6 in parallel (the int-run context)
cat(sprintf("Run time: %.1f min\n", (proc.time()["elapsed"] - t0) / 60))

# ---- summary ---------------------------------------------------------------------------------------------------------------------------------------------
res <- rbindlist(lapply(jobs[file.exists(file), file], fread), fill = TRUE)
fwrite(res, file.path(IN_DIR, "sparse_failure_test.csv"))
cat("\n== Failures per setting: datasets that failed in the int run | that were healthy ==\n")
print(res[, .(fail_in_intfailed = sprintf("%d/%d", sum(failed[int_failed]), sum(int_failed)), fail_in_intok = sprintf("%d/%d", sum(failed[!int_failed]), sum(!int_failed)),
              alpha_failed = round(median(alpha[failed], na.rm = TRUE), 2), alpha_ok = round(median(alpha[!failed], na.rm = TRUE), 2)), by = .(label, setting)])
w <- dcast(res, label + m + int_failed ~ setting, value.var = c("failed", "alpha"))
cat("\n== Determinism with one thread: is t1_a identical to t1_b? ==\n")
print(w[, .(datasets = .N, same_outcome = sum(failed_t1_a == failed_t1_b), alpha_identical = sum(abs(alpha_t1_a - alpha_t1_b) < 1e-8, na.rm = TRUE)), by = label])
cat("\n== Per dataset (TRUE = failed) ==\n")
print(w[, .(label, m, int_failed, t1_a = failed_t1_a, t1_b = failed_t1_b, default = failed_default,
            alpha_t1 = round(alpha_t1_a, 3), alpha_default = round(alpha_default, 3))][order(label, -int_failed, m)])

cat("\n== Tail of INLA's log for failed fits (first 4 failed fits per label) ==\n")
for (lb in LABELS) for (tg in head(jobs[res[label == lb & failed == TRUE, .(m, setting)], on = .(m, setting)][label == lb, tag], 4)) {
  lf <- list.files(file.path(OUT, paste0("inla_", tg)), pattern = "log", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  cat("\n---", tg, "---", if (!length(lf)) "(no log file found)" else "", "\n")
  for (f in lf) { x <- readLines(f, warn = FALSE); cat(tail(x, 15), sep = "\n") } }
