# ==============================================================================
# R/simulation/sim06h_class_tau.R   (2026-11-12)   RUN ON THE MACHINE WITH INLA, from the project root:
#   Rscript R/simulation/sim06h_class_tau.R 2>&1 | tee data-output/simulation/sim06h_run.log      (re-running skips finished fits)
# Option B (user 2026-11-12): every TSD term exponential, exp(-TSD / tau), with tau per disturbance type x transition class
#   (into regen / out of regen / between states) as chosen in sim06e (sim06e_best_tau_by_class.csv). Same 72 transitions x 5 imputations,
#   rows, INLA call and WAIC as sim06 pass 2 (setup sourced from sim06_tsd_tau.R; its reproduction check already passed).
#   Fits saved as fits/p3class_mXX_<label>.rds (fit name "best_all"). Afterwards: R/25_I13e_tsd_shape_B.R (no INLA) for the I13 comparison.
# ==============================================================================
src <- readLines(here::here("R", "simulation", "sim06_tsd_tau.R"))
eval(parse(text = src[1:(grep("^# ---- 0. reproduction check", src) - 1)]))         # setup: inputs, use, ar, rebuild_rows, fit_one, retrans
stopifnot(file.exists(file.path(OUT, "check_ok.txt")))
eval(parse(text = src[grep("^fit_forms <- function", src):(grep("^run_pass <- function", src) - 1)]))   # fit_forms (pass 2 branch used)
BTC <- fread(here("data-output", "simulation", "sim06d_regen_screen", "sim06e_best_tau_by_class.csv"))
tau_of <- function(lb) { cl <- ar[label == lb, class]; b <- BTC[class == cl]; setNames(b$tau_best, b$term) }
logl("==== sim06h start | tau by type x class: %s", paste(sprintf("%s/%s=%g", BTC$term, BTC$class, BTC$tau_best), collapse = " "))
todo <- use[!file.exists(file.path(OUT, "fits", sprintf("p3class_m%02d_%s.rds", m, gsub("->", "_to_", label))))]
prog("[sim06h] %d of %d transition x imputation jobs to run", nrow(todo), nrow(use)); t0 <- Sys.time(); done <- 0L
for (mm in sort(unique(todo$m))) { lb_m <- todo[m == mm, label]; rows <- rebuild_rows(mm); sp <- split(rows[label %in% lb_m], by = "label", keep.by = TRUE); rm(rows); invisible(gc())
  res <- mclapply(lb_m, function(lb) { tms <- terms_of(lb); d <- sp[[lb]][, c("status", "Tstop", "truncation", tms), with = FALSE]
    fits <- fit_forms(d, tms, 2L, tau_of(lb))
    saveRDS(list(label = lb, m = mm, terms = tms, n_events = sum(d$status), n_rows = nrow(d), fits = fits),
            file.path(OUT, "fits", sprintf("p3class_m%02d_%s.rds", mm, gsub("->", "_to_", lb))))
    sum(!is.na(vapply(fits, `[[`, "", "error"))) }, mc.cores = N_WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("worker error: ", paste(res[bad], collapse = " | "))
  done <- done + length(lb_m); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  prog("[sim06h] imputation m%02d done (%d jobs, INLA errors %d) | %d/%d | elapsed %.0f min", mm, length(lb_m), sum(unlist(res)), done, nrow(todo), el) }
F1 <- file.path(OUT, "fits", sprintf("p1_m%02d_%s.rds", use$m, gsub("->", "_to_", use$label)))
F3 <- file.path(OUT, "fits", sprintf("p3class_m%02d_%s.rds", use$m, gsub("->", "_to_", use$label)))
W <- rbindlist(Map(function(f1, f3) { a <- readRDS(f1)$fits$log_all; r <- readRDS(f3); b <- r$fits$best_all; ok <- is.na(a$error) && is.na(b$error)
  data.table(label = r$label, m = r$m, dwaic = if (ok) b$waic - a$waic else NA_real_, error = if (!ok) paste(a$error, b$error) else NA_character_) }, F1, F3))
W <- merge(W, ar, by = "label"); fwrite(W, file.path(OUT, "sim06h_waic.csv"))
cat("\nWAIC, every TSD term exponential at its class tau, minus log (summed over transitions, mean over imputations; negative = better than log):\n")
print(W[, .(transitions = uniqueN(label), dwaic = round(sum(dwaic, na.rm = TRUE) / N_IMP, 1), errors = sum(!is.na(error))), by = class])
logl("==== sim06h done")
