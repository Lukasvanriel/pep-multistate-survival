# ==============================================================================
# R/tests/wp5_C23_real_check.R   (WP5: real-data checks of C2 and C3; read-only; no INLA)
#   Rscript R/tests/wp5_C23_real_check.R [M = 5] [fit_dir] [table_dir]       (project root; defaults FIT_DIR / TABLE_DIR of R/config.R)
# Needs the first imputations of 05 (extras_<tag>_m01..mM.rds) and stage2.rds in fit_dir, and visit_intervals.rds in table_dir.
# PASS / FAIL (plan C2 and C3 checks):
#   C2-1  rule exits per imputation == X->regen change intervals with D_heavy (plan: 1,537 on the full table; printed)
#   C2-2  every rule exit == t_start + D_heavy - year_start (tolerance 1e-9); the rule events are exactly those intervals (none missing, none extra)
#   C3-1  among non-rule events present in all M imputations: number with sd(exit_abs) == 0 (pass = 0)
# INFORMATION (no threshold): exit_source counts per label and imputation, per transition (imputation 1), uniform_nonfinite per transition, floored
#   segments and rule exits at the end visit, mean of the draws vs the Stage 2 moment-matched point.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R"))
args <- commandArgs(TRUE)
M <- if (length(args) >= 1) as.integer(args[1]) else 5L
fit_dir <- if (length(args) >= 2) args[2] else FIT_DIR; table_dir <- if (length(args) >= 3) args[3] else TABLE_DIR
iv <- readRDS(file.path(table_dir, "visit_intervals.rds")); setDT(iv)
p2 <- rbindlist(lapply(seq_len(M), function(m) { e <- readRDS(file.path(fit_dir, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p <- as.data.table(e$p2); p[, m := m]; p }))
stopifnot("extras$p2 has no exit_source (run 05 with the WP5 code)" = "exit_source" %in% names(p2))
p2[, transition := paste(from, to, sep = "->")]
key <- function(d) paste(d$id_pe, d$interval)
xr <- iv[is_change == TRUE & state_end == "regen" & !is.na(D_heavy)]
bad <- character(0)
cat(sprintf("Imputations read: %d | change intervals %d | X->regen intervals with D_heavy %d (plan: 1,537)\n", M, sum(iv$is_change), nrow(xr)))

cat("\nExit source by imputation:\n"); print(dcast(p2[, .N, by = .(m, exit_source)], m ~ exit_source, value.var = "N", fill = 0))
cat("\nFloored segments (exit at the last visit) and rule exits at the end visit, by imputation:\n")
iv_k <- iv[, .(id_pe, interval, t_start, t_end, year_start, D_heavy)]; p2 <- merge(p2, iv_k, by = c("id_pe", "interval"), all.x = TRUE)
print(p2[, .(floored = sum(floored), rule_at_end = sum(exit_source == "rule" & abs(exit_abs - t_end) < 1e-9)), by = m])

# C2-1, C2-2
for (mm in seq_len(M)) {
  r <- p2[m == mm & exit_source == "rule"]
  if (nrow(r) != nrow(xr)) bad <- c(bad, sprintf("C2-1 m=%d: %d rule exits vs %d X->regen intervals with D_heavy", mm, nrow(r), nrow(xr)))
  if (!setequal(key(r), key(xr))) bad <- c(bad, sprintf("C2-2 m=%d: the rule events are not the X->regen intervals with D_heavy (%d missing, %d extra)", mm, length(setdiff(key(xr), key(r))), length(setdiff(key(r), key(xr)))))
  dev <- max(abs(r$exit_abs - (r$t_start + r$D_heavy - r$year_start)))
  if (!is.finite(dev) || dev >= 1e-9) bad <- c(bad, sprintf("C2-2 m=%d: max |exit - (t_start + D_heavy - year_start)| = %.3g", mm, dev))
  if (any(!is.na(p2[m == mm & exit_source != "rule", D_heavy]) & p2[m == mm & exit_source != "rule", to] == "regen")) bad <- c(bad, sprintf("C2-2 m=%d: an X->regen event with D_heavy is not a rule exit", mm))
}
cat(sprintf("\nC2: rule exits per imputation %s | max |exit - rule date| over all imputations %.2e\n", paste(p2[exit_source == "rule", .N, by = m][order(m)]$N, collapse = " / "),
            p2[exit_source == "rule", max(abs(exit_abs - (t_start + D_heavy - year_start)))]))

# C3-1 spread
nr <- p2[exit_source != "rule"]; nr[, k := key(nr)]
sp <- nr[, .(n_m = .N, sd_exit = sd(exit_abs), sources = paste(sort(unique(exit_source)), collapse = "+")), by = .(k, transition)][n_m == M]
n0 <- sp[sd_exit == 0, .N]
cat(sprintf("C3: non-rule events present in all %d imputations: %d | with sd(exit) == 0: %d | median sd %.3f y\n", M, nrow(sp), n0, median(sp$sd_exit)))
if (n0 > 0) { cat("  events with zero spread (first 20):\n"); print(sp[sd_exit == 0][1:min(20, n0)]); bad <- c(bad, sprintf("C3-1: %d non-rule events with sd(exit) == 0", n0)) }

cat("\nExit source by transition, imputation 1 (counts):\n"); print(dcast(p2[m == 1, .N, by = .(transition, exit_source)], transition ~ exit_source, value.var = "N", fill = 0), nrows = 200)
un <- p2[exit_source == "uniform_nonfinite", .N, by = .(transition, m)]
cat("\nuniform_nonfinite by transition and imputation:\n"); if (nrow(un)) print(dcast(un, transition ~ m, value.var = "N", fill = 0)) else cat("  none\n")

# information: mean of the draws vs the Stage 2 moment-matched point
s2 <- tryCatch(as.data.table(readRDS(file.path(fit_dir, "stage2.rds"))$p2), error = function(e) NULL)
if (is.null(s2) || !"exit_source" %in% names(s2)) cat("\n(no Stage 2 event table in stage2.rds: draw-vs-moment comparison skipped)\n") else {
  dr <- p2[exit_source == "draw", .(n_m = .N, mean_draw = mean(exit_abs), width = t_end[1] - t_start[1]), by = .(id_pe, interval)][n_m == M]
  cmp <- merge(dr, s2[exit_source == "moment", .(id_pe, interval, moment = exit_abs)], by = c("id_pe", "interval"))
  cmp[, `:=`(diff = mean_draw - moment, rel = (mean_draw - moment) / width)]
  cat(sprintf("\nMean of the %d draws minus the Stage 2 moment point (information only; %d events with a draw in all imputations and a moment point):\n", M, nrow(cmp)))
  print(rbind(diff_years = quantile(cmp$diff, c(0.05, 0.25, 0.5, 0.75, 0.95)), diff_per_bracket_width = quantile(cmp$rel, c(0.05, 0.25, 0.5, 0.75, 0.95))), digits = 3)
  cat(sprintf("  mean diff %.3f y | mean |diff| %.3f y | mean diff / bracket width %.3f\n", mean(cmp$diff), mean(abs(cmp$diff)), mean(cmp$rel)))
}
if (length(bad)) { cat("\n!! FAILED:\n"); cat(paste0("  - ", bad, "\n"), sep = ""); stop("C2/C3 real-data check failed; report this output to the user, change nothing.", call. = FALSE) }
cat("\nALL CHECKS PASSED\n")
