# ==============================================================================
# R/tests/04t_compare_with_inttest.R  (T9, 2026-10-26) -- compare the 04 output of this run (TABLE_DIR, RUN_TAG in R/config.R) with
# the inttest tables (INTTEST_TABLE_DIR): plots, intervals, state changes, kept transitions and covariate terms per transition.
# Run from the project root after 04:  Rscript R/tests/04t_compare_with_inttest.R
# Output: TABLE_DIR/compare_with_inttest_overview.csv and TABLE_DIR/compare_with_inttest_transitions.csv (only transitions that differ).
# Harvest-term reasons: "none" because < MIN_EVENTS (15) state changes WITHOUT harvest (no reference group) or < 15 WITH harvest.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); options(width = 220)
MIN_EV <- 15L
ld <- function(d) list(tr = as.data.table(readRDS(file.path(d, "transition_set.rds"))), iv = as.data.table(readRDS(file.path(d, "visit_intervals.rds"))),
                       pi = as.data.table(readRDS(file.path(d, "plot_info.rds"))))
A <- ld(INTTEST_TABLE_DIR); B <- ld(TABLE_DIR)
ov <- function(x, tag) { k <- x$tr[kept == TRUE]
  data.table(run = tag, plots = nrow(x$pi), intervals = nrow(x$iv), state_changes = sum(x$iv$is_change), transitions_kept = nrow(k),
             changes_in_kept = sum(k$n_events), fire_terms = sum(k$fire_term), pest_terms = sum(k$pest_term),
             harv_3level = sum(k$harv_term == "3level"), harv_any = sum(k$harv_term == "any"), harv_none = sum(k$harv_term == "none"),
             flagged_fire = sum(k$n_fire), flagged_pest = sum(k$n_pest), flagged_harv_light = sum(k$n_harv1), flagged_harv_heavy = sum(k$n_harv2)) }
overview <- rbind(ov(A, "inttest"), ov(B, RUN_TAG))
cat("Overview (flagged_* = state changes with that flag, in kept transitions):\n"); print(overview, row.names = FALSE)
cat(sprintf("\nPlots only in inttest: %d | only in %s: %d\n", sum(!A$pi$id_pe %in% B$pi$id_pe), RUN_TAG, sum(!B$pi$id_pe %in% A$pi$id_pe)))

m <- merge(A$tr, B$tr, by = c("from", "to"), suffixes = c("_a", "_b"))
h <- function(t) fcase(t == "3level", "3-level", t == "any", "any", default = "-")
d <- m[(kept_a | kept_b) & (kept_a != kept_b | fire_term_a != fire_term_b | pest_term_a != pest_term_b | harv_term_a != harv_term_b)]
out <- d[order(from, to), .(transition = paste(from, "->", to), state_changes = paste(n_events_a, "->", n_events_b),
  kept = fifelse(kept_a != kept_b, paste(kept_a, "->", kept_b), ""),
  fire = fifelse(fire_term_a != fire_term_b, sprintf("%s -> %s (%d -> %d)", fifelse(fire_term_a, "yes", "-"), fifelse(fire_term_b, "yes", "-"), n_fire_a, n_fire_b), ""),
  pest = fifelse(pest_term_a != pest_term_b, sprintf("%s -> %s (%d -> %d)", fifelse(pest_term_a, "yes", "-"), fifelse(pest_term_b, "yes", "-"), n_pest_a, n_pest_b), ""),
  harvest = fifelse(harv_term_a != harv_term_b, sprintf("%s -> %s (none/light/heavy %d/%d/%d -> %d/%d/%d)", h(harv_term_a), h(harv_term_b),
                    n_harv0_a, n_harv1_a, n_harv2_a, n_harv0_b, n_harv1_b, n_harv2_b), ""),
  harvest_reason = fcase(harv_term_a == harv_term_b | !kept_b, "", harv_term_b == "none" & n_harv0_b < MIN_EV, "< 15 without harvest",
                         harv_term_b == "none", "< 15 with harvest", default = "enough with and without harvest"))]
cat(sprintf("\nKept transitions whose status or terms differ (inttest -> %s): %d\n", RUN_TAG, nrow(out))); print(out, row.names = FALSE)
fwrite(overview, file.path(TABLE_DIR, "compare_with_inttest_overview.csv"))
fwrite(out, file.path(TABLE_DIR, "compare_with_inttest_transitions.csv"))
cat("Saved:", file.path(TABLE_DIR, "compare_with_inttest_overview.csv"), "and compare_with_inttest_transitions.csv\n")
