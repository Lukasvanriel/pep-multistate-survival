# ==============================================================================
# R/tests/wp5_recount_stage1.R   (WP5, user decision 2026-11-04; Mac, full data, NO fitting, NO INLA)
#   Rscript R/tests/wp5_recount_stage1.R [table_dir]        (project root; default table_dir = TABLE_DIR of R/config.R)
# run_stage1() skips recount_events() under SUBSAMPLE_N (a subsample's counts cannot equal transition_set), so the smoke test never runs it.
# This script runs it on the full data: prepare_fit_inputs() without subsampling, then EXACTLY run_stage1's row build (stage1_rows(inp) in 05a,
# the function run_stage1 itself calls: same arguments incl. keep_cols and the C2 exit rule), then recount_events(). Stops on a mismatch.
# Run it BEFORE the smoke test. Needs the output of 04 (data/tables/<RUN_TAG>/).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
args <- commandArgs(TRUE); dir <- if (length(args)) args[1] else TABLE_DIR
cat(sprintf("Tables: %s | MIN_EVENTS %d | XREGEN_EXIT_RULE %s\n", dir, MIN_EVENTS, XREGEN_EXIT_RULE))
inp <- prepare_fit_inputs(dir)                        # no subsample_n
stopifnot("full data expected" = !inp$subsampled)
cat(sprintf("Plots %d | intervals %d | kept transitions %d\n", nrow(inp$p1), nrow(inp$iv), nrow(inp$td)))
s <- stage1_rows(inp)                                 # the Stage 1 build_rows call (entry W_j/2, C2 rule, keep_cols)
rc <- recount_events(s$b$rows, inp$tr)                # stops if a count differs from transition_set or a kept term falls below MIN_EVENTS
cat(sprintf("recount_events(): OK, %d kept transitions, counts equal transition_set, all >= %d events\n", nrow(rc), MIN_EVENTS))
cat("Stage 1 exits by source:", p2_source_line(p2_source_counts(s$p2)), "\n")
xr <- inp$iv[is_change == TRUE & state_end == "regen" & !is.na(D_heavy)]
cat(sprintf("X->regen change intervals with D_heavy: %d (plan: 1,537 = heavy harvest 1,223 + fire 315 - 1 both; printed for information) | rule exits %d | rule exits at the end visit %d | floored event segments %d | floored post-rule stubs %d (a stub starting at its end visit = an exit at the end visit)\n",
            nrow(xr), sum(s$p2$exit_source == "rule"), s$dg$n_rule_at_end, s$b$n_floored_event, s$b$n_floored_stub))
cat("ALL CHECKS PASSED\n")
