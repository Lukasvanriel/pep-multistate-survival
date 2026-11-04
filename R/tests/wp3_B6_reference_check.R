# ==============================================================================
# R/tests/wp3_B6_reference_check.R   (WP3, B6 check; run in a SCRATCH tree, NOT from the project root)
# Question: do fit_map_v0() / H_at() as moved into 05a reproduce data-output/diagnostics/cap_test_new_rows_m01_nocut/results.csv (arms A and B)?
# The rows of tests/22 (CAP_CUT=FALSE, main run) need the BASELINE code (commit 2f4382a: config.R with RUN_TAG "main", covariate_spec.R, 04b_build_rows.R,
# 05a_fit_funs.R) and the MAIN tables, because WP2 changed the table columns. Layout used (2026-11-04), all outside the project:
#   <scratch>/ref/  .here, R/{config,covariate_spec,04b_build_rows,05a_fit_funs}.R from `git show 2f4382a:R/<file>`, R/tests/22_cap_test_new_rows.R (unchanged),
#                   R/tests/<this file>, data/tables/main/*, data/build/combined_covariates.rds, data-output/main/fits/extras_main_m01.rds,
#                   data-output/diagnostics/ref_nocut/{results.csv, fits/*.rds}  (copy of cap_test_new_rows_m01_nocut)
#   <scratch>/wp3/R/{config,05a_fit_funs}.R  = the NEW files. Run from <scratch>/ref:  Rscript R/tests/wp3_B6_reference_check.R
# It evaluates tests/22 up to the fitting section (rows only), sources the new config + 05a, and compares maxH_map, alpha_map, n_above_*, maxH_inla
# (H_at at the stored INLA estimates). INLA itself is not rerun.
# ==============================================================================
# WP3 B6 check (scratch tool, not part of the project): rows of tests/22 (CAP_CUT=FALSE, arms A and B) built with the BASELINE code (2f4382a),
# then the MOVED fit_map_v0() / H_at() of the new 05a applied to those rows and compared with cap_test_new_rows_m01_nocut/results.csv and fits/*.rds.
Sys.setenv(CAP_CUT = "FALSE", CAP_ENGINE = "map", CAP_M = "1", CAP_WORKERS = "2")
L <- readLines("R/tests/22_cap_test_new_rows.R")
i_fit <- grep("^# ---- fitting", L); i_cols <- grep('^cols <- c\\("status"', L); i_rm <- grep("^rm\\(A, B\\)", L)
stopifnot(length(i_fit) == 1, length(i_cols) == 1, length(i_rm) == 1)
eval(parse(text = L[1:(i_fit - 1)])); eval(parse(text = L[i_cols:i_rm]))      # builds td, terms_of, DATA (exactly as in 22)
# ---- NEW functions: config + 05a of the work-package copy ----
source("../wp3/R/config.R"); source("../wp3/R/05a_fit_funs.R")
stopifnot(exists("fit_map_v0"), COEF_PRIOR_SD == 5)
ref <- fread("data-output/diagnostics/ref_nocut/results.csv"); fdir <- "data-output/diagnostics/ref_nocut/fits"
jobs <- CJ(arm = c("A", "B"), i = seq_len(nrow(td)))
out <- rbindlist(parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  arm_j <- jobs$arm[j]; i <- jobs$i[j]; lb <- td$label[i]; tms <- terms_of(i); d <- DATA[[arm_j]][[lb]]
  st <- readRDS(file.path(fdir, sprintf("%s_%02d.rds", arm_j, i))); stopifnot(st$label == lb)
  fm <- fit_map_v0(d, tms); hm <- H_at(d, tms, fm$tab)
  hi <- H_at(d, tms, st$inla[, .(parameter, est, sd)])
  r <- ref[which(ref$arm == arm_j & ref$label == lb)]; stopifnot(nrow(r) == 1L)
  cmp <- merge(fm$tab, st$map, by = "parameter", suffixes = c("", "_ref"))
  data.table(arm = arm_j, label = lb,
    alpha_map = fm$tab[parameter == "alpha", est], alpha_map_ref = r$alpha_map,
    maxH_map = hm[["maxH"]], maxH_map_ref = r$maxH_map, n_above_map = hm[["n_above"]], n_above_map_ref = r$n_above_map,
    maxH_inla = hi[["maxH"]], maxH_inla_ref = r$maxH_inla, n_above_inla = hi[["n_above"]], n_above_inla_ref = r$n_above_inla,
    n_above_ev_inla = hi[["n_above_ev"]], n_above_ev_inla_ref = r$n_above_ev_inla,
    coef_maxrel = max(abs(cmp$est - cmp$est_ref) / pmax(abs(cmp$est_ref), 1e-8), na.rm = TRUE),
    se_maxrel = max(abs(cmp$sd - cmp$sd_ref) / pmax(abs(cmp$sd_ref), 1e-8), na.rm = TRUE), conv = fm$conv, conv_ref = st$map_conv) }), fill = TRUE)
rel <- function(a, b) ifelse(is.na(a) & is.na(b), 0, ifelse(is.na(a) | is.na(b), Inf, abs(a - b) / pmax(abs(b), 1e-12)))
out[, `:=`(d_alpha = rel(alpha_map, alpha_map_ref), d_maxH_map = rel(maxH_map, maxH_map_ref), d_maxH_inla = rel(maxH_inla, maxH_inla_ref))]
cnt <- out[, .(n_fits = .N, n_above_map_mismatch = sum(n_above_map != n_above_map_ref, na.rm = TRUE), n_above_inla_mismatch = sum(n_above_inla != n_above_inla_ref, na.rm = TRUE),
               n_above_ev_inla_mismatch = sum(n_above_ev_inla != n_above_ev_inla_ref, na.rm = TRUE), NA_mismatch = sum(is.na(maxH_map) != is.na(maxH_map_ref)),
               max_rel_alpha_map = max(d_alpha), max_rel_maxH_map = max(d_maxH_map), max_rel_maxH_inla = max(d_maxH_inla),
               max_rel_coef = max(coef_maxrel[is.finite(coef_maxrel)]), max_rel_se = max(se_maxrel[is.finite(se_maxrel)]), conv_mismatch = sum(conv != conv_ref)), by = arm]
print(cnt, digits = 4)
fwrite(out, "data-output/diagnostics/ref_nocut/b6_reproduction.csv")
cat("\nWorst 5 fits by relative difference in maxH_map / alpha_map:\n"); print(out[order(-pmax(d_maxH_map, d_alpha))][1:5, .(arm, label, maxH_map, maxH_map_ref, d_maxH_map, alpha_map, alpha_map_ref, d_alpha)], digits = 6)
