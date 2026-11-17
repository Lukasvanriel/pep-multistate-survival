# ==============================================================================
# 25_I11_identifiability.R   (no INLA; <1 min; run from the project root)
# I11: posterior correlations from the pooled total covariance (Rubin T = W + (1 + 1/M) B; pooled_<tag>_final_cov.rds), per transition:
#   (1) disturbance main effect <-> its TSD term; (2) GDD x disturbance <-> TSD term; (3) GDD x disturbance <-> disturbance main effect;
#   (4) alpha <-> TSD terms and alpha <-> GDD (regen exits and all). Strong = |r| > 0.7; very strong = |r| > 0.9.
# ==============================================================================
source(here::here("R", "25_I_common.R"))
cvl <- readRDS(file.path(RESULTS_DIR, sprintf("pooled_%s_final_cov.rds", RUN_TAG)))
ar <- fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[, .(label, class, events)]
DM <- c(fire = "fire_f", pest = "pest_f", harv_light = "harv_light", harv_heavy = "harv_heavy")
pairs <- rbindlist(lapply(names(DM), function(d) data.table(dist = d, p1 = c(DM[[d]], paste0("gdd_x_", d), paste0("gdd_x_", d), "alpha"),
  p2 = paste0("tsd_", d), type = c("main <-> TSD", "GDD x dist <-> TSD", "GDD x dist <-> main", "alpha <-> TSD"))))
pairs[type == "GDD x dist <-> main", p2 := DM[dist]]
pairs <- rbind(pairs, data.table(dist = "-", p1 = "alpha", p2 = "z_gdd_3", type = "alpha <-> GDD"))
R <- rbindlist(lapply(names(cvl), function(lb) { Tm <- cvl[[lb]]$T; nm <- rownames(Tm); C <- cov2cor((Tm + t(Tm)) / 2)
  pp <- pairs[p1 %in% nm & p2 %in% nm]; if (!nrow(pp)) return(NULL); pp[, `:=`(label = lb, r = C[cbind(p1, p2)])] }))
R <- merge(R, ar, by = "label"); fwrite(R, file.path(OUT, "I11_posterior_correlations.csv"))
S <- R[, .(n = .N, median_abs_r = median(abs(r)), n_gt_0.7 = sum(abs(r) > 0.7), n_gt_0.9 = sum(abs(r) > 0.9)), by = .(type, regen_exit = class == "regen_exit")][order(type, regen_exit)]
options(width = 200); print(S)
cat("\n|r| > 0.7:\n"); print(R[abs(r) > 0.7][order(-abs(r)), .(label, class, events, type, dist, r = round(r, 2))])
# (5) the disturbance effect at the median observed TSD (main + f(TSD) x tsd term) vs the main effect alone: posterior SD of each (log HR)
source(here::here("R", "covariate_spec.R")); f_tsd <- tsd_fun(TSD_FORM)
TQ <- c(fire = vi[fire == 1, median(tsd_fire_y, na.rm = TRUE)], pest = vi[pest == 1, median(tsd_pest_y, na.rm = TRUE)],
        harv_light = vi[harv == 1, median(tsd_harv_y, na.rm = TRUE)], harv_heavy = vi[harv == 2, median(tsd_harv_y, na.rm = TRUE)])
E <- rbindlist(lapply(names(cvl), function(lb) { Tm <- cvl[[lb]]$T; nm <- rownames(Tm); rbindlist(lapply(names(DM), function(d) { a <- DM[[d]]; b <- paste0("tsd_", d)
  if (!all(c(a, b) %in% nm)) return(NULL); f <- f_tsd(TQ[[d]])
  data.table(label = lb, dist = d, median_tsd = TQ[[d]], sd_main = sqrt(Tm[a, a]), sd_at_median_tsd = sqrt(Tm[a, a] + f^2 * Tm[b, b] + 2 * f * Tm[a, b])) })) }))
E[, ratio := sd_at_median_tsd / sd_main]; fwrite(E, file.path(OUT, "I11_effect_sd_at_median_tsd.csv"))
cat("\nmedian TSD (y):\n"); print(TQ); cat("\nposterior SD at the median TSD / SD of the main effect, by disturbance:\n")
print(E[, .(n = .N, median_ratio = median(ratio), q10 = quantile(ratio, .1), q90 = quantile(ratio, .9)), by = dist])
