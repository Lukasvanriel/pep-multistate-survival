# ==============================================================================
# 25_I13d_tsd_shape_summary.R   (no INLA; after 25_I13c_tsd_shape.R; 2026-11-12)
# Group shares per GDD bin: observed; main I13 (pooled coefficients, log); log and exp from the same 5 sim06 imputations; undisturbed (none) and
#   disturbance regime (SC). Shift = regime minus undisturbed, per shape.
# ==============================================================================
source(here::here("R", "25_I_common.R"))
O <- fread(file.path(OUT, "I13_group_shares.csv"))[, .(bin, zg, grp, obs, main_none = I8, main_SC = SC)]
SM <- rbindlist(lapply(list.files(file.path(OUT, "cache_I13c"), full.names = TRUE), readRDS))
G <- dcast(SM[, .(share = sum(share)), by = .(bin, grp = GRP[state], version, form)], bin + grp ~ form + version, value.var = "share")
X <- merge(O, G, by = c("bin", "grp"))
X[, `:=`(shift_main = main_SC - main_none, shift_log = log_SC - log_none, shift_exp = exp_SC - exp_none)]
fwrite(X, file.path(OUT, "I13d_tsd_shape_group_shares.csv"))
pc <- function(x) sprintf("%+.0f", 100 * x); p0 <- function(x) sprintf("%.0f", 100 * x)
for (g in c("boreal", "pioneer", "temperate", "regen")) { cat("\n====", g, "(% of stands) ====\n")
  print(X[grp == g, .(bin, observed = p0(obs), undist_log = p0(log_none), undist_exp = p0(exp_none), regime_main = p0(main_SC), regime_log = p0(log_SC), regime_exp = p0(exp_SC),
                      shift_main = pc(shift_main), shift_log = pc(shift_log), shift_exp = pc(shift_exp))]) }
