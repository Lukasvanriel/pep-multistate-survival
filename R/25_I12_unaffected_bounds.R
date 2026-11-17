# ==============================================================================
# 25_I12_unaffected_bounds.R   (no INLA; ~1-2 min; run from the project root)
# I12: bounds for the "unaffected exits" assumption (exits without a fire / pest term keep their undisturbed rate after that disturbance).
# Versions for those exits: x1 = current (no effect); class mean = mean of the point coefficients (main, GDD x dist, TSD) of the transitions of the
#   same class (alpha report: between / to_regen / regen_exit) that have the term; class max = coefficients of the transition of that class with the
#   largest effect at the median observed TSD (main + f(TSD) x TSD term, GDD 0).
# Engine and outcomes as D1 (18): disturbance at t = 0, GDD 0, CMI and soil 0, 0 y in the state at t = 0; point estimates.
# ==============================================================================
src <- readLines(here::here("R", "18_figures_D1_recovery.R")); eval(parse(text = src[1:(grep("^SCD <- CJ", src) - 1)]))
ROWS <- 1L; source(here::here("R", "covariate_spec.R")); OUTI <- file.path(FIG_DIR, "I")
ar <- fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[, .(label, class)]; CLS <- setNames(ar$class, ar$label)
TQ <- c(fire = vi[fire == 1, median(tsd_fire_y, na.rm = TRUE)], pest = vi[pest == 1, median(tsd_pest_y, na.rm = TRUE)])
PAR0 <- PAR; imp <- list()
run_version <- function(d, ver) { PAR <<- PAR0
  if (ver != "x1") { tm <- TERMS[[d]]; has <- names(PAR0)[sapply(PAR0, function(X) tm[1] %in% colnames(X))]
    for (cl in unique(CLS)) { donors <- intersect(has, names(CLS)[CLS == cl]); rec <- setdiff(names(CLS)[CLS == cl], has); if (!length(donors) || !length(rec)) next
      cf <- t(sapply(donors, function(lb) { X <- PAR0[[lb]][1, ]; g <- function(v) if (v %in% names(X)) X[[v]] else 0; c(g(tm[1]), g(tm[2]), g(tm[3])) }))
      v <- if (ver == "class mean") colMeans(cf) else cf[which.max(cf[, 1] + f_tsd(TQ[[d]]) * cf[, 3]), ]
      imp[[length(imp) + 1]] <<- data.table(dist = d, version = ver, class = cl, n_donors = length(donors), n_receivers = length(rec), main = v[1], gdd_x = v[2], tsd = v[3],
        effect_at_median_tsd = v[1] + f_tsd(TQ[[d]]) * v[3], donor = if (ver == "class max") donors[which.max(cf[, 1] + f_tsd(TQ[[d]]) * cf[, 3])] else NA_character_)
      for (lb in rec) { X <- PAR0[[lb]]; add <- matrix(rep(v, each = nrow(X)), nrow(X), dimnames = list(NULL, tm)); PAR[[lb]] <<- cbind(X, add) } } }
  summarise_hist(propagate_hist(zg = 0, dist = d, a0 = 0), list(dist = d, version = ver)) }
R <- rbindlist(lapply(c("fire", "pest"), function(d) rbindlist(lapply(c("x1", "class mean", "class max"), function(v) run_version(d, v)))))
PAR <- PAR0; IMP <- rbindlist(imp); fwrite(IMP, file.path(OUTI, "I12_imputed_terms.csv"))
R[, grp := fcase(outcome %in% c("resistance", "recovery_via_regen", "return_without_regen"), "in start state", outcome == "in_regen", "regen",
                 grepl("temperate", outcome), "temperate", grepl("boreal", outcome), "boreal", default = "pioneer")]
Q <- R[t <= 50, .(est = sum(est)), by = .(dist, version, start, grp, t)]
O <- merge(Q[t == 50, .(dist, version, start, grp, at50 = est)], Q[grp == "in start state", .(years_in_state = sum(est)), by = .(dist, version, start)][, grp := "in start state"], by = c("dist", "version", "start", "grp"), all.x = TRUE)
share <- fread(file.path(FIG_DIR, "B", "B_flags_share_assumed.csv"))[, .(start = origin, dist, share_assumed)]
O <- merge(O, share, by = c("start", "dist"), all.x = TRUE); fwrite(O, file.path(OUTI, "I12_unaffected_bounds.csv"))
options(width = 220); print(IMP[, .(dist, version, class, n_donors, n_receivers, eff = round(effect_at_median_tsd, 2), donor)])
W <- dcast(O[grp == "in start state"], dist + start + round(share_assumed, 2) ~ version, value.var = "years_in_state")
cat("\nyears in the start state out of the first 50:\n"); print(W[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
W2 <- dcast(O, dist + start ~ grp + version, value.var = "at50")
cat("\nshare at 50 y (in start / regen / temperate) x1 | mean | max:\n")
print(O[grp %in% c("in start state", "temperate", "regen"), .(v = paste(sprintf("%.2f", at50[match(c("x1", "class mean", "class max"), version)]), collapse = " | ")), by = .(dist, start, grp)][, dcast(.SD, dist + start ~ grp, value.var = "v")])
