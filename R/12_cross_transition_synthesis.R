# ==============================================================================
# 12_cross_transition_synthesis.R   (no INLA; seconds; run from the project root: Rscript R/12_cross_transition_synthesis.R)
# Needs: 08_final_pool.R and 09_hazard_ratios.R run (pooled_<tag>_final.csv, hazard_ratios_<tag>.csv); package metafor.
# WP7 (D5): no harv_any; TSD block D on the TSD coefficient (_tsd_coef); disturbance block B at each tau = 1, 10, 30, 60 y; no alpha_cov; transitions without a GDD main term
#   (A11) have no GDD rows in block A and no rows in blocks B and C (09 writes none; one-sided bounds are excluded from B).
#
# PURPOSE: replace ~1,300 per-transition contrasts by a few tests of ecological patterns ACROSS transitions.
# UNIT: one pooled estimate per transition (x contrast), log scale, with its Rubin SE (from 09). Analysed by multilevel random-effects
#   meta-analysis/meta-regression (metafor::rma.mv, REML, t-tests):
#     y_i = x_i'b + u_origin + u_destination + u_transition (+ u_obs if a transition contributes >1 estimate) + e_i,  e_i ~ N(0, se_i^2)
#   Crossed random effects for origin and destination state: transitions sharing a state are not independent ecological replicates.
#
# TRANSITION CLASSES:  regen_exit (regen -> X: stand development), to_regen (X -> regen: stand replacement), between (X -> Y, both non-regen).
# STATE THERMAL / MOISTURE / SOIL INDEX (STI, SMI, SSI): mean z-scored GDD / CMI / soil of all visits in that state (data-driven affinity, as the
#   community temperature index of Brice et al. 2019). dSTI = STI(destination) - STI(origin): > 0 means a move toward a warmer-affinity state.
#   CAVEAT for the text: the indices come from the same data, so "climate favours moves toward warm-affinity states" is partly expected;
#   the informative parts are the size/consistency of the slope and how disturbance modifies it.
#
# MODELS (each fitted on two subsets: "all" and "robust" = >= MIN_EVENTS_ROBUST events and no low_events label ("<15" or "15-19" disturbed events of the level, from 09)):
#   A  Climate main effects (per +1 SD; undisturbed plots), for GDD, CMI, soil:
#        A1 class means                         y ~ class - 1
#        A2 between-state: thermophilization    y ~ dIndex          (slope > 0: warmer conditions speed moves toward warmer-affinity states)
#        A3 regen exits: destination affinity   y ~ Index(dest)
#   B  Disturbance effects (vs none, at mean GDD z = 0), run once at each tau = TAU_B (1, 10, 30, 60 y) since the disturbance:
#        B1 cell means disturbance x class      y ~ cell - 1
#        B2 between-state: do disturbances favour moves toward warm-affinity states?   y ~ dist + dSTI_gdd
#      (09's one-sided-bound rows are not used; a level with only a main term has no tau rows.)
#   C  GDD x disturbance interactions (log ratio of the +1 SD GDD hazard ratios, disturbed vs undisturbed):
#        C1 overall mean (is antagonism general?)  y ~ 1
#        C2 by disturbance type                    y ~ dist - 1
#        C3 by class                               y ~ class - 1
#        C4 between-state                          y ~ dSTI_gdd     (does disturbance weaken the climate filter on thermophilic moves?)
#   D  Time since disturbance: the TSD coefficient c (log HR per unit of f(TSD) = log(TSD + 1)), one estimate per transition x type (fire, pest, harv_light, harv_heavy):
#        D1 cell means type x class                y ~ cell - 1
#      Every fixed-time contrast is c times a constant, so the tests are the same for any window; D1 rows also carry the results scaled to the windows 1 vs 10 y, 10 vs 30 y,
#      30 vs 60 y (est = c * (f(t2) - f(t1)), CI scaled the same way, ratio = exp).
#   E  Weibull shape: log(alpha) (pooled on the log scale, 08), origins other than regen:
#        E1 life history                           y ~ shade_tol(origin)      (prediction: negative slope; tolerant late-successional states ~ memoryless)
#        E2 origin means incl. regen (descriptive) y ~ origin - 1   (random: destination, transition)
#        "robust" subset for E = >= MIN_EVENTS_ALPHA events (the alpha prior dominates sparse transitions).
# Cells with < MIN_CELL estimates are dropped before fitting. Sampling covariances between estimates of the same transition are ignored
#   (absorbed in part by u_transition); state this in the text.
# OUTPUTS (RESULTS_DIR = data-output/<tag>/results/): synthesis_inputs_<tag>.csv (every estimate used, with moderators), synthesis_int.csv (all model terms).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
if (!requireNamespace("metafor", quietly = TRUE)) stop("Install metafor first: install.packages('metafor')")
suppressPackageStartupMessages(library(metafor))
options(width = 200)
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, FIT_DIR, RESULTS_DIR
DIR <- RESULTS_DIR
MIN_EVENTS_ROBUST <- 25L; MIN_EVENTS_ALPHA <- 50L; MIN_CELL <- 3L; TAU_B <- c(1, 10, 30, 60)
source(here::here("R", "covariate_spec.R")); f_tsd <- tsd_fun(TSD_FORM)   # tsd_fun(): window constants f(t2) - f(t1)
WINDOWS <- list(w1v10 = c(1, 10), w10v30 = c(10, 30), w30v60 = c(30, 60))

# ---- state traits: VERIFY / EDIT before use -----------------------------------------------------------------------------------------------------------
# shade_tol: Niinemets & Valladares (2006, Ecol. Monogr.) shade-tolerance scale 1 (intolerant) - 5 (very tolerant) of the dominant species.
# Values below are from memory and MUST be checked against the paper's appendix. Composite states: NA (excluded from E1) unless you set a value.
traits <- data.table(
  state     = c("fir", "black_spruce", "white_spruce", "jack_pine", "cedar", "sugar_maple", "red_maple", "yellow_birch",
                "white_birch", "trembling_aspen", "other_pioneer", "other_temperate", "regen"),
  shade_tol = c(5.01,  4.08,           4.15,           1.36,        3.45,    4.76,          3.44,        3.17,
                1.54,          1.21,              NA,              NA,                NA))

# ---- inputs -------------------------------------------------------------------------------------------------------------------------------------------------
iv <- readRDS(file.path(TABLE_DIR, "visit_intervals.rds")); setDT(iv)
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); setDT(sc)
z <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
idx <- iv[, .(STI = mean(z(gdd_3, "gdd_3")), SMI = mean(z(cmi_gs, "cmi_gs")), SSI = mean(z(soil_pc1, "soil_pc1")), n_visits = .N), by = .(state = state_end)]
cat("State indices (mean z-score of visits in the state):\n"); print(idx[order(STI)], digits = 3)
stopifnot(all(traits$state %in% idx$state))

hr <- fread(file.path(DIR, sprintf("hazard_ratios_%s.csv", RUN_TAG)))
hr[, c("from", "to") := tstrsplit(label, "->", fixed = TRUE)]
stopifnot("bound / low_events must hold NA, not empty strings (09 writes na = \"NA\")" = !any(hr$bound %in% "", hr$low_events %in% ""))
hr[, class := fifelse(from == "regen", "regen_exit", fifelse(to == "regen", "to_regen", "between"))]
hr[, robust := n_events >= MIN_EVENTS_ROBUST & is.na(low_events)]   # low_events: "<15" / "15-19" / NA (09)
add_idx <- function(d) {
  d <- merge(d, idx[, .(from = state, STI_from = STI, SMI_from = SMI, SSI_from = SSI)], by = "from")
  d <- merge(d, idx[, .(to = state, STI_to = STI, SMI_to = SMI, SSI_to = SSI)], by = "to")
  d[, `:=`(dSTI = STI_to - STI_from, dSMI = SMI_to - SMI_from, dSSI = SSI_to - SSI_from)][]
}
hr <- add_idx(hr)

clim <- hr[contrast %in% c("per_sd_gdd_3", "per_sd_cmi_gs", "per_sd_soil_pc1")]
clim[, cov := c(per_sd_gdd_3 = "gdd", per_sd_cmi_gs = "cmi", per_sd_soil_pc1 = "soil")[contrast]]
clim[, `:=`(dIndex = fcase(cov == "gdd", dSTI, cov == "cmi", dSMI, cov == "soil", dSSI),
            Index_to = fcase(cov == "gdd", STI_to, cov == "cmi", SMI_to, cov == "soil", SSI_to))]
dmain <- rbindlist(lapply(TAU_B, function(tt) hr[grepl(sprintf("_vs_none z=0 tau=%gy$", tt), contrast) & is.na(bound)][, tau_y := tt]))   # one-sided-bound rows are not estimates for the synthesis
inter <- hr[grepl("_gdd_interaction$", contrast)]
tsd   <- hr[grepl("_tsd_coef$", contrast)]                                            # one estimate per transition x type (no de-duplication)

pl <- fread(file.path(DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
al <- pl[parameter == "alpha", .(label, log_hr = log_est, se = log_se, n_events = n_events_avg)]
al[, c("from", "to") := tstrsplit(label, "->", fixed = TRUE)]
al <- merge(al, traits[, .(from = state, shade_tol)], by = "from", all.x = TRUE)
al[, robust := n_events >= MIN_EVENTS_ALPHA]

fwrite(rbind(clim[, .(block = "A_climate", sub = cov, label, from, to, class, log_hr, se, n_events, robust, dIndex, Index_to)],
             dmain[, .(block = "B_disturbance", sub = dist, tau_y, label, from, to, class, log_hr, se, n_events, robust, dIndex = dSTI, Index_to = STI_to)],
             inter[, .(block = "C_interaction", sub = dist, label, from, to, class, log_hr, se, n_events, robust, dIndex = dSTI, Index_to = STI_to)],
             tsd[,   .(block = "D_tsd", sub = dist, label, from, to, class, log_hr, se, n_events, robust, dIndex = dSTI, Index_to = STI_to)],
             al[,    .(block = "E_alpha", sub = "log_alpha", label, from, to, class = fifelse(from == "regen", "regen_exit", fifelse(to == "regen", "to_regen", "between")),
                       log_hr, se, n_events, robust, dIndex = shade_tol, Index_to = NA_real_)], fill = TRUE),
       file.path(DIR, sprintf("synthesis_inputs_%s.csv", RUN_TAG)))

# ---- meta-analytic fit ---------------------------------------------------------------------------------------------------------------------------------
fit_mv <- function(d, mods, model, subset, random = c("from", "to", "label")) {
  d <- copy(d)[!is.na(log_hr) & !is.na(se) & se > 0]
  if (nrow(d) < 5) return(data.table(model = model, subset = subset, term = NA_character_, k = nrow(d), note = "fewer than 5 estimates"))
  d[, `:=`(yi = log_hr, vi = se^2, obs = .I)]
  rnd <- lapply(random, function(r) as.formula(paste("~ 1 |", r)))
  if ("label" %in% random && anyDuplicated(d$label)) rnd <- c(rnd, list(~ 1 | obs))
  m <- tryCatch(rma.mv(yi, vi, mods = mods, random = rnd, data = d, method = "REML", test = "t"), error = function(e) e)
  if (inherits(m, "error")) return(data.table(model = model, subset = subset, term = NA_character_, k = nrow(d), note = conditionMessage(m)))
  data.table(model = model, subset = subset, term = rownames(m$b), k = m$k, est = m$b[, 1], se = m$se, ci_lo = m$ci.lb, ci_hi = m$ci.ub, p = m$pval,
             ratio = exp(m$b[, 1]), ratio_lo = exp(m$ci.lb), ratio_hi = exp(m$ci.ub),
             sigma2 = paste(sprintf("%s=%.3f", c(random, "obs")[seq_along(m$sigma2)], m$sigma2), collapse = " "),
             QM = m$QM, QMp = m$QMp, note = NA_character_)
}
cells <- function(d, ...) { d <- copy(d); d[, cell := do.call(paste, c(list(...), sep = ":"))]; d[, n_cell := .N, by = cell]; d[n_cell >= MIN_CELL] }
out <- list()
for (sb in c("all", "robust")) {
  S <- function(d) if (sb == "all") d else d[robust == TRUE]
  for (cv in c("gdd", "cmi", "soil")) {
    d <- S(clim[cov == cv])
    out[[length(out) + 1]] <- fit_mv(d, ~ class - 1, sprintf("A1_%s_class_means", cv), sb)
    out[[length(out) + 1]] <- fit_mv(d[class == "between"], ~ dIndex, sprintf("A2_%s_between_vs_dIndex", cv), sb)
    out[[length(out) + 1]] <- fit_mv(d[class == "regen_exit"], ~ Index_to, sprintf("A3_%s_regenexit_vs_destIndex", cv), sb, random = c("to", "label"))
  }
  for (tt in TAU_B) {
    d <- S(dmain[tau_y == tt])
    out[[length(out) + 1]] <- fit_mv(cells(d, d$dist, d$class), ~ cell - 1, sprintf("B1_disturbance_cell_means_tau%gy", tt), sb)
    out[[length(out) + 1]] <- fit_mv(d[class == "between"], ~ dist + dSTI, sprintf("B2_disturbance_between_vs_dSTI_tau%gy", tt), sb)
  }
  d <- S(inter)
  out[[length(out) + 1]] <- fit_mv(d, ~ 1, "C1_interaction_overall", sb)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$dist), ~ cell - 1, "C2_interaction_by_dist", sb)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$class), ~ cell - 1, "C3_interaction_by_class", sb)
  out[[length(out) + 1]] <- fit_mv(d[class == "between"], ~ dSTI, "C4_interaction_between_vs_dSTI", sb)
  d <- S(tsd)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$dist, d$class), ~ cell - 1, "D1_tsd_cell_means", sb)
  d <- if (sb == "all") al else al[robust == TRUE]
  out[[length(out) + 1]] <- fit_mv(d[from != "regen" & !is.na(shade_tol)], ~ shade_tol, "E1_logalpha_vs_shade_tol", sb)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$from), ~ cell - 1, "E2_logalpha_origin_means", sb, random = c("to", "label"))
}
res <- rbindlist(out, fill = TRUE)
# D1 (tsd coefficient c per unit of f(TSD)): results scaled to the fixed-time windows, est * (f(t2) - f(t1)); same sign, same p-value
for (w in names(WINDOWS)) { cst <- f_tsd(WINDOWS[[w]][2]) - f_tsd(WINDOWS[[w]][1]); is_d <- grepl("^D1_", res$model)
  res[is_d, (paste0(w, c("_est", "_lo", "_hi", "_ratio", "_ratio_lo", "_ratio_hi"))) := .(est * cst, ci_lo * cst, ci_hi * cst, exp(est * cst), exp(ci_lo * cst), exp(ci_hi * cst))] }
fwrite(res, file.path(DIR, sprintf("synthesis_%s.csv", RUN_TAG)))

# ---- printed summary -------------------------------------------------------------------------------------------------------------------------------------
cat("\nAll terms (est on the log scale; ratio = exp(est): HR, ratio of HRs, HR per unit of the moderator, or alpha for E2):\n")
print(res[, .(model, subset, term, k, est = round(est, 3), se = round(se, 3), ratio = round(ratio, 2), lo = round(ratio_lo, 2), hi = round(ratio_hi, 2),
              p = signif(p, 2), sigma2, QMp = signif(QMp, 2), note)], nrows = 500)
cat("\nKey tests (robust subset in brackets):\n")
key <- function(mod, trm) { r <- res[model == mod & term == trm]; if (!nrow(r)) return(invisible())
  f <- function(s) { x <- r[subset == s]; if (!nrow(x) || is.na(x$est)) "NA" else sprintf("%.3f (%.3f to %.3f), p = %.2g, k = %d", x$est, x$ci_lo, x$ci_hi, x$p, x$k) }
  cat(sprintf("  %-40s %-10s %s   [%s]\n", mod, trm, f("all"), f("robust"))) }
key("A2_gdd_between_vs_dIndex", "dIndex"); key("A2_cmi_between_vs_dIndex", "dIndex"); key("A3_gdd_regenexit_vs_destIndex", "Index_to")
for (tt in TAU_B) key(sprintf("B2_disturbance_between_vs_dSTI_tau%gy", tt), "dSTI")
key("C1_interaction_overall", "intrcpt"); key("C4_interaction_between_vs_dSTI", "dSTI")
key("E1_logalpha_vs_shade_tol", "shade_tol")
cat("\nD1 cell means scaled to windows (log HR change between the two times; ratio = exp):\n")
print(res[grepl("^D1_", model) & !is.na(est), .(subset, term, k, c_coef = round(est, 3), p = signif(p, 2), w1v10 = round(w1v10_ratio, 2), w10v30 = round(w10v30_ratio, 2), w30v60 = round(w30v60_ratio, 2))], nrows = 100)
# D5 check: every block has k > 0 (number of estimates entering, all / robust subset)
inp_all <- fread(file.path(DIR, sprintf("synthesis_inputs_%s.csv", RUN_TAG)))
kk <- inp_all[, .(k_all = .N, k_robust = sum(robust)), by = .(block, tau_y)][order(block, tau_y)]
cat("\nD5 check, estimates per block (all / robust):\n"); print(kk); stopifnot(nrow(kk) >= 4L + length(TAU_B) - 1L, all(kk$k_all > 0), all(kk$k_robust > 0))
cat("A11 check: GDD / B / C rows of transitions without a GDD main term:\n")
nog <- setdiff(hr$label, hr[contrast == "per_sd_gdd_3", label]); n_a11 <- inp_all[label %in% nog & (block %in% c("B_disturbance", "C_interaction") | (block == "A_climate" & sub == "gdd")), .N]; cat("  transitions without z_gdd_3:", paste(nog, collapse = ", "), "| rows:", n_a11, "\n"); stopifnot(n_a11 == 0L)
