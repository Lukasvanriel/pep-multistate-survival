# ==============================================================================
# 12_cross_transition_synthesis.R   (no INLA; seconds; run from the project root: Rscript R/12_cross_transition_synthesis.R)
# Needs: 08_final_pool.R and 09_hazard_ratios.R run (pooled_int_final.csv, hazard_ratios_int.csv); package metafor.
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
# MODELS (each fitted on two subsets: "all" and "robust" = >= MIN_EVENTS_ROBUST events, not low_events (15-19 disturbed events),
#   alpha_cov_low_share < ALPHA_COV_MAX):
#   A  Climate main effects (per +1 SD; undisturbed plots), for GDD, CMI, soil:
#        A1 class means                         y ~ class - 1
#        A2 between-state: thermophilization    y ~ dIndex          (slope > 0: warmer conditions speed moves toward warmer-affinity states)
#        A3 regen exits: destination affinity   y ~ Index(dest)
#   B  Disturbance main effects (vs none, at mean GDD, TSD = 0):
#        B1 cell means disturbance x class      y ~ cell - 1
#        B2 between-state: do disturbances favour moves toward warm-affinity states?   y ~ dist + dSTI_gdd
#   C  GDD x disturbance interactions (log ratio of the +1 SD GDD hazard ratios, disturbed vs undisturbed):
#        C1 overall mean (is antagonism general?)  y ~ 1
#        C2 by disturbance type                    y ~ dist - 1
#        C3 by class                               y ~ class - 1
#        C4 between-state                          y ~ dSTI_gdd     (does disturbance weaken the climate filter on thermophilic moves?)
#   D  Time since disturbance (log HR per decade), harvest light/heavy/any share one TSD term and are counted once:
#        D1 cell means disturbance family x class  y ~ cell - 1
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
MIN_EVENTS_ROBUST <- 25L; ALPHA_COV_MAX <- 0.2; MIN_EVENTS_ALPHA <- 50L; MIN_CELL <- 3L

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
hr[, class := fifelse(from == "regen", "regen_exit", fifelse(to == "regen", "to_regen", "between"))]
hr[, robust := n_events >= MIN_EVENTS_ROBUST & !(low_events %in% TRUE) & alpha_cov_low_share < ALPHA_COV_MAX]
add_idx <- function(d) {
  d <- merge(d, idx[, .(from = state, STI_from = STI, SMI_from = SMI, SSI_from = SSI)], by = "from")
  d <- merge(d, idx[, .(to = state, STI_to = STI, SMI_to = SMI, SSI_to = SSI)], by = "to")
  d[, `:=`(dSTI = STI_to - STI_from, dSMI = SMI_to - SMI_from, dSSI = SSI_to - SSI_from)][]
}
hr <- add_idx(hr)
hr[, family := fifelse(grepl("^harv", dist), "harvest", dist)]

clim <- hr[contrast %in% c("per_sd_gdd_3", "per_sd_cmi_gs", "per_sd_soil_pc1")]
clim[, cov := c(per_sd_gdd_3 = "gdd", per_sd_cmi_gs = "cmi", per_sd_soil_pc1 = "soil")[contrast]]
clim[, `:=`(dIndex = fcase(cov == "gdd", dSTI, cov == "cmi", dSMI, cov == "soil", dSSI),
            Index_to = fcase(cov == "gdd", STI_to, cov == "cmi", SMI_to, cov == "soil", SSI_to))]
dmain <- hr[grepl("_vs_none z=0 tau=0$", contrast)]
inter <- hr[grepl("_gdd_interaction$", contrast)]
tsd   <- unique(hr[grepl("_tsd_per_decade$", contrast)], by = c("label", "family"))   # light/heavy/any harvest share one TSD term

pl <- fread(file.path(DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
al <- pl[parameter == "alpha", .(label, log_hr = log_est, se = log_se, n_events = n_events_avg)]
al[, c("from", "to") := tstrsplit(label, "->", fixed = TRUE)]
al <- merge(al, traits[, .(from = state, shade_tol)], by = "from", all.x = TRUE)
al[, robust := n_events >= MIN_EVENTS_ALPHA]

fwrite(rbind(clim[, .(block = "A_climate", sub = cov, label, from, to, class, log_hr, se, n_events, robust, dIndex, Index_to)],
             dmain[, .(block = "B_disturbance", sub = dist, label, from, to, class, log_hr, se, n_events, robust, dIndex = dSTI, Index_to = STI_to)],
             inter[, .(block = "C_interaction", sub = dist, label, from, to, class, log_hr, se, n_events, robust, dIndex = dSTI, Index_to = STI_to)],
             tsd[,   .(block = "D_tsd", sub = family, label, from, to, class, log_hr, se, n_events, robust, dIndex = dSTI, Index_to = STI_to)],
             al[,    .(block = "E_alpha", sub = "log_alpha", label, from, to, class = fifelse(from == "regen", "regen_exit", fifelse(to == "regen", "to_regen", "between")),
                       log_hr, se, n_events, robust, dIndex = shade_tol, Index_to = NA_real_)]),
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
  d <- S(dmain)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$dist, d$class), ~ cell - 1, "B1_disturbance_cell_means", sb)
  out[[length(out) + 1]] <- fit_mv(d[class == "between"], ~ dist + dSTI, "B2_disturbance_between_vs_dSTI", sb)
  d <- S(inter)
  out[[length(out) + 1]] <- fit_mv(d, ~ 1, "C1_interaction_overall", sb)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$dist), ~ cell - 1, "C2_interaction_by_dist", sb)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$class), ~ cell - 1, "C3_interaction_by_class", sb)
  out[[length(out) + 1]] <- fit_mv(d[class == "between"], ~ dSTI, "C4_interaction_between_vs_dSTI", sb)
  d <- S(tsd)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$family, d$class), ~ cell - 1, "D1_tsd_cell_means", sb)
  d <- if (sb == "all") al else al[robust == TRUE]
  out[[length(out) + 1]] <- fit_mv(d[from != "regen" & !is.na(shade_tol)], ~ shade_tol, "E1_logalpha_vs_shade_tol", sb)
  out[[length(out) + 1]] <- fit_mv(cells(d, d$from), ~ cell - 1, "E2_logalpha_origin_means", sb, random = c("to", "label"))
}
res <- rbindlist(out, fill = TRUE)
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
key("B2_disturbance_between_vs_dSTI", "dSTI"); key("C1_interaction_overall", "intrcpt"); key("C4_interaction_between_vs_dSTI", "dSTI")
key("E1_logalpha_vs_shade_tol", "shade_tol")
