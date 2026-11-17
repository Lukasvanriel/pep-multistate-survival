# ==============================================================================
# 27_hr_matrices.R   (no INLA; run from the project root; 2026-11-12, user request)
# Matrices (origin x destination) of every raw model coefficient as a hazard ratio, exp(pooled estimate) with its 95% CI (Rubin pooling,
#   pooled_final_final.csv): climate and soil (per SD), disturbance main effects (at TSD = 0, mean GDD), GDD x disturbance, TSD terms.
#   Same heatmap style as A4a (13_figures_A.R). Written to docs/chapter2_proposed_route/figures/hazard_ratios/matrices/ (+ one CSV of all values).
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
OUTM <- here::here("docs", "chapter2_proposed_route", "figures", "hazard_ratios", "matrices"); dir.create(OUTM, recursive = TRUE, showWarnings = FALSE)
PARS <- c(z_gdd_3 = "GDD", z_cmi_gs = "CMI", z_soil_pc1 = "soil PC1",
          fire_f = "fire", pest_f = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest",
          gdd_x_fire = "GDD × fire", gdd_x_pest = "GDD × pest", gdd_x_harv_light = "GDD × light harvest", gdd_x_harv_heavy = "GDD × heavy harvest",
          tsd_fire = "time since fire", tsd_pest = "time since pest", tsd_harv_light = "time since light harvest", tsd_harv_heavy = "time since heavy harvest")
SUB <- c(z_gdd_3 = "climate and soil", fire_f = "disturbance vs none, right after the disturbance, at mean GDD",
         gdd_x_fire = "change in the GDD effect when disturbed", tsd_fire = "change with time since disturbance (log scale of years)")
grpsub <- function(k) SUB[[c(z_gdd_3 = 1, z_cmi_gs = 1, z_soil_pc1 = 1, fire_f = 2, pest_f = 2, harv_light = 2, harv_heavy = 2,
                             gdd_x_fire = 3, gdd_x_pest = 3, gdd_x_harv_light = 3, gdd_x_harv_heavy = 3)[k] %||% 4]]
`%||%` <- function(a, b) if (is.na(a) || is.null(a)) b else a
cc <- pooled[parameter %in% names(PARS), .(label, from, to, parameter, hr = exp(pooled_est), hr_lo = exp(ci_lo), hr_hi = exp(ci_hi))]
stopifnot(nrow(cc) > 0, all(cc$hr_lo <= cc$hr & cc$hr <= cc$hr_hi))
cc[, sig := hr_lo > 1 | hr_hi < 1]
fwrite(cc[, .(transition = label, from, to, coefficient = PARS[parameter], hazard_ratio = hr, ci_low = hr_lo, ci_high = hr_hi, ci_excludes_1 = sig)],
       file.path(OUTM, "all_hazard_ratios.csv"))
i <- 0
for (k in names(PARS)) { i <- i + 1
  x <- cc[parameter == k]; lim <- cap_sym(log(x$hr))
  g <- merge(full_grid(), x, by = c("from", "to"), all.x = TRUE)
  g[, status := fifelse(paste0(from, "->", to) %in% modelled, fifelse(is.na(hr), "noterm", "value"), "notmod")]
  g[, `:=`(fillv = log(hr), txt = fmt_hr(hr), bold = sig %in% TRUE)]
  p <- heat(g, "hazard\nratio", limits = c(-lim, lim), trans_lab = function(v) fmt_hr(exp(v)), breaks = hr_breaks(lim),
            title = paste("Hazard ratio:", PARS[[k]]), subtitle = paste0(grpsub(k), "; bold = 95% CI excludes 1; '–' = term not in this transition"))
  ggsave(file.path(OUTM, sprintf("%02d_%s.png", i, gsub("[^a-z0-9]+", "_", tolower(PARS[[k]])))), p, width = 7.5, height = 7, dpi = 200)
  cat(sprintf("%-28s %3d transitions\n", PARS[[k]], nrow(x))) }
