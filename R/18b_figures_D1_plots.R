# ==============================================================================
# 18b_figures_D1_plots.R   (after 18_figures_D1_recovery.R printed "all D1 scenarios done"; from the project root)
# D1 figures: outcome of a stand over time after a disturbance at t = 0 (definitions: header of 18_figures_D1_recovery.R).
# ==============================================================================
source(here::here("R", "13_fig_helpers.R"))
OUTD <- file.path(FIG_DIR, "D"); OUT <- OUTD; D1 <- rbindlist(lapply(list.files(file.path(OUTD, "cache"), "^d1_.*rds$", full.names = TRUE), readRDS))
share <- fread(file.path(FIG_DIR, "B", "B_flags_share_assumed.csv")); olab <- function(x) gsub("_", " ", x)
DLAB <- c(none = "undisturbed", fire = "fire", pest = "pest", harv_light = "light harvest", harv_heavy = "heavy harvest")
OC <- c(resistance = "stayed (never left)", recovery_via_regen = "back after regen (recovery)", return_without_regen = "back without regen", in_regen = "currently regen",
        reorganized_boreal = "boreal after regen", reorganized_pioneer = "pioneer after regen", reorganized_temperate = "temperate after regen",
        drift_boreal = "boreal, no regen", drift_pioneer = "pioneer, no regen", drift_temperate = "temperate, no regen")
OCOL <- c("stayed (never left)" = "grey35", "back after regen (recovery)" = "#1B7837", "back without regen" = "#A6DBA0", "currently regen" = "grey80",
          "boreal after regen" = "#2166AC", "pioneer after regen" = "#B8860B", "temperate after regen" = "#B2182B",
          "boreal, no regen" = "#92C5DE", "pioneer, no regen" = "#F6E8A6", "temperate, no regen" = "#F4A582")
D1[, outcome_l := factor(OC[outcome], levels = OC)]
ST <- setdiff(ORDER, "regen")
flagtxt <- function(o, d) if (d == "none") "" else { s <- share[origin == o & dist == d]; sprintf("\n%.0f%% assumed%s", 100 * s$share_assumed, if (s$flag) " *" else "") }
fwrite(D1[, .(start, outcome, dist, zg, a0, t, est, lo, hi)], file.path(OUTD, "D1_outcomes.csv"))

for (A0 in c(0, 40)) {   # time in state at t = 0: 0 y (consistent with C) and 40 y (established stand); both shown (user 2026-11-09)
for (d in names(DLAB)) { x <- D1[dist == d & zg == 0 & a0 == A0][, start := factor(start, levels = ST)]
  lf <- setNames(paste0(olab(ST), sapply(ST, flagtxt, d = d)), ST)
  p <- ggplot(x, aes(t, est, fill = outcome_l)) +
    geom_area(colour = "white", linewidth = 0.1) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) + scale_fill_manual(values = OCOL, name = "outcome at t") +
    scale_y_continuous(labels = function(v) paste0(100 * v, "%")) + facet_wrap(~ start, ncol = 4, labeller = as_labeller(lf)) +
    labs(x = sprintf("years since %s", if (d == "none") "start" else DLAB[[d]]), y = "probability",
         title = sprintf("D1a  What a stand becomes after %s (%d y in the state at t = 0)", DLAB[[d]], A0),
         subtitle = "'after regen' = passed through regen at least once; dashed = end of record (55 y)",
         caption = if (d == "none") NULL else "panel: share of the exit rate assumed unaffected by the disturbance (* = > 25%)")
  save_fig(p, paste0("D1a_outcomes_over_time_", d, "_a", A0), 15, 10, NULL) }

# D1b at 50 y: outcome shares by disturbance and GDD
x50 <- D1[t == 50 & a0 == A0][, `:=`(start = factor(start, levels = rev(ST)), dist_l = factor(DLAB[dist], levels = DLAB), gdd = factor(gsd(zg), levels = gsd(c(-1, 0, 1))))]
p <- ggplot(x50, aes(est, start, fill = outcome_l)) + geom_col(width = 0.8) + scale_fill_manual(values = OCOL, name = "outcome at 50 y") +
  scale_x_continuous(labels = function(v) paste0(100 * v, "%")) + scale_y_discrete(labels = olab) + facet_grid(gdd ~ dist_l) +
  labs(x = "share after 50 y", y = "start state", title = sprintf("D1b  Outcome 50 years after disturbance (%d y in the state at t = 0)", A0))
save_fig(p, paste0("D1b_outcomes_50y_by_gdd_a", A0), 17, 12, NULL)

# D1c key quantities at 50 y with CI: recovered (back in X, any route), reorganized toward temperate (after regen), still X
kq <- D1[t == 50 & a0 == A0 & outcome %in% c("resistance", "recovery_via_regen", "reorganized_temperate", "in_regen")]
kq[, `:=`(start = factor(start, levels = rev(ST)), dist_l = factor(DLAB[dist], levels = DLAB), outcome_l = factor(OC[outcome], levels = OC))]
kq <- merge(kq, share[, .(start = origin, dist, flag)], by = c("start", "dist"), all.x = TRUE); kq[is.na(flag), flag := FALSE]
p <- ggplot(kq, aes(est, start, colour = factor(zg))) + geom_errorbar(aes(xmin = lo, xmax = hi, group = factor(zg)), width = 0, orientation = "y", position = position_dodge(width = 0.7)) +
  geom_point(aes(shape = flag, group = factor(zg)), size = 1.5, fill = "white", position = position_dodge(width = 0.7)) +
  scale_colour_manual(values = c(`-1` = "#2166AC", `0` = "grey30", `1` = "#B2182B"), labels = c("GDD -1 SD", "GDD 0", "GDD +1 SD"), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), breaks = "TRUE", labels = c(`TRUE` = "> 25% of exit rate assumed unaffected"), name = NULL) +
  scale_x_continuous(labels = function(v) paste0(100 * v, "%")) + scale_y_discrete(labels = olab) + facet_grid(dist_l ~ outcome_l) +
  labs(x = "probability at 50 y", y = NULL, title = sprintf("D1c  Key outcomes 50 years after disturbance (%d y in the state at t = 0)", A0)) + theme(strip.text.y = element_text(angle = 0))
save_fig(p, paste0("D1c_key_outcomes_50y_CI_a", A0), 15, 16, kq)

}

# D1d check: time in state at the disturbance 0 vs 40 y (GDD 0)
ck <- D1[t == 50][, .(est = sum(est)), by = .(start, dist, zg, a0, grp = fifelse(outcome %in% c("resistance", "recovery_via_regen", "return_without_regen"), "in the start state", fifelse(outcome == "in_regen", "regen", "another state")))]
ck <- dcast(ck, start + dist + zg + grp ~ a0, value.var = "est"); setnames(ck, c("0", "40"), c("a0_0y", "a0_40y")); ck[, diff := a0_40y - a0_0y]
fwrite(ck, file.path(OUTD, "D1d_check_time_in_state_at_disturbance.csv"))
cat(sprintf("D1d: max |difference| (40 y vs 0 y in state at the disturbance) in the 50-y outcome groups = %.3f\n", max(abs(ck$diff))))
print(ck[order(-abs(diff))][1:8])
