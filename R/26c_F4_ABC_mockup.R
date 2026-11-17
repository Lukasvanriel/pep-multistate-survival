# 26c_F4_ABC_mockup.R (no INLA; 2026-11-12; MOCK-UP, user request) -- F4 with three panels, all per pair of states with moves both ways (37 pairs).
# For a pair colder i <-> warmer j with GDD log hazard ratios l_ij (i -> j) and l_ji (j -> i), from 24 (via 26):
#   A  direction: D_obs = l_ij - l_ji (GDD effect on the rate ratio) vs D_null = b_j - b_i (GDD effect on the number ratio, multinomial); detailed balance: D_obs = D_null.
#   B  pace: (l_ij + l_ji) / 2 = GDD effect on the geometric mean of the two rates (turnover); no assumption needed. SE = 0.5 sqrt(s_ij^2 + s_ji^2) (separate fits).
#   C  x = D_obs - D_null (direction beyond balance, SE = se_D), y = pace; one labelled dot per pair.
# Multinomial uncertainty NOT included (draft): CIs and colours are optimistic. Output figures/plan/F4_ABC_mockup.png + F4_ABC_pairs.csv.
src <- readLines(here::here("R", "26_plan_figure_drafts.R")); eval(parse(text = src[1:(grep("^p4b <- ", src) - 1)]))
P <- copy(pa)
P[, `:=`(pace = (lij + lji) / 2, se_pace = 0.5 * sqrt(sij^2 + sji^2), dir = D_obs - D_null)]
P[, type := factor(fifelse(i == "regen" | j == "regen", "pair with regen", "between established states"), c("between established states", "pair with regen"))]
P[, pace_sig := abs(pace / se_pace) > 1.96]
nm <- function(x) gsub("_", " ", x); P[, lab := paste(nm(i), "↔", nm(j))]
setorder(P, -z); P[, id := .I]   # numbers on panel c = rows of F4_ABC_pairs.csv (ggrepel not available)
mb <- lm(pace ~ 0 + type, P, weights = 1 / se_pace^2); cib <- confint(mb)
MB <- data.table(type = factor(levels(P$type), levels(P$type)), est = coef(mb), lo = cib[, 1], hi = cib[, 2])
print(MB[, .(type, hr = round(exp(est), 2), lo = round(exp(lo), 2), hi = round(exp(hi), 2))])
TC <- c(`between established states` = "grey35", `pair with regen` = "#7B3294")
pB <- ggplot(P, aes(type, exp(pace))) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = exp(pace - 1.96 * se_pace), ymax = exp(pace + 1.96 * se_pace)), width = 0, colour = "grey80", position = position_jitter(width = 0.15, seed = 1)) +
  geom_point(aes(colour = type), size = 1.8, position = position_jitter(width = 0.15, seed = 1)) +
  geom_errorbar(data = MB, aes(y = exp(est), ymin = exp(lo), ymax = exp(hi)), width = 0, linewidth = 0.9) + geom_point(data = MB, aes(y = exp(est)), size = 3) +
  scale_y_log10() + scale_colour_manual(values = TC, guide = "none") +
  labs(x = NULL, y = "GDD effect on the geometric mean of the two rates (log scale)", subtitle = "b  turnover between the two states (black = mean, 95% CI)")
pC <- ggplot(P, aes(dir, exp(pace))) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = exp(pace - 1.96 * se_pace), ymax = exp(pace + 1.96 * se_pace)), width = 0, colour = "grey85") +
  geom_errorbarh(aes(xmin = dir - 1.96 * se_D, xmax = dir + 1.96 * se_D), height = 0, colour = "grey85") +
  geom_point(aes(colour = side), size = 2) +
  geom_text(data = P[side != "consistent with detailed balance" | pace_sig], aes(label = id), size = 2.6, nudge_x = 0.12, nudge_y = 0.03) +
  scale_y_log10() + scale_colour_manual(values = SIDE_COL, name = "deviation from detailed balance (95%)", drop = FALSE) +
  labs(x = "direction beyond balance: toward warmer state (+) or colder state (−) (log)", y = "turnover (log scale)",
       subtitle = "c  direction and turnover per pair (numbers: either significant; see list)")
# layout (user 2026-11-16): a and b side by side, c below; larger text
save_fig(((p4a + guides(colour = "none")) + pB + plot_layout(widths = c(1.4, 1))) / pC + plot_layout(heights = c(1, 1.1), guides = "collect") &
           theme(legend.position = "bottom", axis.text = element_text(size = 11), axis.title = element_text(size = 12),
                 plot.subtitle = element_text(size = 12), legend.text = element_text(size = 11), legend.title = element_text(size = 12)),
         "F4_ABC_mockup", 13, 12, NULL)
OUTT <- P[, .(id, pair = lab, type, n_cold_to_warm = nij, n_warm_to_cold = nji, dir = round(dir, 2), dir_z = round(z, 2), side,
              pace_hr = round(exp(pace), 2), pace_lo = round(exp(pace - 1.96 * se_pace), 2), pace_hi = round(exp(pace + 1.96 * se_pace), 2), pace_sig)]
fwrite(OUTT, file.path(OUTP, "F4_ABC_pairs.csv")); print(OUTT, nrows = 50)
