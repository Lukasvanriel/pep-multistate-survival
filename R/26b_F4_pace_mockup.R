# 26b_F4_pace_mockup.R (no INLA; 2026-11-12; MOCK-UP, user request) -- F4 with a new panel b: pace = GDD effect on the exchange speed.
# Symmetric form of a reversible chain (Kelly 1979): q_ij = s_ij * sqrt(pi_j / pi_i), s_ij = s_ji. In logs, per SD GDD:
#   beta_ij = sigma_ij + 1/2 (b_j - b_i), so sigma_ij = beta_ij - 1/2 * null (null = b_to - b_from from 24's multinomial) = GDD effect on the shared exchange speed.
#   SE of sigma = SE of beta (multinomial uncertainty ignored, as in 24). Class means: weighted (1/SE^2) linear model, class intercepts only.
# Panel a = 26's panel a unchanged. Output: figures/plan/F4_pace_mockup.png (does not replace F4_null_model_draft).
src <- readLines(here::here("R", "26_plan_figure_drafts.R")); eval(parse(text = src[1:(grep("^p4b <- ", src) - 1)]))
obs[, sigma := log_hr - null / 2]
mm <- lm(sigma ~ 0 + cls, obs, weights = 1 / se^2); ci <- confint(mm)
CM <- data.table(cls = factor(levels(obs$cls), levels(obs$cls)), est = coef(mm), lo = ci[, 1], hi = ci[, 2])
print(CM[, .(cls, n = obs[, .N, by = cls][match(CM$cls, cls), N], hr = round(exp(est), 2), lo = round(exp(lo), 2), hi = round(exp(hi), 2))])
p4b2 <- ggplot(obs, aes(cls, exp(sigma))) + geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_jitter(aes(colour = cls, size = 1 / se), width = 0.15, height = 0, alpha = 0.55) +
  geom_errorbar(data = CM, aes(y = exp(est), ymin = exp(lo), ymax = exp(hi)), width = 0, linewidth = 0.9) +
  geom_point(data = CM, aes(y = exp(est)), size = 3) +
  scale_y_log10() + scale_colour_manual(values = CC, guide = "none") + scale_size_continuous(range = c(0.6, 3), guide = "none") +
  labs(x = NULL, y = "GDD effect on the exchange speed (hazard ratio, log scale)",
       subtitle = "b  pace: speed-up by warmth shared by both directions (black = mean, 95% CI)")
save_fig(p4a + p4b2 + plot_annotation(title = "F4  Climate direction is sorting; climate pace is regen closure (mock-up)"), "F4_pace_mockup", 12, 5, NULL)
