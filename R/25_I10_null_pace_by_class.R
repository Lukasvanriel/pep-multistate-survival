# ==============================================================================
# 25_I10_null_pace_by_class.R   (no INLA; <1 min; run from the project root)
# I10: A4f Test 1 with a pace intercept per transition class (alpha report: to_regen, regen_exit, between). Same occupancy gradients (multinomial
#   logit of visit state on z-GDD, z-CMI, z-soil) and weights (1/se^2) as 24. (i) class intercepts + common slope; (ii) class intercepts + class slopes.
#   Intercept = GDD log HR at zero expected sorting difference = pace effect of warmth in that class.
# ==============================================================================
src <- readLines(here::here("R", "24_null_model_gdd_directionality.R")); eval(parse(text = src[1:(grep("^t1 <- lm", src) - 1)]))
OUTI <- file.path(FIG_DIR, "I")
ar <- fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[, .(label, class)]
obs <- merge(obs, ar, by = "label"); obs[, class := factor(class, c("between", "to_regen", "regen_exit"))]
m0 <- lm(log_hr ~ null, obs, weights = 1 / se^2); m1 <- lm(log_hr ~ 0 + class + null, obs, weights = 1 / se^2); m2 <- lm(log_hr ~ 0 + class + class:null, obs, weights = 1 / se^2)
tab <- function(m, name) { ci <- confint(m); data.table(model = name, term = names(coef(m)), est = coef(m), lo = ci[, 1], hi = ci[, 2]) }
R <- rbind(tab(m0, "A4f: one intercept"), tab(m1, "class intercepts, common slope"), tab(m2, "class intercepts and slopes"))
R <- merge(R, rbind(data.table(model = c("A4f: one intercept", "class intercepts, common slope", "class intercepts and slopes"),
  R2 = c(summary(m0)$r.squared, NA, NA), AIC = c(AIC(m0), AIC(m1), AIC(m2)))), by = "model")
fwrite(R, file.path(OUTI, "I10_null_pace_by_class.csv"))
print(obs[, .N, by = class]); print(R[, .(model, term, est = round(est, 2), lo = round(lo, 2), hi = round(hi, 2), AIC = round(AIC, 1))])
cat("F-test class intercepts vs one:\n"); print(anova(m0, m1)); cat("F-test class slopes vs common:\n"); print(anova(m1, m2))
