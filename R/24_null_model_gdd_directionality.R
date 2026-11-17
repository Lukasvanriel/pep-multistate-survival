# ==============================================================================
# 24_null_model_gdd_directionality.R   (no INLA; ~1 min; run from the project root)
# Null-model test for the GDD directionality of the hazard ratios (A4; user 2026-11-10).
# Question: is "higher GDD -> faster moves toward warmer states" just a restatement of how the states are sorted along GDD in space?
# Occupancy gradients: multinomial logit of the state of every plot visit on z-GDD, z-CMI, z-soil (same conditioning as the HRs);
#   b_s = GDD slope of log(P(s)/P(ref)); b_j - b_i = GDD slope of log(P(j)/P(i)). Visits: first visits (plot_info) + end visits of all intervals
#   (interval rows carry end-visit covariates).
# Test 1 (all transitions): weighted regression (1/se^2) of the observed GDD log HR on b_to - b_from.
#   slope ~1: HRs mirror the occupancy gradient; ~0.5: equilibrium split equally between both directions; intercept > 0: pace effect of warmth.
# Test 2 (reverse pairs, detailed balance p_i q_ij = p_j q_ji -> logHR_ij - logHR_ji = b_j - b_i): pairs oriented so j = warmer state (b_j > b_i);
#   slope ~1 and intercept ~0 = equilibrium sorting; intercept > 0 = warm direction favoured beyond equilibrium; z = residual / se (se of the observed
#   difference only; occupancy-gradient uncertainty and plot clustering ignored -> |z| optimistic).
# Outputs: figures/A/A4f_null_model_tests.csv (summary), A4f_null_model_pairs.csv, A4f_null_model.png
# ==============================================================================
library(data.table); library(nnet)
source(here::here("R", "config.R"))
OUTA <- file.path(FIG_DIR, "A")
vi   <- readRDS(file.path(TABLE_DIR, "visit_intervals.rds")); pinf <- readRDS(file.path(TABLE_DIR, "plot_info.rds"))
hr   <- fread(file.path(RESULTS_DIR, sprintf("hazard_ratios_%s.csv", RUN_TAG))); sc <- readRDS(file.path(FIT_DIR, "scaling.rds"))
zf   <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
vis <- rbind(pinf[, .(id_pe, state = first_state, gdd = fv_gdd_3, cmi = fv_cmi_gs, soil = fv_soil_pc1)],
             vi[,   .(id_pe, state = state_end,   gdd = gdd_3,    cmi = cmi_gs,    soil = soil_pc1)])
vis <- vis[complete.cases(vis)]
vis[, `:=`(zg = zf(gdd, "gdd_3"), zc = zf(cmi, "cmi_gs"), zs = zf(soil, "soil_pc1"), state = factor(state))]
cat("visits:", nrow(vis), " plots:", uniqueN(vis$id_pe), "\n")
fit <- multinom(state ~ zg + zc + zs, data = vis, maxit = 500, trace = FALSE)
b <- c(setNames(0, levels(vis$state)[1]), coef(fit)[, "zg"]); print(round(sort(b), 2))
obs <- hr[contrast == "per_sd_gdd_3", .(label, log_hr, se, n_events)]
obs[, c("from", "to") := tstrsplit(label, "->", fixed = TRUE)]; obs[, null := b[to] - b[from]]
t1 <- lm(log_hr ~ null, data = obs, weights = 1 / se^2)
pr <- merge(obs[, .(i = from, j = to, lij = log_hr, sij = se, nij = n_events)], obs[, .(i = to, j = from, lji = log_hr, sji = se, nji = n_events)], by = c("i", "j"))
stopifnot(nrow(pr) %% 2 == 0); pr <- pr[b[j] > b[i]]                       # j = warmer state; each pair once
pr[, `:=`(D_obs = lij - lji, D_null = b[j] - b[i], se_D = sqrt(sij^2 + sji^2))]; pr[, `:=`(resid = D_obs - D_null, z = (D_obs - D_null) / se_D)]
t2 <- lm(D_obs ~ D_null, data = pr, weights = 1 / se_D^2); t2b <- lm(D_obs ~ 0 + D_null, data = pr, weights = 1 / se_D^2)
ci <- function(m, k) confint(m)[k, ]
S <- rbind(data.table(test = "Test 1 all transitions", n = nrow(obs), term = c("intercept", "slope"), est = coef(t1), lo = c(ci(t1, 1)[1], ci(t1, 2)[1]), hi = c(ci(t1, 1)[2], ci(t1, 2)[2]), R2 = summary(t1)$r.squared),
           data.table(test = "Test 2 reverse pairs", n = nrow(pr), term = c("intercept (warm direction)", "slope"), est = coef(t2), lo = c(ci(t2, 1)[1], ci(t2, 2)[1]), hi = c(ci(t2, 1)[2], ci(t2, 2)[2]), R2 = summary(t2)$r.squared),
           data.table(test = "Test 2 reverse pairs", n = nrow(pr), term = "slope through origin", est = coef(t2b), lo = ci(t2b, 1)[1], hi = ci(t2b, 1)[2], R2 = NA_real_))
print(S); cat(sprintf("pairs |z| > 2: %d (warm stronger: %d; weaker: %d)\n", sum(abs(pr$z) > 2), sum(pr$z > 2), sum(pr$z < -2)))
fwrite(S, file.path(OUTA, "A4f_null_model_tests.csv"))
fwrite(pr[order(-abs(z)), .(colder = i, warmer = j, D_obs, D_null, resid, z, n_colder_to_warmer = nij, n_warmer_to_colder = nji)], file.path(OUTA, "A4f_null_model_pairs.csv"))
png(file.path(OUTA, "A4f_null_model.png"), width = 2200, height = 1000, res = 200); par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(obs$null, obs$log_hr, xlab = "expected from spatial sorting", ylab = "observed GDD effect (log hazard ratio)", main = "A4f  Test 1: all transitions",
     pch = 19, col = ifelse(abs(obs$log_hr / obs$se) > 1.96, "black", "grey")); abline(0, 1, lty = 2); abline(t1, col = "red")
plot(pr$D_null, pr$D_obs, xlab = "expected from spatial sorting", ylab = "observed: warm move − cold move (log hazard ratio)", main = "A4f  Test 2: reverse pairs",
     pch = 19, col = ifelse(abs(pr$z) > 2, "red", "grey40")); abline(0, 1, lty = 2); abline(t2, col = "red"); dev.off()
cat("24 done\n")
