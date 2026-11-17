# ==============================================================================
# 25_I13_disturbance_regime.R   (no INLA; resumable; run from the project root)
# I13: equilibrium composition under a stationary disturbance regime, per I8 GDD bin, vs the observed composition.
# Rates (new events inside observed intervals / plot-years, by state at the interval start), four kinds as in the model: fire, heavy harvest
#   (stand-origin), light harvest, pest (partial; pest intensities merged as in the model). Two versions (user decision): (S) Poisson per kind,
#   log rate ~ state; (SC) ~ state + z-GDD + z-CMI, evaluated at each bin's medians. Plot-years as offset.
# Microsimulation of individual stands, 1-y steps, point parameters, eta as 22 (climate, soil, disturbance terms with GDD interactions and TSD).
#   Overlap as in the data (03a): a stand-origin event (fire, heavy harvest) ends all older events; same type replaces and restarts TSD (light and
#   heavy harvest are one type); partial types coexist. Same-year events applied pest -> light harvest -> fire -> heavy harvest. Disturbances carry
#   over after a state change. Start: equal numbers in every state; composition = mean over years 801-1000.
# Check: version "none" (no events) vs I8's deterministic equilibrium.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
source(here::here("R", "covariate_spec.R")); f_tsd <- tsd_fun(TSD_FORM)
OUTI <- file.path(FIG_DIR, "I"); CI <- file.path(OUTI, "cache_I13"); dir.create(CI, showWarnings = FALSE)
scl <- readRDS(file.path(FIT_DIR, "scaling.rds")); zz <- function(x, v) (x - scl[var == v, mean]) / scl[var == v, sd]
NPS <- 1500L; TMAX <- 1000L; KEEP <- 801:1000; t0 <- Sys.time()
# ---- bins (as I8) ----
v <- vi[, .(id_pe, visit_end, state = state_end, state_start, L = t_end - t_start, z_g = zz(gdd_3, "gdd_3"), z_c = zz(cmi_gs, "cmi_gs"), z_s = zz(soil_pc1, "soil_pc1"))]
pm <- v[, .(pg = mean(z_g)), by = id_pe]; pm[, bin := cut(pg, quantile(pg, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
v <- merge(v, pm[, .(id_pe, bin)], by = "id_pe")
B <- v[, .(zg = median(z_g), zc = median(z_c), zs = median(z_s, na.rm = TRUE)), by = bin][order(bin)]
# ---- rates ----
ev <- as.data.table(readRDS(file.path(BUILD_DIR, "disturbance_events.rds"))$events)[dropped == FALSE & id_pe %in% vi$id_pe]
ev[, kind := fcase(type == "fire", "fire", type == "harv" & level == 2, "heavy", type == "harv" & level == 1, "light", type == "pest", "pest")]
ne <- ev[!is.na(start_no_mes), .(n = .N), by = .(id_pe, visit_end = start_no_mes, kind)]
KD <- c("fire", "heavy", "light", "pest")
RT <- lapply(setNames(KD, KD), function(k) { d <- merge(v, ne[kind == k, .(id_pe, visit_end, n)], by = c("id_pe", "visit_end"), all.x = TRUE); d[is.na(n), n := 0L]
  d[, state_start := factor(state_start, S)]
  list(S = glm(n ~ 0 + state_start + offset(log(L)), poisson, d), SC = glm(n ~ 0 + state_start + z_g + z_c + offset(log(L)), poisson, d), n = sum(d$n)) })
rate_tab <- function(ver, b) sapply(KD, function(k) { nd <- data.table(state_start = factor(S, S), z_g = B$zg[b], z_c = B$zc[b], L = 1)
  exp(predict(RT[[k]][[ver]], nd)) })                                             # [state, kind] per year
RTAB <- rbindlist(lapply(c("S", "SC"), function(ver) rbindlist(lapply(seq_len(nrow(B)), function(b) data.table(version = ver, bin = B$bin[b], state = S, rate_tab(ver, b))))))
fwrite(RTAB, file.path(OUTI, "I13_rates.csv")); cat("events per kind:", sapply(RT, `[[`, "n"), "\n")
# ---- microsimulation ----
gc_ <- function(X, v) if (v %in% names(X)) X[[v]] else 0
CO <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) PAR[[lb]][1, ])); names(CO) <- S
TT <- list(fire = c("fire_f", "gdd_x_fire", "tsd_fire"), pest = c("pest_f", "gdd_x_pest", "tsd_pest"),
           light = c("harv_light", "gdd_x_harv_light", "tsd_harv_light"), heavy = c("harv_heavy", "gdd_x_harv_heavy", "tsd_harv_heavy"))
simulate <- function(b, ver, seed) { set.seed(seed); zg <- B$zg[b]; zc <- B$zc[b]; zs <- B$zs[b]
  R <- if (ver == "none") matrix(0, nS, 4, dimnames = list(S, KD)) else rate_tab(sub("_(noharv|lightonly|heavyonly|harvonly)$", "", ver), b)
  if (grepl("_harvonly$", ver)) R[, c("fire", "pest")] <- 0                       # harvest only, no natural disturbance (user 2026-11-16)
  if (grepl("_noharv$", ver)) R[, c("light", "heavy")] <- 0                       # I13 no-harvest version (fire and pest only; 2026-11-11)
  if (grepl("_lightonly$", ver)) R[, "heavy"] <- 0                                # fire, pest, light harvest (no heavy harvest; 2026-11-11, Brice 2019 contrast)
  if (grepl("_heavyonly$", ver)) R[, "light"] <- 0                                # fire, pest, heavy harvest (no light harvest)
  P <- 1 - exp(-R)
  n <- NPS * nS; st <- rep(seq_len(nS), each = NPS); age <- numeric(n)
  fire <- logical(n); pest <- logical(n); harv <- integer(n); tf <- tp <- th <- numeric(n)          # active flags + TSD (harv: 1 light, 2 heavy)
  occ <- matrix(0, length(KEEP), nS)
  for (yr in seq_len(TMAX)) {
    if (ver != "none") { U <- matrix(runif(n * 4), n, 4); hit <- U < P[st, , drop = FALSE]
      i <- hit[, "pest"];  pest[i] <- TRUE; tp[i] <- 0
      i <- hit[, "light"]; harv[i] <- 1L; th[i] <- 0
      i <- hit[, "fire"];  fire[i] <- TRUE; tf[i] <- 0; pest[i] <- FALSE; harv[i] <- 0L
      i <- hit[, "heavy"]; harv[i] <- 2L; th[i] <- 0; fire[i] <- FALSE; pest[i] <- FALSE }
    newst <- st
    for (s in seq_len(nS)) { ii <- which(st == s); if (!length(ii) || !length(EX[[s]]$lb)) next
      H <- sapply(CO[[s]], function(X) { e <- X[["b0y"]] + gc_(X, "z_gdd_3") * zg + gc_(X, "z_cmi_gs") * zc + gc_(X, "z_soil_pc1") * zs
        add <- function(on, tm, tsd) ifelse(on, gc_(X, tm[1]) + gc_(X, tm[2]) * zg + gc_(X, tm[3]) * f_tsd(tsd + 0.5), 0)
        e <- e + add(fire[ii], TT$fire, tf[ii]) + add(pest[ii], TT$pest, tp[ii]) + add(harv[ii] == 1L, TT$light, th[ii]) + add(harv[ii] == 2L, TT$heavy, th[ii])
        exp(e) * ((age[ii] + 1)^X[["alpha"]] - age[ii]^X[["alpha"]]) })
      H <- matrix(H, nrow = length(ii)); tot <- rowSums(H); mv <- runif(length(ii)) < 1 - exp(-tot)
      if (any(mv)) { cs <- H[mv, , drop = FALSE] / tot[mv]
        k <- pmin(1L + rowSums(t(apply(cs, 1, cumsum)) < runif(sum(mv))), ncol(cs)); newst[ii[mv]] <- EX[[s]]$to[k] } }
    moved <- newst != st; st <- newst; age <- ifelse(moved, 0, age + 1); tf <- tf + 1; tp <- tp + 1; th <- th + 1
    if (yr %in% KEEP) occ[yr - KEEP[1] + 1, ] <- tabulate(st, nS) / n }
  stopifnot(!anyNA(st), abs(rowSums(occ) - 1) < 1e-9); data.table(bin = B$bin[b], version = ver, state = S, share = colMeans(occ)) }
RUNS <- CJ(b = seq_len(nrow(B)), ver = c("none", "S", "SC", "SC_noharv", "SC_lightonly", "SC_heavyonly", "SC_harvonly"))
for (r in seq_len(nrow(RUNS))) { f <- file.path(CI, sprintf("b%d_%s.rds", RUNS$b[r], RUNS$ver[r])); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > 120) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t1 <- Sys.time(); saveRDS(simulate(RUNS$b[r], RUNS$ver[r], 1300 + r), f); cat(sprintf("bin %d %s (%.0f s)\n", RUNS$b[r], RUNS$ver[r], as.numeric(Sys.time() - t1, units = "secs"))) }
cat("all I13 runs done\n")
