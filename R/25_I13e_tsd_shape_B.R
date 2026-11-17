# ==============================================================================
# 25_I13e_tsd_shape_B.R   (no INLA; resumable; run from the project root, after R/simulation/sim06h_class_tau.R; 2026-11-12)
# Option B for I13: as 25_I13c_tsd_shape.R, but the exponential coefficients come from sim06h (every TSD term exp(-TSD / tau), tau per
#   disturbance type x transition class, sim06e_best_tau_by_class.csv), and each transition uses its own class's tau in the simulation.
#   log = sim06 pass 1 "log_all" (same 5 imputations). Imputations failing the alpha-SD screen in either shape are dropped from both.
#   Same microsimulation, rates (SC), bins and random numbers as 25_I13c. Versions none and SC. Prints the shift table (log, exp A, exp B).
# Self-test (I13E_TEST=1): uses the option-A fits and taus instead; must reproduce 25_I13c's exp results exactly (written to cache_I13e_test).
# ==============================================================================
TEST <- Sys.getenv("I13E_TEST") == "1"
src <- readLines(here::here("R", "25_I13_disturbance_regime.R")); eval(parse(text = src[1:(grep("^# ---- microsimulation", src) - 1)]))
CI3 <- file.path(OUTI, if (TEST) "cache_I13e_test" else "cache_I13e"); dir.create(CI3, showWarnings = FALSE); t0 <- Sys.time()
D6 <- here::here("data-output", "simulation", "sim06_tsd_tau", "fits")
ar <- fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG)))[, .(label, class)]
BTC <- fread(here::here("data-output", "simulation", "sim06d_regen_screen", "sim06e_best_tau_by_class.csv"))
KN <- c(fire = "tsd_fire", pest = "tsd_pest", light = "tsd_harv_light", heavy = "tsd_harv_heavy")
tau_lb <- function(lb) { if (TEST) return(c(fire = 15, pest = 15, light = 25, heavy = 160)); b <- BTC[class == ar[label == lb, class]]; setNames(b$tau_best[match(KN, b$term)], names(KN)) }
FB <- if (TEST) "^p2_fire15_harv_heavy160_harv_light25_pest15_" else "^p3class_"
rd <- function(f, nm) { r <- readRDS(f); x <- r$fits[[nm]]; if (!is.na(x$error)) return(NULL); data.table(label = r$label, m = r$m, ne = r$n_events, x$coef[, .(parameter, mean, sd)]) }
fb <- list.files(D6, FB, full.names = TRUE); stopifnot("sim06h fits missing: run R/simulation/sim06h_class_tau.R first" = length(fb) == 360)
cf <- rbind(rbindlist(lapply(list.files(D6, "^p1_", full.names = TRUE), rd, "log_all"))[, form := "log"], rbindlist(lapply(fb, rd, "best_all"))[, form := "exp"])
cf <- cf[cf[, .(nf = uniqueN(form)), by = .(label, m)][nf == 2, .(label, m)], on = c("label", "m")]          # INLA errors: drop the pair
bad <- unique(cf[parameter == "alpha" & sd / (0.78 * mean / sqrt(ne)) < 0.1, .(label, m)]); cf <- cf[!bad, on = c("label", "m")]
cat("dropped (label, imputation):", paste(bad$label, bad$m, collapse = "; "), "\n"); print(cf[parameter == "alpha", .(n_imp = .N / 2), by = label][, .N, by = n_imp])
K <- cf[, .(mean = mean(mean)), by = .(form, label, parameter)]
make_CO <- function(fm) { o <- lapply(S, function(s) lapply(EX[[s]]$lb, function(lb) { X <- PAR[[lb]][1, ]; k <- K[form == fm & label == lb]
  if (nrow(k)) { Y <- setNames(k$mean, k$parameter); stopifnot(setequal(names(Y), setdiff(names(X), "b0y"))); X <- c(Y, b0y = Y[["(Intercept)"]] - Y[["alpha"]] * log(10)) }
  attr(X, "tau") <- tau_lb(lb); X })); names(o) <- S; o }
COS <- list(log = make_CO("log"), exp = make_CO("exp"))
FT <- list(log = function(x, k, X) log1p(x), exp = function(x, k, X) exp(-x / attr(X, "tau")[[k]]))
gc_ <- function(X, v) if (v %in% names(X)) X[[v]] else 0
TT <- list(fire = c("fire_f", "gdd_x_fire", "tsd_fire"), pest = c("pest_f", "gdd_x_pest", "tsd_pest"),
           light = c("harv_light", "gdd_x_harv_light", "tsd_harv_light"), heavy = c("harv_heavy", "gdd_x_harv_heavy", "tsd_harv_heavy"))
simulate2 <- function(b, ver, fm, seed) { set.seed(seed); zg <- B$zg[b]; zc <- B$zc[b]; zs <- B$zs[b]; CO <- COS[[fm]]; ft <- FT[[fm]]
  R <- if (ver == "none") matrix(0, nS, 4, dimnames = list(S, KD)) else rate_tab(ver, b); P <- 1 - exp(-R)
  n <- NPS * nS; st <- rep(seq_len(nS), each = NPS); age <- numeric(n)
  fire <- logical(n); pest <- logical(n); harv <- integer(n); tf <- tp <- th <- numeric(n); occ <- matrix(0, length(KEEP), nS)
  for (yr in seq_len(TMAX)) {
    if (ver != "none") { U <- matrix(runif(n * 4), n, 4); hit <- U < P[st, , drop = FALSE]
      i <- hit[, "pest"];  pest[i] <- TRUE; tp[i] <- 0
      i <- hit[, "light"]; harv[i] <- 1L; th[i] <- 0
      i <- hit[, "fire"];  fire[i] <- TRUE; tf[i] <- 0; pest[i] <- FALSE; harv[i] <- 0L
      i <- hit[, "heavy"]; harv[i] <- 2L; th[i] <- 0; fire[i] <- FALSE; pest[i] <- FALSE }
    newst <- st
    for (s in seq_len(nS)) { ii <- which(st == s); if (!length(ii) || !length(EX[[s]]$lb)) next
      H <- sapply(CO[[s]], function(X) { e <- X[["b0y"]] + gc_(X, "z_gdd_3") * zg + gc_(X, "z_cmi_gs") * zc + gc_(X, "z_soil_pc1") * zs
        add <- function(on, k, tsd) ifelse(on, gc_(X, TT[[k]][1]) + gc_(X, TT[[k]][2]) * zg + gc_(X, TT[[k]][3]) * ft(tsd + 0.5, k, X), 0)
        e <- e + add(fire[ii], "fire", tf[ii]) + add(pest[ii], "pest", tp[ii]) + add(harv[ii] == 1L, "light", th[ii]) + add(harv[ii] == 2L, "heavy", th[ii])
        exp(e) * ((age[ii] + 1)^X[["alpha"]] - age[ii]^X[["alpha"]]) })
      H <- matrix(H, nrow = length(ii)); tot <- rowSums(H); mv <- runif(length(ii)) < 1 - exp(-tot)
      if (any(mv)) { cs <- H[mv, , drop = FALSE] / tot[mv]
        k <- pmin(1L + rowSums(t(apply(cs, 1, cumsum)) < runif(sum(mv))), ncol(cs)); newst[ii[mv]] <- EX[[s]]$to[k] } }
    moved <- newst != st; st <- newst; age <- ifelse(moved, 0, age + 1); tf <- tf + 1; tp <- tp + 1; th <- th + 1
    if (yr %in% KEEP) occ[yr - KEEP[1] + 1, ] <- tabulate(st, nS) / n }
  stopifnot(!anyNA(st), abs(rowSums(occ) - 1) < 1e-9); data.table(bin = B$bin[b], version = ver, form = fm, state = S, share = colMeans(occ)) }
RUNS <- CJ(b = seq_len(nrow(B)), ver = c("none", "SC"), fm = if (TEST) "exp" else c("log", "exp"))
if (TEST) RUNS <- RUNS[b == 3]
for (r in seq_len(nrow(RUNS))) { f <- file.path(CI3, sprintf("b%d_%s_%s.rds", RUNS$b[r], RUNS$ver[r], RUNS$fm[r])); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > as.numeric(Sys.getenv("I13_BUDGET", "1e9"))) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t1 <- Sys.time(); saveRDS(simulate2(RUNS$b[r], RUNS$ver[r], RUNS$fm[r], 2600 + 10 * RUNS$b[r] + (RUNS$ver[r] == "SC")), f)
  cat(sprintf("bin %d %s %s (%.0f s)\n", RUNS$b[r], RUNS$ver[r], RUNS$fm[r], as.numeric(Sys.time() - t1, units = "secs"))) }
if (TEST) { for (v in c("none", "SC")) { a <- readRDS(file.path(CI3, sprintf("b3_%s_exp.rds", v))); b <- readRDS(file.path(OUTI, "cache_I13c", sprintf("b3_%s_exp.rds", v)))
  cat(sprintf("self-test %s: max |diff| vs 25_I13c = %.2e\n", v, max(abs(a$share - b$share)))) }; quit(save = "no") }
# ---- summary: shift = regime (SC) minus undisturbed (none), group shares, percentage points ----
source(here::here("R", "25_I_common.R"))
rdc <- function(dir, lab) { x <- rbindlist(lapply(list.files(file.path(OUT, dir), full.names = TRUE), readRDS)); x[, form := paste0(form, lab)][] }
SM <- rbind(rdc("cache_I13e", "_B"), rdc("cache_I13c", "_A")[form == "exp_A"])
G <- dcast(SM[, .(share = sum(share)), by = .(bin, grp = GRP[state], version, form)], bin + grp + form ~ version, value.var = "share")
G[, shift := SC - none]; fwrite(G, file.path(OUT, "I13e_tsd_shape_B_group_shares.csv"))
G[, form := factor(form, c("log_B", "exp_A", "exp_B"), c("log", "exp_A", "exp_B"))]
for (g in c("boreal", "pioneer", "temperate", "regen")) { cat("\n====", g, ": regime share (%) | shift vs undisturbed (points) ====\n")
  print(dcast(G[grp == g], bin ~ form, value.var = c("SC", "shift"), fun.aggregate = function(x) round(100 * x))) }
