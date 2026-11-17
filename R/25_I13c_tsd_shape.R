# ==============================================================================
# 25_I13c_tsd_shape.R   (no INLA; resumable; run from the project root; 2026-11-12)
# I13 under the exponential TSD shape (user 2026-11-12, option A): does the equilibrium under the disturbance regime depend on the TSD shape?
#   Same microsimulation, rates (version SC) and bins as 25_I13_disturbance_regime.R; only the coefficients of the 72 transitions with a TSD term
#   change. Both shapes come from the same sim06 fits (5 imputations per transition), so the comparison is not confounded by the imputation set:
#     log: sim06 pass 1 "log_all" (= main model form);  exp: sim06 pass 2 "best_all", every TSD term exp(-TSD/tau), one tau per type pooled over
#     transitions (fire 15, pest 15, light harvest 25, heavy harvest 160 y).
#   Point coefficients = mean over the imputations. Imputations failing the alpha-SD screen (ratio < 0.1) in either shape are dropped from both.
#   The 18 transitions without a TSD term keep the main pooled point estimate (identical in both). Same random numbers for log and exp.
#   Versions: none (no events) and SC. Also compared in 25_I13d: log (5 imputations) vs the main I13 (pooled over all imputations).
# ==============================================================================
src <- readLines(here::here("R", "25_I13_disturbance_regime.R")); eval(parse(text = src[1:(grep("^# ---- microsimulation", src) - 1)]))
CI3 <- file.path(OUTI, "cache_I13c"); dir.create(CI3, showWarnings = FALSE); t0 <- Sys.time()
TAU <- c(fire = 15, pest = 15, light = 25, heavy = 160)
D6 <- here::here("data-output", "simulation", "sim06_tsd_tau", "fits")
rd <- function(f, nm) { r <- readRDS(f); x <- r$fits[[nm]]; stopifnot(is.na(x$error)); data.table(label = r$label, m = r$m, ne = r$n_events, x$coef[, .(parameter, mean, sd)]) }
cf <- rbind(rbindlist(lapply(list.files(D6, "^p1_", full.names = TRUE), rd, "log_all"))[, form := "log"],
            rbindlist(lapply(list.files(D6, "^p2_fire15_harv_heavy160_harv_light25_pest15_", full.names = TRUE), rd, "best_all"))[, form := "exp"])
bad <- unique(cf[parameter == "alpha" & sd / (0.78 * mean / sqrt(ne)) < 0.1, .(label, m)]); cf <- cf[!bad, on = c("label", "m")]
cat("dropped (label, imputation):", paste(bad$label, bad$m, collapse = "; "), "\n")
K <- cf[, .(mean = mean(mean), n_imp = .N), by = .(form, label, parameter)]
stopifnot(K[, uniqueN(n_imp), by = .(label, parameter)]$V1 == 1, K[, uniqueN(form), by = label]$V1 == 2)
make_CO <- function(fm) { o <- lapply(S, function(s) lapply(EX[[s]]$lb, function(lb) { X <- PAR[[lb]][1, ]; k <- K[form == fm & label == lb]
  if (nrow(k)) { Y <- setNames(k$mean, k$parameter); stopifnot(setequal(names(Y), setdiff(names(X), "b0y"))); X <- c(Y, b0y = Y[["(Intercept)"]] - Y[["alpha"]] * log(10)) }; X }))
  names(o) <- S; o }
COS <- list(log = make_CO("log"), exp = make_CO("exp"))
FT <- list(log = function(x, k) log1p(x), exp = function(x, k) exp(-x / TAU[[k]]))
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
        add <- function(on, k, tsd) ifelse(on, gc_(X, TT[[k]][1]) + gc_(X, TT[[k]][2]) * zg + gc_(X, TT[[k]][3]) * ft(tsd + 0.5, k), 0)
        e <- e + add(fire[ii], "fire", tf[ii]) + add(pest[ii], "pest", tp[ii]) + add(harv[ii] == 1L, "light", th[ii]) + add(harv[ii] == 2L, "heavy", th[ii])
        exp(e) * ((age[ii] + 1)^X[["alpha"]] - age[ii]^X[["alpha"]]) })
      H <- matrix(H, nrow = length(ii)); tot <- rowSums(H); mv <- runif(length(ii)) < 1 - exp(-tot)
      if (any(mv)) { cs <- H[mv, , drop = FALSE] / tot[mv]
        k <- pmin(1L + rowSums(t(apply(cs, 1, cumsum)) < runif(sum(mv))), ncol(cs)); newst[ii[mv]] <- EX[[s]]$to[k] } }
    moved <- newst != st; st <- newst; age <- ifelse(moved, 0, age + 1); tf <- tf + 1; tp <- tp + 1; th <- th + 1
    if (yr %in% KEEP) occ[yr - KEEP[1] + 1, ] <- tabulate(st, nS) / n }
  stopifnot(!anyNA(st), abs(rowSums(occ) - 1) < 1e-9); data.table(bin = B$bin[b], version = ver, form = fm, state = S, share = colMeans(occ)) }
RUNS <- CJ(b = seq_len(nrow(B)), ver = c("none", "SC"), fm = c("log", "exp"))
for (r in seq_len(nrow(RUNS))) { f <- file.path(CI3, sprintf("b%d_%s_%s.rds", RUNS$b[r], RUNS$ver[r], RUNS$fm[r])); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > 100) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t1 <- Sys.time(); saveRDS(simulate2(RUNS$b[r], RUNS$ver[r], RUNS$fm[r], 2600 + 10 * RUNS$b[r] + (RUNS$ver[r] == "SC")), f)
  cat(sprintf("bin %d %s %s (%.0f s)\n", RUNS$b[r], RUNS$ver[r], RUNS$fm[r], as.numeric(Sys.time() - t1, units = "secs"))) }
cat("all I13c runs done\n")
