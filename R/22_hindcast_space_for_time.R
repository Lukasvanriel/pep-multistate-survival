# ==============================================================================
# 22_hindcast_space_for_time.R   (no INLA; resumable: re-run until it prints "all draws done"; run from the project root)
# Space-for-time hindcast. The GDD/CMI effects are identified mostly between plots (95% of GDD variance is between plots). Within plots, GDD rose
#   (median +0.49 SD from first to last visit). Question: do the spatially estimated climate effects, applied to the warming each plot actually
#   experienced, predict the change seen in the observed intervals?
# Per observed visit interval (start state, time in state at the start, own covariates of the interval as in the fit: end-visit climate, flags, TSD):
#   predicted end-state distribution over the interval length (first exit on an M-step grid with exact Weibull increments; after a move, the new
#   state's exits from time-in-state 0 with the same interval covariates; a third move is spread by the cumulative-hazard shares of the second state).
# Time in state at the start: first sojourn entered at W_j / 2 (the stage-1 entry rule of 04b); later sojourns entered at the midpoint of the change bracket.
# Scenarios: A = actual climate of the interval; FG = GDD frozen at the plot's first visit (CMI actual); F = GDD and CMI frozen at the first visit.
#   Disturbances, soil and time in state are identical in all three, so A - F is the change the model attributes to the plot's own climate change.
# Test: per bin of within-plot warming dose (and of years since first visit, calendar period, origin group), compare observed - F with A - F.
#   Spatial effects transfer to time if observed - F tracks A - F (ratio ~ 1); ratio ~ 0 = no temporal signal; calibration offset read from the lowest-dose bin.
# Outcomes per interval: change (end != start), dSTI (STI as in 20: mean z-GDD of the visits in each state), end in temperate group, end in boreal group.
# Uncertainty: parameters = point (row 1) + draws 2..(N_DR+1) of 17's PAR; observation noise = plot bootstrap; draw r is paired with bootstrap r.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
OUTH <- file.path(FIG_DIR, "H"); CH <- file.path(OUTH, "cache"); dir.create(CH, recursive = TRUE, showWarnings = FALSE); OUT <- OUTH
M <- 20L; N_DR <- 40L; TIME_BUDGET <- 150; set.seed(22)
source(here::here("R", "00_groups.R"))   # GRP: thermal groups
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
STI <- vi[, .(STI = mean(zf(gdd_3, "gdd_3"))), by = .(state = state_end)]; STI <- setNames(STI$STI, STI$state)[S]

# ---- intervals -----------------------------------------------------------------------------------------------------------------------------------------------
iv <- copy(vi); setorder(iv, id_pe, interval)
pinf <- readRDS(file.path(TABLE_DIR, "plot_info.rds"))
iv[pinf, `:=`(fg = i.fv_gdd_3, fc = i.fv_cmi_gs, yr1 = i.first_year), on = "id_pe"]
iv[, a0 := { E <- W_j[1L] / 2; a <- numeric(.N); for (q in seq_len(.N)) { a[q] <- t_start[q] - E; if (is_change[q]) E <- (t_start[q] + t_end[q]) / 2 }; a }, by = id_pe]
iv[, `:=`(zg = zf(gdd_3, "gdd_3"), zc = zf(cmi_gs, "cmi_gs"), zs = zf(soil_pc1, "soil_pc1"), zg0 = zf(fg, "gdd_3"), zc0 = zf(fc, "cmi_gs"),
          ff = as.integer(fire == 1L), pf = as.integer(pest == 1L), hl = as.integer(harv == 1L), hh = as.integer(harv == 2L), L = t_end - t_start)]
stopifnot(all(iv$state_start %in% S), all(iv$state_end %in% S), !anyNA(iv[, .(a0, zg, zc, zs, zg0, zc0, L)]), all(iv$a0 >= 0), all(iv$L > 0))
iv[, `:=`(dz = zg - zg0, ysf = year_end - yr1, grp0 = GRP[state_start],
          o_chg = as.numeric(is_change), o_dsti = STI[state_end] - STI[state_start],
          o_tmp = as.numeric(GRP[state_end] == "temperate"), o_bor = as.numeric(GRP[state_end] == "boreal"))]
iv[, `:=`(b_dose = cut(dz, c(-Inf, 0.1, 0.3, 0.5, Inf), labels = c("<0.1", "0.1-0.3", "0.3-0.5", ">0.5")),
          b_ysf  = cut(ysf, c(-Inf, 10, 20, 30, Inf), labels = c("0-10", "11-20", "21-30", ">30")),
          b_cal  = cut(year_end, c(-Inf, 1990, 2000, 2010, Inf), labels = c("<=1990", "1991-2000", "2001-2010", ">2010")))]

# ---- engine ----------------------------------------------------------------------------------------------------------------------------------------------------
eta_of <- function(X, d, g, cc) { x <- function(v) if (v %in% names(X)) X[[v]] else 0; ft <- function(y) f_tsd(y)
  x("b0y") + x("z_gdd_3") * g + x("z_cmi_gs") * cc + x("z_soil_pc1") * d$zs +
    d$ff * (x("fire_f") + x("gdd_x_fire") * g + x("tsd_fire") * ft(d$tsd_fire_y)) +
    d$pf * (x("pest_f") + x("gdd_x_pest") * g + x("tsd_pest") * ft(d$tsd_pest_y)) +
    d$hl * (x("harv_light") + x("gdd_x_harv_light") * g + x("tsd_harv_light") * ft(d$tsd_harv_y)) +
    d$hh * (x("harv_heavy") + x("gdd_x_harv_heavy") * g + x("tsd_harv_heavy") * ft(d$tsd_harv_y)) }
pred_end <- function(s, d, r, g, cc) { n <- nrow(d); P <- matrix(0, n, nS); si <- match(s, S); ex <- EX[[s]]
  if (!length(ex$lb)) { P[, si] <- 1; return(P) }
  ed <- outer(d$L, (0:M) / M); um <- (ed[, -1, drop = FALSE] + ed[, -(M + 1), drop = FALSE]) / 2; R <- d$L - um
  dH <- lapply(ex$lb, function(lb) { X <- PAR[[lb]][r, ]; A <- (d$a0 + ed)^X[["alpha"]]; exp(eta_of(X, d, g, cc)) * (A[, -1, drop = FALSE] - A[, -(M + 1), drop = FALSE]) })
  tot <- Reduce(`+`, dH); Hc <- t(apply(tot, 1, cumsum)); Sst <- exp(-cbind(0, Hc[, -M, drop = FALSE])); Sen <- exp(-Hc); lv <- Sst - Sen
  P[, si] <- Sen[, M]
  for (k in seq_along(ex$lb)) { j <- ex$to[k]; mk <- dH[[k]] / pmax(tot, 1e-300) * lv; exj <- EX[[S[j]]]
    if (!length(exj$lb)) { P[, j] <- P[, j] + rowSums(mk); next }
    Hj <- lapply(exj$lb, function(lb) { X <- PAR[[lb]][r, ]; exp(eta_of(X, d, g, cc)) * R^X[["alpha"]] })
    Ht <- Reduce(`+`, Hj); Sj <- exp(-Ht); P[, j] <- P[, j] + rowSums(mk * Sj)
    for (q in seq_along(exj$lb)) P[, exj$to[q]] <- P[, exj$to[q]] + rowSums(mk * (1 - Sj) * Hj[[q]] / pmax(Ht, 1e-300)) }
  P }
SCEN <- list(A = c("zg", "zc"), FG = c("zg0", "zc"), F = c("zg0", "zc0"))
metrics <- function(P, s0) data.table(chg = 1 - P[cbind(seq_len(nrow(P)), match(s0, S))], dsti = drop(P %*% STI) - STI[s0],
                                      tmp = rowSums(P[, GRP[S] == "temperate", drop = FALSE]), bor = rowSums(P[, GRP[S] == "boreal", drop = FALSE]))
predict_all <- function(r) { res <- vector("list", nS); for (s in S) { ii <- which(iv$state_start == s); if (!length(ii)) next; d <- iv[ii]
    out <- data.table(row = ii)
    for (sn in names(SCEN)) { P <- pred_end(s, d, r, d[[SCEN[[sn]][1]]], d[[SCEN[[sn]][2]]])
      if (r == 1L) stopifnot(max(abs(rowSums(P) - 1)) < 1e-8)
      m <- metrics(P, rep(s, nrow(d))); setnames(m, paste0(sn, "_", names(m))); out <- cbind(out, m) }
    res[[match(s, S)]] <- out }
  setorder(rbindlist(res), row)[] }
BINS <- c("b_dose", "b_ysf", "b_cal")
agg <- function(pr, w) { d <- cbind(iv[, c(BINS, "grp0", "o_chg", "o_dsti", "o_tmp", "o_bor"), with = FALSE], pr[, -1], w = w)
  vs <- setdiff(names(d), c(BINS, "grp0", "w"))
  f <- function(by) d[, c(list(n = sum(w)), lapply(.SD, function(v) sum(v * w))), by = by, .SDcols = vs]
  rbindlist(c(lapply(BINS, function(b) { x <- f(b); setnames(x, b, "bin"); x[, `:=`(binning = b, grp0 = "all")] }),
              lapply(BINS, function(b) { x <- f(c(b, "grp0")); setnames(x, b, "bin"); x[, binning := b] }),
              list(f(NULL)[, `:=`(bin = "all", binning = "all", grp0 = "all")])), use.names = TRUE) }
plots <- unique(iv$id_pe); pidx <- split(seq_len(nrow(iv)), iv$id_pe)
for (r in seq_len(N_DR + 1L)) { f <- file.path(CH, sprintf("H_r%03d.rds", r)); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t_start, units = "secs") > TIME_BUDGET) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t0 <- Sys.time(); pr <- predict_all(r)
  if (r == 1L) { w <- rep(1, nrow(iv)); saveRDS(cbind(iv[, .(id_pe, interval, state_start, state_end, year_end, a0, L, dz, ysf, grp0, o_chg, o_dsti, o_tmp, o_bor)], pr[, -1]), file.path(OUTH, "H_interval_point.rds"))
  } else { bs <- sample(plots, length(plots), replace = TRUE); w <- tabulate(unlist(pidx[bs], use.names = FALSE), nbins = nrow(iv)) }
  saveRDS(agg(pr, w)[, draw := r], f); cat(sprintf("draw %d (%.0f s)\n", r, as.numeric(Sys.time() - t0, units = "secs"))) }
cat("all draws done\n")
