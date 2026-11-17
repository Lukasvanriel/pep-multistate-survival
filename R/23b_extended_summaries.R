# ==============================================================================
# 23b_extended_summaries.R   (no INLA; ~1 min; after 23_extended_quantities.R; run from the project root)
# Tables E1-E7 from figures/E/cache (definitions in the header of 23). CI = 2.5/97.5 pct over the 100 draws; point = pooled estimates.
# Mixture = start states weighted by the observed first-visit composition (plot_info), time in state 0 at t = 0.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R")); OUTE <- file.path(FIG_DIR, "E"); CE <- file.path(OUTE, "cache"); OUT <- OUTE
m <- readRDS(file.path(CE, "meta.rds")); S <- m$S; nS <- length(S); STI <- m$STI; TOL <- m$TOL; TGRP <- m$TGRP; W0 <- m$W0
source(here::here("R", "00_groups.R")); GRP <- GRP[S]   # GRP: thermal groups, ordered as S
SCN <- c("none_zm1", "none_z0", "none_zp1", "fire", "pest", "harv_light", "harv_heavy")
SL <- c(none_zm1 = "undisturbed, GDD -1 SD", none_z0 = "undisturbed, GDD 0", none_zp1 = "undisturbed, GDD +1 SD", fire = "fire", pest = "pest",
        harv_light = "light harvest", harv_heavy = "heavy harvest")
PERT <- list(gdd_p1 = c("none_zp1", "none_z0"), gdd_m1 = c("none_zm1", "none_z0"), fire = c("fire", "none_z0"), pest = c("pest", "none_z0"),
             harv_light = c("harv_light", "none_z0"), harv_heavy = c("harv_heavy", "none_z0"))
qq <- function(v) list(est = v[1], lo = unname(quantile(v[-1], 0.025, na.rm = TRUE)), hi = unname(quantile(v[-1], 0.975, na.rm = TRUE)))
base <- setNames(lapply(SCN, function(s) readRDS(file.path(CE, paste0("base_", s, ".rds")))), SCN)

# E1 -------------------------------------------------------------------------------------------------------------------------------------------------------------
e1 <- readRDS(file.path(CE, "E1_E2_first_move.rds")); E1 <- e1$E1; fwrite(E1, file.path(OUTE, "E1_next_state_by_exit_time.csv"))
E1s <- E1[!grepl("^<", dest) & scen == "none_z0" & t %in% c(5, 25, 50)]
E1chg <- E1s[, .(top_t5 = dest[t == 5][which.max(est[t == 5])], share_t5 = max(est[t == 5]), top_t50 = dest[t == 50][which.max(est[t == 50])], share_t50 = max(est[t == 50]),
                 max_abs_change_5_50 = max(abs(est[t == 50] - est[t == 5]))), by = origin]
E1sti <- dcast(E1[dest == "<expected dSTI>" & scen == "none_z0" & t %in% c(5, 50)], origin ~ t, value.var = c("est", "lo", "hi"))
E1chg <- merge(E1chg, E1sti, by = "origin"); fwrite(E1chg, file.path(OUTE, "E1_summary_t5_vs_t50.csv"))

# E2 -------------------------------------------------------------------------------------------------------------------------------------------------------------
FM <- e1$FM; E2 <- list()
for (pn in names(PERT)) { P <- PERT[[pn]][1]; B <- PERT[[pn]][2]
  mx <- list(PB = readRDS(file.path(CE, paste0("mixHP_", pn, ".rds")))$TI[, , 50], BP = readRDS(file.path(CE, paste0("mixJP_", pn, ".rds")))$TI[, , 50])
  PP <- base[[P]]$TI[, , 50]; BB <- base[[B]]$TI[, , 50]
  for (i in seq_len(nS + 1)) { if (i <= nS) { pp <- PP[, i]; bb <- BB[, i]; pb <- mx$PB[, i]; bp <- mx$BP[, i]; who <- S[i] } else {
      pp <- PP %*% W0; bb <- BB %*% W0; pb <- mx$PB %*% W0; bp <- mx$BP %*% W0; who <- "observed composition" }
    tot <- pp - bb; hold <- ((pb - bb) + (pp - bp)) / 2; jump <- tot - hold
    for (q in c("total", "holding", "jump")) E2[[length(E2) + 1]] <- as.data.table(c(list(outcome = "thermal index at 50 y (multi-step)", pert = pn, start = who, part = q), qq(get(c(total = "tot", holding = "hold", jump = "jump")[q])))) }
  for (o in S) { k <- function(h, j) FM[[paste(o, h, j)]]; if (is.null(k(P, P)) || is.null(k(B, B))) next
    pp <- k(P, P); bb <- k(B, B); pb <- k(P, B); bp <- k(B, P); tot <- pp - bb; hold <- ((pb - bb) + (pp - bp)) / 2; jump <- tot - hold
    for (q in c("total", "holding", "jump")) E2[[length(E2) + 1]] <- as.data.table(c(list(outcome = "first-move thermal shift (50 y)", pert = pn, start = o, part = q), qq(get(c(total = "tot", holding = "hold", jump = "jump")[q])))) } }
E2 <- rbindlist(E2); fwrite(E2, file.path(OUTE, "E2_holding_vs_jump.csv"))

# E3 -------------------------------------------------------------------------------------------------------------------------------------------------------------
b0 <- readRDS(file.path(CE, "elas_base.rds"))$MIX[1, , 50]; ti0 <- sum(b0 * STI); g0 <- tapply(b0, GRP, sum)
E3 <- rbindlist(lapply(modelled, function(lb) { r <- readRDS(file.path(CE, paste0("elas_", gsub("->", "__", lb), ".rds")))
  rbindlist(lapply(c("b0", "alpha"), function(ty) { up <- r[[paste(ty, 1.1)]]; dn <- r[[paste(ty, 0.9)]]; gu <- tapply(up, GRP, sum); gd <- tapply(dn, GRP, sum)
    rbind(data.table(label = lb, param = ty, outcome = "thermal index", y0 = ti0, d_plus = sum(up * STI) - ti0, d_minus = sum(dn * STI) - ti0),
          data.table(label = lb, param = ty, outcome = paste("share", names(g0)), y0 = as.numeric(g0), d_plus = as.numeric(gu - g0), d_minus = as.numeric(gd - g0))) })) }))
E3[, `:=`(central_per_10pct = (d_plus - d_minus) / 2, asym = abs(d_plus + d_minus))]
E3[grepl("^share", outcome), elasticity := (d_plus - d_minus) / (0.2 * y0)]
fwrite(E3, file.path(OUTE, "E3_elasticity_50y.csv"))

# E4 -------------------------------------------------------------------------------------------------------------------------------------------------------------
E4 <- rbindlist(lapply(SCN, function(s) { y <- base[[s]]$yrs; rbindlist(lapply(1:2, function(h) { g <- CJ(i = seq_len(nS), j = seq_len(nS))
  rbindlist(lapply(seq_len(nrow(g)), function(r) as.data.table(c(list(scen = s, horizon = c(50, 100)[h], start = S[g$i[r]], state = S[g$j[r]]), qq(y[, g$i[r], g$j[r], h]))))) })) }))
fwrite(E4, file.path(OUTE, "E4_expected_years_in_state.csv"))

# E5 -------------------------------------------------------------------------------------------------------------------------------------------------------------
E5 <- rbindlist(lapply(SCN, function(s) { C <- base[[s]]$CNT; R <- readRDS(file.path(CE, paste0("ret_", s, ".rds")))$RET; cum <- aperm(apply(R, c(1, 2), cumsum), c(2, 3, 1))
  rbindlist(lapply(seq_len(nS + 1), function(i) {
    if (i <= nS) { n50 <- apply(C[, i, , , 1], 1, sum); n100 <- apply(C[, i, , , 2], 1, sum); r50 <- cum[, i, 50]; r100 <- cum[, i, 100]
      med <- apply(cum[, i, ], 1, function(v) { if (v[100] < 1e-6) NA_real_ else which(v >= v[100] / 2)[1] }); who <- S[i] } else {
      n50 <- apply(C[, , , , 1], 1, function(a) sum(apply(a, 1, sum) * W0)); n100 <- apply(C[, , , , 2], 1, function(a) sum(apply(a, 1, sum) * W0))
      r50 <- cum[, , 50] %*% W0; r100 <- cum[, , 100] %*% W0; med <- rep(NA_real_, length(n50)); who <- "observed composition" }
    rbind(as.data.table(c(list(scen = s, start = who, quantity = "state changes within 50 y"), qq(n50))), as.data.table(c(list(scen = s, start = who, quantity = "state changes within 100 y"), qq(n100))),
          as.data.table(c(list(scen = s, start = who, quantity = "P(return to start within 50 y)"), qq(r50))), as.data.table(c(list(scen = s, start = who, quantity = "P(return to start within 100 y)"), qq(r100))),
          as.data.table(c(list(scen = s, start = who, quantity = "median return time among returners (y)"), qq(med)))) })) }))
fwrite(E5, file.path(OUTE, "E5_turnover_and_return.csv"))

# E6 -------------------------------------------------------------------------------------------------------------------------------------------------------------
E6 <- rbindlist(lapply(c("none_zm1", "none_z0", "none_zp1"), function(s) { mx <- readRDS(file.path(CE, paste0("conv_", s, ".rds"))); ref <- mx[, 1000]
  tv <- c(0.5 * sum(abs(W0 - ref)), apply(mx, 2, function(p) 0.5 * sum(abs(p - ref)))); data.table(scen = s, t = 0:1000, tvd = tv) }))
fwrite(E6, file.path(OUTE, "E6_distance_to_long_run.csv"))
E6s <- E6[, .(tvd_t0 = tvd[t == 0], tvd_50 = tvd[t == 50], tvd_100 = tvd[t == 100], years_to_tvd_0.10 = t[which(tvd < 0.10)[1]], years_to_tvd_0.05 = t[which(tvd < 0.05)[1]],
              change_900_1000 = tvd[t == 900]), by = scen]; fwrite(E6s, file.path(OUTE, "E6_summary.csv"))

# E7 -------------------------------------------------------------------------------------------------------------------------------------------------------------
E7 <- list(); E7g <- list()
for (s in SCN) { C <- base[[s]]$CNT
  for (h in 1:2) { M <- apply(C[, , , , h, drop = FALSE], c(1, 3, 4), function(v) sum(v * W0))     # [draw, from, to] expected moves per stand (observed composition)
    dT <- outer(TOL, TOL, function(a, b) b - a); sc <- !is.na(dT); diag(sc) <- FALSE
    tot <- apply(M, 1, sum); scm <- apply(M, 1, function(x) sum(x[sc])); mT <- apply(M, 1, function(x) sum(x[sc] * dT[sc])) / scm
    up <- apply(M, 1, function(x) sum(x[sc & dT > 0])); dn <- apply(M, 1, function(x) sum(x[sc & dT < 0])); toR <- apply(M, 1, function(x) sum(x[, S == "regen"]))
    for (q in list(list("moves per stand (all)", tot), list("share of moves scored for tolerance", scm / tot), list("mean tolerance shift of scored moves", mT),
                   list("scored moves to MORE tolerant", up), list("scored moves to LESS tolerant", dn), list("net (more - less tolerant)", up - dn), list("moves into regen", toR)))
      E7[[length(E7) + 1]] <- as.data.table(c(list(scen = s, horizon = c(50, 100)[h], quantity = q[[1]]), qq(q[[2]])))
    G <- unique(TGRP[!is.na(TGRP)]); for (a in G) for (b in G) if (a != b) { v <- apply(M, 1, function(x) sum(x[which(TGRP == a), which(TGRP == b)])); E7g[[length(E7g) + 1]] <- as.data.table(c(list(scen = s, horizon = c(50, 100)[h], from_g = a, to_g = b), qq(v))) } } }
E7 <- rbindlist(E7); E7g <- rbindlist(E7g); fwrite(E7, file.path(OUTE, "E7_successional_direction.csv")); fwrite(E7g, file.path(OUTE, "E7_tolerance_group_flows.csv"))
saveRDS(list(E1chg = E1chg, E2 = E2, E3 = E3, E4 = E4, E5 = E5, E6s = E6s, E7 = E7, E7g = E7g), file.path(OUTE, "E_summaries.rds")); cat("23b summaries done\n")
