# ==============================================================================
# 21_figures_D4_D7.R   (seconds; after 16 and 17 have run; from the project root). Exploratory set D, parts 4-7. NOT manuscript figures.
# D4 expected composition over time starting from the observed first-visit composition (plot_info first_state; time in state 0 for all, which
#    ignores how long plots had already been in their state): mixture of 17's occupancy runs (undisturbed, GDD -1/0/+1). Point estimates.
# D5 Shannon entropy H(t) = -sum_j P(j, t) log P(j, t) of the expected composition: per start state (17, GDD 0) and for the D4 mixture.
# D6 network of first-move flows within 50 y (16: CIF to each destination at 50 y, undisturbed), edges >= 5%, GDD -1/0/+1.
# D7 HRs in raw units from the pooled coefficients: per 100 GDD degree-days, per 1 unit of cmi_gs, per 1 unit of soil PC1; plus an APPROXIMATE
#    per-degC version: GDD change per +1 degC of growing-season mean temperature (an_meanT_3) from a linear regression across visit intervals.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R")); OUTD <- file.path(FIG_DIR, "D"); OUT <- OUTD; olab <- function(x) gsub("_", " ", x)
OCC <- rbindlist(lapply(list.files(file.path(FIG_DIR, "C", "cache"), "^occ_none_z.*rds$", full.names = TRUE), readRDS))
pinfo <- readRDS(file.path(TABLE_DIR, "plot_info.rds")); w0 <- pinfo[, .N, by = .(start = first_state)][, w := N / sum(N)]
# ---- D4 ----  (a0 = years already spent in the start state at t = 0: 0 (C runs) and 40 (computed here, point estimates; user 2026-11-09))
src17 <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src17[1:(grep("^# ---- propagation engine", src17) - 1)])); OUT <- OUTD   # the eval re-sources the helpers (OUT = figures/A)
occ_a0 <- function(zg, a0, tmax = 100) { nA <- tmax + a0 + 1L; ages <- 0:(nA - 1)
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][1, , drop = FALSE]; a <- X[, "alpha"]
    list(eta = X[, "b0y"] + gcol(X, "z_gdd_3") * zg, inc = (ages + 1)^a - ages^a) })); names(TR) <- S
  M <- array(0, c(nS, nS, nA)); for (i in seq_len(nS)) M[i, i, a0 + 1] <- 1; out <- list()
  for (st in seq_len(tmax)) { New <- array(0, dim(M))
    for (s in seq_len(nS)) { ex <- EX[[s]]; Ms <- matrix(M[, s, ], nS, nA)
      if (!length(ex$lb)) { New[, s, 2:nA] <- New[, s, 2:nA] + Ms[, 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta) * tr$inc); tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      New[, s, 2:nA] <- New[, s, 2:nA] + sweep(Ms, 2, 1 - pl, "*")[, 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; New[, ex$to[k], 1] <- New[, ex$to[k], 1] + drop(Ms %*% pk) } }
    M <- New; out[[st]] <- data.table(start = rep(S, nS), state = rep(S, each = nS), t = st, est = as.vector(apply(M, 1:2, sum))) }
  rbindlist(out)[, zg := zg] }
OCC40 <- rbindlist(lapply(c(-1, 0, 1), occ_a0, a0 = 40))
chk0 <- occ_a0(0, 0)[t == 50]; cmp0 <- merge(chk0, OCC[zg == 0 & t == 50, .(start, state, est0 = est)], by = c("start", "state"))
cat(sprintf("D4 engine check (a0 = 0 vs 17): max |diff| = %.2e\n", max(abs(cmp0$est - cmp0$est0)))); stopifnot(max(abs(cmp0$est - cmp0$est0)) < 1e-8)
mixf <- function(O, a0) { m <- merge(O, w0[, .(start, w)], by = "start")[, .(est = sum(est * w)), by = .(state, t, zg)]
  rbind(m, data.table(state = w0$start, t = 0, est = w0$w)[, .(zz = c(-1, 0, 1)), by = .(state, t, est)][, .(state, t, zg = zz, est)])[, a0 := a0] }
mix0 <- rbind(mixf(OCC[, .(start, state, t, est, zg)], 0), mixf(OCC40, 40)); mix <- mix0[t > 0]
mix0[, `:=`(state = factor(state, levels = ORDER), gdd = factor(gsd(zg), levels = gsd(c(-1, 0, 1))), a0l = factor(sprintf("%d y in state at t = 0", a0), levels = c("0 y in state at t = 0", "40 y in state at t = 0")))]
fwrite(mix0, file.path(OUTD, "D4_composition_from_first_visits.csv"))
p <- ggplot(mix0, aes(t, est, fill = state)) + geom_area(colour = "white", linewidth = 0.1) + geom_vline(xintercept = 55, linetype = 2, linewidth = 0.4) +
  scale_fill_manual(values = PAL, labels = olab, name = NULL) + scale_y_continuous(labels = function(v) paste0(100 * v, "%")) + facet_grid(a0l ~ gdd) +
  labs(x = "years from the first-visit composition", y = "expected share", title = "D4  Expected composition over time from the observed first-visit composition (undisturbed)",
       subtitle = "all plots set to the same GDD; not a projection of the real plots; dashed = end of record (55 y)")
save_fig(p, "D4_composition_from_first_visits", 15, 9, NULL)
# ---- D5 ----
ent <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
E1 <- OCC[zg == 0, .(H = ent(est)), by = .(start, t)][, start := factor(start, levels = ORDER)]
E2 <- mix0[, .(H = ent(est)), by = .(t, zg, a0)]
fwrite(E1, file.path(OUTD, "D5_entropy_by_start.csv")); fwrite(E2, file.path(OUTD, "D5_entropy_mixture.csv"))
p1 <- ggplot(E1, aes(t, H, colour = start)) + geom_line() + scale_colour_manual(values = PAL, labels = olab, name = "start state") + geom_hline(yintercept = log(13), linetype = 3) +
  labs(x = "years since start", y = "entropy of the expected composition", subtitle = "by start state (undisturbed, GDD 0); dotted = maximum")
p2 <- ggplot(E2, aes(t, H, colour = factor(zg), linetype = factor(a0))) + geom_line(linewidth = 0.9) + scale_linetype_manual(values = c(`0` = "solid", `40` = "dashed"), labels = c("0 y in state at t = 0", "40 y"), name = NULL) + scale_colour_manual(values = c(`-1` = "#2166AC", `0` = "grey30", `1` = "#B2182B"), labels = c("GDD -1", "GDD 0", "GDD +1"), name = NULL) +
  geom_hline(yintercept = log(13), linetype = 3) + labs(x = "years from the first-visit composition", y = "entropy", subtitle = "from the first-visit composition (D4)")
save_fig(p1 + p2 + plot_annotation(title = "D5  Entropy of the expected composition", subtitle = "how evenly probability spreads over the 13 states; not within-stand species diversity"),
         "D5_entropy", 14, 5.5, NULL)
# ---- D6 network ----
B <- fread(file.path(FIG_DIR, "B", "B_quantities.csv"))[scen == "none" & quantity == "cif50" & zg %in% c(-1, 0, 1), .(from = origin, to = dest, zg, w = est)]
ang <- seq(90, 90 - 360 + 360 / 13, length.out = 13) * pi / 180; nodes <- data.table(state = ORDER, x = cos(ang), y = sin(ang))
E <- merge(merge(B[w >= 0.05], nodes[, .(from = state, x, y)], by = "from"), nodes[, .(to = state, xend = x, yend = y)], by = "to")
E[, `:=`(x2 = x + 0.88 * (xend - x), y2 = y + 0.88 * (yend - y), x1 = x + 0.12 * (xend - x), y1 = y + 0.12 * (yend - y), gdd = factor(gsd(zg), levels = gsd(c(-1, 0, 1))))]
nn <- nodes[, .(gdd = factor(gsd(c(-1, 0, 1)), levels = gsd(c(-1, 0, 1)))), by = .(state, x, y)]
p <- ggplot() + geom_curve(data = E, aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = w, colour = from), curvature = 0.2, alpha = 0.75,
                           arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_point(data = nn, aes(x, y, fill = state), shape = 21, size = 7) + geom_text(data = nn, aes(x * 1.2, y * 1.2, label = olab(state)), size = 2.6) +
  scale_linewidth(range = c(0.3, 3), labels = function(v) paste0(round(100 * v), "%"), name = "P(moved there\nwithin 50 y)") +
  scale_fill_manual(values = PAL, guide = "none") + scale_colour_manual(values = PAL, guide = "none") + coord_equal(clip = "off") + facet_wrap(~ gdd, nrow = 1) +
  theme_void() + theme(strip.text = element_text(size = 10), plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(10, 10, 10, 10)) + labs(title = "D6  Network of first moves within 50 y (undisturbed)",
  subtitle = "edge width = probability of the move; edges ≥ 5% shown; colour = origin")
save_fig(p, "D6_network", 16, 6.5, E[, .(from, to, zg, w)])
# ---- D7 raw units ----
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); sdv <- setNames(sc$sd, sc$var)
fit <- lm(gdd_3 ~ an_meanT_3, data = vi); dd_per_C <- coef(fit)[[2]]; r2 <- summary(fit)$r.squared
UN <- data.table(parameter = c("z_gdd_3", "z_gdd_3", "z_cmi_gs", "z_soil_pc1"), unit = c("per 100 degree-days GDD", sprintf("per +1 degC (approx.; %.0f degree-days per degC, R2 = %.2f)", dd_per_C, r2), "per 1 unit cmi_gs", "per 1 unit soil PC1"),
                 mult = c(100 / sdv[["gdd_3"]], dd_per_C / sdv[["gdd_3"]], 1 / sdv[["cmi_gs"]], 1 / sdv[["soil_pc1"]]))
R <- merge(pooled[parameter %in% UN$parameter, .(label, from, to, parameter, b = pooled_est, lo = ci_lo, hi = ci_hi)], UN, by = "parameter", allow.cartesian = TRUE)
R[, `:=`(hr = exp(b * mult), hr_lo = exp(lo * mult), hr_hi = exp(hi * mult))]
fwrite(R[, .(label, parameter, unit, hr, hr_lo, hr_hi)], file.path(OUTD, "D7_hazard_ratios_raw_units.csv"))
cat(sprintf("D7: GDD SD = %.0f degree-days; CMI SD = %.2f; soil SD = %.2f; %.0f degree-days per degC (R2 %.2f)\n", sdv[["gdd_3"]], sdv[["cmi_gs"]], sdv[["soil_pc1"]], dd_per_C, r2))
x <- R[unit %like% "per 100"]; g <- merge(full_grid(), x, by = c("from", "to"), all.x = TRUE)
g[, status := fifelse(paste0(from, "->", to) %in% modelled, fifelse(is.na(hr), "noterm", "value"), "notmod")]
g[, `:=`(fillv = log(hr), txt = fmt_hr(hr), bold = !is.na(hr_lo) & (hr_lo > 1 | hr_hi < 1))]
lim <- cap_sym(log(x$hr))
save_fig(heat(g, "hazard\nratio", limits = c(-lim, lim), trans_lab = function(v) fmt_hr(exp(v)), breaks = hr_breaks(lim), title = "D7  GDD hazard ratio per 100 degree-days",
              subtitle = "bold = 95% CI excludes 1"),
         "D7_gdd_hr_per_100dd", 7.5, 7, NULL)
cat("D4-D7 done\n")
