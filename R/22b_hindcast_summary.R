# ==============================================================================
# 22b_hindcast_summary.R   (no INLA; ~2-3 min; run after 22_hindcast_space_for_time.R; run from the project root)
# (1) Binned table (H_summary.csv): per bin of within-plot warming dose / years since first visit / calendar period: observed vs A and F; DiD vs the
#     lowest bin. Compares DIFFERENT plots across bins, so it mixes temporal transfer with between-plot misfit (dose correlates with the plot's GDD level).
# (2) Main estimator (H_within_plot_slopes_controls.csv): interval-level slope of (observed - F) on (A - F) WITHIN plots (plot fixed effects);
#     controls: start state, log(1 + time in state) and its square, interval length, previous interval changed (flicker); optionally calendar period.
#     Slope 1 = the spatially estimated climate effects predict the plot's own change over time; 0 = no temporal signal. CI = plot bootstrap (100),
#     point parameters only. Sojourn fixed effects were tried and dropped: a sojourn is defined by its own ending (outcome-conditioned).
# (3) Figure H1.
# ==============================================================================
source(here::here("R", "13_fig_helpers.R")); OUT <- file.path(FIG_DIR, "H"); set.seed(1)
# ---- (1) binned ----
D <- rbindlist(lapply(list.files(file.path(OUT, "cache"), full.names = TRUE), readRDS)); ms <- c("chg", "dsti", "tmp", "bor")
L <- rbindlist(lapply(ms, function(m) D[, .(draw, binning, bin, grp0, n, metric = m, obs = get(paste0("o_", m)) / n, A = get(paste0("A_", m)) / n, F = get(paste0("F_", m)) / n)]))
L[, `:=`(oF = obs - F, AF = A - F, oA = obs - A)]
ref <- L[bin %in% c("<0.1", "0-10", "<=1990"), .(draw, binning, grp0, metric, oF0 = oF, AF0 = AF)]
L <- merge(L, ref, by = c("draw", "binning", "grp0", "metric"), all.x = TRUE); L[, `:=`(did_o = oF - oF0, did_A = AF - AF0)]; setorder(L, draw)
s <- function(x) if (all(is.na(x))) NA_character_ else sprintf("%.4f [%.4f, %.4f]", x[1], quantile(x[-1], .025, na.rm = TRUE), quantile(x[-1], .975, na.rm = TRUE))
BT <- L[, .(n = n[1], obs = obs[1], A = A[1], F = F[1], obs_minus_A = s(oA), A_minus_F = s(AF), obs_minus_F = s(oF), did_obs = s(did_o), did_A = s(did_A)), by = .(metric, binning, grp0, bin)]
fwrite(BT, file.path(OUT, "H_summary.csv"))
# ---- (2) within-plot slopes ----
d <- readRDS(file.path(OUT, "H_interval_point.rds")); setorder(d, id_pe, interval)
d[, `:=`(per = cut(year_end, c(-Inf, 1990, 2000, 2010, Inf)), prev_chg = shift(o_chg, fill = 0)), by = id_pe]; d[, la := log1p(a0)]
fit <- function(x, m, pfe, gi) { y <- x[[paste0("o_", m)]] - x[[paste0("F_", m)]]; z <- x[[paste0("A_", m)]] - x[[paste0("F_", m)]]
  fo <- paste("~ la + I(la^2) + L + prev_chg", if (uniqueN(x$state_start) > 1) "+ state_start" else "", if (pfe) "+ per" else "")
  M <- cbind(y, z, model.matrix(as.formula(fo), x)[, -1, drop = FALSE]); M <- M - (rowsum(M, gi, reorder = TRUE) / tabulate(gi))[gi, , drop = FALSE]
  C <- M[, -(1:2), drop = FALSE]; C <- C[, colSums(C^2) > 1e-10, drop = FALSE]; Q <- qr(C); yy <- qr.resid(Q, M[, 1]); zz <- qr.resid(Q, M[, 2]); sum(yy * zz) / sum(zz^2) }
res <- rbindlist(lapply(c(FALSE, TRUE), function(pfe) rbindlist(lapply(ms, function(m) rbindlist(lapply(c("all", "boreal", "temperate", "pioneer", "regen"), function(gr) {
  x1 <- if (gr == "all") d else d[grp0 == gr]; x1 <- x1[x1[, .N, by = id_pe][N > 1], on = "id_pe"]; gi <- match(x1$id_pe, unique(x1$id_pe))
  est <- fit(x1, m, pfe, gi); ixg <- split(seq_len(nrow(x1)), gi)
  bt <- replicate(100, { b <- sample(length(ixg), replace = TRUE); ii <- unlist(ixg[b], use.names = FALSE); fit(x1[ii], m, pfe, rep(seq_along(b), lengths(ixg[b]))) })
  data.table(period_FE = pfe, metric = m, origin = gr, n = nrow(x1), slope = est, lo = quantile(bt, .025), hi = quantile(bt, .975)) }))))))
fwrite(res, file.path(OUT, "H_within_plot_slopes_controls.csv"))
# ---- (3) figure ----
ML <- c(chg = "any state change", dsti = "thermal shift (ΔSTI)", tmp = "end in temperate group", bor = "end in boreal group")
res[, `:=`(metric = factor(ML[metric], ML), origin = factor(origin, rev(c("all", "boreal", "pioneer", "regen", "temperate"))),
           version = factor(fifelse(period_FE, "+ calendar-period effects", "plot effects + controls"), c("plot effects + controls", "+ calendar-period effects")))]
p1 <- ggplot(res, aes(slope, origin, colour = version)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") + geom_vline(xintercept = 1, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", position = position_dodge(width = 0.6)) + geom_point(position = position_dodge(width = 0.6)) +
  facet_wrap(~ metric, nrow = 1) + scale_colour_manual(values = c("#1B6E3C", "#7570B3"), name = NULL) + labs(x = "within-plot slope, observed change on predicted change", y = "origin group",
  subtitle = "1 = spatial effects fully predict each plot's change as it warmed; 0 = no temporal signal")
bd <- L[draw == 1 & binning == "b_dose" & grp0 == "all", .(metric, bin, `observed − frozen` = oF, `actual − frozen (model)` = AF)]
bd <- melt(bd, id.vars = c("metric", "bin")); bd[, metric := factor(ML[metric], ML)]
p2 <- ggplot(bd, aes(bin, value, colour = variable, group = variable)) + geom_hline(yintercept = 0, colour = "grey60") + geom_line() + geom_point() +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") + scale_colour_manual(values = c("black", "#B2182B"), name = NULL) +
  labs(x = "GDD increase since first visit (SD)", y = "mean per interval", subtitle = "binned across plots: mixes temporal transfer with between-plot misfit")
save_fig(p1 / p2 + plot_annotation(title = "H1  Space-for-time hindcast: do spatial climate effects predict each plot's change as it warmed?",
         subtitle = "frozen = climate held at the plot's first visit; everything else as observed"), "H1_space_for_time_hindcast", 15, 8, NULL)
cat("22b done\n")
