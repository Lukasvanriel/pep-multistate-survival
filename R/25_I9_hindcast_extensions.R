# ==============================================================================
# 25_I9_hindcast_extensions.R   (no INLA; run after 22_hindcast_space_for_time.R; run from the project root)
# I9 (a): the within-plot slopes of 22b (observed - F on A - F, plot fixed effects, same controls, with calendar-period effects) by individual origin state.
# ==============================================================================
source(here::here("R", "25_I_common.R")); set.seed(9)
d <- readRDS(file.path(FIG_DIR, "H", "H_interval_point.rds")); setorder(d, id_pe, interval)
d[, `:=`(per = cut(year_end, c(-Inf, 1990, 2000, 2010, Inf)), prev_chg = shift(o_chg, fill = 0)), by = id_pe]; d[, la := log1p(a0)]
fit <- function(x, m, pfe, gi) { y <- x[[paste0("o_", m)]] - x[[paste0("F_", m)]]; z <- x[[paste0("A_", m)]] - x[[paste0("F_", m)]]   # as 22b
  fo <- paste("~ la + I(la^2) + L + prev_chg", if (uniqueN(x$state_start) > 1) "+ state_start" else "", if (pfe) "+ per" else "")
  M <- cbind(y, z, model.matrix(as.formula(fo), x)[, -1, drop = FALSE]); M <- M - (rowsum(M, gi, reorder = TRUE) / tabulate(gi))[gi, , drop = FALSE]
  C <- M[, -(1:2), drop = FALSE]; C <- C[, colSums(C^2) > 1e-10, drop = FALSE]; Q <- qr(C); yy <- qr.resid(Q, M[, 1]); zz <- qr.resid(Q, M[, 2])
  c(slope = sum(yy * zz) / sum(zz^2), sd_z = sqrt(mean(zz^2))) }
ms <- c("chg", "dsti", "tmp", "bor")
res <- rbindlist(lapply(ORDER, function(st) rbindlist(lapply(ms, function(m) {
  x1 <- d[state_start == st]; x1 <- x1[x1[, .N, by = id_pe][N > 1], on = "id_pe"]; gi <- match(x1$id_pe, unique(x1$id_pe))
  est <- fit(x1, m, TRUE, gi); ixg <- split(seq_len(nrow(x1)), gi)
  bt <- replicate(100, { b <- sample(length(ixg), replace = TRUE); ii <- unlist(ixg[b], use.names = FALSE); fit(x1[ii], m, TRUE, rep(seq_along(b), lengths(ixg[b])))[1] })
  data.table(origin = st, group = GRP[st], metric = m, n_int = nrow(x1), n_plots = length(ixg), sd_pred = est[2], slope = est[1], lo = quantile(bt, .025), hi = quantile(bt, .975)) }))))
fwrite(res, file.path(OUT, "I9a_slopes_by_origin_state.csv"))
print(dcast(res[, .(origin, n_int, metric, txt = sprintf("%.2f [%.2f, %.2f]", slope, lo, hi))], origin + n_int ~ metric, value.var = "txt")[order(match(origin, ORDER))])
print(res[metric == "chg", .(origin, sd_pred = signif(sd_pred, 2))])
