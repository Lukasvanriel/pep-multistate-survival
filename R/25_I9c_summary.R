# ==============================================================================
# 25_I9c_summary.R   (no INLA; ~2 min; after 25_I9b_placebo.R; run from the project root)
# I9 (b) placebo slopes: (observed - F) on (placebo A - F), 22b estimator with period effects; spread over the 20 permutations.
# I9 (c) extra controls on the real slope: + disturbance flags of the interval (fire, pest, light harvest, heavy harvest); + latitude band x period.
# ==============================================================================
source(here::here("R", "25_I_common.R")); set.seed(9)
d <- readRDS(file.path(FIG_DIR, "H", "H_interval_point.rds"))
vi2 <- vi[, .(id_pe, interval, ff = as.integer(fire == 1L), pf = as.integer(pest == 1L), hl = as.integer(harv == 1L), hh = as.integer(harv == 2L))]
d <- merge(d, vi2, by = c("id_pe", "interval")); setorder(d, id_pe, interval)
con <- DBI::dbConnect(RSQLite::SQLite(), "data/raw/pep/PEP.gpkg", flags = RSQLite::SQLITE_RO)
pl <- unique(as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude FROM placette")), by = "id_pe"); DBI::dbDisconnect(con)
qs <- quantile(pl[id_pe %in% d$id_pe, latitude], c(1/3, 2/3), na.rm = TRUE); pl[, band := cut(latitude, c(-Inf, qs, Inf), labels = c("S", "M", "N"))]
d <- merge(d, pl[, .(id_pe, band)], by = "id_pe", all.x = TRUE); setorder(d, id_pe, interval); stopifnot(!anyNA(d$band))
d[, `:=`(per = cut(year_end, c(-Inf, 1990, 2000, 2010, Inf)), prev_chg = shift(o_chg, fill = 0)), by = id_pe]; d[, la := log1p(a0)]
fit <- function(x, m, gi, pred = "A", extra = "") { y <- x[[paste0("o_", m)]] - x[[paste0("F_", m)]]; z <- x[[paste0(pred, "_", m)]] - x[[paste0("F_", m)]]
  fo <- paste("~ la + I(la^2) + L + prev_chg + per", if (uniqueN(x$state_start) > 1) "+ state_start" else "", extra)
  M <- cbind(y, z, model.matrix(as.formula(fo), x)[, -1, drop = FALSE]); M <- M - (rowsum(M, gi, reorder = TRUE) / tabulate(gi))[gi, , drop = FALSE]
  C <- M[, -(1:2), drop = FALSE]; C <- C[, colSums(C^2) > 1e-10, drop = FALSE]; Q <- qr(C); yy <- qr.resid(Q, M[, 1]); zz <- qr.resid(Q, M[, 2]); sum(yy * zz) / sum(zz^2) }
ms <- c("chg", "dsti", "tmp", "bor"); GR <- c("all", "boreal", "temperate", "pioneer", "regen")
sub <- function(gr) { x1 <- if (gr == "all") d else d[grp0 == gr]; x1 <- x1[x1[, .N, by = id_pe][N > 1], on = "id_pe"]; list(x = x1, gi = match(x1$id_pe, unique(x1$id_pe))) }
# (b) placebo
PF <- list.files(file.path(FIG_DIR, "I", "cache_I9"), "^perm_", full.names = TRUE); stopifnot(length(PF) == 20L)
pb <- rbindlist(lapply(seq_along(PF), function(k) { p <- readRDS(PF[k]); dd <- merge(d, p, by = c("id_pe", "interval")); setorder(dd, id_pe, interval)
  if (k == 1L) { one <- dd[, .N, by = id_pe][N == 1, id_pe]; stopifnot(max(abs(dd[id_pe %in% one, P_chg - A_chg])) < 1e-10) }   # unshuffled plots reproduce A
  rbindlist(lapply(GR, function(gr) { x1 <- if (gr == "all") dd else dd[grp0 == gr]; x1 <- x1[x1[, .N, by = id_pe][N > 1], on = "id_pe"]; gi <- match(x1$id_pe, unique(x1$id_pe))
    rbindlist(lapply(ms, function(m) data.table(perm = k, origin = gr, metric = m, slope = fit(x1, m, gi, "P"),
      cor_PA = cor(x1[[paste0("P_", m)]] - x1[[paste0("F_", m)]], x1[[paste0("A_", m)]] - x1[[paste0("F_", m)]])))) })) }))
PB <- pb[, .(placebo_mean = mean(slope), placebo_min = min(slope), placebo_max = max(slope), cor_placebo_actual = mean(cor_PA)), by = .(origin, metric)]
# (c) extra controls (bootstrap 100 as 22b)
EX <- c(base = "", disturbance = "+ ff + pf + hl + hh", `disturbance + band x period` = "+ ff + pf + hl + hh + band:per")
ct <- rbindlist(lapply(GR, function(gr) { s <- sub(gr); ixg <- split(seq_len(nrow(s$x)), s$gi)
  rbindlist(lapply(names(EX), function(e) rbindlist(lapply(ms, function(m) { est <- fit(s$x, m, s$gi, "A", EX[[e]])
    bt <- replicate(100, { b <- sample(length(ixg), replace = TRUE); ii <- unlist(ixg[b], use.names = FALSE); fit(s$x[ii], m, rep(seq_along(b), lengths(ixg[b])), "A", EX[[e]]) })
    data.table(origin = gr, controls = e, metric = m, slope = est, lo = quantile(bt, .025), hi = quantile(bt, .975)) })))) }))
R <- merge(ct, PB, by = c("origin", "metric")); fwrite(R, file.path(OUT, "I9bc_placebo_controls.csv"))
options(width = 220)
cat("\nPLACEBO (mean [min, max] over 20 shuffles):\n"); print(dcast(PB[, .(origin, metric, t = sprintf("%.2f [%.2f, %.2f] r=%.2f", placebo_mean, placebo_min, placebo_max, cor_placebo_actual))], origin ~ metric, value.var = "t"))
cat("\nCONTROLS:\n"); print(dcast(ct[, .(origin, controls, metric, t = sprintf("%.2f [%.2f, %.2f]", slope, lo, hi))], origin + controls ~ metric, value.var = "t"))
