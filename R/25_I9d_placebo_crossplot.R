# ==============================================================================
# 25_I9d_placebo_crossplot.R   (no INLA; resumable; run from the project root)
# I9 (d) cross-plot placebo: each interval keeps its plot's first-visit climate (zg0, zc0) but gets the warming (zg - zg0, zc - zc0) of a random interval
#   of ANOTHER plot in the same latitude band (tertiles) x period (as 22b); scenario A re-predicted (22's engine, point parameters); F as in 22.
#   Realistic warming magnitudes, no link to the plot's own warming -> slope ~ 0 expected. 20 permutations; cache figures/I/cache_I9d.
# Slopes: 22b estimator with period effects (as I9 b/c), all origins and by group.
# ==============================================================================
src <- readLines(here::here("R", "22_hindcast_space_for_time.R")); eval(parse(text = src[1:(grep("^plots <- unique", src) - 1)]))
CD9 <- file.path(FIG_DIR, "I", "cache_I9d"); dir.create(CD9, recursive = TRUE, showWarnings = FALSE); N_PERM <- 20L
con <- DBI::dbConnect(RSQLite::SQLite(), "data/raw/pep/PEP.gpkg", flags = RSQLite::SQLITE_RO)
pl <- unique(as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude FROM placette")), by = "id_pe"); DBI::dbDisconnect(con)
qs <- quantile(pl[id_pe %in% iv$id_pe, latitude], c(1/3, 2/3), na.rm = TRUE); pl[, band := cut(latitude, c(-Inf, qs, Inf), labels = c("S", "M", "N"))]
iv[pl, band := i.band, on = "id_pe"]; stopifnot(!anyNA(iv$band)); iv[, cell := paste(band, b_cal)]
for (k in seq_len(N_PERM)) { f <- file.path(CD9, sprintf("perm_%02d.rds", k)); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t_start, units = "secs") > 140) { cat("time budget reached; re-run\n"); quit(save = "no") }
  set.seed(990 + k); don <- integer(nrow(iv))
  for (cc in unique(iv$cell)) { ii <- which(iv$cell == cc); pid <- iv$id_pe[ii]
    d <- ii[sample.int(length(ii), length(ii), replace = TRUE)]; bad <- iv$id_pe[d] == pid   # donor from another plot
    while (any(bad)) { d[bad] <- ii[sample.int(length(ii), sum(bad), replace = TRUE)]; bad <- iv$id_pe[d] == pid }
    don[ii] <- d }
  ps <- copy(iv); ps[, `:=`(zg = zg0 + (iv$zg - iv$zg0)[don], zc = zc0 + (iv$zc - iv$zc0)[don])]
  out <- rbindlist(lapply(S, function(s) { ii <- which(ps$state_start == s); if (!length(ii)) return(NULL); d <- ps[ii]
    m <- metrics(pred_end(s, d, 1L, d$zg, d$zc), rep(s, nrow(d))); setnames(m, paste0("P_", names(m))); cbind(row = ii, m) }))
  setorder(out, row); saveRDS(cbind(iv[, .(id_pe, interval)], out[, -1]), f); cat(sprintf("perm %d\n", k)) }
cat("all cross-plot placebos done\n")
# ---- slopes ----
d <- readRDS(file.path(FIG_DIR, "H", "H_interval_point.rds")); setorder(d, id_pe, interval)
d[, `:=`(per = cut(year_end, c(-Inf, 1990, 2000, 2010, Inf)), prev_chg = shift(o_chg, fill = 0)), by = id_pe]; d[, la := log1p(a0)]
fit <- function(x, m, gi) { y <- x[[paste0("o_", m)]] - x[[paste0("F_", m)]]; z <- x[[paste0("P_", m)]] - x[[paste0("F_", m)]]
  fo <- paste("~ la + I(la^2) + L + prev_chg + per", if (uniqueN(x$state_start) > 1) "+ state_start" else "")
  M <- cbind(y, z, model.matrix(as.formula(fo), x)[, -1, drop = FALSE]); M <- M - (rowsum(M, gi, reorder = TRUE) / tabulate(gi))[gi, , drop = FALSE]
  C <- M[, -(1:2), drop = FALSE]; C <- C[, colSums(C^2) > 1e-10, drop = FALSE]; Q <- qr(C); yy <- qr.resid(Q, M[, 1]); zz <- qr.resid(Q, M[, 2]); sum(yy * zz) / sum(zz^2) }
ms <- c("chg", "dsti", "tmp", "bor"); GR <- c("all", "boreal", "temperate", "pioneer", "regen")
PB <- rbindlist(lapply(list.files(CD9, "^perm_", full.names = TRUE), function(f) { dd <- merge(d, readRDS(f), by = c("id_pe", "interval")); setorder(dd, id_pe, interval)
  rbindlist(lapply(GR, function(gr) { x1 <- if (gr == "all") dd else dd[grp0 == gr]; x1 <- x1[x1[, .N, by = id_pe][N > 1], on = "id_pe"]; gi <- match(x1$id_pe, unique(x1$id_pe))
    rbindlist(lapply(ms, function(m) data.table(origin = gr, metric = m, slope = fit(x1, m, gi),
      cor_PA = cor(x1[[paste0("P_", m)]] - x1[[paste0("F_", m)]], x1[[paste0("A_", m)]] - x1[[paste0("F_", m)]])))) })) }))
S9 <- PB[, .(placebo_mean = mean(slope), placebo_min = min(slope), placebo_max = max(slope), cor_placebo_actual = mean(cor_PA)), by = .(origin, metric)]
fwrite(S9, file.path(FIG_DIR, "I", "I9d_placebo_crossplot.csv")); options(width = 220)
print(dcast(S9[, .(origin, metric, t = sprintf("%.2f [%.2f, %.2f] r=%.2f", placebo_mean, placebo_min, placebo_max, cor_placebo_actual))], origin ~ metric, value.var = "t"))
