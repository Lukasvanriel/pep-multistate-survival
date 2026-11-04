# ==============================================================================
# WP1 / A1 check (read-only; writes nothing). Run from the project root:
#   after 03c:  Rscript R/tests/wp1_A1_check.R      -> part A (climate_covariates) + plausibility FLAGS
#   after 03z:  Rscript R/tests/wp1_A1_check.R      -> the same, plus part B (combined_covariates vs its pre-A1 copy)
# Needs the pre-A1 copies made by 03c / 03z: data/build/{climate_by_year_raw,climate_covariates,combined_covariates}.pre_A1.rds
# HARD STOPS (stop() on failure):
#   A1  0 NA in all averaged climate columns | A2  0 differing visits outside the 5 plots (exact, all columns) | A3  every window full
#   B1  combined: same rows/columns, identical outside the 5 plots, 0 NA in gdd_3 / cmi_gs / soil_pc1; in the 5 plots only climate columns change
# FLAGS ONLY (printed, never stop; plausibility): each of the 5 plots vs its 3 nearest plots of the pre-A1 checkpoint (same years), and the pair
#   7209002301 / 7209002302. Bound = 99th percentile of the same statistic over all plots of the pre-A1 checkpoint, same window:
#   neighbours: |value - median(3 nearest)| ; pair: |difference| over nearest-neighbour pairs <= 0.6 km.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(sf); library(here) })
source(here::here("R", "config.R"))
options(width = 220, datatable.print.nrows = 100)
FIVE <- c("7109604001", "7209002301", "7209002302", "7209605702", "7501600102")
must <- function(ok, msg) { cat(sprintf("%s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", msg)); if (!isTRUE(ok)) stop("WP1 check failed: ", msg, call. = FALSE) }
norm <- function(x) { x <- as.data.table(x); x[, id_pe := as.character(id_pe)]; if ("no_mes" %in% names(x)) x[, no_mes := as.integer(no_mes)]; x }
rd   <- function(f) norm(readRDS(file.path(BUILD_DIR, f)))
differs <- function(m, cols) sapply(cols, function(v) { a <- m[[paste0(v, ".o")]]; b <- m[[paste0(v, ".n")]]
                                                          xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b) })

## ---- part A: climate_covariates.rds vs pre-A1 -------------------------------------------------------------------------
new <- rd("climate_covariates.rds"); old <- rd("climate_covariates.pre_A1.rds")
cat(sprintf("rows old %d / new %d | plots old %d / new %d\n", nrow(old), nrow(new), uniqueN(old$id_pe), uniqueN(new$id_pe)))
must(identical(names(old), names(new)), "same columns in old and new climate_covariates")
vcols <- setdiff(names(new), c("id_pe", "no_mes", "year", "win_start", "win_end", "n_years_used"))
na_new <- colSums(is.na(new[, ..vcols])); if (any(na_new > 0)) print(na_new[na_new > 0])
must(all(na_new == 0), sprintf("A1: 0 NA in the %d averaged climate columns of the new file (old file: %d NA in gdd_3 / cmi_gs, expected 14 visits)",
                                length(vcols), old[, sum(is.na(gdd_3) | is.na(cmi_gs))]))
must(new[, all(n_years_used == win_end - win_start + 1L & n_years_used >= 1L)], "A3: every visit has n_years_used == win_end - win_start + 1 (>= 1)")
m <- merge(old, new, by = c("id_pe", "no_mes"), suffixes = c(".o", ".n"))
must(nrow(m) == nrow(old) && nrow(m) == nrow(new), "same visits (id_pe, no_mes) in old and new")
dm <- differs(m, setdiff(names(old), c("id_pe", "no_mes"))); any_d <- rowSums(dm) > 0; outside <- !(m$id_pe %in% FIVE)
if (any(any_d & outside)) { cat("columns that differ outside the 5 plots:\n"); print(colSums(dm[any_d & outside, , drop = FALSE])) }
cat(sprintf("visits that differ: %d (in the 5 plots: %d, expected 14) | outside the 5 plots: %d\n", sum(any_d), sum(any_d & !outside), sum(any_d & outside)))
must(sum(any_d & outside) == 0, "A2: 0 visits differ outside the 5 plots (exact equality, all columns)")

## ---- plausibility FLAGS (no stop) -------------------------------------------------------------------------------------
cat("\n== Plausibility of the 5 plots (FLAGS only) ==\n")
ck <- rd("climate_by_year_raw.rds"); ck[, year := as.integer(year)]
pool_ids <- unique(rd("climate_by_year_raw.pre_A1.rds")$id_pe)
cmi_gs_matched <- function(year, start_gs, end_gs, cmi_mat) {   # COPY of 03c step 4b (validated below against the stored values)
  leap  <- (year %% 4L == 0L & year %% 100L != 0L) | year %% 400L == 0L
  md    <- rbind(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
  mdays <- md[leap + 1L, , drop = FALSE]; mend <- t(apply(md, 1, cumsum))[leap + 1L, , drop = FALSE]; mstart <- mend - mdays + 1
  w <- pmax(pmin(mend, end_gs) - pmax(mstart, start_gs) + 1, 0) / mdays
  cm <- as.matrix(cmi_mat); cm[!is.na(w) & w == 0] <- 0
  rowSums(cm * w)
}
mm <- paste0("cmi_", 1:12)
ck[, cmi_gs := cmi_gs_matched(year, start_gs, end_gs, ck[, ..mm])]
setkey(ck, id_pe, year)
win_mean <- function(pe, y0, y1) { s <- ck[.(pe, as.integer(y0):as.integer(y1)), nomatch = NULL]; c(gdd_3 = mean(s$gdd_3), cmi_gs = mean(s$cmi_gs), n = nrow(s)) }
set.seed(1); smp <- old[n_years_used > 0][sample(.N, 200)]
v <- t(mapply(win_mean, smp$id_pe, smp$win_start, smp$win_end))
copy_ok <- max(abs(v[, "gdd_3"] - smp$gdd_3), abs(v[, "cmi_gs"] - smp$cmi_gs)) < 1e-8
cat(sprintf("copy of cmi_gs_matched reproduces gdd_3 / cmi_gs of 200 random old visits: %s\n", if (copy_ok) "yes" else "NO"))
if (!copy_ok) cat("*** FLAG: the check's own window means do not reproduce the stored values; neighbour / pair comparison SKIPPED\n")
if (copy_ok) {
  pep <- st_read(here("data", "raw", "pep", "PEP.gpkg"), layer = "placette", quiet = TRUE); pep$id_pe <- as.character(pep$id_pe)
  xy <- st_coordinates(st_geometry(pep))[, 1:2]                                  # metres (EPSG:32198)
  pool <- which(pep$id_pe %in% pool_ids); P <- xy[pool, , drop = FALSE]; pool_id <- pep$id_pe[pool]; n <- nrow(P)
  cat(sprintf("pool: %d plots of the pre-A1 checkpoint; computing 3 nearest neighbours of each...\n", n))
  nb <- matrix(NA_integer_, n, 3); dn1 <- numeric(n)
  for (i in seq_len(n)) { d <- (P[, 1] - P[i, 1])^2 + (P[, 2] - P[i, 2])^2; d[i] <- Inf; o <- order(d)[1:3]; nb[i, ] <- o; dn1[i] <- sqrt(d[o[1]]) }
  cache <- list()
  win_stats <- function(ws, we) {
    key <- paste(ws, we); if (!is.null(cache[[key]])) return(cache[[key]])
    s <- ck[id_pe %in% pool_id & year >= ws & year <= we, .(gdd_3 = mean(gdd_3), cmi_gs = mean(cmi_gs)), by = id_pe]
    V <- as.matrix(s[match(pool_id, id_pe), .(gdd_3, cmi_gs)])
    med <- sapply(1:2, function(j) apply(matrix(V[nb, j], ncol = 3), 1, median))
    near <- dn1 <= 600
    r <- list(q_nb = apply(abs(V - med), 2, quantile, probs = 0.99, na.rm = TRUE),
              q_pair = apply(abs(V - V[nb[, 1], , drop = FALSE])[near, , drop = FALSE], 2, quantile, probs = 0.99, na.rm = TRUE), n_near = sum(near))
    cache[[key]] <<- r; r
  }
  tg <- new[id_pe %in% FIVE][order(id_pe, no_mes)]; out <- list()
  for (r in seq_len(nrow(tg))) {
    t <- tg[r]; st <- win_stats(t$win_start, t$win_end)
    ti <- which(pep$id_pe == t$id_pe); d <- sqrt((P[, 1] - xy[ti, 1])^2 + (P[, 2] - xy[ti, 2])^2); o <- order(d)[1:3]
    nbv <- sapply(o, function(k) win_mean(pool_id[k], t$win_start, t$win_end)[c("gdd_3", "cmi_gs")])   # 2 x 3
    med <- apply(nbv, 1, median); dev <- abs(c(t$gdd_3, t$cmi_gs) - med)
    out[[r]] <- data.table(id_pe = t$id_pe, no_mes = t$no_mes, year = t$year, window = sprintf("%d-%d", t$win_start, t$win_end),
      nb = paste(sprintf("%s (%.1f km)", pool_id[o], d[o] / 1000), collapse = ", "),
      gdd_3 = round(t$gdd_3, 1), nb_gdd = paste(round(nbv[1, ], 1), collapse = " / "), med_gdd = round(med[1], 1), dev_gdd = round(dev[1], 1), q99_gdd = round(st$q_nb[1], 1),
      cmi_gs = round(t$cmi_gs, 2), nb_cmi = paste(round(nbv[2, ], 2), collapse = " / "), med_cmi = round(med[2], 2), dev_cmi = round(dev[2], 2), q99_cmi = round(st$q_nb[2], 2),
      flag = dev[1] > st$q_nb[1] || dev[2] > st$q_nb[2])
  }
  out <- rbindlist(out); print(out[, !"nb"]); cat("nearest plots per target visit (pool, EPSG:32198 distance):\n"); print(unique(out[, .(id_pe, nb)]))
  for (r in which(out$flag)) cat(sprintf("*** FLAG neighbours: %s visit %d (%s): gdd_3 %.1f vs median %.1f of 3 nearest (dev %.1f, 99%% bound %.1f) | cmi_gs %.2f vs median %.2f (dev %.2f, bound %.2f)\n",
        out$id_pe[r], out$no_mes[r], out$window[r], out$gdd_3[r], out$med_gdd[r], out$dev_gdd[r], out$q99_gdd[r], out$cmi_gs[r], out$med_cmi[r], out$dev_cmi[r], out$q99_cmi[r]))
  cat(sprintf("neighbour FLAGS: %d of %d target visits\n", sum(out$flag), nrow(out)))
  a <- new[id_pe == "7209002301", .(year, g1 = gdd_3, c1 = cmi_gs, ws = win_start, we = win_end)]; b <- new[id_pe == "7209002302", .(year, g2 = gdd_3, c2 = cmi_gs)]
  pm <- merge(a, b, by = "year")
  if (nrow(pm) == 0) cat("*** FLAG pair 7209002301 / 02: no measurement year in common\n") else {
    pm[, c("dg", "dc") := .(abs(g1 - g2), abs(c1 - c2))]
    pm[, c("qg", "qc") := {s <- win_stats(ws, we); .(s$q_pair[1], s$q_pair[2])}, by = year]
    pm[, flag := dg > qg | dc > qc]
    cat("\npair 7209002301 vs 7209002302 (0.44 km apart, adjacent cells; bound = 99th percentile of |difference| over existing pairs <= 0.6 km):\n")
    print(pm[, .(year, gdd_3_01 = round(g1, 1), gdd_3_02 = round(g2, 1), abs_diff_gdd = round(dg, 1), q99_gdd = round(qg, 1),
                 cmi_gs_01 = round(c1, 2), cmi_gs_02 = round(c2, 2), abs_diff_cmi = round(dc, 2), q99_cmi = round(qc, 2), flag)])
    for (r in which(pm$flag)) cat(sprintf("*** FLAG pair, year %d: |d gdd_3| %.1f (bound %.1f), |d cmi_gs| %.2f (bound %.2f)\n", pm$year[r], pm$dg[r], pm$qg[r], pm$dc[r], pm$qc[r]))
    cat(sprintf("pair FLAGS: %d of %d common years\n", sum(pm$flag), nrow(pm)))
  }
}

## ---- part B: combined_covariates.rds vs pre-A1 (only after 03z) --------------------------------------------------------
bak <- file.path(BUILD_DIR, "combined_covariates.pre_A1.rds")
if (!file.exists(bak)) {
  cat("\n(part B skipped: combined_covariates.pre_A1.rds does not exist yet -- run 03z, then run this script again)\n")
} else {
  cat("\n== Part B: combined_covariates.rds vs pre-A1 ==\n")
  cn <- rd("combined_covariates.rds"); co <- rd("combined_covariates.pre_A1.rds")
  must(identical(names(co), names(cn)), "B1: same columns in old and new combined table")
  must(nrow(co) == nrow(cn) && nrow(merge(co[, .(id_pe, no_mes)], cn[, .(id_pe, no_mes)], by = c("id_pe", "no_mes"))) == nrow(co), "B1: same visits (id_pe, no_mes)")
  must(!anyNA(cn[, .(gdd_3, cmi_gs, soil_pc1)]), sprintf("B1: 0 NA in gdd_3 / cmi_gs / soil_pc1 of the new combined table (old: %d visits with NA)", co[, sum(is.na(gdd_3) | is.na(cmi_gs) | is.na(soil_pc1))]))
  mb <- merge(co, cn, by = c("id_pe", "no_mes"), suffixes = c(".o", ".n")); dmb <- differs(mb, setdiff(names(co), c("id_pe", "no_mes")))
  anyb <- rowSums(dmb) > 0; outb <- !(mb$id_pe %in% FIVE)
  must(sum(anyb & outb) == 0, sprintf("B1: 0 visits differ outside the 5 plots (%d visits differ in the 5 plots)", sum(anyb & !outb)))
  clim_cols <- setdiff(names(new), c("id_pe", "no_mes", "year")); chg <- colnames(dmb)[colSums(dmb[anyb & !outb, , drop = FALSE]) > 0]
  cat("columns that changed in the 5 plots:", paste(chg, collapse = ", "), "\n")
  must(all(chg %in% clim_cols), "B1: in the 5 plots only climate columns changed")
}
cat("\nALL HARD CHECKS PASSED (flags above, if any, are plausibility notes only)\n")
