# ==============================================================================
# 25_I15_elasticity.R   (no INLA; resumable: re-run until "all I15 perturbations done"; run from the project root)
# I15: which transitions control the headline numbers. One transition at a time, point estimates:
#   intercept + log 1.1, GDD effect + log 1.1 (additive on the log-HR per SD), alpha x 1.1 (decade intercept b0 held fixed -> pivot at 10 y; variant 0),
#   heavy-harvest main term + log 1.1 (transitions with one).
# Headlines: (1) D2 boreal <-> temperate balance point (GDD where A = 0; A from 19's group reach within 50 y, first-visit weights; zero crossing by
#   linear interpolation on a 0.1-SD grid); (2) D3 all-plots change in P(leave within 50 y) under heavy harvest at entry; (3) D3 all-plots change in the
#   expected thermal shift of the next state under GDD +1 SD. D3 perturbations only touch the plots of the perturbed transition's origin state.
# 1-SE column by linear scaling: effect x SE / log 1.1 (alpha: SE / (0.1 alpha)); checked by direct 1-SE runs for the top 3 per headline.
# ==============================================================================
src <- readLines(here::here("R", "19_figures_D2_direction.R")); eval(parse(text = src[1:(grep("^RUNS <- CJ", src) - 1)]))
s20 <- readLines(here::here("R", "20_figures_D3_real_plots.R")); eval(parse(text = s20[grep("^sc <- readRDS", s20):(grep("^PERT <- list", s20) + 1)]))
eval(parse(text = s20[grep("^first_move <- function", s20):(grep("^# point estimates, all plots", s20) - 1)]))
OUTI <- file.path(FIG_DIR, "I"); CI <- file.path(OUTI, "cache_I15"); dir.create(CI, showWarnings = FALSE); ROWS <- 1L; t0 <- Sys.time()
pooled <- fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG)))
H <- 50; DTG <- 0.5; TG <- seq(0, H, by = DTG); NG <- length(TG) - 1L   # as 20
PAR0 <- PAR; L11 <- log(1.1)
perturb <- function(lb, kind, size = NULL) { X <- PAR0[[lb]]
  if (kind == "intercept") { d <- if (is.null(size)) L11 else size; X[, "(Intercept)"] <- X[, "(Intercept)"] + d }
  if (kind == "GDD")       { d <- if (is.null(size)) L11 else size; X[, "z_gdd_3"] <- X[, "z_gdd_3"] + d }
  if (kind == "alpha")     { d <- if (is.null(size)) 0.1 * X[1, "alpha"] else size; X[, "alpha"] <- X[, "alpha"] + d }
  if (kind == "heavy")     { d <- if (is.null(size)) L11 else size; X[, "harv_heavy"] <- X[, "harv_heavy"] + d }
  X[, "b0y"] <- X[, "(Intercept)"] - X[, "alpha"] * log(10); X }
# ---- headline functions (use global PAR) ----
gw <- function(g) { x <- w0[state %in% names(GRP)[GRP == g]]; setNames(x$N / sum(x$N), x$state) }   # as 19
TGT <- list(temperate = names(GRP)[GRP == "temperate"], boreal = names(GRP)[GRP == "boreal"]); wB <- gw("boreal"); wT <- gw("temperate")
Az <- function(z) { pT <- drop(reach_set(TGT$temperate, z)[1, match(names(wB), S)] %*% wB); pB <- drop(reach_set(TGT$boreal, z)[1, match(names(wT), S)] %*% wT); (pT - pB) / (pT + pB) }
balance <- function(z0 = 0.3) { zs <- round(z0 + c(-0.1, 0, 0.1), 1); a <- sapply(zs, Az)
  repeat { k <- which(diff(sign(a)) != 0)
    if (length(k)) { k <- k[1]; return(zs[k] - a[k] * (zs[k + 1] - zs[k]) / (a[k + 1] - a[k])) }
    if (a[1] > 0) { zs <- c(round(zs[1] - 0.1, 1), zs); a <- c(Az(zs[1]), a) } else { zs <- c(zs, round(zs[length(zs)] + 0.1, 1)); a <- c(a, Az(zs[length(zs)])) }
    if (length(zs) > 25) return(NA_real_) } }
NPL <- nrow(pl)
d3_state <- function(st) { P <- pl[state == st]; s <- match(st, S); b <- first_move(s, P, 1L, PERT$baseline)
  c(heavy = sum(first_move(s, P, 1L, PERT$`heavy harvest`)[, 1] - b[, 1]), gdd = sum(first_move(s, P, 1L, PERT$`GDD +1 SD`)[, 2] - b[, 2])) / NPL }
# ---- baseline ----
fb <- file.path(CI, "baseline.rds")
if (!file.exists(fb)) { D3B <- sapply(unique(pl$state), d3_state); saveRDS(list(bp = balance(), d3 = D3B), fb) }
BASE <- readRDS(fb); cat(sprintf("baseline: balance point %.3f SD; heavy %.3f; GDD direction %.3f\n", BASE$bp, sum(BASE$d3["heavy", ]), sum(BASE$d3["gdd", ])))
# ---- perturbations ----
PT <- rbind(CJ(label = names(PAR0), kind = c("intercept", "alpha")), CJ(label = names(PAR0)[sapply(PAR0, function(X) "z_gdd_3" %in% colnames(X))], kind = "GDD"),
            CJ(label = names(PAR0)[sapply(PAR0, function(X) "harv_heavy" %in% colnames(X))], kind = "heavy"))
PT[, from := sub("->.*", "", label)]; cat(nrow(PT), "perturbations\n")
run_one <- function(lb, kind, size = NULL) { PAR <<- PAR0; PAR[[lb]] <<- perturb(lb, kind, size); fr <- sub("->.*", "", lb)
  d3 <- if (fr %in% colnames(BASE$d3)) d3_state(fr) - BASE$d3[, fr] else c(heavy = 0, gdd = 0)
  bp <- if (kind == "heavy") BASE$bp else balance(BASE$bp); PAR <<- PAR0
  data.table(label = lb, kind = kind, d_bp = bp - BASE$bp, d_heavy = d3[["heavy"]], d_gdd = d3[["gdd"]]) }
CH <- split(seq_len(nrow(PT)), ceiling(seq_len(nrow(PT)) / 20))
for (ci in seq_along(CH)) { f <- file.path(CI, sprintf("chunk_%02d.rds", ci)); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > 100) { cat("time budget reached; re-run\n"); quit(save = "no") }
  saveRDS(rbindlist(lapply(CH[[ci]], function(i) run_one(PT$label[i], PT$kind[i]))), f); cat(sprintf("chunk %d/%d (%.0f s)\n", ci, length(CH), as.numeric(Sys.time() - t0, units = "secs"))) }
cat("all I15 perturbations done\n")
# ---- summary: rankings, 1-SE column, direct 1-SE check (top 3 per headline) ----
R <- rbindlist(lapply(list.files(CI, "^chunk_", full.names = TRUE), readRDS)); stopifnot(nrow(R) == nrow(PT))
PN <- c(intercept = "(Intercept)", GDD = "z_gdd_3", alpha = "alpha", heavy = "harv_heavy")
R[, parameter := PN[kind]]; R <- merge(R, pooled[, .(label, parameter, est = pooled_est, se)], by = c("label", "parameter"), all.x = TRUE)
R[, f_se := fifelse(kind == "alpha", se / (0.1 * est), se / L11)]
HL <- c(d_bp = "balance point (SD)", d_heavy = "heavy harvest: change in P(leave)", d_gdd = "GDD +1 SD: change in thermal direction")
L <- rbindlist(lapply(names(HL), function(h) R[, .(headline = HL[[h]], label, kind, se, per10 = get(h), per1se = get(h) * f_se)]))
L[, `:=`(share10 = abs(per10) / sum(abs(per10)), share1se = abs(per1se) / sum(abs(per1se))), by = headline]
L[, `:=`(rank10 = frank(-abs(per10), ties.method = "first"), rank1se = frank(-abs(per1se), ties.method = "first")), by = headline]
fc <- file.path(CI, "direct_1se.rds")
if (!file.exists(fc)) { top <- unique(L[rank1se <= 3, .(headline, label, kind, se)])
  DC <- rbindlist(lapply(seq_len(nrow(top)), function(i) { r <- run_one(top$label[i], top$kind[i], size = top$se[i]); h <- names(HL)[HL == top$headline[i]]
    data.table(headline = top$headline[i], label = top$label[i], kind = top$kind[i], direct1se = r[[h]]) })); saveRDS(DC, fc) }
DC <- readRDS(fc); L <- merge(L, DC, by = c("headline", "label", "kind"), all.x = TRUE)
fwrite(L, file.path(OUTI, "I15_elasticity.csv"))
SDG <- sc[var == "gdd_3", sd]; options(width = 230)
for (h in HL) { x <- L[headline == h]; cat("\n====", h, "====\n")
  cat(sprintf("top 10 (x1.1) share of total |effect|: %.2f; top 10 (1 SE): %.2f\n", sum(x[rank10 <= 10, share10]), sum(x[rank1se <= 10, share1se])))
  sc_ <- if (h == HL[["d_bp"]]) SDG else 1
  cat("by x1.1:\n"); print(x[rank10 <= 10][order(rank10), .(label, kind, per10 = signif(per10, 2), dd = if (sc_ > 1) round(per10 * sc_, 1) else NA, share = round(share10, 3))])
  cat("by 1 SE:\n"); print(x[rank1se <= 10][order(rank1se), .(label, kind, per1se = signif(per1se, 2), dd = if (sc_ > 1) round(per1se * sc_, 1) else NA, direct1se = signif(direct1se, 2), share = round(share1se, 3))]) }
