# ==============================================================================
# 19_figures_D2_direction.R   (no INLA; resumable: re-run until "all D2 runs done", then figures; run from the project root)
# Exploratory set D2: directional replacement. NOT manuscript figures. Engine as 17 (1-y steps, competing exits, undisturbed).
# Pair level: asymmetry A_ij = (P_ij - P_ji) / (P_ij + P_ji), P_ij = P(start in i, reached j at least once within 50 y) from 17's first-passage runs
#   (point estimates; only pairs with both directions possible), GDD -1/0/+1.
# Group level: P(start in group G1, reached any state of G2 within 50 y); starts within G1 weighted by the observed first-visit composition
#   (plot_info first_state); all states of G2 absorbing. GDD grid -1.5..1.5 (CMI, soil 0); 50 draws + point -> 95% CI, also for A_G1G2.
# Groups (user 2026-11-09): boreal = black spruce, white spruce, fir, jack pine; temperate = sugar maple, red maple, yellow birch, cedar, other temperate;
#   pioneer = white birch, trembling aspen, other pioneer. Caution: GDD values far from where a state occurs extrapolate in climate space.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
OUTD <- file.path(FIG_DIR, "D"); OUT <- OUTD; CD <- file.path(OUTD, "cache"); dir.create(CD, recursive = TRUE, showWarnings = FALSE)
ROWS <- 1:51; TIME_BUDGET <- 150; t0 <- Sys.time(); H <- 50L
source(here::here("R", "00_groups.R"))   # GRP: thermal groups
G3 <- c("boreal", "temperate", "pioneer")
pinfo <- readRDS(file.path(TABLE_DIR, "plot_info.rds")); w0 <- pinfo[, .N, by = .(state = first_state)]
# probability of having reached any state in 'targets' by H years, from every start state (absorbing set), [draw, start]
reach_set <- function(targets, zg) { D <- length(ROWS); nA <- H + 1L; ages <- 0:(nA - 1); tg <- match(targets, S)
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][ROWS, , drop = FALSE]; a <- X[, "alpha"]
    list(eta = X[, "b0y"] + gcol(X, "z_gdd_3") * zg, inc = outer(a, ages + 1, function(aa, x) x^aa) - outer(a, ages, function(aa, x) x^aa)) })); names(TR) <- S
  M <- array(0, c(D, nS, nS, nA)); for (i in seq_len(nS)) M[, i, i, 1] <- 1; absorbed <- matrix(0, D, nS)
  for (st in seq_len(H)) { New <- array(0, dim(M))
    for (s in setdiff(seq_len(nS), tg)) { ex <- EX[[s]]; Ms <- M[, , s, , drop = FALSE]; dim(Ms) <- c(D, nS, nA)
      if (!length(ex$lb)) { New[, , s, 2:nA] <- New[, , s, 2:nA] + Ms[, , 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta) * tr$inc); tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      stayM <- sweep(Ms, c(1, 3), 1 - pl, "*"); New[, , s, 2:nA] <- New[, , s, 2:nA] + stayM[, , 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; fl <- apply(sweep(Ms, c(1, 3), pk, "*"), 1:2, sum); j <- ex$to[k]
        if (j %in% tg) absorbed <- absorbed + fl else New[, , j, 1] <- New[, , j, 1] + fl } }
    M <- New }
  for (j in tg) absorbed[, j] <- 1                                           # a start inside the target set has reached it at t = 0
  absorbed }
RUNS <- CJ(target = G3, zg = seq(-1.5, 1.5, by = 0.5)); RUNS[, id := sprintf("d2_%s_z%+.1f", target, zg)]
for (i in seq_len(nrow(RUNS))) { f <- file.path(CD, paste0(RUNS$id[i], ".rds")); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > TIME_BUDGET) { cat("time budget reached; re-run to continue\n"); quit(save = "no") }
  saveRDS(reach_set(names(GRP)[GRP == RUNS$target[i]], RUNS$zg[i]), f); cat("  D2", RUNS$id[i], sprintf("(%.0f s)\n", as.numeric(Sys.time() - t0, units = "secs"))) }
cat("all D2 runs done\n")

# ---- group level ----------------------------------------------------------------------------------------------------------------------------------------------
gw <- function(g) { x <- w0[state %in% names(GRP)[GRP == g]]; setNames(x$N / sum(x$N), x$state) }
GP <- rbindlist(lapply(seq_len(nrow(RUNS)), function(i) { A <- readRDS(file.path(CD, paste0(RUNS$id[i], ".rds")))
  rbindlist(lapply(setdiff(G3, RUNS$target[i]), function(g1) { w <- gw(g1); v <- drop(A[, match(names(w), S), drop = FALSE] %*% w)
    data.table(from_g = g1, to_g = RUNS$target[i], zg = RUNS$zg[i], draw = seq_along(v) - 1L, p = v) })) }))
GA <- merge(GP, GP[, .(from_g = to_g, to_g = from_g, zg, draw, p_rev = p)], by = c("from_g", "to_g", "zg", "draw"))
GA[, A := (p - p_rev) / (p + p_rev)]
GS <- GA[order(draw), .(p_est = p[1], p_lo = quantile(p[-1], .025), p_hi = quantile(p[-1], .975), A_est = A[1], A_lo = quantile(A[-1], .025), A_hi = quantile(A[-1], .975)), by = .(from_g, to_g, zg)]
fwrite(GS, file.path(OUTD, "D2_group_reach_and_asymmetry.csv"))
GS[, pair := paste(from_g, "→", to_g)]
p1 <- ggplot(GS, aes(zg, p_est, colour = pair, fill = pair)) + geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = function(v) paste0(100 * v, "%")) + labs(x = "GDD (SD from the mean)", y = "probability of reaching the group within 50 y", colour = NULL, fill = NULL,
  subtitle = "start = observed first-visit composition of the start group")
AS <- GS[from_g < to_g]
p2 <- ggplot(AS, aes(zg, A_est)) + geom_hline(yintercept = 0, linetype = 2) + geom_ribbon(aes(ymin = A_lo, ymax = A_hi), fill = "grey80") + geom_line(linewidth = 0.8) +
  facet_wrap(~ paste(from_g, "vs", to_g)) + labs(x = "GDD (SD from the mean)", y = "asymmetry (> 0 = first group → second dominates)",
  subtitle = "band = 95% CI")
save_fig(p1 / p2 + plot_annotation(title = "D2a  Directional replacement between boreal, temperate and pioneer groups along GDD (undisturbed)"), "D2a_group_direction", 12, 10, NULL)

# ---- pair level (from 17's first-passage runs, point estimates) ---------------------------------------------------------------------------------------------
FPS <- rbindlist(lapply(list.files(CACHE, "^fp_.*rds$", full.names = TRUE), readRDS))[t == H & start != target]
PA <- merge(FPS[, .(i = start, j = target, zg, pij = est)], FPS[, .(i = target, j = start, zg, pji = est)], by = c("i", "j", "zg"))
PA <- PA[pij + pji > 0.01][, A := (pij - pji) / (pij + pji)]; fwrite(PA, file.path(OUTD, "D2_pair_asymmetry.csv"))
for (zz in c(-1, 0, 1)) { g <- CJ(i = ORDER, j = ORDER); g <- merge(g, PA[zg == zz], by = c("i", "j"), all.x = TRUE)
  g[, `:=`(i = factor(i, levels = rev(ORDER)), j = factor(j, levels = ORDER), txt = fifelse(is.na(A), "", sprintf("%+.2f", A)))]
  p <- ggplot(g, aes(j, i)) + geom_tile(aes(fill = A), colour = "white") + geom_text(aes(label = txt), size = 2.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), na.value = "grey95", name = "asymmetry") +
    scale_x_discrete(labels = function(x) gsub("_", " ", x), position = "top") + scale_y_discrete(labels = function(x) gsub("_", " ", x)) + coord_equal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank()) +
    labs(x = "state B", y = "state A", title = sprintf("D2b  Pairwise asymmetry within 50 y (undisturbed, %s)", gsd(zz)),
         subtitle = "red = A → B more likely than B → A; blank = both < 1%")
  save_fig(p, sprintf("D2b_pair_asymmetry_gdd%+d", zz), 9, 8.5, NULL) }
cat("D2 figures done\n")
