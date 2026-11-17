# ==============================================================================
# 25_I14_first_passage_temperate.R   (no INLA; resumable: re-run until "all I14 runs done"; run from the project root)
# I14: first passage into the temperate group (sugar maple, red maple, yellow birch, cedar, other temperate; all absorbing) from each non-temperate start,
#   along GDD (-2 to +2 SD by 0.5; CMI and soil 0), undisturbed, entry into the start state at t = 0. 17's engine (1-y steps, competing exits, moves
#   through any other non-temperate state allowed). Point + 50 draws (first 50 of 17's draws). Outputs: P(reached by 50 / 100 / 200 y), median time.
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
OUTI <- file.path(FIG_DIR, "I"); CI <- file.path(OUTI, "cache_I14"); dir.create(CI, showWarnings = FALSE); t0 <- Sys.time()
source(here::here("R", "00_groups.R"))   # GRP: thermal groups
ABS <- which(GRP[S] == "temperate"); ST <- setdiff(seq_len(nS), ABS); ROWS <- 1:51; TM <- 200L
fp <- function(zg) { D <- length(ROWS); nA <- TM + 1L; ages <- 0:(nA - 1)
  TR <- lapply(S, function(o) lapply(EX[[o]]$lb, function(lb) { X <- PAR[[lb]][ROWS, , drop = FALSE]; a <- X[, "alpha"]
    list(eta = X[, "b0y"] + gcol(X, "z_gdd_3") * zg, inc = outer(a, ages + 1, function(aa, x) x^aa) - outer(a, ages, function(aa, x) x^aa)) })); names(TR) <- S
  M <- array(0, c(D, nS, nS, nA)); for (i in ST) M[, i, i, 1] <- 1; cum <- array(0, c(D, nS, TM)); A <- matrix(0, D, nS)
  for (st in seq_len(TM)) { New <- array(0, dim(M))
    for (s in ST) { ex <- EX[[s]]; Ms <- M[, , s, , drop = FALSE]; dim(Ms) <- c(D, nS, nA)
      if (!length(ex$lb)) { New[, , s, 2:nA] <- New[, , s, 2:nA] + Ms[, , 1:(nA - 1)]; next }
      dH <- lapply(TR[[s]], function(tr) exp(tr$eta) * tr$inc); tot <- Reduce(`+`, dH); pl <- 1 - exp(-tot)
      New[, , s, 2:nA] <- New[, , s, 2:nA] + sweep(Ms, c(1, 3), 1 - pl, "*")[, , 1:(nA - 1)]
      for (k in seq_along(dH)) { pk <- ifelse(tot > 0, dH[[k]] / tot, 0) * pl; flow <- apply(sweep(Ms, c(1, 3), pk, "*"), c(1, 2), sum); j <- ex$to[k]
        if (j %in% ABS) A <- A + flow else New[, , j, 1] <- New[, , j, 1] + flow } }
    M <- New; cum[, , st] <- A }
  cum }
ZG <- seq(-2, 2, 0.5)
for (z in ZG) { f <- file.path(CI, sprintf("z%+.1f.rds", z)); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t0, units = "secs") > 120) { cat("time budget reached; re-run\n"); quit(save = "no") }
  t1 <- Sys.time(); cum <- fp(z)
  res <- rbindlist(lapply(ST, function(i) { C <- cum[, i, ]; med <- apply(C, 1, function(v) if (max(v) >= 0.5) which(v >= 0.5)[1] else NA_real_)
    q <- function(x) c(x[1], quantile(x[-1], c(.025, .975), na.rm = TRUE))
    data.table(start = S[i], zg = z, what = c("P by 50 y", "P by 100 y", "P by 200 y", "median time (y)"),
               rbind(q(C[, 50]), q(C[, 100]), q(C[, 200]), q(med))) }))
  setnames(res, c("start", "zg", "what", "est", "lo", "hi")); saveRDS(res, f); cat(sprintf("z %+.1f (%.0f s)\n", z, as.numeric(Sys.time() - t1, units = "secs"))) }
cat("all I14 runs done\n")
# ---- summary + figure ----
source(here::here("R", "25_I_common.R"))
FP <- rbindlist(lapply(list.files(CI, full.names = TRUE), readRDS)); fwrite(FP, file.path(OUTI, "I14_first_passage_temperate.csv"))
options(width = 220)
cat("P(reached temperate by 100 y) [95% CI]:\n")
print(dcast(FP[what == "P by 100 y" & zg %in% c(-2, -1, 0, 1, 2), .(start, zg, t = sprintf("%.2f [%.2f, %.2f]", est, lo, hi))], start ~ zg, value.var = "t")[order(match(start, ORDER))])
cat("\nmedian time to temperate (y; NA = not reached by 200 y):\n")
print(dcast(FP[what == "median time (y)" & zg %in% c(-2, -1, 0, 1, 2), .(start, zg, t = ifelse(is.na(est), ">200", sprintf("%.0f", est)))], start ~ zg, value.var = "t")[order(match(start, ORDER))])
x <- FP[what == "P by 100 y"][, start := factor(start, ORDER)]
p <- ggplot(x, aes(zg, est, colour = start, fill = start)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() +
  scale_colour_manual(values = PAL, labels = lab, name = "start state") + scale_fill_manual(values = PAL, guide = "none") +
  scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
  labs(x = "GDD (SD)", y = "probability of having entered a temperate state by 100 y", title = "I14  First passage into the temperate group along the GDD gradient",
       subtitle = "undisturbed, CMI and soil at the mean, entry into the start state at t = 0; 95% CI from 50 parameter draws; beyond the 55-y record")
save_fig(p, "I14_first_passage_temperate", 9, 5.5, NULL)
