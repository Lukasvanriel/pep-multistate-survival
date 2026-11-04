# ==============================================================================
# 06x_inla_mode_classic.R   RUN ON THE MACHINE WITH INLA (~ 10-15 min)
# FINDING (06u trace): every run goes through "inla_INLA_preopt_experimental ... Mode [Compact]", i.e. INLA's newer compact/pre-optimised code path.
# In the free m08 fit the SAME theta (the mode, 5.3678) gives maxld -2166.3 during the optimisation but ldens -23417 (relative) during integration ("(interupt)"),
# the Hessian eigenvalue is 5022 (expected ~12) and gradient steps return |grad| 1e6. So the latent-field solve is not reproducible in this code path.
# TEST: switch INLA to its long-standing CLASSIC mode (inla.setOption(inla.mode = "classic")) -- a global, documented switch, no data or model change.
#  Part 1: free fits of all saved collapsed data sets (8 black_spruce->regen + other_temperate->white_birch ones)    -> ratio, alpha, sd, mlik
#  Part 2: m02 fixed-alpha fine grid 1.690-1.710 and cold-start check                                                 -> mlik should be smooth, ~ -2165
# Output: data/processed/fit/inla_mode_classic_<tag>.csv, inla_mode_classic_grid_<tag>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(INLA) })
RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit"); options(width = 200)
ALPHA_NAME <- "alpha parameter for weibullsurv"
cat("inla.mode before:", tryCatch(inla.getOption("inla.mode"), error = function(e) paste("ERROR", conditionMessage(e))), "\n")
MODE <- commandArgs(trailingOnly = TRUE); MODE <- if (length(MODE)) MODE[1] else "classic"
tryCatch(inla.setOption(inla.mode = MODE), error = function(e) stop("inla.setOption(inla.mode) failed: ", conditionMessage(e)))
cat("inla.mode now   :", inla.getOption("inla.mode"), "\n")
run <- function(d, rhs, fixed_alpha = NULL) {
  ci <- list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))
  fam <- list(variant = 1); if (!is.null(fixed_alpha)) fam$hyper <- list(alpha = list(initial = 10 * log(fixed_alpha), fixed = TRUE))
  t0 <- Sys.time(); r <- tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
                control.family = fam, control.inla = ci), error = function(e) conditionMessage(e)); attr(r, "secs") <- as.numeric(difftime(Sys.time(), t0, units = "secs")); r }
# Part 1
out <- list(); files <- sort(list.files(OUT_DIR, pattern = sprintf("^baddata_%s_m[0-9]+\\.rds$", RUN_TAG), full.names = TRUE))
for (f in files) { mm <- as.integer(sub(".*_m([0-9]+)\\.rds$", "\\1", f)); B <- readRDS(f)
  for (lab in names(B)) { b <- B[[lab]]; ev <- sum(b$data$status == 1L); r <- run(b$data, b$rhs)
    out[[length(out) + 1]] <- if (is.character(r)) data.table(m = mm, label = lab, alpha = NA_real_, sd = NA_real_, ratio = NA_real_, mlik = NA_real_, secs = attr(r, "secs"), error = substr(r, 1, 60)) else {
      a <- r$summary.hyperpar[ALPHA_NAME, "mean"]; s <- r$summary.hyperpar[ALPHA_NAME, "sd"]
      data.table(m = mm, label = lab, alpha = a, sd = s, ratio = s / (0.78 * a / sqrt(ev)), mlik = as.numeric(r$mlik[1, 1]), secs = attr(r, "secs"), error = NA_character_) }
    print(tail(rbindlist(out), 1)[, .(m, label, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1), secs = round(secs, 1), error)])
    fwrite(rbindlist(out), file.path(OUT_DIR, sprintf("inla_mode_%s_%s.csv", MODE, RUN_TAG))) } }
# Part 2
b <- readRDS(file.path(OUT_DIR, sprintf("baddata_%s_m02.rds", RUN_TAG)))[["black_spruce->regen"]]; g <- list()
for (a in seq(1.690, 1.710, by = 0.001)) { r <- run(b$data, b$rhs, a); g[[length(g) + 1]] <- data.table(alpha = a, mlik = if (is.character(r)) NA_real_ else as.numeric(r$mlik[1, 1])) }
G <- rbindlist(g); fwrite(G, file.path(OUT_DIR, sprintf("inla_mode_%s_grid_%s.csv", MODE, RUN_TAG)))
G[, d2 := shift(mlik, -1) - 2 * mlik + shift(mlik, 1)]
cat("\n==== Part 1 (healthy ratio ~1.0-1.2; collapse < 0.1) ====\n"); R <- rbindlist(out); print(R[, .(m, label, alpha = round(alpha, 3), sd = signif(sd, 3), ratio = round(ratio, 3), mlik = round(mlik, 1), error)])
cat("\ncollapsed:", sum(R$ratio < 0.1, na.rm = TRUE), " errors:", sum(!is.na(R$error)), " healthy:", sum(R$ratio >= 0.1, na.rm = TRUE), " of ", nrow(R), "\n")
cat("\n==== Part 2: m02 fixed-alpha grid (compact mode gave stuck/partial fits at many points; smooth would rise ~0.06 per step) ====\n")
print(G[, .(alpha, mlik = round(mlik, 3))]); cat("noise = SD(2nd diff)/sqrt(6):", signif(sd(G$d2, na.rm = TRUE) / sqrt(6), 3), " max |mlik + 2165|:", round(max(abs(G$mlik + 2165), na.rm = TRUE), 1), "\n")
