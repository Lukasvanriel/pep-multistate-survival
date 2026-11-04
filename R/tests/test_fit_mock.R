# Mock test of the fitting machinery WITHOUT INLA: fake model objects and a fake posterior sampler.
suppressPackageStartupMessages(library(data.table))
proc_dir <- Sys.getenv("PEP_PROC_DIR"); r_dir <- Sys.getenv("PEP_R_DIR")
source(file.path(r_dir, "covariate_spec.R")); source(file.path(r_dir, "04b_build_rows.R")); source(file.path(r_dir, "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG"); set.seed(7)
N <- as.integer(Sys.getenv("PEP_N", "300"))
inp <- prepare_fit_inputs(proc_dir, subsample_n = N)
cat(sprintf("plots %d, intervals %d, kept transitions %d, terms %d\n", nrow(inp$p1), nrow(inp$iv), nrow(inp$td), length(inp$terms)))

# mock INLA: plausible-ish theta per transition, structure of summary.fixed / summary.hyperpar
set.seed(11); K <- nrow(inp$td)
true_a <- runif(K, 1, 1.6); true_i <- rnorm(K, log(0.03*10), 0.3)   # decade-scale intercept
mock_fit <- function(data_list, rhs_formula, exact = TRUE, ...) lapply(seq_along(data_list), function(k) {
  d <- data_list[[k]]; if (nrow(d) == 0) return(list(error = "no rows"))
  if (k %% 17 == 0) return(list(error = "mock failure"))       # a few failed fits
  tms <- colnames(inp$mask)[inp$mask[k, ]]
  fx <- cbind(mean = c(true_i[k], rnorm(length(tms), 0, .1)), sd = 0.05); rownames(fx) <- c("(Intercept)", tms)
  hy <- cbind(mean = true_a[k] + rnorm(1, 0, .02), sd = 0.03); rownames(hy) <- ALPHA_NAME
  list(summary.fixed = fx, summary.hyperpar = hy, nrow = nrow(d))
})
mock_sampler <- function(n, m, seed = 0L) { fx <- m$summary.fixed
  lapply(seq_len(n), function(i) { lat <- matrix(rnorm(nrow(fx), fx[, "mean"], fx[, "sd"]), ncol = 1, dimnames = list(paste0(rownames(fx), ":1"), NULL))
  list(latent = lat, hyperpar = setNames(rnorm(1, m$summary.hyperpar[1, "mean"], .03), ALPHA_NAME)) }) }

t0 <- proc.time()["elapsed"]
s1 <- run_stage1(inp, mock_fit)
cat(sprintf("Stage 1: %.1fs | theta ok %d/%d | floored segs %d\n", proc.time()["elapsed"] - t0, sum(s1$theta$ok), K, s1$n_floored))
t0 <- proc.time()["elapsed"]
s2 <- run_stage2(inp, s1$theta, mock_fit, M = SIM_M, n_cores = 2L)
cat(sprintf("Stage 2: %.1fs | P1 fallback %d of %d plots (no candidate: %d) | P2 events %d, midpoint %d | floored %d\n",
            proc.time()["elapsed"] - t0, sum(s2$p1$fallback), nrow(s2$p1), sum(s2$p1$n_cand == 0), s2$p2_event, s2$p2_midpoint, s2$n_floored))
cat("entry/W_j (P1 mean) summary among non-fallback plots:\n"); print(summary((s2$p1$entry / inp$p1$W_j)[!s2$p1$fallback & inp$p1$W_j > 0]))
raw <- list(); dg <- list()
for (m in 1:2) { r <- run_one_imputation(m, inp, s2$models, s2$theta, mock_fit, M = SIM_M, sampler = mock_sampler); raw[[m]] <- r$res; dg[[m]] <- r$diag; ex <- r$extras }
N_COV_SAMPLES <- 200L
cat("extras: p1", nrow(ex$p1), "p2", nrow(ex$p2), "theta_m", dim(ex$theta_m), "\n"); e1 <- ex$fits[[1]]; str(e1[c("summary_fixed","cov")], max.level = 2); print(round(e1$cov$cov, 4)); print(sapply(ex$fits, function(e) !is.null(e$error)) |> sum())
print(rbindlist(dg))
pooled <- pool_all(rbindlist(raw)); print(head(pooled, 6)); cat("pooled rows:", nrow(pooled), "\n")
