# ==============================================================================
# 06a_check_sampler.R   RUN ON THE MACHINE WITH INLA (about a minute, no refitting)
#
# Does the seeded posterior sampler reproduce?
#  (a) same seed twice on one seed model -> identical draw? different seed -> different draw?
#  (b) imputation 1 of the "diag" run: set.seed(SEED*1000 + 1), draw_theta_m() on the seed models
#      (stage2.rds) -> equal to theta_m stored in extras_diag_m01.rds?
# Reads only; writes nothing.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) })
SEED <- 2026L; RUN_TAG <- "diag"; OUT_DIR <- here("data", "processed", "fit")
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
RNGkind("L'Ecuyer-CMRG")
inp <- prepare_fit_inputs(here("data", "processed"), subsample_n = NULL, seed = SEED)
s2  <- readRDS(file.path(OUT_DIR, "stage2.rds"))

cat("(a) one seed model, three draws\n")
k <- which(vapply(s2$models, function(m) is.null(m$error), logical(1)))[1]
d1 <- default_sampler(1, s2$models[[k]], seed = 123L)[[1]]
d2 <- default_sampler(1, s2$models[[k]], seed = 123L)[[1]]
d3 <- default_sampler(1, s2$models[[k]], seed = 124L)[[1]]
cat("   same seed, latent identical:   ", isTRUE(all.equal(d1$latent, d2$latent, tolerance = 0)), "\n")
cat("   same seed, hyperpar identical: ", isTRUE(all.equal(unlist(d1$hyperpar), unlist(d2$hyperpar), tolerance = 0)), "\n")
cat("   different seed, latent differs:", !isTRUE(all.equal(d1$latent, d3$latent)), "\n")

cat("\n(b) imputation 1: parameter draw vs stored theta_m\n")
set.seed(SEED * 1000L + 1L)
th <- draw_theta_m(s2$models, inp, s2$theta)
ex <- readRDS(file.path(OUT_DIR, sprintf("extras_%s_m01.rds", RUN_TAG)))$theta_m
bcols <- colnames(th$beta)
cmp <- list(alpha = max(abs(th$alpha - ex$alpha), na.rm = TRUE),
            intercept = max(abs(th$intercept - ex$intercept_years), na.rm = TRUE),
            beta = max(abs(th$beta - as.matrix(ex[, ..bcols])), na.rm = TRUE))
print(unlist(cmp))
cat(if (max(unlist(cmp)) < 1e-8) "   RESULT: reproduces (max abs difference < 1e-8)\n" else "   RESULT: does NOT reproduce\n")
cat("   n_fallback this draw:", th$n_fallback, "\n")
