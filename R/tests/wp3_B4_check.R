# ==============================================================================
# R/tests/wp3_B4_check.R   (WP3, B4 unit test)
# S(t) from the variant-0 parameters, S0 = exp(-exp(b0 + x'beta) (t/10)^alpha), must equal S(t) from the converted P1/P2 form,
# S1 = exp(-(lambda1 t)^alpha), lambda1 from lambda_matrix() on the theta that extract_theta() / draw_theta_m() return. Tolerance 1e-10.
# 3 random covariate rows x t = 5, 20, 50 y x every transition.
#   Rscript R/tests/wp3_B4_check.R          mock fits (no INLA; runs anywhere): extract_theta and draw_theta_m with a fake sampler
#   Rscript R/tests/wp3_B4_check.R inla     ALSO one real INLA variant-0 fit on simulated data and the real inla.posterior.sample
#                                           (needs INLA; run on the Mac). Stops on failure; prints ALL CHECKS PASSED.
# Run from the project root.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("R", "config.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
stopifnot("test is for Weibull variant 0" = WEIBULL_VARIANT == 0L)
RUN_INLA <- "inla" %in% commandArgs(TRUE); TOL <- 1e-10; TT <- c(5, 20, 50)
set.seed(20261104)

# S from raw variant-0 values (theta0) and from the converted theta, for given covariate rows X (n x P) and transition k
S_both <- function(th, X, k, a0, b00, beta0) {
  eta <- b00 + drop(X %*% beta0)                                       # variant-0 linear predictor, decade scale
  S0 <- sapply(TT, function(t) exp(-exp(eta) * (t / 10)^a0))             # n x |TT|
  lam <- lambda_matrix(th, X)[, k]
  S1 <- sapply(TT, function(t) exp(-(lam * t)^a0))
  list(S0 = S0, S1 = S1)
}
compare <- function(th, inp, tag) {
  worst <- 0
  for (k in which(th$ok)) {
    X <- matrix(rnorm(3 * length(inp$terms)), 3, dimnames = list(NULL, inp$terms))
    r <- S_both(th, X, k, th$alpha[k], th$intercept_dec[k], th$beta_v0[k, ])
    worst <- max(worst, abs(r$S0 - r$S1))
  }
  cat(sprintf("%-28s transitions %d | rows 3 x t %s | max |S0 - S1| = %.2e\n", tag, sum(th$ok), paste(TT, collapse = "/"), worst))
  stopifnot(worst < TOL)
}

# ---- mock fits: K transitions with different term sets ---------------------------------------------------------------------------
terms <- c("z_gdd_3", "z_cmi_gs", "z_soil_pc1", "fire_f", "harv_light"); K <- 4L
mask <- rbind(c(T,T,T,F,F), c(T,T,T,T,F), c(T,T,T,F,T), c(F,T,T,T,T)); dimnames(mask) <- list(paste0("s", 1:K, "->r"), terms)
inp <- list(td = data.table(label = rownames(mask)), mask = mask, terms = terms)
truth <- lapply(seq_len(K), function(k) list(a = runif(1, 0.6, 2.5), b0 = rnorm(1, -3, 1), beta = ifelse(mask[k, ], rnorm(length(terms), 0, 0.7), 0)))
mk_model <- function(tr, k) {
  nm <- colnames(mask)[mask[k, ]]
  list(summary.hyperpar = matrix(tr$a, 1, 1, dimnames = list(ALPHA_NAME, "mean")),
       summary.fixed = matrix(c(tr$b0, tr$beta[nm]), ncol = 1, dimnames = list(c("(Intercept)", nm), "mean")))
}
models <- lapply(seq_len(K), function(k) mk_model(truth[[k]], k))
th1 <- extract_theta(models, inp)
stopifnot(all(th1$ok), all(abs(th1$alpha - sapply(truth, `[[`, "a")) < 1e-12), all(abs(th1$intercept_dec - sapply(truth, `[[`, "b0")) < 1e-12))
compare(th1, inp, "extract_theta (means)")
# conversion formula spelled out
for (k in seq_len(K)) stopifnot(abs(th1$intercept[k] - (truth[[k]]$b0 / truth[[k]]$a - log(10))) < 1e-12, all(abs(th1$beta[k, ] - truth[[k]]$beta / truth[[k]]$a) < 1e-12))

fake_sampler <- function(n, m, seed = 0L) {                       # one posterior draw around the model's means; own alpha per draw
  nm <- setdiff(rownames(m$summary.fixed), "(Intercept)")
  lat <- matrix(c(m$summary.fixed["(Intercept)", "mean"] + rnorm(1, 0, 0.2), m$summary.fixed[nm, "mean"] + rnorm(length(nm), 0, 0.2)), ncol = 1,
                dimnames = list(c("(Intercept)", nm), NULL))
  list(list(latent = lat, hyperpar = setNames(list(m$summary.hyperpar[ALPHA_NAME, "mean"] * exp(rnorm(1, 0, 0.1))), ALPHA_NAME)))
}
th2 <- draw_theta_m(models, inp, th1, fake_sampler)
stopifnot(th2$n_fallback == 0L, all(th2$alpha != th1$alpha), all(th2$intercept_dec != th1$intercept_dec))   # really a draw, stored raw
compare(th2, inp, "draw_theta_m (draws)")
# a failed fit falls back to the seed (converted AND raw values from the seed)
models_f <- models; models_f[[2]] <- list(error = "boom")
th3 <- draw_theta_m(models_f, inp, th1, fake_sampler)
stopifnot(th3$n_fallback == 1L, th3$intercept[2] == th1$intercept[2], th3$intercept_dec[2] == th1$intercept_dec[2], identical(th3$beta_v0[2, ], th1$beta_v0[2, ]))
compare(th3, inp, "draw_theta_m (with fallback)")

# ---- real INLA (Mac) ---------------------------------------------------------------------------------------------------------------------
if (RUN_INLA) {
  suppressPackageStartupMessages(library(INLA)); INLA::inla.setOption(num.threads = "1:1")
  n <- 3000; a_true <- 1.4; b0_true <- -3.5; bt <- c(z_gdd_3 = 0.5, z_cmi_gs = -0.3)
  X <- cbind(z_gdd_3 = rnorm(n), z_cmi_gs = rnorm(n)); eta <- b0_true + drop(X %*% bt)
  Tt <- (-log(runif(n)) / exp(eta))^(1 / a_true)                       # variant 0 on the decade scale: S = exp(-e^eta t^alpha), t in decades
  cens <- runif(n, 0.5, 3); d <- data.frame(Tstop = pmin(Tt, cens), status = as.integer(Tt <= cens), truncation = 0, X)
  d$Tstop <- pmax(d$Tstop, 1e-4)
  mods <- fit_inla_models_pc(list(d), "z_gdd_3 + z_cmi_gs", exact = TRUE)
  stopifnot(is.null(mods[[1]]$error))
  cat("INLA call: variant", mods[[1]]$.args$control.family$variant, "| control.fixed prec", mods[[1]]$.args$control.fixed$prec, "\n")
  stopifnot(mods[[1]]$.args$control.family$variant == 0L)
  inp1 <- list(td = data.table(label = "a->b"), mask = matrix(TRUE, 1, 2, dimnames = list("a->b", names(bt))), terms = names(bt))
  th <- extract_theta(mods, inp1)
  cat(sprintf("fitted: alpha %.3f (true %.2f), intercept %.3f (true %.2f), betas %s (true %s)\n", th$alpha, a_true, th$intercept_dec, b0_true,
              paste(round(th$beta_v0, 3), collapse = "/"), paste(bt, collapse = "/")))
  stopifnot(abs(th$alpha - a_true) < 0.15, abs(th$intercept_dec - b0_true) < 0.3)   # the fit recovers the variant-0 parameters (loose)
  compare(th, inp1, "INLA extract_theta")
  thd <- draw_theta_m(mods, inp1, th)                                        # real inla.posterior.sample: names, hyperpar, per-draw alpha
  stopifnot(thd$n_fallback == 0L); compare(thd, inp1, "INLA draw_theta_m")
}
cat("\nALL CHECKS PASSED\n")
