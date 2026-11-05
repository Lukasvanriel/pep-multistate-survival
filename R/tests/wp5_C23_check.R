# ==============================================================================
# R/tests/wp5_C23_check.R   (WP5: C2 exit rule, C3 drawn exit; no INLA, synthetic tables, runs anywhere)
#   Rscript R/tests/wp5_C23_check.R            (project root)
# Checks (a failed check prints what it has and stops; bounds are never changed -- report the output to the user first):
#   1  C3 unit test: 100,000 draws of the production make_p2_draw() vs the truncated Weibull CDF, written here in closed form
#      F(t) = (1 - exp(-(H(t)-H(a)))) / (1 - exp(-(H(b)-H(a)))), H(t) = (lam t)^alpha, independent of the production code.
#      5 cases (own seed each, 20261110 + case), KS p > 0.01. Case 5 has H(a) ~ 800 (S(a) underflows to 0).
#   2  uniform fallback (no usable fit) vs Uniform(a, b) (KS p > 0.01); non-finite closed form -> "uniform_nonfinite"
#   3  C2 rule: Stage 1, Stage 2 (moment) and draw return the same rule date, for ok / failed / non-kept transitions; clamp stop;
#      only X->regen; XREGEN_EXIT_RULE = FALSE switches it off; rule exits at the end visit are counted
#   4  one runif(1) per event: the stream of the other events is identical with the rule on and off
#   5  no D_heavy: make_p2_moment equals the pre-WP5 function, stage1_rows equals build_rows(p2_midpoint)
#   6  synthetic tables through prepare_fit_inputs + stage1_rows: every rule exit = t_start + D_heavy - year_start, rule count = X->regen with D_heavy,
#      exit_source complete; recount_events passes
#   7  run_stage1 / run_stage2 / run_one_imputation plumbing with stubbed fits (counts per source, diag columns, extras$p2$exit_source)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
options(warn = 1)
source(here::here("R", "config.R")); source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
fail <- function(msg, x = NULL) { cat("\n!! CHECK FAILED:", msg, "\n"); if (!is.null(x)) print(x, digits = 6)
  stop("check failed; bounds were NOT changed. Report this output to the user before changing anything.", call. = FALSE) }
ok <- function(cond, msg, x = NULL) if (!isTRUE(cond)) fail(msg, x)
KS_BOUND <- 0.01

# ---- synthetic tables (also used by R/tests/wp5_recount_stage1.R's cloud test through make_synth_dir) ----------------------------------------------
make_synth_dir <- function(n_plots = 400, seed = 7) {
  set.seed(seed); st <- c("regen", "a", "b", "c")
  rows <- list(); pin <- list()
  for (i in seq_len(n_plots)) {
    W <- runif(1, 20, 80); nint <- sample(2:4, 1); gap <- sample(5:12, nint, TRUE); ts <- W + c(0, cumsum(gap))[seq_len(nint)]
    s <- sample(st, 1); s_first <- s; ss <- character(nint); se <- character(nint)
    for (r in seq_len(nint)) { ss[r] <- s; if (runif(1) < 0.4) s <- sample(setdiff(st, s), 1); se[r] <- s }
    xr <- se == "regen" & ss != "regen"
    fire <- rbinom(nint, 1, ifelse(xr, 0.25, 0.08)); pest <- rbinom(nint, 1, 0.1)
    harv <- ifelse(xr, sample(0:2, nint, TRUE, c(.4, .1, .5)), sample(0:2, nint, TRUE, c(.8, .1, .1)))
    d <- data.table(id_pe = i, interval = seq_len(nint), t_start = ts, t_end = ts + gap, state_start = ss, state_end = se,
                    gdd_3 = rnorm(nint, 1500, 150), cmi_gs = rnorm(nint, 5, 3), soil_pc1 = rnorm(nint), fire = fire, harv = harv, pest = pest)
    d[, `:=`(year_start = t_start + 1900, year_end = t_end + 1900, is_change = state_start != state_end)]
    d[, `:=`(tsd_fire_y = fifelse(fire == 1L, runif(.N, 0, 1.6 * (t_end - t_start)), 0), tsd_pest_y = fifelse(pest == 1L, runif(.N, 0, 20), 0),
             tsd_harv_y = fifelse(harv > 0L, runif(.N, 0, 1.6 * (t_end - t_start)), 0))]
    rows[[i]] <- d
    pin[[i]] <- data.table(id_pe = i, W_j = W, first_state = s_first, first_year = 1980L, n_visits = nint + 1L, n_changes = sum(d$is_change),
                           fv_gdd_3 = rnorm(1, 1500, 150), fv_cmi_gs = rnorm(1, 5, 3), fv_soil_pc1 = rnorm(1), fv_fire = 0L, fv_harv = 0L, fv_pest = 0L,
                           fv_tsd_fire_y = 0, fv_tsd_pest_y = 0, fv_tsd_harv_y = 0)
  }
  iv <- rbindlist(rows)
  # A9, as in 04
  iv[, `:=`(D_harv_heavy = fifelse(harv == 2L & tsd_harv_y < year_end - year_start, year_end - tsd_harv_y, NA_real_),
            D_fire = fifelse(fire == 1L & tsd_fire_y < year_end - year_start, year_end - tsd_fire_y, NA_real_))]
  iv[, D_heavy := pmin(D_harv_heavy, D_fire, na.rm = TRUE)]; iv[is.infinite(D_heavy), D_heavy := NA_real_]
  # transition_set from the change intervals (end-visit flags)
  ch <- iv[is_change == TRUE]
  tr <- CJ(from = st, to = st)[from != to]
  cn <- ch[, .(n_events = .N, n_fire = sum(fire == 1L), n_pest = sum(pest == 1L), n_harv0 = sum(harv == 0L), n_harv1 = sum(harv == 1L), n_harv2 = sum(harv == 2L)), by = .(from = state_start, to = state_end)]
  tr <- merge(tr, cn, by = c("from", "to"), all.x = TRUE); for (v in setdiff(names(cn), c("from", "to"))) set(tr, which(is.na(tr[[v]])), v, 0L)
  tr[, `:=`(kept = n_events >= MIN_EVENTS, fire_term = n_fire >= MIN_EVENTS, pest_term = n_pest >= MIN_EVENTS)]
  tr[, harv_term := fifelse(n_harv1 + n_harv2 >= MIN_EVENTS, "3level", "none")]
  tr[kept == FALSE, `:=`(fire_term = FALSE, pest_term = FALSE, harv_term = "none")]
  tr[, n_undist := n_events - pmin(n_events, n_fire * fire_term + n_pest * pest_term + (n_harv1 + n_harv2) * (harv_term == "3level"))]
  d <- tempfile("wp5"); dir.create(d)
  saveRDS(iv, file.path(d, "visit_intervals.rds")); saveRDS(rbindlist(pin), file.path(d, "plot_info.rds")); saveRDS(tr, file.path(d, "transition_set.rds"))
  d
}

# ---- mock inputs for single-event calls ---------------------------------------------------------------------------------------------------------
mock_inp <- function(Dh = NA_real_, t_start = 40, t_end = 52, ys = 1990, n = 1L) {
  st <- c("regen", "a"); km <- matrix(NA_integer_, 2, 2, dimnames = list(st, st)); km["regen", "a"] <- 1L; km["a", "regen"] <- 2L
  iv <- data.table(id_pe = seq_len(n), interval = 1L, t_start = rep(t_start, length.out = n), t_end = rep(t_end, length.out = n), year_start = rep(ys, length.out = n),
                   year_end = rep(ys, length.out = n) + (t_end - t_start), D_heavy = rep(Dh, length.out = n))
  list(states = st, start_idx = 1L, kmat = km, Xiv = matrix(0, n, 1, dimnames = list(NULL, "z")), iv = iv)
}
mock_theta <- function(lam, alpha, ok = TRUE) list(ok = rep(ok, 2), intercept = rep(log(lam), 2), beta = matrix(0, 2, 1), alpha = rep(alpha, 2))

# ---- 1. C3 unit test ---------------------------------------------------------------------------------------------------------------------------
cases <- data.table(case = 1:5, alpha = c(0.7, 1.0, 2.5, 1.4, 2.0), lam = c(0.05, 0.03, 0.04, 0.02, 0.5657), a = c(0, 5, 10, 40, 50), b = c(20, 15, 30, 40.5, 51))
cdf_trunc <- function(t, a, b, lam, al) { H <- function(x) (lam * x)^al; -expm1(-(H(t) - H(a))) / -expm1(-(H(b) - H(a))) }   # written independently of the production code
res1 <- list()
for (i in seq_len(nrow(cases))) {
  cs <- cases[i]; inp <- mock_inp(); f <- make_p2_draw(mock_theta(cs$lam, cs$alpha), inp)
  set.seed(20261110 + cs$case)
  x <- vapply(seq_len(1e5), function(j) f(cs$a, cs$b, "regen", "a", 1L), 0)
  ks <- ks.test(x, function(t) cdf_trunc(t, cs$a, cs$b, cs$lam, cs$alpha))
  Ha <- (cs$lam * cs$a)^cs$alpha
  res1[[i]] <- data.table(case = cs$case, seed = 20261110 + cs$case, alpha = cs$alpha, lam = cs$lam, a = cs$a, b = cs$b, H_a = round(Ha, 1), min_x = min(x), max_x = max(x),
                          n_unique = uniqueN(x), ks_p = ks$p.value)
  ok(all(x >= cs$a & x <= cs$b), sprintf("case %d: draw outside [a, b]", cs$case))
  ok(ks$p.value > KS_BOUND, sprintf("C3 case %d: KS p = %.4g <= %g", cs$case, ks$p.value, KS_BOUND), res1[[i]])
}
res1 <- rbindlist(res1)
cat("1. C3 unit test (production make_p2_draw vs closed-form truncated Weibull CDF, 1e5 draws per case):\n"); print(res1, digits = 4)

# ---- 2. uniform fallback and non-finite ------------------------------------------------------------------------------------------------------
inp <- mock_inp(); dg <- new.env(); f <- make_p2_draw(mock_theta(0.03, 1.5, ok = FALSE), inp, dg)
set.seed(20261120); x <- vapply(1:1e5, function(j) f(8, 31, "regen", "a", 1L), 0)
ksu <- ks.test(x, function(t) (t - 8) / (31 - 8)); ok(ksu$p.value > KS_BOUND, sprintf("uniform fallback: KS p = %.4g", ksu$p.value))
ok(all(x > 8 & x < 31) && dg$n_event == 1e5 && identical(unique(dg$src[1L]), "uniform"), "uniform fallback: range / source")
f <- make_p2_draw(mock_theta(0.03, 1.5), inp, dg); x <- f(8, 31, "a", "regen", 1L)      # to == regen but D_heavy NA: draw (kmat a->regen = 2, theta ok)
ok(dg$src[1L] == "draw", "a->regen without D_heavy must be a draw")
dg <- new.env(); f <- make_p2_draw(mock_theta(0.9, 1000), inp, dg); set.seed(1); x <- f(8, 31, "regen", "a", 1L)
ok(dg$src[1L] == "uniform_nonfinite" && x > 8 && x < 31, "non-finite closed form (H(b) = Inf) must give a uniform_nonfinite draw on (a, b)", list(src = dg$src[1L], x = x))
inp3 <- mock_inp(); inp3$kmat["regen", "a"] <- NA_integer_; dg <- new.env(); f <- make_p2_draw(mock_theta(0.03, 1.5), inp3, dg); invisible(f(8, 31, "regen", "a", 1L))
ok(dg$src[1L] == "uniform", "non-kept transition (kmat NA) must be a uniform draw")
cat(sprintf("2. uniform fallback KS p = %.3f; non-finite -> uniform_nonfinite; non-kept -> uniform: OK\n", ksu$p.value))

# ---- 3. C2 rule ---------------------------------------------------------------------------------------------------------------------------------
# interval r=1: t 40-52, years 1990-2002. D = 1995.5 -> E = 40: a=0, b=12, T = 5.5 ; E = 30: a = 10, b = 22, T = 15.5
inp <- mock_inp(Dh = 1995.5)
thk <- mock_theta(0.03, 1.5); thbad <- mock_theta(0.03, 1.5, ok = FALSE)
for (E in c(40, 30, 12.5)) {
  a <- max(40 - E, 0); b <- 52 - E; want <- 40 + (1995.5 - 1990) - E
  got <- list(stage1 = make_p2_stage1(inp)(a, b, "a", "regen", 1L), moment_ok = make_p2_moment(thk, inp)(a, b, "a", "regen", 1L),
              moment_failed = make_p2_moment(thbad, inp)(a, b, "a", "regen", 1L), draw_ok = make_p2_draw(thk, inp)(a, b, "a", "regen", 1L),
              draw_failed = make_p2_draw(thbad, inp)(a, b, "a", "regen", 1L))
  inpn <- inp; inpn$kmat["a", "regen"] <- NA_integer_            # transition not kept
  got$draw_nonkept <- make_p2_draw(thk, inpn)(a, b, "a", "regen", 1L); got$moment_nonkept <- make_p2_moment(thk, inpn)(a, b, "a", "regen", 1L)
  ok(all(abs(unlist(got) - want) < 1e-12), sprintf("rule date for E = %g (want %.6f)", E, want), unlist(got))
}
dg <- new.env(); invisible(make_p2_moment(thk, inp, dg)(0, 12, "a", "regen", 1L)); ok(dg$src[1L] == "rule" && dg$n_event == 1L, "moment: source rule")
dg <- new.env(); invisible(make_p2_draw(thk, inp, dg)(0, 12, "a", "regen", 1L)); ok(dg$src[1L] == "rule", "draw: source rule")
# only X->regen: the same interval, destination "a" (from regen) gets no rule
dg <- new.env(); f <- make_p2_stage1(inp, dg); ok(f(0, 12, "regen", "a", 1L) == 6 && dg$src[1L] == "midpoint", "rule must only apply to X->regen (destination regen)")
# end visit: D = year_end -> T = b_dur, counted
inp_e <- mock_inp(Dh = 2002); dg <- new.env(); f <- make_p2_stage1(inp_e, dg); ok(abs(f(0, 12, "a", "regen", 1L) - 12) < 1e-12 && dg$n_rule_at_end == 1L, "D_heavy = year_end: T = b_dur, counted at the end")
# clamp: a date outside the bracket stops (A9 would be broken)
for (Dbad in c(2002.001, 1989.99)) { inp_b <- mock_inp(Dh = Dbad); e <- tryCatch(make_p2_stage1(inp_b)(0, 12, "a", "regen", 1L), error = function(e) conditionMessage(e))
  ok(is.character(e) && grepl("clamp moves the rule exit", e), sprintf("D_heavy = %g must stop (clamp)", Dbad), e) }
# switch off
XREGEN_EXIT_RULE <- FALSE
for (nm in c("stage1", "moment", "draw")) { dg <- new.env(); f <- switch(nm, stage1 = make_p2_stage1(inp, dg), moment = make_p2_moment(thk, inp, dg), draw = make_p2_draw(thk, inp, dg))
  invisible(f(0, 12, "a", "regen", 1L)); ok(dg$src[1L] != "rule", sprintf("XREGEN_EXIT_RULE = FALSE: %s still used the rule", nm)) }
XREGEN_EXIT_RULE <- TRUE
# missing D_heavy column stops with a clear message
inp_nc <- inp; inp_nc$iv[, D_heavy := NULL]; e <- tryCatch(make_p2_stage1(inp_nc), error = function(e) conditionMessage(e)); ok(grepl("rerun R/04", e), "missing D_heavy must stop at construction", e)
cat("3. C2 rule: same date in Stage 1 / Stage 2 / draw for ok, failed and non-kept transitions (E = 40, 30, 12.5); X->regen only; end visit counted; clamp stops; switch off: OK\n")

# ---- 4. one runif per event -----------------------------------------------------------------------------------------------------------------------
inp2 <- mock_inp(Dh = c(1995.5, NA), n = 2); th <- mock_theta(0.03, 1.5)
run2 <- function() { set.seed(99); f <- make_p2_draw(th, inp2); c(f(0, 12, "a", "regen", 1L), f(0, 12, "regen", "a", 2L)) }
on <- run2(); XREGEN_EXIT_RULE <- FALSE; off <- run2(); XREGEN_EXIT_RULE <- TRUE
ok(abs(on[1] - 5.5) < 1e-12 && off[1] != on[1] && on[2] == off[2], "rule on/off: the second event must get the same draw", list(on = on, off = off))
set.seed(5); f <- make_p2_draw(th, inp2); for (j in 1:7) invisible(f(0, 12, "regen", "a", 2L)); s1 <- .Random.seed; set.seed(5); invisible(runif(7)); ok(identical(s1, .Random.seed), "exactly one runif(1) per event")
cat("4. one runif(1) per event; the stream of the other events is the same with the rule on and off: OK\n")

# ---- 5/6. synthetic tables, backward compatibility ------------------------------------------------------------------------------------------------
# pre-WP5 make_p2_moment (05a at be794e5), copied literally
old <- new.env(); old$moment_match_weibull_naive <- moment_match_weibull_naive
old$make_p2_moment <- function(theta, inp, diag = new.env()) {
  diag$n_event <- 0L; diag$n_midpoint <- 0L
  function(a_dur, b_dur, from, to, r) {
    diag$n_event <- diag$n_event + 1L
    k <- inp$kmat[from, to]
    if (is.na(k) || !theta$ok[k]) { diag$n_midpoint <- diag$n_midpoint + 1L; return((a_dur + b_dur) / 2) }
    lam <- exp(theta$intercept[k] + sum(theta$beta[k, ] * inp$Xiv[r, ]))
    moment_match_weibull_naive(a_dur, b_dur, theta$alpha[k], lam)
  }
}
d <- make_synth_dir(); inp <- prepare_fit_inputs(d)
K <- nrow(inp$td); P <- length(inp$terms); set.seed(3)
theta <- list(alpha = runif(K, 1.1, 2.1), intercept = log(1 / runif(K, 25, 70)), beta = matrix(rnorm(K * P, 0, 0.25), K, P, dimnames = list(inp$td$label, inp$terms)), ok = rep(TRUE, K))
theta$beta[!inp$mask] <- 0; theta$ok[c(2, 5)] <- FALSE
# 5a  no D_heavy: new moment == old moment
inp_nd <- inp; inp_nd$iv <- copy(inp$iv)[, D_heavy := NA_real_]
fn <- make_p2_moment(theta, inp_nd); fo <- old$make_p2_moment(theta, inp_nd)
chg <- which(inp$iv$is_change); worst <- 0
for (r in chg) { a <- runif(1, 0, 5); b <- a + runif(1, 1, 12); x <- fn(a, b, inp$iv$state_start[r], inp$iv$state_end[r], r); y <- fo(a, b, inp$iv$state_start[r], inp$iv$state_end[r], r); worst <- max(worst, abs(x - y)) }
ok(worst == 0, "make_p2_moment without D_heavy differs from the pre-WP5 function", worst)
# 5b  stage1_rows without D_heavy == build_rows(p2_midpoint)
s1n <- stage1_rows(inp_nd); fe <- setNames(inp$p1$W_j / 2, inp$p1$id_pe); b0 <- build_rows(inp$iv, fe, inp$states, inp$tr, keep_cols = inp$keep_cols)
ok(identical(s1n$b$rows, b0$rows) && s1n$b$n_floored == b0$n_floored, "stage1_rows without D_heavy must equal the old Stage 1 build_rows call")
cat(sprintf("5. without D_heavy: make_p2_moment identical to the pre-WP5 function on %d events; stage1_rows rows identical to build_rows(p2_midpoint) (%d rows)\n", length(chg), nrow(b0$rows)))
# 6  with D_heavy
s1 <- stage1_rows(inp); p2 <- s1$p2
xr <- inp$iv[is_change == TRUE & state_end == "regen" & !is.na(D_heavy)]
rule <- p2[exit_source == "rule"]
ok(nrow(rule) == nrow(xr) && nrow(xr) > 20, "number of rule exits != number of X->regen intervals with D_heavy", c(rule = nrow(rule), expected = nrow(xr)))
mm <- merge(rule, inp$iv[, .(id_pe, interval, t_start, year_start, D_heavy)], by = c("id_pe", "interval"))
dev <- max(abs(mm$exit_abs - (mm$t_start + mm$D_heavy - mm$year_start))); ok(dev < 1e-9, "a rule exit differs from t_start + D_heavy - year_start", dev)
ok(nrow(p2) == sum(inp$iv$is_change) && all(p2$exit_source %in% c("rule", "midpoint")) && sum(p2$exit_source == "midpoint") == nrow(p2) - nrow(rule), "Stage 1 exit_source")
ok(s1$dg$n_rule_at_end == sum(abs(mm$exit_abs - (inp$iv$t_end[match(paste(mm$id_pe, mm$interval), paste(inp$iv$id_pe, inp$iv$interval))])) < 1e-9), "rule_at_end count")
rc <- recount_events(s1$b$rows, inp$tr)
cat(sprintf("6. synthetic tables (%d plots, %d change intervals): %d rule exits = X->regen intervals with D_heavy; max |exit - rule date| %.1e; Stage 1 floored segments %d; recount_events OK (%d kept transitions)\n",
            uniqueN(inp$iv$id_pe), sum(inp$iv$is_change), nrow(rule), dev, s1$b$n_floored, nrow(rc)))
print(dcast(p2_source_counts(p2), transition ~ exit_source, value.var = "n", fill = 0)[1:4])

# ---- 7. run_stage1 / run_stage2 / run_one_imputation plumbing with stubbed fits -----------------------------------------------------------------------
fit_rows <- function(rows, inp, fit_fun, only_k = NULL) { ks <- if (is.null(only_k)) integer(0) else as.integer(only_k); list(models = list(), n_events = integer(0), data = list(), ks = integer(0)) }
extract_theta <- function(models, inp, stat = "mean") theta
tidy_fit <- function(models, inp, ks = integer(0)) data.table(k = integer(0), parameter = character(0), fit_error = logical(0), sd = numeric(0))
add_cap_columns <- function(res, data_list, ks, inp) res
p1_entries <- function(theta, inp, M, mode, n_cores) data.table(id_pe = inp$p1$id_pe, entry = inp$p1$W_j / 2, n_cand = 1L, fallback = FALSE)
draw_theta_m <- function(models, inp, theta_seed, sampler) c(theta, list(intercept_dec = rep(0, K), beta_v0 = matrix(0, K, P, dimnames = list(NULL, inp$terms)), n_fallback = 0L))
t1 <- run_stage1(inp, NULL); ok(t1$p2_src[exit_source == "rule", sum(n)] == nrow(rule) && t1$p2_rule_at_end == s1$dg$n_rule_at_end, "run_stage1 p2_src")
t2 <- run_stage2(inp, theta, NULL)
ok(t2$p2_src[exit_source == "rule", sum(n)] == nrow(rule) && all(t2$p2$exit_source %in% c("rule", "moment", "midpoint")) && t2$p2_midpoint == sum(t2$p2$exit_source == "midpoint") &&
     t2$p2_event == nrow(t2$p2), "run_stage2 sources / counts")
set.seed(11); im <- run_one_imputation(1L, inp, list(), theta, NULL)
ex <- im$extras$p2; ok(all(c("exit_source", "exit_abs", "floored") %in% names(ex)) && nrow(ex) == nrow(rule) + sum(ex$exit_source != "rule"), "extras$p2 columns")
ok(all(ex$exit_source %in% c("rule", "draw", "uniform", "uniform_nonfinite")), "imputation sources", table(ex$exit_source))
dgn <- im$diag; ok(dgn$p2_rule == nrow(rule) && dgn$p2_rule + dgn$p2_draw + dgn$p2_uniform + dgn$p2_uniform_nonfinite == dgn$p2_events && dgn$p2_midpoint == 0L, "mi_diag columns", dgn)
ok(all(c("m", "transition", "exit_source", "n") %in% names(im$p2_src)) && all(im$p2_src$m == 1L), "p2_src per imputation")
# reproducible with the seed
set.seed(11); im2 <- run_one_imputation(1L, inp, list(), theta, NULL); ok(identical(im2$extras$p2, ex), "imputation P2 not reproducible with the seed")
set.seed(12); im3 <- run_one_imputation(1L, inp, list(), theta, NULL); nr <- ex$exit_source != "rule"
ok(all(ex$exit_abs[!nr] == im3$extras$p2$exit_abs[!nr]) && sum(ex$exit_abs[nr] != im3$extras$p2$exit_abs[nr]) > 0.9 * sum(nr), "rule exits fixed, other exits vary between seeds")
cat("7. plumbing: Stage 1 / Stage 2 / imputation counts per source, diag columns, extras$p2$exit_source, reproducible; rule exits identical across seeds, others vary\n")
cat("   Stage 1 :", p2_source_line(t1$p2_src), "\n   Stage 2 :", p2_source_line(t2$p2_src), "\n   imputation:", p2_source_line(im$p2_src), "\n")
cat("\nALL CHECKS PASSED\n")
