# ==============================================================================
# 13_test_oldregen_exclusion.R   RUN ON THE MACHINE WITH INLA (~15-30 min with 6 workers)   -- DIAGNOSTIC ONLY, changes nothing in the pipeline
#
# QUESTION: do the INLA failures of the regen exits (regen->fir 18/128, regen->red_maple 6/128 invalid in the "int" run) disappear when the
#   "old regen" plots are excluded?  Old regen = first state regen with stand age W_j >= WJ_MIN (default 50): BA < 5 m2/ha but an old age class,
#   no stand-origin disturbance recorded (405/458), height class 3-4, density class D -> open mature stands, not regeneration. They create regen
#   stays of 70-140 years (left truncation up to 16 decades), the largest truncation of any regen exit.
#
# DESIGN (re-uses the stored "int" imputations, like script 11; no new imputation):
#   * imputations: every m in which one of LABELS failed the screen in the int run + N_CTRL imputations in which all LABELS were valid
#   * for each m: rows rebuilt from extras_int_mXX.rds (P1 entries, P2 exits); variant "full" = as fitted; variant "excl" = the same rows
#     without the old-regen plots (all their rows dropped)
#   * each LABEL fitted with the pipeline's INLA call (weibullsurv variant 1, adaptive, grid, VB off) in both variants
#   * CORRECTNESS CHECK independent of INLA: maximum-likelihood fit of the same Weibull likelihood (split rows, left truncation) on the same rows.
#     A fit is WRONG if |alpha_INLA - alpha_MLE| > 3 x MLE SE (healthy fits sit at about -1 SE because of the default PC prior on alpha);
#     the pipeline screen (alpha-SD ratio < 0.1 or error) is reported alongside.
#   * step 0: rows of one valid imputation x 2 labels must reproduce mi_raw_int.csv (stops otherwise)
# MUST be run BEFORE 04_prepare_inla.R is re-run for the gs run (it needs the int-version visit_intervals.rds with cmi_sum; checked below).
# OUTPUT: data-output/inttest/fits/test_oldregen_fits.csv (one row per m x label x variant) + printed summary.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 220)
RUN_TAG <- "int"; IN_DIR <- here("data-output", "inttest", "fits")
WJ_MIN <- 50L; N_CTRL <- 10L; N_WORKERS <- 6L
LABELS <- c("regen->fir", "regen->red_maple", "regen->black_spruce", "regen->white_birch", "regen->trembling_aspen", "regen->jack_pine")
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
iv0 <- readRDS(here("data", "tables", "inttest", "visit_intervals.rds"))
if (!"cmi_sum" %in% names(iv0)) stop("visit_intervals.rds is already the gs version: run this test before re-running 04 (or restore the int tables).")
inp <- prepare_fit_inputs(here("data", "tables", "inttest"))
stopifnot(all(LABELS %in% inp$td$label))

# ---- old-regen plots ---------------------------------------------------------------------------------------------------------------------------------
pinfo <- as.data.table(readRDS(here("data", "tables", "inttest", "plot_info.rds")))
old <- pinfo[first_state == "regen" & W_j >= WJ_MIN, as.character(id_pe)]
cat(sprintf("Old-regen plots (first state regen, W_j >= %d): %d of %d regen-first plots\n", WJ_MIN, length(old), pinfo[first_state == "regen", .N]))

# ---- rows of a stored imputation -----------------------------------------------------------------------------------------------------------------------
key_iv <- paste(inp$iv$id_pe, inp$iv$interval); b_abs <- inp$iv$t_end
rebuild_rows <- function(m) {
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(inp$iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  r <- build_rows(inp$iv, setNames(p1$entry, as.character(p1$id_pe)), inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
  r[, label := paste(from, to, sep = "->")][label %in% LABELS] }

# ---- INLA call (identical settings to model_fitting.R, exact = TRUE) and the independent MLE ----------------------------------------------------------
fit_inla <- function(d, rhs) tryCatch(inla(as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", rhs)), family = "weibullsurv", data = d,
  control.family = list(variant = 1), control.compute = list(config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE))), error = function(e) list(error = conditionMessage(e)))
fit_mle <- function(d, tms) {                                   # S(t) = exp(-(lambda t)^alpha), lambda = exp(eta); theta = (log alpha, intercept, betas)
  X <- cbind(1, as.matrix(d[, tms, drop = FALSE])); T <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(T); ltr <- ifelse(tr > 0, log(pmax(tr, 1e-300)), 0)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + a * eta + (a - 1) * lT)) - sum(exp(a * eta) * (T^a - tr^a))) }
  gr <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); e <- exp(a * eta); Ta <- T^a; tra <- tr^a
    gb <- -(colSums(X * (s * a)) - colSums(X * (a * e * (Ta - tra)))); ga <- -(sum(s * (1 / a + eta + lT)) - sum(e * (eta * (Ta - tra) + Ta * lT - tra * ltr))); c(ga * a, gb) }
  o <- optim(c(0, log(sum(s) / sum(T - tr)), rep(0, length(tms))), nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll, gr)), error = function(e) NULL); a <- exp(o$par[1])
  list(alpha = a, alpha_se = if (is.null(V)) NA_real_ else a * sqrt(V[1, 1]), icpt = o$par[2], conv = o$convergence) }
one <- function(d, lb, m, variant) {
  k <- match(lb, inp$td$label); tms <- colnames(inp$mask)[inp$mask[k, ]]; dd <- as.data.frame(d[, c("status", "Tstop", "truncation", tms), with = FALSE])
  mo <- fit_inla(dd, inp$rhs[[k]]); ml <- fit_mle(dd, tms); ev <- sum(dd$status == 1L)
  a <- if (is.null(mo$error)) mo$summary.hyperpar[ALPHA_NAME, ] else NULL
  data.table(m = m, label = lb, variant = variant, rows = nrow(dd), events = ev, inla_error = !is.null(mo$error),
             alpha = if (is.null(a)) NA_real_ else a$mean, alpha_sd = if (is.null(a)) NA_real_ else a$sd,
             icpt = if (is.null(mo$error)) mo$summary.fixed["(Intercept)", "mean"] else NA_real_,
             mle_alpha = ml$alpha, mle_alpha_se = ml$alpha_se, mle_icpt = ml$icpt, mle_conv = ml$conv, trunc_max = max(dd$truncation)) }

# ---- step 0: reproduction ----------------------------------------------------------------------------------------------------------------------------------
fv <- fread(file.path(IN_DIR, sprintf("fit_validity_%s.csv", RUN_TAG)))[label %in% LABELS]
mi <- fread(file.path(IN_DIR, sprintf("mi_raw_%s.csv", RUN_TAG)))
m0 <- fv[m <= M_MI, .(all_ok = all(valid)), by = m][all_ok == TRUE, min(m)]
r0 <- rebuild_rows(m0)
chk <- rbindlist(lapply(LABELS[1:2], function(lb) one(r0[label == lb], lb, m0, "full")))
chk <- merge(chk, mi[m == m0 & parameter == "alpha", .(label, est, sd)], by = "label"); chk[, dz := abs(alpha - est) / sd]
cat(sprintf("Step 0 (m = %d): |alpha - mi_raw| / sd = %s\n", m0, paste(round(chk$dz, 4), collapse = ", ")))
if (any(is.na(chk$dz) | chk$dz > 0.1)) stop("Rows / INLA call do not reproduce the int run; stop.")

# ---- imputations to test ------------------------------------------------------------------------------------------------------------------------------------
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
fail_m <- intersect(fv[valid == FALSE, sort(unique(m))], avail)
ctrl_m <- head(setdiff(intersect(fv[, .(all_ok = all(valid)), by = m][all_ok == TRUE, m], avail), fail_m), N_CTRL)
M_SET <- sort(c(fail_m, ctrl_m))
cat(sprintf("Imputations: %d with an int failure among LABELS (%s) + %d controls\n", length(fail_m), paste(fail_m, collapse = ","), length(ctrl_m)))

res <- rbindlist(mclapply(M_SET, function(m) tryCatch({
  r <- rebuild_rows(m); rx <- r[!(as.character(id_pe) %in% old)]
  rbindlist(lapply(LABELS, function(lb) rbind(one(r[label == lb], lb, m, "full"), one(rx[label == lb], lb, m, "excl")))) },
  error = function(e) data.table(m = m, label = NA_character_, variant = "run_error", inla_error = TRUE)), mc.cores = N_WORKERS), fill = TRUE)
res <- merge(res, fv[, .(m, label, int_valid = valid, int_reason = reason)], by = c("m", "label"), all.x = TRUE)
res[, `:=`(ratio = alpha_sd / (0.78 * alpha / sqrt(events)), z_mle = (alpha - mle_alpha) / mle_alpha_se)]
res[, `:=`(screen_fail = inla_error | is.na(alpha) | ratio < 0.1, wrong = inla_error | is.na(alpha) | abs(z_mle) > 3)]
fwrite(res, file.path(IN_DIR, "test_oldregen_fits.csv"))

# ---- summary ------------------------------------------------------------------------------------------------------------------------------------------------
cat("\n== Per transition and variant: failures (screen) and WRONG fits (INLA alpha > 3 MLE-SE from the MLE) ==\n")
print(res[!is.na(label), .(fits = .N, screen_fail = sum(screen_fail), wrong = sum(wrong), events = round(mean(events)), trunc_max_dec = round(max(trunc_max), 1),
                           mle_alpha = round(median(mle_alpha), 3), inla_alpha_ok = round(median(alpha[!wrong]), 3), median_z_ok = round(median(z_mle[!wrong]), 2)),
          by = .(label, variant)][order(label, variant)])
cat("\n== Imputations that FAILED in the int run: does excluding the old-regen plots repair them? ==\n")
print(dcast(res[int_valid == FALSE, .(label, m, variant, wrong)], label + m ~ variant, value.var = "wrong")[, .(label, m, full_wrong = full, excl_wrong = excl)])
cat("\n== Shift of the estimates (healthy fits), excl vs full: median MLE alpha and events ==\n")
print(res[!is.na(label) & wrong == FALSE, .(mle_alpha = round(median(mle_alpha), 3), mle_icpt = round(median(mle_icpt), 3), events = round(mean(events))), by = .(label, variant)][order(label, variant)])
