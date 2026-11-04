# ==============================================================================
# 14_test_station_disturbances.R   RUN ON THE MACHINE WITH INLA (~20-40 min with 6 workers)   -- DIAGNOSTIC ONLY, changes nothing in the pipeline
#
# QUESTION: are the INLA failures of black_spruce->regen (and the other X->regen exits) caused by disturbances that are recorded ONLY in
#   station_pe? 03_build_covariates.R reads origine/perturb from pee_ori_sond only. Of the 207 black_spruce->regen events the pipeline sees as
#   undisturbed, 121 have a stand-replacing code in station_pe. Station-only codes that are NEW at a visit coincide with >= 25% basal-area loss
#   in 74-91% of cases (9% for visits without any code), so they are real events.
#
# DESIGN (re-uses the stored "int" imputations, like scripts 11 and 13; no new imputation):
#   variant "pipeline": rows as fitted in the int run
#   variant "station" : the visit-interval table gets the station-only, NEW origine codes added (03's code lists):
#        fire (BR, CRB) -> fire = 1; clear-cut family (CT, CPR, ...) -> harv = 2 (heavy); ES -> pest = 1; at the END visit of the interval,
#        carried forward to the plot's later intervals (as 03 does within an epoch); TSD = 0 at that visit, then years since (decades), i.e.
#        the lower-bound rule 03 uses for undated codes (station_pe has no year). CHT / DT (03's truncation codes) -> the plot is truncated
#        before that visit (its intervals from that visit on are dropped).
#     "NEW" = the station code differs from the station code at the plot's previous visit, and pee_ori_sond has no origine code at that visit.
#     INCLUDE_PERTURB = TRUE also adds new station-only perturb codes (BRP -> fire, partial cuts -> harv = 1 if harv is 0, EL -> pest).
#   P1 entries and P2 exits are taken from the stored imputation (exit times stay inside the same brackets); the transition set and each
#   transition's terms are kept as in the int run.
#   Each LABEL is fitted with the pipeline's INLA call in both variants, plus an independent maximum-likelihood fit on the same rows.
#   WRONG = INLA error / NA, or |alpha_INLA - alpha_MLE| > 3 MLE-SE; screen_fail = INLA error / NA or alpha-SD ratio < 0.1 (the pipeline screen).
# MUST be run BEFORE 04_prepare_inla.R is re-run for the gs run (needs the int-version tables; checked below).
# OUTPUT: data/processed/fit_int/test_station_fits.csv + printed summary.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(INLA) }); options(width = 220)
RUN_TAG <- "int"; IN_DIR <- here("data", "processed", "fit_int")
N_CTRL <- 10L; N_FAIL_MAX <- 40L; N_WORKERS <- 6L; INCLUDE_PERTURB <- FALSE
LABELS <- c("black_spruce->regen", "fir->regen", "jack_pine->regen", "white_birch->regen", "trembling_aspen->regen", "fir->sugar_maple", "other_temperate->white_birch")
SIM_M <- 200L; M_MI <- 30L; SD_MIN <- 1e-3; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); COV_CONT <- c("gdd_3", "cmi_sum", "soil_pc1")      # the int run used annual CMI
source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
iv0 <- readRDS(here("data", "processed", "visit_intervals.rds"))
if (!"cmi_sum" %in% names(iv0)) stop("visit_intervals.rds is already the gs version: run this test before re-running 04 (or restore the int tables).")
inp <- prepare_fit_inputs(here("data", "processed"))
stopifnot(all(LABELS %in% inp$td$label))

# ---- station-only disturbance codes (03's lists) --------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(sf))
PEP_PATH <- here("data", "raw", "pep", "PEP.gpkg")
sta <- as.data.table(st_read(PEP_PATH, layer = "station_pe", quiet = TRUE))[, .(id_pe = as.character(id_pe), no_mes = as.integer(no_mes), so = origine, sp = perturb)]
pee <- as.data.table(st_read(PEP_PATH, layer = "pee_ori_sond", quiet = TRUE))[, .(id_pe = as.character(id_pe), no_mes = as.integer(no_mes), po = origine, pp = perturb)]
L <- merge(sta, pee, by = c("id_pe", "no_mes"), all.x = TRUE); setorder(L, id_pe, no_mes)
L[, `:=`(so_prev = shift(so), sp_prev = shift(sp)), by = id_pe]
fire_o <- c("BR", "CRB"); harv_o <- c("CT", "CPR", "CBA", "CBT", "CEF", "CRS", "CPH", "CPT", "CS", "CDV", "CPP"); pest_o <- "ES"; trunc_o <- c("CHT", "DT")
fire_p <- "BRP"; pest_p <- c("EL", "CE")
harv_p <- c("CP", "EPC", "CJ", "CDL", "EC", "CEA", "CAM", "CPI", "CPF", "CJG", "CJP", "CJT", "CPS", "CB", "CTR", "DLD", "CPM", "CPC", "CIP", "CA", "CPX", "ESI", "ECE", "CE", "CPR")
L[, new_o := !is.na(so) & so != "" & (is.na(po) | po == "") & (is.na(so_prev) | so != so_prev)]
L[, new_p := !is.na(sp) & sp != "" & (is.na(pp) | pp == "") & (is.na(sp_prev) | sp != sp_prev)]
add <- L[new_o == TRUE & so %in% c(fire_o, harv_o, pest_o, trunc_o), .(id_pe, visit_end = no_mes, code = so, level = "origine")]
if (INCLUDE_PERTURB) add <- rbind(add, L[new_p == TRUE & sp %in% c(fire_p, harv_p, pest_p), .(id_pe, visit_end = no_mes, code = sp, level = "perturb")])
cat("Station-only NEW codes added:\n"); print(add[, .N, by = .(level, code)][order(-N)])

iv_station <- function(iv) {                                   # returns a modified copy of the interval table (raw flags + TSD), then model columns
  x <- copy(iv); x[, `:=`(id_pe = as.character(id_pe), visit_end = as.integer(visit_end))]
  a <- merge(add, x[, .(id_pe, visit_end, t_end_code = t_end)], by = c("id_pe", "visit_end"))           # codes at an interval end of an included plot
  tr_ <- a[code %in% trunc_o, .(cut_visit = min(visit_end)), by = id_pe]
  x <- merge(x, tr_, by = "id_pe", all.x = TRUE); x <- x[is.na(cut_visit) | visit_end < cut_visit]; x[, cut_visit := NULL]
  a <- a[!code %in% trunc_o]
  for (j in seq_len(nrow(a))) {
    p <- a$id_pe[j]; v <- a$visit_end[j]; t0 <- a$t_end_code[j]; cd <- a$code[j]; lv <- a$level[j]
    w <- which(x$id_pe == p & x$visit_end >= v); if (!length(w)) next
    if (cd %in% c(fire_o, fire_p)) { x[w, tsd_fire := fifelse(fire == 1L, pmin(tsd_fire, (t_end - t0) / 10), (t_end - t0) / 10)]; x[w, fire := 1L] }
    if (cd %in% harv_o) { x[w, tsd_harv := fifelse(harv > 0L, pmin(tsd_harv, (t_end - t0) / 10), (t_end - t0) / 10)]; x[w, harv := 2L] }
    if (lv == "perturb" && cd %in% harv_p) { w2 <- w[x$harv[w] == 0L]; if (length(w2)) { x[w2, harv := 1L]; x[w2, tsd_harv := (t_end - t0) / 10] } }
    if (cd %in% c(pest_o, pest_p)) { x[w, tsd_pest := fifelse(pest == 1L, pmin(tsd_pest, (t_end - t0) / 10), (t_end - t0) / 10)]; x[w, pest := 1L] } }
  setorder(x, id_pe, interval)
  add_model_covariates(x, inp$scaling) }
iv_st <- iv_station(inp$iv)
cat(sprintf("Intervals: pipeline %d, station variant %d (truncation); flagged fire %d -> %d, heavy harvest %d -> %d, pest %d -> %d\n", nrow(inp$iv), nrow(iv_st),
            inp$iv[fire == 1, .N], iv_st[fire == 1, .N], inp$iv[harv == 2, .N], iv_st[harv == 2, .N], inp$iv[pest == 1, .N], iv_st[pest == 1, .N]))

# ---- rows of a stored imputation -----------------------------------------------------------------------------------------------------------------------
rebuild_rows <- function(m, ivx = inp$iv) {
  key_iv <- paste(ivx$id_pe, ivx$interval); b_abs <- ivx$t_end
  ex <- readRDS(file.path(IN_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, m))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
  ex_row <- rep(NA_real_, nrow(ivx)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); ex_row[ii[!is.na(ii)]] <- p2$exit_abs[!is.na(ii)]
  p2f <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
  fe <- setNames(p1$entry, as.character(p1$id_pe)); fe <- fe[names(fe) %in% as.character(unique(ivx$id_pe))]
  r <- build_rows(ivx, fe, inp$states, inp$tr, p2_fun = p2f, keep_cols = inp$keep_cols)$rows[kept == TRUE]
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
m0 <- fv[m <= M_MI & label %in% LABELS[1:2], .(all_ok = all(valid)), by = m][all_ok == TRUE, min(m)]
r0 <- rebuild_rows(m0)
chk <- rbindlist(lapply(LABELS[1:2], function(lb) one(r0[label == lb], lb, m0, "pipeline")))
chk <- merge(chk, mi[m == m0 & parameter == "alpha", .(label, est, sd)], by = "label"); chk[, dz := abs(alpha - est) / sd]
cat(sprintf("Step 0 (m = %d): |alpha - mi_raw| / sd = %s\n", m0, paste(round(chk$dz, 4), collapse = ", ")))
if (any(is.na(chk$dz) | chk$dz > 0.1)) stop("Rows / INLA call do not reproduce the int run; stop.")

# ---- imputations to test ------------------------------------------------------------------------------------------------------------------------------------
avail <- as.integer(sub(".*_m(\\d+)\\.rds$", "\\1", list.files(IN_DIR, pattern = sprintf("^extras_%s_m\\d+\\.rds$", RUN_TAG))))
fail_m <- head(intersect(fv[valid == FALSE, sort(unique(m))], avail), N_FAIL_MAX)
ctrl_m <- head(setdiff(intersect(fv[, .(all_ok = all(valid)), by = m][all_ok == TRUE, m], avail), fail_m), N_CTRL)
M_SET <- sort(c(fail_m, ctrl_m))
cat(sprintf("Imputations: %d with an int failure among LABELS (%s) + %d controls\n", length(fail_m), paste(fail_m, collapse = ","), length(ctrl_m)))

res <- rbindlist(mclapply(M_SET, function(m) tryCatch({
  r <- rebuild_rows(m); rx <- rebuild_rows(m, iv_st)
  rbindlist(lapply(LABELS, function(lb) rbind(one(r[label == lb], lb, m, "pipeline"), one(rx[label == lb], lb, m, "station")))) },
  error = function(e) data.table(m = m, label = NA_character_, variant = "run_error", inla_error = TRUE)), mc.cores = N_WORKERS), fill = TRUE)
res <- merge(res, fv[, .(m, label, int_valid = valid, int_reason = reason)], by = c("m", "label"), all.x = TRUE)
res[, `:=`(ratio = alpha_sd / (0.78 * alpha / sqrt(events)), z_mle = (alpha - mle_alpha) / mle_alpha_se)]
res[, `:=`(screen_fail = inla_error | is.na(alpha) | ratio < 0.1, wrong = inla_error | is.na(alpha) | abs(z_mle) > 3)]
fwrite(res, file.path(IN_DIR, "test_station_fits.csv"))

# ---- summary ------------------------------------------------------------------------------------------------------------------------------------------------
cat("\n== Per transition and variant: failures (screen) and WRONG fits (INLA alpha > 3 MLE-SE from the MLE) ==\n")
print(res[!is.na(label), .(fits = .N, screen_fail = sum(screen_fail), wrong = sum(wrong), events = round(mean(events)), trunc_max_dec = round(max(trunc_max), 1),
                           mle_alpha = round(median(mle_alpha), 3), inla_alpha_ok = round(median(alpha[!wrong]), 3), median_z_ok = round(median(z_mle[!wrong]), 2)),
          by = .(label, variant)][order(label, variant)])
cat("\n== Imputations that FAILED in the int run: repaired by the station codes? (screen_fail / wrong per variant) ==\n")
print(dcast(res[int_valid == FALSE, .(label, m, variant, bad = screen_fail | wrong)], label ~ variant, value.var = "bad", fun.aggregate = function(z) paste0(sum(z), "/", length(z))))
cat("\n== Healthy fits: median MLE alpha / intercept and events, pipeline vs station ==\n")
print(res[!is.na(label) & wrong == FALSE, .(mle_alpha = round(median(mle_alpha), 3), mle_icpt = round(median(mle_icpt), 3), events = round(mean(events))), by = .(label, variant)][order(label, variant)])
