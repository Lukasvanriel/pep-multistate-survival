# ==============================================================================
# R/tests/22_cap_test_new_rows.R   (2026-11-03)   RUN ON THE MACHINE WITH INLA      -- TEST ONLY, changes nothing in the pipeline
#
# QUESTION: once intervals are cut at the disturbance date and an X->regen exit is placed at its stand-replacing
#  disturbance (one-year flagged window before it), do MORE transitions hit INLA's truncation cap (H(tr) > 22.85)
#  than with today's rows? The answer decides between "caveat" and the fixed-alpha Poisson fit (decision list, 2026-11-03).
#
# DATA: ONE stored main-run imputation (CAP_M, default 1): its P1 entries and P2 exits (extras_main_mXX.rds), rows rebuilt
#  with the pipeline's build_rows() as in R/tests/08f. Both arms use the NEW model: weibullsurv variant 0, N(0, 5^2) on all
#  coefficients, flat intercept, pc.alphaw(ALPHA_PRIOR), one thread; harvest = light + heavy whenever light + heavy >= 15
#  events (never "any"); TSD per type and harvest level as log(TSD in years + 1). Term counts from transition_set (end-visit flags).
#   ARM A  today's rows: one row per interval piece, flags and TSD of the END visit, stored P2 exits.
#   ARM B  new rows, built in this order (cutting LAST):
#          1. P1 entries as stored.
#          2. Exits: X->regen event interval with a heavy harvest or fire whose date falls inside the interval -> exit AT the
#             disturbance date (recorded year, or the imputed date when undated; UNDATED_RULE). Several -> the earliest.
#             All other exits: stored P2 exit (the old moment-matched point; the piecewise P2 draw does not exist yet).
#          3. build_rows() (unchanged).
#          4. Cut every row at each disturbance switch time = disturbance date - WINDOW (1 y), not before the start visit.
#             Before the switch: start-visit flag/level of that type, TSD aged from the start visit. After: end-visit values,
#             TSD = end-visit TSD minus the time left to the end visit (0 inside the window). A type whose value changes
#             without its own date inside the interval (ended by another type's stand-origin event) switches at the earliest
#             other date; if there is none it takes the end values for the whole interval (counted as "unexplained").
#          TSD of a piece = value at the piece MIDPOINT; pieces are NOT further split at 1/3/10/30 y (decision still open).
#   ARM B with CAP_CUT=FALSE (added 2026-11-03): today's rows (no cutting; end-visit flags and TSD over the whole interval
#          piece, as arm A) with ONLY step 2 applied: X->regen exits placed at the heavy harvest / fire date. Output folder
#          suffix "_nocut". Question: does the exit rule alone keep the cap at today's level?
#  NOT in this test: pest (ES) as stand-replacing (pest is 0/1 in the tables); piecewise P2; climate NA fix (9 intervals:
#  INLA treats NA as 0, the MAP below does the same).
# PER FIT: INLA error, alpha + alpha-SD ratio screen (>= 0.1), max H(tr) over the rows at INLA's posterior means and the
#  number of rows above the cap, and a MAP reference (penalised maximum likelihood, same likelihood and N(0, 5^2) penalty, no
#  cap) with z = (INLA - MAP) / MAP SE per coefficient. Fits above the cap with large |z| = distorted by the cap.
# OUTPUT data-output/diagnostics/cap_test_new_rows_mXX/: build_stats.txt, fits/ (one rds per arm x transition, resumable),
#  results.csv (arm x transition), coef.csv (arm x transition x parameter), summary printed at the end.
# Run from the project root, nothing else running:
#   Rscript R/tests/22_cap_test_new_rows.R 2>&1 | tee data-output/diagnostics/cap_test_run.log
#   env: CAP_M (imputation, default 1), CAP_WORKERS (default 6), CAP_ENGINE ("inla" default; "map" = MAP only, plumbing check),
#        CAP_CUT (TRUE default; FALSE = exit rule only), CAP_XREGEN_RULE (TRUE), CAP_WINDOW_ALL (TRUE)
#   exit-rule-only run:  CAP_CUT=FALSE Rscript R/tests/22_cap_test_new_rows.R 2>&1 | tee data-output/diagnostics/cap_test_nocut.log
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here); library(parallel); library(tools) })
source(here::here("R", "config.R")); options(width = 220)
FIT_ENGINE   <- Sys.getenv("CAP_ENGINE", "inla")
REAL_M       <- as.integer(Sys.getenv("CAP_M", "1"))
N_WORKERS    <- as.integer(Sys.getenv("CAP_WORKERS", "6"))
SIGMA_B      <- 5; VARIANT <- 0L; MIN_EVENTS <- 15L
CAP          <- -log(0.001 * 2^-23)     # 22.85: INLA clamps F(truncation) at 1 - 0.001 * FLT_EPSILON
WINDOW       <- 1                       # years flagged before the disturbance date (user proposal 2026-11-03)
UNDATED_RULE <- TRUE                    # X->regen exit rule also for undated (imputed-date) heavy harvest / fire
XREGEN_RULE  <- Sys.getenv("CAP_XREGEN_RULE", "TRUE") == "TRUE"   # FALSE: stored P2 exits everywhere (cut only)
WINDOW_ALL   <- Sys.getenv("CAP_WINDOW_ALL", "TRUE") == "TRUE"    # FALSE: the window only for the stand-replacing type of a
                                                                   # rule-placed X->regen exit; every other cut at the date itself
CUT          <- Sys.getenv("CAP_CUT", "TRUE") == "TRUE"           # FALSE (added 2026-11-03): arm B = TODAY's rows (no cutting,
                                                                   # end-visit flags/TSD) with only the X->regen exit rule applied
SIM_M <- 200L; M_MI <- 30L; N_COV_SAMPLES <- 0L
source(here("R", "covariate_spec.R")); source(here("R", "04b_build_rows.R")); source(here("R", "05a_fit_funs.R"))
if (FIT_ENGINE == "inla") suppressPackageStartupMessages(library(INLA))
VTAG <- paste0(if (!XREGEN_RULE) "_norule" else "", if (!WINDOW_ALL) "_windowxregen" else "", if (!CUT) "_nocut" else "")
OUT <- file.path(DIAG_OUT_DIR, sprintf("cap_test_new_rows_m%02d%s%s", REAL_M, VTAG, if (FIT_ENGINE != "inla") paste0("_", FIT_ENGINE) else ""))
dir.create(file.path(OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
statf <- file.path(OUT, "build_stats.txt"); if (file.exists(statf)) file.remove(statf)
logl <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = statf, append = TRUE) }
logl("cap test | %s | engine %s | imputation m%02d | variant %d | sigma %g | pc.alphaw(%g) | cut %s | window %g y (all cuts: %s) | X->regen rule %s (undated %s) | TSD_SOURCE %s",
     format(Sys.time()), FIT_ENGINE, REAL_M, VARIANT, SIGMA_B, ALPHA_PRIOR, CUT, WINDOW, WINDOW_ALL, XREGEN_RULE, UNDATED_RULE, TSD_SOURCE)

# ---- inputs -----------------------------------------------------------------------------------------------------------------
inp <- prepare_fit_inputs(TABLE_DIR)
iv  <- copy(inp$iv)                         # same row order as inp$iv (build_rows indexes it)
iv[, dur := year_end - year_start]
TY <- c("fire", "harv", "pest"); TSDV <- c(fire = "tsd_fire", harv = "tsd_harv", pest = "tsd_pest")
for (T in TY) {                             # end-visit values; TSD back in YEARS (04 stores decades, 0 where unflagged)
  set(iv, j = paste0(T, "_e"), value = iv[[T]])
  set(iv, j = paste0("tsdy_", T, "_e"), value = fifelse(iv[[T]] > 0, 10 * iv[[TSDV[[T]]]], NA_real_))
}
for (T in TY) {                             # start-visit values = previous interval's end values ...
  iv[, (paste0(T, "_s")) := data.table::shift(get(paste0(T, "_e"))), by = id_pe]
  iv[, (paste0("tsdy_", T, "_s")) := data.table::shift(get(paste0("tsdy_", T, "_e"))), by = id_pe]
}
cc <- as.data.table(readRDS(file.path(BUILD_DIR, "combined_covariates.rds")))   # ... and for interval 1 the plot's first visit
cc[, no_mes := as.integer(no_mes)]; setorder(cc, id_pe, no_mes)
sfx <- c(lower = "_imputed", upper = "_imputed_prev", mid = "_imputed_mid")[[TSD_SOURCE]]
fv <- cc[, .SD[1L], by = id_pe][, .(id_pe, y1 = year, f1 = fire, h1 = harv, p1f = pest,
                                    tf1 = get(paste0("tsf", sfx)), th1 = get(paste0("tsh", sfx)), tp1 = get(paste0("sfp", sfx)))]
iv[fv, on = "id_pe", `:=`(y1 = i.y1, f1 = i.f1, h1 = i.h1, p1f = i.p1f, tf1 = i.tf1, th1 = i.th1, tp1 = i.tp1)]
stopifnot(iv[interval == 1L, all(y1 == year_start)])
iv[interval == 1L, `:=`(fire_s = f1, harv_s = h1, pest_s = p1f,
                        tsdy_fire_s = fifelse(f1 > 0, tf1, NA_real_), tsdy_harv_s = fifelse(h1 > 0, th1, NA_real_), tsdy_pest_s = fifelse(p1f > 0, tp1, NA_real_))]
stopifnot("start values differ from 04's *_start flags" = all(iv$fire_s == iv$fire_start) && all(iv$harv_s == iv$harv_start) && all(iv$pest_s == iv$pest_start))
for (T in TY) stopifnot(!anyNA(iv[[paste0("tsdy_", T, "_s")]][iv[[paste0(T, "_s")]] > 0]), !anyNA(iv[[paste0("tsdy_", T, "_e")]][iv[[paste0(T, "_e")]] > 0]))

# ---- disturbance dates inside intervals, switch times, X->regen exit rule ---------------------------------------------------
for (T in TY) {
  te <- iv[[paste0("tsdy_", T, "_e")]]
  on <- iv[[paste0(T, "_e")]] > 0 & !is.na(te) & te < iv$dur - 1e-9           # this type's (newest) event dates inside the interval
  set(iv, j = paste0("on_", T), value = on)
  set(iv, j = paste0("D_", T), value = fifelse(on, iv$year_end - te, NA_real_))
}
iv[, D_any := pmin(D_fire, D_harv, D_pest, na.rm = TRUE)]
iv[, sr_harv := XREGEN_RULE & is_change & state_end == "regen" & on_harv & harv_e == 2L & (UNDATED_RULE | harv_dating %in% "dated")]
iv[, sr_fire := XREGEN_RULE & is_change & state_end == "regen" & on_fire & fire_e == 1L & (UNDATED_RULE | fire_dating %in% "dated")]
iv[, D_sr := pmin(fifelse(sr_harv, D_harv, NA_real_), fifelse(sr_fire, D_fire, NA_real_), na.rm = TRUE)]
for (T in TY) {
  on <- iv[[paste0("on_", T)]]; dif <- iv[[paste0(T, "_s")]] != iv[[paste0(T, "_e")]]
  w_own <- if (WINDOW_ALL) WINDOW else fifelse(!is.na(iv$D_sr) & on & iv[[paste0("D_", T)]] == iv$D_sr, WINDOW, 0)
  w_oth <- if (WINDOW_ALL) WINDOW else 0
  sw <- fifelse(on, pmax(iv[[paste0("D_", T)]] - w_own, iv$year_start),
         fifelse(dif & !is.na(iv$D_any), pmax(iv$D_any - w_oth, iv$year_start), iv$year_start))
  set(iv, j = paste0("sw_", T), value = sw)
  set(iv, j = paste0("unexpl_", T), value = !on & dif & is.na(iv$D_any))
}
logl("intervals %d | with a disturbance date inside: fire %d, harvest %d, pest %d (any %d) | value change without own date: %s",
     nrow(iv), sum(iv$on_fire), sum(iv$on_harv), sum(iv$on_pest), sum(!is.na(iv$D_any)),
     paste(sprintf("%s %d (unexplained %d)", TY, vapply(TY, function(T) sum(iv[[paste0(T, "_s")]] != iv[[paste0(T, "_e")]] & !iv[[paste0("on_", T)]]), integer(1)),
                   vapply(TY, function(T) sum(iv[[paste0("unexpl_", T)]]), integer(1))), collapse = ", "))
logl("X->regen events %d | exit placed at a stand-replacing disturbance: %d (heavy harvest %d, fire %d; dated %d, imputed date %d)",
     iv[is_change & state_end == "regen", .N], sum(!is.na(iv$D_sr)), sum(iv$sr_harv), sum(iv$sr_fire),
     iv[!is.na(D_sr) & ((sr_harv & harv_dating %in% "dated") | (sr_fire & fire_dating %in% "dated")), .N],
     iv[!is.na(D_sr) & !((sr_harv & harv_dating %in% "dated") | (sr_fire & fire_dating %in% "dated")), .N])

# ---- rows ---------------------------------------------------------------------------------------------------------------------
ex <- readRDS(file.path(FIT_DIR, sprintf("extras_%s_m%02d.rds", RUN_TAG, REAL_M))); p1 <- as.data.table(ex$p1); p2 <- as.data.table(ex$p2); rm(ex)
key_iv <- paste(iv$id_pe, iv$interval); b_abs <- iv$t_end
ex_row <- rep(NA_real_, nrow(iv)); ii <- match(paste(p2$id_pe, p2$interval), key_iv); stopifnot(!anyNA(ii)); ex_row[ii] <- p2$exit_abs
first_entry <- setNames(p1$entry, as.character(p1$id_pe))
p2_A <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur; min(max(ex_row[r] - E, a_dur), b_dur) }
p2_B <- function(a_dur, b_dur, from, to, r) { E <- b_abs[r] - b_dur
  T_abs <- if (!is.na(iv$D_sr[r])) iv$t_start[r] + (iv$D_sr[r] - iv$year_start[r]) else ex_row[r]
  min(max(T_abs - E, a_dur), b_dur) }

make_rows <- function(arm) {
  segs <- build_rows(iv, first_entry, inp$states, inp$tr, p2_fun = if (arm == "A") p2_A else p2_B, keep_cols = character(0))$segs
  s <- segs[, .(seg_id, id_pe, iv_row, type, from, event_to, E_abs, start_abs, end_abs)]
  do_cut <- arm == "B" && CUT
  cuts <- if (do_cut) rbindlist(lapply(TY, function(T) data.table(iv_row = seq_len(nrow(iv)), ct = iv$t_start + (iv[[paste0("sw_", T)]] - iv$year_start)))) else NULL
  cc2 <- if (do_cut) merge(s[, .(seg_id, iv_row, start_abs, end_abs)], cuts, by = "iv_row", allow.cartesian = TRUE)[ct > start_abs + 1e-9 & ct < end_abs - 1e-9, .(seg_id, ct)] else NULL
  pts <- unique(rbind(s[, .(seg_id, ct = start_abs)], s[, .(seg_id, ct = end_abs)], cc2))[order(seg_id, ct)]
  # zero-length segments (exit exactly at a visit: build_rows floors their Tstop) stay one piece
  pc <- pts[, if (.N == 1L) .(p_start = ct, p_end = ct) else .(p_start = ct[-.N], p_end = ct[-1L]), by = seg_id][s, on = "seg_id"]
  stopifnot(!anyNA(pc$p_start), nrow(pc[, .N, by = seg_id]) == nrow(s))
  pc[, last := p_end >= end_abs - 1e-12]
  r <- pc$iv_row; cal <- iv$year_start[r] + ((pc$p_start + pc$p_end) / 2 - iv$t_start[r])       # calendar year at the piece midpoint
  for (T in TY) {
    fe <- iv[[paste0(T, "_e")]][r]; te <- iv[[paste0("tsdy_", T, "_e")]][r]
    if (do_cut) {
      fs <- iv[[paste0(T, "_s")]][r]; ts <- iv[[paste0("tsdy_", T, "_s")]][r]; pre <- cal < iv[[paste0("sw_", T)]][r]
      f <- ifelse(pre, fs, fe); tsd <- ifelse(pre, ts + (cal - iv$year_start[r]), pmax(0, te - (iv$year_end[r] - cal)))
    } else { f <- fe; tsd <- te }
    tsd[f == 0] <- 0; stopifnot(!anyNA(tsd), all(tsd >= 0))
    set(pc, j = paste0("f_", T), value = as.integer(f)); set(pc, j = paste0("ty_", T), value = tsd)
  }
  for (v in c("z_gdd_3", "z_cmi_gs", "z_soil_pc1")) set(pc, j = v, value = inp$iv[[v]][r])
  pc[, `:=`(fire_f = as.integer(f_fire == 1L), pest_f = as.integer(f_pest == 1L), harv_light = as.integer(f_harv == 1L), harv_heavy = as.integer(f_harv == 2L))]
  pc[, `:=`(gdd_x_fire = z_gdd_3 * fire_f, gdd_x_pest = z_gdd_3 * pest_f, gdd_x_harv_light = z_gdd_3 * harv_light, gdd_x_harv_heavy = z_gdd_3 * harv_heavy,
            tsd_fire = fire_f * log1p(ty_fire), tsd_pest = pest_f * log1p(ty_pest),
            tsd_harv_light = harv_light * log1p(ty_harv), tsd_harv_heavy = harv_heavy * log1p(ty_harv))]
  pc[, truncation := pmax((p_start - E_abs) / TIME_SCALE, 0)]
  pc[, Tstop := pmax((p_end - E_abs) / TIME_SCALE, truncation + EPSILON)]
  # checks: time at risk per plot = last - first visit; pieces contiguous
  tar <- pc[, .(tar = sum(p_end - p_start)), by = id_pe][iv[, .(span = t_end[.N] - t_start[1L]), by = id_pe], on = "id_pe"]
  stopifnot("time at risk per plot != last - first visit" = max(abs(tar$tar - tar$span)) < 1e-6)
  pairs <- CJ(from = inp$states, to = inp$states)[from != to]
  rows <- merge(pc, pairs, by = "from", allow.cartesian = TRUE)
  rows[, status := as.integer(type == "event" & last & to == event_to)]
  rows[, label := paste(from, to, sep = "->")]
  rows <- rows[label %in% inp$td$label]
  ev <- rows[status == 1L, .N, by = label][inp$td[, .(label, n_events)], on = "label"]
  stopifnot("status-1 rows != events of transition_set" = all(ev$N == ev$n_events))
  list(rows = rows, pieces = pc)
}
A <- make_rows("A"); B <- make_rows("B")
logl("pieces: A %d, B %d (+%d from cuts) | rows (kept transitions): A %d, B %d", nrow(A$pieces), nrow(B$pieces), nrow(B$pieces) - nrow(A$pieces), nrow(A$rows), nrow(B$rows))
mv <- merge(A$pieces[type == "event" & last, .(seg_id_A = seg_id, id_pe, iv_row, TA = p_end)], B$pieces[type == "event" & last, .(id_pe, iv_row, TB = p_end)], by = c("id_pe", "iv_row"))
mv <- mv[!is.na(iv$D_sr[iv_row])]
if (nrow(mv)) logl("X->regen exits moved to the disturbance date: %d | shift (years, B - A): median %.1f, 5-95%% %.1f to %.1f", nrow(mv), median(mv$TB - mv$TA), quantile(mv$TB - mv$TA, .05), quantile(mv$TB - mv$TA, .95))
fl <- function(R) R[status == 1L, .(fire = sum(fire_f), harv_light = sum(harv_light), harv_heavy = sum(harv_heavy), pest = sum(pest_f))]
logl("flagged events (status-1 rows), all kept transitions: transition_set (end visit) fire %d, light %d, heavy %d, pest %d | A %s | B %s",
     sum(inp$td$n_fire), sum(inp$td$n_harv1), sum(inp$td$n_harv2), sum(inp$td$n_pest),
     paste(unlist(fl(A$rows)), collapse = "/"), paste(unlist(fl(B$rows)), collapse = "/"))

# ---- terms (new harvest rule) -----------------------------------------------------------------------------------------------------
td <- copy(inp$td)[, harv_new := n_harv1 + n_harv2 >= MIN_EVENTS]
terms_of <- function(i) c("z_gdd_3", "z_cmi_gs", "z_soil_pc1",
  if (td$fire_term[i]) c("fire_f", "gdd_x_fire", "tsd_fire"), if (td$pest_term[i]) c("pest_f", "gdd_x_pest", "tsd_pest"),
  if (td$harv_new[i]) c("harv_light", "harv_heavy", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_harv_light", "tsd_harv_heavy"))
logl("terms: %d transitions | fire %d | pest %d | harvest light + heavy %d", nrow(td), sum(td$fire_term), sum(td$pest_term), sum(td$harv_new))

# ---- fitting ------------------------------------------------------------------------------------------------------------------
fit_inla_v0 <- function(d, terms) {
  f <- as.formula(sprintf("inla.surv(time = Tstop, event = status, truncation = truncation) ~ %s", paste(c("1", terms), collapse = " + ")))
  m <- tryCatch(INLA::inla(f, family = "weibullsurv", data = d,
                           control.family = list(variant = VARIANT, hyper = list(alpha = list(prior = "pc.alphaw", param = ALPHA_PRIOR))),
                           control.fixed = list(prec = 1 / SIGMA_B^2, prec.intercept = 0),
                           control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)),
                           num.threads = "1:1"), error = function(e) list(error = conditionMessage(e)))
  pars <- c("alpha", "(Intercept)", terms)
  if (!is.null(m$error)) return(list(tab = data.table(parameter = pars, est = NA_real_, sd = NA_real_), error = substr(m$error, 1, 200)))
  fx <- m$summary.fixed; hy <- m$summary.hyperpar
  g <- function(p, col) if (p == "alpha") hy[ALPHA_NAME, col] else if (p %in% rownames(fx)) fx[p, col] else NA_real_
  list(tab = data.table(parameter = pars, est = vapply(pars, g, 1, col = "mean"), sd = vapply(pars, g, 1, col = "sd")), error = NA_character_)
}
fit_map_v0 <- function(d, terms) {      # penalised ML, variant 0: h = alpha exp(eta) t^(alpha-1), H = exp(eta) t^alpha; no truncation cap
  X <- cbind(1, as.matrix(d[, ..terms])); X[is.na(X)] <- 0; Tt <- d$Tstop; tr <- d$truncation; s <- d$status; lT <- log(Tt)
  nll <- function(th) { a <- exp(th[1]); eta <- drop(X %*% th[-1]); -(sum(s * (th[1] + eta + (a - 1) * lT)) - sum(exp(eta) * (Tt^a - tr^a))) + sum(th[-(1:2)]^2) / (2 * SIGMA_B^2) }
  o <- optim(c(0, log(max(sum(s), 0.5) / sum(Tt - tr)), rep(0, length(terms))), nll, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  V <- tryCatch(solve(optimHess(o$par, nll)), error = function(e) NULL); se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  a <- exp(o$par[1])
  list(tab = data.table(parameter = c("alpha", "(Intercept)", terms), est = c(a, o$par[-1]), sd = c(a * se[1], se[-1])), conv = o$convergence)
}
H_at <- function(d, terms, tab) {       # cumulative hazard at the truncation time of every row, at the given estimates
  est <- setNames(tab$est, tab$parameter); if (anyNA(est)) return(c(maxH = NA_real_, n_above = NA_real_, n_above_ev = NA_real_))
  X <- as.matrix(d[, ..terms]); X[is.na(X)] <- 0
  H <- exp(est[["(Intercept)"]] + drop(X %*% est[terms])) * d$truncation^est[["alpha"]]
  c(maxH = max(H), n_above = sum(H > CAP), n_above_ev = sum(H > CAP & d$status == 1L))
}
cols <- c("status", "Tstop", "truncation", "z_gdd_3", "z_cmi_gs", "z_soil_pc1", "fire_f", "pest_f", "harv_light", "harv_heavy",
          "gdd_x_fire", "gdd_x_pest", "gdd_x_harv_light", "gdd_x_harv_heavy", "tsd_fire", "tsd_pest", "tsd_harv_light", "tsd_harv_heavy")
DATA <- list(A = split(A$rows[, c("label", cols), with = FALSE], by = "label", keep.by = FALSE),
             B = split(B$rows[, c("label", cols), with = FALSE], by = "label", keep.by = FALSE))
rm(A, B); invisible(gc())
jobs <- CJ(arm = c("A", "B"), i = seq_len(nrow(td)))
jobs[, file := file.path(OUT, "fits", sprintf("%s_%02d.rds", arm, i))]
todo <- jobs[!file.exists(file)]
cat(sprintf("Fits to run: %d of %d (%d workers)\n", nrow(todo), nrow(jobs), N_WORKERS))
t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(todo)), function(j) {
  arm <- todo$arm[j]; i <- todo$i[j]; lb <- td$label[i]; tms <- terms_of(i); d <- DATA[[arm]][[lb]]
  t1 <- proc.time()[["elapsed"]]
  fi <- if (FIT_ENGINE == "inla") fit_inla_v0(d, tms) else list(tab = data.table(parameter = c("alpha", "(Intercept)", tms), est = NA_real_, sd = NA_real_), error = "not run")
  sec <- proc.time()[["elapsed"]] - t1
  fm <- tryCatch(fit_map_v0(d, tms), error = function(e) list(tab = data.table(parameter = c("alpha", "(Intercept)", tms), est = NA_real_, sd = NA_real_), conv = -1L))
  saveRDS(list(arm = arm, label = lb, events = sum(d$status), rows = nrow(d), error = fi$error, sec = sec, map_conv = fm$conv,
               inla = fi$tab, map = fm$tab, H_inla = H_at(d, tms, fi$tab), H_map = H_at(d, tms, fm$tab)), todo$file[j])
  cat(sprintf("[%s] %s %-35s %s %.0f s\n", format(Sys.time(), "%H:%M:%S"), arm, lb, if (is.na(fi$error)) "ok" else "ERROR", sec))
  NULL }, mc.cores = N_WORKERS, mc.preschedule = FALSE))
cat(sprintf("Fitting done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- results ------------------------------------------------------------------------------------------------------------------
R <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
res <- rbindlist(lapply(R, function(x) { a <- x$inla[parameter == "alpha"]
  data.table(arm = x$arm, label = x$label, events = x$events, rows = x$rows, error = x$error, sec = x$sec, map_conv = x$map_conv,
             alpha = a$est, alpha_sd = a$sd, ratio = a$sd / (0.78 * a$est / sqrt(x$events)), alpha_map = x$map[parameter == "alpha", est],
             maxH_inla = x$H_inla[["maxH"]], n_above_inla = x$H_inla[["n_above"]], n_above_ev_inla = x$H_inla[["n_above_ev"]],
             maxH_map = x$H_map[["maxH"]], n_above_map = x$H_map[["n_above"]]) }))
res[, `:=`(screen_ok = is.na(error) & !is.na(ratio) & ratio >= 0.1, above_cap = !is.na(maxH_inla) & maxH_inla > CAP, above_cap_map = !is.na(maxH_map) & maxH_map > CAP)]
coef <- rbindlist(lapply(R, function(x) merge(x$inla[, .(parameter, inla = est, inla_sd = sd)], x$map[, .(parameter, map = est, map_se = sd)], by = "parameter")[, `:=`(arm = x$arm, label = x$label)]))
coef[, z := (inla - map) / map_se]
fwrite(res, file.path(OUT, "results.csv")); fwrite(coef, file.path(OUT, "coef.csv"))

cat("\n==== SUMMARY (cap = 22.85; H at INLA posterior means unless stated) ====\n")
print(res[, .(transitions = .N, inla_errors = sum(!is.na(error) & error != "not run"), screen_fail = sum(is.na(error) & !screen_ok, na.rm = TRUE),
              above_cap = sum(above_cap), above_cap_at_MAP = sum(above_cap_map)), by = arm])
w <- dcast(res, label ~ arm, value.var = c("above_cap_map", "maxH_map", "maxH_inla", "screen_ok"))
cat("\nTransitions above the cap at the MAP (cap-free reference) in either arm:\n")
print(w[above_cap_map_A | above_cap_map_B, .(label, A = above_cap_map_A, B = above_cap_map_B, maxH_map_A = round(maxH_map_A, 1), maxH_map_B = round(maxH_map_B, 1),
                                            maxH_inla_A = round(maxH_inla_A, 1), maxH_inla_B = round(maxH_inla_B, 1), screen_A = screen_ok_A, screen_B = screen_ok_B)][order(-maxH_map_B)])
cat("\nINLA vs MAP (|z| = |INLA - MAP| / MAP SE), coefficients other than alpha / intercept, by arm and cap status at the MAP:\n")
cz <- merge(coef[!parameter %in% c("alpha", "(Intercept)")], res[, .(arm, label, above_cap_map)], by = c("arm", "label"))
print(cz[, .(fits = uniqueN(label), coefs = .N, median_abs_z = round(median(abs(z), na.rm = TRUE), 2), share_gt1 = round(mean(abs(z) > 1, na.rm = TRUE), 2)), by = .(arm, above_cap_map)][order(arm, above_cap_map)])
cat("\nX->regen: heavy-harvest and fire main effects (log HR), MAP and INLA, arm A vs B:\n")
print(dcast(coef[grepl("->regen$", label) & parameter %in% c("harv_heavy", "fire_f", "alpha")], label + parameter ~ arm, value.var = c("map", "inla"))[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
cat(sprintf("\nOutputs in %s\n", OUT))
