# ==============================================================================
# pep_disturbance_covariates.R
#
# Builds disturbance covariates for PEP permanent plots, starting from the
# full raw PEP geopackage (pee_ori_sond + placette_mes) -- not gated by the
# plot set already filtered in 02_classify_states.R. See discussion: this
# script applies its own complete ignore / drop-all-plot / drop-part-of-plot
# logic, and gets reconciled with state_classification.rds via an inner join
# downstream, with diagnostics on what each side drops that the other keeps.
#
# Planned steps:
#   1. Load pee_ori_sond + placette_mes, merge, compute year, coerce an_* [DONE]
#   2. Diagnostic gate: every non-NA origine/perturb code falls into exactly
#      one list (except CE, deliberately in both harv_p and pest_p)   [DONE]
#   3. Ignore: null origine/perturb (+ an_*) at that row for ignore_* codes [DONE]
#   4. Drop all plot: remove every row for any id_pe ever carrying a
#      dropplot_* code                                                [DONE]
#   5. Drop part of plot: CHT/DT truncation rule                      [DONE]
#   6-8. Disturbance events (rewritten 2026-10-26, R/03a_disturbance_events.R): BA event per visit,
#      records, events, dating (real year or lower/upper bound + midpoint), timeline (stand-origin
#      events end older events; compound events in one interval kept), per-visit flags and TSD
#   9. Assemble output: final covariate shapes (fire/pest binary, harv
#      0/1/2), real + imputed tsf/tsh/sfp, id_pe/no_mes/year/total_ba  [DONE]
#  10. Diagnostics and write final dataset to disk                    [DONE]
#
# Output: data/build/disturbance_covariates.rds
# ==============================================================================

library(sf)
library(data.table)
library(here)

source(here::here("R", "config.R"))   # project paths (data/build etc.), BA_REGEN, BA_LOSS_THRESH
source(here::here("R", "03a_disturbance_events.R"))   # disturbance event logic (steps 6-8)
PEP_PATH <- here("data", "raw", "pep", "PEP.gpkg")

# ==============================================================================
# STEP 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
placette <- as.data.table(st_read(PEP_PATH, layer = "placette_mes", quiet = TRUE))
pee      <- as.data.table(st_read(PEP_PATH, layer = "pee_ori_sond", quiet = TRUE))

# year = calendar year of the plot visit, derived from date_sond (the only
# source for this -- confirmed no separate year/annee field exists anywhere
# in the geopackage; Brice et al.'s own pipeline derives it the same way,
# from the same field, calling it year_measured).
placette[, year := as.integer(format(date_sond, "%Y"))]

dt <- merge(
  pee[, .(id_pe, no_mes, origine, an_origine, perturb, an_perturb)],
  placette[, .(id_pe, no_mes, year)],
  by    = c("id_pe", "no_mes"),
  all.y = TRUE
)

# an_origine/an_perturb come in as text (confirmed: clean 4-digit strings,
# range 1890-2023 / 1933-2022, no sentinel garbage) -- coerce to integer so
# they can be compared numerically against year.
dt[, an_origine := as.integer(an_origine)]
dt[, an_perturb := as.integer(an_perturb)]

setorder(dt, id_pe, no_mes)

# ==============================================================================
# STEP 2. DISTURBANCE CODE LISTS
# ==============================================================================

# codes that are simply not counted at this visit -- plot and observation
# both stay in the dataset, fire/pest/harv just stay 0 for this code
ignore_origine <- c("ETR", "REA", "REC")
ignore_perturb <- c("CHP", "VEP", "DP")

# codes that mean the whole plot gets dropped from the dataset entirely
# (non-natural origin/intervention -- plantation family, reclamation, wasteland)
dropplot_origine <- c("FR", "ENS", "P", "PL", "PLR", "PLN", "PRR", "RPS", "ENM", "EA")
dropplot_perturb <- c("DEG", "NET", "ENR", "DRM", "CD", "RRG", "PL_", "RR")

# EA (water), ENM (seeding w/ mini-greenhouses) added to dropplot_origine, and
# CPR (as a perturb value) added to harv_p -- all three were found missing
# entirely from every list by the diagnostic gate below. EA: both plots
# carrying it also show statut_mes="AB" (Placette abandonnee) at that same
# visit, confirming the plot exits the survey, not a regen state. ENM: same
# family as ENS, already dropplot. CPR (perturb): same code/meaning as the
# CPR already in harv_o for origine, just recorded as a partial event.

# codes that trigger the truncation rule: keep only the plot's observations
# strictly before the first occurrence; drop the plot if <2 remain, or if
# the retained pre-event history itself has an unexplained crash (handled
# separately, not by simple list membership -- see Step 5)
droppart_origine <- c("CHT", "DT")

fire_o <- c("BR", "CRB")
fire_p <- c("BRP")

harv_o <- c("CT", "CPR", "CBA", "CBT", "CEF", "CRS", "CPH",
            "CPT", "CS", "CDV", "CPP")
harv_p <- c("CP", "EPC", "CJ", "CDL", "EC", "CEA", "CAM",
            "CPI", "CPF", "CJG", "CJP", "CJT", "CPS", "CB",
            "CTR", "DLD", "CPM", "CPC", "CIP", "CA", "CPX",
            "ESI", "ECE", "CE", "CPR")

pest_o <- c("ES")
pest_p <- c("EL", "CE")

# Shared basal-area thresholds used by BOTH Step 5's crash detection (an
# unexplained pre-event transition) and Step 6b's ambiguous same-year event
# resolution (ba_before()). Previously Step 5 independently hardcoded 0.5
# (50% loss) with no documented source, while Step 6b used 0.25 -- the one
# documented disturbance-severity threshold found in the PEP dictionary/MRNF
# norm. Unified to 0.25 everywhere (approved 2026-09-05): re-running Step 5's
# truncation with 0.25 instead of 0.5 moves 13 plots from "survives
# truncation" to "dropped entirely" (39 -> 26 of the 48 plots with >=2
# pre-event visits), because those 13 have a pre-event BA transition between
# 25% and 50% loss with no origine/perturb code attached.
# BA_LOSS_THRESH (0.25) and BA_REGEN (5) are set in R/config.R since 2026-10-26 (shared with 02).
stopifnot(exists("BA_LOSS_THRESH"), exists("BA_REGEN"))

# ==============================================================================
# STEP 2b. DIAGNOSTIC GATE
# Every non-NA origine/perturb code appearing in the raw data must fall into
# exactly one of the lists above (CE is the one deliberate exception -- it's
# allowed in both harv_p and pest_p). This runs on dt before any nulling/
# filtering below, on the full raw code vocabulary. Caught real gaps when
# first run: EA, ENM (origine) and CPR (perturb) were entirely unclassified
# -- now fixed above. Stops the script rather than silently letting an
# unclassified code fall through fcase()'s default in Step 6.
# ==============================================================================

cat("Checking code-list coverage...\n")

origine_lists <- list(ignore_origine, dropplot_origine, droppart_origine,
                       fire_o, harv_o, pest_o)
perturb_lists <- list(ignore_perturb, dropplot_perturb, fire_p, harv_p, pest_p)

check_coverage <- function(codes, lists, dual_ok = character(0)) {
  codes <- unique(codes[!is.na(codes)])
  counts <- sapply(codes, function(cd) sum(sapply(lists, function(l) cd %in% l)))
  expected <- ifelse(codes %in% dual_ok, 2L, 1L)
  bad <- codes[counts != expected]
  if (length(bad) > 0) {
    stop(sprintf("Uncovered/miscovered codes: %s", paste(bad, collapse = ", ")))
  }
}

check_coverage(dt$origine, origine_lists)
check_coverage(dt$perturb, perturb_lists, dual_ok = "CE")

cat("  origine and perturb code lists are exhaustive and mutually exclusive")
cat(" (except CE, deliberately dual-coded)\n")

# ==============================================================================
# STEP 3. IGNORE
# codes that don't count toward fire/pest/harv at this visit -- the plot and
# the observation both stay in the dataset. origine/perturb are independent
# columns, so a row can have one nulled while the other survives (confirmed
# real cases exist in the data: ETR/REA + a real harvest perturb code, and
# BR (fire) + a CHP/DP ignore perturb code) -- tested against these exact
# cases in base R before writing this.
# ==============================================================================

cat("Applying ignore codes...\n")

dt[origine %in% ignore_origine, `:=`(origine = NA_character_, an_origine = NA_integer_)]
dt[perturb %in% ignore_perturb, `:=`(perturb = NA_character_, an_perturb = NA_integer_)]

# ==============================================================================
# STEP 4. DROP ALL PLOT
# any plot ever carrying a non-natural-succession code (plantation family,
# reclamation, wasteland...) is removed entirely -- every visit, not just the
# flagged one. Confirmed via a real sample plot (0119600101): two early
# visits carry only CHP (ignore, not dropplot) and would look clean on their
# own, but visit 3 carries DEG (dropplot_perturb), so the whole plot -- all
# 3 visits -- must go. Tested this exact case in base R before writing it.
# Real-data scale: 922 distinct plots / 3,735 of 53,632 total rows.
# ==============================================================================

cat("Applying drop-all-plot filter...\n")

drop_ids <- unique(dt$id_pe[dt$origine %in% dropplot_origine | dt$perturb %in% dropplot_perturb])
cat(sprintf("  %d plots dropped entirely (%d rows)\n",
            length(drop_ids), sum(dt$id_pe %in% drop_ids)))

dt <- dt[!id_pe %in% drop_ids]

# ==============================================================================
# STEP 5. DROP PART OF PLOT (CHT/DT TRUNCATION)
#
# 5a. total_ba is recomputed here directly from dendro_arbres, independent of
#     02_classify_states.R -- confirmed necessary: 8 of the 113 plots that
#     carry CHT/DT are entirely absent from state_classification.rds for
#     reasons unrelated to disturbance (script 2's own exclusions), so
#     joining from that file would silently lose the BA history needed for
#     their truncation check. Verified those exact 8 plots DO have full BA
#     history in dendro_arbres directly.
#
# 5b. For each plot with any origine %in% droppart_origine: truncate to
#     observations strictly before the first occurrence; drop the plot if
#     fewer than 2 remain; among the retained pre-event rows, drop the plot
#     if there's an unexplained (BA loss > BA_LOSS_THRESH, or a drop below
#     the BA_REGEN regen threshold) transition with no origine/perturb code
#     attached. Tested this exact logic in Python against the real
#     geopackage before writing it here: 113 plots carry CHT/DT post steps
#     3-4 (112 CHT, 1 DT) -> 65 dropped for <2 remaining, 22 dropped for an
#     unexplained crash/regen transition, 26 survive. (Re-tested 2026-09-05
#     after unifying the threshold from an undocumented 0.5 to the
#     documented 0.25 shared with Step 6b -- previously 9 dropped for
#     crash/regen and 39 survived.)
# ==============================================================================

cat("Computing live basal area (independent of 02_classify_states.R)...
")

LIVE_ETAT <- c("10", "12", "30", "32", "40", "42", "50", "52")
dendro <- as.data.table(st_read(PEP_PATH, layer = "dendro_arbres", quiet = TRUE))

live_ba <- dendro[etat %in% LIVE_ETAT, .(total_ba = sum(st_ha, na.rm = TRUE)),
                   by = .(id_pe, no_mes)]
dt <- merge(dt, live_ba, by = c("id_pe", "no_mes"), all.x = TRUE)
dt[is.na(total_ba), total_ba := 0]
setorder(dt, id_pe, no_mes)
dt[, ba_event := ba_event_flag(total_ba, BA_REGEN, BA_LOSS_THRESH), by = id_pe]   # E0, shared BA flag (2026-10-26)

cat("Applying drop-part-of-plot (CHT/DT truncation)...
")

droppart_ids <- unique(dt[origine %in% droppart_origine, id_pe])

truncate_one_plot <- function(g) {
  setorder(g, no_mes)
  first_event <- min(g[origine %in% droppart_origine, no_mes])
  pre <- g[no_mes < first_event]
  if (nrow(pre) < 2) return(NULL)
  for (i in seq_len(nrow(pre) - 1)) {
    code_present <- !is.na(pre$origine[i + 1]) || !is.na(pre$perturb[i + 1])
    # ba_event = previous BA >= BA_REGEN and (loss >= BA_LOSS_THRESH or drop below BA_REGEN); since 2026-10-26 evening no BA
    # event when the previous BA was already < BA_REGEN (before: any crash from BA > 0 counted)
    if (pre$ba_event[i + 1] && !code_present) return(NULL)
  }
  pre
}

droppart_result <- rbindlist(
  lapply(droppart_ids, function(pid) truncate_one_plot(dt[id_pe == pid])),
  fill = TRUE
)
cat(sprintf("  %d plots carried CHT/DT; %d survive truncation (%d rows kept)
",
            length(droppart_ids), uniqueN(droppart_result$id_pe), nrow(droppart_result)))

dt <- rbindlist(list(dt[!id_pe %in% droppart_ids], droppart_result), fill = TRUE)
setorder(dt, id_pe, no_mes)

# ==============================================================================
# STEPS 6-8. DISTURBANCE EVENTS (rewritten 2026-10-26; functions in R/03a_disturbance_events.R)
# Replaces the former steps 6 (mark), 6b (salvage pass: forward pass + back-fill + consolidation + reset gate),
# 7 (epochs), 7b (per-row imputed years) and 8/8b (TSD). Old version: R/archive/_archive_pre_events/03_build_covariates.R.
# Fixed order: E0 BA event per visit -> E1 records -> E2 events -> E3 dating -> E4 timeline -> E5 visit values.
# Rules: header of R/03a_disturbance_events.R; tests: R/tests/03t_disturbance_events_test.R (invented plots) and
# R/tests/03u_compare_events_real_plots.R (real plots vs the previous build); decisions: docs/PIPELINE_NOTES.md.
# Output per visit: fire/harv/pest flags; TSD real (dated events), *_imputed (lower bound = min), *_imputed_prev
# (upper bound = max), *_imputed_mid (midpoint); *_dating = "dated" / "bounded" / "first_visit" (stand-origin event
# undated at the plot's first visit: no TSD here, 04 applies the stand-age rule).
# ==============================================================================

cat("Building disturbance events (steps 0-7 of R/03a_disturbance_events.R)...\n")
CODES <- list(fire_o = fire_o, fire_p = fire_p, harv_o = harv_o, harv_p = harv_p, pest_o = pest_o, pest_p = pest_p)
setorder(dt, id_pe, no_mes)
# stand age at each plot's first visit in THIS table (02b gives W_j at the plot's first visit in 02; shifted if 03's first visit
# is a different year -- 3 analysed plots in the 2026-10-26 check)
sa <- as.data.table(readRDS(file.path(BUILD_DIR, "stand_age.rds")))[, .(id_pe, W_j, year_1)]
f03 <- dt[, .(year_first03 = year[1L]), by = id_pe][sa, on = "id_pe", `:=`(W_j = i.W_j, year_1 = i.year_1)]
f03[, W1 := W_j + (year_first03 - year_1)]
cat(sprintf("  stand age at the first visit: %d plots with W_j, %d without; first visit shifted vs 02b: %d plots; negative W1: %d\n",
            f03[!is.na(W_j), .N], f03[is.na(W_j), .N], f03[!is.na(W_j) & year_first03 != year_1, .N], f03[W1 < 0, .N]))
if (f03[W1 < 0, .N] > 0) stop("negative stand age at the first visit -- inspect before continuing")
dt[f03, on = "id_pe", W1 := i.W1]
FV_AGE <- first_visit_partial_ages(dt, CODES, FV_AGE_QUANTILE)
cat(sprintf("  first-visit partial ages (quantile %.2f of dated first-visit partial codes): harvest %.1f y (n %d), pest %.1f y (n %d)\n",
            FV_AGE_QUANTILE, FV_AGE[["harv"]], attr(FV_AGE, "n")[["harv"]], FV_AGE[["pest"]], attr(FV_AGE, "n")[["pest"]]))
if (anyNA(FV_AGE)) stop("first-visit partial age could not be computed for a type")
n_unsorted <- dt[, .(u = is.unsorted(year)), by = id_pe][u == TRUE, .N]
if (n_unsorted > 0) stop(sprintf("%d plots have visit years out of order by no_mes", n_unsorted))

EV_LIST <- list(); REC_LIST <- list()
vis <- dt[, {
  r <- disturbance_plot(year, total_ba, origine, an_origine, perturb, an_perturb, CODES, BA_REGEN, BA_LOSS_THRESH,
                        W1 = W1[1L], fv_age = FV_AGE)
  ev <- r$events
  if (length(ev$type)) {
    EV_LIST[[.GRP]] <<- data.table(id_pe = .BY$id_pe, event = seq_along(ev$type), type = ev$type, level = ev$level, yr = ev$yr,
                                   first_no_mes = no_mes[ev$first], first_year = year[ev$first], first_is_v1 = ev$first == 1L,
                                   start_no_mes = no_mes[ifelse(is.finite(ev$start), as.integer(ev$start), NA_integer_)], dropped = ev$dropped)
    rr <- r$records
    REC_LIST[[.GRP]] <<- data.table(id_pe = .BY$id_pe, no_mes = no_mes[rr$v], type = rr$type, level = rr$level, yr = rr$yr,
                                    yr_src = rr$yr_src, n_raw = rr$n_raw, event = rr$event, action = rr$action)
  }
  r$visits
}, by = id_pe]
stopifnot("event output out of row order" = identical(vis$id_pe, dt$id_pe))
stopifnot("ba_event differs between step 5 and E0" = identical(vis$ba_event, dt$ba_event))
vis_cols <- setdiff(names(vis), c("id_pe", "ba_event"))
dt[, (vis_cols) := vis[, vis_cols, with = FALSE]]
events  <- rbindlist(EV_LIST);  records <- rbindlist(REC_LIST)

cat(sprintf("  %d records (after A1) -> %d events in %d plots; records merged by A1: %d; years from back-fill: %d\n",
            nrow(records), nrow(events), uniqueN(events$id_pe), records[n_raw > 1, .N], records[yr_src == "backfill", .N]))
cat("  Records by action (rule log):\n")
print(dcast(records[, .N, by = .(action, level = fifelse(level == 2L, "stand-origin", "partial"))], action ~ level, value.var = "N", fill = 0))
cat("  Events by type, level and dating:\n")
events[, dating := fcase(!is.na(yr), "dated", first_is_v1, "first_visit", default = "bounded")]
print(dcast(events[, .N, by = .(type, level = fifelse(level == 2L, "stand-origin", "partial"), dating)], type + level ~ dating, value.var = "N", fill = 0))
cat(sprintf("  Stand-origin events dropped by the selection rule: %d (in %d plots)\n", events[dropped == TRUE, .N], events[dropped == TRUE, uniqueN(id_pe)]))

# ---- sanity checks ---------------------------------------------------------------------------------------------
cat("Checking TSD consistency...\n")
for (tp in c("tsf", "tsh", "sfp")) {
  mn <- dt[[paste0(tp, "_imputed")]]; mx <- dt[[paste0(tp, "_imputed_prev")]]; md <- dt[[paste0(tp, "_imputed_mid")]]; re <- dt[[tp]]
  bad <- sum(mn < 0, na.rm = TRUE) + sum(mx < mn, na.rm = TRUE) + sum(md < mn | md > mx, na.rm = TRUE) +
         sum(!is.na(re) & (re != mn | re != mx))
  cat(sprintf("  %s: violations (negative, min > max, mid outside, real != bounds): %d\n", tp, bad))
  if (bad > 0) stop("TSD consistency check failed for ", tp)
}
for (ty in c("fire", "harv", "pest")) {
  fl <- dt[[ty]]; dg <- dt[[paste0(ty, "_dating")]]; ts <- c(fire = "tsf", harv = "tsh", pest = "sfp")[[ty]]
  if (any((fl > 0) != !is.na(dg))) stop("flag and dating marker disagree for ", ty)
  n_na <- dt[get(ty) > 0 & is.na(get(paste0(ts, "_imputed_mid"))), .N]
  cat(sprintf("  %s: flagged visits without a TSD: %d (plots without a stand age only)\n", ty, n_na))
}

# ==============================================================================
# STEP 9. ASSEMBLE OUTPUT
# fire / pest binary, harv 0/1/2 (from E5); TSD in years: real (dated events), *_imputed (lower bound),
# *_imputed_prev (upper bound), *_imputed_mid (midpoint, main analysis since 2026-10-24); *_dating marker for 04;
# ba_event (shared BA flag, for reporting).
# ==============================================================================

cat("Assembling output...
")

out <- dt[, .(
  id_pe, no_mes, year, total_ba, ba_event,
  fire, harv, pest,
  tsf, tsh, sfp,
  tsf_imputed, tsh_imputed, sfp_imputed,
  tsf_imputed_prev, tsh_imputed_prev, sfp_imputed_prev,
  tsf_imputed_mid, tsh_imputed_mid, sfp_imputed_mid,
  fire_dating, harv_dating, pest_dating
)]

# ==============================================================================
# STEP 10. DIAGNOSTICS AND SAVE
# ==============================================================================

cat("Final diagnostics...
")
cat(sprintf("  %d rows, %d distinct plots in final output
",
            nrow(out), uniqueN(out$id_pe)))
cat(sprintf("  fire active: %d rows (%.1f%%)
",
            sum(out$fire == 1), 100 * mean(out$fire == 1)))
cat(sprintf("  harv active: %d rows (%.1f%%) -- of which harv==2 (stand-replacing): %d, harv==1 (partial): %d
",
            sum(out$harv > 0), 100 * mean(out$harv > 0),
            sum(out$harv == 2), sum(out$harv == 1)))
cat(sprintf("  pest active: %d rows (%.1f%%)
",
            sum(out$pest == 1), 100 * mean(out$pest == 1)))

cat(sprintf("  tsf NA rate (of fire==1): %.1f%%  (imputed: %.1f%%, imputed_prev: %.1f%%)
",
            100 * mean(is.na(out[fire == 1, tsf])),
            100 * mean(is.na(out[fire == 1, tsf_imputed])),
            100 * mean(is.na(out[fire == 1, tsf_imputed_prev]))))
cat(sprintf("  tsh NA rate (of harv>0):  %.1f%%  (imputed: %.1f%%, imputed_prev: %.1f%%)
",
            100 * mean(is.na(out[harv > 0, tsh])),
            100 * mean(is.na(out[harv > 0, tsh_imputed])),
            100 * mean(is.na(out[harv > 0, tsh_imputed_prev]))))
cat(sprintf("  sfp NA rate (of pest==1): %.1f%%  (imputed: %.1f%%, imputed_prev: %.1f%%)
",
            100 * mean(is.na(out[pest == 1, sfp])),
            100 * mean(is.na(out[pest == 1, sfp_imputed])),
            100 * mean(is.na(out[pest == 1, sfp_imputed_prev]))))

dir.create(BUILD_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PATH <- file.path(BUILD_DIR, "disturbance_covariates.rds")
dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, OUT_PATH)
saveRDS(list(events = events, records = records,
             visits = dt[, .(id_pe, no_mes, year, origine, an_origine, perturb, an_perturb, fire_event, harv_event, pest_event)]),
        file.path(BUILD_DIR, "disturbance_events.rds"))   # diagnostics / reporting (R/tests/03v_attribute_changes.R)
cat(sprintf("Saved to %s
", OUT_PATH))
