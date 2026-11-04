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
#   6. Mark disturbance events (fire/harv/pest intensity via fcase)   [DONE]
#   6b. Reclassify/salvage disturbance type+year per plot/visit,
#       rebuild fire/harv/pest intensity from the salvage             [DONE]
#   7. Carry forward within epochs (origine resets epoch, perturb doesn't) [DONE]
#   8. Time since disturbance (tsf/tsh/sfp)                           [DONE]
#   9. Assemble output: final covariate shapes (fire/pest binary, harv
#      0/1/2), real + imputed tsf/tsh/sfp, id_pe/no_mes/year/total_ba  [DONE]
#  10. Diagnostics and write final dataset to disk                    [DONE]
#
# Output: data/build/disturbance_covariates.rds
# ==============================================================================

library(sf)
library(data.table)
library(here)

source(here::here("R", "config.R"))   # project paths (data/build etc.)
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
BA_LOSS_THRESH <- 0.25
BA_REGEN <- 5.0

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

cat("Applying drop-part-of-plot (CHT/DT truncation)...
")

droppart_ids <- unique(dt[origine %in% droppart_origine, id_pe])

truncate_one_plot <- function(g) {
  setorder(g, no_mes)
  first_event <- min(g[origine %in% droppart_origine, no_mes])
  pre <- g[no_mes < first_event]
  if (nrow(pre) < 2) return(NULL)
  for (i in seq_len(nrow(pre) - 1)) {
    ba0 <- pre$total_ba[i]
    ba1 <- pre$total_ba[i + 1]
    code_present    <- !is.na(pre$origine[i + 1]) || !is.na(pre$perturb[i + 1])
    crash            <- ba0 > 0 && ba1 <= (1 - BA_LOSS_THRESH) * ba0
    regen_transition <- ba0 >= BA_REGEN && ba1 < BA_REGEN
    if ((crash || regen_transition) && !code_present) return(NULL)
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
# STEP 6. MARK DISTURBANCE EVENTS
# origine match = stand-replacing (severity 2), perturb match = partial
# (severity 1), no match = 0. Three independent fcase() calls, one per
# disturbance type -- this is what lets CE (dual-coded in harv_p and pest_p)
# correctly set both harv_int and pest_int to 1 on the same row, with no
# conflict. Tested against BR/CT/ES/CE/no-code cases before writing this.
# ==============================================================================

cat("Marking disturbance events...
")

dt[, `:=`(
  fire_int = fcase(!is.na(origine) & origine %in% fire_o, 2L,
                    !is.na(perturb) & perturb %in% fire_p, 1L, default = 0L),
  harv_int = fcase(!is.na(origine) & origine %in% harv_o, 2L,
                    !is.na(perturb) & perturb %in% harv_p, 1L, default = 0L),
  pest_int = fcase(!is.na(origine) & origine %in% pest_o, 2L,
                    !is.na(perturb) & perturb %in% pest_p, 1L, default = 0L)
)]

# ==============================================================================
# STEP 6b. RECLASSIFY / SALVAGE MISCLASSIFIED OR MISSING-YEAR DISTURBANCES
#
# The raw origine/perturb classification has real data-quality problems:
# missing (an_origine/an_perturb NA) or future-dated years, and events
# misclassified as origine when they're really a late confirmation of an
# already-tracked perturb-level disturbance (or vice versa). This step
# rebuilds, per plot and per visit, which disturbance TYPE (fire/harv/pest)
# is currently active at the origine level and at the perturb level, and
# what year it started -- salvaging as much as can be soundly inferred
# before Step 7's epoch/carry-forward runs, so Step 8 only has to leave a
# case NA when it's genuinely unresolvable.
#
# KNOWN LIMITATION: this tracks at most ONE type per level per row. The raw
# data lets CE dual-code as both harv_p and pest_p on the same row (Step 6's
# three independent fcase() calls preserve that); this reclassification
# collapses CE to a single type (harv, by list priority order below) and so
# does NOT preserve CE's dual-coding. Flagged for a decision -- not silently
# handled.
#
# Algorithm (tested exhaustively against real data in Python/base-R before
# being written here as data.table; validated on id_pe 7608101001,
# 7000413801, 0119601201, 0301100102, 1219614902, 1419621501, 1419608402,
# 1619609002, 7100600302 -- see chat log for the full before/after tables):
#
#  For each plot, process visits in order, tracking:
#    so_flag/so_year/so_start -- the type/year/start-visit-year currently
#                                 tracked at the ORIGINE (stand-replacing)
#                                 level
#    sp_flag/sp_year/sp_start -- same, for the PERTURB (partial) level
#
#  At each visit, for whichever of origine/perturb has a code this row:
#    - same type already tracked at that level -> continuation: backfill
#      the year if it was unknown and the newly-seen year is <= that
#      level's start-visit-year (i.e. it doesn't postdate when the level
#      was first observed); a genuine year CONFLICT, or a candidate year
#      that postdates the start (temporally impossible as a continuation),
#      falls through to "new event" below.
#    - the type is currently tracked at the OTHER level (e.g. an origine
#      code appears for a type currently tracked at perturb) -> same
#      continuation-or-new-event logic, but checked against that OTHER
#      level's start-visit-year; if it continues, it's recorded as a
#      continuation of THAT level (a perturb-level code confirming an
#      already-origine-tracked type stays origine-level, and vice versa --
#      this is how a stale/misclassified level gets corrected).
#    - neither level currently tracks this type, or the above falls
#      through -> NEW EVENT:
#        - year unknown -> track type with year NA, starting this visit.
#        - year < this visit's own year -> unambiguous, already in the
#          past -> track immediately.
#        - year == this visit's own year -> AMBIGUOUS (did it happen
#          before or after this visit's measurement?) -- resolved via
#          basal area: this_ba < 5 (regeneration threshold) OR basal area
#          loss vs. the previous visit > 25% (the one documented
#          disturbance-severity threshold found in the PEP dictionary/MRNF
#          norm; NOT a general origine/perturb rule -- flagged as a
#          modeling assumption, not a cited rule) => BEFORE => track
#          immediately; otherwise AFTER => this visit keeps whatever was
#          already tracked, and the new event becomes effective starting
#          the NEXT visit (a one-visit-delayed "pending" event).
#      A NEW ORIGINE-level event resets (clears) the perturb-level slot;
#      perturb-level events never reset origine.
#
#  Backward-fill pass (runs after the forward pass, per plot): whenever a
#  slot's year becomes newly known at visit i, walk backward through
#  visits j < i while visit_year[j] >= the newly-known year:
#    - nothing tracked at j -> fill in (type, year) at the level the event
#      was found at (the "source" level).
#    - anything tracked at j (same type or different) with an unknown or
#      older year -> OVERWRITE it entirely with the source's type, level,
#      AND year -- introducing a resolved event at an earlier row always
#      replaces whatever that row had, including which level (origine vs
#      perturb) it's classified at; never downgrade a year that's already
#      as recent or more recent. If j's own type+year is already
#      as-recent-or-newer, stop -- j is already correctly governed by
#      something at least that recent.
#    Stop entirely once visit_year[j] < the newly-known year (the event
#    hadn't happened yet as of that visit). If visit_year[j] == the newly-
#    known year exactly, that visit's own BA trajectory decides (same
#    ba_before() check as the forward pass's ambiguous same-year case): if
#    it shows the event demonstrably HADN'T happened by that visit (BA
#    essentially unchanged from the prior visit), leave that row alone and
#    stop -- it, and everything earlier, predates the disturbance.
# ==============================================================================

cat("Reclassifying disturbance types (salvage pass)...
")

dt[, type_o := fcase(!is.na(origine) & origine %in% fire_o, "fire",
                      !is.na(origine) & origine %in% harv_o, "harv",
                      !is.na(origine) & origine %in% pest_o, "pest",
                      default = NA_character_)]
dt[, type_p := fcase(!is.na(perturb) & perturb %in% fire_p, "fire",
                      !is.na(perturb) & perturb %in% harv_p, "harv",
                      !is.na(perturb) & perturb %in% pest_p, "pest",
                      default = NA_character_)]

ba_before <- function(prev_ba, this_ba) {
  if (!is.na(this_ba) && this_ba < BA_REGEN) return(TRUE)
  if (!is.na(prev_ba) && prev_ba > 0) {
    if ((prev_ba - this_ba) / prev_ba > BA_LOSS_THRESH) return(TRUE)
  }
  FALSE
}

reclassify_plot <- function(year, total_ba, type_o, type_p, an_origine, an_perturb) {
  n <- length(year)
  so_flag <- NA_character_; so_year <- NA_integer_; so_start <- NA_integer_
  sp_flag <- NA_character_; sp_year <- NA_integer_; sp_start <- NA_integer_
  pending <- NULL

  carried_o_type <- rep(NA_character_, n); carried_o_year <- rep(NA_integer_, n)
  carried_p_type <- rep(NA_character_, n); carried_p_year <- rep(NA_integer_, n)

  # same_window(a, b): TRUE if no visit in this plot's own history falls
  # strictly between real years a and b -- i.e. there is no PEP measurement
  # available anywhere that could confirm a second real disturbance occurred
  # between them.
  same_window <- function(a, b) {
    lo <- min(a, b); hi <- max(a, b)
    !any(year > lo & year < hi)
  }

  # resolve_conflict(existing_year, new_year): the unified same-type
  # conflict decision, used for BOTH same-column and cross-column
  # conflicts. same window -> same disturbance (keep the higher year);
  # different window -> genuinely separate disturbance (fork a new,
  # independent track in whichever column just reported it).
  resolve_conflict <- function(existing_year, new_year) {
    if (same_window(existing_year, new_year)) {
      if (new_year > existing_year) "take_new" else "keep"
    } else "fork"
  }

  for (idx in seq_len(n)) {
    Y <- year[idx]
    prev_ba <- if (idx > 1) total_ba[idx - 1] else NA_real_
    this_ba <- total_ba[idx]

    if (!is.null(pending)) {
      if (pending$level == "o") {
        so_flag <- pending$type; so_year <- pending$year; so_start <- pending$year
        sp_flag <- NA_character_; sp_year <- NA_integer_; sp_start <- NA_integer_
      } else {
        sp_flag <- pending$type; sp_year <- pending$year; sp_start <- pending$year
      }
      pending <- NULL
    }

    t_o <- type_o[idx]; t_p <- type_p[idx]
    yr_o <- if (!is.na(t_o)) an_origine[idx] else NA_integer_
    yr_p <- if (!is.na(t_p)) an_perturb[idx] else NA_integer_
    if (!is.na(t_o) && !is.na(yr_o) && yr_o > Y) yr_o <- NA_integer_
    if (!is.na(t_p) && !is.na(yr_p) && yr_p > Y) yr_p <- NA_integer_

    this_o_type <- NA_character_

    new_event <- function(level, T, yr) {
      if (is.na(yr)) {
        if (level == "o") {
          so_flag <<- T; so_year <<- NA_integer_; so_start <<- Y
          sp_flag <<- NA_character_; sp_year <<- NA_integer_; sp_start <<- NA_integer_
        } else {
          sp_flag <<- T; sp_year <<- NA_integer_; sp_start <<- Y
        }
      } else if (yr < Y) {
        if (level == "o") {
          so_flag <<- T; so_year <<- yr; so_start <<- Y
          sp_flag <<- NA_character_; sp_year <<- NA_integer_; sp_start <<- NA_integer_
        } else {
          sp_flag <<- T; sp_year <<- yr; sp_start <<- Y
        }
      } else {
        if (ba_before(prev_ba, this_ba)) {
          if (level == "o") {
            so_flag <<- T; so_year <<- yr; so_start <<- Y
            sp_flag <<- NA_character_; sp_year <<- NA_integer_; sp_start <<- NA_integer_
          } else {
            sp_flag <<- T; sp_year <<- yr; sp_start <<- Y
          }
        } else {
          pending <<- list(level = level, type = T, year = yr)
        }
      }
    }

    if (!is.na(t_o)) {
      if (!is.na(so_flag) && so_flag == t_o) {
        if (is.na(yr_o)) {
        } else if (is.na(so_year)) {
          if (yr_o <= so_start) so_year <- yr_o else new_event("o", t_o, yr_o)
        } else if (so_year != yr_o) {
          action <- resolve_conflict(so_year, yr_o)
          if (action == "keep") {
          } else if (action == "take_new") so_year <- yr_o
          else new_event("o", t_o, yr_o)
        }
      } else if (!is.na(sp_flag) && sp_flag == t_o) {
        if (is.na(yr_o)) {
        } else if (is.na(sp_year)) {
          if (yr_o <= sp_start) sp_year <- yr_o else new_event("o", t_o, yr_o)
        } else if (sp_year != yr_o) {
          if (is.na(t_p)) {
            action <- resolve_conflict(sp_year, yr_o)
            if (action == "keep") {
            } else if (action == "take_new") {
              sp_year <- yr_o  # stays in perturb (already-resolved column), no move
            } else {
              new_event("o", t_o, yr_o)
            }
          } else {
            new_event("o", t_o, yr_o)
          }
        }
      } else new_event("o", t_o, yr_o)
      this_o_type <- so_flag
    }

    if (!is.na(t_p)) {
      if (!is.na(sp_flag) && sp_flag == t_p) {
        if (is.na(yr_p)) {
        } else if (is.na(sp_year)) {
          if (yr_p <= sp_start) sp_year <- yr_p else new_event("p", t_p, yr_p)
        } else if (sp_year != yr_p) {
          action <- resolve_conflict(sp_year, yr_p)
          if (action == "keep") {
          } else if (action == "take_new") sp_year <- yr_p
          else new_event("p", t_p, yr_p)
        }
      } else if (!is.na(so_flag) && so_flag == t_p &&
                 !(!is.na(this_o_type) && this_o_type == t_p)) {
        if (is.na(yr_p)) {
        } else if (is.na(so_year)) {
          if (yr_p <= so_start) so_year <- yr_p else new_event("p", t_p, yr_p)
        } else if (so_year != yr_p) {
          if (is.na(t_o)) {
            action <- resolve_conflict(so_year, yr_p)
            if (action == "keep") {
            } else if (action == "take_new") {
              so_year <- yr_p  # stays in origine (already-resolved column), no move
            } else {
              new_event("p", t_p, yr_p)
            }
          } else {
            new_event("p", t_p, yr_p)
          }
        }
      } else new_event("p", t_p, yr_p)
    }

    carried_o_type[idx] <- so_flag; carried_o_year[idx] <- so_year
    carried_p_type[idx] <- sp_flag; carried_p_year[idx] <- sp_year
  }

  # Backward-fill: two FULLY INDEPENDENT per-level passes. Each only ever
  # looks at its own level's state -- no cross-level magnitude battles --
  # which is what lets two genuinely independent, non-conflicting tracks
  # (e.g. fire in perturb + harv in origine) coexist correctly across every
  # row.
  backfill_one_level <- function(typearr, yrarr) {
    if (n <= 1) return(list(type = typearr, yr = yrarr))
    for (i in 2:n) {
      snap_type <- typearr; snap_year <- yrarr
      if (is.na(snap_type[i]) || is.na(snap_year[i])) next
      newly_known <- is.na(snap_type[i - 1]) || snap_type[i - 1] != snap_type[i] ||
        is.na(snap_year[i - 1]) || snap_year[i - 1] != snap_year[i]
      if (!newly_known) next
      T <- snap_type[i]; Y_new <- snap_year[i]; j <- i - 1
      while (j >= 1) {
        if (year[j] < Y_new) break
        if (year[j] == Y_new) {
          prev_ba_j <- if (j > 1) total_ba[j - 1] else NA_real_
          this_ba_j <- total_ba[j]
          if (!ba_before(prev_ba_j, this_ba_j)) break
        }
        cur_T <- typearr[j]; cur_Y <- yrarr[j]
        if (is.na(cur_T)) {
          typearr[j] <- T; yrarr[j] <- Y_new
        } else if (cur_T == T) {
          if (is.na(cur_Y) || cur_Y < Y_new) { typearr[j] <- T; yrarr[j] <- Y_new }
          else if (cur_Y != Y_new) break
        } else {
          # different type, no year -- can't tell if it predates or postdates
          # the new event, so leave this row alone, but keep walking backward
          # past it (do not break the loop).
          if (is.na(cur_Y)) {
          } else if (cur_Y < Y_new) { typearr[j] <- T; yrarr[j] <- Y_new }
          else break
        }
        j <- j - 1
      }
    }
    list(type = typearr, yr = yrarr)
  }

  res_o <- backfill_one_level(carried_o_type, carried_o_year)
  carried_o_type <- res_o$type; carried_o_year <- res_o$yr
  res_p <- backfill_one_level(carried_p_type, carried_p_year)
  carried_p_type <- res_p$type; carried_p_year <- res_p$yr

  # Consolidation: the two backward-fills run independently and can
  # coincidentally land on the SAME type in both columns at the same row.
  # When that happens it's the same same-type-crossing-columns situation
  # the forward pass already reconciles -- merge instead of leaving two
  # redundant copies, UNLESS the type was ever raw-reported in BOTH
  # columns on the SAME row anywhere in the plot (genuine independent
  # dual-coding, never merged), or the two years are real and fall in
  # DIFFERENT windows (a genuine, already-detected fork -- leave both).
  raw_dual_coded_for_type <- function(T) {
    any(type_o == T & type_p == T, na.rm = TRUE)
  }

  # home_column(T): whichever column tracked this type FIRST (earliest
  # row) across the whole plot -- a merged event always belongs there,
  # not wherever the collision happened to be detected.
  home_column <- function(T) {
    idx_o <- which(carried_o_type == T)
    idx_p <- which(carried_p_type == T)
    first_o <- if (length(idx_o) > 0) idx_o[1] else NA_integer_
    first_p <- if (length(idx_p) > 0) idx_p[1] else NA_integer_
    if (is.na(first_o)) return("p")
    if (is.na(first_p)) return("o")
    if (first_o <= first_p) "o" else "p"
  }

  # track_start_year(typearr, T): the visit year this column's track for
  # type T was first observed (earliest row carrying it) -- an unresolved
  # (NA) year can never truly be LATER than this, since the disturbance
  # was already visible by then. Used to chronology-check a merge against
  # an NA side, mirroring the forward pass's own yr <= start ADOPT test.
  track_start_year <- function(typearr, T) {
    idx <- which(typearr == T)
    if (length(idx) > 0) year[idx[1]] else NA_integer_
  }

  for (i in seq_len(n)) {
    if (!is.na(carried_o_type[i]) && !is.na(carried_p_type[i]) &&
        carried_o_type[i] == carried_p_type[i]) {
      if (raw_dual_coded_for_type(carried_o_type[i])) next
      oy <- carried_o_year[i]; py <- carried_p_year[i]
      T <- carried_o_type[i]
      skip <- FALSE
      if (!is.na(oy) && !is.na(py)) {
        if (!same_window(oy, py)) skip <- TRUE
      } else if (is.na(oy) && !is.na(py)) {
        o_start <- track_start_year(carried_o_type, T)
        if (!is.na(o_start) && py > o_start) skip <- TRUE
        # py postdates origine's own first observation of this type --
        # can't be the same (still-unresolved) event that was already
        # visible at o_start; leave as a genuine fork
      } else if (is.na(py) && !is.na(oy)) {
        p_start <- track_start_year(carried_p_type, T)
        if (!is.na(p_start) && oy > p_start) skip <- TRUE
        # symmetric chronology check
      }
      if (skip) next
      merged_year <- if (!is.na(oy)) oy else NA_integer_
      if (!is.na(py) && (is.na(merged_year) || py > merged_year)) merged_year <- py
      dest <- home_column(T)
      if (dest == "o") {
        carried_o_type[i] <- T; carried_o_year[i] <- merged_year
        carried_p_type[i] <- NA_character_; carried_p_year[i] <- NA_integer_
      } else {
        carried_p_type[i] <- T; carried_p_year[i] <- merged_year
        carried_o_type[i] <- NA_character_; carried_o_year[i] <- NA_integer_
      }
    }
  }

  # Final reset gate: wherever origine's resolved year is strictly later
  # than perturb's resolved year at that row, the origine event
  # supersedes the earlier perturb history -- suppress perturb there.
  # Equal years are left untouched (both survive, same-event case).
  for (i in seq_len(n)) {
    if (is.na(carried_p_type[i]) || is.na(carried_p_year[i])) next
    if (!is.na(carried_o_year[i]) && carried_o_year[i] > carried_p_year[i]) {
      carried_p_type[i] <- NA_character_
      carried_p_year[i] <- NA_integer_
    }
  }

  list(carried_o_type = carried_o_type, carried_o_year = carried_o_year,
       carried_p_type = carried_p_type, carried_p_year = carried_p_year)
}

recl <- dt[, reclassify_plot(year, total_ba, type_o, type_p, an_origine, an_perturb), by = id_pe]
dt[, `:=`(carried_o_type = recl$carried_o_type, carried_o_year = recl$carried_o_year,
          carried_p_type = recl$carried_p_type, carried_p_year = recl$carried_p_year)]

# rebuild fire/harv/pest intensity from the SALVAGED classification, not the
# raw origine/perturb columns -- this is what lets Step 7's epoch/carry-
# forward and Step 8's tsf/tsh/sfp benefit from the corrections above.
dt[, `:=`(
  fire_int = fcase(carried_o_type == "fire", 2L, carried_p_type == "fire", 1L, default = 0L),
  harv_int = fcase(carried_o_type == "harv", 2L, carried_p_type == "harv", 1L, default = 0L),
  pest_int = fcase(carried_o_type == "pest", 2L, carried_p_type == "pest", 1L, default = 0L)
)]

# ==============================================================================
# STEP 7. CARRY FORWARD WITHIN EPOCHS
# epoch = shared clock across fire/harv/pest: any origine-level (stand-
# replacing) event of ANY type starts a new epoch and resets all three
# covariates, not just the type that caused it. has_origine now comes from
# the SALVAGED classification (carried_o_type), not the raw origine column,
# so a promotion/demotion from Step 6b correctly shifts where the epoch
# resets. perturb-level (partial) events never start a new epoch.
# ==============================================================================

cat("Carrying forward within epochs...
")

dt[, has_origine := !is.na(carried_o_type)]
dt[, epoch := cumsum(has_origine), by = id_pe]

dt[, `:=`(
  cov_fire    = nafill(fifelse(fire_int > 0, fire_int, NA_integer_), "locf", fill = 0L),
  cov_harvest = nafill(fifelse(harv_int > 0, harv_int, NA_integer_), "locf", fill = 0L),
  cov_pest    = nafill(fifelse(pest_int > 0, pest_int, NA_integer_), "locf", fill = 0L)
), by = .(id_pe, epoch)]

# ==============================================================================
# STEP 7b. FALL BACK TO OBSERVATION YEAR FOR UNRESOLVED, NON-FIRST-OBSERVATION
# DISTURBANCE YEARS -- IMPUTED VARIANT (kept fully separate from the real
# carried_o_year/carried_p_year; nothing above this point is touched)
#
# Diagnostic run (2026-09-05) on the post-Step-6b classification found:
#   11,289 rows (of 34,002 disturbed-plot rows) have an active carried_o_type
#   or carried_p_type with a NA year; 9,358 of those (83%) are NOT the plot's
#   first observation -- i.e. the gap persists even though an earlier visit
#   existed to (potentially) backfill from, meaning the raw an_origine/
#   an_perturb code itself was simply never dated at any visit where it
#   appears. Split by level/type (rows, not-first-observation count):
#     origine fire: 1526 (1153 not-first)   perturb fire: 32   (29 not-first)
#     origine harv: 3008 (2410 not-first)   perturb harv: 4832 (3952 not-first)
#     origine pest:  299  (258 not-first)   perturb pest: 1755 (1556 not-first)
#
# Rule (approved -- fallback only when NOT the first observation): a row
# where the active disturbance type's year is NA AND this is not the plot's
# first observation falls back to that row's own visit year (`year`). Not
# the first observation means the disturbance was necessarily first
# RECORDED at this visit (not silently carried from an unmeasured prior
# state), so the visit year is a defensible anchor. A first-observation gap
# is left NA -- there's no lower bound on how long the disturbance predates
# the start of monitoring, so no fallback is defensible there.
#
# Kept as separate *_imputed columns (carried_o_year_imputed,
# carried_p_year_imputed, and downstream fire_yr_imputed / harv_yr_imputed /
# pest_yr_imputed / *_carried_imputed / tsf_imputed / tsh_imputed /
# sfp_imputed below) rather than overwriting the real ones, so real vs.
# imputed can always be told apart and compared directly.
#
# A SECOND imputed variant (added 2026-09-05, `_imputed_prev`): the row's own
# visit year is only one defensible anchor -- it's the LATEST the disturbance
# could have started (it was first recorded then), but the disturbance could
# equally have started any time back to the PREVIOUS visit (the last time
# this plot was checked and the type wasn't yet flagged). `_imputed_prev`
# uses that previous visit's year instead, giving the other end of the same
# interval: `_imputed` is the minimum-plausible time-since-disturbance bound,
# `_imputed_prev` is the maximum-plausible bound. Verified against real data
# before writing this: 0 negative time-since-disturbance values, and
# tsf_imputed_prev/tsh_imputed_prev/sfp_imputed_prev >= their `_imputed`
# counterpart in every one of the 5571/15342/4351 rows where both are
# defined (the interval ordering holds, never inverted).
# ==============================================================================

cat("Imputing fallback years for unresolved non-first-observation disturbances...\n")

dt[, first_obs := seq_len(.N) == 1L, by = id_pe]

dt[, carried_o_year_imputed := carried_o_year]
dt[!is.na(carried_o_type) & is.na(carried_o_year) & !first_obs,
   carried_o_year_imputed := year]

dt[, carried_p_year_imputed := carried_p_year]
dt[!is.na(carried_p_type) & is.na(carried_p_year) & !first_obs,
   carried_p_year_imputed := year]

n_imputed_o <- sum(!is.na(dt$carried_o_type) & is.na(dt$carried_o_year) & !dt$first_obs)
n_imputed_p <- sum(!is.na(dt$carried_p_type) & is.na(dt$carried_p_year) & !dt$first_obs)
cat(sprintf("  imputed %d origine-level and %d perturb-level year(s)\n",
            n_imputed_o, n_imputed_p))

dt[, prev_year := shift(year, 1), by = id_pe]

dt[, carried_o_year_imputed_prev := carried_o_year]
dt[!is.na(carried_o_type) & is.na(carried_o_year) & !first_obs,
   carried_o_year_imputed_prev := prev_year]

dt[, carried_p_year_imputed_prev := carried_p_year]
dt[!is.na(carried_p_type) & is.na(carried_p_year) & !first_obs,
   carried_p_year_imputed_prev := prev_year]

# ==============================================================================
# STEP 8. TIME SINCE DISTURBANCE (tsf / tsh / sfp), REAL AND IMPUTED
# fire and pest collapse to binary (any intensity > 0 counts); harvest keeps
# its 0/1/2 intensity per the project's original covariate design. The year
# used for each covariate is whatever Step 6b's salvage resolved for the
# CURRENTLY active event (carried_o_year if the active level is origine,
# carried_p_year if perturb) -- carried forward within the epoch/type the
# same way cov_fire/cov_harvest/cov_pest are, and left NA where Step 6b
# genuinely couldn't resolve a year (kept NA rather than fabricated, per
# INLA's per-covariate NA handling). The *_imputed columns below run the
# identical logic against Step 7b's imputed years instead, side by side.
# ==============================================================================

cat("Computing time-since-disturbance covariates (real and imputed)...
")

# Each type's year comes ONLY from that type's own carried_o_year/carried_p_year --
# no cross-type priority. A single shared `event_year` variable used to sit here,
# resolved by a fire > harv > pest priority order; whenever origine and perturb were
# simultaneously set to TWO DIFFERENT types (403 such rows in the real data), that
# priority order silently leaked the higher-priority type's year into the lower-
# priority type's *_yr column -- e.g. a harv-tracked row got fire's year assigned to
# harv_yr. Confirmed as a real bug (246 contaminated *_yr assignments), not a design
# choice -- nobody approved a cross-type hierarchy, and there's no reason for one:
# fire_yr/harv_yr/pest_yr are three unrelated columns and must each look only at
# their own type's slots. Fixed 2026-09-05.
#
# Separately: when origine AND perturb are BOTH resolved to the SAME type (fire:
# never happens, 0 rows; pest: 3 rows; harv: 639 rows, mostly harv/harv with two
# different years), the MORE RECENT of the two years is used (approved 2026-09-05)
# -- e.g. an older origine-level harvest doesn't hide a genuinely more recent
# perturb-level harvest of the same type. A same-type row where only one side has a
# known year (the other NA) simply uses the known one, rather than the old code's
# type-only match returning NA even when the other slot had a real year.
type_year <- function(otype, oyear, ptype, pyear, target) {
  o_is <- !is.na(otype) & otype == target
  p_is <- !is.na(ptype) & ptype == target
  fcase(
    o_is & p_is & !is.na(oyear) & !is.na(pyear), as.integer(pmax(oyear, pyear)),
    o_is & p_is & !is.na(oyear), oyear,
    o_is & p_is & !is.na(pyear), pyear,
    o_is, oyear,
    p_is, pyear,
    default = NA_integer_
  )
}

dt[, fire_yr := type_year(carried_o_type, carried_o_year, carried_p_type, carried_p_year, "fire")]
dt[, harv_yr := type_year(carried_o_type, carried_o_year, carried_p_type, carried_p_year, "harv")]
dt[, pest_yr := type_year(carried_o_type, carried_o_year, carried_p_type, carried_p_year, "pest")]
dt[, fire_yr_imputed := type_year(carried_o_type, carried_o_year_imputed, carried_p_type, carried_p_year_imputed, "fire")]
dt[, harv_yr_imputed := type_year(carried_o_type, carried_o_year_imputed, carried_p_type, carried_p_year_imputed, "harv")]
dt[, pest_yr_imputed := type_year(carried_o_type, carried_o_year_imputed, carried_p_type, carried_p_year_imputed, "pest")]
dt[, fire_yr_imputed_prev := type_year(carried_o_type, carried_o_year_imputed_prev, carried_p_type, carried_p_year_imputed_prev, "fire")]
dt[, harv_yr_imputed_prev := type_year(carried_o_type, carried_o_year_imputed_prev, carried_p_type, carried_p_year_imputed_prev, "harv")]
dt[, pest_yr_imputed_prev := type_year(carried_o_type, carried_o_year_imputed_prev, carried_p_type, carried_p_year_imputed_prev, "pest")]

dt[, fire_yr := fifelse(fire_int > 0, fire_yr, NA_integer_)]
dt[, harv_yr := fifelse(harv_int > 0, harv_yr, NA_integer_)]
dt[, pest_yr := fifelse(pest_int > 0, pest_yr, NA_integer_)]
dt[, fire_yr_imputed := fifelse(fire_int > 0, fire_yr_imputed, NA_integer_)]
dt[, harv_yr_imputed := fifelse(harv_int > 0, harv_yr_imputed, NA_integer_)]
dt[, pest_yr_imputed := fifelse(pest_int > 0, pest_yr_imputed, NA_integer_)]
dt[, fire_yr_imputed_prev := fifelse(fire_int > 0, fire_yr_imputed_prev, NA_integer_)]
dt[, harv_yr_imputed_prev := fifelse(harv_int > 0, harv_yr_imputed_prev, NA_integer_)]
dt[, pest_yr_imputed_prev := fifelse(pest_int > 0, pest_yr_imputed_prev, NA_integer_)]

# sub-epoch grouping, as before: increments at each row with an active event
# of that type, so a stale known year can't leak forward past a newer
# same-type event whose year is unknown. Depends only on intensity (not on
# year), so it's identical for the real and imputed variants -- computed once.
dt[, fire_grp := cumsum(fire_int > 0), by = .(id_pe, epoch)]
dt[, harv_grp := cumsum(harv_int > 0), by = .(id_pe, epoch)]
dt[, pest_grp := cumsum(pest_int > 0), by = .(id_pe, epoch)]

dt[, fire_yr_carried := fire_yr[1], by = .(id_pe, epoch, fire_grp)]
dt[, harv_yr_carried := harv_yr[1], by = .(id_pe, epoch, harv_grp)]
dt[, pest_yr_carried := pest_yr[1], by = .(id_pe, epoch, pest_grp)]
dt[, fire_yr_carried_imputed := fire_yr_imputed[1], by = .(id_pe, epoch, fire_grp)]
dt[, harv_yr_carried_imputed := harv_yr_imputed[1], by = .(id_pe, epoch, harv_grp)]
dt[, pest_yr_carried_imputed := pest_yr_imputed[1], by = .(id_pe, epoch, pest_grp)]
dt[, fire_yr_carried_imputed_prev := fire_yr_imputed_prev[1], by = .(id_pe, epoch, fire_grp)]
dt[, harv_yr_carried_imputed_prev := harv_yr_imputed_prev[1], by = .(id_pe, epoch, harv_grp)]
dt[, pest_yr_carried_imputed_prev := pest_yr_imputed_prev[1], by = .(id_pe, epoch, pest_grp)]

# cov_fire/cov_pest KEEP their full 0/1/2 epoch-carried severity here --
# do not collapse in place. fire_bin/pest_bin below are an ADDITIONAL view
# (needed now, for tsf/sfp's "is this type currently active at all" check),
# not a replacement. Whether the FINAL saved output keeps cov_fire/cov_pest
# as 0/1/2, exposes only the binary view, or both, is an output-assembly
# decision that belongs in Step 9 -- not decided here by destroying the
# intermediate severity value. (An earlier version of this step overwrote
# cov_fire/cov_pest in place; reverted -- that would have made the severity
# signal unrecoverable for anything written later without redoing the whole
# epoch/locf carry-forward from fire_int/pest_int again.)
dt[, `:=`(
  fire_bin = as.integer(cov_fire > 0),
  pest_bin = as.integer(cov_pest > 0)
)]

dt[, `:=`(
  tsf = fifelse(fire_bin == 1L, year - fire_yr_carried, NA_integer_),
  tsh = fifelse(cov_harvest > 0L,  year - harv_yr_carried, NA_integer_),
  sfp = fifelse(pest_bin == 1L, year - pest_yr_carried, NA_integer_)
)]
dt[, `:=`(
  tsf_imputed = fifelse(fire_bin == 1L, year - fire_yr_carried_imputed, NA_integer_),
  tsh_imputed = fifelse(cov_harvest > 0L,  year - harv_yr_carried_imputed, NA_integer_),
  sfp_imputed = fifelse(pest_bin == 1L, year - pest_yr_carried_imputed, NA_integer_)
)]
dt[, `:=`(
  tsf_imputed_prev = fifelse(fire_bin == 1L, year - fire_yr_carried_imputed_prev, NA_integer_),
  tsh_imputed_prev = fifelse(cov_harvest > 0L,  year - harv_yr_carried_imputed_prev, NA_integer_),
  sfp_imputed_prev = fifelse(pest_bin == 1L, year - pest_yr_carried_imputed_prev, NA_integer_)
)]

cat(sprintf("  tsf NA rate (of fire_bin==1):      %.1f%%  (imputed: %.1f%%, imputed_prev: %.1f%%)\n",
            100 * mean(is.na(dt[fire_bin == 1, tsf])),
            100 * mean(is.na(dt[fire_bin == 1, tsf_imputed])),
            100 * mean(is.na(dt[fire_bin == 1, tsf_imputed_prev]))))
cat(sprintf("  tsh NA rate (of cov_harvest>0):    %.1f%%  (imputed: %.1f%%, imputed_prev: %.1f%%)\n",
            100 * mean(is.na(dt[cov_harvest > 0, tsh])),
            100 * mean(is.na(dt[cov_harvest > 0, tsh_imputed])),
            100 * mean(is.na(dt[cov_harvest > 0, tsh_imputed_prev]))))
cat(sprintf("  sfp NA rate (of pest_bin==1):      %.1f%%  (imputed: %.1f%%, imputed_prev: %.1f%%)\n",
            100 * mean(is.na(dt[pest_bin == 1, sfp])),
            100 * mean(is.na(dt[pest_bin == 1, sfp_imputed])),
            100 * mean(is.na(dt[pest_bin == 1, sfp_imputed_prev]))))

# ==============================================================================
# STEP 8b. SANITY CHECKS -- IMPUTED YEARS MUST NOT CONTRADICT REAL ONES
#
# Check 1: negative time-since-disturbance using the imputed years would mean
# an imputed anchor year is somehow AFTER the visit it's supposed to explain
# -- structurally shouldn't happen (the anchor is always the first row of its
# subgroup, i.e. chronologically earliest), but verified directly rather than
# assumed.
#
# Check 2: within any subgroup (id_pe, epoch, fire/harv/pest_grp) whose
# ANCHOR row's year was imputed (i.e. not real), does any LATER row in that
# same subgroup carry a REAL (non-imputed) year for that type that differs
# from the anchor's imputed value? That would mean the imputed placeholder
# is contradicted by actual evidence later in the same uninterrupted episode
# -- if so, the imputed columns are unsound for those rows and must not be
# used without addressing this. Verified clean against real data
# (2026-09-05, 34,002-row disturbed-plot subset) before this was written:
# 0 negative values and 0 contradicted subgroups across fire/harv/pest.
# ==============================================================================

cat("Checking imputed years for internal consistency...\n")

n_neg <- sum(dt$tsf_imputed < 0, na.rm = TRUE) +
         sum(dt$tsh_imputed < 0, na.rm = TRUE) +
         sum(dt$sfp_imputed < 0, na.rm = TRUE)
cat(sprintf("  negative time-since-disturbance using imputed years: %d\n", n_neg))
if (n_neg > 0) {
  warning("Imputed years produced negative time-since-disturbance -- investigate before using *_imputed columns.")
}

n_neg_prev <- sum(dt$tsf_imputed_prev < 0, na.rm = TRUE) +
              sum(dt$tsh_imputed_prev < 0, na.rm = TRUE) +
              sum(dt$sfp_imputed_prev < 0, na.rm = TRUE)
cat(sprintf("  negative time-since-disturbance using imputed_prev years: %d\n", n_neg_prev))
if (n_neg_prev > 0) {
  warning("Imputed_prev years produced negative time-since-disturbance -- investigate before using *_imputed_prev columns.")
}

check_anchor_conflicts <- function(yr_col, yr_imputed_col, grp_col) {
  dt[, {
    anchor_yr <- get(yr_col)[1]
    anchor_yr_imputed <- get(yr_imputed_col)[1]
    anchor_was_imputed <- is.na(anchor_yr) && !is.na(anchor_yr_imputed)
    if (!anchor_was_imputed || .N < 2) {
      .(conflict = FALSE)
    } else {
      later_real <- get(yr_col)[-1]
      later_real <- later_real[!is.na(later_real)]
      .(conflict = length(later_real) > 0 && any(later_real != anchor_yr_imputed))
    }
  }, by = c("id_pe", "epoch", grp_col)]
}

conf_fire <- check_anchor_conflicts("fire_yr", "fire_yr_imputed", "fire_grp")
conf_harv <- check_anchor_conflicts("harv_yr", "harv_yr_imputed", "harv_grp")
conf_pest <- check_anchor_conflicts("pest_yr", "pest_yr_imputed", "pest_grp")
n_conflict <- sum(conf_fire$conflict) + sum(conf_harv$conflict) + sum(conf_pest$conflict)
cat(sprintf("  subgroups where an imputed anchor conflicts with a later real year: %d\n", n_conflict))
if (n_conflict > 0) {
  warning("An imputed anchor year is contradicted by a real year later in the same episode -- investigate before using *_imputed columns.")
}

conf_fire_prev <- check_anchor_conflicts("fire_yr", "fire_yr_imputed_prev", "fire_grp")
conf_harv_prev <- check_anchor_conflicts("harv_yr", "harv_yr_imputed_prev", "harv_grp")
conf_pest_prev <- check_anchor_conflicts("pest_yr", "pest_yr_imputed_prev", "pest_grp")
n_conflict_prev <- sum(conf_fire_prev$conflict) + sum(conf_harv_prev$conflict) + sum(conf_pest_prev$conflict)
cat(sprintf("  subgroups where an imputed_prev anchor conflicts with a later real year: %d\n", n_conflict_prev))
if (n_conflict_prev > 0) {
  warning("An imputed_prev anchor year is contradicted by a real year later in the same episode -- investigate before using *_imputed_prev columns.")
}

# ==============================================================================
# STEP 9. ASSEMBLE OUTPUT 
#
# Final covariate shapes, per the project's original covariate design and
# the 2026-09-05 decisions above: fire and pest collapse to binary
# (fire_bin/pest_bin), harvest keeps its 0/1/2 intensity (cov_harvest).
# total_ba is kept (2026-09-05: useful diagnostic/modeling reference, even
# though it was only ever a working input for the disturbance classification
# itself, not a disturbance covariate in its own right).
#
# tsf/tsh/sfp is the real value where Step 6b could resolve one.
# tsf_imputed/tsh_imputed/sfp_imputed falls back to this visit's own year
# when unresolved and not the first observation -- the minimum-plausible
# time-since-disturbance bound (the disturbance could not have been recorded
# any later than this visit).
# tsf_imputed_prev/tsh_imputed_prev/sfp_imputed_prev (added 2026-09-05) falls
# back to the PREVIOUS visit's year instead -- the maximum-plausible bound
# (the disturbance could not have started any earlier than the last time the
# plot was checked and the type wasn't yet flagged). Together the two
# imputed variants bracket the true, unrecorded disturbance year.
# All three variants are kept side by side, clearly labeled -- the Step 7b
# imputation approach is still pending your supervisor's sign-off, so
# nothing is decided here by omission; downstream code chooses which to use.
#
# fire/harv/pest are renamed from their working names (fire_bin/cov_harvest/
# pest_bin) to the final, clean covariate names.
# ==============================================================================

cat("Assembling output...
")

out <- dt[, .(
  id_pe, no_mes, year, total_ba,
  fire = fire_bin,
  harv = cov_harvest,
  pest = pest_bin,
  tsf, tsh, sfp,
  tsf_imputed, tsh_imputed, sfp_imputed,
  tsf_imputed_prev, tsh_imputed_prev, sfp_imputed_prev
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
cat(sprintf("Saved to %s
", OUT_PATH))
