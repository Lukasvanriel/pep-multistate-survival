# ==============================================================================
# R/03a_disturbance_events.R -- disturbance event logic (sourced by 03, not run).
# Written 2026-10-26; replaces 03 steps 6b-8 (forward pass + back-fill + consolidation + reset gate, epochs,
# per-row imputed years). Decisions behind every rule: docs/PIPELINE_NOTES.md, "Disturbance event rewrite".
#
# Fixed order, each step only uses the results of the previous one:
#   E0 ba_event   per visit: BA loss >= BA_LOSS_THRESH since the previous visit, or a drop from >= BA_REGEN to < BA_REGEN.
#                 At a plot's first visit (no previous visit): BA < BA_REGEN.
#   E1 records    one record per code per visit: type (fire/harv/pest), level (2 = stand-origin / origine, 1 = partial /
#                 perturb), year (NA if absent or later than the visit). A code in two type lists (CE) gives two records.
#                 Within a visit: origine records first, then perturb records.
#   E2 events     per plot, records in visit order; each record joins the current event of its type or starts a new one:
#                 dated record:   current event dated  -> same window (no visit strictly between the two years) = same event,
#                                                         the latest year sets year and level (tie: the more severe level);
#                                                         other window = new event.
#                                 current event undated (first seen at visit F) -> year <= year(F) = same event, the dated
#                                                         record sets year and level; later year = new event.
#                 undated record: current event first seen at this visit -> same event (level = the more severe);
#                                 current event ended by a stand-origin event seen by now (see E4) -> new event;
#                                 ba_event at this visit -> new event (new year = lower-bound dating, its own level);
#                                 otherwise -> same event (old year, old level).
#                 "Current event" of a type = its event with the latest date (date = year, or year of the first-seen visit).
#   E3 dating     dated event: year. Undated event first seen at visit F: lower bound = year(F) (TSD min), upper bound =
#                 year of the visit before F (TSD max), midpoint. First seen at the plot's first visit: no bounds.
#   E4 timeline   start visit of an event = the first visit at which it is visible: dated: first visit with year(v) > year,
#                 or year(v) == year and ba_event(v) [same-year rule: did the event happen before the measurement?];
#                 undated: F. Back-fill = a dated event starts at earlier visits than the one that recorded it, from its
#                 year on (chronologically possible only). First-visit rule (T3d): an undated PARTIAL event first seen at
#                 the plot's first visit never starts.
#                 Active at v: v >= start, and not ended. Ended at v: a stand-origin event (any type) with a LATER start
#                 visit has started by v. Events with the same start visit (in the same interval between two visits,
#                 incl. "before the first visit") do not end each other (a partial + a stand-origin event: both kept),
#                 EXCEPT stand-origin events of different types with the same start visit (2026-10-26, user: treated as
#                 salvage logging / one stand origin): only ONE is kept, the others are dropped (never active):
#                   all dated -> the latest year (same year: harvest > fire > pest);
#                   any undated -> harvest if present; else (fire / pest) the dated one, else fire.
#   E5 visit      per visit and type, the active event with the latest date (tie: more severe level) sets the flag
#                 (fire/pest 0/1, harv 0/1/2) and the TSD: real (dated only), min, max, mid; dating marker
#                 "dated" / "bounded" / "first_visit" (stand-origin, undated at the first visit -> 04 stand-age rule).
# ==============================================================================

DIST_TYPES <- c("fire", "harv", "pest")

ba_event_flag <- function(total_ba, ba_regen, loss_thresh) {
  n <- length(total_ba)
  if (n == 0L) return(logical(0))
  prev <- c(NA_real_, total_ba[-n])
  out <- (!is.na(prev) & prev > 0 & total_ba <= (1 - loss_thresh) * prev) |
         (!is.na(prev) & prev >= ba_regen & total_ba < ba_regen)
  out[1L] <- total_ba[1L] < ba_regen
  out
}

# E1. codes: named list fire_o, fire_p, harv_o, harv_p, pest_o, pest_p (code vectors from 03 step 2)
make_records <- function(year, origine, an_origine, perturb, an_perturb, codes) {
  v <- integer(0); typ <- character(0); lev <- integer(0); yr <- integer(0)
  yr_ok <- function(y, vy) if (is.na(y) || y > vy) NA_integer_ else as.integer(y)
  for (i in seq_along(year)) {
    if (!is.na(origine[i])) for (T in DIST_TYPES) if (origine[i] %in% codes[[paste0(T, "_o")]]) {
      v <- c(v, i); typ <- c(typ, T); lev <- c(lev, 2L); yr <- c(yr, yr_ok(an_origine[i], year[i])) }
    if (!is.na(perturb[i])) for (T in DIST_TYPES) if (perturb[i] %in% codes[[paste0(T, "_p")]]) {
      v <- c(v, i); typ <- c(typ, T); lev <- c(lev, 1L); yr <- c(yr, yr_ok(an_perturb[i], year[i])) }
  }
  list(v = v, type = typ, level = lev, yr = yr)
}

# E2 + E3. Returns events (type, level, yr, first) and, per record, the event id and the action taken.
# Choose the one stand-origin event kept among events of different types with the same start visit (ids = event indices).
SO_PRIORITY <- c(harv = 3L, fire = 2L, pest = 1L)
choose_stand_origin <- function(ids, type, yr) {
  pr <- SO_PRIORITY[type[ids]]
  if (all(!is.na(yr[ids]))) return(ids[order(-yr[ids], -pr)][1L])
  if (any(type[ids] == "harv")) { h <- ids[type[ids] == "harv"]; return(h[order(-ifelse(is.na(yr[h]), -Inf, yr[h]))][1L]) }
  d <- ids[!is.na(yr[ids])]
  if (length(d)) return(d[order(-yr[d], -SO_PRIORITY[type[d]])][1L])
  ids[order(-pr)][1L]
}

start_visit <- function(yr, level, first, year, ba_ev) {
  if (!is.na(yr)) { s <- which(year > yr | (year == yr & ba_ev))[1L]; return(if (is.na(s)) Inf else s) }
  if (level == 1L && first == 1L) return(Inf)        # T3d: undated partial first seen at the first visit
  first
}

build_events <- function(rec, year, ba_ev) {
  e_type <- character(0); e_level <- integer(0); e_yr <- integer(0); e_first <- integer(0)
  nr <- length(rec$v); r_event <- integer(nr); r_action <- character(nr)
  date_of <- function(e) ifelse(is.na(e_yr[e]), year[e_first[e]], e_yr[e])
  same_window <- function(a, b) { lo <- min(a, b); hi <- max(a, b); !any(year > lo & year < hi) }
  current <- function(T) {
    ids <- which(e_type == T); if (!length(ids)) return(NA_integer_)
    d <- date_of(ids); ids <- ids[d == max(d)]; max(ids)
  }
  st <- function(e) start_visit(e_yr[e], e_level[e], e_first[e], year, ba_ev)
  ended_by_now <- function(e, v) {
    s <- which(e_level == 2L & e_first <= v & seq_along(e_type) != e)
    if (!length(s)) return(FALSE)
    se <- st(e); ss <- vapply(s, st, numeric(1))
    any(ss > se & ss <= v)
  }
  new_event <- function(T, L, Y, v) {
    e_type <<- c(e_type, T); e_level <<- c(e_level, L); e_yr <<- c(e_yr, Y); e_first <<- c(e_first, v)
    length(e_type)
  }
  for (r in seq_len(nr)) {
    v <- rec$v[r]; T <- rec$type[r]; L <- rec$level[r]; Y <- rec$yr[r]
    e <- current(T)
    if (!is.na(Y)) {                                              # ---- dated record
      if (is.na(e)) { e <- new_event(T, L, Y, v); act <- "new_dated" }
      else if (!is.na(e_yr[e])) {
        E <- e_yr[e]
        if (same_window(E, Y)) {
          if (Y > E) { e_yr[e] <- Y; e_level[e] <- L; act <- "merge_later_year" }
          else if (Y == E) { e_level[e] <- max(e_level[e], L); act <- "merge_same_year" }
          else act <- "merge_earlier_year_kept"
        } else { e <- new_event(T, L, Y, v); act <- "new_dated_other_window" }
      } else {
        if (Y <= year[e_first[e]]) { e_yr[e] <- Y; e_level[e] <- L; act <- "date_adopted" }
        else { e <- new_event(T, L, Y, v); act <- "new_dated_after_undated" }
      }
    } else {                                                      # ---- undated record
      if (is.na(e)) { e <- new_event(T, L, NA_integer_, v); act <- "new_undated" }
      else if (e_first[e] == v) { e_level[e] <- max(e_level[e], L); act <- "same_visit" }
      else if (ended_by_now(e, v)) { e <- new_event(T, L, NA_integer_, v); act <- "new_after_reset" }
      else if (isTRUE(ba_ev[v])) { e <- new_event(T, L, NA_integer_, v); act <- "new_ba_event" }
      else act <- "continuation"
    }
    r_event[r] <- e; r_action[r] <- act
  }
  list(type = e_type, level = e_level, yr = e_yr, first = e_first, r_event = r_event, r_action = r_action)
}

# E4. Logical matrix events x visits.
event_activity <- function(ev, year, ba_ev) {
  nE <- length(ev$type); nV <- length(year)
  if (nE == 0L) return(matrix(FALSE, 0L, nV))
  vi <- seq_len(nV)
  sv <- vapply(seq_len(nE), function(e) start_visit(ev$yr[e], ev$level[e], ev$first[e], year, ba_ev), numeric(1))
  so <- which(ev$level == 2L & is.finite(sv))
  for (s0 in unique(sv[so])) {                                   # stand-origin events of different types, same start visit
    g <- so[sv[so] == s0]
    if (length(unique(ev$type[g])) > 1L) { keep <- choose_stand_origin(g, ev$type, ev$yr); sv[setdiff(g, keep)] <- Inf }
  }
  act <- matrix(FALSE, nE, nV)
  so <- which(ev$level == 2L & is.finite(sv))
  for (e in seq_len(nE)) {
    act[e, ] <- vi >= sv[e]
    later <- so[so != e & sv[so] > sv[e] & is.finite(sv[so])]
    if (length(later)) act[e, ] <- act[e, ] & vi < min(sv[later])
  }
  act
}

# E5. Per visit values for one plot.
visit_values <- function(ev, act, year) {
  nV <- length(year); out <- list()
  date_e <- ifelse(is.na(ev$yr), year[ev$first], ev$yr)
  for (T in DIST_TYPES) {
    flag <- integer(nV); real <- mn <- mx <- md <- rep(NA_real_, nV); dating <- rep(NA_character_, nV); eid <- rep(NA_integer_, nV)
    ids <- which(ev$type == T)
    for (v in seq_len(nV)) {
      a <- ids[act[ids, v]]
      if (!length(a)) next
      a <- a[order(-date_e[a], -ev$level[a], -a)][1L]
      eid[v] <- a; flag[v] <- if (T == "harv") ev$level[a] else 1L
      if (!is.na(ev$yr[a])) {
        real[v] <- mn[v] <- mx[v] <- md[v] <- year[v] - ev$yr[a]; dating[v] <- "dated"
      } else if (ev$first[a] > 1L) {
        lo <- year[ev$first[a]]; hi <- year[ev$first[a] - 1L]
        mn[v] <- year[v] - lo; mx[v] <- year[v] - hi; md[v] <- year[v] - (lo + hi) / 2; dating[v] <- "bounded"
      } else dating[v] <- "first_visit"
    }
    out[[T]] <- list(flag = flag, real = real, min = mn, max = mx, mid = md, dating = dating, event = eid)
  }
  out
}

# Wrapper for one plot (rows in visit order). Returns per-visit columns + the plot's events and record actions.
disturbance_plot <- function(year, total_ba, origine, an_origine, perturb, an_perturb, codes, ba_regen, loss_thresh) {
  stopifnot(!is.unsorted(year))
  ba_ev <- ba_event_flag(total_ba, ba_regen, loss_thresh)
  rec <- make_records(year, origine, an_origine, perturb, an_perturb, codes)
  ev <- build_events(rec, year, ba_ev)
  act <- event_activity(ev, year, ba_ev)
  vv <- visit_values(ev, act, year)
  list(visits = list(ba_event = ba_ev,
         fire = vv$fire$flag, harv = vv$harv$flag, pest = vv$pest$flag,
         tsf = vv$fire$real, tsh = vv$harv$real, sfp = vv$pest$real,
         tsf_imputed = vv$fire$min, tsh_imputed = vv$harv$min, sfp_imputed = vv$pest$min,
         tsf_imputed_prev = vv$fire$max, tsh_imputed_prev = vv$harv$max, sfp_imputed_prev = vv$pest$max,
         tsf_imputed_mid = vv$fire$mid, tsh_imputed_mid = vv$harv$mid, sfp_imputed_mid = vv$pest$mid,
         fire_dating = vv$fire$dating, harv_dating = vv$harv$dating, pest_dating = vv$pest$dating,
         fire_event = vv$fire$event, harv_event = vv$harv$event, pest_event = vv$pest$event),
       events = ev, records = rec)
}
