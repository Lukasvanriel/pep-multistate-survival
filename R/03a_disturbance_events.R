# ==============================================================================
# R/03a_disturbance_events.R -- disturbance event logic (sourced by 03, not run).
# Version 2, 2026-10-26 afternoon; A1/A2 revised 2026-10-26 evening (user: "same type, same interval = one event,
# the highest intensity with its year wins; equal intensity -> the higher year"; leftover rule for undated records).
# Rules as decided by the user after the first rewrite broke the back-fill
# (docs/PIPELINE_NOTES.md, "2026-10-26 (afternoon)"). Every rule below was decided by the user; implementation
# details that are not user decisions are marked [impl] and listed in PIPELINE_NOTES.
#
# Per plot, visits in time order. FIRST all years are assigned from the data (steps 1-3), THEN events, order and
# resets (steps 4-7).
#  0 records   one record per code per visit: type (fire/harv/pest), intensity (2 = stand-origin, 1 = partial), year.
#              CE gives two records (partial harvest + partial pest). A year after the visit is kept (event after the
#              measurement).
#  0b station merge (2026-10-26, user; only when station codes are passed): station_pe codes become records the same way
#              (no year). Per visit: a station record is DROPPED if pee has a record of the same type (pee wins, whatever the
#              intensities), or if it is stand-origin and pee has a stand-origin record of another type (one stand origin:
#              pee wins). Otherwise it is ADDED as an undated record (different types are kept together, e.g. salvage).
#              Ignore / drop-plot / CHT-DT codes of either layer are applied in 03 before this (exclusions first).
#  1 BA event  per visit: previous BA >= BA_REGEN and (loss >= BA_LOSS_THRESH or BA now < BA_REGEN); no BA event when the
#              previous BA was already < BA_REGEN (user, 2026-10-26 evening); first visit: BA < BA_REGEN.
#  "Same window" / "same interval" of two years: no visit year strictly between them (all years before the first visit
#  are in one interval).
#  2 A1        same type at the same visit: dated records whose years are in the same interval -> one record with the
#              highest intensity and that record's year; equal intensity -> the higher year. Dated records in
#              different intervals stay separate records (separate events). Undated records -> one record (highest
#              intensity). An undated record next to a dated one of its type at the same visit is a "pair", see 4.
#  3 back-fill an undated record (pairs included) takes the year of the FIRST LATER dated record of its type whose year <=
#              its visit year (incompatible dated records in between are skipped); it keeps its own intensity (it may
#              then be raised by A2). Then A1 again: a back-filled pair record merges with its partner when both
#              years are in one interval [impl, follows "first use the data", 2026-10-26].
#  4 events    pair (undated + dated, same type, same visit): if an event of that type AND the undated record's
#              intensity is already tracked -> the undated record is a leftover of it (continuation, regardless of a
#              BA event) and the dated record is handled on its own; otherwise (incl. the first visit) both are one
#              record: the highest intensity with the dated year (CT - + CP 1987 -> heavy 1987).
#              dated record: same type and same interval as an existing dated event -> A2: the same event, which takes
#              the highest intensity with its year; equal intensity -> the higher year. Else a new event.
#              undated record: if an event of its type is tracked (not ended by a stand-origin event that started at
#              an earlier visit) -> the same event (old year, old intensity), UNLESS there is a BA event at this visit
#              and no stand-origin record of another type at this visit -> new event. No tracked event -> new event.
#  5 dating    dated -> its year. Undated, first seen at visit F > 1 -> between visit F-1 and F: lower (TSD from F),
#              upper (TSD from F-1), mid. Undated at the FIRST visit: stand-origin or partial fire -> stand age W1 at the
#              first visit; partial harvest / pest -> fv_age[type] (median age of dated first-visit partial codes;
#              S19: 10 / 90% quantiles). TSD = that age + years since the first visit (the same in lower/mid/upper).
#  6 timeline  start: dated -> first visit with year > event year, or == event year with a BA event (same-year rule);
#              undated -> first-seen visit. Order of two events: "known older" if the latest possible date of one is
#              before the earliest possible date of the other (dated: the year; bounded: (year[F-1], year[F]];
#              first-visit undated: (-Inf, year[1]]); otherwise unknown. A stand-origin event ends every event (all
#              types) known to be older, from its start on; unknown order -> both stay.
#              Selection: stand-origin events of different types with unknown order (incl. the same year) -> only
#              one kept: all dated -> the latest year (same year: harvest > fire > pest); any undated -> harvest if
#              present, else the dated one, else fire. The others are dropped (never active, end nothing).
#  7 visit     per visit and type, among the active events: the most recent (not known older than another one);
#              unknown order -> the more severe [user]; then [impl] the later latest-possible date, then the later
#              created. Values: flag (fire/pest 0/1, harv 0/1/2), TSD real / lower / upper / mid, dating marker.
# ==============================================================================

DIST_TYPES <- c("fire", "harv", "pest")
same_window <- function(a, b, year) { lo <- min(a, b); hi <- max(a, b); !any(year > lo & year < hi) }
SO_PRIORITY <- c(harv = 3L, fire = 2L, pest = 1L)

ba_event_flag <- function(total_ba, ba_regen, loss_thresh) {
  n <- length(total_ba)
  if (n == 0L) return(logical(0))
  prev <- c(NA_real_, total_ba[-n])
  # 2026-10-26 (user): no BA event when the previous BA was already < ba_regen (applies to all uses: rule B, same-year rule,
  # CHT/DT truncation). Before: loss >= loss_thresh from any previous BA > 0.
  out <- !is.na(prev) & prev >= ba_regen & (total_ba <= (1 - loss_thresh) * prev | total_ba < ba_regen)
  out[1L] <- total_ba[1L] < ba_regen
  out
}

# step 0: raw records (before A1)
make_records <- function(year, origine, an_origine, perturb, an_perturb, codes) {
  v <- integer(0); typ <- character(0); lev <- integer(0); yr <- integer(0)
  for (i in seq_along(year)) {
    if (!is.na(origine[i])) for (T in DIST_TYPES) if (origine[i] %in% codes[[paste0(T, "_o")]]) {
      v <- c(v, i); typ <- c(typ, T); lev <- c(lev, 2L); yr <- c(yr, as.integer(an_origine[i])) }
    if (!is.na(perturb[i])) for (T in DIST_TYPES) if (perturb[i] %in% codes[[paste0(T, "_p")]]) {
      v <- c(v, i); typ <- c(typ, T); lev <- c(lev, 1L); yr <- c(yr, as.integer(an_perturb[i])) }
  }
  data.frame(v = v, type = typ, level = lev, yr = yr, stringsAsFactors = FALSE)
}

# step 0b: station records merged into the pee records (rules in the header)
merge_station_records <- function(rp, rs) {
  if (!nrow(rs)) return(list(r = rp, n_added = 0L, n_same_type = 0L, n_other_so = 0L))
  same <- paste(rs$v, rs$type) %in% paste(rp$v, rp$type)
  so_v <- unique(rp$v[rp$level == 2L])
  other_so <- !same & rs$level == 2L & rs$v %in% so_v
  keep <- !same & !other_so
  r <- rbind(rp, rs[keep, , drop = FALSE]); r <- r[order(r$v), , drop = FALSE]
  list(r = r, n_added = sum(keep), n_same_type = sum(same), n_other_so = sum(other_so))
}

# step 2: A1 -- per (visit, type): dated records merged per interval, undated records merged, pairs marked
collapse_same_visit <- function(r, year) {
  if (is.null(r$n_raw)) r$n_raw <- rep(1L, nrow(r))
  if (is.null(r$yr_src)) r$yr_src <- ifelse(is.na(r$yr), "none", "data")
  if (is.null(r$backfill_from)) r$backfill_from <- rep(NA_integer_, nrow(r))
  if (!nrow(r)) return(cbind(r, pair = logical(0)))
  key <- paste(r$v, r$type)
  out <- lapply(split(seq_len(nrow(r)), factor(key, levels = unique(key))), function(ix) {
    x <- r[ix, , drop = FALSE]; res <- list()
    d <- x[!is.na(x$yr), , drop = FALSE]; u <- x[is.na(x$yr), , drop = FALSE]
    if (nrow(d)) {
      comp <- rep(NA_integer_, nrow(d)); k <- 0L
      for (a in seq_len(nrow(d))) if (is.na(comp[a])) { k <- k + 1L; st <- a
        while (length(st)) { z <- st[1]; st <- st[-1]; if (is.na(comp[z])) { comp[z] <- k
          st <- c(st, which(is.na(comp) & vapply(d$yr, function(y) same_window(y, d$yr[z], year), logical(1)))) } } }
      for (g in unique(comp)) { y <- d[comp == g, , drop = FALSE]; L <- max(y$level)
        w <- which(y$level == L)[which.max(y$yr[y$level == L])]
        res[[length(res) + 1L]] <- data.frame(v = x$v[1], type = x$type[1], level = L, yr = as.integer(y$yr[w]),
                                              n_raw = sum(y$n_raw), yr_src = y$yr_src[w], backfill_from = y$backfill_from[w],
                                              pair = FALSE, stringsAsFactors = FALSE) }
    }
    if (nrow(u)) res[[length(res) + 1L]] <- data.frame(v = x$v[1], type = x$type[1], level = max(u$level), yr = NA_integer_,
                                                      n_raw = sum(u$n_raw), yr_src = "none", backfill_from = NA_integer_,
                                                      pair = nrow(d) > 0, stringsAsFactors = FALSE)
    do.call(rbind, res)
  })
  out <- do.call(rbind, out)
  # within a visit: undated (pair) records first, so the pair decision is taken before the dated partner is handled
  out[order(out$v, match(out$type, DIST_TYPES), !is.na(out$yr), out$yr), , drop = FALSE]
}

# step 3: back-fill (all undated records, pairs included: the data are used first); afterwards A1 is applied again
# (collapse_same_visit), so a back-filled pair record merges with its partner when both years are in one interval
backfill_years <- function(r, year) {
  if (!nrow(r)) return(r)
  for (i in which(is.na(r$yr))) {
    j <- which(r$type == r$type[i] & r$v > r$v[i] & r$yr_src == "data" & r$yr <= year[r$v[i]])
    if (length(j)) { j <- j[which.min(r$v[j])]; r$yr[i] <- r$yr[j]; r$yr_src[i] <- "backfill"; r$backfill_from[i] <- r$v[j] }
  }
  r
}

start_visit_of <- function(yr, first, year, ba_ev) {
  if (!is.na(yr)) { s <- which(year > yr | (year == yr & ba_ev))[1L]; return(if (is.na(s)) Inf else s) }
  first
}
# possible date range of an event: dated -> [yr, yr]; undated first seen F > 1 -> (year[F-1], year[F]]; F == 1 -> (-Inf, year[1]]
date_range <- function(yr, first, year) {
  if (!is.na(yr)) return(c(yr, yr))
  if (first > 1L) return(c(year[first - 1L] + 0.5, year[first]))
  c(-Inf, year[1L])
}
known_older <- function(a, b) a[2] < b[1]          # a, b = date ranges; TRUE if a is certainly before b

# step 4: events
build_events <- function(r, year, ba_ev) {
  e_type <- character(0); e_level <- integer(0); e_yr <- integer(0); e_first <- integer(0)
  r$event <- rep(NA_integer_, nrow(r)); r$action <- rep(NA_character_, nrow(r))
  new_event <- function(T, L, Y, v) { e_type <<- c(e_type, T); e_level <<- c(e_level, L); e_yr <<- c(e_yr, Y); e_first <<- c(e_first, v); length(e_type) }
  so_codes_at <- function(v, T) any(r$v == v & r$type != T & r$level == 2L)
  ended_before <- function(e, v) {        # ended by a stand-origin event (other id) that started at a visit < v and is known later
    s <- which(e_level == 2L & seq_along(e_type) != e)
    if (!length(s)) return(FALSE)
    re <- date_range(e_yr[e], e_first[e], year)
    any(vapply(s, function(k) { sk <- start_visit_of(e_yr[k], e_first[k], year, ba_ev)
      sk < v && known_older(re, date_range(e_yr[k], e_first[k], year)) }, logical(1)))
  }
  merged_into <- rep(NA_integer_, nrow(r))
  for (i in seq_len(nrow(r))) {
    v <- r$v[i]; T <- r$type[i]; L <- r$level[i]; Y <- r$yr[i]
    if (is.na(Y) && r$pair[i]) {                                   # pair: leftover or merged with the dated partner
      tr <- which(e_type == T & e_level == L)
      tr <- tr[!vapply(tr, function(k) ended_before(k, v), logical(1))]
      if (length(tr)) { r$event[i] <- max(tr); r$action[i] <- "pair_leftover" }
      else {
        j <- which(r$v == v & r$type == T & !is.na(r$yr))
        j <- j[order(-r$level[j], -r$yr[j])][1L]                    # [impl] several dated partners: highest, then latest
        r$level[j] <- max(r$level[j], L); r$n_raw[j] <- r$n_raw[j] + r$n_raw[i]; merged_into[i] <- j; r$action[i] <- "pair_merged"
      }
      next
    }
    if (!is.na(Y)) {
      cand <- which(e_type == T & !is.na(e_yr))
      cand <- cand[vapply(cand, function(k) same_window(e_yr[k], Y, year), logical(1))]
      if (length(cand)) { k <- min(cand); r$event[i] <- k
        if (L > e_level[k]) { e_level[k] <- L; e_yr[k] <- Y; r$action[i] <- "A2_higher_intensity" }
        else if (L == e_level[k] && Y > e_yr[k]) { e_yr[k] <- Y; r$action[i] <- "A2_higher_year" }
        else r$action[i] <- "A2_kept"
      }
      else { r$event[i] <- new_event(T, L, Y, v); r$action[i] <- if (r$yr_src[i] == "backfill") "new_dated_backfilled" else "new_dated" }
    } else {
      tracked <- which(e_type == T)
      tracked <- tracked[!vapply(tracked, function(k) ended_before(k, v), logical(1))]
      if (!length(tracked)) { r$event[i] <- new_event(T, L, NA_integer_, v); r$action[i] <- "new_undated" }
      else if (isTRUE(ba_ev[v]) && !so_codes_at(v, T)) { r$event[i] <- new_event(T, L, NA_integer_, v); r$action[i] <- "new_ba_event" }
      else { r$event[i] <- max(tracked); r$action[i] <- if (isTRUE(ba_ev[v])) "continuation_ba_explained" else "continuation" }
    }
  }
  r$event[!is.na(merged_into)] <- r$event[merged_into[!is.na(merged_into)]]
  list(type = e_type, level = e_level, yr = e_yr, first = e_first, records = r)
}

choose_stand_origin <- function(ids, type, yr) {
  pr <- SO_PRIORITY[type[ids]]
  if (all(!is.na(yr[ids]))) return(ids[order(-yr[ids], -pr)][1L])
  if (any(type[ids] == "harv")) { h <- ids[type[ids] == "harv"]; return(h[order(-ifelse(is.na(yr[h]), -Inf, yr[h]))][1L]) }
  d <- ids[!is.na(yr[ids])]
  if (length(d)) return(d[order(-yr[d], -SO_PRIORITY[type[d]])][1L])
  ids[order(-pr)][1L]
}

# step 6: activity matrix (events x visits) + dropped flags
event_activity <- function(ev, year, ba_ev) {
  nE <- length(ev$type); nV <- length(year); vi <- seq_len(nV)
  if (nE == 0L) return(list(act = matrix(FALSE, 0L, nV), dropped = logical(0), start = numeric(0), range = list()))
  sv <- vapply(seq_len(nE), function(e) start_visit_of(ev$yr[e], ev$first[e], year, ba_ev), numeric(1))
  rg <- lapply(seq_len(nE), function(e) date_range(ev$yr[e], ev$first[e], year))
  dropped <- rep(FALSE, nE)
  so <- which(ev$level == 2L & is.finite(sv))
  if (length(so) > 1L) {                                   # selection among stand-origin events of different types, unknown order
    adj <- matrix(FALSE, length(so), length(so))
    for (a in seq_along(so)) for (b in seq_along(so)) if (a != b && ev$type[so[a]] != ev$type[so[b]] &&
        !known_older(rg[[so[a]]], rg[[so[b]]]) && !known_older(rg[[so[b]]], rg[[so[a]]])) adj[a, b] <- TRUE
    comp <- rep(NA_integer_, length(so)); k <- 0L
    for (a in seq_along(so)) if (is.na(comp[a])) { k <- k + 1L; st <- a
      while (length(st)) { x <- st[1]; st <- st[-1]; if (is.na(comp[x])) { comp[x] <- k; st <- c(st, which(adj[x, ] & is.na(comp))) } } }
    for (g in unique(comp)) { ids <- so[comp == g]
      if (length(unique(ev$type[ids])) > 1L) { keep <- choose_stand_origin(ids, ev$type, ev$yr); dropped[setdiff(ids, keep)] <- TRUE } }
  }
  act <- matrix(FALSE, nE, nV)
  live_so <- which(ev$level == 2L & is.finite(sv) & !dropped)
  for (e in which(!dropped)) {
    act[e, ] <- vi >= sv[e]
    later <- live_so[live_so != e & vapply(live_so, function(k) known_older(rg[[e]], rg[[k]]), logical(1))]
    if (length(later)) act[e, ] <- act[e, ] & vi < min(sv[later])
  }
  list(act = act, dropped = dropped, start = sv, range = rg)
}

# step 7: visit values
visit_values <- function(ev, ea, year, W1, fv_age) {
  nV <- length(year); out <- list()
  for (T in DIST_TYPES) {
    flag <- integer(nV); real <- mn <- mx <- md <- rep(NA_real_, nV); dating <- rep(NA_character_, nV); eid <- rep(NA_integer_, nV)
    ids <- which(ev$type == T)
    for (v in seq_len(nV)) {
      a <- ids[ea$act[ids, v]]
      if (!length(a)) next
      if (length(a) > 1L) {
        latest <- a[vapply(a, function(x) !any(vapply(setdiff(a, x), function(y) known_older(ea$range[[x]], ea$range[[y]]), logical(1))), logical(1))]
        a <- latest[order(-ev$level[latest], -vapply(latest, function(x) ea$range[[x]][2], numeric(1)), -latest)][1L]
      }
      eid[v] <- a; flag[v] <- if (T == "harv") ev$level[a] else 1L
      if (!is.na(ev$yr[a])) {
        real[v] <- mn[v] <- mx[v] <- md[v] <- year[v] - ev$yr[a]; dating[v] <- "dated"
      } else if (ev$first[a] > 1L) {
        lo <- year[ev$first[a]]; hi <- year[ev$first[a] - 1L]
        mn[v] <- year[v] - lo; mx[v] <- year[v] - hi; md[v] <- year[v] - (lo + hi) / 2; dating[v] <- "bounded"
      } else {
        age1 <- if (ev$level[a] == 2L || T == "fire") W1 else fv_age[[T]]
        mn[v] <- mx[v] <- md[v] <- age1 + (year[v] - year[1L])
        dating[v] <- if (ev$level[a] == 2L || T == "fire") "first_visit_standage" else "first_visit_age"
      }
    }
    out[[T]] <- list(flag = flag, real = real, min = mn, max = mx, mid = md, dating = dating, event = eid)
  }
  out
}

disturbance_plot <- function(year, total_ba, origine, an_origine, perturb, an_perturb, codes, ba_regen, loss_thresh,
                             W1 = NA_real_, fv_age = c(harv = NA_real_, pest = NA_real_),
                             s_origine = NULL, s_perturb = NULL) {
  stopifnot(!is.unsorted(year))
  ba_ev <- ba_event_flag(total_ba, ba_regen, loss_thresh)
  r0 <- make_records(year, origine, an_origine, perturb, an_perturb, codes)
  st <- list(n_added = 0L, n_same_type = 0L, n_other_so = 0L)
  if (!is.null(s_origine)) {
    nn <- rep(NA_integer_, length(year))
    st <- merge_station_records(r0, make_records(year, s_origine, nn, s_perturb, nn, codes)); r0 <- st$r
  }
  r <- collapse_same_visit(backfill_years(collapse_same_visit(r0, year), year), year)
  ev <- build_events(r, year, ba_ev)
  ea <- event_activity(ev, year, ba_ev)
  vv <- visit_values(ev, ea, year, W1, fv_age)
  list(visits = list(ba_event = ba_ev,
         fire = vv$fire$flag, harv = vv$harv$flag, pest = vv$pest$flag,
         tsf = vv$fire$real, tsh = vv$harv$real, sfp = vv$pest$real,
         tsf_imputed = vv$fire$min, tsh_imputed = vv$harv$min, sfp_imputed = vv$pest$min,
         tsf_imputed_prev = vv$fire$max, tsh_imputed_prev = vv$harv$max, sfp_imputed_prev = vv$pest$max,
         tsf_imputed_mid = vv$fire$mid, tsh_imputed_mid = vv$harv$mid, sfp_imputed_mid = vv$pest$mid,
         fire_dating = vv$fire$dating, harv_dating = vv$harv$dating, pest_dating = vv$pest$dating,
         fire_event = vv$fire$event, harv_event = vv$harv$event, pest_event = vv$pest$event),
       events = list(type = ev$type, level = ev$level, yr = ev$yr, first = ev$first, dropped = ea$dropped, start = ea$start),
       records = ev$records, station = st[c("n_added", "n_same_type", "n_other_so")])
}

# first-visit partial ages from the data: dated partial records at the first visit, age = year[1] - year
first_visit_partial_ages <- function(dt, codes, q) {
  f <- dt[, .SD[1L], by = id_pe]
  out <- c(harv = NA_real_, pest = NA_real_); n <- c(harv = 0L, pest = 0L)
  for (T in c("harv", "pest")) {
    x <- f[!is.na(perturb) & perturb %in% codes[[paste0(T, "_p")]] & !is.na(an_perturb) & an_perturb <= year, year - an_perturb]
    n[T] <- length(x); if (length(x)) out[T] <- as.numeric(quantile(x, q, names = FALSE))
  }
  attr(out, "n") <- n
  out
}
