# ==============================================================================
# 04b_build_rows.R   (B2: row builder; B7: final event recount)
#
# Function definitions only (source() this file). Turns the visit-interval table
# into the split-row dataset that inla.surv() needs, one sequential pass per plot.
# This is the real-data counterpart of build_p1p2_combined_dataset() in the
# simulation folder (which is NOT touched).
#
# ROW RULES (time since entry E, in decades = years / TIME_SCALE; verified in
# docs/PIPELINE_NOTES.md and R/test_split_vs_collapsed.R):
#   * one SEGMENT per visit interval; a plot walks its intervals in order
#   * obs segment   (no state change): Tstop = (end visit - E),  truncation = (start visit - E)
#   * event segment (state change)   : exit T imputed inside the bracket
#                     (start visit, end visit) by p2_fun; Tstop = T - E,
#                     truncation = (start visit - E); then E <- T (clock reset)
#   * stub segment  (emitted right after EVERY event, also when the event interval is the
#                     plot's last one, so the time from the imputed exit to the last visit counts):
#                     from the imputed entry E to the first visit in the new state;
#                     truncation = 0, Tstop = (that visit - E), status 0
#   * truncation of the FIRST segment of the plot = max(W_j - E1, 0)/10 falls out of
#     the same formula because the first visit sits at t = W_j
#   * every segment has rows for ALL destinations (every state != from);
#     status 1 only on the realized destination of an event segment
#   * covariates of a segment = those of its interval (end visit); the stub takes the
#     covariates of the event interval that precedes it (its end visit IS the first
#     visit in the new state)
#   * transitions that are not kept still get their exit imputed (it sets the next
#     entry and the competing censoring times); their rows carry kept = FALSE and are
#     removed only at fitting time.
# ==============================================================================

library(data.table)

TIME_SCALE <- 10
EPSILON    <- 1e-6

# Stage-1 style exit imputation (also the fallback): midpoint of the bracket.
# Signature is the hook for the later moment-matched / joint-draw P2:
#   p2_fun(a_dur, b_dur, from, to, r)  -> T_dur inside [a_dur, b_dur]; r = row index in iv (event interval)
p2_midpoint <- function(a_dur, b_dur, from, to, r) (a_dur + b_dur) / 2

# ------------------------------------------------------------------------------
# build_rows()
#   iv           visit-interval table (04_prepare_inla.R), sorted by id_pe, interval
#   first_entry  numeric, ABSOLUTE entry time of each plot's first sojourn
#                (same clock as t_start/t_end), NAMED by id_pe. Stage 1: W_j / 2.
#   states       all states (destinations = every state != from)
#   tr           transition_set (adds the `kept` flag)
#   keep_cols    interval-table columns carried onto the rows
# Returns list(rows = split rows, segs = one line per segment, n_plots)
# ------------------------------------------------------------------------------
build_rows <- function(iv, first_entry, states, tr, p2_fun = p2_midpoint,
                       keep_cols = c("fire", "harv", "pest", "year_end")) {
  stopifnot(is.data.table(iv))
  keep_cols <- unique(c(keep_cols, intersect(names(iv), COV_CONT)))
  ids <- iv[, .(i1 = .I[1L], n = .N), by = id_pe]
  stopifnot(nrow(ids) == uniqueN(iv$id_pe),                       # each plot is one contiguous block
            iv[, !is.unsorted(interval), by = id_pe][, all(V1)])  # intervals in order within plot
  if (anyNA(first_entry[as.character(ids$id_pe)])) stop("first_entry is NA for some plots: P1 fallback required, plots must not be dropped.")

  a <- iv$t_start; b <- iv$t_end; s0 <- iv$state_start; s1 <- iv$state_end; chg <- iv$is_change
  cap <- nrow(iv) + sum(chg)
  o_row <- integer(cap); o_type <- character(cap); o_from <- character(cap); o_to <- character(cap)
  o_E <- numeric(cap); o_start <- numeric(cap); o_end <- numeric(cap); o_soj <- integer(cap)
  k <- 0L

  for (p in seq_len(nrow(ids))) {
    i1 <- ids$i1[p]; n <- ids$n[p]
    E <- unname(first_entry[as.character(ids$id_pe[p])])
    if (E < 0 || E > a[i1] + 1e-9) stop(sprintf("plot %s: entry %.3f outside [0, W_j=%.3f]", ids$id_pe[p], E, a[i1]))
    soj <- 1L
    for (r in i1:(i1 + n - 1L)) {
      k <- k + 1L
      o_row[k] <- r; o_from[k] <- s0[r]; o_E[k] <- E; o_start[k] <- a[r]; o_soj[k] <- soj
      if (!chg[r]) {
        o_type[k] <- "obs"; o_to[k] <- NA_character_; o_end[k] <- b[r]
      } else {
        a_dur <- max(a[r] - E, 0); b_dur <- b[r] - E
        Td <- p2_fun(a_dur, b_dur, s0[r], s1[r], r)
        if (is.na(Td) || Td < a_dur - 1e-9 || Td > b_dur + 1e-9)
          stop(sprintf("plot %s interval %d: imputed exit %.3f outside bracket (%.3f, %.3f)", ids$id_pe[p], iv$interval[r], Td, a_dur, b_dur))
        o_type[k] <- "event"; o_to[k] <- s1[r]; o_end[k] <- E + Td
        E <- E + Td; soj <- soj + 1L
        k <- k + 1L                                      # stub: new entry -> first visit in the new state
        o_row[k] <- r; o_type[k] <- "stub"; o_from[k] <- s1[r]; o_to[k] <- NA_character_
        o_E[k] <- E; o_start[k] <- E; o_end[k] <- b[r]; o_soj[k] <- soj
      }
    }
  }
  idx <- seq_len(k)
  segs <- data.table(seg_id = idx, id_pe = iv$id_pe[o_row[idx]], interval = iv$interval[o_row[idx]],
                     iv_row = o_row[idx], type = o_type[idx], from = o_from[idx], event_to = o_to[idx],
                     sojourn = o_soj[idx], E_abs = o_E[idx], start_abs = o_start[idx], end_abs = o_end[idx])
  segs[, truncation := pmax((start_abs - E_abs) / TIME_SCALE, 0)]
  segs[, floored := (end_abs - E_abs) / TIME_SCALE < truncation + EPSILON]   # exit imputed AT the last visit: inla.surv needs Tstop > truncation
  segs[, Tstop := pmax((end_abs - E_abs) / TIME_SCALE, truncation + EPSILON)]

  pairs <- CJ(from = states, to = states)[from != to]
  rows <- merge(segs, pairs, by = "from", allow.cartesian = TRUE)
  rows[, status := as.integer(!is.na(event_to) & to == event_to)]
  rows <- cbind(rows, iv[rows$iv_row, ..keep_cols])
  rows[tr, kept := i.kept, on = c("from", "to")]
  setorder(rows, id_pe, seg_id, to)
  list(rows = rows, segs = segs, n_plots = nrow(ids), n_floored = sum(segs$floored))
}

# ------------------------------------------------------------------------------
# recount_events()   (B7: the very last check before INLA)
# Recounts events per kept transition and per disturbance level on the FINAL rows
# and compares them with the transition set. Returns the table; stops if any kept
# transition no longer passes the thresholds or the counts disagree.
# ------------------------------------------------------------------------------
recount_events <- function(rows, tr, min_events = MIN_EVENTS) {   # MIN_EVENTS from R/config.R (source it first)
  ev <- rows[status == 1L & kept == TRUE]
  cnt <- ev[, .(n_events = .N, n_fire = sum(fire == 1L), n_pest = sum(pest == 1L),
                n_harv0 = sum(harv == 0L), n_harv1 = sum(harv == 1L), n_harv2 = sum(harv == 2L)),
            by = .(from, to)]
  chk <- merge(tr[kept == TRUE, .(from, to, fire_term, pest_term, harv_term,
                                  e0 = n_events, f0 = n_fire, p0 = n_pest, h00 = n_harv0, h10 = n_harv1, h20 = n_harv2)],
               cnt, by = c("from", "to"), all = TRUE)
  chk[, ok_counts := n_events == e0 & n_fire == f0 & n_pest == p0 & n_harv0 == h00 & n_harv1 == h10 & n_harv2 == h20]
  chk[, ok_thresholds := n_events >= min_events &
        (!fire_term | n_fire >= min_events) & (!pest_term | n_pest >= min_events) &
        (harv_term != "3level" | (n_harv1 + n_harv2 >= min_events))]    # A2: light + heavy together
  if (anyNA(chk$ok_counts) || !all(chk$ok_counts)) stop("recount_events: event counts differ from transition_set")
  if (!all(chk$ok_thresholds)) stop("recount_events: a kept transition/term no longer passes the >= MIN_EVENTS rule")
  chk[]
}
