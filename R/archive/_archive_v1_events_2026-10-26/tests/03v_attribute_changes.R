# ==============================================================================
# R/tests/03v_attribute_changes.R -- attribute every change between an OLD disturbance build and the CURRENT one
# (made by the event logic of 2026-10-26) to a rule. Diagnostic only; writes nothing.
# Run from the project root AFTER 03:  Rscript R/tests/03v_attribute_changes.R [old_build.rds]
# Default old build: data/build/disturbance_covariates_build_9298c98.rds (copy of the build used by noint/inttest).
# Needs data/build/disturbance_covariates.rds and data/build/disturbance_events.rds from the new 03.
# Compares per visit and type: flag (fire/pest 0/1, harv 0/1/2) and the TSD lower bound (*_imputed).
# Tags are assigned in a fixed priority order; "other: ..." = not explained by a rule (inspect with SHOW below).
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here("R", "config.R")); source(here("R", "03a_disturbance_events.R"))
args <- commandArgs(TRUE)
OLD <- if (length(args)) args[1] else here("data", "build", "disturbance_covariates_build_9298c98.rds")
SHOW <- 3L   # example plots printed per tag

nw  <- as.data.table(readRDS(here("data", "build", "disturbance_covariates.rds")))
old <- as.data.table(readRDS(OLD))
E   <- readRDS(here("data", "build", "disturbance_events.rds")); ev <- E$events; rc <- E$records; vv <- E$visits
setkey(nw, id_pe, no_mes); setkey(old, id_pe, no_mes); setkey(vv, id_pe, no_mes)
stopifnot("old and new builds have different visits" = identical(nw[, .(id_pe, no_mes)], old[, .(id_pe, no_mes)]),
          "events file does not match the build" = identical(nw[, .(id_pe, no_mes)], vv[, .(id_pe, no_mes)]))
nw <- cbind(nw, vv[, .(origine, an_origine, perturb, an_perturb, fire_event, harv_event, pest_event)])
nw[, vi := seq_len(.N), by = id_pe]

# start visit of every event (same function as the pipeline)
pl <- nw[, .(year = list(year), ba = list(ba_event), nm = list(no_mes)), by = id_pe]
ev <- merge(ev, pl, by = "id_pe")
ev[, first := mapply(function(nm, f) match(f, nm), nm, first_no_mes)]
ev[, start := mapply(function(y, l, f, yr, b) start_visit(yr, l, f, y, b), year, level, first, yr, ba)]
ev[, c("year", "ba", "nm") := NULL]
rc <- merge(rc, nw[, .(id_pe, no_mes, perturb)], by = c("id_pe", "no_mes"), all.x = TRUE)
ev <- merge(ev, unique(rc[level == 1L & perturb == "CE", .(id_pe, event)])[, ce := TRUE], by = c("id_pe", "event"), all.x = TRUE)
ev <- merge(ev, rc[, .(acts = paste(sort(unique(action)), collapse = ",")), by = .(id_pe, event)], by = c("id_pe", "event"), all.x = TRUE)
ev[is.na(ce), ce := FALSE]
so_v1 <- ev[level == 2L & start == 1, .(n_so_v1_types = uniqueN(type)), by = id_pe]

rows <- list()
for (T in c("fire", "harv", "pest")) {
  ts <- c(fire = "tsf", harv = "tsh", pest = "sfp")[[T]]
  x <- data.table(id_pe = nw$id_pe, no_mes = nw$no_mes, vi = nw$vi, year = nw$year, type = T, fo = old[[T]], fn = nw[[T]],
                  to = as.numeric(old[[paste0(ts, "_imputed")]]), tn = as.numeric(nw[[paste0(ts, "_imputed")]]),
                  dg = nw[[paste0(T, "_dating")]], e = nw[[paste0(T, "_event")]], ba_event = nw$ba_event)
  x <- x[fo != fn | !mapply(identical, to, tn)]
  EVT <- ev[type == T]
  x <- merge(x, EVT[, .(id_pe, e = event, e_act = acts, e_ce = ce, e_start = start, e_first = first, e_dated = !is.na(yr))],
             by = c("id_pe", "e"), all.x = TRUE)
  x[, n_started := mapply(function(p, v) sum(EVT[id_pe == p, start] <= v), id_pe, vi)]
  x[, has_t3d := id_pe %in% EVT[level == 1L & is.na(yr) & first == 1L, id_pe]]
  erased <- unique(merge(rc[type == T, .(id_pe, no_mes, event)],
                         data.table(id_pe = old$id_pe, no_mes = old$no_mes, of = old[[T]])[of == 0], by = c("id_pe", "no_mes"))[, .(id_pe, e = event)])[, erased := TRUE]
  x <- merge(x, erased, by = c("id_pe", "e"), all.x = TRUE); x[is.na(erased), erased := FALSE]
  x <- merge(x, so_v1, by = "id_pe", all.x = TRUE); x[is.na(n_so_v1_types), n_so_v1_types := 0L]
  x[, tag := fcase(
    fn == 0 & fo > 0 & has_t3d, "T3d: undated partial at the first visit, not flagged",
    fn == 0 & fo > 0 & n_started > 0, "reset: ended by a later stand-origin event (or dated after this visit)",
    fn == 0 & fo > 0, "other: flag lost",
    fo == 0 & fn > 0 & e_ce, "CE: second record (pest or harvest)",
    fo == 0 & fn > 0 & grepl("new_ba_event", e_act), "B: new event after a BA event",
    fo == 0 & fn > 0 & e_start == 1 & n_so_v1_types > 1, "compound: stand-origin events of 2 types before the first visit, both kept",
    fo == 0 & fn > 0 & erased, "restored: recorded code erased by the old back-fill",
    fo == 0 & fn > 0 & e_dated & e_start < e_first, "D/back-fill: dated event active before the visit that recorded it",
    fo == 0 & fn > 0, "other: flag gained",
    fo != fn & n_started > 1, "C: most recent of two events sets the level",
    fo != fn & grepl("new_ba_event", e_act), "B: new event after a BA event",
    fo != fn & grepl("merge_later_year|merge_same_year|date_adopted", e_act), "A/1: level from the latest-dated / dated record",
    fo != fn, "other: level changed",
    dg == "bounded", "F3: undated event, TSD now accumulates",
    dg == "first_visit", "first-visit stand-origin: TSD left to 04 (old: 0 or NA)",
    n_started > 1, "C: most recent of two events sets the TSD",
    grepl("new_ba_event", e_act), "B: new event after a BA event",
    grepl("merge_later_year|merge_earlier|merge_same", e_act), "A: same-window records merged (latest year)",
    grepl("new_dated_other_window|new_dated_after_undated|date_adopted", e_act), "dated: other window / dated record adopted",
    default = "other: TSD changed")]
  rows[[T]] <- x
}
a <- rbindlist(rows, fill = TRUE)
cat(sprintf("Old build: %s\nVisits compared: %d in %d plots. Changed visit x type cells: %d (%d visits, %d plots)\n\n",
            OLD, nrow(nw), uniqueN(nw$id_pe), nrow(a), uniqueN(a[, .(id_pe, no_mes)]), uniqueN(a$id_pe)))
tab <- dcast(a[, .N, by = .(tag, type)], tag ~ type, value.var = "N", fill = 0)
tab <- merge(tab, a[, .(plots = uniqueN(id_pe)), by = tag], by = "tag"); setorder(tab, -plots)
options(width = 200); print(tab, row.names = FALSE)

fmt <- function(f, t) ifelse(f == 0, "0", paste0(f, "|", ifelse(is.na(t), "NA", t)))
show_plot <- function(p) {
  x <- nw[id_pe == p]; o <- old[id_pe == p]
  cat(sprintf("\n--- %s: %s\n", p, paste(unique(a[id_pe == p, paste(type, sub(":.*", "", tag))]), collapse = "; ")))
  print(data.table(v = x$no_mes, year = x$year, BA = round(x$total_ba, 1), ba_ev = as.integer(x$ba_event),
    origine = ifelse(is.na(x$origine), "", paste(x$origine, ifelse(is.na(x$an_origine), "-", x$an_origine))),
    perturb = ifelse(is.na(x$perturb), "", paste(x$perturb, ifelse(is.na(x$an_perturb), "-", x$an_perturb))),
    `old F` = fmt(o$fire, o$tsf_imputed), `old H` = fmt(o$harv, o$tsh_imputed), `old P` = fmt(o$pest, o$sfp_imputed),
    `new F` = fmt(x$fire, x$tsf_imputed), `new H` = fmt(x$harv, x$tsh_imputed), `new P` = fmt(x$pest, x$sfp_imputed)), row.names = FALSE)
}
if (SHOW > 0) { set.seed(5); cat("\n==== Examples (", SHOW, " per tag; flag|TSD lower bound) ====\n", sep = "")
  for (t in tab$tag) { cat("\n#### ", t, "\n", sep = ""); for (p in head(sample(unique(a[tag == t, id_pe])), SHOW)) show_plot(p) } }
