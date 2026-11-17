# ==============================================================================
# 0_01_record_contents.R -- Chapter 2 analysis, section 0 (Data, no model), row 1: plots, intervals, events per transition.
# Run from the project root:  Rscript R/analysis/0_01_record_contents.R
#
# Event = a visit interval whose end state differs from its start state (interval-censored: only start and end state are seen, so an
# A->B->A or A->B->C path inside one interval is not counted as such; events are lower bounds on real changes).
# Choices (Lukas, 2026-11-17): unit = all 156 origin-destination pairs with a "modelled" flag (>= MIN_EVENTS); counts only (no at-risk
# denominator; time at risk and crude rates belong to row 4); record description included (interval length, visits per plot, years);
# exclusion flow included; heatmap too.
# Inputs (read-only): data/tables/<RUN_TAG>/{visit_intervals, plot_info, transition_set, exclusion_flow}
# Outputs: data-output/<RUN_TAG>/analysis/0_data/
#   0_01_exclusion_flow.csv          copy of 04x's table, after a check of its final row against visit_intervals.rds
#   0_01_events_by_transition.csv    from, to, n_events, n_plots, modelled (156 pairs, zeros included)
#   0_01_record_summary.csv          plots, visits, intervals, events, years, visits per plot
#   0_01_visits_per_plot.csv         number of plots by number of visits
#   0_01_interval_length.csv         interval length (years) overall and by start state
#   0_01_events_heatmap.(png|pdf|csv)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, RUN_OUT_DIR, MIN_EVENTS
OUT <- file.path(RUN_OUT_DIR, "analysis", "0_data"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 200)

# state order for display (same order as the old figures, R/13_fig_helpers.R)
ORDER <- c("regen", "jack_pine", "white_spruce", "black_spruce", "fir", "white_birch", "trembling_aspen", "other_pioneer",
           "yellow_birch", "red_maple", "sugar_maple", "cedar", "other_temperate")

iv <- as.data.table(readRDS(file.path(TABLE_DIR, "visit_intervals.rds")))
pi <- as.data.table(readRDS(file.path(TABLE_DIR, "plot_info.rds")))
ts <- fread(file.path(TABLE_DIR, "transition_set.csv"))
ef <- fread(file.path(TABLE_DIR, "exclusion_flow.csv"))
stopifnot("states outside the 13" = all(c(iv$state_start, iv$state_end) %in% ORDER))

# ---- 1. exclusion flow: check the final row against the interval table, then copy ------------------------------------------
fin <- ef[table == "flow" & grepl("^8 ", name)]
stopifnot("exclusion_flow final plots != visit_intervals"   = fin$plots == uniqueN(iv$id_pe),
          "exclusion_flow final intervals != visit_intervals" = fin$intervals == nrow(iv),
          "exclusion_flow final changes != visit_intervals" = fin$state_changes == sum(iv$is_change),
          "plot_info plots != visit_intervals plots"        = nrow(pi) == uniqueN(iv$id_pe))
fwrite(ef, file.path(OUT, "0_01_exclusion_flow.csv"))

# ---- 2. events per transition (all 156 pairs) -----------------------------------------------------------------------------
ev <- iv[is_change == TRUE, .(n_events = .N, n_plots = uniqueN(id_pe)), by = .(from = state_start, to = state_end)]
tr <- merge(CJ(from = ORDER, to = ORDER)[from != to], ev, by = c("from", "to"), all.x = TRUE)
for (v in c("n_events", "n_plots")) set(tr, which(is.na(tr[[v]])), v, 0L)
tr[, modelled := n_events >= MIN_EVENTS]
chk <- merge(tr, ts[, .(from, to, ts_events = n_events, ts_plots = n_plots, ts_kept = kept)], by = c("from", "to"))
stopifnot("pairs: not 156" = nrow(tr) == 156L,
          "events differ from transition_set.csv"   = chk[, all(n_events == ts_events)],
          "plots differ from transition_set.csv"    = chk[, all(n_plots == ts_plots)],
          "modelled flag differs from transition_set.csv" = chk[, all(modelled == ts_kept)],
          "sum of events != state changes" = tr[, sum(n_events)] == sum(iv$is_change))
setorder(tr, -n_events)
fwrite(tr, file.path(OUT, "0_01_events_by_transition.csv"))

# ---- 3. record description ----------------------------------------------------------------------------------------------
iv[, len_y := year_end - year_start]
stopifnot("n_visits in plot_info != intervals + 1" = iv[, .N + 1L, by = id_pe][pi, on = "id_pe"][, all(V1 == n_visits)])
last_year <- iv[, .(last_year = max(year_end)), by = id_pe]
pi <- merge(pi, last_year, by = "id_pe")
q3 <- function(x) as.list(setNames(as.numeric(quantile(x, c(0, 0.25, 0.5, 0.75, 1))), c("min", "q25", "median", "q75", "max")))
rs <- rbind(
  data.table(quantity = c("plots", "visits", "intervals", "events (state changes)", "plots with >= 1 event",
                          "modelled transitions", "events in modelled transitions", "events in unmodelled pairs"),
             value = c(uniqueN(iv$id_pe), pi[, sum(n_visits)], nrow(iv), sum(iv$is_change), uniqueN(iv[is_change == TRUE, id_pe]),
                       tr[modelled == TRUE, .N], tr[modelled == TRUE, sum(n_events)], tr[modelled == FALSE, sum(n_events)]), stat = "count"),
  melt(data.table(`first visit year` = pi$first_year, `last visit year` = pi$last_year, `follow-up per plot (y)` = pi$last_year - pi$first_year,
                  `visits per plot` = pi$n_visits)[, lapply(.SD, function(x) unlist(q3(x))[c(1, 3, 5)])][, stat := c("min", "median", "max")],
       id.vars = "stat", variable.name = "quantity")[, .(quantity = as.character(quantity), value, stat)])
fwrite(rs, file.path(OUT, "0_01_record_summary.csv"))
fwrite(pi[, .(n_plots = .N), by = .(n_visits)][order(n_visits)], file.path(OUT, "0_01_visits_per_plot.csv"))
il <- rbind(iv[, c(list(start_state = "all", n_intervals = .N, mean = mean(len_y)), q3(len_y))],
            iv[, c(list(n_intervals = .N, mean = mean(len_y)), q3(len_y)), by = .(start_state = state_start)][match(ORDER, start_state)])
fwrite(il, file.path(OUT, "0_01_interval_length.csv"))

# ---- 4. heatmap: origin x destination, event counts (log scale), outline = modelled ---------------------------------------
g <- copy(tr)
g[, `:=`(from = factor(from, levels = rev(ORDER)), to = factor(to, levels = ORDER), fillv = log10(pmax(n_events, 1)), seen = n_events > 0)]
p <- ggplot(g, aes(to, from)) +
  geom_tile(data = g[seen == FALSE], fill = "grey92", colour = "white") +
  geom_tile(data = g[seen == TRUE], aes(fill = fillv), colour = "white") +
  geom_tile(data = g[modelled == TRUE], fill = NA, colour = "black", linewidth = 0.4) +
  geom_text(data = g[seen == TRUE], aes(label = n_events, colour = fillv > 0.65 * max(g$fillv)), size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_viridis_c(option = "magma", direction = -1, breaks = log10(c(1, 10, 100, 500)), labels = c(1, 10, 100, 500)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x), position = "top", drop = FALSE) + scale_y_discrete(labels = function(x) gsub("_", " ", x), drop = FALSE) +
  labs(x = "destination", y = "origin", fill = "events\n(log scale)",
       subtitle = sprintf("Black outline: modelled (>= %d events). Grey: none observed.", MIN_EVENTS)) +
  coord_equal() + theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank())
ggsave(file.path(OUT, "0_01_events_heatmap.png"), p, width = 7.5, height = 7, dpi = 300)
ggsave(file.path(OUT, "0_01_events_heatmap.pdf"), p, width = 7.5, height = 7)
fwrite(tr, file.path(OUT, "0_01_events_heatmap.csv"))

# ---- console summary ---------------------------------------------------------------------------------------------------
cat("Exclusion flow (flow rows):\n"); print(ef[table == "flow", .(name, plots, visits, intervals, state_changes)], row.names = FALSE)
cat("\nRecord summary:\n"); print(rs, row.names = FALSE)
cat("\nInterval length (years):\n"); print(il, row.names = FALSE)
cat("\nVisits per plot:\n"); print(pi[, .(n_plots = .N), by = n_visits][order(n_visits)], row.names = FALSE)
cat("\nTop 10 transitions:\n"); print(tr[1:10], row.names = FALSE)
cat("\nChecks passed: exclusion-flow final row, plot_info, 156 pairs, events/plots/modelled flag vs transition_set.csv.\nSaved in", OUT, "\n")
