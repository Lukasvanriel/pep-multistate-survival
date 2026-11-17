# Chapter 2 -- analysis log

One entry per plan row, appended in order. Written by the agent working on the row, checked by Lukas. Keep each entry very concise but complete; numbers in tables.

## Entry format

```
### <section> / <row quantity>   (<date>)
- Question: <the row's question>
- Old work: <script/output> -- correct / partly / wrong / absent; why
- Choices (Lukas): <grouping, climate values, real plots vs scenario stand, other>
- Done: <script in R/analysis/>, output <path>
- Result: <one to three sentences; key numbers in a table if needed>
- Caveat (only if it changes what can be claimed): <...>
```

## Entries

### Data / Plots, intervals, events per transition   (2026-11-17)
- Question: What does the record contain?
- Old work: `data/tables/final/{exclusion_flow,transition_set}.csv`, `figures/A/A1a_observed_events` -- partly: counts correct and consistent (8,770 events; 8,405 in 90 modelled transitions); absent: interval length, visits per plot, years; nothing in `analysis/`.
- Choices (Lukas): all 156 origin-destination pairs with a modelled flag; counts only (no at-risk denominator; time at risk is row 4); record description included; exclusion flow included; heatmap included.
- Done: `R/analysis/0_01_record_contents.R`, outputs `data-output/final/analysis/0_data/0_01_{exclusion_flow,events_by_transition,record_summary,visits_per_plot,interval_length,events_heatmap}`; the script stops unless exclusion-flow final row, plot_info and all 156 pairs match `visit_intervals.rds` / `transition_set.csv`.
- Result:

| Plots | Visits | Intervals | Events | Modelled (90 transitions) | Unmodelled (66 pairs) |
|---|---|---|---|---|---|
| 11,286 | 47,537 | 36,251 | 8,770 (6,025 plots) | 8,405 | 365 |

  Interval length 3-39 y (median 11); visits per plot 2-7 (median 4); first visits 1970-2016, last visits 1978-2025; follow-up per plot 5-54 y (median 39).
- Caveat (what can be claimed): an event is a start/end state difference, so paths hidden inside an interval (A->B->A, A->B->C) are not counted; events are lower bounds on real changes.

### Data / Observed change between visits (thermal index, shade tolerance)   (2026-11-17)
- Question: What did the forest actually do?
- Old work: `R/25_I1_raw_thermophilization.R`, `figures/I/I1_*`, `I1b_*` -- partly: index change per decade (stays = 0, plot-bootstrap CI, groups from `00_groups.R`) correct; absent: observed start -> end occupancy, per-state and tolerance-group resolution; period split overlaps row 4.
- Choices (Lukas): mean of per-interval change per decade; groups = thermal group, 13 origin states, shade-tolerance group; splits = all intervals, period of end visit, latitude tertile; each interval once; 13 x 13 occupancy matrix, all intervals (n and median length per origin); `TOL_SCORE` added to `R/00_groups.R` (vector from `25_I_common.R`, verified identical).
- Done: `R/analysis/0_02_observed_change.R`, outputs `data-output/final/analysis/0_data/0_02_{observed_change,occupancy_matrix,state_index}.csv`, `0_02_change_*` and `0_02_occupancy_heatmap` (png/pdf). Script stops unless its point estimates equal old I1 (25 cells) and I1b (15 cells) and the matrix equals row 1's event counts; CIs: plot bootstrap, 200 draws, seed 25.
- Result: change per decade, all intervals (95% CI):

| Origin | n intervals | Thermal index | n with tolerance | Shade tolerance |
|---|---|---|---|---|
| all | 36,251 | +0.008 (0.005, 0.011) | 26,167 | +0.047 (0.039, 0.055) |
| boreal | 15,183 | +0.054 | 13,532 | -0.158 |
| pioneer | 5,469 | -0.055 | 5,070 | +0.605 |
| temperate | 10,229 | -0.086 | 7,565 | +0.041 |
| regen | 5,370 | +0.121 | 0 | NA |

  Stays: 56% (white spruce) to 90% (sugar maple) of intervals; median interval 10-12 y.
- Caveat (what can be claimed): the thermal index of a state is the mean climate of the plots in it (the same plots), so thermal change measures moves between states' climate positions, not warming within plots; tolerance change exists for 72% of intervals (regen and other temperate unscored).

### Data / Disturbance frequency, co-occurrence, climate position   (2026-11-17)
- Question: What disturbance regime is behind the data?
- Old work: `figures/A/A1f_*`, `R/25_I6_disturbance_regime.R` (`figures/I/I6_*`) -- partly: co-occurrence (I6a), climate position of new events (I6b) and return intervals (I6c) correct; A1f (flags by start state, end-visit flags) correct; absent: new vs legacy split, time since disturbance (TSD) overall, thermal and tolerance groups, period and latitude splits.
- Choices (Lukas): frequency = share of intervals flagged, split new / legacy; new = an event of that type starts at the interval's end visit (event table, as I6b) and its flag is active, legacy = every other flag (incl. events active at a plot's first visit); groups = network, origin thermal group, origin tolerance group, 13 origin states; splits = all intervals, period of end visit, latitude tertile; co-occurrence = flag combinations (new/legacy per flag) and new events starting together; climate position = new events only; return intervals as I6c; TSD as median and percentiles, new vs legacy. Climate position and return intervals are network-wide.
- Done: `R/analysis/0_03_disturbance_regime.R`, outputs `data-output/final/analysis/0_data/0_03_{flag_shares,cooccurrence_flags,cooccurrence_new_events,tsd_flagged,climate_position_new_events,return_intervals,event_accounting,superseded_new_events}.csv` and `0_03_{flag_shares_by_group,flag_shares_by_split,climate_position,tsd_flagged}` (png/pdf). Script stops unless it reproduces old A1f (52 cells), I6a (11 combinations), I6b (15 cells) and I6c (590 gaps); CIs for climate position: plot bootstrap, 200 draws, seed 26; the figure shows the occurrences themselves (median and 5th-95th percentile of the z-score at each new event, and across all intervals); the CSV keeps the mean with its bootstrap interval (`est`, `lo`, `hi`, as old I6b) and adds the occurrence quantiles (`p05`, `p25`, `median`, `p75`, `p95`, `n_obs`). TSD figure and CSV split legacy flags into dated (real event year) and imputed (stand age, typical first-visit age or the midpoint of the bounds).
- Result: share of the 36,251 intervals flagged at the end visit:

| Type | Flagged | New | Legacy | Legacy fraction of flags |
|---|---|---|---|---|
| fire | 13.5% | 1.5% | 12.0% | 89% |
| pest | 18.7% | 5.1% | 13.7% | 73% |
| light harvest | 23.5% | 6.4% | 17.1% | 73% |
| heavy harvest | 21.5% | 5.6% | 15.9% | 74% |

  30% of intervals have no flag; 7.5% have two or more flags active together; new events of two or more types start in the same interval in 134 of 6,581 intervals with a new event. TSD median: new 5-6 y; legacy fire 64 y, pest 31 y, light harvest 22 y, heavy harvest 31 y. Climate position of new events (z-GDD): fire -0.94, pest -0.34, light harvest +0.62, heavy harvest +0.05. Return interval between dated events: median 20 y (n = 590).
- Caveat (what can be claimed): legacy TSD is partly imputed; dated share of legacy flags: fire 72%, pest 33%, light harvest 27%, heavy harvest 58% (`n_dated` in `0_03_tsd_flagged.csv`). Of 12,033 events in retained plots, 5,177 start at the first visit (all legacy), 138 start at a visit where no retained interval ends (134 after the plot's last retained interval, 3 at a visit skipped inside a longer interval, 1 dated after the plot's last visit): the 135 outside any interval never appear in any flag; the 3 inside an interval stay legacy (Lukas: probably back-dated) and 9 new events are superseded by a later stand-replacing event in the same interval (new events, no active flag).
- Events before the first visit: 2 events are dated before the plot's first visit (row "before the first visit" in `0_03_event_accounting.csv`); they start at no interval's end visit, so they are not new events, and where their flag is active it is legacy.

### Data / Raw transition rates over calendar time   (2026-11-17)
- Question: Is the process time-homogeneous over the record?
- Old work: `R/25_I2_rates_by_period.R` (`figures/I/I2_*`) -- partly: estimator correct (observed changes / plot-years at risk in the origin, per 1000), but two periods only, assigned by end-visit year; Poisson Wald CI treats the intervals of a plot as independent; absent: thermal and tolerance groups, overall exit rate per period, a reference for "time-homogeneous".
- Choices (Lukas): four periods (<= 1990, 1991-2000, 2001-2010, > 2010) by interval midpoint year; CI = plot bootstrap (200 draws, seed 26, 2.5-97.5%); resolution = 13 states, 4 thermal groups, 4 tolerance groups, exit rate per origin and for the network; no p-value; heatmap per period (states: the 90 modelled routes; thermal: all 12 routes) and forest plot of the 20 busiest state routes and all thermal routes. Forest reference (my reading of "time-homogeneous", not a separate question): period rate / the route's rate over the whole record (1 = constant rate).
- Done: `R/analysis/0_04_transition_rates.R`, outputs `data-output/final/analysis/0_data/0_04_{rates,period_summary}.csv` (rates: level state / thermal / tolerance / network, route and exit rates, bootstrap CI, ratio to pooled) and `0_04_{rates_state_heatmap,rates_thermal_heatmap,rates_state_forest,rates_thermal_forest,exit_rates_state,exit_rates_thermal}` (png/pdf). Script stops unless events and plot-years reproduce old I2 under the old two-period rule (312 event cells, 26 plot-year cells), route events equal `transition_set.csv` n_events, and 90 routes are modelled.
- Result: network exit rate per 1000 plot-years is almost constant: 22.4 (<= 1990), 23.0, 20.9, 21.0 (> 2010); pooled 21.8. Below the network level it is not constant: thermal exit rates (<= 1990, 1991-2000, 2001-2010, > 2010) boreal 22.8, 13.6, 11.0, 13.5; pioneer 13.7, 26.7, 33.9, 26.9; regen 31.6, 39.7, 33.0, 34.0; temperate 8.0, 9.0, 7.5, 7.2. Largest state shifts: fir 46.1, 22.2, 14.2, 20.1 (pooled 26.0); white birch 14.7, 31.0, 36.8, 34.9 (pooled 28.4). Of the 90 modelled routes, the period rate differs from the route's pooled rate (CI excludes 1) in 39 (<= 1990), 6 (1991-2000), 30 (2001-2010) and 19 (> 2010) routes (descriptive count, no multiplicity adjustment).
- Caveat (what can be claimed): crude rates ignore paths hidden inside an interval, and mean interval length rises from 10.2 y (<= 1990, 1991-2000) to 12.4 and 11.9 y (later periods), which alone lowers observed rates; the plots measured differ between periods (6,587 / 8,386 / 9,722 / 6,396) and so do disturbance mix and climate, so a change in rate is not by itself a change in the process; the bootstrap CI is degenerate for routes with 0-3 events in a period (see `n`).

### Data / Climate trend in the plots over the record   (2026-11-17)
- Question: How much did climate change during the record?
- Old work: absent (row marked new). Only a sentence in the header of `R/22_hindcast_space_for_time.R`: GDD rose by a median of +0.49 SD within plots from first to last visit; no table or figure.
- Choices (Lukas): series = yearly raw climate per plot-year 1960-2022 (GDD = gdd_3; CMI = cmi_gs rebuilt with the 03c step-4b weighting rule, copied unchanged); soil PC1 is constant within plots, so no trend; slope per plot = (last year - first year value) / years between first and last VISIT year x 10 (last visit after 2022 clipped to 2022); main estimate = mean of plot slopes over plots with a span >= 10 y; checks: least-squares slope over the same years, mean over all plots with a span > 0; units raw and z (model scaling); splits all plots, period (of the span midpoint, cuts as rows 3-4), latitude tertile; plot bootstrap (200 draws, seed 26); figures: yearly series with plot sample, slopes by split, network mean with linear trend.
- Done: `R/analysis/0_05_climate_trend.R`, outputs `data-output/final/analysis/0_data/0_05_{plot_slopes,slopes_summary,network_series,network_trend,series_check}.csv` and `0_05_{series_network,slopes_by_split_raw,slopes_by_split_z}` (png/pdf; the network figure draws 1,000 random plots, seed 26). Script stops unless the rebuilt yearly series averaged over each interval's window reproduces gdd_3 and cmi_gs (max abs diff 1.4e-14, 36,251 intervals).
- Result: within-plot change per decade, plots with span >= 10 y (11,000 of 11,286 plots; 2,171 plots clipped to 2022): GDD +53.9 (+0.178 SD; CI 53.2 to 54.7), CMI -0.40 (-0.056 SD; CI -0.45 to -0.35). Checks: least-squares slope GDD +57.4, CMI -0.32; all plots GDD +52.1, CMI -0.49. By period (z per decade) GDD: <= 1990 0.013 (n 277), 1991-2000 0.145, 2001-2010 0.247, > 2010 0.100 (n 122); CMI: 0.151, -0.045, -0.084, -0.162. By latitude GDD: south 0.163, middle 0.194, north 0.178; CMI: -0.070, -0.144, +0.046. Network yearly mean 1960-2022: GDD +44.6 per decade (+0.148 SD), CMI +0.107 (+0.015 SD). Old statement reproduced: median +0.49 SD between a plot's first and last interval covariate (plots with >= 2 intervals).
- Caveat (what can be claimed): GDD rose in every split except the first period, by about 0.15-0.25 SD per decade. The CMI trend is not established: its sign depends on the span (plots' own years -0.40 per decade, 1960-2022 network +0.11), and the north band has the opposite sign to the middle band. The bootstrap CIs cover plot sampling only; plots share the same weather years, so the weather of the endpoint years is not in the CI (the violins show modes set by shared visit years), and the network trend CI (GDD 44.5 to 44.7) is far too narrow as a statement about the trend itself. The first (n = 277) and last (n = 122) periods are small.
- Dropped plot: 1 plot with a span >= 10 y has no yearly value at one endpoint year and is excluded from the slope (11,001 plots with a span >= 10 y -> 11,000).
