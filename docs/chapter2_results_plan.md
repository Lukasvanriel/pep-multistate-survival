---
title: "Chapter 2 — results plan (proposal for supervisors)"
date: "11 November 2026"
---

# 1. Framing

**Working title.** *Disturbance sets the timing and climate the destination of forest reorganisation at a boreal–temperate ecotone*; alternative, closer to the mechanism: *Regeneration is where climate and disturbance meet: timing, direction and memory of forest transitions across an ecotone*.

**The story in one paragraph.** Across 11,300 plots and 50 years, the forest has barely thermophilized (+0.008 SD of GDD per decade, confined to the south before 1990) but has clearly succeeded toward shade tolerance. A semi-Markov transition model explains why. Transition rates carry memory: they rise with time in state, most for stand rebuilding. Their climate effects largely mirror the climatic sorting of the states along the gradient (reverse pairs close to detailed balance), so the direction of change is already written into the gradient; what warmth adds beyond that sorting is one thing, faster closure of regenerating stands. Disturbance decides when a stand leaves, with legacies that change sign over decades; its interactions with climate are antagonistic and confined to regen exits (×0.53) and light harvest (×0.74). Network-wide, transitions moved in step with regional warming as the spatial climate effects predict, except for temperate origins; a response specific to each plot's own warming is shown only for the pace of boreal stands. Under the current disturbance regime, mostly through harvest, the rates point toward more boreal stands in the warm half; without disturbance, toward more temperate stands; the landscape sits between. Reorganisation at this ecotone is gated by timing (disturbance, time in state, regeneration), not by the lack of a direction.

**Contribution.** A transition model across the ecotone that (i) tests Markov vs semi-Markov dynamics per transition, (ii) separates pace from direction for climate and for disturbance, (iii) tests the climate × disturbance interaction directly, (iv) estimates the time course of disturbance legacies, and (v) tests the spatial climate signal against observed within-plot change with two placebos. *Novelty claims ("first …") to be checked against the literature before use.*

# 2. Outline

**Introduction (five paragraphs)**

1. Tree communities lag climate; three candidate gates (slow intrinsic rates, disturbance, dependence on stand history) are rarely tested together.
2. Transition models of forest change are mostly Markov, treat climate and disturbance additively, and read spatial climate effects as temporal ones: three assumptions with contested support.
3. The regeneration phase as the candidate hinge: where disturbance acts, where growth responds to warmth, where memory should matter.
4. Questions: (i) is time-in-state dependence general, and where is it strongest? (ii) do climate and disturbance act on pace or on direction, and do they interact? (iii) what is the time course of disturbance legacies? (iv) do spatially estimated climate effects predict within-plot change?
5. System and approach.

*References for the introduction to be confirmed as primary sources.*

**Results (seven blocks, figures F1–F7; F8 in the Discussion)**

| Block | Content | Figure / table | Sources in the review |
|---|---|---|---|
| R1 What the record shows | raw thermal and successional change per decade, by origin group, period and latitude band | F1; T1 (plots, visits, events per state and disturbance) | I1, I1b, A1 |
| R2 Rates have memory | α by class; hazard curves of emblematic routes | F2 | A3 |
| R3 Disturbance sets the timing; legacies change sign | class-level disturbance effects at 1–60 y; TSD fade by class | F3 | A5, A6B, A6D |
| R4 Climate direction is sorting; climate pace is regen closure | detailed-balance test; class-specific pace intercepts | F4 | A4, A4f, I10 |
| R5 Pace vs direction on the real plots; balance point; first passage | D3 all-plots effects; D2 asymmetry along GDD; I14 first passage within observed range | F5; T2 (home-climate persistence and first destinations, I7) | D2, D3, I14, I7, I15 |
| R6 Climate × disturbance meet in regen | interaction by class; outcome-scale interactions; recovery vs reorganisation by GDD | F6 | A6C, B7, D1 |
| R7 Transfer to time | within-plot hindcast slopes with two placebos and extra controls | F7 | H, I9 |

**Supplementary information.** I6 disturbance climates; I12 fire/pest bounds; I15 elasticity table; I2 rates by period; I4 hysteresis diagnostic (and the refit, when run); group-level states as a sensitivity analysis (S27 in the sensitivity plan); B0 model check; C4 pathways; I8/I13 details; all per-transition hazard ratios (section 4).

**Discussion (seven points)**

1. Reorganisation is gated by timing; direction is already in the gradient: equilibrium in space, lag in time.
2. The regeneration hinge: strongest climate pace effect (I10), strongest memory (A3), the only strong climate × disturbance interaction (A6C), and the weakest data (undisturbed regen intervals are 9% of regen intervals, 1.3% of all; flicker at the regen threshold).
3. Disturbance defers rather than amplifies climate-driven change, by intensity (reset vs release) and by horizon (legacies change sign).
4. Two equilibria and the landscape between them (F8): without disturbance the rates pull warm stands toward temperate states; the current regime, mostly harvest, pulls them back through fir. Labelled as a 1000-y direction, not a forecast.
5. Memory: Markov models mistime recovery rather than direction.
6. Space-for-time: the network moves in step with regional warming as the spatial effects predict (except temperate origins); plot-specific transfer only for boreal pace; regional warming cannot be separated from other regional trends of the period (stand ageing, the budworm cycle).
7. Limitations, by what they change: flicker at the dominance margin (short-term persistence and returns; hysteresis refit), red spruce within "other temperate" (black spruce at +1 SD; I16), fire and pest ranges for assumption-heavy starts (I12), pioneers under-predicted at equilibrium even with disturbance (I13, cause unknown), extrapolation beyond 55 y and beyond each state's climate (I14), unrecorded planting (I5), fire rate in I13 counts only events inside the record.

# 3. What the chapter answers

- **Why do tree communities lag climate?** Not for lack of a climate signal in the rates: the rates reproduce the climatic sorting of the states (A4f, I10). Replacement between established states is slow and nearly memoryless, and the only climate-accelerated step is regen closure. Observed thermophilization is +0.008 SD per decade (I1), while 59% of stands would have to change state to match the long-run composition at +1 SD (E6).
- **Does disturbance accelerate or reset climate-driven change?** Both, depending on intensity and horizon, and it does not amplify the climate signal: interactions are antagonistic and confined to regen exits and light harvest (A6C). Under the current regime, mostly through harvest, the long-run composition in the warm half is more boreal, not more temperate (I13).
- **Does stand history matter?** Yes, by class: stand rebuilding and replacement depend on time in state; lateral replacement between established states is nearly memoryless (A3).
- **Is space-for-time valid for transition rates?** Partly. Transitions moved in step with regional warming as predicted for boreal, pioneer and regenerating origins, but only boreal pace responds to each plot's own warming beyond the regional trend; temperate origins show no temporal signal (H, I9).

# 4. Where the per-transition hazard ratios go

The 90 transitions × up to 14 terms are the model's raw output, but the questions are about classes and derived quantities, so the main text summarises and the SI reports everything:

- **In the main figures as points:** every α (F2a), every disturbance hazard ratio (F3, points behind the class means) and every GDD hazard ratio (F4b) is plotted; the class means, fitted lines and derived quantities are the message.
- **In the SI as tables:** pooled estimates with 95% CI, Rubin SE and fraction of missing information for every transition and term (pooled_final_final.csv, hazard_ratios_final.csv).
- **In the SI as heatmaps:** baselines (A2), α (A3b), climate and soil (A4a), disturbance vs none (A5e), so a reader can look up any route.

# 5. Figures (drafts)

Drafts are montages of existing panels (F1–F3, F5, F6) or new panels (F4, F7, F8); files in data-output/final/figures/plan/. Final versions need a common style and fewer panels per figure.

## F1 The record

![](../data-output/final/figures/plan/F1_record_draft.png)

**Shows.** Observed change per decade in the thermal index of the state (ΔSTI) and in shade tolerance, for every visit interval, by origin group, period of the end visit (top) and latitude band (bottom). **Why we care.** It is the model-free baseline: what the forest actually did before any model. **Answers.** Network-wide thermophilization is near zero (+0.008 SD per decade), positive only before 1990 and in the south; regen and boreal stands move warm, pioneers and temperate stands cool, and succession toward shade tolerance is clear for pioneers.

## F2 Rates have memory

![](../data-output/final/figures/plan/F2_memory_draft.png)

**Shows.** The Weibull shape α of every transition, grouped by class (top), and hazard curves over time in state for selected routes (bottom). α > 1 = the rate rises with time in the state. **Why we care.** Markov models assume α = 1; if it is not, they mistime transitions. **Answers.** Memory is the rule and is strongest for stand rebuilding (out of regen); replacement between established states is close to memoryless.

## F3 Disturbance sets the timing; legacies change sign

![](../data-output/final/figures/plan/F3_disturbance_timing_draft.png)

**Shows.** Top: hazard ratio of each disturbance vs no disturbance, by transition class, at 1, 10, 30 and 60 y since the disturbance (points = transitions, black = class mean). Bottom: how each effect fades with time since disturbance. **Why we care.** Disturbance is the strongest driver of when stands move; whether its effect lasts or reverses decides recovery. **Answers.** Heavy harvest and fire send stands into regen (class means 210× and 75× at 1 y), an effect that fades (heavy harvest 8× at 30 y); by 30 y, harvest and pest instead speed regen exits (2–4×): legacies change sign over decades. Light harvest into regen falls below 1 by 30 y in the model; the observed rates return to about the undisturbed level instead (sim06g), so we treat this as an artefact of the fitted TSD shape, not as a protective effect.

## F4 Climate direction is sorting; climate pace is regen closure

![](../data-output/final/figures/plan/F4_null_model_draft.png)

**Shows.** (a) For each pair of states with moves both ways: the observed difference in the GDD effect of the warm and the cold move against the difference expected if the rates simply mirror where the states occur along GDD (dashed = detailed balance). (b) Every transition's GDD effect against the same expectation, with a separate line per class; the offset at zero is warmth speeding the move beyond sorting. **Why we care.** A warming effect on a transition could be a real response or a restatement of where the states sit; this separates the two. **Answers.** Direction follows sorting (reverse pairs near detailed balance; between-state moves mirror ~60% of the sorting). The only pace effect beyond sorting is out of regen (~2.7× per SD); moves into regen do not follow sorting at all.

## F5 Pace vs direction on the real plots; balance point; first passage

![](../data-output/final/figures/plan/F5_pace_direction_draft.png)

**Shows.** Top: for the real plots (own climate and soil), the change in the probability of leaving within 50 y and in the thermal direction of the next state, for warming and each disturbance. Middle: probability that boreal stands reach a temperate state within 50 y and the reverse, along GDD, and their asymmetry. Bottom: first passage into the temperate group along GDD by start state. **Why we care.** It turns the coefficients into statements about the actual network and locates where net replacement flips. **Answers.** Heavy harvest has the largest effect on pace (+0.35), warming the largest effect on direction (+0.29). The boreal–temperate balance point is at +0.29 SD (~1,355 degree-days; robust to ±11 degree-days for any single parameter, I15). Entry into the temperate group is a threshold between GDD 0 and +1 SD (fir most open; black spruce and jack pine extrapolated there).

## F6 Climate × disturbance meet in regen

![](../data-output/final/figures/plan/F6_climate_x_disturbance_draft.png)

**Shows.** The climate × disturbance interaction by class (top), on the outcome scale (middle), and what stands become 50 y after a disturbance at GDD −1/0/+1 (bottom). **Why we care.** The catalyst hypothesis says disturbance opens the door for climate-driven change; this is the direct test. **Answers.** No amplification: interactions are antagonistic, confined to regen exits (×0.53) and light harvest (×0.74). Disturbance halves the climate signal where it matters most. Final version: bottom panel for fir, white birch and regen only.

## F7 Transfer to time

![](../data-output/final/figures/plan/F7_hindcast_placebo_draft.png)

**Shows.** Within-plot slope of the observed change on the change the spatial model predicts from each plot's own warming (1 = full transfer, 0 = none), by origin group and outcome, with two placebo ranges: own warming in shuffled order (blue) and another plot's warming from the same region and period (orange). **Why we care.** Every climate scenario assumes effects learned across space apply in time. **Answers.** Real slopes exceed the shuffled-order placebo for boreal, pioneer and regen origins, so the response follows the timing of warming. Against the other-plot placebo, only boreal pace (any change, end boreal) is clearly specific to each plot; the direction slopes are close to it, so the network responds to regional warming, not separable from other regional trends. Temperate origins: no signal.

## F8 Two equilibria and the landscape between them (Discussion)

![](../data-output/final/figures/plan/F8_two_equilibria_draft.png)

**Shows.** Share of boreal, regen and temperate stands along GDD: observed now, and the 1000-y composition implied by the rates without disturbance, with the observed disturbance regime, and with fire and pest only. **Why we care.** It shows where the rates point and how much disturbance, and harvest in particular, changes that. **Answers.** Without disturbance the warm half would become more temperate than now (Q5: 86% vs 70%); with the observed regime less (60%); without harvest it returns close to the undisturbed case (80%). The current harvest regime is the main boreal-ward force in the model; regen shares are matched only when disturbance is included.

# 6. Changes from the earlier outline

1. Space-for-time narrowed after the cross-plot placebo (I9d): the response tracks regional warming; plot-specific transfer only for boreal pace; "placebo near zero" removed.
2. Climate × disturbance interaction: light harvest (×0.74) added next to regen exits (×0.53).
3. "Management is a boreal-ward force" tested with a no-harvest regime and with light-only and heavy-only regimes (I13): supported in the model; both harvest types are boreal-ward at equilibrium, heavy about twice as strongly.
4. Novelty claims ("first …") flagged for a literature check.
5. Limitations: pioneer under-prediction at equilibrium (I13) and the I13 fire-rate note added.
6. "Undisturbed regen is 9% of intervals" corrected to 9% of regen intervals (1.3% of all).
7. Hysteresis (I4) treated as a limitation with an SI refit; it does not change the structure of the chapter.
8. Note for the manuscript: TSD shape. Log kept as the main form; an exponential shape (τ per disturbance type × transition class) fits better by WAIC but changes the I13 equilibrium shifts by ≤ 3 percentage points (sim06h, I13e).

# 7. Links to the literature

*Checked against the papers or their abstracts on 11 November 2026, except where marked. Press-release-only items are marked and need the paper before citing.*

## 7.1 Direct predecessors on the same network

**Brice et al. 2019, *Global Ecology and Biogeography* 28:1668–1681.** 6,281 PEP plots, first (1970–80) vs last (2000–16) survey; community temperature index (CTI) and shade index (CSI).

| | Brice et al. 2019 | This chapter | Relation |
|---|---|---|---|
| Thermophilization | +0.03 °C per decade vs +0.14 °C per decade warming | +0.008 SD of GDD per decade (dominant state) | agree: small and lagging |
| Where | positive in sugar maple and fir–yellow birch domains | south band, before 1990 | agree |
| Climate change as driver | climate predictors explain < 1% of ΔCTI | temporal response tracks regional warming; plot-specific only for boreal pace (I9) | agree; ours locates where transfer holds |
| Moderate vs major disturbance | moderate tripled ΔCTI; major → pioneers (CSI −0.037) | first move: light harvest warm-ward, heavy harvest cold-ward and resets to regen (D3) | agree |
| "Amplify" | disturbance speeds thermophilization (main effect); climate × disturbance not tested beyond one term | climate × disturbance interaction antagonistic (A6C) | refines: disturbance sets the pace in a direction the gradient already sets |
| Long run | not addressed | repeated light and heavy harvest are both boreal-ward at equilibrium (I13) | new; leans partly on the unaffected-exit assumption |
| Species | gains of red and sugar maple; losses of fir and black spruce | fir ↔ temperate exchanges control the balance point (I15) | consistent |

Hypothesis for the Discussion (untested): their comparison spans the 1970s–80s budworm outbreak, and fir losses are among the largest contributors to their ΔCTI; our thermophilization is confined to before 1990 (I1) and the fir cycle reverses after 2000 (I2). One text point: their stand-age effect (older stands change less) and our time-in-state effect (α > 1) are different clocks.

**Brice et al. 2020, *Global Change Biology*.** ~11,000 PEP plots, continuous-time Markov multistate model (msm), 4 states (boreal, mixed, pioneer, temperate); moderate disturbances "catalyse transitions to the alternate, temperate-dominated forest state"; major disturbances push toward pioneer; the authors state the model "is not designed to make future range shift projections". The chapter's contribution is stated relative to this paper: semi-Markov rates tested per transition, 13 states, legacy time courses, a direct climate × disturbance test, a sorting null model and a placebo-tested hindcast.

**Vieira et al. 2024, *Ecological Modelling* 497:110813** (with Brice, Gravel). Theoretical state-and-transition model: harvesting boreal stands did not move temperate trees north, partly because harvested boreal stands regenerate boreal. Our rate-based I13 result is the empirical counterpart (fir funnel).

## 7.2 Why communities lag climate (R1, R4, R5)

- Danneyrolles et al. 2019, *Nature Communications* 10:1265: no clear evidence of direct climate impacts on southern Quebec composition since pre-settlement; land-use disturbance dominated.
- Talluto et al. 2017, *Nature Ecology & Evolution* 1:0182: extinction debt and colonization credit for all 21 eastern species.
- Sharma et al. 2022, *PNAS* 119:e2116691118: eastern tree migration lagging; fecundity as a constraint.
- Rosenblad et al. 2023, *PNAS* (press summary read; paper to check): western US thermophilization lags climate about tenfold, driven by mortality.
- This chapter: the rates already carry the direction of the gradient (A4f); the one climate-accelerated step is regen closure (I10). The lag is in the pace of a phase, not in the absence of a direction.

## 7.3 Disturbance and climate-driven change (R3, R6, F8)

- Germain et al. 2025, *Ecosystems* 28:9: disturbance creates recruitment opportunities for warm-adapted trees (NE US).
- Robles et al. 2025, *Landscape Ecology* 40:26: harvest the main driver of composition since 1830 in Témiscamingue (aspen, red maple up).
- Boucher et al. 2017, *Journal of Vegetation Science* 28:57–68: fire and budworm stronger drivers than logging in the central Quebec boreal.
- Donnini & Kross 2026, *Journal of Forestry Research* (press release only): conifers, especially balsam fir, increased across Quebec PEP plots 1985–2021, attributed to budworm recovery; consistent with I2.
- Baltzer et al. 2021, *PNAS*: increasing fire undermines black spruce regeneration; with Quebec's 2023 fire season, a reason the 1970–2020 regime of I13 is not a forecast.

## 7.4 Memory (R2, R3)

- Johnstone et al. 2016, *Frontiers in Ecology and the Environment*: ecological memory and disturbance legacies, conceptual. Forest transition models remain mostly time-homogeneous Markov (Brice et al. 2020 included); A3 and the sign-reversing legacies are a direct test.

## 7.5 Space-for-time (R7)

- Lovell et al. 2023, *Biological Reviews*: "key assumptions remain largely untested". I9 is such a test, with a specific answer (regional temporal response; plot-specific only for boreal pace; none for temperate origins). The cross-plot placebo shows that regional warming is confounded with other regional trends, a methods point for other hindcasts.
- To add after checking: Blois et al. 2013; Damgaard 2019; Svenning & Sandel 2013.

## 7.6 Value beyond Quebec

1. Methods that transfer to any remeasured plot network (FIA, European national inventories): semi-Markov interval-censored transition model with multiple imputation; a detailed-balance null model separating climatic sorting from response; two placebos for space-for-time.
2. Climatic debt: the lag sits in the pace of one phase (regen closure), a testable mechanism.
3. Disturbance is not a universal catalyst: its effect depends on intensity and horizon; under the current regime, harvest opposes the climate signal in the warm half (a management lever).
4. Space-for-time: names the states for which spatially parameterised climate scenarios are and are not supported, relevant to every inventory-based projection model (including LANDIS-II, Chapter 3).
