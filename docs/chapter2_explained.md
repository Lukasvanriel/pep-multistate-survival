---
title: "Chapter 2 explained: the logic, the model, the figures"
date: "11 November 2026"
---

This document walks through the chapter from the question to the conclusion. Read it top to bottom once; the figures are drafts (montages of existing panels or new panels) and will be restyled.

# 1. The problem in one page

**The observation.** Climate at the boreal–temperate ecotone of Quebec has warmed for decades, yet forest composition has barely followed. Brice et al. (2019) measured a community thermophilization of +0.03 °C per decade against +0.14 °C per decade of warming, and climate predictors explained under 1% of it. Danneyrolles et al. (2019) found no climate signal in two centuries of change in southern Quebec; land use dominated. Talluto et al. (2017) found extinction debt and colonization credit for all 21 dominant tree species of eastern North America.

**What the literature leaves unclear.** That forests lag climate is established. *Why* they lag, and what disturbance does to that lag, is not. Four specific gaps:

1. **Where does the lag sit?** Are the rates of change simply insensitive to climate, or do they respond but slowly, and if so in which phase of stand development?
2. **What does disturbance do?** Brice et al. (2019, 2020) found that moderate disturbances speed up thermophilization and transitions toward temperate stands, while major ones reset stands to pioneers. Whether disturbance *amplifies the climate effect itself* (a synergy) or only speeds up change was not tested. Nor is it known what a recurring disturbance regime does in the long run.
3. **Do rates depend on how long a stand has been in its state?** Forest transition models, including Brice et al. (2020), assume they do not (Markov models). Ecological memory and disturbance legacies (Johnstone et al. 2016) suggest they should.
4. **Can climate effects learned across space be used for change in time?** Every climate scenario assumes so (space-for-time). A recent review concludes that its "key assumptions remain largely untested" (Lovell et al. 2023).

**What the chapter does.** It fits a model of *transition rates* between forest states on the full Quebec permanent-plot network, so that each of these four questions becomes a test on the rates.

# 2. Data and states

- **Plots.** 11,286 permanent plots of the Quebec network (PEP), 2–7 visits each (median 4), first visits 1970–2016, intervals 3–39 years (median 11).
- **States.** At each visit a plot is assigned to one of 13 states: the dominant species or species group by basal area (black spruce, white spruce, fir, jack pine, white birch, trembling aspen, other pioneer, yellow birch, red maple, sugar maple, cedar, other temperate), or **regen** when basal area is below 5 m²/ha. States are grouped as boreal, pioneer, temperate and regen.
- **Changes.** 8,770 state changes across 156 possible origin–destination pairs; the 90 transitions with at least 15 events (96% of changes) are modelled.
- **Covariates.** Growing-degree-days (GDD) and growing-season climate moisture index (CMI), each averaged over the 10 years before the visit (NRCan ANUSPLIN), a soil axis (pH and drainage), and four disturbance types from the inventory codes (fire, pest, light harvest, heavy harvest), each with its time since disturbance (TSD).

# 3. The model, and why this one

**What it is.** A *multistate survival model*: for each of the 90 transitions, a hazard (instantaneous rate) of moving from state *i* to state *j*:

> h(t) = α t^(α−1) · exp(intercept + climate, soil and disturbance effects),

where *t* is the time the stand has spent in state *i*. The exits from a state compete: whichever happens first decides where the stand goes. Each transition is fitted separately (cause-specific hazards) in INLA (Bayesian).

**The three features that matter, and why each is needed:**

| Feature | What it does | Why we need it |
|---|---|---|
| **Weibull shape α** (semi-Markov) | the rate can rise (α > 1), stay constant (α = 1) or fall (α < 1) with time in the state | tests gap 3 directly; a Markov model fixes α = 1 |
| **Covariates on each transition** | climate and each disturbance change the rate of each specific route, with GDD × disturbance interactions and TSD legacies | separates **pace** (how fast stands leave) from **direction** (where they go); tests gap 2 |
| **Interval censoring + multiple imputation** | changes happen between visits, at unknown times; exact times are imputed repeatedly from the model and the fits are pooled (Rubin's rules) | uses the real timing information without pretending changes happen at visits |

**Why not the existing approaches.**

- *Comparing composition between two surveys* (Brice et al. 2019) measures how much changed, but not the rates, so it cannot say whether climate acts on pace or direction, and cannot be run forward.
- *A Markov multistate model* (Brice et al. 2020, 4 states) estimates rates but assumes no memory and lumps states (boreal, mixed, pioneer, temperate), which hides routes such as fir ↔ white birch or black spruce → red spruce stands.
- *Process models* (e.g. LANDIS-II, Chapter 3) simulate mechanisms but take their parameters from elsewhere; this chapter estimates the empirical rates such models should reproduce.

**What the model gives.** Two layers of results:

1. **Fitted parameters**: α, hazard ratios for climate, soil and disturbance, legacy (TSD) terms, for every transition.
2. **Derived quantities**: running the fitted rates forward from a starting state gives probabilities of leaving, where stands go, outcomes after a disturbance, and long-run composition. These are scenarios ("what the current rates imply"), not forecasts.

# 4. The logic of the chapter: question → test → result

| Question | Test | Result | Figure |
|---|---|---|---|
| What did the forest actually do? | raw change per decade in the thermal and shade-tolerance index of each plot's state | almost no thermophilization (+0.008 SD of GDD per decade, only in the south and before 1990); clear succession toward shade tolerance | F1 |
| Do rates depend on time in the state? (gap 3) | α per transition | yes: α > 1 wherever the data can tell (median 1.45, never below 1), strongest for stand rebuilding (regen exits) | F2 |
| What does disturbance do to the rates, and for how long? (gap 2) | disturbance hazard ratios by class over time since disturbance | disturbance decides *when* stands move; legacies change sign: first a push into regen, later a faster exit from regen | F3 |
| Does warming change pace or direction, and is that more than spatial sorting? (gap 1) | null model: are the GDD effects just a mirror of where states occur along the gradient? | direction is mostly sorting (near detailed balance); warmth adds pace only to regen exits (~2.7× per SD) | F4 |
| What does this mean for real stands? | effects on the real plots; where net boreal ↔ temperate replacement flips; first passage into temperate | disturbance sets the pace, warming the direction; balance point at +0.29 SD (~1,355 degree-days); temperate entry is a threshold between GDD 0 and +1 | F5 |
| Does disturbance amplify the climate effect? (gap 2) | GDD × disturbance interactions, on rates and outcomes | no: interactions are antagonistic (regen exits ×0.53, light harvest ×0.74); disturbance overrides climate at the first move | F6 |
| Do spatial climate effects predict change in time? (gap 4) | within-plot hindcast with two placebos | transitions follow regional warming as predicted (not temperate origins); a response to each plot's own warming only for boreal pace | F7 |
| What do the rates imply in the long run? | 1000-year composition with and without the observed disturbance regime | without disturbance the warm half would turn more temperate; the current regime, mostly harvest, pulls it back toward boreal | F8 (Discussion) |

# 5. The figures, one by one

For each: what is in it, how to read it, what it shows, and what to keep in mind.

## F1 What the record shows

![](../data-output/final/figures/plan/F1_record_draft.png)

**Contains.** Top: observed change per decade in the thermal index (ΔSTI, left) and shade tolerance (right) of each plot's state, for every visit interval, by origin group and period. Bottom: the same by latitude band. Dots are means with 95% bootstrap intervals; zero = no change.

**How to read.** A point above zero in the thermal panel means stands of that group moved, on average, to states that live in warmer places. "All" (black) is the network.

**Shows.** The network as a whole barely moved warm (black points near zero), and only in the south and before 1990. Regen and boreal stands move warm-ward, pioneers and temperate stands cool-ward: the network churns, but the churn cancels. Shade tolerance rises for pioneers: succession is visible.

**Keep in mind.** This is the model-free baseline; it matches Brice et al. (2019) at a coarser (state) level.

## F2 Rates have memory

![](../data-output/final/figures/plan/F2_memory_draft.png)

**Contains.** Top: the Weibull shape α of every transition with its 95% CI, grouped by class. Bottom: hazard curves over time in the state for selected routes.

**How to read.** α = 1 (dashed line) is a memoryless (Markov) rate. α > 1: the longer a stand has been in a state, the more likely it is to leave.

**Shows.** α is above 1 for 48 of 90 transitions and below 1 for none (median 1.45); the fraction above 1 rises to 77% where there are ≥ 50 events, so the limit is power. Strongest for regen exits (median 1.73) and fir; black spruce and sugar maple are nearly memoryless (≈ 1), consistent with climax-like states.

**Keep in mind.** For regen exits, time in state equals time since disturbance, so α and the TSD terms share the signal.

## F3 Disturbance sets the timing; legacies change sign

![](../data-output/final/figures/plan/F3_disturbance_timing_draft.png)

**Contains.** Top: hazard ratio of each disturbance vs none, by transition class (out of regen, between established states, into regen), at 1, 10, 30 and 60 years since disturbance; coloured points = transitions, black = class mean, log scale. Bottom: how each effect changes with time since disturbance (< 1 = fades).

**How to read.** A value of 100 means the route runs 100× faster than without the disturbance at that time since disturbance.

**Shows.** Heavy harvest and fire send stands into regen (class means 210× and 75× at 1 year); the effect fades (heavy harvest 8× at 30 years). By 30 years, harvest and pest instead speed the *exit* from regen (2–4×). The same disturbance first resets, then accelerates recovery.

**Keep in mind.** 60 years is outside the observed time-since-disturbance range for harvest and pest; fire has terms in only 14 transitions.

## F4 Climate direction is sorting; climate pace is regen closure

![](../data-output/final/figures/plan/F4_null_model_draft.png)

**Contains.** (a) For each pair of states with moves in both directions: the observed difference in the GDD effect (warm move minus cold move) against the difference expected if the rates only reflect where the two states occur along GDD. Dashed line = detailed balance. (b) Every transition's GDD effect against the same expectation, with a fitted line per class.

**How to read.** Points on the dashed line in (a): warming tilts the exchange between two states exactly as their spatial sorting implies, nothing more. In (b), the height of a class's line at zero on the x-axis is the speed-up warmth gives regardless of direction.

**Shows.** Reverse pairs sit close to detailed balance: the climate *direction* of the rates is mostly the spatial sorting of the states. Beyond that, warmth speeds up only the exit from regen (~2.7× per SD of GDD); between established states the extra pace is small, and moves into regen do not follow sorting at all (they are disturbance and mortality).

**Keep in mind.** This is why "warming favours temperate states" is not by itself evidence of a response: it is what sorting predicts.

## F5 Pace vs direction on real plots; balance point; first passage

![](../data-output/final/figures/plan/F5_pace_direction_draft.png)

**Contains.** Top (D3): for the real plots with their own climate and soil, the change in the probability of leaving within 50 years (left) and in the thermal direction of the next state (right), for +1 SD of GDD, CMI, soil, and each disturbance. Middle (D2): probability that boreal stands reach a temperate state within 50 years and the reverse, along GDD, and the asymmetry between the two. Bottom (I14): probability of entering the temperate group within 100 years along GDD, by start state.

**How to read.** Top: a bar far from zero on the left = that perturbation changes how fast stands leave; on the right = where they go. Middle: where the asymmetry crosses zero is the *balance point*, the climate where boreal → temperate and temperate → boreal replacement are equally likely. Bottom: steepness = how sharply the gradient switches stands on.

**Shows.** Heavy harvest has the largest effect on pace (+0.35 probability of leaving), warming the largest on direction (+0.29): **disturbance sets the pace, climate the direction**. The balance point lies at +0.29 SD (~1,355 degree-days), near the middle of the network, and no single parameter moves it by more than ~11 degree-days (I15). Entering a temperate state is a threshold between GDD 0 and +1 SD, most open from fir.

**Keep in mind.** Undisturbed, fixed climate, 50–100-year horizons; black spruce and jack pine at +1 SD and beyond are outside their observed climate.

## F6 Climate × disturbance meet in regen

![](../data-output/final/figures/plan/F6_climate_x_disturbance_draft.png)

**Contains.** Top (A6C): the GDD effect under each disturbance relative to undisturbed, by class (1 = additive). Middle (B7): the same interaction on outcomes (probability of leaving, thermal direction), per state and disturbance; filled = CI excludes 0. Bottom (D1): what a stand becomes 50 years after a disturbance, at GDD −1 / 0 / +1.

**How to read.** Top: below 1 = disturbance weakens the climate effect. Middle: negative = warming matters less when the stand is disturbed.

**Shows.** No amplification. The climate effect is weakened under disturbance in regen exits (×0.53) and light harvest (×0.74); on outcomes every significant interaction for direction is negative. Mechanism: heavy disturbance sends stands to regen whatever the climate; climate gets its say later, in the regen phase.

**Keep in mind.** The bottom panel will be reduced to fir, white birch and regen; fire and pest outcomes for some temperate starts rest on an assumption (exits without a term keep their undisturbed rate; bounds in I12).

## F7 Do climate effects learned across space predict change in time?

![](../data-output/final/figures/plan/F7_hindcast_placebo_draft.png)

**Contains.** For every observed interval, the model predicts the end state twice: with the plot's actual climate and with its climate frozen at the first visit. The slope of (observed − frozen) on (actual − frozen), within plots, is shown as points with 95% CIs, by origin group and outcome. Blue band: placebo with the plot's own warming in shuffled order. Orange band: placebo with another plot's warming from the same region and period.

**How to read.** Slope 1 = stands changed as much as the spatial effects predict from their warming; 0 = no relation. A real slope above the blue band = the timing of warming matters. Above the orange band = the plot's *own* warming matters beyond the regional trend.

**Shows.** For the network and for boreal, pioneer and regen origins the real slopes exceed the shuffled-order placebo: change follows the timing of warming. Against the other-plot placebo, only boreal pace (any change, ending boreal) stands clearly out; the direction slopes are close to it, so the network responds to regional warming, which cannot be separated from other regional trends of the period (stand ageing, the budworm cycle). Temperate origins show no signal at all.

**Keep in mind.** Within-plot warming is small (median +0.49 SD), so this is a weak test, but a placebo-controlled one.

## F8 Two equilibria and the landscape between them (Discussion)

![](../data-output/final/figures/plan/F8_two_equilibria_draft.png)

**Contains.** Share of boreal, regen and temperate stands along GDD (plots in five bins): observed now (black), and the 1000-year composition implied by the rates without disturbance (red), with the observed disturbance regime (dark green), and with fire and pest only (light green).

**How to read.** The distance between a coloured line and the black one is where the rates would take the landscape under that regime. A direction, not a forecast.

**Shows.** Without disturbance the warm half would become more temperate than now (warmest bin 86% vs 70%). With the observed regime it becomes less temperate (60%). Without harvest it returns close to the undisturbed case (80%). Regen is matched only when disturbance is included. Run separately, both light and heavy harvest are boreal-ward at equilibrium, heavy about twice as strongly.

**Keep in mind.** Rates held at their 1970–2020 values; fire under-counted; the light-harvest part leans on the unaffected-exit assumption for other temperate.

# 6. What we tested but do not show in the main text

| Analysis | What it showed | Where it goes |
|---|---|---|
| Model check (B0) and flicker | the model under-predicts change right after a state change, because new states are often near-ties between two species | SI + limitation |
| Hysteresis diagnostic (I4) | a reclassification rule removes 11–24% of changes; reversals persist | SI; refit pending |
| Rates by period (I2) | overall rates stable; the fir ↔ white birch cycle reverses after 2000 (budworm) | SI + Discussion |
| Contents of pooled classes (I3, I16) | "other temperate" is mostly red spruce and white pine; scoring it changes tolerance results slightly | SI + limitation |
| Plantation check (I5) | 25 plots with unrecorded planting, 1 of 584 regen → black spruce events | limitation sentence |
| Disturbance climates (I6) | each disturbance type has its own climate (fire cold north, light harvest warm south) | SI |
| Home-climate scenarios (I7) | at their own climate, temperate states are much more persistent | Table T2 |
| Equilibrium details (I8, I13) | the regen gap of the undisturbed equilibrium is missing disturbance; pioneers stay under-predicted | SI (summary in F8) |
| Null model by class (I10) | the pace effect of warmth sits in regen exits | in F4 |
| Identifiability (I11) | climate terms and α identified; disturbance main effect trades off with its TSD term (coding, not information) | Methods / SI |
| Fire and pest bounds (I12) | results for temperate starts, white spruce and other pioneer are ranges | SI |
| First passage outside observed climate (I14) | black spruce and jack pine at +1 SD and above are extrapolations | marked in F5 |
| Elasticity (I15) | no headline hinges on one transition | SI table |
| Raw per-transition hazard ratios (A2, A4, A5) | all 90 transitions × all terms | SI tables and heatmaps; as points in F2–F4 |
| Scenario sets C1–C5, D1–D6, E1–E7 | occupancy over time, pathways, recovery, turnover, convergence | mostly SI; the key numbers feed F5, F6, F8 |

# 7. Altogether: what does the chapter show?

**The answer in four sentences.** The forest has barely thermophilized not because its transition rates ignore climate, but because the direction climate gives is already written into where the states sit along the gradient, and the only step warmth accelerates is the closure of regenerating stands. Disturbance decides when stands move, with legacies that reverse over decades, and it does not amplify the climate effect; under the current regime, harvest pulls the warm half back toward boreal. Rates carry memory, so models that assume none mistime recovery. Climate effects learned across space predict the network's response to regional warming, but a plot-specific response is shown only for boreal pace, and not at all for temperate stands.

**New or already known?**

| Finding | Status | Relation to earlier work |
|---|---|---|
| Thermophilization small and lagging; disturbance dominates change | **confirmed** | Brice et al. 2019; Danneyrolles et al. 2019; Talluto et al. 2017 |
| Moderate vs major disturbance act differently (warm-ward vs reset) | **confirmed at finer resolution** | Brice et al. 2019, 2020 (4 states → 13 states, per route) |
| Rates rise with time in state; memory strongest in stand rebuilding | **new** for forest transition models | tests the Markov assumption of Brice et al. 2020 directly |
| Disturbance legacies change sign over decades (reset, then faster recovery) | **new** (time course quantified) | conceptual in Johnstone et al. 2016 |
| Climate direction = spatial sorting; warmth's extra pace only in regen closure | **new** | locates the lag in one phase; complements the recruitment explanations of Sharma et al. 2022 |
| Climate × disturbance antagonistic, not synergistic | **new test**; refines "disturbances amplify" | Brice et al. 2019 measured a faster shift (main effect), not a stronger climate sensitivity |
| Balance point and first-passage threshold along GDD | **partly new** | Brice et al. 2020 reported an equilibrium boundary shift; ours is per route with uncertainty and robustness |
| Space-for-time tested with placebos; transfer by origin group | **new** | answers the untested assumption flagged by Lovell et al. 2023 |
| Recurring harvest regime boreal-ward in the long run | **new empirically** | predicted theoretically by Vieira et al. 2024; consistent with the recent fir rebound (Donnini & Kross 2026, press release) |

**Value beyond this system.** The methods (semi-Markov interval-censored transitions with imputation, a sorting null model, placebo hindcasts) apply to any remeasured plot network. The conceptual result, that the lag sits in the pace of one phase and that disturbance can oppose rather than catalyse climate-driven change, is general enough to test elsewhere. And the space-for-time result names which states climate scenarios built from spatial data can and cannot be trusted for, which matters for every projection model parameterised from inventories, including Chapter 3.

**What the chapter does not do.** It does not forecast. It describes what the observed rates imply under stated conditions; the projections (F5 beyond 55 years, F8) are directions, and their main assumptions are fixed climate, a constant 1970–2020 disturbance regime, and the unaffected-exit assumption for fire and pest.
