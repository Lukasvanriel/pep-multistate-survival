---
title: "Chapter 2 - some results"
date: "13 November 2026"
---

There are 2 parts here: Part 1 is a short overview of how I would build the chapter as I understand it now. Part 2 goes through the findings behind it, one figure/result at a time.  The order in which results are presented here is not necessarily the final order in the chapter but simply the order that came most naturally for this overview.

**The model in two lines.** A multistate survival model on the PEP network: 13 forest states (dominant species or group, plus regen), 90 transitions with at least 15 events. Each transition has a Weibull hazard with GDD, CMI, soil, four disturbance types, GDD × disturbance interactions and time since disturbance (TSD). Change times inside visit intervals are imputed (30 imputations, pooled).

**How the results are structured.** (1) We estimate the hazard rates of the 90 transitions: their memory (II), the effect of disturbance and how it fades (III), and the climate effects and their interaction with disturbance (V, top). (2) We combine the rates, with all exits of a state competing, into what happens to stands over 50–100 years: outcomes after disturbance (V, middle and bottom), which factors change today's forest most, and where boreal–temperate replacement balances (VI). (3) We derive the long-run steady state the rates lead to, with and without disturbance, and compare it with today's forest (VIII). Two checks frame this: whether today's forest is at the steady state of its rates, tested from the rates' side (IV), and whether rates estimated across space hold over time (VII). Section I gives the observed record the model has to explain.

**What is new compared to earlier studies.** Many previous studies are based on simulations (e.g. LANDIS, QLDM: Bouchard et al. 2019), this is based on empirical data. It differs from other empirical studies through several points: the rate can depend on time in the state (Weibull α instead of a Markov rate), climate and disturbance can interact, and disturbance effects can change over time. It is also demographically detailed.

# Part 1. Proposed route/manuscript

**Introduction**

1. Tree communities lag climate with significant expected impacts.
2. Three possible brakes (slow intrinsic rates, disturbance, explicit time dependence for sojourns and disturbances) are rarely tested together.
3. Transition models of forest change are either simulations (e.g. LANDIS-II) or, when empirical, mostly Markov, with climate and disturbance additive and spatial climate effects read as temporal ones.
4. Regeneration as the likely hinge: where disturbance sends stands, where growth responds to warmth, where memory should matter most.
5. Questions: (i) do rates depend on time in state, and where most? (ii) do climate and disturbance act on pace or on direction, and do they interact? (iii) how do disturbance effects change over time?
6. System and approach (pure model output, the hazards, and their implications on the system dynamics).

**Results**

| Section | Finding | Figure |
|---|---|---|
| R1 | The forest changes, but hardly toward warmer types | F1 |
| R2 | Transition rates have memory (time dependence), most in regenerating stands | F2 |
| R3 | Disturbance decides when stands move | F3 |
| R4 | Are the climate effects more than today's distribution requires? | F4 |
| R5 | Disturbance seems to weaken rather than amplifies the climate effect | F5 |
| R6 | On the real plots: disturbance sets the pace; warming and heavy harvest push direction in opposite ways; balance point | F6 |
| R7 | Space-for-time (whether it holds up) | F7 |
| R8 | In the long run, harvest pulls the warm half back toward boreal; without disturbance, toward temperate | F8 |

**Overall:** the forest seems to lag climate mainly because of timing: memory, disturbance and slow regrowth.

Raw hazard ratios: as points in F3, and in full in the SI (tables and matrices).

**Discussion**

To come

**Supplementary information.** Data summary (plots, visits, events per state and disturbance); all hazard ratios; model checks; disturbance regimes/ecologies by climate; sensitivity analyses (e.g. group-level states, soil axis, TSD shape).

# Part 2. Findings

## 0. Hazard ratios
Every raw hazard ratio (90 transitions, all terms, interactions included) is in `figures/hazard_ratios/`, as matrices and as one plot per transition. These raw model output could be presented, but given the amount that represents, I think most will immediately have to go into appendices.

## I. The forest changes, but hardly toward warmer types

**Why.** Before any model: what did the forest actually do, seen through our states?

**What.** Per visit interval, the change per decade in the thermal index of the plot's state (mean GDD where that state occurs, like Brice et al. 2019's community index but at state level; only for within Québec) and in its shade tolerance (Niinemets & Valladares 2006). This only is non-zero when a state change occurred.  Origin groups: boreal (black spruce, white spruce, fir, jack pine), pioneer (white birch, trembling aspen, other pioneer), temperate (sugar maple, red maple, yellow birch, cedar, other temperate) and regen.

![](figures/main_figures/F1_record_draft.png)

**How to read.** Above zero = toward warmer states (left) or more shade-tolerant states (right). Top by origin group and period, bottom by latitude.

**Findings.** Overall thermophilization is close to zero (+0.008 SD per decade), positive only before 1990 and in the south. But stands do change: regen and boreal toward warmer states, pioneer and temperate toward cooler ones, pioneers toward shade tolerance. This matches Brice et al. (2019) at the coarser level of our states.

## II. Transition rates have memory

**Why.** Markov models assume the chance of changing does not depend on time in the state. If that is wrong, they get the timing wrong.

**What.** Each transition has a Weibull hazard, h(t) = α t^(α−1) exp(η), where t is the time the stand has spent in its current state and η = intercept + climate, soil and disturbance effects. The shape α per transition: 1 = Markov (constant rate), > 1 = the rate rises with time in the state, < 1 = it falls. Grouped test: log α of every transition analysed by origin group (shade tolerance of the origin state, Niinemets & Valladares 2006: intolerant, intermediate, tolerant; plus regen), testing whether α differs between groups; prediction was that α would follow a succession profile: highest for intolerant (pioneer) origins, near 1 or even below for tolerant ones. 

Grouping: intolerant (jack pine, white birch, trembling aspen, other pioneer), intermediate (cedar, red maple, yellow birch), tolerant (fir, black spruce, white spruce, sugar maple) and regen; other temperate is excluded because its species span the whole tolerance range.

![](figures/main_figures/F2_memory_draft.png)

![](figures/main_figures/F2b_alpha_groups.png)

**How to read.** First figure: α with 95% CI per transition, grouped by origin. The triangles show the fit hitting the INLA numerical cap. I checked, and the fits where INLA hits its cap do not show significantly different estimates for parameters than those that do not, so this is likely minor. Second figure: left = mean α per shade-tolerance group of the origin (95% CI), right = mean α per origin state; grey = transitions (size ~ precision).

**Findings.** Almost all α > 1, many significantly (48 of 90, none below 1; median 1.45): 'memory' of time in state seems to be the rule, and is strongest for stand rebuilding (out of regen, median 1.73); replacement between established states (red) shows less memory (median α 1.37, 27 of 66 significantly above 1; into regen 1.53, 10 of 12; out of regen 1.73, 11 of 12). Where there are at least 50 events, 77% are significantly above 1, so the limit seems to be in the data, not the pattern. For regen exits, time in state equals time since disturbance, so α and the TSD terms share the signal; this confounding needs sorting out. Grouped test: mean α is 1.83 for regen origins (1.40–2.41), 1.47 for intolerant (1.27–1.72), 1.28 for intermediate (1.07–1.52) and 1.29 for tolerant origins (1.11–1.50); all above 1, but the difference between groups is not significant (p = 0.082). Per origin state, α is highest for regen, trembling aspen and fir (≈1.6–1.8) and close to 1 for sugar maple (1.01, 95% CI 0.87–1.16) and black spruce (1.17, 1.05–1.29).

## III. Disturbance decides when stands move

**Why.** Disturbance is an obvious trigger of stand change. How strong is it, and how long does the effect last?

**What.** Hazard ratio of disturbed vs undisturbed stands at 1, 10, 30 and 60 years after the disturbance, per transition (faded dots/crosses) but averaged (black) per type of move. Raw hazard ratios are in the second figure.

![](figures/main_figures/F3_disturbance_timing_draft.png)

**How to read.** Dots = transitions; crosses = transitions outside the observed range (extrapolations, both sides), black = average. Bottom figure: raw hazard ratio for TSD; below 1 = the effect fades.

**Findings.** Heavy harvest and fire send stands into regen (~143× and 87× in year 1), fading over decades (heavy harvest ~8× at 30 y). This also exists for light harvest interestingly, but not for pest. The light-harvest dip below 1 at 60 years is an artefact of the curve's shape, seemingly not seen in the data. Harvest and pest speed up regrowth out of regen (2–4× at 30 years): first reset, then faster recovery.  Visits are about 10 years apart, so few disturbances are observed within 5 years, and 60 years is beyond the observed range for harvest and pest (extrapolation, x's).

## IV. Is the climate effect in the rates more than today's distribution requires?

This part is still in development

**Why.** The raw GDD hazard ratios (matrix below the figure) look like a strong story: almost everywhere, warmth speeds moves toward warm-adapted states; fir → sugar maple is faster in warm plots, sugar maple → fir in cold plots. The obvious headline would be "warming pushes the forest toward temperate states". But a forest that simply stays as it is along the gradient needs exactly such effects. Brice et al. (2020) checked whether the forest is in steady state by deriving the steady state from the rates and comparing it with today's composition and I do the same in VIII. Here we only test the same thing from the side of the raw fitted rates: are the climate effects in the rates (the gdd coefficients) those needed to keep today's distribution, or do they differ?

**What.**

1. *The balance.* At steady state the flows between two states are equal: (number of A) × (rate A → B) = (number of B) × (rate B → A) (detailed balance; Kelly 1979).
2. *Equal ratios.* Rearranged: rate ratio (A → B / B → A) = number ratio (B / A).
3. *Panel a, direction.* Both the rate and number ratios depend on GDD. If the forest is staying as it is along the gradient, the rate ratio must follow the same GDD dependence as the number ratio. The number ratio's GDD dependence was derived from a multinomial model of the state of every plot visit on GDD, CMI and soil (where each state occurs); the rate ratio's from the two GDD hazard ratios of the survival model. Example: red maple / fir grows ×11.5 per SD of GDD; the rate ratio ×8.7 (fir → red maple ×4.2, red maple → fir ×0.48).
4. *Panel b, turnover.* The balance in panel a fixes only the *ratio* of the two hazard rates. They can also rise or fall together without changing the ratio: if warmth made both fir → red maple and red maple → fir 2× faster, panel A would not notice that. So each pair has a second, independent part: the common level of the two rates (their geometric mean), i.e. how much exchange happens in both directions. panel b shows how warmth changes this overall turnover rate per pair. For fir ↔ red maple: √(4.2 × 0.48) ≈ ×1.4.
5. *Panel c.* Both the ratio and mean from a and b as x and y axes, per pair: whether warmth tilts an exchange (x) and whether it speeds its symmetric turnover (y).
6. *Assumptions.* One-on-one reversibility is stricter than steady state alone. The balance is exact for Markov rates and only approximate for our semi-Markov ones (α ≠ 1); an exact version exists (steady state of a semi-Markov process, π_i ∝ ν_i · m_i; Ross 1996) and we will use it if we keep this framing.

![](figures/main_figures/F4_ABC_mockup.png)

![](figures/hazard_ratios/matrices/01_gdd.png)

**How to read.**

- a: one dot per pair of states with moves both ways (37 pairs; both A → B and B → A were modelled in INLA). x = GDD dependence of the number ratio (warmer / colder state of the pair), y = GDD dependence of the rate ratio (move to warmer / move to colder). On the theoretical dashed line y = x, warming keeps today's distribution; red = warm-ward excess, blue = cold-ward excess.
- b: the same pairs; GDD effect on the average turnover of the two rates; 1 = no change in turnover, above 1 = faster exchange in both directions.
- c: x = distance from the line in panel a (measure of how balanced it is), y = b (turnover, dots in panel b); numbers identify pairs significant on either axis.
- Matrix: GDD hazard ratio per transition; red = faster in warmer plots, blue = slower; bold = 95% CI excludes 1.

**Preliminary results.**

- (a) On average and for most pairs, the pairs follow y = x (slope 1.22, CI 0.90–1.54): the climate tilt of the rates matches what today's distribution requires.
- (a) 8 pairs tilt warm-ward beyond balance: 6 toward a temperate state (white birch → other temperate, white birch → red maple, white birch → yellow birch, yellow birch → other temperate, black spruce → other temperate, fir → cedar), 1 toward white birch (black spruce → white birch) and 1 toward fir (regen → fir). 2 pairs tilt cold-ward: black spruce ↔ fir and black spruce ↔ regen.
- (b) Turnover: warmth speeds exchange mainly in pairs involving regen (×1.9 on average, 95% CI 1.56–2.34); between established states on average ×1.1 (0.97–1.28, CI includes 1), with individual pairs faster or slower.

**Open issues.**

1. The rates from the fits are those of undisturbed stands (main GDD terms); the distribution includes the distribution of all plots, which are obviously shaped by fire, budworm and harvest. Does the warm-ward excess of the white birch pairs reflect a disturbance legacy (birch present because of past disturbance, replaced without new disturbance) or a lag behind recent warming?
2. How to read "on the line": the rates maintain today's spatial pattern and would follow it as places warm. What does that imply for why the forest lags, and how does it relate to the long-run composition in VIII and to Brice et al. (2020)?
3. Turnover: the ×2.8 for regen exits is for regen without a recorded disturbance (about 9% of regen intervals); after disturbance the climate effect on regen exits is roughly halved (V). How should turnover be interpreted alongside disturbance?
4. Is detailed balance the right null? We should probably move to the exact semi-Markov version.
5. Not sure how to interpret the average balance (the a slope being near 1 and the averages seemingly in balance) for most while there are clear variations in individual pairs. 

## V. Disturbance does not amplify the climate effect

**Why.** Disturbance is often thought to amplify climate-driven change (e.g. Brice et al. 2019). The interactions test this directly.

**What.** Top: the GDD × disturbance interaction of each transition (177 estimates), with a weighted average per disturbance type and per type of move (more precise estimates weigh more; transitions sharing a state are treated as related); right: the same against the thermal difference between destination and origin. Middle: the interaction on outcomes over 50 years (chance of leaving, years in the state, thermal direction of the next state), per starting state and disturbance. Bottom: where stands are 50 years after a disturbance, at GDD −1, 0 and +1 SD. Middle and bottom are not fitted directly but calculated from the fitted rates, with uncertainty from parameter draws: the middle by combining all competing exits of each starting state over 50 years, the bottom by projecting the rates forward year by year from a disturbance at year 0 (its effect fading with time since disturbance).

![](figures/main_figures/F5_climate_x_disturbance_draft.png)

**How to read.** Top: 1 = GDD and disturbance effects are additive, below 1 = disturbance weakens the warmth effect. Only the top panel shows the interaction itself. The middle panel (effect of +1 SD GDD when disturbed minus when undisturbed, on the chance of leaving by sojourn year 50, number of years in the state and thermal direction of the next state) also reflects the main effects of disturbance and time in state, because these outcomes are not linear in the rates (e.g. if a disturbance already makes leaving very likely, warming cannot add much). The bottom panel is the full model (main effects, interactions, time since disturbance, α) and with competing states: it shows where stands end up, not which term causes it.

**Findings.**

- *Top (interaction):* On average there is no amplification between warming and disturbance. Leaving regen after a disturbance, the warmth effect is roughly halved (×0.53), and weaker under light harvest (×0.74); elsewhere it is additive (around 1). Against the thermal difference between destination and origin (top, right panel), the interaction shows no trend (slope 0.98, 0.82–1.17): disturbance does not selectively weaken warm-ward moves.
- *Middle (outcomes):* After disturbance warming does less to speed exits from regen and fir (stands stay 5–11 years longer), and for black spruce and white birch under heavy harvest it shifts the next state less warm-ward; trembling aspen under light harvest is the exception (it leaves faster).
- *Bottom (where stands end up):* The temperate share after 50 years rises from GDD −1 to +1 SD from 4 to 48% for fir without disturbance, but only from 2 to 28% after heavy harvest; light harvest raises the temperate share at cold and mean GDD (fir 19% vs 4% at −1 SD) without strengthening the warming effect. This is compatible with Brice et al. (2019, 2020), where moderate disturbances favoured temperate states: here light harvest pushes toward temperate states by itself, not by amplifying climate.
- *Overall:* Disturbance overrides climate rather than amplifying it.

**Caveats.**

- Like the main GDD effects, the interaction is a spatial contrast: across the gradient, disturbed stands differ less between warm and cold places than undisturbed stands do. Reading it as "warming matters less for a disturbed stand that warms" is space-for-time. Disturbed and undisturbed stands also differ in where they occur along GDD (each disturbance type has its own climate range) and, for regen, in what they are (undisturbed regen, about 9% of regen intervals, is a different process from regen after harvest or fire); for regen exits after disturbance, time in regen equals time since disturbance. The interaction may partly reflect these differences rather than a dampening of warmth.
- The averages per disturbance type and per type of move are computed separately, so how much of light harvest's ×0.74 comes from regen exits cannot be told.
- Exits without a fire or pest term (too few events to estimate one) are kept at their undisturbed rate in the outcome calculations. For some temperate starts most of the exit rate comes from such exits, so their outcomes after fire and pest rest on that assumption (triangles in the middle panel).

## VI. Full network: disturbances vs warming; boreal–temperate replacement location

**Why.** Which factors matter most for the forest as it actually is right now, and when considering the full network? Hazards don't necessarilty correlate with probabilities due to competing transitions. Part V shows effects for a typical stand per state; here every real plot is used, with its own state, climate and soil, so warming, moisture, soil and each disturbance are compared on one scale for today's network.

**What.** Four parts. (1, top) For every real plot (current state, own climate and soil, undisturbed), the effect of +1 SD GDD, CMI or soil, or of a disturbance, on the chance of changing state within 50 years and on how much warmer or colder the next state is; averaged over all plots, so all factors share one scale. (2) Along the GDD gradient, for an average stand without disturbance: the chance of reaching another group within 50 years, for all six directions between boreal, pioneer and temperate (e.g. boreal → temperate and temperate → boreal). (3) From (2), the balance between each pair of groups (boreal vs pioneer, boreal vs temperate, pioneer vs temperate): where it is zero, both directions are equally likely (the balance point). (4, bottom) For each non-temperate starting state, without disturbance, the chance of becoming temperate within 100 years along GDD.

![](figures/main_figures/F6_pace_direction_draft.png)

**How to read.** (1) Left = change in the chance of changing state within 50 years; right = change in the thermal index of the next state (above 0 = warmer); each row a starting state, top row all plots. (2) One line per direction between groups. (3) Above 0 = the first group's move to the second dominates; the zero crossing is the balance point. (4) One line per starting state.

**Findings.** Heavy harvest has by far the largest effect on the chance of changing state (+0.35 in the probability of leaving within 50 years); on direction, warming and heavy harvest have equally large but opposite effects (+0.29 and −0.29 in the thermal index of the next state). Without disturbance, boreal ↔ temperate replacement balances at +0.29 SD of GDD (~1,355 degree-days); no single parameter moves this by more than about 11 degree-days. Pioneer ↔ temperate balances near the network mean, and between boreal and pioneer the move to boreal dominates everywhere. Entry into temperate states switches on between GDD 0 and +1 SD for every starting state (least for jack pine). Black spruce and jack pine at +1 SD and beyond are outside the climate where they are observed.

**Next steps.** Redo the balance-point analysis (middle panels) with the observed disturbance regime included (stand simulation as in VIII), as Brice et al. (2020) did with their disturbance scenarios; the current version is undisturbed.

## VII. Space for time test?

**Why.** The time for space assumption could impact what we can claim. Variation in GDD comes mainly from spatial variation, not from local warming. Can we claim anything temporal?

**What.** Based on the fitted model, we predict the change based on the GDD data. (1) For every observed visit interval, the model predicts the chance of each end state twice: with the plot's actual climate (A) and with its climate frozen at its first visit (F); everything else is identical, so A − F is the change the model attributes to that plot's own warming. (2) O is what was observed for that interval (1 = the plot changed state, 0 = it did not; likewise for the thermal shift and for ending in a temperate or boreal state), so O − F is how much more happened than expected without the plot's warming. Within each plot (each plot compared only with itself over time), we regress O − F on A − F. The slope says whether intervals where the model predicts more warming-driven change also show more change: 1 = as much as predicted, 0 = no relation. (3) Two null models redo only steps 1–2 with the warming reassigned, 20 times each: the plot's own warming in random order (right amounts, wrong timing), and the warming of another plot from the same latitude band and period (regional trend kept, link to the plot's own warming broken).

![](figures/main_figures/F7_hindcast_placebo_draft.png)

**How to read.** Black dots = real slopes with 95% CI, per origin group (rows) and kind of change (columns: any change, thermal shift, ending temperate, ending boreal). Bands = range of slopes under the two null models. Above the blue band: the timing of the plot's warming matters. Above the orange band: the plot's own warming matters beyond the regional trend. Inside the orange band: any nearby plot's warming predicts as well, so the signal is regional and cannot be separated from other changes in the region over that period.

**Findings.** Over the record, plots changed about as the spatial climate effects predict (slopes 0.7–0.9 network-wide; no signal for temperate origins), so space-for-time holds as far as we can test it. But except for boreal stands (whether they change and whether they stay boreal), this change cannot be attributed to each plot's own warming rather than to other region-wide changes of the same period. So: using the climate effects to describe how the forest would respond to warming (VI, VIII) is supported and not contradicted by the data, but attributing the change observed during the record to warming is only possible for boreal stands. Warming within plots is small (median +0.49 SD), so this needs more work before deciding.

**Further to explore.**
What does this mean exactly for what we can claim?

Maybe some other tests:

- Latitude band × period fixed effects in the within-plot regression, to absorb the regional trend directly (instead of approximating it with the other-plot).
- Lead (falsification) test: does a plot's *future* warming "predict" its current change? If so, the slopes reflect shared trends rather than a response.
- Matched: another plot's warming only from plots with the same start state and similar time in state, within latitude band × period.
- Within–between decomposition in the model: GDD split into each plot's mean (spatial) and its deviation over time (temporal), as two terms (van de Pol & Wright 2009); a refit, so a sensitivity analysis.

## VIII. In the long run, harvesting pulls the warm half back toward boreal

**Why.** Where do today's rates lead, with and without the current disturbances? This is the calculation of a steady state. We cannot simply use linear algebra as with a Markov chain due to the alpha dependence of the hazard. That's why we have to simulate the dynamics for long times.

**What.** The rates run forward 1000 years in five GDD bins for different scenarios: no disturbance, the observed regime, and fire and pest only. A direction, not a forecast.

![](figures/main_figures/F8_two_equilibria_draft.png)

**How to read.** One panel per group (boreal, pioneer, temperate, regen); x = plots binned by mean GDD, coldest (Q1) to warmest (Q5); black = composition now; coloured = the steady state under each regime.

**Findings.**

- *Temperate:* without disturbance, the rates lead the middle and warm bins to more temperate stands than today (Q3: 34% vs 20%; Q5: 86% vs 70%). With the observed regime, fewer than today (Q4: 35% vs 48%; Q5: 60%). Without harvest (fire and pest only), close to the undisturbed result (Q3: 33%; Q5: 80%). With harvest only (no fire or pest), the same as the full regime (Q3: 20%; Q4: 35%; Q5: 61%).
- *Boreal:* the mirror image: with the observed regime the warm half becomes more boreal than today (Q3: 56% vs 44%; Q5: 24% vs 13%); without disturbance less (Q5: 7%); harvest only as the full regime (Q5: 23%).
- *Regen:* only matched when disturbance is included (Q3: 14% observed and with the regime, 15% with harvest only, 8% without disturbance, 3% without harvest); harvest keeps regen in the landscape.
- *Pioneer:* under-predicted in every version (Q3: 22% observed vs 10–12%); the cause is unknown.
- *Cold bin (Q1):* hardly changes under any regime (boreal 69–79%, temperate ≤ 1%).
- *Overall:* today's forest lies between the undisturbed and the disturbed steady states, and harvest is the main force pulling the warm half toward boreal: harvest alone reproduces the full regime, while fire and pest alone leave the forest close to the undisturbed steady state. Heavy harvest alone pulls more strongly (Q5 temperate 43%) than light harvest alone (62%); together 60%, not additive because a light harvest replaces an earlier heavy one in the overlap rules.
- *Robustness:* Since the disturbances have a certain ecology too (where they occur), and the time since disturance can be modelled with an exponential, I also tested taking that into account. It does not matter much; whether disturbance rates depend on climate changes the shares by at most 2 percentage points; the exponential TSD shape by at most 3.
