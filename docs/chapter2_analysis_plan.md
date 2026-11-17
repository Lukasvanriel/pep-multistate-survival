---
title: "Chapter 2 -- analysis plan (frozen list)"
date: "17 November 2026"
---

# Purpose and rules

This is the reference list of analyses for Chapter 2 (multistate Weibull survival model, 13 states, PEP network). It was built from the model's theory, then checked against Brice et al. (2020, *Global Change Biology*) and Putter, Fiocco & Geskus (2007, *Statistics in Medicine*). Work from this file; do not rebuild the list in a new conversation.

- **Completeness** = (a) every quantity the two sources report, plus (b) every new feature of our model (α ≠ 1, time since disturbance, GDD × disturbance, 13 states) followed through every step. Items from either source that we do not do are listed under "Excluded", with the reason.
- **Frozen**: new ideas go to "Later", not into the plan, unless the list is formally reopened.
- **The chain**: hazards (step 1) → one stay with competing exits (step 2) → several changes, i.e. the process (step 3) → the limit (step 4). Each step is a function of the previous one.
- **Status column**: to be filled row by row (exists: output code; missing; partly). Grouping, climate values and use of real plots are decided per row when we reach it.

Notation: h_ij(t | x) = hazard from state i to j at time in state t with covariates x; H_ij(t) = cumulative hazard.

# Conventions (all settled)

| Convention | Decision |
|---|---|
| Starting state | each of the 13 states; the observed mix of states at first visit for whole-forest statements |
| Time in state at start | 0; the model's mean time in state per state as a check (memory) |
| Climate (GDD) | −1 / 0 / +1 SD for comparisons; full gradient where we look for where something changes (decided per row) |
| CMI, soil | 0 (network mean) in scenarios; check at ±1 SD later (lower priority) |
| Disturbance | 10–50 y results: no disturbance, or one event of each type, at the start and after the mean time in state; long run: recurring regimes (observed; no harvest; harvest only; fire and pest only; rates depending on climate) |
| Horizons | 10, 25, 50 y; steady state. 100 y dropped (extrapolation; the convergence time gives the timescale between 50 y and the steady state) |
| Convergence time | time for the distance between the composition and the steady state to halve |
| Resolution | compute per state (or per transition); report per state, thermal group, shade-tolerance group, type of move (type of move only at steps 1–2); chosen per row |
| Real plots vs scenario stand | conditional effects (scenario stand) for mechanism; marginal effects (real plots) for statements about today's forest; chosen per row |
| Uncertainty | posterior draws, including between-imputation variance |
| Extrapolation | results beyond the record are reported as properties of the fitted process (comparisons between scenarios), not forecasts |

# Condensed list (frozen)

*Items added on 17 Nov after checking the sources are marked (new).*

| Section | Quantities |
|---|---|
| Data | record contents; observed change between visits; disturbance regime; raw rates over calendar time; climate trend over the record (new) |
| Model assumptions | fit validity; Weibull shape; proportional hazards; imputation; observed vs predicted; stated untestable assumptions; model comparison by covariate group (new); choice of clock (new) |
| Step 1: hazard | baseline rate; ratio of reciprocal baseline rates (new); α and test of α = 1; GDD, CMI, soil effects; do effects differ between transitions (new); disturbance effects and decay with TSD; GDD × disturbance; α vs TSD in regen exits; previous-state check (new) |
| Step 2: one stay | persistence at 10/25/50 y; mean time in state; chance of leaving to each destination; destination shares; entropy of destinations; pace vs direction split; each driver's effect on these; interaction on this scale |
| Step 3: several changes | composition at 10/25/50 y; composition after a disturbance event; first passage into a group; pathway probabilities (new); balance point along GDD; number of moves and returns; tolerance shift of moves |
| Step 4: long run (along the GDD gradient, per disturbance regime) | steady state; boreal–temperate boundary and its shift under disturbance (new); system turnover (new); system entropy (new); convergence time; steady state vs today's forest; analytic vs simulated steady state |
| Space-for-time | within-plot predicted vs observed change |
| Robustness (headline results only) | Markov refit (α = 1); TSD shape; exits without a disturbance term; CMI and soil at ±1 SD; one transition perturbed at a time; draws include imputation variance; state definition (already in the sensitivity plan) |
| Later | demographic decomposition (subset of transitions); warming over time |

# Full list: theory and question

## Data (no model)

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Plots, intervals, events per transition | interval-censored observation of a multistate process | What does the record contain? | both | |
| Observed change between visits (thermal index, shade tolerance) | empirical counterpart of the occupancy P(t) over one interval | What did the forest actually do? | — | |
| Disturbance frequency, co-occurrence, climate position | distribution of the time-varying covariates | What disturbance regime is behind the data? | Brice 2.3 | |
| Raw transition rates over calendar time | crude hazard estimate (events / time at risk; constant-hazard estimator) | Is the process time-homogeneous over the record? | — | |
| Climate trend in the plots over the record (new) | observed covariate path within plots | How much did climate change during the record? | Brice 2.3 | |

## Model assumptions

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Convergence, failed fits | estimation | Are the fits valid? | — | |
| Fitted vs observed time in state | parametric family: Weibull | Does the Weibull shape fit? | — | |
| Covariate effects over time in state | proportional hazards: a covariate multiplies h by a constant factor | Do effects stay constant over time in state? | Putter 4.4 | |
| Imputation diagnostics | multiple imputation of interval-censored change times; Rubin pooling | Is imputing change times valid? | — | |
| Observed vs predicted transitions | goodness of fit (Putter: model-predicted vs non-parametric cumulative incidence) | Does the model reproduce the data? | Putter 3.3 | |
| Model comparison by covariate group (new) | nested-model comparison (Brice: likelihood ratio and AIC; with INLA: WAIC or LOO) for baseline, climate, soil/CMI, disturbance, full | Does each covariate group improve the model? | Brice 2.4.2 | |
| Choice of clock (new) | clock reset (time in state) vs clock forward (time since an origin); TSD is the model's second clock | Why time in state is the right clock | Putter 4.2, 4.4 | |
| Stated, not testable | competing exits independent given covariates; censoring non-informative | — | Putter 2, 3.2 | |

## Step 1: the hazard, h_ij(t | x) = α_ij t^(α_ij − 1) exp(η_ij(x))

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Baseline rate | intercept of η_ij: hazard at reference covariates | How fast does each transition run by itself? | Brice 3.1 | |
| Ratio of reciprocal baseline rates (new) | q_ij / q_ji at reference covariates | Which direction dominates between two states by itself? | Brice 3.1 | |
| α, test of α = 1 | h ∝ t^(α − 1); α = 1 gives a constant hazard (Markov), nested in the Weibull | Do rates depend on time in state (memory)? | new feature | |
| GDD, CMI, soil effects | hazard ratio exp(β) per SD | How do climate and site change the rates? | Brice 3.2 | |
| Do effects differ between transitions (new) | heterogeneity of β across transitions (test of a common effect) | Is, e.g., the GDD effect the same for all moves, or does it differ by type of move? | Putter 3.3 | |
| Disturbance effects and decay with TSD | time-varying covariate: disturbance effect plus a function of time since disturbance | How does disturbance change the rates, and how fast does that fade? | Brice 3.2; new feature (decay) | |
| GDD × disturbance | β_GDD differs between disturbed and undisturbed stands | Do climate and disturbance interact? | new feature | |
| α vs TSD in regen exits | where time in state = TSD, α and the TSD terms are functions of the same clock (identifiability) | Can memory and recovery be separated? | new feature | |
| Previous-state check (new) | semi-Markov property: given the current state and time in it, the past does not matter. Check: predicted vs observed destinations split by the state before the current one (no refit); single-change intervals only; regen first. Refits with history (previous state as covariate, or split states) are out of scope for this chapter | Does history beyond the current state matter? Supports the memory claims on α | Putter 4.2 | |

## Step 2: one stay in a state (competing risks; built from all h_ij out of state i)

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Chance of still being in the state at 10/25/50 y | S_i(t) = exp(−Σ_j H_ij(t)) | How persistent is each state? | Putter 3.2 | |
| Mean time in state | m_i = ∫ S_i(t) dt (Brice's turnover time −1/q_rr is the Markov case) | Pace | Brice 2.4.4 | |
| Chance of leaving to each destination by 10/25/50 y | cumulative incidence CIF_ij(t) = ∫_0^t S_i(u) h_ij(u) du | Where do stands go, and when? | Putter 3.2 | |
| Destination shares | jump probabilities p_ij = CIF_ij(∞) (the embedded chain) | Direction | Brice 2.4.4 | |
| Entropy of destinations | −Σ_j p_ij log p_ij | How predictable is the next state? | Brice 2.4.4 | |
| Pace vs direction split | h_ij = (total exit hazard h_i·) × (share h_ij / h_i·); a driver's effect split between the two factors | Does a driver act on pace or on direction? | new feature | |
| Each driver's effect on these quantities | change in S_i, m_i, CIF, p when x changes; non-linear in the hazards | Which drivers matter once exits compete? | Putter 3.3 | |
| Interaction on this scale | the same, for GDD × disturbance | Does the interaction matter for outcomes? | new feature | |

## Step 3: several changes (the process; semi-Markov kernel Q_ij(t) = CIF_ij(t))

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Composition at 10/25/50 y | occupancy P_ij(t) from the Markov renewal equation P_ij(t) = δ_ij S_i(t) + Σ_k ∫_0^t dQ_ik(u) P_kj(t − u). With TSD the process has two clocks, so it is computed by simulation | Where is a stand after several moves? | Putter 4.5 | |
| Composition after a disturbance event | the same, disturbance at the start (or after the mean time in state), TSD advancing | Resistance, recovery, reorganisation | new feature | |
| First passage into a group | renewal equation with the target states absorbing | Whether and when, e.g., temperate is reached | — | |
| Pathway probabilities (new) | probability of passing through state k before reaching state j | Which route do stands take (e.g. through regen or pioneer before temperate)? | Putter 4.5 | |
| Balance point along GDD | GDD where first passage A → B equals B → A within a horizon | Where does replacement between two groups balance? | — | |
| Number of moves, returns | expected number of jumps by t (Markov renewal function); return = first passage back to the start state | How dynamic and how cyclic are stands? | — | |
| Tolerance shift of moves | expected jump counts × difference in shade-tolerance score | Which way does succession run? | — | |

## Step 4: the long run (each quantity along the full GDD gradient, per disturbance regime)

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Steady state | undisturbed: π_i = ν_i m_i / Σ_k ν_k m_k, ν = stationary distribution of the embedded chain (ν = ν p); under a recurring regime: second clock, computed by simulation | Where do the rates lead? | Brice Fig. 5 | |
| Boreal–temperate boundary and its shift (new) | GDD where π_boreal = π_temperate, per regime | Where is the front, and does disturbance move it? | Brice Fig. 5 | |
| System turnover (new) | Σ_i π_i m_i | How fast are the dynamics along the gradient? | Brice Fig. 6a,b | |
| System entropy (new) | Σ_i π_i H_i | Where are the dynamics most uncertain? | Brice Fig. 6c,d | |
| Convergence time | time for the distance between P(t) and π to halve (Markov special case: log 2 / spectral gap, Brice's half-life) | How long does it take to get there (lag)? | Brice Fig. 6e,f | |
| Steady state vs today's forest | distance between π and the observed composition; flow balance Σ_i (π_i / m_i) p_ij = π_j / m_j | Is today's forest at the steady state of its rates? | Brice Fig. 5a | |
| Analytic vs simulated steady state | closed form above vs simulation (undisturbed) | Is the simulation right? | — | |

## Space-for-time

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Within-plot predicted vs observed change | β identified mostly between plots; reading it as a response over time assumes between-plot = within-plot effects | Can climate effects estimated across space be read as effects over time? | — | |

## Robustness (headline results only)

| Quantity | Theory | Question | Source | Status |
|---|---|---|---|---|
| Markov refit (α = 1) | nested model, all parameters re-estimated with α fixed at 1 | Does memory change the conclusions? | new feature | |
| TSD shape | alternative functional form for the TSD effect (log vs exponential) | Do conclusions depend on how disturbance effects fade? | new feature | |
| Exits without a disturbance term | bounds on the unestimated effects | Do conclusions depend on that assumption? | — | |
| CMI and soil at ±1 SD | step 2–4 outcomes are non-linear in the hazards, so the GDD contrast can depend on the reference values | Do GDD conclusions hold at other site conditions? | — | |
| One transition perturbed at a time | sensitivity (elasticity) of outcomes to each h_ij | Which transitions drive the results? | Putter 3.3 | |
| Draws include imputation variance | Rubin total variance T = W + (1 + 1/M) B | Is uncertainty fully propagated? | — | |
| State definition | results conditional on the state classification | Do conclusions depend on how states are defined? (already in the sensitivity plan) | Brice 2.2 | |

## Later

| Quantity | Theory | Question | Source |
|---|---|---|---|
| Demographic decomposition (subset of transitions) | basal-area change of each species split into growth, recruitment, mortality, harvest | How do transitions happen? | Brice 2.4.5, 3.5 |
| Warming over time | covariate path x(t) instead of fixed x; relies on space-for-time | Does a warming trajectory change the conclusions? | — |

# Excluded, with reason

| Item | Source | Reason |
|---|---|---|
| Shared or proportional baselines across transitions; reduced-rank models | Putter 3.3, 4.4 | modelling shortcuts for few events per transition; our transitions are fitted separately with their own baselines |
| Fine–Gray subdistribution model | Putter 3.3 | gives effects directly on cumulative incidence, which step 2 already derives from the cause-specific hazards |
| Previous state as covariate; states split by history | Putter 4.2 | require refits with a new state set; out of scope for this chapter (the predicted-vs-observed check stays) |
| Disturbance held permanently "on" as a scenario | Brice 2.4.4 | Markov-specific; we use recurring regimes with fading effects instead (state the difference in the methods: long-run numbers are not directly comparable) |
| 100-year horizon | — | extrapolation; adds nothing beyond the steady state and convergence time |

# Next steps

1. Confirm this list; take it to the supervisor for sign-off.
2. Fill the Status column row by row against the existing outputs.
3. Per row: decide grouping, climate values (−1/0/+1 or gradient) and real plots vs scenario stand.
4. Run what is missing, in an agreed order.

# References

- Brice, M.-H., Vissault, S., Vieira, W., Gravel, D., Legendre, P. & Fortin, M.-J. (2020). Moderate disturbances accelerate forest transition dynamics under climate change in the temperate–boreal ecotone of eastern North America. *Global Change Biology*, 26, 4418–4435. https://doi.org/10.1111/gcb.15143
- Putter, H., Fiocco, M. & Geskus, R. B. (2007). Tutorial in biostatistics: competing risks and multi-state models. *Statistics in Medicine*, 26, 2389–2430. https://doi.org/10.1002/sim.2712
