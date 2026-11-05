# Sensitivity plan (2026-11-05)

Status: **proposal for later review. Nothing here is decided or implemented** unless marked. Written after the final run
(RUN_TAG "final": 05 done, 07 with the z_map screen running, 08-12 pending). Definitions come from `docs/PIPELINE_NOTES.md`
(sections "B. Sensitivity analyses", "Sensitivity analyses", T3/T4 decisions, 2026-10-26 S8 redefinition) and
`docs/PIPELINE_FIX_PLAN_2026-11-03.md`. The ranking reflects the state of the pipeline on 2026-11-05.

## Types

| Type | Meaning | Cost |
|---|---|---|
| post-hoc | re-pooling of existing fits, no refit (`R/sensitivity/10_sensitivity_posthoc.R`) | minutes |
| refit | refit on the stored imputations (`extras_final_mXX.rds`: P1 entries, P2 exits, theta draws), with changed settings, formula or plot subset; same screen, same pooling | minutes to hours (Mac, INLA) |
| full run | own processed tables and own output tag, then 05 -> 07 -> 08 (data or covariate definitions change upstream of 05) | about 6-8 h of Mac time each (from group-C/WP6 timings) |

## Tier 1: High

| S | What | Type | Why high |
|---|---|---|---|
| S1 | Interactions only where the disturbance level has >= 20 disturbed events | post-hoc (script 10, ready) | Addresses the harvest levels with 1-14 events that get main + GDDx + TSD terms under A6 |
| S2 | Pooling without fits whose alpha-SD ratio < 0.5 | post-hoc (script 10, ready) | Tests the 2026-11-05 screen, which still allows fits with ratio 0.1-0.5 |
| S5 | Within- vs between-plot climate (Mundlak decomposition) | refit (script 11, ready) | Core support for the space-for-time framing of the climate effects |
| S7 | Exclude P1-fallback plots (entry W_j/2) | refit (plot subset) | Fallbacks are about 3% of plots per imputation (295-453 of 11,286), not about 1% as originally assumed; concentrated in sugar maple (10-18%) and cedar |
| S9 | Regional random effect f(reg_eco, "iid") | refit (formula) | Spatial confounding; SEs possibly 5-15% too small (design-effect guess). **Catch:** the z_map screen uses fit_map_v0, which has no random effect, so S9 needs its own validity check. Flag transitions with < 5 regions holding events. |
| S21 | Coefficient prior N(0, 2.5^2) instead of N(0, 5^2) | refit (settings) | The main prior was a deliberate choice; reviewers will ask. Same script as S22 and S4. |
| S22 | INLA default coefficient prior (precision 0.001, sigma 31.6) | refit (settings) | As S21. May converge worse in X->regen (near-separation by fire/heavy harvest, hazard ratios about 70-300). |
| S24 | TSD linear in decades instead of log(TSD + 1) | refit (TSD_FORM = "linear") | A key a-priori choice (2026-11-01); disturbance effects at TSD = 0 depend on it. Note: P1's pre-visit TSD also uses TSD_FORM, so a refit on stored imputations is an approximation (P1 not redrawn). |
| **S25 (new)** | X->regen exit rule off (`XREGEN_EXIT_RULE = FALSE`): X->regen exits drawn like all others, not pinned to the disturbance date | refit (rebuild P2 with the switch off from the stored theta draws and seeds) | **Not yet in the plan.** Tests the open question on alpha: sugar_maple->regen alpha 0.81 (main run) -> 1.19 (final; 95% CI 0.93-1.54), 62 of its 78 exits are rule exits; and the rule's effect on all X->regen estimates (rule exits: black_spruce->regen 81%, fir->regen 58%, jack_pine->regen 87%). |

## Tier 2: Medium

| S | What | Type | Note |
|---|---|---|---|
| S3 | All valid fits instead of the first 30 by m | post-hoc (script 10, ready) | Cheap; tests the selection by index |
| S4 | Alpha prior pc.alphaw(5) (INLA default) instead of 0.1 | refit (settings) | Alpha is secondary, but its prior was debated (2026-10-23, 2026-11-01). In the S21/S22 script. |
| S6 | One plot per plot pair | refit (plot subset) | Spatial non-independence |
| S8 | Exclude first-state-regen plots with W_j >= 70 (277 plots on the cloud build) | refit (plot subset) | Redefined by the user 2026-10-26 (replaces "keep old regen / W_j >= 30"); one plot alone moved regen->fir alpha/pest terms by up to 0.7 SE |
| S16 | Without the 177 imputed-soil plots | refit (plot subset) | Soil effects are small |
| S18b | Drop plots with undated partial codes at the first visit | refit (plot subset) | Tests the assumed first-visit ages (harvest 4.5 y, pest 32 y median of dated codes); e.g. 69 of the 72 pest cases among the 121 "flag off" plots rest on the assumed pest age |
| S19 | Assumed first-visit ages at the 10% / 90% quantiles (harvest 3 / 22.5 y, pest 25 / 37 y) | full run (03 ->) | Same assumption as S18b, tested directly |
| S10 | Disturbance flags/TSD at the start visit instead of the end visit | full run | Timing convention (A10 decided end visit); the plan notes it is "relevant again now that rows are not cut" |
| S11 | TSD upper bound (`TSD_SOURCE = imputed_prev`) | full run | Dating assumption for undated later codes (lower-bound dating in the main run) |
| S12 | Hysteresis: state change only if the new dominant leads by >= 10% of BA | full run (from 02) | Classification noise; the user called a hysteresis reclassification acceptable as a sensitivity (2026-10-21) |
| S13 | W_j age-class mapping alternatives (VIN/VIR, JIN/JIR values; +/- one class) | full run | W_j drives P1 |

## Tier 3: Low

| S | What | Type | Note |
|---|---|---|---|
| S14 | Stage-1 entry fraction W/4 and 3W/4 instead of W/2 | full run | Stage 1 only seeds Stage 2; the simulation showed an effect on alpha, which is secondary |
| S15 | Regeneration BA threshold 3-4 m2/ha | full run | Already marked lower priority |
| S17 | Event threshold 20 instead of 15 | refit (fewer terms) | Largely covered by S1 |
| S20 | pee_ori_sond codes only (no station_pe codes) | full run | Essentially the earlier build (`USE_STATION_CODES = FALSE` reproduces it exactly) |
| S23 | Event limit 10 instead of 15 (+17 transitions) | full run | Added transitions are sparse and likely to fail to converge often |
| S26 (new) | Cap-free MAP (fit_map_v0) pooled over all imputations, for the X->regen transitions | ready (`R/tests/24_map_selection_regen.R`) | Already computed for black_spruce, fir, jack_pine and trembling_aspen->regen (log `data-output/diagnostics/24_map_selection.log`); a supplement table if wanted |

## Obsolete

| S | Why |
|---|---|
| S18a | "Flag undated first-visit partial codes with TSD = 0" tested the old main rule (such codes not flagged). The rule changed: they are now flagged with an assumed first-visit age (median of dated codes). S18b and S19 replace S18a. |

## Implementation notes (for when this is taken up)

- **One generic refit script** covers S4, S6, S7, S8, S16, S18b, S21, S22, S24, S25 (and S17): rebuild the rows of each used imputation
  from the stored extras (as script 11 and 08f do, with their identity check), apply the change (settings / formula / plot subset /
  P2 switch), refit with the same INLA call, apply the same screen (no error, alpha-SD ratio >= 0.1, |z_map_alpha| and
  |z_map_intercept| <= 3), replace invalid fits as 07 does, pool as 08 does. Output under its own tag, e.g. `data-output/final_S21/`.
- **S9** needs a separate validity check (no MAP counterpart with a random effect).
- **S24 and S25** change inputs of P1/P2; as refits on stored imputations they redraw only what the switch touches. State this as an
  approximation, or do them as full runs if the difference matters.
- **Full runs** (S10-S13, S19 and the low tier) need their own processed tables (04 overwrites `data/tables/<tag>`), their own RUN_TAG
  and a clean fits folder. They could run on the Mac mini in parallel, once it has a project copy, INLA and the tables.
- After every sensitivity: compare with the main results on the hazard-ratio scale (09 contrasts incl. the TSD grid), not on raw
  coefficients.

## Open decisions (user)

1. Agree with the tiers?
2. Add S25 (exit rule off) and S26 (MAP over all imputations) to the plan?
3. Which tiers to do, and when (e.g. High now, Medium after the paper outline, Low only on reviewer request)?
