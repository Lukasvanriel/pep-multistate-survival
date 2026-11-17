# Agent brief -- Chapter 2 analysis (read this first, every conversation)

You are helping Lukas (PhD, Université de Montréal) carry out the analysis plan of Chapter 2: a multistate Weibull survival model of forest state transitions in Québec's PEP plot network. Lukas supervises; you work on ONE plan section per conversation, row by row, with him.

## 1. Read order, actually read it, do not be lazy (before doing anything)
1. This brief.
2. `docs/chapter2_analysis_plan.md` -- the frozen plan: conventions, condensed list, full list (quantity | theory | question | source | status), exclusions.
3. `docs/chapter2_analysis_log.md` -- what earlier conversations did and found.
4. Old scripts and outputs ONLY as needed for the row you are on.

## 2. Scope
- Work only on the section Lukas names. Do not start the next section.
- Do not add, drop or reword plan rows. New ideas go to the "Later" table of the plan, and only after Lukas agrees.
- Do not rebuild or re-derive the plan. It is frozen.

## 3. Working rules (non-negotiable)
- Diagnose before fixing: inspect the actual data/state before proposing a change.
- No code without Lukas's explicit permission.
- Never make a methodological or data-coding decision yourself. Put every choice to Lukas as an explicit question with an example, then wait.
- Pick one clear position and hold it. Do not flip between recommendations.
- Concise answers; complete tables over narration.
- Caveats proportionate: flag only problems that change what can be claimed, and say which kind they are.
- Never present a list or a check as complete unless you have checked it against the plan row by row.
- Do not invent labels, codes or structures that do not exist in the project. Refer only to files you have opened.
- Say plainly what you have and have not verified ("the script header says X" is not "X is correct").
- References: primary sources only.
- Figures: clean and uncluttered; no technical clutter in labels (no "per +1 SD" when the scale is obvious; say "log scale" once in the axis label).

## 4. Model facts
- 13 states (dominant species/group by basal area, plus regen); 90 modelled transitions (≥ 15 events).
- Hazard per transition: Weibull, variant 0: S(t) = exp(−e^η t^α); β = log hazard ratio. **Fits are on a decade time scale**; the existing engine converts the intercept to years as b0y = intercept − α·log(10) (see `R/17_figures_C.R`, parameter draws).
- Covariates: GDD (gdd_3), CMI (cmi_gs), soil PC1 (all z-scored; scaling in `data-output/final/fits/scaling.rds`); four disturbances (fire, pest, light harvest, heavy harvest); GDD × disturbance; time since disturbance (TSD) as log(TSD + 1).
- Change times inside visit intervals: 30 multiple imputations, Rubin pooling.
- Uncertainty: posterior draws from the pooled estimates with Rubin's total covariance T = W + (1 + 1/M)·B (as in `R/17_figures_C.R`). Several old outputs are point estimates only.
- Paths: `R/config.R` (RUN_TAG = "final"; directories). Key inputs:
  - `data-output/final/results/pooled_final_final.csv` (pooled estimates), `pooled_final_final_cov.rds` (within-imputation covariance), `fit_validity_final.csv`, `hazard_ratios_final.csv`, `alpha_report_final.csv`, `intercept_years_final.csv`
  - `data-output/final/fits/mi_raw_final.csv`, `mi_raw_redraw_final.csv` (per-imputation estimates), `scaling.rds`
  - `data/tables/final/visit_intervals.rds`, `plot_info.rds`, `transition_set.csv`

## 5. Settled decisions
- All conventions: see the plan's "Conventions" table (scenarios, horizons 10/25/50 y + steady state, convergence = halving time, resolution, real plots vs scenario stand, uncertainty, extrapolation).
- Groupings: ONLY from `R/00_groups.R` (`GRP` thermal; `TOL_GROUP` shade tolerance). Never redefine them in a script. Type of move (into regen / out of regen / between states) is a transition property.
- Shade tolerance: other_temperate has NO tolerance group or score; other_pioneer score = 1.18 (Niinemets & Valladares 2006, plot-weighted).

## 6. Folder rules
- Old `R/` scripts and `data-output/final/` (outside `analysis/`) are READ-ONLY references.
- New scripts: `R/analysis/`, named `<section>_<nn>_<short_name>.R` (e.g. `2_01_baseline_rates.R`).
- One script per plan row (nn = the row's number within its section, e.g. `0_01_record_contents.R`). Where rows share a calculation, the shared part goes in the engine file and each row's script only calls it and computes its own quantity. Never merge rows into one script.
- New outputs: `data-output/final/analysis/<section folder>/` (0_data, 1_assumptions, 2_hazard, 3_one_stay, 4_process, 5_long_run, 6_space_for_time, 7_robustness).
- Shared code (draws, hazard, competing-risks, simulation, steady state): one engine file in `R/analysis/`, grown as needed. Borrow code from old scripts after checking it; copy it in properly.
- NEVER run parts of another script by line number (`eval(parse(text = src[a:b]))`); this broke silently once already.
- Before editing any existing file, copy it to `_to_delete/` next to it with a dated name. Deleting files needs Lukas's permission.

## 7. Per-row routine
1. Read the row (quantity, theory, question, conventions that apply).
2. Look for an old script/output that computes it; check it against the row's theory and the conventions. Report: correct / partly / wrong / absent, and why.
3. Explain concisely what that row is about
4. Ask Lukas for the row's open choices (grouping, climate values, real plots vs scenario stand) with an very concise example. Then, with permission: reuse or write the code, run it into the new folder.
5. Fill the row's Status in the plan (output path) and add a log entry (format in the log file).

## 8. Environment
- The shell (`device_bash`) runs in a VM; the project is at `$HOME/mnt/pep-multistate-survival`. Each call is limited to 180 s; background processes are killed when the call ends, so long or resumable scripts are run in foreground chunks.
- **R is not part of the VM. Each conversation installs its own R and removes it again.** The VM's session disk (~9.8 GB) is shared by all conversations and nothing cleans it up, so an install left behind blocks every later conversation. R cannot be installed in the project folder (tested 17 Nov: the Mac file system breaks the install).
- **Before installing:** check `df -h /sessions`. You need about 2 GB free. If there is less, stop and tell Lukas; do not try to free space outside your own session folder.
- **Install** (into your own session folder, never the project folder), then clean the download cache straight away:
  ```
  mkdir -p $HOME/mm && cd $HOME/mm
  curl -sL https://conda.anaconda.org/conda-forge/linux-aarch64/micromamba-2.9.0-0.tar.bz2 -o mm.tar.bz2 && tar -xjf mm.tar.bz2 bin/micromamba && rm mm.tar.bz2
  export MAMBA_ROOT_PREFIX=$HOME/mm/root
  ./bin/micromamba create -y -n r -c conda-forge r-base=4.4 r-data.table r-ggplot2 r-here r-dbi r-rsqlite r-metafor r-scales r-patchwork r-nnet
  ./bin/micromamba clean -a -y
  # Rscript: $HOME/mm/root/envs/r/bin/Rscript ; run from the project root
  ```
  (Check `uname -m` first; use linux-64 instead of linux-aarch64 on x86. The create step can take close to 180 s; if it times out, re-run it.)
- **Remove R immediately after the section commit** (doer), or at the end of the check (checker): `rm -rf $HOME/mm`, then show `df -h /sessions`. Also remove anything else you created in your session folder (scratch files). Deleting needs Lukas's permission (delete-permission tool).
- **Never two installs at the same time:** the doer removes R before the checker starts.

## 9. Known pitfalls
- Caches: resumable scripts reuse caches; when an input changes, the cache can hold stale results (happened with `figures/E/cache/E1_E2_first_move.rds`). Check caches when inputs change.
- Some old outputs are point estimates only (no draws).
- Old scripts used fixed z = 0 scenarios and a 100-y horizon; these do not match the frozen conventions.

## 10. Commits
- One commit per finished plan section, plus one before changing any existing file.
- Ask Lukas to allow deletion first (git needs to delete its lock files), via the delete-permission tool for the project folder.
