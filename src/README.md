# HPAI H5N1 within-herd agent-based model: revised pipeline

## Setup

1. Install the packages in `requirements.txt`.
2. Put `Cases_by_date-4Pranav.xlsx` in this folder.
3. Copy `farm_labels.example.json` to `farm_labels.local.json` and replace the placeholder codes with the herd codes used in the workbook. This file and the workbook are listed in `.gitignore`, so the herd codes never reach GitHub.

## Run order

Open and run the four notebooks in order. All settings are in `config.py`; restart the kernel after editing it.

1. `01_prescreen_optimize.ipynb`: burn-in herds, pre-screen (rejection ABC), Bayesian optimisation
2. `02_post_simulations.ipynb`: `N_SIMULATIONS` runs per farm and metric at the optimum
3. `03_sensitivity_sobol.ipynb`: Sobol global sensitivity
4. `04_sensitivity_local.ipynb`: one-at-a-time local sensitivity around each farm's optimum

Steps 2–4 read the optima from step 1. Steps 3 and 4 do not depend on each other.

## Files

- `config.py`: every setting (iterations, simulations, bounds, pens, seeds).
- `ab_sim.py`: the model, loss functions, burn-in, data input and derived quantities.
- `common.py`: shared helpers (burn-in loading, parallel map, outcome summaries).
- `workers.py`: the four functions that run in parallel worker processes. They must live in a `.py` file: on macOS, Windows and Linux with Python 3.14+, worker processes cannot see functions defined in a notebook.
- `farm_labels.example.json`: template for `farm_labels.local.json`.

## Burn-in

One burn-in is run at `BURN_IN_REFERENCE_SIZE` (10 years of demography, no disease). Every farm herd and every `HS_GRID` herd is a random subsample of it with exactly the target number of lactating cows; other pens are thinned by the same fraction. All burn-ins are cached in `results/burnin/`. Delete that folder after changing anything that affects the herd (pens, culling, herd sizes, reference size).

## What changed from the previous code

| Previous behaviour | Now |
| --- | --- |
| All index cows were new cows added to Pen 8 | Index cows are existing susceptible cows drawn at random from `SEED_PENS`, so herd size is unchanged |
| Fitted incidence counted new infectious cows in every pen, plus index cows on day 1 | Counts only cows becoming infectious in `CASE_PENS`; index cows are excluded unless `COUNT_INDEX_AS_DAY1_CASES = True` |
| Notebook fallback made every farm use the last farm's herd and parameters | Each farm has its own cached burn-in herd; no global state |
| The 305-day rule grouped pens 6 and 7 with the wrong lactation group and summed all past lactations | Uses days in milk of the current lactation |
| Herd size emerged from the burn-in | Each herd is thinned from one reference burn-in to exactly the lactating herd size in sheet 2 |
| One stochastic simulation per loss evaluation | Loss computed on the mean of `N_REPS_PER_EVAL` simulations; best candidates re-evaluated with `OPT_REEVAL_REPS` |
| Pre-screen PLF had its arguments reversed | Fixed |
| βeff from the net change in exposed cows | βH(t) = new exposures × N / (S × I), plus a direct reproduction number from secondary cases |
| Sensitivity: one run per set, T 0.01–0.99, herd size by sub-sampling Farm 1 | Replicates averaged, ranges from `BOUNDS`, herd size from separately built herds, bootstrap CIs, convergence check |
| Hospital pen cows counted as milking only in the figure | Pen 15 is the 2-day post-calving pen for every calving cow, as in the code |

## Decisions to confirm

1. **Transmission bounds.** T is now searched on a log scale between 0.0002 and 0.02. With index cows spread across milking pens, T = 0.01 already infects nearly every lactating cow, and earlier fits sat at the old lower bound (0.01).
2. **Missing dates.** With `MISSING_DAYS_FILL = "zero"` (default), days without a row between a farm's first and last case are counted as zero cases; `"carry_forward"` repeats the last recorded day's count instead, which adds 30 cases for Farm 1, 51 for Farm 2 and 13 for Farm 4. Farm 2 (T) has 33 such days and its window runs from 17 Sep to 10 Nov (55 days). Farm 1 (W) has 14 and Farm 4 (V-A) has 3. The earlier series filled some gaps with 1s and 2s.
3. **Herd size.** Sheet 2 column `Size` is used as the target number of lactating cows (pens 6–11 and 15). The `Date` column, `October sizes` and the unnamed column are not used. Herd V-B has no case data and is skipped.
4. **Case and seed pens.** Both are pens 8–11, matching the manuscript. Add 6, 7 and 15 to both if farm staff also recorded cases in fresh or post-calving cows.
5. **Iterations.** `N_ITERATIONS = 100` is the number of optimiser evaluations per farm and metric. To repeat the whole optimisation 100 times with different seeds (as in the 260210 notebook), set `N_OPT_RUNS = 100`. This multiplies the run time by 100.

## Outputs and the figures they feed

All outputs are CSV files under `results/`. No figures are drawn.

| File | Content | Figure |
| --- | --- | --- |
| `burnin/burnin_summary_*.csv` | Achieved herd size and cows per pen | Methods, Table 1 |
| `optimization/observed_cases.csv` | Observed daily cases (zero-filled) | Fig. 4 |
| `optimization/prescreen_draws.csv`, `prescreen_accepted.csv` | All pre-screen draws with both losses; accepted top 5% | Fig. S1 |
| `optimization/search_bounds.csv` | Prior and narrowed search bounds | Methods |
| `optimization/bo_history.csv` | Every optimiser evaluation: parameters, loss, best so far, relative decrease | Fig. 5A, 5B |
| `optimization/reevaluation.csv` | Best candidates re-evaluated | Supplement |
| `optimization/optimal_parameters.csv`, `*_summary.csv` | Optimum per farm, metric and run; infectious period and duration of infection | Fig. 5C, parameter table |
| `post_simulations/case_incidence*.csv` | Simulated vs observed daily cases | Fig. 4 |
| `post_simulations/herd_daily*.csv` | Daily S, E, I, R, new exposures, βH | Fig. 3, Fig. 2B |
| `post_simulations/pen_daily_summary.csv` | Daily states and new infections per pen | Second-wave diagnostics |
| `post_simulations/infected_by_lactation.csv` | Infections by pen of exposure and lactation group | Second wave |
| `post_simulations/sim_summary*.csv` | Per run: peaks, timing, totals, attack rate, losses, R estimates | Fig. 2A, Results text |
| `sensitivity/sobol_*.csv` | Inputs/outputs, S1, ST and S2 with 95% CIs, at several sample sizes | Fig. 6 |
| `sensitivity/local_*.csv` | Outcomes and paired differences from nominal | Fig. 6, Fig. S3 |

### Reproduction numbers in `sim_summary.csv`

- `R_index_cows`: mean secondary exposures caused by the index cows. This is R0 measured directly in the model, with no homogeneous-mixing formula.
- `R_cohort_<window>`: mean secondary exposures of cows that became infectious during the window (days 1–10, 1–15, and day 1 to the observed peak). Each exposure is credited to a random infectious cow in the same pen; runs last `POST_HORIZON_DAYS`, so cohorts finish their infectious period.
- `R0_H_<window>`: mean βH over the window divided by γ, comparable with the manuscript's Eq. 11. It is an approximation, because the model is pen-structured.

`attack_rate` divides everyone ever infected by the cows at risk: those in transmission pens on day 1 plus all later entrants from the heifer pens.

## Expected run time

With the default settings, a smoke test on 2 cores suggests roughly: step 1 about 20–30 min, step 2 a few minutes, step 3 about 30 min, step 4 about 5 min. Step 03 scales with `SOBOL_N × SOBOL_REPS`.
