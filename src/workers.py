"""
workers.py - functions run inside parallel worker processes.

They live in a .py file rather than in the notebooks because, on macOS, Windows and
Linux with Python 3.14+, worker processes cannot see functions defined in a notebook.
Each worker reads the farm data and burn-in herds from common.FARMS and common.BURN,
which common.init_worker fills in every worker process.
"""
import numpy as np
import pandas as pd

import ab_sim as m
import common as c
import config as cfg


# --------------------------------------------------------------------------
# Step 1: pre-screen
# --------------------------------------------------------------------------
def prescreen_task(task):
    """Simulate one prior draw N_REPS_PER_EVAL times and score the mean curve with every metric."""
    farm, draw, params = task
    f = c.FARMS[farm]
    seed = m.task_seed("prescreen", farm, draw)
    curves = m.simulate_case_curves(c.BURN[farm], params, f["ndays"], cfg.N_REPS_PER_EVAL, seed, cfg)
    mean = curves.mean(axis=0)
    return {"farm": farm, "draw": draw, **params,
            **{f"loss_{k}": m.loss(k, mean, f["observed"], cfg) for k in cfg.METRICS}}


# --------------------------------------------------------------------------
# Step 2: post-optimisation simulations
# --------------------------------------------------------------------------
def r0_windows(f) -> dict:
    """Day windows for the reproduction-number estimates; 'observed_peak' becomes day 1 to the observed peak."""
    out = {}
    for name, w in cfg.R0_WINDOWS.items():
        out[name] = (1, f["observed_peak_day"]) if w == "observed_peak" else w
    return out


def post_simulation_task(task):
    """One post-optimisation simulation: case incidence, herd series, pen arrays, summary and
    infections by lactation."""
    farm, metric, sim, params = task
    f = c.FARMS[farm]
    rng = np.random.default_rng(m.task_seed("post", farm, metric, sim, cfg.POST_SEED))
    ndays = max(cfg.POST_HORIZON_DAYS, f["ndays"])
    res = m.simulate(c.BURN[farm], params["T"], params["R"], params["I"], params["IC"], ndays, rng,
                     cfg.SEED_PENS, cfg.MASKED_PENS, cfg.ANNUAL_CULL_RATE, keep_infected=True)
    key = {"farm": farm, "metric": metric, "sim": sim}

    inc = res.case_incidence(cfg.CASE_PENS, cfg.COUNT_INDEX_AS_DAY1_CASES)
    obs = np.full(ndays, np.nan)
    obs[: f["ndays"]] = f["observed"]
    cases = pd.DataFrame({**key, "day": np.arange(1, ndays + 1), "simulated": inc, "observed": obs})

    herd = m.herd_series(res, cfg.MASKED_PENS)
    for k, v in reversed(list(key.items())):
        herd.insert(0, k, v)

    fit = inc[: f["ndays"]]
    n_index = int(res.index_by_pen.sum())
    at_risk = herd["N"].iloc[0] + herd["entered_from_masked_pens"].iloc[1:].sum()
    ever = herd["new_exposed"].sum() + n_index
    active = np.flatnonzero((herd["E"] + herd["I"]).values > 0)
    s = {**key, **params,
         "peak_cases_fit_window": fit.max(), "peak_day_fit_window": int(np.argmax(fit)) + 1,
         "cumulative_cases_fit_window": fit.sum(),
         "peak_cases_horizon": inc.max(), "peak_day_horizon": int(np.argmax(inc)) + 1,
         "cumulative_cases_horizon": inc.sum(),
         "cumulative_exposed_horizon": herd["new_exposed"].sum(),
         "ever_infected_horizon": ever, "at_risk_horizon": at_risk,
         "attack_rate_horizon": ever / at_risk if at_risk else np.nan,
         "peak_infectious_herd": herd["I"].max(), "peak_infectious_day": int(herd["I"].idxmax()) + 1,
         "final_recovered": herd["R"].iloc[-1],
         "outbreak_end_day": int(active[-1]) + 1 if active.size else 0,
         "loss_DM": m.loss("DM", fit, f["observed"], cfg),
         "loss_PLF": m.loss("PLF", fit, f["observed"], cfg)}
    inf = res.infected
    s["R_index_cows"] = float(inf.loc[inf.is_index, "secondary"].mean()) if n_index else np.nan
    for name, (a, b) in r0_windows(f).items():
        r, n = m.cohort_reproduction(inf, a, b)
        s[f"R_cohort_{name}"] = r
        s[f"n_cohort_{name}"] = n
        bh = herd.loc[(herd.day >= a) & (herd.day <= b), "beta_H"]
        s[f"beta_H_mean_{name}"] = float(bh.mean())
        s[f"R0_H_{name}"] = float(bh.mean() / params["R"])

    lact = (inf.assign(lactation_group=np.where(inf.lactation <= 1, "first", "later"))
            .groupby(["exp_pen", "lactation_group"]).size().reset_index(name="n"))
    for k, v in reversed(list(key.items())):
        lact.insert(0, k, v)
    return cases, herd, res.counts.astype(np.int32), res.new_I.astype(np.int32), s, lact


# --------------------------------------------------------------------------
# Step 3: Sobol
# --------------------------------------------------------------------------
SOBOL_NAMES = ["T", "R", "I", "IC", "HS"]


def sobol_decode(x) -> dict:
    """Convert one Sobol sample row to model parameters (log10 T back-transformed, integers
    rounded, herd size mapped to the nearest HS_GRID herd)."""
    v = dict(zip(SOBOL_NAMES, x))
    out = {}
    for p in SOBOL_NAMES:
        val = v[p]
        if p in cfg.LOG_SCALE_PARAMS:
            val = 10 ** val
        if p in cfg.INTEGER_PARAMS:
            lo, hi = cfg.SOBOL_BOUNDS[p]
            val = int(min(max(round(val), lo), hi))
        out[p] = val
    grid = np.array(cfg.HS_GRID)
    out["HS_herd"] = int(grid[np.argmin(np.abs(grid - v["HS"]))])
    return out


def sobol_task(task):
    """Simulate one Sobol parameter set SOBOL_REPS times and average the outcomes."""
    row, x = task
    p = sobol_decode(x)
    vals = []
    for rep in range(cfg.SOBOL_REPS):
        rng = np.random.default_rng(m.task_seed("sobol", row, rep, cfg.SOBOL_SEED))
        res = m.simulate(c.BURN[p["HS_herd"]], p["T"], p["R"], p["I"], p["IC"], cfg.SA_HORIZON_DAYS,
                         rng, cfg.SEED_PENS, cfg.MASKED_PENS, cfg.ANNUAL_CULL_RATE)
        vals.append(c.outcomes(res))
    mean = pd.DataFrame(vals).mean().to_dict()
    return {"row": row, **{f"x_{k}": xv for k, xv in zip(SOBOL_NAMES, x)}, **p, **mean}


# --------------------------------------------------------------------------
# Step 4: local sensitivity
# --------------------------------------------------------------------------
def local_task(task):
    """One local-sensitivity simulation; the seed depends only on farm and replicate, so all
    values share random numbers (paired comparison with nominal)."""
    farm, param, value, nominal, herd_key, rep = task
    p = dict(nominal)
    if param in p:
        p[param] = value
    rng = np.random.default_rng(m.task_seed("local", farm, rep, cfg.LOCAL_SEED))
    res = m.simulate(c.BURN[herd_key], p["T"], p["R"], int(p["I"]), int(p["IC"]),
                     cfg.SA_HORIZON_DAYS, rng, cfg.SEED_PENS, cfg.MASKED_PENS, cfg.ANNUAL_CULL_RATE)
    return {"farm": farm, "parameter": param, "value": value, "rep": rep, **c.outcomes(res)}
