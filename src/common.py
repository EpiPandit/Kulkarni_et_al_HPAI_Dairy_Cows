"""
common.py - helpers shared by the pipeline scripts (paths, burn-in herds, parameter sampling,
outcome summaries, parallel map).
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import ab_sim as m
import config as cfg

OUT = Path(cfg.OUTPUT_DIR)
DIRS = {k: OUT / k for k in ("burnin", "optimization", "post_simulations", "sensitivity")}
PARAMS = ["T", "R", "I", "IC"]

# worker globals (filled by init_worker)
FARMS: dict = {}
BURN: dict = {}


def ensure_dirs() -> None:
    for d in DIRS.values():
        d.mkdir(parents=True, exist_ok=True)


def n_workers() -> int:
    return cfg.N_WORKERS or os.cpu_count() or 1


# --------------------------------------------------------------------------
# Burn-in herds
# --------------------------------------------------------------------------
def farm_burn_ins(farms: dict, verbose: bool = True) -> dict:
    out = {}
    for label, f in farms.items():
        if verbose:
            print(f"Burn-in {label} (target {f['target_lactating']} lactating cows)")
        out[label] = m.get_burn_in(label, f["target_lactating"], cfg, DIRS["burnin"])
    return out


def grid_burn_ins(verbose: bool = True) -> dict:
    out = {}
    for hs in cfg.HS_GRID:
        if verbose:
            print(f"Burn-in grid herd {hs}")
        out[hs] = m.get_burn_in("grid", hs, cfg, DIRS["burnin"])
    return out


def write_burnin_summary(burns: dict, farms: dict | None = None) -> None:
    rows, calib = [], []
    for key, b in burns.items():
        row = {"herd": key, "target_lactating": b["target_lactating"],
               "lactating_at_end": b["lactating_at_end"], "total_at_end": b["total_at_end"],
               "young_stock_target": b["herd"].target_young}
        row.update({f"pen{p}_{m.PEN_NAMES[p]}": int(c)
                    for p, c in zip(range(1, m.N_PENS + 1), b["pen_counts"])})
        rows.append(row)
        for h in b["calibration"]:
            calib.append({"herd": key, **h})
    name = "farms" if farms else "grid"
    pd.DataFrame(rows).to_csv(DIRS["burnin"] / f"burnin_summary_{name}.csv", index=False)
    pd.DataFrame(calib).to_csv(DIRS["burnin"] / f"burnin_calibration_{name}.csv", index=False)


def init_worker(include_grid: bool = False) -> None:
    """Load farm data and cached burn-in herds in each worker process."""
    global FARMS, BURN
    FARMS = m.load_inputs(cfg)
    BURN = {k: m.get_burn_in(k, f["target_lactating"], cfg, DIRS["burnin"]) for k, f in FARMS.items()}
    if include_grid:
        for hs in cfg.HS_GRID:
            BURN[hs] = m.get_burn_in("grid", hs, cfg, DIRS["burnin"])


def pool_map(func, tasks, include_grid: bool = False, chunksize: int = 4):
    tasks = list(tasks)
    if n_workers() == 1:
        init_worker(include_grid)
        return [func(t) for t in tasks]
    with ProcessPoolExecutor(max_workers=n_workers(), initializer=init_worker,
                             initargs=(include_grid,)) as ex:
        return list(ex.map(func, tasks, chunksize=chunksize))


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
def sample_prior(n: int, rng: np.random.Generator, bounds=None) -> pd.DataFrame:
    bounds = bounds or cfg.BOUNDS
    out = {}
    for p, (lo, hi) in bounds.items():
        if p in cfg.INTEGER_PARAMS:
            out[p] = rng.integers(int(lo), int(hi) + 1, n)
        elif p in cfg.LOG_SCALE_PARAMS:
            out[p] = 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)
        else:
            out[p] = rng.uniform(lo, hi, n)
    return pd.DataFrame(out)


def as_params(row) -> dict:
    return {"T": float(row["T"]), "R": float(row["R"]), "I": int(round(row["I"])),
            "IC": int(round(row["IC"]))}


# --------------------------------------------------------------------------
# Outcomes used in the sensitivity analyses
# --------------------------------------------------------------------------
OUTCOMES = ["peak_daily_cases", "cumulative_cases", "cumulative_exposed", "attack_rate",
            "final_recovered", "outbreak_duration_days"]


def outcomes(res: m.SimResult) -> dict:
    inc = res.case_incidence(cfg.CASE_PENS, cfg.COUNT_INDEX_AS_DAY1_CASES)
    hs = m.herd_series(res, cfg.MASKED_PENS)
    # cows at risk: those in transmission pens on day 1 plus all that entered later
    at_risk = hs["N"].iloc[0] + hs["entered_from_masked_pens"].iloc[1:].sum()
    n_index = int(res.index_by_pen.sum())
    ever = hs["new_exposed"].sum() + n_index
    active = np.flatnonzero((hs["E"] + hs["I"]).values > 0)
    return {
        "peak_daily_cases": float(inc.max()),
        "cumulative_cases": float(inc.sum()),
        "cumulative_exposed": float(hs["new_exposed"].sum()),
        "attack_rate": float(ever / at_risk) if at_risk else np.nan,
        "final_recovered": float(hs["R"].iloc[-1]),
        "outbreak_duration_days": float(active[-1] + 1) if active.size else 0.0,
    }


def summarise(df: pd.DataFrame, by: list, cols: list) -> pd.DataFrame:
    """Mean, SD, median and 2.5/97.5 percentiles of `cols` within groups `by` (long format)."""
    long = df.melt(id_vars=by, value_vars=cols, var_name="variable", value_name="value")
    g = long.groupby(by + ["variable"])["value"]
    return pd.DataFrame({
        "mean": g.mean(), "sd": g.std(), "median": g.median(),
        "p2_5": g.quantile(0.025), "p97_5": g.quantile(0.975), "n": g.count(),
    }).reset_index()
