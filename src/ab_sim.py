"""
ab_sim.py - pen-structured stochastic agent-based SEIR model of HPAI H5N1 in a dairy herd.

Each cow is one agent, stored as one element of a set of numpy arrays (fast, and easy to copy).

Daily order of events (day d = 1, 2, ...):
  1. Culling: every living cow is culled with probability mu = 1 - (1 - annual rate)^(1/365),
     whatever its disease state.
  2. Disease progression: I -> R when the cow's infectious period ends;
     E -> I exactly `latent` days after exposure (fixed latent period).
     Infectious periods are geometric with mean 1/gamma (minimum 1 day).
  3. Movement: each cow's dwell time falls by one; cows whose dwell time ends move to the next pen,
     with a new dwell time drawn uniformly (integers) between the pen's minimum and maximum.
     Calving happens on entry to the post-calving (hospital) pen. A calf is retained only while
     young stock are below a target scaled from the burn-in, which keeps herd size stable.
  4. Transmission: in every pen that is not masked, each susceptible cow is exposed with
     probability 1 - (1 - beta)^I_p, where I_p is the number of infectious cows in that pen.
     There is no transmission between pens. Each new exposure is attributed to one infectious
     cow in the same pen, chosen at random, to count secondary cases.
  5. Recording.
"""
from __future__ import annotations

import pickle
import zlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln

# --------------------------------------------------------------------------
# Herd structure (Konboon et al., 2018; Table S5)
# --------------------------------------------------------------------------
PEN_NAMES = {
    1: "Pre-weaning", 2: "Post-weaned", 3: "Breeding", 4: "Pregnant", 5: "Springer",
    6: "Fresh primiparous", 7: "Fresh multiparous", 8: "High-milk primiparous",
    9: "Low-milk primiparous", 10: "High-milk multiparous", 11: "Low-milk multiparous",
    12: "Dry", 13: "Close-up", 14: "Calving", 15: "Hospital (post-calving)",
}
N_PENS = 15
DWELL_MIN = np.array([60, 396, 60, 234, 18, 5, 30, 180, 100, 90, 200, 48, 7, 1, 2])
DWELL_MAX = np.array([80, 458, 100, 264, 32, 30, 90, 200, 160, 180, 360, 66, 14, 7, 2])
# next pen (1-based); pen 15 goes to 6 (first lactation) or 7 (later lactations)
NEXT_PEN = np.array([2, 3, 4, 5, 15, 8, 10, 9, 12, 11, 12, 13, 14, 15, 0])
YOUNG_PENS = np.array([1, 2, 3, 4, 5])
MILKING_PENS_FOR_DIM = np.array([6, 7, 8, 9, 10, 11, 15])   # days in milk accrue here
MIN_DAYS_IN_MILK = 305      # cows leaving pen 9 or 11 before 305 days in milk stay longer

S, E, I, R = 0, 1, 2, 3
STATE_NAMES = ["S", "E", "I", "R"]


# --------------------------------------------------------------------------
# Herd container
# --------------------------------------------------------------------------
@dataclass
class Herd:
    pen: np.ndarray            # 1..15
    remaining: np.ndarray      # days left in current pen
    lactation: np.ndarray      # 0 = heifer
    dim: np.ndarray            # days in milk, current lactation
    alive: np.ndarray
    state: np.ndarray
    t_infectious: np.ndarray   # day E -> I (-1 if not exposed)
    t_recover: np.ndarray      # day I -> R
    exp_day: np.ndarray        # day exposed
    inf_day: np.ndarray        # day became infectious
    exp_pen: np.ndarray        # pen where exposed
    inf_pen: np.ndarray        # pen where it became infectious
    secondary: np.ndarray      # exposures attributed to this cow
    is_index: np.ndarray
    target_young: int = 0      # herd-size controller

    _FIELDS = ("pen", "remaining", "lactation", "dim", "alive", "state", "t_infectious",
               "t_recover", "exp_day", "inf_day", "exp_pen", "inf_pen", "secondary", "is_index")

    @classmethod
    def empty(cls, n: int = 0) -> "Herd":
        return cls(
            pen=np.zeros(n, np.int16), remaining=np.zeros(n, np.int32),
            lactation=np.zeros(n, np.int16), dim=np.zeros(n, np.int32),
            alive=np.ones(n, bool), state=np.zeros(n, np.int8),
            t_infectious=np.full(n, -1, np.int32), t_recover=np.full(n, -1, np.int32),
            exp_day=np.full(n, -1, np.int32), inf_day=np.full(n, -1, np.int32),
            exp_pen=np.zeros(n, np.int16), inf_pen=np.zeros(n, np.int16),
            secondary=np.zeros(n, np.int32), is_index=np.zeros(n, bool),
        )

    def copy(self) -> "Herd":
        h = Herd(**{f: getattr(self, f).copy() for f in self._FIELDS})
        h.target_young = self.target_young
        return h

    def append(self, other: "Herd") -> None:
        for f in self._FIELDS:
            setattr(self, f, np.concatenate([getattr(self, f), getattr(other, f)]))

    def compact(self) -> None:
        keep = self.alive
        for f in self._FIELDS:
            setattr(self, f, getattr(self, f)[keep])

    def reset_disease(self) -> None:
        n = len(self.pen)
        self.state[:] = S
        for f in ("t_infectious", "t_recover", "exp_day", "inf_day"):
            getattr(self, f)[:] = -1
        self.exp_pen[:] = 0
        self.inf_pen[:] = 0
        self.secondary[:] = 0
        self.is_index[:] = np.zeros(n, bool)

    def n_lactating(self, lactating_pens) -> int:
        return int(np.sum(self.alive & np.isin(self.pen, lactating_pens)))


def _draw_dwell(pens: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    lo, hi = DWELL_MIN[pens - 1], DWELL_MAX[pens - 1]
    return rng.integers(lo, hi + 1)


def _new_calves(n: int, rng) -> Herd:
    h = Herd.empty(n)
    h.pen[:] = 1
    h.remaining[:] = _draw_dwell(h.pen, rng)
    return h


# --------------------------------------------------------------------------
# Daily demographic step (culling + movement + calving)
# --------------------------------------------------------------------------
def _cull(herd: Herd, daily_cull: float, rng) -> None:
    idx = np.flatnonzero(herd.alive)
    culled = idx[rng.random(idx.size) < daily_cull]
    herd.alive[culled] = False


def _move(herd: Herd, rng) -> int:
    """Advance dwell times, move cows, handle calving. Returns the number of calvings."""
    alive = herd.alive
    herd.remaining[alive] -= 1
    in_milk = alive & np.isin(herd.pen, MILKING_PENS_FOR_DIM)
    herd.dim[in_milk] += 1

    due = np.flatnonzero(alive & (herd.remaining <= 0))
    if due.size == 0:
        return 0
    cur = herd.pen[due]

    # 305-day rule: cows due to leave pen 9 or 11 with < 305 days in milk stay longer
    leaving_late = np.isin(cur, (9, 11)) & (herd.dim[due] < MIN_DAYS_IN_MILK)
    if leaving_late.any():
        ext = due[leaving_late]
        herd.remaining[ext] = MIN_DAYS_IN_MILK - herd.dim[ext]
        due, cur = due[~leaving_late], cur[~leaving_late]
        if due.size == 0:
            return 0

    nxt = NEXT_PEN[cur - 1].copy()
    from_hosp = cur == 15
    nxt[from_hosp] = np.where(herd.lactation[due[from_hosp]] <= 1, 6, 7)

    calving = np.isin(cur, (5, 14)) & (nxt == 15)
    calf_idx = due[calving]
    herd.lactation[calf_idx] += 1
    herd.dim[calf_idx] = 0

    herd.pen[due] = nxt
    herd.remaining[due] = _draw_dwell(nxt, rng)
    return int(calving.sum())


def _add_retained_calves(herd: Herd, n_calvings: int, rng) -> None:
    if n_calvings == 0:
        return
    n_young = int(np.sum(herd.alive & np.isin(herd.pen, YOUNG_PENS)))
    keep = min(n_calvings, max(0, herd.target_young - n_young))
    if keep > 0:
        herd.append(_new_calves(keep, rng))


def demography_step(herd: Herd, daily_cull: float, rng) -> None:
    _cull(herd, daily_cull, rng)
    n_calvings = _move(herd, rng)
    _add_retained_calves(herd, n_calvings, rng)


def daily_cull_probability(annual_rate: float) -> float:
    return 1.0 - (1.0 - annual_rate) ** (1.0 / 365.0)


# --------------------------------------------------------------------------
# Burn-in: build a herd at a target lactating herd size
# --------------------------------------------------------------------------
def _initial_herd(target_lactating: int, target_young: int, rng) -> Herd:
    adult_pens = np.arange(5, 16)
    w = (DWELL_MIN[adult_pens - 1] + DWELL_MAX[adult_pens - 1]) / 2.0
    n_adult = int(round(target_lactating * 1.25))
    pens = rng.choice(adult_pens, size=n_adult, p=w / w.sum())
    young_pens = np.arange(1, 5)
    wy = (DWELL_MIN[young_pens - 1] + DWELL_MAX[young_pens - 1]) / 2.0
    ypens = rng.choice(young_pens, size=target_young, p=wy / wy.sum())
    all_pens = np.concatenate([pens, ypens]).astype(np.int16)
    h = Herd.empty(all_pens.size)
    h.pen[:] = all_pens
    h.remaining[:] = rng.integers(1, DWELL_MAX[all_pens - 1] + 1)
    lact = np.where(np.isin(all_pens, (6, 8, 9)), 1,
                    np.where(all_pens >= 6, rng.integers(2, 5, all_pens.size), 0))
    h.lactation[:] = lact
    h.target_young = target_young
    return h


def build_burn_in(target_lactating: int, lactating_pens, annual_cull: float = 0.30,
                  days: int = 3650, seed: int = 1, tol: float = 0.02, max_iter: int = 8,
                  verbose: bool = False) -> dict:
    """Run the demographic model without disease until the herd is at equilibrium, adjusting
    the young-stock target until the mean lactating herd over the last year is within `tol`
    of `target_lactating`. Returns a dict with the herd and diagnostics."""
    daily_cull = daily_cull_probability(annual_cull)
    ratio = 1.0
    history = []
    for it in range(max_iter):
        rng = np.random.default_rng([seed, it])
        target_young = max(1, int(round(ratio * target_lactating)))
        herd = _initial_herd(target_lactating, target_young, rng)
        last_year = []
        for d in range(1, days + 1):
            demography_step(herd, daily_cull, rng)
            if d % 365 == 0:
                herd.compact()
            if d > days - 365:
                last_year.append(herd.n_lactating(lactating_pens))
        achieved = float(np.mean(last_year))
        rel_err = achieved / target_lactating - 1.0
        history.append({"iteration": it, "target_young": target_young,
                        "mean_lactating_last_year": achieved, "relative_error": rel_err})
        if verbose:
            print(f"  burn-in iter {it}: young target {target_young}, lactating {achieved:.0f} "
                  f"(target {target_lactating}, error {rel_err:+.1%})")
        if abs(rel_err) <= tol:
            break
        ratio *= (target_lactating / achieved) ** 0.7   # damped update
    herd.compact()
    herd.reset_disease()
    pen_counts = np.bincount(herd.pen, minlength=N_PENS + 1)[1:]
    return {
        "herd": herd,
        "target_lactating": target_lactating,
        "lactating_at_end": herd.n_lactating(lactating_pens),
        "total_at_end": int(herd.alive.sum()),
        "pen_counts": pen_counts,
        "calibration": history,
    }


def thin_herd(reference: dict, target_lactating: int, lactating_pens, seed: int) -> dict:
    """Random subsample of a larger equilibrium herd with exactly `target_lactating`
    lactating cows. Every other pen is thinned by the same fraction, so the mix of cows
    across pens, lactations and dwell times is kept."""
    rng = np.random.default_rng(seed)
    herd = reference["herd"].copy()
    alive = np.flatnonzero(herd.alive)
    is_lact = np.isin(herd.pen[alive], lactating_pens)
    lact, other = alive[is_lact], alive[~is_lact]
    if target_lactating > lact.size:
        raise ValueError(f"Reference herd has {lact.size} lactating cows; target "
                         f"{target_lactating} is larger. Increase BURN_IN_REFERENCE_SIZE in config.py "
                         f"and delete the burn-in cache folder.")
    frac = target_lactating / lact.size
    keep = [rng.choice(lact, target_lactating, replace=False)]
    for p in np.unique(herd.pen[other]):
        idx = other[herd.pen[other] == p]
        keep.append(rng.choice(idx, int(round(frac * idx.size)), replace=False))
    mask = np.zeros(herd.pen.size, bool)
    mask[np.concatenate(keep)] = True
    herd.alive = mask
    herd.compact()
    herd.target_young = max(1, int(round(reference["herd"].target_young * frac)))
    return {
        "herd": herd,
        "target_lactating": target_lactating,
        "lactating_at_end": herd.n_lactating(lactating_pens),
        "total_at_end": int(herd.alive.sum()),
        "pen_counts": np.bincount(herd.pen, minlength=N_PENS + 1)[1:],
        "calibration": [{"thinned_from_lactating": reference["lactating_at_end"],
                         "fraction_kept": frac}],
    }


def get_reference_burn_in(cfg, cache_dir: Path) -> dict:
    """One burn-in at BURN_IN_REFERENCE_SIZE (built once, then cached)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "burnin_reference.pkl"
    if path.exists():
        with open(path, "rb") as f:
            return pickle.load(f)
    ref = build_burn_in(cfg.BURN_IN_REFERENCE_SIZE, cfg.LACTATING_PENS, cfg.ANNUAL_CULL_RATE,
                        cfg.BURN_IN_DAYS, cfg.BURN_IN_SEED, tol=1.0, max_iter=1)
    with open(path, "wb") as f:
        pickle.dump(ref, f)
    return ref


def get_burn_in(name: str, target_lactating: int, cfg, cache_dir: Path) -> dict:
    """Herd with exactly `target_lactating` cows in this farm's herd-size pens
    (HERD_SIZE_PENS, or LACTATING_PENS if the farm is not listed), thinned from the
    reference burn-in (cached per herd)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    pens = getattr(cfg, "HERD_SIZE_PENS", {}).get(name, cfg.LACTATING_PENS)
    path = cache_dir / f"burnin_{name}_{target_lactating}_p{min(pens)}-{max(pens)}.pkl"
    if path.exists():
        with open(path, "rb") as f:
            return pickle.load(f)
    ref = get_reference_burn_in(cfg, cache_dir)
    seed = (int(cfg.BURN_IN_SEED) + int(target_lactating) + zlib.crc32(name.encode())) % 2**31
    b = thin_herd(ref, target_lactating, pens, seed)
    with open(path, "wb") as f:
        pickle.dump(b, f)
    return b


# --------------------------------------------------------------------------
# Outbreak simulation
# --------------------------------------------------------------------------
@dataclass
class SimResult:
    counts: np.ndarray        # [day, pen, state] living cows at end of day
    new_E: np.ndarray         # [day, pen] new exposures
    new_I: np.ndarray         # [day, pen] cows becoming infectious (index cows excluded)
    new_R: np.ndarray         # [day, pen]
    entered_open: np.ndarray  # [day] cows moving from masked pens into transmission pens
    index_by_pen: np.ndarray = field(default=None)  # [pen] index cows at day 1
    infected: pd.DataFrame = field(default=None)   # one row per cow ever infected

    def case_incidence(self, case_pens, count_index: bool = False) -> np.ndarray:
        cols = np.asarray(case_pens) - 1
        inc = self.new_I[:, cols].sum(axis=1).astype(float)
        if count_index:
            inc[0] += self.index_by_pen[cols].sum()
        return inc


def simulate(burn_in: dict, beta: float, gamma: float, latent: int, n_index: int, ndays: int,
             rng: np.random.Generator, seed_pens, masked_pens, annual_cull: float = 0.30,
             keep_infected: bool = False) -> SimResult:
    herd = burn_in["herd"].copy()
    herd.reset_disease()
    daily_cull = daily_cull_probability(annual_cull)
    latent = int(latent)
    open_pen = np.ones(N_PENS + 1, bool)
    open_pen[0] = False
    open_pen[np.asarray(masked_pens)] = False

    # index cows: existing susceptible cows drawn at random from the seed pens
    cand = np.flatnonzero(herd.alive & np.isin(herd.pen, seed_pens) & (herd.state == S))
    idx = rng.choice(cand, size=min(int(n_index), cand.size), replace=False)
    herd.state[idx] = I
    herd.is_index[idx] = True
    herd.inf_day[idx] = 1
    herd.inf_pen[idx] = herd.pen[idx]
    herd.t_recover[idx] = 1 + rng.geometric(gamma, idx.size)
    index_by_pen = np.bincount(herd.pen[idx], minlength=N_PENS + 1)[1:]

    counts = np.zeros((ndays, N_PENS, 4), np.int32)
    new_E = np.zeros((ndays, N_PENS), np.int32)
    new_I = np.zeros((ndays, N_PENS), np.int32)
    new_R = np.zeros((ndays, N_PENS), np.int32)
    entered_open = np.zeros(ndays, np.int32)

    for d in range(1, ndays + 1):
        t = d - 1
        # 1. culling
        _cull(herd, daily_cull, rng)
        alive = herd.alive
        # 2. progression
        rec = np.flatnonzero(alive & (herd.state == I) & (herd.t_recover <= d))
        if rec.size:
            herd.state[rec] = R
            new_R[t] += np.bincount(herd.pen[rec], minlength=N_PENS + 1)[1:]
        prog = np.flatnonzero(alive & (herd.state == E) & (herd.t_infectious <= d))
        if prog.size:
            herd.state[prog] = I
            herd.inf_day[prog] = d
            herd.inf_pen[prog] = herd.pen[prog]
            herd.t_recover[prog] = d + rng.geometric(gamma, prog.size)
            new_I[t] += np.bincount(herd.pen[prog], minlength=N_PENS + 1)[1:]
        # 3. movement and calving
        was_open = herd.alive & open_pen[herd.pen]
        n_calvings = _move(herd, rng)
        entered_open[t] = int(np.sum(herd.alive & open_pen[herd.pen] & ~was_open))
        _add_retained_calves(herd, n_calvings, rng)
        alive = herd.alive
        # 4. transmission within pens
        inf_mask = alive & (herd.state == I)
        I_p = np.bincount(herd.pen[inf_mask], minlength=N_PENS + 1)
        sus = np.flatnonzero(alive & (herd.state == S) & open_pen[herd.pen])
        if sus.size and I_p[open_pen].any():
            p_exp = 1.0 - (1.0 - beta) ** I_p[herd.pen[sus]]
            hit = sus[rng.random(sus.size) < p_exp]
            if hit.size:
                herd.state[hit] = E
                herd.exp_day[hit] = d
                herd.exp_pen[hit] = herd.pen[hit]
                herd.t_infectious[hit] = d + latent
                by_pen = np.bincount(herd.pen[hit], minlength=N_PENS + 1)
                new_E[t] += by_pen[1:]
                inf_idx = np.flatnonzero(inf_mask)
                inf_pen = herd.pen[inf_idx]
                for p in np.flatnonzero(by_pen):
                    infectors = inf_idx[inf_pen == p]
                    np.add.at(herd.secondary, rng.choice(infectors, size=by_pen[p]), 1)
        # 5. record
        alive = herd.alive
        counts[t] = np.stack([np.bincount(herd.pen[alive & (herd.state == s)],
                                          minlength=N_PENS + 1)[1:] for s in range(4)], axis=1)

    infected = None
    if keep_infected:
        ever = np.flatnonzero(herd.is_index | (herd.exp_day > 0))
        infected = pd.DataFrame({
            "is_index": herd.is_index[ever], "exp_day": herd.exp_day[ever],
            "exp_pen": herd.exp_pen[ever], "inf_day": herd.inf_day[ever],
            "inf_pen": herd.inf_pen[ever], "lactation": herd.lactation[ever],
            "secondary": herd.secondary[ever], "alive_at_end": herd.alive[ever],
            "state_at_end": np.array(STATE_NAMES)[herd.state[ever]],
        })
    return SimResult(counts, new_E, new_I, new_R, entered_open, index_by_pen, infected)


def simulate_case_curves(burn_in, params: dict, ndays: int, n_reps: int, seed, cfg) -> np.ndarray:
    """n_reps simulated daily case curves (rows) for one parameter set."""
    out = np.zeros((n_reps, ndays))
    for r in range(n_reps):
        rng = np.random.default_rng([*np.atleast_1d(seed), r])
        res = simulate(burn_in, params["T"], params["R"], params["I"], params["IC"], ndays, rng,
                       cfg.SEED_PENS, cfg.MASKED_PENS, cfg.ANNUAL_CULL_RATE)
        out[r] = res.case_incidence(cfg.CASE_PENS, cfg.COUNT_INDEX_AS_DAY1_CASES)
    return out


# --------------------------------------------------------------------------
# Loss functions
# --------------------------------------------------------------------------
def dtw_distance(a, b) -> float:
    n, m = len(a), len(b)
    D = np.full((n + 1, m + 1), np.inf)
    D[0, 0] = 0.0
    for i in range(1, n + 1):
        cost = np.abs(a[i - 1] - b)
        for j in range(1, m + 1):
            D[i, j] = cost[j - 1] + min(D[i - 1, j], D[i, j - 1], D[i - 1, j - 1])
    return float(D[n, m])


def crosscorr_alignment_distance(sim, obs, max_shift=None):
    """Smallest MAE between sim and obs over time shifts of up to max_shift days."""
    n = len(obs)
    max_shift = n // 2 if max_shift is None else max_shift
    best, best_shift = float(np.mean(np.abs(sim - obs))), 0
    for s in range(-max_shift, max_shift + 1):
        if s == 0:
            continue
        a, b = (sim[s:], obs[: n - s]) if s > 0 else (sim[: n + s], obs[-s:])
        if a.size:
            mae = float(np.mean(np.abs(a - b)))
            if mae < best:
                best, best_shift = mae, s
    return best, best_shift


def loss_dm(sim, obs, weights) -> float:
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    cc, _ = crosscorr_alignment_distance(sim, obs, len(obs) // 2)
    return (weights["DTW"] * dtw_distance(sim, obs)
            + weights["MAE"] * float(np.mean(np.abs(sim - obs)))
            + weights["CC"] * cc
            + weights.get("TOTAL", 0.0) * abs(sim.sum() - obs.sum()))


def loss_plf(sim, obs, total_weight: float = 0.0) -> float:
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    daily = float(np.sum(sim - obs * np.log(sim + 1e-9) + gammaln(obs + 1)))
    S, Y = sim.sum(), obs.sum()
    total = float(S - Y * np.log(S + 1e-9) - (Y - Y * np.log(Y))) if Y > 0 else float(S)
    return daily + total_weight * total


def loss(metric: str, sim, obs, cfg) -> float:
    if metric == "PLF":
        return loss_plf(sim, obs, getattr(cfg, "PLF_TOTAL_WEIGHT", 0.0))
    if metric == "DM":
        return loss_dm(sim, obs, cfg.DM_WEIGHTS)
    raise ValueError(f"Unknown metric {metric!r}; use 'DM' or 'PLF'")

# --------------------------------------------------------------------------
# Derived quantities
# --------------------------------------------------------------------------
def herd_series(res: SimResult, masked_pens) -> pd.DataFrame:
    """Daily herd-level S, E, I, R (unmasked pens), new exposures and the herd-level
    transmission rate beta_H(t) = X_t * N_t / (S_t * I_t) (corrected Eq. 10)."""
    open_cols = np.setdiff1d(np.arange(1, N_PENS + 1), masked_pens) - 1
    c = res.counts[:, open_cols, :].sum(axis=1)
    X = res.new_E[:, open_cols].sum(axis=1)
    S_, E_, I_, R_ = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
    N_ = c.sum(axis=1)
    # S and I at the start of day t are the end-of-day values of day t-1
    S_prev = np.concatenate([[np.nan], S_[:-1]])
    I_prev = np.concatenate([[np.nan], I_[:-1]])
    N_prev = np.concatenate([[np.nan], N_[:-1]])
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_H = np.where((S_prev > 0) & (I_prev > 0), X * N_prev / (S_prev * I_prev), np.nan)
    return pd.DataFrame({"day": np.arange(1, len(X) + 1), "S": S_, "E": E_, "I": I_, "R": R_,
                         "N": N_, "entered_from_masked_pens": res.entered_open,
                         "new_exposed": X, "new_infectious": res.new_I[:, open_cols].sum(1),
                         "beta_H": beta_H})


def cohort_reproduction(infected: pd.DataFrame, start: int, end: int) -> tuple[float, int]:
    """Mean secondary cases of cows that became infectious on days start..end (index cows included)."""
    sel = infected[(infected.inf_day >= start) & (infected.inf_day <= end)]
    return (float(sel.secondary.mean()) if len(sel) else np.nan), int(len(sel))


# --------------------------------------------------------------------------
# Data input
# --------------------------------------------------------------------------
def load_inputs(cfg) -> dict:
    """Read case series (gaps filled as set by MISSING_DAYS_FILL) and target herd sizes."""
    cases = pd.read_excel(cfg.INPUT_XLSX, sheet_name=cfg.CASES_SHEET)
    herds = pd.read_excel(cfg.INPUT_XLSX, sheet_name=cfg.HERDS_SHEET)
    cases["date"] = pd.to_datetime(cases["date"])
    farms = {}
    for code, label in cfg.FARM_LABELS.items():
        sub = cases[cases["Herd"] == code].groupby("date")["cases"].sum().sort_index()
        if sub.empty:
            continue
        full = pd.date_range(sub.index.min(), sub.index.max(), freq="D")
        missing = full.difference(sub.index)
        if cfg.MISSING_DAYS_FILL == "zero":
            series = sub.reindex(full, fill_value=0)
        elif cfg.MISSING_DAYS_FILL == "carry_forward":
            series = sub.reindex(full).ffill()
        else:
            raise ValueError(f"Unknown MISSING_DAYS_FILL: {cfg.MISSING_DAYS_FILL!r}")
        hs = herds.loc[herds["Herd"] == code, cfg.HERD_SIZE_COLUMN]
        if hs.empty:
            raise ValueError(f"No herd size for {code} in sheet 2")
        farms[label] = {
            "code": code, "dates": full, "observed": series.values.astype(float),
            "ndays": len(full), "n_missing_days_filled": len(missing),
            "target_lactating": int(hs.iloc[0]),
            "observed_peak_day": int(np.argmax(series.values)) + 1,
        }
    return farms


def task_seed(*parts) -> list[int]:
    """Reproducible integer seed list from labels and numbers."""
    out = []
    for p in parts:
        if isinstance(p, str):
            out.append(int.from_bytes(p.encode()[:8].ljust(8, b"\0"), "little") % (2**32))
        else:
            out.append(int(p))
    return out
