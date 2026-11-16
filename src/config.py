"""
config.py - every setting used by the pipeline, in one place.

Edit this file, not the scripts. All scripts read from here.
"""
import json
from pathlib import Path

# --------------------------------------------------------------------------
# Input / output
# --------------------------------------------------------------------------
INPUT_XLSX = Path("Cases_by_date-4Pranav.xlsx")   # sheet 1: cases, sheet 2: herd sizes
CASES_SHEET = 0                                   # "Cases matrix"
HERDS_SHEET = 1                                   # "Dates"
HERD_SIZE_COLUMN = "Size"          # column of sheet 2 used as the target lactating herd size
OUTPUT_DIR = Path("results")

# Herd code in the workbook -> label used in all outputs. The real mapping lives in
# farm_labels.local.json, which is not committed to the repository (see .gitignore).
# To set up: copy farm_labels.example.json to farm_labels.local.json and fill in the codes.
# Herds without case data are skipped automatically.
_LABELS_FILE = Path("farm_labels.local.json")
if not _LABELS_FILE.exists():
    raise FileNotFoundError("farm_labels.local.json not found. Copy farm_labels.example.json "
                            "to farm_labels.local.json and fill in the herd codes.")
FARM_LABELS = json.loads(_LABELS_FILE.read_text())

# Days between a farm's first and last case that have no row in the workbook:
#   "zero"          - no cases that day; cumulative cases stay flat
#   "carry_forward" - repeat the most recent recorded day's count (adds cases)
MISSING_DAYS_FILL = "zero"

# --------------------------------------------------------------------------
# Herd structure (pen numbers 1-15, see PEN_NAMES in ab_sim.py)
# --------------------------------------------------------------------------
MASKED_PENS = [1, 2, 3, 4]            # calves/heifers reared off site: no transmission
LACTATING_PENS = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15] 
HERD_SIZE_PENS = {
    "Farm1": list(range(6, 16)),   # fresh primiparous (6) onwards
    "Farm2": list(range(6, 16)),   # fresh primiparous (6) onwards
    "Farm3": list(range(2, 16)),   # post-weaned (2) onwards
    "Farm4": list(range(5, 16)),   # springers (5) onwards
}
CASE_PENS = [6, 7, 8, 9, 10, 11, 15]   # fresh, high and low milking (primiparous and multiparous) and hospital
SEED_PENS = [6, 7, 8, 9, 10, 11, 15]
ANNUAL_CULL_RATE = 0.30
BURN_IN_DAYS = 3650
# One burn-in is run at this nominal size; every farm and HS_GRID herd is a random subsample
# of it with exactly the target number of lactating cows. Without calibration the burn-in
# reaches roughly 70% of its nominal size, so this must be well above the largest target.
BURN_IN_REFERENCE_SIZE = 5500
BURN_IN_SEED = 2509

# --------------------------------------------------------------------------
# Parameters: T = transmission probability per infectious-susceptible pair per day (beta),
# R = recovery rate (gamma, per day), I = latent/incubation period (days, integer),
# IC = initially infectious (index) cows (integer)
# --------------------------------------------------------------------------
# T is per pair of cows. With index cows spread over milking pens of a few hundred cows,
# T is searched on a log scale.
BOUNDS = {"T": (0.0005, 0.005), "R": (0.04, 0.08), "I": (1, 4), "IC": (1, 10)}
INTEGER_PARAMS = ["I", "IC"]
LOG_SCALE_PARAMS = ["T"]

# Index cows are not new cases. Because the latent period is at least 1 day, the model
# cannot produce cases on day 1 unless they are counted. Set True to count index cows
# in CASE_PENS as day-1 cases (i.e. assume they are the first recorded clinical cows).
COUNT_INDEX_AS_DAY1_CASES = False

# --------------------------------------------------------------------------
# Pre-screening (rejection ABC) and Bayesian optimisation
# --------------------------------------------------------------------------
METRICS = ["DM", "PLF"]
DM_WEIGHTS = {"DTW": 1.0, "MAE": 5.0, "CC": 2.0, "TOTAL": 2.0}
N_PRESCREEN = 1000                    # prior draws per farm
PRESCREEN_ACCEPT_FRACTION = 0.05      # top 5% kept, per metric
N_ITERATIONS = 100                    # Bayesian-optimisation iterations (gp_minimize n_calls) per farm and metric
N_OPT_RUNS = 1                        # independent optimisations per farm and metric (different seeds)
N_SEED_POINTS_FROM_PRESCREEN = 10     # best pre-screen points given to the optimiser as x0/y0
N_RANDOM_START_POINTS = 10            # extra random points before the GP takes over
N_REPS_PER_EVAL = 10                   # simulations averaged per loss evaluation (reduces noise)
OPT_REEVAL_TOP = 5                    # re-evaluate the best candidates ...
OPT_REEVAL_REPS = 20                  # ... with more simulations before choosing the optimum
OPT_SEED = 1
PLF_TOTAL_WEIGHT = 5.0 

# --------------------------------------------------------------------------
# Post-optimisation simulations (descriptives)
# --------------------------------------------------------------------------
N_SIMULATIONS = 100                   # per farm and metric
POST_HORIZON_DAYS = 365               # long run (Fig. 3); R0 cohorts need follow-up beyond the fit window
R0_WINDOWS = {"d1_10": (1, 10), "d1_15": (1, 15), "d1_peak": "observed_peak"}
POST_SEED = 2

# --------------------------------------------------------------------------
# Sensitivity analyses
# --------------------------------------------------------------------------
SA_HORIZON_DAYS = 100
HS_GRID = [500, 1000, 1500, 2000, 2500, 3000] # burn-in herds built at these lactating herd sizes
SOBOL_BOUNDS = {**BOUNDS, "HS": (500, 3000)}
SOBOL_N = 1024                        # base sample size (power of 2); runs = N * (2k + 2)
SOBOL_REPS = 3                        # simulations averaged per parameter set
SOBOL_CONVERGENCE_N = [128, 256, 512, 1024]
SOBOL_SEED = 3
LOCAL_N_POINTS = 20                   # values per parameter across its bound
LOCAL_REPS = 20                       # simulations per value
LOCAL_METRIC = "DM"                   # optimum used as the nominal point
LOCAL_SEED = 4

# --------------------------------------------------------------------------
# Computing
# --------------------------------------------------------------------------
N_WORKERS = None                      # None = all CPU cores
