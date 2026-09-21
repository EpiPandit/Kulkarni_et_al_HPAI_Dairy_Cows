import pandas as pd
import numpy as np
import random
from scipy.special import gammaln
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.stats import skewnorm
from skopt import gp_minimize
from skopt.space import Real
import copy

merge_dwell = pd.DataFrame({
    'modelstate': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    'min_duration': [60, 396, 60, 234, 18, 5, 30, 180, 100, 90, 200, 48, 1, 1, 2],
    'max_duration': [80, 458, 100, 264, 32, 30, 90, 200, 160, 180, 360, 66, 7, 7, 2]
})

merge_dwell['current_pen'] = [f"Pen{i}" for i in range(1, 16)]
transition_rules_dict = {
    "Pen1": "Pen2", "Pen2": "Pen3", "Pen3": "Pen4", "Pen4": "Pen5", "Pen5": "Pen15",
    "Pen6": "Pen8", "Pen7": "Pen10", "Pen8": "Pen9", "Pen9": "Pen12", "Pen10": "Pen11",
    "Pen11": "Pen12", "Pen12": "Pen13", "Pen13": "Pen14", "Pen14": "Pen15", "Pen15": None
}
merge_dwell['next_pen'] = merge_dwell['current_pen'].map(transition_rules_dict)

class Cow:
    def __init__(self, id, pen, lactation=1, state='Susceptible', incubation_period=0):
        self.id = id
        self.current_pen = pen
        self.time_remaining = pen.get_dwell_time()
        self.alive = True
        self.lactation = lactation
        self.birth = None
        self.dim_per_lactation = {}
        self.current_dim_count = 0
        self.dwell_start_absolute_day = 0
        self.completed_dwell_times = []
        self.accumulated_dwell_time_current_lactation = 0
        self.calving_intervals = []
        self.state = state
        self.infection_time = None
        self.recovery_day = None
        self.incubation_period = int(incubation_period if incubation_period is not None else 0)
        if self.state == 'Infectious':
            self.infection_time = 0
            self.recovery_day = None
        self._history_data = []

    def record_history(self, absolute_day, burn_in_days):
        sim_day = absolute_day - burn_in_days if absolute_day > burn_in_days else absolute_day
        self._history_data.append({
            'id': self.id,
            'day': sim_day,
            'absolute_day': absolute_day,
            'pen': self.current_pen.id,
            'alive': self.alive,
            'lactation': self.lactation,
            'state': self.state
        })

    def get_history_df(self):
        return pd.DataFrame(self._history_data)

    def step(self, absolute_day, pens, farm):
        if not self.alive:
            return False
        self.time_remaining -= 1
        self.accumulated_dwell_time_current_lactation += 1
        add_calf = False
        if self.lactation >= 1 and self.current_pen.id in ["Pen15", "Pen6", "Pen7", "Pen8", "Pen9", "Pen10", "Pen11"]:
            self.current_dim_count += 1
        if self.time_remaining <= 0:
            previous_pen = self.current_pen
            if previous_pen.id == "Pen15":
                next_pen_id = "Pen6" if self.lactation < 2 else "Pen7"
            else:
                next_pen_id = previous_pen.get_next_pen_id()
            dwell_days = absolute_day - self.dwell_start_absolute_day
            self.completed_dwell_times.append({'pen': previous_pen.id, 'dwell_days': dwell_days, 'start_day': self.dwell_start_absolute_day})
            min_dwell_sum = 305
            if previous_pen.id == "Pen9":
                group_pens = ["Pen7", "Pen8", "Pen9"]
                current_group_dwell = sum(d['dwell_days'] for d in self.completed_dwell_times if d['pen'] in group_pens)
                if current_group_dwell < min_dwell_sum:
                    extension_days = min_dwell_sum - current_group_dwell
                    self.time_remaining += extension_days
                    self.completed_dwell_times[-1]['dwell_days'] += extension_days
                    return add_calf
            elif previous_pen.id == "Pen11":
                group_pens = ["Pen6", "Pen10", "Pen11"]
                current_group_dwell = sum(d['dwell_days'] for d in self.completed_dwell_times if d['pen'] in group_pens)
                if current_group_dwell < min_dwell_sum:
                    extension_days = min_dwell_sum - current_group_dwell
                    self.time_remaining += extension_days
                    self.completed_dwell_times[-1]['dwell_days'] += extension_days
                    return add_calf
            if next_pen_id is None or next_pen_id not in pens:
                print(f"Warning: Cow {self.id} in {self.current_pen.id} has no valid next pen ({next_pen_id}).")
                self.alive = False
                previous_pen.remove_cow(self)
                return add_calf
            next_pen = pens[next_pen_id]
            farm.record_movement(self.id, absolute_day, self.current_pen.id, next_pen_id)
            if previous_pen.id in ["Pen5", "Pen14"] and next_pen_id == "Pen15":
                completed_lactation_number = self.lactation
                self.calving_intervals.append({
                    'lactation_number': completed_lactation_number,
                    'calv_interval': self.accumulated_dwell_time_current_lactation
                })
                self.accumulated_dwell_time_current_lactation = 0
                self.dim_per_lactation[self.lactation] = self.current_dim_count
                self.current_dim_count = 0
                self.lactation += 1
                add_calf = True
            previous_pen.remove_cow(self)
            self.current_pen = next_pen
            self.time_remaining = next_pen.get_dwell_time()
            self.dwell_start_absolute_day = absolute_day
            self.current_pen.add_cow(self)
        return add_calf

    def calculate_infectious_period(self, current_day):
        if self.infection_time is not None and self.recovery_day is not None:
            return self.recovery_day - self.infection_time
        else:
            return None

    def get_calving_intervals(self):
        return self.calving_intervals

    def get_days_in_milk_data(self):
        return self.dim_per_lactation

    def get_max_lactation(self):
        return self.lactation

    def get_current_days_in_milk(self):
        return self.current_dim_count

    def get_completed_dwell_times(self):
        return self.completed_dwell_times

    def set_dwell_start_day(self, day):
        self.dwell_start_absolute_day = day

class Pen:
    def __init__(self, pen_id, merge_dwell_row):
        self.id = pen_id
        self.dwell_min = merge_dwell_row['min_duration']
        self.dwell_max = merge_dwell_row['max_duration']
        self.next_pen_id = merge_dwell_row['next_pen']
        self.cows_in_pen = []

    def add_cow(self, cow):
        if cow not in self.cows_in_pen:
            self.cows_in_pen.append(cow)

    def remove_cow(self, cow):
        if cow in self.cows_in_pen:
            self.cows_in_pen.remove(cow)

    def get_dwell_time(self):
        dwell_time = np.random.randint(low=self.dwell_min, high=self.dwell_max + 1)
        return max(1, dwell_time)

    def get_next_pen_id(self):
        return self.next_pen_id

class Farm:
    def __init__(self, merge_dwell_df, initial_infected, transmission_rate, recovery_rate, pens_to_mask=None):
        self.pens = {}
        self.cows = {}
        self.cow_counter = 1
        self._cow_history_data = []
        self._pen_counts_data = []
        self._movements_data = []
        self._daily_calf_additions_data = []
        self._days_in_milk_data = []
        self._completed_dwell_times_data = []
        self._calving_intervals_data = []
        self.annual_mortality_rate = 0.3
        self.initial_infected = initial_infected
        self.transmission_rate = transmission_rate
        self.recovery_rate = recovery_rate
        self.burn_in_days = None
        self._infectious_initialized = False
        self.pens_to_mask = pens_to_mask if pens_to_mask is not None else []
        for index, row in merge_dwell_df.iterrows():
            pen_id = row['current_pen']
            self.pens[pen_id] = Pen(pen_id, row)

    def add_cow(self, pen_id, lactation=0, state='Susceptible', current_day=0, incubation_period=0):
        pen = self.pens[pen_id]
        cow = Cow(self.cow_counter, pen, lactation, state=state, incubation_period=incubation_period)
        self.cows[self.cow_counter] = cow
        pen.add_cow(cow)
        self.cow_counter += 1
        cow.set_dwell_start_day(current_day)

    def step_day(self, absolute_day, burn_in_days):
        record = absolute_day > burn_in_days
        day = absolute_day - burn_in_days if record else absolute_day
        cow_ids = list(self.cows.keys())
        new_calves = []
        daily_mortality_prob = 1 - (1 - self.annual_mortality_rate)**(1/365)
        alive_cow_ids = [id for id in cow_ids if id in self.cows and self.cows[id].alive]
        num_alive_cows = len(alive_cow_ids)
        num_to_remove = np.random.binomial(num_alive_cows, daily_mortality_prob)
        if num_alive_cows > 0:
            cows_to_remove_ids = random.sample(alive_cow_ids, min(num_to_remove, num_alive_cows))
            for cow_id in cows_to_remove_ids:
                if cow_id in self.cows:
                    self.cows[cow_id].alive = False
                    self.cows[cow_id].current_pen.remove_cow(self.cows[cow_id])
        for id in list(self.cows.keys()):
            if id in self.cows and self.cows[id].alive:
                cow = self.cows[id]
                if cow.state == 'Exposed':
                    if cow.infection_time is not None and absolute_day - cow.infection_time >= cow.incubation_period:
                        cow.state = 'Infectious'
#                         print(f"New infection: Cow {cow.id} on absolute_day {absolute_day} (day {day})")
                        cow.recovery_day = absolute_day + int(np.random.geometric(self.recovery_rate))
                elif cow.state == 'Infectious':
                    if cow.recovery_day is not None and absolute_day >= cow.recovery_day:
                        cow.state = 'Recovered'
                        cow.recovery_day = None
                previous_pen_id = cow.current_pen.id
                add_calf = cow.step(absolute_day, self.pens, self)
                if cow.current_pen.id == "Pen1" and previous_pen_id != "Pen1" and cow.birth is None:
                    cow.birth = absolute_day
                if add_calf:
                    new_calves.append({'pen': self.pens["Pen1"]})
        for pen in self.pens.values():
            if pen.id not in self.pens_to_mask:
                infectious_cows_in_pen = [cow for cow in pen.cows_in_pen if cow.alive and hasattr(cow, 'state') and cow.state == 'Infectious']
                susceptible_cows_in_pen = [cow for cow in pen.cows_in_pen if cow.alive and hasattr(cow, 'state') and cow.state == 'Susceptible']
                if infectious_cows_in_pen and susceptible_cows_in_pen:
                    num_infectious = len(infectious_cows_in_pen)
                    num_susceptible = len(susceptible_cows_in_pen)
                    prob_exposure = 1 - (1 - self.transmission_rate)**num_infectious
                    num_to_expose = np.random.binomial(num_susceptible, prob_exposure)
                    if num_to_expose > 0:
                        cows_to_expose = random.sample(susceptible_cows_in_pen, num_to_expose)
                        for susc_cow in cows_to_expose:
                            if susc_cow.alive and hasattr(susc_cow, 'state') and susc_cow.state == 'Susceptible':
                                susc_cow.state = 'Exposed'
#                                 print("susceptible to new exposed", " ", susc_cow.id, " ", absolute_day, " ")
                                susc_cow.infection_time = absolute_day
        calves_added_count = 0
        for _ in new_calves:
            self.add_cow("Pen1", lactation=0, state='Susceptible', current_day=absolute_day, incubation_period=0)
            calves_added_count += 1
        if record:
            self._daily_calf_additions_data.append({'day': day, 'count': calves_added_count})
            for id in list(self.cows.keys()):
                if id in self.cows and self.cows[id].alive:
                    cow = self.cows[id]
                    if hasattr(cow, 'state'):
                        cow.record_history(absolute_day, burn_in_days)
                        self._cow_history_data.append({
                            'id': cow.id,
                            'day': day,
                            'absolute_day': absolute_day,
                            'pen': cow.current_pen.id,
                            'alive': cow.alive,
                            'lactation': cow.lactation,
                            'state': cow.state
                        })
            for pen in self.pens.values():
                alive_cows_in_pen = [cow for cow in pen.cows_in_pen if cow.alive]
                self._pen_counts_data.append({
                    'day': day,
                    'pen': pen.id,
                    'count': len(alive_cows_in_pen),
                    'simulation_run': None
                })

    def run_simulation(self, ndays, burn_in_days, incubation_period):
        self.burn_in_days = burn_in_days  # Store burn_in_days as instance variable
        self.initiate_infectious_cows(burn_in_days + 1, incubation_period)
        for day in range(1, ndays + 1):
            absolute_day = burn_in_days + day
            self.step_day(absolute_day, burn_in_days)
        self._collect_lactation_data_from_cows()
        self._collect_completed_dwell_times()
        self._collect_calving_intervals()

    def record_movement(self, cow_id, absolute_day, from_pen, to_pen):
        if absolute_day > self.burn_in_days:
            day = absolute_day - self.burn_in_days

    def initiate_infectious_cows(self, absolute_day, incubation_period=None):
        if self._infectious_initialized:
            print(f"Skipping initiate_infectious_cows: Already initialized on absolute_day {absolute_day}")
            return
        if absolute_day != self.burn_in_days + 1:
            return
        if "Pen8" not in self.pens:
            raise KeyError("Pen8 not found in farm.pens")
        current_infectious = sum(1 for cow in self.cows.values() if cow.alive and cow.state == 'Infectious')
        cows_to_add = max(0, self.initial_infected - current_infectious)

        infected_pen = self.pens["Pen8"]
        infected_ids = []
        for _ in range(cows_to_add):
            new_cow = Cow(self.cow_counter, pen=infected_pen, lactation=1, state="Infectious", incubation_period=incubation_period)
            new_cow.infection_time = absolute_day
            new_cow.recovery_day = absolute_day + int(np.random.geometric(self.recovery_rate))
            if hasattr(infected_pen, "get_dwell_time"):
                new_cow.time_remaining = infected_pen.get_dwell_time()
            new_cow.dwell_start_absolute_day = absolute_day
            self.cows[self.cow_counter] = new_cow
            infected_pen.add_cow(new_cow)
            infected_ids.append(self.cow_counter)
            self.cow_counter += 1
        self._infectious_initialized = True

    def _collect_lactation_data_from_cows(self):
        for cow_id, cow in self.cows.items():
            if cow.current_dim_count > 0:
                self._days_in_milk_data.append({
                    'cow_id': cow_id,
                    'lactation_number': cow.lactation,
                    'days_in_milk': cow.current_dim_count
                })

    def _collect_completed_dwell_times(self):
        for cow_id, cow in self.cows.items():
            if cow.completed_dwell_times:
                for dwell_data in cow.completed_dwell_times:
                    self._completed_dwell_times_data.append({
                        'cow_id': cow_id,
                        'pen': dwell_data['pen'],
                        'dwell_days': dwell_data['dwell_days'],
                        'start_day': dwell_data['start_day']
                    })

    def _collect_calving_intervals(self):
        for cow_id, cow in self.cows.items():
            if cow.calving_intervals:
                for interval_data in cow.calving_intervals:
                    self._calving_intervals_data.append({
                        'cow_id': cow_id,
                        'lactation_number': interval_data['lactation_number'],
                        'calv_interval': interval_data['calv_interval']
                    })

    def get_cow_history_data(self):
        return self._cow_history_data

    def get_pen_counts_data(self):
        return [data for data in self._pen_counts_data if data['pen'] not in self.pens_to_mask]

    def get_pen_counts_masked(self):
        return [data for data in self._pen_counts_data if data['pen'] in self.pens_to_mask]

    def get_movements_data(self):
        return self._movements_data

    def get_daily_calf_additions_data(self):
        return self._daily_calf_additions_data

    def get_days_in_milk_data(self):
        return self._days_in_milk_data

    def get_completed_dwell_times_data(self):
        return self._completed_dwell_times_data

    def get_calving_intervals_data(self):
        return self._calving_intervals_data

    def get_all_cow_lactation_data(self):
        all_cow_lactation_summary = []
        for cow_id, cow in self.cows.items():
            days_in_milk_data = cow.get_days_in_milk_data()
            max_lactation = cow.get_max_lactation()
            all_cow_lactation_summary.append({
                'cow_id': cow_id,
                'max_lactation': max_lactation,
                'days_in_milk_per_lactation': days_in_milk_data
            })
        return all_cow_lactation_summary


def calculate_poisson_loss(simulated_counts, actual_counts):
    if len(simulated_counts) != len(actual_counts):
        raise ValueError("Simulated and actual counts must have the same length.")

    epsilon = 1e-9     # Add a small epsilon to avoid log(0) and log factorial of 0
    log_factorial_actual = gammaln(actual_counts + 1)

    terms = simulated_counts - actual_counts * np.log(simulated_counts + epsilon) + log_factorial_actual

    weights = actual_counts + 1  # Time-varying weights proportional to actual counts +1 to avoid zero

    return np.sum(weights * terms) / np.sum(weights)

def plf(a,b):
    if len(a) != len(b):
        raise ValueError("inputs must have the same length.")
    epsilon = 1e-9
    log_fact = gammaln(b + 1)
    terms = a - b * np.log(a + epsilon) + log_fact
    weights = b + 1 #time varying weights proportional to b + 1 (avoids zero)
    return np.sum(weights * terms)/np.sum(weights)

def run_main_simulation(burn_in_farm_state, merge_dwell_df, ndays, burn_in_days, initial_infected, transmission_rate, recovery_rate, incubation_period=None, actual_infectious_counts=None, pens_to_mask=None, num_simulations=1, plot=False):
    
    if not (0.001 <= transmission_rate <= 0.999):
        raise ValueError(f"Invalid transmission_rate {transmission_rate}. Must be between 0.001 and 0.999.")
    if not (0.001 <= recovery_rate <= 0.999):
        raise ValueError(f"Invalid recovery_rate {recovery_rate}. Must be between 0.001 and 0.999.")
    if transmission_rate is None or recovery_rate is None or incubation_period is None:
        print("Optimal transmission rate, recovery rate, or incubation period not provided. Cannot run simulation.")
        raise ValueError("No parameters to simulate detected.")

#     # Debug: Check infectious cows in burn-in state
#     infectious_cows = [cow for cow in burn_in_farm_state.values() if cow.state == 'Infectious']
#     if len(infectious_cows) > 0:
#         raise ValueError(f"Infectious cows in burn_in_state detected: {len(infectious_cows)}")

    all_simulated_seir_histories = []
    simulated_pen_counts = []
    simulated_daily_calf_additions = []
    simulated_completed_dwell_times = []
    simulated_calving_intervals = []
    simulated_lactation_summary = []
    simulated_days_in_milk = []

    for _ in range(num_simulations):
        # Create a new farm instance for the simulation
        farm = Farm(
            merge_dwell_df,
            initial_infected=initial_infected,
            transmission_rate=transmission_rate,
            recovery_rate=recovery_rate,
            pens_to_mask=pens_to_mask
        )
        start_absolute_day_main_sim = burn_in_days

        # Add cows from the burn-in state, resetting infection-related attributes
        for cow_id, cow in burn_in_farm_state.items():
            new_cow = Cow(
                cow.id,
                farm.pens[cow.current_pen.id],
                cow.lactation,
                state='Susceptible',  # Reset to Susceptible to clear burn-in infections
                incubation_period=incubation_period
            )
            new_cow.time_remaining = cow.time_remaining
            new_cow.birth = cow.birth
            new_cow.dim_per_lactation = cow.dim_per_lactation.copy()
            new_cow.current_dim_count = cow.current_dim_count
            new_cow.dwell_start_absolute_day = start_absolute_day_main_sim
            new_cow.completed_dwell_times = cow.completed_dwell_times.copy()
            new_cow.accumulated_dwell_time_current_lactation = cow.accumulated_dwell_time_current_lactation
            new_cow.calving_intervals = cow.calving_intervals.copy()

            if new_cow.current_pen.id in farm.pens:
                farm.cows[cow_id] = new_cow
                farm.pens[new_cow.current_pen.id].add_cow(new_cow)
            else:
                print(f"Warning: Pen {new_cow.current_pen.id} not found in the new farm instance.")

        farm.cow_counter = max(farm.cows.keys()) + 1 if farm.cows else 1

        # Run the simulation
        farm.run_simulation(ndays, burn_in_days, incubation_period)

        # Collect data from the simulation
        all_simulated_seir_histories.append(farm.get_cow_history_data())
        simulated_pen_counts.append(farm.get_pen_counts_data())
        simulated_daily_calf_additions.append(farm.get_daily_calf_additions_data())
        simulated_completed_dwell_times.append(farm.get_completed_dwell_times_data())
        simulated_calving_intervals.append(farm.get_calving_intervals_data())
        simulated_lactation_summary.append(farm.get_all_cow_lactation_data())

        alive_cows_at_end = [cow for cow in farm.cows.values() if cow.alive]
        final_days_in_milk = [{
            'cow_id': cow.id,
            'lactation_number': cow.lactation,
            'days_in_milk': cow.get_current_days_in_milk(),
        } for cow in alive_cows_at_end]
        simulated_days_in_milk.append(final_days_in_milk)

    # Convert to DataFrames
    simulated_seir_df = [pd.DataFrame(hist) if hist else pd.DataFrame() for hist in all_simulated_seir_histories]
    simulated_pen_counts_df = [pd.DataFrame(counts) if counts else pd.DataFrame() for counts in simulated_pen_counts]
    simulated_daily_calf_additions_df = [pd.DataFrame(calves) if calves else pd.DataFrame() for calves in simulated_daily_calf_additions]
    simulated_dwell_times_df = [pd.DataFrame(dwells) if dwells else pd.DataFrame() for dwells in simulated_completed_dwell_times]
    simulated_calv_interval_df = [pd.DataFrame(intervals) if intervals else pd.DataFrame() for intervals in simulated_calving_intervals]
    simulated_lactation_summary_df = [pd.DataFrame(summary) if summary else pd.DataFrame() for summary in simulated_lactation_summary]
    simulated_days_in_milk_df = [pd.DataFrame(days) if days else pd.DataFrame() for days in simulated_days_in_milk]

    # Plotting
    if plot and num_simulations > 0:
        start_absolute_day_main_sim = burn_in_days
        days_for_plotting = range(1, ndays + 1)

        # Mean and 90% CI Plot for Infectious Counts
        all_infectious_counts_np = []
        for seir_df in simulated_seir_df:
            if not seir_df.empty:
                simulated_infectious_df = seir_df[
                    (seir_df['state'] == 'Infectious') &
                    (seir_df['day'] > 1) &
                    (seir_df['day'] <= ndays) &
                    (~seir_df['pen'].isin(pens_to_mask))
                ].copy()
                simulated_infectious_df = simulated_infectious_df.groupby('id')['day'].min().reset_index()
                simulated_daily_infectious_counts = simulated_infectious_df.groupby('day').size().reindex(
                    days_for_plotting, fill_value=0
                ).values
                all_infectious_counts_np.append(simulated_daily_infectious_counts)

        if all_infectious_counts_np:
            print("\nPlotting Mean Simulated Infectious Counts with 90% CI vs Actual...")
            mean_simulated_counts = np.mean(all_infectious_counts_np, axis=0)
            lower_bound = np.percentile(all_infectious_counts_np, 5, axis=0)  # 90% CI
            upper_bound = np.percentile(all_infectious_counts_np, 95, axis=0)

            plt.figure(figsize=(12, 6))
            plt.plot(days_for_plotting, mean_simulated_counts, label='Mean Simulated Infectious', color='blue', linewidth=2)
            plt.fill_between(days_for_plotting, lower_bound, upper_bound, color='blue', alpha=0.2, label='90% Confidence Interval')
            if actual_infectious_counts is not None and len(actual_infectious_counts) == len(days_for_plotting):
                plt.plot(days_for_plotting, actual_infectious_counts, label='Actual Infectious', color='red', linestyle='--', linewidth=2)
            plt.title(f'Mean Simulated vs Actual Infectious Cow Counts (with 90% CI, {num_simulations} Runs, Excluding Masked Pens)')
            plt.xlabel('Day')
            plt.ylabel('Number of Infectious Cows')
            plt.legend()
            plt.grid(True)
            plt.show()
            
        # Aggregate Pen Populations Across All Runs
        if any(not df.empty for df in simulated_pen_counts_df):
            print("\nPlotting Mean Pen Populations Across All Runs (Excluding Masked Pens)...")
            all_pen_counts = pd.concat([df.assign(run=i) for i, df in enumerate(simulated_pen_counts_df) if not df.empty], ignore_index=True)
            if not all_pen_counts.empty:
                mean_pen_counts = all_pen_counts.groupby(['day', 'pen'])['count'].mean().reset_index()
                plt.figure(figsize=(12, 6))
                sns.lineplot(data=mean_pen_counts, x='day', y='count', hue='pen')
                plt.title(f'Mean Number of Cows per Pen per Day ({num_simulations} Runs)')
                plt.xlabel('Day')
                plt.ylabel('Mean Number of Cows')
                plt.grid(True)
                plt.show()
            else:
                print("\nCannot plot Pen Populations: Aggregated Pen Histories data is empty.")

        # Aggregate SEIR States Across All Runs
        if any(not df.empty for df in simulated_seir_df):
            print("\nPlotting Mean SEIR States over Time with 90% CI Across All Runs (Excluding Masked Pens)...")
            all_seir_counts = []
            for seir_df in simulated_seir_df:
                if not seir_df.empty:
                    seir_counts = seir_df[
                        (~seir_df['pen'].isin(pens_to_mask)) &
                        (seir_df['day'] >= 1) &
                        (seir_df['day'] <= ndays)
                    ].groupby(['day', 'state']).size().reset_index(name='count')
                    # Pivot to ensure all states are present for each day
                    seir_counts_pivot = seir_counts.pivot(index='day', columns='state', values='count').reindex(
                        days_for_plotting, fill_value=0
                    ).fillna(0)
                    all_seir_counts.append(seir_counts_pivot)
            if all_seir_counts:
                # Concatenate and compute mean and percentiles
                all_seir_counts_df = pd.concat(all_seir_counts, axis=0, ignore_index=True)
                mean_seir_counts = all_seir_counts_df.groupby('day').mean().reset_index()
                lower_seir_counts = all_seir_counts_df.groupby('day').apply(lambda x: np.percentile(x, 5, axis=0)).reset_index()
                upper_seir_counts = all_seir_counts_df.groupby('day').apply(lambda x: np.percentile(x, 95, axis=0)).reset_index()
                
                plt.figure(figsize=(12, 6))
                colors = {'Susceptible': 'green', 'Exposed': 'orange', 'Infectious': 'red', 'Recovered': 'blue'}
                for state in ['Susceptible', 'Exposed', 'Infectious', 'Recovered']:
                    if state in mean_seir_counts.columns:
                        plt.plot(days_for_plotting, mean_seir_counts[state], label=f'Mean {state}', color=colors.get(state, 'black'), linewidth=2)
                        plt.fill_between(
                            days_for_plotting,
                            lower_seir_counts[state],
                            upper_seir_counts[state],
                            color=colors.get(state, 'black'),
                            alpha=0.2,
                            label=f'{state} 90% CI'
                        )
                plt.title(f'Mean SEIR States with 90% CI ({num_simulations} Runs, Excluding Masked Pens)')
                plt.xlabel('Day')
                plt.ylabel('Mean Number of Cows')
                plt.legend()
                plt.grid(True)
                plt.show()
            else:
                print("\nCannot plot SEIR States: Aggregated SEIR Histories data is empty.")

        # Aggregate Calving Intervals Across All Runs
        if any(not df.empty for df in simulated_calv_interval_df):
            print("\nPlotting Boxplot of Mean Calving Intervals by Lactation Across All Runs...")
            all_calv_intervals = pd.concat([df for df in simulated_calv_interval_df if not df.empty], ignore_index=True)
            if not all_calv_intervals.empty:
                all_calv_intervals_filtered = all_calv_intervals[all_calv_intervals['calv_interval'] >= 200].copy()
                if not all_calv_intervals_filtered.empty:
                    mean_calv_intervals = all_calv_intervals_filtered.groupby(['cow_id', 'lactation_number'])['calv_interval'].mean().reset_index()
                    plt.figure(figsize=(10, 6))
                    sns.boxplot(data=mean_calv_intervals, x='lactation_number', y='calv_interval')
                    plt.title(f'Boxplot of Mean Calving Intervals by Lactation Number ({num_simulations} Runs)')
                    plt.xlabel('Lactation Number')
                    plt.ylabel('Mean Calving Interval (days)')
                    plt.grid(True)
                    plt.show()
                else:
                    print("No filtered calving interval data (>= 200 days) to plot.")
            else:
                print("\nCannot plot Calving Intervals: Aggregated Calving Intervals data is empty.")

        # Aggregate Days in Milk Across All Runs
        if any(not df.empty for df in simulated_days_in_milk_df):
            print("\nPlotting Boxplot of Mean Days in Milk by Lactation Across All Runs...")
            all_days_in_milk = pd.concat([df for df in simulated_days_in_milk_df if not df.empty], ignore_index=True)
            if not all_days_in_milk.empty:
                mean_days_in_milk = all_days_in_milk.groupby(['cow_id', 'lactation_number'])['days_in_milk'].mean().reset_index()
                plt.figure(figsize=(10, 6))
                sns.boxplot(data=mean_days_in_milk, x='lactation_number', y='days_in_milk')
                plt.title(f'Boxplot of Mean Days in Milk by Lactation ({num_simulations} Runs)')
                plt.xlabel('Lactation')
                plt.ylabel('Mean Days in Milk')
                plt.grid(True)
                plt.show()
            else:
                print("\nCannot plot Days in Milk: Aggregated Days in Milk data is empty.")

#     print("Simulation(s) completed.")
    
    return {
            'seir_history': simulated_seir_df if num_simulations > 1 else simulated_seir_df[0],
            'pen_counts': simulated_pen_counts_df if num_simulations > 1 else simulated_pen_counts_df[0],
            'daily_calf_additions': simulated_daily_calf_additions_df if num_simulations > 1 else simulated_daily_calf_additions_df[0],
            'completed_dwell_times': simulated_dwell_times_df if num_simulations > 1 else simulated_dwell_times_df[0],
            'calving_intervals': simulated_calv_interval_df if num_simulations > 1 else simulated_calv_interval_df[0],
            'lactation_summary': simulated_lactation_summary_df if num_simulations > 1 else simulated_lactation_summary_df[0],
            'days_in_milk': simulated_days_in_milk_df if num_simulations > 1 else simulated_days_in_milk_df[0]
        }

def run_main_simulation2(burn_in_farm_state, merge_dwell_df, ndays, burn_in_days, initial_infected, transmission_rate, recovery_rate, incubation_period=None, actual_infectious_counts=None, pens_to_mask=None, num_simulations=1, plot=False):
    
    if not (0.001 <= transmission_rate <= 0.999):
        raise ValueError(f"Invalid transmission_rate {transmission_rate}. Must be between 0.001 and 0.999.")
    if not (0.001 <= recovery_rate <= 0.999):
        raise ValueError(f"Invalid recovery_rate {recovery_rate}. Must be between 0.001 and 0.999.")
    if transmission_rate is None or recovery_rate is None or incubation_period is None:
        print("Optimal transmission rate, recovery rate, or incubation period not provided. Cannot run simulation.")
        raise ValueError("No parameters to simulate detected.")

#     # Debug: Check infectious cows in burn-in state
#     infectious_cows = [cow for cow in burn_in_farm_state.values() if cow.state == 'Infectious']
#     if len(infectious_cows) > 0:
#         raise ValueError(f"Infectious cows in burn_in_state detected: {len(infectious_cows)}")

    all_simulated_seir_histories = []
    simulated_pen_counts = []
    simulated_daily_calf_additions = []
    simulated_completed_dwell_times = []
    simulated_calving_intervals = []
    simulated_lactation_summary = []
    simulated_days_in_milk = []

    for _ in range(num_simulations):
        # Create a new farm instance for the simulation
        farm = Farm(
            merge_dwell_df,
            initial_infected=initial_infected,
            transmission_rate=transmission_rate,
            recovery_rate=recovery_rate,
            pens_to_mask=pens_to_mask
        )
        start_absolute_day_main_sim = burn_in_days

        # Add cows from the burn-in state, resetting infection-related attributes
        for cow_id, cow in burn_in_farm_state.cows.items():
            new_cow = Cow(
                cow.id,
                farm.pens[cow.current_pen.id],
                cow.lactation,
                state='Susceptible',  # Reset to Susceptible to clear burn-in infections
                incubation_period=incubation_period
            )
            new_cow.time_remaining = cow.time_remaining
            new_cow.birth = cow.birth
            new_cow.dim_per_lactation = cow.dim_per_lactation.copy()
            new_cow.current_dim_count = cow.current_dim_count
            new_cow.dwell_start_absolute_day = start_absolute_day_main_sim
            new_cow.completed_dwell_times = cow.completed_dwell_times.copy()
            new_cow.accumulated_dwell_time_current_lactation = cow.accumulated_dwell_time_current_lactation
            new_cow.calving_intervals = cow.calving_intervals.copy()

            if new_cow.current_pen.id in farm.pens:
                farm.cows[cow_id] = new_cow
                farm.pens[new_cow.current_pen.id].add_cow(new_cow)
            else:
                print(f"Warning: Pen {new_cow.current_pen.id} not found in the new farm instance.")

        farm.cow_counter = max(farm.cows.keys()) + 1 if farm.cows else 1

        # Run the simulation
        farm.run_simulation(ndays, burn_in_days, incubation_period)

        # Collect data from the simulation
        all_simulated_seir_histories.append(farm.get_cow_history_data())
        simulated_pen_counts.append(farm.get_pen_counts_data())
        simulated_daily_calf_additions.append(farm.get_daily_calf_additions_data())
        simulated_completed_dwell_times.append(farm.get_completed_dwell_times_data())
        simulated_calving_intervals.append(farm.get_calving_intervals_data())
        simulated_lactation_summary.append(farm.get_all_cow_lactation_data())

        alive_cows_at_end = [cow for cow in farm.cows.values() if cow.alive]
        final_days_in_milk = [{
            'cow_id': cow.id,
            'lactation_number': cow.lactation,
            'days_in_milk': cow.get_current_days_in_milk(),
        } for cow in alive_cows_at_end]
        simulated_days_in_milk.append(final_days_in_milk)

    # Convert to DataFrames
    simulated_seir_df = [pd.DataFrame(hist) if hist else pd.DataFrame() for hist in all_simulated_seir_histories]
    simulated_pen_counts_df = [pd.DataFrame(counts) if counts else pd.DataFrame() for counts in simulated_pen_counts]
    simulated_daily_calf_additions_df = [pd.DataFrame(calves) if calves else pd.DataFrame() for calves in simulated_daily_calf_additions]
    simulated_dwell_times_df = [pd.DataFrame(dwells) if dwells else pd.DataFrame() for dwells in simulated_completed_dwell_times]
    simulated_calv_interval_df = [pd.DataFrame(intervals) if intervals else pd.DataFrame() for intervals in simulated_calving_intervals]
    simulated_lactation_summary_df = [pd.DataFrame(summary) if summary else pd.DataFrame() for summary in simulated_lactation_summary]
    simulated_days_in_milk_df = [pd.DataFrame(days) if days else pd.DataFrame() for days in simulated_days_in_milk]

    # Plotting
    if plot and num_simulations > 0:
        start_absolute_day_main_sim = burn_in_days
        days_for_plotting = range(1, ndays + 1)

        # Mean and 90% CI Plot for Infectious Counts
        all_infectious_counts_np = []
        for seir_df in simulated_seir_df:
            if not seir_df.empty:
                simulated_infectious_df = seir_df[
                    (seir_df['state'] == 'Infectious') &
                    (seir_df['day'] > 1) &
                    (seir_df['day'] <= ndays) &
                    (~seir_df['pen'].isin(pens_to_mask))
                ].copy()
                simulated_infectious_df = simulated_infectious_df.groupby('id')['day'].min().reset_index()
                simulated_daily_infectious_counts = simulated_infectious_df.groupby('day').size().reindex(
                    days_for_plotting, fill_value=0
                ).values
                all_infectious_counts_np.append(simulated_daily_infectious_counts)

        if all_infectious_counts_np:
            print("\nPlotting Mean Simulated Infectious Counts with 90% CI vs Actual...")
            mean_simulated_counts = np.mean(all_infectious_counts_np, axis=0)
            lower_bound = np.percentile(all_infectious_counts_np, 5, axis=0)  # 90% CI
            upper_bound = np.percentile(all_infectious_counts_np, 95, axis=0)

            plt.figure(figsize=(12, 6))
            plt.plot(days_for_plotting, mean_simulated_counts, label='Mean Simulated Infectious', color='blue', linewidth=2)
            plt.fill_between(days_for_plotting, lower_bound, upper_bound, color='blue', alpha=0.2, label='90% Confidence Interval')
            if actual_infectious_counts is not None and len(actual_infectious_counts) == len(days_for_plotting):
                plt.plot(days_for_plotting, actual_infectious_counts, label='Actual Infectious', color='red', linestyle='--', linewidth=2)
            plt.title(f'Mean Simulated vs Actual Infectious Cow Counts (with 90% CI, {num_simulations} Runs, Excluding Masked Pens)')
            plt.xlabel('Day')
            plt.ylabel('Number of Infectious Cows')
            plt.legend()
            plt.grid(True)
            plt.show()
            
        # Aggregate Pen Populations Across All Runs
        if any(not df.empty for df in simulated_pen_counts_df):
            print("\nPlotting Mean Pen Populations Across All Runs (Excluding Masked Pens)...")
            all_pen_counts = pd.concat([df.assign(run=i) for i, df in enumerate(simulated_pen_counts_df) if not df.empty], ignore_index=True)
            if not all_pen_counts.empty:
                mean_pen_counts = all_pen_counts.groupby(['day', 'pen'])['count'].mean().reset_index()
                plt.figure(figsize=(12, 6))
                sns.lineplot(data=mean_pen_counts, x='day', y='count', hue='pen')
                plt.title(f'Mean Number of Cows per Pen per Day ({num_simulations} Runs)')
                plt.xlabel('Day')
                plt.ylabel('Mean Number of Cows')
                plt.grid(True)
                plt.show()
            else:
                print("\nCannot plot Pen Populations: Aggregated Pen Histories data is empty.")

        # Aggregate SEIR States Across All Runs
        if any(not df.empty for df in simulated_seir_df):
            print("\nPlotting Mean SEIR States over Time with 90% CI Across All Runs (Excluding Masked Pens)...")
            all_seir_counts = []
            for seir_df in simulated_seir_df:
                if not seir_df.empty:
                    seir_counts = seir_df[
                        (~seir_df['pen'].isin(pens_to_mask)) &
                        (seir_df['day'] >= 1) &
                        (seir_df['day'] <= ndays)
                    ].groupby(['day', 'state']).size().reset_index(name='count')
                    # Pivot to ensure all states are present for each day
                    seir_counts_pivot = seir_counts.pivot(index='day', columns='state', values='count').reindex(
                        days_for_plotting, fill_value=0
                    ).fillna(0)
                    all_seir_counts.append(seir_counts_pivot)
            if all_seir_counts:
                # Concatenate and compute mean and percentiles
                all_seir_counts_df = pd.concat(all_seir_counts, axis=0, ignore_index=True)
                mean_seir_counts = all_seir_counts_df.groupby('day').mean().reset_index()
                lower_seir_counts = all_seir_counts_df.groupby('day').apply(lambda x: np.percentile(x, 5, axis=0)).reset_index()
                upper_seir_counts = all_seir_counts_df.groupby('day').apply(lambda x: np.percentile(x, 95, axis=0)).reset_index()
                
                plt.figure(figsize=(12, 6))
                colors = {'Susceptible': 'green', 'Exposed': 'orange', 'Infectious': 'red', 'Recovered': 'blue'}
                for state in ['Susceptible', 'Exposed', 'Infectious', 'Recovered']:
                    if state in mean_seir_counts.columns:
                        plt.plot(days_for_plotting, mean_seir_counts[state], label=f'Mean {state}', color=colors.get(state, 'black'), linewidth=2)
                        plt.fill_between(
                            days_for_plotting,
                            lower_seir_counts[state],
                            upper_seir_counts[state],
                            color=colors.get(state, 'black'),
                            alpha=0.2,
                            label=f'{state} 90% CI'
                        )
                plt.title(f'Mean SEIR States with 90% CI ({num_simulations} Runs, Excluding Masked Pens)')
                plt.xlabel('Day')
                plt.ylabel('Mean Number of Cows')
                plt.legend()
                plt.grid(True)
                plt.show()
            else:
                print("\nCannot plot SEIR States: Aggregated SEIR Histories data is empty.")

        # Aggregate Calving Intervals Across All Runs
        if any(not df.empty for df in simulated_calv_interval_df):
            print("\nPlotting Boxplot of Mean Calving Intervals by Lactation Across All Runs...")
            all_calv_intervals = pd.concat([df for df in simulated_calv_interval_df if not df.empty], ignore_index=True)
            if not all_calv_intervals.empty:
                all_calv_intervals_filtered = all_calv_intervals[all_calv_intervals['calv_interval'] >= 200].copy()
                if not all_calv_intervals_filtered.empty:
                    mean_calv_intervals = all_calv_intervals_filtered.groupby(['cow_id', 'lactation_number'])['calv_interval'].mean().reset_index()
                    plt.figure(figsize=(10, 6))
                    sns.boxplot(data=mean_calv_intervals, x='lactation_number', y='calv_interval')
                    plt.title(f'Boxplot of Mean Calving Intervals by Lactation Number ({num_simulations} Runs)')
                    plt.xlabel('Lactation Number')
                    plt.ylabel('Mean Calving Interval (days)')
                    plt.grid(True)
                    plt.show()
                else:
                    print("No filtered calving interval data (>= 200 days) to plot.")
            else:
                print("\nCannot plot Calving Intervals: Aggregated Calving Intervals data is empty.")

        # Aggregate Days in Milk Across All Runs
        if any(not df.empty for df in simulated_days_in_milk_df):
            print("\nPlotting Boxplot of Mean Days in Milk by Lactation Across All Runs...")
            all_days_in_milk = pd.concat([df for df in simulated_days_in_milk_df if not df.empty], ignore_index=True)
            if not all_days_in_milk.empty:
                mean_days_in_milk = all_days_in_milk.groupby(['cow_id', 'lactation_number'])['days_in_milk'].mean().reset_index()
                plt.figure(figsize=(10, 6))
                sns.boxplot(data=mean_days_in_milk, x='lactation_number', y='days_in_milk')
                plt.title(f'Boxplot of Mean Days in Milk by Lactation ({num_simulations} Runs)')
                plt.xlabel('Lactation')
                plt.ylabel('Mean Days in Milk')
                plt.grid(True)
                plt.show()
            else:
                print("\nCannot plot Days in Milk: Aggregated Days in Milk data is empty.")

#     print("Simulation(s) completed.")
    
    return {
            'seir_history': simulated_seir_df if num_simulations > 1 else simulated_seir_df[0],
            'pen_counts': simulated_pen_counts_df if num_simulations > 1 else simulated_pen_counts_df[0],
            'daily_calf_additions': simulated_daily_calf_additions_df if num_simulations > 1 else simulated_daily_calf_additions_df[0],
            'completed_dwell_times': simulated_dwell_times_df if num_simulations > 1 else simulated_dwell_times_df[0],
            'calving_intervals': simulated_calv_interval_df if num_simulations > 1 else simulated_calv_interval_df[0],
            'lactation_summary': simulated_lactation_summary_df if num_simulations > 1 else simulated_lactation_summary_df[0],
            'days_in_milk': simulated_days_in_milk_df if num_simulations > 1 else simulated_days_in_milk_df[0]
        }

def objective_function_metric(params, burn_in_farm_case, ndays_target_case, burn_in_days, target_curve_case, no_calves, metric = "distance"):

    T, R, I, IC = params
    I = int(round(I)) # gp_minimize works with floats, convert to int
    IC = int(round(IC)) # gp_minimize works with floats, convert to int

    # Ensure seed is at least 1
    IC = max(1, IC)

    sim_curve = simulate_counts(
        burn_in_farm_case,
        ndays_target_case,
        burn_in_days,
        IC,
        float(T),
        float(R),
        I,
        pens_to_mask=no_calves
    )

    if metric == "distance":
        d_dtw = dtw_distance(sim_curve, target_curve_case)
        d_mae = np.mean(np.abs(sim_curve - target_curve_case))
        d_cc, _ = crosscorr_alignment_distance(sim_curve, target_curve_case, max_shift=ndays_target_case//2)
        distance = d_dtw + 5.0 * d_mae + 2.0 * d_cc
    else:
        distance = plf(sim_curve, target_curve_case)
    return distance


def simulate_counts(burn_in_state, ndays, burn_in_days, init_infected, T, R,
                    incubation_period=4, pens_to_mask=None):
    # Convert farm → dict if needed
    if hasattr(burn_in_state, "cows"):
        burn_in_state = initial_cows

    res = run_main_simulation(burn_in_state,
                                 merge_dwell,
                                 ndays, burn_in_days,
                                 init_infected, T, R,
                                 incubation_period, 
                              pens_to_mask = ['Pen1','Pen2','Pen3','Pen4'])

    seir_df = res['seir_history']
    start_abs = burn_in_days
    if seir_df.empty:
        return np.zeros(ndays, dtype=int)

    # choose correct cow identifier column
    cow_col = None
    for candidate in ["cow_id", "id", "animal_id"]:
        if candidate in seir_df.columns:
            cow_col = candidate
            break
    if cow_col is None:
        raise KeyError(f"No cow ID column found in seir_history. Available: {list(seir_df.columns)}")

    # filter to first entry into Infectious
    df = seir_df[(seir_df['state'] == 'Infectious') &
                 (seir_df['absolute_day'] > start_abs) &
                 (seir_df['absolute_day'] <= start_abs + ndays)]
    if df.empty:
        return np.zeros(ndays, dtype=int)

    df_first = df.groupby(cow_col)['absolute_day'].min().reset_index()
    counts = df_first.groupby('absolute_day').size().reindex(
        range(start_abs+1, start_abs+ndays+1), fill_value=0).values

    return counts.astype(int)

def dtw_distance(a, b):
    n, m = len(a), len(b)
    D = np.full((n+1, m+1), np.inf)
    D[0,0] = 0.0
    for i in range(1, n+1):
        for j in range(1, m+1):
            cost = abs(float(a[i-1]) - float(b[j-1]))
            D[i,j] = cost + min(D[i-1,j], D[i,j-1], D[i-1,j-1])
    return D[n,m]

# cross-correlation-based alignment: shift sim to align with obs, return minimal MAE across shifts
def crosscorr_alignment_distance(sim, obs, max_shift=None):
    # compute best shift of sim relative to obs to minimize MAE
    n = len(obs)
    if max_shift is None:
        max_shift = n // 2
    best_mae = np.mean(np.abs(sim - obs))  # zero shift baseline
    best_shift = 0
    for shift in range(-max_shift, max_shift+1):
        if shift == 0:
            continue
        # shift positive => sim delayed -> compare sim[shift:] with obs[:-shift]
        if shift > 0:
            s = sim[shift:]
            o = obs[:len(s)]
        else:
            s = sim[:len(sim)+shift]
            o = obs[-shift:len(obs)]
        if len(s) == 0:
            continue
        mae = np.mean(np.abs(s - o))
        if mae < best_mae:
            best_mae = mae
            best_shift = shift
    return best_mae, best_shift

# additional summary features for ML (for classifier)
def extract_features(curve):
    # raw stats + normalized peak/day + cumulative etc.
    curve = np.asarray(curve, dtype=float)
    total = curve.sum()
    peak = curve.max()
    peak_day = np.argmax(curve) + 1
    mean = curve.mean()
    std = curve.std()
    skew = 0.0
    if std > 0:
        skew = ((np.mean((curve-mean)**3)) / (std**3)) if std>0 else 0.0
    # normalized features
    feat = [total, peak, peak_day, mean, std, skew]
    return np.array(feat, dtype=float)

def posterior_selection(cases, case_inputs, ABC_N_SIMS, T_low, T_high, R_low, R_high, INCUBATION_LOW, INCUBATION_HIGH, SEED_LOW, SEED_HIGH, ABC_ACCEPT_TOP_K, N_post_sims, metric='distance'):
    if metric not in ['distance', 'plf']:
        raise ValueError("Metric must be either 'distance' or 'plf'")

    # Initialize dictionaries to store DataFrames and records
    abc_dfs = {}
    accepted_dfs = {}
    abc_records_by_case = {case: [] for case in cases}

    # Loop over cases to generate ABC records
    for case in cases:
        inputs = case_inputs[case]
        abc_records = []  # List of dicts: will hold T, R, incubation, seed, distance, plf, sim_curve
        
        for i in tqdm(range(ABC_N_SIMS), desc=f"ABC Prior Sampling for case {case}"):
            # Sample priors
            T = float(np.random.uniform(T_low, T_high))
            R = float(np.random.uniform(R_low, R_high))
            incubation = int(np.random.randint(INCUBATION_LOW, INCUBATION_HIGH + 1))
            seed = int(np.random.randint(SEED_LOW, SEED_HIGH + 1))
            no_calves = inputs['no_calves']
            # Run sim with those sampled params
            sim_curve = simulate_counts(
                burn_in_state=inputs['burn_in_farm'],  # Pass Farm object
                ndays=inputs["ndays_target"],
                burn_in_days=inputs["burn_in_days"],
                init_infected=seed,  # Integer for initial infected
                T=T,
                R=R,
                incubation_period=incubation,
                pens_to_mask=inputs["no_calves"]
            )
            
            # Distance: combine DTW, MAE, and CC
            d_dtw = dtw_distance(sim_curve, inputs["target_curve"])
            d_mae = np.mean(np.abs(sim_curve - inputs["target_curve"]))
            d_cc, best_shift = crosscorr_alignment_distance(
                sim_curve, inputs["target_curve"], max_shift=inputs["ndays_target"] // 2
            )
            distance = d_dtw + 5.0 * d_mae + 2.0 * d_cc
            
            # PLF metric
            plf_value = plf(inputs["target_curve"], sim_curve)
            
            abc_records.append({
                "T": T, 
                "R": R,
                "I": incubation, 
                "IC": seed,
                "distance": distance,
                "plf": plf_value,
                "sim_curve": sim_curve
            })
        
        # Store results for this case
        abc_records_by_case[case] = abc_records
        
        # Build DataFrame directly with sim_curves stored separately
        abc_df = pd.DataFrame(
            [{k: v for k, v in r.items() if k != "sim_curve"} for r in abc_records]
        )
        abc_df["sim_curve"] = [r["sim_curve"] for r in abc_records]
        abc_dfs[case] = abc_df

        # Sort by selected metric and keep top K
        accepted_df = abc_df.sort_values(metric, ascending=True).head(int(ABC_ACCEPT_TOP_K)).reset_index(drop=True)
        accepted_dfs[case] = accepted_df

    # Plot scatter plots for all cases
    fig, axes = plt.subplots(1, len(cases), figsize=(6 * len(cases), 5))
    for i, case in enumerate(cases):
        ax = axes[i] if len(cases) > 1 else axes
        abc_df = abc_dfs[case]
        accepted_df = accepted_dfs[case]

        # All prior samples, colored by selected metric
        scatter = ax.scatter(
            abc_df["T"], abc_df["R"],
            c=abc_df[metric], cmap="viridis",
            s=6, alpha=0.6, label="Prior samples"
        )

        # Accepted top K in red
        ax.scatter(
            accepted_df["T"], accepted_df["R"],
            c="red", s=12, label="Accepted top K"
        )

        fig.colorbar(scatter, ax=ax, label=metric)
        ax.set_xlabel("Transmission rate")
        ax.set_ylabel("Recovery rate")
        ax.set_title(f"ABC posterior scatter ({case})")
        ax.legend()

    plt.tight_layout()
    plt.show()
    plt.close()

    # Plot histograms for accepted parameters
    for case in cases:
        accepted_df = accepted_dfs[case]
        plt.figure(figsize=(12, 6))

        plt.subplot(2, 2, 1)
        plt.hist(accepted_df['T'], bins=20)
        plt.title(f"T (accepted) for case {case}")

        plt.subplot(2, 2, 2)
        plt.hist(accepted_df['R'], bins=20)
        plt.title(f"R (accepted) for case {case}")

        plt.subplot(2, 2, 3)
        plt.hist(accepted_df['I'], bins=range(INCUBATION_LOW, INCUBATION_HIGH + 2))
        plt.title(f"Incubation period (accepted) for case {case}")

        plt.subplot(2, 2, 4)
        plt.hist(accepted_df['IC'], bins=range(SEED_LOW, SEED_HIGH + 2))
        plt.title(f"Initial infectious cows (accepted) for case {case}")

        plt.tight_layout()
        plt.show()
        plt.close()

    # Posterior simulations and fan plots
    for case in cases:
        print(f"\nRunning posterior simulations for case {case}...")
        accepted_df = accepted_dfs[case]
        inputs = case_inputs[case]
        ndays_target_case = inputs['ndays_target']
        burn_in_farm_case = inputs['burn_in_farm']
        target_curve_case = inputs['target_curve']
        no_calves_case = inputs['no_calves']
        burn_in_days_case = inputs['burn_in_days']
        # Sample with replacement from accepted posterior dataframe
        posterior_df = accepted_df.sample(n=N_post_sims, replace=True).reset_index(drop=True)

        all_post_curves = []
        for _, row in posterior_df.iterrows():
            T = float(row['T'])
            R = float(row['R'])
            incubation = int(row['I'])
            seed = int(row['IC'])
            sim_curve = simulate_counts(burn_in_farm_case, ndays_target_case, burn_in_days_case,
                                        seed, T, R, incubation, no_calves_case)
            all_post_curves.append(sim_curve)
        all_post_curves = np.array(all_post_curves)

        # Compute posterior predictive median and 95% envelope
        median_curve = np.median(all_post_curves, axis=0)
        low_curve = np.percentile(all_post_curves, 2.5, axis=0)
        high_curve = np.percentile(all_post_curves, 97.5, axis=0)

        # Plot observed vs posterior predictive distribution (fan plot)
        plt.figure(figsize=(10, 6))
        days = np.arange(1, ndays_target_case + 1)

        # Fan plot: all individual posterior simulations
        for curve in all_post_curves:
            plt.plot(days, curve, color='gray', alpha=0.05)

        # Uncertainty band
        plt.fill_between(days, low_curve, high_curve, color='gray', alpha=0.3, label='95% posterior credible interval')

        # Posterior median
        plt.plot(days, median_curve, color='blue', lw=2, label='Posterior median simulation')

        # Observed outbreak
        plt.plot(days, target_curve_case, 'o-', color='red', lw=2, label='Observed outbreak')

        plt.xlabel("Day")
        plt.ylabel("New Infectious")
        plt.title(f"Case {case}: Observed vs Posterior Predictive Simulations (N={N_post_sims})")
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.close()

    return accepted_dfs

def analyze_optimization_results(cases, case_inputs, bayesian_optimization_result, pens_to_mask = ['Pen1', 'Pen2', 'Pen3', 'Pen4']):
    for case in cases:
        print(f"Bayesian Optimization Results Summary for {case}")

    # Plotting Loss vs Parameters from Bayesian Optimization
        if bayesian_optimization_result[case] is not None and hasattr(bayesian_optimization_result[case], 'x_iters') and hasattr(bayesian_optimization_result[case], 'func_vals'):
            print("\nPlotting Loss metric vs Optimization Parameters...")

            # Extract parameter values and corresponding loss values
            transmission_rates_evaluated = [x[0] for x in bayesian_optimization_result[case].x_iters]
            recovery_rates_evaluated = [x[1] for x in bayesian_optimization_result[case].x_iters]
            incubation_periods_evaluated = [int(round(x[2])) for x in bayesian_optimization_result[case].x_iters]
            initial_infected_evaluated = [int(round(x[3])) for x in bayesian_optimization_result[case].x_iters]
            loss_values_evaluated = bayesian_optimization_result[case].func_vals

            # Create scatter plot for Loss vs Transmission Rate
            plt.figure(figsize=(10, 6))
            plt.scatter(transmission_rates_evaluated, loss_values_evaluated, c=loss_values_evaluated, cmap='viridis')
            plt.xlabel('Transmission Rate')
            plt.ylabel('Loss metric')
            plt.title('Loss metric vs Transmission Rate during Bayesian Optimization')
            plt.colorbar(label='Loss metric')
            plt.grid(True)
            plt.show()

            # Create scatter plot for Loss vs Recovery Rate
            plt.figure(figsize=(10, 6))
            plt.scatter(recovery_rates_evaluated, loss_values_evaluated, c=loss_values_evaluated, cmap='viridis')
            plt.xlabel('Recovery Rate')
            plt.ylabel('Loss metric')
            plt.title('Loss metric vs Recovery Rate during Bayesian Optimization')
            plt.colorbar(label='Loss metric')
            plt.grid(True)
            plt.show()

            # Create scatter plot for Loss vs Incubation Period
            plt.figure(figsize=(10, 6))
            plt.scatter(incubation_periods_evaluated, loss_values_evaluated, c=loss_values_evaluated, cmap='viridis')
            plt.xlabel('Incubation Period (days)')
            plt.ylabel('Loss metric')
            plt.title('Loss metric vs Incubation Period during Bayesian Optimization')
            plt.colorbar(label='Loss metric')
            plt.grid(True)
            plt.show()
                        
            # Create scatter plot for Loss vs Initial Infected
            plt.figure(figsize=(10, 6))
            plt.scatter(initial_infected_evaluated, loss_values_evaluated, c=loss_values_evaluated, cmap='viridis')
            plt.xlabel('Incubation Period (days)')
            plt.ylabel('Loss metric')
            plt.title('Loss metric vs Incubation Period during Bayesian Optimization')
            plt.colorbar(label='Loss metric')
            plt.grid(True)
            plt.show()
            
            inputs = case_inputs[case]
            sim = run_main_simulation(burn_in_farm_state = inputs['burn_in_farm'], 
                              merge_dwell_df = merge_dwell, 
                              ndays = inputs['ndays_target'], 
                              burn_in_days = inputs['burn_in_days'], 
                              initial_infected = int(round(bayesian_optimization_result[case].x[3])),
                              transmission_rate = bayesian_optimization_result[case].x[0], 
                              recovery_rate = bayesian_optimization_result[case].x[1], 
                              incubation_period=int(round(bayesian_optimization_result[case].x[2])), 
                              actual_infectious_counts=None, 
                              pens_to_mask=inputs['no_calves'], 
                              num_simulations=500, plot=True)

        else:
            print("\nCannot plot Loss vs Optimization Parameters: Bayesian optimization results are not available or incomplete.")

