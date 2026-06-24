import os
import manifest
import pandas as pd
import numpy as np
from snt.calibration.helpers_seasonality_calibration import get_spline_values4_constantMaxHab, \
    get_spline_values4, update_starting_spline_values_constantMaxHab, update_starting_spline_values

project_path = manifest.project_path

itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
simulation_pop = 10000  # from demographics file
vector_human_scalar = simulation_pop / 1000
max_incidence = 100

# specify name of archetype's representative admin
rep_admin = 'AA'
# burn-in id should be listed for each of the archetypes
burnin_ids = {
    'AA': '7fbc7c53-c44d-ee11-aa0a-b88303911bc1'
}

# set whether this is the first or second set of calibrations for this archetype
round_number = 1  # <-- user needs to specify which round this is to pick up from prior round. Starts at round_number=1
if round_number > 1:
    later_round = True  # if later_round is True, we use the previous best estimate as the starting point
    if round_number >= 2:
        constant_max_LH = True  # if constant_max_LH is True, we fix maxHab at the best estimate and only adjust other params
    else:
        constant_max_LH = False
else:
    later_round = False
    constant_max_LH = False

pull_from_serialization = True  # should be set to True unless run burn-in within each calibration simulation (warning: that takes a long time)
if pull_from_serialization:
    throwaway = 3  # duration of simulation in years, -1.
    burnin_years = 20  # year when serialization occurred in burnin
    burnin_max_habitat_value = 10 + np.log10(
        vector_human_scalar)  # make sure this value matches that used in the burnin
    burnin_id = burnin_ids[rep_admin]
else:
    throwaway = 29  # duration of simulation in years, -1.

# set file describing which intervention inputs are used for this scenario
scenario_fname = os.path.join(manifest.project_path, 'simulation_inputs', '_intervention_file_references',
                              'Interventions_for_calibration.csv')
scen_df = pd.read_csv(scenario_fname)
# row matching this main calibration (not the burnin)
scen_index = scen_df[scen_df['ScenarioName'] == 'seasonality_calibration_noMassITN'].index[0]
expname = '%s_%s_round%i_maxInc%i_v2' % (scen_df.at[scen_index, 'ScenarioName'], rep_admin, round_number, max_incidence)

demographics_file = os.path.join('demographics_and_climate', '_entire_country',
                                 f'demographics_each_admin_{simulation_pop}.json')

# set calibration algorithm parameters (some are different between earlier and later round since we make smaller jumps to match target)
r = 0.3  # radius - determines how close the parameter values of each sample within an iteration is to the center value
r_second = 0.005  # 0.1
r_fourth = 0.001
center_move_scale = 1 / 15
center_move_scale_second = 1 / 25  # 1/20
center_move_scale_fourth = 1 / 30
sim_runs_per_param_set = 2  # outside of testing, generally 5
sim_runs_per_param_set_second = 5  # 10
sim_runs_per_param_set_fourth = 10  # 30
max_iterations = 5  # outside of testing, generally 10
max_iterations_second = 5   # 20
max_iterations_fourth = 5   # 20
samples_per_iteration = 15  # outside of testing, generally 40  # must be at least 8

#############################################################################################################
# set up the parameters that will be used to calibrate vector larval habitat seasonality

# helper files to specify relative abundance of vector species + initialization of monthly habitats
if constant_max_LH:
    spline, fractions = get_spline_values4_constantMaxHab(rep_admin,
                                                          project_path)  # for starting guess on monthly multipliers, use rainfall proportional to maximum rainfall month
else:
    spline, fractions = get_spline_values4(rep_admin,
                                           project_path)  # for starting guess on monthly multipliers, use rainfall proportional to maximum rainfall month

# if the second round, read in the best values from the previous round and use those as the starting point
if later_round:
    round1_fname = os.path.join('%s_%s_round%i_maxInc%i_v2' % (
    scen_df.at[scen_index, 'ScenarioName'], rep_admin, (round_number - 1), max_incidence), '_plots', 'LL_all.csv')

    round1_df = pd.read_csv(round1_fname)
    # use row with best match to reference dataset
    round1_best_df = round1_df[round1_df.total == round1_df.total.max()].reset_index()
    if constant_max_LH:
        spline = update_starting_spline_values_constantMaxHab(spline, round1_best_df)
    else:
        spline = update_starting_spline_values(spline, round1_best_df)
    if round_number < 4:
        r = r_second
        center_move_scale = center_move_scale_second
        sim_runs_per_param_set = sim_runs_per_param_set_second
        max_iterations = max_iterations_second
    else:
        r = r_fourth
        center_move_scale = center_move_scale_fourth
        sim_runs_per_param_set = sim_runs_per_param_set_fourth
        max_iterations = max_iterations_fourth

params = []
# add months to list of params being sampled
for i, row in spline.iterrows():
    d = row.to_dict()
    params.append(d)