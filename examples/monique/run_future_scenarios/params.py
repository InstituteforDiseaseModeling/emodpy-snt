import os
import manifest
import pandas as pd

num_seeds = 3
start_year = 2022
years = 8  # beginning of start_year to beginning of (start_year + year)
serialize = False
pull_from_serialization = True
num_burnin_seeds_calib = 1  # number of seeds run during transmission-intensity calibration simulations (to get xLHs)
num_burnin_seeds = 3  # number of seeds run during "to-present" simulations
ser_date = 12 * 365
population_size = 6000  # needs to match burnin simulation population size
itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
LLIN_decay_2_years = False

burnin_id = '352e2f2f-0c61-f011-aa23-b88303911bc1'  # generated from 2010-2021 run

scenario_fname = os.path.join(manifest.project_path, 'simulation_inputs', '_intervention_file_references',
                              'Interventions_for_projections.csv')
scen_df = pd.read_csv(scenario_fname)
scen_index = scen_df[scen_df['status'] == 'queued'].index[0]
expname = scen_df.at[scen_index, 'ScenarioName']
expname = f'{expname}_v3'

demographics_file = os.path.join('demographics_and_climate', '_entire_country',
                                 f'demographics_each_admin_{population_size}.json')
experiment_type="future_projections"