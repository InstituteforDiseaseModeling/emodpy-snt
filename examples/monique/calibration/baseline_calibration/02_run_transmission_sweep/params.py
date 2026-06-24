import os
import pandas as pd
import manifest

burnin_id = 'eb89ea6e-2f53-ee11-aa0a-b88303911bc1'  # generated from serialize_transmission_sweep
expname = 'PfPR_sweep_main_example'
population_size = 6000  # needs to match burnin simulation population size
num_seeds = 1
num_burnin_seeds = 1  # number of seeds run during previous burnin
ser_date = 10 * 365
years = 8  # beginning of 2010 - end of 2017
start_year = 2010
itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
pull_from_serialization = True

demographics_file = os.path.join('demographics_and_climate', '_entire_country',
                                 f'demographics_each_admin_{population_size}.json')

# set file describing which intervention inputs are used for this scenario
scenario_fname = os.path.join(manifest.project_path, 'simulation_inputs', '_intervention_file_references',
                              'Interventions_for_calibration.csv')
scen_df = pd.read_csv(scenario_fname)
# row matching this main calibration (not the burnin)
scen_index = scen_df[scen_df['ScenarioName'] == 'transmission_calibration_main'].index[0]
