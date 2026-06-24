import os
import pandas as pd
import manifest

# run test
test_run = False

# parameters
expname = 'PfPR_calibration_burnin_sweep_example'
population_size = 6000  # if this population size is different that what was used for seasonality calibration, need to run '00_rescale_demog_vector_files.py'
num_seeds = 3
years = 30
serialize = True
pull_from_serialization = False
itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]

demographics_file = os.path.join('demographics_and_climate', '_entire_country',
                                 f'demographics_each_admin_{population_size}.json')

# set file describing which intervention inputs are used for this scenario
scenario_fname = os.path.join(manifest.project_path, 'simulation_inputs', '_intervention_file_references',
                              'Interventions_for_calibration.csv')
scen_df = pd.read_csv(scenario_fname)
# row matching the burnin
scen_index = scen_df[scen_df['ScenarioName'] == 'transmission_calibration_burnin'].index[0]

