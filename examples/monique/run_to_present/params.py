import os
import manifest
import pandas as pd

num_seeds = 3
start_year = 2010
years = 12  # beginning of 2010 to beginning of 2010+ year
serialize = True
pull_from_serialization = True
num_burnin_seeds = 5  # number of seeds run during previous burnin and habitat multiplier calibration
ser_date = 30 * 365
itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
population_size = 6000  # needs to match burnin simulation population size

demographics_file = os.path.join('demographics_and_climate', '_entire_country',
                                 f'demographics_each_admin_{population_size}.json')

use_arch_burnin = True
burnin_id = 'c2a74e6d-a199-f011-aa26-b88303911bc1' # <- Calculon  #     'eb89ea6e-2f53-ee11-aa0a-b88303911bc1'  # 'cab2ccbc-f70a-ee11-aa07-b88303911bc1'  # <- on IDMCloud  # generated from serialize_transmission_sweep (1960-2010)

scenario_fname = os.path.join(manifest.project_path, 'simulation_inputs', '_intervention_file_references',
                              'Interventions_to_present.csv')
scen_df = pd.read_csv(scenario_fname)
scen_index = scen_df[scen_df['status'] == 'queued'].index[0]  # submitter changed this to queued (means queued for submitting to comps)
expname = scen_df.at[scen_index, 'ScenarioName']
experiment_type="to_present"
