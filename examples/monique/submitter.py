import os
import pandas as pd
from snt.load_paths import load_box_paths

# Specify data location
USER_PATH = None
# USER_PATH = r'C:\Projects\emodpy-snt'
_, project_path = load_box_paths(user_path=USER_PATH, country_name='Example')

wdir = os.path.join(project_path, 'simulation_inputs')
future_projections = True

while True:
    if future_projections:
        scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references',
                                      'Interventions_for_projections.csv')
    else:
        scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references',
                                      'Interventions_to_present.csv')
    df = pd.read_csv(scenario_fname)

    if len(df[df['status'] == 'run']) == 0:
        break

    if future_projections:
        os.system('python ./run_future_scenarios/run_simulations.py')
    else:
        os.system('python ./run_to_present/run_simulations.py')

    i = df[df['status'] == 'run'].index[0]
    df.loc[i, 'status'] = 'queued'
    df.to_csv(scenario_fname, index=False)
