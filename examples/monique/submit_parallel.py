import os
import sys
import subprocess
import pandas as pd
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform
from idmtools.entities import Suite
from snt.load_paths import load_box_paths
from datetime import datetime

# Configuration
USER_PATH = None
CURRENT_DIRECTORY = os.path.dirname(__file__)
FUTURE_PROJECTIONS = False
ANALYZER_STATUS = None

experiment_type = "future_projections" if FUTURE_PROJECTIONS else "to_present"
SUITE_TRACK_FILE = os.path.join(CURRENT_DIRECTORY, f'suite_tracking_{experiment_type}.csv')
# Format: YYYYMMDD_HHMMSS
time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
background_processes = []

def ensure_tracking_file():
    """
    Create the tracking file if it doesn't exist or is empty.
    """
    # Create a tracking file with the title row if it doesn't exist or is empty.
    if not os.path.exists(SUITE_TRACK_FILE) or os.path.getsize(SUITE_TRACK_FILE) == 0:
        with open(SUITE_TRACK_FILE, 'w') as f:
            f.write("suite_name,suite_id,experiment_name,experiment_id,experiment_type,start_time_stamp,end_time_stamp\n")


def create_suite(platform, experiment_type):
    """
    Create a new suite in COMPS.
    Args
        platform: COMPS platform object.
        experiment_type: experiment type, either 'future_projections' or 'to_present'
    Return:
        suite object.
    """
    # Create a new suite in COMPS.
    suite = Suite(name='SNT Suite')
    suite.update_tags({'name': f'suite_{experiment_type}'})
    platform.create_items([suite])
    print(f"Created new suite: {suite.name} ({suite.id})")
    return suite

def get_scenario_file_and_script(project_path, experiment_type):
    if experiment_type=='future_projections':
        scenario_file = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_for_projections.csv')
        sim_script = './run_future_scenarios/run_simulations.py'
    elif experiment_type=='to_present':
        scenario_file = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_to_present.csv')
        sim_script = './run_to_present/run_simulations.py'
    else:
        raise Exception("no scenario file.")
    return scenario_file, sim_script


def submit_experiment_nonblocking(sim_script, suite_id, idx, scenario_file, log_file_path):
    """
    Submit an experiment to COMPS in non-blocking mode
    Args:
        sim_script: str, path to the simulation script.
        suite_id: str, ID of the suite to which the experiment belongs.
        idx: int, index of the experiment in the scenario file.
        scenario_file: str, path to the scenario file.
        log_file_path: str, path to the log file.
    """
    os.environ["NO_COLOR"] = "1"  # disable color in subprocess's log
    with open(log_file_path, "a", encoding="utf-8") as f:
        proc = subprocess.Popen(
            [
                sys.executable,
                sim_script,
                "--suite-id", str(suite_id),
                "--scen-index", str(idx),
                "--scenario-fname", str(scenario_file)
            ],
            stdout=f,
            stderr=subprocess.STDOUT,
            start_new_session=True,  # fully detached background process
            encoding='utf-8',
            errors='replace',
            text=True
        )
    print(f"Launched scenario {idx} → PID {proc.pid}, log: {log_file_path}")

def update_tracking_file(suite_name, suite_id, experiment_type):
    """
    Write the tracking file with new suite_id.
    Args:
        suite_name: suite name
        suite_id: suite id
        experiment_type: experiment type, either 'future_projections' or 'to_present'
    """
    df = pd.read_csv(SUITE_TRACK_FILE)
    if not (df['suite_id'] == suite_id).any():  # Only add new suite_id to tracking file
        with open(SUITE_TRACK_FILE, 'a') as f:
            f.write(f"{suite_name},{suite_id},None,None,{experiment_type},{time_stamp},None\n")


def main():
    # Create suite tracking file if it doesn't exist.'
    ensure_tracking_file()
    # Load Box path
    _, project_path = load_box_paths(user_path=USER_PATH, country_name='Example')
    # Initialize platform.
    platform = Platform("CALCULON")
    # Create a new suite in COMPS.
    suite = create_suite(platform, experiment_type)
    # To use existing suite,comment out above line and uncomment following line
    #suite = platform.get_item("your suite_id", item_type=ItemType.SUITE)
    suite_id, suite_name = str(suite.id), suite.name
    print(f"\nThe created suite can be viewed at {platform.endpoint}/#explore/"
                             f"Suites?filters=Id={suite_id}\n")
    # Append a new suite_id to a tracking file.
    update_tracking_file(suite_name, suite_id, experiment_type)

    # Select a simulation run script based on an experiment type.
    scenario_file, sim_script = get_scenario_file_and_script(project_path, experiment_type)

    # Read a scenario file and filter for rows with status 'run' to submit simulations to COMPS.
    df = pd.read_csv(scenario_file)
    df_run = df[df['status'] == 'run'].copy()

    # Create a log directory if it doesn't exist.'
    log_dir = os.path.join(CURRENT_DIRECTORY, "logs")
    os.makedirs(log_dir, exist_ok=True)
    print(f"Log directory created at {log_dir}\n")

    # Submit simulations to COMPS in non-blocking mode.
    # Iterate over rows with the status 'run' in the scenario file. And submit them to COMPS.
    for idx in df_run.index:
        print(f"Before submit: Experiment at row {idx} status: '{df.loc[idx, 'status']}' in CSV file.")
        # log file path:
        log_path = os.path.join(log_dir, f"scenario_{experiment_type}_{idx}.log")
        # submit the experiment to COMPS in non-blocking mode with background processes.
        # to see progress, open the log file
        submit_experiment_nonblocking(sim_script, suite_id, idx, scenario_file, log_path)
        df.loc[idx, 'status'] = 'queued'
        print(f"After submit: Experiment at row {idx} status: '{df.loc[idx, 'status']}' in CSV file.")
        print(f"Dispatched experiment at row {idx} to COMPS\n")

    df.to_csv(scenario_file, index=False)
    print("\nAll experiments dispatched in non-blocking mode.")


if __name__ == "__main__":
    main()
