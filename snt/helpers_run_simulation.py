import os
import subprocess
import sys
from importlib import import_module

import pandas as pd
from filelock import FileLock
from datetime import datetime

from idmtools.core import ItemType
from idmtools.entities.experiment import Experiment
from idmtools.entities.templated_simulation import TemplatedSimulations

def log(msg):
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S,%f')[:-3]} - {msg}")

def _pre_run(experiment, **kwargs):
    from snt.utility.plugins import initialize_plugins
    initialize_plugins(**kwargs)


def _post_run(params, manifest, experiment, suite_id, tracking_file, **kwargs):
    """
    Add extra work after run experiment.
    In this case, we launch the SSMT analyzer after the experiment is done with subprocess in background.
    Logs are written to logs/ssmt_analyzer_{experiment_type}.log.
    Args:
        params:
        manifest:
        experiment: idmtools Experiment
        suite_id: idmtools Suite ID
        tracking_file: suite/experiment tracking file path, e.g. suite_tracking_future_projections.csv or suite_tracking_to_present.csv
        kwargs: additional parameters including platform
    Return:
        None
    """
    platform = kwargs['platform']
    if experiment.succeeded:
        suite = platform.get_item(suite_id, ItemType.SUITE, raw=True)
        time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        update_suite_tracking(experiment, suite, params.experiment_type, tracking_file, time_stamp)


        analyzer_script = os.path.join(manifest.CURRENT_DIR, "..", "analyzers", "post_ssmt.py")
        # replace with local analyzer script if needed
        #analyzer_script = os.path.join(manifest.CURRENT_DIR, "..", "analyzers", "post_analysis.py")

        ssmt_log = os.path.abspath(os.path.join(manifest.TRACKING_DIR, "..", "logs", f"ssmt_analyzer_{experiment.name}_{experiment.id}.log"))

        os.environ["NO_COLOR"] = "1"  # disable color in subprocess's log
        # run ssmt (or local analyzer) subprocess in background
        with open(ssmt_log, "a+", encoding="utf-8") as log_file:
            subprocess.Popen([
                sys.executable,
                analyzer_script,
                "--exp-id", str(experiment.id),
                "--type", str(params.experiment_type),  # or hardcode 'future_projections' etc.
                "--name", str(experiment.name)
            ],
                stdout=log_file,
                stderr=subprocess.STDOUT,
                encoding='utf-8'
            )
        log(f"Launched SSMT analyzer for experiment: {experiment.name}--{experiment.id}")
        log(f"Analyzer_log: {ssmt_log}")


def post_run(params, manifest, experiment, suite_id, tracking_file, **kwargs):
    """
    post_run is called after the experiment is done. It writes experiment name and id to a file for the analyzer manager to pick up.
    Args:
        params:
        manifest:
        experiment: idmtools Experiment
        suite_id: idmtools Suite ID
        tracking_file: suite/experiment tracking file path, e.g. suite_tracking_future_projections.csv or suite_tracking_to_present.csv
        kwargs: additional parameters including platform

    Returns:

    """
    platform = kwargs['platform']
    if experiment.succeeded:
        suite = platform.get_item(suite_id, ItemType.SUITE, raw=True)
        time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        update_suite_tracking(experiment, suite, params.experiment_type, tracking_file, time_stamp)
        # Signal to analyzer manager
        queue_dir = os.path.join(manifest.TRACKING_DIR, "analyzer_queue")
        os.makedirs(queue_dir, exist_ok=True)
        flag_file = os.path.abspath(os.path.join(queue_dir, f"exp_{experiment.id}.ready"))

        with open(flag_file, "w") as f:
            f.write(f"{experiment.id},{experiment.name},{params.experiment_type}\n")

        log(f"Enqueued experiment for analyzer: {flag_file}")

        # uncomment out the following if you want to update suite tracking file with experiment status
        # suite = platform.get_item(suite_id, ItemType.SUITE, raw=True)
        # time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        # update_suite_tracking(experiment, suite, params.experiment_type, tracking_file, time_stamp)

def _config_experiment(params, manifest, **kwargs):
    """
    Build experiment from task and builder. task is EMODTask. builder is SimulationBuilder used for config parameter sweeping.
    Args:
        kwargs: additional parameters
    Return:
        experiment
    """
    # Dynamically import config_task and config_sweep_builders from current working directory
    get_task = import_module("config_task").get_task
    get_sweep_builders = import_module("config_sweep_builders").get_sweep_builders

    builders = get_sweep_builders(**kwargs)
    task = get_task(**kwargs)

    if manifest.sif_id:
        task.set_sif(manifest.sif_id)

    ts = TemplatedSimulations(base_task=task, builders=builders)

    df = pd.read_csv(params.scenario_fname)
    scen_name = df.at[params.scen_index, 'ScenarioName']
    experiment_name = f"{params.experiment_type}_{scen_name}"

    experiment = Experiment.from_template(ts, name=experiment_name)
    return experiment

def run_experiment(platform, params, manifest, print_params_fn, suite_id, tracking_file, **kwargs):
    """
    Get configured experiment and run.
    Args:
        kwargs: user inputs
    Returns:
        None
    """
    kwargs['platform'] = platform  # This is need to _config_experiment's get_sweep_builders(**kwargs)
    print_params_fn()  # call the injected function

    experiment = _config_experiment(params, manifest, **kwargs)  # use the injected builder
    _pre_run(experiment, **kwargs)

    if suite_id:
        experiment.parent_id = suite_id
    try:
        experiment.run(wait_until_done=True, platform=platform)
        # experiment = platform.get_item("cc510fbd-7b67-f011-9f17-b88303912b51", ItemType.EXPERIMENT)
        update_scenario_status_to_done(params.scenario_fname, params.scen_index)

        # Note:
        # Choose one of the following two lines (and comment out no needed one), depending on your preference or requirement
        # _post_run: Analyzer is launched in background automatically after experiment is done.
        # post_run: For signaling to analyzer queue manager (another script runs analyzer in the background).
       # _post_run(params, manifest, experiment, suite_id, tracking_file, **kwargs)
        post_run(params, manifest, experiment, suite_id, tracking_file, **kwargs)

    except Exception as e:
        log(f"Experiment run failed: {e}")
        with open(os.path.join(manifest.TRACKING_DIR, "failed_experiments.log"), "a") as f:
            f.write(f"{params.scenario_fname}, index {params.scen_index} failed for exp {experiment}: {e}\n")


def update_suite_tracking(experiment, suite, experiment_type, tracking_file, time_stamp):
    """
    Update suite tracking file.
    Args:
        experiment:
        suite:
        experiment_type:
        tracking_file:
        time_stamp:

    Returns:

    """
    lock_path = tracking_file + ".lock"
    with FileLock(lock_path, timeout=30):
        # Ensure the file exists with headers
        if not os.path.exists(tracking_file):
            with open(tracking_file, "w") as f:
                f.write("suite_name,suite_id,experiment_name,experiment_id,experiment_type,start_time_stamp,end_time_stamp\n")

        df = pd.read_csv(tracking_file)

        match = df['suite_id'] == str(suite.id)

        if match.any():
            # Update existing row(s)
            df.loc[match, 'experiment_name'] = experiment.name
            df.loc[match, 'experiment_id'] = experiment.id
            df.loc[match, 'end_time_stamp'] = time_stamp
        else:
            # Add new row if suite_id not found
            new_row = {
                'suite_name': suite.name,
                'suite_id': suite.id,
                'experiment_name': experiment.name,
                'experiment_id': experiment.id,
                'experiment_type': experiment_type,
                'start_time_stamp': time_stamp,
                'end_time_stamp': None
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

        # Save changes
        df.to_csv(tracking_file, index=False)


def update_scenario_status_to_done(scenario_fname, scen_index):
    """
    Update scenario status to 'done' in scenario_fname. This is set to 'done' when the experiment is done.
    Args:
        scenario_fname: scenario file path
        scen_index: scenario index
    """
    lockfile = scenario_fname + ".lock"
    with FileLock(lockfile, timeout=30):
        df = pd.read_csv(scenario_fname)
        df.loc[scen_index, 'status'] = 'done'
        df.to_csv(scenario_fname, index=False)
        log(f"Updated status in scenario index {scen_index} to 'done' in {scenario_fname}")
