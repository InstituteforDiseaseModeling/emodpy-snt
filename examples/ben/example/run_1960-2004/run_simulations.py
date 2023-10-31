import manifest
import params
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from idmtools.entities.templated_simulation import TemplatedSimulations


def _print_params():
    """
    Just a useful convenient function for the user.
    """
    print("test_run: ", params.test_run)
    print("iopath: ", params.iopath)
    print("expname: ", params.expname)
    print("homepath: ", params.homepath)
    print("user: ", params.user)
    print("serialize: ", params.serialize)
    print("num_seeds: ", params.num_seeds)
    print("years: ", params.years)
    print("pull_from_serialization: ", params.pull_from_serialization)


def _pre_run(experiment: Experiment, **kwargs):
    """
    Add extra work before run experiment.
    Args:
        experiment: idmtools Experiment
        kwargs: additional parameters
    Return:
        None
    """
    from snt.utility.plugins import initialize_plugins
    initialize_plugins(**kwargs)


def _post_run(experiment: Experiment, **kwargs):
    """
    Add extra work after run experiment.
    Args:
        experiment: idmtools Experiment
        kwargs: additional parameters
    Return:
        None
    """
    wait_until_done = kwargs.get('wait_until_done', None)
    if wait_until_done and experiment.succeeded:
        with open("ben/example/run_1960-2004/experiment_id.txt", "w") as fd:
            fd.write(experiment.uid.hex)


def _config_experiment(**kwargs):
    """
    Build experiment from task and builder. task is EMODTask. builder is SimulationBuilder used for config parameter sweeping.
    Args:
        kwargs: additional parameters
    Return:
        experiment
    """
    from config_sweep_builders import get_sweep_builders
    from config_task import get_task

    builders = get_sweep_builders(**kwargs)

    task = get_task(**kwargs)

    if manifest.sif_path:
        task.set_sif(manifest.sif_path)

    ts = TemplatedSimulations(base_task=task, builders=builders)

    experiment = Experiment.from_template(ts, name=params.expname)

    return experiment


def run_experiment(**kwargs):
    """
    Get configured experiment and run.
    Args:
        kwargs: user inputs
    Returns:
        None
    """
    # make sure pass platform through
    kwargs['platform'] = platform

    _print_params()

    experiment = _config_experiment(**kwargs)
    _pre_run(experiment, **kwargs)
    experiment.run(**kwargs)
    _post_run(experiment, **kwargs)


if __name__ == "__main__":
    """
    - show_warnings_once=True:  show api warnings for only one simulation
    - show_warnings_once=False: show api warnings for all simulations
    - show_warnings_once=None:  not show api warnings
    """
    platform = Platform('CALCULON', node_group='idm_48cores')
    # platform = Platform('IDMCLOUD', node_group='emod_abcd')

    # To use Slurm Platform: Specify job directory
    # job_directory = r'C:\Projects\emodpy-snt\data\TEST_DEST'
    # platform = Platform('SLURM_LOCAL', job_directory=job_directory)

    # If you don't have Eradication, un-comment out the following to download Eradication
    # import emod_malaria.bootstrap as dtk
    # import pathlib
    # import os
    # dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    # os.chdir(os.path.dirname(__file__))
    # print("...done.")
    run_experiment(show_warnings_once=True, wait_until_done=True)
