import copy
import numpy as np
from idmtools.core.platform_factory import Platform
from idmtools_calibra.calib_manager import CalibManager
from idmtools_calibra.plotters.likelihood_plotter import LikelihoodPlotter
from idmtools_calibra.plotters.site_data_plotter import SiteDataPlotter
import emod_api.config.default_from_schema_no_validation as dfs
from emodpy_malaria.malaria_config import configure_linear_spline, set_species_param
from snt.algorithms.optim_tool import OptimTool
from snt.calibration.SeasonalityCalibSite import SeasonalityCalibSite

import manifest
import params


def print_params():
    """
    Just a useful convenience function for the user.
    """
    # Display exp_name and nSims
    print("expname: ", params.expname)
    print("max_iterations: ", params.max_iterations)
    print("samples_per_iteration: ", params.samples_per_iteration)
    print("pull_from_serialization: ", params.pull_from_serialization)
    print("later_round: ", params.later_round)
    print("round_number: ", params.round_number)
    print("burnin_ids: ", params.burnin_ids)


def _pre_run(**kwargs):
    """
    Add extra work before run calibration.
    Args:
        kwargs: additional parameters
    Return:
        None
    """
    from snt.utility.plugins import initialize_plugins
    initialize_plugins(**kwargs)


def constrain_sample(sample):
    """
    This function is called on every sample and allow the user to edit them before they are passed
    to the map_sample_to_model_input function.
    It is useful to round some parameters as demonstrated below.
    Can do much more here, e.g. for
    # Clinical Fever Threshold High <  MSP1 Merozoite Kill Fraction
    if 'Clinical Fever Threshold High' and 'MSP1 Merozoite Kill Fraction' in sample:
        sample['Clinical Fever Threshold High'] = \
            min( sample['Clinical Fever Threshold High'], sample['MSP1 Merozoite Kill Fraction'] )
    You can omit this function by not specifying it in the OptimTool constructor.

    Args:
        sample: represents a sample used to generate simulations

    Returns:
        sample
    """
    return sample


def map_sample_to_model_input(simulation, sample):
    # set up some variables to hold the values we'll need to update the monthly habitat (hab)
    tags = {}
    sample = copy.deepcopy(sample)

    month_vals = []
    # pull out the variable values at the current step of calibration
    for ii in range(1, 13):
        val_name = 'MonthVal%d' % ii
        cur_value = sample.pop(val_name)
        month_vals.append(cur_value)
        tags.update({val_name: cur_value})
    my_spline = month_vals

    max_habitat_name = 'MaxHab'
    maxvalue = sample.pop(max_habitat_name)
    tags.update({max_habitat_name: maxvalue})
    # constant habitat is currently set as a fixed parameter, not fitted
    const_name = 'Constant'
    const = 3 + np.log10(params.vector_human_scalar)
    tags.update({const_name: const})

    habitat = dfs.schema_to_config_subnode(manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"])
    for (s, sp) in zip(params.fractions, ['arabiensis', 'funestus', 'gambiae']):
        linear_spline_habitat = configure_linear_spline(manifest,
                                                        max_larval_capacity=pow(10, maxvalue) * s,
                                                        capacity_distribution_number_of_years=1,
                                                        capacity_distribution_over_time={
                                                            "Times": [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
                                                                      334],
                                                            "Values": list(my_spline)
                                                        }
                                                        )
        set_species_param(simulation.task.config, sp, "Habitats", linear_spline_habitat, overwrite=True)

        new_habitat = copy.deepcopy(habitat)
        new_habitat.parameters.Habitat_Type = "CONSTANT"
        new_habitat.parameters.Max_Larval_Capacity = pow(10, const) * s
        set_species_param(simulation.task.config, sp, "Habitats", new_habitat.parameters)

    # tags.update({'Pop_Scale' : 1})

    return tags


def _config_manager(task):
    """
    This function is designed to create and config calibration.
    Args:
        task: represents EMODTask
    Returns:
        Calib CalibManager object
    """
    sites = [SeasonalityCalibSite(hfca=params.rep_admin, throwaway=params.throwaway, project_path=params.project_path)]
    plotters = [LikelihoodPlotter(combine_sites=True), SiteDataPlotter(num_to_plot=0, combine_sites=True)]

    optimtool = OptimTool(params.params,
                          mu_r=params.r,
                          # <-- radius for numerical derivatve.  CAREFUL not to go too small with integer parameters
                          sigma_r=params.r / 10.,  # <-- stdev of radius
                          center_repeats=1,
                          # <-- Number of times to replicate the center (current guess).  Nice to compare intrinsic to extrinsic noise
                          samples_per_iteration=params.samples_per_iteration,
                          # <-- Samples per iteration, includes center repeats.  Actual number of sims run is this number times number of sites.
                          center_move_scale=params.center_move_scale  # controls
                          )

    calib_manager = CalibManager(name=params.expname,
                                 task=task,
                                 map_sample_to_model_input_fn=map_sample_to_model_input,
                                 sites=sites,
                                 next_point=optimtool,
                                 sim_runs_per_param_set=params.sim_runs_per_param_set,  # <-- Replicates
                                 max_iterations=params.max_iterations,  # <-- Iterations
                                 plotters=plotters)

    return calib_manager


def get_manager(**kwargs):
    """
    Config a task and then config a CalibManager.
    Args:
        kwargs: user inputs
    Returns:
        CalibManager
    """
    from config_task import _config_task
    task = _config_task(**kwargs)

    if manifest.sif_id:
        task.set_sif(manifest.sif_id)

    calib_manager = _config_manager(task)

    return calib_manager


def _post_run(suite_id, **kwargs):
    """
    Add extra work after run calibration.
    Args:
        suite_id:
        kwargs: additional parameters
    Return:
        None
    """
    from idmtools.core import ItemType
    suite = platform.get_item(suite_id, ItemType.SUITE)

    if all([exp.succeeded for exp in suite.experiments]):
        with open("monique/calibration/seasonality_calibration/02_seasonality_calibration/experiment_id.txt", "w") as fd:
            fd.write("Calibration good.")


def run_calibration(directory: str = '.', **kwargs):
    """
    Get configured calibration and run.
    Args:
        directory: str, where to keep calibration results
        kwargs: user inputs
    Returns:
        None
    """
    # make sure pass platform through
    kwargs['platform'] = platform
    kwargs['directory'] = directory

    print_params()

    calib_manager = get_manager(**kwargs)
    _pre_run(**kwargs)
    calib_manager.run_calibration(**kwargs)
    _post_run(calib_manager.suite_id, **kwargs)


if __name__ == "__main__":
    """
    - show_warnings_once=True:  show api warnings for only one simulation
    - show_warnings_once=False: show api warnings for all simulations
    - show_warnings_once=None:  not show api warnings
    """
    platform = Platform('CALCULON', node_group='idm_48cores')
    # platform = Platform('IDMCLOUD', node_group='emod_abcd')

    # If you don't have Eradication, un-comment out the following to download Eradication
    # import emod_malaria.bootstrap as dtk
    # import pathlib
    # import os
    # dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    # os.chdir(os.path.dirname(__file__))
    # print("...done.")
    run_calibration(directory=manifest.directory, show_warnings_once=True)
