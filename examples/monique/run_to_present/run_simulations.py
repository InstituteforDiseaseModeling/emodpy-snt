import os
import sys
from pathlib import Path
import manifest
import params
import argparse
from idmtools.core.platform_factory import Platform

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snt.helpers_run_simulation import  run_experiment, log

tracking_file = os.path.join(manifest.CURRENT_DIR, "..", f"suite_tracking_{params.experiment_type}.csv")

def _print_params():
    """
    Just a useful convenient function for the user.
    """
    print("expname: ", params.expname)
    print("population_size: ", params.population_size)
    print("serialize: ", params.serialize)
    print("num_seeds: ", params.num_seeds)
    print("years: ", params.years)
    print("pull_from_serialization: ", params.pull_from_serialization)


if __name__ == "__main__":
    """
    - show_warnings_once=True:  show api warnings for only one simulation
    - show_warnings_once=False: show api warnings for all simulations
    - show_warnings_once=None:  not show api warnings
    """
    parser = argparse.ArgumentParser(description="Run experiment optionally using an existing suite_id")
    parser.add_argument('--suite-id', type=str, help='Optional suite ID to reuse or track')
    parser.add_argument('--scen-index', type=int, required=True, help='Index of scenario to run')
    parser.add_argument('--scenario-fname', type=str, required=True)
    parser.add_argument(
        '--show-warnings-once',
        type=str,
        choices=['True', 'False', 'None'],
        default='True',
        help='True: show warning once, False: for all, None: suppress warnings'
    )
    args = parser.parse_args()
    params.scen_index = args.scen_index
    params.scenario_fname = args.scenario_fname

    def str_to_bool_none(val):
        return {'True': True, 'False': False, 'None': None}[val]

    # Determine warning level
    show_warnings_once = str_to_bool_none(args.show_warnings_once)

    print("\n\n")
    print("======================================New Experiment========================================")
    log(f'Start experiment with suite_id: {args.suite_id}')
    platform = Platform('CALCULON')
    # If you don't have Eradication, un-comment out the following to download Eradication
    # import emod_malaria.bootstrap as dtk
    # import pathlib
    # import os
    # dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    # os.chdir(os.path.dirname(__file__))
    # print("...done.")
    run_experiment(
        platform=platform,
        params=params,
        manifest=manifest,
        print_params_fn=_print_params,
        suite_id=args.suite_id,
        tracking_file = tracking_file,
        show_warnings_once=show_warnings_once
    )
