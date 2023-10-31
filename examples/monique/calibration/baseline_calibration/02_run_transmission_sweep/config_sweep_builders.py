import os
import pandas as pd
import numpy as np
import manifest
import params
from typing import List
from functools import partial
from idmtools.builders import SimulationBuilder
from idmtools.entities.simulation import Simulation
from snt.helpers_add_interventions import add_all_interventions
from snt.helpers_sim_setup import load_master_csv, habitat_scales, set_up_hfca
from snt.utility.sweeping import set_param, ItvFn, CfgFn
from snt.helpers_sim_setup import get_burnin_exp
from emodpy_malaria.reporters.builtin import add_report_event_counter

platform = None


def sweep_interventions(simulation: Simulation, func_list: List):
    """
    Sweeping on simulation.
    Args:
        simulation: idmtools Simulation
        func_list: a list of functions

    Returns:
        dict of parameters
    """
    tags_updated = {}
    for func in func_list:
        tags = func(simulation)
        if tags:
            if isinstance(func, ItvFn):
                fname = func.func.__name__
                if fname == 'add_all_interventions':
                    # add counter report
                    add_report_event_counter(simulation.task, manifest, event_trigger_list=tags["events"])
            else:
                tags_updated.update(tags)
    return tags_updated


###################################
# Common interface
###################################
def get_sweep_builders(**kwargs):
    global platform
    platform = kwargs.get('platform', None)
    builder = SimulationBuilder()

    burnin_id = params.burnin_id
    scen_df = params.scen_df
    project_path = manifest.project_path
    scen_index = params.scen_index

    ser_df = get_burnin_exp(burnin_id=burnin_id, platform=platform)

    # INTERVENTIONS
    # health-seeking
    if (not pd.isna(scen_df.at[scen_index, 'CM_filename'])) and (not (scen_df.at[scen_index, 'CM_filename'] == '')):
        hs_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'CM_filename']))
    else:
        hs_df = pd.DataFrame()
    # NMFs
    if (not pd.isna(scen_df.at[scen_index, 'NMF_filename'])) and (not (scen_df.at[scen_index, 'NMF_filename'] == '')):
        nmf_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'NMF_filename']))
    else:
        nmf_df = pd.DataFrame()
    # ITNs
    if (not pd.isna(scen_df.at[scen_index, 'ITN_filename'])) and (not (scen_df.at[scen_index, 'ITN_filename'] == '')):
        itn_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ITN_filename']))
    else:
        itn_df = pd.DataFrame()
    if (not pd.isna(scen_df.at[scen_index, 'ANC_ITN_filename'])) and (
            not (scen_df.at[scen_index, 'ANC_ITN_filename'] == '')):
        itn_anc_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ANC_ITN_filename']))
    else:
        itn_anc_df = pd.DataFrame()
    if (not pd.isna(scen_df.at[scen_index, 'EPI_ITN_filename'])) and (
            not (scen_df.at[scen_index, 'EPI_ITN_filename'] == '')):
        itn_epi_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'EPI_ITN_filename']))
    else:
        itn_epi_df = pd.DataFrame()

    itn_decay_params = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'itn_discard_decay_params.csv'))
    itn_use_seasonality = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'ITN_use_seasonality.csv'))

    # IRS
    if (not pd.isna(scen_df.at[scen_index, 'IRS_filename'])) and (not (scen_df.at[scen_index, 'IRS_filename'] == '')):
        irs_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'IRS_filename']))
    else:
        irs_df = pd.DataFrame()
    # SMC
    if (not pd.isna(scen_df.at[scen_index, 'SMC_filename'])) and (not (scen_df.at[scen_index, 'SMC_filename'] == '')):
        smc_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'SMC_filename']))
    else:
        smc_df = pd.DataFrame()

    df = load_master_csv(project_path=project_path)
    rel_abundance_df = habitat_scales(project_path=project_path)
    lhdf = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'larval_habitats',
                                    'monthly_habitats_1_%i.csv' % params.population_size))
    hfca_list = list(df.index.values)

    sweepspace = {
        rep_admin: [np.logspace(-1.6, 2, 25)[ii] for ii in [0, 4, 8, 12, 16, 20]] for rep_admin in
        lhdf['archetype'].values  # make sure this matches the values from serialize_transmission_sweep
    }

    ser_df.to_csv("ser_df.csv")
    funcs_list = [[CfgFn(set_up_hfca, manifest=manifest,
                         hfca=my_admin,
                         archetype_hfca=df.at[my_admin, 'seasonality_archetype'],
                         pull_from_serialization=params.pull_from_serialization,
                         ser_date=params.ser_date,
                         hdf=rel_abundance_df,
                         lhdf=lhdf,
                         population_size=params.population_size,
                         hab_multiplier=hab_scale,
                         run_number=(x % params.num_burnin_seeds),
                         use_arch_burnin=True,
                         ser_df=ser_df),
                   partial(set_param, param='Run_Number', value=x),
                   partial(set_param, param='Habitat_Multiplier', value=hab_scale),
                   ItvFn(add_all_interventions,
                         hfca=my_admin,
                         itn_anc_adult_birthday_years=params.itn_anc_adult_birthday_years,
                         itn_decay_params=itn_decay_params,
                         itn_use_seasonality=itn_use_seasonality,
                         seed_index=x + 1,
                         hs_df=hs_df,
                         nmf_df=nmf_df,
                         itn_df=itn_df,
                         itn_anc_df=itn_anc_df,
                         itn_epi_df=itn_epi_df,
                         irs_df=irs_df,
                         smc_df=smc_df)
                   ]
                  for x in range(params.num_seeds)
                  for my_admin in hfca_list
                  # for my_admin in ['Bukinanyana', 'Gisuru', 'Rutovu']
                  for hab_scale in sweepspace[df.at[my_admin, 'seasonality_archetype']]
                  ]

    builder.add_sweep_definition(sweep_interventions, funcs_list)
    print(f"builder.count: {builder.count}")

    return [builder]
