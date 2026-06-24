import warnings
import os
import math
import pandas as pd
import numpy as np
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.usage_dependent_bednet import add_scheduled_usage_dependent_bednet, \
    add_triggered_usage_dependent_bednet
from emodpy_malaria.interventions.irs import add_scheduled_irs_housing_modification
from emodpy_malaria.interventions.larvicide import add_larvicide
from emod_api.interventions.common import change_individual_property_triggered
from emodpy_malaria.interventions.drug_campaign import add_drug_campaign
from emodpy_malaria.interventions.diag_survey import add_diagnostic_survey
from snt.support_files.malaria_vaccdrug_campaigns import add_vaccdrug_campaign
from emodpy_malaria.interventions.adherentdrug import adherent_drug
from snt.helpers_sim_setup import update_smc_access_ips, update_spaq_access_ips
from emodpy_malaria.interventions.vaccine import add_scheduled_vaccine, add_triggered_vaccine
from emodpy_malaria.interventions.common import add_triggered_campaign_delay_event, add_campaign_event
from emod_api.interventions.common import BroadcastEvent, PropertyValueChanger, DelayedIntervention, change_individual_property_scheduled

def tryread_intervention_csv_from_scen_df(project_path, scen_df, scen_index, intervention_colname):
    if ((intervention_colname not in scen_df.columns.values) or pd.isna(scen_df.at[scen_index, intervention_colname])) or ((scen_df.at[scen_index, intervention_colname] == '')):
        df = pd.DataFrame()
    else:
        try:
            df = pd.read_csv(os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, intervention_colname]))
        except IOError:
            print(f"WARNING: Cannot read intervention file for: {intervention_colname}.")
            df = pd.DataFrame()
    return df


def add_hfca_hs(campaign, hs_df, hfca, seed_index=0):
    # df = hs_df[hs_df['repDS'] == hfca]
    if 'LGA' in hs_df.columns.values:
        df = hs_df[hs_df['LGA'] == hfca]
    elif 'repDS' in hs_df.columns.values:
        df = hs_df[hs_df['repDS'] == hfca]
    else:
        df = hs_df[hs_df['admin_name'] == hfca]
    if 'seed' in df.columns.values:
        df = df[df['seed'] == seed_index]
    for r, row in df.iterrows():
        add_hs_from_file(campaign, row)

    return len(df)


def add_hs_from_file(campaign, row):
    hs_child = row['U5_coverage']
    hs_adult = row['adult_coverage']
    start_day = row['simday']
    severe_cases = row['severe_coverage']

    add_treatment_seeking(campaign, start_day=start_day,
                          targets=[{'trigger': 'NewClinicalCase', 'coverage': hs_child, 'agemin': 0, 'agemax': 5,
                                    'rate': 0.33},
                                   {'trigger': 'NewClinicalCase', 'coverage': hs_adult, 'agemin': 5, 'agemax': 100,
                                    'rate': 0.33},
                                   ],
                          drug=['Artemether', 'Lumefantrine'], duration=row['duration'])
    add_treatment_seeking(campaign, start_day=start_day,
                          targets=[{'trigger': 'NewSevereCase', 'coverage': severe_cases, 'rate': 0.5}],
                          # change by adding column and reviewing literature
                          drug=['Artemether', 'Lumefantrine'], duration=row['duration'],
                          broadcast_event_name='Received_Severe_Treatment')


def add_nmf_hs(campaign, hs_df, nmf_df, hfca, seed_index=0):
    # if no NMF rate is specified, assume all age groups have 0.0038 probability each day
    if nmf_df.empty:
        nmf_df = pd.DataFrame({'U5_nmf': [0.0038], 'adult_nmf': [0.0038]})
    elif nmf_df.shape[0] != 1:
        warnings.warn('The NMF dataframe has more than one row. Only values in the first row will be used.')
    nmf_row = nmf_df.iloc[0]

    # apply the health-seeking rate for clinical malaria to NMFs
    df = hs_df[hs_df['admin_name'] == hfca]
    if 'seed' in df.columns.values:
        df = df[df['seed'] == seed_index]
    for r, row in df.iterrows():
        add_nmf_hs_from_file(campaign, row, nmf_row)

    return len(df)


def add_nmf_hs_from_file(campaign, row, nmf_row):
    hs_child = row['U5_coverage']
    hs_adult = row['adult_coverage']
    start_day = row['simday']
    duration = row['duration']
    if start_day == 0:  # due to dtk diagnosis/treatment configuration, a start day of 0 is not supported
        start_day = 1  # start looking for NMFs on day 1 (not day 0) of simulation
        if duration > 1:
            duration = duration - 1
    nmf_child = nmf_row['U5_nmf']
    nmf_adult = nmf_row['adult_nmf']
    if duration < 1000:
        start_day_list = [start_day]
        duration_list = [duration]
    else:
        interval = 900  # should be less than the EMOD max (which is currently 1000)
        start_day_list_0 = list(range(0, duration, interval))
        start_day_list = [start_day + ss for ss in start_day_list_0]
        duration_list = [interval if ((ss + interval) < duration) else (duration - ss) for ss in start_day_list_0]

    for ii in range(len(start_day_list)):
        add_drug_campaign(campaign, 'MSAT', 'AL', start_days=[start_day_list[ii]],
                          target_group={'agemin': 0, 'agemax': 5},
                          coverage=nmf_child * hs_child,
                          repetitions=duration_list[ii], tsteps_btwn_repetitions=1,
                          diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                          receiving_drugs_event_name='Received_NMF_Treatment')
        add_drug_campaign(campaign, 'MSAT', 'AL', start_days=[start_day_list[ii]],
                          target_group={'agemin': 5, 'agemax': 120},
                          coverage=nmf_adult * hs_adult,
                          repetitions=duration_list[ii], tsteps_btwn_repetitions=1,
                          diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                          receiving_drugs_event_name='Received_NMF_Treatment')


def add_hfca_irs(campaign, irs_df, hfca, seed_index=0):
    irs_df = irs_df[irs_df['admin_name'].str.upper() == hfca.upper()]
    if 'seed' in irs_df.columns.values:
        irs_df = irs_df[irs_df['seed'] == seed_index]
    for r, row in irs_df.iterrows():
        add_scheduled_irs_housing_modification(campaign, start_day=row['simday'],
                                               demographic_coverage=row['effective_coverage'],
                                               target_age_min=0,
                                               target_age_max=100,
                                               killing_initial_effect=row['initial_kill'],
                                               killing_decay_time_constant=row['mean_duration'])  # !! this was 'initial_kill' before, but I think that was an error in the translation from dtktools?

    return len(irs_df)


# larviciding
def add_hfca_lsm(campaign, lsm_df, hfca, seed_index=0):
    lsm_df = lsm_df[lsm_df['admin_name'].str.upper() == hfca.upper()]
    if 'seed' in lsm_df.columns.values:
        lsm_df = lsm_df[lsm_df['seed'] == seed_index]
    for r, row in lsm_df.iterrows():
        # note that the add_larvicide() function has an error such that no matter what you specify for num_repetitions and timesteps_between_reps, it only delivers the intervention once,
        #   so will add each campaign separately here until that is fixed
        for ii in range(row['num_repetitions']):
            add_larvicide(campaign, start_day=row['simday'] + ii * row['timesteps_between_reps'],
                          spray_coverage=row['effective_coverage'],
                          killing_effect=row['killing_effect'],
                          habitat_target="ALL_HABITATS",
                          box_duration=row['box_duration'],
                          decay_time_constant=row['decay_time_constant'])

    return len(lsm_df)


# all itn
def add_hfca_itns(campaign, itn_df, itn_anc_df, itn_anc_adult_birthday_years, itn_epi_df, itn_cont_df, itn_chw_df, itn_chw_annual_df,
                  hfca,
                  itn_use_seasonality, itn_decay_params, seed_index=0):  #
    if not itn_anc_adult_birthday_years:
        itn_anc_adult_birthday_years = []
    if not itn_df.empty:
        df = itn_df[itn_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        nets = len(df)
        for r, row in df.iterrows():
            add_itn_from_file(campaign, row, itn_use_seasonality, itn_decay_params)
    else:
        nets = 0

    if not itn_anc_df.empty:
        df = itn_anc_df[itn_anc_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        add_itn_anc(campaign, df, itn_anc_adult_birthday_years, itn_use_seasonality, itn_decay_params)
        nets += len(df)

    if not itn_epi_df.empty:
        df = itn_epi_df[itn_epi_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        add_birthday_routine_itn_from_file(campaign, df, itn_use_seasonality, itn_decay_params)
        nets += len(df)

    if not itn_cont_df.empty:
        df = itn_cont_df[itn_cont_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        add_birthday_routine_itn_from_file(campaign, df, itn_use_seasonality, itn_decay_params)
        nets += len(df)

    if not itn_chw_df.empty:
        df = itn_chw_df[itn_chw_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        nets += len(df)
        for r, row in df.iterrows():
            add_monthly_chw_dist(campaign, row, itn_use_seasonality, itn_decay_params)

    if not itn_chw_annual_df.empty:
        df = itn_chw_annual_df[itn_chw_annual_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        nets += len(df)
        for r, row in df.iterrows():
            add_annual_chw_dist(campaign, row, itn_use_seasonality, itn_decay_params)

    return nets


def add_itn_from_file(campaign, row, itn_use_seasonality, itn_decay_params):
    # net retention time, decay rates of killing and blocking
    itn_discard_config = {"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                          "Expiration_Period_Log_Normal_Mu": row['net_life_lognormal_mu'],
                          "Expiration_Period_Log_Normal_Sigma": row['net_life_lognormal_sigma']}
    itn_decay_kill = itn_decay_params['kill_decay_time'][0].item()
    itn_decay_block = itn_decay_params['block_decay_time'][0].item()

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']
    seasonal_offset = row['simday'] % 365
    seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

    zipped_lists = zip(seasonal_times, seasonal_scales)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
    if seasonal_times[0] > 0:
        seasonal_times.insert(0, 0)
        seasonal_scales.insert(0, seasonal_scales[-1])

    # use-coverage by age
    itn_u5 = row['itn_u5']
    itn_5_10 = row['itn_5_10']
    itn_10_15 = row['itn_10_15']
    itn_15_20 = row['itn_15_20']
    itn_o20 = row['itn_o20']
    # a single distribution coverage is given to all age groups, but use differs by age.
    #   set the distribution coverage to be the maximum use coverage across all age groups, then use age_dependence
    #   uses to end up with the appropriate use for that age group
    coverage_all = max([itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20])
    if coverage_all > 0:
        age_dep = [x / coverage_all for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]
    else:
        age_dep = [1 for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]

    # fraction of indoor time protected by net
    indoor_net_protection = row['indoor_net_protection']

    add_scheduled_usage_dependent_bednet(campaign, start_day=row['simday'], demographic_coverage=coverage_all,
                                         killing_initial_effect=row['kill_initial'],
                                         killing_decay_time_constant=itn_decay_kill,
                                         blocking_initial_effect=row['block_initial'],
                                         blocking_decay_time_constant=itn_decay_block,
                                         discard_config=itn_discard_config,
                                         age_dependence={'Times': [0, 5, 10, 15, 20],
                                                         'Values': [x * indoor_net_protection for x in age_dep]},
                                         seasonal_dependence={"Times": seasonal_times,
                                                              "Values": seasonal_scales})
    # it says simday == 895, but intervention definitely starts doing its thing on start_day 0, is it
    # just that the data is in simday 895?
    if 'ITN coverage severe treatment' in row:
        if row['simday'] == 895:
            add_triggered_usage_dependent_bednet(campaign, start_day=0,
                                                 demographic_coverage=row['ITN coverage severe treatment'],
                                                 trigger_condition_list=['Received_Severe_Treatment'],
                                                 listening_duration=-1,
                                                 killing_initial_effect=row['kill_initial'],
                                                 killing_decay_time_constant=itn_decay_kill,
                                                 blocking_initial_effect=row['block_initial'],
                                                 blocking_decay_time_constant=itn_decay_block,
                                                 discard_config=itn_discard_config,
                                                 age_dependence={'Times': [0, 5, 10, 15, 20],
                                                                 'Values': [x * indoor_net_protection for x in
                                                                            age_dep]},
                                                 seasonal_dependence={"Times": seasonal_times,
                                                                      "Values": seasonal_scales})
    # st question as above
    if 'ITN coverage uncomplicated treatment' in row:
        if row['simday'] == 895:
            add_triggered_usage_dependent_bednet(campaign, start_day=0,
                                                 demographic_coverage=row['ITN coverage uncomplicated treatment'],
                                                 trigger_condition_list=['Received_Treatment'],
                                                 listening_duration=-1,
                                                 killing_initial_effect=row['kill_initial'],
                                                 killing_decay_time_constant=itn_decay_kill,
                                                 blocking_initial_effect=row['block_initial'],
                                                 blocking_decay_time_constant=itn_decay_block,
                                                 discard_config=itn_discard_config,
                                                 age_dependence={'Times': [0, 5, 10, 15, 20],
                                                                 'Values': [x * indoor_net_protection for x in
                                                                            age_dep]},
                                                 seasonal_dependence={"Times": seasonal_times,
                                                                      "Values": seasonal_scales})


def add_birthday_routine_itn_from_file(campaign, itn_epi_df, itn_use_seasonality, itn_decay_params):
    # decay rates of killing and blocking
    itn_decay_kill = itn_decay_params['kill_decay_time'][0].item()
    itn_decay_block = itn_decay_params['block_decay_time'][0].item()

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']

    for r, row in itn_epi_df.iterrows():
        # net retention time
        itn_discard_config = {"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                              "Expiration_Period_Log_Normal_Mu": row['net_life_lognormal_mu'],
                              "Expiration_Period_Log_Normal_Sigma": row['net_life_lognormal_sigma']}
        indoor_net_protection = row['indoor_net_protection']

        # age of individuals receiving nets
        birthday_age = row['birthday_age']

        seasonal_offset = row['simday'] % 365
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0:
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        add_triggered_usage_dependent_bednet(campaign, start_day=row['simday'],
                                             demographic_coverage=row['coverage'],
                                             trigger_condition_list=['HappyBirthday'],
                                             listening_duration=row['duration'],
                                             target_age_min=(birthday_age - 0.5),
                                             target_age_max=(birthday_age + 0.5),
                                             killing_initial_effect=row['kill_initial'],
                                             killing_decay_time_constant=itn_decay_kill,
                                             blocking_initial_effect=row['block_initial'],
                                             blocking_decay_time_constant=itn_decay_block,
                                             discard_config=itn_discard_config,
                                             age_dependence={'Times': [0, 120],
                                                             'Values': [indoor_net_protection, indoor_net_protection]},
                                             seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales})


def add_itn_anc(campaign, itn_anc_df, itn_anc_adult_birthday_years, itn_use_seasonality, itn_decay_params):
    # decay rates of killing and blocking
    itn_decay_kill = itn_decay_params['kill_decay_time'][0].item()
    itn_decay_block = itn_decay_params['block_decay_time'][0].item()

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']

    for r, row in itn_anc_df.iterrows():
        # net retention time
        itn_discard_config = {"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                              "Expiration_Period_Log_Normal_Mu": row['net_life_lognormal_mu'],
                              "Expiration_Period_Log_Normal_Sigma": row['net_life_lognormal_sigma']}
        indoor_net_protection = row['indoor_net_protection']

        seasonal_offset = row['simday'] % 365
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0:
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        # ITNs protecting infants
        add_triggered_usage_dependent_bednet(campaign, start_day=row['simday'],
                                             demographic_coverage=row['coverage'],
                                             trigger_condition_list=['Births'],
                                             listening_duration=row['duration'],
                                             killing_initial_effect=row['kill_initial'],
                                             killing_decay_time_constant=itn_decay_kill,
                                             blocking_initial_effect=row['block_initial'],
                                             blocking_decay_time_constant=itn_decay_block,
                                             discard_config=itn_discard_config,
                                             age_dependence={'Times': [0, 5],
                                                             'Values': [indoor_net_protection, indoor_net_protection]},
                                             seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales})

        # ITNs protecting mothers (since we don't simulate pregnancy explicitly, we approximate ANC ITN coverage among
        #    women by giving them ITNs at certain birthdays, using the ANC ITN coverage for that admin)
        for birthday_year in itn_anc_adult_birthday_years:
            add_triggered_usage_dependent_bednet(campaign, start_day=row['simday'],
                                                 demographic_coverage=row['coverage'],
                                                 target_age_min=(birthday_year - 0.5),
                                                 target_age_max=(birthday_year + 0.5),
                                                 target_gender="Female",
                                                 trigger_condition_list=["HappyBirthday"],
                                                 killing_initial_effect=row['kill_initial'],
                                                 killing_decay_time_constant=itn_decay_kill,
                                                 blocking_initial_effect=row['block_initial'],
                                                 blocking_decay_time_constant=itn_decay_block,
                                                 discard_config=itn_discard_config,
                                                 age_dependence={'Times': [0, 120],
                                                                 'Values': [indoor_net_protection,
                                                                            indoor_net_protection]},
                                                 seasonal_dependence={"Times": seasonal_times,
                                                                      "Values": seasonal_scales})


def add_monthly_chw_dist(campaign, row, itn_use_seasonality, itn_decay_params):
    # net retention time, decay rates of killing and blocking
    itn_discard_config = {"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                          "Expiration_Period_Log_Normal_Mu": row['net_life_lognormal_mu'],
                          "Expiration_Period_Log_Normal_Sigma": row['net_life_lognormal_sigma']}
    itn_decay_kill = itn_decay_params['kill_decay_time'][0].item()
    itn_decay_block = itn_decay_params['block_decay_time'][0].item()

    # use-coverage by age
    itn_u5 = row['itn_u5']
    itn_5_10 = row['itn_5_10']
    itn_10_15 = row['itn_10_15']
    itn_15_20 = row['itn_15_20']
    itn_o20 = row['itn_o20']
    # a single distribution coverage is given to all age groups, but use differs by age.
    #   set the distribution coverage to be the maximum use coverage across all age groups, then use age_dependence
    #   uses to end up with the appropriate use for that age group
    coverage_all = max([itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20])
    age_dep = [x / coverage_all for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]

    # fraction of indoor time protected by net
    indoor_net_protection = row['indoor_net_protection']

    # add distribution every month for a year (every 30 days starting from simday)
    simday = row['simday']
    for mm in range(12):
        chw_check_day = simday + 30 * mm

        # seasonality in ITN use
        seasonal_scales = itn_use_seasonality['itn_use_scalar']
        seasonal_days = itn_use_seasonality['day']
        seasonal_offset = chw_check_day
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0:
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        add_scheduled_usage_dependent_bednet(campaign, start_day=chw_check_day, demographic_coverage=coverage_all,
                                             killing_initial_effect=row['kill_initial'],
                                             killing_decay_time_constant=itn_decay_kill,
                                             blocking_initial_effect=row['block_initial'],
                                             blocking_decay_time_constant=itn_decay_block,
                                             discard_config=itn_discard_config,
                                             age_dependence={'Times': [0, 5, 10, 15, 20],
                                                             'Values': [x * indoor_net_protection for x in age_dep]},
                                             seasonal_dependence={"Times": seasonal_times,
                                                                  "Values": seasonal_scales},
                                             dont_allow_duplicates=True
                                             )


def add_annual_chw_dist(campaign, row, itn_use_seasonality, itn_decay_params):
    itn_decay_kill = itn_decay_params['kill_decay_time'][0].item()
    itn_decay_block = itn_decay_params['block_decay_time'][0].item()

    # use-coverage by age
    itn_u5 = row['itn_u5']
    itn_5_10 = row['itn_5_10']
    itn_10_15 = row['itn_10_15']
    itn_15_20 = row['itn_15_20']
    itn_o20 = row['itn_o20']
    # a single distribution coverage is given to all age groups, but use differs by age.
    #   set the distribution coverage to be the maximum use coverage across all age groups, then use age_dependence
    #   uses to end up with the appropriate use for that age group
    coverage_all = max([itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20])
    age_dep = [x / coverage_all for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]

    # fraction of indoor time protected by net
    indoor_net_protection = row['indoor_net_protection']

    # add distribution every month for a year (every 30 days starting from simday)
    simday = row['simday']

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']
    seasonal_offset = simday % 365
    seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

    zipped_lists = zip(seasonal_times, seasonal_scales)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
    if seasonal_times[0] > 0:
        seasonal_times.insert(0, 0)
        seasonal_scales.insert(0, seasonal_scales[-1])

    add_scheduled_usage_dependent_bednet(campaign, start_day=simday, demographic_coverage=coverage_all,
                                         killing_initial_effect=row['kill_initial'],
                                         killing_decay_time_constant=itn_decay_kill,
                                         blocking_initial_effect=row['block_initial'],
                                         blocking_decay_time_constant=itn_decay_block,
                                         discard_config={"Expiration_Period_Distribution": "CONSTANT_DISTRIBUTION",
                                                         "Expiration_Period_Constant": 365},
                                         age_dependence={'Times': [0, 5, 10, 15, 20],
                                                         'Values': [x * indoor_net_protection for x in age_dep]},
                                         seasonal_dependence={"Times": seasonal_times,
                                                              "Values": seasonal_scales})


def add_hfca_vaccsmc(campaign, smc_df, hfca, effective_coverage_resistance_multiplier=1, seed_index=0):
    df = smc_df[smc_df['admin_name'] == hfca]
    if df.shape[0] > 0:
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]

    if len(df) == 0:
        return len(df)

    for r, row in df.iterrows():
        if 'max_age' in df.columns.values:
            max_smc_age = row['max_age']
        else:
            max_smc_age = 5
        if 'TAT' in df.columns.values:
            tat = row['TAT']
        else:
            tat = 0
        if tat:
            raise ValueError('TAT not yet implemented for vaccSMC')

        # get coverage and ages for SMC in different groups
        cov_high_u5 = row['coverage_high_access_U5'] * effective_coverage_resistance_multiplier
        cov_low_u5 = row['coverage_low_access_U5'] * effective_coverage_resistance_multiplier
        if 'coverage_high_access_5_10' in df.columns.values:
            cov_high_5_10 = row['coverage_high_access_5_10'] * effective_coverage_resistance_multiplier
            cov_low_5_10 = row['coverage_low_access_5_10'] * effective_coverage_resistance_multiplier
            # create list of values for all groups receiving SMC with order: [U5-high, U5-low, O5-high, O5-low]
            age_min_list = [0.25, 0.25, 5, 5]
            age_max_list = [5, 5, max(10, max_smc_age), max(10, max_smc_age)]
            cov_list = [cov_high_u5, cov_low_u5, cov_high_5_10, cov_low_5_10]
            access_list = ['High', 'Low', 'High', 'Low']
        else:
            # create list of values for all groups receiving SMC with order: [U5-high, U5-low, O5-high, O5-low]
            age_min_list = [0.25, 0.25]
            age_max_list = [5, 5]
            cov_list = [cov_high_u5, cov_low_u5]
            access_list = ['High', 'Low']

        for ii in range(len(cov_list)):
            add_vaccdrug_campaign(campaign, campaign_type='SMC', start_days=[row['simday']],
                                  coverages=cov_list[ii],
                                  target_group={'agemin': age_min_list[ii],
                                                'agemax': age_max_list[ii]},
                                  ind_property_restrictions=[{'SMCAccess': access_list[ii]}],
                                  receiving_drugs_event=False)  ## If False uses vaccSMC with automatic offset of 17 days, if True, uses vaccDrugSMC

    return len(df)


def smc_adherent_configuration(campaign, adherence, sp_resist_day1_multiply):
    smc_adherent_config = adherent_drug(campaign,
                                        doses=[["SulfadoxinePyrimethamine", 'Amodiaquine'],
                                               ['Amodiaquine'],
                                               ['Amodiaquine']],
                                        dose_interval=1,
                                        non_adherence_options=['Stop'],
                                        non_adherence_distribution=[1],
                                        adherence_values=[
                                            sp_resist_day1_multiply,  # for day 1
                                            adherence,  # day 2
                                            adherence  # day 3
                                        ]
                                        )
    return smc_adherent_config


def add_hfca_smc(campaign, smc_df, hfca, adherence_multiplier=1, sp_resist_day1_multiply=1, seed_index=0):
    df = smc_df[smc_df['admin_name'] == hfca]
    if df.shape[0] > 0:
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        if 'adherence' in smc_df.columns.values:
            adherent_drug_configs = smc_adherent_configuration(campaign,
                                                               adherence=df['adherence'].values[
                                                                             0] * adherence_multiplier,
                                                               sp_resist_day1_multiply=sp_resist_day1_multiply)
        else:
            default_adherence = 0.8
            adherent_drug_configs = smc_adherent_configuration(campaign,
                                                               adherence=default_adherence * adherence_multiplier,
                                                               sp_resist_day1_multiply=sp_resist_day1_multiply)
        for r, row in df.iterrows():
            cov_high_u5 = row['coverage_high_access_U5']
            cov_low_u5 = row['coverage_low_access_U5']
            cov_high_5_10 = row['coverage_high_access_5_10']
            cov_low_5_10 = row['coverage_low_access_5_10']
            if 'max_age' in smc_df.columns.values:
                max_smc_age = row['max_age']
            else:
                max_smc_age = 5
            if 'TAT' in smc_df.columns.values:
                tat = row['TAT']
            else:
                tat = 0

            # set diagnostic so that no one should be excluded from SMC due to fever
            add_diagnostic_survey(campaign, start_day=row['simday'],
                                  coverage=1,
                                  target={"agemin": 0.25, "agemax": max(10, max_smc_age)},
                                  diagnostic_type='FEVER',
                                  diagnostic_threshold=0.5,
                                  # SVET - do these configs works? seems like we should make an actual broadcast
                                  # event, but should give it a shot
                                  negative_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever",
                                      "class": "BroadcastEvent"}],
                                  positive_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever",
                                      # currently set so that no one should be excluded from SMC due to fever,
                                      # to include fever discrimination, change to "Has_SMC_Fever"
                                      "class": "BroadcastEvent"}]
                                  )
            # U5 - high access
            add_drug_campaign(campaign, 'SMC', start_days=[row['simday']],
                              coverage=cov_high_u5, target_group={'agemin': 0.25, 'agemax': 5},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'High'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # U5 - low access
            add_drug_campaign(campaign, 'SMC', start_days=[row['simday']],
                              coverage=cov_low_u5, target_group={'agemin': 0.25, 'agemax': 5},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'Low'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # U10 (or older) - high access  (
            add_drug_campaign(campaign, 'SMC', start_days=[row['simday']],
                              coverage=cov_high_5_10, target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'High'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # U10 (or older) - low access
            add_drug_campaign(campaign, 'SMC', start_days=[row['simday']],
                              coverage=cov_low_5_10, target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'Low'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )

            if tat:
                add_drug_campaign(campaign, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_high_u5,
                                  target_group={'agemin': 0.25, 'agemax': 5},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess': 'High'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
                add_drug_campaign(campaign, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_low_u5,
                                  target_group={'agemin': 0.25, 'agemax': 5},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess': 'Low'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
                add_drug_campaign(campaign, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_high_u5,
                                  target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess': 'High'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
                add_drug_campaign(campaign, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_low_u5,
                                  target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess': 'Low'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
    return len(df)




def add_hfca_spaq(campaign, spaq_df, hfca, adherence_multiplier=1, sp_resist_day1_multiply=1, seed_index=0):
    df = spaq_df[spaq_df['admin_name'] == hfca]
    if df.shape[0] > 0:
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        if 'adherence' in spaq_df.columns.values:
            adherent_drug_configs = smc_adherent_configuration(campaign,
                                                               adherence=df['adherence'].values[
                                                                             0] * adherence_multiplier,
                                                               sp_resist_day1_multiply=sp_resist_day1_multiply)
        else:
            default_adherence = 0.8
            adherent_drug_configs = smc_adherent_configuration(campaign,
                                                               adherence=default_adherence * adherence_multiplier,
                                                               sp_resist_day1_multiply=sp_resist_day1_multiply)
        for r, row in df.iterrows():
            cov_high = row['coverage_high_access']
            cov_low = row['coverage_low_access']
            min_age = row['min_age']
            max_age = row['max_age']

            # set diagnostic so that no one should be excluded from SMC due to fever
            add_diagnostic_survey(campaign, start_day=row['simday'],
                                  coverage=1,
                                  target={"agemin": min_age, "agemax": max_age},
                                  diagnostic_type='FEVER',
                                  diagnostic_threshold=0.5,
                                  # SVET - do these configs works? seems like we should make an actual broadcast
                                  # event, but should give it a shot
                                  negative_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever",
                                      "class": "BroadcastEvent"}],
                                  positive_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever",
                                      # currently set so that no one should be excluded from SMC due to fever,
                                      # to include fever discrimination, change to "Has_SMC_Fever"
                                      "class": "BroadcastEvent"}]
                                  )
            # High access
            add_drug_campaign(campaign, 'SMC', start_days=[row['simday']],
                              coverage=cov_high, target_group={'agemin': min_age, 'agemax': max_age},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'High'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # Low access
            add_drug_campaign(campaign, 'SMC', start_days=[row['simday']],
                              coverage=cov_low, target_group={'agemin': min_age, 'agemax': max_age},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'Low'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )

    return len(df)





def calc_high_low_access_coverages(coverage_all, high_access_frac):
    if (high_access_frac < 1) & (coverage_all >= high_access_frac):
        coverage_high = 1
        coverage_low = (coverage_all - high_access_frac) / (1 - high_access_frac)
    else:
        coverage_high = coverage_all / high_access_frac
        coverage_low = 0
    return [coverage_high, coverage_low]


def change_rtss_ips(campaign):
    # SVET - blackout is not used, but default within EMOD is to not have blackout period, so person will
    # be getting multiple of this intervention a day if they received multiple vaccines
    # made issue: https://github.com/InstituteforDiseaseModeling/emod-api/issues/630
    change_individual_property_triggered(campaign,
                                         triggers=['Received_Vaccine'],
                                         new_ip_key='VaccineStatus',
                                         new_ip_value='GotVaccine',
                                         ip_restrictions=[{'VaccineStatus': 'None'}],
                                         blackout=False)

    change_individual_property_triggered(campaign,
                                         triggers=['Received_Vaccine'],
                                         new_ip_key='VaccineStatus',
                                         new_ip_value='GotBooster1',
                                         ip_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                         blackout=False)

    change_individual_property_triggered(campaign,
                                         triggers=['Received_Vaccine'],
                                         new_ip_key='VaccineStatus',
                                         new_ip_value='GotBooster2',
                                         ip_restrictions=[{'VaccineStatus': 'GotBooster1'}],
                                         blackout=False)


def add_epi_rtss(campaign, rtss_df):
    start_days = list(rtss_df['RTSS_day'])
    coverage_levels = list(rtss_df['coverage'])
    rtss_types = list(rtss_df['vaccine'])
    rtss_touchpoints = list(rtss_df['rtss_touchpoints'])
    rtss_event_names = [f'RTSS_{x + 1}_eligible' for x in range(len(rtss_touchpoints))]  # okay to have more than one each for primary or booster, just distinguishes them for broadcast events

    delay_distribution_name = list(rtss_df['distribution_name'])[0]
    std_dev_list = list(rtss_df['distribution_std'])

    initial_effect_list = list(rtss_df['initial_effect'])
    decay_time_constant_list = list(rtss_df['decay_time_constant'])
    duration_list = list(rtss_df['duration']) if 'duration' in rtss_df.columns else [-1] * len(rtss_df)

    for tp_time_trigger, coverage, vtype, event_name, std, init_eff, decay_t, start_day, duration in \
            zip(rtss_touchpoints, coverage_levels, rtss_types, rtss_event_names, std_dev_list,
                initial_effect_list, decay_time_constant_list, start_days, duration_list):

        if delay_distribution_name == "LOG_NORMAL_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                  "Delay_Period_Log_Normal_Mu": (
                                          np.log(tp_time_trigger + 14) - ((1 / 2) * std ** 2)),
                                  "Delay_Period_Log_Normal_Sigma": std}
        elif delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                  "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                  "Delay_Period_Gaussian_Std_Dev": std}
        else:  # no delay
            delay_distribution = {"Delay_Period_Distribution": "CONSTANT_DISTRIBUTION",
                                  "Delay_Period_Constant": tp_time_trigger}

        # triggered_campaign_delay_event only has option for constant delay, but we need different
        # distributions, so we're manually creating a delayed intervention that broadcasts an event
        # and slipping it into the triggered intervention
        broadcast_event = BroadcastEvent(campaign, event_name)
        if delay_distribution:
            broadcast_event = DelayedIntervention(campaign, Configs=[broadcast_event],
                                                  Delay_Dict=delay_distribution)
        add_triggered_campaign_delay_event(campaign, start_day=start_day,
                                           trigger_condition_list=['Births'],
                                           demographic_coverage=coverage,
                                           listening_duration=duration,
                                           individual_intervention=broadcast_event)

        # TODO: Make EPI support booster1 and booster2
        if not vtype == 'booster':
            add_triggered_vaccine(campaign,
                                  start_day=start_day,
                                  trigger_condition_list=[event_name],
                                  listening_duration=duration,
                                  intervention_name='RTSS',
                                  broadcast_event='Received_Vaccine',
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=init_eff,
                                  vaccine_box_duration=0,
                                  vaccine_decay_time_constant=decay_t,
                                  efficacy_is_multiplicative=False)
        else:
            add_triggered_vaccine(campaign,
                                  start_day=start_day,
                                  trigger_condition_list=[event_name],
                                  listening_duration=duration,
                                  ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                  intervention_name='RTSS',
                                  broadcast_event='Received_Vaccine',
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=init_eff,
                                  vaccine_box_duration=0,
                                  vaccine_decay_time_constant=decay_t,
                                  efficacy_is_multiplicative=False
                                  )


def add_campaign_rtss(campaign, rtss_df):
    for r, row in rtss_df.iterrows():
        vtype = row['vaccine']
        broadcast = 'Received_Vaccine'
        if vtype == 'booster1':
            add_scheduled_vaccine(campaign,
                                  start_day=row['RTSS_day'],
                                  demographic_coverage=row['coverage'],
                                  ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                  repetitions=row['repetitions'],
                                  timesteps_between_repetitions=row['tsteps_btwn_repetitions'],
                                  target_age_min=row['agemin'],
                                  target_age_max=row['agemax'],
                                  broadcast_event=broadcast,
                                  intervention_name='RTSS',
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=row['initial_effect'],
                                  vaccine_box_duration=0,
                                  vaccine_decay_time_constant=row['decay_time_constant'],
                                  efficacy_is_multiplicative=False)

        elif vtype == 'booster2':
            add_scheduled_vaccine(campaign,
                                  start_day=row['RTSS_day'],
                                  demographic_coverage=row['coverage'],
                                  ind_property_restrictions=[{'VaccineStatus': 'GotBooster1'}],
                                  repetitions=row['repetitions'],
                                  timesteps_between_repetitions=row['tsteps_btwn_repetitions'],
                                  target_age_min=row['agemin'],
                                  target_age_max=row['agemax'],
                                  broadcast_event=broadcast,
                                  intervention_name='RTSS',
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=row['initial_effect'],
                                  vaccine_box_duration=0,
                                  vaccine_decay_time_constant=row['decay_time_constant'],
                                  efficacy_is_multiplicative=False)
        elif vtype == 'booster3':
            add_scheduled_vaccine(campaign,
                                  start_day=row['RTSS_day'],
                                  demographic_coverage=row['coverage'],
                                  ind_property_restrictions=[{'VaccineStatus': 'GotBooster2'}],
                                  repetitions=row['repetitions'],
                                  timesteps_between_repetitions=row['tsteps_btwn_repetitions'],
                                  target_age_min=row['agemin'],
                                  target_age_max=row['agemax'],
                                  broadcast_event=broadcast,
                                  intervention_name='RTSS',
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=row['initial_effect'],
                                  vaccine_box_duration=0,
                                  vaccine_decay_time_constant=row['decay_time_constant'],
                                  efficacy_is_multiplicative=False)
        else:
            add_scheduled_vaccine(campaign,
                                  start_day=row['RTSS_day'],
                                  demographic_coverage=row['coverage'],
                                  repetitions=row['repetitions'],
                                  timesteps_between_repetitions=row['tsteps_btwn_repetitions'],
                                  target_age_min=row['agemin'],
                                  target_age_max=row['agemax'],
                                  broadcast_event=broadcast,
                                  intervention_name='RTSS',
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=row['initial_effect'],
                                  vaccine_box_duration=0,
                                  vaccine_decay_time_constant=row['decay_time_constant'],
                                  efficacy_is_multiplicative=False)


def add_ds_rtss(campaign, rtss_df, hfca):
    change_ips = False
    rtss_df = rtss_df[rtss_df['admin_name'].str.upper() == hfca.upper()]
    # First, process EPI style distribution
    rtss_df1 = rtss_df[rtss_df['deploy_type'] == 'EPI']
    if len(rtss_df1) > 0:
        add_epi_rtss(campaign, rtss_df1)
        change_ips = True

    # Second, process campaign style distribution
    rtss_df2 = rtss_df[rtss_df['deploy_type'] == 'campaign']
    if len(rtss_df2) > 0:
        add_campaign_rtss(campaign, rtss_df2)
        change_ips = True

    if change_ips:
        change_rtss_ips(campaign)

    return len(rtss_df)


def add_ds_vaccpmc(campaign, pmc_df, hfca):
    df = pmc_df[pmc_df['admin_name'].str.upper() == hfca.upper()]
    if len(df) == 0:
        return 0
    if "num_IIV_groups" in pmc_df:
        num_iiv_groups = pmc_df['num_IIV_groups'].unique()[0]
    else:
        num_iiv_groups = 1

    pmc_touchpoints_dict = {}
    for i, tp in enumerate(df['pmc_touchpoints']):
        pmc_touchpoints_dict[f'{i}'] = tp

    listening_durations = list(df['duration']) if 'duration' in df.columns else [-1] * len(df)
    add_vaccdrug_campaign(campaign, campaign_type='PMC', start_days=list(df['simday']),
                          coverages=df['coverage'],
                          listening_durations=listening_durations,
                          target_group=pmc_touchpoints_dict,
                          vaccine_param_dict={'vacc_initial_effect': 0.80, 'vacc_box_duration': 20,  # old: 32, Manuela-updated for NGA: 20
                                              'vacc_decay_duration': 9},  # old: 10, Manuela-updated for NGA: 9
                          delay_distribution_dict={'delay_distribution_name': df['distribution_name'],
                                                   'delay_distribution_mean': df['distribution_mean'],
                                                   'delay_distribution_std': df['distribution_std']},
                          num_iiv_groups=num_iiv_groups,
                          receiving_drugs_event=False)  ## use vaccine effects only with default offset of -10 days

    return len(pmc_df)









#####################
# ################### vaccine campaigns created with triggered event  ############
#####################

def get_concentration_at_time(tt, initial_concentration, fast_frac, k1, k2):
    # calculate the concentration at a specified time
    concentration_at_tt = initial_concentration * (fast_frac * math.exp(-1 * tt / k1) + (1 - fast_frac) * math.exp(-1 * tt / k2))
    return concentration_at_tt


def get_time_efficacy_values(initial_concentration, max_efficacy, fast_frac, k1, k2, hh, nn, total_time):
    # get concentration and efficacy through time
    concentration_through_time = [get_concentration_at_time(tt, initial_concentration, fast_frac, k1, k2) for tt in range(total_time)]
    # efficacy_through_time = [max_efficacy * (1 - math.exp(mm * cc)) for cc in concentration_through_time]
    efficacy_through_time = [max_efficacy / (1 + math.pow((hh / cc), nn)) for cc in concentration_through_time]
    return [[i for i in range(total_time)], efficacy_through_time]


def get_vacc_params_from_pkpd_df(row):
    # extract the parameters describing PKPD from the row of an input file
    initial_concentration = row['initial_concentration']
    max_efficacy = row['max_efficacy']
    fast_frac = row['fast_frac']
    k1 = row['k1']
    k2 = row['k2']
    hh = row['hh']
    nn = row['nn']
    total_time = row['total_time']
    time_efficacy_values = get_time_efficacy_values(initial_concentration, max_efficacy, fast_frac, k1, k2, hh, nn, total_time)
    return time_efficacy_values


def add_triggered_vacc(campaign, vacc_char_df, my_ds=''):
    # set up vaccines to be distributed when triggered by the 'event_add_new_vaccine' event
    if my_ds != '':
        if 'admin_name' in vacc_char_df.columns:
            vacc_char_df = vacc_char_df[vacc_char_df['admin_name'].str.upper() == my_ds.upper()]
    if 'vacc_type' in vacc_char_df.columns:
        row_initial = vacc_char_df[vacc_char_df['vacc_type'] == 'initial'].iloc[0]
        row_boost = vacc_char_df[vacc_char_df['vacc_type'] == 'booster'].iloc[0]
    else:
        row_initial = vacc_char_df.iloc[0]
        row_boost = vacc_char_df.iloc[0]

    """Set vaccine properties (e.g., initial efficacy, waning)"""
    # TODO: currently, assumes that boost has same waning type as initial dose. Should probably change this.
    try:
        waning_type = row_initial['vacc_waning_type']
    except:
        waning_type = 'exponential'

    if waning_type == 'pkpd':
        time_efficacy_values_initial = get_vacc_params_from_pkpd_df(row_initial)
        time_efficacy_values_boost = get_vacc_params_from_pkpd_df(row_boost)
        time_efficacy_initial = time_efficacy_values_initial[1][0]
        time_efficacy_multipliers = [time_efficacy_values_initial[1][yy] / time_efficacy_initial for yy in
                                     range(len(time_efficacy_values_initial[1]))]
        time_efficacy_boost_initial = time_efficacy_values_boost[1][0]
        time_efficacy_boost_multipliers = [time_efficacy_values_boost[1][yy] / time_efficacy_boost_initial for yy in
                                     range(len(time_efficacy_values_boost[1]))]
    else:
        raise ValueError("Unknown vaccine decay type. Only 'pkpd' currently supported.")

    # vaccine is added in response to broadcast event
    # initial vaccine
    add_triggered_vaccine(campaign,
                          start_day=1,
                          vaccine_type="AcquisitionBlocking",
                          trigger_condition_list=['event_add_new_vaccine'],
                          listening_duration=-1,
                          demographic_coverage=1,
                          repetitions=1,
                          timesteps_between_repetitions=-1,
                          ind_property_restrictions=[{'VaccineStatus': 'None'}],
                          vaccine_initial_effect=time_efficacy_initial,
                          vaccine_linear_times=time_efficacy_values_initial[0],
                          vaccine_linear_values=time_efficacy_boost_multipliers,
                          vaccine_expire_at_end=True,
                          disqualifying_properties=[{"vaccine_selected": "Yes"}],
                          efficacy_is_multiplicative=True)


    # booster vaccines
    add_triggered_vaccine(campaign,
                          start_day=1,
                          vaccine_type="AcquisitionBlocking",
                          trigger_condition_list=['event_add_new_vaccine'],
                          listening_duration=-1,
                          demographic_coverage=1,
                          repetitions=1,
                          timesteps_between_repetitions=-1,
                          ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                          vaccine_initial_effect=time_efficacy_boost_initial,
                          vaccine_linear_times=time_efficacy_values_boost[0],
                          vaccine_linear_values=time_efficacy_multipliers,
                          vaccine_expire_at_end=True,
                          disqualifying_properties=[{"vaccine_selected": "Yes"}],
                          efficacy_is_multiplicative=True)


def change_vacc_ips(campaign):
    # update VaccineStatus IP to 'Received_Vaccine' if individual receives a vaccine (note: it looks like the add_triggered_vaccine function may now support doing this internally)
    change_individual_property_triggered(campaign,
                                         new_ip_key='VaccineStatus',
                                         new_ip_value='GotVaccine',
                                         ip_restrictions=[{'VaccineStatus': 'None'}],
                                         triggers=['Received_Vaccine'],
                                         blackout=False)


def add_pkpd_vacc(campaign, vacc_df, my_ds='', cohort_month_shift=0):
    # add delivery of vaccine  on a particular day, which works for seasonal distribution and for age-based ONLY in cohort simulations

    # Sequence of vaccine events:
    #  - on start_day, select people to receive vaccine (change their IP vaccine_selected to True). This will remove old vaccines
    #  - on start_day+1, create a campaign which changes IP vaccine_selected to False (allowing new vaccines to be given) and broadcasts an event that triggers a node-level intervention where the new vaccine is distributed to these same individuals. The vaccine should have Disqualifying_Properties set to {“vaccine_selected”: “True”}. Also change VaccineStatus IP to ReceivedVaccine.
    if my_ds != '':
        if 'admin_name' in vacc_df.columns:
            vacc_df = vacc_df[vacc_df['admin_name'].str.upper() == my_ds.upper()]

    """Note: for cohort-EPI model, we use campaign-style deployment targeted to specific ages (since births disabled in cohort simulation)"""
    for r, row in vacc_df.iterrows():
        # calculate vaccine delivery day, given cohort month shift. If EPI type, don't adjust for cohort month
        #    (because vaccine is given according to individual's age instead of in a mass campaign)
        start_day0 = row['simday']
        if row['deploy_type'] == 'EPI_cohort':
            start_day = start_day0
        elif 'season_cohort' in row['deploy_type']:
            start_day = start_day0 - round(30.4 * cohort_month_shift)
        elif 'season' in row['deploy_type']:
            start_day = start_day0
        else:
            print('WARNING: vaccine delivery name not recognized, age-based vaccination.')
            start_day = start_day0
        if start_day == 0:
            start_day = 1  # avoid issue with vaccines not being given if set to begin on day 0
        if start_day > 0:
            cov_high = row['coverage_high_access']
            cov_low = row['coverage_low_access']
            # Select people to receive vaccine (change their IP vaccine_selected to True)
            """Set group of individuals to receive vaccine and change IPs accordingly"""

            # high-access coverage
            change_individual_property_scheduled(campaign,
                                                 start_day=start_day,
                                                 coverage=cov_high,
                                                 new_ip_key='vaccine_selected',
                                                 new_ip_value='Yes',
                                                 target_age_min=row['agemin'],
                                                 target_age_max=row['agemax'],
                                                 ip_restrictions=[{'SMCAccess': 'High'}])

            # low-access coverage
            change_individual_property_scheduled(campaign,
                                                 start_day=start_day,
                                                 coverage=cov_low,
                                                 new_ip_key='vaccine_selected',
                                                 new_ip_value='Yes',
                                                 target_age_min=row['agemin'],
                                                 target_age_max=row['agemax'],
                                                 ip_restrictions=[{'SMCAccess': 'Low'}])

            # On start_day+1, create a campaign which changes IP vaccine_selected to No (allowing new vaccines to be given) and broadcasts an event that triggers a node-level intervention where the new vaccine is distributed to these same individuals. The vaccine should have Disqualifying_Properties set to {“vaccine_selected”: “Yes”}. Also change VaccineStatus IP to ReceivedVaccine.
            broadcast = BroadcastEvent(campaign, "event_add_new_vaccine")
            change = PropertyValueChanger(campaign, Target_Property_Key="vaccine_selected", Target_Property_Value="No")
            add_campaign_event(campaign=campaign,
                               start_day=start_day+1,
                               demographic_coverage=1,
                               node_ids=None,  # all nodes will get intervention
                               repetitions=1,
                               timesteps_between_repetitions=-1,
                               ind_property_restrictions=[{'vaccine_selected': 'Yes'}],
                               individual_intervention=[change, broadcast])

    change_vacc_ips(campaign)
    return len(vacc_df)


def add_pkpd_epi_vacc(campaign, epi_vacc_df, hfca):
    # add delivery of vaccine a particular number of days after birth, which works for age-based distribution in simulations with vital dynamics

    # Sequence of vaccine events:
    #  - When someone is born, it begins a countdown until the specified touchpoint(s), when an 'epi_touchpoint' event is broadcast (each will have a unique event name).
    #  - Any time an individual has an 'epi_touchpoint' event, it causes them to be eligible for an IP change to 'vaccine_selected:Yes' (whether or not this occurs depends on their access group and the coverage associated with that access group)
    #  - Also have a daily campaign with 100% coverage that applies only to individuals with vaccine_selected:Yes. It should broadcast the 'event_add_new_vaccine' event and change the 'vaccine_selected' IP to 'No'
    #  - The 'event_add_new_vaccine' will trigger a node-level intervention, created in add_triggered_vacc(), where the new vaccine is distributed to these individuals. The vaccine should have Disqualifying_Properties set to {“vaccine_selected”: “True”} so that old vaccines are removed when new ones are given. It also changes VaccineStatus IP to ReceivedVaccine.

    if hfca != '':
        if 'admin_name' in epi_vacc_df.columns:
            epi_vacc_df = epi_vacc_df[epi_vacc_df['admin_name'].str.upper() == hfca.upper()]
    if len(epi_vacc_df) == 0:
        return 0

    # iterate through touchpoints in dataframe, since each one may have a different coverage
    for r, row in epi_vacc_df.iterrows():
        cur_touchpoint = row['epi_touchpoint']
        cur_touchpoint_event = f'epi_touchpoint{cur_touchpoint}'
        cov_high = row['coverage_high_access']
        cov_low = row['coverage_low_access']
        start_epi = row['simday']  # EPI vaccines are only given after this day of the simulation

        # have a campaign to broadcast the current EPI touchpoint event
        broadcast_event = BroadcastEvent(campaign, cur_touchpoint_event)
        add_triggered_campaign_delay_event(campaign, start_day=0,
                                           trigger_condition_list=['Births'],
                                           delay_period_constant=cur_touchpoint,
                                           demographic_coverage=1,
                                           individual_intervention=broadcast_event)

        # Select people to receive vaccine on this EPI touchpoint day (change their IP vaccine_selected to Yes), with coverage depending on access group
        # high-access coverage
        change_individual_property_triggered(campaign,
                                             start_day=start_epi,
                                             coverage=cov_high,
                                             new_ip_key='vaccine_selected',
                                             new_ip_value='Yes',
                                             ip_restrictions=[{'SMCAccess': 'High'}],
                                             triggers=[cur_touchpoint_event],
                                             blackout=False)
        # low-access coverage
        change_individual_property_triggered(campaign,
                                             start_day=start_epi,
                                             coverage=cov_low,
                                             new_ip_key='vaccine_selected',
                                             new_ip_value='Yes',
                                             ip_restrictions=[{'SMCAccess': 'Low'}],
                                             triggers=[cur_touchpoint_event],
                                             blackout=False)

    # Whenever an individual has the 'vaccine_selected:Yes' IP, change the IP vaccine_selected to No (allowing new vaccines to be given) and broadcast an event that triggers a node-level intervention where the new vaccine is distributed to these same individuals. The vaccine should have Disqualifying_Properties set to {“vaccine_selected”: “Yes”}. Also change VaccineStatus IP to ReceivedVaccine.
    broadcast = BroadcastEvent(campaign, "event_add_new_vaccine")
    change = PropertyValueChanger(campaign, Target_Property_Key="vaccine_selected", Target_Property_Value="No")
    add_campaign_event(campaign=campaign,
                       start_day=0,
                       demographic_coverage=1,
                       repetitions=-1,
                       timesteps_between_repetitions=1,
                       ind_property_restrictions=[{'vaccine_selected': 'Yes'}],
                       individual_intervention=[change, broadcast])

    change_vacc_ips(campaign)
    return len(epi_vacc_df)


#####################
# ################### main function to coordinate adding all interventions  ############
#####################


def add_all_interventions(campaign, hfca, seed_index=1, hs_df=pd.DataFrame(), nmf_df=pd.DataFrame(),
                          itn_df=pd.DataFrame(),
                          itn_anc_df=pd.DataFrame(), itn_use_seasonality=pd.DataFrame(),
                          itn_decay_params=pd.DataFrame(),
                          itn_anc_adult_birthday_years=None, itn_epi_df=pd.DataFrame(), itn_cont_df=pd.DataFrame(),
                          itn_chw_df=pd.DataFrame(), itn_chw_annual_df=pd.DataFrame(),
                          irs_df=pd.DataFrame(), lsm_df=pd.DataFrame(), spaq_df=pd.DataFrame(), smc_df=pd.DataFrame(), pmc_df=pd.DataFrame(), vacc_df=pd.DataFrame(),
                          vacc_char_df=pd.DataFrame(), vacc_df_2=pd.DataFrame(), epi_vacc_df=pd.DataFrame(), use_same_access_ips_all_ages=False,
                          sp_resist_day1_multiply=1, adherence_multiplier=1, use_smc_vaccine_proxy=False):
    event_list = []

    if not irs_df.empty:
        has_irs = add_hfca_irs(campaign, irs_df, hfca, seed_index=seed_index)
        if has_irs > 0:
            event_list.append('Received_IRS')
    if not lsm_df.empty:
        has_lsm = add_hfca_lsm(campaign, lsm_df, hfca, seed_index=seed_index)
    if not smc_df.empty:
        has_smc = update_smc_access_ips(campaign, smc_df=smc_df, hfca=hfca, use_same_access_ips_all_ages=use_same_access_ips_all_ages)
        if use_smc_vaccine_proxy:
            has_smc = add_hfca_vaccsmc(campaign, smc_df, hfca,
                                       effective_coverage_resistance_multiplier=sp_resist_day1_multiply,
                                       seed_index=seed_index)
        else:
            has_smc = add_hfca_smc(campaign, smc_df=smc_df, hfca=hfca, adherence_multiplier=adherence_multiplier,
                                   sp_resist_day1_multiply=sp_resist_day1_multiply, seed_index=seed_index)
        if has_smc > 0:
            event_list.append('Received_Campaign_Drugs')
    if not spaq_df.empty:
        has_spaq = update_spaq_access_ips(campaign, spaq_df=spaq_df, hfca=hfca,
                                        use_same_access_ips_all_ages=use_same_access_ips_all_ages)
        has_spaq = add_hfca_spaq(campaign, spaq_df=spaq_df, hfca=hfca, adherence_multiplier=adherence_multiplier,
                               sp_resist_day1_multiply=sp_resist_day1_multiply, seed_index=seed_index)
        if has_spaq > 0:
            event_list.append('Received_Campaign_Drugs')

    if not pmc_df.empty:
        has_pmc = add_ds_vaccpmc(campaign, pmc_df=pmc_df, hfca=hfca)  # per default use vaccpmc
        if has_pmc > 0:
            event_list = event_list + [
                'Received_PMC_VaccDrug']  # 'Received_Vehicle_1','Received_Vehicle_2','Received_Vehicle_3'
    if not vacc_df.empty:
        has_vacc = add_ds_rtss(campaign, rtss_df=vacc_df, hfca=hfca)
        if has_vacc > 0:
            event_list = event_list + ['Received_Vaccine']
    if not vacc_df_2.empty:
        has_pkpd_vacc = update_smc_access_ips(campaign, smc_df=vacc_df_2, hfca=hfca, use_same_access_ips_all_ages=use_same_access_ips_all_ages)
        add_triggered_vacc(campaign, vacc_char_df, hfca)
        has_vacc = add_pkpd_vacc(campaign, vacc_df_2, hfca)
        if has_vacc > 0:
            event_list = event_list + ['Received_Vaccine']
            event_list = event_list + ['event_add_new_vaccine']
    if not epi_vacc_df.empty:
        has_pkpd_vacc = update_smc_access_ips(campaign, smc_df=epi_vacc_df, hfca=hfca,
                                              use_same_access_ips_all_ages=use_same_access_ips_all_ages)
        add_triggered_vacc(campaign, vacc_char_df, hfca)
        has_vacc = add_pkpd_epi_vacc(campaign, epi_vacc_df, hfca)
    if not (itn_df.empty and itn_anc_df.empty and itn_epi_df.empty and itn_cont_df.empty and itn_chw_df.empty and itn_chw_annual_df.empty):
        has_itn = add_hfca_itns(campaign=campaign, itn_df=itn_df, itn_anc_df=itn_anc_df,
                                itn_anc_adult_birthday_years=itn_anc_adult_birthday_years, itn_epi_df=itn_epi_df,
                                itn_cont_df=itn_cont_df, itn_chw_df=itn_chw_df, itn_chw_annual_df=itn_chw_annual_df, hfca=hfca,
                                itn_use_seasonality=itn_use_seasonality, itn_decay_params=itn_decay_params,
                                seed_index=seed_index)
        if has_itn > 0:
            event_list.append('Bednet_Got_New_One')
            event_list.append('Bednet_Using')
            event_list.append('Bednet_Discarded')
    if not hs_df.empty:
        # case management for malaria
        has_cm = add_hfca_hs(campaign, hs_df, hfca, seed_index=seed_index)
        if has_cm:
            event_list.append('Received_Treatment')
            event_list.append('Received_Severe_Treatment')
        # case management for NMFs
        add_nmf_hs(campaign, hs_df, nmf_df, hfca, seed_index=seed_index)
        event_list.append('Received_NMF_Treatment')

    return {"events": event_list}


