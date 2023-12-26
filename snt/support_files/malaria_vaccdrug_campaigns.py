from emodpy_malaria.interventions.drug_campaign import add_drug_campaign
from emodpy_malaria.interventions.vaccine import add_scheduled_vaccine, add_triggered_vaccine
from emodpy_malaria.malaria_config import set_drug_param
from emodpy_malaria.interventions.common import add_triggered_campaign_delay_event
import math
import scipy.stats as stats
import numpy as np


def add_vaccdrug_campaign(campaign, campaign_type: str = 'SMC', start_days: list = None,
                          coverages: list = None, target_group: dict = None,
                          vaccine_param_dict: dict = None, drug_param_dict: dict = None,
                          delay_distribution_dict: dict = None, node_ids: list = None,
                          ind_property_restrictions: list = None,
                          listening_duration: int = -1,
                          target_residents_only: bool = True,
                          trigger_condition_list: list = None,
                          check_eligibility_at_trigger: bool = False,
                          receiving_vaccine_event_name: str = None,
                          receiving_drugs_event_name: str = None,
                          num_iiv_groups: int = 1,
                          receiving_drugs_event: bool = True):
    """
        Add a vaccine + vehicle drug intervention to approximate efficacy of SMC or PMC, as specified in **campaign_type**,
        to the campaign. This intervention uses default parameters corresponding to SMC with SPAQ and PMC with SP, if not
        otherwise specified via vacc_param_dic, drug_param_dic. The vehicle drug instantly clear parasites (blood stage +
        liver stage) and the prophylactic effect is added by the vaccine event.

        Campaign type specifications:
        For SMC, the drug event (MDA drug campaign) is initiated at the specified simdays and generates a broadcast event
        that is used to trigger the vaccine event without delay.
        For PMC, the drug event (MDA drug campaign) is triggered when the individual reaches the eligible age, which is
        broadcast using a birth-triggered triggered_campaign_delay_event, and similar to SMC, the vaccine campaign is
        triggered by the drug event broadcast without delay.
        PMC supports inter-individual-variation (IIV) in initial vaccine efficacy, and delay_distributions around
        specified age to receive the event.

    Args:
        campaign: campaign object to which the intervention will be added, and schema_path container
        campaign_type: Type of drug campaign. Default is **SMC**.
            Available options are:
            * SMC
            * PMC
        start_days: List of start days (integers) when the drug regimen will
            be distributed. Due to diagnostic/treatment configuration,
            the earliest start day is 1. When trigger_condition_list is used
            then the first entry of start_days is the day to start listening
            for the trigger(s).
        coverages: List of effective coverage levels per round of the campaign (fraction of people who receive the
            campaign at each round).

            Examples::

                [0.8, 0.8, 0.8, 0.8] or  [0.8, 0.7, 0.6, 0.5]

        target_group: A dictionary of to specify age range for SMC or age-touchpoints for PMC. Default is Everyone.

            Examples::

                 {'agemin': x, 'agemax': y} , if campaign_type = SMC
                 {'1': 76, '2': 106, '3': 274, ...} , if campaign_type = PMC

        vaccine_param_dict: dictionary of parameters for vaccine to use with this intervention.
            keys expected: 'vaccine_initial_effect', 'vaccine_box_duration', 'vaccine_decay_duration'
            'vaccine_decay_duration' will be converted to vaccine's Decay_Time_Constant by
            vaccine_decay_duration/log(2)

            Examples::

                {'vaccine_initial_effect': 0.598, 'vaccine_box_duration': 21.7, 'vaccine_decay_duration': 1.18}

        drug_param_dict: Dictionary of three parameters for initial parasite clearing efficacy (drug campaign).
            Default: Fitted parameters for SMC-SPAQ and PMC-SP.

            Examples::

                {'drug_irbc_killing': 18.6, 'drug_hep_killing': 1.5, 'drug_box_day': 2.0}

        delay_distribution_dict: Dictionary of lists of distribution parameters. Distributions allowed:
            "GAUSSIAN_DISTRIBUTION","LOG_NORMAL_DISTRIBUTION", "CONSTANT_DISTRIBUTION".
            for age-touchpoints, there is also an internal offset (please see code)
            for "CONSTANT_DISTRIBUTION", only 'distribution_name' is used, age-touchpoint is used as the constant
            for "GAUSSIAN_DISTRIBUTION", 'delay_distribution_mean' is ignored, age-touchpoint is used as mean,
            std is used directly.
            for "LOG_NORMAL_DISTRIBUTION", 'delay_distribution_mean' is converted to mu with
            np.log(tp_time_trigger + mean) - ((1 / 2) * std ** 2)) and std is used directly as sigma

            Examples::

                {'delay_distribution_name': df['distribution_name'],
                 'delay_distribution_mean': df['distribution_mean'],
                 'delay_distribution_std': df['distribution_std']}

        node_ids: The list of nodes to apply this intervention to (**Node_List**
            parameter). If not provided, set value of NodeSetAll.
        ind_property_restrictions:  List of IndividualProperty key:value pairs that
            individuals must have to receive the diagnostic intervention.
            For example, ``[{"IndividualProperty1":"PropertyValue1"},
            {"IndividualProperty2":"PropertyValue2"}]``. Default is no restrictions.
        listening_duration: Length of time, in days, for which the triggered event will be listening for the triggers
        target_residents_only: When set to True the intervention is only distributed to individuals that began the
            simulation in that node.
        trigger_condition_list: List of events that will begin a triggerable
            campaign if campaign_type is SMC. If campaign_type is PMC, campaign is triggered by birth per default.
        check_eligibility_at_trigger: If triggered event is delayed, you have an
            option to check individual/node's eligibility at the initial trigger
            or when the event is actually distributed after delay. (for example, a person might've aged out of the
            intervention before the initial trigger and the intervention distribution)
        receiving_vaccine_event_name:  Event to send out when person received vaccine.
            Default: 'Received_<campaign_type>_VaccDrug'
        receiving_drugs_event_name: Event to send out when person received drugs.
            Event name needs to include 'Received_Vehicle' in it, as otherwise overwritten in drug_campaigns function
            (see drug_campaigns.py L247)
            Default: SMC: 'Received_Vehicle' ; PMC: 'Received_Vehicle_X' with X being number of PMC dose
        num_iiv_groups: Number of individual drug response groups.
            If >1, ind_property_restrictions is set to {'DrugResponseGroup': val} if campaign_type = PMC, not used for
            SMC. Default: SMC: Not used. PMC: 1, IIV only acts on the vaccine event, not the drug event
        receiving_drugs_event: Specify whether to deploy the parasite clearing drug event or the vaccine event only.
            Default: True
            Exception: for PMC, set to False

    Returns:
        Dictionary with drug campaign parameters
    """
    if not receiving_drugs_event_name:
        receiving_drugs_event_name = 'Received_Vehicle'
    if not receiving_vaccine_event_name:
        receiving_vaccine_event_name = f'Received_{campaign_type}_VaccDrug'

    if campaign_type == 'SMC':
        if receiving_drugs_event:
            add_vaccdrug_smc(campaign, start_days=start_days, coverages=coverages,
                             vaccine_param_dict=vaccine_param_dict, drug_param_dict=drug_param_dict,
                             target_group=target_group,
                             receiving_drugs_event=receiving_drugs_event_name,
                             receiving_vaccine_event=receiving_vaccine_event_name,
                             node_ids=node_ids,
                             ind_property_restrictions=ind_property_restrictions,
                             listening_duration=listening_duration,
                             trigger_condition_list=trigger_condition_list,
                             target_residents_only=target_residents_only,
                             check_eligibility_at_trigger=check_eligibility_at_trigger)
        else:
            add_vacc_smc(campaign, start_days=start_days, coverages=coverages,
                         vaccine_param_dict=vaccine_param_dict,
                         target_group=target_group,
                         receiving_vaccine_event=receiving_vaccine_event_name,
                         node_ids=node_ids,
                         ind_property_restrictions=ind_property_restrictions,
                         target_residents_only=target_residents_only)
    elif campaign_type == 'PMC':
        if trigger_condition_list:
            raise ValueError(
                'You passed in a trigger_condition for campaign_type PMC that is per default triggered by "Births". '
                'trigger_condition_list needs to be empty.\n')
        if receiving_drugs_event:
            add_vaccdrug_pmc(campaign, start_days=start_days, coverages=coverages,
                             target_group=target_group, num_iiv_groups=num_iiv_groups,
                             vaccine_param_dict=vaccine_param_dict,
                             drug_param_dict=drug_param_dict,
                             receiving_drugs_event=receiving_drugs_event_name,
                             receiving_vaccine_event=receiving_vaccine_event_name,
                             delay_distribution_dict=delay_distribution_dict,
                             node_ids=node_ids,
                             ind_property_restrictions=ind_property_restrictions,
                             listening_duration=listening_duration,
                             target_residents_only=target_residents_only,
                             check_eligibility_at_trigger=check_eligibility_at_trigger)
        else:
            add_vacc_pmc(campaign, start_days=start_days, coverages=coverages,
                         target_group=target_group, num_iiv_groups=num_iiv_groups,
                         vaccine_param_dict=vaccine_param_dict,
                         receiving_vaccine_event=receiving_vaccine_event_name,
                         delay_distribution_dict=delay_distribution_dict,
                         node_ids=node_ids,
                         ind_property_restrictions=ind_property_restrictions,
                         listening_duration=listening_duration,
                         target_residents_only=target_residents_only)
    else:
        raise ValueError('Invalid campaign_type specified, valid options: "SMC" or "PMC"')

def make_vehicle_drug(config, decay_time, Cmax, C50, irbc_killing, hep_killing=0.0):

    if not (decay_time == 0):
        set_drug_param(config, "Vehicle", "Drug_Decay_T1", decay_time)
        set_drug_param(config, "Vehicle", "Drug_Decay_T2", decay_time)

    if not (Cmax == 0):
        set_drug_param(config, "Vehicle", "Drug_Cmax", Cmax)

    if not (C50 == 0):
        set_drug_param(config, "Vehicle", "Drug_PKPD_C50", C50)

    if not (irbc_killing == 0):
        set_drug_param(config, "Vehicle", "Max_Drug_IRBC_Kill", irbc_killing)

    if not (hep_killing == 0):
        set_drug_param(config, "Vehicle", "Drug_Hepatocyte_Killrate", hep_killing)

    return {'drug_decay_time': decay_time,
            'drug_cmax': Cmax,
            'drug_c50': C50,
            'drug_irbc_killing': irbc_killing,
            'drug_hep_killing': hep_killing}

#def make_vehicle_drug(campaign, drug_box_day: float = 0, drug_irbc_killing: float = 0, drug_hep_killing: float = 0):
#    if drug_box_day:
#        set_drug_param(campaign, "Vehicle", "Drug_Decay_T1", drug_box_day)
#        set_drug_param(campaign, "Vehicle", "Drug_Decay_T2", drug_box_day)
#    if drug_irbc_killing:
#        set_drug_param(campaign, "Vehicle", "Max_Drug_IRBC_Kill", drug_irbc_killing)
#    if drug_hep_killing:
#        set_drug_param(campaign, "Vehicle", "Drug_Hepatocyte_Killrate", drug_hep_killing)
#
#    return {'drug_box_day': drug_box_day,
#            'drug_irbc_killing': drug_irbc_killing,
#            'drug_hep_killing': drug_hep_killing}


def add_vaccdrug_smc(campaign, start_days: list, coverages: list,
                     target_group: dict = None,
                     node_ids: list = None,
                     vaccine_param_dict: dict = None, drug_param_dict: dict = None,
                     receiving_vaccine_event: str = None, receiving_drugs_event: str = None,
                     listening_duration: int = -1, trigger_condition_list: list = None,
                     ind_property_restrictions: dict = None,
                     target_residents_only: int = 1,
                     check_eligibility_at_trigger: bool = False):
    """
        Add a vaccine + vehicle drug intervention to approximate efficacy of SMC. This intervention uses default
        parameters corresponding to SMC with SPAQ, if not otherwise specified via vaccine_param_dict and drug_param_dict.
        The vehicle drug instantly clear parasites (blood stage + liver stage) and the prophylactic effect is added
        by the vaccine event.

        Campaign type specifications:
        For SMC, the drug event (MDA drug campaign) is initiated at the specified simdays and generates a broadcast event
        that is used to trigger the vaccine event without delay.

    Args:
        campaign: campaign object to which the intervention will be added, and schema_path container
        start_days: list of days on which to run the interventions
        coverages: list of coverages for each day in start_days
        target_group: A dictionary of to specify age range for SMC
            Default is Everyone.

             Example::

                 {'agemin': x, 'agemax': y} for campaign_type = SMC

        node_ids: The list of nodes to apply this intervention to (**Node_List**
            parameter). If not provided, set value of NodeSetAll.
        vaccine_param_dict: dictionary of parameters for vaccine to use with this intervention.
            keys expected: 'vaccine_initial_effect', 'vaccine_box_duration', 'vaccine_decay_duration'
            'vaccine_decay_duration' will be converted to vaccine's Decay_Time_Constant by
            vaccine_decay_duration/log(2)
            Example and Default: {'vaccine_initial_effect': 0.598, 'vaccine_box_duration': 21.7,
            'vaccine_decay_duration': 1.18}
        drug_param_dict: dictionary of parameters for a drug to use with this intervention, these will be assigned
            to the 'Vehicle' drug. Example and Default: {'drug_box_day': 2.0, 'drug_irbc_killing': 10.8,
            'drug_hep_killing': 3.64}
        receiving_vaccine_event:  Event to send out when person received vaccine.
            Default: 'Received_<campaign_type>_VaccDrug'
        receiving_drugs_event: Event to send out when person received drugs.
            Event name needs to include 'Received_Vehicle' in it, as otherwise overwritten in drug_campaigns function
            (see drug_campaigns.py L247)
            Default: SMC: 'Received_Vehicle'
        listening_duration: Length of time, in days, for which the triggered event will be listening for the triggers
        trigger_condition_list: List of events that will begin a triggerable
            campaign if campaign_type is SMC. If campaign_type is PMC, campaign is triggered by birth per default.
        ind_property_restrictions:  List of IndividualProperty key:value pairs that
            individuals must have to receive the diagnostic intervention.
            For example, ``[{"IndividualProperty1":"PropertyValue1"},
            {"IndividualProperty2":"PropertyValue2"}]``. Default is no restrictions.
        target_residents_only: When set to True the intervention is only distributed to individuals that began the
            simulation in that node.
        check_eligibility_at_trigger: If triggered event is delayed, you have an
            option to check individual/node's eligibility at the initial trigger
            or when the event is actually distributed after delay. (for example, a person might've aged out of the
            intervention before the initial trigger and the intervention distribution)
        receiving_drugs_event_name: Event to send out when person received drugs.
            Event name needs to include 'Received_Vehicle' in it, as otherwise overwritten in drug_campaigns function
            (see drug_campaigns.py L247)
            Default: SMC: 'Received_Vehicle'; PMC: 'Received_Vehicle_X' with X being number of PMC dose
        num_iiv_groups: Number of individual drug response groups.
            If >1, ind_property_restrictions is set to {'DrugResponseGroup': val} if campaign_type = PMC, not used for
            SMC. Default: SMC: Not used; PMC: 1, IIV only acts on the vaccine event, not the drug event
        receiving_drugs_event: Specify whether to deploy the parasite clearing drug event or the vaccine event only.
            Default: True
            Exception: for PMC, set to False

    Returns:
        dictionary of tags
    """
    if len(start_days) != len(coverages):
        raise ValueError(f"Length of start_days - {len(start_days)}, should be equal to length of coverages - "
                         f"{len(coverages)}, but it's not.\n")
    target_age_min = 0
    target_age_max = 125
    if target_group:
        target_age_min = target_group['agemin']
        target_age_max = target_group['agemax']
    if not vaccine_param_dict:
        vaccine_param_dict = {'vacc_initial_effect': 0.598, 'vacc_box_duration': 21.7, 'vacc_decay_duration': 1.18}
    if not drug_param_dict:
        # drug_param_dic = {'drug_irbc_killing': 18.6, 'drug_hep_killing': 1.5, 'drug_box_day': 2.0}
        drug_param_dict = {'drug_box_day': 2.0, 'drug_irbc_killing': 10.8, 'drug_hep_killing': 3.64}

    vaccine_initial_effect = vaccine_param_dict['vacc_initial_effect']
    vaccine_box_duration = vaccine_param_dict['vacc_box_duration']
    vaccine_decay_duration = vaccine_param_dict['vacc_decay_duration']

    #drug_decay_time = drug_param_dict['drug_decay_time']
    #drug_Cmax = drug_param_dict['drug_Cmax']
    #drug_C50 = drug_param_dict['drug_C50']
    #drug_irbc_killing = drug_param_dict['drug_irbc_killing']
    #drug_hep_killing = drug_param_dict['drug_hep_killing']
    #make_vehicle_drug(campaign, 
    #                  decay_time=drug_decay_time,
    #                  Cmax=drug_Cmax,
    #                  C50=drug_C50,
    #                  maxeff=drug_irbc_killing,
    #                  hep=drug_hep_killing)

    for (d, cov) in zip(start_days, coverages):
        add_drug_campaign(campaign, campaign_type='MDA',
                          drug_code='Vehicle',
                          start_days=[d + 1],  # when triggered and triggering interventions are deployed on the
                          # same day, that-day events are not "heard" by the triggered intervention, so I'm
                          # off-setting them by one day
                          coverage=cov,
                          repetitions=-1,
                          tsteps_btwn_repetitions=-1,
                          listening_duration=listening_duration,
                          target_group=target_group,
                          ind_property_restrictions=ind_property_restrictions,
                          receiving_drugs_event_name=receiving_drugs_event,
                          trigger_condition_list=trigger_condition_list,
                          target_residents_only=target_residents_only,
                          node_ids=node_ids,
                          check_eligibility_at_trigger=check_eligibility_at_trigger)

    add_triggered_vaccine(campaign,
                          start_day=start_days[0],  # otherwise it won't "hear" the first round of drugs
                          demographic_coverage=1,
                          trigger_condition_list=[receiving_drugs_event],
                          listening_duration=listening_duration,
                          target_age_min=target_age_min,
                          target_age_max=target_age_max,
                          node_ids=node_ids,
                          ind_property_restrictions=ind_property_restrictions,
                          intervention_name="RTSS",
                          vaccine_type="AcquisitionBlocking",
                          vaccine_initial_effect=vaccine_initial_effect,
                          vaccine_box_duration=vaccine_box_duration,
                          vaccine_decay_time_constant=vaccine_decay_duration / math.log(2),
                          efficacy_is_multiplicative=True,
                          broadcast_event=receiving_vaccine_event)

    return {'smc_cov': sum(coverages) / len(coverages),
            'total_smc_rounds': len(coverages)}


def add_vacc_smc(campaign, start_days, coverages, target_group: dict = None,
                 vaccine_param_dict: dict = None, receiving_vaccine_event: str = None,
                 node_ids: list = None,
                 ind_property_restrictions: dict = None,
                 target_residents_only: bool = True):
    """

    Args:
        campaign:
        start_days:
        coverages:
        target_group:
        vaccine_param_dict:
        receiving_vaccine_event:
        node_ids:
        ind_property_restrictions:
        target_residents_only:

    Returns:

    """
    if not vaccine_param_dict:
        vaccine_param_dict = {'vacc_initial_effect': 0.98, 'vacc_box_duration': 30,
                              'vacc_decay_duration': 28}
        # vaccine_param_dict = {'vacc_initial_effect': 0.85, 'vacc_box_duration': 13,  # !30,
        #                       'vacc_decay_duration': 20}
    vacc_initial_effect = vaccine_param_dict['vacc_initial_effect']
    vacc_box_duration = vaccine_param_dict['vacc_box_duration']
    vacc_decay_duration = vaccine_param_dict['vacc_decay_duration']
    if len(start_days) != len(coverages):
        raise ValueError(f"Length of start_days - {len(start_days)}, should be equal to length of coverages - "
                         f"{len(coverages)}, but it's not.\n")
    target_age_min = 0
    target_age_max = 125
    if target_group:
        target_age_min = target_group['agemin']
        target_age_max = target_group['agemax']
    vacc_smc_offset = 17  # - x days since no drug clearing event, hence delayed efficacy

    for (d, cov) in zip(start_days, coverages):
        add_scheduled_vaccine(campaign,
                              start_day=d - vacc_smc_offset,
                              demographic_coverage=cov,
                              target_age_min=target_age_min,
                              target_age_max=target_age_max,
                              target_residents_only=target_residents_only,
                              node_ids=node_ids,
                              ind_property_restrictions=ind_property_restrictions,
                              intervention_name="RTSS",
                              vaccine_type="AcquisitionBlocking",
                              vaccine_initial_effect=vacc_initial_effect,
                              vaccine_box_duration=vacc_box_duration,
                              vaccine_decay_time_constant=vacc_decay_duration / math.log(2),
                              efficacy_is_multiplicative=True,
                              broadcast_event=receiving_vaccine_event)

    return {'smc_cov': sum(coverages) / len(coverages),
            'total_smc_rounds': len(coverages)}


def add_vaccdrug_pmc(campaign, start_days: list, coverages: list,
                     target_group: dict = None,
                     num_iiv_groups: int = 1,
                     vaccine_param_dict: dict = None, drug_param_dict: dict = None,
                     receiving_vaccine_event: str = None, receiving_drugs_event: str = None,
                     listening_duration: int = -1, node_ids: list = None,
                     delay_distribution_dict: dict = None,
                     ind_property_restrictions: dict = None, target_residents_only: int = 1,
                     check_eligibility_at_trigger: bool = False):
    """

    Args:
        campaign:
        start_days:
        coverages:
        target_group:
        num_iiv_groups:
        vaccine_param_dict:
        drug_param_dict:
        receiving_vaccine_event:
        receiving_drugs_event:
        listening_duration:
        node_ids:
        delay_distribution_dict:
        ind_property_restrictions:
        target_residents_only:
        check_eligibility_at_trigger:

    Returns:

    """

    if not vaccine_param_dict:
        vaccine_param_dict = {'vacc_initial_effect': 0.85, 'vacc_box_duration': 13.20, 'vacc_decay_duration': 11.53}
    if not drug_param_dict:
        drug_param_dict = {'drug_irbc_killing': 16.03, 'drug_hep_killing': 0.92, 'drug_box_day': 2.0}

    vaccine_initial_effect = vaccine_param_dict['vacc_initial_effect']
    vaccine_box_duration = vaccine_param_dict['vacc_box_duration']
    vaccine_decay_duration = vaccine_param_dict['vacc_decay_duration']

    drug_box_day = drug_param_dict['drug_box_day']
    drug_irbc_killing = drug_param_dict['drug_irbc_killing']
    drug_hep_killing = drug_param_dict['drug_hep_killing']
    make_vehicle_drug(campaign, drug_box_day=drug_box_day, drug_irbc_killing=drug_irbc_killing,
                      drug_hep_killing=drug_hep_killing)

    pmc_touchpoints = list(target_group.values())
    pmc_event_names = [f'PMC_{x + 1}' for x in range(len(pmc_touchpoints))]
    if len(pmc_touchpoints) != len(coverages) or len(pmc_touchpoints) != len(
            delay_distribution_dict['delay_distribution_name']):
        raise ValueError(f"Length of target_groups's values - {len(pmc_touchpoints)}, should be equal to "
                         f"length of coverages - "
                         f"{len(coverages)} and length of the list of delay distributions "
                         f"{len(delay_distribution_dict['delay_distribution_name'])} but it's not.\n")

    for i, (tp_time_trigger, cov, event_name) in enumerate(zip(pmc_touchpoints, coverages, pmc_event_names)):

        delay_distribution_name = list(delay_distribution_dict['delay_distribution_name'])[i]
        mean = list(delay_distribution_dict['delay_distribution_mean'])[i]
        std = list(delay_distribution_dict['delay_distribution_std'])[i]

        if delay_distribution_name == "LOG_NORMAL_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                  "Delay_Period_Log_Normal_Mu": (
                                          np.log(tp_time_trigger + mean) - ((1 / 2) * std ** 2)),
                                  "Delay_Period_Log_Normal_Sigma": std}
        elif delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                  "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                  "Delay_Period_Gaussian_Std_Dev": std}
        else:  # Assume we want CONSTANT_DISTRIBUTION for our delay with tp_time_trigger being the delay
            delay_distribution = {"Delay_Period_Distribution": "CONSTANT_DISTRIBUTION",
                                  "Delay_Period_Constant": tp_time_trigger}

        """When running with IIV birth triggered event needs to happen separate from drug campaigns"""
        from emod_api.interventions.common import BroadcastEvent, DelayedIntervention
        # triggered_campaign_delay_event only has option for constant delay, but we need different
        # distributions, so we're manually creating a delayed intervention that broadcasts an event
        # and slipping it into the triggered intervention
        broadcast_event = BroadcastEvent(campaign, event_name)
        delayed_intervention = DelayedIntervention(campaign, Configs=[broadcast_event], Delay_Dict=delay_distribution)
        add_triggered_campaign_delay_event(campaign, start_day=start_days[0],
                                           trigger_condition_list=['Births'],
                                           demographic_coverage=1,
                                           listening_duration=listening_duration,
                                           individual_intervention=delayed_intervention)

        add_drug_campaign(campaign, campaign_type='MDA',
                          drug_code='Vehicle',
                          start_days=[start_days[0]],
                          coverage=cov,
                          repetitions=-1,
                          tsteps_btwn_repetitions=-1,
                          listening_duration=listening_duration,
                          trigger_condition_list=[event_name],
                          ind_property_restrictions=ind_property_restrictions,
                          receiving_drugs_event_name=f'{receiving_drugs_event}_{i + 1}',
                          node_ids=node_ids,
                          target_residents_only=target_residents_only,
                          check_eligibility_at_trigger=check_eligibility_at_trigger
                          )

        if num_iiv_groups > 1:
            iiv_groups = ["Group%d" % x for x in range(num_iiv_groups)]
            for index, val in enumerate(iiv_groups):
                eff_lower = vaccine_initial_effect - (
                        vaccine_initial_effect * 0.25)  # FIXME, arbitrary default parameters
                eff_upper = vaccine_initial_effect + (vaccine_initial_effect * 0.25)
                if eff_lower <= 0:
                    eff_lower = vaccine_initial_effect
                if eff_upper > 1:
                    eff_upper = 0.98
                eff_sd = 0.025
                vaccine_initial_effect = stats.truncnorm.rvs((eff_lower - vaccine_initial_effect) / eff_sd,
                                                             (eff_upper - vaccine_initial_effect) / eff_sd,
                                                             loc=vaccine_initial_effect,
                                                             scale=eff_sd)

                add_triggered_vaccine(campaign,
                                      start_day=start_days[0],
                                      trigger_condition_list=[f'{receiving_drugs_event}_{i + 1}'],
                                      demographic_coverage=1,
                                      ind_property_restrictions=[{'DrugResponseGroup': val}],
                                      broadcast_event=receiving_vaccine_event,
                                      intervention_name="RTSS",
                                      vaccine_type="AcquisitionBlocking",
                                      vaccine_initial_effect=vaccine_initial_effect,
                                      vaccine_box_duration=vaccine_box_duration,
                                      vaccine_decay_time_constant=vaccine_decay_duration / math.log(2),
                                      efficacy_is_multiplicative=True
                                      )
        else:
            add_triggered_vaccine(campaign,
                                  start_day=start_days[0],
                                  trigger_condition_list=[f'{receiving_drugs_event}_{i + 1}'],
                                  demographic_coverage=1,
                                  ind_property_restrictions=ind_property_restrictions,
                                  broadcast_event=receiving_vaccine_event,
                                  intervention_name="RTSS",
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=vaccine_initial_effect,
                                  vaccine_box_duration=vaccine_box_duration,
                                  vaccine_decay_time_constant=vaccine_decay_duration / math.log(2),
                                  efficacy_is_multiplicative=True
                                  )

    return {'pmc_cov': sum(coverages) / len(coverages),
            'pmc_coverages': '-'.join([str(round(x, 2)) for x in coverages]),
            'n_pmc_touchpoints': len(pmc_touchpoints),
            'pmc_touchpoints': '-'.join([str(x) for x in pmc_touchpoints]),
            }


def add_vacc_pmc(campaign, start_days: list, coverages: list, target_group: dict = None,
                 num_iiv_groups: int = 1,
                 vaccine_param_dict: dict = None,
                 receiving_vaccine_event: str = None,
                 listening_duration: int = -1,
                 node_ids: list = None,
                 delay_distribution_dict: dict = None,
                 ind_property_restrictions: dict = None, target_residents_only: bool = True):
    if not vaccine_param_dict:
        vaccine_param_dict = {'vacc_initial_effect': 0.80, 'vacc_box_duration': 32, 'vacc_decay_duration': 10}

    vaccine_initial_effect = vaccine_param_dict['vacc_initial_effect']
    vaccine_box_duration = vaccine_param_dict['vacc_box_duration']
    vaccine_decay_duration = vaccine_param_dict['vacc_decay_duration']

    pmc_touchpoints = list(target_group.values())
    pmc_event_names = [f'PMC_{x + 1}' for x in range(len(pmc_touchpoints))]
    vacc_pmc_offset = 17  # - x days since no drug clearing event, hence delayed efficacy
    if len(pmc_touchpoints) != len(coverages):
        raise ValueError(f"Length of target_groups's values - {len(pmc_touchpoints)}, should be equal to "
                         f"length of coverages - "
                         f"{len(coverages)}, but it's not.\n")

    for i, (tp_time_trigger, cov, event_name) in enumerate(zip(pmc_touchpoints, coverages, pmc_event_names)):
        tp_time_trigger = tp_time_trigger - vacc_pmc_offset
        delay_distribution_name = list(delay_distribution_dict['delay_distribution_name'])[i]
        mean = list(delay_distribution_dict['delay_distribution_mean'])[i]
        std = list(delay_distribution_dict['delay_distribution_std'])[i]

        if delay_distribution_name == "LOG_NORMAL_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                  "Delay_Period_Log_Normal_Mu": (
                                          np.log(tp_time_trigger + mean) - ((1 / 2) * std ** 2)),
                                  "Delay_Period_Log_Normal_Sigma": std}
        elif delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                  "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                  "Delay_Period_Gaussian_Std_Dev": std}
        elif delay_distribution_name == "CONSTANT_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "CONSTANT_DISTRIBUTION",
                                  "Delay_Period_Constant": tp_time_trigger}

        """When running with IIV birth triggered event needs to happen separate from drug campaigns"""
        from emod_api.interventions.common import BroadcastEvent, DelayedIntervention
        # triggered_campaign_delay_event only has option for constant delay, but we need different
        # distributions, so we're manually creating a delayed intervention that broadcasts an event
        # and slipping it into the triggered intervention
        broadcast_event = BroadcastEvent(campaign, event_name)
        delayed_intervention = DelayedIntervention(campaign, Configs=[broadcast_event], Delay_Dict=delay_distribution)
        add_triggered_campaign_delay_event(campaign, start_day=start_days[0],
                                           trigger_condition_list=['Births'],
                                           demographic_coverage=1,
                                           listening_duration=listening_duration,
                                           node_ids=node_ids,
                                           individual_intervention=delayed_intervention)

        if num_iiv_groups > 1:
            iiv_groups = ["Group%d" % x for x in range(num_iiv_groups)]
            for index, val in enumerate(iiv_groups):
                eff_lower = vaccine_initial_effect - (
                        vaccine_initial_effect * 0.25)  # FIXME, arbitrary default parameters
                eff_upper = vaccine_initial_effect + (vaccine_initial_effect * 0.25)
                if eff_lower <= 0:
                    eff_lower = vaccine_initial_effect
                if eff_upper > 1:
                    eff_upper = 0.98
                eff_sd = 0.025
                vaccine_initial_effect = stats.truncnorm.rvs((eff_lower - vaccine_initial_effect) / eff_sd,
                                                             (eff_upper - vaccine_initial_effect) / eff_sd,
                                                             loc=vaccine_initial_effect,
                                                             scale=eff_sd)
                add_triggered_vaccine(campaign,
                                      start_day=start_days[0],
                                      trigger_condition_list=[event_name],
                                      demographic_coverage=cov,
                                      ind_property_restrictions=[{'DrugResponseGroup': val}],
                                      node_ids=node_ids,
                                      broadcast_event=receiving_vaccine_event,
                                      intervention_name="RTSS",
                                      vaccine_type="AcquisitionBlocking",
                                      vaccine_initial_effect=vaccine_initial_effect,
                                      vaccine_box_duration=vaccine_box_duration,
                                      vaccine_decay_time_constant=vaccine_decay_duration / math.log(2),
                                      efficacy_is_multiplicative=True
                                      )

        else:
            add_triggered_vaccine(campaign,
                                  start_day=start_days[0],
                                  trigger_condition_list=[event_name],
                                  demographic_coverage=cov,
                                  ind_property_restrictions=ind_property_restrictions,
                                  node_ids=node_ids,
                                  target_residents_only=target_residents_only,
                                  broadcast_event=receiving_vaccine_event,
                                  intervention_name="RTSS",
                                  vaccine_type="AcquisitionBlocking",
                                  vaccine_initial_effect=vaccine_initial_effect,
                                  vaccine_box_duration=vaccine_box_duration,
                                  vaccine_decay_time_constant=vaccine_decay_duration / math.log(2),
                                  efficacy_is_multiplicative=True
                                  )

    return {'pmc_cov': sum(coverages) / len(coverages),
            'pmc_coverages': '-'.join([str(round(x, 2)) for x in coverages]),
            'n_pmc_touchpoints': len(pmc_touchpoints),
            'pmc_touchpoints': '-'.join([str(x) for x in pmc_touchpoints]),
            }
