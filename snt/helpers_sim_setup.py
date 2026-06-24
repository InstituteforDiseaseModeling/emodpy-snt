import os
import warnings
import pandas as pd
import numpy as np
import emodpy_malaria.malaria_config as malaria_config
from emodpy_malaria.malaria_config import configure_linear_spline, set_species_param, add_species
from emod_api.interventions.common import change_individual_property_scheduled, change_individual_property_at_age
import emod_api.config.default_from_schema_no_validation as dfs


def update_basic_params(config, manifest, project_path):
    vector_species = ['arabiensis', 'funestus', 'gambiae']
    add_species(config, manifest, vector_species)
    config.parameters.Birth_Rate_Dependence = "FIXED_BIRTH_RATE"
    config.parameters.Age_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'
    config.parameters['logLevel_JsonConfigurable'] = 'ERROR'
    config.parameters['logLevel_VectorHabitat'] = 'ERROR'
    config.parameters['logLevel_StandardEventCoordinator'] = 'ERROR'
    config.parameters['logLevel_InterventionsContainer'] = 'ERROR'
    config.parameters['Habitat_Multiplier'] = 1
    config.parameters.Enable_Default_Reporting = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Property_Output = 0
    config.parameters.Enable_Demographics_Risk = 1
    config.parameters.Enable_Demographics_Birth = 1
    config.parameters.Report_Detection_Threshold_Blood_Smear_Parasites = 0
    config.parameters.Report_Parasite_Smear_Sensitivity = 0.01  # number of microliters of blood examined
    config.parameters.Incubation_Period_Distribution = 'CONSTANT_DISTRIBUTION'
    config.parameters.Incubation_Period_Constant = 7  # parameter shortened from 7 to yield a 12-13 day incubation period

    # read in and set vector bionomics
    vector_bionomics = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'vector_bionomics.csv'))
    for vector in vector_species:
        bionomics_row = vector_bionomics[vector_bionomics['species'] == vector].reset_index()
        # overwrite=True is only needed for parameters that are lists (as you can either overwrite or add to them)
        set_species_param(config, vector, 'Anthropophily', bionomics_row['Anthropophily'][0])
        set_species_param(config, vector, 'Indoor_Feeding_Fraction', bionomics_row['Indoor_Feeding_Fraction'][0])
        set_species_param(config, vector, 'Days_Between_Feeds', bionomics_row['Days_Between_Feeds'][0])

    # added to match the old config to the new config
    config.parameters.Enable_Natural_Mortality = 1
    config.parameters.Enable_Initial_Prevalence = 1
    config.parameters.Base_Air_Temperature = 22
    # config.parameters.Enable_Vector_Migration = 0


def habitat_scales(project_path):
    # done
    rel_abundance_fname = os.path.join(project_path, 'ento', 'DS_vector_rel_abundance.csv')
    rdf = pd.read_csv(rel_abundance_fname)
    rdf = rdf.rename(columns={'Anopheles_arabiensis': 'arabiensis_scale_factor',
                              'Anopheles_coluzzii_gambiae': 'gambiae_scale_factor',
                              'Anopheles_funestus_subgroup': 'funestus_scale_factor'})
    rdf = rdf.set_index('DS')
    return rdf


def set_input_files(config, hfca, archetype_hfca, population_size=1000):
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Base_Air_Temperature = 25
    config.parameters.Base_Land_Temperature = 25
    config.parameters["DS"] = hfca
    config.parameters["Archetype"] = archetype_hfca
    return {'DS': hfca}


def get_burnin_exp(platform, burnin_id='', burnin_fname=''):
    if burnin_id:
        print("building from pickup: identifying burnin addresses from %s" % burnin_id)
        ser_df = platform.create_sim_directory_df(burnin_id)
    else:
        print("building from pickup: identifying burnin addresses from %s" % burnin_fname)
        ser_df = pd.read_csv(burnin_fname)
    return ser_df


def set_up_hfca(config, manifest, hfca, archetype_hfca=None,
                pull_from_serialization=False, ser_date=50 * 365,
                hdf=None, lhdf=None, population_size=1000,
                hab_multiplier=-1, run_number=-1, use_arch_burnin=False, ser_df=pd.DataFrame()):
    set_input_files(config, hfca, archetype_hfca, population_size)
    if not archetype_hfca:
        archetype_hfca = hfca

    set_habitats(config, manifest, hfca, hdf, lhdf, archetype_hfca, abs(hab_multiplier))

    if pull_from_serialization:
        hab_scale_factor_param_name = 'Habitat_Multiplier'
        if 'Run_Number' in ser_df.columns.values:
            ser_df['Run_Number'] = ser_df['Run_Number'].astype(int)
        if use_arch_burnin:
            if hab_multiplier >= 0 and run_number >= 0:
                ser_df[hab_scale_factor_param_name] = ser_df[hab_scale_factor_param_name].apply(
                    lambda x: np.round(float(x), 5))
                sdf = ser_df[(ser_df[hab_scale_factor_param_name] >= (np.round(hab_multiplier, 5) - 0.00001)) &
                             (ser_df[hab_scale_factor_param_name] <= (np.round(hab_multiplier, 5) + 0.00001)) &
                             (ser_df['Run_Number'] == run_number) &
                             (ser_df['admin_name'] == archetype_hfca)]
                config.parameters[hab_scale_factor_param_name] = hab_multiplier
                ser_path = sdf['outpath'].values[0]
            else:
                ser_path = ser_df['outpath'].values[0]
        elif 'admin_name' in ser_df.columns.values:
            ser_df = ser_df[ser_df['Run_Number'] == run_number]
            ser_df = ser_df.set_index('admin_name')
            ser_path = ser_df.at[hfca, 'outpath']
        else:
            ser_path = ser_df['outpath'].values[0]

        config.parameters.Serialized_Population_Reading_Type = 'READ'
        config.parameters.Serialized_Population_Path = os.path.join(ser_path, 'output')
        config.parameters.Serialized_Population_Filenames = ['state-%05d.dtk' % ser_date]
        config.parameters.Enable_Random_Generator_From_Serialized_Population = 0
        config.parameters.Serialization_Mask_Node_Read = 0
        # 0 corresponds to the previous version default: the same larval habitat parameters will be used
        # in the burnin and pickup (from the burnin config)

    return {'admin_name': hfca}


def set_habitats(config, manifest, hfca, hdf, lhdf, archetype_hfca, hab_multiplier):
    # done
    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0.00001, a), max(0.00001, f), max(0.00001, g))

    my_spline, maxvalue, const = load_spline_and_scale_factors(lhdf, archetype_hfca)
    const_mult = 1 if hab_multiplier >= 1 else hab_multiplier

    for (s, sp) in zip(fraction, ['arabiensis', 'funestus', 'gambiae']):
        linear_spline_habitat = configure_linear_spline(manifest,
                                                        max_larval_capacity=pow(10, maxvalue) * s * hab_multiplier,
                                                        capacity_distribution_number_of_years=1,
                                                        capacity_distribution_over_time={
                                                            "Times": [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
                                                                      334],
                                                            "Values": my_spline
                                                        }
                                                        )
        set_species_param(config, sp, "Habitats", linear_spline_habitat, overwrite=True)
        habitat = dfs.schema_to_config_subnode(manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"])
        habitat.parameters.Habitat_Type = "CONSTANT"
        habitat.parameters.Max_Larval_Capacity = pow(10, const) * s * const_mult
        set_species_param(config, sp, "Habitats", habitat.parameters, overwrite=False)


def load_spline_and_scale_factors(lhdf, archetype_hfca):
    lhdf = lhdf.set_index('archetype')
    my_spline = [lhdf.at[archetype_hfca, 'MonthVal%d' % x] for x in range(1, 13)]
    maxvalue = lhdf.at[archetype_hfca, 'MaxHab']
    const = lhdf.at[archetype_hfca, 'Constant']
    # pop_scale = lhdf.at[archetype_hfca, 'pop_scale']

    return my_spline, maxvalue, const  # , pop_scale


def load_master_csv(project_path):
    master_csv = os.path.join(project_path, 'admin_pop_archetype.csv')
    df = pd.read_csv(master_csv, encoding='latin')
    df['admin_name'] = df['admin_name'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
    df = df.set_index('admin_name')
    return df


def update_smc_access_ips(campaign, hfca, smc_df, use_same_access_ips_all_ages=False):
    df = smc_df[smc_df['admin_name'] == hfca]

    if use_same_access_ips_all_ages:
        # Set IPs at the beginning of the simulation and then update individuals' IPs at their birthday.
        # However, individuals born before the start of the simulation are missed by the birthday-triggered IP change
        #  (it is birth-triggered with a delay, so it doesn't get applied to anyone who was already born when the simulation begins).
        # I believe this approach only works if we use the same access IPs for individuals of all ages (therefore, it does not support different high-access fractions for U5 versus 5-10
        if 'high_access' in df.columns:
            high_access_fraction = df['high_access'].values[0]
        else:
            high_access_fraction = df['high_access_U5'].values[0]

        # set IPs at beginning of simulation (does not use IPs from burnin since SMC coverage may change)
        change_individual_property_scheduled(campaign, coverage=1,
                                             new_ip_key='SMCAccess', new_ip_value="Low")
        change_individual_property_scheduled(campaign, coverage=df['high_access_U5'].values[0],
                                             new_ip_key='SMCAccess', new_ip_value="High")

        # set how IPs change as individuals are born or age
        change_individual_property_at_age(campaign, 'SMCAccess', 'Low', 1, coverage=1)
        change_individual_property_at_age(campaign, 'SMCAccess', 'High', 2, coverage=high_access_fraction)

    else:  # change the SMCAccess each year, but this does mean scrambling who is in low versus high access group between years
        # change SMCAccess property of newborns (important for those born after the first SMC round in a year)
        if df.shape[0] > 0:
            change_individual_property_scheduled(campaign, coverage=1,
                                                 new_ip_key='SMCAccess', new_ip_value="Low",
                                                 target_age_min=0, target_age_max=5)
            change_individual_property_scheduled(campaign, coverage=df['high_access_U5'].values[0],
                                                 new_ip_key='SMCAccess', new_ip_value="High",
                                                 target_age_min=0, target_age_max=5)

            # the week before the first SMC round in each year, change the SMCAccess IP for the U5 and O5 age groups
            first_round_df = df[df["round"] == 1]
            for idx, row in first_round_df.iterrows():
                change_individual_property_scheduled(campaign, start_day=(row['simday'] - 7), coverage=1,
                                                     new_ip_key='SMCAccess', new_ip_value="Low",
                                                     target_age_min=0, target_age_max=5)
                change_individual_property_scheduled(campaign, start_day=(row['simday'] - 7), coverage=row['high_access_U5'],
                                                     new_ip_key='SMCAccess', new_ip_value="High",
                                                     target_age_min=0, target_age_max=5)
                change_individual_property_scheduled(campaign, start_day=(row['simday'] - 7), coverage=1,
                                                     new_ip_key='SMCAccess', new_ip_value="Low",
                                                     target_age_min=5, target_age_max=120)
                change_individual_property_scheduled(campaign, start_day=(row['simday'] - 7), coverage=row['high_access_5_10'],
                                                     new_ip_key='SMCAccess', new_ip_value="High",
                                                     target_age_min=5, target_age_max=120)

    return {'admin_name': hfca}




def update_spaq_access_ips(campaign, hfca, spaq_df, use_same_access_ips_all_ages=False):
    # change the SMCAccess each year, but this does mean scrambling who is in low versus high access group between years
    df = spaq_df[spaq_df['admin_name'] == hfca]

    # ---- derive access age-groups for this admin from the flexible df ----
    # One access group per unique (min_age, max_age); high-access fraction can differ by age group
    age_groups = (df[['min_age', 'max_age', 'high_access_frac']]
                  .drop_duplicates()
                  .sort_values('min_age')
                  .reset_index(drop=True))

    # guard: within each year, age bands must be disjoint (touching endpoints OK).
    # Overlap ACROSS years is fine — the IP is re-assigned each year.
    for yr, yr_df in df.groupby('year'):
        bands = (yr_df[['min_age', 'max_age']].drop_duplicates()
                 .sort_values('min_age').reset_index(drop=True))
        overlaps = [(tuple(bands.loc[i, ['min_age', 'max_age']]),
                     tuple(bands.loc[j, ['min_age', 'max_age']]))
                    for i in range(len(bands)) for j in range(i + 1, len(bands))
                    if (bands.loc[i, 'min_age'] < bands.loc[j, 'max_age']) and
                    (bands.loc[j, 'min_age'] < bands.loc[i, 'max_age'])]
        if overlaps:
            msg = ("update_spaq_access_ips: overlapping age groups within year %s for "
                   "admin '%s': %s. Age groups must be disjoint within a year so each "
                   "person is in exactly one access group." % (yr, hfca, overlaps))
            warnings.warn('\n' + '!' * 70 + '\n' + msg + '\n' + '!' * 70)
            raise ValueError(msg)


    if df.shape[0] > 0:
        # Set SMCAccess the week before each year's first round, per age group.
        # The youngest group's assignment is extended down to age 0 so children
        # below its min_age already hold the IP by the time they age into the
        # eligible band (replaces the removed change_individual_property_at_age).
        youngest_min_by_year = df.groupby('year')['min_age'].min().to_dict()
        first_round_df = df.loc[df.groupby(['year', 'min_age', 'max_age'])['simday'].idxmin()]
        for idx, row in first_round_df.iterrows():
            max_age = row['max_age']
            is_youngest = np.isclose(row['min_age'], youngest_min_by_year[row['year']])
            ip_age_min = 0 if is_youngest else row['min_age']
            change_individual_property_scheduled(campaign, start_day=(row['simday'] - 7), coverage=1,
                                                 new_ip_key='SMCAccess', new_ip_value="Low",
                                                 target_age_min=ip_age_min, target_age_max=max_age)
            change_individual_property_scheduled(campaign, start_day=(row['simday'] - 7),
                                                 coverage=row['high_access_frac'],
                                                 new_ip_key='SMCAccess', new_ip_value="High",
                                                 target_age_min=ip_age_min, target_age_max=max_age)
    return {'admin_name': hfca}



def set_drug_params(config):
    # Amodaquine
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Cmax", 270)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Decay_T1", 0.7)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Decay_T2", 15.9)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_PKPD_C50", 150)  # param2: 55)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Vd", 1)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Max_Drug_IRBC_Kill", 0.23)  # param2: 0.2)

    # SulfadoxinePyrimethamine
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_Decay_T1", 11.5)
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_Decay_T2", 11.5)
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_PKPD_C50", 3)  # param2: 0.9)
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Max_Drug_IRBC_Kill", 0.28)  # param2: 0.28)


def update_smc_drug_params(config, row):
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_PKPD_C50", row['sp_c50'])
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Max_Drug_IRBC_Kill", row['sp_kill'])
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_PKPD_C50", row['aq_c50'])
    malaria_config.set_drug_param(config, 'Amodiaquine', "Max_Drug_IRBC_Kill", row['sp_kill'])


if __name__ == "__main__":
    pass
