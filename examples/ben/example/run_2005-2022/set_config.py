import manifest
import params
from emodpy_malaria.vector_config import set_species_param
from snt.hbhi.set_up_general import initialize_config


def set_config(config):
    # BASIC SETUP
    initialize_config(config, manifest, params.years, serialize=params.serialize)

    config.parameters.x_Temporary_Larval_Habitat = 1
    config.parameters.x_Base_Population = 1
    config.parameters.x_Birth = 1

    config.parameters.Report_Event_Recorder = 1
    config.parameters.Report_Event_Recorder_Individual_Properties = []
    config.parameters.Report_Event_Recorder_Events = ['Received_NMF_Treatment', 'Received_Severe_Treatment',
                                                      'Received_Treatment']
    config.parameters.Report_Event_Recorder_Ignore_Events_In_List = 0
    config.parameters.Custom_Individual_Events = ['Received_NMF_Treatment']

    set_species_param(config, 'arabiensis', 'Anthropophily', 0.88, overwrite=True)
    set_species_param(config, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
    set_species_param(config, 'funestus', 'Anthropophily', 0.5, overwrite=True)
    set_species_param(config, 'funestus', 'Indoor_Feeding_Fraction', 0.86, overwrite=True)
    set_species_param(config, 'gambiae', 'Anthropophily', 0.74, overwrite=True)
    set_species_param(config, 'gambiae', 'Indoor_Feeding_Fraction', 0.9, overwrite=True)

    return config
