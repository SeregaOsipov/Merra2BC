import re
from emac2wrf.lexical_utils import parse_mapping_rules

# mapping between MERRA2 species and WRF species
pipe_to_process = []  # set of mapping rules from merra to wrf


def initialise(merra2wrf_config):
    pipe_to_process = parse_mapping_rules(merra2wrf_config.spc_map)

    print("\nConversion MAP:")
    for item in pipe_to_process:
        print(str(item) + ":\t" )


def get_list_of_wrf_spec_by_merra_var(name):
    return pipe_to_process.get(name)

