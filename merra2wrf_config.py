import numpy as np
from emac2wrf.bc_ic_utils import IcBcMappingItem

class Config(object):
    pass


# old implementation for varios chem opts preserved for reference purposes
# def get_merra2wrf_config():
#     config = Config()
#
#     config.root_path = '/project/k1090/osipovs'  # SHAHEEN
#     config.root_path = '/work/mm0062/b302074'  # MISTRAL
#
#     config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/"
#     config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/2017-09-01"
#     config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/2017-10-31"
#     # config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/2018-03-08"
#     config.wrf_input_file = "wrfinput_d01"
#     config.wrf_bdy_file = "wrfbdy_d01"
#
#     config.wrf_met_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/met_em/"
#     config.wrf_met_files = "met_em.d01.2017-0*"
#
#     config.mera_dir = config.root_path + "/Data/NASA/MERRA2/inst3_3d_aer_Nv/"
#     config.mera_files = "MERRA2_400.inst3_3d_aer_Nv.201*"  # "MERRA2_400.inst3_3d_aer_Nv.20170*"
#
#     config.do_IC = False  # True
#     config.do_BC = True
#
#     # GOCART DUST ONLY
#     config.spc_map = ['DUST_1 -> 1.0*[DU001];1.e9',
#                'DUST_2 -> 1.0*[DU002];1.e9',
#                'DUST_3 -> 1.0*[DU003];1.e9',
#                'DUST_4 -> 1.0*[DU004];1.e9',
#                'DUST_5 -> 1.0*[DU005];1.e9']
#
#     # GOCART FULL
#     config.spc_map = ['DUST_1 -> 1.0*[DU001];1.e9',
#                'DUST_2 -> 1.0*[DU002];1.e9',
#                'DUST_3 -> 1.0*[DU003];1.e9',
#                'DUST_4 -> 1.0*[DU004];1.e9',
#                'DUST_5 -> 1.0*[DU005];1.e9',
#                'SEAS_1 -> 1.0*[SS002];1.e9',
#                'SEAS_2 -> 1.0*[SS003];1.e9',
#                'SEAS_3 -> 1.0*[SS004];1.e9',
#                'SEAS_4 -> 1.0*[SS005];1.e9',
#                'so2 -> 0.453*[SO2];1.e6',
#                'sulf -> 0.302*[SO4];1.e6',
#                'BC1 -> 0.4143*[BCPHOBIC];1.e9', 'BC2 -> 0.4143*[BCPHILIC];1.e9',
#                'OC1 -> 0.4143*[OCPHOBIC];1.e9', 'OC2 -> 0.4143*[OCPHILIC];1.e9',
#                'dms -> 0.467*[DMS];1.e6']
#     # ,'msa -> 0.302*[MSA];1.e6'
#
#     # config.spc_map = ['so2 -> 0.453*[SO2];1.e6', 'sulf -> 0.302*[SO4];1.e6']
#     # config.spc_map = ['o3 -> 0.604*[O3];1.e6', 'co -> 1*[CO];1.e6']
#
#     '''
#     #CBMZ-MOSAIC_8bins  DUST only
#     config.spc_map =['oin_a02->4.92e-4*[DU001];1.e9',
#               'oin_a03->8.94e-3*[DU001];1.e9',
#               'oin_a04->0.14300*[DU001];1.e9',
#               'oin_a05->0.84740*[DU001]+0.1520*[DU002];1.e9',
#               'oin_a06->0.84800*[DU002]+0.4055*[DU003];1.e9',
#               'oin_a07->0.59450*[DU003]+0.4480*[DU004];1.e9',
#               'oin_a08->0.55200*[DU004]+1*[DU005];1.e9',
#
#               'num_a02->0.213e14*[DU001];1',
#               'num_a03->0.598e14*[DU001];1',
#               'num_a04->1.192e14*[DU001];1',
#               'num_a05->1.434e14*[DU001]+0.948e13*[DU002];1',
#               'num_a06->2.086e13*[DU002]+0.358e13*[DU003];1',
#               'num_a07->0.255e13*[DU003]+0.592e12*[DU004];1',
#               'num_a08->0.296e12*[DU004]+1.655e11*[DU005];1']
#     '''
#
#     return config


def get_merra2wrf_mapping():
    '''
    chem 100/106 case (MADE SORGAM)

    To prescribe aerosol in the MADE scheme we need to specify two parameters:
    mass concentration (linked to 3rd moment) and number concentration (0th moment).
    3rd moment (mass) can be derived from MERRA2. For 0th moment (number concentration) I assumed fixed size.

    See aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient
    '''


    gas_stream = 'inst3_3d_chm_Nv'
    aerosol_stream = 'inst3_3d_aer_Nv'
    mappings = [
        # chem 100: Gases
        IcBcMappingItem('so2 -> 0.453*[SO2];1.e6', aerosol_stream),
        IcBcMappingItem('o3 -> 0.604*[O3];1.e6', gas_stream),
        IcBcMappingItem('co -> 0.966*[CO];1.e6', gas_stream),
        # TODO: setup CH4

        # chem 100: Aerosols
        # TODO: BC & SO4 should likely be assigned to the i-mode. Check emissions (speciation) break down.
        # TODO: Here currently everything goes to j mode.
        # put all SO4 into accumulation mode
        IcBcMappingItem('ac0 -> 3.56*[SO4];1.e11', aerosol_stream),  # MERRA2_400.inst3_3d_aer_Nv.20*
        IcBcMappingItem('so4aj -> 1.0*[SO4];1.e9', aerosol_stream),
        # put all dust into the coarse mode
        IcBcMappingItem('corn -> 1.68*[DU001]+1.68*[DU002]+1.68*[DU003]+1.68*[DU004]+1.68*[DU005];1.e7', aerosol_stream),
        IcBcMappingItem('soila -> 1.0*[DU001]+1.0*[DU002]+1.0*[DU003]+1.0*[DU004]+1.0*[DU005];1.e9', aerosol_stream),
        # put all sea salt into the coarse mode.
        # TODO: Check what should be the SD, for now assume the same
        IcBcMappingItem('corn -> 1.98*[SS001]+1.98*[SS002]+1.98*[SS003]+1.98*[SS004]+1.98*[SS005];1.e7', aerosol_stream),
        IcBcMappingItem('seas -> 1.0*[SS001]+1.0*[SS002]+1.0*[SS003]+1.0*[SS004]+1.0*[SS005];1.e9', aerosol_stream),
        # put all organic aerosols into acc mode
        # put BC into elemental carbon, consistent with emissions
        IcBcMappingItem('ac0 -> 6.41*[BCPHOBIC]+6.41*[BCPHILIC];1.e11', aerosol_stream),
        IcBcMappingItem('ecj -> 1.0*[BCPHOBIC]+1.0*[BCPHILIC];1.e9', aerosol_stream),
        # assign OC into primary orgpaj in WRF, similar to the emissions
        IcBcMappingItem('ac0 -> 6.41*[OCPHOBIC]+6.41*[OCPHILIC];1.e11', aerosol_stream),
        IcBcMappingItem('orgpaj -> 1.0*[OCPHOBIC]+1.0*[OCPHILIC];1.e9', aerosol_stream),
    ]
    return mappings