import numpy as np

# SHAHEEN
root_path = '/project/k1090/osipovs'
# MISTRAL
root_path = '/work/mm0062/b302074'

wrf_dir = root_path + "/Data/AirQuality/AQABA/chem_106"
wrf_input_file = "wrfinput_d01"
wrf_bdy_file = "wrfbdy_d01"

wrf_met_dir = root_path + "/Data/AirQuality/AQABA/met_em/"
wrf_met_files = "met_em.d01.2017-0*"

# mera_files="svc_MERRA2_300.inst3_3d_aer_Nv.2010*"

mera_dir = root_path + "/Data/NASA/MERRA2/inst3_3d_aer_Nv/"
mera_files = "MERRA2_400.inst3_3d_aer_Nv.20170*"

# mera_dir = root_path + "/Data/NASA/MERRA2/inst3_3d_chm_Nv/"
# mera_files = "MERRA2_400.inst3_3d_chm_Nv.20170*"

do_IC = True
do_BC = True

# GOCART DUST ONLY
spc_map = ['DUST_1 -> 1.0*[DU001];1.e9',
           'DUST_2 -> 1.0*[DU002];1.e9',
           'DUST_3 -> 1.0*[DU003];1.e9',
           'DUST_4 -> 1.0*[DU004];1.e9',
           'DUST_5 -> 1.0*[DU005];1.e9']

# GOCART FULL
spc_map = ['DUST_1 -> 1.0*[DU001];1.e9',
           'DUST_2 -> 1.0*[DU002];1.e9',
           'DUST_3 -> 1.0*[DU003];1.e9',
           'DUST_4 -> 1.0*[DU004];1.e9',
           'DUST_5 -> 1.0*[DU005];1.e9',
           'SEAS_1 -> 1.0*[SS002];1.e9',
           'SEAS_2 -> 1.0*[SS003];1.e9',
           'SEAS_3 -> 1.0*[SS004];1.e9',
           'SEAS_4 -> 1.0*[SS005];1.e9',
           'so2 -> 0.453*[SO2];1.e6',
           'sulf -> 0.302*[SO4];1.e6',
           'BC1 -> 0.4143*[BCPHOBIC];1.e9', 'BC2 -> 0.4143*[BCPHILIC];1.e9',
           'OC1 -> 0.4143*[OCPHOBIC];1.e9', 'OC2 -> 0.4143*[OCPHILIC];1.e9',
           'dms -> 0.467*[DMS];1.e6']
# ,'msa -> 0.302*[MSA];1.e6'

# spc_map = ['so2 -> 0.453*[SO2];1.e6', 'sulf -> 0.302*[SO4];1.e6']
# spc_map = ['o3 -> 0.604*[O3];1.e6', 'co -> 1*[CO];1.e6']

'''
chem_106 case (MADE SORGAM)

To prescribe aerosol in the MADE scheme we need to specify two parameters:
mass concentration (linked to 3rd moment) and number concentration (0th moment).
3rd moment (mass) can be derived from MERRA2. For 0th moment (number concentration) I assumed fixed size.

See aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient
'''

# put all sulfate into accumulation mode
spc_map = ['ac0 -> 1.22*[SO4];1.e11', 'so4aj -> 1.0*[SO4];1.e9']

# put all dust into the coarse mode
spc_map += ['corn -> 1.68*[DU001]+1.68*[DU002]+1.68*[DU003]+1.68*[DU004]+1.68*[DU005];1.e7',
            'soila -> 1.0*[DU001]+1.0*[DU002]+1.0*[DU003]+1.0*[DU004]+1.0*[DU005];1.e9']

# put all sea salt into the coarse mode
# TODO: Check what should be the SD, for now assume the same
spc_map += ['corn -> 1.98*[SS001]+1.98*[SS002]+1.98*[SS003]+1.98*[SS004]+1.98*[SS005];1.e7',
            'seas -> 1.0*[SS001]+1.0*[SS002]+1.0*[SS003]+1.0*[SS004]+1.0*[SS005];1.e9']

# put all organic aerosols into acc mode
# TODO: emissions (speciation) breaks down carbon into i&j modes (20/80). Here I put everything to j mode.
# put BC into elemental carbon, consistent with emissions
spc_map += ['ac0 -> 6.41*[BCPHOBIC]+6.41*[BCPHILIC];1.e11', 'ecj -> 1.0*[BCPHOBIC]+1.0*[BCPHILIC];1.e9']

# assign OC into primary orgpaj in WRF, similar to the emissions
spc_map += ['ac0 -> 6.41*[OCPHOBIC]+6.41*[OCPHILIC];1.e11', 'orgpaj -> 1.0*[OCPHOBIC]+1.0*[OCPHILIC];1.e9']

# spc_map += ['so2 -> 0.453*[SO2];1.e6', ]
# spc_map = ['o3 -> 0.604*[O3];1.e6', 'co -> 1*[CO];1.e6']

# TODO: add DMS
#'dms -> 0.467*[DMS];1.e6']


'''
#CBMZ-MOSAIC_8bins  DUST only
spc_map =['oin_a02->4.92e-4*[DU001];1.e9',
          'oin_a03->8.94e-3*[DU001];1.e9',
          'oin_a04->0.14300*[DU001];1.e9',
          'oin_a05->0.84740*[DU001]+0.1520*[DU002];1.e9',
          'oin_a06->0.84800*[DU002]+0.4055*[DU003];1.e9',
          'oin_a07->0.59450*[DU003]+0.4480*[DU004];1.e9',
          'oin_a08->0.55200*[DU004]+1*[DU005];1.e9',

          'num_a02->0.213e14*[DU001];1',
          'num_a03->0.598e14*[DU001];1',
          'num_a04->1.192e14*[DU001];1',
          'num_a05->1.434e14*[DU001]+0.948e13*[DU002];1',
          'num_a06->2.086e13*[DU002]+0.358e13*[DU003];1',
          'num_a07->0.255e13*[DU003]+0.592e12*[DU004];1',
          'num_a08->0.296e12*[DU004]+1.655e11*[DU005];1']
'''
