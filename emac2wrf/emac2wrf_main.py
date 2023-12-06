import argparse
import xarray as xr
import time
import wrf_module
import netCDF4
from netCDF4 import Dataset
import numpy as np
import datetime as dt
from dateutil import rrule
from emac2wrf import emac_module
from emac2wrf.emac2wrf_config import get_emac2wrf_mapping
from emac2wrf.lexical_utils import parse_mapping_rule, get_unique_wrf_keys_from_mappings
import wrf as wrf
import pandas as pd
from distutils.util import strtobool

"""
EMAC2WRF was developed from MERRA2WRF alike
 
Steps:
1. This script will increment the IC and BC according to config
2. Change config and run again to include other trace gases and aerosols (if necessary)

Notes:
1. It is best to have both WRF and EMAC output as daily per file. MFDataset is possible, but probably less reliable
2. Either set start & end dates & hourly interval or None. Dates to process then will be deduced from the wrfbdy.

How to run:
gogomamba
python -u ${MERRA2BC}/emac2wrf/emac2wrf_main.py --start_date=2017-06-15_00:00:00 --end_date=2017-09-01_00:00:00 --hourly_interval=3


EMME 2017 / 2050 examples:
year=2017
scenario=

year=2050
scenario=CLE

data_dir=/work/mm0062/b302074/Data/AirQuality/EMME/${year}/${scenario}/IC_BC/
python -u ${MERRA2BC}/emac2wrf/emac2wrf_main.py --hourly_interval=3 --do_IC --do_BC --zero_out_first --emac_dir=${data_dir}/emac/ --wrf_dir=${data_dir} --wrf_met_dir=${data_dir}/met_em/ --wrf_met_files=met_em.d01.${year}-* >& log.emac2wrf
"""

#%%
root_path = '/project/k1090/osipovs'  # SHAHEEN
root_path = '/work/mm0062/b302074'  # MISTRAL

parser = argparse.ArgumentParser()
parser.add_argument("--start_date", help="YYYY-MM-DD_HH:MM:SS format", default=None)  # applies filter on the existing nc dates
parser.add_argument("--end_date", help="YYYY-MM-DD_HH:MM:SS format", default=None)  # applies filter on the existing nc dates
parser.add_argument("--hourly_interval", help="dt in dates to process", default=3)
parser.add_argument("--do_IC", help="process initial conditions?", action='store_true')
parser.add_argument("--do_BC", help="process boundary conditions?", action='store_true')
parser.add_argument("--zero_out_first", help="zero out fields in IC and BC first?", action='store_true')
# setup parent model (EMAC)
parser.add_argument("--emac_dir", help="folder containing emac output")  # /work/mm0062/b302011/script/Osipov/simulations/AQABA_2050
parser.add_argument("--emac_file_name_template", help="folder containing emac output", default='test01_________{date_time}_{stream}.nc')  # sim label has fixed width and then filled with ___
# setup downscaling model (WRF)
parser.add_argument("--wrf_dir", help="folder containing WRF IC & BC files")   # , default=data_dir+'/1-week-icbc/')
parser.add_argument("--wrf_input", help="use default wrfinput_d01", default='wrfinput_d01')
parser.add_argument("--wrf_bdy_file", help="use default wrfbdy_d01", default='wrfbdy_d01')
parser.add_argument("--wrf_met_dir", help="use default wrfbdy_d01", default=root_path + '/Data/AirQuality/EMME/IC_BC/met_em')
parser.add_argument("--wrf_met_files", help="met_em file names template")  # , default='met_em.d01.2017-0*')

parser.add_argument("--mode", "--port", help="the are only to support pycharm debugging")
args = parser.parse_args()

# debug section
# print("DEBUG SETTINGS ARE ON")
# year = 2017
# data_dir = '/work/mm0062/b302074/Data/AirQuality/EMME/2017/IC_BC/'
#
# args.do_IC = True
# args.do_BC = True
# args.zero_out_first = True
# args.emac_dir = data_dir + '/emac/'
# args.wrf_dir = data_dir
# args.wrf_met_dir = data_dir + '/met_em/'
# args.wrf_met_files = 'met_em.d01.{}-*'.format(year)

if args.start_date:
    start_date = dt.datetime.strptime(args.start_date, '%Y-%m-%d_%H:%M:%S')
if args.end_date:
    end_date = dt.datetime.strptime(args.end_date, '%Y-%m-%d_%H:%M:%S')
config = args
#%%
start_time = time.time()

# modules initialisation
mappings = get_emac2wrf_mapping()
wrf_module.initialise(args)

# prepare the mapping to process
print("\nConversion MAP:")
for mapping in mappings:
    print(mapping.mapping_rule_str + ":\t")

#%%
print('\nDeriving dates to process from wrf_bdy file\n')
wrfbdy_f = Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')
wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)
dates_to_process = wrf_bdy_dates

if args.start_date is not None:
    dates_to_process = dates_to_process[dates_to_process >= start_date]
if args.end_date is not None:
    dates_to_process = dates_to_process[dates_to_process <= end_date]

# dates_to_process = list(rrule.rrule(rrule.HOURLY, interval=int(args.hourly_interval), dtstart=start_date, until=end_date))  # has to match exactly dates in EMAC output

#%%


def get_emac_df(fp):  # aux function, which shifts long for EMAC specific case
    emac_df = xr.open_mfdataset(fp)
    emac_df = emac_df.assign_coords(lon=(((emac_df.lon + 180) % 360) - 180))  # EMAC lon spans [0; 360], switch to [-180; 180]
    emac_df = emac_df.sortby('lon')
    return emac_df


#%% IC
if config.do_IC:
    date = dates_to_process[0]
    print("INITIAL CONDITIONS: {}".format(date))

    # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_0000'), stream='ECHAM5')  # '%Y%m%d_%H%M'
    fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_0000'), stream='WRF_bc_met')  # '%Y%m%d_%H%M'
    print("Opening file: {}".format(fp))
    emac_df = get_emac_df(fp)
    emac_df = emac_df.sel(time=date)

    emac_pressure_rho_3d = emac_df['press'].load()  # press_ave
    emac_pressure_rho_3d_hi = emac_module.hor_interpolate_3d_field_on_wrf_grid(emac_pressure_rho_3d, wrf_module.ny, wrf_module.nx, wrf_module.xlon, wrf_module.xlat)

    met_file_name = wrf_module.get_met_file_by_time(date.strftime('%Y-%m-%d_%H_%M_%S'))
    met_fp = config.wrf_met_dir + "/" + met_file_name
    print("Opening metfile: " + met_fp)
    wrf_met_nc = Dataset(met_fp, 'r')
    WRF_PRES = wrf_module.get_pressure_from_metfile(wrf_met_nc)

    print("Opening wrfintput: " + config.wrf_input)
    wrfinput_f = Dataset(config.wrf_dir + "/" + config.wrf_input, 'r+')

    if config.zero_out_first:  # zero out the fields in wrfinput before increments
        unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
        print("Following fields will be zeroed out FIRST in ICs: {}".format(unique_wrf_keys))
        for wrf_key in unique_wrf_keys:
            wrfinput_f.variables[wrf_key][0, :] = 0  # minuscule might be better than 0

    print("INITIAL CONDITIONS: Increments")
    for mapping in mappings:
        print('Processing mapping: {}'.format(mapping.mapping_rule_str))
        pipe_to_process = parse_mapping_rule(mapping.mapping_rule_str)
        for rule_vo in pipe_to_process:
            fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_0000'), stream=mapping.output_stream)
            print("\t\t - Reading " + rule_vo['merra_key'] + " from {}".format(fp))
            emac_stream_df = get_emac_df(fp)
            emac_stream_df = emac_stream_df.sel(time=date)
            parent_var_da = emac_stream_df[rule_vo['merra_key']].load()
            emac_stream_df.close()

            print("\t\t - Horizontal interpolation of " + rule_vo['merra_key'] + " on WRF horizontal grid")
            parent_var_hi = emac_module.hor_interpolate_3d_field_on_wrf_grid(parent_var_da, wrf_module.ny, wrf_module.nx, wrf_module.xlon, wrf_module.xlat)
            print("\t\t - Vertical interpolation of " + rule_vo['merra_key'] + " on WRF vertical grid")
            # casting xarray to numpy accelerates things a lot
            parent_var_hvi = emac_module.ver_interpolate_3d_field_on_wrf_grid(parent_var_hi.to_numpy(), emac_pressure_rho_3d_hi.to_numpy(), WRF_PRES, wrf_module.nz, wrf_module.ny, wrf_module.nx)
            parent_var_hvi = np.flipud(parent_var_hvi)

            wrf_key = rule_vo['wrf_key']
            print("\t\t - Updating wrfinput field {}[0] += {} * {} * {:.1e}".format(wrf_key, rule_vo['merra_key'], rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
            wrfinput_f.variables[wrf_key][0, :] = wrfinput_f.variables[wrf_key][0, :] + parent_var_hvi * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent']

    wrfinput_f.close()
    emac_df.close()
    wrf_met_nc.close()
    print("DONE: INITIAL CONDITIONS")

#%%
if config.do_BC:
    print("\n\nSTART BOUNDARY CONDITIONS")

    print("Opening " + config.wrf_bdy_file)
    wrfbdy_f = Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')

    wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
    wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)
    dt = (dates_to_process[1]-dates_to_process[0]).total_seconds()

    for date in dates_to_process:
        print("\n\tBCs, next date: {}".format(date))

        # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_0000'), stream='ECHAM5')  # file with the pressure
        # emac_nc = Dataset(fp, 'r')
        # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m*'), stream='ECHAM5')  # MF
        # emac_nc = netCDF4.MFDataset(fp, 'r')
        # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_*'), stream='ECHAM5')  # daily MF. Due to EMAC restarts, files can split sub-daily
        # emac_nc = netCDF4.MFDataset(fp, 'r')
        # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_*'), stream='ECHAM5')  # daily MF. Due to EMAC restarts, files can split sub-daily
        fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_*'), stream='WRF_bc_met')  # daily MF. Due to EMAC restarts, files can split sub-daily
        emac_df = get_emac_df(fp)
        emac_df = emac_df.sel(time=date)

        time_index_in_wrfbdy = wrf_bdy_dates.get_loc(date)

        print("\tReading EMAC Pressure at {} from {}".format(date, fp))
        emac_pressure_rho_3d = emac_df['press']  # press_ave
        print("\tHorizontal interpolation of EMAC Pressure on WRF boundary")
        emac_pressure_rho_3d_hi = emac_module.hor_interpolate_3d_field_on_wrf_boubdary(emac_pressure_rho_3d, len(wrf_module.boundary_lons), wrf_module.boundary_lons, wrf_module.boundary_lats)

        met_file_name = wrf_module.get_met_file_by_time(date.strftime('%Y-%m-%d_%H_%M_%S'))
        met_fp = config.wrf_met_dir + "/" + met_file_name
        print("\tReading WRF Pressure from: {}".format(met_fp))
        wrf_met_nc = Dataset(met_fp, 'r')
        WRF_PRES = wrf_module.get_pressure_from_metfile(wrf_met_nc)
        WRF_PRES_BND = np.concatenate((WRF_PRES[:, :, 0], WRF_PRES[:, wrf_module.ny - 1, :],
                                       WRF_PRES[:, :, wrf_module.nx - 1], WRF_PRES[:, 0, :]), axis=1)
        wrf_met_nc.close()

        if config.zero_out_first:  # zero out the fields in wrfinput before increments
            unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
            print("Following fields will be zeroed out FIRST in BCs: {}".format(unique_wrf_keys))
            for wrf_key in unique_wrf_keys:
                # wrf_module.update_boundaries(0, wrfbdy_f, wrf_key, time_index)
                wrfbdy_f.variables[wrf_key + "_BXS"][time_index_in_wrfbdy, :] = 0  # BCs
                wrfbdy_f.variables[wrf_key + "_BXE"][time_index_in_wrfbdy, :] = 0
                wrfbdy_f.variables[wrf_key + "_BYS"][time_index_in_wrfbdy, :] = 0
                wrfbdy_f.variables[wrf_key + "_BYE"][time_index_in_wrfbdy, :] = 0
                wrfbdy_f.variables[wrf_key + "_BTXS"][time_index_in_wrfbdy, :] = 0  # BCs tendencies
                wrfbdy_f.variables[wrf_key + "_BTXE"][time_index_in_wrfbdy, :] = 0
                wrfbdy_f.variables[wrf_key + "_BTYS"][time_index_in_wrfbdy, :] = 0
                wrfbdy_f.variables[wrf_key + "_BTYE"][time_index_in_wrfbdy, :] = 0

        for mapping in mappings:
            print('Processing mapping: {}'.format(mapping.mapping_rule_str))
            pipe_to_process = parse_mapping_rule(mapping.mapping_rule_str)
            for rule_vo in pipe_to_process:
                print("\t\t - Reading " + rule_vo['merra_key'])

                merra_key = rule_vo['merra_key']
                wrf_key = rule_vo['wrf_key']

                fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_*'), stream=mapping.output_stream)  # daily MF
                emac_stream_df = get_emac_df(fp)
                emac_stream_df = emac_stream_df.sel(time=date)
                parent_var_da = emac_stream_df[rule_vo['merra_key']].load()
                emac_stream_df.close()

                print("\t\tHorizontal interpolation of " + merra_key + " on WRF boundary")
                parent_var_hi = emac_module.hor_interpolate_3d_field_on_wrf_boubdary(parent_var_da, len(wrf_module.boundary_lons), wrf_module.boundary_lons, wrf_module.boundary_lats)
                print("\t\tVertical interpolation of " + merra_key + " on WRF boundary")
                # casting to numpy() speeds up things a lot
                parent_var_hvi = emac_module.ver_interpolate_3dfield_on_wrf_boubdary(parent_var_hi.to_numpy(), emac_pressure_rho_3d_hi.to_numpy(), WRF_PRES_BND, wrf_module.nz, len(wrf_module.boundary_lons))
                parent_var_hvi = np.flipud(parent_var_hvi)

                print("\t\t - Updating wrfbdy: {}[{}] += {} * {} * {:.1e}".format(wrf_key, time_index_in_wrfbdy, merra_key, rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
                wrf_module.update_boundaries(parent_var_hvi * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent'], wrfbdy_f, wrf_key, time_index_in_wrfbdy)

        unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
        print('Prescribing uniform tendency in BC for these fields: {}'.format(unique_wrf_keys))
        for wrf_key in unique_wrf_keys:
            wrf_module.update_tendency_boundaries(wrfbdy_f, wrf_key, time_index_in_wrfbdy, dt)

        print("--- %s seconds ---" % (time.time() - start_time))
        emac_df.close()

    print("Closing " + config.wrf_bdy_file)
    wrfbdy_f.close()

    print("FINISH BOUNDARY CONDITIONS")

print("--- %s seconds ---" % (time.time() - start_time))
