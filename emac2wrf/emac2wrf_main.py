import argparse
import merra2wrf_config
import time
import wrf_module
import merra2wrf_mapper
import netCDF4
from netCDF4 import Dataset
import numpy as np
from datetime import datetime
import datetime as dt
from dateutil import rrule
from emac2wrf import emac2wrf_config, emac_module
from emac2wrf.bc_ic_utils import is_child_domain_covered_by_parent_domain
from emac2wrf.emac2wrf_config import get_emac2wrf_mapping
from emac2wrf.lexical_utils import parse_mapping_rule, get_unique_wrf_keys_from_mappings
import wrf as wrf
import pandas as pd
from distutils.util import strtobool

"""
EMAC2WRF was developed from MERRA2WRF alike
 
Steps:
1. Run zero_fields.py first to zero out the IC and BC for specified chem species
2. Run this file, which will increment the IC and BC according to config
3. Change config and run again to include other trace gases and aerosols

Notes:
1. It is best, when both WRF and EMAC output are in daily per file. MFDataset is possible, but probably less reliable
2. Either set start & end dates & hourly interval or None. Dates to process then will be deduced from the wrfbdy.

How to run (to not specify dates to derive them from wrf files):
gogomamba
data_dir=/work/mm0062/b302074/Data/AirQuality/EMME/IC_BC/2017
python -u ${MERRA2BC}/emac2wrf/emac2wrf_main.py  --hourly_interval=3  --do_IC=True --zero_out_first=True --emac_dir=${data_dir}/emac/ --wrf_dir=${data_dir}/1-week-icbc/ --wrf_met_dir=${data_dir}/met_em/
"""

#%%
root_path = '/project/k1090/osipovs'  # SHAHEEN
root_path = '/work/mm0062/b302074'  # MISTRAL

# debug
# data_dir='/work/mm0062/b302074/Data/AirQuality/EMME/IC_BC/2017'

parser = argparse.ArgumentParser()
parser.add_argument("--start_date", help="YYYY-MM-DD_HH:MM:SS format", default=None)  # None means derive  # default='2017-01-01_03:00:00')  #
parser.add_argument("--end_date", help="YYYY-MM-DD_HH:MM:SS format", default=None)  # default='2017-09-01_00:00:00')  #
parser.add_argument("--hourly_interval", type=int, help="dt in dates to process", default=3)
parser.add_argument("--do_IC", type=strtobool, help="process initial conditions?", default=False)
parser.add_argument("--do_BC", type=strtobool, help="process boundary conditions?", default=True)
parser.add_argument("--zero_out_first", help="zero out fields in IC and BC first?", default=True)
# setup parent model (EMAC)
parser.add_argument("--emac_dir", help="folder containing emac output")  # , default=data_dir+'/emac/')  # /work/mm0062/b302011/script/Osipov/simulations/AQABA_2050
parser.add_argument("--emac_file_name_template", help="folder containing emac output", default='MIM_STD________{date_time}_{stream}.nc')  # sim label has fixed width and then filled with ___
# setup downscaling model (WRF)
parser.add_argument("--wrf_dir", help="folder containing WRF IC & BC files")   # , default=data_dir+'/1-week-icbc/')
parser.add_argument("--wrf_input", help="use default wrfinput_d01", default='wrfinput_d01')
parser.add_argument("--wrf_bdy_file", help="use default wrfbdy_d01", default='wrfbdy_d01')
parser.add_argument("--wrf_met_dir", help="use default wrfbdy_d01")  #, default=data_dir+'/met_em')
parser.add_argument("--wrf_met_files", help="met_em file names template", default='met_em.d01.20*')

parser.add_argument("--mode", "--port", help="the are only to support pycharm debugging")
args = parser.parse_args()

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
emac_module.initialise(args, mappings)

# prepare the mapping to process
# pipe_to_process = parse_mapping_rules(config.spc_map)
print("\nConversion MAP:")
for mapping in mappings:
    print(mapping.mapping_rule_str + ":\t")

#%% -----------------------------
# Sanity checks:
# check species availability in wrf and in merra files
# for rule_vo in pipe_to_process:
#     if rule_vo['merra_key'] not in emac_module.vars:
#         raise Exception("Could not find variable " + rule_vo['merra_key'] + " in parent model file. Exiting...")

# for rule_vo in pipe_to_process:
#     if rule_vo['wrf_key'] not in wrf_module.wrf_vars:
#         raise Exception("Could not find variable " + rule_vo['wrf_key'] + " in WRF input file. Exiting...")

# check domain coverage
# if not is_child_domain_covered_by_parent_domain(emac_module, wrf_module):
#     raise Exception("Child domain (WRF) is not fully covered by Parent domain area. Exiting...")

# dates_to_process = wrf_module.dates.keys() & emac_module.dates.keys()
# Or get all dates from wrf_bdy
if args.start_date is None or args.end_date is None:
    print('\nDates to process will be derived from wrf_bdy file\n')
    wrfbdy_f = Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')
    wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
    wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)
    dates_to_process = wrf_bdy_dates
else:
    print('Dates to process will be generated between start and end dates')
    # manually setup dates to process. Check the coverage in WRF
    dates_to_process = list(rrule.rrule(rrule.HOURLY, interval=args.hourly_interval, dtstart=start_date, until=end_date))  # has to match exactly dates in EMAC output

# check that merra2 time is covered by wrf time
# if len(time_intersection) != len(wrf_module.wrf_times):
#     print('These datetimes are missing in MERRA2 dataset:')
#     time_intersection = dict.fromkeys(time_intersection, 0)
#     for key in wrf_module.wrf_times.keys():
#         if not key in time_intersection:
#             print(key)
#     raise Exception("WRF time range is not fully covered by MERRA2 time range. Exiting...")

# -----------------------------
# sorting times for processing
# time_intersection = sorted(time_intersection, key=lambda x: time.mktime(time.strptime(x, "%Y-%m-%d_%H:%M:%S")))
# print("\nTimes for processing: {}".format(time_intersection))
# exit()

if config.do_IC:
    date = dates_to_process[0]
    print("INITIAL CONDITIONS: {}".format(date))

    fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_%H%M'), stream='ECHAM5')
    print("Opening file: {}".format(fp))
    emac_nc = Dataset(fp, 'r')

    emac_nc_dates = netCDF4.num2date(emac_nc['time'][:], emac_nc['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    emac_nc_dates = pd.to_datetime(emac_nc_dates)
    time_index_in_emac = emac_nc_dates.get_loc(date)

    emac_pressure_rho_3d = emac_module.get_3d_field_by_time_index(time_index_in_emac, emac_nc, 'press')  # press_ave
    emac_pressure_rho_3d_hi = emac_module.hor_interpolate_3d_field_on_wrf_grid(emac_pressure_rho_3d, wrf_module.ny, wrf_module.nx, wrf_module.xlon, wrf_module.xlat)

    # visual debug
    # import matplotlib.pyplot as plt
    # plt.figure(num=0)
    # plt.clf()
    # plt.contourf(emac_module.lon, emac_module.lat, emac_pressure_rho_3d[-1])
    # plt.figure(num=1)
    # plt.clf()
    # plt.contourf(wrf_module.xlon, wrf_module.xlat, emac_pressure_rho_3d_hi[-1])

    met_file_name = wrf_module.get_met_file_by_time(date.strftime('%Y-%m-%d_%H:%M:%S'))
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
            print("\t\t - Reading " + rule_vo['merra_key'])
            fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_%H%M'), stream=mapping.emac_stream)
            emac_stream_nc = Dataset(fp)
            # parent_var = emac_module.get_3d_field_by_time(date, emac_stream_nc, rule_vo['merra_key'])
            parent_var = emac_module.get_3d_field_by_time_index(0, emac_stream_nc, rule_vo['merra_key'])
            emac_stream_nc.close()

            print("\t\t - Horizontal interpolation of " + rule_vo['merra_key'] + " on WRF horizontal grid")
            parent_var_hi = emac_module.hor_interpolate_3d_field_on_wrf_grid(parent_var, wrf_module.ny, wrf_module.nx, wrf_module.xlon, wrf_module.xlat)

            print("\t\t - Vertical interpolation of " + rule_vo['merra_key'] + " on WRF vertical grid")
            parent_var_hvi = emac_module.ver_interpolate_3dfield_on_wrf_grid(parent_var_hi, emac_pressure_rho_3d_hi, WRF_PRES, wrf_module.nz, wrf_module.ny, wrf_module.nx)
            parent_var_hvi = np.flipud(parent_var_hvi)

            wrf_key = rule_vo['wrf_key']
            print("\t\t - Updating wrfinput field {}[0] += {} * {} * {:.1e}".format(wrf_key, rule_vo['merra_key'], rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
            wrfinput_f.variables[wrf_key][0, :] = wrfinput_f.variables[wrf_key][0, :] + parent_var_hvi * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent']

    wrfinput_f.close()
    emac_nc.close()
    wrf_met_nc.close()
    print("DONE: INITIAL CONDITIONS")

if config.do_BC:
    print("\n\nSTART BOUNDARY CONDITIONS")

    print("Opening " + config.wrf_bdy_file)
    wrfbdy_f = Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')

    wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
    wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)

    # difference between two given times
    dt = (dates_to_process[1]-dates_to_process[0]).total_seconds()

    for date in dates_to_process:
        print("\n\tBCs, next date: {}".format(date))

        # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_0000'), stream='ECHAM5')  # file with the pressure
        # emac_nc = Dataset(fp, 'r')
        # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m*'), stream='ECHAM5')  # MF
        # emac_nc = netCDF4.MFDataset(fp, 'r')
        fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_*'), stream='ECHAM5')  # daily MF. Due to EMAC restarts, files can split sub-daily
        emac_nc = netCDF4.MFDataset(fp, 'r')

        emac_nc_dates = netCDF4.num2date(emac_nc['time'][:], emac_nc['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        emac_nc_dates = pd.to_datetime(emac_nc_dates)

        time_index_in_emac = emac_nc_dates.get_loc(date)
        time_index_in_wrfbdy = wrf_bdy_dates.get_loc(date)

        print("\tReading EMAC Pressure at time index {} from {}".format(time_index_in_emac, fp))
        # emac_pressure_rho_3d = emac_module.get_3d_field_by_time(date, emac_nc, 'press')  # press_ave
        emac_pressure_rho_3d = emac_module.get_3d_field_by_time_index(time_index_in_emac, emac_nc, 'press')  # press_ave
        print("\tHorizontal interpolation of EMAC Pressure on WRF boundary")
        emac_pressure_rho_3d_hi = emac_module.hor_interpolate_3dfield_on_wrf_boubdary(emac_pressure_rho_3d, len(wrf_module.boundary_lons), wrf_module.boundary_lons, wrf_module.boundary_lats)

        met_file_name = wrf_module.get_met_file_by_time(date.strftime('%Y-%m-%d_%H:%M:%S'))
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

                # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_0000'), stream=mapping.emac_stream)
                # emac_stream_nc = Dataset(fp)
                # fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m*'), stream=mapping.emac_stream)  # MF
                # emac_stream_nc = netCDF4.MFDataset(fp)
                fp = config.emac_dir + config.emac_file_name_template.format(date_time=date.strftime('%Y%m%d_*'), stream=mapping.emac_stream)  # daily MF
                emac_stream_nc = netCDF4.MFDataset(fp)
                parent_var = emac_module.get_3d_field_by_time_index(time_index_in_emac, emac_stream_nc, rule_vo['merra_key'])
                emac_stream_nc.close()

                print("\t\tHorizontal interpolation of " + merra_key + " on WRF boundary")
                parent_var_hi = emac_module.hor_interpolate_3dfield_on_wrf_boubdary(parent_var, len(wrf_module.boundary_lons), wrf_module.boundary_lons, wrf_module.boundary_lats)
                print("\t\tVertical interpolation of " + merra_key + " on WRF boundary")
                parent_var_hvi = emac_module.ver_interpolate_3dfield_on_wrf_boubdary(parent_var_hi, emac_pressure_rho_3d_hi, WRF_PRES_BND, wrf_module.nz, len(wrf_module.boundary_lons))
                parent_var_hvi = np.flipud(parent_var_hvi)

                print("\t\t - Updating wrfbdy: {}[{}] += {} * {} * {:.1e}".format(wrf_key, time_index_in_wrfbdy, merra_key, rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
                wrf_module.update_boundaries(parent_var_hvi * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent'], wrfbdy_f, wrf_key, time_index_in_wrfbdy)

        unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
        print('Prescribing uniform tendency in BC for these fields: {}'.format(unique_wrf_keys))
        for wrf_key in unique_wrf_keys:
            wrf_module.update_tendency_boundaries(wrfbdy_f, wrf_key, time_index_in_wrfbdy, dt)

        print("--- %s seconds ---" % (time.time() - start_time))
        emac_nc.close()

    print("Closing " + config.wrf_bdy_file)
    wrfbdy_f.close()

    print("FINISH BOUNDARY CONDITIONS")

print("--- %s seconds ---" % (time.time() - start_time))
