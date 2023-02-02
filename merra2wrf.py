import argparse
import time
import merra2_module
import wrf_module
from dateutil import rrule
import numpy as np
from datetime import datetime
from emac2wrf.bc_ic_utils import is_child_domain_covered_by_parent_domain
from emac2wrf.lexical_utils import parse_mapping_rule, get_unique_wrf_keys_from_mappings
from merra2wrf_config import get_merra2wrf_mapping
import wrf as wrf
import pandas as pd
import datetime as dt
import netCDF4 as nc
import xarray as xr
from climpy.utils.merra_utils import derive_merra2_pressure_profile

"""
How to run (AREAD 2022 example):
gogomamba
year=2022
data_dir=/work/mm0062/b302074/Data/AirQuality/AREAD/IC_BC/
python -u ${MERRA2BC}/merra2wrf.py --hourly_interval=3 --do_IC --do_BC --zero_out_first --wrf_dir=${data_dir} --wrf_met_dir=${data_dir}/met_em/ --wrf_met_files=met_em.d01.${year}-* >& log.emac2wrf
"""

#%% TODO: try to merge with emac2wrf
root_path = '/project/k1090/osipovs'  # SHAHEEN
root_path = '/work/mm0062/b302074'  # MISTRAL

parser = argparse.ArgumentParser()
parser.add_argument("--start_date", help="YYYY-MM-DD_HH:MM:SS format", default=None)  # default='2017-08-31_21:00:00')  #
parser.add_argument("--end_date", help="YYYY-MM-DD_HH:MM:SS format", default=None)  # default='2017-09-01_00:00:00')  #
parser.add_argument("--hourly_interval", help="dt in dates to process", default=3)
parser.add_argument("--do_IC", help="process initial conditions?", action='store_true')
parser.add_argument("--do_BC", help="process boundary conditions?", action='store_true')
parser.add_argument("--zero_out_first", help="zero out fields in IC and BC first?", action='store_true')
# setup parent model
parser.add_argument("--merra2_dir", default='/work/mm0062/b302074/Data/NASA/MERRA2/')  # stream
parser.add_argument("--merra2_file_name_template", default='{stream}/MERRA2_400.{stream}.{date_time}.nc4')
# setup downscaling model (WRF)
parser.add_argument("--wrf_dir", help="folder containing WRF IC & BC files")   # , default=/Data/AirQuality/AQABA/IC_BC')
parser.add_argument("--wrf_input", help="use default wrfinput_d01", default='wrfinput_d01')
parser.add_argument("--wrf_bdy_file", help="use default wrfbdy_d01", default='wrfbdy_d01')
parser.add_argument("--wrf_met_dir", help="use default wrfbdy_d01")  # default=root_path + '/Data/AirQuality/EMME/IC_BC/met_em'
parser.add_argument("--wrf_met_files", help="met_em file names template")  # , default='met_em.d01.2017-0*')

parser.add_argument("--mode", "--port", "--host", help="the are only to support pycharm debugging")
args = parser.parse_args()

if args.start_date:
    start_date = dt.datetime.strptime(args.start_date, '%Y-%m-%d_%H:%M:%S')
if args.end_date:
    end_date = dt.datetime.strptime(args.end_date, '%Y-%m-%d_%H:%M:%S')
config = args

#%% DEBUG settings
#args.wrf_dir = '/work/mm0062/b302074/Data/AirQuality/AREAD/IC_BC/debugICBC/'
#args.wrf_met_dir = '/work/mm0062/b302074/Data/AirQuality/AREAD/IC_BC/met_em/'
#args.wrf_met_files = 'met_em.d01.2022-0*'
# config.do_IC = True
#config.do_BC = True
# config.zero_out_first = False  # just speed up
#%%
start_time = time.time()

mappings = get_merra2wrf_mapping()
wrf_module.initialise(args)
# merra2_module.initialise(args, mappings)

print("\nConversion MAP:")
for mapping in mappings:
    print(mapping.mapping_rule_str + ":\t")

#%%
if args.start_date is None or args.end_date is None:
    print('\nDates to process will be derived from wrf_bdy file\n')
    wrfbdy_f = nc.Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')
    wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
    wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)
    dates_to_process = wrf_bdy_dates
else:
    print('Dates to process will be generated between start and end dates')
    # manually setup dates to process. Check the coverage in WRF
    dates_to_process = list(rrule.rrule(rrule.HOURLY, interval=int(args.hourly_interval), dtstart=start_date, until=end_date))  # has to match exactly dates in EMAC output

#%% IC
if config.do_IC:
    date = dates_to_process[0]
    print("INITIAL CONDITIONS: {}".format(date))

    fp = config.merra2_dir + config.merra2_file_name_template.format(date_time=date.strftime('%Y%m%d'), stream='inst3_3d_aer_Nv')
    print("Opening file: {}".format(fp))
    merra_df = xr.open_dataset(fp)
    merra_df = merra_df.sel(time=date)
    merra_df = merra_df.isel(lev=slice(None, None, -1))  # flip vertical dimension

    pressure_stag, pressure_rho = derive_merra2_pressure_profile(merra_df)
    merra_pressure_rho_3d = pressure_rho
    merra_pressure_rho_3d_hi = merra2_module.hor_interpolate_3d_field_on_wrf_grid(merra_pressure_rho_3d, wrf_module.ny, wrf_module.nx, wrf_module.xlon, wrf_module.xlat)

    met_file_name = wrf_module.get_met_file_by_time(date.strftime('%Y-%m-%d_%H:%M:%S'))
    met_fp = config.wrf_met_dir + "/" + met_file_name
    print("Opening metfile: " + met_fp)
    wrf_met_nc = nc.Dataset(met_fp, 'r')
    WRF_PRES = wrf_module.get_pressure_from_metfile(wrf_met_nc)

    print("Opening wrfintput: " + config.wrf_input)
    wrfinput_f = nc.Dataset(config.wrf_dir + "/" + config.wrf_input, 'r+')

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
            fp = config.merra2_dir + config.merra2_file_name_template.format(date_time=date.strftime('%Y%m%d'), stream=mapping.output_stream)
            merra_stream_df = xr.open_dataset(fp)
            merra_stream_df = merra_stream_df.sel(time=date)
            merra_stream_df = merra_stream_df.isel(lev=slice(None, None, -1))  # flip vertical dimension

            parent_var = merra_stream_df[rule_vo['merra_key']]
            print("\t\t - Horizontal interpolation of " + rule_vo['merra_key'] + " on WRF horizontal grid")
            parent_var_hi = merra2_module.hor_interpolate_3d_field_on_wrf_grid(parent_var, wrf_module.ny, wrf_module.nx, wrf_module.xlon, wrf_module.xlat)

            print("\t\t - Vertical interpolation of " + rule_vo['merra_key'] + " on WRF vertical grid")
            parent_var_hvi = merra2_module.ver_interpolate_3d_field_on_wrf_grid(parent_var_hi, merra_pressure_rho_3d_hi, WRF_PRES, wrf_module.nz, wrf_module.ny, wrf_module.nx)
            parent_var_hvi = np.flipud(parent_var_hvi)

            wrf_key = rule_vo['wrf_key']
            print("\t\t - Updating wrfinput field {}[0] += {} * {} * {:.1e}".format(wrf_key, rule_vo['merra_key'], rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
            wrfinput_f.variables[wrf_key][0, :] = wrfinput_f.variables[wrf_key][0, :] + parent_var_hvi * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent']

            merra_stream_df.close()

    wrfinput_f.close()
    merra_df.close()
    wrf_met_nc.close()
    print("DONE: INITIAL CONDITIONS")

#%% BC
if config.do_BC:
    print("\n\nSTART BOUNDARY CONDITIONS")

    print("Opening " + config.wrf_bdy_file)
    wrfbdy_f = nc.Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')

    wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
    wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)
    dt = (dates_to_process[1] - dates_to_process[0]).total_seconds()

    for date in dates_to_process:
        print("\n\tBCs, next date: {}".format(date))

        fp = config.merra2_dir + config.merra2_file_name_template.format(date_time=date.strftime('%Y%m%d'), stream='inst3_3d_aer_Nv')
        print("Opening file: {}".format(fp))
        merra_df = xr.open_dataset(fp)
        merra_df = merra_df.sel(time=date)
        merra_df = merra_df.isel(lev=slice(None, None, -1))  # flip vertical dimension

        pressure_stag, pressure_rho = derive_merra2_pressure_profile(merra_df)
        merra_pressure_rho_3d = pressure_rho
        print("\tHorizontal interpolation of MERRA Pressure on WRF boundary")
        merra_pressure_rho_3d_hi = merra2_module.hor_interpolate_3dfield_on_wrf_boubdary(merra_pressure_rho_3d, len(wrf_module.boundary_lons), wrf_module.boundary_lons, wrf_module.boundary_lats)

        met_file_name = wrf_module.get_met_file_by_time(date.strftime('%Y-%m-%d_%H:%M:%S'))
        met_fp = config.wrf_met_dir + "/" + met_file_name
        print("\tReading WRF Pressure from: {}".format(met_fp))
        wrf_met_nc = nc.Dataset(met_fp, 'r')
        WRF_PRES = wrf_module.get_pressure_from_metfile(wrf_met_nc)
        WRF_PRES_BND = np.concatenate((WRF_PRES[:, :, 0], WRF_PRES[:, wrf_module.ny - 1, :],
                                       WRF_PRES[:, :, wrf_module.nx - 1], WRF_PRES[:, 0, :]), axis=1)
        wrf_met_nc.close()

        time_index_in_wrfbdy = wrf_bdy_dates.get_loc(date)

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

            fp = config.merra2_dir + config.merra2_file_name_template.format(date_time=date.strftime('%Y%m%d'), stream=mapping.output_stream)
            merra_stream_df = xr.open_dataset(fp)
            merra_stream_df = merra_stream_df.sel(time=date)
            merra_stream_df = merra_stream_df.isel(lev=slice(None, None, -1))  # flip vertical dimension

            for rule_vo in pipe_to_process:
                print("\t\t - Reading " + rule_vo['merra_key'])

                merra_key = rule_vo['merra_key']
                wrf_key = rule_vo['wrf_key']

                parent_var = merra_stream_df[rule_vo['merra_key']]
                print("\t\tHorizontal interpolation of " + merra_key + " on WRF boundary")
                parent_var_hi = merra2_module.hor_interpolate_3dfield_on_wrf_boubdary(parent_var, len(wrf_module.boundary_lons), wrf_module.boundary_lons, wrf_module.boundary_lats)
                print("\t\tVertical interpolation of " + merra_key + " on WRF boundary")
                parent_var_hvi = merra2_module.ver_interpolate_3dfield_on_wrf_boubdary(parent_var_hi, merra_pressure_rho_3d_hi, WRF_PRES_BND, wrf_module.nz, len(wrf_module.boundary_lons))
                parent_var_hvi = np.flipud(parent_var_hvi)

                print("\t\t - Updating wrfbdy: {}[{}] += {} * {} * {:.1e}".format(wrf_key, time_index_in_wrfbdy, merra_key, rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
                wrf_module.update_boundaries(parent_var_hvi * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent'], wrfbdy_f, wrf_key, time_index_in_wrfbdy)

            merra_stream_df.close()

        unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
        print('Prescribing uniform tendency in BC for these fields: {}'.format(unique_wrf_keys))
        for wrf_key in unique_wrf_keys:
            wrf_module.update_tendency_boundaries(wrfbdy_f, wrf_key, time_index_in_wrfbdy, dt)

        print("--- %s seconds ---" % (time.time() - start_time))
        merra_df.close()

    print("Closing " + config.wrf_bdy_file)
    wrfbdy_f.close()

    print("FINISH BOUNDARY CONDITIONS")

#%% time stats
print("--- %s seconds ---" % (time.time() - start_time))
