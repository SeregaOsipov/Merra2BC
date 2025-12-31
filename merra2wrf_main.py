import argparse
import time

import regridding_utils
from dateutil import rrule
import numpy as np
from datetime import datetime

from climpy.utils.time_utils import argsparse_datetime
from climpy.utils.wrf_utils import derive_wrf_pressure_from_met_and_input_files
from lexical_utils import parse_mapping_rule, get_unique_wrf_keys_from_mappings
from mapping_utils import get_merra2wrf_mapping
import wrf as wrf
import pandas as pd
import datetime as dt
import netCDF4 as nc
import xarray as xr
from climpy.utils.merra_utils import derive_merra2_pressure_profile, generate_merra2_file_name
from regridding_utils import flatten_and_line_up_boundaries

"""

TODO: I SHOULD BE ABLE TO MERGE ALL THREE (CAMS, MERRA2, EMAC) into one


How to run AREAD 2022:
source /project/k10066/osipovs/.commonrc; gogomamba; mamba activate py311
year=2022
data_dir=/work/mm0062/b302074/Data/AirQuality/AREAD/IC_BC/
python -u ${MERRA2BC}/merra2wrf_main.py --hourly_interval=3 --do_IC --do_BC --zero_out_first --wrf_dir=${data_dir}/debugICBC/ --wrf_met_dir=${data_dir}/met_em/ --wrf_met_files=met_em.d01.${year}-* >& log.emac2wrf

How to run NCEC 2021:
source /project/k10066/osipovs/.commonrc; gogomamba; mamba activate py311
data_dir=/work/mm0062/b302074/Data/AirQuality/NCEC/2021/IC_BC/
python -u ${MERRA2BC}/merra2wrf_main.py --hourly_interval=3 --do_IC --do_BC --zero_out_first --wrf_dir=${data_dir} --wrf_met_dir=${data_dir}/met_em/ --wrf_met_files=met_em.d01.* >& log.merra2wrf

source /project/k10066/osipovs/.commonrc; gogomamba; mamba activate py311
data_dir=/scratch/osipovs/Models/WRF_run/run_real
python -u ${MERRA2BC}/merra2wrf_main.py --hourly_interval=3 --do_IC --do_BC --zero_out_first --wrf_dir=${data_dir} --wrf_met_dir=${data_dir}/ --wrf_met_files=met_em.d01.* --wrf_met_files_date_format=%Y-%m-%d_%H:%M:%S >& log.merra2wrf
"""

#%%
parser = argparse.ArgumentParser()
parser.add_argument("--start_date", help="YYYY-MM-DD_HH:MM:SS format", default=None, type=argsparse_datetime)  # default='2017-08-31_21:00:00')  #
parser.add_argument("--end_date", help="YYYY-MM-DD_HH:MM:SS format", default=None, type=argsparse_datetime)  # default='2017-09-01_00:00:00')  #
parser.add_argument("--hourly_interval", help="dt in dates to process", default=3)
parser.add_argument("--do_IC", help="process initial conditions?", action='store_true')
parser.add_argument("--do_BC", help="process boundary conditions?", action='store_true')
parser.add_argument("--zero_out_first", help="zero out fields in IC and BC first?", action='store_true')
# setup parent model
parser.add_argument("--reanalysis_dir", default='/project/k10048/Data/NASA/MERRA2/')  # MERRA2
parser.add_argument("--reanalysis_file_name_template", default='{stream}/MERRA2_{VVV}.{stream}.{date_time}.nc4')  # {VVV} is 400, 401 etc
# setup downscaling model (WRF)
parser.add_argument("--wrf_dir", help="folder containing WRF IC & BC files")   # , default=/Data/AirQuality/AQABA/IC_BC')
parser.add_argument("--wrf_input", help="use default wrfinput_d01", default='wrfinput_d01')
parser.add_argument("--wrf_bdy_file", help="use default wrfbdy_d01", default='wrfbdy_d01')
parser.add_argument("--wrf_met_dir", help="use default wrfbdy_d01")
parser.add_argument("--wrf_met_files", help="met_em file names template")  # , default='met_em.d01.2017-0*')
parser.add_argument("--wrf_met_files_date_format", help="met_em file names template", default='%Y-%m-%d_%H_%M_%S')  # or '%Y-%m-%d_%H:%M:%S'

parser.add_argument("--mode", "--port", "--host", help="the are only to support pycharm debugging")
args = parser.parse_args()
#%% DEBUG settings
# args.wrf_dir = '/work/mm0062/b302074/Data/AirQuality/AREAD/IC_BC/debugICBC/'
# args.wrf_met_dir = '/work/mm0062/b302074/Data/AirQuality/AREAD/IC_BC/met_em/'
# args.wrf_met_files = 'met_em.d01.2022-0*'
# args.do_IC = True
# args.do_BC = True
# args.zero_out_first = False  # just speed up
#%%
wrf_bdy_fp = args.wrf_dir + "/" + args.wrf_bdy_file
wrf_ic_fp = args.wrf_dir + "/" + args.wrf_input
print('WRF bdy fp: {}'.format(wrf_bdy_fp))
print('WRF ic fp: {}'.format(wrf_ic_fp))

wrf_bdy_ds = xr.open_dataset(wrf_bdy_fp)
wrf_ic_ds = xr.open_dataset(wrf_ic_fp)

start_time = time.time()

mappings = get_merra2wrf_mapping()

print("\nConversion MAP:")
for mapping in mappings:
    print(mapping.mapping_rule_str + ":\t")
#%%
print('\nDeriving dates to process from wrf_bdy file\n')
wrfbdy_f = nc.Dataset(wrf_bdy_fp, 'r')

wrf_bdy_dates = wrf.extract_times(wrfbdy_f, wrf.ALL_TIMES)
wrf_bdy_dates = pd.to_datetime(wrf_bdy_dates)
dates_to_process = wrf_bdy_dates
delta_t = (dates_to_process[1] - dates_to_process[0]).total_seconds()

if args.start_date is not None:
    dates_to_process = dates_to_process[dates_to_process >= args.start_date]
if args.end_date is not None:
    dates_to_process = dates_to_process[dates_to_process <= args.end_date]
#%% Aux
def xr_open_merra2_dataset(fp):  # standardize file opening
    print("Opening file: {}".format(fp))
    ds = xr.open_dataset(fp)
    # ds = ds.isel(lat=slice(None, None, -1))  # flip lat dimension, to make ascending, satisfying interpolation requirements

    if 'lev' in ds.dims:
        ds = ds.rename({'lev':'level'})
        # ds = ds.isel(level=slice(None, None, -1))  # flip dimensions

    return ds


xr_open_dataset_impl = xr_open_merra2_dataset
#%% IC
if args.do_IC:
    date = dates_to_process[0]
    print("INITIAL CONDITIONS: {}".format(date))

    fp = args.reanalysis_dir + generate_merra2_file_name(args.reanalysis_file_name_template, date, stream='inst3_3d_aer_Nv')
    donor_ml_ds = xr_open_dataset_impl(fp).sel(time=date)
    # grab pressure from another stream instead of deriving it
    # derive_merra2_pressure_profile(donor_ml_ds)
    fp = args.reanalysis_dir + generate_merra2_file_name(args.reanalysis_file_name_template, date, stream='inst3_3d_asm_Nv')
    donor_asm_ds = xr_open_dataset_impl(fp).sel(time=date)
    donor_ml_ds['p_rho'] = donor_asm_ds.PL

    # alternative way to sort
    # merra_df = merra_df.sortby('level', ascending=False) # set TOA as last layer
    # pressure_stag = pressure_stag.sortby('level', ascending=False)
    # pressure_rho = pressure_rho.sortby('level', ascending=False)

    donor_p_rho_hi_da = regridding_utils.hor_interpolate_3d_field_on_wrf_grid(donor_ml_ds.p_rho, wrf_ic_ds)

    met_file_name = regridding_utils.get_met_file_by_time(date.strftime(args.wrf_met_files_date_format))
    met_fp = args.wrf_met_dir + "/" + met_file_name
    print("Opening metfile: " + met_fp)
    wrf_met_ds = xr.open_dataset(met_fp).set_coords(['XLAT_M', 'XLONG_M'])
    derive_wrf_pressure_from_met_and_input_files(wrf_met_ds, wrf_ic_ds)
    wrf_met_ds = wrf_met_ds.squeeze()  # squeeze Time, otherwise causes size issues with pure netcdf

    print("Opening wrfintput: " + wrf_ic_fp)
    wrfinput_f = nc.Dataset(wrf_ic_fp, 'r+')

    if args.zero_out_first:  # zero out the fields in wrfinput before increments
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
            fp = args.reanalysis_dir + generate_merra2_file_name(args.reanalysis_file_name_template, date, stream=mapping.output_stream)
            donor_stream_ds = xr_open_dataset_impl(fp).sel(time=date)

            donor_var_da = donor_stream_ds[rule_vo['merra_key']].load()

            print("\t\t - Horizontal interpolation of " + rule_vo['merra_key'] + " on WRF horizontal grid")
            donor_var_hi_da = regridding_utils.hor_interpolate_3d_field_on_wrf_grid(donor_var_da, wrf_ic_ds)

            print("\t\t - Vertical interpolation of " + rule_vo['merra_key'] + " on WRF vertical grid")
            donor_var_hvi_da = regridding_utils.ver_interpolate_3d_field_on_wrf_grid(donor_var_hi_da, donor_p_rho_hi_da, wrf_met_ds.p_rho)

            wrf_key = rule_vo['wrf_key']
            print("\t\t - Updating wrfinput field {}[0] += {} * {} * {:.1e}".format(wrf_key, rule_vo['merra_key'], rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
            increment_da = donor_var_hvi_da * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent']
            wrfinput_f.variables[wrf_key][0, :] = wrfinput_f.variables[wrf_key][0, :] + increment_da.transpose('bottom_top', 'south_north', 'west_east')

            donor_stream_ds.close()

    wrfinput_f.close()
    donor_ml_ds.close()
    print("DONE: INITIAL CONDITIONS")

#%% BC
if args.do_BC:
    print("\n\nSTART BOUNDARY CONDITIONS")

    print("Opening WRF bdy: " + wrf_bdy_fp)
    wrfbdy_f = nc.Dataset(wrf_bdy_fp, 'r+')

    for date in dates_to_process:
        print("\n\tBCs, next date: {}".format(date))

        fp = args.reanalysis_dir + generate_merra2_file_name(args.reanalysis_file_name_template, date, stream='inst3_3d_aer_Nv')
        donor_ml_ds = xr_open_dataset_impl(fp).sel(time=date)
        # read pressure instead of deriving it
        # derive_merra2_pressure_profile(donor_ml_ds)
        fp = args.reanalysis_dir + generate_merra2_file_name(args.reanalysis_file_name_template, date, stream='inst3_3d_asm_Nv')
        donor_asm_ds = xr_open_dataset_impl(fp).sel(time=date)
        donor_ml_ds['p_rho'] = donor_asm_ds.PL

        met_file_name = regridding_utils.get_met_file_by_time(date.strftime(args.wrf_met_files_date_format))
        met_fp = args.wrf_met_dir + "/" + met_file_name
        print("\tReading WRF Pressure from: {}".format(met_fp))
        wrf_met_ds = xr.open_dataset(met_fp).set_coords(['XLAT_M', 'XLONG_M'])
        derive_wrf_pressure_from_met_and_input_files(wrf_met_ds, wrf_ic_ds)
        wrf_met_ds = wrf_met_ds.squeeze()  # squeeze Time, otherwise causes size issues with pure netcdf

        # prep boundary arrays for interpolation
        wrf_met_ds['p_rho_on_bnd'] = flatten_and_line_up_boundaries(wrf_met_ds.p_rho)  # concatenate east, west, nort and south boundaries
        wrf_met_ds['boundary_lons'] = flatten_and_line_up_boundaries(wrf_met_ds.XLONG_M)
        wrf_met_ds['boundary_lats'] = flatten_and_line_up_boundaries(wrf_met_ds.XLAT_M)

        print("\tHorizontal interpolation of Pressure on WRF boundary")
        donor_p_rho_hi_da = regridding_utils.hor_interpolate_3d_field_on_wrf_boundary(donor_ml_ds.p_rho, wrf_met_ds)

        time_index_in_wrfbdy = wrf_bdy_dates.get_loc(date)

        if args.zero_out_first:  # zero out the fields before increments
            unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
            print("Following fields will be zeroed out FIRST in BCs: {}".format(unique_wrf_keys))
            for wrf_key in unique_wrf_keys:
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

            fp = args.reanalysis_dir + generate_merra2_file_name(args.reanalysis_file_name_template, date, mapping.output_stream)
            donor_stream_ds = xr_open_dataset_impl(fp).sel(time=date)

            for rule_vo in pipe_to_process:
                print("\t\t - Reading " + rule_vo['merra_key'])

                merra_key = rule_vo['merra_key']
                wrf_key = rule_vo['wrf_key']

                donor_var_da = donor_stream_ds[rule_vo['merra_key']].load()
                print("\t\tHorizontal interpolation of " + merra_key + " on WRF boundary")
                donor_var_hi_da = regridding_utils.hor_interpolate_3d_field_on_wrf_boundary(donor_var_da, wrf_met_ds)
                print("\t\tVertical interpolation of " + merra_key + " on WRF boundary")
                donor_var_hvi_da = regridding_utils.ver_interpolate_3d_field_on_wrf_boundary(donor_var_hi_da, donor_p_rho_hi_da, wrf_met_ds.p_rho_on_bnd)

                print("\t\t - Updating wrfbdy: {}[{}] += {} * {} * {:.1e}".format(wrf_key, time_index_in_wrfbdy, merra_key, rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
                increment_on_bnd_da = donor_var_hvi_da * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent']
                regridding_utils.update_boundaries(increment_on_bnd_da, wrf_met_ds, wrfbdy_f, wrf_key, time_index_in_wrfbdy)

            donor_stream_ds.close()

        unique_wrf_keys = get_unique_wrf_keys_from_mappings(mappings)
        print('Prescribing uniform tendency in BC for these fields: {}'.format(unique_wrf_keys))
        for wrf_key in unique_wrf_keys:
            regridding_utils.update_tendency_boundaries(wrfbdy_f, wrf_key, time_index_in_wrfbdy, delta_t)

        print("--- %s seconds ---" % (time.time() - start_time))
        donor_ml_ds.close()

    print("Closing " + args.wrf_bdy_file)
    wrfbdy_f.close()

    print("FINISH BOUNDARY CONDITIONS")

#%% time stats
print("--- %s seconds ---" % (time.time() - start_time))
