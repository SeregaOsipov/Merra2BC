# It is OK to have duplicates in spc_map

import config
import time
import merra2_module
import wrf_module
import merra2wrf_mapper
import utils
from netCDF4 import Dataset
import numpy as np
from datetime import datetime

"""
Steps:
1. Run zero_fields.py first to zero out the IC and BC for specified chem species
2. Run this file, which will increment the IC and BC by the values obtained from the MERRA2
3. Change config and run this file again to include SO2 and sulfate
"""

start_time = time.time()

# modules initialisation
wrf_module.initialise()
merra2wrf_mapper.initialise()
merra2_module.initialise()

# -----------------------------
# Sanity checks:
# check species availability in wrf and in merra files
for rule_vo in merra2wrf_mapper.pipe_to_process:
    if rule_vo['merra_key'] not in merra2_module.merra_vars:
        utils.error_message("Could not find variable " + rule_vo + " in MERRA2 file. Exiting...")

for rule_vo in merra2wrf_mapper.pipe_to_process:
    if rule_vo['wrf_key'] not in wrf_module.wrf_vars:
        utils.error_message("Could not find variable " + rule_vo + " in WRF input file. Exiting...")

# check that spatial dimensions are covered
if ((min(wrf_module.wrf_bnd_lons) < min(merra2_module.mera_lon)) | (
        max(wrf_module.wrf_bnd_lons) > max(merra2_module.mera_lon)) | (
        min(wrf_module.wrf_bnd_lats) < min(merra2_module.mera_lat)) | (
        max(wrf_module.wrf_bnd_lats) > max(merra2_module.mera_lat))):
    utils.error_message("WRF area is not fully covered by MERRA2 area. Exiting...")

time_intersection = wrf_module.wrf_times.keys() & merra2_module.mera_times.keys()

# check that merra2 time is covered by wrf time
if len(time_intersection) != len(wrf_module.wrf_times):
    print('These datetimes are missing in MERRA2 dataset:')
    time_intersection = dict.fromkeys(time_intersection, 0)
    for key in wrf_module.wrf_times.keys():
        if not key in time_intersection:
            print(key)
    utils.error_message("WRF time range is not fully covered by MERRA2 time range. Exiting...")

# -----------------------------
# sorting times for processing
time_intersection = sorted(time_intersection, key=lambda x: time.mktime(time.strptime(x, "%Y-%m-%d_%H:%M:%S")))

print("\nTimes for processing: {}".format(time_intersection))
# exit()

if config.do_IC:
    print("START INITIAL CONDITIONS")
    cur_time = time_intersection[0]
    index_of_opened_merra_file = merra2_module.get_file_index_by_time(cur_time)
    print("Opening merra file: " + merra2_module.get_file_name_by_index(
        index_of_opened_merra_file) + " with initial time: " + cur_time)
    print(config.mera_dir + "/" + merra2_module.get_file_name_by_index(index_of_opened_merra_file))
    merra_f = Dataset(config.mera_dir + "/" + merra2_module.get_file_name_by_index(index_of_opened_merra_file), 'r')
    MERRA_PRES = merra2_module.get_pressure_by_time(cur_time, merra_f)

    # Horizontal interpolation of Merra pressure on WRF horizontal grid
    MER_HOR_PRES = merra2_module.hor_interpolate_3dfield_on_wrf_grid(MERRA_PRES, wrf_module.ny, wrf_module.nx,
                                                                     wrf_module.xlon, wrf_module.xlat)

    print("Opening metfile: " + wrf_module.get_met_file_by_time(cur_time))
    metfile = Dataset(config.wrf_met_dir + "/" + wrf_module.get_met_file_by_time(cur_time), 'r')
    WRF_PRES = wrf_module.get_pressure_from_metfile(metfile)

    print("Opening wrfintput: " + config.wrf_input_file)
    wrfinput_f = Dataset(config.wrf_dir + "/" + config.wrf_input_file, 'r+')

    for rule_vo in merra2wrf_mapper.pipe_to_process:
        print("\t\t - Reading " + rule_vo['merra_key'] + " field from Merra2.")
        MER_SPECIE = merra2_module.get_3dfield_by_time(cur_time, merra_f, rule_vo['merra_key'])

        print("\t\t - Horisontal interpolation of " + rule_vo['merra_key'] + " on WRF horizontal grid")
        MER_HOR_SPECIE = merra2_module.hor_interpolate_3dfield_on_wrf_grid(MER_SPECIE, wrf_module.ny, wrf_module.nx,
                                                                           wrf_module.xlon, wrf_module.xlat)

        print("\t\t - Vertical interpolation of " + rule_vo['merra_key'] + " on WRF vertical grid")
        WRF_SPECIE = merra2_module.ver_interpolate_3dfield_on_wrf_grid(MER_HOR_SPECIE, MER_HOR_PRES, WRF_PRES,
                                                                       wrf_module.nz, wrf_module.ny, wrf_module.nx)
        WRF_SPECIE = np.flipud(WRF_SPECIE)

        wrf_key = rule_vo['wrf_key']
        print("\t\t - Updating wrfinput field {}[0] += {} * {} * {:.1e}".format(wrf_key, rule_vo['merra_key'], rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
        wrfinput_f.variables[wrf_key][0, :] = wrfinput_f.variables[wrf_key][0, :] + WRF_SPECIE * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent']

    print("Closing wrfintput: " + config.wrf_input_file)
    wrfinput_f.close()

    print("Closing mera file " + merra2_module.get_file_name_by_index(index_of_opened_merra_file))
    merra_f.close()

    print("Closing metfile " + wrf_module.get_met_file_by_time(cur_time))
    metfile.close()

    print("FINISH INITIAL CONDITIONS")

if config.do_BC:
    print("\n\nSTART BOUNDARY CONDITIONS")

    print("Opening " + config.wrf_bdy_file)
    wrfbdy_f = Dataset(config.wrf_dir + "/" + config.wrf_bdy_file, 'r+')

    # difference betweeen two given times
    dt = (datetime.strptime(time_intersection[1], '%Y-%m-%d_%H:%M:%S') - datetime.strptime(time_intersection[0],
                                                                                           '%Y-%m-%d_%H:%M:%S')).total_seconds()

    cur_time = time_intersection[0]
    index_of_opened_merra_file = merra2_module.get_file_index_by_time(cur_time)
    print("\nOpening MERRA2 file: " + merra2_module.get_file_name_by_index(
        index_of_opened_merra_file) + " file which has index " + str(index_of_opened_merra_file))
    merra_f = Dataset(config.mera_dir + "/" + merra2_module.get_file_name_by_index(index_of_opened_merra_file), 'r')

    for cur_time in time_intersection:
        if merra2_module.get_file_index_by_time(cur_time) != index_of_opened_merra_file:
            print("Closing prev. opened MERRA2 file with index " + str(index_of_opened_merra_file))
            merra_f.close()

            index_of_opened_merra_file = merra2_module.get_file_index_by_time(cur_time)
            print("\nOpening MERRA2 file: " + merra2_module.get_file_name_by_index(
                index_of_opened_merra_file) + " file which has index " + str(index_of_opened_merra_file))
            merra_f = Dataset(config.mera_dir + "/" + merra2_module.get_file_name_by_index(index_of_opened_merra_file),
                              'r')

        print("\n\tCur_time=" + cur_time)
        print("\tReading MERRA Pressure at index " + str(merra2_module.get_index_in_file_by_time(cur_time)))
        MERRA_PRES = merra2_module.get_pressure_by_time(cur_time, merra_f)

        print("\tHorizontal interpolation of MERRA Pressure on WRF boundary")
        MER_HOR_PRES_BND = merra2_module.hor_interpolate_3dfield_on_wrf_boubdary(MERRA_PRES,
                                                                                 len(wrf_module.wrf_bnd_lons),
                                                                                 wrf_module.wrf_bnd_lons,
                                                                                 wrf_module.wrf_bnd_lats)

        print("\tReading WRF Pressure from: " + wrf_module.get_met_file_by_time(cur_time))
        metfile = Dataset(config.wrf_met_dir + "/" + wrf_module.get_met_file_by_time(cur_time), 'r')
        WRF_PRES = wrf_module.get_pressure_from_metfile(metfile)
        WRF_PRES_BND = np.concatenate((WRF_PRES[:, :, 0], WRF_PRES[:, wrf_module.ny - 1, :],
                                       WRF_PRES[:, :, wrf_module.nx - 1], WRF_PRES[:, 0, :]), axis=1)
        metfile.close()

        time_index = wrf_module.get_index_in_file_by_time(cur_time)
        for rule_vo in merra2wrf_mapper.pipe_to_process:
            merra_key = rule_vo['merra_key']
            wrf_key = rule_vo['wrf_key']
            print("\n\t\t - Reading " + merra_key + " field from MERRA.")
            MER_SPECIE = merra2_module.get_3dfield_by_time(cur_time, merra_f, merra_key)

            print("\t\tHorizontal interpolation of " + merra_key + " on WRF boundary")
            MER_HOR_SPECIE_BND = merra2_module.hor_interpolate_3dfield_on_wrf_boubdary(MER_SPECIE,
                                                                                       len(wrf_module.wrf_bnd_lons),
                                                                                       wrf_module.wrf_bnd_lons,
                                                                                       wrf_module.wrf_bnd_lats)

            print("\t\tVertical interpolation of " + merra_key + " on WRF boundary")
            WRF_SPECIE_BND = merra2_module.ver_interpolate_3dfield_on_wrf_boubdary(MER_HOR_SPECIE_BND, MER_HOR_PRES_BND,
                                                                                   WRF_PRES_BND, wrf_module.nz,
                                                                                   len(wrf_module.wrf_bnd_lons))
            WRF_SPECIE_BND = np.flipud(WRF_SPECIE_BND)

            print("\t\t - Updating wrfbdy field: {}[{}] += {} * {} * {:.1e}".format(wrf_key, time_index, merra_key, rule_vo['wrf_multiplier'], rule_vo['wrf_exponent']))
            wrf_module.update_boundaries(WRF_SPECIE_BND * rule_vo['wrf_multiplier'] * rule_vo['wrf_exponent'], wrfbdy_f, wrf_key, time_index)

        unique_wrf_keys = []
        for vo in merra2wrf_mapper.pipe_to_process:
            if vo['wrf_key'] not in unique_wrf_keys:
                unique_wrf_keys.append(vo['wrf_key'])

        wrf_sp_index = 0
        for wrf_key in unique_wrf_keys:
            wrf_module.update_tendency_boundaries(wrfbdy_f, wrf_key, time_index, dt, wrf_sp_index)
            wrf_sp_index = wrf_sp_index + 1

        print("--- %s seconds ---" % (time.time() - start_time))

    print("Closing prev. opened MERRA2 file with index " + str(index_of_opened_merra_file))
    merra_f.close()

    print("Closing " + config.wrf_bdy_file)
    wrfbdy_f.close()

    print("FINISH BOUNDARY CONDITIONS")

print("--- %s seconds ---" % (time.time() - start_time))
