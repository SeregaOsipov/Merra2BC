import merra2wrf_config
import re
import os
from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
from scipy import interpolate

from multiprocessing import Pool
from functools import partial

p_top = 1  # Pa   (=0.01 hPa)
lat = 0
lon = 0

shifted_lons = False
shift_index = 0

files = []
dates = {}  # map between time and index in file
mera_times_files = {}  # map between time and file index
vars = []
n_x_points = 0
n_y_points = 0
n_z = 0

numbers = re.compile(r'(\d+)')


# ********************************
# Horizontal interpolation of 3d Merra field on WRF boundary
def hor_interpolate_3dfield_on_wrf_boubdary(FIELD, wrf_length, wrf_lon, wrf_lat):
    FIELD_BND = np.zeros([n_z, wrf_length])
    for z_level in range(n_z):
        f = interpolate.RectBivariateSpline(np.flipud(lat), lon, np.flipud(FIELD[z_level, :, :]))  # flip to produce monotonic increasing lat
        FIELD_BND[z_level, :] = f(wrf_lat, wrf_lon, grid=False)
    return FIELD_BND


# Vertical interpolation of Merra boundary on WRF boundary
def ver_interpolate_3dfield_on_wrf_boubdary(MER_HOR_SPECIE_BND, MER_HOR_PRES_BND, WRF_PRES_BND, wrf_nz, wrf_length):
    WRF_SPECIE_BND = np.zeros([wrf_nz, wrf_length])  # Required SPEC on WRF boundary
    for i in range(0, wrf_length):
        f = interpolate.interp1d(MER_HOR_PRES_BND[:, i], MER_HOR_SPECIE_BND[:, i], kind='linear', bounds_error=False,
                                 fill_value=0)
        WRF_SPECIE_BND[:, i] = f(WRF_PRES_BND[:, i])
    return WRF_SPECIE_BND


def hor_interpolate_3d_field_on_wrf_grid(FIELD, wrf_ny, wrf_nx, wrf_lon, wrf_lat):
    FIELD_HOR = np.zeros([n_z, wrf_ny, wrf_nx])
    for z_level in range(n_z):
        # interpolation requires monotonically increasing values, have to flip lat dimension
        f = interpolate.RectBivariateSpline(np.flipud(lat), lon, np.flipud(FIELD[z_level, :, :]))
        FIELD_HOR[z_level, :, :] = f(wrf_lat, wrf_lon, grid=False).reshape(wrf_ny, wrf_nx)
    return FIELD_HOR


# Vertical interpolation on WRF grid
def ver_interpolate_3dfield_on_wrf_grid(MER_HOR_SPECIE, MER_HOR_PRES, WRF_PRES, wrf_nz, wrf_ny, wrf_nx):
    WRF_SPECIE = np.zeros([wrf_nz, wrf_ny, wrf_nx])  # Required SPEC on WRF grid
    for x in range(0, wrf_nx, 1):
        for y in range(0, wrf_ny, 1):
            f = interpolate.interp1d(MER_HOR_PRES[:, y, x], MER_HOR_SPECIE[:, y, x], kind='linear', bounds_error=False, fill_value=0)
            WRF_SPECIE[:, y, x] = f(WRF_PRES[:, y, x])
    return WRF_SPECIE


# def get_3d_field_by_time(date, emac_nc, var_key):
#     time_index = date.day - 1
#     field = emac_nc.variables[var_key][time_index]
#     if shifted_lons:
#         field = np.roll(field, shift_index, axis=2)
#     return field  # np.flipud


def get_3d_field_by_time_index(time_index, emac_nc, var_key):
    field = emac_nc.variables[var_key][time_index]
    if shifted_lons:
        field = np.roll(field, shift_index, axis=2)
    return field  # np.flipud


# def get_pressure_by_time(date, emac_nc):
#
#     # global p_top
#     # MER_Pres will be restored on 73 edges
#     MER_Pres = np.zeros([n_z + 1, n_y_points, n_x_points])
#     # filling top edge with Ptop_mera
#     MER_Pres[0, :, :] = p_top
#
#     # Extract deltaP from NetCDF file at index defined by time
#     mera_time_idx = get_index_in_file_by_time(date)
#     DELP = emac_nc.variables['DELP'][mera_time_idx, :]  # Pa
#
#     for z_level in range(n_z):
#         MER_Pres[z_level + 1] = MER_Pres[z_level] + DELP[z_level]
#
#     # BUT! we need pressure on 72 levels
#     # => averaging pressure values on adjacent edges
#     MER_Pres = (MER_Pres[0:n_z:1][:, :] + MER_Pres[1::1][:, :]) / 2
#
#     if shifted_lons:
#         MER_Pres = np.roll(MER_Pres, shift_index, axis=2)
#
#     MER_Pres = np.flipud(MER_Pres)
#     return pressure_rho


def initialise(config, mappings):  # TODO: get rid of this
    global files, n_x, n_y, n_z, lon, lat, vars, shifted_lons, shift_index

    # get the aux info from any mapping rule
    mapping = mappings[0]

    file_name = config.emac_file_name_template.format(date_time='\w+', stream=mapping.output_stream)  # only works for each mapping separatelly
    file_names = [f for f in os.listdir(config.emac_dir) if re.match(file_name, f) is not None]
    # files = sorted([f for f in os.listdir(config.dir) if re.match(mapping.files, f) is not None], key=numericalSort)
    # print("Open "+config.dir+"/"+files[0])
    nc = Dataset(config.emac_dir + "/" + file_names[0], 'r')
    n_x = nc.variables['lon'].size
    n_y = nc.variables['lat'].size
    try:
        n_z = nc.variables['lev'].size
    except Exception:
        pass

    print("MERRA2 dimensions: [bottom_top]=" + str(n_z) + " [south_north]=" + str(n_y) + " [west_east]=" + str(n_x))

    vars = [var for var in nc.variables]

    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]

    # if data is given in range of 0_360, then we need to shift lons and data to the -180_180
    if max(lon) > 180:
        print("SHIFTING LONGITUDES")
        index = 0
        for l_item in lon:
            if l_item > 180:
                lon[index] = lon[index] - 360.0
            index = index + 1
        shift_index = int(len(lon) / 2) - 1
        lon = np.roll(lon, shift_index)
        shifted_lons = True
        print("###########################")

    print("Lower left corner: lat=" + str(min(lat)) + " long=" + str(min(lon)))
    print("Upper right corner: lat=" + str(max(lat)) + " long=" + str(max(lon)))

    # number of times in mera file
    times_per_file = nc.variables['time'].size
    nc.close()

    index = 0
    for merra_file in files:
        date = numbers.split(merra_file)[5]
        for i in range(0, times_per_file, 1):
            # t = datetime.strptime(date, '%Y%m%d') + timedelta(minutes=(i * (24 / times_per_file) * 60))
            # mera_times_files.update({t.strftime("%Y-%m-%d_%H:%M:%S"): index})
            # times.update({t.strftime("%Y-%m-%d_%H:%M:%S"): i})
            t = datetime.strptime(date, '%Y%m') + timedelta(minutes=(i * (24 / times_per_file) * 60))
            mera_times_files.update({t.strftime("%Y-%m-%d_%H:%M:%S"): index})
            dates.update({t.strftime("%Y-%m-%d_%H:%M:%S"): i})
        index = index + 1
