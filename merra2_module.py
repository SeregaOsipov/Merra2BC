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
mera_times = {}  # map between time and index in file
mera_times_files = {}  # map between time and file index
vars = []
n_x_points = 0
n_y_points = 0
n_z_points = 0

numbers = re.compile(r'(\d+)')


def numericalSort(value):
    parts = numbers.split(value)
    return parts[9]


def get_file_index_by_time(time):
    return mera_times_files.get(time)


def get_index_in_file_by_time(time):
    return mera_times.get(time)


def get_file_name_by_index(index):
    return files[index]


# ********************************
# Horizontal interpolation of 3d Merra field on WRF boundary
def hor_interpolate_3dfield_on_wrf_boubdary(field_df, wrf_length, wrf_lon, wrf_lat):
    FIELD_BND = np.zeros([field_df.sizes['lev'], wrf_length])
    for z_level in range(field_df.sizes['lev']):
        f = interpolate.RectBivariateSpline(field_df.lat, field_df.lon, field_df.isel(lev=z_level))
        FIELD_BND[z_level, :] = f(wrf_lat, wrf_lon, grid=False)
    return FIELD_BND


# Vertical interpolation of Merra boundary on WRF boundary
def ver_interpolate_3dfield_on_wrf_boubdary(var_df, merra_p_df, WRF_PRES_BND, wrf_nz, wrf_length):
    WRF_SPECIE_BND = np.zeros([wrf_nz, wrf_length])  # Required SPEC on WRF boundary
    for i in range(0, wrf_length):
        f = interpolate.interp1d(merra_p_df[:, i], var_df[:, i], kind='linear', bounds_error=False, fill_value=0)
        WRF_SPECIE_BND[:, i] = f(WRF_PRES_BND[:, i])
    return WRF_SPECIE_BND


# Horizontal interpolation of 3d Merra field on WRF horizontal grid
def hor_interpolate_3d_field_on_wrf_grid(field_df, wrf_ny, wrf_nx, wrf_lon, wrf_lat):
    FIELD_HOR = np.zeros([field_df.sizes['lev'], wrf_ny, wrf_nx])

    for z_level in range(field_df.sizes['lev']):
        f = interpolate.RectBivariateSpline(field_df.lat, field_df.lon, field_df.isel(lev=z_level))
        FIELD_HOR[z_level, :, :] = f(wrf_lat, wrf_lon, grid=False).reshape(wrf_ny, wrf_nx)

    return FIELD_HOR


# Vertical interpolation on WRF grid
def ver_interpolate_3d_field_on_wrf_grid(MER_HOR_SPECIE, MER_HOR_PRES, WRF_PRES, wrf_nz, wrf_ny, wrf_nx):
    WRF_SPECIE = np.zeros([wrf_nz, wrf_ny, wrf_nx])  # Required SPEC on WRF grid
    for x in range(0, wrf_nx, 1):
        for y in range(0, wrf_ny, 1):
            f = interpolate.interp1d(MER_HOR_PRES[:, y, x], MER_HOR_SPECIE[:, y, x], kind='linear', bounds_error=False,
                                     fill_value=0)
            WRF_SPECIE[:, y, x] = f(WRF_PRES[:, y, x])
    return WRF_SPECIE


# ********************************

# extracts 3d field from merra2 file from given time
def get_3dfield_by_time(time, merra_file, field_name):
    mera_time_idx = get_index_in_file_by_time(time)
    field = merra_file.variables[field_name][mera_time_idx, :]

    if shifted_lons:
        field = np.roll(field, shift_index, axis=2)

    return np.flipud(field)


def get_pressure_by_time(time, merra_file):
    # global p_top
    # MER_Pres will be restored on 73 edges
    MER_Pres = np.zeros([n_z_points + 1, n_y_points, n_x_points])
    # filling top edge with Ptop_mera
    MER_Pres[0, :, :] = p_top

    # Extract deltaP from NetCDF file at index defined by time
    mera_time_idx = get_index_in_file_by_time(time)
    DELP = merra_file.variables['DELP'][mera_time_idx, :]  # Pa

    for z_level in range(n_z_points):
        MER_Pres[z_level + 1] = MER_Pres[z_level] + DELP[z_level]

    # BUT! we need pressure on 72 levels
    # => averaging pressure values on adjacent edges
    MER_Pres = (MER_Pres[0:n_z_points:1][:, :] + MER_Pres[1::1][:, :]) / 2

    if shifted_lons:
        MER_Pres = np.roll(MER_Pres, shift_index, axis=2)

    MER_Pres = np.flipud(MER_Pres)
    return MER_Pres


def initialise(config, mappings):
    # global files, n_x_points, n_y_points, n_z_points, lon, lat, vars, shifted_lons, shift_index

    # get the aux info from any mapping rule
    mapping = mappings[0]
    fp = config.merra2_dir + config.merra2_file_name_template.format(date_time='\w+', stream=mapping.output_stream)
    file_name = os.path.basename(fp)
    folder_path = os.path.dirname(fp)
    files = sorted([f for f in os.listdir(folder_path) if re.match(file_name, f)], key=numericalSort)

    merra_f = Dataset('{}/{}'.format(folder_path, files[0]), 'r')
    n_x_points = merra_f.variables['lon'].size
    n_y_points = merra_f.variables['lat'].size
    try:  # not all merra2 files (loading diagnostic) have 'lev' variable
        n_z_points = merra_f.variables['lev'].size
    except Exception:
        pass

    print("MERRA2 dimensions: [bottom_top]=" + str(n_z_points) + " [south_north]=" + str(
        n_y_points) + " [west_east]=" + str(n_x_points))

    vars = [var for var in merra_f.variables]

    lon = merra_f.variables['lon'][:]
    lat = merra_f.variables['lat'][:]

    # if data is given in range of 0_360, then we need to shift lons and data to the -180_180
    if max(lon) > 180:
        print("###########################")
        print("ATTENTION!!!:")
        print("SHIFTING LONGITUDES")
        index = 0
        for lon in lon:
            if lon > 180:
                lon[index] = lon[index] - 360.0
            index = index + 1
        shift_index = len(lon) / 2
        lon = np.roll(lon, shift_index)
        shifted_lons = True
        print("###########################")

    print("Lower left corner: lat=" + str(min(lat)) + " long=" + str(min(lon)))
    print("Upper right corner: lat=" + str(max(lat)) + " long=" + str(max(lon)))

    # number of times in  mera file
    times_per_file = merra_f.variables['time'].size
    merra_f.close()

    index = 0
    for merra_file in files:
        date = numbers.split(merra_file)[9]
        for i in range(0, times_per_file, 1):
            t = datetime.strptime(date, '%Y%m%d') + timedelta(minutes=(i * (24 / times_per_file) * 60))
            mera_times_files.update({t.strftime("%Y-%m-%d_%H:%M:%S"): index})
            mera_times.update({t.strftime("%Y-%m-%d_%H:%M:%S"): i})
        index = index + 1
