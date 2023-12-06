from datetime import datetime
import re
import os
from netCDF4 import Dataset
import numpy as np

met_files = []
met_times_files = {}
dates = {}

dx = dy = true_lat1 = true_lat2 = 0
cen_lat = 0
cen_lon = 0
projection = ""

nx = ny = nz = nw = 0
wrf_p_top = 0
znu = []
xlon = [[]]
xlat = [[]]

wrf_vars = []

# mapping for basemap projections
projections = {"Lambert Conformal": "lcc", "Mercator": "merc"}

boundary_lons = []
boundary_lats = []


def get_pressure_from_metfile(metfile):
    PSFC = metfile.variables['PSFC'][:]
    WRF_Pres = np.zeros([nz, ny, nx])
    for z_level in reversed(range(nz)):
        WRF_Pres[nz - 1 - z_level, :] = PSFC * znu[0, z_level] + (1.0 - znu[0, z_level]) * wrf_p_top
    return WRF_Pres


def get_met_file_by_time(time):
    return "met_em.d01." + time + ".nc"


def get_index_in_file_by_time(time):
    return dates.get(time)


def get_BaseMapProjectionByWrfProjection():
    return projections.get(projection)


numbers = re.compile(r'(\d+)')


def numericalSort1(value):
    parts = numbers.split(value)
    return int(float(parts[3]) * 1e6 + float(parts[5]) * 1e4 + float(parts[7]) * 1e2 + float(parts[9]))


def get_ordered_met_files():
    return met_files


def update_boundaries(WRF_SPECIE_BND, wrfbdy_f, name, index):
    WRF_SPECIE_LEFT_BND = WRF_SPECIE_BND[:, 0:ny]
    WRF_SPECIE_TOP_BND = WRF_SPECIE_BND[:, ny:ny + nx]
    WRF_SPECIE_RIGHT_BND = WRF_SPECIE_BND[:, ny + nx:2 * ny + nx]
    WRF_SPECIE_BOT_BND = WRF_SPECIE_BND[:, 2 * ny + nx:2 * ny + 2 * nx]

    wrfbxs = np.repeat(WRF_SPECIE_LEFT_BND[np.newaxis, :, :], nw, axis=0)
    wrfbxe = np.repeat(WRF_SPECIE_RIGHT_BND[np.newaxis, :, :], nw, axis=0)
    wrfbys = np.repeat(WRF_SPECIE_BOT_BND[np.newaxis, :, :], nw, axis=0)
    wrfbye = np.repeat(WRF_SPECIE_TOP_BND[np.newaxis, :, :], nw, axis=0)

    # print("\t\t\tUpdating BC for "+name)
    wrfbdy_f.variables[name + "_BXS"][index, :] = wrfbdy_f.variables[name + "_BXS"][index, :] + wrfbxs
    wrfbdy_f.variables[name + "_BXE"][index, :] = wrfbdy_f.variables[name + "_BXE"][index, :] + wrfbxe
    wrfbdy_f.variables[name + "_BYS"][index, :] = wrfbdy_f.variables[name + "_BYS"][index, :] + wrfbys
    wrfbdy_f.variables[name + "_BYE"][index, :] = wrfbdy_f.variables[name + "_BYE"][index, :] + wrfbye


def update_tendency_boundaries(wrfbdy_f, name, index, dt):
    if index > 0:
        print("\t\t\tUpdating Tendency BC for " + name)
        wrfbdy_f.variables[name + "_BTXS"][index - 1, :] = (wrfbdy_f.variables[name + "_BXS"][index, :] -
                                                            wrfbdy_f.variables[name + "_BXS"][index - 1, :]) / dt
        wrfbdy_f.variables[name + "_BTXE"][index - 1, :] = (wrfbdy_f.variables[name + "_BXE"][index, :] -
                                                            wrfbdy_f.variables[name + "_BXE"][index - 1, :]) / dt
        wrfbdy_f.variables[name + "_BTYS"][index - 1, :] = (wrfbdy_f.variables[name + "_BYS"][index, :] -
                                                            wrfbdy_f.variables[name + "_BYS"][index - 1, :]) / dt
        wrfbdy_f.variables[name + "_BTYE"][index - 1, :] = (wrfbdy_f.variables[name + "_BYE"][index, :] -
                                                            wrfbdy_f.variables[name + "_BYE"][index - 1, :]) / dt


def initialise(config):
    global met_files, dates, wrf_p_top, znu, xlon, xlat, nx, ny, nz, nw, boundary_lons, boundary_lats, vars, cen_lat, cen_lon, projection, dx, dy, true_lat2, true_lat1

    met_files = sorted([f for f in os.listdir(config.wrf_met_dir) if re.match(config.wrf_met_files, f)], key=numericalSort1)
    wrf_bdy_fp = config.wrf_dir + "/" + config.wrf_bdy_file
    print('wrf_bdy_fp is {}'.format(wrf_bdy_fp))
    wrf_bdy_nc = Dataset(wrf_bdy_fp, 'r')
    for i in range(0, len(wrf_bdy_nc.variables['Times'][:]), 1):
        # wrf_times.update({''.join(wrfbddy.variables['Times'][i]):i})
        date_string = ''.join([char.decode("utf-8") for char in wrf_bdy_nc.variables['Times'][i]])
        dates.update({date_string: i})
        met_times_files.update({date_string: met_files[i]})

    nx = len(wrf_bdy_nc.dimensions['west_east'])
    ny = len(wrf_bdy_nc.dimensions['south_north'])
    nz = len(wrf_bdy_nc.dimensions['bottom_top'])
    nw = len(wrf_bdy_nc.dimensions['bdy_width'])
    wrf_bdy_nc.close()

    # Reading "PRESSURE TOP OF THE MODEL, PA" and "eta values on half (mass) levels"
    wrfinput_nc = Dataset(config.wrf_dir + "/" + config.wrf_input, 'r')
    wrf_p_top = wrfinput_nc.variables['P_TOP'][:]
    znu = wrfinput_nc.variables['ZNU'][:]
    xlon = wrfinput_nc.variables['XLONG'][0]
    xlat = wrfinput_nc.variables['XLAT'][0]
    vars = [var for var in wrfinput_nc.variables]

    projection = wrfinput_nc.getncattr('MAP_PROJ_CHAR')
    cen_lat = wrfinput_nc.getncattr('CEN_LAT')
    cen_lon = wrfinput_nc.getncattr('CEN_LON')
    dy = wrfinput_nc.getncattr('DY')
    dx = wrfinput_nc.getncattr('DX')

    true_lat1 = wrfinput_nc.getncattr('TRUELAT1')
    true_lat2 = wrfinput_nc.getncattr('TRUELAT2')

    wrfinput_nc.close()

    boundary_lons = np.concatenate((xlon[:, 0], xlon[ny - 1, :], xlon[:, nx - 1], xlon[0, :]), axis=0)
    boundary_lats = np.concatenate((xlat[:, 0], xlat[ny - 1, :], xlat[:, nx - 1], xlat[0, :]), axis=0)

    print("\nWRF dimensions: [bottom_top]=" + str(nz) + " [south_north]=" + str(ny) + " [west_east]=" + str(
        nx) + " [bdy_width]=" + str(nw))
    print("P_TOP=" + str(wrf_p_top) + " Pa")

    print("Lower left corner: lat=" + str(min(boundary_lats)) + " long=" + str(min(boundary_lons)))
    print("Upper right corner: lat=" + str(max(boundary_lats)) + " long=" + str(max(boundary_lons)))