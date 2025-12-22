import xarray as xr
import re
import numpy as np
from scipy import interpolate
import netCDF4 as nc
from types import SimpleNamespace
import itertools


def hor_interpolate_3d_field_on_wrf_boundary(donor_var_da, wrf_met_ds):
    field_bnd_da = xr.DataArray(  # wrap into xarray structure
        data=np.zeros([donor_var_da.level.size, wrf_met_ds.boundary_points.size]),
        dims=["level", "boundary_points"],
        coords=dict(level=donor_var_da.level, ),
    )

    for level_index, level in enumerate(donor_var_da.level):
        f = interpolate.RectBivariateSpline(donor_var_da.lat, donor_var_da.lon, donor_var_da.isel(level=level_index))
        field_bnd_da.isel(level=level_index).values[:] = f(wrf_met_ds.boundary_lats, wrf_met_ds.boundary_lons, grid=False)

    return field_bnd_da


def ver_interpolate_3d_field_on_wrf_boundary(donor_var_da, donor_p_rho_da, wrf_p_rho_da):
    wrf_hvi_da = xr.full_like(wrf_p_rho_da, fill_value=0).rename('diag')
    wrf_hvi_da.attrs.clear()

    for boundary_point_index, boundary_point in enumerate(wrf_p_rho_da.boundary_points):
        f = interpolate.interp1d(donor_p_rho_da.isel(boundary_points=boundary_point_index), donor_var_da.isel(boundary_points=boundary_point_index), kind='linear', bounds_error=False, fill_value=0)
        wrf_hvi_da.isel(boundary_points=boundary_point_index)[:] = f(wrf_p_rho_da.isel(boundary_points=boundary_point_index))

    return wrf_hvi_da


def hor_interpolate_3d_field_on_wrf_grid(donor_var_da, wrf_grid_ds):  #, wrf_ny, wrf_nx, wrf_lon, wrf_lat):
    field_hor_da = xr.DataArray(  # wrap into xarray structure
        data=np.zeros([donor_var_da.level.size, wrf_grid_ds.south_north.size, wrf_grid_ds.west_east.size]),
        dims=["level", "south_north", "west_east"],
        coords=dict(level=donor_var_da.level, ),
    )

    for level_index, level in enumerate(donor_var_da.level):
        f = interpolate.RectBivariateSpline(donor_var_da.lat, donor_var_da.lon, donor_var_da.isel(level=level_index))
        field_hor_da.isel(level=level_index).values[:] = f(wrf_grid_ds.XLAT, wrf_grid_ds.XLONG, grid=False).reshape(wrf_grid_ds.south_north.size, wrf_grid_ds.west_east.size)

    return field_hor_da


def ver_interpolate_3d_field_on_wrf_grid(donor_var_da, donor_p_rho_da, wrf_p_rho_da):
    wrf_hvi_da = xr.full_like(wrf_p_rho_da, fill_value=0).rename('diag')
    wrf_hvi_da.attrs.clear()

    for lat_index, lon_index in itertools.product(range(wrf_p_rho_da.south_north.size), range(wrf_p_rho_da.west_east.size)):  # This creates a single loop that visits every grid cell
        f = interpolate.interp1d(donor_p_rho_da.isel(south_north=lat_index, west_east=lon_index), donor_var_da.isel(south_north=lat_index, west_east=lon_index), kind='linear', bounds_error=False, fill_value=0)
        wrf_hvi_da.isel(south_north=lat_index, west_east=lon_index)[:] = f(wrf_p_rho_da.isel(south_north=lat_index, west_east=lon_index))

    return wrf_hvi_da


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


def update_boundaries(increment_on_bnd_da, wrf_met_ds, wrfbdy_f, name, index):
    '''
    wrf_met_ds is wrf_gid_ds, needs to provide grid information
    '''
    bdy_width = wrfbdy_f.dimensions['bdy_width'].size

    west, north, east, south = unravel_1d_boundary(increment_on_bnd_da, wrf_met_ds)
    west = west.transpose('bottom_top', 'south_north')
    east = east.transpose('bottom_top', 'south_north')
    south = south.transpose('bottom_top', 'west_east')
    north = north.transpose('bottom_top', 'west_east')

    west = west.expand_dims({'bdy_width': np.arange(bdy_width)})#, axis=1)  #
    east = east.expand_dims({'bdy_width': np.arange(bdy_width)})
    south = south.expand_dims({'bdy_width': np.arange(bdy_width)})
    north = north.expand_dims({'bdy_width': np.arange(bdy_width)})

    # print("\t\t\tUpdating BC for "+name)
    # make sure this is the dims order in da: bdy_width, bottom_top, south_north
    wrfbdy_f.variables[name + "_BXS"][index, :] = wrfbdy_f.variables[name + "_BXS"][index, :] + west
    wrfbdy_f.variables[name + "_BXE"][index, :] = wrfbdy_f.variables[name + "_BXE"][index, :] + east
    wrfbdy_f.variables[name + "_BYS"][index, :] = wrfbdy_f.variables[name + "_BYS"][index, :] + south
    wrfbdy_f.variables[name + "_BYE"][index, :] = wrfbdy_f.variables[name + "_BYE"][index, :] + north


def get_met_file_by_time(time):
    return "met_em.d01." + time + ".nc"


def is_child_domain_covered_by_parent_domain(parent_domain, child_domain):
    '''
    This check does not account for the distortions
    :param parent_domain:
    :param child_domain:
    :return:
    '''

    result = (min(child_domain.boundary_lons) > min(parent_domain.lon)) | (
            max(child_domain.boundary_lons) < max(parent_domain.lon)) | (
                     min(child_domain.boundary_lats) > min(parent_domain.lat)) | (
                     max(child_domain.boundary_lats) < max(parent_domain.lat))

    return result


def flatten_and_line_up_boundaries(da):
    '''
    Extract the four boundaries and align in 1D for interpolation
    '''
    north = da.isel(south_north=-1).rename({'west_east': 'boundary_points'})  # WRF_PRES[:, ny-1, :]
    south = da.isel(south_north=0).rename({'west_east': 'boundary_points'})  # WRF_PRES[:, 0, :]
    east = da.isel(west_east=-1).rename({'south_north': 'boundary_points'})  # WRF_PRES[:, :, nx-1]
    west = da.isel(west_east=0).rename({'south_north': 'boundary_points'})  # WRF_PRES[:, :, 0]
    diag_on_bnd_da = xr.concat([west, north, east, south], dim='boundary_points')
    return diag_on_bnd_da


def unravel_1d_boundary(diag_on_bnd_da, wrf_grid_ds):
    ny = wrf_grid_ds.sizes['south_north']
    nx = wrf_grid_ds.sizes['west_east']

    # 2. Slice the 1D array based on the original concat order:
    # Order was: [west (ny), north (nx), east (ny), south (nx)]

    west = diag_on_bnd_da.isel(boundary_points=slice(0, ny)).rename({'boundary_points': 'south_north'})
    north = diag_on_bnd_da.isel(boundary_points=slice(ny, ny + nx)).rename({'boundary_points': 'west_east'})
    east = diag_on_bnd_da.isel(boundary_points=slice(ny + nx, 2 * ny + nx)).rename({'boundary_points': 'south_north'})
    south = diag_on_bnd_da.isel(boundary_points=slice(2 * ny + nx, 2 * ny + 2 * nx)).rename({'boundary_points': 'west_east'})

    return west, north, east, south

