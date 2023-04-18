import xarray as xr
import re
import numpy as np
from scipy import interpolate

# p_top = 1  # Pa   (=0.01 hPa)
numbers = re.compile(r'(\d+)')


def hor_interpolate_3d_field_on_wrf_boubdary(emac_field_da, wrf_length, wrf_lon, wrf_lat):
    '''
    Horizontal interpolation of 3d Merra field on WRF boundary

    :param emac_field_da:
    :param wrf_length:
    :param wrf_lon:
    :param wrf_lat:
    :return:
    '''
    field_bnd_da = xr.DataArray(      # wrap into xarray structure
        data=np.zeros([emac_field_da.lev.size, wrf_length]),
        dims=["lev", "wrf_length"],
        coords=dict(
            lev=emac_field_da.lev,
        ),
    )

    for z_level in range(emac_field_da.lev.size):
        f = interpolate.RectBivariateSpline(np.flipud(emac_field_da.lat), emac_field_da.lon, np.flipud(emac_field_da[z_level, :, :]))  # flip to produce monotonic increasing lat
        field_bnd_da[z_level, :] = f(wrf_lat, wrf_lon, grid=False)

    return field_bnd_da


def ver_interpolate_3dfield_on_wrf_boubdary(MER_HOR_SPECIE_BND, MER_HOR_PRES_BND, WRF_PRES_BND, wrf_nz, wrf_length):
    '''
    Vertical interpolation of Merra boundary on WRF boundary

    :param MER_HOR_SPECIE_BND:
    :param MER_HOR_PRES_BND:
    :param WRF_PRES_BND:
    :param wrf_nz:
    :param wrf_length:
    :return:
    '''
    WRF_SPECIE_BND = np.zeros([wrf_nz, wrf_length])  # Required SPEC on WRF boundary

    for i in range(0, wrf_length):
        f = interpolate.interp1d(MER_HOR_PRES_BND[:, i], MER_HOR_SPECIE_BND[:, i], kind='linear', bounds_error=False,
                                 fill_value=0)
        WRF_SPECIE_BND[:, i] = f(WRF_PRES_BND[:, i])
    return WRF_SPECIE_BND


def hor_interpolate_3d_field_on_wrf_grid(emac_field_da, wrf_ny, wrf_nx, wrf_lon, wrf_lat):
    field_hor_da = xr.DataArray(  # wrap into xarray structure
        data=np.zeros([emac_field_da.lev.size, wrf_ny, wrf_nx]),
        dims=["lev", "wrf_ny", "wrf_nx"],
        coords=dict(
            lev=emac_field_da.lev,
        ),
    )

    for z_level in range(emac_field_da.lev.size):
        # interpolation requires monotonically increasing values, have to flip lat dimension
        f = interpolate.RectBivariateSpline(np.flipud(emac_field_da.lat), emac_field_da.lon, np.flipud(emac_field_da[z_level, :, :]))
        field_hor_da[z_level, :, :] = f(wrf_lat, wrf_lon, grid=False).reshape(wrf_ny, wrf_nx)

    return field_hor_da


def ver_interpolate_3d_field_on_wrf_grid(MER_HOR_SPECIE, MER_HOR_PRES, WRF_PRES, wrf_nz, wrf_ny, wrf_nx):
    WRF_SPECIE = np.zeros([wrf_nz, wrf_ny, wrf_nx])  # Required SPEC on WRF grid
    for x in range(0, wrf_nx, 1):
        for y in range(0, wrf_ny, 1):
            f = interpolate.interp1d(MER_HOR_PRES[:, y, x], MER_HOR_SPECIE[:, y, x], kind='linear', bounds_error=False, fill_value=0)
            WRF_SPECIE[:, y, x] = f(WRF_PRES[:, y, x])
    return WRF_SPECIE



