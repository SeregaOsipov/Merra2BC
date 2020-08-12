# -*- coding: utf-8 -*-
import config
import time
from netCDF4 import Dataset

start_time = time.time()

zero = 1e-16
fields_to_zero = ['o3', 'co', 'so2', 'sulf']

# ---------------------------------------
# INITIAL CONDITIONS
nc_fp = config.wrf_dir + "/" + config.wrf_input_file
print("SETTING TO ZERO INITIAL CONDITIONS in {}".format(nc_fp))
wrfinput = Dataset(nc_fp, 'r+')
for field in fields_to_zero:
    print("Setting to zero IC for ", field)
    wrfinput.variables[field][:] = zero
wrfinput.close()

# BOUNDARY CONDITIONS
nc_fp = config.wrf_dir + "/" + config.wrf_bdy_file
print("\n\nSETTING TO ZERO BOUNDARY CONDITIONS AND TENDENCIES in {}".format(nc_fp))
wrfbddy = Dataset(nc_fp, 'r+')
for field in fields_to_zero:
    print("Setting to zero BC for ", field)
    wrfbddy.variables[field + "_BXS"][:] = zero
    wrfbddy.variables[field + "_BXE"][:] = zero
    wrfbddy.variables[field + "_BYS"][:] = zero
    wrfbddy.variables[field + "_BYE"][:] = zero

    print("Setting to zero Tendency BC for ", field)
    wrfbddy.variables[field + "_BTXS"][:] = zero
    wrfbddy.variables[field + "_BTXE"][:] = zero
    wrfbddy.variables[field + "_BTYS"][:] = zero
    wrfbddy.variables[field + "_BTYE"][:] = zero
wrfbddy.close()
print("--- %s seconds ---" % (time.time() - start_time))
