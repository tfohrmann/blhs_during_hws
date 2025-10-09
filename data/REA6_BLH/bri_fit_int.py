#!/usr/bin/env python
# coding: utf-8

import numpy as np
from pandas import date_range
import xarray as xr
import dask.array as da

from os import listdir, path
import sys

step = 1

args = sys.argv
year = int(args[1])
month = int(args[2])
day = int(args[3])

path_rea = f"./{year}/"


# Load static information:
rea_static = xr.open_dataset("../rea_static_rot.nc")
if step > 1:
    rea_static = rea_static.isel(rlat=slice(0, -1, step), rlon=slice(0, None, step))
rlats = rea_static.rlat.values
rlons = rea_static.rlon.values

sea_mask = (rea_static.LSM.values >= 0.5)

heights = (rea_static.HHL - rea_static.h0) #height levels above ground


# Functions:
def get_fnames(directory, var, dts):
    """
    This function matches all required datetimes 'dts' to matching files in the path "directory"
    for the variable "var".
    """
    fnames = []
    all_files = listdir(directory)

    for y, m, d, h in zip(dts.year, dts.month, dts.day, dts.hour):
        match = [item for item in all_files if item.startswith(f"{var}.3D.{y}{m:02}{d:02}{h:02}")]
        if len(match) != 1:
            print(f"Warning: No file or multiple files found for variable {var} @ {y}.{m}.{d} {h}h")
            continue
        fnames.append(directory + match[0])
    
    return(fnames)


def bulk_richardson_number(h, t, q, v, varname_lvl=None):
    """
    Computes the bulk Richardson number from the four input vertical profiles: Height 'h' (m),
    temperature 't' (K), specific humidity 'q' (no unit), wind speed 'v' (m/s).
    The profiles are expected to be in descending order, i.e., height decreases with model level,
    i.e., h[n_lvl] is the lowest height. If the variables are of type numpy.ndarray, the height
    dimension needs to be indexed by the last axis. If the variables are stored in xarray.DataArrays,
    'varname_lvl' refers to the name of the height coordinate and the level with the heighest index
    needs to be closest to the surface.

    The formula is taken from:
    'https://www.ecmwf.int/sites/default/files/elibrary/2017/17736-part-iv-physical-processes.pdf#section.3.10'

    :param h: numpy.ndarray or xarray.DataArray
    :param t: numpy.ndarray or xarray.DataArray
    :param q: numpy.ndarray or xarray.DataArray
    :param v: numpy.ndarray or xarray.DataArray
    :return Ri_b: numpy.ndarray or xarray.DataArray
    """
    assert (type(t) == type(q)) and (type(t) == type(v))
    assert type(t) == xr.core.dataarray.DataArray
    
    g     = 9.81   #gravitational acceleration in m/s^2
    c_p   = 1005.  #speficif heat capacity in J/kgK
    R_dry = 287.1  #gas constant for dry air in J/kgK
    R_vap = 461.5  #gas constant for water vapor in J/kgK
    p_0   = 1000.  #reference surface pressure in hPa
    
    eps = R_vap / R_dry - 1.
    s = c_p * t*(1 + eps*q) + g*h

    Ri_b = h * 2. * g * (s - s.sel({varname_lvl: s[varname_lvl].max()}))
    Ri_b /= ( (s + s.sel({varname_lvl: s[varname_lvl].max()}) - g*h - g*h[-1]) * v)

    return Ri_b


# Execution:
dts = date_range(start=np.datetime64(f"{year}-{month:02}-{day:02}T01"), periods=24, freq="1h")

fnames_qv = get_fnames("./QV/", "QV", dts)
fnames_t = get_fnames("./T/", "T", dts)
fnames_pp = get_fnames("./PP/", "PP", dts)
fnames_U = get_fnames("./U/", "U", dts)
fnames_V = get_fnames("./V/", "V", dts)


ds_q = xr.open_mfdataset(fnames_qv)["var51"]
ds_q["rlon"], ds_q["rlat"] = rea_static["rlon"], rea_static["rlat"]

ds_t = xr.open_mfdataset(fnames_t)["var11"]
ds_t["rlon"], ds_t["rlat"] = rea_static["rlon"], rea_static["rlat"]

ds_p = rea_static["p0"].values + xr.open_mfdataset(fnames_pp)["var139"]/100.
ds_p["rlon"], ds_p["rlat"] = rea_static["rlon"], rea_static["rlat"]

ds_U = xr.open_mfdataset(fnames_U)["var33"]
ds_U["rlon"], ds_U["rlat"] = rea_static["rlon"], rea_static["rlat"]

ds_V = xr.open_mfdataset(fnames_V)["var34"]
ds_V["rlon"], ds_V["rlat"] = rea_static["rlon"], rea_static["rlat"]

ds_v = ds_U**2 + ds_V**2


if step > 1:
    ds_q = ds_q.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))
    ds_t = ds_t.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))
    ds_p = ds_p.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))
    ds_v = ds_v.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))


Ri_b = bulk_richardson_number(heights, ds_t, ds_q, ds_v, varname_lvl="lev")
upper = (Ri_b >= 0.25).argmin(dim="lev").compute()
lower = upper - 1
blh_bri = (heights.isel(lev=upper) - heights.isel(lev=lower)) / (Ri_b.isel(lev=upper) - Ri_b.isel(lev=lower)) * (0.25 - Ri_b.isel(lev=upper)) + heights.isel(lev=upper)



###### Save
data_vars = {"blh": (["time", "rlat", "rlon"], blh_bri.transpose("time", "rlat", "rlon").data, 
                     {"units": "m", "long_name": "planetary boundary layer height above ground"})}
coords = {"time": (["time"], dts),
          "rlat": (["rlat"], rlats),
          "rlon": (["rlon"], rlons)}
attrs = {"creation_date": str(np.datetime64("now"))[:10], 
         "author": "Till Fohrmann", 
         "email": "tfohrmann@uni-bonn.de"}

blhs = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)

if path.exists("./blhs.nc"):
    blhs_prev = xr.open_dataset("./blhs.nc")
    (xr.concat([blhs_prev, blhs], dim="time")).to_netcdf("./blhs_merged.nc", "w")
else:
    blhs.to_netcdf("./blhs.nc", "w")
