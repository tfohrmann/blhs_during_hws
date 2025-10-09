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

def pot_t(T, p):
    '''Computes the potential temperature given a vertical temperature 
    profile "T" (K) and a pressure profile "p" (hPa).'''
    
    #Constants
    R_L = 287.
    c_p = 1005. #J/kgK
    p_0 = 1000. #hPa
    
    #Compute the potential temperature
    theta = T * (p_0/p)**(R_L/c_p)
    
    return(theta)

def virtual_temp(T, p, q):
    """
    Computes the virtual temperature from T and p.
    """
    r_l = 0
    #e = magnus(T - 273.15)
    #r_v = 0.622 * e / (p - e)
    # these old lines would mean saturation all the time... but q ~ r, so:
    r_v = q
    Tv = T * (1 + 0.61 * r_v - r_l)
    return Tv

def bri_alt(h, t, p, v, q, theta_v_sfc):
    """
    BRI ALT: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JD018143
    """
    g     = 9.81   #gravitational acceleration in m/s^2
    
    theta_v = virtual_temp(pot_t(t, p), p, q)

    Ri_b = g / theta_v_sfc * (theta_v - theta_v_sfc) * h / v

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

theta_v_sfc = xr.open_dataset(f"theta_v2m.{year}{month:02}.nc")["theta_v_2m"]

if step > 1:
    ds_q = ds_q.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))
    ds_t = ds_t.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))
    ds_p = ds_p.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))
    ds_v = ds_v.isel(rlat=slice(0, None, step), rlon=slice(0, None, step))


Ri_b_alt = bri_alt(heights, ds_t, ds_p, ds_v, ds_q, theta_v_sfc).compute()
condition = Ri_b_alt.sel({"lev": slice(None, None, -1)}) >= 0.25
first_index_from_back = condition.argmax(dim="lev")
upper = (Ri_b_alt["lev"].size - 1 - first_index_from_back)
lower = xr.where(upper < 39, upper + 1, 39)
blh_bri_alt = (heights.isel(lev=upper) - heights.isel(lev=lower)) / (Ri_b_alt.isel(lev=upper) - Ri_b_alt.isel(lev=lower)) * (0.25 - Ri_b_alt.isel(lev=upper)) + heights.isel(lev=upper)
# Correct instances where NaNs get produced at upper=lower
blh_bri_alt = xr.where(upper == 39, (heights.isel(lev=upper) - 0) / (Ri_b_alt.isel(lev=upper) - 0) * (0.25 - Ri_b_alt.isel(lev=upper)) + heights.isel(lev=upper), blh_bri_alt)

###### Save
data_vars = {"blh": (["time", "rlat", "rlon"], blh_bri_alt.transpose("time", "rlat", "rlon").data, 
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
