import numpy as np
import pandas as pd


### TOOLS FOR COMPUTING THE BLH ###

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

def magnus(t):
    """
    Magnus' formula for computing the saturated water vapour pressure.
    Input temperature 't' in degree Celsius.
    """
    return 6.112 * np.exp(17.62 * t / (243.14 + t))

def virtual_temp(T, p, q):
    """
    Computes the virtual temperature from T and p.
    """
    r_l = 0
    r_v = q # q ~ r
    Tv = T * (1 + 0.61 * r_v - r_l)
    return Tv

def sh_from_rh(rh, t, p):
    """
    Approximates specific humidity from relative humidity 'rh'.
    Takes temperature 't' in Kelvin and pressure 'p' in hPa.
    """
    e = rh * magnus(t - 273.15)
    return 0.622 * e / (p - 0.378 * e)

def bulk_richardson_number(df_synop, df_temp):
    """
    Computes the bulk Richardson number from four variables contained in the input data frames of the TEMP type
    and the corresponding SYNOP observations. The four variables are: Height 'Z' (gpm),
    temperature 'T' (K), specific humidity 'Q' (kg/kg), squared wind speed 'FF2' (m^2/s^2).
    The profiles are expected to be in descending order, i.e., the data frame index is monotonically decreasing
    pressure 'P' in Pa.

    The formula is taken from:
    'https://www.ecmwf.int/sites/default/files/elibrary/2017/17736-part-iv-physical-processes.pdf#section.3.10'
    """   
    g     = 9.81   #gravitational acceleration in m/s^2
    c_p   = 1005.  #speficif heat capacity in J/kgK
    R_dry = 287.1  #gas constant for dry air in J/kgK
    R_vap = 461.5  #gas constant for water vapor in J/kgK
    eps = R_vap / R_dry - 1.

    s_sfc = (c_p * df_synop["T"] * ( 1 + eps * df_synop["Q"]) + g * df_synop["Z"]).values
    s = c_p * df_temp["T"] * ( 1 + eps * df_temp["Q"] ) + g * df_temp["Z"]
    
    Ri_b = df_temp["Z"] * 2. * g * ( s - s_sfc )
    Ri_b /= ( (s + s_sfc - g * df_temp["Z"] - g * df_synop["Z"].values) * df_temp["FF2"])

    return Ri_b

def bri_alt(df_synop, df_temp):
    """
    Computes the bulk Richardson number from the commonly used formular in literature.
    For example: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JD018143
    """
    g     = 9.81   #gravitational acceleration in m/s^2
    
    theta_v_sfc = virtual_temp(pot_t(df_synop["T"], df_synop.index/100.), df_synop.index/100., df_synop["Q"]).item()
    theta_v = virtual_temp(pot_t(df_temp["T"], df_temp.index/100.), df_temp.index/100., df_temp["Q"])

    return g/theta_v_sfc * (theta_v - theta_v_sfc) * (df_temp["Z"] - df_synop["Z"].item()) / df_temp["FF2"]

def compute_blh(df_temp, Z_0):
    """
    Computes the boundary layer height by finding the level at which the bulk Richardson number 
    exceeds the critical value of 0.25. This level and the next lower one (or the ground level)
    are used to linearly interpolate to the height corresponding to exactly 0.25.
    Value is with respect to station height and in gpm.
    """
    # Scan for exceedance:
    for plvl in reversed(df_temp.index):
        if pd.notna(df_temp.loc[plvl, "Ri_b"]) and df_temp.loc[plvl, "Ri_b"] >= 0.25:
            Z_upper = df_temp.loc[plvl, "Z"]
            Rib_upper = df_temp.loc[plvl, "Ri_b"]
            break
    
    # Scan back down for lower bound:
    for plvl in df_temp.loc[plvl:].index:
        if pd.notna(df_temp.loc[plvl, "Ri_b"]) and df_temp.loc[plvl, "Ri_b"] < 0.25:
            Z_lower = df_temp.loc[plvl, "Z"]
            Rib_lower= df_temp.loc[plvl, "Ri_b"]
            break

    # If the last level:
    if plvl == df_temp.index[-1]:
        Z_lower = Z_0
        Rib_lower = 0.

    # Return linearly interpolated value w.r.t. station height
    return (0.25 - Rib_lower) / (Rib_upper - Rib_lower) * (Z_upper - Z_lower) + Z_lower - Z_0