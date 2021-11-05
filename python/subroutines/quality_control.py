"""
# all quality controls for data in ESMAC_Diags
# examples of quality controls:
# 1. qc flags in ARM data
# 2. minimum and maximum cutoffs
# 3. systematic data corrections
# 4. additional data masking
"""
import numpy as np

#%%
def qc_mask_cloudflag(data, cflag):
    """
    mask data with cloud flag, typically remove cloud_flag = 1 
    """
    if len(data.shape) == 1:
        data[cflag == 1] = np.nan
    elif len(data.shape) == 2:
        data[cflag == 1, :] = np.nan
    else: 
        raise ValueError("Error: input data should be 1-d or 2-d")
    return(data)

#%%
def qc_mask_qcflag(data, qc):
    """
    mask data with qc_flag, typically remove all qc_flag  !=  0 
    """
    if len(data.shape) == 1:
        data[qc != 0] = np.nan
    elif len(data.shape) == 2:
        data[qc != 0, :] = np.nan
    else: 
        raise ValueError("Error: input data should be 1-d or 2-d")
    return(data)

#%%
def qc_mask_qcflag_cpc(data, qc):
    """
    mask data with qc_flag, with some exceptions
    maximum value flag for CPC is too low (8000 cm-3) that some reasonable large values at SGP are flagged 
    retain these data by keeping qc flag >20
    """
    data[np.logical_and(qc>0, qc <= 20)] = np.nan
    return(data)

#%%
def qc_mask_takeoff_landing(time, data):
    """
    mask aircraft takeoff/landing time
    mask time is set as 30min to exclude possible land contamination for ocean measurements

    Parameters
    ----------
    time : measurement time in seconds
    data : input data
    time is only 1-dimensional while data can be 1- or 2-dimensional

    Returns
    -------
    masked data

    """
    # time should be in unit of seconds
    idx = np.logical_or(time<(time[0]+1800), time>(time[-1]-1800))
    if len(data.shape) == 1:
        data[idx] = np.nan
    elif len(data.shape) == 2:
        data[:, idx] = np.nan
    else:
        raise ValueError("Error: input data should be 1-d or 2-d")
    return(data)

#%%
def qc_remove_neg(data):
    """
    remove negative values
    """
    data[data<0] = np.nan
    return(data)

#%%
def qc_ccn_max(ccn, SS):
    """
    set a max value for CCN measurement
    different maximum threshold for different supersaturations
    SS can be a fixed value of changable for different ccn values
    
    Parameters
    ----------
    ccn : ccn number concentration
    SS : supersaturation
    ccn and SS should be the same size

    Returns
    -------
    ccn with quality controls

    """
    ccn[np.logical_and(SS<0.2, ccn>2000)] = np.nan
    ccn[np.logical_and(SS<0.6, ccn>4000)] = np.nan
    ccn[ccn>5000] = np.nan
    return(ccn)
    
#%%
def qc_cn_max(data, size):
    """
    set a max value for CN measurements
    different threshold for different minimal cufoff size
    
    Parameters
    ----------
    data : input aerosol number concentration
    size : minimum detected aerosol size (nm) in CPC

    Returns
    -------
    data with value greater than maximum removed

    """
    if size == 3:
        data[data>5e4] = np.nan
    elif size == 10:
        data[data>2.5e4] = np.nan
    elif size == 100:
        data[data>0.5e4] = np.nan
    return(data)

#%%
def qc_cpc_air(cpc3, cpc10):
    """
    set a minimum threshold for ARM G1 flight CPC measurements
    """
    cpc3[cpc3<20] = np.nan
    cpc10[cpc10<10] = np.nan
    return(cpc3, cpc10)

#%%
def qc_correction_nanosmps(data):
    """
    nanoSMPS at SGP is known to overcount particle number. make a correction to ensure smooth transition with SMPS
    the fraction of 3.8 is chosen so that nanoSMPS and SMPS is smoothly transit at ~18 nm
    """
    data = data/3.8
    return(data)

#%%
def qc_acsm_org_max(data):
    """
    set max value for ACSM measured total organic matters
    """
    data[data>10] = np.nan # max value set as 10 ug/m3
    return(data)

#%%
def qc_uhsas_RF_NCAR(data):
    """
    set a maximum value for NCAR research flight UHSAS measurements
    """
    
    data[data>500] = np.nan
    return(data)