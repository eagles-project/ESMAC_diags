# all quality controls for data in ESMAC_Diags
# examples of quality controls:
# 1. qc flags in ARM data
# 2. minimum and maximum cutoffs
# 3. systematic data corrections
# 4. additional data masking

import numpy as np


#%% mask data with cloud flag, typically remove cloud_flag=1 
def qc_mask_cloudflag(data,cflag):
    if len(data.shape)==1:
        data[cflag==1]=np.nan
    elif len(data.shape)==2:
        data[cflag==1,:]=np.nan
    else: 
        error
    return(data)

#%% mask data with qc_flag, typically remove all qc_flag != 0 
def qc_mask_qcflag(data,qc):
    if len(data.shape)==1:
        data[qc!=0]=np.nan
    elif len(data.shape)==2:
        data[qc!=0,:]=np.nan
    else: 
        error
    return(data)

#%% mask data with qc_flag, with some exceptions
# maximum value flag for CPC is too low (8000 cm-3) that some reasonable large values at SGP are flagged 
# retain these data by keeping qc flag >20
def qc_mask_qcflag_cpc(data,qc):
    data[np.logical_and(qc>0, qc<=20)]=np.nan
    return(data)

#%% mask aircraft takeoff/landing time
# mask time is set as 30min to exclude possible land contamination for ocean measurements
def qc_mask_takeoff_landing(time,data):
    # time should be in unit of seconds
    idx=np.logical_or(time<(time[0]+1800), time>(time[-1]-1800))
    if len(data.shape)==1:
        data[idx]=np.nan
    elif len(data.shape)==2:
        data[:,idx]=np.nan
    else:
        error
    return(data)

#%% remove negative values
def qc_remove_neg(data):
    data[data<0]=np.nan
    return(data)

#%% set a max value for CCN measurement
# different maximum threshold for different supersaturations
# SS can be a fixed value of changable for different ccn values
def qc_ccn_max(ccn,SS):
    ccn[np.logical_and(SS<0.2,ccn>2000)]=np.nan
    ccn[np.logical_and(SS<0.6,ccn>4000)]=np.nan
    ccn[ccn>5000]=np.nan
    return(ccn)
    
#%% set a max value for CN measurements
# different threshold for different minimal cufoff size
def qc_cn_max(data,size):
    if size==3:
        data[data>5e4]=np.nan
    elif size==10:
        data[data>2.5e4]=np.nan
    elif size==100:
        data[data>0.5e4]=np.nan
    return(data)

#%% set a minimum threshold for ARM G1 flight CPC measurements
# for HISCALE and ACEENA only
def qc_cpc_air(cpc3, cpc10):
    cpc3[cpc3<20]=np.nan
    cpc10[cpc10<10]=np.nan
    return(cpc3, cpc10)

#%% correction on nanoSMPS
# nanoSMPS at SGP is known to overcount particle number. make a correction to ensure smooth transition with SMPS
def qc_correction_nanosmps(data):
    # the fraction of 3.8 is chosen so that nanoSMPS and SMPS is smoothly transit at ~18 nm
    data=data/3.8
    return(data)

#%% set max value for ACSM measured total organic matters
def qc_acsm_org_max(data):
    data[data>10]=np.nan # max value set as 10 ug/m3
    return(data)

#%% set a maximum value for NCAR research flight UHSAS measurements
def qc_uhsas_RF_NCAR(data):
    data[data>500]=np.nan
    return(data)