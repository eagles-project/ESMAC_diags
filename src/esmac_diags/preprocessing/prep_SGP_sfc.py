"""
prepare long-term surface data for SGP site
options of output data into coarser resolution
"""

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d,\
                avg_time_1d, avg_time_2d, interp_time_1d, avg_height_2d
from esmac_diags.subroutines.read_surface import read_smpsb_pnnl,read_smps_bin
from esmac_diags.subroutines.quality_control import qc_remove_neg, qc_mask_qcflag, \
                qc_mask_qcflag_cpc, qc_correction_nanosmps
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_ARM

# acsmpath = '../../../raw_data/obs/SGP/sgpaosacsm/'
# ccnpath = '../../../raw_data/obs/SGP/sgpccn/'
# cpcpath = '../../../raw_data/obs/SGP/sgpcpc/'
# tdmapath = '../../../raw_data/obs/SGP/sgptdmasizeC1.b1/'
# uhsaspath = '../../../raw_data/obs/SGP/sgpaosuhsasE13.b1/'
# smpspath = '../../../raw_data/obs/SGP/sgpaossmpsE13.b1/'
# nanosmpspath = '../../../raw_data/obs/SGP/sgpaosnanosmpsE13.b1/'
# predatapath = '../../../prep_data/SGP/'
# year=2020
# dt=3600

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_ACSM(acsmpath, predatapath, year, dt=3600):
    """
    prepare acsm data

    Parameters
    ----------
    acsmpath : str
        input datapath
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                           
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(acsmpath, '*.b2.'+year+'*.nc'))
    lst.sort()
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    org = obsdata['total_organics'].load()
    qc_org = obsdata['qc_total_organics'].load()
    so4 = obsdata['sulfate'].load()
    qc_so4 = obsdata['qc_sulfate'].load()
    nh4 = obsdata['ammonium'].load()
    qc_nh4 = obsdata['qc_ammonium'].load()
    no3 = obsdata['nitrate'].load()
    qc_no3 = obsdata['qc_nitrate'].load()
    chl = obsdata['chloride'].load()
    qc_chl = obsdata['qc_chloride'].load()
    obsdata.close()
    
    # quality controls
    org = qc_mask_qcflag(org,qc_org)
    no3 = qc_mask_qcflag(no3,qc_no3)
    so4 = qc_mask_qcflag(so4,qc_so4)
    nh4 = qc_mask_qcflag(nh4,qc_nh4)
    chl = qc_mask_qcflag(chl,qc_chl)
    org = qc_remove_neg(org)
    no3 = qc_remove_neg(no3)
    so4 = qc_remove_neg(so4)
    nh4 = qc_remove_neg(nh4)
    chl = qc_remove_neg(chl)
    
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    org_new = median_time_1d(time, org, time_new)
    no3_new = median_time_1d(time, no3, time_new)
    so4_new = median_time_1d(time, so4, time_new)
    nh4_new = median_time_1d(time, nh4, time_new)
    chl_new = median_time_1d(time, chl, time_new)
    
    #%% output file
    outfile = predatapath + 'sfc_ACSM_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'org': (['time'], np.float32(org_new)),
                    'so4': (['time'], np.float32(so4_new)),
                    'nh4': (['time'], np.float32(nh4_new)),
                    'no3': (['time'], np.float32(no3_new)),
                    'chl': (['time'], np.float32(chl_new)),
                    },
                      coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['org'].attrs["long_name"] = org.long_name
    ds['org'].attrs["units"] = org.units
    ds['so4'].attrs["long_name"] = so4.long_name
    ds['so4'].attrs["units"] = so4.units
    ds['nh4'].attrs["long_name"] = nh4.long_name
    ds['nh4'].attrs["units"] = nh4.units
    ds['no3'].attrs["long_name"] = no3.long_name
    ds['no3'].attrs["units"] = no3.units
    ds['chl'].attrs["long_name"] = chl.long_name
    ds['chl'].attrs["units"] = chl.units
    
    ds.attrs["title"] = 'Aerosol composition from surface ACSM'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_ccn(ccnpath, predatapath, year, dt=3600):
    """
    prepare surface CCN data. 

    Parameters
    ----------
    ccnpath : char
        input datapath of CCN data
    predatapath : char
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
    
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)

    #%%
    lst = glob.glob(os.path.join(ccnpath, '*.b1.'+year+'*.nc'))
    lst.sort()
    # first data
    ccndata = xr.open_dataset(lst[0])
    ccntime = ccndata['time']
    coefs = ccndata['N_CCN_fit_coefs']
    ss_m = ccndata['supersaturation_calculated'].load().data
    ccn_m = ccndata['N_CCN'].load().data
    idx1 = np.nanargmin(np.abs(ss_m-0.1), axis=1)
    idx2 = np.nanargmin(np.abs(ss_m-0.2), axis=1)
    idx5 = np.nanargmin(np.abs(ss_m-0.5), axis=1)
    ss1 = np.array([ss_m[i,idx1[i]] for i in range(len(idx1))])
    ss2 = np.array([ss_m[i,idx2[i]] for i in range(len(idx2))])
    ss5 = np.array([ss_m[i,idx5[i]] for i in range(len(idx5))])
    ccn1 = np.array([ccn_m[i,idx1[i]] for i in range(len(idx1))])
    ccn2 = np.array([ccn_m[i,idx2[i]] for i in range(len(idx2))])
    ccn5 = np.array([ccn_m[i,idx5[i]] for i in range(len(idx5))])
    qc_ccn_tmp = ccndata['qc_N_CCN'].load().data
    qc_ccns = np.array([qc_ccn_tmp[i, [idx1[i], idx2[i], idx5[i]]] for i in range(len(idx5))])
    ccndata.close()
    for file in lst[1:]:
        ccndata = xr.open_dataset(file)
        ccntime = xr.concat([ccntime, ccndata['time']], dim="time")
        coefs = xr.concat([coefs, ccndata['N_CCN_fit_coefs']], dim="time")
        ss_m = ccndata['supersaturation_calculated'].load().data
        ccn_m = ccndata['N_CCN'].load().data
        idx1 = np.nanargmin(np.abs(ss_m-0.1), axis=1)
        idx2 = np.nanargmin(np.abs(ss_m-0.2), axis=1)
        idx5 = np.nanargmin(np.abs(ss_m-0.5), axis=1)
        ss1 = np.hstack((ss1, np.array([ss_m[i,idx1[i]] for i in range(len(idx1))])))
        ccn1 = np.hstack((ccn1, np.array([ccn_m[i,idx1[i]] for i in range(len(idx1))])))
        ss2 = np.hstack((ss2, np.array([ss_m[i,idx2[i]] for i in range(len(idx2))])))
        ccn2 = np.hstack((ccn2, np.array([ccn_m[i,idx2[i]] for i in range(len(idx2))])))
        ss5 = np.hstack((ss5, np.array([ss_m[i,idx5[i]] for i in range(len(idx5))])))
        ccn5 = np.hstack((ccn5, np.array([ccn_m[i,idx5[i]] for i in range(len(idx5))])))
        qc_ccn_tmp = ccndata['qc_N_CCN'].load().data
        qc_ccn_125 = np.array([qc_ccn_tmp[i, [idx1[i], idx2[i], idx5[i]]] for i in range(len(idx5))])
        qc_ccns = np.vstack((qc_ccns, qc_ccn_125))
        ccndata.close()
    
    #%% these are computed from CCN spectra polynomial fits
    #this accounts for fluctuations in supersaturation that are different than the target supersaturation
    #but the fits do not always work, so the the sample size is less than the measured CCN
    ccn1_fit = coefs[:,0] + coefs[:,1]*0.1 + coefs[:,2]*(0.1**2)
    ccn2_fit = coefs[:,0] + coefs[:,1]*0.2 + coefs[:,2]*(0.2**2)
    ccn5_fit = coefs[:,0] + coefs[:,1]*0.5 + coefs[:,2]*(0.5**2)
    
    #apply basic QC flags
    ccn1 = qc_mask_qcflag(ccn1, qc_ccns[:,0])
    ccn2 = qc_mask_qcflag(ccn2, qc_ccns[:,1])
    ccn5 = qc_mask_qcflag(ccn5, qc_ccns[:,2])
          
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    ccn1_fit_i = median_time_1d(ccntime, ccn1_fit, time_new)
    ccn2_fit_i = median_time_1d(ccntime, ccn2_fit, time_new)
    ccn5_fit_i = median_time_1d(ccntime, ccn5_fit, time_new)
    ccn1_measure = median_time_1d(ccntime, ccn1, time_new)
    ccn2_measure = median_time_1d(ccntime, ccn2, time_new)
    ccn5_measure = median_time_1d(ccntime, ccn5, time_new)
    ss1_i = median_time_1d(ccntime, ss1, time_new)
    ss2_i = median_time_1d(ccntime, ss2, time_new)
    ss5_i = median_time_1d(ccntime, ss5, time_new)
    
    #%% output file
    outfile = predatapath + 'sfc_CCN_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'ccn1_fit': ('time', np.float32(ccn1_fit_i)),
                     'ccn2_fit': ('time', np.float32(ccn2_fit_i)),
                     'ccn5_fit': ('time', np.float32(ccn5_fit_i)),
                     'ccn1_m': ('time', np.float32(ccn1_measure)),
                     'ccn2_m': ('time', np.float32(ccn2_measure)),
                     'ccn5_m': ('time', np.float32(ccn5_measure)),
                     'ss1': ('time', np.float32(ss1_i)),
                     'ss2': ('time', np.float32(ss2_i)),
                     'ss5': ('time', np.float32(ss5_i)),},
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['ccn1_fit'].attrs["long_name"] = "0.1% Cloud Condensation Nuclei"
    ds['ccn1_fit'].attrs["units"] = "cm-3"
    ds['ccn1_fit'].attrs["description"] = "Interpolated hourly values calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['ccn2_fit'].attrs["long_name"] = "0.2% Cloud Condensation Nuclei"
    ds['ccn2_fit'].attrs["units"] = "cm-3"
    ds['ccn2_fit'].attrs["description"] = "Interpolated hourly values calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['ccn5_fit'].attrs["long_name"] = "0.5% Cloud Condensation Nuclei"
    ds['ccn5_fit'].attrs["units"] = "cm-3"
    ds['ccn5_fit'].attrs["description"] = "Interpolated hourly values calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['ccn1_m'].attrs["long_name"] = "0.1% Cloud Condensation Nuclei - measured"
    ds['ccn1_m'].attrs["units"] = "cm-3"
    ds['ccn1_m'].attrs["description"] = "Interpolated hourly values ARM-measured CCN targetted to 0.1% SS. see SS1 for actual measured SS"
    ds['ss1'].attrs["long_name"] = "Actual Supersaturation targetted to 0.1%"
    ds['ss1'].attrs["units"] = "%"
    ds['ss1'].attrs["description"] = "measured SS that is closest to 0.1%. Interpolated into hourly. ccn1_m is measured at this SS"
    ds['ccn2_m'].attrs["long_name"] = "0.2% Cloud Condensation Nuclei"
    ds['ccn2_m'].attrs["units"] = "cm-3"
    ds['ccn2_m'].attrs["description"] = "Interpolated hourly values ARM-measured CCN targetted to 0.2% SS. see SS2 for actual measured SS"
    ds['ss2'].attrs["long_name"] = "Actual Supersaturation targetted to 0.2%"
    ds['ss2'].attrs["units"] = "%"
    ds['ss2'].attrs["description"] = "measured SS that is closest to 0.2%. Interpolated into hourly. ccn2_m is measured at this SS"
    ds['ccn5_m'].attrs["long_name"] = "0.5% Cloud Condensation Nuclei"
    ds['ccn5_m'].attrs["units"] = "cm-3"
    ds['ccn5_m'].attrs["description"] = "Interpolated hourly values ARM-measured CCN targetted to 0.5% SS. see SS5 for actual measured SS"
    ds['ss5'].attrs["long_name"] = "Actual Supersaturation targetted to 0.5%"
    ds['ss5'].attrs["units"] = "%"
    ds['ss5'].attrs["description"] = "measured SS that is closest to 0.5%. Interpolated into hourly. ccn5_m is measured at this SS"
    
    ds.attrs["title"] = 'Surface CCN number concentration'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_cloud_2d(armbepath, predatapath, height_out, year, dt=3600):
    """
    prepare cloud fraction data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    height_out : numpy array
        vertical dimension of output data. will average the original ARSCL resolution to it.
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                           
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbecldradC1.c1.'+year+'*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    height = obsdata['height'].load()
    cloud = obsdata['cld_frac'].load()
    qc_cloud = obsdata['qc_cld_frac'].load()
    obsdata.close()    
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cloud_i = np.full((len(time_new),len(height)), np.nan)
    for kk in range(len(height)):
        # quality controls. For ARMBE cloud fraction, remove data with <30% valid points within 1-hr window 
        cl = cloud[:,kk]
        cl[qc_cloud[:,kk]>=2] = np.nan
        # interpolate into standard time
        cloud_i[:,kk] = np.interp(time_new, time, cl)
        
    cloud_o = avg_time_2d(height,cloud_i.T,height_out).T
        
    
    #%% output file
    outfile = predatapath + 'cloud_2d_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cloud': (['time','height'], np.float32(cloud_o))
                    },
                     coords={'time': ('time', time_new), 'height':('height',np.float32(height_out))})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cloud'].attrs["long_name"] = "Cloud Fraction"
    ds['cloud'].attrs["units"] = "%"
    ds['cloud'].attrs["description"] = "Cloud Fraction based on radar and MPL"
    
    ds.attrs["title"] = 'ARMBE hourly 2-d cloud fraction data derived from ARSCL data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'interpolated into each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_cloudheight_ARSCL(arsclpath, predatapath, year, dt=3600):
    """
    prepare cloud base and top height data at ARM sites from ARSCL
    include multi-layer clouds
    
    Parameters
    ----------
    arsclpath : char
        input datapath.  
    predatapath : char
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string

    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
            
    #%% read in data
    lst1 = glob.glob(os.path.join(arsclpath, 'sgparscl*.'+year+'*.nc'))
    lst1.sort()
    arscldata = xr.open_mfdataset(lst1, combine='by_coords')
    arscltime = arscldata['time'].load()
    cbh = arscldata['cloud_base_best_estimate'].load()
    cbhs = arscldata['cloud_layer_base_height'].load()
    cths = arscldata['cloud_layer_top_height'].load()
    arscldata.close()
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # data treatments
    cbh = qc_remove_neg(cbh, remove_zero='True')
    for ll in range(cbhs.shape[1]):
        cbhs[:,ll] = qc_remove_neg(cbhs[:,ll], remove_zero='True')
        cths[:,ll] = qc_remove_neg(cths[:,ll], remove_zero='True')
    cth = np.nanmax(cths,axis=1)  # cloud top height for all clouds
        
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cbh_new = avg_time_1d(arscltime, cbh, time_new)
    cth_new = avg_time_1d(arscltime, cth, time_new)
    cbhs_new = avg_time_2d(arscltime, cbhs, time_new)
    cths_new = avg_time_2d(arscltime, cths, time_new)
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # output file
    outfile = predatapath + 'cloudheight_ARSCL_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cbh': ('time', np.float32(cbh_new)),
                     'cth': ('time', np.float32(cth_new)),
                     'cbhs': (['time', 'layer'], np.float32(cbhs_new)),
                     'cths': (['time', 'layer'], np.float32(cths_new)),
                    },
                     coords={'time': ('time', time_new), 'layer': ('layer', np.arange(cbhs.shape[1]))})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cbh'].attrs["long_name"] = 'cloud base height, best estimate from ceilometer and MPL'
    ds['cbh'].attrs["units"] = "m"
    ds['cth'].attrs["long_name"] = 'cloud top height, the highest level of cths'
    ds['cth'].attrs["units"] = "m"
    ds['cbhs'].attrs["long_name"] = 'Base height of hydrometeor layers for up to 10 layers'
    ds['cbhs'].attrs["units"] = "m"
    ds['cths'].attrs["long_name"] = 'Top height of hydrometeor layers for up to 10 layers'
    ds['cths'].attrs["units"] = "m"
    
    
    ds.attrs["title"] = "cloud heights for up to 10 cloud layers from ARSCL data"
    ds.attrs["description"] = 'average value of each time window'
    ds.attrs["inputfile_sample"] = lst1[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CPC(cpcpath, predatapath, year, dt=3600):
    """
    prepare CPC and CPCu data

    Parameters
    ----------
    cpcpath : str
        input datapath for CPC (10nm)
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                           
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst1 = glob.glob(os.path.join(cpcpath, '*.b1.'+year+'*.nc'))
    lst1.sort()
    obsdata = xr.open_mfdataset(lst1, combine='by_coords')
    time10 = obsdata['time']
    cpc10 = obsdata['concentration'].load()
    qc_cpc10 = obsdata['qc_concentration'].load()
    obsdata.close()
    
    # quality controls
    cpc10 = qc_mask_qcflag_cpc(cpc10, qc_cpc10.astype(int).data)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cpc10_new = median_time_1d(time10, cpc10, time_new)
    
    #%% output file
    outfile = predatapath + 'sfc_CPC_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cpc10': (['time'], np.float32(cpc10_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cpc10'].attrs["long_name"] = 'CPC aerosol number concentration (>10nm)'
    ds['cpc10'].attrs["units"] = '1/cm3'
    
    ds.attrs["title"] = 'Aerosol number concentration from surface CPC and CPCu'
    ds.attrs["inputfile_sample"] = lst1[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_UHSAS(uhsaspath, predatapath, year, dt=3600):
    """
    prepare UHSAS data

    Parameters
    ----------
    uhsaspath : str
        input datapath for UHSAS
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                           
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% read in data
    lst = glob.glob(os.path.join(uhsaspath, '*.b1.'+year+'*.nc'))
    if len(lst)==0:
        print('UHSAS is not available in year '+year)
        return
    lst.sort()
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time'].data
    uhsas = obsdata['dN_dlogDp'].load().data
    qc_uhsas = obsdata['qc_total_N_conc'].load()
    dmin = obsdata['diameter_optical_lower_bound'].load().data
    dmax = obsdata['diameter_optical_upper_bound'].load().data
    obsdata.close()
    
    size = (dmin[0,:]+dmax[0,:])/2
    dlogDp = np.log10(dmax[0,:]/dmin[0,:])
    uhsas = qc_mask_qcflag(uhsas,qc_uhsas)
    uhsas = qc_remove_neg(uhsas)
    uhsas = uhsas*np.tile(dlogDp,[len(time),1])
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    # startdate = np.datetime_as_string(np.datetime64(time[0]))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1]))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    uhsas_new = median_time_2d(time, uhsas, time_new)
        
    idx100 = dmin[0,:]>=100
    uhsas100_new = np.nansum(uhsas_new[:,idx100], 1)
    uhsas100_new[uhsas100_new==0] = np.nan
    
    #%% output file
    outfile = predatapath + 'sfc_UHSAS_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'uhsas_all': (['time', 'size'], uhsas_new),
                    'uhsas100': (['time'], uhsas100_new),      
                    'dlogDp': (['size'], np.float32(dlogDp)),   
                    },
                      coords={'time': ('time', time_new), 'size': ('size', size)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['uhsas_all'].attrs["long_name"] = 'aerosol number size distribution'
    ds['uhsas_all'].attrs["units"] = '1/cm3'
    ds['uhsas100'].attrs["long_name"] = 'aerosol number concentration for size >100nm'
    ds['uhsas100'].attrs["units"] = '1/cm3'
    ds['dlogDp'].attrs["units"] = "N/A"
    ds['dlogDp'].attrs["description"] = "to calculate dN/dlogDp"
    
    ds.attrs["title"] = 'Aerosol number concentration and size distribution from UHSAS'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_SMPS(smpspath, nanosmpspath, predatapath, year, dt=3600):
    """
    prepare SMPS data

    Parameters
    ----------
    smpspath : str
        input datapath for SMPS
    nanosmpspath : str
        input datapath for nanoSMPS
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                           
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% read in data
    lst = glob.glob(os.path.join(smpspath, '*.b1.'+year+'*.nc'))
    if len(lst)==0:
        print('SMPS is not available in year '+year)
        return
    lst.sort()
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time1 = obsdata['time'].data
    smps = obsdata['dN_dlogDp'].load()
    qc_smps = obsdata['qc_total_N_conc'].load()
    size = obsdata['diameter_mobility'].load()
    bounds = obsdata['diameter_mobility_bounds'].load().data
    obsdata.close()
    
    dmin = bounds[0,:,0]
    dmax = bounds[0,:,1]
    dlogDp = np.log10(dmax/dmin)
    smps = qc_mask_qcflag(smps,qc_smps)
    smps = smps*np.tile(dlogDp,[len(time1),1])
    
    
    lst2 = glob.glob(os.path.join(nanosmpspath, '*.b1.'+year+'*.nc'))
    lst2.sort()
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].data
    nanosmps = obsdata['dN_dlogDp'].load()
    qc_nanosmps = obsdata['qc_total_N_conc'].load()
    size = obsdata['diameter_mobility'].load()
    bounds = obsdata['diameter_mobility_bounds'].load().data
    obsdata.close()
    
    dmin2 = bounds[0,:,0]
    dmax2 = bounds[0,:,1]
    dlogDp2 = np.log10(dmax2/dmin2)
    if any(dlogDp!=dlogDp2):
        print (dmin, dmin2)
        raise ValueError(' ')
    nanosmps = qc_mask_qcflag(nanosmps,qc_nanosmps)
    nanosmps = nanosmps*np.tile(dlogDp,[len(time2),1])
    nanosmps = qc_correction_nanosmps(nanosmps)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    # startdate = np.datetime_as_string(np.datetime64(time1[0]))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time1[-1]))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    smps_new = median_time_2d(time1, smps, time_new)
    nanosmps_new = median_time_2d(time2, nanosmps, time_new)
        
    # combine SMPS and nanoSMPS
    smps_new[:,0:80]=nanosmps_new[:,0:80]
    
    idx100 = dmin[:]>=100
    smps100_new = np.nansum(smps_new[:,idx100], 1)
    smps100_new[smps100_new==0] = np.nan
    
    
    #%% output file
    outfile = predatapath + 'sfc_SMPS_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'smps_all': (['time', 'size'], smps_new),
                    'smps100': (['time'], smps100_new),      
                    'dlogDp': (['size'], np.float32(dlogDp)),   
                    },
                      coords={'time': ('time', time_new), 'size': ('size', size.data)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['smps_all'].attrs["long_name"] = 'aerosol number size distribution'
    ds['smps_all'].attrs["units"] = '1/cm3'
    ds['smps100'].attrs["long_name"] = 'aerosol number concentration for size >100nm'
    ds['smps100'].attrs["units"] = '1/cm3'
    ds['dlogDp'].attrs["long_name"] = "d log10 Dp"
    ds['dlogDp'].attrs["units"] = "N/A"
    ds['dlogDp'].attrs["description"] = "to calculate dN/dlogDp"
    
    ds.attrs["title"] = 'Aerosol number concentration and size distribution from SMPS'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_TDMA(tdmapath, predatapath, year, dt=3600):
    """
    prepare aerosol number size distribution and concentration for size>100nm from TDMA
    
    Parameters
    ----------
    tdmapath : char
        input datapath
    predatapath : char
        output datapath
    filename_out : char
        output filename.
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
        
    #%% settings
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% read in data
    lst = glob.glob(os.path.join(tdmapath, '*.b1.'+year+'*.cdf'))
    if len(lst)==0:
        print('TDMA is not available in year '+year)
        return
    lst.sort()
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time'].data
    tdma = obsdata['number_concentration'].load()
    qc_tdma = obsdata['qc_number_concentration'].load()
    size = obsdata['diameter'].load()
    obsdata.close()
    
    size = size.data*1000.   # um to nm
    # find dlogDp
    d = np.log10(size[1:]/size[0:-1])
    if np.std(d)>0.01:
        print(d)
        raise ValueError('dlogDp should be fairly consistent for TDMA. The standard deviation is too large. ' \
                         + 'check dlogDp for each size above:')
    dlogDp = np.ones_like(size) * np.mean(d)
    
    tdma = qc_mask_qcflag(tdma,qc_tdma)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    # startdate = np.datetime_as_string(np.datetime64(time[0]))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1]))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    tdma_new = median_time_2d(time, tdma, time_new)
    tdma_new = tdma_new * np.mean(d)
    
    idx100 = size>=100
    tdma100 = np.nansum(tdma_new[:,idx100], 1)
    tdma100[tdma100==0] = np.nan
    
    #%% output file
    outfile = predatapath + 'sfc_TDMA_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                'tdma_all': (['time', 'size'], np.float32(tdma_new)),
                'tdma100': (['time'], np.float32(tdma100)),        
                'dlogDp': (['size'], np.float32(dlogDp)),         
                },
                 coords={'time': ('time', time_new), 'size': ('size', np.float32(size))})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['size'].attrs["description"] = "Midpoint of geometric mean mobility diameter"
    ds['tdma_all'].attrs["long_name"] = 'aerosol number size distribution'
    ds['tdma_all'].attrs["units"] = '1/cm3'
    ds['tdma_all'].attrs["description"] = "from tdma, average into hourly"
    ds['tdma100'].attrs["long_name"] = 'aerosol number concentration for size >100nm'
    ds['tdma100'].attrs["units"] = '1/cm3'
    ds['tdma100'].attrs["description"] = "from tdma, average into hourly"
    ds['dlogDp'].attrs["long_name"] = "d log10 Dp"
    ds['dlogDp'].attrs["units"] = "N/A"
    ds['dlogDp'].attrs["description"] = "to calculate dN/dlogDp"
    
    ds.attrs["title"] = 'Aerosol number concentration and size distribution from TDMA'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mfrsr_cod(mfrsrpath,  predatapath, year, dt=3600):
    """
    prepare cloud optical depth data from MFRSR

    Parameters
    ----------
    mfrsrpath : str
        input datapath for MFRSR
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                   
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(mfrsrpath, '*.c1.'+year+'*.cdf'))
    lst.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst[0])
    time = mfrsrdata['time']
    cod = mfrsrdata['optical_depth_instantaneous']
    qc_cod = mfrsrdata['qc_optical_depth_instantaneous']
    cos_d = mfrsrdata['cosine_solar_zenith_angle']
    mfrsrdata.close()
    for file in lst[1:]:
        mfrsrdata = xr.open_dataset(file)
        time = xr.concat([time, mfrsrdata['time']], dim="time")
        cod = xr.concat([cod, mfrsrdata['optical_depth_instantaneous']], dim="time")
        qc_cod = xr.concat([qc_cod, mfrsrdata['qc_optical_depth_instantaneous']], dim="time")
        cos_d = xr.concat([cos_d, mfrsrdata['cosine_solar_zenith_angle']], dim="time")
        
        mfrsrdata.close()
    
    # quality controls
    cod.load()
    qc_cod.load()
    cod = qc_mask_qcflag(cod, qc_cod)
    
    # treat all missing COD with zero and keep nighttime (cosine_solar_zenith_angle <0.2) removed
    cod[np.isnan(cod)] = 0
    cod[cos_d<=0.2] = np.nan
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cod_new = avg_time_1d(time, cod, time_new)
    
    #%% output file
    outfile = predatapath + 'cod_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cod': (['time'], np.float32(cod_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cod'].attrs["long_name"] = 'cloud optical depth'
    ds['cod'].attrs["units"] = cod.units
    
    ds.attrs["title"] = 'cloud optical depth from surface MFRSR'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean value of each time window. treat missing as zero during daytime (cosine_zenith_angle>0.2)'
    ds.attrs['note'] = 'From MFRSRCLDDOD technical report (DOE/SC-ARM-TR-047): the results are only valid for '+\
            '“horizontally homogeneous” stratiform clouds with optical depths larger than approximately 7. '+\
            'The retrieval assumes a single cloud layer consisting solely of liquid water drops.'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mfrsr_Reff(mfrsrpath,  predatapath, year, dt=3600):
    """
    prepare cloud effective radius data from MFRSR

    Parameters
    ----------
    mfrsrpath : str
        input datapath for MFRSR
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(mfrsrpath, '*.c1.'+year+'*.cdf'))
    lst.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst[0])
    time = mfrsrdata['time']
    reff = mfrsrdata['effective_radius_instantaneous']
    qc_reff = mfrsrdata['qc_effective_radius_instantaneous']
    mfrsrdata.close()
    for file in lst[1:]:
        mfrsrdata = xr.open_dataset(file)
        time = xr.concat([time, mfrsrdata['time']], dim="time")
        reff = xr.concat([reff, mfrsrdata['effective_radius_instantaneous']], dim="time")
        qc_reff = xr.concat([qc_reff, mfrsrdata['qc_effective_radius_instantaneous']], dim="time")
        mfrsrdata.close()
    
    # quality controls
    reff.load()
    qc_reff.load()
    reff = qc_mask_qcflag(reff, qc_reff)
    
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    reff_new = median_time_1d(time, reff, time_new)
    
    #%% output file
    outfile = predatapath + 'reff_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'reff': (['time'], np.float32(reff_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['reff'].attrs["long_name"] = 'cloud effective radius'
    ds['reff'].attrs["units"] = reff.units
    
    ds.attrs["title"] = 'cloud effective radius from surface MFRSR'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_LWP(armbepath, mfrsrpath, predatapath, year, dt=3600):
    """
    prepare liquid water path
    Although LWP is measured by microwave radiometer (MWR), it is processed in 
    different ways and stored in different datastreams.
    Here we get LWP from two datasets: ARMBE and MFRSR

    Parameters
    ----------
    armbepath : str
        input datapath for ARMBE (MWR data)
    mfrsrpath : str
        input datapath for MFRSR
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string

    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% read in data
    lst1 = glob.glob(os.path.join(armbepath, '*armbecldradC1.c1.'+year+'*.nc'))
    obsdata = xr.open_mfdataset(lst1, combine='by_coords')
    
    time1 = obsdata['time']
    lwp = obsdata['lwp']
    qc_lwp = obsdata['qc_lwp']
    obsdata.close()
    
    lwp.load()
    qc_lwp.load()
    
    # quality controls. For ARMBE lwp, remove data with <30% valid points within 1-hr window 
    lwp[qc_lwp>=2] = np.nan
    
    #%% read in MFRSR LWP for comparison
    lst2 = glob.glob(os.path.join(mfrsrpath, '*.c1.'+year+'*.cdf'))
    lst2.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst2[0])
    time2 = mfrsrdata['time']
    lwp2 = mfrsrdata['lwp']
    qc_lwp2 = mfrsrdata['qc_lwp']
    mfrsrdata.close()
    for file in lst2[1:]:
        mfrsrdata = xr.open_dataset(file)
        time2 = xr.concat([time2, mfrsrdata['time']], dim="time")
        lwp2 = xr.concat([lwp2, mfrsrdata['lwp']], dim="time")
        qc_lwp2 = xr.concat([qc_lwp2, mfrsrdata['qc_lwp']], dim="time")
        mfrsrdata.close()
        
    lwp2.load()
    qc_lwp2.load()
    # lwp2 = qc_mask_qcflag(lwp2, qc_lwp2)
    lwp2 = lwp2*1000. # change unit from kg/m2 (mm) to g/m2
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    lwp_new = avg_time_1d(time1, lwp, time_new)
    lwp2_new = avg_time_1d(time2, lwp2, time_new)
    
    #%% output file
    outfile = predatapath + 'LWP_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwp_armbe': ('time', np.float32(lwp_new)),
                    'lwp_mfrsr': ('time', np.float32(lwp2_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwp_armbe'].attrs["long_name"] = "liquid water path"
    ds['lwp_armbe'].attrs["units"] = "g/m2"
    ds['lwp_armbe'].attrs["description"] = "liquid water path from ARMBE data based on MWR measurements"
    ds['lwp_mfrsr'].attrs["long_name"] = "liquid water path"
    ds['lwp_mfrsr'].attrs["units"] = "g/m2"
    ds['lwp_mfrsr'].attrs["description"] = "liquid water path from MFRSR data based on MWR measurements"
    
    ds.attrs["title"] = 'surface-retrieved cloud liquid water path'
    ds.attrs["inputfile_sample"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["description"] = 'mean value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_LTS(armbepath, arsclpath, predatapath, year, dt=3600):
    """
    prepare lower tropospheric stability (potential temperature difference between 700hPa and surface) from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    arsclpath : str
        input datapath for cloud base information
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
               
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbeatmC1.c1.'+year+'*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    pres = obsdata['pressure'].load()
    ht = obsdata['height'].load()
    T = obsdata['temperature_p'].load()
    Tz = obsdata['temperature_h'].load()
    Ts = obsdata['temperature_sfc'].load()
    ps = obsdata['pressure_sfc'].load()
    obsdata.close()    
    
    # in ARMBEATM, pres[12] is 700hPa. check if not
    p700_idx = 12
    if np.abs(pres[p700_idx]-700)>0.1:
        raise ValueError('this level is not 700hPa, check: '+str(pres[p700_idx]))
    p850_idx = 6
    if np.abs(pres[p850_idx]-850)>0.1:
        raise ValueError('this level is not 850hPa, check: '+str(pres[p850_idx]))
    
    # remove one timestep that temperature is incorrect
    idx = T[:,30]>300
    T[idx.data,:]=np.nan
    Tz[idx.data,:]=np.nan
    
    # interpolate to fill NaNs 
    Tz_interp = Tz.interpolate_na(dim='time')
    
    #%% read in data
    lst1 = glob.glob(os.path.join(arsclpath, 'sgparscl*.'+year+'*.nc'))
    lst1.sort()
    arscldata = xr.open_mfdataset(lst1, combine='by_coords')
    arscltime = arscldata['time'].load()
    cbh = arscldata['cloud_base_best_estimate'].load()
    arscldata.close()
    
    # data treatments
    cbh = qc_remove_neg(cbh, remove_zero='True')
    cbh_armbe = avg_time_1d(arscltime, cbh, time)
    
    # use dry static energy as an approxy to calcalate potential temperature (theta = T + gz/Cp) 
    thetadiff_cb = np.empty((0))
    for tt in range(len(time)):
        if np.isnan(cbh_armbe[tt]):
            thetadiff_cb = np.append(thetadiff_cb, np.nan)
        else:
            T_cb = np.interp(cbh_armbe[tt], ht, Tz_interp[tt,:])
            thetadiff = T_cb - Ts[tt] + 9.8/1005.7*cbh_armbe[tt]
            thetadiff_cb = np.append(thetadiff_cb, thetadiff)
            
    # only extract valid data when sounding is available
    LTS700_valid  = np.empty((0))
    LTS850_valid  = np.empty((0))
    time700_valid = np.empty(0,dtype='datetime64')
    time850_valid = np.empty(0,dtype='datetime64')
    
    for tt in range(len(time)):
        if ~np.isnan(T[tt,p700_idx]):
            time700_valid = np.hstack((time700_valid,time[tt]))
            theta_s = Ts[tt].data * (1000/ps[tt].data)**0.286
            theta_700 = T[tt,p700_idx].data * (1000/700)**0.286       
            LTS700_valid = np.append(LTS700_valid, theta_700-theta_s )
        if ~np.isnan(T[tt,p850_idx]):
            time850_valid = np.hstack((time850_valid,time[tt]))
            theta_s = Ts[tt].data * (1000/ps[tt].data)**0.286
            theta_850 = T[tt,p850_idx].data * (1000/850)**0.286       
            LTS850_valid = np.append(LTS850_valid, theta_850-theta_s )
        
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    LTS700_new = avg_time_1d(time700_valid, LTS700_valid, time_new)
    LTS850_new = avg_time_1d(time850_valid, LTS850_valid, time_new)
    thetadiff_new = avg_time_1d(time, thetadiff_cb, time_new)
        
    #%% output file
    outfile = predatapath + 'LTS_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'LTS700': (['time'], np.float32(LTS700_new)),
                    'LTS850': (['time'], np.float32(LTS850_new)),
                    'thetadiff_cb': (['time'], np.float32(thetadiff_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['LTS700'].attrs["long_name"] = "Lower Tropospheric Stability"
    ds['LTS700'].attrs["units"] = "K"
    ds['LTS700'].attrs["description"] = "potential temperature difference between 700hPa and surface, " + \
        "calculated from radiosonde  and surface data"
    ds['LTS850'].attrs["long_name"] = "Lower Tropospheric Stability"
    ds['LTS850'].attrs["units"] = "K"
    ds['LTS850'].attrs["description"] = "potential temperature difference between 850hPa and surface, " + \
        "calculated from radiosonde  and surface data"
    ds['thetadiff_cb'].attrs["long_name"] = "potential temperature difference between cloud base and surface"
    ds['thetadiff_cb'].attrs["units"] = "K"
    ds['thetadiff_cb'].attrs["description"] = "potential temperature differencebetween cloud base and surface"
    
    ds.attrs["title"] = 'Lower Tropospheric Stability from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_precip(armbepath, predatapath, year, dt=3600):
    """
    prepare surface precipitation data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
               
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbeatmC1.c1.'+year+'*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    precip = obsdata['precip_rate_sfc'].load()
    obsdata.close()    
                
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    precip_new = avg_time_1d(time, precip, time_new)
        
    
    #%% output file
    outfile = predatapath + 'precip_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'precip': (['time'], np.float32(precip_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['precip'].attrs["long_name"] = "Surface Precipitation Rate"
    ds['precip'].attrs["units"] = "mm/hr"
    
    ds.attrs["title"] = 'surface precipitation data from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'precipitation at SGP is measured by Tipping Bucket Rain Gauge (TBRG). mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_radiation(armbepath, predatapath, year, dt=3600):
    """
    prepare surface radiation data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                   
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbecldradC1.c1.'+year+'*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    lwdn = obsdata['lwdn'].load()
    lwup = obsdata['lwup'].load()
    swdn = obsdata['swdn'].load()
    swup = obsdata['swup'].load()
    obsdata.close()    
                
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    lwdn_new = avg_time_1d(time, lwdn, time_new)
    swdn_new = avg_time_1d(time, swdn, time_new)
    lwup_new = avg_time_1d(time, lwup, time_new)
    swup_new = avg_time_1d(time, swup, time_new)
        
    #%% output file
    outfile = predatapath + 'sfc_radiation_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwdn': (['time'], np.float32(lwdn_new)),
                    'swdn': (['time'], np.float32(swdn_new)),
                    'lwup': (['time'], np.float32(lwup_new)),
                    'swup': (['time'], np.float32(swup_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwdn'].attrs["long_name"] = "Surface downward longwave flux"
    ds['lwdn'].attrs["units"] = "W/m2"
    ds['swdn'].attrs["long_name"] = "Surface downward shortwave flux"
    ds['swdn'].attrs["units"] = "W/m2"
    ds['lwup'].attrs["long_name"] = "Surface upward longwave flux"
    ds['lwup'].attrs["units"] = "W/m2"
    ds['swup'].attrs["long_name"] = "Surface upward shortwave flux"
    ds['swup'].attrs["units"] = "W/m2"
    
    ds.attrs["title"] = 'surface radiative flux data from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_totcld(armbepath, predatapath, year, dt=3600):
    """
    prepare total cloud fraction data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbecldradC1.c1.'+year+'*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    cf_arscl = obsdata['tot_cld']
    qc_cf_arscl = obsdata['qc_tot_cld']
    cf_tsi = obsdata['tot_cld_tsi']
    qc_cf_tsi = obsdata['qc_tot_cld_tsi']
    cf_visst = obsdata['cld_tot']
    obsdata.close()
    
    cf_arscl.load()
    cf_tsi.load()
    cf_visst.load()
    qc_cf_arscl.load()
    qc_cf_tsi.load()
    
    # quality controls. For ARMBE cloud fraction, <30% valid points within 1-hr window are flagged and removed
    cf_arscl[qc_cf_arscl>=2] = np.nan
    cf_tsi[qc_cf_tsi>=2] = np.nan
                
    # change unit from 1 to %
    cf_arscl = cf_arscl*100
    cf_tsi = cf_tsi*100
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cf_arscl_new = avg_time_1d(time, cf_arscl, time_new)
    cf_tsi_new = avg_time_1d(time, cf_tsi, time_new)
    cf_visst_new = avg_time_1d(time, cf_visst, time_new)
    
    #%% output file
    outfile = predatapath + 'totcld_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'tot_cld_arscl': ('time', np.float32(cf_arscl_new)),
                    'tot_cld_tsi': ('time', np.float32(cf_tsi_new)),
                    'tot_cld_visst': ('time', np.float32(cf_visst_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['tot_cld_arscl'].attrs["long_name"] = "cloud_area_fraction"
    ds['tot_cld_arscl'].attrs["units"] = "%"
    ds['tot_cld_arscl'].attrs["description"] = "Total cloud fraction based on ARSCL radar/MPL measurements"
    ds['tot_cld_tsi'].attrs["long_name"] = "cloud_area_fraction"
    ds['tot_cld_tsi'].attrs["units"] = "%"
    ds['tot_cld_tsi'].attrs["description"] = "Total cloud fraction based on total sky imager, 100 degree FOV"
    ds['tot_cld_visst'].attrs["long_name"] = "cloud_area_fraction"
    ds['tot_cld_visst'].attrs["units"] = "%"
    ds['tot_cld_visst'].attrs["description"] = "Total cloud fraction based on VISST satellite product"
    
    
    ds.attrs["title"] = 'total cloud fraction from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Ndrop(ndroppath, predatapath, year, dt=3600):
    """
    prepare cloud deoplet number concentration from Ndrop data
    
    Parameters
    ----------
    ndroppath : char
        input datapath.  
    predatapath : char
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
        
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(ndroppath, '*.c1.'+year+'*.nc'))
    lst.sort()
    
    obsdata = xr.open_dataset(lst[0])
    time = obsdata['time']
    nd = obsdata['drop_number_conc']
    qc_nd = obsdata['qc_drop_number_conc']
    ctype = obsdata['cloud_base_type']
    obsdata.close()
    for file in lst[1:]:
        obsdata = xr.open_dataset(file)
        time = xr.concat([time, obsdata['time']], dim="time")
        nd = xr.concat([nd, obsdata['drop_number_conc']], dim="time")
        qc_nd = xr.concat([qc_nd, obsdata['qc_drop_number_conc']], dim="time")
        ctype = xr.concat([ctype, obsdata['cloud_base_type']], dim="time")
        obsdata.close()
        
    # quality control
    nd = qc_mask_qcflag(nd,qc_nd)
    nd = nd*1e-6   # m-3 to cm-3
    
    # # exclude ice clouds or multi-layer clouds
    # nd[ctype!=1] = np.nan
    
    # exclude small values
    nd[nd<10] = np.nan
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, nd, time_new)
    # nd_new = avg_time_1d(time, nd, time_new)
    
    #%% output file
    outfile = predatapath + 'Ndrop_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cdnc': ('time', np.float32(nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cdnc'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['cdnc'].attrs["units"] = "cm-3"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from ARM retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd_ARMretrieval(mfrsrpath, arsclpath, predatapath, year, dt=3600):
    """
    prepare cloud deoplet number concentration (Nd) data at ARM sites
    input data is cloud optical depth from MFRSR, LWP from MWR (in the MFRSR data), 
    and cloud top/base height from ARSCL
    using  ARM surface-based retrievals algorithm
    https://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-140.pdf 
    assuming cloud is homogeneous between cloud top and cloud base in the column
    
    Parameters
    ----------
    mfrsrpath : char
        input datapath. 
    arsclpath : char
        input datapath.  
    predatapath : char
        output datapath
    year : int
        specify the year of data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    year = str(year)   # change to string
    
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst1 = glob.glob(os.path.join(arsclpath, 'sgparscl*.'+year+'*.nc'))
    lst1.sort()
    arscldata = xr.open_mfdataset(lst1, combine='by_coords')
    arscltime = arscldata['time']
    cbh = arscldata['cloud_base_best_estimate']
    cbhs = arscldata['cloud_layer_base_height']
    cths = arscldata['cloud_layer_top_height']
    arscldata.close()
    
    lst2 = glob.glob(os.path.join(mfrsrpath, '*.c1.'+year+'*.cdf'))
    lst2.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst2[0])
    mfrsrtime = mfrsrdata['time']
    lwp = mfrsrdata['lwp']
    qc_lwp = mfrsrdata['qc_lwp']
    cod = mfrsrdata['optical_depth_instantaneous']
    qc_cod = mfrsrdata['qc_optical_depth_instantaneous']
    mfrsrdata.close()
    for file in lst2[1:]:
        mfrsrdata = xr.open_dataset(file)
        mfrsrtime = xr.concat([mfrsrtime, mfrsrdata['time']], dim="time")
        lwp = xr.concat([lwp, mfrsrdata['lwp']], dim="time")
        qc_lwp = xr.concat([qc_lwp, mfrsrdata['qc_lwp']], dim="time")
        cod = xr.concat([cod, mfrsrdata['optical_depth_instantaneous']], dim="time")
        qc_cod = xr.concat([qc_cod, mfrsrdata['qc_optical_depth_instantaneous']], dim="time")
        mfrsrdata.close()
    
    lwp.load()
    qc_lwp.load()
    cod.load()
    qc_cod.load()
    cbh.load()
    cbhs.load()
    cths.load()
    
    # lwp = qc_mask_qcflag(lwp,qc_lwp)  # do not mask LWP since clearsky (LWP=0) is flagged
    cod = qc_mask_qcflag(cod,qc_cod)
    cbh = qc_remove_neg(cbh, remove_zero='True')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # calculate CDNC
    # calculate cloud depth for the lowest cloud layer in ARSCL
    cth = cths[:,0]
    H = cths[:,0] - cbh
    H.load()
    H = qc_remove_neg(H, remove_zero='true')
    H[cbhs[:,1] > 0] = np.nan
    
    # filter some data samples
    # lwp = qc_remove_neg(lwp, remove_zero='true')
    # cod = qc_remove_neg(cod, remove_zero='true')
    # cod[cod>500] = np.nan
    # lwp[lwp>1] = np.nan
    # lwp[lwp<0.01] = np.nan
    H[cth>5000.] = np.nan   # remove deep clouds with cloud top >5km
    lwp[np.logical_or(lwp<0.02, lwp>0.3)] = np.nan
    cod[np.logical_or(cod<2, cod>60)] = np.nan
    
    # calculate CDNC first then average into 1hr
    time = mfrsrtime.data
    H_tmp = np.interp(np.int64(time), np.int64(arscltime), H)
    nd = calc_cdnc_ARM(lwp, cod, H_tmp)
    
    # exclude small values
    nd[nd<10] = np.nan
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, nd, time_new)
        
    #%% output file
    outfile = predatapath + 'Nd_ARMretrieval_SGP_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cdnc': ('time', np.float32(nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cdnc'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['cdnc'].attrs["units"] = "cm-3"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from ARM retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["reference"] = 'https://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-140.pdf '
    ds.attrs["inputfile_sample"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
