"""
prepare long-term surface data for ENA site
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
                avg_time_1d, avg_time_2d
from esmac_diags.subroutines.read_ARMdata import read_uhsas
from esmac_diags.subroutines.quality_control import qc_remove_neg, qc_mask_qcflag, \
                qc_mask_qcflag_cpc
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_ARM


# acsmpath = '../../../data/ENA/enaaosacsm/'
# armbepath = '../../../data/ENA/enaarmbe/'
# arsclpath = '../../../data/ENA/enaarscl/'
# arsclbndpath = '../../../data/ENA/enaarsclkazrbnl1kolliasC1.c0/'
# ccnpath = '../../../data/ENA/enaaosccn1colspectraC1.b1/'
# cpcpath = '../../../data/ENA/enaaoscpcf/'
# uhsaspath = '../../../data/ENA/enaaosuhsasC1.b1/'
# mfrsrpath = '../../../data/ENA/enamfrsrcldod1minC1.c1/'
# metpath = '../../../data/ENA/enametC1.b1/'
# parspath =  '../../../data/ENA/enapars/'
# radfluxpath =  '../../../data/ENA/enaradflux/'
# Nd_WUpath = '../../../data/ENA/Wu_etal_retrieval/'
# ndroppath = '../../../data/ENA/enandrop/'
# # predatapath = 'C:/Users/tang357/Downloads/ENA/'
# dt=3600
# year='2017'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_ACSM(acsmpath, predatapath, year, dt=300):
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

    # data resolution is 30-min, so interpolate for finer resolution; a mean could also be used for coarser resolution is warranted
    if dt >= 1800:
        org_new = median_time_1d(time, org, time_new)
        no3_new = median_time_1d(time, no3, time_new)
        so4_new = median_time_1d(time, so4, time_new)
        nh4_new = median_time_1d(time, nh4, time_new)
        chl_new = median_time_1d(time, chl, time_new)
    if dt < 1800:
        org_new = interp_time_1d(time, org, time_new)
        no3_new = interp_time_1d(time, no3, time_new)
        so4_new = interp_time_1d(time, so4, time_new)
        nh4_new = interp_time_1d(time, nh4, time_new)
        chl_new = interp_time_1d(time, chl, time_new)
    
    #%% output file
    outfile = predatapath + 'sfc_ACSM_ENA_'+year+'.nc'
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
    if dt >= 1800:
        ds.attrs["description"] = 'median value of each time window'
    if dt < 1800:
        ds.attrs["description"] = 'interpolated value from 30-min resolution data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_ccn(ccnpath, predatapath, year, dt=300):
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
    
    lst = glob.glob(os.path.join(ccnpath, '*.b1.'+year+'*.nc'))
    ccndata = xr.open_mfdataset(lst, combine='by_coords')
    ccntime = ccndata['time'].load()
    # base_time = ccndata['base_time']
    coefs = ccndata['N_CCN_fit_coefs'].load()
    ss_m = ccndata['supersaturation_calculated'].load().data
    idx1 = np.nanargmin(np.abs(ss_m-0.1), axis=1)
    idx2 = np.nanargmin(np.abs(ss_m-0.2), axis=1)
    idx5 = np.nanargmin(np.abs(ss_m-0.5), axis=1)
  
    ccn_m = ccndata['N_CCN'].load().data
    ss1 = np.array([ss_m[i,idx1[i]] for i in range(len(idx1))])
    ss2 = np.array([ss_m[i,idx2[i]] for i in range(len(idx2))])
    ss5 = np.array([ss_m[i,idx5[i]] for i in range(len(idx5))])
    ccn1 = np.array([ccn_m[i,idx1[i]] for i in range(len(idx1))])
    ccn2 = np.array([ccn_m[i,idx2[i]] for i in range(len(idx2))])
    ccn5 = np.array([ccn_m[i,idx5[i]] for i in range(len(idx5))])
    qc_ccns = ccndata['qc_N_CCN'].load()
    qc_ccns = np.array([qc_ccns[i, [idx1[i], idx2[i], idx5[i]]] for i in range(len(idx5))])
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

    #apply to ccn fits
    ccn1_fit = qc_mask_qcflag(ccn1_fit, qc_ccns[:,0])
    ccn2_fit = qc_mask_qcflag(ccn2_fit, qc_ccns[:,1])
    ccn5_fit = qc_mask_qcflag(ccn5_fit, qc_ccns[:,2])

    #remove known bad data ranges
    if year == '2020':
      ccn5[:] = np.nan
      ccn5_fit[:] = np.nan
    if year == '2019':
      ind1 = np.where(ccntime.dt.month <= 3)
      ind2 = np.where(np.logical_and(ccntime.dt.month <= 4, ccntime.dt.day <= 2))
      ind3 = np.where(np.logical_and(ccntime.dt.month >= 9, ccntime.dt.day >= 22))
      ind4 = np.where(ccntime.dt.month >= 10)
      ccn5[ind1] = np.nan
      ccn5[ind2] = np.nan
      ccn5[ind3] = np.nan
      ccn5[ind4] = np.nan
      ccn5_fit[ind1] = np.nan
      ccn5_fit[ind2] = np.nan
      ccn5_fit[ind3] = np.nan
      ccn5_fit[ind4] = np.nan
    if year == '2018':
      ind1 = np.where(np.logical_and(ccntime.dt.month >= 12, ccntime.dt.day >= 8))
      ccn5[ind1] = np.nan
      ccn5_fit[ind1] = np.nan
          
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    #below differs from ACE-ENA code that uses median_time_1d and interp_time_1d functions
    ccn1_fit_i = np.interp(np.int64(time_new), np.int64(ccntime), ccn1_fit, left=np.nan, right=np.nan)
    ccn2_fit_i = np.interp(np.int64(time_new), np.int64(ccntime), ccn2_fit, left=np.nan, right=np.nan)
    ccn5_fit_i = np.interp(np.int64(time_new), np.int64(ccntime), ccn5_fit, left=np.nan, right=np.nan)
    ccn1_measure = np.interp(np.int64(time_new), np.int64(ccntime), ccn1, left=np.nan, right=np.nan)
    ccn2_measure = np.interp(np.int64(time_new), np.int64(ccntime), ccn2, left=np.nan, right=np.nan)
    ccn5_measure = np.interp(np.int64(time_new), np.int64(ccntime), ccn5, left=np.nan, right=np.nan)
    ss1_i = np.interp(np.int64(time_new), np.int64(ccntime), ss1, left=np.nan, right=np.nan)
    ss2_i = np.interp(np.int64(time_new), np.int64(ccntime), ss2, left=np.nan, right=np.nan)
    ss5_i = np.interp(np.int64(time_new), np.int64(ccntime), ss5, left=np.nan, right=np.nan)
    
    
    #%% output file
    outfile = predatapath + 'sfc_CCN_ENA_'+year+'.nc'
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
def prep_cloud_2d(armbepath, predatapath, height_out, year, dt=300):
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

    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    #%% read in data
    if dt >= 3600:
        lst = glob.glob(os.path.join(armbepath, '*armbecldradC1.c1.'+year+'*.nc'))
        obsdata = xr.open_mfdataset(lst, combine='by_coords')
        time = obsdata['time']
        height = obsdata['height'].load()
        cloud = obsdata['cld_frac'].load()
        qc_cloud = obsdata['qc_cld_frac'].load()
        obsdata.close()    
        
        cloud_i = np.full((len(time_new),len(height)), np.nan)
        for kk in range(len(height)):
            # quality controls. For ARMBE cloud fraction, remove data with <30% valid points within 1-hr window 
            cl = cloud[:,kk]
            cl[qc_cloud[:,kk]>=2] = np.nan
            # interpolate into standard time
            cloud_i[:,kk] = np.interp(time_new, time, cl)
            
        cloud_o = avg_time_2d(height,cloud_i.T,height_out).T

    if dt < 3600:
        lst = glob.glob(os.path.join(arsclpath, 'enaarsclkazr1kolliasC1.c0*.nc'))
        lst.sort()
        arscldata = xr.open_mfdataset(lst, combine='by_coords')
        time = arscldata['time'].load()
        height = arscldata['height'].load()
        tmpcloud_flag = arscldata['cloud_source_flag'].load() #0=missing; 1=clear, 2+=cloud
        arscldata.close()

        cloud_flag = xr.where(tmpcloud_flag >= 2, 1, 0) #sets cloud to 1 and no cloud to 0
        cloud_flag = xr.where(tmpcloud_flag == 0, np.nan, cloud_flag) #sets missing to NaN
        dt_new = time_new[1]-time_new[0]
        #%% count the number of cloudy points at each height in the new time interval and divide by all points at that height to get cloud fraction
        #%% do a half time interval offset so that the time arrays don't shift
        cloud_i = cloud_flag.resample(time = dt_new, offset = dt_new/2).sum()/cloud_flag.resample(time = dt_new, offset = dt_new/2).count()
        cloud_i['time'] = cloud_i['time'] + dt_new/2
        
        cloud_o = avg_height_2d(height,cloud_i,height_out)
    
    #%% output file
    outfile = predatapath + 'cloud_2d_ENA_'+year+'.nc'
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
    
    if dt >= 3600:
        ds.attrs["title"] = 'ARMBE hourly 2-d cloud fraction data derived from ARSCL data'
        ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
        ds.attrs["description"] = 'interpolated into each time window'
    if dt < 3600:
        ds.attrs["title"] = 'ARSCL 2-d cloud fraction data'
        ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
        ds.attrs["description"] = 'accumulated into each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_cloudheight_ARSCL(arsclpath, predatapath, year, dt=300):
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
    lst1 = glob.glob(os.path.join(arsclpath, 'enaarsclkazrbnd1kolliasC1.c0.'+year+'*.nc'))
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
    outfile = predatapath + 'cloudheight_ARSCL_ENA_'+year+'.nc'
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
def prep_CPC(cpcpath,  predatapath, year, dt=300):
    """
    prepare CPC data

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
    outfile = predatapath + 'sfc_CPC_ENA_'+year+'.nc'
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
    
    ds.attrs["title"] = 'Aerosol number concentration from surface CPC'
    ds.attrs["inputfile_sample"] = lst1[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CPC_withENAmask(aerosolmaskpath, predatapath, year, dt=300):
    """
    prepare aerosol number concentration data (>10nm) from CPC
    apply aerosol mask data specifically for ENA site
    aerosol mask data provided by Francesca Gallo and Allison Aiken 
    reference: Gallo et al., 2020 ACP. https://doi.org/10.5194/acp-20-7553-2020

    Parameters
    ----------
    aerosolmaskpath : char
        input datapath of aerosol number concentration with mask
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
        
    #%% settings
    
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data with mask
    # data from Francesca Gallo <effegal@gmail.com>. 
    # reference: https://acp.copernicus.org/articles/20/7553/2020/
    lst = glob.glob(os.path.join(aerosolmaskpath, 'CPC_ENA_AM_'+year+'*.txt'))
    lst.sort()
    dateall = []
    timeall = []
    cpcvalue = []
    maskflag = []
    for ll in range(len(lst)):
        filename=lst[ll]
        f = open(filename, 'r')
        for i in range(13):
            h = f.readline()
        h = f.readline()
        varname = h.strip()
        h = f.readline()
        varunit = h.strip()
        print(filename)
        print(varname)
        print(varunit)
        f.readline()
        for line in f:
            line = line.strip()  # remove \n
            columns = line.split()
            dateall.append(columns[0])
            timeall.append(columns[1])
            if columns[2][0] == 'n':
                cpcvalue.append(-9999.0)
            else:
                cpcvalue.append(float(columns[2]))
            maskflag.append(int(columns[3]))
        f.close()
    
    cpcvalue = np.array(cpcvalue)
    maskflag = np.array(maskflag,dtype='float')
    cpcvalue[cpcvalue<0] = np.nan
    maskflag[maskflag<-999] = np.nan
    
    mask_starttime = dateall[0] + ' ' + timeall[0]
    mask_endtime = dateall[-1] + ' ' + timeall[-1]
    mask_time = pd.date_range(start=mask_starttime,end=mask_endtime,freq="min")    
    
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cpcvalue_masked = np.array(cpcvalue)
    cpcvalue_masked[maskflag!=0] = np.nan
    maskcount = np.zeros(len(cpcvalue))
    maskcount[maskflag==0] = 1.
    
    cpc_1hr_nomask = median_time_1d(mask_time, cpcvalue, time_new)
    cpc_1hr_masked = median_time_1d(mask_time, cpcvalue_masked, time_new)
    valid_fraction = avg_time_1d(mask_time, maskcount, time_new)
    
    #%% output file
    outfile = predatapath + 'sfc_CPC_ENA_withmask_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cpc_origin': (['time'], np.float32(cpc_1hr_nomask)), 
                    'cpc_masked': (['time'], np.float32(cpc_1hr_masked)), 
                    'f_valid': (['time'], np.float32(valid_fraction)),              
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['f_valid'].attrs["long_name"] = 'fraction of valid data'
    ds['f_valid'].attrs["units"] = '1'
    ds['f_valid'].attrs["description"] = "number of valid data points divided by total points within the hour"
    ds['cpc_origin'].attrs["long_name"] = 'CPC aerosol number concentration (>10nm)'
    ds['cpc_origin'].attrs["units"] = '1/cm3'
    ds['cpc_origin'].attrs["description"] = "average into hourly without applying mask"
    ds['cpc_masked'].attrs["long_name"] = 'CPC aerosol number concentration (>10nm)'
    ds['cpc_masked'].attrs["units"] = '1/cm3'
    ds['cpc_masked'].attrs["description"] = "average into hourly, applying 1-min aerosol mask"
    
    ds.attrs["title"] = 'Aerosol number concentration from surface CPC with contamination mask from Gallo et al., 2020'
    ds.attrs["reference"] = 'https://doi.org/10.5194/acp-20-7553-2020'
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_UHSAS(uhsaspath, predatapath, year, dt=300):
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
    outfile = predatapath + 'sfc_UHSAS_ENA_'+year+'.nc'
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
def prep_LWP(armbepath, mfrsrpath, predatapath, year, dt=300):
    """
    prepare liquid water path
    Although LWP is measured by microwave radiometer (MWR), it is processed in 
    different ways and stored in different datastreams.
    Here we get LWP from two datasets: ARMBE and MFRSR

    Parameters
    ----------
    armbepath : str
        input datapath for ARMBE (MWR data)
    mwrpath : str
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
    
    # #%% read in MFRSR LWP for comparison
    # lst2 = glob.glob(os.path.join(mfrsrpath, '*.c1.'+year+'*.cdf'))
    # lst2.sort()
    # # first data
    # mfrsrdata = xr.open_dataset(lst2[0])
    # time2 = mfrsrdata['time']
    # lwp2 = mfrsrdata['lwp']
    # qc_lwp2 = mfrsrdata['qc_lwp']
    # mfrsrdata.close()
    # for file in lst2[1:]:
    #     mfrsrdata = xr.open_dataset(file)
    #     time2 = xr.concat([time2, mfrsrdata['time']], dim="time")
    #     lwp2 = xr.concat([lwp2, mfrsrdata['lwp']], dim="time")
    #     qc_lwp2 = xr.concat([qc_lwp2, mfrsrdata['qc_lwp']], dim="time")
    #     mfrsrdata.close()
        
    # lwp2.load()
    # qc_lwp2.load()
    # lwp2 = qc_mask_qcflag(lwp2, qc_lwp2)
    # lwp2 = lwp2*1000. # change unit from kg/m2 (mm) to g/m2

    #%% read in MWR LWP
    lst2 = glob.glob(os.path.join(mwrpath, '*.nc'))
    mwrdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = mwrdata['time']
    lwp2 = mwrdata['phys_lwp'] #units are g/m2
    qc_lwp2 = mwrdata['qc_phys_lwp']
    mwrdata.close()

    mwr_lwp.load()
    mwr_lwp[qc_mwr_lwp > 0] = np.nan
  
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    lwp_new = avg_time_1d(time1, lwp, time_new)
    lwp2_new = avg_time_1d(time2, lwp2, time_new)
    
    #%% output file
    outfile = predatapath + 'LWP_ENA_'+year+'.nc'
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
    ds['lwp_armbe'].attrs["description"] = "liquid water path from hourly ARMBE data based on MWR measurements"
    ds['lwp_mwr'].attrs["long_name"] = "liquid water path"
    ds['lwp_mwr'].attrs["units"] = "g/m2"
    ds['lwp_mwr'].attrs["description"] = "liquid water path from MWR retrievals"
    
    ds.attrs["title"] = 'surface-retrieved cloud liquid water path'
    ds.attrs["inputfile_sample"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["description"] = 'mean value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_LTS(armbepath, arsclpath, predatapath, year, dt=300):
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
    
    #%% read in cloud base data
    lst1 = glob.glob(os.path.join(arsclpath, 'enaarsclkazrbnd1kolliasC1.c0.'+year+'*.nc'))
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

    if dt >= 3600:
        LTS700_new = avg_time_1d(time700_valid, LTS700_valid, time_new)
        LTS850_new = avg_time_1d(time850_valid, LTS850_valid, time_new)
        thetadiff_new = avg_time_1d(time, thetadiff_cb, time_new)
    if dt < 3600:
        LTS700_new = interp_time_1d(time700_valid, LTS700_valid, time_new)
        LTS850_new = interp_time_1d(time850_valid, LTS850_valid, time_new)
        thetadiff_new = interp_time_1d(time, thetadiff_cb, time_new)
        
    #%% output file
    outfile = predatapath + 'LTS_ENA_'+year+'.nc'
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
    ds.attrs["description"] = 'mean of each time window (hourly or longer time); interpolated from hourly (shorter than hourly time)'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mfrsr_cod(mfrsrpath,  predatapath, year, dt=300):
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
    mfrsrdata.close()
    for file in lst[1:]:
        mfrsrdata = xr.open_dataset(file)
        time = xr.concat([time, mfrsrdata['time']], dim="time")
        cod = xr.concat([cod, mfrsrdata['optical_depth_instantaneous']], dim="time")
        qc_cod = xr.concat([qc_cod, mfrsrdata['qc_optical_depth_instantaneous']], dim="time")
        mfrsrdata.close()
    
    # quality controls
    cod.load()
    qc_cod.load()
    cod = qc_mask_qcflag(cod, qc_cod)
    
    
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    cod_new = avg_time_1d(time, cod, time_new)
    
    #%% output file
    outfile = predatapath + 'cod_ENA_'+year+'.nc'
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
    ds.attrs["description"] = 'mean value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mfrsr_Reff(mfrsrpath,  predatapath, year, dt=300):
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
    lwp_source = mfrsrdata['lwp_aource']
    mfrsrdata.close()
    for file in lst[1:]:
        mfrsrdata = xr.open_dataset(file)
        time = xr.concat([time, mfrsrdata['time']], dim="time")
        reff = xr.concat([reff, mfrsrdata['effective_radius_instantaneous']], dim="time")
        qc_reff = xr.concat([qc_reff, mfrsrdata['qc_effective_radius_instantaneous']], dim="time")
        lwp_source = xr.concat([lwp_soource, mfrsrdata['lwp_source']], dim="time")
        mfrsrdata.close()
    
    # quality controls
    reff.load()
    qc_reff.load()
    lwp_source.load()
    reff = qc_mask_qcflag(reff, qc_reff)
    # remove effective radii values not derived from LWP retrieved from MWR
    reff[lwp_source == 2] = np.nan
    reff[lwp_source < 0] = np.nan
    
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    reff_new = median_time_1d(time, reff, time_new)
    
    #%% output file
    outfile = predatapath + 'reff_ENA_'+year+'.nc'
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
def prep_precip(armbepath, metpath, parspath, predatapath, year, dt=300):
    """
    prepare surface precipitation data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    metpath : str
        input datapath for ORG data
    parspath : str
        input datapath for disdrometer data
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
        
    # #%% read in data (old way used ARMBE hourly rainfalls from MET, potentially from the PWD
    # lst = glob.glob(os.path.join(armbepath, '*armbeatmC1.c1.'+year+'*.nc'))
    # obsdata = xr.open_mfdataset(lst, combine='by_coords')
    # time = obsdata['time']
    # precip = obsdata['precip_rate_sfc'].load()
    # obsdata.close()    

    #%% read in data (new way uses the optical rain gauge (ORG) and Parsivel disdrometer rates that are better for light
    lst = glob.glob(os.path.join(metpath, '*.cdf'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    precip_org = obsdata['org_precip_rate_mean'].load()
    qc_precip_org = obsdata['qc_org_precip_rate_mean'].load()
    precip_pwd = obsdata['pwd_precip_rate_mean_1min'].load()
    qc_precip_pwd = obsdata['qc_pwd_precip_rate_mean_1min'].load()
    obsdata.close()

    precip_org = qc_mask_qcflag(precip_org, qc_precip_org)
    precip_pwd = qc_mask_qcflag(precip_pwd, qc_precip_pwd)

    lst = glob.glob(os.path.join(parspath, '*.nc'))
    parsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    precip_pars = obsdata['rain_rate'].load()
    parsdata.close()
  
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    precip_org_new = avg_time_1d(time, precip_org, time_new)
    precip_pwd_new = avg_time_1d(time, precip_pwd, time_new)
    precip_pars_new = avg_time_1d(time, precip_pars, time_new)
        
    #%% output file
    outfile = predatapath + 'precip_ENA_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'precip_org': (['time'], np.float32(precip_org_new)),
                    'precip_pwd': (['time'], np.float32(precip_pwd_new)),
                    'precip_pars': (['time'], np.float32(precip_pars_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['precip_org'].attrs["long_name"] = "ORG Surface Precipitation Rate"
    ds['precip_org'].attrs["units"] = "mm/hr"
    ds['precip_pwd'].attrs["long_name"] = "PWD Surface Precipitation Rate"
    ds['precip_pwd'].attrs["units"] = "mm/hr"
    ds['precip_pars'].attrs["long_name"] = "Parsivel Surface Precipitation Rate"
    ds['precip_pars'].attrs["units"] = "mm/hr"

    ds.attrs["title"] = 'surface precipitation data from MET data and Parsivel disdrometer'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'precipitation at ENA is measured by ARM optical rain gauge (ORG), present weather detector (PWD), and Parsivel2 disdrometer. mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_radiation(armbepath, radfluxpath, predatapath, year, dt=300):
    """
    prepare surface radiation data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    radfluxpath : str
        input datapath.
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
    if dt >= 3600:
        lst = glob.glob(os.path.join(armbepath, '*armbecldradC1*.nc'))
        obsdata = xr.open_mfdataset(lst, combine='by_coords')
        time = obsdata['time']
        lwdn = obsdata['lwdn'].load()
        lwup = obsdata['lwup'].load()
        swdn = obsdata['swdn'].load()
        swup = obsdata['swup'].load()
        obsdata.close()
    if dt < 3600:
        lst = glob.glob(os.path.join(radfluxpath, '*.nc'))
        obsdata = xr.open_mfdataset(lst, combine='by_coords')
        time = obsdata['time']
        lwdn = obsdata['downwelling_longwave'].load()
        lwup = obsdata['upwelling_longwave'].load()
        swdn = obsdata['downwelling_shortwave'].load()
        swup = obsdata['upwelling_shortwave'].load()
        obsdata.close()    
                
    #%% re-shape the data into coarser resolution
    
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    lwdn_new = avg_time_1d(time, lwdn, time_new)
    swdn_new = avg_time_1d(time, swdn, time_new)
    lwup_new = avg_time_1d(time, lwup, time_new)
    swup_new = avg_time_1d(time, swup, time_new)
        
    #%% output file
    outfile = predatapath + 'sfc_radiation_ENA_'+year+'.nc'
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
    
    if dt >= 3600:
        ds.attrs["title"] = 'surface radiative flux data from ARMBE hourly data'
    if dt < 3600:
        ds.attrs["title"] = 'surface radiative flux data from RADFLUX data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_totcld(armbepath, arsclbndpath, tsipath, predatapath, year, dt=300):
    """
    prepare total cloud fraction data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    arsclbndpath : char
        input datapath.
    tsipath : char
        input datapath.
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

    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
  
    #%% read in data
    if dt >= 3600:
        lst = glob.glob(os.path.join(armbepath, '*armbecldrad*.nc'))
        obsdata = xr.open_mfdataset(lst, combine='by_coords')
        time = obsdata['time']
        cf_arscl = obsdata['tot_cld']
        qc_cf_arscl = obsdata['qc_tot_cld']
        cf_tsi = obsdata['tot_cld_tsi']
        qc_cf_tsi = obsdata['qc_tot_cld_tsi']
        # cf_visst = obsdata['cld_tot'] #this data is already preprocessed in the satellite program
        obsdata.close()

        cf_arscl.load()
        cf_tsi.load()
        # cf_visst.load()
        qc_cf_arscl.load()
        qc_cf_tsi.load()

        # quality controls. For ARMBE cloud fraction, <30% valid points within 1-hr window are flagged and removed
        cf_arscl[qc_cf_arscl>=2] = np.nan
        cf_tsi[qc_cf_tsi>=2] = np.nan

        # change unit from 1 to %
        cf_arscl = cf_arscl*100
        cf_tsi = cf_tsi*100

        cf_arscl_new = avg_time_1d(time, cf_arscl, time_new)
        cf_tsi_new = avg_time_1d(time, cf_tsi, time_new)
        # cf_visst_new = avg_time_1d(time, cf_visst, time_new)

    if dt < 3600:
        lst = glob.glob(os.path.join(arsclbndpath, '*.nc'))
        arscldata = xr.open_mfdataset(lst, combine='by_coords')
        arscltime = arscldata['time']
        cloud_base = arscldata['cloud_layer_base_height'][:,0] #first cloud layer base
        arscldata.close()

        lst = glob.glob(os.path.join(tsipath, '*.nc'))
        tsidata = xr.open_mfdataset(lst, combine='by_coords')
        cf_opaque_tsi = obsdata['percent_opaque']
        qc_cf_opaque_tsi = obsdata['qc_percent_opaque']
        cf_thin_tsi = obsdata['percent_thin']
        qc_cf_thin_tsi = obsdata['qc_percent_thin']
        tsidata.close()

        # compute ARSCL cloud fraction over dt
        cloud.load()
        arscl_dt = np.abs(arscltime[1].dt.second - arscltime[0].dt.second) # arscl time step
        arscl_nt = dt/arscl_dt # number of arscl timesteps in the time windoe
        cf_arscl_new = np.size(np.where(cloud > 0))/arscl_nt # fraction of times in window with cloud bases present
        # change unit from 1 to %
        cf_arscl = cf_arscl*100
      
        # average TSI cloud fraction over dt
        cf_opaque_tsi.load()
        qc_cf_opaque_tsi.load()
        cf_thin_tsi.load()
        qc_cf_thin_tsi.load()

        cf_tsi = cf_opaque_tsi + cf_thin_tsi # add opaque and thin cloud fractions to get total (should be in units of %)
        cf_tsi[qc_cf_opaque_tsi > 0] = np.nan
        cf_tsi[qc_cf_thin_tsi > 0] = np.nan
        cf_tsi_new = avg_time_1d(tsitime, cf_tsi, time_new)
        
    #%% output file
    outfile = predatapath + 'totcld_ENA_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'tot_cld_arscl': ('time', np.float32(cf_arscl_new)),
                    'tot_cld_tsi': ('time', np.float32(cf_tsi_new)),
                    # 'tot_cld_visst': ('time', np.float32(cf_visst_new))
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
    # ds['tot_cld_visst'].attrs["long_name"] = "cloud_area_fraction"
    # ds['tot_cld_visst'].attrs["units"] = "%"
    # ds['tot_cld_visst'].attrs["description"] = "Total cloud fraction based on VISST satellite product"
    
    if dt >= 3600:
        ds.attrs["title"] = 'total cloud fraction from ARMBE hourly data (from ARSCL and TSI)'
    if dt < 3600:
        ds.attrs["title"] = 'total cloud fraction from ARSCL and TSI'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Ndrop(ndroppath, predatapath, year, dt=300):
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
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time'].load()
    nd = obsdata['drop_number_conc'].load()
    qc_nd = obsdata['qc_drop_number_conc'].load()
    ctype = obsdata['cloud_base_type']
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
    outfile = predatapath + 'Ndrop_ENA_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'nd': ('time', np.float32(nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['nd'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['nd'].attrs["units"] = "cm-3"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from ARM retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd_ARMretrieval(mfrsrpath, arsclbndpath, mwrpath, predatapath, year, dt=300):
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
    arsclbndpath : char
        input datapath.
    mwrpath : char
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
    lst1 = glob.glob(os.path.join(arsclpath, 'enaarsclkazrbnd1kolliasC1.c0.'+year+'*.nc'))
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

    lst3 = glob.glob(os.path.join(mwrpath, '*.nc'))
    mwrdata = xr.open_mfdataset(files, combine='by_coords')
    
    mwrtime = mwrdata['time']
    lwp = mwrdata['phys_lwp']*1e-3 #kg/m**2
    qc_lwp = mwrdata['qc_phys_lwp']
  
    lwp.load()
    qc_lwp.load()
    cod.load()
    qc_cod.load()
    cbh.load()
    cbhs.load()
    cths.load()
    
    lwp = qc_mask_qcflag(lwp, qc_lwp)  # do not mask LWP with MFRSR input since clearsky (LWP=0) is flagged (should be okay in MWR datastream, but need to check, AV 6/24/2024)
    lwp[lwp < 0] = 0 #can be small negative values retrieved that should be set to 0
    cod = qc_mask_qcflag(cod,qc_cod)
    cbh = qc_remove_neg(cbh, remove_zero='True')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # calculate CDNC
    # calculate cloud depth for the lowest cloud layer in ARSCL
    cth = cths[:,0]
    H = cths[:,0] - cbh
    H.load()
    H = qc_remove_neg(H, remove_zero='true')
    H[cbhs[:,1] > 0] = np.nan #remove multiple cloud layers
    
    # filter some data samples
    H[cth>5000.] = np.nan   # remove deep clouds with cloud top > 5 km (AV: a better way to do this is by cloud top temperature because this will include cloud tops < 0C, especially in winter, but cloud top temperature requires interpolating from ARMBE or interpsonde to cloud top heights)
    lwp[np.logical_or(lwp < 0.02, lwp > 0.3)] = np.nan # may want to revisit upper limit on LWP
    cod[np.logical_or(cod < 4, cod > 60)] = np.nan # changed lower COD limit from 2 to 4 (AV 6/24/2024)
    
    # calculate CDNC first then average into 1hr
    # use 5-min inputs rather than 20 s, which was originally used; 5-min is closer to resolution of coarsest input and 20-s is noisy
    # time = mfrsrtime.data
    # H_tmp = np.interp(np.int64(time), np.int64(arscltime), H)
    time_5min = pd.date_range(start='2017-06-20', end='2018-02-21', freq=str(int(300))+"s") # make inputs every 5 min to avoid high frequency noise in nd retrieval (COD looks like it doesn't vary at timescales < ~5 min)
    H_5min = avg_time_1d(arscltime, H, time_5min)
    cod_5min = avg_time_1d(mfrsrtime, cod, time_5min)
    lwp_5min = avg_time_1d(mwrtime, lwp, time_5min
    nd = calc_cdnc_ARM(lwp_5min, cod_5min, H_5min)
    
    # exclude small values (AV removed this filter 6/24/2024)
    # nd[nd<10] = np.nan
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, nd, time_new)
        
    #%% output file
    outfile = predatapath + 'Nd_ARMretrieval_ENA_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'nd': ('time', np.float32(nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['nd'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['nd'].attrs["units"] = "cm-3"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from ARM retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["reference"] = 'https://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-140.pdf '
    ds.attrs["inputfile_sample"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd_WU(Nd_WUpath, predatapath, year, dt=300):
    """
    prepare cloud deoplet number concentration (Nd) and effective radius data for ENA
    retrieval algorithm developed by Peng Wu, et al., 2020
    reference: https://doi.org/10.1029/2019JD032205
    
    Parameters
    ----------
    Nd_WUpath : char
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
    time = np.array([],dtype='datetime64[s]')
    ndall = np.array([])
    reall = np.array([])
    lst = glob.glob(os.path.join(Nd_WUpath, 'ENA_micro_sfc_retrieval_'+year+'*.nc'))
    lst.sort()
    for fname in lst:
        f = Dataset(fname,'r')
        hr = f.variables['Time'][:]
        nc = f.variables['NC'][:]
        re0 = f.variables['RC'][:]
        re = np.nanmedian(re0, axis=0)
        day = fname[-11:-7]+'-'+fname[-7:-5]+'-'+fname[-5:-3]
        f.close()
        
        time1 = np.array([np.datetime64(day) + np.timedelta64(int(x*3600),'s') for x in hr])
        time = np.hstack((time, time1))
        ndall = np.hstack((ndall, nc))
        reall = np.hstack((reall, re))
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start=year+'-01-01', end=year+'-12-31 23:59:00', freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, ndall, time_new)
    re_new = median_time_1d(time, reall, time_new)
    
    #%% output file
    outfile = predatapath + 'Nd_Reff_Wu_etal_ENA_'+year+'.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'nd': ('time', np.float32(nd_new)),
                     'reff': ('time', np.float32(re_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['nd'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['nd'].attrs["units"] = "cm-3"
    ds['reff'].attrs["long_name"] = "cloud droplet effective radius"
    ds['reff'].attrs["units"] = "um"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from Peng Wu et al's retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["reference"] = 'https://doi.org/10.1029/2019JD032205'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
