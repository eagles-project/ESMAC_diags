"""
prepare ship data from MARCUS
options of output data into coarser resolution
"""

import glob
import os
import xarray as xr
import pandas as pd
import numpy as np
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, median_time_1d, median_time_2d, interp_time_1d
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag
from esmac_diags.subroutines.specific_data_treatment import calc_cldfrac_from_highres

# shipmetpath = '../../../../ESMAC_Diags_Tool/data/MARCUS/obs/ship/maraadmetX1.b1/'
# mwrpath = '../../../../ESMAC_Diags_Tool/data/MARCUS/obs/ship/marmwrret1liljclouM1.s2/'
# cpcpath = '../../../../ESMAC_Diags_Tool/data/MARCUS/obs/ship/maraoscpcf1mM1.b1/'
# ccnpath = '../../../../ESMAC_Diags_Tool/data/MARCUS/obs/ship/maraosccn1colavgM1.b1/'
# uhsaspath = '../../../../ESMAC_Diags_Tool/data/MARCUS/obs/ship/maraosuhsasM1.a1/'
# exhaustfreepath = '../../../../ESMAC_Diags_Tool/data/MARCUS/obs/ship/ship_exhaustfree/'
# prep_data_path = '../../../prep_data/MARCUS/ship/'

# dt=3600


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CCN(shipmetpath, ccnpath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    ccnpath : str
        input path of CCN data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    # lst2 = glob.glob(ccnpath+'maraosccn1colavgM1.b1.*')
    lst2 = glob.glob(ccnpath+'maraosccn1colspectraM1.b1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    ccn = obsdata['N_CCN'].load()
    qc_ccns = obsdata['qc_N_CCN'].load()
    ss = obsdata['supersaturation_calculated'].load()
    coefs = obsdata['N_CCN_fit_coefs'].load()
    obsdata.close()
    
    # ccn = qc_mask_qcflag(ccn, qc_ccn)
    
    # ccn_1s = np.array(ccn)
    # ccn_2s = np.array(ccn)
    # ccn_3s = np.array(ccn)
    # ccn_5s = np.array(ccn)
    # ccn_6s = np.array(ccn)
    # ccn_1s[np.logical_or(SS<0.05, SS>0.15)] = np.nan
    # ccn_2s[np.logical_or(SS<0.15, SS>0.25)] = np.nan
    # ccn_3s[np.logical_or(SS<0.25, SS>0.35)] = np.nan
    # ccn_5s[np.logical_or(SS<0.45, SS>0.55)] = np.nan
    # ccn_6s[np.logical_or(SS<0.55, SS>0.65)] = np.nan

    ss1 = ss[:,1]
    ss2 = ss[:,2]
    ss5 = ss[:,3]
    ss8 = ss[:,4]
    ss10 = ss[:,5]
    ccn1 = ccn[:,1]
    ccn2 = ccn[:,2]
    ccn5 = ccn[:,3]
    ccn8 = ccn[:,4]
    ccn10 = ccn[:,5]

    #%% these are computed from CCN spectra polynomial fits
    #this accounts for fluctuations in supersaturation that are different than the target supersaturation
    #but the fits do not always work, so the the sample size is less than the measured CCN
    ccn1_fit = coefs[:,0] + coefs[:,1]*0.1 + coefs[:,2]*(0.1**2)
    ccn2_fit = coefs[:,0] + coefs[:,1]*0.2 + coefs[:,2]*(0.2**2)
    ccn5_fit = coefs[:,0] + coefs[:,1]*0.5 + coefs[:,2]*(0.5**2)
    ccn8_fit = coefs[:,0] + coefs[:,1]*0.8 + coefs[:,2]*(0.8**2)
    ccn10_fit = coefs[:,0] + coefs[:,1]*1.0 + coefs[:,2]*(1.0**2)
  
    #apply basic QC flags
    ccn1 = qc_mask_qcflag(ccn1, qc_ccns[:,0])
    ccn2 = qc_mask_qcflag(ccn2, qc_ccns[:,1])
    ccn5 = qc_mask_qcflag(ccn5, qc_ccns[:,2])
    ccn8 = qc_mask_qcflag(ccn5, qc_ccns[:,3])
    ccn10 = qc_mask_qcflag(ccn5, qc_ccns[:,4])

    #apply to ccn fits
    ccn1_fit = qc_mask_qcflag(ccn1_fit, qc_ccns[:,0])
    ccn2_fit = qc_mask_qcflag(ccn2_fit, qc_ccns[:,1])
    ccn5_fit = qc_mask_qcflag(ccn5_fit, qc_ccns[:,2])
    ccn8_fit = qc_mask_qcflag(ccn8_fit, qc_ccns[:,3])
    ccn10_fit = qc_mask_qcflag(ccn10_fit, qc_ccns[:,4])
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period

    if dt >= 3600:
        lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
        lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
        ccn1_measure = median_time_1d(time2, ccn1, time_new, arraytype='xarray')
        ss1_i = median_time_1d(time2, ss1, time_new, arraytype='xarray')
        ccn2_measure = median_time_1d(time2, ccn2, time_new, arraytype='xarray')
        ss2_i = median_time_1d(time2, ss2, time_new, arraytype='xarray')
        ccn5_measure = median_time_1d(time2, ccn5, time_new, arraytype='xarray')
        ss5_i = median_time_1d(time2, ss5, time_new, arraytype='xarray')
        ccn8_measure = median_time_1d(time2, ccn8, time_new, arraytype='xarray')
        ss8_i = median_time_1d(time2, ss8, time_new, arraytype='xarray')
        ccn10_measure = median_time_1d(time2, ccn10, time_new, arraytype='xarray')
        ss10_i = median_time_1d(time2, ss10, time_new, arraytype='xarray')
        ccn1_fit_i = median_time_1d(time2, ccn1_fit, time_new, arraytype='xarray')
        ccn2_fit_i = median_time_1d(time2, ccn2_fit, time_new, arraytype='xarray')
        ccn5_fit_i = median_time_1d(time2, ccn5_fit, time_new, arraytype='xarray')   
        ccn8_fit_i = median_time_1d(time2, ccn8_fit, time_new, arraytype='xarray')  
        ccn10_fit_i = median_time_1d(time2, ccn10_fit, time_new, arraytype='xarray')  
    if dt < 3660:
        lon1 = interp_time_1d(time, lon, time_new, arraytype='xarray')
        lat1 = interp_time_1d(time, lat, time_new, arraytype='xarray')
        ccn1_measure = interp_time_1d(time2, ccn1, time_new, arraytype='xarray')
        ss1_i = interp_time_1d(time2, ss1, time_new, arraytype='xarray')
        ccn2_measure = interp_time_1d(time2, ccn2, time_new, arraytype='xarray')
        ss2_i = interp_time_1d(time2, ss2, time_new, arraytype='xarray')
        ccn5_measure = interp_time_1d(time2, ccn5, time_new, arraytype='xarray')
        ss5_i = interp_time_1d(time2, ss5, time_new, arraytype='xarray')
        ccn8_measure = interp_time_1d(time2, ccn8, time_new, arraytype='xarray')
        ss8_i = interp_time_1d(time2, ss8, time_new, arraytype='xarray')
        ccn10_measure = interp_time_1d(time2, ccn10, time_new, arraytype='xarray')
        ss10_i = interp_time_1d(time2, ss10, time_new, arraytype='xarray')
        ccn1_fit_i = interp_time_1d(time2, ccn1_fit, time_new, arraytype='xarray')
        ccn2_fit_i = interp_time_1d(time2, ccn2_fit, time_new, arraytype='xarray')
        ccn5_fit_i = interp_time_1d(time2, ccn5_fit, time_new, arraytype='xarray')   
        ccn8_fit_i = interp_time_1d(time2, ccn8_fit, time_new, arraytype='xarray')  
        ccn10_fit_i = interp_time_1d(time2, ccn10_fit, time_new, arraytype='xarray')  
    
    #%% output file
    outfile = prep_data_path + 'CCN_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                     'CCN1_fit': ('time', np.float32(ccn1_fit_i)),
                     'CCN2_fit': ('time', np.float32(ccn2_fit_i)),
                     'CCN5_fit': ('time', np.float32(ccn5_fit_i)),
                     'CCN8_fit': ('time', np.float32(ccn8_fit_i)),
                     'CCN10_fit': ('time', np.float32(ccn10_fit_i)),
                     'CCN1': ('time', np.float32(ccn1_measure)),
                     'CCN2': ('time', np.float32(ccn2_measure)),
                     'CCN5': ('time', np.float32(ccn5_measure)),
                     'CCN8': ('time', np.float32(ccn8_measure)),
                     'CCN10': ('time', np.float32(ccn10_measure)),
                     'ss1': ('time', np.float32(ss1_i)),
                     'ss2': ('time', np.float32(ss2_i)),
                     'ss5': ('time', np.float32(ss5_i)),
                     'ss8': ('time', np.float32(ss8_i)),
                     'ss10': ('time', np.float32(ss10_i)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CCN1_fit'].attrs["long_name"] = "0.1% Cloud Condensation Nuclei"
    ds['CCN1_fit'].attrs["units"] = "cm-3"
    ds['CCN1_fit'].attrs["description"] = "Calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['CCN2_fit'].attrs["long_name"] = "0.2% Cloud Condensation Nuclei"
    ds['CCN2_fit'].attrs["units"] = "cm-3"
    ds['CCN2_fit'].attrs["description"] = "Calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['CCN5_fit'].attrs["long_name"] = "0.5% Cloud Condensation Nuclei"
    ds['CCN5_fit'].attrs["units"] = "cm-3"
    ds['CCN5_fit'].attrs["description"] = "Calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['CCN8_fit'].attrs["long_name"] = "0.8% Cloud Condensation Nuclei"
    ds['CCN8_fit'].attrs["units"] = "cm-3"
    ds['CCN8_fit'].attrs["description"] = "Calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['CCN10_fit'].attrs["long_name"] = "1.0% Cloud Condensation Nuclei"
    ds['CCN10_fit'].attrs["units"] = "cm-3"
    ds['CCN10_fit'].attrs["description"] = "Calculated using a polynomial fit to ARM-measured CCN spectra"
    ds['CCN1'].attrs["long_name"] = "0.1% Cloud Condensation Nuclei - measured"
    ds['CCN1'].attrs["units"] = "cm-3"
    ds['CCN1'].attrs["description"] = "ARM-measured CCN targeted to 0.1% SS. see SS1 for actual measured SS"
    ds['ss1'].attrs["long_name"] = "Actual Supersaturation targeted to 0.1%"
    ds['ss1'].attrs["units"] = "%"
    ds['ss1'].attrs["description"] = "measured SS that is closest to 0.1%. ccn1_m is measured at this SS"
    ds['CCN2'].attrs["long_name"] = "0.2% Cloud Condensation Nuclei"
    ds['CCN2'].attrs["units"] = "cm-3"
    ds['CCN2'].attrs["description"] = "ARM-measured CCN targeted to 0.2% SS. see SS2 for actual measured SS"
    ds['ss2'].attrs["long_name"] = "Actual Supersaturation targeted to 0.2%"
    ds['ss2'].attrs["units"] = "%"
    ds['ss2'].attrs["description"] = "measured SS that is closest to 0.2%. ccn2_m is measured at this SS"
    ds['CCN5'].attrs["long_name"] = "0.5% Cloud Condensation Nuclei"
    ds['CCN5'].attrs["units"] = "cm-3"
    ds['CCN5'].attrs["description"] = "ARM-measured CCN targeted to 0.5% SS. see SS5 for actual measured SS"
    ds['ss5'].attrs["long_name"] = "Actual Supersaturation targeted to 0.5%"
    ds['ss5'].attrs["units"] = "%"
    ds['ss5'].attrs["description"] = "measured SS that is closest to 0.5%. ccn5_m is measured at this SS"
    ds['CCN8'].attrs["long_name"] = "0.8% Cloud Condensation Nuclei"
    ds['CCN8'].attrs["units"] = "cm-3"
    ds['CCN8'].attrs["description"] = "ARM-measured CCN targeted to 0.8% SS. see SS8 for actual measured SS"
    ds['ss8'].attrs["long_name"] = "Actual Supersaturation targeted to 0.8%"
    ds['ss8'].attrs["units"] = "%"
    ds['ss8'].attrs["description"] = "measured SS that is closest to 0.8%. ccn5_m is measured at this SS"
    ds['CCN10'].attrs["long_name"] = "1.0% Cloud Condensation Nuclei"
    ds['CCN10'].attrs["units"] = "cm-3"
    ds['CCN10'].attrs["description"] = "ARM-measured CCN targeted to 1.0% SS. see SS10 for actual measured SS"
    ds['ss10'].attrs["long_name"] = "Actual Supersaturation targeted to 1.0%"
    ds['ss10'].attrs["units"] = "%"
    ds['ss10'].attrs["description"] = "measured SS that is closest to 1.0%. ccn5_m is measured at this SS"
    
    ds.attrs["title"] = 'Surface CCN number concentration'
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    if dt >= 3600:
        ds.attrs["description"] = 'median value of each time window'
    if dt < 3600:
        ds.attrs["description"] = 'interpolated value from ~hourly resolution data'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CCN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
exhaustfreepath : str
    input path of exhaust-free data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    

    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    filename_exhaustfree = exhaustfreepath+'CCN_exhaustfree_1hr.nc'
    obsdata = xr.open_dataset(filename_exhaustfree)
    time2 = obsdata['time'].load()
    ccn1s = obsdata['CCN1'].load()
    ccn2s = obsdata['CCN2'].load()
    ccn5s = obsdata['CCN5'].load()
    obsdata.close()
    
    ccn1s = qc_remove_neg(ccn1s)
    ccn2s = qc_remove_neg(ccn2s)
    ccn5s = qc_remove_neg(ccn5s)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period

    if dt >= 3600:
        lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
        lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
        ccn1 = median_time_1d(time2, ccn1s, time_new, arraytype='xarray')
        ccn2 = median_time_1d(time2, ccn2s, time_new, arraytype='xarray')
        ccn5 = median_time_1d(time2, ccn5s, time_new, arraytype='xarray')
    if dt < 3600:
        lon1 = interp_time_1d(time, lon, time_new, arraytype='xarray')
        lat1 = interp_time_1d(time, lat, time_new, arraytype='xarray')
        ccn1 = interp_time_1d(time2, ccn1s, time_new, arraytype='xarray')
        ccn2 = interp_time_1d(time2, ccn2s, time_new, arraytype='xarray')
        ccn5 = interp_time_1d(time2, ccn5s, time_new, arraytype='xarray')

    #%% output file
    outfile = prep_data_path + 'CCN_MARCUS_exhaustfree.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1),
                    'lon': (['time'], lon1),
                    'CCN1': (['time'], ccn1.data),
                    'CCN2': (['time'], ccn2.data),
                    'CCN5': (['time'], ccn5.data),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CCN1'].attrs["long_name"] = ccn1s.long_name
    ds['CCN1'].attrs["description"] = ccn1s.description
    ds['CCN1'].attrs["units"] = "1/cm3"
    ds['CCN2'].attrs["long_name"] = ccn2s.long_name
    ds['CCN2'].attrs["description"] = ccn2s.description
    ds['CCN2'].attrs["units"] = "1/cm3"
    ds['CCN5'].attrs["long_name"] = ccn5s.long_name
    ds['CCN5'].attrs["description"] = ccn5s.description
    ds['CCN5'].attrs["units"] = "1/cm3"

    ds.attrs["title"] = 'Surface CCN number concentration'
    ds.attrs["input data_example"] = filename_exhaustfree.split('/')[-1]
    if dt >= 3600:
        ds.attrs["description"] = 'median value of each time window'
    if dt < 3600:
        ds.attrs["description"] = 'interpolated value from ~hourly resolution data'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CN(shipmetpath, cpcpath, uhsaspath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    cpcpath : str
        input path of CPC data
    uhsaspath : str
        input path of UHSAS data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    lst1 = sorted(glob.glob(cpcpath+'maraoscpcf1mM1.b1.*'))
    obsdata = xr.open_mfdataset(lst1, combine="nested", concat_dim="time")
    time1 = obsdata['time'].load()
    cpc = obsdata['concentration'].load()
    qc_cpc = obsdata['qc_concentration'].load()
    obsdata.close()

    cpc = cpc.drop_duplicates(dim="time")
    time1 = time1.drop_duplicates(dim="time")
    qc_cpc = qc_cpc.drop_duplicates(dim="time")

    cpc = qc_remove_neg(cpc.data)
    cpc = qc_mask_qcflag(cpc, qc_cpc)
    
    lst2 = sorted(glob.glob(uhsaspath+'maraosuhsasM1.a1.*'))
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    uhsas = obsdata['concentration'].load()
    dmin = obsdata['lower_size_limit'][0,:].load()
    dmax = obsdata['upper_size_limit'][0,:].load()
    obsdata.close()
    
    uhsas = qc_remove_neg(uhsas.data)
    uhsas100 = np.nansum(uhsas[:, dmin>100], axis=1)
    uhsas100 = qc_remove_neg(uhsas100, remove_zero='True')
    
    #%% re-shape the data into coarser resolution (ship position starts at 3Z on 10-21)
    time_new = pd.date_range(start='2017-10-21 03:00:00', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period

    tmpcpc = xr.DataArray(data=np.array(cpc), dims=["time"], coords=dict(time=time1))
    tmpuhsas100 = xr.DataArray(data=np.array(uhsas100), dims=["time"], coords=dict(time=time2))

    lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
    lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
    cpc1 = median_time_1d(time1, tmpcpc, time_new, arraytype='xarray')
    uhsas1 = median_time_1d(time2, tmpuhsas100, time_new, arraytype='xarray')
    
    #%% output file
    outfile = prep_data_path + 'CN_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'CPC10': (['time'], cpc1.data),
                    'UHSAS100': (['time'], uhsas1.data),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CPC10'].attrs["long_name"] = 'CPC measured aerosol number (size>10nm)'
    ds['CPC10'].attrs["units"] = "1/cm3"
    ds['UHSAS100'].attrs["long_name"] = 'UHSAS measured aerosol number (size>100nm)'
    ds['UHSAS100'].attrs["units"] = "1/cm3"
    
    ds.attrs["input data_example"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    exhaustfreepath : str
        input path of exhaust-free data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    

    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    filename_exhaustfree = exhaustfreepath + 'CPC_UHSAS_exhaustfree_1hr.nc'
    obsdata = xr.open_dataset(filename_exhaustfree)
    time2 = obsdata['time'].load()
    cpc = obsdata['CPC'].load()
    uhsas = obsdata['UHSAS100'].load()
    obsdata.close()
    
    cpc = qc_remove_neg(cpc.data)
    uhsas = qc_remove_neg(uhsas.data)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period
    
    tmpcpc = xr.DataArray(data=np.array(cpc), dims=["time"], coords=dict(time=time2))
    tmpuhsas = xr.DataArray(data=np.array(uhsas), dims=["time"], coords=dict(time=time2))

    lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
    lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
    cpc1 = median_time_1d(time1, tmpcpc, time_new, arraytype='xarray')
    uhsas1 = median_time_1d(time2, tmpuhsas, time_new, arraytype='xarray')
    
    #%% output file
    outfile = prep_data_path + 'CN_MARCUS_exhaustfree.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'CPC10': (['time'], cpc1),
                    'UHSAS100': (['time'], uhsas1),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CPC10'].attrs["long_name"] = 'CPC measured aerosol number (size>10nm)'
    ds['CPC10'].attrs["units"] = "1/cm3"
    ds['UHSAS100'].attrs["long_name"] = 'UHSAS measured aerosol number (size>100nm)'
    ds['UHSAS100'].attrs["units"] = "1/cm3"
    
    ds.attrs["input data_example"] = filename_exhaustfree.split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    exhaustfreepath : str
        input path of exhaust-free UHSAS size distribution data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    

    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in UHSAS data
    filename_exhaustfree = exhaustfreepath + 'CPC_UHSAS_exhaustfree_1hr.nc'
    obsdata = xr.open_dataset(filename_exhaustfree)
    time2 = obsdata['time'].load()
    size = obsdata['size'].load()
    dmin = obsdata['size_low'].load()
    dmax = obsdata['size_high'].load()
    uhsas = obsdata['UHSAS'].load()
    obsdata.close()
    
    uhsas = qc_remove_neg(uhsas.data)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period
    
    tmpuhsas = xr.DataArray(data=np.array(uhsas), dims=["time"], coords=dict(time=time2))
    
    lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
    lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
    uhsas1 = median_time_2d(time2, tmpuhsas, time_new, arraytype='xarray')
    
    #%% output file
    outfile = prep_data_path + 'CNsize_UHSAS_MARCUS_exhaustfree.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'size_low': (['size'], dmin.data),
                    'size_high': (['size'], dmax.data),
                    'size_distribution_uhsas': (['time', 'size'], uhsas1),
                    },
                     coords={'time': ('time', time_new), 'size': ('size', size.data)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['size_low'].attrs["long_name"] = "lower bound of size bin"
    ds['size_low'].attrs["units"] = "nm"
    ds['size_high'].attrs["long_name"] = "upper bound of size bin"
    ds['size_high'].attrs["units"] = "nm"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['size'].attrs["description"] = "middle of bin: 0.5*(dmin+dmax)"
    ds['size_distribution_uhsas'].attrs["long_name"] = 'aerosol number size distribution'
    ds['size_distribution_uhsas'].attrs["units"] = '1/cm3'
    
    ds.attrs["input data_example"] = filename_exhaustfree.split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize(shipmetpath, uhsaspath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    uhsaspath : str
        input path for UHSAS size distribution data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    

    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in UHSAS data
    lst2 = glob.glob(uhsaspath+'maraosuhsasM1.a1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    uhsas = obsdata['concentration'].load()
    dmin = obsdata['lower_size_limit'][0,:].load()
    dmax = obsdata['upper_size_limit'][0,:].load()
    obsdata.close()
    
    uhsas = qc_remove_neg(uhsas.data)
    size = (dmin+dmax)/2
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period
    
    tmpuhsas = xr.DataArray(data=np.array(uhsas), dims=["time", "bin_num"], coords=dict(time=time2, bin_num=size))
    
    lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
    lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
    uhsas1 = median_time_2d(time2, tmpuhsas, time_new, arraytype='xarray')
    
    #%% output file
    outfile = prep_data_path + 'CNsize_UHSAS_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'size_low': (['size'], dmin.data),
                    'size_high': (['size'], dmax.data),
                    'size_distribution_uhsas': (['time', 'size'], uhsas1),
                    },
                     coords={'time': ('time', time_new), 'size': ('size', size.data)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['size_low'].attrs["long_name"] = "lower bound of size bin"
    ds['size_low'].attrs["units"] = "nm"
    ds['size_high'].attrs["long_name"] = "upper bound of size bin"
    ds['size_high'].attrs["units"] = "nm"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['size'].attrs["description"] = "middle of bin: 0.5*(dmin+dmax)"
    ds['size_distribution_uhsas'].attrs["long_name"] = 'aerosol number size distribution'
    ds['size_distribution_uhsas'].attrs["units"] = '1/cm3'
    
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_MET(shipmetpath, prep_data_path, dt=3600):
    """
    prepare basic meteorological fields
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)

    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    # T1 = shipdata['air_temperature_starboard'].load()
    # qc_T1 = shipdata['qc_air_temperature_starboard'].load()
    # RH1 = shipdata['relative_humidity_starboard'].load()
    # qc_RH1 = shipdata['qc_relative_humidity_starboard'].load()
    T = shipdata['air_temperature_port'].load()
    qc_T = shipdata['qc_air_temperature_port'].load()
    RH = shipdata['relative_humidity_port'].load()
    qc_RH = shipdata['qc_relative_humidity_port'].load()
    ps = shipdata['atmospheric_pressure'].load()
    qc_ps = shipdata['qc_atmospheric_pressure'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    T = qc_mask_qcflag(T, qc_T)
    RH = qc_mask_qcflag(RH, qc_RH)
    ps = qc_mask_qcflag(ps, qc_ps)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period
    
    lon1 = median_time_1d(time, lon, time_new, arraytype='xarray')
    lat1 = median_time_1d(time, lat, time_new, arraytype='xarray')
    T1 = median_time_1d(time, T, time_new, arraytype='xarray')
    RH1 = median_time_1d(time, RH, time_new, arraytype='xarray')
    ps1 = median_time_1d(time, ps, time_new, arraytype='xarray')
    
    #%% output file
    outfile = prep_data_path + 'T_RH_Ps_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'T': (['time'], T1),
                    'RH': (['time'], RH1),
                    'Ps': (['time'], ps1),
                    },
                      coords={'time': ('time', time_new), })
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['T'].attrs["long_name"] = T.long_name
    ds['T'].attrs["units"] = T.units
    ds['RH'].attrs["long_name"] = RH.long_name
    ds['RH'].attrs["units"] = RH.units
    ds['Ps'].attrs["long_name"] = ps.long_name
    ds['Ps'].attrs["units"] = ps.units
    
    ds.attrs["input data_example"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_MWR(shipmetpath, mwrpath, prep_data_path, dt=3600):
    """
    prepare LWP retrieved from MWR
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    mwrpath : str
        input path for MWR data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    

    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)

    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in MWR
    lst2 = glob.glob(mwrpath+'marmwrret1liljclouM1.s2.*')
    lst2.sort()
    
    # first data
    obsdata = xr.open_dataset(lst2[0])
    time2 = obsdata['time'].load()
    lwp = obsdata['be_lwp'].load()
    qc_lwp = obsdata['qc_be_lwp'].load()
    obsdata.close()
    for file in lst2[1:]:
        obsdata = xr.open_dataset(file)
        if file[-18:] == '20171223.000018.nc':
            t1 = obsdata['time'].load()
            lwp1 = obsdata['be_lwp'].load()
            qc_lwp1 = obsdata['qc_be_lwp'].load()
            t2 = xr.concat([t1[:358], t1[4082:], t1[366:4082]], dim="time")
            lwp2 = xr.concat([lwp1[:358], lwp1[4082:], lwp1[366:4082]], dim="time")
            qc_lwp2 = xr.concat([qc_lwp1[:358], qc_lwp1[4082:], qc_lwp1[366:4082]], dim="time")
            time2 = xr.concat([time2, t2], dim="time")
            lwp = xr.concat([lwp, lwp2], dim="time")
            qc_lwp = xr.concat([qc_lwp, qc_lwp2], dim="time")
        elif file[-18:] == '20171227.000000.nc':
            t1 = obsdata['time'].load()
            lwp1 = obsdata['be_lwp'].load()
            qc_lwp1 = obsdata['qc_be_lwp'].load()
            time2 = xr.concat([time2, t1[0:3417]], dim="time")
            lwp = xr.concat([lwp, lwp1[0:3417]], dim="time")
            qc_lwp = xr.concat([qc_lwp, qc_lwp1[0:3417]], dim="time")
        else:
            time2 = xr.concat([time2, obsdata['time']], dim="time")
            lwp = xr.concat([lwp, obsdata['be_lwp']], dim="time")
            qc_lwp = xr.concat([qc_lwp, obsdata['qc_be_lwp']], dim="time")
        obsdata.close()
    lwp = qc_mask_qcflag(lwp, qc_lwp)
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period
    
    lon1 = avg_time_1d(time, lon, time_new, arraytype='xarray')
    lat1 = avg_time_1d(time, lat, time_new, arraytype='xarray')
    lwp1 = avg_time_1d(time2, lwp, time_new, arraytype='xarray')
    lwp1 = qc_remove_neg(lwp1)
    
    #%% calculate cloud fraction from LWP
    # from MWR handbook: a value of LWP that is +/- 0.03 mm of zero could be clear sky
    lwp_thres = 30
    cf_out = calc_cldfrac_from_highres(lwp, time2, time_new, thres=lwp_thres)
    
    cf_5 = calc_cldfrac_from_highres(lwp, time2, time_new, thres=5)
    cf_10 = calc_cldfrac_from_highres(lwp, time2, time_new, thres=10)
    cf_20 = calc_cldfrac_from_highres(lwp, time2, time_new, thres=20)
    cf_30 = calc_cldfrac_from_highres(lwp, time2, time_new, thres=30)
    
    #%% output file
    outfile = prep_data_path + 'LWP_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'lwp': (['time'], lwp1),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['lwp'].attrs["long_name"] = lwp.long_name
    ds['lwp'].attrs["units"] = lwp.units
    
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'average into '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    outfile = prep_data_path + 'totcld_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'cldfrac': (['time'], cf_out),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['cldfrac'].attrs["long_name"] = 'total cloud fraction'
    ds['cldfrac'].attrs["units"] = '%'
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'calculated from LWP with threshold of '+str(lwp_thres)+' g/m2'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    ds.to_netcdf(outfile, mode='w')

    outfile = prep_data_path + 'totcld_sensitivity_MARCUS.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1.data),
                    'lon': (['time'], lon1.data),
                    'cldfrac_5': (['time'], cf_5),
                    'cldfrac_10': (['time'], cf_10),
                    'cldfrac_20': (['time'], cf_20),
                    'cldfrac_30': (['time'], cf_30),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['cldfrac_5'].attrs["long_name"] = 'total cloud fraction'
    ds['cldfrac_5'].attrs["units"] = '%'
    ds['cldfrac_10'].attrs["long_name"] = 'total cloud fraction'
    ds['cldfrac_10'].attrs["units"] = '%'
    ds['cldfrac_20'].attrs["long_name"] = 'total cloud fraction'
    ds['cldfrac_20'].attrs["units"] = '%'
    ds['cldfrac_30'].attrs["long_name"] = 'total cloud fraction'
    ds['cldfrac_30'].attrs["units"] = '%'
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'calculated from LWP with different LWP threshold (5/10/20/30 g/m2)'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    ds.to_netcdf(outfile, mode='w')

        
