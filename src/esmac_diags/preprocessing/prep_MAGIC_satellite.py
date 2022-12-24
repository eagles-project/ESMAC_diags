"""
prepare satellite data from MAGIC
only pixel-level data are available, average into 0.5x0.5 degree at the ship location
options of output data into coarser resolution
"""

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag
from esmac_diags.subroutines.time_format_change import datetime2cday
from esmac_diags.subroutines.specific_data_treatment import find_nearest, insolation, \
                calc_cdnc_VISST, calc_clouddepth_VISST

#%% test settings
# shipmetpath = '../../../data/MAGIC/obs/ship/magmarinemetM1.b1/'
# visstpixpath = '../../../data/MAGIC/obs/visst/pix/'
# predatapath = 'C:/Users/tang357/Downloads/prep_data/MAGIC/satellite/'
# shipmetpath = '../../../raw_data/obs/MAGIC/ship/magmarinemetM1.b1/'
# visstpixpath = '../../../raw_data/obs/MAGIC/visst/grid/'
# predatapath = '../../../prep_data/MAGIC/satellite/'
# dt=3600

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_VISST_pixel(shipmetpath, visstpixpath, predatapath, dt=3600):
    """
    prepare VISST-satellite data, average pixel level (4km) into 0.5x0.5 grid

    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    visstpixpath : str
        input datapath
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                               
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in location data
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    shiptime = shipdata['time'].load()
    lon0 = shipdata['lon_mean_gps'].load().data
    qc_lon = shipdata['qc_lon_mean_gps'].load().data
    lat0 = shipdata['lat_mean_gps'].load().data
    qc_lat = shipdata['qc_lat_mean_gps'].load().data
    shipdata.close()
    
    lat0 = qc_mask_qcflag(lat0, qc_lat)
    lon0 = qc_mask_qcflag(lon0, qc_lon)
    
    idx = np.logical_and(~np.isnan(lat0), ~np.isnan(lon0))
    shiptime = shiptime[idx]
    lon0 = lon0[idx]
    lat0 = lat0[idx]
    
    
    #%% read in satellite data
    lst = glob.glob(os.path.join(visstpixpath, 'magvisstpxg15*.cdf'))
    lst.sort()
    # first data
    visstdata = xr.open_dataset(lst[0])
    vissttime = visstdata['time_offset'].load()
    lat = visstdata['latitude'].data
    lon = visstdata['longitude'].data
    
    # get location of ship for VISST time
    lon_ship = np.interp(vissttime, shiptime, lon0)
    lat_ship = np.interp(vissttime, shiptime, lat0)
    lon_ship_all = np.array(lon_ship)
    lat_ship_all = np.array(lat_ship)
    
    # choose 0.5x0.5 degree around ship location
    idx = np.logical_and(np.abs(lon-lon_ship)<=0.25, np.abs(lat-lat_ship)<=0.25)
    idx = idx.flatten()   # change 2d array to 1-dimensional for index operation
    
    vis_reflectance0 = visstdata['reflectance_vis'].data.flatten()
    wp0 = visstdata['cloud_lwp_iwp'].data.flatten()
    phase0 = visstdata['cloud_phase'].data.flatten() #0:snow, 1:water, 2:ice, 3:no retrieval, 4:clear, 5:bad data, 6:suspected water, 7:suspected ice, 13:cleaned data
    particle_size0 = visstdata['cloud_particle_size'].data.flatten()
    cod0 = visstdata['cloud_visible_optical_depth'].data.flatten()
    ctt0 = visstdata['cloud_top_temperature'].data.flatten()
    ctp0 = visstdata['cloud_top_pressure'].data.flatten()
    cth0 = visstdata['cloud_top_height'].data.flatten()
    bb_lw0 = visstdata['broadband_longwave_flux'].data.flatten()
    bb_sw_albedo0 = visstdata['broadband_shortwave_albedo'].data.flatten()
    visstdata.close()
    
    # cloud mask
    cld_mask = phase0 * np.nan
    cld_mask[phase0==4] = 0  # clear
    cld_mask[np.logical_or(phase0==1, phase0==2)] = 1 # water or ice clouds
    # cld_mask[np.logical_or(phase0==6, phase0==7)] = 1 # suspected water or ice clouds
    
    # grid-mean lwp and iwp
    lwp0 = np.array(wp0)
    iwp0 = np.array(wp0)
    lwp0[phase0!=1] = np.nan
    iwp0[phase0!=2] = np.nan
    # count clear-sky pixel for gridbox-mean LWP and IWP
    lwp0[np.logical_or(phase0==2, phase0==4)] = 0
    iwp0[np.logical_or(phase0==1, phase0==4)] = 0
    
    # calculate grid-mean value
    vis_reflectance = np.nanmean(vis_reflectance0[idx])
    lwp = np.nanmean(lwp0[idx])
    iwp = np.nanmean(iwp0[idx])
    particle_size = np.nanmean(particle_size0[idx])
    cod = np.nanmean(cod0[idx])
    ctt = np.nanmean(ctt0[idx])
    ctp = np.nanmean(ctp0[idx])
    cth = np.nanmean(cth0[idx])
    bb_lw = np.nanmean(bb_lw0[idx])
    bb_sw_albedo = np.nanmean(bb_sw_albedo0[idx])
    cldfrac = np.nanmean(cld_mask[idx])
    
    vissttime = vissttime.data
    
    #%% all other files
    for filename in lst[1:]:
        print(filename)
        visstdata = xr.open_dataset(filename)
        vissttime0 = visstdata['time_offset'].load()
        lat = visstdata['latitude'].data
        lon = visstdata['longitude'].data
        
        # get location of ship for VISST time
        lon_ship = np.interp(vissttime0, shiptime, lon0)
        lat_ship = np.interp(vissttime0, shiptime, lat0)
        
        # choose 0.5x0.5 degree around ship location
        idx = np.logical_and(np.abs(lon-lon_ship)<=0.25, np.abs(lat-lat_ship)<=0.25)
        idx = idx.flatten()   # change 2d array to 1-dimensional for index operation
        
        vis_reflectance0 = visstdata['reflectance_vis'].data.flatten()
        wp0 = visstdata['cloud_lwp_iwp'].data.flatten()
        phase0 = visstdata['cloud_phase'].data.flatten() #0:snow, 1:water, 2:ice, 3:no retrieval, 4:clear, 5:bad data, 6:suspected water, 7:suspected ice, 13:cleaned data
        particle_size0 = visstdata['cloud_particle_size'].data.flatten()
        cod0 = visstdata['cloud_visible_optical_depth'].data.flatten()
        ctt0 = visstdata['cloud_top_temperature'].data.flatten()
        ctp0 = visstdata['cloud_top_pressure'].data.flatten()
        cth0 = visstdata['cloud_top_height'].data.flatten()
        bb_lw0 = visstdata['broadband_longwave_flux'].data.flatten()
        bb_sw_albedo0 = visstdata['broadband_shortwave_albedo'].data.flatten()
        visstdata.close()
        
        # cloud mask
        cld_mask = phase0 * np.nan
        cld_mask[phase0==4] = 0  # clear
        cld_mask[np.logical_or(phase0==1, phase0==2)] = 1 # water or ice clouds
        # cld_mask[np.logical_or(phase0==6, phase0==7)] = 1 # suspected water or ice clouds
        
        # grid-mean lwp and iwp
        lwp0 = np.array(wp0)
        iwp0 = np.array(wp0)
        lwp0[phase0!=1] = np.nan
        iwp0[phase0!=2] = np.nan
        # count clear-sky pixel for gridbox-mean LWP and IWP
        lwp0[np.logical_or(phase0==2, phase0==4)] = 0
        iwp0[np.logical_or(phase0==1, phase0==4)] = 0
        
        # combine data with grid-mean value
        vissttime = np.hstack((vissttime, vissttime0.data))
        lon_ship_all = np.hstack((lon_ship_all, lon_ship))
        lat_ship_all = np.hstack((lat_ship_all, lat_ship))
        vis_reflectance = np.hstack((vis_reflectance, np.nanmean(vis_reflectance0[idx])))
        lwp = np.hstack((lwp, np.nanmean(lwp0[idx])))
        iwp = np.hstack((iwp, np.nanmean(iwp0[idx])))
        particle_size = np.hstack((particle_size, np.nanmean(particle_size0[idx])))
        cod = np.hstack((cod, np.nanmean(cod0[idx])))
        ctt = np.hstack((ctt, np.nanmean(ctt0[idx])))
        ctp = np.hstack((ctp, np.nanmean(ctp0[idx])))
        cth = np.hstack((cth, np.nanmean(cth0[idx])))
        bb_lw = np.hstack((bb_lw, np.nanmean(bb_lw0[idx])))
        bb_sw_albedo = np.hstack((bb_sw_albedo, np.nanmean(bb_sw_albedo0[idx])))
        cldfrac = np.hstack((cldfrac, np.nanmean(cld_mask[idx])))
            
    #%% calculate TOA SW flux from albedo
    print('calculate some variables:')
    # change time to calendar day
    calday = datetime2cday(vissttime)
    # calculate insolation
    ins = insolation(calday, lon_ship_all, lat_ship_all, leap_year='leap')
    
    # calculate net SW flux
    bb_sw = ins * (1 - bb_sw_albedo*0.01)
    
    cldfrac = cldfrac * 100 # unit: %
    
    #%% retrieve CDNC
    H = calc_clouddepth_VISST(lwp*0.001, ctt, adiabaticity=0.8)
    H_ad = calc_clouddepth_VISST(lwp*0.001, ctt, adiabaticity=1.0)
    Nd = calc_cdnc_VISST(lwp*0.001, ctt, cod, adiabaticity=0.8)
    Nd_ad = calc_cdnc_VISST(lwp*0.001, ctt, cod, adiabaticity=1.0)
    
    #filter out columns with ice and bad retrievals
    H_array = np.array(H)
    H_ad_array = np.array(H_ad)
    Nd_array = np.array(Nd)
    Nd_ad_array = np.array(Nd_ad)
    
    ind = np.isinf(Nd_array)
    H_array[ind] = np.nan
    H_ad_array[ind] = np.nan
    Nd_array[ind] = np.nan
    Nd_ad_array[ind] = np.nan
    
    ind = np.array(iwp >= 1)
    H_array[ind] = np.nan
    H_ad_array[ind] = np.nan
    Nd_array[ind] = np.nan
    Nd_ad_array[ind] = np.nan
    
    # effective radius
    reff = np.array(particle_size)
    reff[ind] = np.nan
    
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2012-10-05', end='2013-10-09 23:59:00', freq=str(int(dt))+"s")  # MAGIC time period
    
    Nd_new = avg_time_1d(vissttime, Nd_array, time_new)
    H_new = avg_time_1d(vissttime, H, time_new)
    lwp_new = avg_time_1d(vissttime, lwp, time_new)
    iwp_new = avg_time_1d(vissttime, iwp, time_new)
    reff_new = avg_time_1d(vissttime, reff, time_new)
    cod_new = avg_time_1d(vissttime, cod, time_new)
    ctt_new = avg_time_1d(vissttime, ctt, time_new)
    ctp_new = avg_time_1d(vissttime, ctp, time_new)
    cth_new = avg_time_1d(vissttime, cth, time_new)
    lw_new = avg_time_1d(vissttime, bb_lw, time_new)
    sw_new = avg_time_1d(vissttime, bb_sw, time_new)
    albedo_new = avg_time_1d(vissttime, bb_sw_albedo, time_new)
    cldfrac_new = avg_time_1d(vissttime, cldfrac, time_new)
    
    #%% output file
    outfile = predatapath + 'Nd_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'Nd': (['time'], np.float32(Nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['Nd'].attrs["long_name"] = 'cloud droplet number concentration'
    ds['Nd'].attrs["units"] = '#/cm3'
    
    ds.attrs["title"] = 'cloud droplet number concentration retrieved average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["description"] = 'retrieved following Bennartz 2007, assuming adiabaticity = 0.8'
    ds.attrs["reference"] = 'https://doi.org/10.1029/2006JD007547'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'Hcld_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'Hcld': (['time'], np.float32(H_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['Hcld'].attrs["long_name"] = 'cloud depth for liquid cloud only'
    ds['Hcld'].attrs["units"] = 'm'
    
    ds.attrs["title"] = 'liquid cloud depth retrieved average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["description"] = 'retrieved from LWP and cloud top temperature assuming adiabaticity = 0.8'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'LWP_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwp': (['time'], np.float32(lwp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwp'].attrs["long_name"] = 'liquid water path'
    ds['lwp'].attrs["units"] = 'g/m2'
    
    ds.attrs["title"] = 'liquid water path average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'IWP_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'iwp': (['time'], np.float32(iwp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['iwp'].attrs["long_name"] = 'ice water path'
    ds['iwp'].attrs["units"] = 'g/m2'
    
    ds.attrs["title"] = 'ice water path average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'Reff_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'reff': (['time'], np.float32(reff_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['reff'].attrs["long_name"] = 'effective radius for liquid clouds'
    ds['reff'].attrs["units"] = 'um'
    
    ds.attrs["title"] = 'liquid clouds effective radius average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'cod_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cod': (['time'], np.float32(cod_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cod'].attrs["long_name"] = 'cloud optical depth for liquid clouds'
    ds['cod'].attrs["units"] = 'N/A'
    
    ds.attrs["title"] = 'liquid clouds opical depth average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    outfile = predatapath + 'cloudfraction_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cldtot': (['time'], np.float32(cldfrac_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cldtot'].attrs["long_name"] = 'cloud fraction for all heights'
    ds['cldtot'].attrs["units"] = '%'
    
    ds.attrs["title"] = 'cloud fraction calculated from VISST pixel-level (4km) data in 0.5x0.5 grid'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    ds.to_netcdf(outfile, mode='w')
        
    #
    outfile = predatapath + 'cloudtop_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'ctt': (['time'], np.float32(ctt_new)),
                    'cth': (['time'], np.float32(cth_new)),
                    'ctp': (['time'], np.float32(ctp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['ctt'].attrs["long_name"] = 'cloud top temperature'
    ds['ctt'].attrs["units"] = 'K'
    ds['ctp'].attrs["long_name"] = 'cloud top pressure'
    ds['ctp'].attrs["units"] = 'hPa'
    ds['cth'].attrs["long_name"] = 'cloud top height'
    ds['cth'].attrs["units"] = 'km'
    
    ds.attrs["title"] = 'cloud top temperature, pressure and height average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["description"] = 'for any cloud type'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'lwflx_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwnettoa': (['time'], np.float32(lw_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwnettoa'].attrs["long_name"] = 'net LW flux at TOA'
    ds['lwnettoa'].attrs["units"] = 'W/m2'
    
    ds.attrs["title"] = 'net longwave flux at TOA average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["description"] = 'upward positive'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'swflx_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'swnettoa': (['time'], np.float32(sw_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['swnettoa'].attrs["long_name"] = 'net SW flux at TOA'
    ds['swnettoa'].attrs["units"] = 'W/m2'
    
    ds.attrs["title"] = 'net shortwave flux at TOA average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["description"] = 'calculate from insolation and SW albedo, downward positive'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    outfile = predatapath + 'albedo_VISSTgrid_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'albedo': (['time'], np.float32(albedo_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['albedo'].attrs["long_name"] = 'broadband_shortwave_albedo at TOA'
    ds['albedo'].attrs["units"] = '%'
    
    ds.attrs["title"] = 'broadband_shortwave_albedo at TOA average from VISST pixel-level (4km) data to grid (0.5x0.5 degree)'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
