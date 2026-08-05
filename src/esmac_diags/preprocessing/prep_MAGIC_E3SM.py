"""
prepare E3SM data for MAGIC ship tracks
options of average data into coarser resolution
"""

import glob
import os
import re
import fnmatch
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import time as ttt
from scipy.interpolate import interp1d
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d
from esmac_diags.subroutines.quality_control import  qc_mask_qcflag
from esmac_diags.subroutines.specific_data_treatment import find_nearest, calc_Reff_from_REL, \
                            calc_cdnc_ARM, calc_cdnc_VISST
from esmac_diags.subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm

# input_path = '../../../data/ACEENA/model/E3SMv1_out/'
# output_path = 'C:/Users/tang357/Downloads/MAGIC/'
# input_filehead = 'E3SMv1_MAGIC'
# output_filehead = 'E3SMv1_MAGIC'

# lev_out=np.arange(25.,1001,25.)
# height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
#                     1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
#                     3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
#                     10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

# dt = 3600

# shipmetpath = '../../../data/MAGIC/obs/ship/magmarinemetM1.b1/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_sfc(input_path, input2d_filehead, input3d_filehead, input3d_dryaerosol_filehead, input3d_cloudaerosol_filehead, input_cosp_filehead, output_path, output_filehead, dt=3600, config=config):
    """
    prepare surface (include TOA and vertical integrated) variables from E3SM output along ship tracks
    choose the grid nearest to the ship track location
    interpolate into new time with resolution of dt

    Parameters
    ----------
    input_path : str
        data path of E3SM output data.
    input_filehead : str
        filehead of the E3SM data, end at ".cam.h*" (E3SMv1) or ".eam.h*" (E3SMv2)
    output_path : str
        data path of output data..
    output_filehead : str
        filehead of the preprocessed E3SM data
    shipmetpath : str
        input path for ship location data
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    #%% settings specific for each site
    E3SMdomain_range = config['E3SMdomain_range']    # domain range in E3SM regional output
  
    # output time range and resolution
    time_new = pd.date_range(start='2012-10-05', end='2013-10-09 23:59:00', freq=str(int(dt))+"s")  # MAGIC time period
    
    #%% read in ship location data
    print(ttt.ctime(ttt.time())+': read in ship location data:')
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load().data
    lon = shipdata['lon_mean_gps'].load().data
    qc_lon = shipdata['qc_lon_mean_gps'].load().data
    lat = shipdata['lat_mean_gps'].load().data
    qc_lat = shipdata['qc_lat_mean_gps'].load().data
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    lon[lon<0] += 360
    
    lon_new = median_time_1d(time, lon, time_new)
    lat_new = median_time_1d(time, lat, time_new)
    
    #%% read in E3SM data
    print(ttt.ctime(ttt.time())+': read in E3SM data:')
    variable_names = list()
    variables = list()
    
    # lst = glob.glob(input_path + input_filehead+'.*.h?.*-00000.nc') 
    # lst.sort()
    lst3d = glob.glob(input_path + input3d_filehead+'.*2012-1?-??-*.nc') + \
            glob.glob(input_path + input3d_filehead+'.*2013-0?-??-*.nc') + \
            glob.glob(input_path + input3d_filehead+'.*2013-10-0?-*.nc') + \
            glob.glob(input_path + input3d_filehead+'.*2013-10-10-*.nc')
    lst3d_dryaer = glob.glob(input_path + input3d_dryaerosol_filehead+'.*2012-1?-??-*.nc') + \
            glob.glob(input_path + input3d_dryaerosol_filehead+'.*2013-0?-??-*.nc') + \
            glob.glob(input_path + input3d_dryaerosol_filehead+'.*2013-10-0?-*.nc') + \
            glob.glob(input_path + input3d_dryaerosol_filehead+'.*2013-10-10-*.nc')
    lst3d_wetaer = glob.glob(input_path + input3d_cloudaerosol_filehead+'.*2012-1?-??-*.nc') + \
            glob.glob(input_path + input3d_cloudaerosol_filehead+'.*2013-0?-??-*.nc') + \
            glob.glob(input_path + input3d_cloudaerosol_filehead+'.*2013-10-0?-*.nc') + \
            glob.glob(input_path + input3d_cloudaerosol_filehead+'.*2013-10-10-*.nc')
    lst2d = glob.glob(input_path + input2d_filehead+'.*2012-1?-??-*.nc') + \
            glob.glob(input_path + input2d_filehead+'.*2013-0?-??-*.nc') + \
            glob.glob(input_path + input2d_filehead+'.*2013-10-0?-*.nc') + \
            glob.glob(input_path + input2d_filehead+'.*2013-10-10-*.nc')
    lstcosp = glob.glob(input_path + input_cosp_filehead+'.*2012-1?-??-*.nc') + \
            glob.glob(input_path + input_cosp_filehead+'.*2013-0?-??-*.nc') + \
            glob.glob(input_path + input_cosp_filehead+'.*2013-10-0?-*.nc') + \
            glob.glob(input_path + input_cosp_filehead+'.*2013-10-10-*.nc')
    lst3d.sort()
    lst3d_dryaer.sort()
    lst3d_wetaer.sort()
    lst2d.sort()
    lstcosp.sort()
  
    # first data
    # e3smdata = xr.open_mfdataset(lst)
    # e3smtime = e3smdata.indexes['time'].to_datetimeindex()
    # lonm = e3smdata['lon'+'_'+E3SMdomain_range].load().data
    # latm = e3smdata['lat'+'_'+E3SMdomain_range].load().data
    # if len(lonm.shape)>1:   # when use xr.open_mfdataset, lat,lon will add a time dimension
    #     lonm = lonm[0,:]
    #     latm = latm[0,:]
    e3smdata3d = xr.open_dataset(lst3d[0])
    e3smdata3d = e3smdata3d.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
    dimensions3d = (config['time_dim'], config['vert_dim'], config['latlon_dim']+E3SMdomain_range)
    variables3d = []
    for var_name, var_data in e3smdata3d.data_vars.items():
        if set(var_data.dims) == set(dimensions3d):
            variables3d.append(var_name)
          
    e3smdata2d = xr.open_dataset(lst2d[0])
    e3smdata2d = e3smdata2d.transpose(config['time_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time and location
    dimensions2d = (config['time_dim'], config['latlon_dim']+E3SMdomain_range)
    variables2d = []
    for var_name, var_data in e3smdata2d.data_vars.items():
        if set(var_data.dims) == set(dimensions2d):
            variables2d.append(var_name)

    e3smdatacosp = xr.open_dataset(lstcosp[0])
    e3smdatacosp = e3smdatacosp.transpose(config['time_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time and location
    dimensions2d = (config['time_dim'], config['latlon_dim']+E3SMdomain_range)
    variablescosp = []
    for var_name, var_data in e3smdatacosp.data_vars.items():
        if set(var_data.dims) == set(dimensions2d):
            variablescosp.append(var_name)
          
    e3smtime = e3smdata3d.indexes[config['time_dim']].to_datetimeindex()
    lonm = e3smdata3d[config['LON']+E3SMdomain_range].load()
    latm = e3smdata3d[config['LAT']+E3SMdomain_range].load()
    
    lon_0 = median_time_1d(time, lon, e3smtime)
    lat_0 = median_time_1d(time, lat, e3smtime)
    x_idx_all = []
    for tt in range(len(e3smtime)):
        if np.isnan(lon_0[tt]):
            x_idx_all.append(-1)
        else:
            x_idx_tmp = find_nearest(lonm, latm, lon_0[tt], lat_0[tt])
            x_idx_all.append(x_idx_tmp)

    # getting column and levels
    len_ncol = len(e3smdata3d[config['latlon_dim']+E3SMdomain_range])
    len_lev = len(e3smdata3d[config['vert_dim']])
    
    # Get all simulated variables
    # vlist3d = list(e3smdata3d.variables.keys())
    vlist3d = variables3d
    av_vars3d = fnmatch.filter(vlist3d,'*'+E3SMdomain_range)

    # vlist2d = list(e3smdata2d.variables.keys())
    vlist2d = variables2d
    av_vars2d = fnmatch.filter(vlist2d,'*'+E3SMdomain_range)

    vlistcosp = variablescosp
    av_varscosp = fnmatch.filter(vlistcosp,'*'+E3SMdomain_range)

    # variables to calculate Reff and Nd
    print(ttt.ctime(ttt.time())+': Reff and Nd:')
    # z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
    # cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
    # z3 = xr.concat([z3[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # cloud = xr.concat([cloud[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # dz = (z3[:,:-2] - z3[:,2:])/2
    # dz = np.append(dz, (z3[:,-2:-1]+z3[:,-1:])/2, axis=1)
    # dz = np.insert(dz,0,dz[:,0],axis=1)
    # weight = cloud*dz
    # # mid-level T, P, z
    # Tmid = 0.5*(T[:,0:-1] + T[:,1:])
    # Pmid = 0.5*(Pres[:,0:-1] + Pres[:,1:])
    # Zmid = 0.5*(z3[:,0:-1] + z3[:,1:])
    # # cloud top
    # cf_sum_top = np.cumsum(cloud, axis=1)
    # cf_sum_top[cf_sum_top > 1] = 1
    # cf_sum_top_diff = cf_sum_top[:,1:] - cf_sum_top[:,0:-1]
    # z_cldtop = np.sum(Zmid[:,:]*cf_sum_top_diff[:,:], axis=1)
    # T_cldtop = np.sum(Tmid[:,:]*cf_sum_top_diff[:,:], axis=1)
    # # cloud base
    # cf_sum_base = np.cumsum(cloud[:, ::-1], axis=1)
    # cf_sum_base[cf_sum_base > 1] = 1
    # cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
    # z_cldbase = np.sum(Zmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # # normalize cloud fraction when not 100%
    # z_cldtop = z_cldtop / cf_sum_top[:,-1]
    # T_cldtop = T_cldtop / cf_sum_top[:,-1]
    # z_cldbase = z_cldbase / cf_sum_base[:,-1]
    # T_cldbase = T_cldbase / cf_sum_base[:,-1]
    # e3sm_cloud_depth = z_cldtop - z_cldbase
    # # # alternatively, one can assume a minimum cloud depth in which the full mass levels are used instead of half levels
    # # e3sm_z2 = z3[:,:-1]
    # # e3sm_z3 = z3[:,1:]
    # # e3sm_t2 = T[:,:-1]
    # # e3sm_t3 = T[:,1:]
    # # #cloud base and top calculations for minimum cloud depth
    # # z_cldbase2 = np.sum(e3sm_z2[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # # T_cldbase2 = np.sum(e3sm_t2[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # # z_cldtop2 = np.sum(e3sm_z3[:,:]*cf_sum_top_diff[:,:], axis=1)
    # # T_cldtop2 = np.sum(e3sm_t3[:,:]*cf_sum_top_diff[:,:], axis=1)
    # # e3sm_cloud_depth2 = z_cldtop2 - z_cldbase2
    # cloud_depth = xr.DataArray(data=e3sm_cloud_depth,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="cloud depth",units="m"),)
    # cbh = xr.DataArray(data=z_cldbase,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="cloud base height",units="m"),)
    # cth = xr.DataArray(data=z_cldtop,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="cloud top height",units="m"),)
    # cbt = xr.DataArray(data=T_cldbase,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="cloud base temperature",units="K"),)
    # ctt = xr.DataArray(data=T_cldtop,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="cloud top temperature",units="K"),)

    req_vlist = [config['Z'], config['CF']]
    req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars3d).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing cloud base, top, and depth')
        z3 = e3smdata3d[config['Z']+E3SMdomain_range][:,:,x_idx_all].load()
        if config['cloudfraction3d_legit'] == True:
            cloud = e3smdata3d[config['CF']+E3SMdomain_range][:,:,x_idx_all].load()
        else:
            cloud = e3smdata3d[config['CFLIQ']+E3SMdomain_range][:,:,x_idx_all].load() + e3smdata3d[config['CFICE']+E3SMdomain_range][:,:,x_idx_all].load()
            cloud = cloud.where(cloud <= 1, 1) # make sure CF doesn't exceed 100%
        dz = (z3[:,:-2].data - z3[:,2:].data)/2
        dz = np.append(dz, (z3[:,-2:-1].data+z3[:,-1:].data)/2, axis=1)
        dz = np.insert(dz,0,dz[:,0],axis=1)
        
        # mid-level T, P, z
        Tmid = 0.5*(T[:,0:-1].data + T[:,1:].data)
        Pmid = 0.5*(Pres[:,0:-1].data + Pres[:,1:].data)
        Zmid = 0.5*(z3[:,0:-1].data + z3[:,1:].data)
        # cloud top
        cf_sum_top = np.cumsum(cloud.data, axis=1)
        cf_sum_top[cf_sum_top > 1] = 1
        cf_sum_top_diff = cf_sum_top[:,1:] - cf_sum_top[:,0:-1]
        z_cldtop = np.sum(Zmid*cf_sum_top_diff, axis=1)
        T_cldtop = np.sum(Tmid*cf_sum_top_diff, axis=1)
        # cloud base
        cf_sum_base = np.cumsum(cloud[:, ::-1].data, axis=1)
        cf_sum_base[cf_sum_base > 1] = 1
        cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
        z_cldbase = np.sum(Zmid[:,::-1]*cf_sum_base_diff, axis=1)
        T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff, axis=1)
        # normalize cloud fraction when not 100%
        z_cldtop = z_cldtop / cf_sum_top[:,-1]
        T_cldtop = T_cldtop / cf_sum_top[:,-1]
        z_cldbase = z_cldbase / cf_sum_base[:,-1]
        T_cldbase = T_cldbase / cf_sum_base[:,-1]
        e3sm_cloud_depth = z_cldtop - z_cldbase

        cloud_depth = xr.DataArray(data=e3sm_cloud_depth,  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="cloud depth",units="m"),)
        cbh = xr.DataArray(data=z_cldbase,  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="cloud base height",units="m"),)
        cth = xr.DataArray(data=z_cldtop,  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="cloud top height",units="m"),)
        cbt = xr.DataArray(data=T_cldbase,  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="cloud base temperature",units="K"),)
        ctt = xr.DataArray(data=T_cldtop,  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="cloud top temperature",units="K"),)
    else:
        cloud_depth = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cbh = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cth = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cbt = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        ctt = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

    # # cloud optical depth and effective radius
    # rel = e3smdata['REL'+'_'+E3SMdomain_range].load()
    # freql = e3smdata['FREQL'+'_'+E3SMdomain_range].load()
    # icwnc = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
    # cod_a = e3smdata['TOT_CLD_VISTAU'+'_'+E3SMdomain_range].load()
    # cod_m = e3smdata['TAUWMODIS'+'_'+E3SMdomain_range].load()*0.01 #  cod from MODIS simulator, cloud fraction is treated as 1 but is 100
    # solin = e3smdata['SOLIN'+'_'+E3SMdomain_range].load()
    # rel = xr.concat([rel[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # freql = xr.concat([freql[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # icwnc = xr.concat([icwnc[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # cod_a = xr.concat([cod_a[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # cod_m = xr.concat([cod_m[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # solin = xr.concat([solin[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # # calculate mean effective radius. 
    # reff = calc_Reff_from_REL(rel, dz, freql, icwnc)
    # reff[reff==0] = np.nan
    # # calculate mean optical depth
    # cod = np.sum(cod_a,axis=1)
    # # cod[solin==0] = np.nan
    # reff_mean = xr.DataArray(data=reff,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
    # cod_mean = xr.DataArray(data=cod,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
    
    # # mean cloud droplet number concentration
    # cdnc_col = e3smdata['CDNUMC'+'_'+E3SMdomain_range].load()
    # cdnc_col = xr.concat([cdnc_col[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # cdnc_mean = cdnc_col/np.sum(weight,axis=1)
    # cdnc_mean[cdnc_mean >2e9] = np.nan
    # cdnc_mean = xr.DataArray(data=cdnc_mean,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
    
    # # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
    # lwp = e3smdata['TGCLDLWP'+'_'+E3SMdomain_range]
    # lwp = xr.concat([lwp[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
    # T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
    # nd_arm = calc_cdnc_ARM(lwp, cod_m, e3sm_cloud_depth)
    # nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m, adiabaticity=0.8)
    # cdnc_arm = xr.DataArray(data=nd_arm*1e6,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
    #                 description='Retrieved using ARM Ndrop algorithm'),)
    # cdnc_sat = xr.DataArray(data=nd_sat*1e6,  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
    #                 description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)

    if config['reff_output'] == True:
        req_vlist = [config['REL'], config['CFLIQ'], config['NC']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars3d).intersection(req_vlist))
  
        if len(matched_vlist) == len(req_vlist):
          print('\nAnalyzing effective radius')
          rel = e3smdata3d[config['REL']+E3SMdomain_range][:,:,x_idx_all].load()
          freql = e3smdata3d[config['CFLIQ']+E3SMdomain_range][:,:,x_idx_all].load()
          icwnc = e3smdata3d[config['NC']+E3SMdomain_range][:,:,x_idx_all].load()
  
          # calculate mean effective radius. 
          reff = calc_Reff_from_REL(rel.data, dz, freql.data, icwnc.data)
          reff[reff==0] = np.nan
    
          reff_mean = xr.DataArray(data=reff,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime)),
                attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
        else:
            reff_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='reff',attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    if config['tau3d_output'] == True:
        req_vlist3d = [config['TAU3D']]
        req_vlist2d = [config['SWDOWNTOA']]
        req_vlist3d = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist3d]
        req_vlist2d = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist2d]
        matched_vlist3d = list(set(av_vars3d).intersection(req_vlist3d))
        matched_vlist2d = list(set(av_vars2d).intersection(req_vlist2d))
        
        if np.logical_and(len(matched_vlist3d) == len(req_vlist3d), len(matched_vlist2d) == len(req_vlist2d)):
            print('\nAnalyzing cloud optical depth')
            cod_a = e3smdata3d[config['TAU3D']+E3SMdomain_range][:,:,x_idx_all].load()
            solin = e3smdata2d[config['SWDOWNTOA']+E3SMdomain_range][:,x_idx_all].load()
          
            # calculate mean optical depth
            cod = np.sum(cod_a.data,axis=1)
            # cod[solin==0] = np.nan        
            
            cod_mean = xr.DataArray(data=cod,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime)),
                attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
        else:
            cod_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cod',attrs={'units':'dummy_unit','long_name':'Dummy'})

    if config['cosp_output'] == True:   
        req_vlist = [config['TAULIQMODIS']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_varscosp).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing MODIS simulator cloud optical depth')
            cod_m = e3smdatacosp[config['TAULIQMODIS']+E3SMdomain_range][:,x_idx_all]
            cod_m[:] = e3smdatacosp[config['TAULIQMODIS']+E3SMdomain_range][:,x_idx_all].load()/e3smdatacosp[config['SUNLIT']+E3SMdomain_range][:,x_idx_all].values
        else:
            cod_m = xr.DataArray(np.zeros(len(e3smtime))*np.nan)
    
    # mean cloud droplet number concentration
    if config['colnc_output'] == True:
        req_vlist = [config['NC2D']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars2d).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing mean cloud droplet number concentration')
            cdnc_col = e3smdata3d[config['NC2D']+E3SMdomain_range][:,x_idx_all].load()
            weight = cloud*dz
            cdnc_mean = cdnc_col/weight.sum(dim=config['vert_dim'])
            cdnc_mean[cdnc_mean >2e9] = np.nan
            cdnc_mean = xr.DataArray(data=cdnc_mean,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
        else:
            cdnc_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    else:
        req_vlist = [config['NC']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars3d).intersection(req_vlist))
      
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing mean cloud droplet number concentration')
            #compute cloud layer mean CDNC from 3D NC (note that if NC is not in-cloud only, one needs to divide by cloud fraction)
            nc3d = e3smdata3d[config['NC']+E3SMdomain_range][:,:,x_idx_all].load()      
            if nc3d.attrs['units'] == '1/kg':
              rho = np.array(Pres/T/287.06)
              cdnc_rel = nc3d*rho/cloud
            if nc3d.attrs['units'] == 'm-3':
              cdnc_rel = nc3d/cloud
            cdnc_rel = cdnc_rel.where(cloud > 0, other = 0)
            weight = cloud*dz
            weight_column = weight.sum(dim=config['vert_dim'])
            cdnc_rel_avg = cdnc_rel.dot(weight, dims=config['vert_dim'])
            cdnc_mean = np.divide(cdnc_rel_avg, weight_column)
            
            cdnc_mean = xr.DataArray(data=cdnc_mean,  dims=[config['time_dim']],
                    coords=dict(time=([config['time_dim']], e3smtime)),
                    attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
        else:
          cdnc_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
    if config['tau3d_output'] == True:
        req_vlist = [config['LWP']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars2d).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing cloud droplet number concentration retrieved like Ndrop')
            lwp = e3smdata2d[config['LWP']+E3SMdomain_range][:,x_idx_all].data
            e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
            nd_arm = calc_cdnc_ARM(lwp, cod, e3sm_cloud_depth)
          
            cdnc_arm = xr.DataArray(data=nd_arm*1e6,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                            description='Retrieved using ARM Ndrop algorithm'),)
        else:
            cdnc_arm = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    elif config['tau3d_output'] == False and config['cosp_output'] == True:
        req_vlist = [config['LWP']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars2d).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing cloud droplet number concentration retrieved like Ndrop')
            lwp = e3smdata2d[config['LWP']+E3SMdomain_range][:,x_idx_all].data
            e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
            nd_arm = calc_cdnc_ARM(lwp, cod_m, e3sm_cloud_depth)
          
            cdnc_arm = xr.DataArray(data=nd_arm*1e6,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                            description='Retrieved using ARM Ndrop algorithm'),)
        else:
            cdnc_arm = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

    if config['cosp_output'] == True:
        # req_vlist = [config['LWP']]
        req_vlist = [config['LWPMODIS']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_varscosp).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing cloud droplet number concentration retrieved like MODIS')
            # lwp = e3smdata2d[config['LWP']+E3SMdomain_range][:,x_idx_all].data
            lwp = e3smdatacosp[config['LWPMODIS']+E3SMdomain_range][:,x_idx_all]
            lwp[:] = e3smdatacosp[config['LWPMODIS']+E3SMdomain_range][:,x_idx_all].load()/e3smdatacosp[config['SUNLIT']+E3SMdomain_range][:,x_idx_all].values
            T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
            nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m, adiabaticity=0.8)
          
            cdnc_sat = xr.DataArray(data=nd_sat*1e6,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                            description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)
        else:
            cdnc_sat = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

    # all other 2D (surface and vertical integrated) variables
    print(ttt.ctime(ttt.time())+': Other variables:')
    # variable2d_names = ['AODABS', 'AODALL', 'CDNUMC', 'CLDHGH', 'CLDMED', 'CLDLOW', 'CLDTOT', 
    #                     'FLDS', 'FLNS', 'FLNT', 'FLUT', 'FSDS', 'FSNS', 'FSNT', 'FSUTOA', 'SOLIN', 
    #                     'LHFLX', 'SHFLX', 'LWCF', 'SWCF', "TGCLDLWP", "TGCLDIWP", 
    #                     'IWPMODIS', 'LWPMODIS', 'REFFCLWMODIS', 'TAUIMODIS', 'TAUTMODIS', 'TAUWMODIS', 
    #                     #'MEANPTOP_ISCCP', 'MEANCLDALB_ISCCP', 'MEANTAU_ISCCP',
    #                     'PBLH', 'PRECT', 'PRECL', 'PRECC', 'PS', 'TREFHT', ]
    # for varname in variable2d_names:
    #     var = e3smdata[varname + '_'+E3SMdomain_range].load()
    #     var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
    #     if varname=='AODABS' or varname=='AODALL':
    #         var.attrs['units']='N/A'
    #     variable_names.append(varname)
    #     var2 = xr.concat([var[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    #     variables.append(var2)
    
    # # all other 3D (with vertical level) variables at the lowest model level
    # variable3d_names = ['CCN1', 'CCN3', 'CCN4', 'CCN5', 'Q', 'T', 'RELHUM', 'U', 'V'] 
    # for varname in variable3d_names:
    #     var = e3smdata[varname + '_'+E3SMdomain_range].load()
    #     var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
    #     var3 = xr.concat([var[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    #     variables.append(var3)
    #     variable_names.append(varname)
  
    variable2d_names = [config['AOD'], config['CLDTOT'], 
                        config['LWDOWNSFC'], config['LWUPTOA'],
                        config['SWDOWNSFC'], config['SWUPTOA'], config['SWDOWNTOA'], 
                        config['LWCF'], config['SWCF'], config['LWP'], config['IWP'], 
                        config['PBLH'], config['PS'], config['TS']]

    if config['netradiation_output'] == True:
        variable2d_names.append(config['LWNETSFC'])
        variable2d_names.append(config['LWNETTOA'])
        variable2d_names.append(config['SWNETSFC'])
        variable2d_names.append(config['SWNETTOA'])
    else:
        variable2d_names.append(config['LWUPSFC'])
        variable2d_names.append(config['LWDOWNTOA'])
        variable2d_names.append(config['SWUPSFC'])

    if config['clrskyradiation_output'] == True:
        variable2d_names.append(config['SWUPTOACLR'])
        variable2d_names.append(config['LWUPTOACLR'])
        variable2d_names.append(config['SWDOWNSFCCLR'])
        variable2d_names.append(config['LWDOWNSFCLR'])

    if config['convectiveparam'] == True:
        variable2d_names.append(config['PRECIPSFCTOT'])
        variable2d_names.append(config['PRECIPSFCSTRAT'])
        variable2d_names.append(config['PRECIPSFCCONV'])
    else:
        variable2d_names.append(config['PRECIPSFCLIQ'])
        variable2d_names.append(config['PRECIPSFCICE'])
    
    if config['colnc_output'] == True:
        variable2d_names.append(config['NC2D'])

    if config['aodabs_output'] == True:
        variable2d_names.append(config['AODABS'])

    if config['cosp_output'] == True:
        variablecosp_names = []
        variablecosp_names.append(config['IWPMODIS'])
        variablecosp_names.append(config['LWPMODIS'])
        variablecosp_names.append(config['REFFLIQMODIS'])
        variablecosp_names.append(config['REFFICEMODIS'])
        variablecosp_names.append(config['TAUICEMODIS'])
        variablecosp_names.append(config['TAUTOTMODIS'])
        variablecosp_names.append(config['TAULIQMODIS'])
        variablecosp_names.append(config['NDMODIS'])
        variablecosp_names.append(config['CFLIQMODIS'])
        variablecosp_names.append(config['CFICEMODIS'])
        variablecosp_names.append(config['CFTOTMODIS'])

    if config['cloudfraction_layers_output'] == True:
        variable2d_names.append(config['CLDHGH'])
        variable2d_names.append(config['CLDMED'])
        variable2d_names.append(config['CLDLOW'])

    if config['surfaceflux_output'] == True:
        variable2d_names.append(config['LHFLX'])
        variable2d_names.append(config['SHFLX'])

    for varname in variable2d_names:
        try:
            var = e3smdata2d[varname + E3SMdomain_range][:,x_idx_all].load()
            var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
        except:
            # var = xr.DataArray(np.zeros((len(e3smtime),len_ncol))*np.nan,name=varname,\
            #                    dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime,config['latlon_dim']+E3SMdomain_range:np.arange(len_ncol)},\
            #                    attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            var = xr.DataArray(np.zeros((len(e3smtime)))*np.nan,name=varname,\
                               dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
        if varname==config['AODABS'] + E3SMdomain_range or varname==config['AOD'] + E3SMdomain_range:
            var.attrs['units']='N/A'
        variable_names.append(varname)
        variables.append(var)

    for varname in variablecosp_names:
        try:
            var = e3smdatacosp[varname + E3SMdomain_range][:,x_idx_all]
            var[:] = e3smdatacosp[varname + E3SMdomain_range][:,x_idx_all].load()/e3smdatacosp[config['SUNLIT']+E3SMdomain_range][:,x_idx_all].values
            if varname == config['IWPMODIS'] + E3SMdomain_range or varname == config['REFFICEMODIS'] + E3SMdomain_range or varname == config['TAUICEMODIS'] + E3SMdomain_range:
                var[:] = var/e3smdatacosp[config['CFICEMODIS'] + E3SMdomain_range][:,x_idx_all].values
            if varname == config['TAUTOTMODIS'] + E3SMdomain_range:
                var[:] = var/e3smdatacosp[config['CFTOTMODIS'] + E3SMdomain_range][:,x_idx_all].values
            var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
        except:
            # var = xr.DataArray(np.zeros((len(e3smtime),len_ncol))*np.nan,name=varname,\
            #                    dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime,config['latlon_dim']+E3SMdomain_range:np.arange(len_ncol)},\
            #                    attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            var = xr.DataArray(np.zeros((len(e3smtime)))*np.nan,name=varname,\
                               dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
        variable_names.append(varname)
        variables.append(var)
    
    # all other 3D (with vertical level) variables at the lowest model level
    # variable3d_names = [config['Q'], config['T'], config['RH'], config['U'], config['V']] 
    variable3d_names = [config['T']]
    # if config['ccn_output'] == True:
    #     variable3d_names.append(config['CCN1'])
    #     variable3d_names.append(config['CCN2'])
    #     variable3d_names.append(config['CCN5'])
    #     variable3d_names.append(config['CCN10'])

    for varname in variable3d_names:
        try:
            var = e3smdata3d[varname + E3SMdomain_range][:,-1,x_idx_all].load()
            var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
        except:
            # var = xr.DataArray(np.zeros((len(e3smtime),len_lev,len_ncol))*np.nan,name=varname,\
            #                    dims=[config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime,config['vert_dim']:e3smdata3d[config['vert_dim']],config['latlon_dim']+E3SMdomain_range:e3smdata3d[config['latlon_dim']+E3SMdomain_range]},\
            #                    attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            var = xr.DataArray(np.zeros((len(e3smtime),len_lev))*np.nan,name=varname,\
                               dims=[config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime,config['vert_dim']:e3smdata3d[config['vert_dim']]},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
        variables.append(var)
        variable_names.append(varname)
    
    e3smdata3d.close()
    e3smdata2d.close()

    # aerosol composition
    print(ttt.ctime(ttt.time())+': aerosol composition:')
    # bc_a1 = e3smdata['bc_a1'+'_'+E3SMdomain_range].load()
    # bc_a3 = e3smdata['bc_a3'+'_'+E3SMdomain_range].load()
    # bc_a4 = e3smdata['bc_a4'+'_'+E3SMdomain_range].load()
    # dst_a1 = e3smdata['dst_a1'+'_'+E3SMdomain_range].load()
    # dst_a3 = e3smdata['dst_a3'+'_'+E3SMdomain_range].load()
    # mom_a1 = e3smdata['mom_a1'+'_'+E3SMdomain_range].load()
    # mom_a2 = e3smdata['mom_a2'+'_'+E3SMdomain_range].load()
    # mom_a3 = e3smdata['mom_a3'+'_'+E3SMdomain_range].load()
    # mom_a4 = e3smdata['mom_a4'+'_'+E3SMdomain_range].load()
    # ncl_a1 = e3smdata['ncl_a1'+'_'+E3SMdomain_range].load()
    # ncl_a2 = e3smdata['ncl_a2'+'_'+E3SMdomain_range].load()
    # ncl_a3 = e3smdata['ncl_a3'+'_'+E3SMdomain_range].load()
    # pom_a1 = e3smdata['pom_a1'+'_'+E3SMdomain_range].load()
    # pom_a3 = e3smdata['pom_a3'+'_'+E3SMdomain_range].load()
    # pom_a4 = e3smdata['pom_a4'+'_'+E3SMdomain_range].load()
    # so4_a1 = e3smdata['so4_a1'+'_'+E3SMdomain_range].load()
    # so4_a2 = e3smdata['so4_a2'+'_'+E3SMdomain_range].load()
    # so4_a3 = e3smdata['so4_a3'+'_'+E3SMdomain_range].load()
    # soa_a1 = e3smdata['soa_a1'+'_'+E3SMdomain_range].load()
    # soa_a2 = e3smdata['soa_a2'+'_'+E3SMdomain_range].load()
    # soa_a3 = e3smdata['soa_a3'+'_'+E3SMdomain_range].load()
    # bc_all = xr.concat([bc_a1[i,-1,x_idx_all[i]] + bc_a3[i,-1,x_idx_all[i]] + bc_a4[i,-1,x_idx_all[i]]
    #              for i in range(len(e3smtime))],dim='time')
    # dst_all = xr.concat([dst_a1[i,-1,x_idx_all[i]] + dst_a3[i,-1,x_idx_all[i]]
    #              for i in range(len(e3smtime))],dim='time')
    # mom_all = xr.concat([mom_a1[i,-1,x_idx_all[i]] + mom_a2[i,-1,x_idx_all[i]] + mom_a3[i,-1,x_idx_all[i]] + mom_a4[i,-1,x_idx_all[i]]
    #              for i in range(len(e3smtime))],dim='time')
    # ncl_all = xr.concat([ncl_a1[i,-1,x_idx_all[i]] + ncl_a2[i,-1,x_idx_all[i]] + ncl_a3[i,-1,x_idx_all[i]]
    #              for i in range(len(e3smtime))],dim='time')
    # pom_all = xr.concat([pom_a1[i,-1,x_idx_all[i]] + pom_a3[i,-1,x_idx_all[i]] + pom_a4[i,-1,x_idx_all[i]]
    #              for i in range(len(e3smtime))],dim='time')
    # so4_all = xr.concat([so4_a1[i,-1,x_idx_all[i]] + so4_a2[i,-1,x_idx_all[i]] + so4_a3[i,-1,x_idx_all[i]] 
    #              for i in range(len(e3smtime))],dim='time')
    # soa_all = xr.concat([soa_a1[i,-1,x_idx_all[i]] + soa_a2[i,-1,x_idx_all[i]] + soa_a3[i,-1,x_idx_all[i]] 
    #              for i in range(len(e3smtime))],dim='time')
    # bc_all.attrs['units'] = bc_a1.units
    # bc_all.attrs['long_name'] = 'black carbon concentration'
    # dst_all.attrs['units'] = dst_a1.units
    # dst_all.attrs['long_name'] = 'dust concentration'
    # mom_all.attrs['units'] = mom_a1.units
    # mom_all.attrs['long_name'] = 'marine organic matter concentration'
    # ncl_all.attrs['units'] = ncl_a1.units
    # ncl_all.attrs['long_name'] = 'sea salt concentration'
    # pom_all.attrs['units'] = pom_a1.units
    # pom_all.attrs['long_name'] = 'primary organic matter concentration'
    # so4_all.attrs['units'] = so4_a1.units
    # so4_all.attrs['long_name'] = 'sulfate concentration'
    # soa_all.attrs['units'] = soa_a1.units
    # soa_all.attrs['long_name'] = 'secondary organic aerosol concentration'
    # # change time to standard datetime64 format
    # bc_all.coords['time'] = bc_all.indexes['time'].to_datetimeindex()
    # dst_all.coords['time'] = dst_all.indexes['time'].to_datetimeindex()
    # mom_all.coords['time'] = mom_all.indexes['time'].to_datetimeindex()
    # ncl_all.coords['time'] = ncl_all.indexes['time'].to_datetimeindex()
    # pom_all.coords['time'] = pom_all.indexes['time'].to_datetimeindex()
    # so4_all.coords['time'] = so4_all.indexes['time'].to_datetimeindex()
    # soa_all.coords['time'] = soa_all.indexes['time'].to_datetimeindex()
    
    # # aerosol size
    # print(ttt.ctime(ttt.time())+': aerosol size:')
    # num_a1 = e3smdata['num_a1'+'_'+E3SMdomain_range].load()
    # num_a2 = e3smdata['num_a2'+'_'+E3SMdomain_range].load()
    # num_a3 = e3smdata['num_a3'+'_'+E3SMdomain_range].load()
    # num_a4 = e3smdata['num_a4'+'_'+E3SMdomain_range].load()
    # dn1 = e3smdata['dgnd_a01'+'_'+E3SMdomain_range].load()
    # dn2 = e3smdata['dgnd_a02'+'_'+E3SMdomain_range].load()
    # dn3 = e3smdata['dgnd_a03'+'_'+E3SMdomain_range].load()
    # dn4 = e3smdata['dgnd_a04'+'_'+E3SMdomain_range].load()
    # P0 = e3smdata['P0'].load()
    # hyam = e3smdata['hyam'].load()
    # hybm = e3smdata['hybm'].load()
    # if len(hyam.shape)>1:   # when use xr.open_mfdataset, hyam will add a time dimension
    #     hyam = hyam[0,:]
    #     hybm = hybm[0,:]
    #     P0 = P0[0]
    # T = e3smdata['T'+'_'+E3SMdomain_range].load()
    # PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
    # Pres = np.nan*T
    # zlen = T.shape[1]
    # for kk in range(zlen):
    #     Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
    # T = xr.concat([T[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # Pres = xr.concat([Pres[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data
    # numall = [xr.concat([num_a1[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #           xr.concat([num_a2[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #           xr.concat([num_a3[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #           xr.concat([num_a4[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #           ]
    # dnall = [xr.concat([dn1[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #          xr.concat([dn2[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #          xr.concat([dn3[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #          xr.concat([dn4[i,-1,x_idx_all[i]] for i in range(len(e3smtime))],dim='time').data,
    #          ]
    # NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:,-1], Pres[:,-1])
    # NCNall = xr.DataArray(data=NCN,  dims=["size", "time"],
    #     coords=dict(time=(["time"], e3smtime),size=(["size"], np.arange(1,3001))),
    #     attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)

    # e3smdata.close()

    if config['aerosol_output'] == True:
      print(lst3d_dryaer[0])
      e3smdata3d_dryaer = xr.open_dataset(lst3d_dryaer[0])
      e3smdata3d_dryaer = e3smdata3d_dryaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location

      vlist = list(e3smdata3d_dryaer.variables.keys())
      av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)

      if config['composition_output'] == True:
          req_vlist = [config['bc_a1'], config['bc_a3'], config['bc_a4'], config['dst_a1'], config['dst_a3'], config['mom_a1'], \
                       config['mom_a2'], config['mom_a3'], config['mom_a4'], config['ncl_a1'], config['ncl_a2'], config['ncl_a3'], \
                       config['pom_a1'], config['pom_a3'], config['pom_a4'], config['so4_a1'], config['so4_a2'], config['so4_a3'], \
                       config['soa_a1'], config['soa_a2'], config['soa_a3']]
          req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
          matched_vlist = list(set(av_vars).intersection(req_vlist))
      
          if len(matched_vlist) == len(req_vlist):
              print('\nAnalyzing aerosol composition')
              bc_a1 = e3smdata3d_dryaer[config['bc_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              bc_a3 = e3smdata3d_dryaer[config['bc_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              bc_a4 = e3smdata3d_dryaer[config['bc_a4']+E3SMdomain_range][:,-1,x_idx_all].load()
              dst_a1 = e3smdata3d_dryaer[config['dst_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              dst_a3 = e3smdata3d_dryaer[config['dst_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              mom_a1 = e3smdata3d_dryaer[config['mom_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              mom_a2 = e3smdata3d_dryaer[config['mom_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
              mom_a3 = e3smdata3d_dryaer[config['mom_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              mom_a4 = e3smdata3d_dryaer[config['mom_a4']+E3SMdomain_range][:,-1,x_idx_all].load()
              ncl_a1 = e3smdata3d_dryaer[config['ncl_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              ncl_a2 = e3smdata3d_dryaer[config['ncl_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
              ncl_a3 = e3smdata3d_dryaer[config['ncl_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              pom_a1 = e3smdata3d_dryaer[config['pom_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              pom_a3 = e3smdata3d_dryaer[config['pom_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              pom_a4 = e3smdata3d_dryaer[config['pom_a4']+E3SMdomain_range][:,-1,x_idx_all].load()
              so4_a1 = e3smdata3d_dryaer[config['so4_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              so4_a2 = e3smdata3d_dryaer[config['so4_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
              so4_a3 = e3smdata3d_dryaer[config['so4_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              soa_a1 = e3smdata3d_dryaer[config['soa_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
              soa_a2 = e3smdata3d_dryaer[config['soa_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
              soa_a3 = e3smdata3d_dryaer[config['soa_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
              bc_all  = bc_a1 +  bc_a3 + bc_a4
              dst_all = dst_a1 +  dst_a3
              mom_all = mom_a1 + mom_a2 + mom_a3 + mom_a4
              ncl_all = ncl_a1 + ncl_a2 + ncl_a3
              pom_all = pom_a1 + pom_a3 + pom_a4
              so4_all = so4_a1 + so4_a2 + so4_a3
              soa_all = soa_a1 + soa_a2 + soa_a3
              bc_all.attrs['units'] = bc_a1.units
              bc_all.attrs['long_name'] = 'black carbon concentration'
              dst_all.attrs['units'] = dst_a1.units
              dst_all.attrs['long_name'] = 'dust concentration'
              mom_all.attrs['units'] = mom_a1.units
              mom_all.attrs['long_name'] = 'marine organic matter concentration'
              ncl_all.attrs['units'] = ncl_a1.units
              ncl_all.attrs['long_name'] = 'sea salt concentration'
              pom_all.attrs['units'] = pom_a1.units
              pom_all.attrs['long_name'] = 'primary organic matter concentration'
              so4_all.attrs['units'] = so4_a1.units
              so4_all.attrs['long_name'] = 'sulfate concentration'
              soa_all.attrs['units'] = soa_a1.units
              soa_all.attrs['long_name'] = 'secondary organic aerosol concentration'
              # change time to standard datetime64 format
              bc_all.coords[config['time_dim']] = bc_all.indexes[config['time_dim']].to_datetimeindex()
              dst_all.coords[config['time_dim']] = dst_all.indexes[config['time_dim']].to_datetimeindex()
              mom_all.coords[config['time_dim']] = mom_all.indexes[config['time_dim']].to_datetimeindex()
              ncl_all.coords[config['time_dim']] = ncl_all.indexes[config['time_dim']].to_datetimeindex()
              pom_all.coords[config['time_dim']] = pom_all.indexes[config['time_dim']].to_datetimeindex()
              so4_all.coords[config['time_dim']] = so4_all.indexes[config['time_dim']].to_datetimeindex()
              soa_all.coords[config['time_dim']] = soa_all.indexes[config['time_dim']].to_datetimeindex()
          else:
              bc_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              dst_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              mom_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              ncl_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              pom_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              so4_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              soa_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
      # aerosol size
      if config['dgnum_output_combined'] == False:
          req_vlist = [config['num_a1'], config['num_a2'], config['num_a3'], config['num_a4'], \
                       config['dgnd_a01'], config['dgnd_a02'], config['dgnd_a03'], config['dgnd_a04']]
          req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
          matched_vlist = list(set(av_vars).intersection(req_vlist))
          if len(matched_vlist) == len(req_vlist):
              print('\nAnalyzing aerosol size')
              num_a1 = e3smdata3d_dryaer[config['num_a1']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_a2 = e3smdata3d_dryaer[config['num_a2']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_a3 = e3smdata3d_dryaer[config['num_a3']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_a4 = e3smdata3d_dryaer[config['num_a4']+E3SMdomain_range][:, -1, x_idx_all].load()
              numall = [num_a1.data, num_a2.data, num_a3.data, num_a4.data]
              num_dn1 = e3smdata3d_dryaer[config['dgnd_a01']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_dn2 = e3smdata3d_dryaer[config['dgnd_a02']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_dn3 = e3smdata3d_dryaer[config['dgnd_a03']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_dn4 = e3smdata3d_dryaer[config['dgnd_a04']+E3SMdomain_range][:, -1, x_idx_all].load()
              dnall  = [num_dn1.data,    num_dn2.data,    num_dn3.data,    num_dn4.data]
              NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1].data, Pres[:, -1].data)
              NCNall = xr.DataArray(data=NCN,  dims=["size", "time"],
                  coords=dict(time=([config['time_dim']], e3smtime),size=(["size"], np.arange(1,3001))),
                  attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
          else:
              NCNall = xr.DataArray(np.zeros((3000,len(e3smtime)))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
      else:
          req_vlist = [config['num_a1'], config['num_a2'], config['num_a3'], config['num_a4'], \
                       config['dgnum']]
          req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
          matched_vlist = list(set(av_vars).intersection(req_vlist))
          if len(matched_vlist) == len(req_vlist):
              print('\nAnalyzing for aerosol size')
              num_a1 = e3smdata3d_dryaer[config['num_a1']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_a2 = e3smdata3d_dryaer[config['num_a2']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_a3 = e3smdata3d_dryaer[config['num_a3']+E3SMdomain_range][:, -1, x_idx_all].load()
              num_a4 = e3smdata3d_dryaer[config['num_a4']+E3SMdomain_range][:, -1, x_idx_all].load()
              numall = [num_a1.data, num_a2.data, num_a3.data, num_a4.data]
              num_dn1 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 0].load()
              num_dn2 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 1].load()
              num_dn3 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 2].load()
              num_dn4 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 3].load()
              dnall  = [num_dn1.data,    num_dn2.data,    num_dn3.data,    num_dn4.data]
              NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1].data, Pres[:, -1].data)
              NCNall = xr.DataArray(data=NCN,  dims=["size", "time"],
                  coords=dict(time=([config['time_dim']], e3smtime),size=(["size"], np.arange(1,3001))),
                  attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
          else:
              NCNall = xr.DataArray(np.zeros((3000,len(e3smtime)))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            
    #CCN
    if config['ccn_output'] == True:
        print(lst3d_wetaer[0])
        e3smdata3d_wetaer = xr.open_dataset(lst3d_wetaer[0])
        e3smdata3d_wetaer = e3smdata3d_wetaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location

        vlist = list(e3smdata3d_wetaer.variables.keys())
        av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)

        req_vlist = [config['CCN1'], config['CCN2']]#, config['CCN5'], config['CCN10']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
      
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing CCN')
            ccn1_all = e3smdata3d_wetaer[config['CCN1']+E3SMdomain_range][:,-1,x_idx_all].load()
            ccn2_all = e3smdata3d_wetaer[config['CCN2']+E3SMdomain_range][:,-1,x_idx_all].load()
            # ccn5_all = e3smdata3d_wetaer[config['CCN5']+E3SMdomain_range][:,-1,x_idx_all].load()
            # ccn10_all = e3smdata3d_wetaer[config['CCN10']+E3SMdomain_range][:,-1,x_idx_all].load()
        else:
            ccn1_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            ccn2_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            # ccn5_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            # ccn10_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

    #%%  add data for each day
    for ii in np.arange(len(lst3d[1:]))+1:
        print(lst3d[ii])
        # e3smdata = xr.open_dataset(file)
        e3smdata3d = xr.open_dataset(lst3d[ii])
        e3smdata3d = e3smdata3d.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
        dimensions3d = (config['time_dim'], config['vert_dim'], config['latlon_dim']+E3SMdomain_range)
        variables3d = []
        for var_name, var_data in e3smdata3d.data_vars.items():
            if set(var_data.dims) == set(dimensions3d):
                variables3d.append(var_name)
      
        e3smdata2d = xr.open_dataset(lst2d[ii])
        e3smdata2d = e3smdata2d.transpose(config['time_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time and location
        dimensions2d = (config['time_dim'], config['latlon_dim']+E3SMdomain_range)
        variables2d = []
        for var_name, var_data in e3smdata2d.data_vars.items():
            if set(var_data.dims) == set(dimensions2d):
                variables2d.append(var_name)

        e3smdatacosp = xr.open_dataset(lstcosp[ii])
        e3smdatacosp = e3smdatacosp.transpose(config['time_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time and location
        dimensions2d = (config['time_dim'], config['latlon_dim']+E3SMdomain_range)
        variablescosp = []
        for var_name, var_data in e3smdatacosp.data_vars.items():
            if set(var_data.dims) == set(dimensions2d):
                variablescosp.append(var_name)
          
        e3smtime_i = e3smdata3d.indexes[config['time_dim']].to_datetimeindex()
        e3smtime = np.hstack((e3smtime, e3smtime_i))

        lon_0 = median_time_1d(time, lon, e3smtime_i)
        lat_0 = median_time_1d(time, lat, e3smtime_i)
        x_idx_all = []
        for tt in range(len(e3smtime_i)):
            if np.isnan(lon_0[tt]):
                x_idx_all.append(-1)
            else:
                x_idx_tmp = find_nearest(lonm, latm, lon_0[tt], lat_0[tt])
                x_idx_all.append(x_idx_tmp)   

        T = e3smdata3d[config['T']+E3SMdomain_range][:,:,x_idx_all].load()

        if config['pres_output'] == False:
            PS = e3smdata3d[config['PS']+E3SMdomain_range][:,x_idx_all].load()
            Pres = np.nan*T
            zlen = T.shape[1]
            for kk in range(zlen):
              Pres[:, kk] = hyam[kk]*P0  +  hybm[kk]*PS
        else:
            Pres = e3smdata3d[config['PRES']+E3SMdomain_range][:,:,x_idx_all].load()

        # vlist3d = list(e3smdata3d.variables.keys())
        vlist3d = variables3d
        av_vars3d = fnmatch.filter(vlist3d,'*'+E3SMdomain_range)
    
        # vlist2d = list(e3smdata2d.variables.keys())
        vlist2d = variables2d
        av_vars2d = fnmatch.filter(vlist2d,'*'+E3SMdomain_range)

        vlistcosp = variablescosp
        av_varscosp = fnmatch.filter(vlistcosp,'*'+E3SMdomain_range)

        # variables to calculate cloud heights and depth
        print('\nAnalyzing cloud base, top, and depth')
        req_vlist = [config['Z'], config['CF']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars3d).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing cloud base, top, and depth')
            z3 = e3smdata3d[config['Z']+E3SMdomain_range][:,:,x_idx_all].load()
            if config['cloudfraction3d_legit'] == True:
                cloud = e3smdata3d[config['CF']+E3SMdomain_range][:,:,x_idx_all].load()
            else:
                cloud = e3smdata3d[config['CFLIQ']+E3SMdomain_range][:,:,x_idx_all].load() + e3smdata3d[config['CFICE']+E3SMdomain_range][:,:,x_idx_all].load()
                cloud = cloud.where(cloud <= 1, 1) # make sure CF doesn't exceed 100%
            dz = (z3[:,:-2].data - z3[:,2:].data)/2
            dz = np.append(dz, (z3[:,-2:-1].data+z3[:,-1:].data)/2, axis=1)
            dz = np.insert(dz,0,dz[:,0],axis=1)
            
            # mid-level T, P, z
            Tmid = 0.5*(T[:,0:-1].data + T[:,1:].data)
            Pmid = 0.5*(Pres[:,0:-1].data + Pres[:,1:].data)
            Zmid = 0.5*(z3[:,0:-1].data + z3[:,1:].data)
            # cloud top
            cf_sum_top = np.cumsum(cloud.data, axis=1)
            cf_sum_top[cf_sum_top > 1] = 1
            cf_sum_top_diff = cf_sum_top[:,1:] - cf_sum_top[:,0:-1]
            z_cldtop = np.sum(Zmid*cf_sum_top_diff, axis=1)
            T_cldtop = np.sum(Tmid*cf_sum_top_diff, axis=1)
            # cloud base
            cf_sum_base = np.cumsum(cloud[:, ::-1].data, axis=1)
            cf_sum_base[cf_sum_base > 1] = 1
            cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
            z_cldbase = np.sum(Zmid[:,::-1]*cf_sum_base_diff, axis=1)
            T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff, axis=1)
            # normalize cloud fraction when not 100%
            z_cldtop = z_cldtop / cf_sum_top[:,-1]
            T_cldtop = T_cldtop / cf_sum_top[:,-1]
            z_cldbase = z_cldbase / cf_sum_base[:,-1]
            T_cldbase = T_cldbase / cf_sum_base[:,-1]
            e3sm_cloud_depth = z_cldtop - z_cldbase
            cloud_depth_2 = xr.DataArray(data=e3sm_cloud_depth,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime_i)),
                attrs=dict(long_name="cloud depth",units="m"),)
            cbh_2 = xr.DataArray(data=z_cldbase,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime_i)),
                attrs=dict(long_name="cloud base height",units="m"),)
            cth_2 = xr.DataArray(data=z_cldtop,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime_i)),
                attrs=dict(long_name="cloud top height",units="m"),)
            cbt_2 = xr.DataArray(data=T_cldbase,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime_i)),
                attrs=dict(long_name="cloud base temperature",units="K"),)
            ctt_2 = xr.DataArray(data=T_cldtop,  dims=[config['time_dim']],
                coords=dict(time=([config['time_dim']], e3smtime_i)),
                attrs=dict(long_name="cloud top temperature",units="K"),)
            cbh = xr.concat([cbh, cbh_2], dim=config['time_dim'])
            cth = xr.concat([cth, cth_2], dim=config['time_dim'])
            cbt = xr.concat([cbt, cbt_2], dim=config['time_dim'])
            ctt = xr.concat([ctt, ctt_2], dim=config['time_dim'])
            cloud_depth = xr.concat([cloud_depth, cloud_depth_2], dim=config['time_dim'])
        else:
            cloud_depth = xr.concat([cloud_depth, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cloud_depth',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            cbh = xr.concat([cbh, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cbh',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            cth = xr.concat([cth, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cth',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            cbt = xr.concat([cbt, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cbt',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            ctt = xr.concat([ctt, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='ctt',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
        
    
        # cloud optical depth and effective radius
        if config['reff_output'] == True:
            req_vlist = [config['REL'], config['CFLIQ'], config['NC']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars3d).intersection(req_vlist))
      
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing effective radius')
                rel = e3smdata3d[config['REL']+E3SMdomain_range][:,:,x_idx_all].load()
                freql = e3smdata3d[config['CFLIQ']+E3SMdomain_range][:,:,x_idx_all].load()
                icwnc = e3smdata3d[config['NC']+E3SMdomain_range][:,:,x_idx_all].load()
        
                # calculate mean effective radius. 
                reff = calc_Reff_from_REL(rel.data, dz, freql.data, icwnc.data)
                reff[reff==0] = np.nan
          
                reff_2 = xr.DataArray(data=reff,  dims=[config['time_dim']],
                    coords=dict(time=([config['time_dim']], e3smtime_i)),
                    attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
                reff_mean = xr.concat([reff_mean, reff_2], dim=config['time_dim'])
            else:
                reff_mean = xr.concat([reff_mean, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='reff',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
        
        if config['tau3d_output'] == True:
            req_vlist3d = [config['TAU3D']]
            req_vlist2d = [config['SWDOWNTOA']]
            req_vlist3d = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist3d]
            req_vlist2d = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist2d]
            matched_vlist3d = list(set(av_vars3d).intersection(req_vlist3d))
            matched_vlist2d = list(set(av_vars2d).intersection(req_vlist2d))
        
            if np.logical_and(len(matched_vlist3d) == len(req_vlist3d), len(matched_vlist2d) == len(req_vlist2d)):
                print('\nAnalyzing cloud optical depth')
                cod_a = e3smdata3d[config['TAU3D']+E3SMdomain_range][:,:,x_idx_all].load()
                solin = e3smdata2d[config['SWDOWNTOA']+E3SMdomain_range][:,x_idx_all].load()
              
                # calculate mean optical depth
                cod = np.sum(cod_a.data,axis=1)
                # cod[solin==0] = np.nan        
                
                cod_2 = xr.DataArray(data=cod,  dims=[config['time_dim']],
                  coords=dict(time=([config['time_dim']], e3smtime_i)),
                  attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
                cod_mean = xr.concat([cod_mean, cod_2], dim=config['time_dim'])
            else:
                cod_mean = xr.concat([cod_mean, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cod',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])

        if config['cosp_output'] == True:   
            req_vlist = [config['TAULIQMODIS']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_varscosp).intersection(req_vlist))
            
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing MODIS simulator cloud optical depth')
                cod_m2 = e3smdatacosp[config['TAULIQMODIS']+E3SMdomain_range][:,x_idx_all]
                cod_m2[:] = e3smdatacosp[config['TAULIQMODIS']+E3SMdomain_range][:,x_idx_all].load()/e3smdatacosp[config['SUNLIT']+E3SMdomain_range][:,x_idx_all].values
                cod_m = xr.concat([cod_m, cod_m2], dim=config['time_dim'])
            else:
                cod_m = xr.concat([cod_m, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])

        # mean cloud droplet number concentration
        if config['colnc_output'] == True:
            req_vlist = [config['NC2D']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars2d).intersection(req_vlist))
            
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing mean cloud droplet number concentration')
                cdnc_col = e3smdata3d[config['NC2D']+E3SMdomain_range][:,x_idx_all].load()
                weight = cloud*dz
                cdnc = cdnc_col/np.sum(weight,axis=1)
                cdnc[cdnc >2e9] = np.nan
                cdnc = xr.DataArray(data=cdnc,  dims=[config['time_dim']],
                    coords=dict(time=([config['time_dim']], e3smtime_i)),
                    attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
                cdnc_mean = xr.concat([cdnc_mean, cdnc], dim=config['time_dim'])
            else:
                cdnc_mean = xr.concat([cdnc_mean, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cdnc_mean',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
        else:
            req_vlist = [config['NC']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars3d).intersection(req_vlist))
          
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing mean cloud droplet number concentration')
                #compute cloud layer mean CDNC from 3D NC (note that if NC is not in-cloud only, one needs to divide by cloud fraction)
                nc3d = e3smdata3d[config['NC']+E3SMdomain_range][:,:,x_idx_all].load()   
                if nc3d.attrs['units'] == '1/kg':
                  rho = np.array(Pres/T/287.06)
                  cdnc_rel = nc3d*rho/cloud
                if nc3d.attrs['units'] == 'm-3':
                  cdnc_rel = nc3d/cloud
                cdnc_rel = cdnc_rel.where(cloud > 0, other = 0)
                weight = cloud*dz
                weight_column = weight.sum(dim=config['vert_dim'])
                cdnc_rel_avg = cdnc_rel.dot(weight, dims=config['vert_dim'])
                cdnc = np.divide(cdnc_rel_avg, weight_column)
                cdnc = xr.DataArray(data=cdnc,  dims=[config['time_dim']],
                        coords=dict(time=([config['time_dim']], e3smtime_i)),
                        attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
                cdnc_mean = xr.concat([cdnc_mean, cdnc], dim=config['time_dim'])
            else:
                cdnc_mean = r.concat([cdnc_mean, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            
        # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
        if config['tau3d_output'] == True:
            req_vlist = [config['LWP']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars2d).intersection(req_vlist))
            
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing for cloud droplet number concentration retrieved like Ndrop')
                lwp = e3smdata2d[config['LWP']+E3SMdomain_range][:,x_idx_all].data
                e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
                nd_arm = calc_cdnc_ARM(lwp, cod_2, e3sm_cloud_depth)
                nd_arm = xr.DataArray(data=nd_arm*1e6,  dims=[config['time_dim']],
                    coords=dict(time=([config['time_dim']], e3smtime_i)),
                    attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                                description='Retrieved using ARM Ndrop algorithm'),)
                cdnc_arm = xr.concat([cdnc_arm, nd_arm], dim=config['time_dim'])
            else:
                cdnc_arm = xr.concat([cdnc_arm, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cdnc_arm',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
          
        elif config['tau3d_output'] == False and config['cosp_output'] == True:
            req_vlist = [config['LWP']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars2d).intersection(req_vlist))
            
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing for cloud droplet number concentration retrieved like Ndrop')
                lwp = e3smdata2d[config['LWP']+E3SMdomain_range][:,x_idx_all].data
                e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
                nd_arm = calc_cdnc_ARM(lwp, cod_m2, e3sm_cloud_depth)
                nd_arm = xr.DataArray(data=nd_arm*1e6,  dims=[config['time_dim']],
                    coords=dict(time=([config['time_dim']], e3smtime_i)),
                    attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                                description='Retrieved using ARM Ndrop algorithm'),)
                cdnc_arm = xr.concat([cdnc_arm, nd_arm], dim=config['time_dim'])
            else:
                cdnc_arm = xr.concat([cdnc_arm, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cdnc_arm',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
    
        if config['cosp_output'] == True:
            # req_vlist = [config['LWP']]
            req_vlist = [config['LWPMODIS']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_varscosp).intersection(req_vlist))
            
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing cloud droplet number concentration retrieved like MODIS')
                # lwp = e3smdata2d[config['LWP']+E3SMdomain_range][:,x_idx_all].data
                lwp = e3smdatacosp[config['LWPMODIS']+E3SMdomain_range][:,x_idx_all]
                lwp[:] = e3smdatacosp[config['LWPMODIS']+E3SMdomain_range][:,x_idx_all].load()/e3smdatacosp[config['SUNLIT']+E3SMdomain_range][:,x_idx_all].values
                T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
                nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m2, adiabaticity=0.8)
                nd_sat = xr.DataArray(data=nd_sat*1e6,  dims=[config['time_dim']],
                    coords=dict(time=([config['time_dim']], e3smtime_i)),
                    attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                                description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)
                cdnc_sat = xr.concat([cdnc_sat, nd_sat], dim=config['time_dim'])
            else:
                cdnc_sat = xr.concat([cdnc_sat, xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cdnc_sat',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
        
        # all other 2D (surface and vertical integrated) variables
        for varname in variable2d_names:
            try:
                var = e3smdata2d[varname + E3SMdomain_range][:,x_idx_all].load()
                var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
            except:
                # var = xr.DataArray(np.zeros((len(e3smtime_i),len_ncol))*np.nan,name=varname,\
                #                dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime_i,config['latlon_dim']+E3SMdomain_range:np.arange(len_ncol)},\
                #                attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
                var = xr.DataArray(np.zeros((len(e3smtime_i)))*np.nan,name=varname,\
                               dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime_i},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            vv = variable_names.index(varname)
            variables[vv] = xr.concat([variables[vv], var],dim=config['time_dim'])

        for varname in variablecosp_names:
          try:
              var = e3smdatacosp[varname + E3SMdomain_range][:,x_idx_all]
              var[:] = e3smdatacosp[varname + E3SMdomain_range][:,x_idx_all].load()/e3smdatacosp[config['SUNLIT']+E3SMdomain_range][:,x_idx_all].values
              if varname == config['IWPMODIS'] + E3SMdomain_range or varname == config['REFFICEMODIS'] + E3SMdomain_range or varname == config['TAUICEMODIS'] + E3SMdomain_range:
                  var[:] = var/e3smdatacosp[config['CFICEMODIS'] + E3SMdomain_range][:,x_idx_all].values
              if varname == config['TAUTOTMODIS'] + E3SMdomain_range:
                  var[:] = var/e3smdatacosp[config['CFTOTMODIS'] + E3SMdomain_range][:,x_idx_all].values
              var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
          except:
              # var = xr.DataArray(np.zeros((len(e3smtime),len_ncol))*np.nan,name=varname,\
              #                    dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime,config['latlon_dim']+E3SMdomain_range:np.arange(len_ncol)},\
              #                    attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
              var = xr.DataArray(np.zeros((len(e3smtime)))*np.nan,name=varname,\
                                 dims=[config['time_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime},\
                                 attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
          vv = variable_names.index(varname)
          variables[vv] = xr.concat([variables[vv], var],dim=config['time_dim'])
        
        # all other 3D (with vertical level) variables at the lowest model level
        for varname in variable3d_names:
            try:
                var = e3smdata3d[varname + E3SMdomain_range][:,-1,x_idx_all].load()
                var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
            except:
                # var = xr.DataArray(np.zeros((len(e3smtime_i),len_lev,len_ncol))*np.nan,name=varname,\
                #                dims=[config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime_i,config['vert_dim']:e3smdata3d[config['vert_dim']],config['latlon_dim']+E3SMdomain_range:e3smdata3d[config['latlon_dim']+E3SMdomain_range]},\
                #                attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
                var = xr.DataArray(np.zeros((len(e3smtime_i),len_lev))*np.nan,name=varname,\
                               dims=[config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime_i,config['vert_dim']:e3smdata3d[config['vert_dim']]},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            vv = variable_names.index(varname)
            variables[vv] = xr.concat([variables[vv], var],dim=config['time_dim'])
    
        e3smdata3d.close()
        e3smdata2d.close()
      
        if config['aerosol_output'] == True:
            print(lst3d_dryaer[ii])
            e3smdata3d_dryaer = xr.open_dataset(lst3d_dryaer[ii])
            e3smdata3d_dryaer = e3smdata3d_dryaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location

            vlist = list(e3smdata3d_dryaer.variables.keys())
            av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)

            if config['composition_output'] == True:
                # aerosol composition
                req_vlist = [config['bc_a1'], config['bc_a3'], config['bc_a4'], config['dst_a1'], config['dst_a3'], config['mom_a1'], \
                           config['mom_a2'], config['mom_a3'], config['mom_a4'], config['ncl_a1'], config['ncl_a2'], config['ncl_a3'], \
                           config['pom_a1'], config['pom_a3'], config['pom_a4'], config['so4_a1'], config['so4_a2'], config['so4_a3'], \
                           config['soa_a1'], config['soa_a2'], config['soa_a3']]
                req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                matched_vlist = list(set(av_vars).intersection(req_vlist))
              
                if len(matched_vlist) == len(req_vlist):
                    print('\nAnalyzing aerosol composition')
                    bc_a1 = e3smdata3d_dryaer[config['bc_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    bc_a3 = e3smdata3d_dryaer[config['bc_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    bc_a4 = e3smdata3d_dryaer[config['bc_a4']+E3SMdomain_range][:,-1,x_idx_all].load()
                    dst_a1 = e3smdata3d_dryaer[config['dst_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    dst_a3 = e3smdata3d_dryaer[config['dst_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    mom_a1 = e3smdata3d_dryaer[config['mom_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    mom_a2 = e3smdata3d_dryaer[config['mom_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
                    mom_a3 = e3smdata3d_dryaer[config['mom_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    mom_a4 = e3smdata3d_dryaer[config['mom_a4']+E3SMdomain_range][:,-1,x_idx_all].load()
                    ncl_a1 = e3smdata3d_dryaer[config['ncl_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    ncl_a2 = e3smdata3d_dryaer[config['ncl_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
                    ncl_a3 = e3smdata3d_dryaer[config['ncl_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    pom_a1 = e3smdata3d_dryaer[config['pom_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    pom_a3 = e3smdata3d_dryaer[config['pom_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    pom_a4 = e3smdata3d_dryaer[config['pom_a4']+E3SMdomain_range][:,-1,x_idx_all].load()
                    so4_a1 = e3smdata3d_dryaer[config['so4_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    so4_a2 = e3smdata3d_dryaer[config['so4_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
                    so4_a3 = e3smdata3d_dryaer[config['so4_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    soa_a1 = e3smdata3d_dryaer[config['soa_a1']+E3SMdomain_range][:,-1,x_idx_all].load()
                    soa_a2 = e3smdata3d_dryaer[config['soa_a2']+E3SMdomain_range][:,-1,x_idx_all].load()
                    soa_a3 = e3smdata3d_dryaer[config['soa_a3']+E3SMdomain_range][:,-1,x_idx_all].load()
                    bc  = bc_a1 + bc_a3 + bc_a4
                    dst = dst_a1 + dst_a3
                    mom = mom_a1 + mom_a2 + mom_a3 + mom_a4
                    ncl = ncl_a1 + ncl_a2 + ncl_a3
                    pom = pom_a1 + pom_a3 + pom_a4
                    so4 = so4_a1 + so4_a2 + so4_a3
                    soa = soa_a1 + soa_a2 + soa_a3
                    # change time to standard datetime64 format
                    bc.coords[config['time_dim']] = bc.indexes[config['time_dim']].to_datetimeindex()
                    dst.coords[config['time_dim']] = dst.indexes[config['time_dim']].to_datetimeindex()
                    mom.coords[config['time_dim']] = mom.indexes[config['time_dim']].to_datetimeindex()
                    ncl.coords[config['time_dim']] = ncl.indexes[config['time_dim']].to_datetimeindex()
                    pom.coords[config['time_dim']] = pom.indexes[config['time_dim']].to_datetimeindex()
                    so4.coords[config['time_dim']] = so4.indexes[config['time_dim']].to_datetimeindex()
                    soa.coords[config['time_dim']] = soa.indexes[config['time_dim']].to_datetimeindex()
                    bc_all = xr.concat([bc_all, bc], dim=config['time_dim'])
                    dst_all = xr.concat([dst_all, dst], dim=config['time_dim'])
                    mom_all = xr.concat([mom_all, mom], dim=config['time_dim'])
                    ncl_all = xr.concat([ncl_all, ncl], dim=config['time_dim'])
                    pom_all = xr.concat([pom_all, pom], dim=config['time_dim'])
                    so4_all = xr.concat([so4_all, so4], dim=config['time_dim'])
                    soa_all = xr.concat([soa_all, soa], dim=config['time_dim'])
                else:
                    bc_all  = xr.concat([bc_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='bc_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                    dst_all = xr.concat([dst_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='dst_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                    mom_all = xr.concat([mom_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='mom_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                    ncl_all = xr.concat([ncl_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='ncl_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                    pom_all = xr.concat([pom_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='pom_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                    so4_all = xr.concat([so4_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='so4_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                    soa_all = xr.concat([soa_all, xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='soa_all',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            
            # aerosol size
            if config['dgnum_output_combined'] == False:
                req_vlist = [config['num_a1'], config['num_a2'], config['num_a3'], config['num_a4'], \
                             config['dgnd_a01'], config['dgnd_a02'], config['dgnd_a03'], config['dgnd_a04']]
                req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                matched_vlist = list(set(av_vars).intersection(req_vlist))
            
                if len(matched_vlist) == len(req_vlist):
                    print('\nAnalyzing aerosol size')
                    num_a1 = e3smdata3d_dryaer[config['num_a1']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_a2 = e3smdata3d_dryaer[config['num_a2']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_a3 = e3smdata3d_dryaer[config['num_a3']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_a4 = e3smdata3d_dryaer[config['num_a4']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_dn1 = e3smdata3d_dryaer[config['dgnd_a01']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_dn2 = e3smdata3d_dryaer[config['dgnd_a02']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_dn3 = e3smdata3d_dryaer[config['dgnd_a03']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_dn4 = e3smdata3d_dryaer[config['dgnd_a04']+E3SMdomain_range][:, -1, x_idx_all].load()
                    numall = [num_a1.data, num_a2.data, num_a3.data, num_a4.data]
                    dnall  = [num_dn1.data,    num_dn2.data,    num_dn3.data,    num_dn4.data]
                    NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1].data, Pres[:, -1].data)
                    NCN2 = xr.DataArray(data=NCN,  dims=["size", config['time_dim']],
                        coords=dict(time=([config['time_dim']], e3smtime_i),size=(["size"], np.arange(1,3001))),
                        attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
                    NCNall = xr.concat([NCNall, NCN2], dim=config['time_dim'])
                else:
                    NCNall = xr.concat([NCNall, xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name='NCNall',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
            else:
                req_vlist = [config['num_a1'], config['num_a2'], config['num_a3'], config['num_a4'], \
                             config['dgnum']]
                req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                matched_vlist = list(set(av_vars).intersection(req_vlist))
            
                if len(matched_vlist) == len(req_vlist):
                    print('\nAnalyzing aerosol size')
                    num_a1 = e3smdata3d_dryaer[config['num_a1']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_a2 = e3smdata3d_dryaer[config['num_a2']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_a3 = e3smdata3d_dryaer[config['num_a3']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_a4 = e3smdata3d_dryaer[config['num_a4']+E3SMdomain_range][:, -1, x_idx_all].load()
                    num_dn1 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 0].load()
                    num_dn2 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 1].load()
                    num_dn3 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 2].load()
                    num_dn4 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:, -1, x_idx_all, 3].load()
                    numall = [num_a1.data, num_a2.data, num_a3.data, num_a4.data]
                    dnall  = [num_dn1.data,    num_dn2.data,    num_dn3.data,    num_dn4.data]
                    NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1].data, Pres[:, -1].data)
                    NCN2 = xr.DataArray(data=NCN,  dims=["size", config['time_dim']],
                        coords=dict(time=([config['time_dim']], e3smtime_i),size=(["size"], np.arange(1,3001))),
                        attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
                    NCNall = xr.concat([NCNall, NCN2], dim=config['time_dim'])
                else:
                    NCNall = xr.concat([NCNall, xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name='NCNall',attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
      
        #CCN
        if config['ccn_output'] == True:
            print(lst3d_wetaer[ii])
            e3smdata3d_wetaer = xr.open_dataset(lst3d_wetaer[ii])
            e3smdata3d_wetaer = e3smdata3d_wetaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
    
            vlist = list(e3smdata3d_wetaer.variables.keys())
            av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
          
            req_vlist = [config['CCN1'], config['CCN2']]#, config['CCN5'], config['CCN10']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars).intersection(req_vlist))
          
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing CCN')
                ccn1 = e3smdata3d_wetaer[config['CCN1']+E3SMdomain_range][:,-1,x_idx_all].load()
                ccn2 = e3smdata3d_wetaer[config['CCN2']+E3SMdomain_range][:,-1,x_idx_all].load()
                # ccn5 = e3smdata3d_wetaer[config['CCN5']+E3SMdomain_range][:,-1,x_idx_all].load()
                # ccn10 = e3smdata3d_wetaer[config['CCN10']+E3SMdomain_range][:,-1,x_idx_all].load()
                ccn1_all = xr.concat([ccn1_all, ccn1], dim=config['time_dim'])
                ccn2_all = xr.concat([ccn2_all, ccn2], dim=config['time_dim'])
                # ccn5_all = xr.concat([ccn5_all, ccn5], dim=config['time_dim'])
                # ccn10_all = xr.concat([ccn10_all, ccn10], dim=config['time_dim'])
            else:
                ccn1_all = xr.concat([ccn1_all, xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name=config['CCN1'],attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                ccn2_all = xr.concat([ccn2_all, xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name=config['CCN2'],attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                # ccn5_all = xr.concat([ccn5_all, xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name=config['CCN5'],attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
                # ccn10_all = xr.concat([ccn10_all, xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name=config['CCN10'],attrs={'units':'dummy_unit','long_name':'Dummy'})], dim=config['time_dim'])
          
           
    # put all variables into the list
    # # aerosol composition    
    # variable_names = variable_names + ['bc','dst','mom','ncl','pom','so4','soa']
    # variables = variables + [bc_all,dst_all,mom_all,ncl_all,pom_all,so4_all,soa_all]
    # # aerosol size and number
    # NCN3 = xr.DataArray(data=np.nansum(NCNall[3:, :], 0),  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="Aerosol number concentration for size >3nm",units="#/m3"),)    # >3nm
    # NCN10 = xr.DataArray(data=np.nansum(NCNall[10:, :], 0),  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="Aerosol number concentration for size >10nm",units="#/m3"),)     # >10nm
    # NCN100 = xr.DataArray(data=np.nansum(NCNall[100:, :], 0),  dims=["time"],
    #     coords=dict(time=(["time"], e3smtime)),
    #     attrs=dict(long_name="Aerosol number concentration for size >100nm",units="#/m3"),)     # >100nm
    # variable_names = variable_names + [ 'NCN3', 'NCN10', 'NCN100']
    # variables = variables + [NCN3, NCN10, NCN100]  # size distribution data will be added later
    # # mean cloud droplet number concentration
    # variable_names = variable_names + ['Nd_mean', 'Nd_ARM', 'Nd_VISST']
    # variables = variables + [cdnc_mean, cdnc_arm, cdnc_sat]
    # # mean cloud optical depth and effective radius
    # variable_names = variable_names + ['reff','cod']
    # variables = variables + [reff_mean, cod_mean]
    # # cloud depth
    # variable_names = variable_names + ['cbt','ctt','cbh','cth','clddepth']
    # variables = variables + [cbt, ctt, cbh, cth, cloud_depth]

    if config['composition_output'] == True:
        # aerosol composition    
        variable_names = variable_names + ['bc', 'dst', 'mom', 'ncl', 'pom', 'so4', 'soa']
        variables = variables + [bc_all, dst_all, mom_all, ncl_all, pom_all, so4_all, soa_all]
    if config['aerosol_output'] == True:
        # aerosol size and number
        NCN3 = xr.DataArray(data=np.nansum(NCNall[3:, :], 0),  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="Aerosol number concentration for size >3nm",units="#/m3"),)    # >3nm
        NCN10 = xr.DataArray(data=np.nansum(NCNall[10:, :], 0),  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="Aerosol number concentration for size >10nm",units="#/m3"),)     # >10nm
        NCN100 = xr.DataArray(data=np.nansum(NCNall[100:, :], 0),  dims=[config['time_dim']],
            coords=dict(time=([config['time_dim']], e3smtime)),
            attrs=dict(long_name="Aerosol number concentration for size >100nm",units="#/m3"),)     # >100nm
        variable_names = variable_names + [ 'NCN3', 'NCN10', 'NCN100']
        variables = variables + [NCN3, NCN10, NCN100]  # size distribution data will be added later
    if config['ccn_output'] == True:
        variable_names = variable_names + [ 'CCN1', 'CCN2']#, 'CCN5', 'CCN10']
        variables = variables + [ccn1_all, ccn2_all]#, ccn5_all, ccn10_all]
    # mean cloud droplet number concentration
    variable_names = variable_names + ['Nd_mean']
    variables = variables + [cdnc_mean]
    if config['tau3d_output'] == True:
        variable_names = variable_names + ['cod']
        variables = variables + [cod_mean]
        variable_names = variable_names + ['Nd_ARM']
        variables = variables + [cdnc_arm]
    elif config['tau3d_output'] == False and config['cosp_output'] == True:
        variable_names = variable_names + ['Nd_ARM']
        variables = variables + [cdnc_arm]
    if config['cosp_output'] == True:
        # variable_names = variable_names + ['cod_VISST']
        # variables = variables + [cod_m]
        variable_names = variable_names + ['Nd_VISST']
        variables = variables + [cdnc_sat]
    if config['reff_output'] == True:
        variable_names = variable_names + ['reff']
        variables = variables + [reff_mean]
        
    # cloud depth
    variable_names = variable_names + ['cbt','ctt','cbh','cth','clddepth']
    variables = variables + [cbt, ctt, cbh, cth, cloud_depth]

    #%% data treatments
    print(ttt.ctime(ttt.time())+': some data treatments:')
    lon_0 = median_time_1d(time, lon, e3smtime)
    lat_0 = median_time_1d(time, lat, e3smtime)
    x_idx_all = []
    for tt in range(len(e3smtime)):
        if np.isnan(lon_0[tt]):
            x_idx_all.append(-1)
        else:
            x_idx_tmp = find_nearest(lonm, latm, lon_0[tt], lat_0[tt])
            x_idx_all.append(x_idx_tmp)              
    # mask data with invalid ship location (x_idx_all == -1)
    for var in variables:
        var[np.array(x_idx_all) == -1] = np.nan
    
    # # change some units
    # # composition
    # T = variables[variable_names.index('T')]
    # ps = variables[variable_names.index('PS')]
    # rho = np.array(ps/T/287.06)
    # for vv in ['bc','dst','mom','ncl','pom','so4','soa']:
    #     variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1e9*rho
    #     variables[variable_names.index(vv)].attrs['units']='ug/m3'
    # # aerosol number
    # NCNall.data = NCNall.data * 1e-6
    # NCNall.attrs['units']='#/cm3'
    # for vv in ['NCN3', 'NCN10', 'NCN100', 'Nd_mean', 'Nd_ARM', 'Nd_VISST']:
    #     variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data * 1e-6
    #     variables[variable_names.index(vv)].attrs['units']='#/cm3'
    # # LWC and IWC
    # for vv in ['TGCLDIWP','TGCLDLWP']:
    #     variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1000
    #     variables[variable_names.index(vv)].attrs['units']='g/m2'
    # # cloud fraction
    # for vv in ['CLDTOT','CLDLOW','CLDMED','CLDHGH']:
    #     variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *100
    #     variables[variable_names.index(vv)].attrs['units']='%'

    #%% change some units
    for vv in ['Nd_mean']:
        variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data * 1e-6
        variables[variable_names.index(vv)].attrs['units']='#/cm3'
    if config['aerosol_output'] == True:
        # composition
        T = variables[variable_names.index(config['T'])]
        ps = variables[variable_names.index(config['PS'])]
        rho = np.array(ps/T/287.06)
        if config['composition_output'] == True:
            for vv in ['bc','dst','mom','ncl','pom','so4','soa']:
                variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1e9*rho
                variables[variable_names.index(vv)].attrs['units']='ug/m3'
        # aerosol number
        NCNall.data = NCNall.data * 1e-6
        NCNall.attrs['units']='#/cm3'
        for vv in ['NCN3', 'NCN10', 'NCN100']:
            variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data * 1e-6
            variables[variable_names.index(vv)].attrs['units']='#/cm3'
    if config['cosp_output'] == True:
        for vv in ['Nd_VISST']:
            variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data * 1e-6
            variables[variable_names.index(vv)].attrs['units']='#/cm3'
    if config['tau3d_output'] == True:
        for vv in ['Nd_ARM']:
            variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data * 1e-6
            variables[variable_names.index(vv)].attrs['units']='#/cm3'
    # LWC and IWC
    varbls = [config['IWP'],config['LWP']]
    varbls = ["{}{}".format(i,E3SMdomain_range) for i in varbls]
    varbls = list(set(av_vars).intersection(varbls))
    if len(varbls) > 0:
        for vv in [config['IWP'],config['LWP']]:
            variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1000
            variables[variable_names.index(vv)].attrs['units']='g/m2'
    # cloud fraction
    varbls = [config['CLDTOT']]#,config['CLDLOW'],config['CLDMED'],config['CLDHGH']]
    varbls = ["{}{}".format(i,E3SMdomain_range) for i in varbls]
    varbls = list(set(av_vars).intersection(varbls))
    if len(varbls) > 0:
        for vv in [config['CLDTOT']]:#,config['CLDLOW'],config['CLDMED'],config['CLDHGH']]:
            if variables[variable_names.index(vv)].ndim == 1:
                variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data*100
                variables[variable_names.index(vv)].attrs['long_name']='max cloud fraction in column'
                variables[variable_names.index(vv)].attrs['standard_name']='cloud fraction'
                variables[variable_names.index(vv)].attrs['units']='%'
            if variables[variable_names.index(vv)].ndim == 2:
                var2dto1d = variables[variable_names.index(vv)].max(dim=config['vert_dim'])*100 #max cloud fraction in the column
                variables[variable_names.index(vv)] = variables[variable_names.index(vv)][:,0] #retain only time dimension (time is first dimension based on dimension reordering earlier)
                variables[variable_names.index(vv)] = var2dto1d
                variables[variable_names.index(vv)].attrs['long_name']='max cloud fraction in column'
                variables[variable_names.index(vv)].attrs['standard_name']='cloud fraction'
                variables[variable_names.index(vv)].attrs['units']='%'
              
    # #%% re-shape the data into pre-defined resolution
    # print(ttt.ctime(ttt.time())+': re-shape to final time resolution:')
    # variables_new = list()
    # #1d variable. only numpy.interp can keep some single-point values (see Nd_mean)
    # for var in variables:
    #     # var_new = np.interp(np.int64(time_new), np.int64(e3smtime), var, left=np.nan, right=np.nan)
    #     var_new = np.interp(time_new, e3smtime, var, left=np.nan, right=np.nan)
    #     variables_new.append(var_new)
    # # treat variables with other dimensions (e.g., size distribution)    
    # f = interp1d(np.int64(e3smtime), NCNall, bounds_error=False)
    # NCNall_new = f(np.int64(time_new))

    variables_new = list()
    #1d variable. only numpy.interp can keep some single-point values (see Nd_mean)
    for var in variables:
        var_new = np.interp(np.int64(time_new), np.int64(e3smtime), var, left=np.nan, right=np.nan)
        variables_new.append(var_new)
    # treat variables with other dimensions (e.g., size distribution)
    if config['aerosol_output'] == True:
        f = interp1d(np.int64(e3smtime), NCNall, bounds_error=False)
        NCNall_new = f(np.int64(time_new))
      
    # add longtitude and latitude
    lon_new = xr.DataArray(data=lon_new,  dims=[config['time_dim']],
        coords=dict(time=([config['time_dim']], time_new)),
        attrs=dict(long_name="longitude",units="degree_east"),)
    lat_new = xr.DataArray(data=lat_new,  dims=[config['time_dim']],
        coords=dict(time=([config['time_dim']], time_new)),
        attrs=dict(long_name="latitude",units="degree_north"),)
    variable_names = variable_names + ['lon','lat']
    variables_new = variables_new + [lon_new, lat_new]
    variables = variables + [lon_new, lat_new]   # need this for attrs
    
    # %% output extacted file
    print('\nWriting file')
    varall_1d = {
            variable_names[vv]: (config['time_dim'], np.float32(variables_new[vv])) for vv in range(len(variable_names))
    }
    if config['aerosol_output'] == True:
        varall_2d = {
                'NCNall': (['size','time',], np.float32(NCNall_new))
        }
        varall_1d.update(varall_2d)
    outfile = output_path + output_filehead + '_ship.nc'
    print(ttt.ctime(ttt.time())+': output file '+outfile)
    ds = xr.Dataset( varall_1d,
                      coords={config['time_dim']: (config['time_dim'], time_new), 'size':('size', np.arange(1,3001))})
    #assign attributes
    ds[config['time_dim']].attrs["long_name"] = "Time"
    ds[config['time_dim']].attrs["standard_name"] = "time"
    for vv in range(len(variable_names)):
        ds[variable_names[vv]].attrs["long_name"] = variables[vv].long_name
        ds[variable_names[vv]].attrs["units"] = variables[vv].units
        ds[variable_names[vv]].attrs["description"] = "variables at surface, TOA or the lowest model level"
    if config['aerosol_output'] == True:
        ds['NCNall'].attrs["long_name"] = NCNall.long_name
        ds['NCNall'].attrs["units"] = NCNall.units
        ds['NCNall'].attrs["description"] = 'calculated from modal information into 1nm increment'
    
    ds.attrs["title"] = 'preprocessed E3SM data at surface, TOA or the lowest model level, along the ship track'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_profiles(input_path, input_filehead, output_path, output_filehead, shipmetpath, 
                      height_out, lev_out=np.arange(25.,1001,25.), dt=3600):
    """
    prepare vertical profile (to p or to z) variables from E3SM output along ship tracks
    choose the grid nearest to the ship location
    interpolate into new time with resolution of dt

    Parameters
    ----------
    input_path : str
        data path of E3SM output data.
    input_filehead : str
        filehead of the E3SM data, end at ".cam.h*" (E3SMv1) or ".eam.h*" (E3SMv2)
    output_path : str
        data path of output data..
    output_filehead : str
        filehead of the preprocessed E3SM data
    shipmetpath : str
        input path for ship location data
    lev_out : numpy array
        vertical coordinates of pressure, upside down
        a default value is given to be consistent with armbeatm
    height_out : numpy array
        vertical coordinates of height, upside down
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    #%% settings specific for each site
    E3SMdomain_range = '202e_to_243e_20n_to_35n'    # domain range in E3SM regional output
    
    # output time range and resolution
    time_new = pd.date_range(start='2012-10-05', end='2013-10-09 23:59:00', freq=str(int(dt))+"s")  # MAGIC time period
    
    #%% read in ship location data
    print(ttt.ctime(ttt.time())+': read in ship location data:')
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load().data
    lon = shipdata['lon_mean_gps'].load().data
    qc_lon = shipdata['qc_lon_mean_gps'].load().data
    lat = shipdata['lat_mean_gps'].load().data
    qc_lat = shipdata['qc_lat_mean_gps'].load().data
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    lon[lon<0] += 360
    
    lon_new = median_time_1d(time, lon, time_new)
    lat_new = median_time_1d(time, lat, time_new)
    
    #%% read in E3SM data
    print(ttt.ctime(ttt.time())+': read in E3SM data:')
    
    lst = glob.glob(input_path + input_filehead+'.*.h?.*-00000.nc') 
    lst.sort()
    
    e3smdata = xr.open_mfdataset(lst)
    e3smtime = e3smdata.indexes['time'].to_datetimeindex()
    lonm = e3smdata['lon'+'_'+E3SMdomain_range].load()
    latm = e3smdata['lat'+'_'+E3SMdomain_range].load()
    z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
    hyam = e3smdata['hyam'].load()
    hybm = e3smdata['hybm'].load()
    p0 = e3smdata['P0'].load()
    ps = e3smdata['PS'+'_'+E3SMdomain_range].load()
    Ts = e3smdata['TREFHT'+'_'+E3SMdomain_range].load()
    T = e3smdata['T'+'_'+E3SMdomain_range].load()
    Q = e3smdata['Q'+'_'+E3SMdomain_range].load()
    U = e3smdata['U'+'_'+E3SMdomain_range].load()
    V = e3smdata['V'+'_'+E3SMdomain_range].load()
    RH = e3smdata['RELHUM'+'_'+E3SMdomain_range].load()
    cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
    e3smdata.close()
    
    if len(lonm.shape)>1:   # when use xr.open_mfdataset, lat,lon will add a time dimension
        lonm = lonm[0,:].data
        latm = latm[0,:].data
        hyam = hyam[0,:].data
        hybm = hybm[0,:].data
        p0 = p0[0].data
    lon_0 = median_time_1d(time, lon, e3smtime)
    lat_0 = median_time_1d(time, lat, e3smtime)
    x_idx_all = []
    for tt in range(len(e3smtime)):
        if np.isnan(lon_0[tt]):
            x_idx_all.append(-1)
        else:
            x_idx_tmp = find_nearest(lonm, latm, lon_0[tt], lat_0[tt])
            x_idx_all.append(x_idx_tmp)
            
    ps = xr.concat([ps[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    Ts = xr.concat([Ts[i,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    z3 = xr.concat([z3[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    T = xr.concat([T[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    Q = xr.concat([Q[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    U = xr.concat([U[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    V = xr.concat([V[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    RH = xr.concat([RH[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    cloud = xr.concat([cloud[i,:,x_idx_all[i]] for i in range(len(e3smtime))],dim='time')
    
    levm = np.nan*T
    zlen = T.shape[1]
    for kk in range(zlen):
        levm[:, kk] = hyam[kk]*p0  +  hybm[kk]*ps
    levm = 0.01*levm  # hPa
    # calculate theta
    theta = T * (1000./levm)**0.286
    theta_s = Ts * (100000./ps)**0.286
    
    #%% interpolate data into pressure coordinate
    print(ttt.ctime(ttt.time())+': interpolate into final vertical coordinates:')
    cloud_p = np.empty((len(e3smtime),len(lev_out)))
    T_p = np.empty((len(e3smtime),len(lev_out)))
    Q_p = np.empty((len(e3smtime),len(lev_out)))
    RH_p = np.empty((len(e3smtime),len(lev_out)))
    theta_p = np.empty((len(e3smtime),len(lev_out)))
    z_p = np.empty((len(e3smtime),len(lev_out)))
    for tt in range(len(e3smtime)):
        cloud_p[tt,:] = np.interp(lev_out,levm[tt,:],cloud[tt,:])
        T_p[tt,:] = np.interp(lev_out,levm[tt,:],T[tt,:])
        Q_p[tt,:] = np.interp(lev_out,levm[tt,:],Q[tt,:])
        RH_p[tt,:] = np.interp(lev_out,levm[tt,:],RH[tt,:])
        theta_p[tt,:] = np.interp(lev_out,levm[tt,:],theta[tt,:])
        z_p[tt,:] = np.interp(lev_out,levm[tt,:],z3[tt,:])
            
    # interpolate data into height coordinate. flip model data since numpy.interp only works for increasing dimension
    cloud_z = np.empty((len(e3smtime),len(height_out)))
    T_z = np.empty((len(e3smtime),len(height_out)))
    RH_z = np.empty((len(e3smtime),len(height_out)))
    theta_z = np.empty((len(e3smtime),len(height_out)))
    Q_z = np.empty((len(e3smtime),len(height_out)))
    p_z = np.empty((len(e3smtime),len(height_out)))
    for tt in range(len(e3smtime)):
        cloud_z[tt,:] = np.interp(height_out,np.flip(z3[tt,:]),np.flip(cloud[tt,:]))
        T_z[tt,:] = np.interp(height_out,np.flip(z3[tt,:]),np.flip(T[tt,:]))
        RH_z[tt,:] = np.interp(height_out,np.flip(z3[tt,:]),np.flip(RH[tt,:]))
        theta_z[tt,:] = np.interp(height_out,np.flip(z3[tt,:]),np.flip(theta[tt,:]))
        Q_z[tt,:] = np.interp(height_out,np.flip(z3[tt,:]),np.flip(Q[tt,:]))
        p_z[tt,:] = np.interp(height_out,np.flip(z3[tt,:]),np.flip(levm[tt,:]))
        
    # lower tropospheric stability (theta diff between sfc and 700hPa)
    idx700 = 27
    if lev_out[idx700]!=700:
        print(lev_out[idx700])
        raise ValueError('the index of 700hPa has changed! check idx')
    LTS700 = theta_p[:, idx700] - theta_s    
    idx850 = 33
    if lev_out[idx850]!=850:
        print(lev_out[idx850])
        raise ValueError('the index of 850hPa has changed! check idx')
    LTS850 = theta_p[:, idx850] - theta_s  
    
        
    #%% re-shape the data into pre-defined resolution
    print(ttt.ctime(ttt.time())+': re-shape to final time resolution:')
    
    # treat variables with other dimensions (e.g., size distribution)    
    f = interp1d(np.int64(e3smtime), cloud_p, axis=0, bounds_error=False)
    cloud_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), T_p, axis=0, bounds_error=False)
    T_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), Q_p, axis=0, bounds_error=False)
    Q_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), RH_p, axis=0, bounds_error=False)
    RH_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), theta_p, axis=0, bounds_error=False)
    theta_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), z_p, axis=0, bounds_error=False)
    z_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), cloud_z, axis=0, bounds_error=False)
    cloud_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), T_z, axis=0, bounds_error=False)
    T_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), Q_z, axis=0, bounds_error=False)
    Q_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), RH_z, axis=0, bounds_error=False)
    RH_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), theta_z, axis=0, bounds_error=False)
    theta_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), p_z, axis=0, bounds_error=False)
    p_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), LTS700, axis=0, bounds_error=False)
    LTS700_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), LTS850, axis=0, bounds_error=False)
    LTS850_new = f(np.int64(time_new))
    
        
    # put all variables into the list
    # p-level
    cp = xr.DataArray(data=cloud_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=cloud.long_name, units=cloud.units),)
    tp = xr.DataArray(data=T_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=T.long_name, units=T.units),)
    qp = xr.DataArray(data=Q_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=Q.long_name, units=Q.units),)
    rhp = xr.DataArray(data=RH_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=RH.long_name, units=RH.units),)
    thp = xr.DataArray(data=theta_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name='Potential Temperature', units='K'),)
    zp = xr.DataArray(data=z_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=z3.long_name, units=z3.units),)
    varnames_p = [ 'cloud_p', 'T_p', 'Q_p', 'RH_p', 'theta_p', 'Z_p']
    variables_p = [   cp,      tp,     qp,    rhp,    thp,    zp]
    # z-level
    cz = xr.DataArray(data=cloud_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name=cloud.long_name, units=cloud.units),)
    tz = xr.DataArray(data=T_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name=T.long_name, units=T.units),)
    qz = xr.DataArray(data=Q_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name=Q.long_name, units=Q.units),)
    rhz = xr.DataArray(data=RH_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name=RH.long_name, units=RH.units),)
    thz = xr.DataArray(data=theta_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name='Potential Temperature', units='K'),)
    pz = xr.DataArray(data=p_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name='Pressure', units='hPa'),)
    varnames_z = [ 'cloud_z', 'T_z', 'Q_z', 'RH_z', 'theta_z', 'P_z']
    variables_z = [   cz,      tz,     qz,    rhz,    thz,    pz]
    #
    l700 = xr.DataArray(data=LTS700_new,  dims=["time"],
        coords=dict(time=(["time"], time_new), ),
        attrs=dict(long_name='lower troposphere stability (700hPa theta - surface theta)', units='K'),)
    l850 = xr.DataArray(data=LTS850_new,  dims=["time"],
        coords=dict(time=(["time"], time_new), ),
        attrs=dict(long_name='lower troposphere stability (850hPa theta - surface theta)', units='K'),)
    varnames_1d = [ 'LTS700', 'LTS850']
    variables_1d = [l700, l850]
    # add longtitude and latitude
    lon_new = xr.DataArray(data=lon_new,  dims=["time"],
        coords=dict(time=(["time"], time_new)),
        attrs=dict(long_name="longitude",units="degree_east"),)
    lat_new = xr.DataArray(data=lat_new,  dims=["time"],
        coords=dict(time=(["time"], time_new)),
        attrs=dict(long_name="latitude",units="degree_north"),)
    varnames_1d = varnames_1d + ['lon','lat']
    variables_1d = variables_1d + [lon_new, lat_new]
    
    # %% output extacted file
    varall_p = {
            varnames_p[vv]: (['time','lev',], np.float32(variables_p[vv])) for vv in range(len(varnames_p))
    }
    varall_z = {
            varnames_z[vv]: (['time','height',], np.float32(variables_z[vv])) for vv in range(len(varnames_z))
    }
    varall_1d = {
            varnames_1d[vv]: (['time',], np.float32(variables_1d[vv])) for vv in range(len(varnames_1d))
    }
    varall_1d.update(varall_p)
    varall_1d.update(varall_z)
    outfile = output_path + output_filehead + '_profiles.nc'
    print(ttt.ctime(ttt.time())+': output file '+outfile)
    ds = xr.Dataset( varall_1d,
                     coords={'time' : ('time', time_new), 
                             'lev' : ('lev', lev_out), 
                             'height' : ('height', height_out),
                             })
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    for vv in range(len(varnames_p)):
        ds[varnames_p[vv]].attrs["long_name"] = variables_p[vv].long_name
        ds[varnames_p[vv]].attrs["units"] = variables_p[vv].units
        ds[varnames_p[vv]].attrs["description"] = "interpolate into pressure level"
    for vv in range(len(varnames_z)):
        ds[varnames_z[vv]].attrs["long_name"] = variables_z[vv].long_name
        ds[varnames_z[vv]].attrs["units"] = variables_z[vv].units
        ds[varnames_z[vv]].attrs["description"] = "interpolate into height level"
    for vv in range(len(varnames_1d)):
        ds[varnames_1d[vv]].attrs["long_name"] = variables_1d[vv].long_name
        ds[varnames_1d[vv]].attrs["units"] = variables_1d[vv].units
        
    ds.attrs["title"] = 'preprocessed E3SM vertical profiles at ARM site'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
