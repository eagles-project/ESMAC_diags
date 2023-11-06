"""
prepare E3SM output data 
for surface (include remote sensing and satellite) measurements at SGP site
input data is regional hourly output data in E3SM
all variables has a region appendix similar to "lat_260e_to_265e_34n_to_39n"
"""

import glob
import os
import numpy as np
from scipy.interpolate import interp1d
import xarray as xr
import pandas as pd
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import  avg_time_2d
from esmac_diags.subroutines.specific_data_treatment import find_nearest, calc_Reff_from_REL, \
                            calc_cdnc_ARM, calc_cdnc_VISST
from esmac_diags.subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm


import warnings
warnings.filterwarnings("ignore")
# site = 'ENA'
# input_path = '../../../raw_data/model/'
# output_path = '../../../prep_data/ENA/model/'
# input_filehead = 'E3SMv2_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv2_ENA'
# dt=3600

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600):
    """
    prepare surface (include TOA and vertical integrated) variables from E3SM output
    choose the grid nearest to the ARM site
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
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    #%% settings specific for each site
    if site=='SGP':
        # SGP
        lat0 = 36.6059
        lon0 = -97.48792
        E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
        
        # output time range and resolution
        time_new = pd.date_range(start='2011-01-01', end='2020-12-31 23:59:00', freq=str(int(dt))+"s")   # SGP time period
    elif site=='ENA':
        # ENA
        lat0 = 39.09527
        lon0 = -28.0339
        E3SMdomain_range = '330e_to_335e_37n_to_42n'    # domain range in E3SM regional output
        
        # output time range and resolution
        time_new = pd.date_range(start='2016-01-01', end='2018-12-31 23:59:00', freq=str(int(dt))+"s")  # ENA time period
        
    
    #%% read in data
    variable_names = list()
    variables = list()
    
    if site=='SGP':
        lst = glob.glob(input_path + input_filehead+'.*.h?.*-00000.nc') 
    elif site=='ENA':
        lst = glob.glob(input_path + input_filehead+'.*.h?.*2016-*-00000.nc') + \
                glob.glob(input_path + input_filehead+'.*.h?.*2017-*-00000.nc') + \
                glob.glob(input_path + input_filehead+'.*.h?.*2018-*-00000.nc')
    lst.sort()
    # first data
    e3smdata = xr.open_dataset(lst[0])
    e3smtime = e3smdata.indexes['time'].to_datetimeindex()
    lonm = e3smdata['lon'+'_'+E3SMdomain_range].load()
    latm = e3smdata['lat'+'_'+E3SMdomain_range].load()
    ncols = len(lonm)
    
    # aerosol composition
    bc_a1 = e3smdata['bc_a1'+'_'+E3SMdomain_range].load()
    bc_a3 = e3smdata['bc_a3'+'_'+E3SMdomain_range].load()
    bc_a4 = e3smdata['bc_a4'+'_'+E3SMdomain_range].load()
    dst_a1 = e3smdata['dst_a1'+'_'+E3SMdomain_range].load()
    dst_a3 = e3smdata['dst_a3'+'_'+E3SMdomain_range].load()
    mom_a1 = e3smdata['mom_a1'+'_'+E3SMdomain_range].load()
    mom_a2 = e3smdata['mom_a2'+'_'+E3SMdomain_range].load()
    mom_a3 = e3smdata['mom_a3'+'_'+E3SMdomain_range].load()
    mom_a4 = e3smdata['mom_a4'+'_'+E3SMdomain_range].load()
    ncl_a1 = e3smdata['ncl_a1'+'_'+E3SMdomain_range].load()
    ncl_a2 = e3smdata['ncl_a2'+'_'+E3SMdomain_range].load()
    ncl_a3 = e3smdata['ncl_a3'+'_'+E3SMdomain_range].load()
    pom_a1 = e3smdata['pom_a1'+'_'+E3SMdomain_range].load()
    pom_a3 = e3smdata['pom_a3'+'_'+E3SMdomain_range].load()
    pom_a4 = e3smdata['pom_a4'+'_'+E3SMdomain_range].load()
    so4_a1 = e3smdata['so4_a1'+'_'+E3SMdomain_range].load()
    so4_a2 = e3smdata['so4_a2'+'_'+E3SMdomain_range].load()
    so4_a3 = e3smdata['so4_a3'+'_'+E3SMdomain_range].load()
    soa_a1 = e3smdata['soa_a1'+'_'+E3SMdomain_range].load()
    soa_a2 = e3smdata['soa_a2'+'_'+E3SMdomain_range].load()
    soa_a3 = e3smdata['soa_a3'+'_'+E3SMdomain_range].load()
    bc_all  = bc_a1[:,-1,:] +                       bc_a3[:,-1,:] + bc_a4[:,-1,:]
    dst_all = dst_a1[:,-1,:] +                      dst_a3[:,-1,:]
    mom_all = mom_a1[:,-1,:] + mom_a2[:,-1,:] + mom_a3[:,-1,:] + mom_a4[:,-1,:]
    ncl_all = ncl_a1[:,-1,:] + ncl_a2[:,-1,:] + ncl_a3[:,-1,:]
    pom_all = pom_a1[:,-1,:] +                    + pom_a3[:,-1,:] + pom_a4[:,-1,:]
    so4_all = so4_a1[:,-1,:] + so4_a2[:,-1,:] + so4_a3[:,-1,:]
    soa_all = soa_a1[:,-1,:] + soa_a2[:,-1,:] + soa_a3[:,-1,:]
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
    bc_all.coords['time'] = bc_all.indexes['time'].to_datetimeindex()
    dst_all.coords['time'] = dst_all.indexes['time'].to_datetimeindex()
    mom_all.coords['time'] = mom_all.indexes['time'].to_datetimeindex()
    ncl_all.coords['time'] = ncl_all.indexes['time'].to_datetimeindex()
    pom_all.coords['time'] = pom_all.indexes['time'].to_datetimeindex()
    so4_all.coords['time'] = so4_all.indexes['time'].to_datetimeindex()
    soa_all.coords['time'] = soa_all.indexes['time'].to_datetimeindex()
    
    # aerosol size
    num_a1 = e3smdata['num_a1'+'_'+E3SMdomain_range].load()
    num_a2 = e3smdata['num_a2'+'_'+E3SMdomain_range].load()
    num_a3 = e3smdata['num_a3'+'_'+E3SMdomain_range].load()
    num_a4 = e3smdata['num_a4'+'_'+E3SMdomain_range].load()
    dn1 = e3smdata['dgnd_a01'+'_'+E3SMdomain_range].load()
    dn2 = e3smdata['dgnd_a02'+'_'+E3SMdomain_range].load()
    dn3 = e3smdata['dgnd_a03'+'_'+E3SMdomain_range].load()
    dn4 = e3smdata['dgnd_a04'+'_'+E3SMdomain_range].load()
    P0 = e3smdata['P0'].load()
    hyam = e3smdata['hyam'].load()
    hybm = e3smdata['hybm'].load()
    T = e3smdata['T'+'_'+E3SMdomain_range].load()
    Ts = e3smdata['TREFHT'+'_'+E3SMdomain_range].load()
    PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
    Pres = np.nan*T
    zlen = T.shape[1]
    for kk in range(zlen):
        Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
    numall = [num_a1[:, -1, :].data, num_a2[:, -1, :].data, num_a3[:, -1, :].data, num_a4[:, -1, :].data]
    dnall  = [dn1[:, -1, :].data,    dn2[:, -1, :].data,    dn3[:, -1, :].data,    dn4[:, -1, :].data]
    NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1, :].data, Pres[:, -1, :].data)
    NCNall = xr.DataArray(data=NCN,  dims=["size", "time", "ncol"],
        coords=dict(time=(["time"], e3smtime),size=(["size"], np.arange(1,3001)),ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
    
    # variables to calculate Reff and Nd
    z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
    cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
    z3 = z3[:,:,:]
    cloud = cloud[:,:,:]
    dz = (z3[:,:-2,:].data - z3[:,2:,:].data)/2
    dz = np.append(dz, (z3[:,-2:-1,:].data+z3[:,-1:,:].data)/2, axis=1)
    dz = np.insert(dz,0,dz[:,0,:],axis=1)
    weight = cloud.data*dz
    # mid-level T, P, z
    Tmid = 0.5*(T[:,0:-1,:].data + T[:,1:,:].data)
    Pmid = 0.5*(Pres[:,0:-1,:].data + Pres[:,1:,:].data)
    Zmid = 0.5*(z3[:,0:-1,:].data + z3[:,1:,:].data)
    # cloud top
    cf_sum_top = np.cumsum(cloud.data, axis=1)
    cf_sum_top[cf_sum_top > 1] = 1
    cf_sum_top_diff = cf_sum_top[:,1:,:] - cf_sum_top[:,0:-1,:]
    z_cldtop = np.sum(Zmid[:,:,:]*cf_sum_top_diff[:,:,:], axis=1)
    T_cldtop = np.sum(Tmid[:,:,:]*cf_sum_top_diff[:,:,:], axis=1)
    # cloud base
    cf_sum_base = np.cumsum(cloud[:, ::-1,:].data, axis=1)
    cf_sum_base[cf_sum_base > 1] = 1
    cf_sum_base_diff = cf_sum_base[:,1:,:] - cf_sum_base[:,0:-1,:]
    z_cldbase = np.sum(Zmid[:,::-1,:]*cf_sum_base_diff[:,:,:], axis=1)
    p_cldbase = np.sum(Pmid[:,::-1,:]*cf_sum_base_diff[:,:,:], axis=1)
    T_cldbase = np.sum(Tmid[:,::-1,:]*cf_sum_base_diff[:,:,:], axis=1)
    # normalize cloud fraction when not 100%
    z_cldtop = z_cldtop / cf_sum_top[:,-1,:]
    T_cldtop = T_cldtop / cf_sum_top[:,-1,:]
    z_cldbase = z_cldbase / cf_sum_base[:,-1,:]
    T_cldbase = T_cldbase / cf_sum_base[:,-1,:]
    e3sm_cloud_depth = z_cldtop - z_cldbase
    # # alternatively, one can assume a minimum cloud depth in which the full mass levels are used instead of half levels
    # e3sm_z2 = z3[:,:-1].data
    # e3sm_z3 = z3[:,1:].data
    # e3sm_t2 = T[:,:-1,:].data
    # e3sm_t3 = T[:,1:,:].data
    # #cloud base and top calculations for minimum cloud depth
    # z_cldbase2 = np.sum(e3sm_z2[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # T_cldbase2 = np.sum(e3sm_t2[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # z_cldtop2 = np.sum(e3sm_z3[:,:]*cf_sum_top_diff[:,:], axis=1)
    # T_cldtop2 = np.sum(e3sm_t3[:,:]*cf_sum_top_diff[:,:], axis=1)
    # e3sm_cloud_depth2 = z_cldtop2 - z_cldbase2
    cloud_depth = xr.DataArray(data=e3sm_cloud_depth,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="cloud depth",units="m"),)
    cbh = xr.DataArray(data=z_cldbase,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="cloud base height",units="m"),)
    cth = xr.DataArray(data=z_cldtop,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="cloud top height",units="m"),)
    cbt = xr.DataArray(data=T_cldbase,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="cloud base temperature",units="K"),)
    ctt = xr.DataArray(data=T_cldtop,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="cloud top temperature",units="K"),)
    thetadiff_cb = xr.DataArray(data = T_cldbase - Ts + 9.8/1005.7*z_cldbase,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="Theta diff between cloud base and surface",units="K"),)
    
    # cloud optical depth and effective radius
    rel = e3smdata['REL'+'_'+E3SMdomain_range].load()
    freql = e3smdata['FREQL'+'_'+E3SMdomain_range].load()
    icwnc = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
    cod_a = e3smdata['TOT_CLD_VISTAU'+'_'+E3SMdomain_range].load()
    cod_m = e3smdata['TAUWMODIS'+'_'+E3SMdomain_range].load()*0.01   # need a correction due to CF problem
    solin = e3smdata['SOLIN'+'_'+E3SMdomain_range].load()
    # calculate mean effective radius. 
    cloud_tmp = np.array(freql)
    # filter out small cloud fraction
    for jj in range(len(lonm)):
        for ii in range(freql.shape[1]):
            # cloud fraction <= 1% or in cloud Nd < 1 cm-3
            idx = np.logical_or(freql[:,ii,jj]<=0.01, icwnc[:,ii,jj]<1e-6)
            cloud_tmp[idx, ii,jj] = 0
    weight2 = cloud_tmp * dz
    reff = np.nansum(rel.data*weight2,axis=1)/np.sum(weight2,axis=1)
    reff[reff==0] = np.nan
    # calculate mean optical depth
    cod = np.sum(cod_a.data,axis=1)
    cod[solin.data==0] = np.nan
    # cod from MODIS simulator
    cod_m = cod_m.data
    reff_mean = xr.DataArray(data=reff,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
    cod_mean = xr.DataArray(data=cod,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
    
    # mean cloud droplet number concentration
    cdnc_col = e3smdata['CDNUMC'+'_'+E3SMdomain_range].load()
    cdnc_mean = cdnc_col.data/np.sum(weight,axis=1)
    cdnc_mean[cdnc_mean >1e10] = np.nan
    cdnc_mean = xr.DataArray(data=cdnc_mean,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
    
    # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
    lwp = e3smdata['TGCLDLWP'+'_'+E3SMdomain_range][:,:].data
    e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
    T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
    nd_arm = calc_cdnc_ARM(lwp, cod_m, e3sm_cloud_depth)
    nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m, adiabaticity=0.8)
    cdnc_arm = xr.DataArray(data=nd_arm*1e6,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                    description='Retrieved using ARM Ndrop algorithm'),)
    cdnc_sat = xr.DataArray(data=nd_sat*1e6,  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                    description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)
    
    # all other 2D (surface and vertical integrated) variables
    variable2d_names = ['AODABS', 'AODALL', 'CDNUMC', 'CLDHGH', 'CLDMED', 'CLDLOW', 'CLDTOT', 
                        'FLDS', 'FLNS', 'FLNT', 'FLUT', 'FSDS', 'FSNS', 'FSNT', 'FSUTOA', 'SOLIN', 
                        'LHFLX', 'SHFLX', 'LWCF', 'SWCF', "TGCLDLWP", "TGCLDIWP", 
                        'IWPMODIS', 'LWPMODIS', 'REFFCLWMODIS', 'TAUIMODIS', 'TAUTMODIS', 'TAUWMODIS', 
                       # 'MEANPTOP_ISCCP', 'MEANCLDALB_ISCCP', 'MEANTAU_ISCCP',
                        'PBLH', 'PRECT', 'PRECL', 'PRECC', 'PS', 'TREFHT', ]
    for varname in variable2d_names:
        var = e3smdata[varname + '_'+E3SMdomain_range].load()
        var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
        if varname=='AODABS' or varname=='AODALL':
            var.attrs['units']='N/A'
        variable_names.append(varname)
        variables.append(var)
    
    # all other 3D (with vertical level) variables at the lowest model level
    variable3d_names = ['CCN1', 'CCN3', 'CCN4', 'CCN5', 'Q', 'T', 'RELHUM', 'U', 'V'] 
    for varname in variable3d_names:
        var = e3smdata[varname + '_'+E3SMdomain_range].load()
        var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
        variables.append(var[:,-1,:])
        variable_names.append(varname)
    
    e3smdata.close()
    
    #%%  add data for each day
    for file in lst[1:]:
        print(file)
        e3smdata = xr.open_dataset(file)
        e3smtime_i = e3smdata.indexes['time'].to_datetimeindex()
        e3smtime = np.hstack((e3smtime, e3smtime_i))
        
        # aerosol composition
        bc_a1 = e3smdata['bc_a1'+'_'+E3SMdomain_range].load()
        bc_a3 = e3smdata['bc_a3'+'_'+E3SMdomain_range].load()
        bc_a4 = e3smdata['bc_a4'+'_'+E3SMdomain_range].load()
        dst_a1 = e3smdata['dst_a1'+'_'+E3SMdomain_range].load()
        dst_a3 = e3smdata['dst_a3'+'_'+E3SMdomain_range].load()
        mom_a1 = e3smdata['mom_a1'+'_'+E3SMdomain_range].load()
        mom_a2 = e3smdata['mom_a2'+'_'+E3SMdomain_range].load()
        mom_a3 = e3smdata['mom_a3'+'_'+E3SMdomain_range].load()
        mom_a4 = e3smdata['mom_a4'+'_'+E3SMdomain_range].load()
        ncl_a1 = e3smdata['ncl_a1'+'_'+E3SMdomain_range].load()
        ncl_a2 = e3smdata['ncl_a2'+'_'+E3SMdomain_range].load()
        ncl_a3 = e3smdata['ncl_a3'+'_'+E3SMdomain_range].load()
        pom_a1 = e3smdata['pom_a1'+'_'+E3SMdomain_range].load()
        pom_a3 = e3smdata['pom_a3'+'_'+E3SMdomain_range].load()
        pom_a4 = e3smdata['pom_a4'+'_'+E3SMdomain_range].load()
        so4_a1 = e3smdata['so4_a1'+'_'+E3SMdomain_range].load()
        so4_a2 = e3smdata['so4_a2'+'_'+E3SMdomain_range].load()
        so4_a3 = e3smdata['so4_a3'+'_'+E3SMdomain_range].load()
        soa_a1 = e3smdata['soa_a1'+'_'+E3SMdomain_range].load()
        soa_a2 = e3smdata['soa_a2'+'_'+E3SMdomain_range].load()
        soa_a3 = e3smdata['soa_a3'+'_'+E3SMdomain_range].load()
        bc  = bc_a1[:,-1,:] +                       bc_a3[:,-1,:] + bc_a4[:,-1,:]
        dst = dst_a1[:,-1,:] +                      dst_a3[:,-1,:]
        mom = mom_a1[:,-1,:] + mom_a2[:,-1,:] + mom_a3[:,-1,:] + mom_a4[:,-1,:]
        ncl = ncl_a1[:,-1,:] + ncl_a2[:,-1,:] + ncl_a3[:,-1,:]
        pom = pom_a1[:,-1,:] +                    + pom_a3[:,-1,:] + pom_a4[:,-1,:]
        so4 = so4_a1[:,-1,:] + so4_a2[:,-1,:] + so4_a3[:,-1,:]
        soa = soa_a1[:,-1,:] + soa_a2[:,-1,:] + soa_a3[:,-1,:]
        # change time to standard datetime64 format
        bc.coords['time'] = bc.indexes['time'].to_datetimeindex()
        dst.coords['time'] = dst.indexes['time'].to_datetimeindex()
        mom.coords['time'] = mom.indexes['time'].to_datetimeindex()
        ncl.coords['time'] = ncl.indexes['time'].to_datetimeindex()
        pom.coords['time'] = pom.indexes['time'].to_datetimeindex()
        so4.coords['time'] = so4.indexes['time'].to_datetimeindex()
        soa.coords['time'] = soa.indexes['time'].to_datetimeindex()
        bc_all = xr.concat([bc_all, bc], dim="time")
        dst_all = xr.concat([dst_all, dst], dim="time")
        mom_all = xr.concat([mom_all, mom], dim="time")
        ncl_all = xr.concat([ncl_all, ncl], dim="time")
        pom_all = xr.concat([pom_all, pom], dim="time")
        so4_all = xr.concat([so4_all, so4], dim="time")
        soa_all = xr.concat([soa_all, soa], dim="time")
        
        # aerosol size
        num_a1 = e3smdata['num_a1'+'_'+E3SMdomain_range].load()
        num_a2 = e3smdata['num_a2'+'_'+E3SMdomain_range].load()
        num_a3 = e3smdata['num_a3'+'_'+E3SMdomain_range].load()
        num_a4 = e3smdata['num_a4'+'_'+E3SMdomain_range].load()
        dn1 = e3smdata['dgnd_a01'+'_'+E3SMdomain_range].load()
        dn2 = e3smdata['dgnd_a02'+'_'+E3SMdomain_range].load()
        dn3 = e3smdata['dgnd_a03'+'_'+E3SMdomain_range].load()
        dn4 = e3smdata['dgnd_a04'+'_'+E3SMdomain_range].load()
        P0 = e3smdata['P0'].load()
        hyam = e3smdata['hyam'].load()
        hybm = e3smdata['hybm'].load()
        T = e3smdata['T'+'_'+E3SMdomain_range].load()
        Ts = e3smdata['TREFHT'+'_'+E3SMdomain_range].load()
        PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
        Pres = np.nan*T
        zlen = T.shape[1]
        for kk in range(zlen):
            Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
        numall = [num_a1[:, -1, :].data, num_a2[:, -1, :].data, num_a3[:, -1, :].data, num_a4[:, -1, :].data]
        dnall  = [dn1[:, -1, :].data,    dn2[:, -1, :].data,    dn3[:, -1, :].data,    dn4[:, -1, :].data]
        NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1, :].data, Pres[:, -1, :].data)
        NCN2 = xr.DataArray(data=NCN,  dims=["size", "time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i),size=(["size"], np.arange(1,3001)), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
        NCNall = xr.concat([NCNall, NCN2], dim="time")
        
        # variables to calculate Reff and Nd
        z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
        cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
        z3 = z3[:,:,:]
        cloud = cloud[:,:,:]
        dz = (z3[:,:-2,:].data - z3[:,2:,:].data)/2
        dz = np.append(dz, (z3[:,-2:-1,:].data+z3[:,-1:,:].data)/2, axis=1)
        dz = np.insert(dz,0,dz[:,0,:],axis=1)
        weight = cloud.data*dz
        # mid-level T, P, z
        Tmid = 0.5*(T[:,0:-1,:].data + T[:,1:,:].data)
        Pmid = 0.5*(Pres[:,0:-1,:].data + Pres[:,1:,:].data)
        Zmid = 0.5*(z3[:,0:-1,:].data + z3[:,1:,:].data)
        # cloud top
        cf_sum_top = np.cumsum(cloud.data, axis=1)
        cf_sum_top[cf_sum_top > 1] = 1
        cf_sum_top_diff = cf_sum_top[:,1:,:] - cf_sum_top[:,0:-1,:]
        z_cldtop = np.sum(Zmid[:,:,:]*cf_sum_top_diff[:,:,:], axis=1)
        T_cldtop = np.sum(Tmid[:,:,:]*cf_sum_top_diff[:,:,:], axis=1)
        # cloud base
        cf_sum_base = np.cumsum(cloud[:, ::-1,:].data, axis=1)
        cf_sum_base[cf_sum_base > 1] = 1
        cf_sum_base_diff = cf_sum_base[:,1:,:] - cf_sum_base[:,0:-1,:]
        z_cldbase = np.sum(Zmid[:,::-1,:]*cf_sum_base_diff[:,:,:], axis=1)
        T_cldbase = np.sum(Tmid[:,::-1,:]*cf_sum_base_diff[:,:,:], axis=1)
        # normalize cloud fraction when not 100%
        z_cldtop = z_cldtop / cf_sum_top[:,-1,:]
        T_cldtop = T_cldtop / cf_sum_top[:,-1,:]
        z_cldbase = z_cldbase / cf_sum_base[:,-1,:]
        T_cldbase = T_cldbase / cf_sum_base[:,-1,:]
        e3sm_cloud_depth = z_cldtop - z_cldbase
        cloud_depth_2 = xr.DataArray(data=e3sm_cloud_depth,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="cloud depth",units="m"),)
        cbh_2 = xr.DataArray(data=z_cldbase,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="cloud base height",units="m"),)
        cth_2 = xr.DataArray(data=z_cldtop,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="cloud top height",units="m"),)
        cbt_2 = xr.DataArray(data=T_cldbase,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="cloud base temperature",units="K"),)
        ctt_2 = xr.DataArray(data=T_cldtop,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="cloud top temperature",units="K"),)
        thetadiff_cb_2 = xr.DataArray(data = T_cldbase - Ts + 9.8/1005.7*z_cldbase,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="Theta diff between cloud base and surface",units="K"),)
        cloud_depth = xr.concat([cloud_depth, cloud_depth_2], dim="time")
        cbh = xr.concat([cbh, cbh_2], dim="time")
        cth = xr.concat([cth, cth_2], dim="time")
        cbt = xr.concat([cbt, cbt_2], dim="time")
        ctt = xr.concat([ctt, ctt_2], dim="time")
        thetadiff_cb = xr.concat([thetadiff_cb, thetadiff_cb_2], dim="time")
        
        # cloud optical depth and effective radius
        rel = e3smdata['REL'+'_'+E3SMdomain_range].load()
        freql = e3smdata['FREQL'+'_'+E3SMdomain_range].load()
        icwnc = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
        cod_a = e3smdata['TOT_CLD_VISTAU'+'_'+E3SMdomain_range].load()
        cod_m = e3smdata['TAUWMODIS'+'_'+E3SMdomain_range].load()*0.01   # correction for CF problem
        solin = e3smdata['SOLIN'+'_'+E3SMdomain_range].load()
        # calculate mean effective radius. 
        cloud_tmp = np.array(freql)
        # filter out small cloud fraction
        for jj in range(len(lonm)):
            for ii in range(freql.shape[1]):
                # cloud fraction <= 1% or in cloud Nd < 1 cm-3
                idx = np.logical_or(freql[:,ii,jj]<=0.01, icwnc[:,ii,jj]<1e-6)
                cloud_tmp[idx, ii,jj] = 0
        weight2 = cloud_tmp * dz
        reff = np.nansum(rel.data*weight2,axis=1)/np.sum(weight2,axis=1)
        reff[reff==0] = np.nan
        # calculate mean optical depth
        cod = np.sum(cod_a.data,axis=1)
        cod[solin.data==0] = np.nan
        # cod from MODIS simulator
        cod_m = cod_m.data
        reff_2 = xr.DataArray(data=reff,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
        cod_2 = xr.DataArray(data=cod,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
        reff_mean = xr.concat([reff_mean, reff_2], dim="time")
        cod_mean = xr.concat([cod_mean, cod_2], dim="time")
        
        # mean cloud droplet number concentration
        cdnc_col = e3smdata['CDNUMC'+'_'+E3SMdomain_range].load()
        cdnc = cdnc_col.data/np.sum(weight,axis=1)
        cdnc[cdnc >1e10] = np.nan
        cdnc = xr.DataArray(data=cdnc,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
        cdnc_mean = xr.concat([cdnc_mean, cdnc], dim="time")
        
        # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
        lwp = e3smdata['TGCLDLWP'+'_'+E3SMdomain_range][:,:].data
        e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
        T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
        nd_arm = calc_cdnc_ARM(lwp, cod_m, e3sm_cloud_depth)
        nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m, adiabaticity=0.8)
        nd_arm = xr.DataArray(data=nd_arm*1e6,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                        description='Retrieved using ARM Ndrop algorithm'),)
        nd_sat = xr.DataArray(data=nd_sat*1e6,  dims=["time", "ncol"],
            coords=dict(time=(["time"], e3smtime_i), ncol=(["ncol"], np.arange(1,ncols+1))),
            attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                        description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)
        cdnc_arm = xr.concat([cdnc_arm, nd_arm], dim="time")
        cdnc_sat = xr.concat([cdnc_sat, nd_sat], dim="time")
        
        # all other 2D (surface and vertical integrated) variables
        for varname in variable2d_names:
            var = e3smdata[varname + '_'+E3SMdomain_range].load()
            var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
            vv = variable_names.index(varname)
            variables[vv] = xr.concat([variables[vv], var],dim='time')
        
        # all other 3D (with vertical level) variables at the lowest model level
        for varname in variable3d_names:
            var = e3smdata[varname + '_'+E3SMdomain_range].load()
            var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
            vv = variable_names.index(varname)
            variables[vv] = xr.concat([variables[vv], var[:,-1,:]],dim='time')
            
        e3smdata.close()
    
    # put all variables into the list
    # aerosol composition    
    variable_names = variable_names + ['bc','dst','mom','ncl','pom','so4','soa']
    variables = variables + [bc_all,dst_all,mom_all,ncl_all,pom_all,so4_all,soa_all]
    # aerosol size and number
    NCN3 = xr.DataArray(data=np.nansum(NCNall[3:, :, :], 0),  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="Aerosol number concentration for size >3nm",units="#/m3"),)    # >3nm
    NCN10 = xr.DataArray(data=np.nansum(NCNall[10:, :, :], 0),  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="Aerosol number concentration for size >10nm",units="#/m3"),)     # >10nm
    NCN100 = xr.DataArray(data=np.nansum(NCNall[100:, :, :], 0),  dims=["time", "ncol"],
        coords=dict(time=(["time"], e3smtime), ncol=(["ncol"], np.arange(1,ncols+1))),
        attrs=dict(long_name="Aerosol number concentration for size >100nm",units="#/m3"),)     # >100nm
    variable_names = variable_names + [ 'NCN3', 'NCN10', 'NCN100']
    variables = variables + [NCN3, NCN10, NCN100]  # size distribution data will be added later
    # mean cloud droplet number concentration
    variable_names = variable_names + ['Nd_mean', 'Nd_ARM', 'Nd_VISST']
    variables = variables + [cdnc_mean, cdnc_arm, cdnc_sat]
    # mean cloud optical depth and effective radius
    variable_names = variable_names + ['reff','cod']
    variables = variables + [reff_mean, cod_mean]
    # cloud depth
    variable_names = variable_names + ['cbt','ctt','cbh','cth','clddepth','thetadiff_cb']
    variables = variables + [cbt,ctt,cbh,cth,cloud_depth,thetadiff_cb]
    
    #%% change some units
    # composition
    T = variables[variable_names.index('T')]
    ps = variables[variable_names.index('PS')]
    rho = np.array(ps/T/287.06)
    for vv in ['bc','dst','mom','ncl','pom','so4','soa']:
        variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1e9*rho
        variables[variable_names.index(vv)].attrs['units']='ug/m3'
    # aerosol number
    NCNall.data = NCNall.data * 1e-6
    NCNall.attrs['units']='#/cm3'
    for vv in ['NCN3', 'NCN10', 'NCN100', 'Nd_mean', 'Nd_ARM', 'Nd_VISST']:
        variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data * 1e-6
        variables[variable_names.index(vv)].attrs['units']='#/cm3'
    # LWC and IWC
    for vv in ['TGCLDIWP','TGCLDLWP']:
        variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1000
        variables[variable_names.index(vv)].attrs['units']='g/m2'
    # cloud fraction
    for vv in ['CLDTOT','CLDLOW','CLDMED','CLDHGH']:
        variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *100
        variables[variable_names.index(vv)].attrs['units']='%'
    
    #%% re-shape the data into pre-defined resolution
    variables_new = list()
    #1d variable. only numpy.interp can keep some single-point values (see Nd_mean)
    for var in variables:
        var_new = avg_time_2d(e3smtime, var, time_new)
        variables_new.append(var_new)
        
    
    # %% output extacted file
    varall_1d = {
            variable_names[vv]: (['time','ncol',], np.float32(variables_new[vv])) for vv in range(len(variable_names))
    }
    lonlat = {'lon' : ('ncol', np.float32(lonm)), 'lat' : ('ncol', np.float32(latm))}
    varall_1d.update(lonlat)
    variable_names = variable_names + ['lon','lat']
    variables = variables + [lonm, latm]
     
    outfile = output_path + output_filehead + '_sfc_allgrids.nc'
    print('output file '+outfile)
    ds = xr.Dataset( varall_1d,
                     coords={'time' : ('time', time_new), 
                             'ncol':('ncol', np.arange(1,ncols+1)),
                             })
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    for vv in range(len(variable_names)):
        ds[variable_names[vv]].attrs["long_name"] = variables[vv].long_name
        ds[variable_names[vv]].attrs["units"] = variables[vv].units
        ds[variable_names[vv]].attrs["description"] = "variables at surface, TOA or the lowest model level"
    
    ds.attrs["title"] = 'preprocessed E3SM data at surface, TOA or the lowest model level'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

# site = 'ENA'
# input_path = '../../../raw_data/model/'
# output_path = '../../../prep_data/ENA/model/'
# input_filehead = 'E3SMv1_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv1_ENA'
# dt=3600
# prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)
# input_filehead = 'E3SMv2_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv2_ENA'
# prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)

# site = 'SGP'
# input_path = '../../../raw_data/model/'
# output_path = '../../../prep_data/SGP/model/'
# input_filehead = 'E3SMv1_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv1_SGP'
# dt=3600
# prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)
# input_filehead = 'E3SMv2_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv2_SGP'
# prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)
