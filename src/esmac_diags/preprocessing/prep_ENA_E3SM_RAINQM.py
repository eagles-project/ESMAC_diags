"""
prepare E3SM output data 
for surface (include remote sensing and satellite) measurements at ENA site
input data is regional hourly output data in E3SM
all variables has a region appendix similar to "lat_330e_to_335e_37n_to_42n"
"""

import glob
import os
import numpy as np
from scipy.interpolate import interp1d
import xarray as xr
import pandas as pd
import time as ttt
import esmac_diags
from esmac_diags.subroutines.specific_data_treatment import find_nearest, calc_Reff_from_REL, \
                            calc_cdnc_ARM, calc_cdnc_VISST
from esmac_diags.subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm

# input_path = '../../../raw_data/model/'
# output_path = '../../../prep_data/ENA/model/'
# input_filehead = 'E3SMv1_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv1_ENA'
# dt=3600
# height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
#                     1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
#                     3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
#                     10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])
# lev_out=np.arange(25.,1001,25.)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_profiles(input_path, input_filehead, output_path, output_filehead, 
                      height_out, lev_out=np.arange(25.,1001,25.), dt=3600):
    """
    prepare vertical profile (to p or to z) variables from E3SM output
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
    # ENA
    lat0 = 39.09527
    lon0 = -28.0339
    E3SMdomain_range = '330e_to_335e_37n_to_42n'    # domain range in E3SM regional output
    
    # output time range and resolution
    time_new = pd.date_range(start='2016-01-01', end='2018-12-31 23:59:00', freq=str(int(dt))+"s")  # ENA time period
    
    #%% read in data
    
    lst = glob.glob(input_path + input_filehead+'.*2016-*-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*2017-*-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*2018-*-00000.nc')
    lst.sort()
    # first data
    e3smdata = xr.open_dataset(lst[0])
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
    lwc = e3smdata['RAINQM'+'_'+E3SMdomain_range].load()
    iwc = e3smdata['SNOWQM'+'_'+E3SMdomain_range].load()
    e3smdata.close()
    # only extract the model column at the site
    if lon0<0:
        lon0=lon0+360   # make longitude consistent with E3SM from 0 to 360
    x_idx = find_nearest(lonm,latm,lon0,lat0)
    
    Pres = np.nan*T
    zlen = T.shape[1]
    for kk in range(zlen):
        Pres[:, kk, :] = hyam[kk]*p0  +  hybm[kk]*ps
    levm = 0.01* (ps[:,x_idx]*hybm + hyam*p0)  # hPa
    # calculate theta
    theta = T * (1000./levm)**0.286
    theta_s = Ts[:, x_idx] * (100000./ps[:, x_idx])**0.286
    
    # interpolate data into pressure coordinate
    cloud_p = np.empty((cloud.shape[0],len(lev_out)))
    T_p = np.empty((T.shape[0],len(lev_out)))
    Q_p = np.empty((Q.shape[0],len(lev_out)))
    RH_p = np.empty((RH.shape[0],len(lev_out)))
    theta_p = np.empty((T.shape[0],len(lev_out)))
    z_p = np.empty((T.shape[0],len(lev_out)))
    lwc_p = np.empty((lwc.shape[0],len(lev_out)))
    iwc_p = np.empty((iwc.shape[0],len(lev_out)))
    for i in range(len(e3smtime)):
        cloud_p[i,:] = np.interp(lev_out,levm[i,:],cloud[i,:,x_idx])
        T_p[i,:] = np.interp(lev_out,levm[i,:],T[i,:,x_idx])
        Q_p[i,:] = np.interp(lev_out,levm[i,:],Q[i,:,x_idx])
        RH_p[i,:] = np.interp(lev_out,levm[i,:],RH[i,:,x_idx])
        theta_p[i,:] = np.interp(lev_out,levm[i,:],theta[i,:,x_idx])
        z_p[i,:] = np.interp(lev_out,levm[i,:],z3[i,:,x_idx])
        lwc_p[i,:] = np.interp(lev_out,levm[i,:],lwc[i,:,x_idx])
        iwc_p[i,:] = np.interp(lev_out,levm[i,:],iwc[i,:,x_idx])
            
    # interpolate data into height coordinate. flip model data since numpy.interp only works for increasing dimension
    cloud_z = np.empty((len(e3smtime),len(height_out)))
    T_z = np.empty((len(e3smtime),len(height_out)))
    RH_z = np.empty((len(e3smtime),len(height_out)))
    theta_z = np.empty((len(e3smtime),len(height_out)))
    Q_z = np.empty((len(e3smtime),len(height_out)))
    p_z = np.empty((len(e3smtime),len(height_out)))
    lwc_z = np.empty((len(e3smtime),len(height_out)))
    iwc_z = np.empty((len(e3smtime),len(height_out)))
    for i in range(len(e3smtime)):
        cloud_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(cloud[i,:,x_idx]))
        T_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(T[i,:,x_idx]))
        RH_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(RH[i,:,x_idx]))
        theta_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(theta[i,:,x_idx]))
        Q_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(Q[i,:,x_idx]))
        p_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(levm[i,:]))
        lwc_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(lwc[i,:,x_idx]))
        iwc_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(iwc[i,:,x_idx]))
        
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
    
    # mid-level T, P, z
    Tmid = 0.5*(T[:,0:-1,x_idx].data + T[:,1:,x_idx].data)
    Pmid = 0.5*(Pres[:,0:-1,x_idx].data + Pres[:,1:,x_idx].data)
    # cloud base
    cf_sum_base = np.cumsum(cloud[:, ::-1,x_idx].data, axis=1)
    cf_sum_base[cf_sum_base > 1] = 1
    cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
    p_cldbase = np.sum(Pmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
    # normalize cloud fraction when not 100%
    p_cldbase = p_cldbase / cf_sum_base[:,-1]
    T_cldbase = T_cldbase / cf_sum_base[:,-1]
    theta_cb = T_cldbase * (100000./p_cldbase)**0.286
    thetadiff_cb = theta_cb - theta_s
    
    #%%  add data for each day
    for file in lst[1:]:
        print(file)
        e3smdata = xr.open_dataset(file)
        e3smtime_i = e3smdata.indexes['time'].to_datetimeindex()
        e3smtime = np.hstack((e3smtime, e3smtime_i))
        
        z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
        ps = e3smdata['PS'+'_'+E3SMdomain_range].load()
        Ts = e3smdata['TREFHT'+'_'+E3SMdomain_range].load()
        T = e3smdata['T'+'_'+E3SMdomain_range].load()
        Q = e3smdata['Q'+'_'+E3SMdomain_range].load()
        U = e3smdata['U'+'_'+E3SMdomain_range].load()
        V = e3smdata['V'+'_'+E3SMdomain_range].load()
        RH = e3smdata['RELHUM'+'_'+E3SMdomain_range].load()
        cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
        lwc = e3smdata['RAINQM'+'_'+E3SMdomain_range].load()
        iwc = e3smdata['SNOWQM'+'_'+E3SMdomain_range].load()
        e3smdata.close()
        
        Pres = np.nan*T
        zlen = T.shape[1]
        for kk in range(zlen):
            Pres[:, kk, :] = hyam[kk]*p0  +  hybm[kk]*ps
        levm = 0.01* (ps[:,x_idx]*hybm + hyam*p0)  # hPa
        # calculate theta
        theta = T * (1000./levm)**0.286
        theta_s2 = Ts[:, x_idx] * (100000./ps[:, x_idx])**0.286
    
        # interpolate data into pressure coordinate
        cloud_p2 = np.empty((len(e3smtime_i),len(lev_out)))
        T_p2 = np.empty((len(e3smtime_i),len(lev_out)))
        Q_p2 = np.empty((len(e3smtime_i),len(lev_out)))
        RH_p2 = np.empty((len(e3smtime_i),len(lev_out)))
        theta_p2 = np.empty((len(e3smtime_i),len(lev_out)))
        z_p2 = np.empty((T.shape[0],len(lev_out)))
        lwc_p2 = np.empty((lwc.shape[0],len(lev_out)))
        iwc_p2 = np.empty((iwc.shape[0],len(lev_out)))
        for i in range(len(e3smtime_i)):
            cloud_p2[i,:] = np.interp(lev_out,levm[i,:],cloud[i,:,x_idx])
            T_p2[i,:] = np.interp(lev_out,levm[i,:],T[i,:,x_idx])
            Q_p2[i,:] = np.interp(lev_out,levm[i,:],Q[i,:,x_idx])
            RH_p2[i,:] = np.interp(lev_out,levm[i,:],RH[i,:,x_idx])
            theta_p2[i,:] = np.interp(lev_out,levm[i,:],theta[i,:,x_idx])
            z_p2[i,:] = np.interp(lev_out,levm[i,:],z3[i,:,x_idx])
            lwc_p2[i,:] = np.interp(lev_out,levm[i,:],lwc[i,:,x_idx])
            iwc_p2[i,:] = np.interp(lev_out,levm[i,:],iwc[i,:,x_idx])
        
        # interpolate data into height coordinate
        cloud_z2 = np.empty((len(e3smtime_i),len(height_out)))
        T_z2 = np.empty((len(e3smtime_i),len(height_out)))
        RH_z2 = np.empty((len(e3smtime_i),len(height_out)))
        theta_z2 = np.empty((len(e3smtime_i),len(height_out)))
        Q_z2 = np.empty((len(e3smtime_i),len(height_out)))
        p_z2 = np.empty((len(e3smtime_i),len(height_out)))
        lwc_z2 = np.empty((len(e3smtime_i),len(height_out)))
        iwc_z2 = np.empty((len(e3smtime_i),len(height_out)))
        for i in range(len(e3smtime_i)):
            cloud_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(cloud[i,:,x_idx]))
            T_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(T[i,:,x_idx]))
            RH_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(RH[i,:,x_idx]))
            theta_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(theta[i,:,x_idx]))
            Q_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(Q[i,:,x_idx]))
            p_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(levm[i,:]))
            lwc_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(lwc[i,:,x_idx]))
            iwc_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(iwc[i,:,x_idx]))
            
        # lower tropospheric stability (theta diff between sfc and 700hPa)
        LTS700_2 = theta_p2[:, idx700] - theta_s2   
        LTS850_2 = theta_p2[:, idx850] - theta_s2  
           
        # mid-level T, P, z
        Tmid = 0.5*(T[:,0:-1,x_idx].data + T[:,1:,x_idx].data)
        Pmid = 0.5*(Pres[:,0:-1,x_idx].data + Pres[:,1:,x_idx].data)
        # cloud base
        cf_sum_base = np.cumsum(cloud[:, ::-1,x_idx].data, axis=1)
        cf_sum_base[cf_sum_base > 1] = 1
        cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
        p_cldbase = np.sum(Pmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
        T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
        # normalize cloud fraction when not 100%
        p_cldbase = p_cldbase / cf_sum_base[:,-1]
        T_cldbase = T_cldbase / cf_sum_base[:,-1]
        theta_cb_2 = T_cldbase * (100000./p_cldbase)**0.286 
        
        # combine data
        cloud_p = np.vstack((cloud_p,cloud_p2))
        T_p = np.vstack((T_p,T_p2))
        Q_p = np.vstack((Q_p,Q_p2))
        RH_p = np.vstack((RH_p,RH_p2))
        theta_p = np.vstack((theta_p,theta_p2))
        z_p = np.vstack((z_p,z_p2))
        lwc_p = np.vstack((lwc_p,lwc_p2))
        iwc_p = np.vstack((iwc_p,iwc_p2))
        cloud_z = np.vstack((cloud_z,cloud_z2))
        T_z = np.vstack((T_z,T_z2))
        RH_z = np.vstack((RH_z,RH_z2))
        theta_z = np.vstack((theta_z,theta_z2))
        Q_z = np.vstack((Q_z,Q_z2))
        p_z = np.vstack((p_z,p_z2))
        lwc_z = np.vstack((lwc_z,lwc_z2))
        iwc_z = np.vstack((iwc_z,iwc_z2))
        LTS700 = np.hstack((LTS700,LTS700_2))
        LTS850 = np.hstack((LTS850,LTS850_2))
        thetadiff_cb = np.hstack((thetadiff_cb, theta_cb_2-theta_s2))
        
    #%% re-shape the data into pre-defined resolution
    
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
    f = interp1d(np.int64(e3smtime), lwc_p, axis=0, bounds_error=False)
    lwc_p_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), iwc_p, axis=0, bounds_error=False)
    iwc_p_new = f(np.int64(time_new))
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
    f = interp1d(np.int64(e3smtime), lwc_z, axis=0, bounds_error=False)
    lwc_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), iwc_z, axis=0, bounds_error=False)
    iwc_z_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), LTS700, axis=0, bounds_error=False)
    LTS700_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), LTS850, axis=0, bounds_error=False)
    LTS850_new = f(np.int64(time_new))
    f = interp1d(np.int64(e3smtime), thetadiff_cb, axis=0, bounds_error=False)
    thetadiff_cb_new = f(np.int64(time_new))
        
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
    lp = xr.DataArray(data=lwc_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=lwc.long_name, units=lwc.units),)
    ip = xr.DataArray(data=iwc_p_new,  dims=["time","lev"],
        coords=dict(lev=(["lev"], lev_out), time=(["time"], time_new), ),
        attrs=dict(long_name=iwc.long_name, units=iwc.units),)
    varnames_p = [ 'cloud_p', 'T_p', 'Q_p', 'RH_p', 'theta_p', 'Z_p', 'RAINQM_p', 'SNOWQM_p']
    variables_p = [   cp,      tp,     qp,    rhp,    thp,    zp,       lp,     ip]
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
    lz = xr.DataArray(data=lwc_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name=lwc.long_name, units=lwc.units),)
    iz = xr.DataArray(data=iwc_z_new,  dims=["time","height"],
        coords=dict(height=(["height"], height_out), time=(["time"], time_new), ),
        attrs=dict(long_name=iwc.long_name, units=iwc.units),)
    varnames_z = [ 'cloud_z', 'T_z', 'Q_z', 'RH_z', 'theta_z', 'P_z', 'RAINQM_z', 'SNOWQM_z']
    variables_z = [   cz,      tz,     qz,    rhz,    thz,    pz,       lz,     iz]
    #
    l700 = xr.DataArray(data=LTS700_new,  dims=["time"],
        coords=dict(time=(["time"], time_new), ),
        attrs=dict(long_name='lower troposphere stability (700hPa theta - surface theta)', units='K'),)
    l850 = xr.DataArray(data=LTS850_new,  dims=["time"],
        coords=dict(time=(["time"], time_new), ),
        attrs=dict(long_name='lower troposphere stability (850hPa theta - surface theta)', units='K'),)
    thetadiff_cb = xr.DataArray(data=thetadiff_cb_new,  dims=["time"],
        coords=dict(time=(["time"], time_new), ),
        attrs=dict(long_name='Theta difference between cloud base and surface', units='K'),)
    varnames_1d = [ 'LTS700', 'LTS850', 'thetadiff_cb']
    variables_1d = [l700, l850, thetadiff_cb]
    
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
    print('output file '+outfile)
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

     