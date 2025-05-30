"""
prepare E3SM output data 
for surface (include remote sensing and satellite) and aircraft measurements in HISCALE campaign
input data is regional hourly output data in E3SM
all variables has a region appendix similar to "lat_260e_to_265e_34n_to_39n"
"""

import glob
import os
import re
import fnmatch
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import gamma
import xarray as xr
import pandas as pd
import time as ttt
import esmac_diags
from esmac_diags.subroutines.read_aircraft import read_iwg1
from esmac_diags.subroutines.time_format_change import hhmmss2sec
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d
from esmac_diags.subroutines.specific_data_treatment import find_nearest, calc_Reff_from_REL, \
                            calc_cdnc_ARM, calc_cdnc_VISST
from esmac_diags.subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm
from netCDF4 import Dataset

#%% test settings
# # input_path = '/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_hourly/'
# # output_path = '../../../data/HISCALE/model/'
# input_path = '../../../data/HISCALE/model/E3SMv2_out/'
# output_path = 'C:/Users/tang357/Downloads/HISCALE/'
# input_filehead = 'E3SMv1_HISCALE_test'
# output_filehead = 'E3SMv1_HISCALE'

# lev_out=np.arange(25.,1001,25.)
# height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
#                     1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
#                     3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
#                     10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

# dt = 3600
    
# # dt = 60
# # iwgpath = '../../../data/HISCALE/obs/aircraft/mei-iwg1/'

# import warnings
# warnings.filterwarnings("ignore")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, 
                      iwgpath, dt=60):
    """
    prepare E3SM output along flight tracks
    choose the nearest grid and level of the aircraft location
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
    iwgpath : str
        data path of aircraft location
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
        
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    #%% settings specific for each site
    # HISCALE
    E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    
    #%% find all data
    lst0 = glob.glob(iwgpath + '*.a2.txt')
    lst0.sort()
    
    for filename in lst0[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
        
        # re-shape the data into coarser resolution for output
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        lon_new = median_time_1d(time, lon, time_new)
        lat_new = median_time_1d(time, lat, time_new)
        height_new = median_time_1d(time, height, time_new)
        lon_new=lon_new+360   # make longitude consistent with E3SM from 0 to 360
    
        #%% read in E3SM data
        variable3d_names = ['T', 'Q', 'U', 'V', 'Z3', 'REI', 'REL', 'CCN1', 'CCN3', 'CCN4', 'CCN5', 
                            'CLDICE', 'CLDLIQ', 'CLOUD',  'FREQL',
                            'IWC', 'LWC', 'ICWNC', 'ICINC', ]
        variables = list()
        variables_new = list()
        for varname in variable3d_names:
            variables_new.append([])
        NCNall = np.empty((3000,0))
        p = list()      # pressure
        bc_all  = list()
        dst_all = list()
        mom_all = list()
        ncl_all = list()
        pom_all = list()
        so4_all = list()
        soa_all = list()
        phi_all = np.empty((999,0))
        
        lst = glob.glob(input_path + input_filehead+'.*'+timestr[0]+'-00000.nc')
        if len(lst)!=1:
            raise ValueError('Should only contain one file: '+lst)
        e3smdata = xr.open_dataset(lst[0])
        e3smtime = e3smdata.indexes['time'].to_datetimeindex()
        lonm = e3smdata['lon'+'_'+E3SMdomain_range].load()
        latm = e3smdata['lat'+'_'+E3SMdomain_range].load()
        z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
        
        P0 = e3smdata['P0'].load()
        hyam = e3smdata['hyam'].load()
        hybm = e3smdata['hybm'].load()
        T = e3smdata['T'+'_'+E3SMdomain_range].load()
        PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
        Pres = np.nan*T
        zlen = T.shape[1]
        for kk in range(zlen):
            Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
        
        # change time format into seconds of the day
        timem = np.float64((e3smtime - e3smtime[0]).seconds)
        
        # Get all simulated variables
        vlist = list(e3smdata.variables.keys())
        av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
        
        # variables for calculating aerosol size
        req_vlist = ['num_a1', 'num_a2', 'num_a3', 'num_a4', 'dgnd_a01', 'dgnd_a02', \
                     'dgnd_a03', 'dgnd_a04']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing for aerosol size')
            num_a1 = e3smdata['num_a1'+'_'+E3SMdomain_range].load()
            num_a2 = e3smdata['num_a2'+'_'+E3SMdomain_range].load()
            num_a3 = e3smdata['num_a3'+'_'+E3SMdomain_range].load()
            num_a4 = e3smdata['num_a4'+'_'+E3SMdomain_range].load()
            dn1 = e3smdata['dgnd_a01'+'_'+E3SMdomain_range].load()
            dn2 = e3smdata['dgnd_a02'+'_'+E3SMdomain_range].load()
            dn3 = e3smdata['dgnd_a03'+'_'+E3SMdomain_range].load()
            dn4 = e3smdata['dgnd_a04'+'_'+E3SMdomain_range].load()
        else:
            num_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            num_a2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            num_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            num_a4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            dn1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            dn2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            dn3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            dn4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            
        # aerosol composition
        req_vlist = ['bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', \
                     'mom_a2', 'mom_a3', 'mom_a4', 'ncl_a1', 'ncl_a2', 'ncl_a3', \
                     'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3', \
                     'soa_a1', 'soa_a2', 'soa_a3']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
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
        else:
            bc_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            bc_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            bc_a4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            dst_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            dst_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            mom_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            mom_a2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            mom_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            mom_a4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            ncl_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            ncl_a2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            ncl_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            pom_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            pom_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            pom_a4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            so4_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            so4_a2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            so4_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            soa_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            soa_a2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            soa_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        
        # droplet size distribution
        req_vlist = ['ICWNC', 'lambda_cloud', 'mu_cloud']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        if len(matched_vlist) == len(req_vlist):
            nd_cld = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
            lmda = e3smdata['lambda_cloud'+'_'+E3SMdomain_range].load()
            mu = e3smdata['mu_cloud'+'_'+E3SMdomain_range].load()
        else:
            nd_cld = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            lmda = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            mu = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        
        
        # other variables
        for varname in variable3d_names:
            try:
                var = e3smdata[varname + '_'+E3SMdomain_range].load()
            except:
                var = xr.DataArray(np.zeros(z3.shape)*np.nan,name=varname,\
                                   dims=["time","lev","ncol_"+E3SMdomain_range],coords={"time":e3smtime,"lev":e3smdata['lev'],"ncol_"+E3SMdomain_range:e3smdata['ncol_'+E3SMdomain_range]},\
                                   attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            variables.append(var)
        e3smdata.close()
        
        #%% find the flight track grid
        for tt in range(len(time_new)):
            t_idx = np.abs(timem-time_new[tt]).argmin()
            x_idx = find_nearest(lonm, latm, lon_new[tt], lat_new[tt])
            z_idx = np.abs(z3[t_idx, :, x_idx]-height_new[tt]).argmin()
            for vv in range(len(variable3d_names)):
                variables_new[vv].append(float(variables[vv][t_idx, z_idx, x_idx]))
            p.append(Pres[t_idx, z_idx, x_idx].data)
            # calculate aerosol size
            numall = [num_a1[t_idx, z_idx, x_idx].data, num_a2[t_idx, z_idx, x_idx].data, 
                      num_a3[t_idx, z_idx, x_idx].data, num_a4[t_idx, z_idx, x_idx].data]
            dnall  = [dn1[t_idx, z_idx, x_idx].data,    dn2[t_idx, z_idx, x_idx].data,    
                      dn3[t_idx, z_idx, x_idx].data,    dn4[t_idx, z_idx, x_idx].data]
            NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[t_idx, z_idx, x_idx].data, Pres[t_idx, z_idx, x_idx].data)
            NCNall = np.hstack((NCNall, np.reshape(NCN,(3000,1))))
            # calculate aerosol composition
            bc_all.append(bc_a1[t_idx, z_idx, x_idx].data +                       
                    bc_a3[t_idx, z_idx, x_idx].data + bc_a4[t_idx, z_idx, x_idx].data)
            dst_all.append(dst_a1[t_idx, z_idx, x_idx].data +                      
                    dst_a3[t_idx, z_idx, x_idx].data)
            mom_all.append(mom_a1[t_idx, z_idx, x_idx].data + mom_a2[t_idx, z_idx, x_idx].data + 
                    mom_a3[t_idx, z_idx, x_idx].data + mom_a4[t_idx, z_idx, x_idx].data)
            ncl_all.append(ncl_a1[t_idx, z_idx, x_idx].data + ncl_a2[t_idx, z_idx, x_idx].data + 
                    ncl_a3[t_idx, z_idx, x_idx].data)
            pom_all.append(pom_a1[t_idx, z_idx, x_idx].data +                    
                    pom_a3[t_idx, z_idx, x_idx].data + pom_a4[t_idx, z_idx, x_idx].data)
            so4_all.append(so4_a1[t_idx, z_idx, x_idx].data + so4_a2[t_idx, z_idx, x_idx].data + 
                    so4_a3[t_idx, z_idx, x_idx].data)
            soa_all.append(soa_a1[t_idx, z_idx, x_idx].data + soa_a2[t_idx, z_idx, x_idx].data + 
                    soa_a3[t_idx, z_idx, x_idx].data)
            # calculate droplet size distribution
            N0 = nd_cld[t_idx, z_idx, x_idx].data * (lmda[t_idx, z_idx, x_idx].data ** (mu[t_idx, z_idx, x_idx].data+1)) / \
                    gamma(mu[t_idx, z_idx, x_idx].data+1)    # parameter N0
            D_cld = np.arange(1, 1000) * 1e-6  # in m
            phi = N0 * (D_cld**mu[t_idx, z_idx, x_idx].data) * np.exp(- lmda[t_idx, z_idx, x_idx].data * D_cld)
            phi_all = np.hstack((phi_all, np.reshape(phi,(len(D_cld),1))))
            nd_bin = phi_all * (D_cld[1] - D_cld[0])   # droplet number concentration in each size bin
            
        NCN3 = np.nansum(NCNall[3:, :], 0)   # >3nm
        NCN10 = np.nansum(NCNall[10:, :], 0)    # >10nm
        NCN100 = np.nansum(NCNall[100:, :], 0)    # >100nm
        
        # #%%
        # nd_bin[nd_bin<1e-6]=np.nan
        
        # import matplotlib.pyplot as plt
        # fig,(ax2) = plt.subplots(1,1,figsize=(8,4))
        # ax2.plot(np.nanmean(nd_bin*1e-6,axis=1))
        # ax2.set_xscale('log')
        # ax2.set_yscale('log')
        # ax2.set_ylim(1e-8, 1e-2)
        
        # fig,(ax2) = plt.subplots(1,1,figsize=(8,4))
        # ax2.contourf(time_new, np.arange(1,1000), nd_bin*1e-6)
        # # ax2.set_xscale('log')
        # ax2.plot(time_new, height_new)
        # ax2.set_yscale('log')
        # e
        #%% change some units
        # composition
        T = variables_new[variable3d_names.index('T')]
        rho = np.array(p)/T/287.06
        bc_all = np.array(bc_all)*1e9*rho
        dst_all = np.array(dst_all)*1e9*rho
        mom_all = np.array(mom_all)*1e9*rho
        ncl_all = np.array(ncl_all)*1e9*rho
        pom_all = np.array(pom_all)*1e9*rho
        so4_all = np.array(so4_all)*1e9*rho
        soa_all = np.array(soa_all)*1e9*rho
        composition_units = 'ug/m3'
        # aerosol number
        NCNall = NCNall * 1e-6
        NCN3 = NCN3 * 1e-6
        NCN10 = NCN10 * 1e-6
        NCN100 = NCN100 * 1e-6
        ncn_units = '#/cm3'
        # cloud number size distribution
        nd_bin = nd_bin * 1e-6
        nd_units = '#/cm3'
        # LWC and IWC
        idx = variable3d_names.index('LWC')
        variables_new[idx] = np.array(variables_new[idx])*1000
        variables[idx].attrs['units']='g/m3'
        idx = variable3d_names.index('IWC')
        variables_new[idx] = np.array(variables_new[idx])*1000
        variables[idx].attrs['units']='g/m3'
        # droplet number
        idx = variable3d_names.index('ICWNC')
        variables_new[idx] = np.array(variables_new[idx])*1e-6
        variables[idx].attrs['units']='#/cm3'
        idx = variable3d_names.index('ICINC')
        variables_new[idx] = np.array(variables_new[idx])*1e-6
        variables[idx].attrs['units']='#/cm3'
    
        #%% output 
        
        outfile = output_path + output_filehead + '_flight_'+date+'.nc'
        print('output file '+outfile)
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        t = f.createDimension('time', None)  # unlimited
        s = f.createDimension('CNsize', 3000)  # unlimited
        s = f.createDimension('Ndsize', 999)  # unlimited
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time",))
        height_o = f.createVariable("height", 'f8', ("time",))
        var_o = list()
        for vv in range(len(variable3d_names)):
            var_o.append (f.createVariable(variable3d_names[vv], 'f8', ("time", )))
        p_o = f.createVariable('pres', 'f8', ("time",))
        bc_o = f.createVariable('bc', 'f8', ("time",))
        dst_o = f.createVariable('dst', 'f8', ("time",))
        pom_o = f.createVariable('pom', 'f8', ("time",))
        mom_o = f.createVariable('mom', 'f8', ("time",))
        ncl_o = f.createVariable('ncl', 'f8', ("time",))
        so4_o = f.createVariable('so4', 'f8', ("time",))
        soa_o = f.createVariable('soa', 'f8', ("time",))
        ncn_o = f.createVariable('NCNall', 'f8', ("CNsize", "time",))
        ncn3_o = f.createVariable('NCN3', 'f8', ("time",))
        ncn10_o = f.createVariable('NCN10', 'f8', ("time",))
        ncn100_o = f.createVariable('NCN100', 'f8', ("time",))
        nd_o = f.createVariable('Nd_bin', 'f8', ("Ndsize", "time",))
        
        # write data
        time_o[:] = time_new
        height_o[:] = height_new
        for vv in range(len(variable3d_names)):
            var_o[vv][:] = np.array(variables_new[vv])
        p_o[:] = np.array(p)
        bc_o[:] = np.array(bc_all)
        dst_o[:] = np.array(dst_all)
        pom_o[:] = np.array(pom_all)
        mom_o[:] = np.array(mom_all)
        ncl_o[:] = np.array(ncl_all)
        so4_o[:] = np.array(so4_all)
        soa_o[:] = np.array(soa_all)
        ncn_o[:,:] = NCNall
        ncn3_o[:] = NCN3
        ncn10_o[:] = NCN10
        ncn100_o[:] = NCN100
        nd_o[:,:] = nd_bin
        
        # attributes
        time_o.units = "Seconds since " + timestr[0] + ' 00:00 UTC'
        height_o.units = 'm MSL'
        for vv in range(len(variable3d_names)):
            var_o[vv].units = variables[vv].units
            var_o[vv].long_name = variables[vv].long_name
        p_o.units = 'Pa'
        p_o.long_name = 'Pressure'
        bc_o.units = composition_units
        bc_o.long_name = 'total black carbon aerosol concentration'
        dst_o.units = composition_units
        dst_o.long_name = 'total dust aerosol concentration'
        ncl_o.units = composition_units
        ncl_o.long_name = 'total sea salt aerosol concentration'
        pom_o.units = composition_units
        pom_o.long_name = 'total primary organic aerosol concentration'
        mom_o.units = composition_units
        mom_o.long_name = 'total marine organic aerosol concentration'
        so4_o.units = composition_units
        so4_o.long_name = 'total sulfate aerosol concentration'
        soa_o.units = composition_units
        soa_o.long_name = 'total secondary organic aerosol concentration'
        ncn_o.units = ncn_units
        ncn_o.long_name = 'Aerosol number size distribution'
        ncn_o.description = 'calculated from modal information into 1nm increment'
        ncn3_o.units = ncn_units
        ncn3_o.long_name = 'Aerosol number concentration for size >3nm'
        ncn10_o.units = ncn_units
        ncn10_o.long_name = 'Aerosol number concentration for size >10nm'
        ncn100_o.units = ncn_units
        ncn100_o.long_name = 'Aerosol number concentration for size >100nm'
        nd_o.units = nd_units
        nd_o.long_name = 'cloud droplet number size distribution'
        nd_o.description = 'calculated from MG2 cloud spectrum output into 1um increment from 1um to 1000um'
        
        # global attributes
        f.title = 'preprocessed E3SM data along aircraft track at the nearest time, grid, and vertical level'
        f.aircraftfile = filename.split('/')[-1]
        f.modelfile = lst[0].split('/')[-1]
        f.date = ttt.ctime(ttt.time())
        
        f.close()
        
    
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
    # HISCALE
    lat0 = 36.6059
    lon0 = -97.48792
    E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    
    # output time range and resolution
    time_new = pd.date_range(start='2016-04-25', end='2016-09-23', freq=str(int(dt))+"s")  # HISCALE time period
    
    #%% read in data
    
    lst = glob.glob(input_path + input_filehead+'.*.h?.*2016-04-2?-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-05-??-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-06-0?-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-08-2?-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-09-??-00000.nc') 
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
    e3smdata.close()
    # only extract the model column at the site
    if lon0<0:
        lon0=lon0+360   # make longitude consistent with E3SM from 0 to 360
    x_idx = find_nearest(lonm,latm,lon0,lat0)
    
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
    for i in range(len(e3smtime)):
        cloud_p[i,:] = np.interp(lev_out,levm[i,:],cloud[i,:,x_idx])
        T_p[i,:] = np.interp(lev_out,levm[i,:],T[i,:,x_idx])
        Q_p[i,:] = np.interp(lev_out,levm[i,:],Q[i,:,x_idx])
        RH_p[i,:] = np.interp(lev_out,levm[i,:],RH[i,:,x_idx])
        theta_p[i,:] = np.interp(lev_out,levm[i,:],theta[i,:,x_idx])
        z_p[i,:] = np.interp(lev_out,levm[i,:],z3[i,:,x_idx])
            
    # interpolate data into height coordinate. flip model data since numpy.interp only works for increasing dimension
    cloud_z = np.empty((len(e3smtime),len(height_out)))
    T_z = np.empty((len(e3smtime),len(height_out)))
    RH_z = np.empty((len(e3smtime),len(height_out)))
    theta_z = np.empty((len(e3smtime),len(height_out)))
    Q_z = np.empty((len(e3smtime),len(height_out)))
    p_z = np.empty((len(e3smtime),len(height_out)))
    for i in range(len(e3smtime)):
        cloud_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(cloud[i,:,x_idx]))
        T_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(T[i,:,x_idx]))
        RH_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(RH[i,:,x_idx]))
        theta_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(theta[i,:,x_idx]))
        Q_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(Q[i,:,x_idx]))
        p_z[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(levm[i,:]))
        
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
        e3smdata.close()
        
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
        for i in range(len(e3smtime_i)):
            cloud_p2[i,:] = np.interp(lev_out,levm[i,:],cloud[i,:,x_idx])
            T_p2[i,:] = np.interp(lev_out,levm[i,:],T[i,:,x_idx])
            Q_p2[i,:] = np.interp(lev_out,levm[i,:],Q[i,:,x_idx])
            RH_p2[i,:] = np.interp(lev_out,levm[i,:],RH[i,:,x_idx])
            theta_p2[i,:] = np.interp(lev_out,levm[i,:],theta[i,:,x_idx])
            z_p2[i,:] = np.interp(lev_out,levm[i,:],z3[i,:,x_idx])
        
        # interpolate data into height coordinate
        cloud_z2 = np.empty((len(e3smtime_i),len(height_out)))
        T_z2 = np.empty((len(e3smtime_i),len(height_out)))
        RH_z2 = np.empty((len(e3smtime_i),len(height_out)))
        theta_z2 = np.empty((len(e3smtime_i),len(height_out)))
        Q_z2 = np.empty((len(e3smtime_i),len(height_out)))
        p_z2 = np.empty((len(e3smtime_i),len(height_out)))
        for i in range(len(e3smtime_i)):
            cloud_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(cloud[i,:,x_idx]))
            T_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(T[i,:,x_idx]))
            RH_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(RH[i,:,x_idx]))
            theta_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(theta[i,:,x_idx]))
            Q_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(Q[i,:,x_idx]))
            p_z2[i,:] = np.interp(height_out,np.flip(z3[i,:,x_idx]),np.flip(levm[i,:]))
            
        # lower tropospheric stability (theta diff between sfc and 700hPa)
        LTS700_2 = theta_p2[:, idx700] - theta_s2   
        LTS850_2 = theta_p2[:, idx850] - theta_s2  
           
        # combine data
        cloud_p = np.vstack((cloud_p,cloud_p2))
        T_p = np.vstack((T_p,T_p2))
        Q_p = np.vstack((Q_p,Q_p2))
        RH_p = np.vstack((RH_p,RH_p2))
        theta_p = np.vstack((theta_p,theta_p2))
        z_p = np.vstack((z_p,z_p2))
        cloud_z = np.vstack((cloud_z,cloud_z2))
        T_z = np.vstack((T_z,T_z2))
        RH_z = np.vstack((RH_z,RH_z2))
        theta_z = np.vstack((theta_z,theta_z2))
        Q_z = np.vstack((Q_z,Q_z2))
        p_z = np.vstack((p_z,p_z2))
        LTS700 = np.hstack((LTS700,LTS700_2))
        LTS850 = np.hstack((LTS850,LTS850_2))
        
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

     
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, dt=3600):
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
    # HISCALE
    lat0 = 36.6059
    lon0 = -97.48792
    E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    
    # output time range and resolution
    time_new = pd.date_range(start='2016-04-25', end='2016-09-23', freq=str(int(dt))+"s")  # HISCALE time period
    print(len(time_new))
    
    #%% read in data
    variable_names = list()
    variables = list()
    
    lst = glob.glob(input_path + input_filehead+'.*.h?.*2016-04-2?-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-05-??-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-06-0?-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-08-2?-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*.h?.*2016-09-??-00000.nc') 
    lst.sort()
    # first data
    e3smdata = xr.open_dataset(lst[0])
    e3smtime = e3smdata.indexes['time'].to_datetimeindex()
    lonm = e3smdata['lon'+'_'+E3SMdomain_range].load()
    latm = e3smdata['lat'+'_'+E3SMdomain_range].load()
    P0 = e3smdata['P0'].load()
    hyam = e3smdata['hyam'].load()
    hybm = e3smdata['hybm'].load()
    T = e3smdata['T'+'_'+E3SMdomain_range].load()
    PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
    # getting column and levels
    len_ncol = len(e3smdata['ncol_'+E3SMdomain_range])
    len_lev = len(e3smdata['lev'])
    # only extract the model column at the site
    if lon0<0:
        lon0=lon0+360   # make longitude consistent with E3SM from 0 to 360
    x_idx = find_nearest(lonm,latm,lon0,lat0)
    
    # Get all simulated variables
    vlist = list(e3smdata.variables.keys())
    av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
    
    # aerosol composition
    req_vlist = ['bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', \
                 'mom_a2', 'mom_a3', 'mom_a4', 'ncl_a1', 'ncl_a2', 'ncl_a3', \
                 'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3', \
                 'soa_a1', 'soa_a2', 'soa_a3']
    req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing for aerosol composition')
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
        bc_all  = bc_a1[:,-1,x_idx] +                       bc_a3[:,-1,x_idx] + bc_a4[:,-1,x_idx]
        dst_all = dst_a1[:,-1,x_idx] +                      dst_a3[:,-1,x_idx]
        mom_all = mom_a1[:,-1,x_idx] + mom_a2[:,-1,x_idx] + mom_a3[:,-1,x_idx] + mom_a4[:,-1,x_idx]
        ncl_all = ncl_a1[:,-1,x_idx] + ncl_a2[:,-1,x_idx] + ncl_a3[:,-1,x_idx]
        pom_all = pom_a1[:,-1,x_idx] +                    + pom_a3[:,-1,x_idx] + pom_a4[:,-1,x_idx]
        so4_all = so4_a1[:,-1,x_idx] + so4_a2[:,-1,x_idx] + so4_a3[:,-1,x_idx]
        soa_all = soa_a1[:,-1,x_idx] + soa_a2[:,-1,x_idx] + soa_a3[:,-1,x_idx]
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
    else:
        bc_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        dst_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        mom_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        ncl_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        pom_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        so4_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        soa_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    # aerosol size
    req_vlist = ['num_a1', 'num_a2', 'num_a3', 'num_a4', 'dgnd_a01', 'dgnd_a02', \
                 'dgnd_a03', 'dgnd_a04']
    req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing for aerosol size')
        num_a1 = e3smdata['num_a1'+'_'+E3SMdomain_range].load()
        num_a2 = e3smdata['num_a2'+'_'+E3SMdomain_range].load()
        num_a3 = e3smdata['num_a3'+'_'+E3SMdomain_range].load()
        num_a4 = e3smdata['num_a4'+'_'+E3SMdomain_range].load()
        dn1 = e3smdata['dgnd_a01'+'_'+E3SMdomain_range].load()
        dn2 = e3smdata['dgnd_a02'+'_'+E3SMdomain_range].load()
        dn3 = e3smdata['dgnd_a03'+'_'+E3SMdomain_range].load()
        dn4 = e3smdata['dgnd_a04'+'_'+E3SMdomain_range].load()
        Pres = np.nan*T
        zlen = T.shape[1]
        for kk in range(zlen):
            Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
        numall = [num_a1[:, -1, x_idx].data, num_a2[:, -1, x_idx].data, num_a3[:, -1, x_idx].data, num_a4[:, -1, x_idx].data]
        dnall  = [dn1[:, -1, x_idx].data,    dn2[:, -1, x_idx].data,    dn3[:, -1, x_idx].data,    dn4[:, -1, x_idx].data]
        NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1, x_idx].data, Pres[:, -1, x_idx].data)
        NCNall = xr.DataArray(data=NCN,  dims=["size", "time"],
            coords=dict(time=(["time"], e3smtime),size=(["size"], np.arange(1,3001))),
            attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
    else:
        NCNall = xr.DataArray(np.zeros((3000,len(e3smtime)))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    # variables to calculate Reff and Nd
    req_vlist = ['Z3', 'CLOUD']
    req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing for variables to calculate Reff and Nd')
        z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
        cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
        z3 = z3[:,:,x_idx]
        cloud = cloud[:,:,x_idx]
        dz = (z3[:,:-2].data - z3[:,2:].data)/2
        dz = np.append(dz, (z3[:,-2:-1].data+z3[:,-1:].data)/2, axis=1)
        dz = np.insert(dz,0,dz[:,0],axis=1)
        weight = cloud.data*dz
        # mid-level T, P, z
        Tmid = 0.5*(T[:,0:-1,x_idx].data + T[:,1:,x_idx].data)
        Pmid = 0.5*(Pres[:,0:-1,x_idx].data + Pres[:,1:,x_idx].data)
        Zmid = 0.5*(z3[:,0:-1].data + z3[:,1:].data)
        # cloud top
        cf_sum_top = np.cumsum(cloud.data, axis=1)
        cf_sum_top[cf_sum_top > 1] = 1
        cf_sum_top_diff = cf_sum_top[:,1:] - cf_sum_top[:,0:-1]
        z_cldtop = np.sum(Zmid[:,:]*cf_sum_top_diff[:,:], axis=1)
        T_cldtop = np.sum(Tmid[:,:]*cf_sum_top_diff[:,:], axis=1)
        # cloud base
        cf_sum_base = np.cumsum(cloud[:, ::-1].data, axis=1)
        cf_sum_base[cf_sum_base > 1] = 1
        cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
        z_cldbase = np.sum(Zmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
        T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
        # normalize cloud fraction when not 100%
        z_cldtop = z_cldtop / cf_sum_top[:,-1]
        T_cldtop = T_cldtop / cf_sum_top[:,-1]
        z_cldbase = z_cldbase / cf_sum_base[:,-1]
        T_cldbase = T_cldbase / cf_sum_base[:,-1]
        e3sm_cloud_depth = z_cldtop - z_cldbase
        # # alternatively, one can assume a minimum cloud depth in which the full mass levels are used instead of half levels
        # e3sm_z2 = z3[:,:-1].data
        # e3sm_z3 = z3[:,1:].data
        # e3sm_t2 = T[:,:-1,x_idx].data
        # e3sm_t3 = T[:,1:,x_idx].data
        # #cloud base and top calculations for minimum cloud depth
        # z_cldbase2 = np.sum(e3sm_z2[:,::-1]*cf_sum_base_diff[:,:], axis=1)
        # T_cldbase2 = np.sum(e3sm_t2[:,::-1]*cf_sum_base_diff[:,:], axis=1)
        # z_cldtop2 = np.sum(e3sm_z3[:,:]*cf_sum_top_diff[:,:], axis=1)
        # T_cldtop2 = np.sum(e3sm_t3[:,:]*cf_sum_top_diff[:,:], axis=1)
        # e3sm_cloud_depth2 = z_cldtop2 - z_cldbase2
        cloud_depth = xr.DataArray(data=e3sm_cloud_depth,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="cloud depth",units="m"),)
        cbh = xr.DataArray(data=z_cldbase,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="cloud base height",units="m"),)
        cth = xr.DataArray(data=z_cldtop,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="cloud top height",units="m"),)
        cbt = xr.DataArray(data=T_cldbase,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="cloud base temperature",units="K"),)
        ctt = xr.DataArray(data=T_cldtop,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="cloud top temperature",units="K"),)
    else:
        cloud_depth = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cbh = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cth = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cbt = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        ctt = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    # cloud optical depth and effective radius
    req_vlist = ['REL', 'FREQL', 'ICWNC', 'TOT_CLD_VISTAU', 'TAUWMODIS', 'SOLIN']
    req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing for cloud optical depth and effective radius')
        rel = e3smdata['REL'+'_'+E3SMdomain_range].load()
        freql = e3smdata['FREQL'+'_'+E3SMdomain_range].load()
        icwnc = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
        cod_a = e3smdata['TOT_CLD_VISTAU'+'_'+E3SMdomain_range].load()
        cod_m = e3smdata['TAUWMODIS'+'_'+E3SMdomain_range].load()*0.01   # cloud fraction is treated as 1 but is 100
        solin = e3smdata['SOLIN'+'_'+E3SMdomain_range].load()
        rel = rel[:,:,x_idx]
        freql = freql[:,:,x_idx]
        icwnc = icwnc[:,:,x_idx]
        cod_a = cod_a[:,:,x_idx]
        solin = solin[:,x_idx]
        # calculate mean effective radius. 
        reff = calc_Reff_from_REL(rel.data, dz, freql.data, icwnc.data)
        reff[reff==0] = np.nan
        # calculate mean optical depth
        cod = np.sum(cod_a.data,axis=1)
        # cod[solin==0] = np.nan
        # cod from MODIS simulator
        cod_m = cod_m[:,x_idx]
        reff_mean = xr.DataArray(data=reff,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
        cod_mean = xr.DataArray(data=cod,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
    else:
        cod_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cod',attrs={'units':'dummy_unit','long_name':'Dummy'})
        reff_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='reff',attrs={'units':'dummy_unit','long_name':'Dummy'})
        cod_m = xr.DataArray(np.zeros(len(e3smtime))*np.nan)
    
    # mean cloud droplet number concentration
    req_vlist = ['CDNUMC']
    req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing for mean cloud droplet number concentration')
        cdnc_col = e3smdata['CDNUMC'+'_'+E3SMdomain_range].load()
        cdnc_col = cdnc_col[:,x_idx]
        cdnc_mean = cdnc_col/np.sum(weight,axis=1)
        cdnc_mean[cdnc_mean >2e9] = np.nan
        cdnc_mean = xr.DataArray(data=cdnc_mean,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
    else:
        cdnc_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
    req_vlist = ['TGCLDLWP']
    req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
    matched_vlist = list(set(av_vars).intersection(req_vlist))
    
    if len(matched_vlist) == len(req_vlist):
        print('\nAnalyzing for cloud droplet number concentration retrieved like Ndrop and Bennartz 2007')
        lwp = e3smdata['TGCLDLWP'+'_'+E3SMdomain_range][:,x_idx].data
        e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
        T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
        nd_arm = calc_cdnc_ARM(lwp, cod_m, e3sm_cloud_depth)
        nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m, adiabaticity=0.8)
        cdnc_arm = xr.DataArray(data=nd_arm*1e6,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                        description='Retrieved using ARM Ndrop algorithm'),)
        cdnc_sat = xr.DataArray(data=nd_sat*1e6,  dims=["time"],
            coords=dict(time=(["time"], e3smtime)),
            attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                        description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)
    else:
        cdnc_arm = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        cdnc_sat = xr.DataArray(np.zeros(len(e3smtime))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
    
    # all other 2D (surface and vertical integrated) variables
    variable2d_names = ['AODABS', 'AODALL', 'CDNUMC', 'CLDHGH', 'CLDMED', 'CLDLOW', 'CLDTOT', 
                        'FLDS', 'FLNS', 'FLNT', 'FLUT', 'FSDS', 'FSNS', 'FSNT', 'FSUTOA', 'SOLIN', 
                        'LHFLX', 'SHFLX', 'LWCF', 'SWCF', "TGCLDLWP", "TGCLDIWP", 
                        'IWPMODIS', 'LWPMODIS', 'REFFCLWMODIS', 'TAUIMODIS', 'TAUTMODIS', 'TAUWMODIS', 
                        #'MEANPTOP_ISCCP', 'MEANCLDALB_ISCCP', 'MEANTAU_ISCCP',
                        'PBLH', 'PRECT', 'PRECL', 'PRECC', 'PS', 'TREFHT', ]
    
    for varname in variable2d_names:
        try:
            var = e3smdata[varname + '_'+E3SMdomain_range].load()
            var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
        except:
            var = xr.DataArray(np.zeros((len(e3smtime),len_ncol))*np.nan,name=varname,\
                               dims=["time","ncol_"+E3SMdomain_range],coords={"time":e3smtime,"ncol_"+E3SMdomain_range:np.arange(len_ncol)},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
        if varname=='AODABS' or varname=='AODALL':
            var.attrs['units']='N/A'
        variable_names.append(varname)
        variables.append(var[:,x_idx])
    
    # all other 3D (with vertical level) variables at the lowest model level
    variable3d_names = ['CCN1', 'CCN3', 'CCN4', 'CCN5', 'Q', 'T', 'RELHUM', 'U', 'V'] 
    for varname in variable3d_names:
        try:
            var = e3smdata[varname + '_'+E3SMdomain_range].load()
            var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
        except:
            var = xr.DataArray(np.zeros((len(e3smtime),len_lev,len_ncol))*np.nan,name=varname,\
                               dims=["time","lev","ncol_"+E3SMdomain_range],coords={"time":e3smtime,"lev":e3smdata['lev'],"ncol_"+E3SMdomain_range:e3smdata['ncol_'+E3SMdomain_range]},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
        variables.append(var[:,-1,x_idx])
        variable_names.append(varname)
    
    e3smdata.close()
    
    #%%  add data for each day
    for file in lst[1:]:
        print(file)
        e3smdata = xr.open_dataset(file)
        e3smtime_i = e3smdata.indexes['time'].to_datetimeindex()
        e3smtime = np.hstack((e3smtime, e3smtime_i))
        
        # aerosol composition
        req_vlist = ['bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', \
                     'mom_a2', 'mom_a3', 'mom_a4', 'ncl_a1', 'ncl_a2', 'ncl_a3', \
                     'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3', \
                     'soa_a1', 'soa_a2', 'soa_a3']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
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
            bc  = bc_a1[:,-1,x_idx] +                       bc_a3[:,-1,x_idx] + bc_a4[:,-1,x_idx]
            dst = dst_a1[:,-1,x_idx] +                      dst_a3[:,-1,x_idx]
            mom = mom_a1[:,-1,x_idx] + mom_a2[:,-1,x_idx] + mom_a3[:,-1,x_idx] + mom_a4[:,-1,x_idx]
            ncl = ncl_a1[:,-1,x_idx] + ncl_a2[:,-1,x_idx] + ncl_a3[:,-1,x_idx]
            pom = pom_a1[:,-1,x_idx] +                    + pom_a3[:,-1,x_idx] + pom_a4[:,-1,x_idx]
            so4 = so4_a1[:,-1,x_idx] + so4_a2[:,-1,x_idx] + so4_a3[:,-1,x_idx]
            soa = soa_a1[:,-1,x_idx] + soa_a2[:,-1,x_idx] + soa_a3[:,-1,x_idx]
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
        else:
            bc_all  = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='bc_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
            dst_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='dst_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
            mom_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='mom_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
            ncl_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='ncl_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
            pom_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='pom_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
            so4_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='so4_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
            soa_all = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='soa_all',attrs={'units':'dummy_unit','long_name':'Dummy'})
        
        # aerosol size
        req_vlist = ['num_a1', 'num_a2', 'num_a3', 'num_a4', 'dgnd_a01', 'dgnd_a02', \
                     'dgnd_a03', 'dgnd_a04']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            num_a1 = e3smdata['num_a1'+'_'+E3SMdomain_range].load()
            num_a2 = e3smdata['num_a2'+'_'+E3SMdomain_range].load()
            num_a3 = e3smdata['num_a3'+'_'+E3SMdomain_range].load()
            num_a4 = e3smdata['num_a4'+'_'+E3SMdomain_range].load()
            dn1 = e3smdata['dgnd_a01'+'_'+E3SMdomain_range].load()
            dn2 = e3smdata['dgnd_a02'+'_'+E3SMdomain_range].load()
            dn3 = e3smdata['dgnd_a03'+'_'+E3SMdomain_range].load()
            dn4 = e3smdata['dgnd_a04'+'_'+E3SMdomain_range].load()
            Pres = np.nan*T
            zlen = T.shape[1]
            for kk in range(zlen):
                Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
            numall = [num_a1[:, -1, x_idx].data, num_a2[:, -1, x_idx].data, num_a3[:, -1, x_idx].data, num_a4[:, -1, x_idx].data]
            dnall  = [dn1[:, -1, x_idx].data,    dn2[:, -1, x_idx].data,    dn3[:, -1, x_idx].data,    dn4[:, -1, x_idx].data]
            NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1, x_idx].data, Pres[:, -1, x_idx].data)
            NCN2 = xr.DataArray(data=NCN,  dims=["size", "time"],
                coords=dict(time=(["time"], e3smtime_i),size=(["size"], np.arange(1,3001))),
                attrs=dict(long_name="Aerosol number size distribution",units="#/m3"),)
            NCNall = xr.concat([NCNall, NCN2], dim="time")
        else:
            NCNall = xr.DataArray(np.zeros((3000,len(e3smtime_i)))*np.nan,name='NCNall',attrs={'units':'dummy_unit','long_name':'Dummy'})
    
        # variables to calculate cloud heights and depth
        req_vlist = ['Z3', 'CLOUD']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
            cloud = e3smdata['CLOUD'+'_'+E3SMdomain_range].load()
            z3 = z3[:,:,x_idx]
            cloud = cloud[:,:,x_idx]
            dz = (z3[:,:-2].data - z3[:,2:].data)/2
            dz = np.append(dz, (z3[:,-2:-1].data+z3[:,-1:].data)/2, axis=1)
            dz = np.insert(dz,0,dz[:,0],axis=1)
            weight = cloud.data*dz
            # mid-level T, P, z
            Tmid = 0.5*(T[:,0:-1,x_idx].data + T[:,1:,x_idx].data)
            Pmid = 0.5*(Pres[:,0:-1,x_idx].data + Pres[:,1:,x_idx].data)
            Zmid = 0.5*(z3[:,0:-1].data + z3[:,1:].data)
            # cloud top
            cf_sum_top = np.cumsum(cloud.data, axis=1)
            cf_sum_top[cf_sum_top > 1] = 1
            cf_sum_top_diff = cf_sum_top[:,1:] - cf_sum_top[:,0:-1]
            z_cldtop = np.sum(Zmid[:,:]*cf_sum_top_diff[:,:], axis=1)
            T_cldtop = np.sum(Tmid[:,:]*cf_sum_top_diff[:,:], axis=1)
            # cloud base
            cf_sum_base = np.cumsum(cloud[:, ::-1].data, axis=1)
            cf_sum_base[cf_sum_base > 1] = 1
            cf_sum_base_diff = cf_sum_base[:,1:] - cf_sum_base[:,0:-1]
            z_cldbase = np.sum(Zmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
            T_cldbase = np.sum(Tmid[:,::-1]*cf_sum_base_diff[:,:], axis=1)
            # normalize cloud fraction when not 100%
            z_cldtop = z_cldtop / cf_sum_top[:,-1]
            T_cldtop = T_cldtop / cf_sum_top[:,-1]
            z_cldbase = z_cldbase / cf_sum_base[:,-1]
            T_cldbase = T_cldbase / cf_sum_base[:,-1]
            e3sm_cloud_depth = z_cldtop - z_cldbase
            cloud_depth_2 = xr.DataArray(data=e3sm_cloud_depth,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="cloud depth",units="m"),)
            cbh_2 = xr.DataArray(data=z_cldbase,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="cloud base height",units="m"),)
            cth_2 = xr.DataArray(data=z_cldtop,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="cloud top height",units="m"),)
            cbt_2 = xr.DataArray(data=T_cldbase,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="cloud base temperature",units="K"),)
            ctt_2 = xr.DataArray(data=T_cldtop,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="cloud top temperature",units="K"),)
            cbh = xr.concat([cbh, cbh_2], dim="time")
            cth = xr.concat([cth, cth_2], dim="time")
            cbt = xr.concat([cbt, cbt_2], dim="time")
            ctt = xr.concat([ctt, ctt_2], dim="time")
            cloud_depth = xr.concat([cloud_depth, cloud_depth_2], dim="time")
        else:
            cloud_depth = xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cloud_depth',attrs={'units':'dummy_unit','long_name':'Dummy'})
            cbh = xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cbh',attrs={'units':'dummy_unit','long_name':'Dummy'})
            cth = xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cth',attrs={'units':'dummy_unit','long_name':'Dummy'})
            cbt = xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='cbt',attrs={'units':'dummy_unit','long_name':'Dummy'})
            ctt = xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,name='ctt',attrs={'units':'dummy_unit','long_name':'Dummy'})
        
    
        # cloud optical depth and effective radius
        req_vlist = ['REL', 'FREQL', 'ICWNC', 'TOT_CLD_VISTAU', 'TAUWMODIS', 'SOLIN']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            rel = e3smdata['REL'+'_'+E3SMdomain_range].load()
            freql = e3smdata['FREQL'+'_'+E3SMdomain_range].load()
            icwnc = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
            cod_a = e3smdata['TOT_CLD_VISTAU'+'_'+E3SMdomain_range].load()
            cod_m = e3smdata['TAUWMODIS'+'_'+E3SMdomain_range].load()*0.01   # cloud fraction is treated as 1 but is 100
            solin = e3smdata['SOLIN'+'_'+E3SMdomain_range].load()
            rel = rel[:,:,x_idx]
            freql = freql[:,:,x_idx]
            icwnc = icwnc[:,:,x_idx]
            cod_a = cod_a[:,:,x_idx]
            solin = solin[:,x_idx]
            # calculate mean effective radius. 
            reff = calc_Reff_from_REL(rel.data, dz, freql.data, icwnc.data)
            reff[reff==0] = np.nan
            # calculate mean optical depth
            cod = np.sum(cod_a.data,axis=1)
            # cod[solin==0] = np.nan
            # cod from MODIS simulator
            cod_m = cod_m[:,x_idx]
            
            reff_2 = xr.DataArray(data=reff,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="mean cloud liquid effective radius",units="um"),)
            cod_2 = xr.DataArray(data=cod,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="column-total cloud optical depth",units="N/A"),)
            reff_mean = xr.concat([reff_mean, reff_2], dim="time")
            cod_mean = xr.concat([cod_mean, cod_2], dim="time")
        else:
            cod_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cod',attrs={'units':'dummy_unit','long_name':'Dummy'})
            reff_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='reff',attrs={'units':'dummy_unit','long_name':'Dummy'})
            cod_m = xr.DataArray(np.zeros(len(e3smtime_i))*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        
        # mean cloud droplet number concentration
        req_vlist = ['CDNUMC']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            cdnc_col = e3smdata['CDNUMC'+'_'+E3SMdomain_range].load()
            cdnc_col = cdnc_col[:,x_idx]
            cdnc = cdnc_col/np.sum(weight,axis=1)
            cdnc[cdnc >2e9] = np.nan
            cdnc = xr.DataArray(data=cdnc,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3"),)
            cdnc_mean = xr.concat([cdnc_mean, cdnc], dim="time")
        else:
            cdnc_mean = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cdnc_mean',attrs={'units':'dummy_unit','long_name':'Dummy'})
            
        # cloud droplet number concentration retrieved like Ndrop and Bennartz 2007
        req_vlist = ['TGCLDLWP']
        req_vlist = ["{}_{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        
        if len(matched_vlist) == len(req_vlist):
            lwp = e3smdata['TGCLDLWP'+'_'+E3SMdomain_range][:,x_idx].data
            e3sm_cloud_depth[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
            T_cldtop[z_cldtop>5000] = np.nan  # remove deep clouds with cloud top >5km
            nd_arm = calc_cdnc_ARM(lwp, cod_m, e3sm_cloud_depth)
            nd_sat = calc_cdnc_VISST(lwp, T_cldtop, cod_m, adiabaticity=0.8)
            nd_arm = xr.DataArray(data=nd_arm*1e6,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                            description='Retrieved using ARM Ndrop algorithm'),)
            nd_sat = xr.DataArray(data=nd_sat*1e6,  dims=["time"],
                coords=dict(time=(["time"], e3smtime_i)),
                attrs=dict(long_name="mean cloud water number concentration",units="#/m3",\
                            description='Retrieved using Bennartz(2007) algorithm, also used for VISST data'),)
            cdnc_arm = xr.concat([cdnc_arm, nd_arm], dim="time")
            cdnc_sat = xr.concat([cdnc_sat, nd_sat], dim="time")
        else:
            cdnc_arm = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cdnc_arm',attrs={'units':'dummy_unit','long_name':'Dummy'})
            cdnc_sat = xr.DataArray(np.zeros(len(e3smtime))*np.nan,name='cdnc_sat',attrs={'units':'dummy_unit','long_name':'Dummy'})
    
        
        # all other 2D (surface and vertical integrated) variables
        for varname in variable2d_names:
            try:
                var = e3smdata[varname + '_'+E3SMdomain_range].load()
                var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
            except:
                var = xr.DataArray(np.zeros((len(e3smtime_i),len_ncol))*np.nan,name=varname,\
                               dims=["time","ncol_"+E3SMdomain_range],coords={"time":e3smtime_i,"ncol_"+E3SMdomain_range:np.arange(len_ncol)},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            vv = variable_names.index(varname)
            variables[vv] = xr.concat([variables[vv], var[:,x_idx]],dim='time')
        
        # all other 3D (with vertical level) variables at the lowest model level
        for varname in variable3d_names:
            try:
                var = e3smdata[varname + '_'+E3SMdomain_range].load()
                var.coords['time'] = var.indexes['time'].to_datetimeindex() # change time to standard datetime64 format
            except:
                var = xr.DataArray(np.zeros((len(e3smtime_i),len_lev,len_ncol))*np.nan,name=varname,\
                               dims=["time","lev","ncol_"+E3SMdomain_range],coords={"time":e3smtime_i,"lev":e3smdata['lev'],"ncol_"+E3SMdomain_range:e3smdata['ncol_'+E3SMdomain_range]},\
                               attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            vv = variable_names.index(varname)
            variables[vv] = xr.concat([variables[vv], var[:,-1,x_idx]],dim='time')
    
        e3smdata.close()
    # put all variables into the list
    # aerosol composition    
    variable_names = variable_names + ['bc','dst','mom','ncl','pom','so4','soa']
    variables = variables + [bc_all,dst_all,mom_all,ncl_all,pom_all,so4_all,soa_all]
    # aerosol size and number
    NCN3 = xr.DataArray(data=np.nansum(NCNall[3:, :], 0),  dims=["time"],
        coords=dict(time=(["time"], e3smtime)),
        attrs=dict(long_name="Aerosol number concentration for size >3nm",units="#/m3"),)    # >3nm
    NCN10 = xr.DataArray(data=np.nansum(NCNall[10:, :], 0),  dims=["time"],
        coords=dict(time=(["time"], e3smtime)),
        attrs=dict(long_name="Aerosol number concentration for size >10nm",units="#/m3"),)     # >10nm
    NCN100 = xr.DataArray(data=np.nansum(NCNall[100:, :], 0),  dims=["time"],
        coords=dict(time=(["time"], e3smtime)),
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
    variable_names = variable_names + ['cbt','ctt','cbh','cth','clddepth']
    variables = variables + [cbt, ctt, cbh, cth, cloud_depth]
    
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
    varbls = ['TGCLDIWP','TGCLDLWP']
    varbls = ["{}_{}".format(i,E3SMdomain_range) for i in varbls]
    varbls = list(set(av_vars).intersection(varbls))
    if len(varbls) > 0:
        for vv in ['TGCLDIWP','TGCLDLWP']:
            variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *1000
            variables[variable_names.index(vv)].attrs['units']='g/m2'
    # cloud fraction
    varbls = ['CLDTOT','CLDLOW','CLDMED','CLDHGH']
    varbls = ["{}_{}".format(i,E3SMdomain_range) for i in varbls]
    varbls = list(set(av_vars).intersection(varbls))
    if len(varbls) > 0:
        for vv in ['CLDTOT','CLDLOW','CLDMED','CLDHGH']:
            variables[variable_names.index(vv)].data = variables[variable_names.index(vv)].data *100
            variables[variable_names.index(vv)].attrs['units']='%'
    
    #%% re-shape the data into pre-defined resolution
    variables_new = list()
    #1d variable. only numpy.interp can keep some single-point values (see Nd_mean)
    for var in variables:
        var_new = np.interp(np.int64(time_new), np.int64(e3smtime), var, left=np.nan, right=np.nan)
        variables_new.append(var_new)
    # treat variables with other dimensions (e.g., size distribution)    
    f = interp1d(np.int64(e3smtime), NCNall, bounds_error=False)
    NCNall_new = f(np.int64(time_new))
    
    # %% output extacted file
    varall_1d = {
            variable_names[vv]: ('time', np.float32(variables_new[vv])) for vv in range(len(variable_names))
    }
    varall_2d = {
            'NCNall': (['size','time',], np.float32(NCNall_new))
    }
    varall_1d.update(varall_2d)
    outfile = output_path + output_filehead + '_sfc.nc'
    print('output file '+outfile)
    ds = xr.Dataset( varall_1d,
                      coords={'time': ('time', time_new), 'size':('size', np.arange(1,3001))})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    for vv in range(len(variable_names)):
        ds[variable_names[vv]].attrs["long_name"] = variables[vv].long_name
        ds[variable_names[vv]].attrs["units"] = variables[vv].units
        ds[variable_names[vv]].attrs["description"] = "variables at surface, TOA or the lowest model level"
    ds['NCNall'].attrs["long_name"] = NCNall.long_name
    ds['NCNall'].attrs["units"] = NCNall.units
    ds['NCNall'].attrs["description"] = 'calculated from modal information into 1nm increment'
    
    ds.attrs["title"] = 'preprocessed E3SM data at surface, TOA or the lowest model level'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
