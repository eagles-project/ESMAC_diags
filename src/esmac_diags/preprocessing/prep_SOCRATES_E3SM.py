"""
prepare E3SM data for CSET aircraft tracks
options of average data into coarser resolution
"""

import glob
import os
import numpy as np
import pandas as pd
import xarray as xr
import time as ttt
from scipy.special import gamma
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d
from esmac_diags.subroutines.read_aircraft import read_RF_NCAR
from esmac_diags.subroutines.specific_data_treatment import find_nearest 
from esmac_diags.subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm
from netCDF4 import Dataset

# input_path = '../../../data/ACEENA/model/E3SMv1_out/'
# output_path = 'C:/Users/tang357/Downloads/SOCRATES/'
# input_filehead = 'E3SMv1_SO'
# output_filehead = 'E3SMv1_SOCRATES'

# lev_out=np.arange(25.,1001,25.)
# height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
#                     1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
#                     3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
#                     10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

# dt = 3600
    
# dt = 60
# RFpath = '../../../data/SOCRATES/obs/aircraft/aircraft_lowrate/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, 
                      RFpath, dt=60):
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
    RFpath : str
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
    # SOCRATES
    E3SMdomain_range = '133e_to_164e_42s_to_63s'    # domain range in E3SM regional output
    
    #%% find all data
    lst = glob.glob(RFpath + 'RF*.PNI.nc')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = filename.split('.')
        date = fname[-4]
        print(date)
        
        #%% read in aircraft data
        (time, height, timeunit, hunit, hlongname, cellsize, cellunit) = read_RF_NCAR(filename, 'ALT')
        (time, lat, timeunit, latunit, latlongname, cellsize, cellunit) = read_RF_NCAR(filename, 'LAT')
        (time, lon, timeunit, lonunit, lonlongname, cellsize, cellunit) = read_RF_NCAR(filename, 'LON')
        lon[lon<0] = lon[lon<0] + 360
        timestr = timeunit.split(' ')[2]
        
        # re-shape the data into coarser resolution for output
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        lon_new = median_time_1d(time, lon, time_new)
        lat_new = median_time_1d(time, lat, time_new)
        height_new = median_time_1d(time, height, time_new)
        if sum(np.isnan(lat_new))>0.1*len(lat_new):
            raise ValueError('Too many (>10%) NaNs in RF data, check the file: '+filename)
        else:  # interpolate to fill NaNs
            lon_new = pd.Series(lon_new).interpolate().tolist()
            lat_new = pd.Series(lat_new).interpolate().tolist()
            height_new = np.array(pd.Series(height_new).interpolate().tolist())
        
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
        
        # SOCRATES flight extend to the next day
        lstm = glob.glob(input_path + input_filehead+'.*'+timestr+'-00000.nc') + \
            glob.glob(input_path + input_filehead+'.*'+str(np.datetime64(timestr)+1)+'-00000.nc')
        if len(lstm)!=2:
            print(lstm)
            raise ValueError('Should only contain two files')
        e3smdata = xr.open_mfdataset(lstm)
        e3smtime = e3smdata.indexes['time'].to_datetimeindex()
        lonm = e3smdata['lon'+'_'+E3SMdomain_range].load()
        latm = e3smdata['lat'+'_'+E3SMdomain_range].load()
        z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
        # change time format into seconds of the day
        timem = np.float64((e3smtime - e3smtime[0]).seconds)
        
        # variables for calculating aerosol size
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
        if len(hyam.shape)>1:   # when use xr.open_mfdataset, hyam will add a time dimension
            hyam = hyam[0,:]
            hybm = hybm[0,:]
        T = e3smdata['T'+'_'+E3SMdomain_range].load()
        PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
        Pres = np.nan*T
        zlen = T.shape[1]
        for kk in range(zlen):
            Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
            
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
        
        # droplet size distribution
        nd_cld = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
        lmda = e3smdata['lambda_cloud'+'_'+E3SMdomain_range].load()
        mu = e3smdata['mu_cloud'+'_'+E3SMdomain_range].load()
        
        
        # other variables
        for varname in variable3d_names:
            var = e3smdata[varname+'_'+E3SMdomain_range].load()
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
        time_o.units = "Seconds since " + timestr + ' 00:00 UTC'
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
        f.modelfile = lstm[0].split('/')[-1]
        f.date = ttt.ctime(ttt.time())
        
        f.close()
    
