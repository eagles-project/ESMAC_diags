"""
prepare E3SM data for CSET aircraft tracks
options of average data into coarser resolution
"""

import glob
import os
import yaml
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
def prep_E3SM_flight(input_path, input2d_filehead, input3d_filehead, input3d_dryaerosol_filehead, input3d_cloudaerosol_filehead, output_path, output_filehead, 
                      iwgpath, dt=60, config=config):
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
    E3SMdomain_range = config['E3SMdomain_range']    # domain range in E3SM regional output
    
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
        time_new_int = np.array([int(i) for i in time_new])
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
        # variable3d_names = [config['T'], config['Q'], config['U'], config['V'], config['Z'], 
        #                     config['CF'], config['CFLIQ'], config['NC'], config['NI']] #,config['QI'], config['QC']
        variable3d_names = [config['T'], config['Z'], 
                            config['CF'], config['CFLIQ'], config['CFICE'], config['NC']]

        if config['rain_output'] == True:
            variable3d_names.append(config['QR'])
            variable3d_names.append(config['NR'])
        if config['reff_output'] == True:
            variable3d_names.append(config['REL'])
        if config['ccn_output'] == True:
            variable3d_names.append(config['CCN1'])
            variable3d_names.append(config['CCN2'])
        #     variable3d_names.append(config['CCN5'])
        #     variable3d_names.append(config['CCN10'])
        
        variables = list()
        variables_new = list()
        for varname in variable3d_names:
            variables_new.append([])
        NCNall = np.empty((3000,0))
        p = list()      # pressure
        cwc = list()
        iwc = list()
        if config['rain_output'] == True:
            rwc = list()
        if config['aerosol_output'] == True:
            NCNall = np.empty((3000,0))
        if config['composition_output'] == True:
            bc_all  = list()
            dst_all = list()
            mom_all = list()
            ncl_all = list()
            pom_all = list()
            so4_all = list()
            soa_all = list()
            phi_all = np.empty((999,0))
        if config['ccn_output'] == True:
            ccn1_all = list()
            ccn2_all = list()
            # ccn5_all = list()
            # ccn10_all = list()
        
        # SOCRATES flight extend to the next day
        lst3d = glob.glob(input_path + input3d_filehead+'.*'+timestr[0]+'-*.nc')
        lst3d.sort()
        if len(lst3d)==0:
            raise ValueError('No model output files on flight day')

        #string split model file names to get sssss arrays
        timesecstr_int = [int(filename.split('.')[-2].split('-')[-1]) for filename in lst3d]
        #find last file time before first iwg flight time
        firstfile_ind = np.where(timesecstr_int <= np.min(time_new_int))[-1][-1]
        #Add any files to the list with times during the flight
        lastfile_ind = np.where(np.logical_and(timesecstr_int > np.min(time_new_int), timesecstr_int <= np.max(time_new_int)))[0]
        file_inds = np.concatenate([np.array(firstfile_ind).reshape(1), lastfile_ind])
        
        # lst3d = glob.glob(input_path + input_filehead+'.*'+timestr+'-00000.nc') + \
        #     glob.glob(input_path + input_filehead+'.*'+str(np.datetime64(timestr)+1)+'-00000.nc')
        # if len(lst3d)!=2:
        #     print(lst3d)
        #     raise ValueError('Should only contain two files')
        # e3smdata = xr.open_mfdataset(lst3d)
        # e3smtime = e3smdata.indexes['time'].to_datetimeindex()
        # lonm = e3smdata['lon'+'_'+E3SMdomain_range].load()
        # latm = e3smdata['lat'+'_'+E3SMdomain_range].load()
        # z3 = e3smdata['Z3'+'_'+E3SMdomain_range].load()
        # # change time format into seconds of the day
        # timem = np.float64((e3smtime - e3smtime[0]).seconds)

        #define model variable arrays with first model file
        print(lst3d[file_inds[0]])
        e3smdata3d = xr.open_dataset(lst3d[file_inds[0]])
        e3smdata3d = e3smdata3d.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
        e3smtime = e3smdata3d.indexes[config['time_dim']].to_datetimeindex()
        tmplonm = e3smdata3d[config['LON']+E3SMdomain_range].load()
        tmplatm = e3smdata3d[config['LAT']+E3SMdomain_range].load()

        #also only include portions of the domain that include the flight track by using the flight track min/max lat and lon +/- delta defined in config file
        lon[lon < -180] = np.nan
        lon[lon > 180] = np.nan
        lat[lat < -90] = np.nan
        lat[lat > 90] = np.nan
        minlon = np.nanmin(lon)+360 - config['model_grid_deg']
        maxlon = np.nanmax(lon)+360 + config['model_grid_deg']
        minlat = np.nanmin(lat) - config['model_grid_deg']
        maxlat = np.nanmax(lat) + config['model_grid_deg']
        print(str(minlon), str(maxlon), str(minlat),  str(maxlat))
        latlon_ind = np.where(np.logical_and(np.logical_and(np.logical_and(tmplonm >= minlon, tmplonm <= maxlon), tmplatm >= minlat), tmplatm <= maxlat))[0]
        lonm = tmplonm[latlon_ind]
        latm = tmplatm[latlon_ind]
        z3 = e3smdata3d[config['Z']+E3SMdomain_range][:,:,latlon_ind,...].load()
        T = e3smdata3d[config['T']+E3SMdomain_range][:,:,latlon_ind,...].load()

        if config['pres_output'] == False:
            P0 = e3smdata3d[config['P0']].load()
            hyam = e3smdata3d[config['HYAM']].load()
            hybm = e3smdata3d[config['HYBM']].load()
            PS = e3smdata3d[config['PS']+E3SMdomain_range][:,latlon_ind,...].load()
            # Pres = np.nan*T
            # zlen = T.shape[1]
            Pres = xr.full_like(T, np.nan)
            Pres = Pres.assign_attrs(units='Pa',long_name='Pressure',standard_name='air_pressure')
            zlen = T.sizes[config['vert_dim']]
            for kk in range(zlen):
                Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
        else:
            Pres = e3smdata3d[config['PRES']+E3SMdomain_range][:,:,latlon_ind,...].load()
      
        # change time format into seconds of the day
        # timem = np.float64((e3smtime - e3smtime[0]).seconds) # this only works for the first time being 0Z
        timem = np.float64(e3smtime.hour + e3smtime.minute + e3smtime.second)*3600

        # # variables for calculating aerosol size
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
        # T = e3smdata['T'+'_'+E3SMdomain_range].load()
        # PS = e3smdata['PS'+'_'+E3SMdomain_range].load()
        # Pres = np.nan*T
        # zlen = T.shape[1]
        # for kk in range(zlen):
        #     Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS

        # Get all simulated variables
        vlist = list(e3smdata3d.variables.keys())
        av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)

        # condensate mass and number
        # req_vlist = [config['QC'], config['QI']]
        req_vlist = [config['QC']]
        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
        matched_vlist = list(set(av_vars).intersection(req_vlist))
        if len(matched_vlist) == len(req_vlist):
            print('\nAnalyzing Condensate')
            qc = e3smdata3d[config['QC']+E3SMdomain_range][:,:,latlon_ind,...].load()
            # qi = e3smdata3d[config['QI']+E3SMdomain_range][:,:,latlon_ind,...].load()
        else:
            qc = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            # qi = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

        if config['rain_output'] == True:
            req_vlist = [config['QR']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars).intersection(req_vlist))
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing Rain')
                qr = e3smdata3d[config['QR']+E3SMdomain_range][:,:,latlon_ind,...].load()
            else:
                qr = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
          
        # droplet size distribution
        if config['dsd_output'] == True:
          req_vlist = [config['LAMBDA_CLOUD'], config['MU_CLOUD']]
          req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
          matched_vlist = list(set(av_vars).intersection(req_vlist))
          if len(matched_vlist) == len(req_vlist):
              print('\nAnalyzing DSD')
              nd_cld = e3smdata3d[config['NC']+E3SMdomain_range][:,:,latlon_ind,...].load()
              lmda = e3smdata3d[config['LAMBDA_CLOUD']+E3SMdomain_range][:,:,latlon_ind,...].load()
              mu = e3smdata3d[config['MU_CLOUD']+E3SMdomain_range][:,:,latlon_ind,...].load()
          else:
              nd_cld = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              lmda = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
              mu = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
        
        # other variables
        for varname in variable3d_names:
            try:
                var = e3smdata3d[varname + E3SMdomain_range][:,:,latlon_ind,...].load()
                var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
            except:
                var = xr.DataArray(np.zeros(z3.shape)*np.nan,name=varname,\
                                   dims=[config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime,config['vert_dim']:e3smdata3d[config['vert_dim']],config['latlon_dim']+E3SMdomain_range:e3smdata3d[config['latlon_dim']+E3SMdomain_range]},\
                                   attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
            variables.append(var)
        
        e3smdata3d.close()

        # variables for calculating aerosol size
        if config['aerosol_output'] == True:
            lst3d_dryaer = glob.glob(input_path + input3d_dryaerosol_filehead+'.*'+timestr[0]+'-*.nc')
            lst3d_dryaer.sort()
            if len(lst3d_dryaer)==0:
                raise ValueError('No model output files on flight day')
  
            #define model variable arrays with first model file
            print(lst3d_dryaer[file_inds[0]])
            e3smdata3d_dryaer = xr.open_dataset(lst3d_dryaer[file_inds[0]])
            e3smdata3d_dryaer = e3smdata3d_dryaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location

            vlist = list(e3smdata3d_dryaer.variables.keys())
            av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)

            req_vlist = [config['num_a1'], config['num_a2'], config['num_a3'], config['num_a4']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars).intersection(req_vlist))
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing for aerosol size (1)')
                num_a1 = e3smdata3d_dryaer[config['num_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                num_a2 = e3smdata3d_dryaer[config['num_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                num_a3 = e3smdata3d_dryaer[config['num_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                num_a4 = e3smdata3d_dryaer[config['num_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
            else:
                num_a1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                num_a2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                num_a3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                num_a4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

            if config['dgnum_output_combined'] == False:
                req_vlist = [config['dgnd_a01'], config['dgnd_a02'], config['dgnd_a03'], config['dgnd_a04']]
                req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                matched_vlist = list(set(av_vars).intersection(req_vlist))
                if len(matched_vlist) == len(req_vlist):
                    print('\nAnalyzing for aerosol size (2)')
                    num_dn1 = e3smdata3d_dryaer[config['dgnd_a01']+E3SMdomain_range][:,:,latlon_ind,...].load()
                    num_dn2 = e3smdata3d_dryaer[config['dgnd_a02']+E3SMdomain_range][:,:,latlon_ind,...].load()
                    num_dn3 = e3smdata3d_dryaer[config['dgnd_a03']+E3SMdomain_range][:,:,latlon_ind,...].load()
                    num_dn4 = e3smdata3d_dryaer[config['dgnd_a04']+E3SMdomain_range][:,:,latlon_ind,...].load()
                else:
                    num_dn1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    num_dn2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    num_dn3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    num_dn4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
            else:
                req_vlist = [config['dgnum']]
                req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                matched_vlist = list(set(av_vars).intersection(req_vlist))
                if len(matched_vlist) == len(req_vlist):
                    print('\nAnalyzing for aerosol size (2)')
                    num_dn1 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,0].load()
                    num_dn2 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,1].load()
                    num_dn3 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,2].load()
                    num_dn4 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,3].load()
                else:
                    num_dn1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    num_dn2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    num_dn3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    num_dn4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})

            if config['composition_output'] == True:
              # aerosol composition
              req_vlist = [config['bc_a1'], config['bc_a3'], config['bc_a4'], config['dst_a1'], config['dst_a3'], config['mom_a1'], \
                           config['mom_a2'], config['mom_a3'], config['mom_a4'], config['ncl_a1'], config['ncl_a2'], config['ncl_a3'], \
                           config['pom_a1'], config['pom_a3'], config['pom_a4'], config['so4_a1'], config['so4_a2'], config['so4_a3'], \
                           config['soa_a1'], config['soa_a2'], config['soa_a3']]
              req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
              matched_vlist = list(set(av_vars).intersection(req_vlist))
          
              if len(matched_vlist) == len(req_vlist):
                  print('\nAnalyzing for aerosol composition')
                  bc_a1 = e3smdata3d_dryaer[config['bc_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  bc_a3 = e3smdata3d_dryaer[config['bc_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  bc_a4 = e3smdata3d_dryaer[config['bc_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  dst_a1 = e3smdata3d_dryaer[config['dst_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  dst_a3 = e3smdata3d_dryaer[config['dst_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  mom_a1 = e3smdata3d_dryaer[config['mom_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  mom_a2 = e3smdata3d_dryaer[config['mom_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  mom_a3 = e3smdata3d_dryaer[config['mom_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  mom_a4 = e3smdata3d_dryaer[config['mom_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  ncl_a1 = e3smdata3d_dryaer[config['ncl_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  ncl_a2 = e3smdata3d_dryaer[config['ncl_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  ncl_a3 = e3smdata3d_dryaer[config['ncl_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  pom_a1 = e3smdata3d_dryaer[config['pom_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  pom_a3 = e3smdata3d_dryaer[config['pom_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  pom_a4 = e3smdata3d_dryaer[config['pom_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  so4_a1 = e3smdata3d_dryaer[config['so4_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  so4_a2 = e3smdata3d_dryaer[config['so4_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  so4_a3 = e3smdata3d_dryaer[config['so4_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  soa_a1 = e3smdata3d_dryaer[config['soa_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  soa_a2 = e3smdata3d_dryaer[config['soa_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                  soa_a3 = e3smdata3d_dryaer[config['soa_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
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

        if config['ccn_output'] == True:
            lst3d_wetaer = glob.glob(input_path + input3d_cloudaerosol_filehead+'.*'+timestr[0]+'-*.nc')
            lst3d_wetaer.sort()
            if len(lst3d_wetaer)==0:
                raise ValueError('No model output files on flight day')
    
            #define model variable arrays with first model file
            print(lst3d_wetaer[file_inds[0]])
            e3smdata3d_wetaer = xr.open_dataset(lst3d_wetaer[file_inds[0]])
            e3smdata3d_wetaer = e3smdata3d_wetaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location

            vlist = list(e3smdata3d_wetaer.variables.keys())
            av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
      
            req_vlist = [config['CCN1'], config['CCN2']]#, config['CCN5'], config['CCN10']]
            req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
            matched_vlist = list(set(av_vars).intersection(req_vlist))
            if len(matched_vlist) == len(req_vlist):
                print('\nAnalyzing CCN')
                ccn1 = e3smdata3d_wetaer[config['CCN1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                ccn2 = e3smdata3d_wetaer[config['CCN2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                # ccn5 = e3smdata3d_wetaer[config['CCN5']+E3SMdomain_range][:,:,latlon_ind,...].load()
                # ccn10 = e3smdata3d_wetaer[config['CCN10']+E3SMdomain_range][:,:,latlon_ind,...].load()
            else:
                ccn1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                ccn2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                # ccn5 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                # ccn10 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})


        #add model variable arrays at additional times if there are additional model output files within the flight period
        if len(file_inds) > 1:
            for ii in np.arange(len(file_inds)-1):
                print(lst3d[file_inds[ii+1]])
                e3smdata3d = xr.open_dataset(lst3d[file_inds[ii+1]])
                e3smdata3d = e3smdata3d.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
                e3smtime_i = e3smdata3d.indexes[config['time_dim']].to_datetimeindex()
                e3smtime = np.hstack((e3smtime, e3smtime_i))
                
                newz3 = e3smdata3d[config['Z']+E3SMdomain_range][:,:,latlon_ind,...].load()
                z3 = xr.concat([z3, newz3], dim=config['time_dim'])
    
                newT = e3smdata3d[config['T']+E3SMdomain_range][:,:,latlon_ind,...].load()
                T = xr.concat([T, newT], dim=config['time_dim'])
    
                if config['pres_output'] == False:
                    PS = e3smdata3d[config['PS']+E3SMdomain_range][:,latlon_ind,...].load()
                    newPres = xr.full_like(T, np.nan)
                    newPres = Pres.assign_attrs(units='Pa',long_name='Pressure',standard_name='air_pressure')
                    zlen = T.sizes[config['vert_dim']]
                    for kk in range(zlen):
                        newPres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
                else:
                    newPres = e3smdata3d[config['PRES']+E3SMdomain_range][:,:,latlon_ind,...].load()
    
                Pres = xr.concat([Pres, newPres], dim=config['time_dim'])
          
                # change time format into seconds of the day
                timem = np.concatenate([timem, np.float64(e3smtime_i.hour + e3smtime_i.minute + e3smtime_i.second)*3600])

                vlist = list(e3smdata3d.variables.keys())
                av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
      
                # condensate mass and number
                req_vlist = [config['QC']]#, config['QI']]
                req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                matched_vlist = list(set(av_vars).intersection(req_vlist))
                if len(matched_vlist) == len(req_vlist):
                    print('\nAnalyzing Condensate')
                    new_qc = e3smdata3d[config['QC']+E3SMdomain_range][:,:,latlon_ind,...].load()
                    # new_qi = e3smdata3d[config['QI']+E3SMdomain_range][:,:,latlon_ind,...].load()
                    qc = xr.concat([qc, new_qc], dim=config['time_dim'])
                    # qi = xr.concat([qi, new_qi], dim=config['time_dim'])
                else:
                    new_qc = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    # new_qi = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                    qc = xr.concat([qc, new_qc], dim=config['time_dim'])
                    # qi = xr.concat([qi, new_qi], dim=config['time_dim'])
        
                if config['rain_output'] == True:
                    req_vlist = [config['QR']]
                    req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                    matched_vlist = list(set(av_vars).intersection(req_vlist))
                    if len(matched_vlist) == len(req_vlist):
                        print('\nAnalyzing Rain')
                        new_qr = e3smdata3d[config['QC']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        qr = xr.concat([qr, new_qr], dim=config['time_dim'])
                    else:
                        new_qr = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        qr = xr.concat([qr, new_qr], dim=config['time_dim'])
                  
                # droplet size distribution
                if config['dsd_output'] == True:
                  req_vlist = [config['LAMBDA_CLOUD'], config['MU_CLOUD']]
                  req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                  matched_vlist = list(set(av_vars).intersection(req_vlist))
                  if len(matched_vlist) == len(req_vlist):
                      print('\nAnalyzing DSD')
                      new_nd_cld = e3smdata3d[config['NC']+E3SMdomain_range][:,:,latlon_ind,...].load()
                      new_lmda = e3smdata3d[config['LAMBDA_CLOUD']+E3SMdomain_range][:,:,latlon_ind,...].load()
                      new_mu = e3smdata3d[config['MU_CLOUD']+E3SMdomain_range][:,:,latlon_ind,...].load()
                      nd_cld = xr.concat([nd_cld, new_nd_cld], dim=config['time_dim'])
                      lmda = xr.concat([lmda, new_lmda], dim=config['time_dim'])
                      mu = xr.concat([mu, new_mu], dim=config['time_dim'])
                  else:
                      new_nd_cld = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                      new_lmda = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                      new_mu = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                      nd_cld = xr.concat([nd_cld, new_nd_cld], dim=config['time_dim'])
                      lmda = xr.concat([lmda, new_lmda], dim=config['time_dim'])
                      mu = xr.concat([mu, new_mu], dim=config['time_dim'])
                
                # other variables
                for varname in variable3d_names:
                    try:
                        var = e3smdata3d[varname + E3SMdomain_range][:,:,latlon_ind,...].load()
                        var.coords[config['time_dim']] = var.indexes[config['time_dim']].to_datetimeindex() # change time to standard datetime64 format
                    except:
                        var = xr.DataArray(np.zeros(z3.shape)*np.nan,name=varname,\
                                           dims=[config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range],coords={config['time_dim']:e3smtime_i,config['vert_dim']:e3smdata3d[config['vert_dim']],config['latlon_dim']+E3SMdomain_range:e3smdata3d[config['latlon_dim']+E3SMdomain_range]},\
                                           attrs={'units':'dummy_unit','long_name':'dummy_long_name'})
                    
                    vv = variable3d_names.index(varname)
                    variables[vv] = xr.concat([variables[vv], var],dim=config['time_dim'])            
              
                e3smdata3d.close()  
              
                # variables for calculating aerosol size
                if config['aerosol_output'] == True:
                    print(lst3d_dryaer[file_inds[ii+1]])
                    e3smdata3d_dryaer = xr.open_dataset(lst3d_dryaer[file_inds[ii+1]])
                    e3smdata3d_dryaer = e3smdata3d_dryaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
                    
                    vlist = list(e3smdata3d_dryaer.variables.keys())
                    av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
              
                    req_vlist = [config['num_a1'], config['num_a2'], config['num_a3'], config['num_a4']]
                    req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                    matched_vlist = list(set(av_vars).intersection(req_vlist))
                    if len(matched_vlist) == len(req_vlist):
                        print('\nAnalyzing for aerosol size (1)')
                        new_num_a1 = e3smdata3d_dryaer[config['num_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        new_num_a2 = e3smdata3d_dryaer[config['num_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        new_num_a3 = e3smdata3d_dryaer[config['num_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        new_num_a4 = e3smdata3d_dryaer[config['num_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        num_a1 = xr.concat([num_a1, new_num_a1], dim=config['time_dim'])
                        num_a2 = xr.concat([num_a2, new_num_a2], dim=config['time_dim'])
                        num_a3 = xr.concat([num_a3, new_num_a3], dim=config['time_dim'])
                        num_a4 = xr.concat([num_a4, new_num_a4], dim=config['time_dim'])
                    else:
                        new_num_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        new_num_a2 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        new_num_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        new_num_a4 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        num_a1 = xr.concat([num_a1, new_num_a1], dim=config['time_dim'])
                        num_a2 = xr.concat([num_a2, new_num_a2], dim=config['time_dim'])
                        num_a3 = xr.concat([num_a3, new_num_a3], dim=config['time_dim'])
                        num_a4 = xr.concat([num_a4, new_num_a4], dim=config['time_dim'])

                    if config['dgnum_output_combined'] == False:
                        req_vlist = [config['dgnd_a01'], config['dgnd_a02'], config['dgnd_a03'], config['dgnd_a04']]
                        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                        matched_vlist = list(set(av_vars).intersection(req_vlist))
                        if len(matched_vlist) == len(req_vlist):
                            print('\nAnalyzing for aerosol size (2)')
                            new_num_dn1 = e3smdata3d_dryaer[config['dgnd_a01']+E3SMdomain_range][:,:,latlon_ind,...].load()
                            new_num_dn2 = e3smdata3d_dryaer[config['dgnd_a02']+E3SMdomain_range][:,:,latlon_ind,...].load()
                            new_num_dn3 = e3smdata3d_dryaer[config['dgnd_a03']+E3SMdomain_range][:,:,latlon_ind,...].load()
                            new_num_dn4 = e3smdata3d_dryaer[config['dgnd_a04']+E3SMdomain_range][:,:,latlon_ind,...].load()
                            num_dn1 = xr.concat([num_dn1, new_num_dn1], dim=config['time_dim'])
                            num_dn2 = xr.concat([num_dn2, new_num_dn2], dim=config['time_dim'])
                            num_dn3 = xr.concat([num_dn3, new_num_dn3], dim=config['time_dim'])
                            num_dn4 = xr.concat([num_dn4, new_num_dn4], dim=config['time_dim'])
                        else:
                            new_num_dn1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            new_num_dn2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            new_num_dn3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            new_num_dn4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            num_dn1 = xr.concat([num_dn1, new_num_dn1], dim=config['time_dim'])
                            num_dn2 = xr.concat([num_dn2, new_num_dn2], dim=config['time_dim'])
                            num_dn3 = xr.concat([num_dn3, new_num_dn3], dim=config['time_dim'])
                            num_dn4 = xr.concat([num_dn4, new_num_dn4], dim=config['time_dim'])
                    else:
                        req_vlist = [config['dgnum']]
                        req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                        matched_vlist = list(set(av_vars).intersection(req_vlist))
                        if len(matched_vlist) == len(req_vlist):
                            print('\nAnalyzing for aerosol size (2)')
                            new_num_dn1 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,0].load()
                            new_num_dn2 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,1].load()
                            new_num_dn3 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,2].load()
                            new_num_dn4 = e3smdata3d_dryaer[config['dgnum']+E3SMdomain_range][:,:,latlon_ind,3].load()
                            num_dn1 = xr.concat([num_dn1, new_num_dn1], dim=config['time_dim'])
                            num_dn2 = xr.concat([num_dn2, new_num_dn2], dim=config['time_dim'])
                            num_dn3 = xr.concat([num_dn3, new_num_dn3], dim=config['time_dim'])
                            num_dn4 = xr.concat([num_dn4, new_num_dn4], dim=config['time_dim'])
                        else:
                            new_num_dn1 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            new_num_dn2 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            new_num_dn3 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            new_num_dn4 = xr.DataArray(np.zeros(z3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                            num_dn1 = xr.concat([num_dn1, new_num_dn1], dim=config['time_dim'])
                            num_dn2 = xr.concat([num_dn2, new_num_dn2], dim=config['time_dim'])
                            num_dn3 = xr.concat([num_dn3, new_num_dn3], dim=config['time_dim'])
                            num_dn4 = xr.concat([num_dn4, new_num_dn4], dim=config['time_dim'])

                    if config['composition_output'] == True:
                      # aerosol composition
                      req_vlist = [config['bc_a1'], config['bc_a3'], config['bc_a4'], config['dst_a1'], config['dst_a3'], config['mom_a1'], \
                                   config['mom_a2'], config['mom_a3'], config['mom_a4'], config['ncl_a1'], config['ncl_a2'], config['ncl_a3'], \
                                   config['pom_a1'], config['pom_a3'], config['pom_a4'], config['so4_a1'], config['so4_a2'], config['so4_a3'], \
                                   config['soa_a1'], config['soa_a2'], config['soa_a3']]
                      req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                      matched_vlist = list(set(av_vars).intersection(req_vlist))
                  
                      if len(matched_vlist) == len(req_vlist):
                          print('\nAnalyzing for aerosol composition')
                          new_bc_a1 = e3smdata3d_dryaer[config['bc_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_bc_a3 = e3smdata3d_dryaer[config['bc_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_bc_a4 = e3smdata3d_dryaer[config['bc_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_dst_a1 = e3smdata3d_dryaer[config['dst_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_dst_a3 = e3smdata3d_dryaer[config['dst_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_mom_a1 = e3smdata3d_dryaer[config['mom_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_mom_a2 = e3smdata3d_dryaer[config['mom_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_mom_a3 = e3smdata3d_dryaer[config['mom_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_mom_a4 = e3smdata3d_dryaer[config['mom_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_ncl_a1 = e3smdata3d_dryaer[config['ncl_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_ncl_a2 = e3smdata3d_dryaer[config['ncl_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_ncl_a3 = e3smdata3d_dryaer[config['ncl_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_pom_a1 = e3smdata3d_dryaer[config['pom_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_pom_a3 = e3smdata3d_dryaer[config['pom_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_pom_a4 = e3smdata3d_dryaer[config['pom_a4']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_so4_a1 = e3smdata3d_dryaer[config['so4_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_so4_a2 = e3smdata3d_dryaer[config['so4_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_so4_a3 = e3smdata3d_dryaer[config['so4_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_soa_a1 = e3smdata3d_dryaer[config['soa_a1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_soa_a2 = e3smdata3d_dryaer[config['soa_a2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          new_soa_a3 = e3smdata3d_dryaer[config['soa_a3']+E3SMdomain_range][:,:,latlon_ind,...].load()
                          bc_a1 = xr.concat([bc_a1, new_bc_a1], dim=config['time_dim'])
                          bc_a3 = xr.concat([bc_a3, new_bc_a3], dim=config['time_dim'])
                          bc_a4 = xr.concat([bc_a4, new_bc_a4], dim=config['time_dim'])
                          dst_a1 = xr.concat([dst_a1, new_dst_a1], dim=config['time_dim'])
                          dst_a3 = xr.concat([dst_a3, new_dst_a3], dim=config['time_dim'])
                          mom_a1 = xr.concat([mom_a1, new_mom_a1], dim=config['time_dim'])
                          mom_a2 = xr.concat([mom_a2, new_mom_a2], dim=config['time_dim'])
                          mom_a3 = xr.concat([mom_a3, new_mom_a3], dim=config['time_dim'])
                          mom_a4 = xr.concat([mom_a4, new_mom_a4], dim=config['time_dim'])
                          ncl_a1 = xr.concat([ncl_a1, new_ncl_a1], dim=config['time_dim'])
                          ncl_a2 = xr.concat([ncl_a2, new_ncl_a2], dim=config['time_dim'])
                          ncl_a3 = xr.concat([ncl_a3, new_ncl_a3], dim=config['time_dim'])
                          pom_a1 = xr.concat([pom_a1, new_pom_a1], dim=config['time_dim'])
                          pom_a3 = xr.concat([pom_a3, new_pom_a3], dim=config['time_dim'])
                          pom_a4 = xr.concat([pom_a4, new_pom_a4], dim=config['time_dim'])
                          so4_a1 = xr.concat([so4_a1, new_so4_a1], dim=config['time_dim'])
                          so4_a2 = xr.concat([so4_a2, new_so4_a2], dim=config['time_dim'])
                          so4_a3 = xr.concat([so4_a3, new_so4_a3], dim=config['time_dim'])
                          soa_a1 = xr.concat([soa_a1, new_soa_a1], dim=config['time_dim'])
                          soa_a2 = xr.concat([soa_a2, new_soa_a2], dim=config['time_dim'])
                          soa_a3 = xr.concat([soa_a3, new_soa_a3], dim=config['time_dim'])
                      else:
                          new_bc_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_bc_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_bc_a4 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_dst_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_dst_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_mom_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_mom_a2 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_mom_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_mom_a4 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_ncl_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_ncl_a2 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_ncl_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_pom_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_pom_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_pom_a4 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_so4_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_so4_a2 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_so4_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_soa_a1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_soa_a2 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          new_soa_a3 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                          bc_a1 = xr.concat([bc_a1, new_bc_a1], dim=config['time_dim'])
                          bc_a3 = xr.concat([bc_a3, new_bc_a3], dim=config['time_dim'])
                          bc_a4 = xr.concat([bc_a4, new_bc_a4], dim=config['time_dim'])
                          dst_a1 = xr.concat([dst_a1, new_dst_a1], dim=config['time_dim'])
                          dst_a3 = xr.concat([dst_a3, new_dst_a3], dim=config['time_dim'])
                          mom_a1 = xr.concat([mom_a1, new_mom_a1], dim=config['time_dim'])
                          mom_a2 = xr.concat([mom_a2, new_mom_a2], dim=config['time_dim'])
                          mom_a3 = xr.concat([mom_a3, new_mom_a3], dim=config['time_dim'])
                          mom_a4 = xr.concat([mom_a4, new_mom_a4], dim=config['time_dim'])
                          ncl_a1 = xr.concat([ncl_a1, new_ncl_a1], dim=config['time_dim'])
                          ncl_a2 = xr.concat([ncl_a2, new_ncl_a2], dim=config['time_dim'])
                          ncl_a3 = xr.concat([ncl_a3, new_ncl_a3], dim=config['time_dim'])
                          pom_a1 = xr.concat([pom_a1, new_pom_a1], dim=config['time_dim'])
                          pom_a3 = xr.concat([pom_a3, new_pom_a3], dim=config['time_dim'])
                          pom_a4 = xr.concat([pom_a4, new_pom_a4], dim=config['time_dim'])
                          so4_a1 = xr.concat([so4_a1, new_so4_a1], dim=config['time_dim'])
                          so4_a2 = xr.concat([so4_a2, new_so4_a2], dim=config['time_dim'])
                          so4_a3 = xr.concat([so4_a3, new_so4_a3], dim=config['time_dim'])
                          soa_a1 = xr.concat([soa_a1, new_soa_a1], dim=config['time_dim'])
                          soa_a2 = xr.concat([soa_a2, new_soa_a2], dim=config['time_dim'])
                          soa_a3 = xr.concat([soa_a3, new_soa_a3], dim=config['time_dim'])

                if config['ccn_output'] == True:
                    print(lst3d_wetaer[file_inds[ii+1]])
                    e3smdata3d_wetaer = xr.open_dataset(lst3d_wetaer[file_inds[ii+1]])
                    e3smdata3d_wetaer = e3smdata3d_wetaer.transpose(config['time_dim'],config['vert_dim'],config['latlon_dim']+E3SMdomain_range,...) # ensure ordering of time, height, and location
                    
                    vlist = list(e3smdata3d_wetaer.variables.keys())
                    av_vars = fnmatch.filter(vlist,'*'+E3SMdomain_range)
                  
                    req_vlist = [config['CCN1'], config['CCN2']]#, config['CCN5'], config['CCN10']]
                    req_vlist = ["{}{}".format(i,E3SMdomain_range) for i in req_vlist]
                    matched_vlist = list(set(av_vars).intersection(req_vlist))
                    if len(matched_vlist) == len(req_vlist):
                        print('\nAnalyzing CCN')
                        new_ccn1 = e3smdata3d_wetaer[config['CCN1']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        new_ccn2 = e3smdata3d_wetaer[config['CCN2']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        # new_ccn5 = e3smdata3d_wetaer[config['CCN5']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        # new_ccn10 = e3smdata3d_wetaer[config['CCN10']+E3SMdomain_range][:,:,latlon_ind,...].load()
                        ccn1 = xr.concat([ccn1, new_ccn1], dim=config['time_dim'])
                        ccn2 = xr.concat([ccn2, new_ccn2], dim=config['time_dim'])
                        # ccn5 = xr.concat([ccn5, new_ccn5], dim=config['time_dim'])
                        # ccn10 = xr.concat([ccn10, new_ccn10], dim=config['time_dim'])
                    else:
                        new_ccn1 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        new_ccn2 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        # new_ccn5 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        # new_ccn10 = xr.DataArray(np.zeros(newz3.shape)*np.nan,attrs={'units':'dummy_unit','long_name':'Dummy'})
                        ccn1 = xr.concat([ccn1, new_ccn1], dim=config['time_dim'])
                        ccn2 = xr.concat([ccn2, new_ccn2], dim=config['time_dim'])
                        # ccn5 = xr.concat([ccn5, new_ccn5], dim=config['time_dim'])
                        # ccn10 = xr.concat([ccn10, new_ccn10], dim=config['time_dim'])
        
        # # aerosol composition
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
        
        # # droplet size distribution
        # nd_cld = e3smdata['ICWNC'+'_'+E3SMdomain_range].load()
        # lmda = e3smdata['lambda_cloud'+'_'+E3SMdomain_range].load()
        # mu = e3smdata['mu_cloud'+'_'+E3SMdomain_range].load()
        
        # # other variables
        # for varname in variable3d_names:
        #     var = e3smdata[varname+'_'+E3SMdomain_range].load()
        #     variables.append(var)
        # e3smdata.close()

             
        #%% find the flight track grid
        print('\nInterpolating to Flight Track')
        for tt in range(len(time_new)):
            t_idx = np.abs(timem-time_new[tt]).argmin()
            x_idx = find_nearest(lonm, latm, lon_new[tt], lat_new[tt])
            z_idx = np.abs(z3[t_idx, :, x_idx]-height_new[tt]).argmin()
            for vv in range(len(variable3d_names)):
                variables_new[vv].append(float(variables[vv][t_idx, z_idx, x_idx]))
            p.append(Pres[t_idx, z_idx, x_idx].data)

            #Cloud LWC and IWC
            rho = np.array(Pres[t_idx, z_idx, x_idx].data)/T[t_idx, z_idx, x_idx].data/287.06
            if qc.attrs['units'] == 'kg/kg':
                cwc.append(qc[t_idx, z_idx, x_idx].data * rho * 1000)
            if qc.attrs['units'] == 'kg/m3':
                cwc.append(qc[t_idx, z_idx, x_idx].data * 1000)
            # if qi.attrs['units'] == 'kg/kg':
            #     iwc.append(qi[t_idx, z_idx, x_idx].data * rho * 1000)
            # if qi.attrs['units'] == 'kg/m3':
            #     iwc.append(qi[t_idx, z_idx, x_idx].data * 1000  )     

            #Rain Water Content
            if config['rain_output'] == True:
                if qr.attrs['units'] == 'kg/kg':
                    rwc.append(qr[t_idx, z_idx, x_idx].data * rho * 1000)
                if qr.attrs['units'] == 'kg/m3':
                    rwc.append(qr[t_idx, z_idx, x_idx].data * 1000)
                    
            if config['aerosol_output'] == True:
                # calculate aerosol size
                numall = [num_a1[t_idx, z_idx, x_idx].data, num_a2[t_idx, z_idx, x_idx].data, 
                          num_a3[t_idx, z_idx, x_idx].data, num_a4[t_idx, z_idx, x_idx].data]
                dnall  = [dn1[t_idx, z_idx, x_idx].data,    dn2[t_idx, z_idx, x_idx].data,    
                          dn3[t_idx, z_idx, x_idx].data,    dn4[t_idx, z_idx, x_idx].data]
                NCN = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[t_idx, z_idx, x_idx].data, Pres[t_idx, z_idx, x_idx].data)
                NCNall = np.hstack((NCNall, np.reshape(NCN,(3000,1))))
                if config['composition_output'] == True:
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

            if config['ccn_output'] == True:
                ccn1_all.append(ccn1[t_idx, z_idx, x_idx].data)
                ccn2_all.append(ccn2[t_idx, z_idx, x_idx].data)
                # ccn5_all.append(ccn5[t_idx, z_idx, x_idx].data)
                # ccn10_all.append(ccn10[t_idx, z_idx, x_idx].data)

            if config['dsd_output'] == True:
                # calculate droplet size distribution
                N0 = nd_cld[t_idx, z_idx, x_idx].data * (lmda[t_idx, z_idx, x_idx].data ** (mu[t_idx, z_idx, x_idx].data+1)) / \
                        gamma(mu[t_idx, z_idx, x_idx].data+1)    # parameter N0
                D_cld = np.arange(1, 1000) * 1e-6  # in m
                phi = N0 * (D_cld**mu[t_idx, z_idx, x_idx].data) * np.exp(- lmda[t_idx, z_idx, x_idx].data * D_cld)
                phi_all = np.hstack((phi_all, np.reshape(phi,(len(D_cld),1))))
                nd_bin = phi_all * (D_cld[1] - D_cld[0])   # droplet number concentration in each size bin
            
        if config['aerosol_output'] == True:
            NCN3 = np.nansum(NCNall[3:, :], 0)   # >3nm
            NCN10 = np.nansum(NCNall[10:, :], 0)    # >10nm
            NCN100 = np.nansum(NCNall[100:, :], 0)    # >100nm
        
      
        #%% change some units
        if config['aerosol_output'] == True:
            # composition
            T = variables_new[variable3d_names.index('T')]
            rho = np.array(p)/T/287.06
            # aerosol number
            NCNall = NCNall * 1e-6
            NCN3 = NCN3 * 1e-6
            NCN10 = NCN10 * 1e-6
            NCN100 = NCN100 * 1e-6
            ncn_units = '#/cm3'
        if config['composition_output'] == True:
            bc_all = np.array(bc_all)*1e9*rho
            dst_all = np.array(dst_all)*1e9*rho
            mom_all = np.array(mom_all)*1e9*rho
            ncl_all = np.array(ncl_all)*1e9*rho
            pom_all = np.array(pom_all)*1e9*rho
            so4_all = np.array(so4_all)*1e9*rho
            soa_all = np.array(soa_all)*1e9*rho
            composition_units = 'ug/m3'
        if config['dsd_output'] == True:
            # cloud number size distribution
            nd_bin = nd_bin * 1e-6
            nd_units = '#/cm3'
            
        # # LWC and IWC
        # idx = variable3d_names.index('LWC')
        # variables_new[idx] = np.array(variables_new[idx])*1000
        # variables[idx].attrs['units']='g/m3'
        # idx = variable3d_names.index('IWC')
        # variables_new[idx] = np.array(variables_new[idx])*1000
        # variables[idx].attrs['units']='g/m3'
        # # droplet number
        # idx = variable3d_names.index('ICWNC')
        # variables_new[idx] = np.array(variables_new[idx])*1e-6
        # variables[idx].attrs['units']='#/cm3'
        # idx = variable3d_names.index('ICINC')
        # variables_new[idx] = np.array(variables_new[idx])*1e-6
        # variables[idx].attrs['units']='#/cm3'
        
        #%% output  
        outfile = output_path + output_filehead + '_flight_'+date+'.nc'
        print('output file '+outfile)
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        t = f.createDimension('time', None)  # unlimited
        if config['aerosol_output'] == True:
            s = f.createDimension('CNsize', 3000)  # unlimited
        if config['dsd_output'] == True:
            s = f.createDimension('Ndsize', 999)  # unlimited
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time",))
        height_o = f.createVariable("height", 'f8', ("time",))
        var_o = list()
        for vv in range(len(variable3d_names)):
            var_o.append (f.createVariable(variable3d_names[vv], 'f8', ("time", )))
        p_o = f.createVariable('pres', 'f8', ("time",))
        cwc_o = f.createVariable('cwc', 'f8', ("time",))
        # iwc_o = f.createVariable('iwc', 'f8', ("time",))
        if config['rain_output'] == True:
            rwc_o = f.createVariable('rwc', 'f8', ("time",))
        if config['aerosol_output'] == True:
            ncn_o = f.createVariable('NCNall', 'f8', ("CNsize", "time",))
            ncn3_o = f.createVariable('NCN3', 'f8', ("time",))
            ncn10_o = f.createVariable('NCN10', 'f8', ("time",))
            ncn100_o = f.createVariable('NCN100', 'f8', ("time",))
        if config['composition_output'] == True:
            bc_o = f.createVariable('bc', 'f8', ("time",))
            dst_o = f.createVariable('dst', 'f8', ("time",))
            pom_o = f.createVariable('pom', 'f8', ("time",))
            mom_o = f.createVariable('mom', 'f8', ("time",))
            ncl_o = f.createVariable('ncl', 'f8', ("time",))
            so4_o = f.createVariable('so4', 'f8', ("time",))
            soa_o = f.createVariable('soa', 'f8', ("time",))
        if config['ccn_output'] == True:
            ccn1_o = f.createVariable('CCN1', 'f8', ("time",))
            ccn2_o = f.createVariable('CCN2', 'f8', ("time",))
            # ccn5_o = f.createVariable('CCN5', 'f8', ("time",))
            # ccn10_o = f.createVariable('CCN10', 'f8', ("time",))
        if config['dsd_output'] == True:
            nd_o = f.createVariable('Nd_bin', 'f8', ("Ndsize", "time",))
        
        # write data
        time_o[:] = time_new
        height_o[:] = height_new
        for vv in range(len(variable3d_names)):
            var_o[vv][:] = np.array(variables_new[vv])
        p_o[:] = np.array(p)
        cwc_o[:] = cwc
        # iwc_o[:] = iwc
        if config['rain_output'] == True:
            rwc_o[:] = rwc
        if config['aerosol_output'] == True:
            ncn_o[:,:] = NCNall
            ncn3_o[:] = NCN3
            ncn10_o[:] = NCN10
            ncn100_o[:] = NCN100
        if config['composition_output'] == True:
            bc_o[:] = np.array(bc_all)
            dst_o[:] = np.array(dst_all)
            pom_o[:] = np.array(pom_all)
            mom_o[:] = np.array(mom_all)
            ncl_o[:] = np.array(ncl_all)
            so4_o[:] = np.array(so4_all)
            soa_o[:] = np.array(soa_all)
        if config['ccn_output'] == True:
            ccn1_o[:] = ccn1_all
            ccn2_o[:] = ccn2_all
            # ccn5_o[:] = ccn5_all
            # ccn10_o[:] = ccn10_all
        if config['dsd_output'] == True:
            nd_o[:,:] = nd_bin
        
        # attributes
        time_o.units = "Seconds since " + timestr + ' 00:00 UTC'
        height_o.units = 'm MSL'
        for vv in range(len(variable3d_names)):
            var_o[vv].units = variables[vv].units
            var_o[vv].long_name = variables[vv].long_name
        p_o.units = 'Pa'
        p_o.long_name = 'Pressure'
        cwc_o.units = 'g/m3'
        cwc_o.long_name = 'Cloud Water Cotent'
        # iwc_o.units = 'g/m3'
        # iwc_o.long_name = 'Ice Water Cotent'
        if config['rain_output'] == True:
            rwc_o.units = 'g/m3'
            rwc_o.long_name = 'Rain Water Cotent'
        if config['aerosol_output'] == True:
            ncn_o.units = ncn_units
            ncn_o.long_name = 'Aerosol number size distribution'
            ncn_o.description = 'calculated from modal information into 1nm increment'
            ncn3_o.units = ncn_units
            ncn3_o.long_name = 'Aerosol number concentration for size >3nm'
            ncn10_o.units = ncn_units
            ncn10_o.long_name = 'Aerosol number concentration for size >10nm'
            ncn100_o.units = ncn_units
            ncn100_o.long_name = 'Aerosol number concentration for size >100nm'
        if config['composition_output'] == True:
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
        if config['ccn_output'] == True:
            ccn1_o.units = ccn1.units
            ccn1_o.long_name = ccn1.long_name
            ccn2_o.units = ccn2.units
            ccn2_o.long_name = ccn2.long_name
            # ccn5_o.units = ccn5.units
            # ccn5_o.long_name = ccn5.long_name
            # ccn10_o.units = ccn10.units
            # ccn10_o.long_name = ccn10.long_name
        if config['dsd_output'] == True:
            nd_o.units = nd_units
            nd_o.long_name = 'cloud droplet number size distribution'
            nd_o.description = 'calculated from microphysics scheme output into 1um increment from 1um to 1000um'
        
        # global attributes
        f.title = 'preprocessed E3SM data along aircraft track at the nearest time, grid, and vertical level'
        f.aircraftfile = filename.split('/')[-1]
        f.modelfile = lst3d[0].split('/')[-1]
        f.date = ttt.ctime(ttt.time())
        
        f.close()
    
