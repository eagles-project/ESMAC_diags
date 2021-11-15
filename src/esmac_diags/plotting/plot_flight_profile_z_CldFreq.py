"""
# plot vertical profile of cloud fraction for all flights in each IOP
# compare models and aircraft measurements
"""

import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import  read_RF_NCAR
from ..subroutines.specific_data_treatment import lwc2cflag
from ..subroutines.read_netcdf import read_extractflight,read_merged_size

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    height_bin = settings['height_bin']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    if campaign in ['HISCALE', 'ACEENA']:
        IOP = settings.get('IOP', None)
        merged_size_path = settings.get('merged_size_path', None)
    elif campaign in ['CSET', 'SOCRATES']:
        RFpath = settings.get('RFpath', None)
    else:
        raise ValueError('campaign name is not recognized: '+campaign)

    #%% other settings
    
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
        
    #%%
    z=height_bin
    dz = z[1]-z[0]
    zmin=z-np.insert((z[1:]-z[0:-1])/2,0,dz)
    zmax=z+np.append((z[1:]-z[0:-1])/2,dz)
    
    zlen=len(z)   
    
    
    #%% find files for flight information
    
    lst = glob.glob(E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[0]+'_*.nc')
    lst.sort()
    if len(lst)==0:
        raise ValueError('cannot find any file')
    # choose files for specific IOP
    if campaign=='HISCALE':
        if IOP=='IOP1':
            lst=lst[0:17]
        elif IOP=='IOP2':
            lst=lst[17:]
        elif IOP[0:4]=='2016':
            a=lst[0].split('_'+Model_List[0]+'_')
            lst = glob.glob(a[0]+'_'+Model_List[0]+'_'+IOP+'*')
            lst.sort()
    elif campaign=='ACEENA':
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+Model_List[0]+'_')
            lst = glob.glob(a[0]+'_'+Model_List[0]+'_'+IOP+'*')
            lst.sort()
            
    alldates = [x.split('_')[-1].split('.')[0] for x in lst]
        
    #%% read all data
    
    cflagall=[]
    heightall=[]
    cldmall=[]
    
    nmodels=len(Model_List)
    for mm in range(nmodels):
        cldmall.append([])
        
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        
        #%% read in models
        
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            
            (timem,heightm,cloud,timeunit,cldunit,cldname)=read_extractflight(filename_m,'CLOUD')
                
            cldmall[mm].append(cloud)
            
        #%% read in obs
        if campaign in ['HISCALE', 'ACEENA']:

            #% read in flight information
            if campaign=='HISCALE':
                filename = merged_size_path+'merged_bin_fims_pcasp_'+campaign+'_'+date+'.nc'
            elif campaign=='ACEENA':
                filename = merged_size_path+'merged_bin_fims_pcasp_opc_'+campaign+'_'+date+'.nc'
            (time,size,height,timeunit,cunit,long_name)=read_merged_size(filename,'height')
            (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
            time=np.ma.compressed(time)
        
        elif campaign in ['CSET', 'SOCRATES']:
            filename = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
            if len(filename)==1 or len(filename)==2:  # SOCRATES has two flights in 20180217, choose the later one
                (time,lwc,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename[-1],'PLWCC')
            # calculate cloud flag based on LWC
            cflag=lwc2cflag(lwc,lwcunit)
            
        heightall.append(heightm)
        cflagall.append(cflag)
        
    
    #%% calculate percentiles for each height bin
    
    cflag_z = list()
    cldm_z = []
    nmodels=len(Model_List)
    for mm in range(nmodels):
        cldm_z.append([])
    for zz in range(zlen):
        cflag_z.append(np.empty(0))
        for mm in range(nmodels):
            cldm_z[mm].append(np.empty(0))
        
    ndays=len(heightall)
    # ndays=1;
    for dd in range(ndays):
        height = heightall[dd]
        cflag  = cflagall[dd]
        for zz in range(zlen):
            idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
            cflag_z[zz]=np.append(cflag_z[zz],cflag[idx])
            
        for mm in range(nmodels):
            cldm = cldmall[mm][dd]
            for zz in range(zlen):
                idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
                cldm_z[mm][zz]=np.append(cldm_z[mm][zz],cldm[idx])
            
    #%% remove all NANs and calculate cloud frequency
    cldfreq_flag = np.full(zlen,np.nan)
    cldfreq_m = []
    for mm in range(nmodels):
        cldfreq_m.append(np.full(zlen,np.nan))
        
    for zz in range(zlen):
        data = cflag_z[zz]
        data = data[data>=0]
        if len(data)>0:
            cldfreq_flag[zz] = sum(data==1)/len(data)
        for mm in range(nmodels):
            data = cldm_z[mm][zz]
            data = data[~np.isnan(data)]
            if len(data)>0:
                cldfreq_m[mm][zz] = np.mean(data)
      
    #%% plot frequency  
    if campaign in ['HISCALE', 'ACEENA']:
        figname = figpath_aircraft_statistics+'profile_height_CldFreq_'+campaign+'_'+IOP+'.png'
    else:
        figname = figpath_aircraft_statistics+'profile_height_CldFreq_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,ax = plt.subplots(figsize=(4,8))
    
    ax.plot(cldfreq_flag,z,color='k',linewidth=1,linestyle='-',label='Obs')
    for mm in range(nmodels):
        ax.plot(cldfreq_m[mm],z,color=color_model[mm],linewidth=1,label=Model_List[mm])
    
    ax.tick_params(color='k',labelsize=12)
    # ax.set_ylim(-1,zlen)
    # ax.set_yticks(range(zlen))
    # ax.set_yticks(z[0:-1:2])
    ax.set_ylabel('Height (m MSL)',fontsize=12)
    ax.legend(loc='upper right', fontsize='large')
    ax.set_xlabel('Cloud Frequency',fontsize=12)
    if campaign in ['HISCALE', 'ACEENA']:
        ax.set_title(IOP,fontsize=15)
    else:
        ax.set_title(campaign,fontsize=15)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)