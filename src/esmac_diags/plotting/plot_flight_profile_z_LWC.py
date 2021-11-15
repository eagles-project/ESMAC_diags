"""
# plot vertical profile of cloud fraction for all flights in each IOP
# compare models and aircraft measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_wcm, read_RF_NCAR
from ..subroutines.read_netcdf import read_extractflight
from ..subroutines.quality_control import qc_mask_qcflag,qc_remove_neg

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
        wcmpath = settings.get('wcmpath', None)
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
    
    heightall=[]
    lwcobsall=[]
    lwcmall=[]
    
    nmodels=len(Model_List)
    for mm in range(nmodels):
        lwcmall.append([])
        
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
            
        #%% read in obs
        if campaign in ['HISCALE', 'ACEENA']:
            if date[-1]=='a':
                flightidx=1
            else:
                flightidx=2
            
            filename_wcm = glob.glob(wcmpath+'WCM_G1_'+date[0:8]+'*')
            filename_wcm.sort()
            if len(filename_wcm)==0:
                print('skip this date: '+date)
                continue
            (wcm,wcmlist)=read_wcm(filename_wcm[flightidx-1])
            time0=wcm[0,:]
            flag=wcm[-1,:]
            lwcobs=wcm[2,:]
            lwcobs=qc_remove_neg(lwcobs)
            lwcobs=qc_mask_qcflag(lwcobs,flag)
        
        elif campaign in ['CSET', 'SOCRATES']:
            filename = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
            if len(filename)==1 or len(filename)==2:  # SOCRATES has two flights in 20180217, choose the later one
                (time,lwcobs,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename[-1],'PLWCC')
            lwcobs=qc_remove_neg(lwcobs)
            
        lwcobsall.append(lwcobs)
        
        #%% read in models
        
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            
            (timem,heightm,lwc,timeunit,cldunit,cldname)=read_extractflight(filename_m,'LWC')
                
             # change E3SM unit from kg/m3 to g/m3 
            lwcmall[mm].append(lwc*1000)
            
        heightall.append(heightm)
        
    #%% calculate percentiles for each height bin
    
    lwcobs_z = list()
    lwcm_z = []
    for mm in range(nmodels):
        lwcm_z.append([])
    for zz in range(zlen):
        lwcobs_z.append(np.empty(0))
        for mm in range(nmodels):
            lwcm_z[mm].append(np.empty(0))
        
    ndays=len(heightall)
    # ndays=1;
    for dd in range(ndays):
        height = heightall[dd]
        lwcobs  = lwcobsall[dd]
        for zz in range(zlen):
            idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
            lwcobs_z[zz]=np.append(lwcobs_z[zz],lwcobs[idx])
            
        for mm in range(nmodels):
            lwcm = lwcmall[mm][dd]
            for zz in range(zlen):
                idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
                lwcm_z[mm][zz]=np.append(lwcm_z[mm][zz],lwcm[idx])
          
    #%% remove all NANs and calculate cloud frequency
    lwcmean_o = np.full(zlen,np.nan)
    std_lwc_o = np.full(zlen,np.nan)
    lwcmean_m = []
    for mm in range(nmodels):
        lwcmean_m.append(np.full(zlen,np.nan))
        
    for zz in range(zlen):
        data = lwcobs_z[zz]
        data = data[~np.isnan(data)]
        if len(data)>0:
            lwcmean_o[zz] = np.mean(data)
            std_lwc_o[zz] = np.std(data)/np.sqrt(len(data))
        for mm in range(nmodels):
            data = lwcm_z[mm][zz]
            data = data[~np.isnan(data)]
            if len(data)>0:
                lwcmean_m[mm][zz] = np.mean(data)
                
    #%% plot frequency  
    if campaign in ['HISCALE', 'ACEENA']:
        figname = figpath_aircraft_statistics+'profile_height_LWC_'+campaign+'_'+IOP+'.png'
    else:
        figname = figpath_aircraft_statistics+'profile_height_LWC_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,ax = plt.subplots(figsize=(3,8))
    
    ax.plot(lwcmean_o,z,color='k',linewidth=1,linestyle='-',label='Obs')
    ax.fill_betweenx(z,lwcmean_o-std_lwc_o,lwcmean_o+std_lwc_o,facecolor='k',alpha=0.2)
    
    for mm in range(nmodels):
        ax.plot(lwcmean_m[mm],z,color=color_model[mm],linewidth=1,label=Model_List[mm])
    
    ax.tick_params(color='k',labelsize=16)
    # ax.set_ylim(-1,zlen)
    # ax.set_yticks(range(zlen))
    if campaign=='HISCALE':
        ax.set_ylim(0,4500)
    ax.set_yticks(z)
    ax.set_ylabel('Height (m MSL)',fontsize=16)
    ax.legend(loc='upper right', fontsize='large')
    ax.set_xlabel('LWC (g/m3)',fontsize=16)
    if campaign in ['HISCALE', 'ACEENA']:
        ax.set_title(IOP,fontsize=18)
    else:
        ax.set_title(campaign,fontsize=18)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
                
                
