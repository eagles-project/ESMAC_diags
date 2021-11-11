"""
# plot vertical profile of cloud fraction
# for each day of selected IOP
# compare models and surface measurements
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import timeunit2cday,yyyymmdd2cday,cday2mmdd
from ..subroutines.read_ARMdata import read_armbe
from ..subroutines.read_netcdf import read_E3SM_z

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    lon0 = settings['lon0']
    Model_List = settings['Model_List']
    armbepath = settings['armbepath']
    start_date = settings['start_date']
    end_date = settings['end_date']
    E3SM_profile_path = settings['E3SM_profile_path']
    figpath_profile_timeseries = settings['figpath_profile_timeseries']

    IOP = settings.get('IOP', None)

    #%% other settings
        
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date,'noleap')
    cday2 = yyyymmdd2cday(end_date,'noleap')
    if start_date[0:4]!=end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
        
    if not os.path.exists(figpath_profile_timeseries):
        os.makedirs(figpath_profile_timeseries)
    
    
    #%% read in obs data
    if campaign=='ACEENA':
        if IOP=='IOP1':
            filename_armbe = armbepath+'enaarmbecldradC1.c1.20170101.003000.nc'
            year='2017'
        elif IOP=='IOP2':
            filename_armbe = armbepath+'enaarmbecldradC1.c1.20180101.003000.nc'
            year='2018'
    elif campaign=='HISCALE':  
        filename_armbe = armbepath+'sgparmbecldradC1.c1.20160101.003000.nc'
        year='2016'
            
    (time0,height0,cld0,time0unit,cld0unit) = read_armbe(filename_armbe,'cld_frac')
    
    
    time0=time0/86400.+timeunit2cday(time0unit)
    if campaign=='HISCALE':
        # observation is leap year. change the time for comparison with noleap model output. 
        # note that it is not suitable for January and February
        time0=time0-1   
        
    #%% read in model
    
    cldm = []
    nmodels = len(Model_List)
    for mm in range(nmodels):
        timem=np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            
            filename_input = E3SM_profile_path+'Profile_vars_'+campaign+'_'+Model_List[mm]+'.'+date+'.nc'
            (time,height,data,timemunit,dataunit,long_name)=read_E3SM_z(filename_input,'CLOUD')
            
            timem = np.hstack((timem,time))
            if cday==cday1:
                datam=data*100
            else:
                datam = np.vstack((datam,data*100))
        
        data=data*100.
        dataunit='%'
        cldm.append(datam)
    
    # change to local solar time
    timeshift = lon0/360*24
    if timeshift>12:
        timeshift=timeshift-24
    time0 = time0+timeshift/24.
    timem = timem+timeshift/24.
    
    #%% plot cloud for each day in time_range
    
    
    for cday in range(cday1,cday2+1):
        idxo=np.logical_and(time0>cday-0.1, time0<cday+1.1)
        idxm=np.logical_and(timem>cday-0.1, timem<cday+1.1)
        
        mmdd = cday2mmdd(cday)
        figname = figpath_profile_timeseries+'cloudfraction_'+campaign+'_'+year+mmdd+'.png'
        print('plotting figures to '+figname)
        
        fig,ax = plt.subplots(nmodels+1,1,figsize=(6,2*nmodels+1))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
        h0=ax[0].contourf(time0[idxo],height0/1000,cld0[idxo,:].T,np.arange(0,101,10),cmap=plt.get_cmap('jet'))
        ax[0].set_title('Cloud Fraction (%) Obs')
        ax[0].set_ylim(0,15)
        ax[0].set_xlim(cday,cday+1)
        ax[0].set_xticks(np.arange(cday,cday+1.01,0.125))
        ax[0].set_xticklabels([])
        ax[0].set_ylabel('Height (km)')
        for mm in range(nmodels):
            ax[mm+1].contourf(timem[idxm],height/1000,cldm[mm][idxm,:].T,np.arange(0,101,10),cmap=plt.get_cmap('jet'))
            ax[mm+1].set_title(Model_List[mm])
            ax[mm+1].set_ylim(0,15)
            ax[mm+1].set_xlim(cday,cday+1)
            ax[mm+1].set_xticks(np.arange(cday,cday+1.01,0.125))
            ax[mm+1].set_xticklabels([])
            ax[mm+1].set_ylabel('Height (km)')
            
        ax[-1].set_xticklabels(['00','03','06','09','12','15','18','21','24'])
        ax[-1].set_xlabel('Local Solar Time for ' + campaign +' day '+year+mmdd)
        
        cax = plt.axes([1.01, 0.2, 0.02, 0.6])
        cbar=fig.colorbar(h0, cax=cax)
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()
