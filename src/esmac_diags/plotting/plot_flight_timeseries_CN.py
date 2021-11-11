"""
# plot aircraft track data
# timeseries of aerosol number concentration (CN)
# compare models and CPC measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_cpc, read_RF_NCAR
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
from ..subroutines.quality_control import qc_cpc_air

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_timeseries = settings['figpath_aircraft_timeseries']
    
    if campaign in ['HISCALE', 'ACEENA']:
        IOP = settings.get('IOP', None)
        merged_size_path = settings.get('merged_size_path', None)
        cpcpath = settings.get('cpcpath', None)
    elif campaign in ['CSET', 'SOCRATES']:
        RFpath = settings.get('RFpath', None)
    else:
        raise ValueError('campaign name is not recognized: '+campaign)

    #%% other settings
    
    if not os.path.exists(figpath_aircraft_timeseries):
        os.makedirs(figpath_aircraft_timeseries)
       
    
    #%% find files for flight information
    if campaign in ['HISCALE', 'ACEENA']:
        lst = glob.glob(merged_size_path+'merged_bin_*'+campaign+'*.nc')
    elif campaign in ['CSET', 'SOCRATES']:
        lst = glob.glob(RFpath+'RF*.PNI.nc')
    else:
        raise ValueError('campaign name is not recognized: '+campaign)
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
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
    elif campaign=='ACEENA':
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
        
    # for each flight
    for filename in lst:
        
        #%% read in flight data (for HISCALE and ACEENA)
        if campaign in ['HISCALE', 'ACEENA']:
            # get date info:        
            date=filename[-12:-3]
            if date[-1]=='a':
                flightidx=1
            else:
                flightidx=2
        
            #% read in flight information
            (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
            time=np.ma.compressed(time)
            if campaign=='HISCALE':
                filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
            elif campaign=='ACEENA':
                filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')    
            filename_c.sort()
            # read in data
            if len(filename_c)==1 or len(filename_c)==2: # some days have two flights
                (cpc,cpclist)=read_cpc(filename_c[flightidx-1])
                if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                    cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1)
                time_cpc = cpc[0,:]
                cpc10 = cpc[1,:]
                cpc3 = cpc[2,:]
            elif len(filename_c)==0:
                time_cpc=time
                cpc10=np.nan*np.empty([len(time)])
                cpc3=np.nan*np.empty([len(time)])
            else:
                raise ValueError('find too many files')
            # some quality checks
            (cpc3,cpc10) = qc_cpc_air(cpc3,cpc10)
            
        #%% read in flight data (for CSET and SOCRATES)
        elif campaign in ['CSET', 'SOCRATES']:
            fname=filename.split('.')
            date=fname[-4]
            (time_cpc,cpc10,timeunit,cpc10unit,cpc10longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCN')
            if campaign=='CSET':
                (time_cpc,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_RWOOU')
            elif campaign=='SOCRATES':
                # there are two variables: CONCU100_CVIU and CONCU100_LWII
                (time_cpc,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_LWII')
            
        #%% read in Models
        nmodels=len(Model_List)
        cpc100_m = []
        cpc10_m = []
        cpc3_m = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,cpc_m,timeunitm,ncn_unit,ncn_longname)=read_extractflight(filename_m,'NCN')
            (timem,heightm,cpcu_m,timeunitm,ncnu_unit,ncnu_longname)=read_extractflight(filename_m,'NUCN')
            (timem,heightm,ncnall,timeunitm,ncnall_unit,ncnall_longname)=read_extractflight(filename_m,'NCNall')
            # if len(cpc_m)!=cpc.shape[1]:
            #     print('CPC and MAM have different dimensions! check')
            #     print(cpc.shape,cpc_m.shape)
            #     errors
            cpc100_m.append(np.sum(ncnall[100:,:],0)*1e-6) # #/m3 to #/cm3
            cpc10_m.append(cpc_m*1e-6) # #/m3 to #/cm3
            cpc3_m.append(cpcu_m*1e-6) # #/m3 to #/cm3
        
        timem2 = timem/3600
        
        #%% make plot
            
        figname = figpath_aircraft_timeseries+'CN_'+campaign+'_'+date+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
        ax1.plot(time_cpc/3600,cpc10,color='k',linewidth=1,label='CPC(>10nm)')
        for mm in range(nmodels):
            ax1.plot(timem2, cpc10_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        ax1.set_yscale('log')
        ax1.tick_params(color='k',labelsize=12)
        ylim1 = ax1.get_ylim()
        
        if campaign in ['HISCALE', 'ACEENA']:
            ax2.plot(time_cpc/3600,cpc3,color='k',linewidth=1,label='CPC(>3nm)')
            for mm in range(nmodels):
                ax2.plot(timem2, cpc3_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        elif campaign in ['CSET', 'SOCRATES']:
            ax2.plot(time_cpc/3600,uhsas100,color='k',linewidth=1,label='UHSAS(>100nm)')
            for mm in range(nmodels):
                ax2.plot(timem2, cpc100_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        ax2.set_yscale('log')
        ax2.tick_params(color='k',labelsize=12)
        ylim2 = ax2.get_ylim()
        
        # set ylimit consistent in subplots
        ax1.set_ylim([max(1,min(ylim1[0],ylim2[0])), max(ylim1[1],ylim2[1])])
        ax2.set_ylim([max(1,min(ylim1[0],ylim2[0])), max(ylim1[1],ylim2[1])])
        
        ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
        ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
        
        ax2.set_xlabel('time (hour UTC) '+date,fontsize=14)
        ax1.set_title('Aerosol Number Concentration (cm$^{-3}$)',fontsize=15)
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()