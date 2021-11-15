"""
# plot aircraft track data
# timeseries of aerosol size distribution
# compare models and aircraft measurements
"""


import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_RF_NCAR
from ..subroutines.specific_data_treatment import lwc2cflag, avg_time_2d
# from time_format_change import yyyymmdd2cday, hhmmss2sec
from ..subroutines.read_netcdf import read_merged_size,read_extractflight

def run_plot(settings):
    
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_timeseries = settings['figpath_aircraft_timeseries']

    IOP = settings.get('IOP', None)
    merged_size_path = settings.get('merged_size_path', None)
    RFpath = settings.get('RFpath', None)

    #%% other settings
    if not os.path.exists(figpath_aircraft_timeseries):
        os.makedirs(figpath_aircraft_timeseries)
        
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
        
    # dN/dlnDp for model
    dlnDp_m = np.empty((3000))
    for bb in range(3000):
        dlnDp_m[bb]=np.log((bb+2)/(bb+1))
    
    for date in alldates:
        
        #%% read in Models
        nmodels=len(Model_List)
        data_m = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (timem,heightm,datam,timeunitm,datamunit,datamlongname)=read_extractflight(filename_m,'NCNall')
            datam=datam*1e-6    # #/m3 to #/cm3
            # average in time for quicker plot
            time2 = np.arange(timem[0],timem[-1],60)
            data2 = avg_time_2d(timem,datam.T,time2)
            datam = data2.T
            # change to dN/dlnDp
            for tt in range(len(time2)):
                datam[:,tt]=datam[:,tt]/dlnDp_m
            data_m.append(datam) 
            
        # timem = (np.array(timem)-int(timem[0]))*24
        timem = time2/3600.
        
        
        #%% read observation        
        if campaign in ['HISCALE', 'ACEENA']:

            if campaign=='HISCALE':
                filename = merged_size_path+'merged_bin_fims_pcasp_'+campaign+'_'+date+'.nc'
            elif campaign=='ACEENA':
                filename = merged_size_path+'merged_bin_fims_pcasp_opc_'+campaign+'_'+date+'.nc'
            #% read in flight information
            (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
            (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
            (time,size,height,timeunit,zunit,long_name)=read_merged_size(filename,'height')
            (time,size,sizeh,timeunit,dataunit,long_name)=read_merged_size(filename,'size_high')
            (time,size,sizel,timeunit,dataunit,long_name)=read_merged_size(filename,'size_low')
            (time,size,merge,timeunit,dataunit,long_name)=read_merged_size(filename,'size_distribution_merged')
            time=np.ma.compressed(time)
            size=size*1000.
        
        elif campaign in ['CSET', 'SOCRATES']:
            filename = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
            # cloud flag
            (time,lwc,timeunit,lwcunit,lwclongname,size,cellunit)=read_RF_NCAR(filename[-1],'PLWCC')
            # calculate cloud flag based on LWC
            cflag=lwc2cflag(lwc,lwcunit)
            if campaign=='CSET':
                (time,uhsas,timeunit,dataunit,long_name,size,cellunit)=read_RF_NCAR(filename[-1],'CUHSAS_RWOOU')
            elif campaign=='SOCRATES':
                # there are two variables: CUHSAS_CVIU and CUHSAS_LWII
                (time,uhsas,timeunit,dataunit,long_name,size,cellunit)=read_RF_NCAR(filename[-1],'CUHSAS_LWII')
            merge = uhsas[:,0,:]
            size=size*1000.
            sizeh = size
            sizel = np.hstack((2*size[0]-size[1],  size[0:-1]))
            
        # merge=merge.T
        # time=time/3600.
        ## average in time for quicker plot
        time2=np.arange(time[0],time[-1],60)
        data2 = avg_time_2d(time,merge,time2)
        merge = data2.T
        time=time2/3600.
    
        # change to dN/dlnDp
        for bb in range(len(size)):
            dlnDp=np.log(sizeh[bb]/sizel[bb])
            merge[bb,:]=merge[bb,:]/dlnDp
            
            
        
        #%% make plot
        
        figname = figpath_aircraft_timeseries+'AerosolSize_'+campaign+'_'+date+'.png'
        print('plotting figures to '+figname)
        
        #fig = plt.figure()
        fig,ax = plt.subplots(nmodels+1,1,figsize=(8,2*(nmodels+1)))   # figsize in inches
        plt.tight_layout(h_pad=1.1)   #pad=0.4, w_pad=0.5, h_pad=1.0
        plt.subplots_adjust(right=0.9,bottom=0.1)
        
        leveltick=[0.1,1,10,100,1000,10000]
        levellist=np.arange(np.log(leveltick[0]),11,.5)
        
        merge[merge<0.01]=0.01
        h1 = ax[0].contourf(time,size,np.log(merge),levellist,cmap=plt.get_cmap('jet'))
        
        d_mam=np.arange(1,3001)
        h2=[]
        for mm in range(nmodels):
            datam = data_m[mm]
            datam[datam<0.01]=0.01
            h_m = ax[mm+1].contourf(timem,d_mam,np.log(datam),levellist,cmap=plt.get_cmap('jet'))
            h2.append(h_m)
    
        # colorbar
        cax = plt.axes([0.95, 0.2, 0.02, 0.6])
        cbar=fig.colorbar(h2[0], cax=cax, ticks=np.log(leveltick))
        cbar.ax.set_yticklabels(leveltick, fontsize=14)
        
        # set axis
        for ii in range(nmodels+1):
            ax[ii].set_xlim(timem[0],timem[-1])
            ax[ii].set_yscale('log')
            ax[ii].set_ylim(3, 5000)
            ax[ii].set_yticks([10,100,1000])
            ax[ii].tick_params(color='k',labelsize=14)
            if ii==0:
                ax[ii].text(0.01, 0.94, 'OBS', fontsize=14,transform=ax[ii].transAxes, verticalalignment='top')
            else:
                ax[ii].text(0.01, 0.94, Model_List[ii-1], fontsize=14,transform=ax[ii].transAxes, verticalalignment='top')
            
        ax[1].set_ylabel('Diameter (nm)',fontsize=14)
        ax[0].set_title('Size Distribution (#/dlnDp, cm-3)',fontsize=15)
        ax[nmodels].set_xlabel('time (hour UTC) in '+date,fontsize=14)
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()
