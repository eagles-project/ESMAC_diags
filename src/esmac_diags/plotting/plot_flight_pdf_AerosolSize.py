"""
# plot mean aerosol size ditribution for aircraft track data
# average for each IOP
# compare models and aircraft measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_RF_NCAR
from ..subroutines.specific_data_treatment import lwc2cflag
# from time_format_change import yyyymmdd2cday, hhmmss2sec
from ..subroutines.read_netcdf import read_merged_size,read_extractflight

from ..subroutines.specific_data_treatment import  avg_time_2d
from ..subroutines.quality_control import qc_mask_cloudflag, qc_uhsas_RF_NCAR,qc_remove_neg,qc_mask_takeoff_landing

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    IOP = settings.get('IOP', None)
    merged_size_path = settings.get('merged_size_path', None)

    #%% other settings

    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
        
        
    #%% find files for flight information

    lst = sorted(glob.glob(E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[0]+'_*.nc'))
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
        
    print('reading '+format(len(alldates))+' files to calculate mean aerosol pdf: ')

    nmodels=len(Model_List)
    pdfall_m = [np.empty((3000,0)) for mm in range(nmodels)]
    size_m = np.zeros(3000)
    pdf_model = [size_m for mm in range(nmodels)]
    if 'pdf_obs' in locals():
        del pdf_obs

    # number of valid timesteps
    n_o = 0
    n_m = [0 for mm in range(nmodels)]
        

    # dN/dlnDp for model
    dlnDp_m = np.empty((3000))
    for bb in range(3000):
        dlnDp_m[bb]=np.log((bb+2)/(bb+1))

    for date in alldates[:]:
        print(date)
        
        #%% read in Models
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (timem,heightm,datam,timeunitm,datamunit,datamlongname)=read_extractflight(filename_m,'NCNall')
            datam=datam*1e-6    # #/m3 to #/cm3
            
            # average in time for quicker plot
            time2=np.arange(300,86400,600)
            data2 = avg_time_2d(timem,datam.T,time2)
            datam=data2.T
            timem=time2
            
            for tt in range(len(timem)):
                datam[:,tt]=datam[:,tt]/dlnDp_m
                
            pdfall_m[mm] = np.column_stack((pdfall_m[mm],datam))
            for tt in range(len(timem)):
                if ~np.isnan(datam[0,tt]):
                    pdf_model[mm] = pdf_model[mm]+datam[:,tt]
                    n_m[mm]=n_m[mm]+1
            
        #%% read observation        
        if campaign in ['HISCALE', 'ACEENA']:
            if date[-1]=='a':
                flightidx=1
            else:
                flightidx=2
                
            if campaign=='HISCALE':
                filename = merged_size_path+'merged_bin_fims_pcasp_'+campaign+'_'+date+'.nc'
            elif campaign=='ACEENA':
                filename = merged_size_path+'merged_bin_fims_pcasp_opc_'+campaign+'_'+date+'.nc'
        
            (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
            (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
            (time,size,legnum,timeunit,zunit,long_name)=read_merged_size(filename,'leg_number')
            (time,size,sizeh,timeunit,dataunit,long_name)=read_merged_size(filename,'size_high')
            (time,size,sizel,timeunit,dataunit,long_name)=read_merged_size(filename,'size_low')
            (time,size,merge,timeunit,dataunit,long_name)=read_merged_size(filename,'size_distribution_merged')
            time=np.ma.compressed(time)
            size=size*1000.
            merge = qc_mask_cloudflag(merge,cflag)
            
            # average in time for quicker plot
            time2=np.arange(300,86400,600)
            data2 = avg_time_2d(time,merge,time2)
            merge = data2.T
            time=time2/3600.
            

        elif campaign in ['CSET', 'SOCRATES']:
            filename = glob.glob(settings['RFpath']+'RF*'+date+'*.PNI.nc')
            # cloud flag
            (time,lwc,timeunit,lwcunit,lwclongname,size,cellunit)=read_RF_NCAR(filename[-1],'PLWCC')
            if campaign=='CSET':
                (time,uhsas,timeunit,dataunit,long_name,size,cellunit)=read_RF_NCAR(filename[-1],'CUHSAS_RWOOU')
            elif campaign=='SOCRATES':
                # there are two variables: CUHSAS_CVIU and CUHSAS_LWII
                (time,uhsas,timeunit,dataunit,long_name,size,cellunit)=read_RF_NCAR(filename[-1],'CUHSAS_LWII')
            uhsas=uhsas[:,0,:]
            # calculate cloud flag based on LWC
            cflag=lwc2cflag(lwc,lwcunit)
            uhsas = qc_mask_cloudflag(uhsas,cflag)
            uhsas= qc_uhsas_RF_NCAR(uhsas)
            
            # average in time for quicker plot
            time2=np.arange(300,86400,600)
            data2 = avg_time_2d(time,uhsas,time2)
            merge = data2.T
            time0 = np.array(time)
            time=time2/3600.
            
            size=size*1000.
            sizeh = size
            sizel = np.hstack((2*size[0]-size[1],  size[0:-1]))
        
        # change to dN/dlnDp
        for bb in range(len(size)):
            dlnDp=np.log(sizeh[bb]/sizel[bb])
            merge[bb,:]=merge[bb,:]/dlnDp
        
        merge=qc_remove_neg(merge)
        
        # exclude 30min after takeoff and before landing
        merge = qc_mask_takeoff_landing(time2,merge)
        
        # fig,ax=plt.subplots()
        # ax.plot(merge[9,:])
        # ax.set_title(date)
        # error
        
        if ('pdf_obs' in locals()) == False:
            pdf_obs = np.zeros(len(size)) 
            pdfall_o = np.empty((len(size),0))
        idx_valid = ~np.isnan(np.mean(merge,0))
        pdf_obs = pdf_obs + np.sum(merge[:,idx_valid],1)
        pdfall_o = np.hstack((pdfall_o,np.array(merge[:,idx_valid])))
        n_o = n_o + np.sum(idx_valid)
        

    #%% calculate mean pdf

    pdf_obs[pdf_obs<1e-3]=np.nan
    pdf_obs=pdf_obs/n_o
    for mm in range(nmodels):
        pdf_model[mm]=pdf_model[mm]/n_m[mm]

    #%%
    pdfall_o[pdfall_o<0]=np.nan
    pct1_o = [np.nanpercentile(pdfall_o[i,:],10) for i in range(len(size))]
    pct2_o = [np.nanpercentile(pdfall_o[i,:],90) for i in range(len(size))]
    pct1_m = [[] for mm in range(nmodels)]
    pct2_m = [[] for mm in range(nmodels)]
    for mm in range(nmodels):
        pct1_m[mm] = [np.nanpercentile(pdfall_m[mm][i,:],10) for i in range(3000)]
        pct2_m[mm] = [np.nanpercentile(pdfall_m[mm][i,:],90) for i in range(3000)]

    #%% make plot

    if campaign in ['HISCALE', 'ACEENA']:
        figname = figpath_aircraft_statistics+'pdf_AerosolSize_'+campaign+'_'+IOP+'.png'
    else:
        figname = figpath_aircraft_statistics+'pdf_AerosolSize_'+campaign+'.png'

    print('plotting figures to '+figname)

    #fig = plt.figure()
    fig,ax = plt.subplots(figsize=(4,2.5))   # figsize in inches

    ax.plot(size,pdf_obs,color='k',label='Obs')
    for mm in range(nmodels):
        ax.plot(np.arange(1,3001),pdf_model[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])

    ax.fill_between(size,pct1_o,pct2_o, alpha=0.5, facecolor='gray')
    for mm in range(nmodels):
        ax.fill_between(np.arange(1,3001),pct1_m[mm],pct2_m[mm], alpha=0.2, facecolor=color_model[mm])

    ax.legend(loc='upper right', shadow=False, fontsize='medium')
    ax.tick_params(color='k',labelsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.01,1e4)
    ax.set_xlim(0.67,4500)
    ax.set_xlabel('Diameter (nm)',fontsize=13)
    ax.set_ylabel('#/dlnDp (cm$^{-3}$)',fontsize=13)

    if campaign in ['HISCALE', 'ACEENA']:
        ax.set_title(campaign+' '+IOP,fontsize=14)
    else:
        ax.set_title(campaign,fontsize=14)

    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()

