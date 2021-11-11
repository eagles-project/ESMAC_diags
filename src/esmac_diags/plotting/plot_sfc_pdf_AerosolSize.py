"""
# plot mean aerosol size ditribution for surface data
# compare models and surface measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import yyyymmdd2cday, cday2mmdd
from ..subroutines.read_surface import read_smpsb_pnnl,read_smps_bin
from ..subroutines.read_ARMdata import read_uhsas, read_smps_bnl
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.specific_data_treatment import  avg_time_2d
from ..subroutines.quality_control import qc_mask_qcflag,qc_correction_nanosmps

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    start_date = settings['start_date']
    end_date = settings['end_date']
    E3SM_sfc_path = settings['E3SM_sfc_path']
    figpath_sfc_statistics = settings['figpath_sfc_statistics']

    if campaign=='ACEENA':
        IOP = settings['IOP']
        uhsassfcpath = settings['uhsassfcpath']
    elif campaign=='HISCALE':
        IOP = settings['IOP']
        if IOP=='IOP1':
            smps_bnl_path = settings['smps_bnl_path']
            nanosmps_bnl_path = settings['nanosmps_bnl_path']
        elif IOP=='IOP2':
            smps_pnnl_path = settings['smps_pnnl_path']
    else:
        raise ValueError('campaign name is not recognized: '+campaign)
        
    #%% other settings
        
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date,'noleap')
    cday2 = yyyymmdd2cday(end_date,'noleap')
    if start_date[0:4]!=end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
    
    if not os.path.exists(figpath_sfc_statistics):
        os.makedirs(figpath_sfc_statistics)
            
    #%% read in obs data
    if campaign=='ACEENA':
        if IOP=='IOP1':
            lst = glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.2017062*')+glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.201707*')
        elif IOP=='IOP2':
            lst = glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.201801*')+glob.glob(uhsassfcpath+'enaaosuhsasC1.a1.201802*')
        lst.sort()
        t_uhsas=np.empty(0)
        uhsas=np.empty((0,99))
        for filename in lst:
            (time,dmin,dmax,data,timeunit,dataunit,long_name) = read_uhsas(filename)
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date,'noleap')
            # average in time for quicker plot
            time2=np.arange(1800,86400,3600)
            data2 = avg_time_2d(time,data,time2)
            t_uhsas=np.hstack((t_uhsas, cday+time2/86400))
            uhsas=np.vstack((uhsas, data2))
        size_u = (dmin+dmax)/2
        # change to dN/dlogDp
        dlnDp_u=np.empty(99)
        for bb in range(len(size_u)):
            dlnDp_u[bb]=np.log10(dmax[bb]/dmin[bb])
            uhsas[:,bb]=uhsas[:,bb]/dlnDp_u[bb]
        
        time = np.array(t_uhsas)
        size = np.array(size_u)
        obs = np.array(uhsas.T)
        
    elif campaign=='HISCALE':    
        if IOP=='IOP1':
            lst = glob.glob(smps_bnl_path+'*.nc')
            lst.sort()
            t_smps=np.empty(0)
            smps=np.empty((0,192))
            for filename in lst:
                (time,size,flag,timeunit,dataunit,smps_longname)=read_smps_bnl(filename,'status_flag')
                (time,size,data,timeunit,smpsunit,smps_longname)=read_smps_bnl(filename,'number_size_distribution')
                data=qc_mask_qcflag(data,flag)
                timestr=timeunit.split(' ')
                date=timestr[2]
                cday=yyyymmdd2cday(date,'noleap')
                # average in time for quicker plot
                time2=np.arange(1800,86400,3600)
                data2 = avg_time_2d(time,data,time2)
                t_smps=np.hstack((t_smps, cday+time2/86400))
                smps=np.vstack((smps, data2))
            smps=smps.T
            # combine with nanoSMPS
            lst2 = glob.glob(nanosmps_bnl_path+'*.nc')
            lst2.sort()
            t_nano=np.empty(0)
            nanosmps=np.empty((0,192))
            for filename2 in lst2:
                (timen,sizen,flagn,timenunit,datanunit,long_name)=read_smps_bnl(filename2,'status_flag')
                (timen,sizen,datan,timenunit,nanounit,nanoname)=read_smps_bnl(filename2,'number_size_distribution')
                datan=qc_mask_qcflag(datan,flagn)
                timestr=timenunit.split(' ')
                date=timestr[2]
                cday=yyyymmdd2cday(date,'noleap')
                # average in time for quicker plot
                time2=np.arange(1800,86400,3600)
                data2 = avg_time_2d(timen,datan,time2)
                t_nano=np.hstack((t_nano, cday+time2/86400))
                nanosmps=np.vstack((nanosmps, data2))
            # nanosmps is overcounting, adjust nanosmps value for smooth transition to SMPS
            nanosmps=qc_correction_nanosmps(nanosmps.T)
            for tt in range(smps.shape[1]):
                if any(t_nano==t_smps[tt]):
                    smps[0:80,tt]=nanosmps[0:80,t_nano==t_smps[tt]].reshape(80)
            
        elif IOP=='IOP2':
            data=read_smpsb_pnnl(smps_pnnl_path+'HiScaleSMPSb_SGP_20160827_R1.ict')
            size=read_smps_bin(smps_pnnl_path+'NSD_column_size_chart.txt')
            time=data[0,:]
            smps=data[1:-1,:]
            flag=data[-1,:]
            smps=qc_mask_qcflag(smps.T,flag).T
            cday=yyyymmdd2cday('2016-08-27')
            # average in time for quicker plot
            time2=np.arange(time[0],time[-1]+1800,3600)
            data2 = avg_time_2d(time,smps.T,time2)
            t_smps=cday+time2/86400
            smps=data2.T
            
        time = np.array(t_smps)
        size = np.array(size)
        obs = np.array(smps)  
        
        # SMPS is already divided by log10
        
    else:
        raise ValueError('does not recognize this campaign: '+campaign)
        
    #%% read in models
    model = []
    nmodels = len(Model_List)
    for mm in range(nmodels):
        tmp_data=np.empty((3000,0))
        timem=np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            
            filename_input = E3SM_sfc_path+'SFC_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (time,ncn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NCNall')
            
            timem = np.hstack((timem,time))
            tmp_data = np.hstack((tmp_data,ncn*1e-6))
        
        # change to dN/dlog10Dp
        for bb in range(3000):
            dlnDp=np.log10((bb+2)/(bb+1))
            tmp_data[bb,:]=tmp_data[bb,:]/dlnDp
        
        model.append(tmp_data)
            
    #%% calculate mean pdf
    pdf_obs=np.nanmean(obs,1)
    pdf_model=[None]*nmodels
    for mm in range(nmodels):
        pdf_model[mm]=np.nanmean(model[mm],1)
    
    #%%
    pct1_o = [np.nanpercentile(obs[i,:],10) for i in range(len(size))]
    pct2_o = [np.nanpercentile(obs[i,:],90) for i in range(len(size))]
    pct1_m = [[] for mm in range(nmodels)]
    pct2_m = [[] for mm in range(nmodels)]
    for mm in range(nmodels):
        pct1_m[mm] = [np.nanpercentile(model[mm][i,:],10) for i in range(3000)]
        pct2_m[mm] = [np.nanpercentile(model[mm][i,:],90) for i in range(3000)]
    
    # import scipy.stats as stats
    # sem_o = np.ma.filled(stats.sem(obs,1,nan_policy='omit'),np.nan)
    # sem_m = [[] for mm in range(nmodels)]
    # for mm in range(nmodels):
    #     sem_m[mm] = np.ma.filled(stats.sem(model[mm],1,nan_policy='omit'),np.nan)
    
    #%% make plot
    # not plotting data if the mean value is 0
    pdf_obs[pdf_obs==0] = np.nan
    
    figname = figpath_sfc_statistics+'pdf_AerosolSize_'+campaign+'_'+IOP+'.png'
    
    print('plotting figures to '+figname)
    
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
    ax.set_title(campaign+' '+IOP,fontsize=14)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()
