"""
# plot surface timeseries of aerosol number concentration
# compare models and surface measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import yyyymmdd2cday,cday2mmdd
from ..subroutines.read_ARMdata import read_cpc
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.specific_data_treatment import  avg_time_1d
from ..subroutines.quality_control import qc_remove_neg,qc_mask_qcflag_cpc

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    cpcsfcpath = settings['cpcsfcpath']
    cpcusfcpath = settings['cpcusfcpath']
    start_date = settings['start_date']
    end_date = settings['end_date']
    E3SM_sfc_path = settings['E3SM_sfc_path']
    figpath_sfc_timeseries = settings['figpath_sfc_timeseries']

    IOP = settings.get('IOP', None)

    #%% other settings
    
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date,'noleap')
    cday2 = yyyymmdd2cday(end_date,'noleap')
    if start_date[0:4]!=end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
        
    if not os.path.exists(figpath_sfc_timeseries):
        os.makedirs(figpath_sfc_timeseries)
        
    #%% read in obs data
    if campaign=='ACEENA':
        # cpc
        if IOP=='IOP1':
            lst = glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.2017062*')+glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201707*')
        elif IOP=='IOP2':
            lst = glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201801*')+glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201802*')
        lst.sort()
        t_cpc=np.empty(0)
        cpc=np.empty(0)
        for filename in lst:
            (time,data,qc,timeunit,cpcunit)=read_cpc(filename)
            data=qc_mask_qcflag_cpc(data,qc)
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date,'noleap')
            # average in time for quicker and clearer plot
            time2=np.arange(1800,86400,3600)
            data2 = avg_time_1d(np.array(time),np.array(data),time2)
            t_cpc=np.hstack((t_cpc, cday+time2/86400))
            cpc=np.hstack((cpc, data2))
        # no cpcu
        t_cpcu = np.array(np.nan)
        cpcu = np.array(np.nan)
        
    elif campaign=='HISCALE':  
        # cpc
        if IOP=='IOP1':
            lst = glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201604*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201605*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201606*')
        elif IOP=='IOP2':
            lst = glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201608*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201609*')
        lst.sort()
        t_cpc=np.empty(0)
        cpc=np.empty(0)
        if len(lst)==0:
            t_cpc = np.array(np.nan)
            cpc = np.array(np.nan)
        else:
            for filename in lst:
                (time,data,qc,timeunit,cpcunit)=read_cpc(filename)
                data=qc_mask_qcflag_cpc(data,qc)
                timestr=timeunit.split(' ')
                date=timestr[2]
                cday=yyyymmdd2cday(date,'noleap')
                t_cpc=np.hstack((t_cpc, cday+time/86400))
                cpc=np.hstack((cpc, data))
            cpc=qc_remove_neg(cpc)
      
        # cpcu
        if IOP=='IOP1':
            lst = glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201604*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201605*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201606*')
        elif IOP=='IOP2':
            lst = glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201608*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201609*')
        lst.sort()
        t_cpcu=np.empty(0)
        cpcu=np.empty(0)
        if len(lst)==0:
            t_cpcu = np.array(np.nan)
            cpcu = np.array(np.nan)
        else:
            for filename in lst:
                (time,data,qc,timeunit,cpcuunit)=read_cpc(filename)
                data=qc_mask_qcflag_cpc(data,qc)
                timestr=timeunit.split(' ')
                date=timestr[2]
                cday=yyyymmdd2cday(date,'noleap')
                t_cpcu=np.hstack((t_cpcu, cday+time/86400))
                cpcu=np.hstack((cpcu, data))
            cpcu=qc_remove_neg(cpcu)
        
    #%% read in models
    ncn_m = []
    nucn_m = []
    nmodels = len(Model_List)
    for mm in range(nmodels):
        tmp_ncn=np.empty(0)
        tmp_nucn=np.empty(0)
        timem=np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            
            filename_input = E3SM_sfc_path+'SFC_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (time,ncn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NCN')
            (time,nucn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NUCN')
            
            timem = np.hstack((timem,time))
            tmp_ncn = np.hstack((tmp_ncn,ncn*1e-6))
            tmp_nucn = np.hstack((tmp_nucn,nucn*1e-6))
        
        ncn_m.append(tmp_ncn)
        nucn_m.append(tmp_nucn)
        
    
    #%% make plot
        
    figname = figpath_sfc_timeseries+'timeseries_CN_'+campaign+'_'+IOP+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax1.plot(t_cpc,cpc,color='k',linewidth=1,label='CPC(>10nm)')
    for mm in range(nmodels):
        ax1.plot(timem, ncn_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax1.set_yscale('log')
    ax1.tick_params(color='k',labelsize=12)
    # ylim1 = ax1.get_ylim()
    
    ax2.plot(t_cpcu,cpcu,color='k',linewidth=1,label='CPC(>3nm)')
    for mm in range(nmodels):
        ax2.plot(timem, nucn_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax2.set_yscale('log')
    ax2.tick_params(color='k',labelsize=12)
    # ylim2 = ax2.get_ylim()
    
    # ax1.set_yticks([10,100,1000,10000,100000])
    # ax2.set_yticks([10,100,1000,10000,100000])
    ax1.set_xlim(cday1,cday2)
    ax2.set_xlim(cday1,cday2)
    
    # # set ylimit consistent in subplots
    # ax1.set_ylim([ylim1[0], ylim2[1]])
    # ax2.set_ylim([ylim1[0], ylim2[1]])
    
    
    ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    
    ax2.set_xlabel('Calendar Day',fontsize=14)
    ax1.set_title('Aerosol Number Concentration (cm$^{-3}$) '+campaign+' '+IOP,fontsize=15)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
