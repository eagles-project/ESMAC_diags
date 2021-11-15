"""
# plot surface timeseries of CCN size distribution
# compare models and surface measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import yyyymmdd2cday,cday2mmdd
from ..subroutines.read_ARMdata import read_ccn
from ..subroutines.read_surface import read_CCN_hiscale_IOP1, read_CCN_hiscale_IOP2
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.quality_control import qc_remove_neg,qc_mask_qcflag,qc_ccn_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    ccnsfcpath = settings['ccnsfcpath']
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
            lst = glob.glob(ccnsfcpath+'enaaosccn1colavgC1.b1.201706*')+glob.glob(ccnsfcpath+'enaaosccn1colavgC1.b1.201707*')
        elif IOP=='IOP2':
            lst = glob.glob(ccnsfcpath+'enaaosccn1colavgC1.b1.201801*')+glob.glob(ccnsfcpath+'enaaosccn1colavgC1.b1.201802*')
        lst.sort()
        t_ccn=np.empty(0)
        ccn=np.empty(0)
        SS=np.empty(0)
        for filename in lst:
            (time,timeunit,data,qc,dataunit,SS0)=read_ccn(filename)
            data=qc_mask_qcflag(data,qc)
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date,'noleap')
            t_ccn=np.hstack((t_ccn, cday+time/86400))
            ccn=np.hstack((ccn, data))
            SS=np.hstack((SS, SS0))
        ccn=qc_remove_neg(ccn)
        ccn=qc_ccn_max(ccn,SS)
        # SS=0.1%
        idx = np.logical_and(SS>0.05, SS<0.15)
        t_ccna = t_ccn[idx]
        ccna = ccn[idx]
        SSa = 0.1
        # SS=0.5%
        idx = np.logical_and(SS>0.4, SS<0.6)
        t_ccnb = t_ccn[idx]
        ccnb = ccn[idx]
        SSb = 0.5
    
    elif campaign=='HISCALE':  
        if IOP=='IOP1':
            (times_ccn,ccnsfc,sssfc,timeunit)=read_CCN_hiscale_IOP1(ccnsfcpath)
            sssfc=[int(x*10) for x in sssfc]
            sssfc=np.array(sssfc)/10.
            times_ccn=np.array(times_ccn)
            ccnsfc=np.array(ccnsfc)
        elif IOP=='IOP2':
            (times_ccn,ccnsfc,sssfc,timeunit)=read_CCN_hiscale_IOP2(ccnsfcpath)
            sssfc=[int(x*10) for x in sssfc]
            sssfc=np.array(sssfc)/10.
            times_ccn=np.array(times_ccn)
            ccnsfc=np.array(ccnsfc)
        # find the nearest Supersaturation in Obs comparing to model
        # 0.1%
        idx = sssfc==0.1
        ccna = ccnsfc[idx]
        t_ccna = times_ccn[idx]
        SSa = 0.1
        # 0.5%
        idx = sssfc==0.5
        ccnb = ccnsfc[idx]
        t_ccnb = times_ccn[idx]
        SSb = 0.5
        
    #%% read in models
    ccna_m = []
    ccnb_m = []
    nmodels = len(Model_List)
    for mm in range(nmodels):
        tmp_CCN3=np.empty(0)
        tmp_CCN5=np.empty(0)
        timem=np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            
            filename_input = E3SM_sfc_path+'SFC_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            (time,ccn3,timemunit,dataunit,ccn3_longname)=read_E3SM(filename_input,'CCN3')
            (time,ccn5,timemunit,dataunit,ccn5_longname)=read_E3SM(filename_input,'CCN5')
            
            timem = np.hstack((timem,time))
            tmp_CCN3 = np.hstack((tmp_CCN3,ccn3))
            tmp_CCN5 = np.hstack((tmp_CCN5,ccn5))
        
        ccna_m.append(tmp_CCN3)
        ccnb_m.append(tmp_CCN5)
        
        # get supersaturation
        SS3 = ccn3_longname.split('=')[-1]
        SS5 = ccn5_longname.split('=')[-1]
        
    #%% make plot
        
    figname = figpath_sfc_timeseries+'timeseries_CCN_'+campaign+'_'+IOP+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax1.plot(t_ccna,ccna,color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax1.plot(timem, ccna_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax1.set_yscale('log')
    ax1.tick_params(color='k',labelsize=12)
    ylim1 = ax1.get_ylim()
    
    ax2.plot(t_ccnb,ccnb,color='k',linewidth=1,label='Obs')
    for mm in range(nmodels):
        ax2.plot(timem, ccnb_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax2.set_yscale('log')
    ax2.tick_params(color='k',labelsize=12)
    ylim2 = ax2.get_ylim()
    
    # ax1.set_yticks([10,100,1000,10000,100000])
    # ax2.set_yticks([10,100,1000,10000,100000])
    ax1.set_xlim(cday1,cday2)
    ax2.set_xlim(cday1,cday2)
    
    # set ylimit consistent in subplots
    ax1.set_ylim([ylim1[0], ylim2[1]])
    ax2.set_ylim([ylim1[0], ylim2[1]])
    
    # supersaturation
    fig.text(0.67,0.9,'SS_obs='+format(np.nanmean(SSa),'.2f')+'%, SS_model='+SS3)
    fig.text(0.67,0.4,'SS_obs='+format(np.nanmean(SSb),'.2f')+'%, SS_model='+SS5)
    
    ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    
    ax2.set_xlabel('Calendar Day',fontsize=14)
    ax1.set_title('CCN Number Concentration (cm$^{-3}$) '+campaign+' '+IOP,fontsize=15)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
