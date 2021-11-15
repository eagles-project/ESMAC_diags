"""
# plot surface diurnal cycle of aerosol composition
# compare models and surface measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import yyyymmdd2cday,cday2mmdd
from ..subroutines.read_ARMdata import read_acsm
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.specific_data_treatment import  avg_time_1d
from ..subroutines.quality_control import qc_remove_neg,qc_acsm_org_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    acsmpath = settings['acsmpath']
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
        if IOP=='IOP1':
            lst = glob.glob(acsmpath+'enaaosacsmC1.a1.201706*') + glob.glob(acsmpath+'enaaosacsmC1.a1.201707*')
        elif IOP=='IOP2':
            lst = glob.glob(acsmpath+'enaaosacsmC1.a1.201801*') + glob.glob(acsmpath+'enaaosacsmC1.a1.201802*')
        lst.sort()
    elif campaign=='HISCALE':  
        if IOP=='IOP1':
            lst = glob.glob(acsmpath+'sgpaosacsmC1.b1.201604*') + glob.glob(acsmpath+'sgpaosacsmC1.b1.201605*') + glob.glob(acsmpath+'sgpaosacsmC1.b1.201606*')
        elif IOP=='IOP2':
            lst = glob.glob(acsmpath+'sgpaosacsmC1.b1.201608*.cdf') + glob.glob(acsmpath+'sgpaosacsmC1.b1.201609*.cdf')
        lst.sort()
        
    t_obs=np.empty(0)
    so4_obs=np.empty(0)
    org_obs=np.empty(0)
    for filename in lst:
        (times_obs,so4sfc,timeunit,so4sfcunit)=read_acsm(filename,'sulfate')
        (times_obs,orgsfc,timeunit,orgsfcunit)=read_acsm(filename,'total_organics')
        timestr=timeunit.split(' ')
        date=timestr[2]
        cday=yyyymmdd2cday(date,'noleap')
        # average in time for quicker plot
        time2=np.arange(0,86400,3600)
        so42 = avg_time_1d(np.array(times_obs),np.array(so4sfc),time2)
        org2 = avg_time_1d(np.array(times_obs),np.array(orgsfc),time2)
        t_obs=np.hstack((t_obs, cday+time2/86400))
        so4_obs=np.hstack((so4_obs, so42))
        org_obs=np.hstack((org_obs, org2))
    
    so4_obs=qc_remove_neg(so4_obs)
    org_obs=qc_remove_neg(org_obs)
    org_obs=qc_acsm_org_max(org_obs)
        
        
    
    #%% read in models
    nmodels = len(Model_List)
    model_org = list()
    model_so4 = list()
    
    for mm in range(nmodels):
        so4varname=['so4_a1','so4_a2','so4_a3']
        orgvarname=['soa_a1','soa_a2','soa_a3','pom_a1','pom_a3','pom_a4',\
                    'mom_a1','mom_a2','mom_a3','mom_a4']
        if Model_List[mm]=='NucSoaCond':
            so4varname.append('so4_a5')
            orgvarname.append('soa_a5')
    
        timem2 = np.array([])
        tmp_so4 = np.empty(0)
        tmp_org = np.empty(0)
        ps = np.empty(0)
        ts = np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            filename_input = E3SM_sfc_path+'SFC_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            
            (timem,so4all,timeunitm,so4unit,so4name)=read_E3SM(filename_input,so4varname)
            (timem,orgall,timeunitm,orgunit,orgname)=read_E3SM(filename_input,orgvarname)
            (timem,[psm,tsm],timeunitm,varunit,varlongname)=read_E3SM(filename_input,['PS','T'])
            
            tmp_so4 = np.hstack((tmp_so4,sum(so4all)))
            tmp_org = np.hstack((tmp_org,sum(orgall)))
            ps = np.hstack((ps,psm))
            ts = np.hstack((ts,tsm))
            timem2 = np.hstack((timem2,timem))
        
        model_org.append(tmp_org)
        model_so4.append(tmp_so4)
    
    # change E3SM unit from kg/kg to ug/m3 
    rho = ps/287.06/ts
    
    for mm in range(nmodels):
        model_so4[mm]=model_so4[mm]*1e9*rho
        model_org[mm]=model_org[mm]*1e9*rho
        
    #%% calculate diurnal cycle
    days = np.arange(cday1, cday2+1)
    
    time_dc = np.arange(30,1440.,60)
    so4_o_dc = np.full((len(time_dc),len(days)),np.nan)
    org_o_dc = np.full((len(time_dc),len(days)),np.nan)
    n_valid = list()
    if len(so4_obs)>1: # not NaN
        for dd in range(len(days)):
            for tt in range(len(time_dc)):
                time_tmp = days[dd]+time_dc[tt]/1440.
                idx = np.abs(t_obs-time_tmp).argmin()
                if (t_obs[idx]-time_tmp)*1440 <= 30:    
                    so4_o_dc[tt,dd] = so4_obs[idx]
    if len(org_obs)>1:
        for dd in range(len(days)):
            for tt in range(len(time_dc)):
                time_tmp = days[dd]+time_dc[tt]/1440.
                idx = np.abs(t_obs-time_tmp).argmin()
                if (t_obs[idx]-time_tmp)*1440 <= 30:    
                    org_o_dc[tt,dd] = org_obs[idx]
    so4_o_dc = np.nanmean(so4_o_dc,1)
    org_o_dc = np.nanmean(org_o_dc,1)
    
    # for E3SM data
    so4_m_dc = []
    org_m_dc = []
    for mm in range(nmodels):
        tmp_so4 = np.full((24,len(days)),np.nan)
        tmp_org = np.full((24,len(days)),np.nan)
        for dd in range(len(days)):
            idx=np.logical_and(timem2>=days[dd], timem2<days[dd]+1)
            tmp_so4[:,dd] = model_so4[mm][idx]
            tmp_org[:,dd] = model_org[mm][idx]
        so4_m_dc.append(np.nanmean(tmp_so4,1))
        org_m_dc.append(np.nanmean(tmp_org,1))
        
    #%% make plot
        
    figname = figpath_sfc_timeseries+'diurnalcycle_AerosolComposition_'+campaign+'_'+IOP+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(6,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax1.plot(time_dc/60,so4_o_dc,color='k',linewidth=1,label='OBS (SO4)')
    for mm in range(nmodels):
        ax1.plot(time_dc/60, so4_m_dc[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax1.tick_params(color='k',labelsize=12)
    # ylim1 = ax1.get_ylim()
    
    ax2.plot(time_dc/60,org_o_dc,color='k',linewidth=1,label='OBS (ORG)')
    for mm in range(nmodels):
        ax2.plot(time_dc/60, org_m_dc[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax2.tick_params(color='k',labelsize=12)
    # ylim2 = ax2.get_ylim()
    
    ax1.set_xlim(0,24)
    ax2.set_xlim(0,24)
    ax1.set_xticks(np.arange(0,24,3))
    ax2.set_xticks(np.arange(0,24,3))
    
    # set ylimit consistent in subplots
    # ax1.set_yticks([10,100,1000,10000,100000])
    # ax2.set_yticks([10,100,1000,10000,100000])
    # ax1.set_yscale('log')
    # ax2.set_yscale('log')
    # ax1.set_ylim([ylim1[0], ylim2[1]])
    # ax2.set_ylim([ylim1[0], ylim2[1]])
    
    
    ax1.legend(loc='center right', shadow=False, fontsize='medium',bbox_to_anchor=(1.3, .5))
    ax2.legend(loc='center right', shadow=False, fontsize='medium',bbox_to_anchor=(1.3, .5))
    
    ax2.set_xlabel('Hour (UTC)',fontsize=12)
    ax1.set_title('Aerosol Sulfate and Organic Concentration ($\mu$g/m$^3$) '+campaign+' '+IOP,fontsize=14)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    # plt.close()