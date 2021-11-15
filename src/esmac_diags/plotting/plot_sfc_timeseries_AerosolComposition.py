"""
# plot surface timeseries of aerosol composition
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
        time2=np.arange(1800,86400,3600)
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
    
    #%% make plot
        
    figname = figpath_sfc_timeseries+'timeseries_AerosolComposition_'+campaign+'_'+IOP+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax1.plot(t_obs,so4_obs,color='k',linewidth=1,label='OBS (SO4)')
    for mm in range(nmodels):
        ax1.plot(timem2, model_so4[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax1.set_yscale('log')
    ax1.tick_params(color='k',labelsize=12)
    
    ax2.plot(t_obs,org_obs,color='k',linewidth=1,label='OBS (ORG)')
    for mm in range(nmodels):
        ax2.plot(timem2, model_org[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax2.set_yscale('log')
    ax2.tick_params(color='k',labelsize=12)
    
    ax1.set_xlim(cday1,cday2)
    ax2.set_xlim(cday1,cday2)
    
    ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    
    ax2.set_xlabel('Calendar Day',fontsize=14)
    
    ax1.set_title('Aerosol Sulfate and Organic Concentration ($\mu$g/m$^3$) '+campaign+' '+IOP,fontsize=14)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
