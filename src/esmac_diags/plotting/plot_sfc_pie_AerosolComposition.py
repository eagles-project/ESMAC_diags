"""
# plot surface aerosol composition in a pie plot
# plot models and surface measurements separately
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import yyyymmdd2cday,cday2mmdd
from ..subroutines.read_ARMdata import read_acsm
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.quality_control import qc_remove_neg,qc_acsm_org_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    acsmpath = settings['acsmpath']
    start_date = settings['start_date']
    end_date = settings['end_date']
    E3SM_sfc_path = settings['E3SM_sfc_path']
    figpath_sfc_statistics = settings['figpath_sfc_statistics']

    IOP = settings.get('IOP', None)

    #%% other settings
    
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date,'noleap')
    cday2 = yyyymmdd2cday(end_date,'noleap')
    if start_date[0:4]!=end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
        
    import os
    if not os.path.exists(figpath_sfc_statistics):
        os.makedirs(figpath_sfc_statistics)
        
                
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
    else:
        raise ValueError('surface aerosol composition is only available in HISCALE or ACEENA. check: '+campaign)
        
    t_obs=np.empty(0)
    so4_obs=np.empty(0)
    org_obs=np.empty(0)
    nh4_obs=np.empty(0)
    no3_obs=np.empty(0)
    chl_obs=np.empty(0)
    for filename in lst:
        (times_obs,so4sfc,timeunit,so4sfcunit)=read_acsm(filename,'sulfate')
        (times_obs,orgsfc,timeunit,orgsfcunit)=read_acsm(filename,'total_organics')
        (times_obs,nh4sfc,timeunit,nh4sfcunit)=read_acsm(filename,'ammonium')
        (times_obs,no3sfc,timeunit,no3sfcunit)=read_acsm(filename,'nitrate')
        (times_obs,chlsfc,timeunit,chlsfcunit)=read_acsm(filename,'chloride')
        timestr=timeunit.split(' ')
        date=timestr[2]
        cday=yyyymmdd2cday(date,'noleap')
        so4_obs=np.hstack((so4_obs, so4sfc))
        org_obs=np.hstack((org_obs, orgsfc))
        nh4_obs=np.hstack((nh4_obs, nh4sfc))
        no3_obs=np.hstack((no3_obs, no3sfc))
        chl_obs=np.hstack((chl_obs, chlsfc))
    so4_obs=qc_remove_neg(so4_obs)
    nh4_obs=qc_remove_neg(nh4_obs)
    no3_obs=qc_remove_neg(no3_obs)
    chl_obs=qc_remove_neg(chl_obs)
    org_obs=qc_remove_neg(org_obs)
    org_obs=qc_acsm_org_max(org_obs)
    
    #%% read in models
    nmodels = len(Model_List)
    model_org = list()
    model_so4 = list()
    model_bc = list()
    model_dst = list()
    model_ncl = list()
    
    for mm in range(nmodels):
        bcvarname=['bc_a1','bc_a3','bc_a4']
        dstvarname=['dst_a1','dst_a3']
        nclvarname=['ncl_a1','ncl_a2','ncl_a3']
        so4varname=['so4_a1','so4_a2','so4_a3']
        orgvarname=['soa_a1','soa_a2','soa_a3','pom_a1','pom_a3','pom_a4',\
                    'mom_a1','mom_a2','mom_a3','mom_a4']
        if Model_List[mm]=='NucSoaCond':
            so4varname.append('so4_a5')
            orgvarname.append('soa_a5')
    
        timem2 = np.array([])
        tmp_so4 = np.empty(0)
        tmp_org = np.empty(0)
        tmp_bc = np.empty(0)
        tmp_dst = np.empty(0)
        tmp_ncl = np.empty(0)
        ps = np.empty(0)
        ts = np.empty(0)
        for cday in range(cday1,cday2+1):
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
            filename_input = E3SM_sfc_path+'SFC_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
            
            (timem,so4all,timeunitm,so4unit,so4name)=read_E3SM(filename_input,so4varname)
            (timem,orgall,timeunitm,orgunit,orgname)=read_E3SM(filename_input,orgvarname)
            (timem,bcall,timeunitm,bcunit,bcname)=read_E3SM(filename_input,bcvarname)
            (timem,dstall,timeunitm,dstunit,dstname)=read_E3SM(filename_input,dstvarname)
            (timem,nclall,timeunitm,nclunit,nclname)=read_E3SM(filename_input,nclvarname)
            (timem,[psm,tsm],timeunitm,varunit,varlongname)=read_E3SM(filename_input,['PS','T'])
            
            tmp_so4 = np.hstack((tmp_so4,sum(so4all)))
            tmp_org = np.hstack((tmp_org,sum(orgall)))
            tmp_bc = np.hstack((tmp_bc,sum(bcall)))
            tmp_dst = np.hstack((tmp_dst,sum(dstall)))
            tmp_ncl = np.hstack((tmp_ncl,sum(nclall)))
            ps = np.hstack((ps,psm))
            ts = np.hstack((ts,tsm))
            timem2 = np.hstack((timem2,timem))
        
        model_so4.append(tmp_so4)
        model_org.append(tmp_org)
        model_bc.append(tmp_bc)
        model_dst.append(tmp_dst)
        model_ncl.append(tmp_ncl)
    
    # change E3SM unit from kg/kg to ug/m3 
    rho = ps/287.06/ts
    
    for mm in range(nmodels):
        model_so4[mm]=model_so4[mm]*1e9*rho
        model_org[mm]=model_org[mm]*1e9*rho
        model_bc[mm]=model_bc[mm]*1e9*rho
        model_dst[mm]=model_dst[mm]*1e9*rho
        model_ncl[mm]=model_ncl[mm]*1e9*rho
        
    #%% Pie plot
    
    figname = figpath_sfc_statistics+'Pieplot_AerosolComposition_'+campaign+'_'+IOP+'.png'
    print('plotting figures to '+figname)
    
    fig,ax = plt.subplots(1,nmodels+1,figsize=((nmodels+1)*3.5,3.5))   # figsize in inches
    # colors = ['limegreen', 'red', 'b', 'y', 'orange' ]
    
    colors_o = ['limegreen', 'red', 'orange', 'lightblue', 'yellow']
    labels_o = ['ORG', 'SO4', 'NO3', 'NH4', 'CHL']
    sizeo = [np.nanmean(org_obs),np.nanmean(so4_obs),np.nanmean(no3_obs),np.nanmean(nh4_obs),np.nanmean(chl_obs)]
    
    colors_m = ['limegreen', 'red', 'k', 'silver','gray']
    labels_m = ['ORG', 'SO4', 'BC', 'DST', 'NCL']
    sizem = []
    for mm in range(nmodels):
        sizem.append([np.mean(model_org[mm]),np.mean(model_so4[mm]),np.mean(model_bc[mm]),np.mean(model_dst[mm]),np.mean(model_ncl[mm])])
    
    def absolute_value(val):
        a=np.round(val*sum(sizeo))/100
        return a
    ax[0].pie(sizeo/sum(sizeo),labels=labels_o,colors=colors_o, autopct=absolute_value)  # autopct='%1.1f%%'
    for mm in range(nmodels):
        def absolute_valuemm(val):
            a=np.round(val*sum(sizem[mm]))/100
            return a
        ax[mm+1].pie(sizem[mm]/sum(sizem[mm]),labels=labels_m, colors=colors_m, autopct=absolute_valuemm)
        
    ax[0].set_title('Obs',fontsize=14)
    for mm in range(nmodels):
        ax[mm+1].set_title(Model_List[mm],fontsize=14)
    fig.text(.5,.15,'unit: $\mu$g/m$^3$')
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)