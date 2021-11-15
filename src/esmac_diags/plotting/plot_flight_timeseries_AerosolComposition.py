"""
# plot aircraft track data
# timeseries of aerosol composition (SO4 and total organic) concentration 
# compare models and aircraft measurements
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.time_format_change import hhmmss2sec
from ..subroutines.read_aircraft import read_ams,read_iwg1
from ..subroutines.read_netcdf import read_merged_size,read_extractflight
from ..subroutines.quality_control import qc_mask_qcflag,qc_remove_neg

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_timeseries = settings['figpath_aircraft_timeseries']

    IOP = settings.get('IOP', None)
    merged_size_path = settings.get('merged_size_path', None)
    amspath = settings.get('amspath', None)
    iwgpath = settings.get('iwgpath', None)
    
    if campaign in ['CSET', 'SOCRATES']:
        raise ValueError('CSET and SOCRATES do not have composition data')
    elif campaign not in ['HISCALE', 'ACEENA']:
        raise ValueError('campaign name is not recognized: '+campaign)

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
        
    # for each flight
    for date in alldates:
        
        if date[-1]=='a':
            flightidx=1
        else:
            flightidx=2
    
        #% read in flight information
        if campaign=='HISCALE':
            filename = merged_size_path+'merged_bin_fims_pcasp_'+campaign+'_'+date+'.nc'
        elif campaign=='ACEENA':
            filename = merged_size_path+'merged_bin_fims_pcasp_opc_'+campaign+'_'+date+'.nc'
        (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
        (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
        (time,size,height,timeunit,zunit,long_name)=read_merged_size(filename,'height')
        time=np.ma.compressed(time)
        
        #%% read T and P from iwg
        filename_i=glob.glob(iwgpath+'aaf.iwg*.'+date+'*txt')
        filename_i.sort()
        # read in data
        if len(filename_i)==1: 
            (iwg,iwgvars)=read_iwg1(filename_i[0])
            timelen = len(iwg)
            if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                iwg.insert(1403,list(iwg[1403]))
                tstr=iwg[1403][1]
                tstr=tstr[0:-1]+str(int(tstr[-1])-1)
                iwg[1403][1]=tstr
                del iwg[-1]
            # get variables
            time_iwg=np.empty(timelen)
            T_iwg=np.empty(timelen)
            P_iwg=np.empty(timelen)
            for t in range(timelen):
                T_iwg[t]=float(iwg[t][20])+273.15
                P_iwg[t]=float(iwg[t][23])*100
                timestr=iwg[t][1].split(' ')
                time_iwg[t]=hhmmss2sec(timestr[1])
        else:
            raise ValueError('find no file or multiple files: ' + filename_i)
        
        #%% read aerosol composition in AMS
        
        filename_ams=glob.glob(amspath+'*'+date[0:8]+'*')
        filename_ams.sort()
        
        if len(filename_ams)==1 or len(filename_ams)==2:
            (ams,amslist)=read_ams(filename_ams[flightidx-1])
            time_ams=ams[0,:]
            flag=ams[-1,:]
            orgaaf=ams[1,:]
            so4aaf=ams[5,:]
            # flag=1 is also good data but behind CVI inlet. currently only use good data behind isokinetic inlet (flag=0)
            orgaaf=qc_mask_qcflag(orgaaf,flag)
            so4aaf=qc_mask_qcflag(so4aaf,flag)
        elif len(filename_ams)==0:
            time_ams = time_iwg
            orgaaf = np.full(len(time_ams),np.nan)
            so4aaf = np.full(len(time_ams),np.nan)
        else:
            raise ValueError('find too many files')
        
        # change values from standardize condition to ambient condition
        T_ams = np.interp(time_ams,time,T_iwg)
        P_ams = np.interp(time_ams,time,P_iwg)
        so4aaf = so4aaf * (296.15/T_ams) * (P_ams/101325.)
        orgaaf = orgaaf * (296.15/T_ams) * (P_ams/101325.)
        
        # some quality check:
        orgaaf=qc_remove_neg(orgaaf)
        so4aaf=qc_remove_neg(so4aaf)
        
        
        #%% read in Models
        nmodels=len(Model_List)
        so4_m = []
        org_m = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,soa_a1,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a1')
            (timem,heightm,soa_a2,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a2')
            (timem,heightm,soa_a3,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a3')
            (timem,heightm,so4_a1,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a1')
            (timem,heightm,so4_a2,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a2')
            (timem,heightm,so4_a3,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a3')
            (timem,heightm,pom_a1,timeunitm,pomunit,pomname)=read_extractflight(filename_m,'pom_a1')
            (timem,heightm,pom_a3,timeunitm,pomunit,pomname)=read_extractflight(filename_m,'pom_a3')
            (timem,heightm,pom_a4,timeunitm,pomunit,pomname)=read_extractflight(filename_m,'pom_a4')
            (timem,heightm,mom_a1,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a1')
            (timem,heightm,mom_a2,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a2')
            (timem,heightm,mom_a3,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a3')
            (timem,heightm,mom_a4,timeunitm,momunit,momname)=read_extractflight(filename_m,'mom_a4')
            
            # add nucleation mode if available
            try:
                (timem,heightm,soa_a5,timeunitm,soaunit,soaname)=read_extractflight(filename_m,'soa_a5')
                model_org = soa_a1+soa_a2+soa_a3+soa_a5 + pom_a1+pom_a3+pom_a4 + mom_a1+mom_a2+mom_a3+mom_a4
            except:
                model_org = soa_a1+soa_a2+soa_a3 + pom_a1+pom_a3+pom_a4 + mom_a1+mom_a2+mom_a3+mom_a4
            try:
                (timem,heightm,so4_a5,timeunitm,so4unit,so4name)=read_extractflight(filename_m,'so4_a5')
                model_so4 = so4_a1+so4_a2+so4_a3+so4_a5
            except:
                model_so4 = so4_a1+so4_a2+so4_a3
            
            # change E3SM unit from kg/kg to ug/m3 
            rho = P_iwg/T_iwg/287.06
            model_so4=model_so4*1e9*rho
            model_org=model_org*1e9*rho
            
            so4_m.append(model_so4) 
            org_m.append(model_org) 
        
        timem2 = timem/3600
           
        #%% make plot
            
        figname = figpath_aircraft_timeseries+'AerosolComposition_'+campaign+'_'+date+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=2.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
        ax1.plot(time_ams/3600,so4aaf,color='k',linewidth=1,label='OBS')
        for mm in range(nmodels):
            ax1.plot(timem2, so4_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        # ax1.set_yscale('log')
        ax1.tick_params(color='k',labelsize=12)
        # ylim1 = ax1.get_ylim()
        
        ax2.plot(time_ams/3600,orgaaf,color='k',linewidth=1,label='OBS')
        for mm in range(nmodels):
            ax2.plot(timem2, org_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        # ax2.set_yscale('log')
        ax2.tick_params(color='k',labelsize=12)
        # ylim2 = ax2.get_ylim()
        
        ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
        ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
        
        ax2.set_xlabel('time (hour UTC) '+date,fontsize=14)
    
        ax1.set_title('Aerosol Sulfate Concentration ($\mu$g/m$^3$)',fontsize=13)
        ax2.set_title('Aerosol Organic Concentration ($\mu$g/m$^3$)',fontsize=13)
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()   