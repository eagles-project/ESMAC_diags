"""
# plot timeseries of surface CCN number concentration along each ship leg
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_ccn_magic, read_ccn
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import cday2mmdd
from ..subroutines.specific_data_treatment import mask_model_ps
from ..subroutines.quality_control import qc_mask_qcflag,qc_ccn_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    shipccnpath = settings['shipccnpath']
    shipmetpath = settings['shipmetpath']
    E3SM_ship_path = settings['E3SM_ship_path']
    figpath_ship_timeseries = settings['figpath_ship_timeseries']

    #%% other settings
    
    if not os.path.exists(figpath_ship_timeseries):
        os.makedirs(figpath_ship_timeseries)
    
    lst = glob.glob(E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[0]+'_shipleg*.nc')
    lst.sort()
    
    for ll in range(len(lst)):
        
        if campaign=='MAGIC':
            legnum=lst[ll][-5:-3]
        elif campaign=='MARCUS':
            legnum=lst[ll][-4]
        
         
        #%% read in model
        nmodels=len(Model_List)
        ccn1m = list()
        ccn5m = list()
        for mm in range(nmodels):
            filenamem = E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
        
            (timem,ccn1,timeunitm,ccn1unit,ccn1longname)=read_E3SM(filenamem,'CCN3')
            (timem,ccn5,timeunitm,ccn5unit,ccn5longname)=read_E3SM(filenamem,'CCN5')
        
            ccn1m.append(ccn1)
            ccn5m.append(ccn5)
            
        # mask data where model grid is not at ocean surface (Ps is too different than obs)
        (timem,psm,timeunitm,psmunit,psmlongname)=read_E3SM(filenamem,'PS')
        datamask = mask_model_ps(timem,0.01*psm,legnum,campaign,shipmetpath)
        # for mm in range(nmodels):
        #     ccn1m[mm][datamask]=np.nan
        #     ccn5m[mm][datamask]=np.nan
        
        year0 = str(int(timeunitm.split()[2][0:4])+1)
        
        #%% read in observations
        # find the days related to the ship leg
        day = [int(a) for a in timem]
        day = list(set(day))
        day.sort()
        
        t_ccn=np.empty(0)
        ccn=np.empty(0)
        SS=np.empty(0)
        for dd in day:
            
            if campaign=='MAGIC':
                if int(legnum)<=9:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipccnpath+'magaosccn100M1.a1.2012'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
                    else:
                        filenameo = glob.glob(shipccnpath+'magaosccn100M1.a1.2013'+cday2mmdd(dd-365,calendar='noleap')+'.*.cdf')
                else:
                    filenameo = glob.glob(shipccnpath+'magaosccn100M1.a1.2013'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
                if len(filenameo)==0:
                    continue  # some days may be missing
                (time,timeunit,obs,dataunit,SS0)=read_ccn_magic(filenameo[0])
            elif campaign=='MARCUS':
                if int(legnum)<=2:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipccnpath+'maraosccn1colavgM1.b1.2017'+cday2mmdd(dd,calendar='noleap')+'.*')
                    else:
                        filenameo = glob.glob(shipccnpath+'maraosccn1colavgM1.b1.2018'+cday2mmdd(dd-365,calendar='noleap')+'.*')
                else:
                    filenameo = glob.glob(shipccnpath+'maraosccn1colavgM1.b1.2018'+cday2mmdd(dd,calendar='noleap')+'.*')
                if len(filenameo)==0:
                    continue  # some days may be missing
                (time,timeunit,obs,qc,dataunit,SS0)=read_ccn(filenameo[0])            
                obs=qc_mask_qcflag(obs,qc)
                
            t_ccn=np.hstack((t_ccn, dd+time/86400))
            ccn=np.hstack((ccn, obs))
            SS=np.hstack((SS, SS0))
            ccn=qc_ccn_max(ccn,SS)
            
        # if time expands two years, add 365 days to the second year
        if t_ccn[0]>t_ccn[-1]:
            t_ccn[t_ccn<=t_ccn[-1]]=t_ccn[t_ccn<=t_ccn[-1]]+365
            
        # ccn[np.logical_or(ccn<0,ccn>1500)]=np.nan
        # SS=0.1%
        idx = np.logical_and(SS>0.05, SS<0.15)
        t_ccn1 = t_ccn[idx]
        ccn1o = ccn[idx]
        SS1 = 0.1
        # SS=0.5%
        idx = np.logical_and(SS>0.4, SS<0.6)
        t_ccn5 = t_ccn[idx]
        ccn5o = ccn[idx]
        SS5 = 0.5
        
        #%% make plot
            
        figname = figpath_ship_timeseries+'timeseries_CCN_'+campaign+'_ship'+legnum+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
        ax1.plot(t_ccn1,ccn1o,color='k',linewidth=1,label='OBS')
        for mm in range(nmodels):
            ax1.plot(timem, ccn1m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        # ax1.set_yscale('log')
        ax1.tick_params(color='k',labelsize=12)
        ylim1 = ax1.get_ylim()
        
        ax2.plot(t_ccn5,ccn5o,color='k',linewidth=1,label='OBS')
        for mm in range(nmodels):
            ax2.plot(timem, ccn5m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        # ax2.set_yscale('log')
        ax2.tick_params(color='k',labelsize=12)
        ylim2 = ax2.get_ylim()
        
        # set ylimit consistent in subplots
        # ax1.set_ylim([ylim1[0], ylim2[1]])
        # ax2.set_ylim([ylim1[0], ylim2[1]])
        
        ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.2, .5))
        ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.2, .5))
        
        # supersaturation
        fig.text(0.08,0.9,'SS='+str(SS1)+'%')
        fig.text(0.08,0.4,'SS='+str(SS5)+'%')
        
        ax2.set_xlabel('Calendar Day in '+year0,fontsize=14)
            
        ax1.set_title('CCN Number Concentration (cm$^{-3}$)',fontsize=15)
        
        fig.text(.08, .999,'trip # '+legnum, fontsize=12)
        
        # mask non-ocean model grid (ps is inconsistent with obs)
        ax1.vlines(timem[datamask],ylim1[0],ylim1[1],color='lightgray')
        ax2.vlines(timem[datamask],ylim2[0],ylim2[1],color='lightgray')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()