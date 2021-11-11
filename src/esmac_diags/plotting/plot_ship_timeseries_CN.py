"""
# plot timeseries of surface aerosol number concentration along each ship leg
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_cpc, read_uhsas
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import  cday2mmdd
from ..subroutines.specific_data_treatment import mask_model_ps
from ..subroutines.quality_control import qc_mask_qcflag,qc_remove_neg,qc_cn_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    shipcpcpath = settings['shipcpcpath']
    shipuhsaspath = settings['shipuhsaspath']
    shipmetpath = settings['shipmetpath']
    E3SM_ship_path = settings['E3SM_ship_path']
    figpath_ship_timeseries = settings['figpath_ship_timeseries']

    #%% other settings
    
    if not os.path.exists(figpath_ship_timeseries):
        os.makedirs(figpath_ship_timeseries)
    
    lst = glob.glob(E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[0]+'_shipleg*.nc')
    lst.sort()
    
    for ll in range(len(lst)):
        
        if campaign=='MAGIC':
            legnum=lst[ll][-5:-3]
        elif campaign=='MARCUS':
            legnum=lst[ll][-4]
        
        #%% read in model
        nmodels=len(Model_List)
        datam = list()
        databins = list()
        for mm in range(nmodels):
            filenamem = E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
        
            (timem,NCNall,timeunitm,datamunit,datamlongname)=read_E3SM(filenamem,'NCNall')
            (timem,data,timeunitm,datamunit,datamlongname)=read_E3SM(filenamem,'NCN')
        
            datam.append(data*1e-6)    # change unit from 1/m3 to 1/cm3
            databins.append(NCNall*1e-6)    # change unit from 1/m3 to 1/cm3
            
            # mask data where model grid is not at ocean surface (Ps is too different than obs)
            filenamem = E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
            (timem,psm,timeunitx,psmunit,psmlongname)=read_E3SM(filenamem,'PS')
            datamask = mask_model_ps(timem,0.01*psm,legnum,campaign,shipmetpath)
            # for mm in range(nmodels):
            #     datam[mm][datamask]=np.nan
            
        year0 = str(int(timeunitm.split()[2][0:4])+1)
        
        #%% read in observations
        # find the days related to the ship leg
        day = [int(a) for a in timem]
        day = list(set(day))
        day.sort()
    
        # CPC    
        t_cpc=np.empty(0)
        cpc=np.empty(0)
        for dd in day:
            
            if campaign=='MAGIC':
                if int(legnum)<=9:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipcpcpath+'magaoscpcfM1.a1.2012'+cday2mmdd(dd,calendar='noleap')+'.*')
                    else:
                        filenameo = glob.glob(shipcpcpath+'magaoscpcfM1.a1.2013'+cday2mmdd(dd-365,calendar='noleap')+'.*')
                else:
                    filenameo = glob.glob(shipcpcpath+'magaoscpcfM1.a1.2013'+cday2mmdd(dd,calendar='noleap')+'.*')
                if len(filenameo)==0:
                    continue  # some days may be missing
            elif campaign=='MARCUS':
                if int(legnum)<=2:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipcpcpath+'maraoscpcf1mM1.b1.2017'+cday2mmdd(dd,calendar='noleap')+'.*')
                    else:
                        filenameo = glob.glob(shipcpcpath+'maraoscpcf1mM1.b1.2018'+cday2mmdd(dd-365,calendar='noleap')+'.*')
                else:
                    filenameo = glob.glob(shipcpcpath+'maraoscpcf1mM1.b1.2018'+cday2mmdd(dd,calendar='noleap')+'.*')
                if len(filenameo)==0:
                    continue  # some days may be missing
                    
                    
            (time,obs,qc,timeunit,dataunit)=read_cpc(filenameo[0])
            obs=qc_mask_qcflag(obs,qc)
            t_cpc=np.hstack((t_cpc, dd+time/86400))
            cpc=np.hstack((cpc, obs))
            
        cpc=qc_remove_neg(cpc)
        cpc=qc_cn_max(cpc,10)
        # if time expands two years, add 365 days to the second year
        if t_cpc[0]>t_cpc[-1]:
            t_cpc[t_cpc<=t_cpc[-1]]=t_cpc[t_cpc<=t_cpc[-1]]+365
    
        # UHSAS
        t_uh=np.empty(0)
        uhsas=np.empty(0)
        for dd in day:
            
            if campaign=='MAGIC':
                if int(legnum)<=9:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2012'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
                    else:
                        filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2013'+cday2mmdd(dd-365,calendar='noleap')+'.*.cdf')
                else:
                    filenameo = glob.glob(shipuhsaspath+'magaosuhsasM1.a1.2013'+cday2mmdd(dd,calendar='noleap')+'.*.cdf')
            elif campaign=='MARCUS':
                if int(legnum)<=2:
                    if dd<=365:  # year 2012
                        filenameo = glob.glob(shipuhsaspath+'maraosuhsasM1.a1.2017'+cday2mmdd(dd,calendar='noleap')+'.*')
                    else:
                        filenameo = glob.glob(shipuhsaspath+'maraosuhsasM1.a1.2018'+cday2mmdd(dd-365,calendar='noleap')+'.*')
                else:
                    filenameo = glob.glob(shipuhsaspath+'maraosuhsasM1.a1.2018'+cday2mmdd(dd,calendar='noleap')+'.*')
            
            if len(filenameo)==0:
                continue  # some days may be missing
            if len(filenameo)>1:
                raise ValueError('find too many files')
                
            (time,dmin,dmax,obs,timeunit,uhunit,uhlongname)=read_uhsas(filenameo[0])
            obs=np.ma.filled(obs)
            obs=qc_remove_neg(obs)
            uhsas=np.hstack((uhsas, np.nansum(obs,1)))
            t_uh = np.hstack((t_uh,time/86400+dd))
            
        uhsas=qc_cn_max(uhsas,100)
        # if no obs available, fill one data with NaN
        if len(t_uh)==0:
            t_uh=[timem[0],timem[1]]
            uhsas=np.full((2),np.nan)
            
        # if time expands two years, add 365 days to the second year
        if t_uh[0]>t_uh[-1]:
            t_uh[t_uh<=t_uh[-1]]=t_uh[t_uh<=t_uh[-1]]+365
                
        #%% Calculate model aerosol number concentration for UHSAS size range
        b1 = int(dmin[0])
        b2 = int(dmax[-1])
        datam2=list()
        for mm in range(nmodels):
            datam2.append(np.nansum(databins[mm][b1-1:b2,:],0))
            # datam2[mm][datamask]=np.nan
        
        #%% make plot
            
        figname = figpath_ship_timeseries+'timeseries_CN_'+campaign+'_ship'+legnum+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
        ax1.plot(t_cpc,cpc,color='k',linewidth=1,label='CPC')
        for mm in range(nmodels):
            ax1.plot(timem, datam[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        # ax1.set_yscale('log')
        ax1.tick_params(color='k',labelsize=12)
        ylim1 = ax1.get_ylim()
        
        ax2.plot(t_uh,uhsas,color='k',linewidth=1,label='UHSAS')
        # ax2.plot(t_uh,uhsas,color='k',linewidth=1,label='UHSAS ('+str(b1)+'-'+str(b2)+'nm)')
        for mm in range(nmodels):
            ax2.plot(timem, datam2[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
        # ax1.set_yscale('log')
        ax2.tick_params(color='k',labelsize=12)
        ylim2 = ax2.get_ylim()
        
        ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.2, .5))
        ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.2, .5))
        
        ax2.set_xlabel('Calendar Day in '+year0,fontsize=14)
            
        ax1.set_title('Aerosol Number Concentration (cm$^{-3}$)',fontsize=15)
        
        fig.text(.08, .999,'trip # '+legnum, fontsize=12)
        
        # mask non-ocean model grid (ps is inconsistent with obs)
        ax1.vlines(timem[datamask],ylim1[0],ylim1[1],color='lightgray')
        ax2.vlines(timem[datamask],ylim2[0],ylim2[1],color='lightgray')
        
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()