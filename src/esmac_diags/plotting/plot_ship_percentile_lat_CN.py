"""
# plot ship-track aerosol number concentration binned by different latitudes
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_cpc, read_uhsas
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import  cday2mmdd
from ..subroutines.quality_control import qc_mask_qcflag,qc_remove_neg,qc_cn_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    shipcpcpath = settings['shipcpcpath']
    shipuhsaspath = settings['shipuhsaspath']
    E3SM_ship_path = settings['E3SM_ship_path']
    figpath_ship_statistics = settings['figpath_ship_statistics']

    #%% other settings
    
    dlat = latbin[1]-latbin[0]
    latmin = latbin-dlat/2
    latmax = latbin+dlat/2
    latlen = len(latbin)
    
    if not os.path.exists(figpath_ship_statistics):
        os.makedirs(figpath_ship_statistics)
    
    
    #%% read in model
    nmodels=len(Model_List)
    cpc_m = list()
    uhsas_m = list()
    for mm in range(nmodels):
        # initialize variables by latitude bins
        cpc_tmp = list()
        uh_tmp = list()
        for bb in range(latlen):
            cpc_tmp.append(np.empty(0))
            uh_tmp.append(np.empty(0))
        
        lst = glob.glob(E3SM_ship_path+'Ship_CNsize_'+campaign+'_'+Model_List[mm]+'_shipleg*.nc')
            
        for ll in range(len(lst)):        
            filenamem = lst[ll]
            (timem,varm,timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['NCN','NCNall','lat','lon'])
            for ii in range(len(varm)):
                varm[ii][varm[ii]<-9000] = np.nan
            
            lat0=varm[2]
            NCN=varm[0]        
            NCN_uh = np.nansum(varm[1][54:1000,:],0)   # UHSAS CN size range: 55-1000nm
            
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat0>=latmin[bb], lat0<latmax[bb])
                cpc_tmp[bb]=np.hstack((cpc_tmp[bb],NCN[idx]*1e-6))   # change unit from 1/m3 to 1/cm3
                uh_tmp[bb]=np.hstack((uh_tmp[bb],NCN_uh[idx]*1e-6))
            
        cpc_m.append(cpc_tmp)
        uhsas_m.append(uh_tmp)
            
    #%% read in observation
    cpc_o = list()
    uhsas_o = list()
    for bb in range(latlen):
        cpc_o.append(np.empty(0))
        uhsas_o.append(np.empty(0))
        
    for ll in range(len(lst)):       
        # use lat/lon from extracted model data
        filenamem = lst[ll]
        if campaign=='MAGIC':
            legnum=lst[ll][-5:-3]
        elif campaign=='MARCUS':
            legnum=lst[ll][-4]
        (timem,[lat1,lon1],timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['lat','lon'])
        lat1[lat1<-9000]=np.nan
        lon1[lon1<-9000]=np.nan
            
        # find the days related to the ship leg
        day = [int(a) for a in timem]
        day = list(set(day))
        day.sort()
        
        # read in CPC    
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
        lat2=np.interp(t_cpc,timem,lat1)
        # separate into latitude bins
        for bb in range(latlen):
            idx = np.logical_and(lat2>=latmin[bb], lat2<latmax[bb])
            cpc_o[bb]=np.hstack((cpc_o[bb],cpc[idx]))  
        
        # read in UHSAS
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
        lat3=np.interp(t_uh,timem,lat1)
        # separate into latitude bins
        for bb in range(latlen):
            idx = np.logical_and(lat3>=latmin[bb], lat3<latmax[bb])
            uhsas_o[bb]=np.hstack((uhsas_o[bb],uhsas[idx]))  
            
    #%% 
    for bb in range(latlen):
        cpc_o[bb] = cpc_o[bb][~np.isnan(cpc_o[bb])]
        uhsas_o[bb] = uhsas_o[bb][~np.isnan(uhsas_o[bb])]
    
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
        
    figname = figpath_ship_statistics+'percentile_lat_CN_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(cpc_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(cpc_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax1.tick_params(color='k',labelsize=15)
    ax1.set_yscale('log')
    ax1.set_xlim(-1,latlen)
    ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    ax1.set_xticklabels([])
    # plot temporal lines for label
    ax1.plot([],c='k',label='CPC')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
    ax1.legend(loc='upper left', fontsize='x-large')
        
    ax2.boxplot(uhsas_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax2.boxplot(uhsas_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax2.tick_params(color='k',labelsize=15)
    ax2.set_yscale('log')
    ax2.set_ylim(10,3000)
    ax2.set_yticks([10,100,1000])
    ax2.set_xlim(-1,latlen)
    ax2.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    ax2.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
    # plot temporal lines for label
    ax2.plot([],c='k',label='UHSAS')
    for mm in range(nmodels):
        ax2.plot([],c=color_model[mm],label=Model_List[mm])
    ax2.legend(loc='upper left', fontsize='x-large')
        
    # # set xlimit consistent in subplots
    # ylim1 = ax1.get_ylim()
    # ylim2 = ax2.get_ylim()
    # ax1.set_ylim([min(ylim1[0],ylim2[0]), max(ylim1[1],ylim2[1])])
    # ax2.set_ylim([min(ylim1[0],ylim2[0]), max(ylim1[1],ylim2[1])])
    
    ax2.set_xlabel('Latitude',fontsize=16)
    ax1.set_title('Aerosol Number Concentration (cm$^{-3}$)',fontsize=17)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
