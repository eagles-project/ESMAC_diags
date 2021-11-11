"""
# plot ship-track CCN number concentration binned by different latitudes
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_ccn_magic, read_ccn
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import  cday2mmdd
from ..subroutines.quality_control import qc_mask_qcflag,qc_ccn_max

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    shipccnpath = settings['shipccnpath']
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
    ccn3_m = list()
    ccn5_m = list()
    for mm in range(nmodels):
        # initialize variables by latitude bins
        ccn3_tmp = list()
        ccn5_tmp = list()
        for bb in range(latlen):
            ccn3_tmp.append(np.empty(0))
            ccn5_tmp.append(np.empty(0))
        
        lst = glob.glob(E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg*.nc')
            
        for ll in range(len(lst)):        
            filenamem = lst[ll]
            (timem,varm,timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['CCN3','CCN5','lat','lon'])
            for ii in range(len(varm)):
                varm[ii][varm[ii]<-9000] = np.nan
            
            lat0=varm[2]
            lon0=varm[3]
            ccn3=varm[0]        
            ccn5=varm[1]
            
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat0>=latmin[bb], lat0<latmax[bb])
                ccn3_tmp[bb]=np.hstack((ccn3_tmp[bb],ccn3[idx]))  
                ccn5_tmp[bb]=np.hstack((ccn5_tmp[bb],ccn5[idx]))
            
        ccn3_m.append(ccn3_tmp)
        ccn5_m.append(ccn5_tmp)
            
    #%% read in observation
    ccn3_o = list()
    ccn5_o = list()
    for bb in range(latlen):
        ccn3_o.append(np.empty(0))
        ccn5_o.append(np.empty(0))
        
    for ll in range(len(lst)):       
        # use lat/lon from extracted model data
        filenamem = lst[ll]
        if campaign=='MAGIC':
            legnum=lst[ll][-5:-3]
        elif campaign=='MARCUS':
            legnum=lst[ll][-4]
            
        (timem,[lat0,lon0],timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['lat','lon'])
        lat0[lat0<-9000]=np.nan
        lon0[lon0<-9000]=np.nan
        # if time expands two years, add 365 days to the second year
        if timem[0]>timem[-1]:
            timem[timem<=timem[-1]]=timem[timem<=timem[-1]]+365
            
        # find the days related to the ship leg
        day = [int(a) for a in timem]
        day = list(set(day))
        day.sort()
        
        # read in CCN
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
        
        lat1=np.interp(t_ccn1,timem,lat0)
        lon1=np.interp(t_ccn1,timem,lon0)
        lat5=np.interp(t_ccn5,timem,lat0)
        lon5=np.interp(t_ccn5,timem,lon0)
        # separate into latitude bins
        for bb in range(latlen):
            idx = np.logical_and(lat1>=latmin[bb], lat1<latmax[bb])
            ccn3_o[bb]=np.hstack((ccn3_o[bb],ccn1o[idx]))  
            idx = np.logical_and(lat5>=latmin[bb], lat5<latmax[bb])
            ccn5_o[bb]=np.hstack((ccn5_o[bb],ccn5o[idx]))  
    
    #%%
    for bb in range(latlen):
        ccn3_o[bb] = ccn3_o[bb][~np.isnan(ccn3_o[bb])]
        ccn5_o[bb] = ccn5_o[bb][~np.isnan(ccn5_o[bb])]
    
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
        
    figname = figpath_ship_statistics+'percentile_lat_CCN_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(ccn3_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(ccn3_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax1.tick_params(color='k',labelsize=15)
    #ax1.set_yscale('log')
    ax1.set_xlim(-1,latlen)
    ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    ax1.set_xticklabels([])
    # plot temporal lines for label
    ax1.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
        
    ax2.boxplot(ccn5_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax2.boxplot(ccn5_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax2.tick_params(color='k',labelsize=15)
    #ax2.set_yscale('log')
    ax2.set_xlim(-1,latlen)
    ax2.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    ax2.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
    # plot temporal lines for label
    ax2.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax2.plot([],c=color_model[mm],label=Model_List[mm])
        
    # ax1.legend(loc='upper right', fontsize='large')
    ax2.legend(loc='upper right', fontsize='x-large')
    
    # supersaturation
    fig.text(0.08,0.98,'SS='+str(SS1)+'%')
    fig.text(0.08,0.47,'SS='+str(SS5)+'%')
        
    ax2.set_xlabel('Latitude',fontsize=16)
    ax1.set_title('CCN Number Concentration (cm$^{-3}$)',fontsize=17)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
    
