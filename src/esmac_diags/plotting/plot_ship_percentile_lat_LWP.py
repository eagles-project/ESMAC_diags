"""
# plot ship-track liquid water path binned by different latitudes
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_mwr
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import yyyymmdd2cday, cday2mmdd

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    shipmwrpath = settings['shipmwrpath']
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
    lwp_m = list()
    rain_m = list()
    lonm = np.empty(0)
    latm = np.empty(0)
    timem = np.empty(0)
    
    for mm in range(nmodels):
        
        # initialize variables by latitude bins
        lwp_tmp = list()
        rain_tmp = list()
        for bb in range(latlen):
            lwp_tmp.append(np.empty(0))
            rain_tmp.append(np.empty(0))
            
        lst = glob.glob(E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg*.nc')
        lst.sort()   
        for ll in range(len(lst)):
            filenamem = lst[ll]
            (timem0,varm,timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['TGCLDLWP','PRECT','lat','lon'])
            for ii in range(len(varm)):
                varm[ii][varm[ii]<-9000] = np.nan
                
            lat0=varm[2]
            lon0=varm[3]
            
            # save all time and lon/lat to interpolate observation
            if mm==0:
                timem=np.hstack((timem,timem0))
                lonm=np.hstack((lonm,lon0))
                latm=np.hstack((latm,lat0))
            
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat0>=latmin[bb], lat0<latmax[bb])
                lwp_tmp[bb]=np.hstack((lwp_tmp[bb],varm[0][idx]*1000))
            
        lwp_m.append(lwp_tmp)
        
    # change the unit
    varmunit[0]='g/m2'
    varmunit[1]='mm/hr'
    varmlongname[0]='LWP'
    varmlongname[1]='Rainrate'
    
    idx_neg = np.where((timem[1:]-timem[:-1])<0)  # find the day from Dec 31 to Jan 1.
    if len(idx_neg)!=1:
        raise ValueError('this code is only designed for two continuous years. If more than two years please edit this part')
    else:
        timem[idx_neg[0][0]+1:] = timem[idx_neg[0][0]+1:]+365
    
    #%% read in observation
    
    # initialize variables by latitude bins
    lwp_o = list()
    rain_o = list()
    for bb in range(latlen):
        lwp_o.append(np.empty(0))
        rain_o.append(np.empty(0))
    
    timeo = np.empty(0)
    lwpo = np.empty(0)
    
    if campaign=='MAGIC':
        startdate='2012-10-05'
        enddate='2013-09-26'
        # startdate='2013-07-01'
        # enddate='2013-08-31'
        mwrdatastream='magmwrret1liljclouM1'
    elif campaign=='MARCUS':
        startdate='2017-10-30'
        enddate='2018-03-22'
        mwrdatastream='marmwrret1liljclouM1'
        
    cday1=yyyymmdd2cday(startdate,'noleap')
    cday2=yyyymmdd2cday(enddate,'noleap')
    if startdate[0:4]!=enddate[0:4]:
        cday2=cday2+365  # cover two years
            
    for cc in range(cday1,cday2+1):
        if cc<=365:
            yyyymmdd=startdate[0:4]+cday2mmdd(cc)
        else:
            yyyymmdd=enddate[0:4]+cday2mmdd(cc-365)
                
        # read in MWR
        filenameo = glob.glob(shipmwrpath+mwrdatastream+'.s2.'+yyyymmdd+'.*')
        if len(filenameo)==0:
            continue  # some days may be missing
        (time,lwp,timeunit,lwpunit,lwpflag)=read_mwr(filenameo[0],'be_lwp')
        lwp[lwp<-9000]=np.nan
        # if no obs available, fill one data with NaN
        if len(time)==0:
            continue
        
        timeo=np.hstack((timeo,cc+time/86400.))
        lwpo = np.hstack((lwpo,lwp))
        
    lat1=np.interp(timeo,timem,latm)
    # separate into latitude bins
    for bb in range(latlen):
        idx = np.logical_and(lat1>=latmin[bb], lat1<latmax[bb])
        lwp_o[bb]=np.hstack((lwp_o[bb],lwpo[idx]))  
                
    #%% 
    for bb in range(latlen):
        lwp_o[bb] = lwp_o[bb][~np.isnan(lwp_o[bb])]
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    figname = figpath_ship_statistics+'percentile_lat_LWP_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1) = plt.subplots(figsize=(8,2))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(lwp_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(lwp_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax1.tick_params(color='k',labelsize=15)
    # ax1.set_yscale('log')
    ax1.set_xlim(-1,latlen)
    ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    ax1.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
    # plot temporal lines for label
    ax1.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
    ax1.legend(loc='upper left', fontsize='x-large')
        
    
    ax1.set_xlabel('Latitude',fontsize=16)
    ax1.set_title('LWP (g/m$^{2}$)',fontsize=17)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
