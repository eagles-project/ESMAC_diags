"""
# plot ship-track meteorological variables binned by different latitudes
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from ..subroutines.read_ship import read_marmet
from ..subroutines.read_ARMdata import read_mwr, read_met
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import yyyymmdd2cday, cday2mmdd

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    shipmetpath = settings['shipmetpath']
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
    
    
    #%% read in observation
    
    # initialize variables by latitude bins
    lwp_o = list()
    rain_o = list()
    for bb in range(latlen):
        lwp_o.append(np.empty(0))
        rain_o.append(np.empty(0))
    
    if campaign=='MAGIC':
        lst = glob.glob(shipmetpath+'marmet*.txt')
        lst.sort()
    
        for ll in range(len(lst)):
            legnum=lst[ll][-6:-4]
            
            # read in MET
            filenameo = shipmetpath+'marmet'+legnum+'.txt'
            (shipdata,shipvarlist) = read_marmet(filenameo)        
            # get variables
            lat=np.array([float(a[shipvarlist.index('lat')]) for a in shipdata]) 
            lon=np.array([float(a[shipvarlist.index('lon')]) for a in shipdata]) 
            rain=np.array([float(a[shipvarlist.index('org')]) for a in shipdata]) 
            lat[lat==-999]=np.nan
            lon[lon==-999]=np.nan
            rain[rain==-999]=np.nan
            # rain rate in leg 19 are unrealistic. mask all data
            if legnum=='19':
                rain=rain*np.nan
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat>=latmin[bb], lat<latmax[bb])
                rain_o[bb]=np.hstack((rain_o[bb],rain[idx]))
                
            # read in MWR
            t_lwp=np.empty(0)
            lwp=np.empty(0)
            # find the days related to the ship leg
            year=[a[1] for a in shipdata]
            month=[a[2] for a in shipdata]
            day=[a[3] for a in shipdata]
            hh=[int(a[4]) for a in shipdata]
            mm=[int(a[5]) for a in shipdata]
            ss=[int(a[6]) for a in shipdata]
            
            yyyymmdd = [year[i]+month[i]+day[i] for i in range(len(year))]   # yyyymmdd
            time0 = np.array(hh)/24. + np.array(mm)/1440. + np.array(ss)/86400.
            time0 = np.array([time0[i] + yyyymmdd2cday(yyyymmdd[i],'noleap') for i in range(len(time0))])
        
            sday = [year[a]+month[a]+day[a] for a in range(len(mm))]
            sday = list(set(sday))
            sday.sort()
            for dd in sday:
                filenameo = glob.glob(shipmwrpath+'magmwrret1liljclouM1.s2.'+dd+'.*')
                if len(filenameo)==0:
                    continue  # some days may be missing
                (time,obs,timeunit,lwpunit,lwpflag)=read_mwr(filenameo[0],'be_lwp')
                t_lwp=np.hstack((t_lwp, yyyymmdd2cday(dd,'noleap')+time/86400))
                obs[obs<-9000]=np.nan
                lwp=np.hstack((lwp, obs))
            # if no obs available, fill one data with NaN
            if len(t_lwp)==0:
                t_lwp=[time0[0],time0[1]]
                lwp=np.full((2),np.nan)
            # if time expands two years, add 365 days to the second year
            if t_lwp[0]>t_lwp[-1]:
                t_lwp[t_lwp<=t_lwp[-1]]=t_lwp[t_lwp<=t_lwp[-1]]+365
            lat1=np.interp(t_lwp,time0,lat)
            lon1=np.interp(t_lwp,time0,lon)
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat1>=latmin[bb], lat1<latmax[bb])
                lwp_o[bb]=np.hstack((lwp_o[bb],lwp[idx]))  
            
    elif campaign=='MARCUS':
        
        startdate='2017-10-30'
        enddate='2018-03-22'
        cday1=yyyymmdd2cday(startdate,'noleap')
        cday2=yyyymmdd2cday(enddate,'noleap')
        if startdate[0:4]!=enddate[0:4]:
            cday2=cday2+365  # cover two years
                
        for cc in range(cday1,cday2+1):
            if cc<=365:
                yyyymmdd=startdate[0:4]+cday2mmdd(cc)
            else:
                yyyymmdd=enddate[0:4]+cday2mmdd(cc-365)
                
            lst0 = glob.glob(shipmetpath+'maraadmetX1.b1.'+yyyymmdd+'*')
            if len(lst0)==0:
                continue
            (time0,lon,timeunit,lonunit,lon_long_name)=read_met(lst0[0],'lon')
            (time0,lat,timeunit,lonunit,lon_long_name)=read_met(lst0[0],'lat')
            
            lat[lat==-999]=np.nan
            lon[lon==-999]=np.nan
            rain=np.array(lon)*np.nan
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat>=latmin[bb], lat<latmax[bb])
                rain_o[bb]=np.hstack((rain_o[bb],rain[idx]))
            
            # read in MWR
            filenameo = glob.glob(shipmwrpath+'marmwrret1liljclouM1.s2.'+yyyymmdd+'.*')
            if len(filenameo)==0:
                continue  # some days may be missing
            (time,lwp,timeunit,lwpunit,lwpflag)=read_mwr(filenameo[0],'be_lwp')
            lwp[lwp<-9000]=np.nan
            # if no obs available, fill one data with NaN
            if len(time)==0:
                continue
            
            lat1=np.interp(time,time0,lat)
            lon1=np.interp(time,time0,lon)
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat1>=latmin[bb], lat1<latmax[bb])
                lwp_o[bb]=np.hstack((lwp_o[bb],lwp[idx]))  
            
                
    #%% read in model
    nmodels=len(Model_List)
    lwp_m = list()
    rain_m = list()
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
            (timem,varm,timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['TGCLDLWP','PRECT','lat','lon'])
            for ii in range(len(varm)):
                varm[ii][varm[ii]<-9000] = np.nan
                
            lat0=varm[2]
            lon0=varm[3]
            
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat0>=latmin[bb], lat0<latmax[bb])
                lwp_tmp[bb]=np.hstack((lwp_tmp[bb],varm[0][idx]*1000))
                rain_tmp[bb]=np.hstack((rain_tmp[bb],varm[1][idx]*3600*1000))
            
        lwp_m.append(lwp_tmp)
        rain_m.append(rain_tmp)
        
    # change the unit
    varmunit[0]='g/m2'
    varmunit[1]='mm/hr'
    varmlongname[0]='LWP'
    varmlongname[1]='Rainrate'
    
    #%% calculate the mean and standard error for each bin
    mean_lwp_o = np.array([np.nanmean(a) for a in lwp_o])
    mean_rain_o = np.array([np.nanmean(a) for a in rain_o])
    sem_lwp_o = np.array([scipy.stats.sem(a,nan_policy='omit') for a in lwp_o])
    sem_rain_o = np.array([scipy.stats.sem(a,nan_policy='omit') for a in rain_o])
    
    mean_lwp_m = list()
    mean_rain_m = list()
    sem_lwp_m = list()
    sem_rain_m = list()
    for mm in range(nmodels):
        mean_lwp_m.append(np.array([np.nanmean(a) for a in lwp_m[mm]]))
        mean_rain_m.append(np.array([np.nanmean(a) for a in rain_m[mm]]))
        sem_lwp_m.append(np.array([scipy.stats.sem(a,nan_policy='omit') for a in lwp_m[mm]]))
        sem_rain_m.append(np.array([scipy.stats.sem(a,nan_policy='omit') for a in rain_m[mm]]))
    
    #%% make plot
        
    figname = figpath_ship_statistics+'composite_LWPrain_bylat_'+campaign+'_all.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax1.plot(latbin,mean_lwp_o,color='k',linewidth=1,label='OBS')
    ax1.fill_between(latbin, mean_lwp_o-sem_lwp_o, mean_lwp_o+sem_lwp_o, alpha=0.5, facecolor='gray')
    for mm in range(nmodels):
        ax1.plot(latbin, mean_lwp_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax1.tick_params(color='k',labelsize=15)
    
    ax2.plot(latbin,mean_rain_o,color='k',linewidth=1,label='OBS') 
    ax2.fill_between(latbin,mean_rain_o-sem_rain_o, mean_rain_o+sem_rain_o, alpha=0.5, facecolor='gray')
    for mm in range(nmodels):
        ax2.plot(latbin, mean_rain_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax2.tick_params(color='k',labelsize=15)
    
    ax1.set_xticks(np.arange(np.floor(latbin[0]),int(latbin[-1])+1, dlat*2))
    ax2.set_xticks(np.arange(np.floor(latbin[0]),int(latbin[-1])+1, dlat*2))
    ax1.set_xticklabels([])
    ax2.set_xlabel('Latitude',fontsize=16)
    ax1.set_title(varmlongname[0]+' ('+varmunit[0]+')',fontsize=17)
    ax2.set_title(varmlongname[1]+' ('+varmunit[1]+')',fontsize=17)
    
    ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
