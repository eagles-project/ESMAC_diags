# plot timeseries of basic meteorological variables along ship track

import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from read_ship import read_marmet
from read_ARMdata import read_met
from read_netcdf import read_E3SM
from time_format_change import yyyymmdd2cday,  cday2mmdd

# define function of averaging in time for faster plotting
def avg_time(time0,data0,time):
    data0[data0<-900]=np.nan
    if len(data0)!=len(time0):
        error
    data = np.full((len(time)),np.nan)
    dt=(time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0>=time[tt]-dt,time0<=time[tt]+dt)
        data[tt]=np.nanmean(data0[idx],axis=0)
    return(data)

#%% settings

from settings import campaign, Model_List, color_model, shipmetpath, E3SM_ship_path, figpath_ship_timeseries

import os
if not os.path.exists(figpath_ship_timeseries):
    os.makedirs(figpath_ship_timeseries)


lst = glob.glob(E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[0]+'_shipleg*.nc')
lst.sort()

for ll in range(len(lst)):
    
    #%% for MAGIC, read each ship leg
    if campaign=='MAGIC':
        legnum=lst[ll][-5:-3]
        filenameo = shipmetpath+'marmet'+legnum+'.txt'
        (shipdata,shipvarlist) = read_marmet(filenameo)
        year=[a[1] for a in shipdata]
        month=[a[2] for a in shipdata]
        day=[a[3] for a in shipdata]
        hh=[int(a[4]) for a in shipdata]
        mm=[int(a[5]) for a in shipdata]
        ss=[int(a[6]) for a in shipdata]
        yyyymmdd = [year[i]+month[i]+day[i] for i in range(len(year))]   # yyyymmdd
        # get time in calendar day
        time = np.array(hh)/24. + np.array(mm)/1440. + np.array(ss)/86400. 
        time = np.array([time[i] + yyyymmdd2cday(yyyymmdd[i],'noleap') for i in range(len(time))])
        if time[-1]<time[0]:
            time[time<=time[-1]]=time[time<=time[-1]]+365
        
        
        # get variables
        ps=np.array([float(a[shipvarlist.index('bp')]) for a in shipdata])    
        rh=np.array([float(a[shipvarlist.index('rh')]) for a in shipdata]) 
        ta=np.array([float(a[shipvarlist.index('ta')]) for a in shipdata]) 
        rain=np.array([float(a[shipvarlist.index('org')]) for a in shipdata]) 
        lat=np.array([float(a[7]) for a in shipdata])
    
        lat[lat==-999]=np.nan
        ps[ps==-999]=np.nan
        rh[rh==-999]=np.nan
        ta[ta==-999]=np.nan
        rain[rain==-999]=np.nan
        
        # rain needs to be averaged into 1-hr timewindow for comparison with model
        time1hr = np.arange(int(time[0]*24),int(time[-1]*24)+1)/24.
        rain2=avg_time(time,rain,time1hr)
        
        # rain rate in leg 19 are unrealistic. mask all data
        if legnum=='19':
            rain=rain*np.nan
            rain2=rain2*np.nan
    
    #%% for MARCUS, read files for each vessel trip
    elif campaign=='MARCUS':
        legnum=lst[ll][-4]
        if legnum=='1':
            startdate='2017-10-30'
            enddate='2017-12-02'
        elif legnum=='2':
            startdate='2017-12-13'
            enddate='2018-01-11'
        elif legnum=='3':
            startdate='2018-01-16'
            enddate='2018-03-04'
        elif legnum=='4':
            startdate='2018-03-09'
            enddate='2018-03-22'
            
        year=[startdate[0:4]]
        
        cday1=yyyymmdd2cday(startdate,'noleap')
        cday2=yyyymmdd2cday(enddate,'noleap')
        if startdate[0:4]!=enddate[0:4]:
            cday2=cday2+365  # cover two years

        time=np.empty(0)
        lon=np.empty(0)
        lat=np.empty(0)
        ta=np.empty(0)
        rh=np.empty(0)
        ps=np.empty(0)
        for cc in range(cday1,cday2+1):
            if cc<=365:
                yyyymmdd=startdate[0:4]+cday2mmdd(cc)
            else:
                yyyymmdd=enddate[0:4]+cday2mmdd(cc-365)
                
            lst0 = glob.glob(shipmetpath+'maraadmetX1.b1.'+yyyymmdd+'*')
            (time0,lon0,timeunit,lonunit,lon_long_name)=read_met(lst0[0],'lon')
            (time0,lat0,timeunit,lonunit,lon_long_name)=read_met(lst0[0],'lat')
            (time0,ta1,timeunit,taunit,ta_long_name)=read_met(lst0[0],'air_temperature_port')
            (time0,ta2,timeunit,taunit,ta_long_name)=read_met(lst0[0],'air_temperature_starboard')
            (time0,rh1,timeunit,rhunit,rh_long_name)=read_met(lst0[0],'relative_humidity_port')
            (time0,rh2,timeunit,rhunit,rh_long_name)=read_met(lst0[0],'relative_humidity_starboard')
            (time0,ps0,timeunit,psunit,ps_long_name)=read_met(lst0[0],'atmospheric_pressure')
            
            time = np.hstack((time, time0/86400. + cc))
            lat = np.hstack((lat,lat0))
            lon = np.hstack((lon,lon0))
            ta = np.hstack((ta,(ta1+ta2)/2))
            rh = np.hstack((rh,(rh1+rh2)/2))
            ps = np.hstack((ps,ps0))
            
        ps[ps<=-999]=np.nan
        rh[rh<=-999]=np.nan
        ta[ta<=-999]=np.nan
        lat[lat<=-999]=np.nan
        
        
    #%% read in model
    nmodels=len(Model_List)
    T_m = list()
    RH_m = list()
    ps_m = list()
    rain_m = list()
    for mm in range(nmodels):
        filenamem = E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg'+legnum+'.nc'
    
        (timem,varm,timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['T','RELHUM','PS','PRECT','TREFHT'])
    
        T_m.append(varm[0]-273.16)
        RH_m.append(varm[1])
        ps_m.append(varm[2]*0.01)
        rain_m.append(varm[3]*3600*1000)
        
    # change the unit
    varmunit[0]='C'
    varmunit[2]='hPa'
    varmunit[3]='mm/hr'
    varmlongname[3]='Rainrate'
        
    if len(time)!=len(timem):
        print('ERROR: the model shiptrack data should have the same time as in obs. check!')
        error
    
    #%% make plot
        
    figname = figpath_ship_timeseries+'timeseries_met_'+campaign+'_ship'+legnum+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax0,ax1,ax2,ax3) = plt.subplots(4,1,figsize=(8,7))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax0.plot(time,lat,color='k')
    ax0.tick_params(color='k',labelsize=12)
    
    ax1.plot(time,ta,color='k',linewidth=1,label='OBS')
    for mm in range(nmodels):
        ax1.plot(time, T_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax1.tick_params(color='k',labelsize=12)
    
    ax2.plot(time,rh,color='k',linewidth=1,label='OBS')
    for mm in range(nmodels):
        ax2.plot(time, RH_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax2.tick_params(color='k',labelsize=12)
    
    ax3.plot(time,ps,color='k',linewidth=1,label='OBS')
    for mm in range(nmodels):
        ax3.plot(time, ps_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    ax3.tick_params(color='k',labelsize=12)
    ax3.set_ylim(int(min(np.nanmin(ps), np.nanmin(ps_m[0]-2))),int(max(np.nanmax(ps),np.nanmax(ps_m[0])))+2)
    
    # # ax4.plot(time,rain,color='k',linewidth=1,label='OBS')
    # ax4.plot(time1hr,rain2,color='k',linewidth=1,label='OBS')  # 1hr averaged rain. 
    # for mm in range(nmodels):
    #     ax4.plot(time, rain_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax4.tick_params(color='k',labelsize=12)
    
    ax0.set_xticklabels([])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xlabel('Calenday Day in '+year[0],fontsize=14)
    ax0.set_title('latitude',fontsize=15)
    ax1.set_title(varmlongname[0]+' ('+varmunit[0]+')',fontsize=15)
    ax2.set_title(varmlongname[1]+' ('+varmunit[1]+')',fontsize=15)
    ax3.set_title(varmlongname[2]+' ('+varmunit[2]+')',fontsize=15)
    # ax4.set_title(varmlongname[3]+' ('+varmunit[3]+')',fontsize=15)
    
    ax3.legend(loc='lower right', shadow=False, fontsize='large')
    
    fig.text(.1, .999,'trip # '+legnum, fontsize=15)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    