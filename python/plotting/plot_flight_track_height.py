# plot aircraft information
# plot 1: plot flight track location (lat/lon) with height in color
# plot 2: plot timeseries of flight height with cloud and CVI flags


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from read_ARMdata import read_cvi_aceena
from time_format_change import yyyymmdd2cday, hhmmss2sec
from read_aircraft import read_iwg1, read_cvi_hiscale


#%% settings

from settings import campaign, iwgpath, cvipath, lat0, lon0, IOP, figpath_aircraft_timeseries

import os
if not os.path.exists(figpath_aircraft_timeseries):
    os.makedirs(figpath_aircraft_timeseries)
    

#%% find all IWG data

lst = glob.glob(iwgpath+'aaf.iwg*.txt')
lst.sort()

if len(lst)==0:
    print('ERROR: cannot find any file at '+iwgpath)
    error
    
# choose files for specific IOP
if campaign=='HISCALE':
    if IOP=='IOP1':
        lst=lst[0:17]
    elif IOP=='IOP2':
        lst=lst[17:]
    elif IOP[0:4]=='2016':
        a=lst[0].split('.'+campaign.lower()+'.')
        lst = glob.glob(a[0]+'*'+IOP+'*')
        lst.sort()
elif campaign=='ACEENA':
    if IOP=='IOP1':
        lst=lst[0:20]
    elif IOP=='IOP2':
        lst=lst[20:]
    elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
        a=lst[0].split('.'+campaign.lower()+'.')
        lst = glob.glob(a[0]+'*'+IOP+'*')
        lst.sort()
else:
    print('ERROR: campaign name is not recognized: '+campaign)
    error

for filename in lst:
    
    # get date info:        
    fname=filename.split('.')
    date=fname[-3]
    if date[-1]=='a':
        flightidx=1
    else:
        flightidx=2
    
    # read in IWG data
    (iwg,iwgvars)=read_iwg1(filename)
    timelen = len(iwg)
    if np.logical_and(campaign=='ACEENA', date=='20180216a'):
        iwg.insert(1403,list(iwg[1403]))
        tstr=iwg[1403][1]
        tstr=tstr[0:-1]+str(int(tstr[-1])-1)
        iwg[1403][1]=tstr
        del iwg[-1]
    # get lat, lon, height, time
    lon=np.empty(timelen)
    lat=np.empty(timelen)
    height=np.empty(timelen)
    time=np.empty(timelen)
    cldflag=np.empty(timelen)
    legnum=np.full(timelen,0)
    for t in range(timelen):
        lat[t]=float(iwg[t][2])
        lon[t]=float(iwg[t][3])
        height[t]=float(iwg[t][4])
        cldflag[t]=int(iwg[t][35])
        legnum[t]=int(iwg[t][-1])
        timestr=iwg[t][1].split(' ')
        time[t]=hhmmss2sec(timestr[1])
    datestr=timestr[0]
    
    # read in CVI
    if campaign=='HISCALE':
        filename_c=glob.glob(cvipath+'CVI_G1_'+date[0:8]+'*R4_HISCALE_001s.ict.txt')
        filename_c.sort()
        # read in data
        if len(filename_c)==1 or len(filename_c)==2:
            (cvi,cvilist)=read_cvi_hiscale(filename_c[flightidx-1])
            time_cvi = cvi[0,:]
            cvi_inlet=cvi[-1,:]
            if all(time_cvi==time)==False:
                print('time is not consistent for CVI')
                error
        elif len(filename_c)==0:
            time_cvi=time
            cvi_inlet=np.nan*np.empty([len(time)])
        else:
            print('find too many files, check: ')
            print(filename_c)
            error
            
    elif campaign=='ACEENA':
        filename_c=glob.glob(cvipath+'enaaafinletcviF1.c1.'+date[0:8]+'*.nc')
        filename_c.sort()
        # read in data
        if len(filename_c)==1:
            (time_c,lon_c,lat_c,alt_c,timeunit_c,cvimode,cvi_inlet)=read_cvi_aceena(filename_c[0])
            if date=='20180216a':
                time_c=np.insert(time_c,1403,(time_c[1402]+time_c[1403])/2)
                cvi_inlet=np.insert(cvi_inlet,1403,cvi_inlet[1403])
            if all(time_c==time)==False:
                print('time is not consistent for CVI')
                error
        elif len(filename_c)==0:
            time_cvi=time
            cvi_inlet=np.nan*np.empty([len(time)])
        else:
            print('find too many files, check: ')
            print(filename_c)
            error
        # cvi_inlet[cvi_inlet==-9]=1  # if cvi_inlet is unfunctional, use fims as good data
    
    else:
        print('ERROR: not recognize this campaign: '+campaign)
        error
    
    #%% plot flight tracks:
    # change longitude to [-180, 180]
    if lon0>180:
        lon0=lon0-360
        
    figname = figpath_aircraft_timeseries + 'flighttrack_'+campaign+'_'+date+'.png'
    print('plot flight track to '+figname)
    fig,ax = plt.subplots(figsize=(8,5))   # figsize in inches
    # plot the location of the campaign site:
    ax.plot([lon0,lon0],[lat0-5, lat0+5],':',color=[.8,.8,.8])
    ax.plot([lon0-5, lon0+5],[lat0,lat0],':',color=[.8,.8,.8])
    # plot flight track
    h=ax.scatter(lon,lat,s=1,c=height,cmap='jet',vmin=0,vmax=4000)  #vmin/vmax: color range
    ax.set_xlim(min(np.floor(min(lon)),np.floor(lon0)), max(np.ceil(max(lon)),np.ceil(lon0)))
    ax.set_ylim(min(np.floor(min(lat)),np.floor(lat0)), max(np.ceil(max(lat)),np.ceil(lat0)))
    ax.tick_params(color='k',labelsize=14)
    ax.set_xlabel('longitude',fontsize=14)
    ax.set_ylabel('latitude',fontsize=14)
    ax.set_title('Flight track '+date,fontsize=15)
    cbar=fig.colorbar(h)
    fig.text(0.81,0.91, 'm MSL')
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    plt.close()
        
    #%% plot flight height and flag/leg timeseries
    figname = figpath_aircraft_timeseries + 'flightheight_timeseries_'+campaign+'_'+date+'.png'
    print('plot flight height timeseries to '+figname)
    
    fig,ax1 = plt.subplots(figsize=(8,2))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    h11=ax1.plot(time/3600,height/1000,color='k',linewidth=1)
    for ll in range(1,max(legnum)+1):
        idx=legnum==ll
        ax1.plot(time[idx]/3600,height[idx]/1000,color='b',linewidth=2)
    h12=ax1.plot(time/3600,time*0+4.4,color='k',linewidth=.2)
    cvi2=0.0*cvi_inlet
    cvi2[cvi_inlet==1]=np.nan
    cvi2=cvi2+4.4
    h13=ax1.plot(time/3600,cvi2,color='k',linewidth=2)
    h14=ax1.vlines(time[cldflag==1]/3600,0,5,color='silver',linewidth=0.1)
    # ax1.set_xlim(time[0]/3600-0.3, time[-1]/3600+0.3)
    ax1.set_ylim(0,4.5)
    ax1.set_ylabel('height (km)',fontsize=12)
    ax1.set_xlabel('time (hour UTC) '+date,fontsize=12)
    ax1.set_title('thin black: flight track. blue: flight legs. gray vertical lines: cloud flag. thick black: CVI mode', fontsize=10)
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    plt.close()