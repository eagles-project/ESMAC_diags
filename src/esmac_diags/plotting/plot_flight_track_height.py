"""
# plot aircraft information
# plot 1: plot flight track location (lat/lon) with height in color
# plot 2: plot timeseries of flight height with cloud and CVI flags
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ARMdata import read_cvi_aceena
from ..subroutines.specific_data_treatment import lwc2cflag
from ..subroutines.time_format_change import hhmmss2sec
from ..subroutines.read_aircraft import read_iwg1, read_cvi_hiscale, read_RF_NCAR


def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    lat0 = settings['lat0']
    lon0 = settings['lon0']
    figpath_aircraft_timeseries = settings['figpath_aircraft_timeseries']

    IOP = settings.get('IOP', None)
    iwgpath = settings.get('iwgpath', None)
    cvipath = settings.get('cvipath', None)
    RFpath = settings.get('RFpath', None)

    #%% other settings

    if not os.path.exists(figpath_aircraft_timeseries):
        os.makedirs(figpath_aircraft_timeseries)
        

    #%% find all flight data

    if campaign=='HISCALE':
        lst = glob.glob(iwgpath+'*a2.txt')
        lst.sort()
        if IOP=='IOP1':
            lst=lst[0:17]
        elif IOP=='IOP2':
            lst=lst[17:]
        elif IOP[0:4]=='2016':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
    elif campaign=='ACEENA':
        lst = glob.glob(iwgpath+'*a2.txt')
        lst.sort()
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
    elif campaign in ['CSET', 'SOCRATES']:
        lst = glob.glob(RFpath+'RF*.PNI.nc')
    else:
        raise ValueError('campaign name is not recognized: '+campaign)
    lst.sort()

    #%% read in data and make plot
    for filename in lst:
        
        # get date info:        
        fname=filename.split('.')
        
        #%% read in flight data (for HISCALE and ACEENA)
        if campaign in ['HISCALE', 'ACEENA']:
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
                        raise ValueError('time dimension is incosistent')
                elif len(filename_c)==0:
                    time_cvi=time
                    cvi_inlet=np.nan*np.empty([len(time)])
                else:
                    raise ValueError('find too many files: '+filename_c)
                    
            elif campaign=='ACEENA':
                filename_c=glob.glob(cvipath+'enaaafinletcviF1.c1.'+date[0:8]+'*.nc')
                filename_c.sort()
                # read in data
                if len(filename_c)==1:
                    (time_c,lon_c,lat_c,alt_c,timeunit_c,cvimode,cvi_inlet,enhance_factor,dilution_factor)=read_cvi_aceena(filename_c[0])
                    if date=='20180216a':
                        time_c=np.insert(time_c,1403,(time_c[1402]+time_c[1403])/2)
                        cvi_inlet=np.insert(cvi_inlet,1403,cvi_inlet[1403])
                    if all(time_c==time)==False:
                        raise ValueError('time dimension is incosistent')
                elif len(filename_c)==0:
                    time_cvi=time
                    cvi_inlet=np.nan*np.empty([len(time)])
                else:
                    raise ValueError('find too many files: '+filename_c)
                # cvi_inlet[cvi_inlet==-9]=1  # if cvi_inlet is unfunctional, use fims as good data
            
            else:
                raise ValueError('do not recognize this campaign: '+campaign)
        
        #%% read in flight data (for CSET and SOCRATES)
        elif campaign in ['CSET', 'SOCRATES']:
            date=fname[-4]
            print('input data for '+date)
            (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'ALT')
            (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LAT')
            (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LON')
            (time,lwc,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename,'PLWCC')
            # lon[lon<0]=lon[lon<0]+360
            # calculate cloud flag based on LWC
            cldflag=lwc2cflag(lwc,lwcunit)
            if campaign=='SOCRATES':
                (time,cvi_inlet,timeunit,cviunit,cvilongname,cellsize,cellunit)=read_RF_NCAR(filename,'CVINLET')
            else:
                cvi_inlet=np.nan*np.empty([len(time)])
        
        
        #%% plot flight tracks:
        lat[lat<-9000]=np.nan
        lon[lon<-9000]=np.nan
        height[height<-9000]=np.nan
        
        # change longitude to [-180, 180]
        if lon0>180:
            lon0=lon0-360
            
        try:
        #     os.environ['PROJ_LIB'] = r'c:\Users\tang357\Anaconda3\pkgs\basemap-1.3.0-py38ha7665c8_0\Library\share'
        #     from mpl_toolkits.basemap import Basemap
        #     figname = figpath_aircraft_timeseries + 'flighttrack_'+campaign+'_'+date+'.png'
        #     print('plot flight track to '+figname)
        #     fig,ax = plt.subplots(figsize=(8,5))   # figsize in inches
        #     plt.tight_layout(pad=0.1, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        #     if campaign in ['CSET', 'SOCRATES']:
        #         m = Basemap(llcrnrlon=min(np.floor(min(lon)),np.floor(lon0))-2,llcrnrlat=min(np.floor(min(lat)),np.floor(lat0))-2,\
        #             urcrnrlon=max(np.ceil(max(lon)),np.ceil(lon0))+2,urcrnrlat=max(np.ceil(max(lat)),np.ceil(lat0))+2,\
        #             resolution='l',rsphere=(6378137.00,6356752.3142),projection='lcc',lat_0=np.min(lat),lon_0=np.min(lon)) #,lat_ts=5.)
        #         m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0])
        #         m.drawmeridians(np.arange(-180,180,5),labels=[0,0,0,1])
        #         m.drawcoastlines()
        #         m.fillcontinents()
        #     elif campaign=='HISCALE':
        #         m = Basemap(llcrnrlon=-99,llcrnrlat=35,urcrnrlon=-95,urcrnrlat=38,\
        #             resolution='l',rsphere=(6378137.00,6356752.3142),projection='lcc',lat_0=lat0,lon_0=lon0) #,lat_ts=5.)
        #         m.drawparallels(np.arange(30,40,1),labels=[1,0,0,0])
        #         m.drawmeridians(np.arange(-110,-90,1),labels=[0,0,0,1])
        #         m.drawstates()
        #         x2,y2=m(lon0,lat0)
        #         m.scatter(x2,y2,s=100,marker='*',color='k')
        #     elif campaign=='ACEENA':
        #         m = Basemap(llcrnrlon=-30,llcrnrlat=37,urcrnrlon=-25,urcrnrlat=41,\
        #             resolution='l',rsphere=(6378137.00,6356752.3142),projection='lcc',lat_0=lat0,lon_0=lon0) #,lat_ts=5.)
        #         m.drawparallels(np.arange(30,42,1),labels=[1,0,0,0])
        #         m.drawmeridians(np.arange(-30,-20,1),labels=[0,0,0,1])
        #         m.drawcoastlines()
        #         m.fillcontinents()
        #         x2,y2=m(lon0,lat0)
        #         m.scatter(x2,y2,s=100,marker='*',color='k')
        #     x, y = m(lon,lat)
        #     h=m.scatter(x,y,s=1,c=height,cmap='jet')
        #     ax.set_title('Flight track '+date,fontsize=15)
        #     cbar=fig.colorbar(h)
        # except:
            figname = figpath_aircraft_timeseries + 'flighttrack_'+campaign+'_'+date+'.png'
            print('plot flight track to '+figname)
            fig,ax = plt.subplots(figsize=(8,5))   # figsize in inches
            # plot the location of the campaign site:
            ax.plot([lon0,lon0],[lat0-50, lat0+50],':',color=[.8,.8,.8])
            ax.plot([lon0-50, lon0+50],[lat0,lat0],':',color=[.8,.8,.8])
            # plot flight track
            h=ax.scatter(lon,lat,s=1,c=height,cmap='jet',vmin=0,vmax=max(height))  #vmin/vmax: color range
            ax.set_xlim(min(np.floor(min(lon)),np.floor(lon0)), max(np.ceil(max(lon)),np.ceil(lon0)))
            ax.set_ylim(min(np.floor(min(lat)),np.floor(lat0)), max(np.ceil(max(lat)),np.ceil(lat0)))
            ax.tick_params(color='k',labelsize=14)
            ax.set_xlabel('longitude',fontsize=14)
            ax.set_ylabel('latitude',fontsize=14)
            ax.set_title('Flight track '+date,fontsize=15)
            cbar=fig.colorbar(h)
            fig.text(0.81,0.91, 'm MSL')
        except:
            raise ValueError("cannot make flight track plot")
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()
            
        #%% plot flight height and flag/leg timeseries
        figname = figpath_aircraft_timeseries + 'flightheight_'+campaign+'_'+date+'.png'
        print('plot flight height timeseries to '+figname)
        
        fig,ax1 = plt.subplots(figsize=(8,2))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
        h11=ax1.plot(time/3600,height/1000,color='k',linewidth=1)
        if campaign in ['HISCALE', 'ACEENA']:
            for ll in range(1,max(legnum)+1):
                idx=legnum==ll
                ax1.plot(time[idx]/3600,height[idx]/1000,color='b',linewidth=2)
        h12=ax1.plot(time/3600,time*0+max(height)*0.00105,color='k',linewidth=.2)
        cvi2=0.0*cvi_inlet
        cvi2[cvi_inlet==1]=np.nan
        cvi2=cvi2+max(height)*0.00105
        h13=ax1.plot(time/3600,cvi2,color='k',linewidth=2)
        h14=ax1.vlines(time[cldflag==1]/3600,0,max(height)*0.0011,color='silver',linewidth=0.1)
        # ax1.set_xlim(time[0]/3600-0.3, time[-1]/3600+0.3)
        ax1.set_ylim(0,max(height)*0.0011)
        ax1.set_ylabel('height (km)',fontsize=12)
        ax1.set_xlabel('time (hour UTC) '+date,fontsize=12)
        ax1.set_title('thin black: flight track. blue: flight legs. gray vertical lines: cloud flag. thick black: CVI mode', fontsize=10)
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        plt.close()
    
