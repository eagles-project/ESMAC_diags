"""
# plot ship-track meteorological variables binned by different latitudes
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_ship import read_marmet
from ..subroutines.read_ARMdata import read_met
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.time_format_change import yyyymmdd2cday,  cday2mmdd

def run_plot(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    shipmetpath = settings['shipmetpath']
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
    T_o = list()
    RH_o = list()
    ps_o = list()
    # rain_o = list()
    for bb in range(latlen):
        T_o.append(np.empty(0))
        RH_o.append(np.empty(0))
        ps_o.append(np.empty(0))
        # rain_o.append(np.empty(0))
    
    if campaign=='MAGIC':
        lst = glob.glob(shipmetpath+'marmet*.txt')
        lst.sort()
        for ll in range(len(lst)):
            legnum=lst[ll][-6:-4]
            
            filenameo = shipmetpath+'marmet'+legnum+'.txt'
            (shipdata,shipvarlist) = read_marmet(filenameo)        
            
            # get variables
            lat=np.array([float(a[shipvarlist.index('lat')]) for a in shipdata]) 
            lon=np.array([float(a[shipvarlist.index('lon')]) for a in shipdata]) 
            sst=np.array([float(a[shipvarlist.index('ssst')]) for a in shipdata]) 
            ps=np.array([float(a[shipvarlist.index('bp')]) for a in shipdata])    
            rh=np.array([float(a[shipvarlist.index('rh')]) for a in shipdata]) 
            ta=np.array([float(a[shipvarlist.index('ta')]) for a in shipdata]) 
            rain=np.array([float(a[shipvarlist.index('org')]) for a in shipdata]) 
        
            lat[lat==-999]=np.nan
            lon[lon==-999]=np.nan
            sst[sst==-999]=np.nan
            ps[ps==-999]=np.nan
            rh[rh==-999]=np.nan
            ta[ta==-999]=np.nan
            rain[rain==-999]=np.nan
            
            # rain rate in leg 19 are unrealistic. mask all data
            if legnum=='19':
                rain=rain*np.nan
        
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat>=latmin[bb], lat<latmax[bb])
                T_o[bb]=np.hstack((T_o[bb],ta[idx]))
                RH_o[bb]=np.hstack((RH_o[bb],rh[idx]))
                ps_o[bb]=np.hstack((ps_o[bb],ps[idx]))
                # rain_o[bb]=np.hstack((rain_o[bb],rain[idx]))
                
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
            (time0,ta1,timeunit,taunit,ta_long_name)=read_met(lst0[0],'air_temperature_port')
            (time0,ta2,timeunit,taunit,ta_long_name)=read_met(lst0[0],'air_temperature_starboard')
            (time0,rh1,timeunit,rhunit,rh_long_name)=read_met(lst0[0],'relative_humidity_port')
            (time0,rh2,timeunit,rhunit,rh_long_name)=read_met(lst0[0],'relative_humidity_starboard')
            (time0,ps,timeunit,psunit,ps_long_name)=read_met(lst0[0],'atmospheric_pressure')
            
            ta = (ta1+ta2)/2
            rh = (rh1+rh2)/2
            
            ps[ps<=-999]=np.nan
            rh[rh<=-999]=np.nan
            ta[ta<=-999]=np.nan
            
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat>=latmin[bb], lat<latmax[bb])
                T_o[bb]=np.hstack((T_o[bb],ta[idx]))
                RH_o[bb]=np.hstack((RH_o[bb],rh[idx]))
                ps_o[bb]=np.hstack((ps_o[bb],ps[idx]))
            
    
        
    #%% read in model
    nmodels=len(Model_List)
    T_m = list()
    RH_m = list()
    ps_m = list()
    rain_m = list()
    for mm in range(nmodels):
        
        # initialize variables by latitude bins
        T_tmp = list()
        RH_tmp = list()
        ps_tmp = list()
        rain_tmp = list()
        for bb in range(latlen):
            T_tmp.append(np.empty(0))
            RH_tmp.append(np.empty(0))
            ps_tmp.append(np.empty(0))
            rain_tmp.append(np.empty(0))
        
        lst = glob.glob(E3SM_ship_path+'Ship_vars_'+campaign+'_'+Model_List[mm]+'_shipleg*.nc')
        lst.sort()
        for filenamem in lst:    
            (timem,varm,timeunitm,varmunit,varmlongname)=read_E3SM(filenamem,['T','RELHUM','PS','PRECT','lat','lon'])
            for ii in range(len(varm)):
                varm[ii][varm[ii]<-9000] = np.nan
                
            lat0=varm[4]
            lon0=varm[5]
            
            # separate into latitude bins
            for bb in range(latlen):
                idx = np.logical_and(lat0>=latmin[bb], lat0<latmax[bb])
                T_tmp[bb]=np.hstack((T_tmp[bb],varm[0][idx]-273.16))
                RH_tmp[bb]=np.hstack((RH_tmp[bb],varm[1][idx]))
                ps_tmp[bb]=np.hstack((ps_tmp[bb],varm[2][idx]*0.01))
                rain_tmp[bb]=np.hstack((rain_tmp[bb],varm[3][idx]*3600*1000))
            
        T_m.append(T_tmp)
        RH_m.append(RH_tmp)
        ps_m.append(ps_tmp)
        rain_m.append(rain_tmp)
        
    # change the unit
    varmunit[0]='C'
    varmunit[2]='hPa'
    varmunit[3]='mm/hr'
    varmlongname[3]='Rainrate'
    
    #%% calculate the mean and standard error for each bin
    for bb in range(latlen):
        T_o[bb] = T_o[bb][~np.isnan(T_o[bb])]
        RH_o[bb] = RH_o[bb][~np.isnan(RH_o[bb])]
        ps_o[bb] = ps_o[bb][~np.isnan(ps_o[bb])]
        # rain_o[bb] = rain_o[bb][~np.isnan(rain_o[bb])]
        
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
        
    figname = figpath_ship_statistics+'percentile_lat_met_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(8,6))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
    ax1.boxplot(T_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax1.boxplot(T_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax1.tick_params(color='k',labelsize=15)
    # ax1.set_yscale('log')
    ax1.set_xlim(-1,latlen)
    ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    # plot temporal lines for label
    ax1.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax1.plot([],c=color_model[mm],label=Model_List[mm])
        
    ax2.boxplot(RH_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax2.boxplot(RH_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax2.tick_params(color='k',labelsize=15)
    # ax2.set_yscale('log')
    ax2.set_xlim(-1,latlen)
    ax2.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    # plot temporal lines for label
    ax2.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax2.plot([],c=color_model[mm],label=Model_List[mm])
        
    ax3.boxplot(ps_o,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax3.boxplot(ps_m[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax3.tick_params(color='k',labelsize=15)
    # ax3.set_yscale('log')
    ax3.set_xlim(-1,latlen)
    ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    # plot temporal lines for label
    ax3.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
    # ax4.boxplot(rain_o,whis=(5,95),showmeans=False,showfliers=False,
    #             positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
    #             boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
    #             medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
    #             vert=True, patch_artist=True)    # need patch_artist to fill color in box
    # for mm in range(nmodels):
    #     c = color_model[mm]
    #     ax4.boxplot(rain_m[mm],whis=(5,95),showmeans=False,showfliers=False,
    #             positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
    #             boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
    #             medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
    #             vert=True, patch_artist=True)    # need patch_artist to fill color in box
    # ax4.tick_params(color='k',labelsize=12)
    # # ax4.set_yscale('log')
    # ax4.set_xlim(-1,latlen)
    # ax4.set_xticks(np.arange(-.25,latlen,2))
    # # plot temporal lines for label
    # ax4.plot([],c='k',label='OBS')
    # for mm in range(nmodels):
    #     ax4.plot([],c=color_model[mm],label=Model_List[mm])
        
    
    ax3.set_xlabel('Latitude',fontsize=16)
    
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
    ax1.set_title(varmlongname[0]+' ('+varmunit[0]+')',fontsize=17)
    ax2.set_title(varmlongname[1]+' ('+varmunit[1]+')',fontsize=17)
    ax3.set_title(varmlongname[2]+' ('+varmunit[2]+')',fontsize=17)
    # ax4.set_title(varmlongname[3]+' ('+varmunit[3]+')',fontsize=15)
    
    ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
