"""
# plot percentile of aerosol number concentration binned by different latitudes
# separated by below-cloud, near-cloud and above-cloud
# for aircraft measurements in CSET or SOCRATES
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_RF_NCAR
from ..subroutines.read_netcdf import read_extractflight
from ..subroutines.specific_data_treatment import lwc2cflag
from ..subroutines.quality_control import qc_mask_takeoff_landing,qc_remove_neg,qc_mask_cloudflag

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    if campaign in ['CSET', 'SOCRATES']:
        RFpath = settings['RFpath']
    else:
        raise ValueError('This code is only for CSET or SOCRATES. check campaign setting: '+campaign)
    
    #%% other settings
    
    plot_method = 'all'  # 'cloud': separate by cloud. 'height': separate by height. 'all': all heights below 5km
    
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
        
    dlat = latbin[1]-latbin[0]
    latmin = latbin-dlat/2
    latmax = latbin+dlat/2
    latlen = len(latbin)
        
    nmodels=len(Model_List)
    
    #%% find files for flight information
    
    lst = glob.glob(E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[0]+'_*.nc')
    lst.sort()
    if len(lst)==0:
        raise ValueError('cannot find any file')
    alldates = [x.split('_')[-1].split('.')[0] for x in lst]
    
    
    #%% define variables by latitude bins below, near and above clouds
        
    cpc_below_lat = []
    cpc_near_lat = []
    cpc_above_lat = []
    uhsas_below_lat = []
    uhsas_near_lat = []
    uhsas_above_lat = []
    for bb in range(latlen):
        cpc_below_lat.append(np.empty(0))
        cpc_near_lat.append(np.empty(0))
        cpc_above_lat.append(np.empty(0))
        uhsas_below_lat.append(np.empty(0))
        uhsas_near_lat.append(np.empty(0))
        uhsas_above_lat.append(np.empty(0))
        
    ncn10_below_lat = []
    ncn10_near_lat = []
    ncn10_above_lat = []
    ncn100_below_lat = []
    ncn100_near_lat = []
    ncn100_above_lat = []
    for mm in range(nmodels):
        ncn10_below_lat.append([np.empty(0) for bb in range(latlen)])
        ncn10_near_lat.append([np.empty(0) for bb in range(latlen)])
        ncn10_above_lat.append([np.empty(0) for bb in range(latlen)])
        ncn100_below_lat.append([np.empty(0) for bb in range(latlen)])
        ncn100_near_lat.append([np.empty(0) for bb in range(latlen)])
        ncn100_above_lat.append([np.empty(0) for bb in range(latlen)])
    
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
        
        #%% read in Models
        cpc100_m = []
        cpc10_m = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,cpc_m,timeunitm,ncn_unit,ncn_longname)=read_extractflight(filename_m,'NCN')
            (timem,heightm,ncnall,timeunitm,ncnall_unit,ncnall_longname)=read_extractflight(filename_m,'NCNall')
            
            cpc100_m.append(np.sum(ncnall[100:,:],0)*1e-6) # #/m3 to #/cm3
            cpc10_m.append(cpc_m*1e-6)    # #/m3 to #/cm3
        
        #%% read in observations
        # note that it is only for CSET and SOCRATES
        lst = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
        if len(lst)==1 or len(lst)==2:  # SOCRATES has two flights in 20180217, choose the later one
            filename=lst[-1]
        else:
            raise ValueError('find no file or too many files: '+lst)
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'ALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LAT')
        (time,lwc,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename,'PLWCC')
        (time,cpc10,timeunit,cpc10unit,cpc10longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCN')
        if campaign=='CSET':
            (time,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_RWOOU')
            # (time,uhsas500,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU500_RWOOU')
        elif campaign=='SOCRATES':
            # there are two variables: CONCU100_CVIU and CONCU100_LWII
            (time,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_LWII')
            # (time,uhsas500,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU500_LWII')
        
        # exclude 30min after takeoff and before landing
        height=qc_mask_takeoff_landing(time,height)
        lat=qc_mask_takeoff_landing(time,lat)
        lwc=qc_mask_takeoff_landing(time,lwc)
        cpc10=qc_mask_takeoff_landing(time,cpc10)
        uhsas100=qc_mask_takeoff_landing(time,uhsas100)
        timem=qc_mask_takeoff_landing(time,timem)
        for mm in range(nmodels):
            cpc10_m[mm]=qc_mask_takeoff_landing(time,cpc10_m[mm])
            cpc100_m[mm]=qc_mask_takeoff_landing(time,cpc100_m[mm])
        
        # calculate cloud flag based on LWC
        cldflag=lwc2cflag(lwc,lwcunit)
        
        cpc10 = qc_mask_cloudflag(cpc10,cldflag)
        cpc10 = qc_remove_neg(cpc10)
        uhsas100 = qc_mask_cloudflag(uhsas100,cldflag)
        uhsas100 = qc_remove_neg(uhsas100)
        
        # if min(lat)<28:
        #     print(np.nanmax(uhsas100[np.logical_and(lat>25,lat<28)]))
        
        #%% separate data by cloud or height
        flag_below = np.zeros(len(time))
        flag_near = np.zeros(len(time))
        flag_above = np.zeros(len(time))
        
        # option 1: separate data by cloud and put in each latitude bin
        if plot_method == 'cloud':
            for ii in range(len(time)):
                if height[ii]>5000:
                    continue   # exclude measurements above 5km
                # check if  there is cloud within 1hr window
                i_start = max(ii-1800, 0)
                i_end = min(ii+1800, len(time))
                if any(cldflag[i_start:i_end]==1):
                    cheight=height[i_start:i_end][cldflag[i_start:i_end]==1]
                    cldmax = np.max(cheight)
                    cldmin = np.min(cheight)
                    if height[ii]<min(cldmin,2000):
                        flag_below[ii]=1
                    elif height[ii]>=cldmin and height[ii]<=cldmax:
                        flag_near[ii]=1
                    elif height[ii]>max(cldmax,1000):
                        flag_above[ii]=1
        
        # option 2: separate data by height
        elif plot_method == 'height':
            for ii in range(len(time)):
                if height[ii]>5000:
                    continue   # exclude measurements above 5km
                # check if  there is cloud within 1hr window
                i_start = max(ii-1800, 0)
                i_end = min(ii+1800, len(time))
                if any(cldflag[i_start:i_end]==1):
                    cheight=height[i_start:i_end][cldflag[i_start:i_end]==1]
                    cldmax = np.max(cheight)
                    cldmin = np.min(cheight)
                    if height[ii]<min(cldmin,2000):
                        flag_below[ii]=1
                    elif height[ii]>max(cldmax,2000):
                        flag_above[ii]=1
                else:
                    if height[ii]<2000:
                        flag_below[ii]=1
                    elif height[ii]>=2000:
                        flag_above[ii]=1
                        
        # option 3: use all heights below 5km
        elif plot_method == 'all':
            for ii in range(len(time)):
                if height[ii]<=5000: # exclude measurements above 5km
                    flag_below[ii]=1
                        
        for bb in range(latlen):
            idx = np.logical_and(lat>=latmin[bb], lat<latmax[bb])
            if any(flag_below[idx]==1):
                cpc_below_lat[bb] = np.hstack((cpc_below_lat[bb], cpc10[idx][flag_below[idx]==1]))
                uhsas_below_lat[bb] = np.hstack((uhsas_below_lat[bb], uhsas100[idx][flag_below[idx]==1]))
                for mm in range(nmodels):
                    ncn10_below_lat[mm][bb] = np.hstack((ncn10_below_lat[mm][bb], cpc10_m[mm][idx][flag_below[idx]==1]))
                    ncn100_below_lat[mm][bb] = np.hstack((ncn100_below_lat[mm][bb], cpc100_m[mm][idx][flag_below[idx]==1]))
            if any(flag_near[idx]==1):
                cpc_near_lat[bb] = np.hstack((cpc_near_lat[bb], cpc10[idx][flag_near[idx]==1]))
                uhsas_near_lat[bb] = np.hstack((uhsas_near_lat[bb], uhsas100[idx][flag_near[idx]==1]))
                for mm in range(nmodels):
                    ncn10_near_lat[mm][bb] = np.hstack((ncn10_near_lat[mm][bb], cpc10_m[mm][idx][flag_near[idx]==1]))
                    ncn100_near_lat[mm][bb] = np.hstack((ncn100_near_lat[mm][bb], cpc100_m[mm][idx][flag_near[idx]==1]))
            if any(flag_above[idx]==1):
                cpc_above_lat[bb] = np.hstack((cpc_above_lat[bb], cpc10[idx][flag_above[idx]==1]))
                uhsas_above_lat[bb] = np.hstack((uhsas_above_lat[bb], uhsas100[idx][flag_above[idx]==1]))
                for mm in range(nmodels):
                    ncn10_above_lat[mm][bb] = np.hstack((ncn10_above_lat[mm][bb], cpc10_m[mm][idx][flag_above[idx]==1]))
                    ncn100_above_lat[mm][bb] = np.hstack((ncn100_above_lat[mm][bb], cpc100_m[mm][idx][flag_above[idx]==1]))
              
    #%% remove nan elements in the observations
    for bb in range(latlen):
        cpc=cpc_below_lat[bb]
        cpc_below_lat[bb] = cpc[~np.isnan(cpc)]
        uhsas=uhsas_below_lat[bb]
        uhsas_below_lat[bb] = uhsas[~np.isnan(uhsas)]
        cpc=cpc_near_lat[bb]
        cpc_near_lat[bb] = cpc[~np.isnan(cpc)]
        uhsas=uhsas_near_lat[bb]
        uhsas_near_lat[bb] = uhsas[~np.isnan(uhsas)]
        cpc=cpc_above_lat[bb]
        cpc_above_lat[bb] = cpc[~np.isnan(cpc)]
        uhsas=uhsas_above_lat[bb]
        uhsas_above_lat[bb] = uhsas[~np.isnan(uhsas)]
    
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    #%% plot separate by cloud
    if plot_method == 'cloud':
        #%% for CPC (>10nm)
        figname = figpath_aircraft_statistics+'percentile_lat_CN10nm_bycldheight_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(8,6))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(cpc_above_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ncn10_above_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax1.plot([],c='k',label='CPC')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax2.boxplot(cpc_near_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax2.boxplot(ncn10_near_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax2.tick_params(color='k',labelsize=15)
        # ax2.set_yscale('log')
        ax2.set_xlim(-1,latlen)
        ax2.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax2.plot([],c='k',label='CPC')
        for mm in range(nmodels):
            ax2.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(cpc_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ncn10_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        # plot temporal lines for label
        ax3.plot([],c='k',label='CPC')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CN (>10nm) # (cm$^{-3}$) '+campaign,fontsize=17)
        fig.text(0.1,0.95,'Above Clouds',fontsize=15)
        fig.text(0.1,0.6,'Near Clouds',fontsize=15)
        fig.text(0.1,0.25,'Below Cloud',fontsize=15)
        
        ax2.legend(loc='upper right', shadow=False, fontsize='x-large')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%% plot for UHSAS (>100nm)
        figname = figpath_aircraft_statistics+'percentile_lat_CN100nm_bycldheight_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(8,6))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(uhsas_above_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ncn100_above_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax1.plot([],c='k',label='UHSAS100')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax2.boxplot(uhsas_near_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax2.boxplot(ncn100_near_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax2.tick_params(color='k',labelsize=15)
        # ax2.set_yscale('log')
        ax2.set_xlim(-1,latlen)
        ax2.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax2.plot([],c='k',label='UHSAS100')
        for mm in range(nmodels):
            ax2.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(uhsas_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ncn100_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        # plot temporal lines for label
        ax3.plot([],c='k',label='UHSAS100')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CN (>100nm) # (cm$^{-3}$) '+campaign,fontsize=17)
        fig.text(0.1,0.95,'Above Clouds',fontsize=15)
        fig.text(0.1,0.6,'Near Clouds',fontsize=15)
        fig.text(0.1,0.25,'Below Cloud',fontsize=15)
        
        ax2.legend(loc='upper right', shadow=False, fontsize='x-large')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
    elif plot_method == 'height':
        #%% for CPC (>10nm)
        figname = figpath_aircraft_statistics+'percentile_lat_CN10nm_byheight_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(cpc_above_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ncn10_above_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        if campaign=='SOCRATES':
            ax1.set_ylim(-100,4000)
        elif campaign=='CSET':
            ax1.set_ylim(-20,1200)
        ax1.plot([],c='k',label='CPC')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(cpc_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ncn10_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        if campaign=='SOCRATES':
            ax3.set_ylim(-100,2000)
        elif campaign=='CSET':
            ax3.set_ylim(-50,4000)
        # plot temporal lines for label
        ax3.plot([],c='k',label='CPC')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CN (>10nm) # (cm$^{-3}$) '+campaign,fontsize=17)
        fig.text(0.1,0.9,'2-5km',fontsize=15)
        fig.text(0.1,0.4,'0-2km',fontsize=15)
        
        ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        #%% plot for UHSAS (>100nm)
        figname = figpath_aircraft_statistics+'percentile_lat_CN100nm_byheight_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(uhsas_above_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ncn100_above_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax1.plot([],c='k',label='UHSAS100')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(uhsas_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ncn100_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax3.set_ylim(-10,400)
        # plot temporal lines for label
        ax3.plot([],c='k',label='UHSAS100')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CN (>100nm) # (cm$^{-3}$) '+campaign,fontsize=17)
        fig.text(0.1,0.9,'2-5km',fontsize=15)
        fig.text(0.1,0.4,'0-2km',fontsize=15)
        
        ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    #%%    
    
    elif plot_method == 'all':
        #%% for CPC (>10nm)
        figname = figpath_aircraft_statistics+'percentile_lat_CN_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(cpc_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ncn10_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        if campaign=='SOCRATES':
            ax1.set_ylim(-100,4000)
        elif campaign=='CSET':
            ax1.set_ylim(-20,2500)
        ax1.plot([],c='k',label='CPC')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(uhsas_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ncn100_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        if campaign=='SOCRATES':
            ax3.set_ylim(-10,400)
        elif campaign=='CSET':
            ax3.set_ylim(-10,1000)
        # plot temporal lines for label
        ax3.plot([],c='k',label='UHSAS100')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CN # (cm$^{-3}$) '+campaign,fontsize=17)
        
        ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
        ax3.legend(loc='upper right', shadow=False, fontsize='x-large')
        fig.text(0.1,0.9,'>10nm',fontsize=15)
        fig.text(0.1,0.4,'>100nm',fontsize=15)
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    #%%
    else:
        raise ValueError('does not recognize plot_method: '+plot_method)