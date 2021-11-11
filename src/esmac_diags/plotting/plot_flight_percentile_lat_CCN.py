"""
# plot percentile of aerosol number concentration binned by different latitudes
# separated by below-cloud, near-cloud and above-cloud
# for aircraft measurements in CSET or SOCRATES
"""


import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_ccn_socrates, read_RF_NCAR
from ..subroutines.read_netcdf import read_extractflight
from ..subroutines.quality_control import qc_remove_neg

def run_plot(settings):

    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    color_model = settings['color_model']
    latbin = settings['latbin']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']
    figpath_aircraft_statistics = settings['figpath_aircraft_statistics']
    
    if campaign in ['CSET', 'SOCRATES']:
        ccnpath = settings['ccnpath']
        RFpath = settings['RFpath']
    else:
        raise ValueError('This code is only for CSET or SOCRATES. check campaign setting: '+campaign)
    
    #%% other settings

    plot_method = 'all'  # 'height': separate by height. 'all': all heights below 5km
    
    if not os.path.exists(figpath_aircraft_statistics):
        os.makedirs(figpath_aircraft_statistics)
       
    dlat = latbin[1]-latbin[0]
    latmin = latbin-dlat/2
    latmax = latbin+dlat/2
    latlen = len(latbin)
        
    nmodels=len(Model_List)
    
    #%% find files for flight information
    
    lst = glob.glob(E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[0]+'_*.nc')
    lst.sort()
    if len(lst)==0:
        raise ValueError('cannot find any file')
    alldates = [x.split('_')[-1].split('.')[0] for x in lst]
    
    #%% define variables by latitude bins below, near and above clouds
        
    ccna_below_lat = []
    ccna_near_lat = []
    ccna_above_lat = []
    ccnb_below_lat = []
    ccnb_near_lat = []
    ccnb_above_lat = []
    for bb in range(latlen):
        ccna_below_lat.append(np.empty(0))
        ccna_near_lat.append(np.empty(0))
        ccna_above_lat.append(np.empty(0))
        ccnb_below_lat.append(np.empty(0))
        ccnb_near_lat.append(np.empty(0))
        ccnb_above_lat.append(np.empty(0))
        
    ccn3_below_lat = []
    ccn3_near_lat = []
    ccn3_above_lat = []
    ccn5_below_lat = []
    ccn5_near_lat = []
    ccn5_above_lat = []
    for mm in range(nmodels):
        ccn3_below_lat.append([np.empty(0) for bb in range(latlen)])
        ccn3_near_lat.append([np.empty(0) for bb in range(latlen)])
        ccn3_above_lat.append([np.empty(0) for bb in range(latlen)])
        ccn5_below_lat.append([np.empty(0) for bb in range(latlen)])
        ccn5_near_lat.append([np.empty(0) for bb in range(latlen)])
        ccn5_above_lat.append([np.empty(0) for bb in range(latlen)])
    
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
        
        #%% read in Models
        
        ccn3=[]
        ccn5=[]
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,ccn3_tmp,timeunitm,ccn3_unit,ccn3_longname)=read_extractflight(filename_m,'CCN3')
            (timem,heightm,ccn5_tmp,timeunitm,ccn5_unit,ccn5_longname)=read_extractflight(filename_m,'CCN5')
            ccn3.append(ccn3_tmp)
            ccn5.append(ccn5_tmp)
            
        # get supersaturation
        SS3 = ccn3_longname.split('=')[-1]
        SS5 = ccn5_longname.split('=')[-1]
           
        #%% read in observations for CSET and SOCRATES
        # CSET does not have observed CCN
        if campaign=='CSET':
            timea=timem
            SSa=np.nan*np.empty([len(timem)])
            ccna=np.nan*np.empty([len(timem)])
            timeb=timem
            SSb=np.nan*np.empty([len(timem)])
            ccnb=np.nan*np.empty([len(timem)])
            
        # SOCRATES
        elif campaign=='SOCRATES':
            filename_ccn=glob.glob(ccnpath+'CCNscanning_SOCRATES_GV_RF*'+date[0:8]+'_R0.ict')
            if len(filename_ccn)==1:
                (data0,ccnlist)=read_ccn_socrates(filename_ccn[0])
                time_ccn = data0[0,:]
                ccn = data0[1,:]
                SS = data0[3,:]
                ccn=qc_remove_neg(ccn)
                timea=time_ccn
                timeb=time_ccn
                ccna=np.array(ccn)
                ccnb=np.array(ccn)
                idxa=np.logical_and(SS>0.05, SS<0.15)
                ccna[idxa==False]=np.nan
                SSa=np.full((len(timea)),0.1)
                idxb=np.logical_and(SS>0.45, SS<0.55)
                ccnb[idxb==False]=np.nan
                SSb=np.full((len(timeb)),0.5)
            elif len(filename_ccn)==0:
                timea=timem
                SSa=np.nan*np.empty([len(timem)])
                ccna=np.nan*np.empty([len(timem)])
                timeb=timem
                SSb=np.nan*np.empty([len(timem)])
                ccnb=np.nan*np.empty([len(timem)])
            else:
                raise ValueError('find too many files: ' + filename_ccn)
                
        if any(timea!=timeb):
            raise ValueError('inconsitent time dimension')
            
        
        # need latitude from RF file
        lst = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
        if len(lst)==1 or len(lst)==2:  # SOCRATES has two flights in 20180217, choose the later one
            filename=lst[-1]
        else:
            raise ValueError('find no file or too many files: ' + lst)
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LAT')
        
        # exclude NaNs
        idx = np.logical_or(~np.isnan(ccna), ~np.isnan(ccnb))
        ccna=ccna[idx]
        ccnb=ccnb[idx]
        SSa=SSa[idx]
        SSb=SSb[idx]
        
        # for interpolation of model results
        timea=timea[idx]
        timeb=timeb[idx]
        time=timea
        # interpolate model results into observational time
        for mm in range(nmodels):
            ccn3[mm] = (np.interp(timea,timem,ccn3[mm])) 
            ccn5[mm] = (np.interp(timeb,timem,ccn5[mm])) 
        height = np.interp(timeb,timem,heightm)
        lat = np.interp(timeb,timem,lat)
            
            
        #%% separate data by cloud or height
        flag_below = np.zeros(len(time))
        flag_near = np.zeros(len(time))
        flag_above = np.zeros(len(time))
        
        if plot_method == 'height':
            for ii in range(len(time)):
                if height[ii]>5000:
                    continue   # exclude measurements above 5km
                elif height[ii]<2000:
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
                ccna_below_lat[bb] = np.hstack((ccna_below_lat[bb], ccna[idx][flag_below[idx]==1]))
                ccnb_below_lat[bb] = np.hstack((ccnb_below_lat[bb], ccnb[idx][flag_below[idx]==1]))
                for mm in range(nmodels):
                    ccn3_below_lat[mm][bb] = np.hstack((ccn3_below_lat[mm][bb], ccn3[mm][idx][flag_below[idx]==1]))
                    ccn5_below_lat[mm][bb] = np.hstack((ccn5_below_lat[mm][bb], ccn5[mm][idx][flag_below[idx]==1]))
            if any(flag_near[idx]==1):
                ccna_near_lat[bb] = np.hstack((ccna_near_lat[bb], ccna[idx][flag_near[idx]==1]))
                ccnb_near_lat[bb] = np.hstack((ccnb_near_lat[bb], ccnb[idx][flag_near[idx]==1]))
                for mm in range(nmodels):
                    ccn3_near_lat[mm][bb] = np.hstack((ccn3_near_lat[mm][bb], ccn3[mm][idx][flag_near[idx]==1]))
                    ccn5_near_lat[mm][bb] = np.hstack((ccn5_near_lat[mm][bb], ccn5[mm][idx][flag_near[idx]==1]))
            if any(flag_above[idx]==1):
                ccna_above_lat[bb] = np.hstack((ccna_above_lat[bb], ccna[idx][flag_above[idx]==1]))
                ccnb_above_lat[bb] = np.hstack((ccnb_above_lat[bb], ccnb[idx][flag_above[idx]==1]))
                for mm in range(nmodels):
                    ccn3_above_lat[mm][bb] = np.hstack((ccn3_above_lat[mm][bb], ccn3[mm][idx][flag_above[idx]==1]))
                    ccn5_above_lat[mm][bb] = np.hstack((ccn5_above_lat[mm][bb], ccn5[mm][idx][flag_above[idx]==1]))
              
    #%% remove nan elements in the observations
    for bb in range(latlen):
        ccna=ccna_below_lat[bb]
        ccna_below_lat[bb] = ccna[~np.isnan(ccna)]
        ccnb=ccnb_below_lat[bb]
        ccnb_below_lat[bb] = ccnb[~np.isnan(ccnb)]
        ccna=ccna_near_lat[bb]
        ccna_near_lat[bb] = ccna[~np.isnan(ccna)]
        ccnb=ccnb_near_lat[bb]
        ccnb_near_lat[bb] = ccnb[~np.isnan(ccnb)]
        ccna=ccna_above_lat[bb]
        ccna_above_lat[bb] = ccna[~np.isnan(ccna)]
        ccnb=ccnb_above_lat[bb]
        ccnb_above_lat[bb] = ccnb[~np.isnan(ccnb)]
    
    
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
    #%% plot separate by height
    if plot_method == 'height':
        # for ccna 
        figname = figpath_aircraft_statistics+'percentile_lat_CCN3_byheight_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(ccna_above_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ccn3_above_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax1.set_ylim(-10,50)
        ax1.plot([],c='k',label='OBS')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(ccna_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ccn3_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax3.set_ylim(-10,100)
        # plot temporal lines for label
        ax3.plot([],c='k',label='OBS')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CCN (SS='+format(np.nanmean(SSa),'.2f')+'%) # (cm$^{-3}$) '+campaign,fontsize=17)
        fig.text(0.1,0.9,'2-5km',fontsize=15)
        fig.text(0.1,0.4,'0-2km',fontsize=15)
        
        ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
        
        # plot for ccnb (SS=0.5%)
        figname = figpath_aircraft_statistics+'percentile_lat_CCN5_byheight_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(ccnb_above_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ccn5_above_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        ax1.plot([],c='k',label='OBS')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(ccnb_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ccn5_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        # ax3.set_ylim(-10,400)
        # plot temporal lines for label
        ax3.plot([],c='k',label='OBS')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CCN (SS='+format(np.nanmean(SSb),'.2f')+'%) # (cm$^{-3}$) '+campaign,fontsize=17)
        fig.text(0.1,0.9,'2-5km',fontsize=15)
        fig.text(0.1,0.4,'0-2km',fontsize=15)
        
        ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    #%%    
    
    elif plot_method == 'all':
        #%% for OBS
        figname = figpath_aircraft_statistics+'percentile_lat_CCN_'+campaign+'.png'
        print('plotting figures to '+figname)
        
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
            
        ax1.boxplot(ccna_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax1.boxplot(ccn3_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax1.tick_params(color='k',labelsize=15)
        # ax1.set_yscale('log')
        ax1.set_xlim(-1,latlen)
        ax1.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        if campaign=='SOCRATES':
            ax1.set_ylim(-10,200)
        elif campaign=='CSET':
            ax1.set_ylim(-10,200)
        ax1.plot([],c='k',label='OBS')
        for mm in range(nmodels):
            ax1.plot([],c=color_model[mm],label=Model_List[mm])
            
        ax3.boxplot(ccnb_below_lat,whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                    boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        for mm in range(nmodels):
            c = color_model[mm]
            ax3.boxplot(ccn5_below_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                    positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                    boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                    medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                    vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax3.tick_params(color='k',labelsize=15)
        # ax3.set_yscale('log')
        ax3.set_xlim(-1,latlen)
        ax3.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
        if campaign=='SOCRATES':
            ax3.set_ylim(-10,500)
        elif campaign=='CSET':
            ax3.set_ylim(-10,500)
        # plot temporal lines for label
        ax3.plot([],c='k',label='OBS')
        for mm in range(nmodels):
            ax3.plot([],c=color_model[mm],label=Model_List[mm])
        
        ax3.set_xlabel('Latitude',fontsize=16)
        
        ax1.set_xticklabels([])
        ax3.set_xticklabels([])
        ax3.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
        ax1.set_title('Percentile of CCN # (cm$^{-3}$) '+campaign,fontsize=17)
        
        ax1.legend(loc='upper right', shadow=False, fontsize='x-large')
        ax3.legend(loc='upper right', shadow=False, fontsize='x-large')
        fig.text(0.1,0.9,'SS='+format(np.nanmean(SSa),'.2f')+'%',fontsize=15)
        fig.text(0.1,0.4,'SS='+format(np.nanmean(SSb),'.2f')+'%',fontsize=15)
        
        fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    #%%
    else:
        raise ValueError('does not recognize plot_method: '+plot_method)
