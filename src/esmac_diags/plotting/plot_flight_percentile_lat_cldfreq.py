"""# plot percentile of meteorological variables binned by different latitudes
# for aircraft measurements in CSET or SOCRATES
# only select a certain height ranges for warm clouds (the height range needs to be further tuned)
"""

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from ..subroutines.read_aircraft import read_RF_NCAR
from ..subroutines.read_netcdf import read_extractflight
from ..subroutines.specific_data_treatment import lwc2cflag
from ..subroutines.quality_control import qc_mask_takeoff_landing

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
        raise ValueError('cannot fine any file')
    alldates = [x.split('_')[-1].split('.')[0] for x in lst]
    
    
    #%% define variables by latitude bins
        
    height_lat = []
    cbheight = []         # cloud base height
    cflag_lat = []
    cloudo_lat = []        # cloud fraction by flag
    
    for bb in range(latlen):
        height_lat.append(np.empty(0))
        cbheight.append(np.empty(0))
        cflag_lat.append(np.empty(0))
        cloudo_lat.append(np.empty(0))
        
    cloudm_lat = []
    for mm in range(nmodels):
        cloudm_lat.append(list(cloudo_lat))
    
    print('reading '+format(len(alldates))+' files to calculate the statistics: ')
    
    for date in alldates:
        print(date)
        
        #%% read in Models
        cloudm = []
        for mm in range(nmodels):
            filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
            (timem,heightm,cloud,timeunitm,clunit,cllongname)=read_extractflight(filename_m,'CLOUD')
            cloudm.append(cloud)   
        
        #%% read in observations
        # note that it is only for CSET and SOCRATES
        lst = glob.glob(RFpath+'RF*'+date+'*.PNI.nc')
        if len(lst)==1 or len(lst)==2:  # SOCRATES has two flights in 20180217, choose the later one
            filename=lst[-1]
        else:
            raise ValueError('find no file or too many files: '+lst)
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'ALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LAT')
        (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LON')
        (time,lwc,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename,'PLWCC')
        
        # exclude 30min after takeoff and before landing
        height=qc_mask_takeoff_landing(time,height)
        lat=qc_mask_takeoff_landing(time,lat)
        lon=qc_mask_takeoff_landing(time,lon)
        lwc=qc_mask_takeoff_landing(time,lwc)
        timem=qc_mask_takeoff_landing(time,timem)
        for mm in range(nmodels):
            cloudm[mm]=qc_mask_takeoff_landing(time,cloudm[mm])
        
        # calculate cloud flag based on LWC
        cldflag=lwc2cflag(lwc,lwcunit)
        
            
        #%% put data in each latitude bin
        for bb in range(latlen):
            idx = np.logical_and(lat>=latmin[bb], lat<latmax[bb])
            height2 = height[idx]
            cldflag2 = cldflag[idx]
            
            # set specific height range
            idx2 = height2<5000
            
            if any(cldflag2==1):
                cbheight[bb] = np.hstack((cbheight[bb],min(height2[cldflag2==1])))
            
            if len(idx2)!=0:
                height_lat[bb] = np.hstack((height_lat[bb],height2[idx2]))
                cflag_lat[bb] = np.hstack((cflag_lat[bb],cldflag2[idx2]))
                cloudo_lat[bb] = np.hstack((cloudo_lat[bb],sum(cldflag2)/len(cldflag2)))
                for mm in range(nmodels):
                    cld_temp = np.nanmean(cloudm[mm][idx][idx2])
                    if ~np.isnan(cld_temp):
                        cloudm_lat[mm][bb] = np.hstack((cloudm_lat[mm][bb], cld_temp))
        
    #%% make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(nmodels+1)
    p_shift = (p_shift - p_shift.mean())*0.2
    
        
    figname = figpath_aircraft_statistics+'percentile_lat_CldFreq_'+campaign+'.png'
    print('plotting figures to '+figname)
    
    fig,ax = plt.subplots(1,1,figsize=(8,2))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
        
        
    ax.boxplot(cloudo_lat,whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[-1],widths=0.15,
                boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    for mm in range(nmodels):
        c = color_model[mm]
        ax.boxplot(cloudm_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
    ax.tick_params(color='k',labelsize=15)
    # ax.set_yscale('log')
    ax.set_xlim(-1,latlen)
    ax.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    # plot temporal lines for label
    ax.plot([],c='k',label='OBS')
    for mm in range(nmodels):
        ax.plot([],c=color_model[mm],label=Model_List[mm])
    
    ax.set_xlabel('Latitude',fontsize=16)
    
    ax.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
    ax.set_title('Cloud Fraction (Fraction)',fontsize=17)
    # ax4.set_title(varmlongname[3]+' ('+varmunit[3]+')',fontsize=15)
    
    ax.legend(loc='upper right', shadow=False, fontsize='x-large')
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)