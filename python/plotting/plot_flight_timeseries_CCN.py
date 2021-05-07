# plot aircraft track data
# timeseries of CCN number concentration 
# compare models and aircraft measurements


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from read_aircraft import read_ccn_hiscale
from read_ARMdata import read_ccn
from read_netcdf import read_merged_size,read_extractflight

#%% settings

from settings import campaign, ccnpath,merged_size_path, Model_List, color_model, \
    IOP, E3SM_aircraft_path, figpath_aircraft_timeseries

import os
if not os.path.exists(figpath_aircraft_timeseries):
    os.makedirs(figpath_aircraft_timeseries)
    

#%% find files for flight information
lst = glob.glob(merged_size_path+'merged_bin_*'+campaign+'*.nc')
lst.sort()


if len(lst)==0:
    print('ERROR: cannot find any file at '+merged_size_path)
    error
    
# choose files for specific IOP
if campaign=='HISCALE':
    if IOP=='IOP1':
        lst=lst[0:17]
    elif IOP=='IOP2':
        lst=lst[17:]
    elif IOP[0:4]=='2016':
        a=lst[0].split('_'+campaign+'_')
        lst = glob.glob(a[0]+'*'+IOP+'*')
        lst.sort()
elif campaign=='ACEENA':
    if IOP=='IOP1':
        lst=lst[0:20]
    elif IOP=='IOP2':
        lst=lst[20:]
    elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
        a=lst[0].split('_'+campaign+'_')
        lst = glob.glob(a[0]+'*'+IOP+'*')
        lst.sort()
else:
    print('ERROR: campaign name is not recognized: '+campaign)
    error
    
# for each flight
for filename in lst:
    
    # get date info:        
    date=filename[-12:-3]
    if date[-1]=='a':
        flightidx=1
    else:
        flightidx=2

    #% read in flight information
    (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
    (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
    (time,size,height,timeunit,zunit,long_name)=read_merged_size(filename,'height')
    time=np.ma.compressed(time)
    
    #%% read in ccn
    if campaign=='HISCALE':
        filename_ccn=glob.glob(ccnpath+'CCN_G1_'+date[0:8]+'*R2_HiScale001s.*')
        filename_ccn.sort()
        # read in data
        if len(filename_ccn)==1 or len(filename_ccn)==2:
            (data0,ccnlist)=read_ccn_hiscale(filename_ccn[flightidx-1])
            # only choose data quality is good (flag=0)
            flag = data0[7,:]
            time_ccn = data0[0,:]
            ccna = data0[10,:]
            ccnb = data0[11,:]
            SSa = data0[2,:]
            SSb = data0[5,:]
            idx = flag==0.0
            ccna[flag!=0]=np.nan
            ccnb[flag!=0]=np.nan
        elif len(filename_ccn)==0:
            time_ccn=time
            ccna=np.nan*np.empty([len(time)])
            ccnb=np.nan*np.empty([len(time)])
            SSa=0.24
            SSb=0.46
        else:
            print('find too many files, check: ')
            print(filename_ccn)
            error
        timea=time_ccn
        timeb=time_ccn
        
    elif campaign=='ACEENA':
        filename_ccna=glob.glob(ccnpath+'enaaafccn2colaF1.b1.'+date[0:8]+'*.nc')
        filename_ccnb=glob.glob(ccnpath+'enaaafccn2colbF1.b1.'+date[0:8]+'*.nc')
        # read in data
        if len(filename_ccna)==1:
            (timea,timeunita,ccna,ccnunit,SSa)=read_ccn(filename_ccna[0])
            ccna[ccna<0]=np.nan
        elif len(filename_ccna)==0:
            # print('no CCN data found. set as NaN')
            timea=time
            SSa=np.nan*np.empty([len(time)])
            ccna=np.nan*np.empty([len(time)])
        else:
            print('find too many files, check: ')
            print(filename_ccna)
            error
        if len(filename_ccnb)==1:
            (timeb,timeunitb,ccnb,ccnunit,SSb)=read_ccn(filename_ccnb[0])
            ccnb[ccnb<0]=np.nan
        elif len(filename_ccnb)==0:
            # print('no CCN data found. set as NaN')
            timeb=time
            SSb=np.nan*np.empty([len(time)])
            ccnb=np.nan*np.empty([len(time)])
        else:
            print('find too many files, check: ')
            print(filename_ccnb)
            error
        
    #%% read in Models
    nmodels=len(Model_List)
    ccn3_m = []
    ccn5_m = []
    for mm in range(nmodels):
        filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
    
        (timem,heightm,ccn3,timeunitm,ccn3_unit,ccn3_longname)=read_extractflight(filename_m,'CCN3')
        (timem,heightm,ccn5,timeunitm,ccn5_unit,ccn5_longname)=read_extractflight(filename_m,'CCN5')
        
        ccn3_m.append(ccn3) 
        ccn5_m.append(ccn5) 
        
    # get supersaturation
    SS3 = ccn3_longname.split('=')[-1]
    SS5 = ccn5_longname.split('=')[-1]
    timem2 = timem/3600
    
    #%% make plot
        
    figname = figpath_aircraft_timeseries+'timeseries_CCN_'+campaign+'_'+date+'.png'
    print('plotting figures to '+figname)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))   # figsize in inches
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
    ax1.plot(timea/3600,ccna,color='k',linewidth=1,label='OBS')
    for mm in range(nmodels):
        ax1.plot(timem2, ccn3_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax1.set_yscale('log')
    ax1.tick_params(color='k',labelsize=12)
    ylim1 = ax1.get_ylim()
    
    ax2.plot(timeb/3600,ccnb,color='k',linewidth=1,label='OBS')
    for mm in range(nmodels):
        ax2.plot(timem2, ccn5_m[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
    # ax2.set_yscale('log')
    ax2.tick_params(color='k',labelsize=12)
    ylim2 = ax2.get_ylim()
    
    # set ylimit consistent in subplots
    ax1.set_ylim([ylim1[0], ylim2[1]])
    ax2.set_ylim([ylim1[0], ylim2[1]])
    
    ax1.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    ax2.legend(loc='center right', shadow=False, fontsize='large',bbox_to_anchor=(1.25, .5))
    
    # supersaturation
    fig.text(0.67,0.9,'SS_obs='+format(np.nanmean(SSa),'.2f')+'%, SS_model='+SS3)
    fig.text(0.67,0.4,'SS_obs='+format(np.nanmean(SSb),'.2f')+'%, SS_model='+SS5)
    
    ax2.set_xlabel('time (hour UTC) '+date,fontsize=14)
    ax1.set_title('CCN Number Concentration (cm$^{-3}$)',fontsize=15)
    
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    plt.close()
    