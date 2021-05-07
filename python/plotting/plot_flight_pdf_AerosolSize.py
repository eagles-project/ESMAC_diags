# plot mean aerosol size ditribution for aircraft track data
# average for each IOP
# compare models and aircraft measurements


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
# matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
# from time_format_change import yyyymmdd2cday, hhmmss2sec
from read_netcdf import read_merged_size,read_extractflight

# define function of averaging in time for faster plotting
def avg_time(time0,data0,time):
    data0[data0<0]=np.nan
    if data0.shape[0]!=len(time0):
        error
    data = np.full((len(time),data0.shape[1]),np.nan)
    dt=(time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0>=time[tt]-dt,time0<=time[tt]+dt)
        data[tt,:]=np.nanmean(data0[idx,:],axis=0)
    return(data)

#%% settings

from settings import campaign, merged_size_path, Model_List, color_model, IOP, \
    E3SM_aircraft_path, figpath_aircraft_statistics

import os
if not os.path.exists(figpath_aircraft_statistics):
    os.makedirs(figpath_aircraft_statistics)
    
# option of exclude measurements not in a flight leg
is_flight_leg = 'False'
    
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
    
    
#%% read all data

print('reading '+format(len(lst))+' files to calculate mean aerosol pdf: ')

nmodels=len(Model_List)
size_m = np.zeros(3000)
pdf_model = [size_m for mm in range(nmodels)]
(time,size,cvi,timeunit,cunit,long_name)=read_merged_size(lst[0],'CVI_inlet')
pdf_obs = np.zeros(len(size)) 

n_o = 0
n_m = [0 for mm in range(nmodels)]
    
for filename in lst:
    
    # get date info:        
    date=filename[-12:-3]
    if date[-1]=='a':
        flightidx=1
    else:
        flightidx=2
    print(date)

    #% read in observation
    (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
    (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
    (time,size,legnum,timeunit,zunit,long_name)=read_merged_size(filename,'leg_number')
    (time,size,sizeh,timeunit,dataunit,long_name)=read_merged_size(filename,'size_high')
    (time,size,sizel,timeunit,dataunit,long_name)=read_merged_size(filename,'size_low')
    (time,size,merge,timeunit,dataunit,long_name)=read_merged_size(filename,'size_distribution_merged')
    time=np.ma.compressed(time)
    size=size*1000.
    merge=np.array(merge.T)
    time=time/3600.

    # change to dN/dlnDp
    for bb in range(len(size)):
        dlnDp=np.log(sizeh[bb]/sizel[bb])
        merge[bb,:]=merge[bb,:]/dlnDp
    
    merge[merge<0]=np.nan
    merge[:,cflag==1]=np.nan
    if is_flight_leg=='True':
        merge[:,np.logical_or(legnum<=0, legnum>=100)]=np.nan
    
    idx_valid = ~np.isnan(np.mean(merge,0))
    pdf_obs = pdf_obs + np.sum(merge[:,idx_valid],1)
    n_o = n_o + np.sum(idx_valid)
    
    #%% read in Models
    for mm in range(nmodels):
        filename_m = E3SM_aircraft_path+'Aircraft_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        (timem,heightm,datam,timeunitm,datamunit,datamlongname)=read_extractflight(filename_m,'NCNall')
        datam=datam*1e-6    # #/m3 to #/cm3
        # change to dN/dlnDp
        for bb in range(3000):
            dlnDp=np.log((bb+2)/(bb+1))
            datam[bb,:]=datam[bb,:]/dlnDp
            
        if is_flight_leg=='True':
            datam[:,np.logical_or(legnum<=0, legnum>=100)]=np.nan
        
        for tt in range(len(timem)):
            if ~np.isnan(datam[0,tt]):
                pdf_model[mm] = pdf_model[mm]+datam[:,tt]
                n_m[mm]=n_m[mm]+1
                

#%% calculate mean pdf

pdf_obs=pdf_obs/n_o
for mm in range(nmodels):
    pdf_model[mm]=pdf_model[mm]/n_m[mm]

#%% make plot

if is_flight_leg=='True':
    figname = figpath_aircraft_statistics+'pdf_AerosolSize_'+campaign+'_'+IOP+'_legsonly.png'
else:
    figname = figpath_aircraft_statistics+'pdf_AerosolSize_'+campaign+'_'+IOP+'.png'

print('plotting figures to '+figname)

#fig = plt.figure()
fig,ax = plt.subplots(figsize=(6,3))   # figsize in inches

ax.plot(size,pdf_obs,color='k',label='Obs')
for mm in range(nmodels):
    ax.plot(np.arange(1,3001),pdf_model[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])

ax.legend(loc='upper right', shadow=False, fontsize='medium')
ax.tick_params(color='k',labelsize=12)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.01,100000)
ax.set_xlabel('Diameter (nm)',fontsize=12)
ax.set_title('Aerosol Size Distribution (#/dlnDp, cm$^{-3}$) '+IOP,fontsize=13)

fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# plt.close()