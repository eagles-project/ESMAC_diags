# plot vertical profile of cloud fraction for all flights in each IOP
# compare models and aircraft measurements


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from read_aircraft import read_wcm
from read_netcdf import read_extractflight,read_merged_size

#%% settings

from settings import campaign, merged_size_path, wcmpath, Model_List, color_model, IOP,  \
    height_bin, E3SM_aircraft_path, figpath_aircraft_statistics
    
import os
if not os.path.exists(figpath_aircraft_statistics):
    os.makedirs(figpath_aircraft_statistics)
    
    
#%%
z=height_bin
dz = z[1]-z[0]
zmin=z-np.insert((z[1:]-z[0:-1])/2,0,dz)
zmax=z+np.append((z[1:]-z[0:-1])/2,dz)

zlen=len(z)   

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

if len(lst)==0:
    print('ERROR: cannot find any file for '+IOP)
    error
    
#%% read all data

heightall=[]
lwc083all=[]
lwc021all=[]
lwcmall=[]

nmodels=len(Model_List)
for mm in range(nmodels):
    lwcmall.append([])
    
print('reading '+format(len(lst))+' files to calculate the statistics: ')

for filename in lst:
    
    # get date info:        
    date=filename[-12:-3]
    if date[-1]=='a':
        flightidx=1
    else:
        flightidx=2
    print(date)
    
    #% read in flight information
    (time,size,height,timeunit,cunit,long_name)=read_merged_size(filename,'height')
    time=np.ma.compressed(time)
    
    
     #%% read in WCM
    filename_wcm = glob.glob(wcmpath+'WCM_G1_'+date[0:8]+'*')
    if len(filename_wcm)==0:
        print('skip this date: '+date)
        continue
    (wcm,wcmlist)=read_wcm(filename_wcm[flightidx-1])
    time0=wcm[0,:]
    flag=wcm[-1,:]
    lwc083=wcm[2,:]
    lwc021=wcm[3,:]
    lwc083[lwc083<=0]=0.
    lwc021[lwc021<=0]=0.
    lwc083[flag!=0]=np.nan
    lwc021[flag!=0]=np.nan
    
    lwc083all.append(lwc083)
    lwc021all.append(lwc021)
    heightall.append(height)
    
    #%% read in models
    
    for mm in range(nmodels):
        filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        
        (timem,heightm,lwc,timeunit,cldunit,cldname)=read_extractflight(filename_m,'LWC')
            
         # change E3SM unit from kg/m3 to g/m3 
        lwcmall[mm].append(lwc*1000)
        
#%% calculate percentiles for each height bin

lwc083_z = list()
lwc021_z = list()
lwcm_z = []
for mm in range(nmodels):
    lwcm_z.append([])
for zz in range(zlen):
    lwc083_z.append(np.empty(0))
    lwc021_z.append(np.empty(0))
    for mm in range(nmodels):
        lwcm_z[mm].append(np.empty(0))
    
ndays=len(heightall)
# ndays=1;
for dd in range(ndays):
    height = heightall[dd]
    lwc083  = lwc083all[dd]
    lwc021  = lwc021all[dd]
    for zz in range(zlen):
        idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
        lwc083_z[zz]=np.append(lwc083_z[zz],lwc083[idx])
        lwc021_z[zz]=np.append(lwc021_z[zz],lwc021[idx])
        
    for mm in range(nmodels):
        lwcm = lwcmall[mm][dd]
        for zz in range(zlen):
            idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
            lwcm_z[mm][zz]=np.append(lwcm_z[mm][zz],lwcm[idx])
      
#%% remove all NANs and calculate cloud frequency
lwcmean_083 = np.full(zlen,np.nan)
lwcmean_021 = np.full(zlen,np.nan)
std_lwc_083 = np.full(zlen,np.nan)
lwcmean_m = []
for mm in range(nmodels):
    lwcmean_m.append(np.full(zlen,np.nan))
    
for zz in range(zlen):
    data = lwc083_z[zz]
    data = data[~np.isnan(data)]
    if len(data)>0:
        lwcmean_083[zz] = np.mean(data)
        std_lwc_083[zz] = np.std(data)/np.sqrt(len(data))
    data = lwc021_z[zz]
    data = data[~np.isnan(data)]
    if len(data)>0:
        lwcmean_021[zz] = np.mean(data)
    for mm in range(nmodels):
        data = lwcm_z[mm][zz]
        data = data[~np.isnan(data)]
        if len(data)>0:
            lwcmean_m[mm][zz] = np.mean(data)
            
#%% plot frequency  
figname = figpath_aircraft_statistics+'profile_height_LWC_'+campaign+'_'+IOP+'.png'
print('plotting figures to '+figname)

fig,ax = plt.subplots(figsize=(4,8))

ax.plot(lwcmean_083,z,color='k',linewidth=1,linestyle='-',label='Obs')
ax.fill_betweenx(z,lwcmean_083-std_lwc_083,lwcmean_083+std_lwc_083,facecolor='k',alpha=0.2)

for mm in range(nmodels):
    ax.plot(lwcmean_m[mm],z,color=color_model[mm],linewidth=1,label=Model_List[mm])

ax.tick_params(color='k',labelsize=12)
# ax.set_ylim(-1,zlen)
# ax.set_yticks(range(zlen))
# ax.set_yticks(z[0:-1:2])
ax.set_ylabel('Height (m MSL)',fontsize=12)
ax.legend(loc='upper right', fontsize='large')
ax.set_xlabel('LWC (g/m3)',fontsize=12)
ax.set_title(IOP,fontsize=15)

fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            