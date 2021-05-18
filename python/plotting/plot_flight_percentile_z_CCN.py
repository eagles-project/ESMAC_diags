# plot percentile of CCN number concentration with height
# for flight data in IOPs
# compare models and CCN measurements


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from read_aircraft import read_ccn_hiscale, read_ccn_socrates
from read_ARMdata import read_ccn
from read_netcdf import read_merged_size,read_extractflight

#%% settings

from settings import campaign, Model_List, color_model, \
    height_bin, E3SM_aircraft_path, figpath_aircraft_statistics

if campaign=='HISCALE' or campaign=='ACEENA':
    from settings import IOP, ccnpath, merged_size_path
elif campaign=='CSET' or campaign=='SOCRATES':
    from settings import ccnpath
else:
    print('ERROR: campaign name is not recognized: '+campaign)
    error
    
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

lst = glob.glob(E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[0]+'*.nc')
lst.sort()
if len(lst)==0:
    print('ERROR: cannot find any file at '+E3SM_aircraft_path)
    error
# choose files for specific IOP
if campaign=='HISCALE':
    if IOP=='IOP1':
        lst=lst[0:17]
    elif IOP=='IOP2':
        lst=lst[17:]
    elif IOP[0:4]=='2016':
        a=lst[0].split('_'+Model_List[0]+'_')
        lst = glob.glob(a[0]+'_'+Model_List[0]+'_'+IOP+'*')
        lst.sort()
elif campaign=='ACEENA':
    if IOP=='IOP1':
        lst=lst[0:20]
    elif IOP=='IOP2':
        lst=lst[20:]
    elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
        a=lst[0].split('_'+Model_List[0]+'_')
        lst = glob.glob(a[0]+'_'+Model_List[0]+'_'+IOP+'*')
        lst.sort()
        
alldates = [x.split('_')[-1].split('.')[0] for x in lst]
    
#%% read all data

height_all = []
ccna_all = []
ccnb_all = []
SSa_all = np.array([])
SSb_all = np.array([])
ccn3_all = []
ccn5_all = []
nmodels=len(Model_List)
for mm in range(nmodels):
    ccn3_all.append([])
    ccn5_all.append([])
    
print('reading '+format(len(alldates))+' files to calculate the statistics: ')

for date in alldates:
    
    #%% read in Models
    
    for mm in range(nmodels):
        filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
    
        (timem,heightm,ccn3,timeunitm,ccn3_unit,ccn3_longname)=read_extractflight(filename_m,'CCN3')
        (timem,heightm,ccn5,timeunitm,ccn5_unit,ccn5_longname)=read_extractflight(filename_m,'CCN5')
        
        
    # get supersaturation
    SS3 = ccn3_longname.split('=')[-1]
    SS5 = ccn5_longname.split('=')[-1]
    
    #%% read in flight data (for HISCALE)
    if campaign=='HISCALE':
        filename_ccn=glob.glob(ccnpath+'CCN_G1_'+date[0:8]+'*R2_HiScale001s.*')
        filename_ccn.sort()
        if date[-1]=='a':
            flightidx=1
        else:
            flightidx=2
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
            SSa[SSa<0]=np.nan
            SSb[SSb<0]=np.nan
        elif len(filename_ccn)==0:
            time_ccn=timem
            ccna=np.nan*np.empty([len(timem)])
            ccnb=np.nan*np.empty([len(timem)])
            SSa=0.24*np.full(len(timem),1)
            SSb=0.46*np.full(len(timem),1)
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
            SSa[SSa<0]=np.nan
        elif len(filename_ccna)==0:
            # print('no CCN data found. set as NaN')
            timea=timem
            SSa=np.nan*np.empty([len(timem)])
            ccna=np.nan*np.empty([len(timem)])
        else:
            print('find too many files, check: ')
            print(filename_ccna)
            error
        if len(filename_ccnb)==1:
            (timeb,timeunitb,ccnb,ccnunit,SSb)=read_ccn(filename_ccnb[0])
            ccnb[ccnb<0]=np.nan
            SSb[SSb<0]=np.nan
        elif len(filename_ccnb)==0:
            # print('no CCN data found. set as NaN')
            timeb=timem
            SSb=np.nan*np.empty([len(timem)])
            ccnb=np.nan*np.empty([len(timem)])
        else:
            print('find too many files, check: ')
            print(filename_ccnb)
            error
         
    # CSET does not have observed CCN
    elif campaign=='CSET':
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
            ccn[ccn<-9000]=np.nan
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
            print('find too many files, check: ')
            print(filename_ccn)
            error
            
    if any(timea!=timeb):
        error   
        
    # exclude NaNs
    idx = np.logical_or(~np.isnan(ccna), ~np.isnan(ccnb))
    ccna_all.append(ccna[idx])
    ccnb_all.append(ccnb[idx])
    SSa_all=np.append(SSa_all,SSa[idx])
    SSb_all=np.append(SSb_all,SSb[idx])
    
    height2=np.interp(timea,timem,heightm)
    height_all.append(height2[idx])
    
    # for interpolation of model results
    timea=timea[idx]
    timeb=timeb[idx]
    
    # interpolate model results into observational time
    timem2 = (timem-int(timem[0]))*86400
    for mm in range(nmodels):
        ccn3_all[mm].append(np.interp(timea,timem2,ccn3)) 
        ccn5_all[mm].append(np.interp(timeb,timem2,ccn5)) 
         
#%% calculate percentiles for each height bin

ccna_z = list()
ccnb_z = list()
ccn3_z = []
ccn5_z = []
nmodels=len(Model_List)
for mm in range(nmodels):
    ccn3_z.append([])
    ccn5_z.append([])
for zz in range(zlen):
    ccna_z.append(np.empty(0))
    ccnb_z.append(np.empty(0))
    for mm in range(nmodels):
        ccn3_z[mm].append(np.empty(0))
        ccn5_z[mm].append(np.empty(0))
    
ndays=len(height_all)
for dd in range(ndays):
    height = height_all[dd]
    ccna = ccna_all[dd]
    ccnb = ccnb_all[dd]
    for zz in range(zlen):
        idx = np.logical_and(height>=zmin[zz], height<zmax[zz])
        ccna_z[zz]=np.append(ccna_z[zz],ccna[np.logical_and(idx,~np.isnan(ccna))])
        ccnb_z[zz]=np.append(ccnb_z[zz],ccnb[np.logical_and(idx,~np.isnan(ccnb))])
        for mm in range(nmodels):
            ccn3 = ccn3_all[mm][dd]
            ccn5 = ccn5_all[mm][dd]
            ccn3_z[mm][zz]=np.append(ccn3_z[mm][zz],ccn3[idx])
            ccn5_z[mm][zz]=np.append(ccn5_z[mm][zz],ccn5[idx])
        

#%% make plot
# set position shift so that models and obs are not overlapped
p_shift = np.arange(nmodels+1)
p_shift = (p_shift - p_shift.mean())*0.2

if campaign=='HISCALE' or campaign=='ACEENA':
    figname = figpath_aircraft_statistics+'percentile_height_CCN_'+campaign+'_'+IOP+'.png'
else:
    figname = figpath_aircraft_statistics+'percentile_height_CCN_'+campaign+'.png'
print('plotting figures to '+figname)

fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,8))   # figsize in inches
# plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0
    
ax1.boxplot(ccna_z,whis=(5,95),showmeans=False,showfliers=False,
            positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
            boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
            medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
            vert=False, patch_artist=True)    # need patch_artist to fill color in box
for mm in range(nmodels):
    c = color_model[mm]
    ax1.boxplot(ccn3_z[mm],whis=(5,95),showmeans=False,showfliers=False,
            positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
            boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
            medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
            vert=False, patch_artist=True)    # need patch_artist to fill color in box
ax1.tick_params(color='k',labelsize=12)
ax1.set_xscale('log')
ax1.set_ylim(-1,zlen)
ax1.set_yticks(range(zlen))
ax1.set_yticklabels(z)
# ax1.set_yticks([1,3,5,7,9,11,12,13,14,15,16])
# ax1.set_yticklabels(range(400,4100,400))
# plot temporal lines for label
ax1.plot([],c='k',label='Obs ('+format(np.nanmean(SSa_all),'.2f')+'%)')
for mm in range(nmodels):
    ax1.plot([],c=color_model[mm],label=Model_List[mm])
ax1.legend(loc='upper right', fontsize='large')
    
ax2.boxplot(ccnb_z,whis=(5,95),showmeans=False,showfliers=False,
            positions=np.array(range(zlen))+p_shift[-1],widths=0.15,
            boxprops=dict(facecolor='k', color='k'),whiskerprops=dict(color='k'),
            medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color='k'),
            vert=False, patch_artist=True)    # need patch_artist to fill color in box
for mm in range(nmodels):
    c = color_model[mm]
    ax2.boxplot(ccn5_z[mm],whis=(5,95),showmeans=False,showfliers=False,
            positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
            boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
            medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
            vert=False, patch_artist=True)    # need patch_artist to fill color in box
ax2.tick_params(color='k',labelsize=12)
ax2.set_xscale('log')
ax2.set_ylim(-1,zlen)
ax2.set_yticks(range(zlen))
ax2.set_yticklabels([])
# ax1.set_yticks(np.arange(0,20,2))
# ax1.set_yticklabels(range(400,4100,400))
# plot temporal lines for label
ax2.plot([],c='k',label='Obs ('+format(np.nanmean(SSb_all),'.2f')+'%)')
for mm in range(nmodels):
    ax2.plot([],c=color_model[mm],label=Model_List[mm])
ax2.legend(loc='upper right', fontsize='large')
    
# set xlimit consistent in subplots
xlim1 = ax1.get_xlim()
xlim2 = ax2.get_xlim()
ax1.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])
ax2.set_xlim([min(xlim1[0],xlim2[0]), max(xlim1[1],xlim2[1])])

ax1.set_ylabel('Height (m MSL)',fontsize=14)
fig.text(0.4,0.06, 'CCN number (cm$^{-3}$)', fontsize=14)
ax1.set_title('SS = '+SS3,fontsize=14)
ax2.set_title('SS = '+SS5,fontsize=14)
if campaign=='HISCALE' or campaign=='ACEENA':
    fig.text(0.48,0.92, IOP, fontsize=16)

fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
plt.close()
    