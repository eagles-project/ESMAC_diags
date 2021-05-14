# plot surface timeseries of aerosol size distribution
# compare models and surface measurements


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from time_format_change import yyyymmdd2cday, hhmmss2sec,cday2mmdd
from read_surface import read_smpsb_pnnl,read_smps_bin
from read_ARMdata import read_uhsas, read_smps_bnl
from read_netcdf import read_E3SM

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

from settings import campaign, Model_List, IOP, start_date, end_date, E3SM_sfc_path, figpath_sfc_timeseries
if campaign=='ACEENA':
    from settings import uhsaspath
elif campaign=='HISCALE':
    if IOP=='IOP1':
        from settings import smps_bnl_path, nanosmps_bnl_path
    elif IOP=='IOP2':
        from settings import smps_pnnl_path

import os
if not os.path.exists(figpath_sfc_timeseries):
    os.makedirs(figpath_sfc_timeseries)
    
# change start date into calendar day
cday1 = yyyymmdd2cday(start_date,'noleap')
cday2 = yyyymmdd2cday(end_date,'noleap')
if start_date[0:4]!=end_date[0:4]:
    print('ERROR: currently not support multiple years. please set start_date and end_date in the same year')
    error
year0 = start_date[0:4]

# set time resolution for plotting. longer time needs coarser resolution to prevent memory error
dt_res = 3600  # in sec


#%% read in obs data
if campaign=='ACEENA':
    if IOP=='IOP1':
        lst = glob.glob(uhsaspath+'enaaosuhsasC1.a1.2017062*')+glob.glob(uhsaspath+'enaaosuhsasC1.a1.201707*')
    elif IOP=='IOP2':
        lst = glob.glob(uhsaspath+'enaaosuhsasC1.a1.201801*')+glob.glob(uhsaspath+'enaaosuhsasC1.a1.201802*')
    lst.sort()
    t_uhsas=np.empty(0)
    uhsas=np.empty((0,99))
    for filename in lst:
        (time,dmin,dmax,data,timeunit,dataunit,long_name) = read_uhsas(filename,'concentration')
        timestr=timeunit.split(' ')
        date=timestr[2]
        cday=yyyymmdd2cday(date)
        # average in time for quicker plot
        time2=np.arange(0,86400,dt_res)
        data2 = avg_time(time,data,time2)
        t_uhsas=np.hstack((t_uhsas, cday+time2/86400))
        uhsas=np.vstack((uhsas, data2))
    size_u = (dmin+dmax)/2
    uhsas[uhsas<0]=np.nan
    # change to dN/dlogDp
    dlnDp_u=np.empty(99)
    for bb in range(len(size_u)):
        dlnDp_u[bb]=np.log10(dmax[bb]/dmin[bb])
        uhsas[:,bb]=uhsas[:,bb]/dlnDp_u[bb]
    
    timeo = np.array(t_uhsas)
    size = np.array(size_u)
    obs = np.array(uhsas.T)
    
elif campaign=='HISCALE':    
    if IOP=='IOP1':
        lst = glob.glob(smps_bnl_path+'*.nc')
        lst.sort()
        t_smps=np.empty(0)
        smps=np.empty((0,192))
        for filename in lst:
            (time,size,flag,timeunit,dataunit,smps_longname)=read_smps_bnl(filename,'status_flag')
            (time,size,data,timeunit,smpsunit,smps_longname)=read_smps_bnl(filename,'number_size_distribution')
            d=data<0
            data[d]=np.nan
            data[flag!=0,:]=np.nan
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date)
            # average in time for quicker plot
            time2=np.arange(0,86400,dt_res)
            data2 = avg_time(time,data,time2)
            t_smps=np.hstack((t_smps, cday+time2/86400))
            smps=np.vstack((smps, data2))
        smps=smps.T
        # combine with nanoSMPS
        lst2 = glob.glob(nanosmps_bnl_path+'*.nc')
        lst2.sort()
        t_nano=np.empty(0)
        nanosmps=np.empty((0,192))
        for filename2 in lst2:
            (timen,sizen,flagn,timenunit,datanunit,long_name)=read_smps_bnl(filename2,'status_flag')
            (timen,sizen,datan,timenunit,nanounit,nanoname)=read_smps_bnl(filename2,'number_size_distribution')
            datan[datan<0]=np.nan
            datan[flagn!=0,:]=np.nan
            timestr=timenunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date)
            # average in time for quicker plot
            time2=np.arange(0,86400,dt_res)
            data2 = avg_time(timen,datan,time2)
            t_nano=np.hstack((t_nano, cday+time2/86400))
            nanosmps=np.vstack((nanosmps, data2))
        nanosmps=nanosmps.T
        for tt in range(smps.shape[1]):
            if any(t_nano==t_smps[tt]):
                smps[0:80,tt]=nanosmps[0:80,t_nano==t_smps[tt]].reshape(80)/3.8
        
    elif IOP=='IOP2':
        data=read_smpsb_pnnl(smps_pnnl_path+'HiScaleSMPSb_SGP_20160827_R1.ict')
        size=read_smps_bin(smps_pnnl_path+'NSD_column_size_chart.txt')
        time=data[0,:]
        smps=data[1:-1,:]
        flag=data[-1,:]
        smps[:,flag!=0]=np.nan
        cday=yyyymmdd2cday('2016-08-27')
        # average in time for quicker plot
        time2 = np.arange(time[0],time[-1],dt_res)
        smps = avg_time(time,smps.T,time2)
        smps = smps.T
        t_smps=cday+time2/86400
        
    timeo = np.array(t_smps)
    size = np.array(size)
    obs = np.array(smps)  
    
    # SMPS is already divided by log10
    
else:
    print('ERROR: does NOT recognize this campaign: '+campaign)
    error
    
#%% read in models
model = []
nmodels = len(Model_List)
for mm in range(nmodels):
    tmp_data=np.empty((3000,0))
    timem=np.empty(0)
    for cday in range(cday1,cday2+1):
        mmdd=cday2mmdd(cday)
        date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
        
        filename_input = E3SM_sfc_path+'SFC_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        (time,ncn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NCNall')
        
        timem = np.hstack((timem,time))
        tmp_data = np.hstack((tmp_data,ncn*1e-6))
    
    # change to dN/dlog10Dp
    for bb in range(3000):
        dlnDp=np.log10((bb+2)/(bb+1))
        tmp_data[bb,:]=tmp_data[bb,:]/dlnDp
    
    model.append(tmp_data)

sizem = np.arange(1,3001)

#%% make plot

figname = figpath_sfc_timeseries+'timeseries_AerosolSize_'+campaign+'_'+IOP+'.png'
print('plotting figures to '+figname)

#fig = plt.figure()
fig,ax = plt.subplots(nmodels+1,1,figsize=(8,2*(nmodels+1)))   # figsize in inches
plt.tight_layout(h_pad=1.1)   #pad=0.4, w_pad=0.5, h_pad=1.0
plt.subplots_adjust(right=0.9,bottom=0.1)

leveltick=[0.1,1,10,100,1000,10000,100000]
levellist=np.arange(np.log(leveltick[0]),12,.5)

obs[obs<0.01]=0.01
h1 = ax[0].contourf(timeo,size,np.log(obs),levellist,cmap=plt.get_cmap('jet'))

# d_mam=np.arange(1,3001)
h2=[]
for mm in range(nmodels):
    datam = model[mm]
    datam[datam<0.01]=0.01
    h_m = ax[mm+1].contourf(timem,sizem,np.log(datam),levellist,cmap=plt.get_cmap('jet'))
    h2.append(h_m)

# colorbar
cax = plt.axes([0.95, 0.2, 0.02, 0.6])
cbar=fig.colorbar(h2[0], cax=cax, ticks=np.log(leveltick))
cbar.ax.set_yticklabels(leveltick, fontsize=14)

# set axis
for ii in range(nmodels+1):
    ax[ii].set_xlim(cday1,cday2)
    ax[ii].set_yscale('log')
    ax[ii].set_ylim(1, 5000)
    ax[ii].set_yticks([1,10,100,1000])
    ax[ii].tick_params(color='k',labelsize=14)
    if ii==0:
        ax[ii].text(0.01, 0.94, 'OBS', fontsize=14,transform=ax[ii].transAxes, verticalalignment='top')
    else:
        ax[ii].text(0.01, 0.94, Model_List[ii-1], fontsize=14,transform=ax[ii].transAxes, verticalalignment='top')
    
ax[1].set_ylabel('Diameter (nm)',fontsize=14)
ax[0].set_title('Size Distribution (#/dlogDp, cm-3) '+campaign+' '+IOP,fontsize=15)
ax[nmodels].set_xlabel('Calendar Day',fontsize=14)

fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# plt.close()



