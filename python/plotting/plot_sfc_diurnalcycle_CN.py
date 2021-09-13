# plot surface diurnal cycle of aerosol size distribution
# compare models and surface measurements


import sys
sys.path.insert(1,'../subroutines/')

import matplotlib
matplotlib.use('AGG') # plot without needing X-display setting
import matplotlib.pyplot as plt
import numpy as np
import glob
from time_format_change import yyyymmdd2cday,cday2mmdd
from read_ARMdata import read_cpc
from read_netcdf import read_E3SM

# define function used in the code
def avg_time(time0,data0,time):
    data0[data0<0]=np.nan
    if data0.shape[0]!=len(time0):
        error
    data = np.full((len(time)),np.nan)
    dt=(time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0>=time[tt]-dt,time0<=time[tt]+dt)
        data[tt]=np.nanmean(data0[idx],axis=0)
    return(data)

#%% settings

from settings import campaign, cpcsfcpath, cpcusfcpath, Model_List, color_model, \
    IOP, start_date, end_date, E3SM_sfc_path, figpath_sfc_timeseries


# change start date into calendar day
cday1 = yyyymmdd2cday(start_date,'noleap')
cday2 = yyyymmdd2cday(end_date,'noleap')
if start_date[0:4]!=end_date[0:4]:
    print('ERROR: currently not support multiple years. please set start_date and end_date in the same year')
    error
year0 = start_date[0:4]
    
import os
if not os.path.exists(figpath_sfc_timeseries):
    os.makedirs(figpath_sfc_timeseries)
    
#%% read in obs data
if campaign=='ACEENA':
    # cpc
    if IOP=='IOP1':
        lst = glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.2017062*')+glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201707*')
    elif IOP=='IOP2':
        lst = glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201801*')+glob.glob(cpcsfcpath+'enaaoscpcfC1.b1.201802*')
    lst.sort()
    t_cpc=np.empty(0)
    cpc=np.empty(0)
    for filename in lst:
        (time,data,timeunit,cpcunit)=read_cpc(filename)
        timestr=timeunit.split(' ')
        date=timestr[2]
        cday=yyyymmdd2cday(date,'noleap')
        # average in time for better comparison with obs
        time2=np.arange(0,86400,3600)
        data2 = avg_time(np.array(time),np.array(data),time2)
        t_cpc=np.hstack((t_cpc, cday+time2/86400))
        cpc=np.hstack((cpc, data2))
    cpc[cpc<0]=np.nan
    # no cpcu
    t_cpcu = np.array([np.nan])
    cpcu = np.array([np.nan])
    
elif campaign=='HISCALE':  
    # cpc
    if IOP=='IOP1':
        lst = glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201604*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201605*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201606*')
    elif IOP=='IOP2':
        lst = glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201608*')+glob.glob(cpcsfcpath+'sgpaoscpcC1.b1.201609*')
    lst.sort()
    t_cpc=np.empty(0)
    cpc=np.empty(0)
    if len(lst)==0:
        t_cpc = np.array([np.nan])
        cpc = np.array([np.nan])
    else:
        for filename in lst:
            (time,data,timeunit,cpcunit)=read_cpc(filename)
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date,'noleap')
            # average in time for better comparison with obs
            time2=np.arange(0,86400,3600)
            data2 = avg_time(np.array(time),np.array(data),time2)
            t_cpc=np.hstack((t_cpc, cday+time2/86400))
            cpc=np.hstack((cpc, data2))
        cpc[cpc<0]=np.nan
  
    # cpcu
    if IOP=='IOP1':
        lst = glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201604*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201605*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201606*')
    elif IOP=='IOP2':
        lst = glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201608*')+glob.glob(cpcusfcpath+'sgpaoscpcuS01.b1.201609*')
    lst.sort()
    t_cpcu=np.empty(0)
    cpcu=np.empty(0)
    if len(lst)==0:
        t_cpcu = np.array([np.nan])
        cpcu = np.array([np.nan])
    else:
        for filename in lst:
            (time,data,timeunit,cpcuunit)=read_cpc(filename)
            timestr=timeunit.split(' ')
            date=timestr[2]
            cday=yyyymmdd2cday(date,'noleap')
            # average in time for better comparison with obs
            time2=np.arange(0,86400,3600)
            data2 = avg_time(np.array(time),np.array(data),time2)
            t_cpcu=np.hstack((t_cpcu, cday+time2/86400))
            cpcu=np.hstack((cpcu, data2))
        cpcu[cpcu<0]=np.nan
    
#%% read in models
ncn_m = []
nucn_m = []
nmodels = len(Model_List)
for mm in range(nmodels):
    tmp_ncn=np.empty(0)
    tmp_nucn=np.empty(0)
    timem=np.empty(0)
    for cday in range(cday1,cday2+1):
        mmdd=cday2mmdd(cday)
        date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
        
        filename_input = E3SM_sfc_path+'SFC_CNsize_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
        (time,ncn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NCN')
        (time,nucn,timemunit,dataunit,long_name)=read_E3SM(filename_input,'NUCN')
        
        timem = np.hstack((timem,time))
        tmp_ncn = np.hstack((tmp_ncn,ncn*1e-6))
        tmp_nucn = np.hstack((tmp_nucn,nucn*1e-6))
    
    ncn_m.append(tmp_ncn)
    nucn_m.append(tmp_nucn)
    
#%% calculate diurnal cycle
days = np.arange(cday1,cday2+1)

time_dc = np.arange(30,1440.,60)
cpc_o_dc = np.full((len(time_dc),len(days)),np.nan)
cpcu_o_dc = np.full((len(time_dc),len(days)),np.nan)
n_valid = list()
if len(cpc)>1: # not NaN
    for dd in range(len(days)):
        for tt in range(len(time_dc)):
            time_tmp = days[dd]+time_dc[tt]/1440.
            idx = np.abs(t_cpc-time_tmp).argmin()
            if (t_cpc[idx]-time_tmp)*1440 <= 30:    
                cpc_o_dc[tt,dd] = cpc[idx]
if len(cpcu)>1:
    for dd in range(len(days)):
        for tt in range(len(time_dc)):
            time_tmp = days[dd]+time_dc[tt]/1440.
            idx = np.abs(t_cpcu-time_tmp).argmin()
            if (t_cpcu[idx]-time_tmp)*1440 <= 30:    
                cpcu_o_dc[tt,dd] = cpcu[idx]
cpc_o_dc = np.nanmean(cpc_o_dc,1)
cpcu_o_dc = np.nanmean(cpcu_o_dc,1)

# for E3SM data
ncn_m_dc = []
nucn_m_dc = []
for mm in range(nmodels):
    tmp_ncn = np.full((24,len(days)),np.nan)
    tmp_nucn = np.full((24,len(days)),np.nan)
    for dd in range(len(days)):
        idx=np.logical_and(timem>=days[dd], timem<days[dd]+1)
        tmp_ncn[:,dd] = ncn_m[mm][idx]
        tmp_nucn[:,dd] = nucn_m[mm][idx]
    ncn_m_dc.append(np.nanmean(tmp_ncn,1))
    nucn_m_dc.append(np.nanmean(tmp_nucn,1))
    
#%% make plot
    
figname = figpath_sfc_timeseries+'diurnalcycle_CN_'+campaign+'_'+IOP+'.png'
print('plotting figures to '+figname)

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(6,4))   # figsize in inches
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)   #pad=0.4, w_pad=0.5, h_pad=1.0

ax1.plot(time_dc/60,cpc_o_dc,color='k',linewidth=1,label='CPC(>10nm)')
for mm in range(nmodels):
    ax1.plot(time_dc/60, ncn_m_dc[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
ax1.tick_params(color='k',labelsize=12)
ylim1 = ax1.get_ylim()

ax2.plot(time_dc/60,cpcu_o_dc,color='k',linewidth=1,label='CPC(>3nm)')
for mm in range(nmodels):
    ax2.plot(time_dc/60, nucn_m_dc[mm],color=color_model[mm],linewidth=1, label=Model_List[mm])
ax2.tick_params(color='k',labelsize=12)
ylim2 = ax2.get_ylim()

ax1.set_xlim(0,24)
ax2.set_xlim(0,24)
ax1.set_xticks(np.arange(0,24,3))
ax2.set_xticks(np.arange(0,24,3))

# set ylimit consistent in subplots
# ax1.set_yticks([10,100,1000,10000,100000])
# ax2.set_yticks([10,100,1000,10000,100000])
# ax1.set_yscale('log')
# ax2.set_yscale('log')
# ax1.set_ylim([ylim1[0], ylim2[1]])
# ax2.set_ylim([ylim1[0], ylim2[1]])


ax1.legend(loc='center right', shadow=False, fontsize='medium',bbox_to_anchor=(1.3, .5))
ax2.legend(loc='center right', shadow=False, fontsize='medium',bbox_to_anchor=(1.3, .5))

ax2.set_xlabel('Hour (UTC)',fontsize=12)
ax1.set_title('Aerosol Number Concentration (cm$^{-3}$) '+campaign+' '+IOP,fontsize=14)

fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# plt.close()