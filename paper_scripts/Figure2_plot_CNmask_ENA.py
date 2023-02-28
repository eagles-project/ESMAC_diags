"""
script to generate ENA aerosol mask data

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from CBcolors import CB_color_cycle

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings
# set site name and datapath

# set site name.
site = 'ENA'

predatapath = '../prep_data/'+site+'/surface/'
rawdatapath = 'C:/Users/tang357/OneDrive - PNNL/EAGLES/python_diag_pkg/ESMAC_Diags_Tool/test_raw_data/obs/ENA/ENA_AerosolMask/'

#%% read in obs
lst = glob.glob(os.path.join(rawdatapath, 'CPC_ENA_AM_*.txt'))
lst.sort()
dateall = []
timeall = []
cpcvalue = []
maskflag = []
for ll in range(len(lst)):
    filename=lst[ll]
    f = open(filename, 'r')
    for i in range(13):
        h = f.readline()
    h = f.readline()
    varname = h.strip()
    h = f.readline()
    varunit = h.strip()
    print(filename)
    print(varname)
    print(varunit)
    f.readline()
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split()
        dateall.append(columns[0])
        timeall.append(columns[1])
        if columns[2][0] == 'n':
            cpcvalue.append(-9999.0)
        else:
            cpcvalue.append(float(columns[2]))
        maskflag.append(int(columns[3]))
    f.close()

cpcvalue = np.array(cpcvalue)
maskflag = np.array(maskflag,dtype='float')
cpcvalue[cpcvalue<0] = np.nan
maskflag[maskflag<-999] = np.nan
mask_starttime = dateall[0] + ' ' + timeall[0]
mask_endtime = dateall[-1] + ' ' + timeall[-1]
mask_time = pd.date_range(start=mask_starttime,end=mask_endtime,freq="min")  
cpc_o = np.array(cpcvalue)
cpc_m = np.array(cpcvalue)
cpc_m[maskflag!=0] = np.nan

#%%
lst = sorted(glob.glob(predatapath + 'sfc_ACSM_ENA_*.nc'))
obsdata = xr.open_mfdataset(lst)
acsmtime = obsdata['time'].load()
so4_o = obsdata['so4'].load().data
org_o = obsdata['org'].load().data
nh4_o = obsdata['nh4'].load().data
no3_o = obsdata['no3'].load().data
chl_o = obsdata['chl'].load().data
obsdata.close()

lst = sorted(glob.glob(predatapath + 'sfc_CPC_ENA_withmask_*.nc'))
cpc2data = xr.open_mfdataset(lst)
cpc2time = cpc2data['time']
f_valid = cpc2data['f_valid'].load().data
cpc2data.close()

org_m = np.array(org_o)
org_m[f_valid<0.5] = np.nan


#%% plot timeseries and PDFs in one figure
fig = plt.figure(figsize=(15,7))
plt.rcParams.update({'font.size': 13})
ax1 = fig.add_subplot(2, 4, (1,3))
ax1.plot(mask_time,cpc_o,color=CB_color_cycle[0])
ax1.plot(mask_time,cpc_m,color=CB_color_cycle[1])
ax1.set_title('(a) Aerosol number (cm$^{-3}$) from CPC', fontsize=15)
ax1.legend(['contamined','good'])
ax1.grid()
ax1.set_xlim(np.datetime64('2017-10-10'), np.datetime64('2017-10-15'))
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=1))

ax2 = fig.add_subplot(2,4,4)
wt_o = np.ones_like(cpc_o)/np.count_nonzero(~np.isnan(cpc_o))
wt_m = np.ones_like(cpc_m)/np.count_nonzero(~np.isnan(cpc_m))
ax2.hist([cpc_o, cpc_m], bins=np.arange(21)*100, histtype='bar', weights=[wt_o, wt_m],color=CB_color_cycle[0:2])
ax2.set_xlim(-20,2000)
# ax2.set_xlabel('Aerosol number from CPC (cm$^{-3}$)', fontsize=15)
ax2.set_xlabel('cm$^{-3}$')
ax2.set_ylabel('Fraction')
ax2.set_title('(b) Aerosol Number')
ax2.legend(['contamined','good'])
ax2.grid()

ax3 = fig.add_subplot(2,4, (5,7))
ax3.plot(acsmtime,org_o,color=CB_color_cycle[0],marker='.')
ax3.plot(acsmtime,org_m,color=CB_color_cycle[1],marker='.')
ax3.set_title('(c) Total Organic ($\mu$m/m$^3$) from ACSM', fontsize=15)
ax3.legend(['contamined','good'])
ax3.grid()
ax3.set_xlim(np.datetime64('2017-10-10'), np.datetime64('2017-10-15'))
ax3.set_ylim(-0.2, 4)
ax3.xaxis.set_major_locator(mdates.DayLocator(interval=1))

ax4 = fig.add_subplot(2,4,8)
wt_o = np.ones_like(org_o)/np.count_nonzero(~np.isnan(org_o))
wt_m = np.ones_like(org_m)/np.count_nonzero(~np.isnan(org_m))
ax4.hist([org_o, org_m], bins=np.arange(21)*0.075, histtype='bar', weights=[wt_o, wt_m],color=CB_color_cycle[0:2])
ax4.set_xlim(-0.02,1.52)
# ax4.set_xlabel('Total Organic from ACSM ($\mu$m/m$^3$)')
ax4.set_xlabel('$\mu$m/m$^3$')
ax4.set_ylabel('Fraction')
ax4.set_title('(d) Total Organic')
ax4.legend(['contamined','good'])
ax4.grid()

plt.tight_layout()
