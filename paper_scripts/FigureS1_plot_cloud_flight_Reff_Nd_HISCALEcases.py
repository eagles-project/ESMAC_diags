
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import esmac_diags
from esmac_diags.subroutines.time_format_change import yyyymmdd2cday,hhmmss2sec
from esmac_diags.subroutines.time_resolution_change import avg_time_1d
from esmac_diags.subroutines.read_aircraft import read_iwg1, read_wcm, read_mergedSD
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_VISST, calc_clouddepth_VISST
import matplotlib.dates as mdates
from CBcolors import CB_color_cycle

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'
path_raw = 'C:/Users/tang357/OneDrive - PNNL/EAGLES/python_diag_pkg/ESMAC_Diags_Tool/data/'

# raw data path
iwgpath = path_raw + 'HISCALE/obs/aircraft/mei-iwg1/'
mergeSDpath = path_raw + 'HISCALE/obs/aircraft/mergedSD/'
wcmpath = path_raw + 'HISCALE/obs/aircraft/matthews-wcm/'
mfrsrpath = path_raw + 'HISCALE/obs/surface/arm_mfrsr/'
ndroppath = path_raw + 'HISCALE/obs/surface/sgpndrop/'
visstgridpath = path_raw + 'HISCALE/obs/satellite/visst/grid/'
visstpixpath = path_raw + 'HISCALE/obs/satellite/visst/pix_3x3/'
# path of prepared files
prep_model_path = '../prep_data/' +site+'/model/'
prep_flight_path = '../prep_data/'+site+'/flight/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
modis_path = 'C:/Users/tang357/OneDrive - PNNL/EAGLES/Hi-Scale/data/MODIS/'
# set output path for plots
figpath= '../figures/'

# date8 = '20160502'
# date8 = '20160516'
# date8 = '20160519'
date8 = '20160520'
# date8 = '20160829'
# date8 = '20160901'
# date8 = '20160907'
# date8 = '20160915'
date10 = date8[0:4]+'-'+date8[4:6]+'-'+date8[6:8]
if date8[4:6]=='05':
    time_start = date10+'T16'
    time_end = date10+'T22'
else:
    time_start = date10+'T13'
    time_end = date10+'T20'


#%% read preprocessed cloud fraction and aircraft location data
# lst = glob.glob(prep_flight_path + 'WCM_'+site+'_'+date8+'*.nc')
# obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
# time_air = obsdata['time'].load()
# z_air = obsdata['height'].load()
# lwc = obsdata['LWC'].load()
# obsdata.close()

filename = prep_sfc_path + 'cloud_2d_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cf = obsdata['time'].load()
height_cf = obsdata['height'].load()
cf_obs = obsdata['cloud'].load()
obsdata.close()

#%% read aircraft data and calculate Reff

# cloud size
time_mergesd = np.array([])
nd_mergeSD = np.array([])
reff_air = np.array([])
lst1 = glob.glob(mergeSDpath + 'aaf.g1.hiscale.mergedSD.'+date8+'*.txt')
lst1.sort()
for filename in lst1:
    fname = filename.split('.')
    date = fname[-2]
    print(filename)
    (time1, Nd, Ndsize, dmean, dmin, dmax) = read_mergedSD(filename)
    if date=='20160425a' or date=='20160427a' or date=='20160428a':
        dmax[11] = 24.
        dmin[11] = 21.
        dmean[11] = 22.5
    elif date=='20160503a' or date=='20160506a' or date=='20160507a':
        Ndsize = np.vstack((Ndsize[0:12,:], Ndsize[13:,:]))
        dmin = np.hstack((dmin[:12], dmin[13:]))
        dmax = np.hstack((dmax[:12], dmax[13:]))
        dmean = np.hstack((dmean[:12], dmean[13:]))
    print(dmean[8:14])
        
    Nd = qc_remove_neg(Nd, remove_zero='True')
    Ndsize = qc_remove_neg(Ndsize)
    # unit change from dNd/Dp to Nd
    Dp = dmax-dmin
    for tt in range(Ndsize.shape[1]):   
        Ndsize[:,tt] = Ndsize[:,tt]*Dp
    time_mergesd = np.hstack((time_mergesd, time1/86400+yyyymmdd2cday(date[0:8])))
    nd_mergeSD = np.hstack((nd_mergeSD, Nd))
    # calculate effective radius from observed size distribution
    reff_a = 0.5 * np.nansum(Ndsize.T * (dmean**3), axis=1) / np.nansum(Ndsize.T * (dmean**2), axis=1)
    reff_air = np.hstack((reff_air, reff_a))
# unit change
nd_mergeSD = nd_mergeSD/1000 # #/L to #/cm3

# LWC
time_lwc = np.empty((0))
lwc_all = np.empty((0))
lst2 = glob.glob(wcmpath+'WCM_G1_'+date8+'*')
lst2.sort()
for filename in lst2:
    fname = filename.split('_')
    date = fname[-3]
    (wcm,wcmlist)=read_wcm(filename)
    time2=wcm[0,:]
    flag=wcm[-1,:]
    twcobs=wcm[1,:]
    lwcobs=wcm[2,:]
    # quality controls
    twcobs=qc_mask_qcflag(twcobs,flag)
    lwcobs=qc_mask_qcflag(lwcobs,flag)
    twcobs=qc_remove_neg(twcobs)
    lwcobs=qc_remove_neg(lwcobs)
    time_lwc = np.hstack((time_lwc, time2/86400+yyyymmdd2cday(date[0:8])))
    lwc_all = np.hstack((lwc_all, lwcobs))

# location of flight
time_iwg = []
lon = []
lat = []
height = []
cldflag = []
lst = glob.glob(iwgpath + 'aaf.iwg1001s.g1.hiscale.'+date8+'*.a2.txt')
lst.sort()
for filename in lst:
    fname = filename.split('.')
    date = fname[-3]
    cday = yyyymmdd2cday(date[0:8])
    (iwg, iwgvars) = read_iwg1(filename)
    timelen = len(iwg)
    # get lat, lon, height, time
    for t in range(timelen):
        lat1 = float(iwg[t][2])
        lon1 = float(iwg[t][3])
        height1 = float(iwg[t][4])
        cldflag1 = int(iwg[t][35])
        timestr = iwg[t][1].split(' ')
        time1 = hhmmss2sec(timestr[1])
        # time_iwg = np.hstack((time_iwg, time1/86400+cday))
        # lon = np.hstack((lon, lon1))
        # lat = np.hstack((lat, lat1))
        # height = np.hstack((height, height1))
        # cldflag = np.hstack((cldflag, cldflag1))
        time_iwg.append(time1/86400+cday)
        lon.append(lon1)
        lat.append(lat1)
        height.append(height1)
        cldflag.append(cldflag1)
    
time_iwg = np.array(time_iwg)
lon = np.array(lon)
lat = np.array(lat)
height = np.array(height)
cldflag = np.array(cldflag)

    
#%% surface retrievals
# surface Nd
lst = glob.glob(os.path.join(ndroppath, '*ndropmfrsrC1.c1.'+date8+'*.nc'))
lst.sort()
obsdata = xr.open_mfdataset(lst, combine='by_coords')
time_ndrop = obsdata['time'].load()
nd_sfc = obsdata['drop_number_conc'].load()
qc_nd = obsdata['qc_drop_number_conc'].load()
cth_sfc = obsdata['cloud_top_height'].load()
lwp_sfc = obsdata['lwp_meas'].load()
obsdata.close()
# quality control
nd_sfc = qc_mask_qcflag(nd_sfc,qc_nd)
nd_sfc = nd_sfc*1e-6   # m-3 to cm-3

# Reff surface
lst = glob.glob(os.path.join(mfrsrpath, '*.c1.'+date8+'*.cdf'))
lst.sort()
# first data
mfrsrdata = xr.open_dataset(lst[0])
mfrsrtime = mfrsrdata['time']
reff = mfrsrdata['effective_radius_instantaneous']
qc_reff = mfrsrdata['qc_effective_radius_instantaneous']
cf_sfc = mfrsrdata['cloudfraction']
mfrsrdata.close()
for file in lst[1:]:
    mfrsrdata = xr.open_dataset(file)
    mfrsrtime = xr.concat([mfrsrtime, mfrsrdata['time']], dim="time")
    reff = xr.concat([reff, mfrsrdata['effective_radius_instantaneous']], dim="time")
    qc_reff = xr.concat([qc_reff, mfrsrdata['qc_effective_radius_instantaneous']], dim="time")
    cf_sfc = xr.concat([cf_sfc, mfrsrdata['cloudfraction']], dim="time")
    mfrsrdata.close()
# quality controls
reff.load()
qc_reff.load()
reff_sfc = qc_mask_qcflag(reff, qc_reff)


#%% satellite retrievals
lst = glob.glob(os.path.join(visstgridpath, 'sgpvisstgrid*c1.'+date8+'*.cdf'))
filetime = [a.split('.c1.')[1] for a in lst]
sortidx = np.argsort(filetime)
# first data
# site=='SGP':
x_idx = 15
y_idx = 10
visstdata = xr.open_dataset(lst[sortidx[0]])
vissttime = visstdata['time']
reff_sat = visstdata['particle_size'][:,y_idx,x_idx,1]
lwp = visstdata['water_path'][:,y_idx,x_idx,1]
ctt_liq = visstdata['cloud_temperature'][:,y_idx,x_idx,2]
cod_liq_linavg = visstdata['optical_depth_linear'][:,y_idx,x_idx,2]
cf_sat = visstdata['cloud_percentage'][:,y_idx,x_idx,0]
cth_sat = visstdata['cloud_height_top'][:,y_idx,x_idx,0]
visstdata.close()
for ii in range(1,len(lst)):
    file = lst[sortidx[ii]]
    visstdata = xr.open_dataset(file)
    vissttime = xr.concat([vissttime, visstdata['time']], dim="time")
    reff_sat = xr.concat([reff_sat, visstdata['particle_size'][:,y_idx,x_idx,1]], dim="time")
    lwp = xr.concat([lwp, visstdata['water_path'][:,y_idx,x_idx,1]], dim="time")
    cod_liq_linavg = xr.concat([cod_liq_linavg, visstdata['optical_depth_linear'][:,y_idx,x_idx,2]], dim="time")
    ctt_liq = xr.concat([ctt_liq, visstdata['cloud_temperature'][:,y_idx,x_idx,2]], dim="time")
    cf_sat = xr.concat([cf_sat, visstdata['cloud_percentage'][:,y_idx,x_idx,0]], dim="time")
    cth_sat = xr.concat([cth_sat, visstdata['cloud_height_top'][:,y_idx,x_idx,0]], dim="time")
    visstdata.close()

# calculate Nd_sat
lwp_sat = lwp.data * 0.001 # to kg/m2
ctt = ctt_liq.data
cod = cod_liq_linavg.data
nd_sat = calc_cdnc_VISST(lwp_sat, ctt, cod, adiabaticity=0.8)


# MODIS data

# modisfile = modis_path + 'Reff_CTH_HISCALE_0.5x0.5.nc'
# modisdata = xr.open_dataset(modisfile)
# modistime = modisdata['time']
# reff_21 = modisdata['reff']
# reff_16 = modisdata['reff_16']
# reff_37 = modisdata['reff_37']
# reff_1621 = modisdata['reff_1621']
# reff_pcl = modisdata['reff_pcl']
# modisdata.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

# height_air = np.interp(time_mergesd, time_iwg, height)
# lwc_nd = np.interp(time_mergesd, time_lwc, lwc_all)
# nd_air = np.array(nd_mergeSD)
lwc_nd = np.interp(time_iwg, time_lwc, lwc_all)
nd_air = np.interp(time_iwg, time_mergesd, nd_mergeSD,left=np.nan,right=np.nan)
reff_air = np.interp(time_iwg, time_mergesd, reff_air,left=np.nan,right=np.nan)
# only select sample in reasonable Nd and LWC range
idx_air = np.logical_or(lwc_nd<0.01, nd_air<1)
reff_air[idx_air] = np.nan
reff_air[reff_air>50] = np.nan

nd_air[nd_air<1]=np.nan
nd_sfc[nd_sfc<1]=np.nan
nd_sat[nd_sat<1]=np.nan
nd_air[lwc_nd<0.01]=np.nan

# # exclude flight far from SGP site
# radius = 0.5   # threshold distance from SGP, in degree lat/lon
# dist = np.sqrt((lon+97.48792)**2 + (lat-36.6059)**2)
# dist_nd = np.interp(time_mergesd, time_iwg, dist)
# reff_air[dist_nd > radius] = np.nan
# nd_air[dist_nd > radius] = np.nan

cloudfraction = cf_obs.data
cloudfraction[cloudfraction<5]=np.nan

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot

# time2 = time_mergesd  * 86400
# mergesdtime = np.datetime64('2015-12-31T00') + time2.astype('timedelta64[s]')
# time2 = time_iwg  * 86400
# iwgtime = np.datetime64('2015-12-31T00') + time2.astype('timedelta64[s]')

# plt.rcParams.update({'font.size': 16})
# fig = plt.figure(figsize=(8,8))

# ax1 = fig.add_subplot(3,1,1)
# h0=ax1.contourf(time_cf, height_cf*0.001, cloudfraction.T,np.arange(0,101,10),cmap='jet')
# ax1.plot(iwgtime,height*0.001,'k.')
# ax1.set_ylim(0,4)
# ax1.set_ylabel('Height (km)')
# ax1.set_title('Cloud Fraction (%) and Flight Height')
# cax = plt.axes([0.96, 0.72, 0.02, 0.2])
# fig.colorbar(h0, cax=cax) 
# ax1.grid(True,linestyle=':')

# ax2 = fig.add_subplot(3,1,2)
# ax2.plot(iwgtime,reff_air,'k.',label='Flight')
# ax2.plot(mfrsrtime,reff_sfc,'r.',label='MFRSR')
# ax2.plot(vissttime,reff_sat,'b-',marker='o',label='VISST')
# # ax2.scatter(modistime,reff_21,300,color='g',linewidth=6,marker='x',label='MODIS')
# # ax2.scatter(modistime,reff_16,100,color='g',linewidth=1,marker='x')
# # ax2.scatter(modistime,reff_37,100,color='g',linewidth=1,marker='x')
# # ax2.scatter(modistime,reff_1621,100,color='g',linewidth=1,marker='x')

# ax2.set_ylim(0,40)
# ax2.set_title('Reff ($\mu$m)')
# # ax2.legend(loc='upper left')
# ax2.grid(True,linestyle=':')

# ax3 = fig.add_subplot(3,1,3)
# ax3.plot(iwgtime,nd_air,'k.',label='Flight')
# ax3.plot(mfrsrtime,nd_sfc,'r.',label='MFRSR')
# ax3.plot(vissttime,nd_sat,'b-',marker='o',label='VISST')
# # ax3.scatter([],[],300,color='g',linewidth=6,marker='x',label='MODIS')
# ax3.set_ylim(-20,1000)
# ax3.set_title('Nd (cm$^{-3}$)')
# ax3.legend(loc='upper left')
# ax3.grid(True,linestyle=':')

# ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# ax1.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
# ax2.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
# ax3.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
# ax1.set_xticklabels([])
# ax2.set_xticklabels([])
# ax3.set_xlabel('Time (UTC) in '+date10,fontsize=18)

# plt.tight_layout()

# ax2.set_xlim(np.datetime64('2016-05-20T00'), np.datetime64('2016-05-20T23'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot averaged into vissttime

time2 = time_mergesd  * 86400
mergesdtime = np.datetime64('2015-12-31T00') + time2.astype('timedelta64[s]')
time2 = time_iwg  * 86400
iwgtime = np.datetime64('2015-12-31T00') + time2.astype('timedelta64[s]')


# mask small Nd then average
nd_sat_3 = np.array(nd_sat)
nd_air_3 = np.array(nd_air)
nd_sfc_3 = np.array(nd_sfc)
reff_sat_3 = np.array(reff_sat)
reff_air_3 = np.array(reff_air)
reff_sfc_3 = np.array(reff_sfc)
nd_sat_3[nd_sat_3<10] = np.nan
nd_air_3[nd_air_3<10] = np.nan
nd_sfc_3[nd_sfc_3<10] = np.nan
reff_sat_3[nd_sat<10] = np.nan
reff_air_3[nd_air<10] = np.nan
reff_sfc_3[nd_sfc<10] = np.nan
# nd_sfc_3[nd_sfc_3>800] = np.nan
# with overcast clouds
nd_air_xr = xr.DataArray(
    data=nd_air_3,dims=["time"],
    coords=dict(time=time_iwg),
    attrs=dict(description="Nd",units="cm-3",),
)
reff_air_xr = xr.DataArray(
    data=reff_air_3,dims=["time"],
    coords=dict(time=time_iwg),
    attrs=dict(description="reff",units="um",),
)
nd_air_incld = nd_air_xr.rolling(time=9, center=True).median()
reff_air_incld = reff_air_xr.rolling(time=9, center=True).median()
# nd_sfc_3[np.logical_or(cf_sfc<0.9, cth_sfc>4000)] = np.nan
# nd_sat_3[np.logical_or(cf_sat<90, cth_sat>4)] = np.nan
# reff_sfc_3[np.logical_or(cf_sfc<0.9, cth_sfc>4000)] = np.nan
# reff_sat_3[np.logical_or(cf_sat<90, cth_sat>4)] = np.nan
nd_sat_5 = np.array(nd_sat_3)
reff_sat_5 = np.array(reff_sat_3)
nd_air_5 = avg_time_1d(iwgtime, nd_air_incld, vissttime.data)
nd_sfc_5 = avg_time_1d(mfrsrtime.data, nd_sfc_3, vissttime.data)
reff_air_5 = avg_time_1d(iwgtime, reff_air_incld, vissttime.data)
reff_sfc_5 = avg_time_1d(mfrsrtime.data, reff_sfc_3, vissttime.data)
# cl_sfc_5 = avg_time_1d(mfrsrtime.data, np.int32(~np.isnan(nd_sfc_3)), vissttime.data)
# nd_sfc_5[cl_sfc_5<0.9] = np.nan
# cl_sfc_5 = avg_time_1d(mfrsrtime.data, np.int32(~np.isnan(reff_sfc_3)), vissttime.data)
# reff_sfc_5[cl_sfc_5<0.9] = np.nan

time_start = iwgtime[0] - np.timedelta64(30,'m')
time_end = iwgtime[-1] + np.timedelta64(30,'m')

#%%
plt.rcParams.update({'font.size': 16})
fig = plt.figure(figsize=(8,10))

ax1 = fig.add_subplot(4,1,1)
h0=ax1.contourf(time_cf, height_cf*0.001, cloudfraction.T,np.arange(0,101,10),cmap='viridis')
ax1.plot(iwgtime,height*0.001,'k.')
ax1.set_ylim(0,4)
ax1.set_ylabel('Height (km)')
ax1.set_title('Cloud Fraction (%) and Flight Height')
cax = plt.axes([0.99, 0.79, 0.02, 0.16])
fig.colorbar(h0, cax=cax) 
ax1.grid(True,linestyle=':')

ax2 = fig.add_subplot(4,1,2)
ax2.plot(mfrsrtime,lwp_sfc,color=CB_color_cycle[1],label='Surface')
ax2.plot(vissttime,lwp_sat,color=CB_color_cycle[0],marker='^',markersize=8,label='Satellite')
ax2.set_ylim(-0.02, 0.41)
ax2.set_title('LWP (kg/m${^2}$)')
ax2.legend(loc='upper left',fontsize=14)
ax2.grid(True,linestyle=':')

ax3 = fig.add_subplot(4,1,3)
ax3.plot(iwgtime,nd_air,color='k',linestyle='none',marker='.',alpha=0.2,markersize=2)
ax3.plot(mfrsrtime,nd_sfc,color=CB_color_cycle[1],linestyle='none',marker='.',alpha=0.4,markersize=4)
ax3.plot(vissttime,nd_air_5,color='k',marker='o',label='Aircraft')
ax3.plot(vissttime,nd_sfc_5,color=CB_color_cycle[1],marker='s',markersize=6,label='MFRSR/Ndrop')
ax3.plot(vissttime,nd_sat_5,color=CB_color_cycle[0],marker='^',markersize=8,label='VISST')
# ax3.scatter([],[],300,color='g',linewidth=6,marker='x',label='MODIS')
ax3.set_ylim(-20,1000)
ax3.set_title('Nd (cm$^{-3}$)')
ax3.grid(True,linestyle=':')

ax4 = fig.add_subplot(4,1,4)
ax4.plot(iwgtime,reff_air,color='k',linestyle='none',marker='.',alpha=0.2,markersize=2)
ax4.plot(mfrsrtime,reff_sfc,color=CB_color_cycle[1],linestyle='none',marker='.',alpha=0.2,markersize=4)
ax4.plot(vissttime,reff_air_5,color='k',marker='o',label='Aircraft')
ax4.plot(vissttime,reff_sfc_5,color=CB_color_cycle[1],marker='s',markersize=6,label='MFRSR/Ndrop')
ax4.plot(vissttime,reff_sat_5,color=CB_color_cycle[0],marker='^',markersize=8,label='VISST')
ax4.set_ylim(0,30)
ax4.set_title('Reff ($\mu$m)')
ax4.legend(loc='upper left',fontsize=14)
ax4.grid(True,linestyle=':')

ax1.xaxis.set_major_locator(mdates.HourLocator(interval=1))
ax2.xaxis.set_major_locator(mdates.HourLocator(interval=1))
ax3.xaxis.set_major_locator(mdates.HourLocator(interval=1))
ax4.xaxis.set_major_locator(mdates.HourLocator(interval=1))
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax4.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax1.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
ax2.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
ax3.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
ax4.set_xlim(np.datetime64(time_start), np.datetime64(time_end))
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax4.set_xlabel('Time (UTC) in '+date10,fontsize=18)

plt.tight_layout()