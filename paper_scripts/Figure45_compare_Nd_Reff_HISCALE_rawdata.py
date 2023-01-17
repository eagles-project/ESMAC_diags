
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
from esmac_diags.subroutines.time_format_change import yyyymmdd2cday,hhmmss2sec
from esmac_diags.subroutines.time_resolution_change import avg_time_1d
from esmac_diags.subroutines.read_aircraft import read_iwg1, read_wcm, read_mergedSD
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_VISST, calc_clouddepth_VISST
import matplotlib.dates as mdates


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
arsclpath = path_raw + 'HISCALE/obs/profile/arscl/'
visstgridpath = path_raw + 'HISCALE/obs/satellite/visst/grid/'
visstpixpath = path_raw + 'HISCALE/obs/satellite/visst/pix_3x3/'
# path of prepared files
prep_model_path = 'C:/Users/tang357/Downloads/prep_data/' +site+'/model/'
prep_flight_path = 'C:/Users/tang357/Downloads/prep_data/'+site+'/flight/'
prep_sfc_path = 'C:/Users/tang357/Downloads/prep_data/'+site+'/surface/'
prep_sat_path = 'C:/Users/tang357/Downloads/prep_data/'+site+'/satellite/'
# set output path for plots
figpath= '../figures/'+site+'/'

#%% read aircraft data 
# time_iwg = np.empty(0)
# lon = np.empty(0)
# lat = np.empty(0)
# height = np.empty(0)
# cldflag = np.empty(0)
time_iwg = []
lon = []
lat = []
height = []
cldflag = []
lst = glob.glob(iwgpath + '*.a2.txt')
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

#%% read aircraft data and calculate Reff
# time_mergesd = np.array([])
# nd_mergeSD = np.empty((0))
# reff_air = np.array([])
# lst1 = glob.glob(mergeSDpath + 'aaf.g1.hiscale.mergedSD.*.txt')
# lst1.sort()
# for filename in lst1:
#     fname = filename.split('.')
#     date = fname[-2]
#     print(filename)
#     (time1, Nd, Ndsize, dmean, dmin, dmax) = read_mergedSD(filename)
#     if date=='20160425a' or date=='20160427a' or date=='20160428a':
#         dmax[11] = 24.
#         dmin[11] = 21.
#         dmean[11] = 22.5
#     elif date=='20160503a' or date=='20160506a' or date=='20160507a':
#         Ndsize = np.vstack((Ndsize[0:12,:], Ndsize[13:,:]))
#         dmin = np.hstack((dmin[:12], dmin[13:]))
#         dmax = np.hstack((dmax[:12], dmax[13:]))
#         dmean = np.hstack((dmean[:12], dmean[13:]))        
#     Nd = qc_remove_neg(Nd, remove_zero='True')
#     Ndsize = qc_remove_neg(Ndsize)
#     # unit change from dNd/Dp to Nd
#     Dp = dmax-dmin
#     for tt in range(Ndsize.shape[1]):   
#         Ndsize[:,tt] = Ndsize[:,tt]*Dp
#     time_mergesd = np.hstack((time_mergesd, time1/86400+yyyymmdd2cday(date[0:8])))
#     nd_mergeSD = np.hstack((nd_mergeSD, Nd))
#     # Nd_bin = np.hstack((Nd_bin, Ndsize))
#     # calculate effective radius from observed size distribution
#     reff_a = 0.5 * np.nansum(Ndsize.T * (dmean**3), axis=1) / np.nansum(Ndsize.T * (dmean**2), axis=1)
#     reff_air = np.hstack((reff_air, reff_a))
    
# # unit change
# nd_mergeSD = nd_mergeSD/1000 # #/L to #/cm3
# # Nd_bin = Nd_bin/1000 # #/L to #/cm3

# #%% save time_mergesd, nd_mergeSD and reff_air for quick read
# np.savetxt('time_HISCALE.txt',time_mergesd)
# np.savetxt('nd_HISCALE.txt',nd_mergeSD)
# np.savetxt('reff_HISCALE.txt',reff_air)

time_mergesd = np.loadtxt('time_HISCALE.txt')
nd_mergeSD = np.loadtxt('nd_HISCALE.txt')
reff_air = np.loadtxt('reff_HISCALE.txt')

#%% LWC
time_lwc = np.empty((0))
lwc_all = np.empty((0))
lst2 = glob.glob(wcmpath+'WCM_G1_*')
lst2.sort()
for filename in lst2:
    fname = filename.split('_')
    date = fname[-3]
    (wcm,wcmlist)=read_wcm(filename)
    time=wcm[0,:]
    flag=wcm[-1,:]
    twcobs=wcm[1,:]
    lwcobs=wcm[2,:]
    # quality controls
    twcobs=qc_mask_qcflag(twcobs,flag)
    lwcobs=qc_mask_qcflag(lwcobs,flag)
    twcobs=qc_remove_neg(twcobs)
    lwcobs=qc_remove_neg(lwcobs)
    time_lwc = np.hstack((time_lwc, time/86400+yyyymmdd2cday(date[0:8])))
    lwc_all = np.hstack((lwc_all, lwcobs))
    
    
#%% surface Ndrop
lst = glob.glob(os.path.join(ndroppath, '*ndropmfrsrC1.c1*.nc'))
lst.sort()
obsdata = xr.open_mfdataset(lst, combine='by_coords')
time_sfc = obsdata['time'].load()
nd_sfc = obsdata['drop_number_conc'].load()
qc_nd = obsdata['qc_drop_number_conc'].load()
cth_sfc = obsdata['cloud_top_height'].load()
dH_nd = obsdata['cloud_thickness'].load()
obsdata.close()
# quality control
nd_sfc = qc_mask_qcflag(nd_sfc,qc_nd)
nd_sfc = nd_sfc*1e-6   # m-3 to cm-3

# Reff surface
lst = glob.glob(os.path.join(mfrsrpath, '*.cdf'))
lst.sort()
# first data
mfrsrdata = xr.open_dataset(lst[0])
mfrsrtime = mfrsrdata['time']
reff_sfc = mfrsrdata['effective_radius_instantaneous']
qc_reff = mfrsrdata['qc_effective_radius_instantaneous']
reffa_sfc = mfrsrdata['effective_radius_average']
qc_reffa = mfrsrdata['qc_effective_radius_average']
cf_sfc = mfrsrdata['cloudfraction']
mfrsrdata.close()
for file in lst[1:]:
    mfrsrdata = xr.open_dataset(file)
    mfrsrtime = xr.concat([mfrsrtime, mfrsrdata['time']], dim="time")
    reff_sfc = xr.concat([reff_sfc, mfrsrdata['effective_radius_instantaneous']], dim="time")
    qc_reff = xr.concat([qc_reff, mfrsrdata['qc_effective_radius_instantaneous']], dim="time")
    reffa_sfc = xr.concat([reff_sfc, mfrsrdata['effective_radius_average']], dim="time")
    qc_reffa = xr.concat([qc_reff, mfrsrdata['qc_effective_radius_average']], dim="time")
    cf_sfc = xr.concat([cf_sfc, mfrsrdata['cloudfraction']], dim="time")
    mfrsrdata.close()
# quality controls
reff_sfc.load()
qc_reff.load()
reff_sfc = qc_mask_qcflag(reff_sfc, qc_reff)
reffa_sfc = qc_mask_qcflag(reffa_sfc, qc_reffa)
cf_sfc[cf_sfc<0] = np.nan

#%% satellite Ndrop
lst = glob.glob(os.path.join(visstgridpath, '*visstgrid*.cdf'))
filetime = [a.split('.c1.')[1] for a in lst]
sortidx = np.argsort(filetime)
# first data
# site=='SGP':
x_idx = 15
y_idx = 10
visstdata = xr.open_dataset(lst[sortidx[0]])
vissttime = visstdata['time']
reff_sat = visstdata['particle_size'][:,y_idx,x_idx,1]
lwp_sat = visstdata['water_path'][:,y_idx,x_idx,1]
ctt_liq = visstdata['cloud_temperature'][:,y_idx,x_idx,2]
cod_liq_linavg = visstdata['optical_depth_linear'][:,y_idx,x_idx,2]
cf_sat = visstdata['cloud_percentage'][:,y_idx,x_idx,0]
cth_sat = visstdata['cloud_height_top'][:,y_idx,x_idx,0]
solar_angle = visstdata['solar_zenith_angle'][:,y_idx,x_idx]
visstdata.close()
for ii in range(1,len(lst)):
    file = lst[sortidx[ii]]
    visstdata = xr.open_dataset(file)
    vissttime = xr.concat([vissttime, visstdata['time']], dim="time")
    reff_sat = xr.concat([reff_sat, visstdata['particle_size'][:,y_idx,x_idx,1]], dim="time")
    lwp_sat = xr.concat([lwp_sat, visstdata['water_path'][:,y_idx,x_idx,1]], dim="time")
    cod_liq_linavg = xr.concat([cod_liq_linavg, visstdata['optical_depth_linear'][:,y_idx,x_idx,2]], dim="time")
    ctt_liq = xr.concat([ctt_liq, visstdata['cloud_temperature'][:,y_idx,x_idx,2]], dim="time")
    cf_sat = xr.concat([cf_sat, visstdata['cloud_percentage'][:,y_idx,x_idx,0]], dim="time")
    cth_sat = xr.concat([cth_sat, visstdata['cloud_height_top'][:,y_idx,x_idx,0]], dim="time")
    solar_angle = xr.concat([solar_angle, visstdata['solar_zenith_angle'][:,y_idx,x_idx]], dim="time")
    visstdata.close()

# calculate Nd_sat
lwp = lwp_sat.data * 0.001 # to kg/m2
ctt = ctt_liq.data
cod = cod_liq_linavg.data
nd_sat_0 = calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=0.8)
# Nd_ad = calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=1.0)

# change nd_sat to xarray Dataset
nd_sat = xr.DataArray(
    data=nd_sat_0,dims=["time"],
    coords=dict(time=vissttime.data),
    attrs=dict(description="Nd",units="cm-3",),
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

# remove satellite measurements for large solar zenith angle
nd_sat[solar_angle>70] = np.nan
reff_sat[solar_angle>70] = np.nan

# change to iwg time
nd_air = np.interp(time_iwg, time_mergesd, nd_mergeSD)
reff_air = np.interp(time_iwg, time_mergesd, reff_air)
lwc_nd = np.interp(time_iwg, time_lwc, lwc_all)

reff_air[np.logical_or(lwc_nd<0.01, nd_air<1)] = np.nan
reff_air[reff_air>50] = np.nan
nd_air[np.logical_or(lwc_nd<0.01, nd_air<1)]=np.nan

# nd_sfc[nd_sfc<1]=np.nan
# nd_sat[nd_sat<1]=np.nan

time2 = time_iwg  * 86400
time_air = np.datetime64('2015-12-31T00') + time2.astype('timedelta64[s]')


#%% all data
w0 = np.ones_like(nd_air)/sum(~np.isnan(nd_air.data))
w1 = np.ones_like(nd_sfc)/sum(~np.isnan(nd_sfc.data))
w2 = np.ones_like(nd_sat)/sum(~np.isnan(nd_sat.data))
fig,ax = plot.hist([nd_air, nd_sfc, nd_sat], weights=[w0, w1, w2], bins=np.arange(0,810,40), 
                    legend = ['aircraft ('+format(1/w0[0],'.0f')+')','Ndrop ('+format(1/w1[0],'.0f')+')','VISST ('+format(1/w2[0],'.0f')+')'], 
                    xlabel='cm$^{-3}$', ylabel='Fraction', color=['k','r','b'],  ylimit=(0,0.4),
                     title = '(a) All data in original resolution' )
                    # title = 'Cloud Droplet Number Concentration '+site, )

#%% select only for low-level overcast clouds
nd_air_xr = xr.DataArray(
    data=nd_air,dims=["time"],
    coords=dict(time=time_air),
    attrs=dict(description="Nd",units="cm-3",),
)
nd_air_incld = nd_air_xr.rolling(time=9, center=True).median()
# nd_sfc_st = nd_sfc[cth_sfc<4000].rolling(time=9, center=True).median()
nd_sfc_st = nd_sfc[np.logical_and(cf_sfc>0.90, cth_sfc<4000)]
nd_sat_st = nd_sat[np.logical_and(cf_sat>90, cth_sat<4)]

w0 = np.ones_like(nd_air_incld)/sum(~np.isnan(nd_air_incld.data))
w1 = np.ones_like(nd_sfc_st)/sum(~np.isnan(nd_sfc_st.data))
w2 = np.ones_like(nd_sat_st)/sum(~np.isnan(nd_sat_st.data))
fig,ax = plot.hist([nd_air_incld, nd_sfc_st, nd_sat_st], weights=[w0, w1, w2], bins=np.arange(0,810,40), 
                    legend = ['aircraft ('+format(1/w0[0],'.0f')+')','Ndrop ('+format(1/w1[0],'.0f')+')','VISST ('+format(1/w2[0],'.0f')+')'], 
                    xlabel='cm$^{-3}$', ylabel='Fraction', color=['k','r','b'],  ylimit=(0,0.4),
                     title = '(b) Overcast low clouds in original resolution' )
                    # title = 'Cloud Droplet Number Concentration '+site, )

#%% remove small Nd and average in VISST resolution, overcast condition
nd_sat_5 = np.array(nd_sat_st)
nd_sfc_1 = np.array(nd_sfc)
nd_air_1 = np.array(nd_air_incld)
# nd_sat_5[nd_sat_5<10] = np.nan
# nd_air_1[nd_air_1<10] = np.nan
# nd_sfc_1[nd_sfc_1<10] = np.nan
nd_sfc_1[np.logical_or(cf_sfc<0.90, cth_sfc>4000)] = np.nan
nd_air_5 = avg_time_1d(time_air, nd_air_1, vissttime.data)
nd_sfc_5 = avg_time_1d(time_sfc.data, nd_sfc_1, vissttime.data)
cf_sfc_5 = avg_time_1d(time_sfc.data, cf_sfc, vissttime.data)
nd_sfc_5[cf_sfc_5<0.90] = np.nan

w0 = np.ones_like(nd_air_5)/sum(~np.isnan(nd_air_5))
w1 = np.ones_like(nd_sfc_5)/sum(~np.isnan(nd_sfc_5.data))
w2 = np.ones_like(nd_sat_5)/sum(~np.isnan(nd_sat_5.data))
# nd_bins = np.array([0,2,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,30,35,40,45,55])
fig,ax = plot.hist([nd_air_5, nd_sfc_5, nd_sat_5], weights=[w0, w1, w2], bins=np.arange(0,810,40), 
                   legend = ['aircraft ('+format(1/w0[0],'.0f')+')','Ndrop ('+format(1/w1[0],'.0f')+')','VISST ('+format(1/w2[0],'.0f')+')'], 
                   ylabel='Fraction', xlabel='cm$^{-3}$', color=['k','r','b'], ylimit=(0,0.4), 
                    title = '(c) Overcast low clouds averaged into 30-min' )
                     # title = 'Cloud Droplet Number Concentration '+site, )

#%% all data
w0 = np.ones_like(reff_air)/sum(~np.isnan(reff_air.data))
w1 = np.ones_like(reff_sfc)/sum(~np.isnan(reff_sfc.data))
w2 = np.ones_like(reff_sat)/sum(~np.isnan(reff_sat.data))
fig,ax = plot.hist([reff_air, reff_sfc, reff_sat], weights=[w0, w1, w2], 
                    bins=np.arange(0,41), xlimit=(0,30), ylimit=(0,0.3),
                    legend = ['aircraft ('+format(1/w0[0],'.0f')+')','MFRSR ('+format(1/w1[0],'.0f')+')',\
                              'VISST ('+format(1/w2[0],'.0f')+')'], 
                    color=['k','r','b'], title = '(a) All data in original resolution', #title='Cloud Effective Radius '+site, 
                    ylabel='Fraction', xlabel='$\mu$m')
    
#%% select only for low-level overcast clouds
reff_air_xr = xr.DataArray(
    data=reff_air,dims=["time"],
    coords=dict(time=time_air),
    attrs=dict(description="reff",units="um",),
)
reff_air_incld = reff_air_xr.rolling(time=9, center=True).median()
reff_sfc_st = reff_sfc[np.logical_and(cf_sfc>0.90, cth_sfc<4000)]
reff_sat_st = reff_sat[np.logical_and(cf_sat>90, cth_sat<4)]

w0 = np.ones_like(reff_air_incld)/sum(~np.isnan(reff_air_incld.data))
w1 = np.ones_like(reff_sfc_st)/sum(~np.isnan(reff_sfc_st.data))
w3 = np.ones_like(reff_sat_st)/sum(~np.isnan(reff_sat_st.data))
fig,ax = plot.hist([reff_air_incld, reff_sfc_st,  reff_sat_st], weights=[w0, w1,  w3], 
                    bins=np.arange(0,41), xlimit=(0,30), ylimit=(0,0.3),
                    legend = ['aircraft ('+format(1/w0[0],'.0f')+')','MFRSR ('+format(1/w1[0],'.0f')+')',\
                              'VISST ('+format(1/w3[0],'.0f')+')'], 
                    color=['k','r','b'], title = '(b) Overcast low clouds in original resolution', #title='Cloud Effective Radius '+site, 
                    ylabel='Fraction', xlabel='$\mu$m')
    
#%% remove small Nd and average in VISST resolution, in overcast condition
reff_sat_5 = np.array(reff_sat)
reff_sfc_1 = np.array(reff_sfc)
reff_air_1 = np.array(reff_air)
# reff_sat_5[nd_sat<10] = np.nan
# reff_air_1[nd_air<10] = np.nan
# reff_sfc_1[nd_sfc<10] = np.nan
reff_sfc_1[np.logical_or(cf_sfc<0.90, cth_sfc>4000)] = np.nan
reff_sat_5[np.logical_or(cf_sat<90, cth_sat>4)] = np.nan
reff_air_2 = xr.DataArray(
    data=reff_air_1,dims=["time"],
    coords=dict(time=time_air),
    attrs=dict(description="reff",units="um",),
)
reff_air_3 = reff_air_2.rolling(time=9, center=True).median()
reff_air_5 = avg_time_1d(time_air, reff_air_3, vissttime.data)
reff_sfc_5 = avg_time_1d(time_sfc.data, reff_sfc_1, vissttime.data)
cf_sfc_5 = avg_time_1d(time_sfc.data, cf_sfc, vissttime.data)
reff_sfc_5[cf_sfc_5<0.90] = np.nan

w0 = np.ones_like(reff_air_5)/sum(~np.isnan(reff_air_5))
w1 = np.ones_like(reff_sfc_5)/sum(~np.isnan(reff_sfc_5.data))
w3 = np.ones_like(reff_sat_5)/sum(~np.isnan(reff_sat_5.data))
fig,ax = plot.hist([reff_air_5, reff_sfc_5, reff_sat_5], weights=[w0, w1, w3], 
                    bins=np.arange(0,35), xlimit=(0,30), ylimit=(0,0.3),
                    legend = ['aircraft ('+format(1/w0[0],'.0f')+')','MFRSR ('+format(1/w1[0],'.0f')+')',\
                              'VISST ('+format(1/w3[0],'.0f')+')'], 
                    color=['k','r','b'], title = '(c) Overcast low clouds averaged into 30-min',#title='Cloud Effective Radius '+site, 
                    ylabel='Fraction', xlabel='$\mu$m')


