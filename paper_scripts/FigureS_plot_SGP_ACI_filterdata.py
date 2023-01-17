"""
script to generate all plots for HISCALE surface data

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
import esmac_diags.plotting.calc_statistics as calc
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings
# set site name and datapath

# set site name.
site = 'SGP'

prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
   
# path of output figures
figpath= '../figures/'+site+'/'

if not os.path.exists(figpath):
    os.makedirs(figpath)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data

lst = sorted(glob.glob(prep_sfc_path + 'sfc_CPC_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_cpc = obsdata['time'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()
obsdata = xr.open_mfdataset(prep_sfc_path + 'sfc_UHSAS_'+site+'_*.nc')
time_uhsas = obsdata['time']
uhsas100 = obsdata['uhsas100'].load()
obsdata.close()
obsdata = xr.open_mfdataset(prep_sfc_path + 'sfc_SMPS_'+site+'_*.nc')
time_smps = obsdata['time']
smps100 = obsdata['smps100'].load()
obsdata.close()
obsdata = xr.open_mfdataset(prep_sfc_path + 'sfc_TDMA_'+site+'_*.nc')
time_tdma = obsdata['time']
tdma100 = obsdata['tdma100'].load()
obsdata.close()
a=xr.combine_by_coords([uhsas100,smps100,tdma100])
time_cn100 = a.time
cn100 = np.nanmean([a.uhsas100, a.smps100, a.tdma100],axis=0)
cn100_all = np.interp(time_cpc, time_cn100, cn100)

lst = sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn2 = obsdata['ccn2_fit'].load()
ccn2m = obsdata['ccn2_m'].load()
ss2 = obsdata['ss2'].load()
ccn1 = obsdata['ccn1_fit'].load()
ccn1m = obsdata['ccn1_m'].load()
ss1 = obsdata['ss1'].load()
obsdata.close()
dss = 0.05
ccn2m[np.abs(ss2-0.2)>dss] = np.nan
ccn2[np.logical_or(ccn2<0, ccn2>1e4)]=np.nan
ccn2_all = ccn2.interp(time=time_cpc)
ccn2m_all = ccn2m.interp(time=time_cpc)
ccn1m[np.abs(ss1-0.1)>dss] = np.nan
ccn1[np.logical_or(ccn1<0, ccn1>1e4)]=np.nan
ccn1_all = ccn1.interp(time=time_cpc)
ccn1m_all = ccn1m.interp(time=time_cpc)

obsdata = xr.open_mfdataset(prep_sfc_path + 'Ndrop_'+site+'*.nc')
time_ndrop = obsdata['time']
ndrop = obsdata['cdnc'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'totcld_'+site+'*.nc')
time_cld = obsdata['time']
cld_armbe = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_sat = obsdata['tot_cld_visst'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'reff_'+site+'*.nc')
time_reff = obsdata['time']
reff = obsdata['reff'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'LWP_'+site+'*.nc')
time_lwp = obsdata['time']
lwp = obsdata['lwp_armbe'].load()
lwp_mfrsr = obsdata['lwp_mfrsr'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'LTS_'+site+'*.nc')
time_lst = obsdata['time']
LTS850 = obsdata['LTS850'].load()
LTS700 = obsdata['LTS700'].load()
thetadiff_cb = obsdata['thetadiff_cb'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'cod_'+site+'*.nc')
time_cod = obsdata['time']
cod = obsdata['cod'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'cloudheight_ARSCL_'+site+'*.nc')
time_reff = obsdata['time']
cth = obsdata['cth'].load()
cths = obsdata['cths'].load()
obsdata.close()

# satellite
satdata = xr.open_mfdataset(prep_sat_path + 'albedo_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
albedo_sat = satdata['albedo'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'cloudfraction_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
cfall_sat = satdata['cldtot'].load()
cflow_sat = satdata['cldlow'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'cloudtop_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
cth_sat = satdata['cth'].load()
ctt_sat = satdata['ctt'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'cod_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
cod_sat = satdata['cod'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'Hcld_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
Hcld_sat = satdata['Hcld'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'LWP_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
lwp_sat = satdata['lwp'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'IWP_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
iwp_sat = satdata['iwp'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'Nd_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
nd_sat = satdata['Nd'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'Reff_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
reff_sat = satdata['reff'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'albedo_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
albedo_sat = satdata['albedo'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'solarzenith_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
zenith_sat = satdata['solar_zenith_angle'].load()
satdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
ccn1_m = modeldata['CCN3'].load()
ccn2_m = modeldata['CCN4'].load()
ncn10_m = modeldata['NCN10'].load()
ncn100_m = modeldata['NCN100'].load()
cod_m = modeldata['cod'].load()
reff_m = modeldata['reff'].load()
lwp_m = modeldata['TGCLDLWP'].load()
iwp_m = modeldata['TGCLDIWP'].load()
nd_m = modeldata['Nd_mean'].load()
nd_arm_m = modeldata['Nd_ARM'].load()
nd_visst_m = modeldata['Nd_VISST'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
cbh_m = modeldata['cbh'].load()
cth_m = modeldata['cth'].load()
lwdnsfc_m = modeldata['FLDS'].load()
lwnetsfc_m = modeldata['FLNS'].load()
lwnettoa_m = modeldata['FLNT'].load()
lwuptoa_m = modeldata['FLUT'].load()
swdnsfc_m = modeldata['FSDS'].load()
swnetsfc_m = modeldata['FSNS'].load()
swdntoa_m = modeldata['SOLIN'].load()
swnettoa_m = modeldata['FSNT'].load()
swuptoa_m = modeldata['FSUTOA'].load()
modeldata.close()
lwupsfc_m = lwnetsfc_m + lwdnsfc_m
swupsfc_m = swdnsfc_m - swnetsfc_m
albedo_m = swuptoa_m/swdntoa_m*100
filename = prep_model_path + 'E3SMv2_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
thetadiff_cb_m = modeldata['thetadiff_cb'].load()
modeldata.close()


#%% specific data treatments
nd_sat.data[zenith_sat.data>70]=np.nan
lwp_sat.data[zenith_sat.data>70]=np.nan

# select low-level ovarcasting clouds
idx_sfc = np.logical_and(np.logical_and(cth<4000, cld_armbe>90), thetadiff_cb<2)
idx_sat = np.logical_and(np.logical_and(cth_sat<4, cfall_sat>90), thetadiff_cb<2)
idx_m1 = np.logical_and(np.logical_and(cth_m<4000, cld_m>90), thetadiff_cb_m<2)

idx_sat[zenith_sat.data>70] = False


idx_sfc[iwp_sat>10] = False
idx_sat[iwp_sat>10] = False
idx_m1[iwp_m.data>10] = False

ndrop[ndrop<10] = np.nan
ndrop[ndrop>800] = np.nan
nd_sat[nd_sat<10] = np.nan
nd_sat[nd_sat>800] = np.nan
nd_m[nd_m<10] = np.nan
nd_m[nd_m>800] = np.nan
nd_arm_m[nd_arm_m<10] = np.nan
nd_arm_m[nd_arm_m>800] = np.nan
nd_visst_m[nd_visst_m<10] = np.nan
nd_visst_m[nd_visst_m>800] = np.nan

lwp[lwp<10] = np.nan
lwp_sat[lwp_sat<10] = np.nan
lwp_m[lwp_m<10] = np.nan

print(sum(~np.isnan(ndrop.data)), sum(~np.isnan(nd_sat.data)), sum(~np.isnan(nd_m.data)))
print(np.sum(idx_sfc.data),np.sum(idx_sat.data),np.sum(idx_m1.data))

#%%
# # Nd vs surface CCN
fig,ax = plot.jointhist([ccn1m_all[:].data, ccn1m_all[:].data, ccn1_m[:].data], 
                        [ndrop[:].data, nd_sat[:].data, nd_m[:].data,], 
                    # xedges=np.arange(0,410,20), yedges=np.arange(0,310,20), 
                    xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
                    normalize_x=True, vmin=0, vmax=0.15,
                    xlabel=['0.1%CCN (cm$^{-3}$)','0.1%CCN (cm$^{-3}$)','0.1%CCN (cm$^{-3}$)'], 
                    ylabel='Nd (cm$^{-3}$)', title=['Surface','Satellite','E3SMv2'], )
for ax0 in ax[0,:]:
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xticks([10,30,100,300,1000])
    ax0.set_yticks([10,30,100,300,1000])
    ax0.set_yticklabels([10,30,100,300,1000])
    ax0.set_xlim(10,999)
for ax0 in ax[1,:]:
    ax0.set_xscale('log')
    ax0.set_xticks([10,30,100,300,1000])
    ax0.set_xticklabels([10,30,100,300,1000])
    ax0.set_xlim(10,999)
regress,sample = calc.linear_regress([np.log(ccn1m_all[:].data), np.log(ccn1m_all[:].data), np.log(ccn1_m[:].data),], 
                        [np.log(ndrop[:].data), np.log(nd_sat[:].data), np.log(nd_m[:].data),],
                        figpath+'linearfit_Nd_CCN_SGPall.txt',legend=['Surface','Satellite','E3SMv2'],
                        labelx='ln(CCN1)',labely='ln(Nd)')
xdata = [ccn1m_all[~np.isnan(ndrop)].data, ccn1m_all[~np.isnan(nd_sat)].data, ccn1_m[~np.isnan(nd_m)].data]
for nn in range(3):
    x = [max(np.nanmin(xdata[nn]),10), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    ax[0,nn].plot(x, y, color='r')

# %% LWP vs Nd
fig,ax = plot.jointhist([ndrop.data[:], nd_sat[:].data, nd_m[:].data], 
                        [lwp[:].data, lwp_sat[:].data, lwp_m[:].data], 
                        #[ndrop.data, nd_sat.data, nd_visst_m.data, nd_visst_m2.data], 
                        xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
                        normalize_x=True, vmin=0, vmax=0.15,
                    xlabel=['Nd (cm$^{-3}$)','Nd (cm$^{-3}$)','Nd (cm$^{-3}$)'], 
                    ylabel='LWP (g/m$^2$)', title=['Surface','Satellite','E3SMv2',], )
for ax0 in ax[0,:]:
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xticks([10,30,100,300])
    ax0.set_yticks([10,30,100,300])
    ax0.set_yticklabels([10,30,100,300])
for ax0 in ax[1,:]:
    ax0.set_xscale('log')
    ax0.set_xticks([10,30,100,300])
    ax0.set_xticklabels([10,30,100,300])
regress,sample = calc.linear_regress([np.log(ndrop[:].data), np.log(nd_sat[:].data), np.log(nd_m[:].data)],
                        [np.log(lwp[:].data), np.log(lwp_sat[:].data), np.log(lwp_m[:].data)], 
                        figpath+'linearfit_LWP_Nd_SGPall.txt',legend=['Surface','Satellite','E3SMv2',],
                        labelx='ln(Nd)',labely='ln(LWP)')
xdata = [ndrop[~np.isnan(lwp)].data, nd_sat[~np.isnan(lwp_sat)].data, nd_m[~np.isnan(lwp_m)].data]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    ax[0,nn].plot(x, y, color='r')
    
# linear regression in different Nd bins
nd_1 = ndrop.data
nd_2 = nd_sat.data.flatten()
nd_3 = nd_m.data.flatten()
lwp_1 = lwp.data
lwp_2 = lwp_sat.data.flatten()
lwp_3 = lwp_m.data.flatten()
nd_1[np.isnan(lwp_1)] = np.nan
nd_2[np.isnan(lwp_2)] = np.nan
nd_3[np.isnan(lwp_3)] = np.nan

nd_thres = 50
regress,sample = calc.linear_regress(
        [np.log(nd_1[nd_1<nd_thres]), np.log(nd_2[nd_2<nd_thres]), np.log(nd_3[nd_3<nd_thres]), ],
        [np.log(lwp_1[nd_1<nd_thres]), np.log(lwp_2[nd_2<nd_thres]), np.log(lwp_3[nd_3<nd_thres])], 
        figpath+'linearfit_LWP_Nd_SGPall.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(Nd)',labely='ln(LWP)')
xdata=[nd_1[nd_1<nd_thres], nd_2[nd_2<nd_thres], nd_3[nd_3<nd_thres],]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    if regress[nn][3]<0.01:
        ax[0,nn].plot(x, y, color='k',linewidth=5,linestyle='--')
    print(regress[nn][2:4])

regress,sample = calc.linear_regress(
        [np.log(nd_1[nd_1>nd_thres]), np.log(nd_2[nd_2>nd_thres]), np.log(nd_3[nd_3>nd_thres]), ],
        [np.log(lwp_1[nd_1>nd_thres]), np.log(lwp_2[nd_2>nd_thres]), np.log(lwp_3[nd_3>nd_thres])], 
        figpath+'linearfit_LWP_Nd_SGPall.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(Nd)',labely='ln(LWP)')
xdata=[nd_1[nd_1>nd_thres], nd_2[nd_2>nd_thres], nd_3[nd_3>nd_thres],]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    if regress[nn][3]<0.01:
        ax[0,nn].plot(x, y, color='k',linewidth=5,linestyle='--')
    print(regress[nn][2:4])

#%% heatmap of albedo to LWP and Nd
xedges=np.exp(np.arange(np.log(10),7,0.5))
yedges=np.exp(np.arange(np.log(10),7,0.5))
fig,ax = plot.heatmap([nd_sat[:].data, nd_visst_m[:].data,], 
                      [lwp_sat[:].data, lwp_m[:].data, ], 
                      [albedo_sat[:].data, albedo_m[:].data, ], 
                    # xedges=np.arange(0,310,20), yedges=np.arange(10,300,20), 
                    xedges=xedges, yedges=yedges, 
                    vmin=20, vmax=60,cmap='plasma',min_sample=4,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', 
                    title=['Satellite','E3SMv2'])
fig.text(.91, .85,'TOA Albedo (%)')      # unit of color variable
xticklabels = (xedges[0:-1] + xedges[1:])/2
yticklabels = (yedges[0:-1] + yedges[1:])/2
xticks = np.interp(np.array([10,30,100,300]),xticklabels,np.arange(len(xticklabels)))
yticks = np.interp(np.array([10,30,100,300]),yticklabels,np.arange(len(yticklabels)))
xticks[0] = -0.5
yticks[0] = -0.5
for ax0 in ax:
    ax0.set_xticks(xticks)
    ax0.set_yticks(yticks)
    ax0.set_xticklabels([10,30,100,300])
    ax0.set_yticklabels([10,30,100,300])


#%%
# # Nd vs surface CCN
fig,ax = plot.jointhist([ccn1m_all[idx_sfc].data, ccn1m_all[idx_sat].data, ccn1_m[idx_m1].data], 
                        [ndrop[idx_sfc].data, nd_sat[idx_sat].data, nd_m[idx_m1].data,], 
                    # xedges=np.arange(0,410,20), yedges=np.arange(0,310,20), 
                    xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
                    normalize_x=True, vmin=0, vmax=0.201,
                    xlabel=['0.1%CCN (cm$^{-3}$)','0.1%CCN (cm$^{-3}$)','0.1%CCN (cm$^{-3}$)'], 
                    ylabel='Nd (cm$^{-3}$)', title=['Surface','Satellite','E3SMv2'], )
for ax0 in ax[0,:]:
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xticks([10,30,100,300,1000])
    ax0.set_yticks([10,30,100,300,1000])
    ax0.set_yticklabels([10,30,100,300,1000])
    ax0.set_xlim(10,999)
for ax0 in ax[1,:]:
    ax0.set_xscale('log')
    ax0.set_xticks([10,30,100,300,1000])
    ax0.set_xticklabels([10,30,100,300,1000])
    ax0.set_xlim(10,999)
regress,sample = calc.linear_regress([np.log(ccn1m_all[idx_sfc].data), np.log(ccn1m_all[idx_sat].data), np.log(ccn1_m[idx_m1].data),], 
                        [np.log(ndrop[idx_sfc].data), np.log(nd_sat[idx_sat].data), np.log(nd_m[idx_m1].data),],
                        figpath+'linearfit_Nd_CCN_SGPfilter.txt',legend=['Surface','Satellite','E3SMv2'],
                        labelx='ln(CCN1)',labely='ln(Nd)')
xdata = [ccn1m_all[idx_sfc].data, ccn1m_all[idx_sat].data, ccn1_m[idx_m1].data]
for nn in range(3):
    x = [max(np.nanmin(xdata[nn]),10), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    ax[0,nn].plot(x, y, color='r')

# %% LWP vs Nd
fig,ax = plot.jointhist([ndrop.data[idx_sfc], nd_sat[idx_sat].data, nd_m[idx_m1].data], 
                        [lwp[idx_sfc].data, lwp_sat[idx_sat].data, lwp_m[idx_m1].data], 
                        xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
                        normalize_x=True, vmin=0, vmax=0.15,
                    xlabel=['Nd (cm$^{-3}$)','Nd (cm$^{-3}$)','Nd (cm$^{-3}$)'], 
                    ylabel='LWP (g/m$^2$)', title=['Surface','Satellite','E3SMv2',], )
for ax0 in ax[0,:]:
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xticks([10,30,100,300])
    ax0.set_yticks([10,30,100,300])
    ax0.set_yticklabels([10,30,100,300])
for ax0 in ax[1,:]:
    ax0.set_xscale('log')
    ax0.set_xticks([10,30,100,300])
    ax0.set_xticklabels([10,30,100,300])
regress,sample = calc.linear_regress([np.log(ndrop[idx_sfc].data), np.log(nd_sat[idx_sat].data), np.log(nd_m[idx_m1].data)],
                        [np.log(lwp[idx_sfc].data), np.log(lwp_sat[idx_sat].data), np.log(lwp_m[idx_m1].data),], 
                        figpath+'linearfit_LWP_Nd_SGPfilter.txt',legend=['Surface','Satellite','E3SMv2',],
                        labelx='ln(Nd)',labely='ln(LWP)')
xdata = [ndrop[idx_sfc].data, nd_sat[idx_sat].data, nd_m[idx_m1].data]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    ax[0,nn].plot(x, y, color='r')
    
# linear regression in different Nd bins
nd_1 = ndrop[idx_sfc].data
nd_2 = nd_sat[idx_sat].data.flatten()
nd_3 = nd_m[idx_m1].data.flatten()
lwp_1 = lwp[idx_sfc].data
lwp_2 = lwp_sat[idx_sat].data.flatten()
lwp_3 = lwp_m[idx_m1].data.flatten()
nd_1[np.isnan(lwp_1)] = np.nan
nd_2[np.isnan(lwp_2)] = np.nan
nd_3[np.isnan(lwp_3)] = np.nan

nd_thres = 50
regress,sample = calc.linear_regress(
        [np.log(nd_1[nd_1<nd_thres]), np.log(nd_2[nd_2<nd_thres]), np.log(nd_3[nd_3<nd_thres]), ],
        [np.log(lwp_1[nd_1<nd_thres]), np.log(lwp_2[nd_2<nd_thres]), np.log(lwp_3[nd_3<nd_thres])], 
        figpath+'linearfit_LWP_Nd_SGPfilter.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(Nd)',labely='ln(LWP)')
xdata=[nd_1[nd_1<nd_thres], nd_2[nd_2<nd_thres], nd_3[nd_3<nd_thres],]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    if regress[nn][3]<0.01:
        ax[0,nn].plot(x, y, color='k',linewidth=5,linestyle='--')
    print(regress[nn][2:4])

regress,sample = calc.linear_regress(
        [np.log(nd_1[nd_1>nd_thres]), np.log(nd_2[nd_2>nd_thres]), np.log(nd_3[nd_3>nd_thres]), ],
        [np.log(lwp_1[nd_1>nd_thres]), np.log(lwp_2[nd_2>nd_thres]), np.log(lwp_3[nd_3>nd_thres])], 
        figpath+'linearfit_LWP_Nd_SGPfilter.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(Nd)',labely='ln(LWP)')
xdata=[nd_1[nd_1>nd_thres], nd_2[nd_2>nd_thres], nd_3[nd_3>nd_thres],]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    if regress[nn][3]<0.01:
        ax[0,nn].plot(x, y, color='k',linewidth=5,linestyle='--')
    print(regress[nn][2:4])

#%% heatmap of albedo to LWP and Nd
xedges=np.exp(np.arange(np.log(10),7,0.5))
yedges=np.exp(np.arange(np.log(10),7,0.5))
fig,ax = plot.heatmap([nd_sat[idx_sat].data, nd_visst_m[idx_m1].data,], 
                      [lwp_sat[idx_sat].data, lwp_m[idx_m1].data, ], 
                      [albedo_sat[idx_sat].data, albedo_m[idx_m1].data, ], 
                    # xedges=np.arange(0,310,20), yedges=np.arange(10,300,20), 
                    xedges=xedges, yedges=yedges, 
                    vmin=20, vmax=60,cmap='plasma',min_sample=4,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', 
                    title=['Satellite','E3SMv2'])
fig.text(.91, .85,'TOA Albedo (%)')      # unit of color variable
xticklabels = (xedges[0:-1] + xedges[1:])/2
yticklabels = (yedges[0:-1] + yedges[1:])/2
xticks = np.interp(np.array([10,30,100,300]),xticklabels,np.arange(len(xticklabels)))
yticks = np.interp(np.array([10,30,100,300]),yticklabels,np.arange(len(yticklabels)))
xticks[0] = -0.5
yticks[0] = -0.5
for ax0 in ax:
    ax0.set_xticks(xticks)
    ax0.set_yticks(yticks)
    ax0.set_xticklabels([10,30,100,300])
    ax0.set_yticklabels([10,30,100,300])
