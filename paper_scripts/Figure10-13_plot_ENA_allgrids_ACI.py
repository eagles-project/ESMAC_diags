"""
script to generate all plots for ENA data - all grids around ENA site

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
site = 'ENA'
if site=='ENA':
    lat0 = 39.09527
    lon0 = -28.0339
elif site=='SGP':
    lat0 = 36.6059
    lon0 = -97.48792
    
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

obsdata = xr.open_mfdataset(prep_sfc_path + 'Ndrop_'+site+'*.nc')
time_ndrop = obsdata['time']
ndrop = obsdata['nd'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'Nd_Reff_Wu_etal_'+site+'*.nc')
time_wu = obsdata['time']
nd_wu = obsdata['nd'].load()
reff_wu = obsdata['reff'].load()
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
cbh = obsdata['cbh'].load()
cths = obsdata['cths'].load()
obsdata.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data

# satellite
satdata = xr.open_mfdataset(prep_sat_path + 'albedo_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
albedo_sat = satdata['albedo'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'cloudfraction_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
cfall_sat = satdata['cldtot'].load()
cflow_sat = satdata['cldlow'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'cloudtop_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
cth_sat = satdata['cth'].load()
ctt_sat = satdata['ctt'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'cod_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
cod_sat = satdata['cod'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'Hcld_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
Hcld_sat = satdata['Hcld'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'LWP_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
lwp_sat = satdata['lwp'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'IWP_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
iwp_sat = satdata['iwp'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'Nd_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
nd_sat = satdata['Nd'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'Reff_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
reff_sat = satdata['reff'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'albedo_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
albedo_sat = satdata['albedo'].load()
satdata.close()

satdata = xr.open_mfdataset(prep_sat_path + 'solarzenith_VISSTallgrid_'+site+'*.nc')
time_sat = satdata['time']
zenith_sat = satdata['solar_zenith_angle'].load()
satdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv2_'+site+'_sfc_allgrids.nc'
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
cld_m = modeldata['CLDTOT'].load()
cbh_m = modeldata['cbh'].load()
cth_m = modeldata['cth'].load()
ctt_m = modeldata['ctt'].load()
Hcld_m = modeldata['clddepth'].load()
cldlow_m = modeldata['CLDLOW'].load()
cldmid_m = modeldata['CLDMED'].load()
cldhgh_m = modeldata['CLDHGH'].load()
swdntoa_m = modeldata['SOLIN'].load()
swnettoa_m = modeldata['FSNT'].load()
swuptoa_m = modeldata['FSUTOA'].load()
thetadiff_cb_m = modeldata['thetadiff_cb'].load()
modeldata.close()
albedo_m = swuptoa_m/swdntoa_m*100

cldlow_m.data[cldlow_m.data<0] = np.nan
cldlow_m.data[cldlow_m.data>100] = np.nan
cldmid_m.data[cldmid_m.data<0] = np.nan
cldmid_m.data[cldmid_m.data>100] = np.nan
cldhgh_m.data[cldhgh_m.data<0] = np.nan
cldhgh_m.data[cldhgh_m.data>100] = np.nan

#%% specific data treatments

# select satellite data within 5x5 degree
lon = satdata['lon'].data
lat = satdata['lat'].data
x_idx = np.logical_and(lon>(lon0-2.5), lon<(lon0+2.5))
y_idx = np.logical_and(lat>(lat0-2.5), lat<(lat0+2.5))
x0 = np.abs(lon[x_idx]-lon0).argmin()
y0 = np.abs(lat[y_idx]-lat0).argmin()

albedo_sat = albedo_sat[:,y_idx,x_idx]
cfall_sat = cfall_sat[:,y_idx,x_idx]
cflow_sat = cflow_sat[:,y_idx,x_idx]
cth_sat = cth_sat[:,y_idx,x_idx]
ctt_sat = ctt_sat[:,y_idx,x_idx]
cod_sat = cod_sat[:,y_idx,x_idx]
Hcld_sat = Hcld_sat[:,y_idx,x_idx]
lwp_sat = lwp_sat[:,y_idx,x_idx]
iwp_sat = iwp_sat[:,y_idx,x_idx]
nd_sat = nd_sat[:,y_idx,x_idx]
reff_sat = reff_sat[:,y_idx,x_idx]
zenith_sat = zenith_sat[:,y_idx,x_idx]

nd_sat.data[zenith_sat.data>65]=np.nan
lwp_sat.data[zenith_sat.data>65]=np.nan

# select low-level ovarcasting clouds
idx_sfc = np.logical_and(cth<4000, cld_armbe>90)
idx_sat = np.logical_and(cth_sat.data<4, cfall_sat.data>90)
idx_m1 = np.logical_and(cth_m.data<4000, cldlow_m.data>90)
# some data has NaNs so have to use the following
idx_sfc[thetadiff_cb.data>2] = False
idx_sat[thetadiff_cb.data>2, :,:] = False
idx_m1[thetadiff_cb_m.data>2] = False
idx_sfc[iwp_sat[:,y0,x0]>10] = False
idx_sat[iwp_sat>10] = False
idx_m1[iwp_m.data>10] = False
idx_sat[zenith_sat.data>65] = False
print(np.sum(idx_sfc.data),np.sum(idx_sat.data),np.sum(idx_m1.data))

ndrop[ndrop<10] = np.nan
ndrop[ndrop>800] = np.nan
nd_wu[nd_wu<10] = np.nan
nd_wu[nd_wu>800] = np.nan
nd_sat.data[nd_sat.data<10] = np.nan
nd_sat.data[nd_sat.data>800] = np.nan
nd_m.data[nd_m.data<10] = np.nan
nd_m.data[nd_m.data>800] = np.nan
nd_arm_m.data[nd_arm_m.data<10] = np.nan
nd_arm_m.data[nd_arm_m.data>800] = np.nan
nd_visst_m.data[nd_visst_m.data<10] = np.nan
nd_visst_m.data[nd_visst_m.data>800] = np.nan

lwp.data[lwp.data<10] = np.nan
lwp_sat.data[lwp_sat.data<10] = np.nan
lwp_m.data[lwp_m.data<10] = np.nan
ccn1_all.data[ccn1_all.data<20] = np.nan
ccn1_m.data[ccn1_m.data<20] = np.nan


#%%
ccn_sat = np.ones_like(nd_sat)
for zz in range(nd_sat.shape[2]):
    for yy in range(nd_sat.shape[1]):
        ccn_sat[:,yy,zz] = ccn1_all[:]
# # Nd vs surface CCN
fig,ax = plot.jointhist([ccn1_all[idx_sfc].data, ccn_sat[idx_sat].flatten(), ccn1_m.data[idx_m1].flatten()], 
                        [ndrop[idx_sfc].data, nd_sat.data[idx_sat].flatten(), nd_m.data[idx_m1].flatten()], 
                    xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
                    normalize_x=True, vmin=0, vmax=0.16,
                    xlabel=['0.1%CCN (cm$^{-3}$)','0.1%CCN (cm$^{-3}$)','0.1%CCN (cm$^{-3}$)'], 
                    ylabel='Nd (cm$^{-3}$)', title=['Ground','Satellite','E3SMv2',], )
for ax0 in ax[0,:]:
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xticks([10,30,100,300,1000])
    ax0.set_yticks([10,30,100,300,1000])
    ax0.set_yticklabels([10,30,100,300,1000])
for ax0 in ax[1,:]:
    ax0.set_xscale('log')
    ax0.set_xticks([10,30,100,300,1000])
    ax0.set_xticklabels([10,30,100,300,1000])
regress,sample = calc.linear_regress(
        [np.log(ccn1_all[idx_sfc].data), np.log(ccn_sat[idx_sat].flatten()), np.log(ccn1_m.data[idx_m1].flatten()),], 
        [np.log(ndrop[idx_sfc].data), np.log(nd_sat.data[idx_sat].flatten()), np.log(nd_m.data[idx_m1].flatten()), ],
        figpath+'linearfit_Nd_CCN_sfccouple.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(CCN1)',labely='ln(Nd)')
xdata = [ccn1_all[idx_sfc].data, ccn_sat[idx_sat].flatten(), ccn1_m.data[idx_m1].flatten()]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    ax[0,nn].plot(x, y, color='r')
    a = regress[nn][0]
    b = regress[nn][1]
    if b<0:
        ax[0,nn].text(50,600, 'y = '+format(a,'3.2f')+'x - '+format(-b,'3.2f'),color='r')
    else:
        ax[0,nn].text(50,600, 'y = '+format(a,'3.2f')+'x + '+format(b,'3.2f'),color='r')
    ax[0,nn].text(100,400, 'R = '+format(regress[nn][2],'3.2f'),color='r')


# %% LWP vs Nd

fig,ax = plot.jointhist([ndrop.data[idx_sfc], nd_sat.data[idx_sat].flatten(), nd_m.data[idx_m1].flatten(),], 
                        [lwp[idx_sfc].data, lwp_sat.data[idx_sat].flatten(), lwp_m.data[idx_m1].flatten(),], 
                    xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
                    normalize_x=True, vmin=0, vmax=0.15,
                    xlabel=['Nd (cm$^{-3}$)','Nd (cm$^{-3}$)','Nd (cm$^{-3}$)',], 
                    ylabel='LWP (g/m$^2$)', title=['Ground','Satellite','E3SMv2',], )
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
    
regress,sample = calc.linear_regress(
        [np.log(ndrop[idx_sfc].data), np.log(nd_sat.data[idx_sat].flatten()), np.log(nd_m.data[idx_m1].flatten()), ],
        [np.log(lwp.data[idx_sfc]), np.log(lwp_sat.data[idx_sat].flatten()), np.log(lwp_m.data[idx_m1].flatten())], 
        figpath+'linearfit_LWP_Nd_sfccouple_allgrids.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(Nd)',labely='ln(LWP)')
xdata=[ndrop.data[idx_sfc], nd_sat.data[idx_sat].flatten(), nd_m.data[idx_m1].flatten(),]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    ax[0,nn].plot(x, y, color='r')
    a = regress[nn][0]
    b = regress[nn][1]
    if b<0:
        ax[0,nn].text(50,600, 'y = '+format(a,'3.2f')+'x - '+format(-b,'3.2f'),color='r')
    else:
        ax[0,nn].text(50,600, 'y = '+format(a,'3.2f')+'x + '+format(b,'3.2f'),color='r')
    ax[0,nn].text(100,400, 'R = '+format(regress[nn][2],'3.2f'),color='r')

# linear regression in different Nd bins
nd_1 = ndrop[idx_sfc].data
nd_2 = nd_sat.data[idx_sat].flatten()
nd_3 = nd_m.data[idx_m1].flatten()
lwp_1 = lwp.data[idx_sfc]
lwp_2 = lwp_sat.data[idx_sat].flatten()
lwp_3 = lwp_m.data[idx_m1].flatten()

nd_thres = 50
regress,sample = calc.linear_regress(
        [np.log(nd_1[nd_1<nd_thres]), np.log(nd_2[nd_2<nd_thres]), np.log(nd_3[nd_3<nd_thres]), ],
        [np.log(lwp_1[nd_1<nd_thres]), np.log(lwp_2[nd_2<nd_thres]), np.log(lwp_3[nd_3<nd_thres])], 
        figpath+'linearfit_LWP_Nd_sfccouple_allgrids.txt',legend=['Surface','Satellite','E3SMv2',],
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
        figpath+'linearfit_LWP_Nd_sfccouple_allgrids.txt',legend=['Surface','Satellite','E3SMv2',],
        labelx='ln(Nd)',labely='ln(LWP)')
xdata=[nd_1[nd_1>nd_thres], nd_2[nd_2>nd_thres], nd_3[nd_3>nd_thres],]
for nn in range(3):
    x = [np.nanmin(xdata[nn]), min(np.nanmax(xdata[nn]),1000)]
    y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
    if regress[nn][3]<0.01:
        ax[0,nn].plot(x, y, color='k',linewidth=5,linestyle='--')
    print(regress[nn][2:4])

    
#%% heatmap of albedo to LWP and Nd

# idx_sat[cod_sat<4] = False
# idx_m1[cod_m<4] = False

xedges=np.exp(np.arange(np.log(10),7,0.5))
yedges=np.exp(np.arange(np.log(10),7,0.5))
fig,ax = plot.heatmap([nd_sat.data[idx_sat].flatten(), nd_m.data[idx_m1].flatten(), ], 
                      [lwp_sat.data[idx_sat].flatten(), lwp_m.data[idx_m1].flatten(),], 
                      [albedo_sat.data[idx_sat].flatten(), albedo_m.data[idx_m1].flatten(),], 
                    xedges=xedges, yedges=yedges, 
                    vmin=20, vmax=60, cmap='plasma',min_sample=10,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', 
                    # title=['Ground','Satellite','E3SMv2','E3SMv2'])
                    title=['Satellite','E3SMv2',])
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

# add sample number PDF
X=np.arange(0,len(xedges)-1)
Y=np.arange(0,len(yedges)-1)
H_count, x, y = np.histogram2d(nd_sat.data[idx_sat].flatten(), lwp_sat.data[idx_sat].flatten(), bins=(xedges, yedges))
cs = ax[0].contour(X,Y,H_count.T,np.array([10, 100, 300, 1000, 3000]), colors=['k'])
ax[0].clabel(cs, cs.levels, inline=True, fontsize=10)
H_count, x, y = np.histogram2d(nd_m.data[idx_m1].flatten(), lwp_m.data[idx_m1].flatten(), bins=(xedges, yedges))
cs1 = ax[1].contour(X,Y,H_count.T,np.array([10, 100, 300, 1000]), colors=['k'])
ax[1].clabel(cs1, cs1.levels, inline=True, fontsize=10)
