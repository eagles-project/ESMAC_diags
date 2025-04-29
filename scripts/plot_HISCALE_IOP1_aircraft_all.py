"""
script to generate all plots for HISCALE aircraft data

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
import esmac_diags.plotting.calc_statistics as calc
import matplotlib.dates as mdates

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings
# set site name and datapath

# set site name.
site = 'HISCALE'

prep_model_path = '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/prep_data/'+site+'/model/flight/60s/'
prep_obs_path = '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/prep_data/'+site+'/flight/60s/'

height_bin = np.arange(300,3200,200)

time1 = np.datetime64('2016-04-25')
time2 = np.datetime64('2016-05-21')
IOP = 'IOP1'

figpath= '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/figures/'+site+'/flight/60s/'
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = sorted(glob.glob(prep_obs_path + 'CPC_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_cpc = obsdata['time'].load()
cpc3 = obsdata['cpc3'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()
cpc3 = cpc3.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
cpc10 = cpc10.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

lst = sorted(glob.glob(prep_obs_path + 'CCN_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn2 = obsdata['CCN2'].load()
ccn5 = obsdata['CCN5'].load()
obsdata.close()
ccn2 = ccn2.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
ccn5 = ccn5.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

lst = sorted(glob.glob(prep_obs_path + 'AMS_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ams = obsdata['time'].load()
org = obsdata['ORG'].load()
so4 = obsdata['SO4'].load()
no3 = obsdata['NO3'].load()
nh4 = obsdata['NH4'].load()
chl = obsdata['CHL'].load()
obsdata.close()
org = org.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
so4 = so4.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
nh4 = nh4.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
no3 = no3.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
chl = chl.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

lst = sorted(glob.glob(prep_obs_path + 'WCM_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_wcm = obsdata['time'].load()
lwc = obsdata['LWC'].load()
obsdata.close()
lwc = lwc.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

lst = sorted(glob.glob(prep_obs_path + 'PCASP100_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_pcasp = obsdata['time'].load()
pcasp100 = obsdata['pcasp100'].load()
obsdata.close()
pcasp100 = pcasp100.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

lst = sorted(glob.glob(prep_obs_path + 'beasd_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_beasd = obsdata['time'].load()
size_beasd = obsdata['size'].load()
sizeh_beasd = obsdata['size_high'].load()
sizel_beasd = obsdata['size_low'].load()
beasd = obsdata['size_distribution_merged'].load()
obsdata.close()
beasd = beasd.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

lst = sorted(glob.glob(prep_obs_path + 'mergedSD_HISCALE_*.nc'))
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_nd = obsdata['time'].load()
size_nd = obsdata['size'].load()
sizeh_nd = obsdata['size_high'].load()
sizel_nd = obsdata['size_low'].load()
nd = obsdata['Nd'].load()
nd_size = obsdata['Nd_bin'].load()
height = obsdata['height'].load()
obsdata.close()
nd = nd.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
nd_size = nd_size.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)
height = height.where(np.logical_and(obsdata.time>time1, obsdata.time<time2), drop=True)

# read in E3SM data
lst = sorted(glob.glob(prep_model_path + 'E3SMv1_HISCALE_flight_*.nc'))
modeldata = xr.open_mfdataset(lst)
time_m = modeldata['time'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn3_m = modeldata['NCN3'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn10_m = modeldata['NCN10'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn100_m = modeldata['NCN100'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn_m = modeldata['NCNall'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ccn1_m = modeldata['CCN3'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ccn2_m = modeldata['CCN4'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ccn5_m = modeldata['CCN5'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
bc_m = modeldata['bc'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
dst_m = modeldata['dst'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncl_m = modeldata['ncl'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
pom_m = modeldata['pom'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
mom_m = modeldata['mom'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
so4_m = modeldata['so4'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
soa_m = modeldata['soa'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
lwc_m = modeldata['LWC'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
reff_m = modeldata['REL'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
nd_m = modeldata['ICWNC'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
nd_bin_m = modeldata['Nd_bin'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
modeldata.close()

lst = sorted(glob.glob(prep_model_path + 'E3SMv2_HISCALE_flight_*.nc'))
modeldata = xr.open_mfdataset(lst)
time_m2 = modeldata['time'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn3_m2 = modeldata['NCN3'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn10_m2 = modeldata['NCN10'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn100_m2 = modeldata['NCN100'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncn_m2 = modeldata['NCNall'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ccn1_m2 = modeldata['CCN3'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ccn2_m2 = modeldata['CCN4'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ccn5_m2 = modeldata['CCN5'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
bc_m2 = modeldata['bc'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
dst_m2 = modeldata['dst'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
ncl_m2 = modeldata['ncl'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
pom_m2 = modeldata['pom'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
mom_m2 = modeldata['mom'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
so4_m2 = modeldata['so4'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
soa_m2 = modeldata['soa'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
lwc_m2 = modeldata['LWC'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
reff_m2 = modeldata['REL'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
nd_m2 = modeldata['ICWNC'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
nd_bin_m2 = modeldata['Nd_bin'].load().where(np.logical_and(modeldata.time>time1, modeldata.time<time2), drop=True)
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# total organic in E3SM
org_m = pom_m+mom_m+soa_m
org_m2 = pom_m2+mom_m2+soa_m2

# unit change
nd = nd/1000 # #/L to #/cm3
nd_size = nd_size/1000 # #/L to #/cm3

# calculate effective radius from observed size distribution
reff = np.sum(nd_size * (size_nd**3), axis=1) / np.sum(nd_size * (size_nd**2), axis=1)
reff[reff>50] = np.nan
    
# remove Nd less than 10 cm-3
nd[nd<10] = np.nan
nd_m[nd_m<10] = np.nan
nd_m2[nd_m2<10] = np.nan
nd_bin_m[:,np.nansum(nd_bin_m,axis=0)<10] = np.nan
nd_bin_m2[:,np.nansum(nd_bin_m2,axis=0)<10] = np.nan
reff_m[np.isnan(nd_m)] = np.nan
reff_m2[np.isnan(nd_m2)] = np.nan

# remove lwc less than 0.02 g/m3
lwc[lwc<0.02] = np.nan
lwc_m[lwc_m<0.02] = np.nan
lwc_m2[lwc_m2<0.02] = np.nan

# size distribution /dlogDp
dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
ncn_m = ncn_m.T/dlogDp_e3sm
ncn_m2 = ncn_m2.T/dlogDp_e3sm
dlogDp_beasd = np.log10(sizeh_beasd/sizel_beasd)
beasd = beasd/dlogDp_beasd

dlogDp_nd = np.log10(sizeh_nd/sizel_nd)
nd_size = nd_size/dlogDp_nd
dlogDp_nd_e3sm = np.log10(np.arange(2,1001)/np.arange(1,1000))
nd_bin_m = nd_bin_m.T/dlogDp_nd_e3sm
nd_bin_m2 = nd_bin_m2.T/dlogDp_nd_e3sm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% bar plot
datagroup0 = [org,so4,nh4,no3,chl, [], []]
datagroup1 = [org_m, so4_m, [], [], [], bc_m, dst_m]
datagroup2 = [org_m2, so4_m2, [], [], [], bc_m2, dst_m2]
dataall=[datagroup0, datagroup1, datagroup2,]
labelall = ['Organic', 'SO$_4$', 'NH$_4$', 'NO$_3$', 'Chl', 'BC', 'Dust']
colorall = ['limegreen', 'red', 'lightblue', 'orange', 'cyan', 'k', 'silver']
fig,ax = plot.bar(dataall, datalabel=['Obs','E3SMv1','E3SMv2',], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
                  title='Aerosol Composition  '+site+' '+IOP, varlabel= labelall, colorall=colorall)
fig.savefig(figpath+'bar_composition_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% mean size distribution
fig,ax = plot.mean_size_witherror([size_beasd, np.arange(1,3001), np.arange(1,3001)], 
                        [beasd,ncn_m,ncn_m2], 
                  legend = ['BEASD','E3SMv1', 'E3SMv2'],color=['k','r','b'],
                  marker=['.',None,None], linestyles=['none','-','-'],
                  xlimit=(10, 1e4), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                    title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
fig.savefig(figpath+'mean_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.mean_size_witherror([size_nd, np.arange(1,1000), np.arange(1,1000)], 
                        [nd_size,nd_bin_m,nd_bin_m2], 
                        marker=['.',None,None], linestyles=['none','-','-'],
                        # marker=['.','.','.'], linestyles=['none','none','none'],
                  legend = ['MergedSD','E3SMv1', 'E3SMv2'],color=['k','r','b'], xlimit=(1, 1e3), ylimit=(1e-5,1e3), 
                  xlabel='Diameter ($\mu$m)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                    title = 'Mean Cloud Droplet Size Distribution '+site+' '+IOP)
fig.savefig(figpath+'mean_cloud_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% 1d histogram
w0 = np.ones_like(org)/sum(~np.isnan(org.data))
w1 = np.ones_like(org_m)/sum(~np.isnan(org_m.data))
w2 = np.ones_like(org_m2)/sum(~np.isnan(org_m2.data))
fig,ax = plot.hist([org,org_m,org_m2], weights=[w0,w1,w2], bins=np.arange(0,6.2,0.2),color=['k','r','b'],
                          legend = ['Obs','E3SMv1','E3SMv2'], title='Total Organic '+site+' '+IOP, 
                          ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
w0 = np.ones_like(so4)/sum(~np.isnan(so4.data))
w1 = np.ones_like(so4_m)/sum(~np.isnan(so4_m.data))
w2 = np.ones_like(so4_m2)/sum(~np.isnan(so4_m2.data))
fig,ax = plot.hist([so4,so4_m,so4_m2], weights=[w0,w1,w2], bins=np.arange(0,5.2,0.2),color=['k','r','b'],
                          legend = ['Obs','E3SMv1','E3SMv2'],title='Sulfate '+site+' '+IOP, 
                          ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(ccn2)/sum(~np.isnan(ccn2.data))
w1 = np.ones_like(ccn2_m)/sum(~np.isnan(ccn2_m.data))
w2 = np.ones_like(ccn2_m2)/sum(~np.isnan(ccn2_m2.data))
fig,ax = plot.hist([ccn2,ccn2_m,ccn2_m2], weights=[w0,w1,w2], legend = ['Obs','E3SMv1','E3SMv2'], 
                    color=['k','r','b'], title='CCN (SS=0.2%) '+site+' '+IOP, 
                    bins=np.arange(0,1500,50), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
w0 = np.ones_like(ccn5)/sum(~np.isnan(ccn5.data))
w1 = np.ones_like(ccn5_m)/sum(~np.isnan(ccn5_m.data))
w2 = np.ones_like(ccn5_m2)/sum(~np.isnan(ccn5_m2.data))
fig,ax = plot.hist([ccn5,ccn5_m, ccn5_m2], weights=[w0,w1,w2], legend = ['Obs','E3SMv1','E3SMv2'], title='CCN (SS=0.5%) '+site+' '+IOP, 
                    color=['k','r','b'], bins=np.arange(0,2200,100), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CCN5_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cpc3)/sum(~np.isnan(cpc3.data))
w1 = np.ones_like(ncn3_m)/sum(~np.isnan(ncn3_m.data))
w2 = np.ones_like(ncn3_m2)/sum(~np.isnan(ncn3_m2.data))
fig,ax = plot.hist([cpc3,ncn3_m,ncn3_m2], weights=[w0,w1,w2], bins=np.arange(0,29000,1000), legend = ['Obs','E3SMv1','E3SMv2'], 
                        color=['k','r','b'],title='Aerosol number (>3nm) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CPC3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
w0 = np.ones_like(cpc10)/sum(~np.isnan(cpc10.data))
w1 = np.ones_like(ncn10_m)/sum(~np.isnan(ncn10_m.data))
w2 = np.ones_like(ncn10_m2)/sum(~np.isnan(ncn10_m2.data))
fig,ax = plot.hist([cpc10,ncn10_m,ncn10_m2], weights=[w0,w1,w2], bins=np.arange(0,22000,1000), legend = ['Obs','E3SMv1','E3SMv2'], 
                        color=['k','r','b'], title='Aerosol number (>10nm) '+site+' '+IOP,ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
w0 = np.ones_like(pcasp100)/sum(~np.isnan(pcasp100.data))
w1 = np.ones_like(ncn100_m)/sum(~np.isnan(ncn100_m.data))
w2 = np.ones_like(ncn100_m2)/sum(~np.isnan(ncn100_m2.data))
fig,ax = plot.hist([pcasp100, ncn100_m, ncn100_m2], weights=[w0,w1,w2], bins=np.arange(0,1500,50), legend = ['Obs','E3SMv1','E3SMv2'], 
                        color=['k','r','b'],title='Aerosol number (>100nm) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(lwc)/sum(~np.isnan(lwc.data))
w1 = np.ones_like(lwc_m)/sum(~np.isnan(lwc_m.data))
w2 = np.ones_like(lwc_m2)/sum(~np.isnan(lwc_m2.data))
fig,ax = plot.hist([lwc, lwc_m, lwc_m2], weights=[w0,w1,w2], bins=np.arange(0.01,0.5,0.02), legend = ['Obs','E3SMv1','E3SMv2'], 
                        color=['k','r','b'], title='lwc '+site+' '+IOP, ylabel='Fraction', xlabel=lwc.units)
# ax.set_yscale('log')
# ax.set_ylim(1e-4,1)
fig.savefig(figpath+'hist_lwc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
w0 = np.ones_like(nd)/sum(~np.isnan(nd.data))
w1 = np.ones_like(nd_m)/sum(~np.isnan(nd_m.data))
w2 = np.ones_like(nd_m2)/sum(~np.isnan(nd_m2.data))
fig,ax = plot.hist([nd, nd_m, nd_m2], weights=[w0, w1, w2], bins=np.arange(1,510,20), legend = ['Obs','E3SMv1','E3SMv2'], 
                        color=['k','r','b'], title='Nd '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
# ax.set_yscale('log')
# ax.set_ylim(1e-4,.1)
fig.savefig(figpath+'hist_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(reff)/sum(~np.isnan(reff.data))
w1 = np.ones_like(reff_m)/sum(~np.isnan(reff_m.data))
w2 = np.ones_like(reff_m2)/sum(~np.isnan(reff_m2.data))
fig,ax = plot.hist([reff, reff_m, reff_m2], weights=[w0, w1, w2], bins=np.arange(3,30,1), legend = ['Obs','E3SMv1','E3SMv2'], 
                        color=['k','r','b'], title='Cloud Effective Radius '+site+' '+IOP, ylabel='Fraction', xlabel='$\mu$m')
fig.savefig(figpath+'hist_Reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% percentiles with height
fig,ax = plot.percentile_z([org,org_m,org_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='Organic',color=['k','r','b'],
                      xlabel='$\mu$g/m$^3$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.percentile_z([so4, so4_m, so4_m2], [height,height,height,], 
                      height_bin, figsize=(3,8), title='Sulfate', color=['k','r','b'],
                      xlabel='$\mu$g/m$^3$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([ccn2,ccn2_m,ccn2_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CCN (0.2%)',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.percentile_z([ccn5,ccn5_m,ccn5_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CCN (0.5%)',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CCN5_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([cpc3,ncn3_m,ncn3_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CN (>3nm)',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CN3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.percentile_z([cpc10,ncn10_m,ncn10_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CN (>10nm)',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CN10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.percentile_z([pcasp100, ncn100_m, ncn100_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CN (>100nm)',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([lwc, lwc_m, lwc_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='LWC',color=['k','r','b'],
                      xlabel='g/m$^{3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_lwc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([nd, nd_m, nd_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='Nd',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([reff, reff_m, reff_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='Reff',color=['k','r','b'],
                      xlabel='$\mu$m', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_Reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statistics
# calc.mean_std_percentiles([org,org_m,org_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_ORG_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([so4, so4_m, so4_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_SO4_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn2,ccn2_m,ccn2_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_CCN2_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn5,ccn5_m,ccn5_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_CCN5_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cpc3,ncn3_m,ncn3_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_CPC3_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cpc10,ncn10_m,ncn10_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([pcasp100, ncn100_m, ncn100_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_PCASP100_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([lwc, lwc_m, lwc_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_lwc_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([nd, nd_m, nd_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_Nd_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([reff, reff_m, reff_m2], legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_Reff_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(org,org_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_ORG_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(org,org_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_ORG_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(so4, so4_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_SO4_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(so4, so4_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_SO4_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn2,ccn2_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn2,ccn2_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn5,ccn5_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN5_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn5,ccn5_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN5_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(cpc3,ncn3_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN3nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(cpc3,ncn3_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN3nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(cpc10,ncn10_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(cpc10,ncn10_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(pcasp100, ncn100_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN100_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(pcasp100, ncn100_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN100_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(lwc, lwc_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_LWC_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(lwc, lwc_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_LWC_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(nd, nd_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(nd, nd_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(reff, reff_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(reff, reff_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

#%% joint histogram
fig,ax = plot.jointhist([pcasp100, ncn100_m, ncn100_m2],[ccn2.data, ccn2_m.data, ccn2_m2.data],
                    xedges=np.arange(0,500,20), yedges=np.arange(0,500,20), vmax=0.3,
                    normalize_x=True, #xlimit=(0,1000), ylimit=(0,1000),
                    xlabel=['CN (>100nm) (cm$^{-3}$)','CN (>100nm) (cm$^{-3}$)','CN (>100nm) (cm$^{-3}$)'], 
                    ylabel='0.2%CCN (cm$^{-3}$)', title=['Obs','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_CN_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ccn2.data, ccn2_m.data, ccn2_m2.data], [nd.data, nd_m.data, nd_m2.data],
                    xedges=np.arange(0,1200,100), yedges=np.arange(0,500,50), 
                    normalize_x=True, #xlimit=(0,1000), ylimit=(0,1000),
                    xlabel=['0.2%CCN (cm$^{-3}$)','0.2%CCN (cm$^{-3}$)','0.2%CCN (cm$^{-3}$)'], 
                    ylabel='Nd (cm$^{-3}$)', title=['Obs','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_Nd_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist( [nd.data, nd_m.data, nd_m2.data],[lwc.data*1000, lwc_m.data*1000, lwc_m2.data*1000], 
                    xedges=np.arange(0,500,50), yedges=np.arange(0,500,50), 
                    normalize_x=True, #xlimit=(0,1000), ylimit=(0,1000),
                    xlabel=['Nd (cm$^{-3}$)','Nd (cm$^{-3}$)','Nd (cm$^{-3}$)'], 
                    ylabel='LWC (mg/m$^{3}$)', title=['Obs','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_LWC_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([reff, reff_m, reff_m2],[lwc, lwc_m, lwc_m2],
                    xedges=np.arange(0,30,2), yedges=np.arange(0,0.6,0.05), vmax=0.3,
                    normalize_x=True, #xlimit=(0,1000), ylimit=(0,1000),
                    xlabel=['Reff ($\mu$m)','Reff ($\mu$m)','Reff ($\mu$m)'], 
                    ylabel='LWC (g/m3)', title=['Obs','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_LWC_Reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter
fig,ax = plot.scatter([ccn2.data, ccn2_m.data, ccn2_m2.data], [nd.data, nd_m.data, nd_m2.data],
                      title=['Obs','E3SMv1','E3SMv2'], xlimit=(0,800), ylimit=(0,500),
                    xlabel='0.2%CCN (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', 
                    linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_CCN_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.scatter([nd.data, nd_m.data, nd_m2.data], [lwc.data*1000, lwc_m.data*1000, lwc_m2.data*1000], 
                      title=['Obs','E3SMv1','E3SMv2'], xlimit=(0,500), ylimit=(0,500),
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWC (mg/m$^{3}$)', 
                    linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_LWC_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.scatter([nd.data, nd_m.data, nd_m2.data], [reff.data, reff_m.data, reff_m2.data], 
                      title=['Obs','E3SMv1','E3SMv2'], xlimit=(0,500), ylimit=(0,50),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', 
                    linear_fit=False, intercept=True)
fig.savefig(figpath+'scatter_Reff_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
