"""
script to generate all plots for ACEENA aircraft data

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
site = 'ACEENA'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/flight/'

# path of output figures
figpath= '../figures/'+site+'/flight/'

height_bin = np.arange(100,4300,300)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = glob.glob(prep_obs_path + 'CPC_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst)
time_cpc = obsdata['time'].load()
cpc3 = obsdata['cpc3'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()

lst = glob.glob(prep_obs_path + 'CCN_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn1 = obsdata['CCN1'].load()
ccn3 = obsdata['CCN3'].load()
obsdata.close()

lst = glob.glob(prep_obs_path + 'AMS_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst)
time_ams = obsdata['time'].load()
org = obsdata['ORG'].load()
so4 = obsdata['SO4'].load()
no3 = obsdata['NO3'].load()
nh4 = obsdata['NH4'].load()
chl = obsdata['CHL'].load()
obsdata.close()

lst = glob.glob(prep_obs_path + 'WCM_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst)
time_wcm = obsdata['time'].load()
LWC = obsdata['LWC'].load()
obsdata.close()

lst = glob.glob(prep_obs_path + 'PCASP100_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst)
time_pcasp = obsdata['time'].load()
pcasp100 = obsdata['pcasp100'].load()
obsdata.close()

lst = sorted(glob.glob(prep_obs_path + 'beasd_ACEENA_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_beasd = obsdata['time'].load()
size_beasd = obsdata['size'].load()
sizeh_beasd = obsdata['size_high'].load()
sizel_beasd = obsdata['size_low'].load()
beasd = obsdata['size_distribution_merged'].load()
obsdata.close()
dlogDp_beasd = np.log10(sizeh_beasd/sizel_beasd)
beasd = beasd/dlogDp_beasd

lst = glob.glob(prep_obs_path + 'merged_bin_cpc_fims_pcasp_opc_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst)
time_mergeCN = obsdata['time'].load()
size_mergeCN = obsdata['size'].load()
sizeh = obsdata['size_high'].load()
sizel = obsdata['size_low'].load()
mergeCN = obsdata['size_distribution_merged'].load()
obsdata.close()
dlogDp_mergeCN = np.log10(sizeh/sizel)
mergeCN = mergeCN/dlogDp_mergeCN

lst = glob.glob(prep_obs_path + 'mergedSD_ACEENA_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_nd = obsdata['time'].load()
size_nd = obsdata['size'].load()
sizeh = obsdata['size_high'].load()
sizel = obsdata['size_low'].load()
nd = obsdata['Nd'].load()
nd_size = obsdata['Nd_bin'].load()
height = obsdata['height'].load()
obsdata.close()
dlogDp_nd = np.log10(sizeh/sizel)
nd_size = nd_size/dlogDp_nd

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# trim for the same time period
# IOP = 'IOP1'
# time1 = np.datetime64('2017-06-21')
# time2 = np.datetime64('2017-07-21')
IOP = 'IOP2'
time1 = np.datetime64('2018-01-19')
time2 = np.datetime64('2018-02-20')

height = height[np.logical_and(time_nd>=time1, time_nd<=time2)]
org = org[np.logical_and(time_ams>=time1, time_ams<=time2)]
so4 = so4[np.logical_and(time_ams>=time1, time_ams<=time2)]
nh4 = nh4[np.logical_and(time_ams>=time1, time_ams<=time2)]
no3 = no3[np.logical_and(time_ams>=time1, time_ams<=time2)]
chl = chl[np.logical_and(time_ams>=time1, time_ams<=time2)]
ccn1 = ccn1[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
ccn3 = ccn3[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
cpc3 = cpc3[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
cpc10 = cpc10[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
pcasp100 = pcasp100[np.logical_and(time_pcasp>=time1, time_pcasp<=time2)]
LWC = LWC[np.logical_and(time_wcm>=time1, time_wcm<=time2)]
nd = nd[np.logical_and(time_nd>=time1, time_nd<=time2)]

nd_size = nd_size[np.logical_and(time_nd>=time1, time_nd<=time2), :]
mergeCN = mergeCN[np.logical_and(time_mergeCN>=time1, time_mergeCN<=time2), :]
beasd = beasd[np.logical_and(time_beasd>=time1, time_beasd<=time2), :]

# remove Nd retrieval less than 1 cm-3
nd[nd<1000] = np.nan

# change time to standard time
ccn1 = xr.DataArray(data=np.interp(nd.time,ccn1.time, ccn1), coords=dict(time=nd.time))
ccn3 = xr.DataArray(data=np.interp(nd.time,ccn3.time, ccn3), coords=dict(time=nd.time))
pcasp100 = xr.DataArray(data=np.interp(nd.time, pcasp100.time, pcasp100), coords=dict(time=nd.time))
LWC = xr.DataArray(data=np.interp(nd.time,LWC.time, LWC), coords=dict(time=nd.time))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not os.path.exists(figpath):
    os.makedirs(figpath)

#%% 1d histogram
w1 = np.ones_like(org)/sum(~np.isnan(org.data))
w2 = np.ones_like(so4)/sum(~np.isnan(so4.data))
w3 = np.ones_like(nh4)/sum(~np.isnan(nh4.data))
w4 = np.ones_like(no3)/sum(~np.isnan(no3.data))
w5 = np.ones_like(chl)/sum(~np.isnan(chl.data))
fig,ax = plot.hist([org,so4], weights=[w1,w2], bins=np.arange(0,1.2,0.1),
                          legend = ['organic','SO4'], color=['limegreen', 'red'],
                          ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_org_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.hist([nh4,no3,chl], weights=[w3,w4,w5], bins=np.arange(0,0.14,0.01),
                          legend = ['NH4','NO3','Chl'], color=['b', 'orange', 'magenta'],
                          ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_NH4_NO3_Chl_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w2 = np.ones_like(ccn1)/sum(~np.isnan(ccn1.data))
w3 = np.ones_like(ccn3)/sum(~np.isnan(ccn3.data))
fig,ax = plot.hist([ccn1,ccn3], weights=[w2,w3], legend = ['0.1%CCN','0.3%CCN'], 
                    bins=np.arange(0,410,20), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CCN_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w1 = np.ones_like(cpc3)/sum(~np.isnan(cpc3.data))
# w2 = np.ones_like(cpc10)/sum(~np.isnan(cpc10.data))
# fig,ax = plot.hist([cpc3,cpc10], weights=[w1,w2], bins=np.arange(0,2100,100), legend = ['CPC(>3nm)','CPC(>10nm)'], 
#                         ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath+'hist_CPC_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w1 = np.ones_like(pcasp100)/sum(~np.isnan(pcasp100.data))
# fig,ax = plot.hist([pcasp100], weights=[w1], legend = ['PCASP100'], 
#                     bins=np.arange(0,400,20), ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath+'hist_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w1 = np.ones_like(LWC)/sum(~np.isnan(LWC.data))
# fig,ax = plot.hist([LWC], weights=[w1], legend = ['LWC'], bins=np.arange(0,0.6,0.03), 
#                     ylabel='Fraction', xlabel=LWC.units)
# ax.set_yscale('log')
# ax.set_ylim(1e-3,1)
# fig.savefig(figpath+'hist_LWC_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
# w1 = np.ones_like(nd)/sum(~np.isnan(nd.data))
# fig,ax = plot.hist([nd/1000], weights=[w1], legend = ['Nd'], bins=np.arange(1,510,20), 
#                     ylabel='Fraction', xlabel='cm$^{-3}$')
# # ax.set_yscale('log')
# # ax.set_ylim(1e-4,.1)
# fig.savefig(figpath+'hist_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% percentiles with height
# fig,ax = plot.percentile_z([org,so4], [height,height], 
#                       height_bin, figsize=(3,8), 
#                       xlabel='$\mu$g/m$^3$', ylabel='height (m)', legend=['ORG','SO4'])
# fig.savefig(figpath+'percentile_z_org_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.percentile_z([no3,nh4,chl], [height,height,height,], 
#                       height_bin, figsize=(3,8), color=['g','c','r'],
#                       xlabel='$\mu$g/m$^3$', ylabel='height (m)', legend=['NO3','NH4','Chl'])
# fig.savefig(figpath+'percentile_z_no3_nh4_chl_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.percentile_z([ccn1,ccn3], [height,height], 
#                       height_bin, figsize=(3,8), 
#                       xlabel='cm$^{-3}$', ylabel='height (m)', legend=['0.1%CCN','0.3%CCN'])
# fig.savefig(figpath+'percentile_z_CCN_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.percentile_z([cpc3,cpc10,pcasp100], [height,height,height], 
#                       height_bin, figsize=(3,8), 
#                       xlabel='cm$^{-3}$', ylabel='height (m)', legend=['CPC3','CPC10','PCASP100'])
# fig.savefig(figpath+'percentile_z_CN_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.percentile_z([LWC], [height], 
#                       height_bin, figsize=(3,8), 
#                       xlabel='g/m$^{3}$', ylabel='height (m)', legend=['LWC'])
# fig.savefig(figpath+'percentile_z_LWC_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.percentile_z([nd/1000], [height,], 
#                       height_bin, figsize=(3,8), 
#                       xlabel='cm$^{-3}$', ylabel='height (m)', legend=['Nd'])
# fig.savefig(figpath+'percentile_z_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% mean size distribution
fig,ax = plot.mean_size([size_mergeCN, size_beasd], [np.nanmean(mergeCN,axis=0),np.nanmean(beasd,axis=0)], 
                  legend = ['mergeCN','BEASD'],color=['k','b'],
                  xlimit=(5, 2e4), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                    title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
# fig.savefig(figpath+'mean_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.mean_size([size_nd], [np.nanmean(nd_size/1000,axis=0)], 
                  legend = ['Nd'],color=['r'], xlimit=(1, 1e3), ylimit=(1e-5,1e2), 
                  xlabel='Diameter (um)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                    title = 'Mean Cloud Droplet Size Distribution '+site+' '+IOP)
# fig.savefig(figpath+'mean_cloud_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


#%% scatter
fig,ax,fit_line = plot.scatter(pcasp100.data, ccn1.data, xlimit=(0,500), ylimit=(0,500),
                    xlabel='PCASP100 (cm$^{-3}$)', ylabel='0.1%CCN (cm$^{-3}$)', 
                    linear_fit=True, intercept=True, title=None)
fig.savefig(figpath+'scatter_ccn1_PCASP100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax,fit_line = plot.scatter(pcasp100.data, ccn3.data, xlimit=(0,500), ylimit=(0,500),
                    xlabel='PCASP100 (cm$^{-3}$)', ylabel='0.3%CCN (cm$^{-3}$)', 
                    linear_fit=True, intercept=True, title=None)
fig.savefig(figpath+'scatter_ccn3_PCASP100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax,fit_line = plot.scatter(cpc10.data, ccn3.data, xlimit=(0,5000), ylimit=(0,5000),
#                     xlabel='CPC10 (cm$^{-3}$)', ylabel='0.3%CCN (cm$^{-3}$)', 
#                     linear_fit=True, intercept=False, title=None)
# fig.savefig(figpath+'scatter_ccn3_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax,fit_line = plot.scatter(ccn3.data, nd.data/1000, xlimit=(0,300), ylimit=(0,300),
                    xlabel='0.3%CCN (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', 
                    linear_fit=True, intercept=False, title=None)
fig.savefig(figpath+'scatter_Nd_CCN3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% joint histogram
# fig,ax = plot.jointhist([ccn1.data, ccn3.data], [nd.data/1000, nd.data/1000], 
#                     xedges=np.arange(0,300,10), yedges=np.arange(0,300,10), 
#                     normalize_x=False, xlimit=(0,200), ylimit=(0,200),
#                     xlabel=['0.1%CCN (cm$^{-3}$)','0.3%CCN (cm$^{-3}$)'], ylabel='Nd (cm$^{-3}$)', title=None)
# fig.savefig(figpath+'jointhist_Nd_CCN_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statistics
# calc.mean_std_percentiles(org,outfile=figpath+'statistics_1var_ORG_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(so4,outfile=figpath+'statistics_1var_SO4_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(no3,outfile=figpath+'statistics_1var_NO3_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(nh4,outfile=figpath+'statistics_1var_NH4_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(chl,outfile=figpath+'statistics_1var_Chl_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(ccn1,outfile=figpath+'statistics_1var_CCN1_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(ccn3,outfile=figpath+'statistics_1var_CCN3_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(cpc3,outfile=figpath+'statistics_1var_CPC3_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(cpc10,outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(pcasp100,outfile=figpath+'statistics_1var_PCASP100_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(LWC,outfile=figpath+'statistics_1var_LWC_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles(nd/1000,outfile=figpath+'statistics_1var_Nd_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn1, ccn3,outfile=figpath+'statistics_2vars_CCN1vsCCN3_'+site+'_'+IOP+'.txt',
#                         label1='0.1%CCN',label2='0.3%CCN')
