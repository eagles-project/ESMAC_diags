"""
script to generate all plots for CSET aircraft data

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
site = 'CSET'

prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/flight/'

height_bin = np.arange(100,4300,300)

figpath= '../figures/'+site+'/flight/'
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = sorted(glob.glob(prep_obs_path + 'CN_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
uhsas100 = obsdata['UHSAS100'].load()
cpc10 = obsdata['CPC10'].load()
obsdata.close()

lst = sorted(glob.glob(prep_obs_path + 'LWC_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
lwc = obsdata['LWC'].load()
obsdata.close()

lst = sorted(glob.glob(prep_obs_path + 'Nd_size_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
size_1dc = obsdata['size_1dc'].load()
size_2dc = obsdata['size_2dc'].load()
size_cdp = obsdata['size_cdp'].load()
nd_bin_1dc = obsdata['Nd_1dc'].load()
nd_bin_2dc = obsdata['Nd_2dc'].load()
nd_bin_cdp = obsdata['Nd_cdp'].load()
nd_all = obsdata['Nd'].load()
obsdata.close()

lst = sorted(glob.glob(prep_obs_path + 'UHSASsize_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
size_u = obsdata['size'].load()
sizeh_u = obsdata['size_high'].load()
sizel_u = obsdata['size_low'].load()
uhsas = obsdata['size_distribution_uhsas'].load()
height = obsdata['height'].load()
obsdata.close()

# read in E3SM data
lst = sorted(glob.glob(prep_model_path + 'E3SMv1_'+site+'_flight_*.nc'))
modeldata = xr.open_mfdataset(lst)
time_m = modeldata['time'].load()
ncn10_m = modeldata['NCN10'].load()
ncn100_m = modeldata['NCN100'].load()
ncn_m = modeldata['NCNall'].load()
ccn1_m = modeldata['CCN3'].load()
ccn2_m = modeldata['CCN4'].load()
ccn5_m = modeldata['CCN5'].load()
lwc_m = modeldata['LWC'].load()
reff_m = modeldata['REL'].load()
nd_m = modeldata['ICWNC'].load()
nd_bin_m = modeldata['Nd_bin'].load()
modeldata.close()

lst = sorted(glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc'))
modeldata = xr.open_mfdataset(lst)
time_m2 = modeldata['time'].load()
ncn10_m2 = modeldata['NCN10'].load()
ncn100_m2 = modeldata['NCN100'].load()
ncn_m2 = modeldata['NCNall'].load()
ccn1_m2 = modeldata['CCN3'].load()
ccn2_m2 = modeldata['CCN4'].load()
ccn5_m2 = modeldata['CCN5'].load()
lwc_m2 = modeldata['LWC'].load()
reff_m2 = modeldata['REL'].load()
nd_m2 = modeldata['ICWNC'].load()
nd_bin_m2 = modeldata['Nd_bin'].load()
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# unit change
nd_all = nd_all/1000 # #/L to #/cm3
nd_bin_cdp = nd_bin_cdp/1000 # #/L to #/cm3
nd_bin_1dc = nd_bin_1dc/1000 # #/L to #/cm3
nd_bin_2dc = nd_bin_2dc/1000 # #/L to #/cm3

# calculate effective radius from observed size distribution
reff = (np.sum(nd_bin_cdp * (size_cdp**3), axis=1) + np.sum(nd_bin_1dc * (size_1dc**3), axis=1) + np.sum(nd_bin_2dc * (size_2dc**3), axis=1)) / \
        (np.sum(nd_bin_cdp * (size_cdp**2), axis=1) + np.sum(nd_bin_1dc * (size_1dc**2), axis=1) + np.sum(nd_bin_2dc * (size_2dc**2), axis=1))
reff[reff>50] = np.nan
    
# remove Nd less than 10 cm-3
nd_all[nd_all<10] = np.nan
nd_bin_cdp[nd_all<10,:] = np.nan
nd_bin_2dc[nd_all<10,:] = np.nan
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
dlogDp_u = np.log10(sizeh_u/sizel_u)
uhsas = uhsas/dlogDp_u

size_cdp_low = np.hstack((size_cdp[0:1]*0.5, size_cdp[0:-1]))
dlogDp_cdp = np.log10(size_cdp/size_cdp_low)
nd_bin_cdp = nd_bin_cdp/dlogDp_cdp
size_1dc_low = np.hstack((size_1dc[0:1]*0.5, size_1dc[0:-1]))
dlogDp_1dc = np.log10(size_1dc/size_1dc_low)
nd_bin_1dc = nd_bin_1dc/dlogDp_1dc
size_2dc_low = np.hstack((size_2dc[0:1]*0.5, size_2dc[0:-1]))
dlogDp_2dc = np.log10(size_2dc/size_2dc_low)
nd_bin_2dc = nd_bin_2dc/dlogDp_2dc
dlogDp_nd_e3sm = np.log10(np.arange(2,1001)/np.arange(1,1000))
nd_bin_m = nd_bin_m.T/dlogDp_nd_e3sm
nd_bin_m2 = nd_bin_m2.T/dlogDp_nd_e3sm

# combine different instruments
instruments = 'CDP-2DC'
size_nd = np.hstack((size_cdp,  size_2dc))
nd_bin = np.hstack((nd_bin_cdp,  nd_bin_2dc))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% mean size distribution
fig,ax = plot.mean_size_witherror([size_u, np.arange(1,3001), np.arange(1,3001)], 
                        [uhsas,ncn_m,ncn_m2], 
                  legend = ['UHSAS','E3SMv1', 'E3SMv2'],color=['k','r','b'],
                  marker=['.',None,None], linestyles=['none','-','-'],
                  xlimit=(10, 3e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                    title = 'Mean Aerosol Size Distribution '+site)
fig.savefig(figpath+'mean_aerosol_size_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.mean_size_witherror([size_nd, np.arange(1,1000), np.arange(1,1000)], 
                        [nd_bin,nd_bin_m,nd_bin_m2], 
                        marker=['.',None,None], linestyles=['none','-','-'],
                        # marker=['.','.','.'], linestyles=['none','none','none'],
                  legend = [instruments,'E3SMv1', 'E3SMv2'],color=['k','r','b'], xlimit=(2, 1e3), ylimit=(1e-5,1e3), 
                  xlabel='Diameter ($\mu$m)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                    title = 'Mean Cloud Droplet Size Distribution '+site)
fig.savefig(figpath+'mean_cloud_size_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% histogram
w0 = np.ones_like(cpc10)/sum(~np.isnan(cpc10.data))
w1 = np.ones_like(ncn10_m)/sum(~np.isnan(ncn10_m.data))
w2 = np.ones_like(ncn10_m2)/sum(~np.isnan(ncn10_m2.data))
fig,ax = plot.hist([cpc10,ncn10_m,ncn10_m2],  weights=[w0,w1,w2], 
                    legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,1500,50), 
                    title = 'CN (>10nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_CN10_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(uhsas100)/sum(~np.isnan(uhsas100.data))
w1 = np.ones_like(ncn100_m)/sum(~np.isnan(ncn100_m.data))
w2 = np.ones_like(ncn100_m2)/sum(~np.isnan(ncn100_m2.data))
fig,ax = plot.hist([uhsas100,ncn100_m,ncn100_m2],  weights=[w0,w1,w2], 
                    legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,410,20), 
                    title = 'CN (>100nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_CN100_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(lwc)/sum(~np.isnan(lwc.data))
w1 = np.ones_like(lwc_m)/sum(~np.isnan(lwc_m.data))
w2 = np.ones_like(lwc_m2)/sum(~np.isnan(lwc_m2.data))
fig,ax = plot.hist([lwc,lwc_m,lwc_m2],  weights=[w0,w1,w2], 
                    legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,0.6,0.03), 
                    title = 'LWC '+site, ylabel='Fraction', xlabel='g/m$^3$')
fig.savefig(figpath + 'histogram_LWC_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(nd_all)/sum(~np.isnan(nd_all.data))
w1 = np.ones_like(nd_m)/sum(~np.isnan(nd_m.data))
w2 = np.ones_like(nd_m2)/sum(~np.isnan(nd_m2.data))
fig,ax = plot.hist([nd_all,nd_m,nd_m2],  weights=[w0,w1,w2], 
                    legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(10,250,10), 
                    title = 'Nd '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_Nd_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(reff)/sum(~np.isnan(reff.data))
w1 = np.ones_like(reff_m)/sum(~np.isnan(reff_m.data))
w2 = np.ones_like(reff_m2)/sum(~np.isnan(reff_m2.data))
fig,ax = plot.hist([reff,reff_m,reff_m2],  weights=[w0,w1,w2], 
                    legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(2,36,2), 
                    title = 'Reff '+site, ylabel='Fraction', xlabel='$\mu$m')
fig.savefig(figpath + 'histogram_Reff_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% percentiles with height
fig,ax = plot.percentile_z([cpc10,ncn10_m,ncn10_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CN (>10nm)',color=['k','r','b'], xlimit=(-100,3000),
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CN10_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.percentile_z([uhsas100, ncn100_m, ncn100_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='CN (>100nm)',color=['k','r','b'],
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_CN100_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([lwc, lwc_m, lwc_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='LWC',color=['k','r','b'], xlimit=(-.1,1),
                      xlabel='g/m$^{3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_lwc_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([nd_all, nd_m, nd_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='Nd',color=['k','r','b'], xlimit=(-20,300),
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentile_z([reff, reff_m, reff_m2], [height,height,height], 
                      height_bin, figsize=(3,8), title='Reff',color=['k','r','b'],
                      xlabel='$\mu$m', ylabel='height (m)', legend = ['Obs', 'E3SMv1', 'E3SMv2'], )
fig.savefig(figpath+'percentile_z_Reff_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statisticscalc.mean_std_percentiles([ccn5,ccn5_m,ccn5_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_ccn5_'+site+'.txt')
# calc.mean_std_percentiles([cpc10,ncn10_m,ncn10_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_CPC10_'+site+'.txt')
# calc.mean_std_percentiles([uhsas100, ncn100_m, ncn100_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_CN100_'+site+'.txt')
# calc.mean_std_percentiles([lwc, lwc_m, lwc_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_lwc_'+site+'.txt')
# calc.mean_std_percentiles([nd_all, nd_m, nd_m2],legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_Nd_'+site+'.txt')
# calc.mean_std_percentiles([reff, reff_m, reff_m2], legend=['Obs','E3SMv1','E3SMv2'], outfile=figpath+'statistics_1var_Reff_'+site+'.txt')

# calc.bias_corrcoef_RMSE(cpc10,ncn10_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cpc10,ncn10_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(uhsas100, ncn100_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN100_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(uhsas100, ncn100_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN100_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(lwc, lwc_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_LWC_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwc, lwc_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_LWC_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(nd_all, nd_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_all, nd_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(reff, reff_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff, reff_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsOBS_'+site+'.txt')

#%% joint histogram

fig,ax = plot.jointhist([uhsas100, ncn100_m, ncn100_m2], [nd_all, nd_m, nd_m2], 
                    xedges=np.arange(0,400,30), yedges=np.arange(0,400,30), 
                    normalize_x=True, #xlimit=(0,1000), ylimit=(0,1000),
                    xlabel=['CN (>100nm) (cm$^{-3}$)','CN (>100nm) (cm$^{-3}$)','CN (>100nm) (cm$^{-3}$)'], 
                    ylabel='Nd (cm$^{-3}$)', title=['Obs','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_CN_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([reff, reff_m, reff_m2],[lwc, lwc_m, lwc_m2],
                    xedges=np.arange(0,36,3), yedges=np.arange(0,0.9,0.1), vmax=0.3,
                    normalize_x=True, #xlimit=(0,1000), ylimit=(0,1000),
                    xlabel=['Reff ($\mu$m)','Reff ($\mu$m)','Reff ($\mu$m)'], 
                    ylabel='LWC (g/m3)', title=['Obs','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_LWC_Reff_'+site+'_'+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter
fig,ax = plot.scatter([nd_all.data, nd_m.data, nd_m2.data], [reff.data, reff_m.data, reff_m2.data], 
                      title=['Obs','E3SMv1','E3SMv2'], xlimit=(0,300), ylimit=(0,50),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', 
                    linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_Reff_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)