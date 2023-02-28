import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
from CBcolors import CB_color_cycle

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

lst = sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn2 = obsdata['ccn2_fit'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_sfc_path + 'totcld_'+site+'*.nc')
time_cld = obsdata['time']
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_visst = obsdata['tot_cld_visst'].load()
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

lst = sorted(glob.glob(prep_sfc_path + 'Ndrop_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ndrop = obsdata['time'].load()
ndrop = obsdata['cdnc'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'reff_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_reff = obsdata['time'].load()
reff = obsdata['reff'].load()
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

obsdata = xr.open_mfdataset(prep_sfc_path + 'precip_'+site+'*.nc')
time_pr = obsdata['time']
precip = obsdata['precip'].load()
obsdata.close()

# satellite data
satdata = xr.open_mfdataset(prep_sat_path + 'cod_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
cod_sat = satdata['cod'].load()
satdata.close()
satdata = xr.open_mfdataset(prep_sat_path + 'LWP_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
lwp_sat = satdata['lwp'].load()
satdata.close()
lst = sorted(glob.glob(prep_sat_path + 'Reff_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
reff_sat = obsdata['reff'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'Nd_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
nd_sat = obsdata['Nd'].load()
obsdata.close()
satdata = xr.open_mfdataset(prep_sat_path + 'IWP_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
iwp_sat = satdata['iwp'].load()
satdata.close()
satdata = xr.open_mfdataset(prep_sat_path + 'solarzenith_VISSTgrid_'+site+'*.nc')
time_sat = satdata['time']
zenith_sat = satdata['solar_zenith_angle'].load()
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

# E3SM data
filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
ccn2_m = modeldata['CCN4'].load()
ncn10_m = modeldata['NCN10'].load()
ncn100_m = modeldata['NCN100'].load()
cod_m = modeldata['cod'].load()
reff_m = modeldata['reff'].load()
lwp_m = modeldata['TGCLDLWP'].load()
nd_m = modeldata['Nd_mean'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
iwp_m = modeldata['TGCLDIWP'].load()
cth_m = modeldata['cth'].load()
modeldata.close()
filename = prep_model_path + 'E3SMv2_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
time_mp = modeldata['time'].load()
LTS700_m = modeldata['LTS700'].load()
LTS850_m = modeldata['LTS850'].load()
thetadiff_cb_m = modeldata['thetadiff_cb'].load()
modeldata.close()

precip_m = precip_m * 3600 * 1000   # m/s to mm/hr
lwp[lwp<20] = np.nan
lwp_sat[lwp_sat<20] = np.nan
lwp_m[lwp_m<20] = np.nan
cod[cod<2] = np.nan
cod_sat[cod_sat<2] = np.nan
cod_m[cod_m<2] = np.nan
cod[cod>91] = np.nan
cod_sat[cod_sat>91] = np.nan
cod_m[cod_m>91] = np.nan

ndrop[ndrop<10] = np.nan
nd_sat[nd_sat<10] = np.nan
nd_m[nd_m<10] = np.nan
ndrop[ndrop>800] = np.nan
nd_sat[nd_sat>800] = np.nan
nd_m[nd_m>800] = np.nan

ccn2.data[ccn2.data<20] = np.nan
ccn2_m.data[ccn2_m.data<20] = np.nan

#%%

# select low-level ovarcasting clouds
# select surface coupling samples
idx_sfc = np.logical_and(np.logical_and(cth<4000, cld_arscl>90), thetadiff_cb<2)
idx_sat = np.logical_and(np.logical_and(cth_sat<4, cfall_sat>90), thetadiff_cb<2)
idx_m = np.logical_and(np.logical_and(cth_m<4000, cld_m>90), thetadiff_cb_m<2)

idx_sfc[iwp_sat>10] = False
idx_sat[iwp_sat>10] = False

idx_sat[zenith_sat.data>70] = False

#%% 1d histogram
w1 = np.ones_like(ccn2)/sum(~np.isnan(ccn2.data))
w2 = np.ones_like(ccn2_m)/sum(~np.isnan(ccn2_m.data))
fig,ax = plot.hist([ccn2,ccn2_m], weights=[w1,w2], bins=np.arange(0,1500,50),
                    legend =['CCN Counter', 'E3SMv2',], color=[CB_color_cycle[1],'k'],
                    xlabel = '(a) Surface CCN (SS=0.2%) (cm$^{-3}$)', ylabel='PDF')
                    # title = '(a) Surface CCN (SS=0.2%) '+site, ylabel='PDF', xlabel='cm$^{-3}$')

w0 = np.ones_like(ndrop[idx_sfc])/sum(~np.isnan(ndrop[idx_sfc].data))
w000 = np.ones_like(nd_sat[idx_sat])/sum(~np.isnan(nd_sat[idx_sat].data))
w1 = np.ones_like(nd_m[idx_m])/sum(~np.isnan(nd_m[idx_m].data))
fig,ax = plot.hist([ndrop[idx_sfc],nd_sat[idx_sat],nd_m[idx_m]],  weights=[w0,w000,w1], 
                   legend = ['Ndrop','Satellite', 'E3SMv2'], color=[CB_color_cycle[1],CB_color_cycle[0],'k'],bins=np.arange(0,410,20), 
                   xlabel = '(c) Nd (cm$^{-3}$)', ylabel='PDF')
                    # title = '(c) Nd '+site, ylabel='PDF', xlabel='cm$^{-3}$')

w0 = np.ones_like(reff[idx_sfc])/sum(~np.isnan(reff[idx_sfc].data))
w000 = np.ones_like(reff_sat[idx_sat])/sum(~np.isnan(reff_sat[idx_sat].data))
w1 = np.ones_like(reff_m[idx_m])/sum(~np.isnan(reff_m[idx_m].data))
fig,ax = plot.hist([reff[idx_sfc],reff_sat[idx_sat],reff_m[idx_m]], weights=[w0,w000,w1], 
                    legend = ['MFRSR','Satellite', 'E3SMv2'], color=[CB_color_cycle[1],CB_color_cycle[0],'k'], bins=np.arange(3,22,1), 
                    xlabel = '(e) Cloud Effective Radius ($\mu$m)', ylabel='PDF')
                    # title = '(e) Cloud Effective Radius '+site,ylabel='PDF', xlabel='$\mu$m')

w0 = np.ones_like(cod[idx_sfc])/sum(~np.isnan(cod[idx_sfc].data))
w00 = np.ones_like(cod_sat[idx_sat])/sum(~np.isnan(cod_sat[idx_sat].data))
w1 = np.ones_like(cod_m[idx_m])/sum(~np.isnan(cod_m[idx_m].data))
fig,ax = plot.hist( [cod[idx_sfc], cod_sat[idx_sat], cod_m[idx_m], ], weights=[w0,w00,w1,], 
                    legend = ['MFRSR','Satellite','E3SMv2',], color=[CB_color_cycle[1],CB_color_cycle[0],'k',],bins=np.arange(2,81,3), 
                    xlabel = '(g) Cloud Optical Depth (unitless)', ylabel='PDF')
                    # title='(g) Cloud Optical Depth '+site, ylabel='PDF', xlabel='unitless')
# fig.savefig(figpath+'hist_cod_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cld_arscl)/sum(~np.isnan(cld_arscl.data))
w00 = np.ones_like(cld_visst)/sum(~np.isnan(cld_visst.data))
w1 = np.ones_like(cld_m)/sum(~np.isnan(cld_m.data))
fig,ax = plot.hist([cld_arscl,cld_visst,cld_m], 
                    weights=[w0,w00,w1],  bins=np.arange(0,101,5), 
                    legend = ['ARMBE','Satellite','E3SMv2'], color=[CB_color_cycle[1],CB_color_cycle[0],'k'],
                    xlabel = '(i) Cloud Fraction (%)', ylabel='PDF')
                      # title = '(i) Cloud PDF '+site, ylabel='PDF', xlabel="%")
# fig.savefig(figpath+'hist_totcld_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

