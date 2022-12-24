
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
import esmac_diags.plotting.calc_statistics as calc

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'MARCUS'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/ship/'
prep_sat_path = '../prep_data/'+site+'/satellite/'

# set output path for plots
figpath= '../figures/'+site+'/sfc_toa/'
if not os.path.exists(figpath):
    os.makedirs(figpath)

filename = prep_sfc_path + 'T_RH_Ps_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
T_marcus = obsdata['T'].load()
RH_marcus = obsdata['RH'].load()
Ps_marcus = obsdata['Ps'].load()
lon_marcus = obsdata['lon']
lat_marcus = obsdata['lat']
obsdata.close()
filename = prep_sfc_path + 'CCN_'+site+'_exhaustfree.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
ccn2_marcus = obsdata['CCN2'].load()
lon_marcus = obsdata['lon']
lat_marcus = obsdata['lat']
obsdata.close()
filename = prep_sfc_path + 'CN_'+site+'_exhaustfree.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
cpc10_marcus = obsdata['CPC10'].load()
uhsas100_marcus = obsdata['UHSAS100'].load()
obsdata.close()
filename = prep_sfc_path + 'CNsize_UHSAS_'+site+'_exhaustfree.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
CNsize_marcus = obsdata['size_distribution_uhsas'].load()
dmean_marcus = obsdata['size'].load()
dmin_marcus = obsdata['size_low'].load()
dmax_marcus = obsdata['size_high'].load()
obsdata.close()
filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
lwp_marcus = obsdata['lwp'].load()
obsdata.close()
# filename = prep_sfc_path + 'sfc_radiation_'+site+'.nc'
# obsdata = xr.open_dataset(filename)
# time_rad = obsdata['time'].load()
# lwdnsfc = obsdata['lwdn'].load()
# swdnsfc = obsdata['swdn'].load()
# lwupsfc = obsdata['lwup'].load()
# swupsfc = obsdata['swup'].load()
# obsdata.close()
# lwdnsfc_aceena = lwdnsfc.sel(time=time_aceena)
# swdnsfc_aceena = swdnsfc.sel(time=time_aceena)
# lwupsfc_aceena = lwupsfc.sel(time=time_aceena)
# swupsfc_aceena = swupsfc.sel(time=time_aceena)

filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
cld_marcus = obsdata['cldfrac'].load()
obsdata.close()

filename = prep_sfc_path + 'totcld_sensitivity_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
cld_5 = obsdata['cldfrac_5'].load()
cld_10 = obsdata['cldfrac_10'].load()
cld_20 = obsdata['cldfrac_20'].load()
cld_30 = obsdata['cldfrac_30'].load()
obsdata.close()

# satellite data
filename = prep_sat_path + 'cod_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cod_sat_marcus = obsdata['cod'].load()
obsdata.close()
filename = prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
reff_sat_marcus = obsdata['reff'].load()
obsdata.close()
filename = prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwp_sat_marcus = obsdata['lwp'].load()
obsdata.close()
filename = prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
nd_sat_marcus = obsdata['Nd'].load()
obsdata.close()
filename = prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwnettoa_marcus = obsdata['lwnettoa'].load()
obsdata.close()
filename = prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
swnettoa_marcus = obsdata['swnettoa'].load()
obsdata.close()
filename = prep_sat_path + 'cloudfraction_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cld_sat_marcus = obsdata['cldtot'].load()
obsdata.close()
filename = prep_sat_path + 'albedo_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
albedo_marcus = obsdata['albedo'].load()
obsdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv1_'+site+'_ship.nc'
modeldata = xr.open_dataset(filename)
time_m_marcus = modeldata['time'].load()
lat_m_marcus = modeldata['lat'].load()
lon_m_marcus = modeldata['lon'].load()
T_m_marcus = modeldata['TREFHT'].load()
RH_m_marcus = modeldata['RELHUM'].load()
Ps_m_marcus = modeldata['PS'].load()
ccn2_m_marcus = modeldata['CCN4'].load()
ncn10_m_marcus = modeldata['NCN10'].load()
ncn100_m_marcus = modeldata['NCN100'].load()
CNsize_m_marcus = modeldata['NCNall'].load()
cod_m_marcus = modeldata['cod'].load()
reff_m_marcus = modeldata['reff'].load()
lwp_m_marcus = modeldata['TGCLDLWP'].load()
nd_m_marcus = modeldata['Nd_mean'].load()
precip_m_marcus = modeldata['PRECT'].load()
cld_m_marcus = modeldata['CLDTOT'].load()
lwdnsfc_m_marcus = modeldata['FLDS'].load()
lwnetsfc_m_marcus = modeldata['FLNS'].load()
lwnettoa_m_marcus = modeldata['FLNT'].load()
lwuptoa_m_marcus = modeldata['FLUT'].load()
swdnsfc_m_marcus = modeldata['FSDS'].load()
swnetsfc_m_marcus = modeldata['FSNS'].load()
swdntoa_m_marcus = modeldata['SOLIN'].load()
swnettoa_m_marcus = modeldata['FSNT'].load()
swuptoa_m_marcus = modeldata['FSUTOA'].load()
modeldata.close()
lwupsfc_m_marcus = lwnetsfc_m_marcus + lwdnsfc_m_marcus
swupsfc_m_marcus = swdnsfc_m_marcus - swnetsfc_m_marcus
albedo_m_marcus = swuptoa_m_marcus/swdntoa_m_marcus*100

filename = prep_model_path + 'E3SMv2_'+site+'_ship.nc'
modeldata = xr.open_dataset(filename)
time_m2_marcus = modeldata['time'].load()
lat_m2_marcus = modeldata['lat'].load()
lon_m2_marcus = modeldata['lon'].load()
T_m2_marcus = modeldata['TREFHT'].load()
RH_m2_marcus = modeldata['RELHUM'].load()
Ps_m2_marcus = modeldata['PS'].load()
ccn2_m2_marcus = modeldata['CCN4'].load()
ncn10_m2_marcus = modeldata['NCN10'].load()
ncn100_m2_marcus = modeldata['NCN100'].load()
CNsize_m2_marcus = modeldata['NCNall'].load()
cod_m2_marcus = modeldata['cod'].load()
reff_m2_marcus = modeldata['reff'].load()
lwp_m2_marcus = modeldata['TGCLDLWP'].load()
nd_m2_marcus = modeldata['Nd_mean'].load()
precip_m2_marcus = modeldata['PRECT'].load()
cld_m2_marcus = modeldata['CLDTOT'].load()
lwdnsfc_m2_marcus = modeldata['FLDS'].load()
lwnetsfc_m2_marcus = modeldata['FLNS'].load()
lwnettoa_m2_marcus = modeldata['FLNT'].load()
lwuptoa_m2_marcus = modeldata['FLUT'].load()
swdnsfc_m2_marcus = modeldata['FSDS'].load()
swnetsfc_m2_marcus = modeldata['FSNS'].load()
swdntoa_m2_marcus = modeldata['SOLIN'].load()
swnettoa_m2_marcus = modeldata['FSNT'].load()
swuptoa_m2_marcus = modeldata['FSUTOA'].load()
modeldata.close()
lwupsfc_m2_marcus = lwnetsfc_m2_marcus + lwdnsfc_m2_marcus
swupsfc_m2_marcus = swdnsfc_m2_marcus - swnetsfc_m2_marcus
albedo_m2_marcus = swuptoa_m2_marcus/swdntoa_m2_marcus*100

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

nd_sat_marcus[nd_sat_marcus<10] = np.nan
nd_m_marcus[nd_m_marcus<10] = np.nan
nd_m2_marcus[nd_m2_marcus<10] = np.nan
nd_sat_marcus[nd_sat_marcus>500] = np.nan
nd_m_marcus[nd_m_marcus>500] = np.nan
nd_m2_marcus[nd_m2_marcus>500] = np.nan

lwp_sat_marcus[lwp_sat_marcus<20] = np.nan
lwp_marcus[lwp_marcus<20] = np.nan
lwp_m_marcus[lwp_m_marcus<20] = np.nan
lwp_m2_marcus[lwp_m2_marcus<20] = np.nan

cod_sat_marcus[cod_sat_marcus<2] = np.nan
cod_m_marcus[cod_m_marcus<2] = np.nan
cod_m2_marcus[cod_m2_marcus<2] = np.nan
cod_sat_marcus[cod_sat_marcus>100] = np.nan
cod_m_marcus[cod_m_marcus>100] = np.nan
cod_m2_marcus[cod_m2_marcus>100] = np.nan

# remove MARCUS data near Antarctic
T_m_marcus[lat_m_marcus<-65] = np.nan
RH_m_marcus[lat_m_marcus<-65] = np.nan
Ps_m_marcus[lat_m_marcus<-65] = np.nan
ccn2_m_marcus[lat_m_marcus<-65] = np.nan
ncn10_m_marcus[lat_m_marcus<-65] = np.nan
ncn100_m_marcus[lat_m_marcus<-65] = np.nan
CNsize_m_marcus[:,lat_m_marcus<-65] = np.nan
cod_m_marcus[lat_m_marcus<-65] = np.nan
reff_m_marcus[lat_m_marcus<-65] = np.nan
lwp_m_marcus[lat_m_marcus<-65] = np.nan
nd_m_marcus[lat_m_marcus<-65] = np.nan
precip_m_marcus[lat_m_marcus<-65] = np.nan
cld_m_marcus[lat_m_marcus<-65] = np.nan
albedo_m_marcus[lat_m_marcus<-65] = np.nan
T_m2_marcus[lat_m2_marcus<-65] = np.nan
RH_m2_marcus[lat_m2_marcus<-65] = np.nan
Ps_m2_marcus[lat_m2_marcus<-65] = np.nan
ccn2_m2_marcus[lat_m2_marcus<-65] = np.nan
ncn10_m2_marcus[lat_m2_marcus<-65] = np.nan
ncn100_m2_marcus[lat_m2_marcus<-65] = np.nan
CNsize_m2_marcus[:,lat_m2_marcus<-65] = np.nan
cod_m2_marcus[lat_m2_marcus<-65] = np.nan
reff_m2_marcus[lat_m2_marcus<-65] = np.nan
lwp_m2_marcus[lat_m2_marcus<-65] = np.nan
nd_m2_marcus[lat_m2_marcus<-65] = np.nan
precip_m2_marcus[lat_m2_marcus<-65] = np.nan
cld_m2_marcus[lat_m2_marcus<-65] = np.nan
albedo_m2_marcus[lat_m2_marcus<-65] = np.nan

# unit change
T_m_marcus = T_m_marcus - 273.15
T_m2_marcus = T_m2_marcus - 273.15
Ps_m_marcus = Ps_m_marcus * 0.01
Ps_m2_marcus = Ps_m2_marcus * 0.01

#%% change to dN/dlnDp
dlogDp_u=np.empty(len(dmean_marcus))
for bb in range(len(dmean_marcus)):
    dlogDp_u[bb]=np.log10(dmax_marcus[bb]/dmin_marcus[bb])
dlogDp=np.empty(3000)
for bb in range(3000):
    dlogDp[bb]=np.log10((bb+2)/(bb+1))
CNsize_m_marcus = CNsize_m_marcus.T    
CNsize_m2_marcus = CNsize_m2_marcus.T    
dN_dlogDp_o = np.ones_like(CNsize_marcus)
dN_dlogDp_m = np.ones_like(CNsize_m_marcus)
dN_dlogDp_m2 = np.ones_like(CNsize_m2_marcus)
for tt in range(len(time_marcus)):
    dN_dlogDp_o[tt,:] = CNsize_marcus[tt,:]/dlogDp_u
    dN_dlogDp_m[tt,:] = CNsize_m_marcus[tt,:]/dlogDp
    dN_dlogDp_m2[tt,:] = CNsize_m2_marcus[tt,:]/dlogDp
pdf_obs = np.nanmean(dN_dlogDp_o,axis=0)
pdf_m = np.nanmean(dN_dlogDp_m,axis=0)
pdf_m2 = np.nanmean(dN_dlogDp_m2,axis=0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# #%% timeseries
# if site=='MAGIC':
#     legnum = np.arange(2,20)
#     startdate = [np.datetime64('2012-09-22'), np.datetime64('2012-10-06'), np.datetime64('2012-10-20'), np.datetime64('2012-11-03'), 
#                   np.datetime64('2012-11-17'), np.datetime64('2012-12-01'), np.datetime64('2012-12-15'), np.datetime64('2012-12-29'), 
#                   np.datetime64('2013-05-11'), np.datetime64('2013-05-25'), np.datetime64('2013-06-08'), np.datetime64('2013-06-22'), 
#                   np.datetime64('2013-07-07'), np.datetime64('2013-07-20'), np.datetime64('2013-08-03'), np.datetime64('2013-08-17'), 
#                   np.datetime64('2013-08-31'), np.datetime64('2013-09-14'), ]
#     enddate = [np.datetime64('2012-10-04'), np.datetime64('2012-10-18'), np.datetime64('2012-11-01'), np.datetime64('2012-11-15'), 
#                 np.datetime64('2012-11-30'), np.datetime64('2012-12-13'), np.datetime64('2012-12-27'), np.datetime64('2013-01-06'), 
#                 np.datetime64('2013-05-23'), np.datetime64('2013-06-06'), np.datetime64('2013-06-20'), np.datetime64('2013-07-03'), 
#                 np.datetime64('2013-07-18'), np.datetime64('2013-08-01'), np.datetime64('2013-08-15'), np.datetime64('2013-08-29'), 
#                 np.datetime64('2013-09-12'), np.datetime64('2013-09-26'), ]
# elif site=='MARCUS':
#     legnum = np.arange(1,5)
#     startdate = [np.datetime64('2017-10-30'), np.datetime64('2017-12-13'), np.datetime64('2018-01-16'), np.datetime64('2018-03-09')]
#     enddate = [np.datetime64('2017-12-02'), np.datetime64('2018-01-11'), np.datetime64('2018-03-04'), np.datetime64('2018-03-22')]
        
# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[T_marcus, T_m_marcus, T_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='Temperature (C) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   color=['k','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_T_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[RH_marcus, RH_m_marcus, RH_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='RH (%) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   color=['k','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_RH_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[Ps_marcus, Ps_m_marcus, Ps_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='Ps (hPa) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   color=['k','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_Ps_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[cpc10_marcus, ncn10_m_marcus, ncn10_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='CN (>10nm) (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   color=['k','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_CN10_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[uhsas100_marcus,ncn100_m_marcus,ncn100_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='CN (>100nm) (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   color=['k','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_CN100_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[ccn2_marcus,ccn2_m_marcus,ccn2_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='CCN (SS=0.2%) (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   color=['k','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_CCN2_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus, time_marcus],[lwp_marcus,lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='LWP (g/m$^2$) '+site, legend=['Ship','Satellite','E3SMv1','E3SMv2'], 
#                   color=['k','gray','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_LWP_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus, time_marcus],[cld_marcus,cld_sat_marcus, cld_m_marcus,cld_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='Cloud Fraction (%) '+site, legend=['Ship','Satellite','E3SMv1','E3SMv2'], 
#                   color=['k','gray','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_CF_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[nd_sat_marcus, nd_m_marcus,nd_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='Nd (cm$^{-3}$) '+site, legend=['Satellite','E3SMv1','E3SMv2'], 
#                   color=['gray','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_Nd_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[reff_sat_marcus,reff_m_marcus,reff_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='Reff ($\mu$m) '+site, legend=['Satellite','E3SMv1','E3SMv2'], 
#                   color=['gray','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_Reff_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.timeseries([time_marcus, time_marcus, time_marcus],[cod_sat_marcus,cod_m_marcus,cod_m2_marcus], figsize=(10,3),
#                   xlabel='Time', ylabel=None, title='Cloud Optical Depth (N/A) '+site, legend=['Satellite','E3SMv1','E3SMv2'], 
#                   color=['gray','r','b'])
# for ll in range(len(legnum)):
#     ax.set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.86, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_COD_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# ax.set_xlim(time_marcus[0],time_marcus[-1])

# #%% aerosol size distribution, timeseries and mean
# fig,ax = plot.timeseries_size([time_marcus, time_marcus, time_marcus],[dmean_marcus, np.arange(1,3001), np.arange(1,3001)], 
#                           [dN_dlogDp_o.T,dN_dlogDp_m.T,dN_dlogDp_m2.T], figsize=(10,6),
#                   xlabel='Time', ylabel=None, title='dN/dlog10Dp (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'])
# for ll in range(len(legnum)):
#     for ii in range(len(ax)):
#         ax[ii].set_xlim(startdate[ll],enddate[ll])
#     txt = fig.text(0.1, 0.9, 'trip # '+format(legnum[ll]))
#     figname = figpath + 'timeseries_CNsize_'+site+'_leg'+format(legnum[ll])+'.png'
#     fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     txt.remove()
# for ii in range(len(ax)):
#     ax[ii].set_xlim(time_marcus[0],time_marcus[-1])

# fig,ax = plot.mean_size([dmean_marcus, np.arange(1,3001), np.arange(1,3001)], [pdf_obs, pdf_m, pdf_m2],
#                         figsize=(7,5), xlimit=None, ylimit=None, xscale='log',yscale='log',
#                   xlabel='Size (nm)', ylabel=None, title='dN/dlog10Dp (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   linestyles=['none','-','-'], marker=['.',None,None], color=['k','r','b'])
# figname = figpath + 'CNsize_mean_'+site+'.png'
# fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% histogram
# w0 = np.ones_like(ccn2_marcus)/sum(~np.isnan(ccn2_marcus.data))
# w1 = np.ones_like(ccn2_m_marcus)/sum(~np.isnan(ccn2_m_marcus.data))
# w2 = np.ones_like(ccn2_m2_marcus)/sum(~np.isnan(ccn2_m2_marcus.data))
# fig,ax = plot.hist([ccn2_marcus,ccn2_m_marcus,ccn2_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,410,20), 
#                     title = 'CCN (SS=0.2%) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath + 'histogram_CCN2_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(cpc10_marcus)/sum(~np.isnan(cpc10_marcus.data))
# w1 = np.ones_like(ncn10_m_marcus)/sum(~np.isnan(ncn10_m_marcus.data))
# w2 = np.ones_like(ncn10_m2_marcus)/sum(~np.isnan(ncn10_m2_marcus.data))
# fig,ax = plot.hist([cpc10_marcus,ncn10_m_marcus,ncn10_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,1500,50), 
#                     title = 'CN (>10nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath + 'histogram_CN10_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(uhsas100_marcus)/sum(~np.isnan(uhsas100_marcus.data))
# w1 = np.ones_like(ncn100_m_marcus)/sum(~np.isnan(ncn100_m_marcus.data))
# w2 = np.ones_like(ncn100_m2_marcus)/sum(~np.isnan(ncn100_m2_marcus.data))
# fig,ax = plot.hist([uhsas100_marcus,ncn100_m_marcus,ncn100_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,410,20), 
#                     title = 'CN (>100nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath + 'histogram_CN100_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(lwp_marcus)/sum(~np.isnan(lwp_marcus.data))
# w00 = np.ones_like(lwp_sat_marcus)/sum(~np.isnan(lwp_sat_marcus.data))
# w1 = np.ones_like(lwp_m_marcus)/sum(~np.isnan(lwp_m_marcus.data))
# w2 = np.ones_like(lwp_m2_marcus)/sum(~np.isnan(lwp_m2_marcus.data))
# fig,ax = plot.hist([lwp_marcus,lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus],  weights=[w0,w00,w1,w2], 
#                     legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], bins=np.arange(10,510,25), 
#                     title = 'LWP '+site, ylabel='Fraction', xlabel='g/m$^2$')
# fig.savefig(figpath + 'histogram_LWP_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(cld_marcus)/sum(~np.isnan(cld_marcus.data))
# w00 = np.ones_like(cld_sat_marcus)/sum(~np.isnan(cld_sat_marcus.data))
# w1 = np.ones_like(cld_m_marcus)/sum(~np.isnan(cld_m_marcus.data))
# w2 = np.ones_like(cld_m2_marcus)/sum(~np.isnan(cld_m2_marcus.data))
# fig,ax = plot.hist([cld_marcus,cld_sat_marcus, cld_m_marcus,cld_m2_marcus],  weights=[w0,w00,w1,w2], 
#                     legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],  bins=np.arange(0,101,5), 
#                     title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')
# fig.savefig(figpath + 'histogram_CF_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w00 = np.ones_like(nd_sat_marcus)/sum(~np.isnan(nd_sat_marcus.data))
# w1 = np.ones_like(nd_m_marcus)/sum(~np.isnan(nd_m_marcus.data))
# w2 = np.ones_like(nd_m2_marcus)/sum(~np.isnan(nd_m2_marcus.data))
# fig,ax = plot.hist([ nd_sat_marcus, nd_m_marcus,nd_m2_marcus],  weights=[w00,w1,w2], 
#                     legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'], bins=np.arange(10,160,5), 
#                     title = 'Nd '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath + 'histogram_Nd_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w00 = np.ones_like(reff_sat_marcus)/sum(~np.isnan(reff_sat_marcus.data))
# w1 = np.ones_like(reff_m_marcus)/sum(~np.isnan(reff_m_marcus.data))
# w2 = np.ones_like(reff_m2_marcus)/sum(~np.isnan(reff_m2_marcus.data))
# fig,ax = plot.hist([ reff_sat_marcus,reff_m_marcus,reff_m2_marcus], weights=[w00,w1,w2], 
#                     legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'], bins=np.arange(4,31,1), 
#                     title = 'Reff '+site, ylabel='Fraction', xlabel='$\mu$m')
# fig.savefig(figpath + 'histogram_Reff_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(cod_sat_marcus)/sum(~np.isnan(cod_sat_marcus.data))
# w1 = np.ones_like(cod_m_marcus)/sum(~np.isnan(cod_m_marcus.data))
# w2 = np.ones_like(cod_m2_marcus)/sum(~np.isnan(cod_m2_marcus.data))
# fig,ax = plot.hist([cod_sat_marcus,cod_m_marcus,cod_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'], bins=np.arange(0,61,3), 
#                     title = 'Cloud Optical Depth '+site, ylabel='Fraction', xlabel='N/A')
# fig.savefig(figpath + 'histogram_COD_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% sensitivity to different LWP threshold
# w0 = np.ones_like(cld_5)/sum(~np.isnan(cld_5.data))
# w1 = np.ones_like(cld_m_marcus)/sum(~np.isnan(cld_m_marcus.data))
# w2 = np.ones_like(cld_m2_marcus)/sum(~np.isnan(cld_m2_marcus.data))
# fig,ax = plot.hist([cld_5,cld_m_marcus,cld_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship (thres_lwp=5)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
#                     title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

# w0 = np.ones_like(cld_10)/sum(~np.isnan(cld_10.data))
# w1 = np.ones_like(cld_m_marcus)/sum(~np.isnan(cld_m_marcus.data))
# w2 = np.ones_like(cld_m2_marcus)/sum(~np.isnan(cld_m2_marcus.data))
# fig,ax = plot.hist([cld_10,cld_m_marcus,cld_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship (thres_lwp=10)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
#                     title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

# w0 = np.ones_like(cld_20)/sum(~np.isnan(cld_20.data))
# w1 = np.ones_like(cld_m_marcus)/sum(~np.isnan(cld_m_marcus.data))
# w2 = np.ones_like(cld_m2_marcus)/sum(~np.isnan(cld_m2_marcus.data))
# fig,ax = plot.hist([cld_20,cld_m_marcus,cld_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship (thres_lwp=20)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
#                     title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

# w0 = np.ones_like(cld_30)/sum(~np.isnan(cld_30.data))
# w1 = np.ones_like(cld_m_marcus)/sum(~np.isnan(cld_m_marcus.data))
# w2 = np.ones_like(cld_m2_marcus)/sum(~np.isnan(cld_m2_marcus.data))
# fig,ax = plot.hist([cld_30,cld_m_marcus,cld_m2_marcus],  weights=[w0,w1,w2], 
#                     legend = ['Ship (thres_lwp=30)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
#                     title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

# #%% precentile with latitude
# latbin = np.arange(-65.5,-42,1)

# plot.percentile_lat([T_marcus,T_m_marcus,T_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, #ylimit=(0,3000), 
#                   xlabel='Latitude', ylabel=None, title = 'Temperature (C) '+site,
#                   legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
# fig.savefig(figpath + 'percentile_lat_T_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([RH_marcus,RH_m_marcus,RH_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, #ylimit=(0,3000), 
#                   xlabel='Latitude', ylabel=None, title = 'RH (%) '+site,
#                   legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
# fig.savefig(figpath + 'percentile_lat_RH_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([Ps_marcus,Ps_m_marcus,Ps_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, #ylimit=(900,1030), 
#                   xlabel='Latitude', ylabel=None, title = 'Ps (hPa) '+site,
#                   legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
# fig.savefig(figpath + 'percentile_lat_Ps_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([cpc10_marcus,ncn10_m_marcus,ncn10_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=(0,2000), 
#                   xlabel='Latitude', ylabel=None, title = 'CN (>10nm) (cm$^{-3}$) '+site,
#                   legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
# fig.savefig(figpath + 'percentile_lat_CN10_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([ccn2_marcus,ccn2_m_marcus,ccn2_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=(0,400), 
#                   xlabel='Latitude', ylabel=None, title = 'CCN (SS=0.2%) (cm$^{-3}$) '+site,
#                   legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
# fig.savefig(figpath + 'percentile_lat_CCN2_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([uhsas100_marcus,ncn100_m_marcus,ncn100_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=None, 
#                   xlabel='Latitude', ylabel=None, title = 'CN (>100nm) (cm$^{-3}$) '+site,
#                   legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
# fig.savefig(figpath + 'percentile_lat_CN100_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([lwp_marcus,lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus], [lat_marcus,lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=None, 
#                   xlabel='Latitude', ylabel=None, title = 'LWP (g/m$^2$) '+site,
#                   legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'])
# fig.savefig(figpath + 'percentile_lat_LWP_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([cld_marcus,cld_sat_marcus,cld_m_marcus,cld_m2_marcus], [lat_marcus,lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=None, 
#                   xlabel='Latitude', ylabel=None, title = 'Cloud Fraction (%) '+site,
#                   legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'])
# fig.savefig(figpath + 'percentile_lat_CF_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([nd_sat_marcus,nd_m_marcus,nd_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=(0,200), 
#                   xlabel='Latitude', ylabel=None, title = 'Nd (cm$^{-3}$) '+site,
#                   legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'])
# fig.savefig(figpath + 'percentile_lat_Nd_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([reff_sat_marcus,reff_m_marcus,reff_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=None, 
#                   xlabel='Latitude', ylabel=None, title = 'Cloud Effective Radius ($\mu$m) '+site,
#                   legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'])
# fig.savefig(figpath + 'percentile_lat_Reff_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# plot.percentile_lat([cod_sat_marcus,cod_m_marcus,cod_m2_marcus], [lat_marcus,lat_marcus,lat_marcus], latbin, 
#                     figsize=(8,2), xlimit=None, ylimit=None, 
#                   xlabel='Latitude', ylabel=None, title = 'Cloud Optical Depth (N/A) '+site,
#                   legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'])
# fig.savefig(figpath + 'percentile_lat_COD_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% mean statistics

# calc.bias_corrcoef_RMSE(T_m_marcus, T_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_T_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(RH_m_marcus, RH_marcus,label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_RH_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(Ps_m_marcus,Ps_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_Ps_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn10_m_marcus,cpc10_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CN10_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn100_m_marcus,uhsas100_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CN100_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ccn2_m_marcus,ccn2_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m_marcus,lwp_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_LWP_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m_marcus,cld_marcus, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CF_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cod_m_marcus,cod_sat_marcus, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_COD_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m_marcus,lwp_sat_marcus, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_LWP_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m_marcus,cld_sat_marcus, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_CF_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_m_marcus,nd_sat_marcus, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_Nd_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_m_marcus,reff_sat_marcus, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_Reff_E3SMv1vsSat_'+site+'.txt')

# calc.bias_corrcoef_RMSE(T_m2_marcus, T_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_T_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(RH_m2_marcus, RH_marcus,label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_RH_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(Ps_m2_marcus,Ps_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_Ps_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn10_m2_marcus,cpc10_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CN10_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn100_m2_marcus,uhsas100_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CN100_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ccn2_m2_marcus,ccn2_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m2_marcus,lwp_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_LWP_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m2_marcus,cld_marcus, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CF_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cod_m2_marcus,cod_sat_marcus, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_COD_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m2_marcus,lwp_sat_marcus, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_LWP_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m2_marcus,cld_sat_marcus, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_CF_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_m2_marcus,nd_sat_marcus, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_Nd_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_m2_marcus,reff_sat_marcus, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_Reff_E3SMv2vsSat_'+site+'.txt')

# calc.mean_std_percentiles([T_marcus,T_m_marcus,T_m2_marcus],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_T_'+site+'.txt')
# calc.mean_std_percentiles([RH_marcus,RH_m_marcus,RH_m2_marcus],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_RH_'+site+'.txt')
# calc.mean_std_percentiles([Ps_marcus,Ps_m_marcus,Ps_m2_marcus],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Ps_'+site+'.txt')
# calc.mean_std_percentiles([cpc10_marcus,ncn10_m_marcus,ncn10_m2_marcus],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CN10_'+site+'.txt')
# calc.mean_std_percentiles([uhsas100_marcus,ncn100_m_marcus,ncn100_m2_marcus],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CN100_'+site+'.txt')
# calc.mean_std_percentiles([ccn2_marcus,ccn2_m_marcus,ccn2_m2_marcus],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CCN2_'+site+'.txt')
# calc.mean_std_percentiles([lwp_marcus,lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus],['Ship','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_LWP_'+site+'.txt')
# calc.mean_std_percentiles([cld_marcus,cld_sat_marcus,cld_m_marcus,cld_m2_marcus],['Ship','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CF_'+site+'.txt')
# calc.mean_std_percentiles([nd_sat_marcus,nd_m_marcus,nd_m2_marcus],['Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Nd_'+site+'.txt')
# calc.mean_std_percentiles([reff_sat_marcus,reff_m_marcus,reff_m2_marcus],['Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Reff_'+site+'.txt')
# calc.mean_std_percentiles([cod_sat_marcus,cod_m_marcus,cod_m2_marcus],['Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_COD_'+site+'.txt')


#%% joint histogram
fig,ax = plot.jointhist([uhsas100_marcus,ncn100_m_marcus,ncn100_m2_marcus], [ccn2_marcus,ccn2_m_marcus,ccn2_m2_marcus], 
                    xedges=np.arange(0,400,20),yedges=np.arange(0,400,20), normalize_x=True,
                    xlabel='CN (>100nm) (cm$^{-3}$)', ylabel='CCN (SS=0.2%) (cm$^{-3}$)', vmax=0.5,
                    title=['Ship','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CN100_CCN2_ship_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ccn2_marcus,ccn2_m_marcus,ccn2_m2_marcus], [nd_sat_marcus,nd_m_marcus,nd_m2_marcus],
                    xedges=np.arange(0,350,20),yedges=np.arange(0,300,15), normalize_x=True,
                    xlabel='CCN (SS=0.2%) (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', vmax=0.4,
                    title=['Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CCN2_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([nd_sat_marcus,nd_m_marcus,nd_m2_marcus],[lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus], 
                    xedges=np.arange(0,300,15),yedges=np.arange(0,300,10), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', vmax=0.4,
                    title=['Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_LWP_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([nd_sat_marcus,nd_m_marcus,nd_m2_marcus],[reff_sat_marcus,reff_m_marcus,reff_m2_marcus],
                    xedges=np.arange(0,300,15),yedges=np.arange(4,25,1), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', vmax=0.25,
                    title=['Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_Reff_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([cod_sat_marcus,cod_sat_marcus,cod_m_marcus,cod_m2_marcus],[lwp_marcus,lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus], 
                    xedges=np.arange(0,40,2),yedges=np.arange(0,300,10), normalize_x=True,
                    xlabel='Cloud Optical Depth (N/A)', ylabel='LWP (g/m$^2$)', vmax=0.25,
                    title=['Ship','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_COD_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter plot

fig,ax = plot.scatter([nd_sat_marcus.data,nd_m_marcus.data,nd_m2_marcus.data], 
                      [ccn2_marcus.data,ccn2_m_marcus.data,ccn2_m2_marcus.data],
                      xlimit=(0,300), ylimit=(0,350),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', title=['Satellite','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=True)
# fig.savefig(figpath+'scatter_Nd_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.scatter([uhsas100_marcus.data,ncn100_m_marcus.data,ncn100_m2_marcus.data], 
                      [ccn2_marcus.data,ccn2_m_marcus.data,ccn2_m2_marcus.data],
                      xlimit=(0,500), ylimit=(0,500),
                    xlabel='Surface CN (>100nm) (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', title=['Ship','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=True)
# fig.savefig(figpath+'scatter_CN100_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% heatmaps

# xedges=np.exp(np.arange(np.log(10),6.5,0.5))
# yedges=np.exp(np.arange(np.log(10),6.5,0.5))
fig,ax = plot.heatmap([nd_sat_marcus.data,nd_m_marcus.data,nd_m2_marcus.data],
                      [lwp_sat_marcus,lwp_m_marcus,lwp_m2_marcus],
                      [albedo_marcus,albedo_m_marcus,albedo_m2_marcus],vmax=60,
                    xedges=np.arange(0,300,20), yedges=np.arange(10,300,20),
                    # xedges=xedges, yedges=yedges, 
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', zlabel='TOA Albedo (%)',
                    title=['Satellite','E3SMv1','E3SMv2'])
# fig.savefig(figpath+'heatmap_CCN2_LWP_Nd_sfc_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


