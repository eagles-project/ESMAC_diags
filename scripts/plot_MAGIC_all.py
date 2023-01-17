
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
site = 'MAGIC'
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
time_magic = obsdata['time'].load()
T_magic = obsdata['T'].load()
RH_magic = obsdata['RH'].load()
Ps_magic = obsdata['Ps'].load()
lon_magic = obsdata['lon']
lat_magic = obsdata['lat']
obsdata.close()
filename = prep_sfc_path + 'CCN_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
ccn2_magic = obsdata['CCN2'].load()
obsdata.close()
filename = prep_sfc_path + 'CN_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
cpc10_magic = obsdata['CPC10'].load()
uhsas100_magic = obsdata['UHSAS100'].load()
obsdata.close()
filename = prep_sfc_path + 'CNsize_UHSAS_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
CNsize_magic = obsdata['size_distribution_uhsas'].load()
dmean_magic = obsdata['size'].load()
dmin_magic = obsdata['size_low'].load()
dmax_magic = obsdata['size_high'].load()
obsdata.close()
filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
lwp_magic = obsdata['lwp'].load()
obsdata.close()

filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
cld_magic = obsdata['cldfrac'].load()
obsdata.close()

filename = prep_sfc_path + 'totcld_sensitivity_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
cld_5 = obsdata['cldfrac_5'].load()
cld_10 = obsdata['cldfrac_10'].load()
cld_20 = obsdata['cldfrac_20'].load()
cld_30 = obsdata['cldfrac_30'].load()
obsdata.close()

filename = prep_sfc_path + 'Nd_Reff_Wu_etal_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
reff_wu_magic = obsdata['reff'].load()
nd_wu_magic = obsdata['cdnc'].load()
obsdata.close()

# satellite data
filename = prep_sat_path + 'cod_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cod_sat_magic = obsdata['cod'].load()
obsdata.close()
filename = prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
reff_sat_magic = obsdata['reff'].load()
obsdata.close()
filename = prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwp_sat_magic = obsdata['lwp'].load()
obsdata.close()
filename = prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
nd_sat_magic = obsdata['Nd'].load()
obsdata.close()
filename = prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwnettoa_magic = obsdata['lwnettoa'].load()
obsdata.close()
filename = prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
swnettoa_magic = obsdata['swnettoa'].load()
obsdata.close()
filename = prep_sat_path + 'cloudfraction_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cld_sat_magic = obsdata['cldtot'].load()
obsdata.close()
filename = prep_sat_path + 'albedo_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
albedo_magic = obsdata['albedo'].load()
obsdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv1_'+site+'_ship.nc'
modeldata = xr.open_dataset(filename)
time_m_magic = modeldata['time'].load()
lat_m_magic = modeldata['lat'].load()
lon_m_magic = modeldata['lon'].load()
T_m_magic = modeldata['TREFHT'].load()
RH_m_magic = modeldata['RELHUM'].load()
Ps_m_magic = modeldata['PS'].load()
ccn2_m_magic = modeldata['CCN4'].load()
ncn10_m_magic = modeldata['NCN10'].load()
ncn100_m_magic = modeldata['NCN100'].load()
CNsize_m_magic = modeldata['NCNall'].load()
cod_m_magic = modeldata['cod'].load()
reff_m_magic = modeldata['reff'].load()
lwp_m_magic = modeldata['TGCLDLWP'].load()
nd_m_magic = modeldata['Nd_mean'].load()
precip_m_magic = modeldata['PRECT'].load()
cld_m_magic = modeldata['CLDTOT'].load()
lwdnsfc_m_magic = modeldata['FLDS'].load()
lwnetsfc_m_magic = modeldata['FLNS'].load()
lwnettoa_m_magic = modeldata['FLNT'].load()
lwuptoa_m_magic = modeldata['FLUT'].load()
swdnsfc_m_magic = modeldata['FSDS'].load()
swnetsfc_m_magic = modeldata['FSNS'].load()
swdntoa_m_magic = modeldata['SOLIN'].load()
swnettoa_m_magic = modeldata['FSNT'].load()
swuptoa_m_magic = modeldata['FSUTOA'].load()
modeldata.close()
lwupsfc_m_magic = lwnetsfc_m_magic + lwdnsfc_m_magic
swupsfc_m_magic = swdnsfc_m_magic - swnetsfc_m_magic
albedo_m_magic = swuptoa_m_magic/swdntoa_m_magic*100

filename = prep_model_path + 'E3SMv2_'+site+'_ship.nc'
modeldata = xr.open_dataset(filename)
time_m2_magic = modeldata['time'].load()
lat_m2_magic = modeldata['lat'].load()
lon_m2_magic = modeldata['lon'].load()
T_m2_magic = modeldata['TREFHT'].load()
RH_m2_magic = modeldata['RELHUM'].load()
Ps_m2_magic = modeldata['PS'].load()
ccn2_m2_magic = modeldata['CCN4'].load()
ncn10_m2_magic = modeldata['NCN10'].load()
ncn100_m2_magic = modeldata['NCN100'].load()
CNsize_m2_magic = modeldata['NCNall'].load()
cod_m2_magic = modeldata['cod'].load()
reff_m2_magic = modeldata['reff'].load()
lwp_m2_magic = modeldata['TGCLDLWP'].load()
nd_m2_magic = modeldata['Nd_mean'].load()
precip_m2_magic = modeldata['PRECT'].load()
cld_m2_magic = modeldata['CLDTOT'].load()
lwdnsfc_m2_magic = modeldata['FLDS'].load()
lwnetsfc_m2_magic = modeldata['FLNS'].load()
lwnettoa_m2_magic = modeldata['FLNT'].load()
lwuptoa_m2_magic = modeldata['FLUT'].load()
swdnsfc_m2_magic = modeldata['FSDS'].load()
swnetsfc_m2_magic = modeldata['FSNS'].load()
swdntoa_m2_magic = modeldata['SOLIN'].load()
swnettoa_m2_magic = modeldata['FSNT'].load()
swuptoa_m2_magic = modeldata['FSUTOA'].load()
modeldata.close()
lwupsfc_m2_magic = lwnetsfc_m2_magic + lwdnsfc_m2_magic
swupsfc_m2_magic = swdnsfc_m2_magic - swnetsfc_m2_magic
albedo_m2_magic = swuptoa_m2_magic/swdntoa_m2_magic*100

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

nd_wu_magic[nd_wu_magic<10] = np.nan
nd_sat_magic[nd_sat_magic<10] = np.nan
nd_m_magic[nd_m_magic<10] = np.nan
nd_m2_magic[nd_m2_magic<10] = np.nan
nd_wu_magic[nd_wu_magic>500] = np.nan
nd_sat_magic[nd_sat_magic>500] = np.nan
nd_m_magic[nd_m_magic>500] = np.nan
nd_m2_magic[nd_m2_magic>500] = np.nan

lwp_magic[lwp_magic<20] = np.nan
lwp_sat_magic[lwp_sat_magic<20] = np.nan
lwp_m_magic[lwp_m_magic<20] = np.nan
lwp_m2_magic[lwp_m2_magic<20] = np.nan

cod_sat_magic[cod_sat_magic<2] = np.nan
cod_m_magic[cod_m_magic<2] = np.nan
cod_m2_magic[cod_m2_magic<2] = np.nan
cod_sat_magic[cod_sat_magic>100] = np.nan
cod_m_magic[cod_m_magic>100] = np.nan
cod_m2_magic[cod_m2_magic>100] = np.nan

# remove magic data near California coast
T_m_magic[lon_m_magic>240] = np.nan
RH_m_magic[lon_m_magic>240] = np.nan
Ps_m_magic[lon_m_magic>240] = np.nan
ccn2_m_magic[lon_m_magic>240] = np.nan
ncn10_m_magic[lon_m_magic>240] = np.nan
ncn100_m_magic[lon_m_magic>240] = np.nan
CNsize_m_magic[:,lon_m_magic>240] = np.nan
cod_m_magic[lon_m_magic>240] = np.nan
reff_m_magic[lon_m_magic>240] = np.nan
lwp_m_magic[lon_m_magic>240] = np.nan
nd_m_magic[lon_m_magic>240] = np.nan
precip_m_magic[lon_m_magic>240] = np.nan
cld_m_magic[lon_m_magic>240] = np.nan
albedo_m_magic[lon_m_magic>240] = np.nan
T_m2_magic[lon_m2_magic>240] = np.nan
RH_m2_magic[lon_m2_magic>240] = np.nan
Ps_m2_magic[lon_m2_magic>240] = np.nan
ccn2_m2_magic[lon_m2_magic>240] = np.nan
ncn10_m2_magic[lon_m2_magic>240] = np.nan
ncn100_m2_magic[lon_m2_magic>240] = np.nan
CNsize_m2_magic[:,lon_m2_magic>240] = np.nan
cod_m2_magic[lon_m2_magic>240] = np.nan
reff_m2_magic[lon_m2_magic>240] = np.nan
lwp_m2_magic[lon_m2_magic>240] = np.nan
nd_m2_magic[lon_m2_magic>240] = np.nan
precip_m2_magic[lon_m2_magic>240] = np.nan
cld_m2_magic[lon_m2_magic>240] = np.nan
albedo_m2_magic[lon_m2_magic>240] = np.nan

# unit change
T_m_magic = T_m_magic - 273.15
T_m2_magic = T_m2_magic - 273.15
Ps_m_magic = Ps_m_magic * 0.01
Ps_m2_magic = Ps_m2_magic * 0.01

#%% change to dN/dlnDp
dlogDp_u=np.empty(len(dmean_magic))
for bb in range(len(dmean_magic)):
    dlogDp_u[bb]=np.log10(dmax_magic[bb]/dmin_magic[bb])
dlogDp=np.empty(3000)
for bb in range(3000):
    dlogDp[bb]=np.log10((bb+2)/(bb+1))
CNsize_m_magic = CNsize_m_magic.T    
CNsize_m2_magic = CNsize_m2_magic.T    
dN_dlogDp_o = np.ones_like(CNsize_magic)
dN_dlogDp_m = np.ones_like(CNsize_m_magic)
dN_dlogDp_m2 = np.ones_like(CNsize_m2_magic)
for tt in range(len(time_magic)):
    dN_dlogDp_o[tt,:] = CNsize_magic[tt,:]/dlogDp_u
    dN_dlogDp_m[tt,:] = CNsize_m_magic[tt,:]/dlogDp
    dN_dlogDp_m2[tt,:] = CNsize_m2_magic[tt,:]/dlogDp
pdf_obs = np.nanmean(dN_dlogDp_o,axis=0)
pdf_m = np.nanmean(dN_dlogDp_m,axis=0)
pdf_m2 = np.nanmean(dN_dlogDp_m2,axis=0)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%% timeseries
if site=='MAGIC':
    legnum = np.arange(2,20)
    startdate = [np.datetime64('2012-09-22'), np.datetime64('2012-10-06'), np.datetime64('2012-10-20'), np.datetime64('2012-11-03'), 
                  np.datetime64('2012-11-17'), np.datetime64('2012-12-01'), np.datetime64('2012-12-15'), np.datetime64('2012-12-29'), 
                  np.datetime64('2013-05-11'), np.datetime64('2013-05-25'), np.datetime64('2013-06-08'), np.datetime64('2013-06-22'), 
                  np.datetime64('2013-07-07'), np.datetime64('2013-07-20'), np.datetime64('2013-08-03'), np.datetime64('2013-08-17'), 
                  np.datetime64('2013-08-31'), np.datetime64('2013-09-14'), ]
    enddate = [np.datetime64('2012-10-04'), np.datetime64('2012-10-18'), np.datetime64('2012-11-01'), np.datetime64('2012-11-15'), 
                np.datetime64('2012-11-30'), np.datetime64('2012-12-13'), np.datetime64('2012-12-27'), np.datetime64('2013-01-06'), 
                np.datetime64('2013-05-23'), np.datetime64('2013-06-06'), np.datetime64('2013-06-20'), np.datetime64('2013-07-03'), 
                np.datetime64('2013-07-18'), np.datetime64('2013-08-01'), np.datetime64('2013-08-15'), np.datetime64('2013-08-29'), 
                np.datetime64('2013-09-12'), np.datetime64('2013-09-26'), ]
elif site=='MARCUS':
    legnum = np.arange(1,5)
    startdate = [np.datetime64('2017-10-30'), np.datetime64('2017-12-13'), np.datetime64('2018-01-16'), np.datetime64('2018-03-09')]
    enddate = [np.datetime64('2017-12-02'), np.datetime64('2018-01-11'), np.datetime64('2018-03-04'), np.datetime64('2018-03-22')]
        
fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[T_magic, T_m_magic, T_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='Temperature (C) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
                  color=['k','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_T_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[RH_magic, RH_m_magic, RH_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='RH (%) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
                  color=['k','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_RH_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[Ps_magic, Ps_m_magic, Ps_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='Ps (hPa) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
                  color=['k','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_Ps_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[cpc10_magic, ncn10_m_magic, ncn10_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='CN (>10nm) (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
                  color=['k','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_CN10_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[uhsas100_magic,ncn100_m_magic,ncn100_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='CN (>100nm) (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
                  color=['k','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_CN100_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[ccn2_magic,ccn2_m_magic,ccn2_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='CCN (SS=0.2%) (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
                  color=['k','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_CCN2_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic, time_magic],[lwp_magic,lwp_sat_magic,lwp_m_magic,lwp_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='LWP (g/m$^2$) '+site, legend=['Ship','Satellite','E3SMv1','E3SMv2'], 
                  color=['k','gray','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_LWP_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic, time_magic],[cld_magic,cld_sat_magic, cld_m_magic,cld_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='Cloud Fraction (%) '+site, legend=['Ship','Satellite','E3SMv1','E3SMv2'], 
                  color=['k','gray','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_CF_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic, time_magic],[nd_wu_magic, nd_sat_magic, nd_m_magic,nd_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='Nd (cm$^{-3}$) '+site, legend=['Ship','Satellite','E3SMv1','E3SMv2'], 
                  color=['k','gray','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_Nd_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic, time_magic],[reff_wu_magic, reff_sat_magic,reff_m_magic,reff_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='Reff ($\mu$m) '+site, legend=['Ship','Satellite','E3SMv1','E3SMv2'], 
                  color=['k','gray','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_Reff_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

fig,ax = plot.timeseries([time_magic, time_magic, time_magic],[cod_sat_magic,cod_m_magic,cod_m2_magic], figsize=(10,3),
                  xlabel='Time', ylabel=None, title='Cloud Optical Depth (N/A) '+site, legend=['Satellite','E3SMv1','E3SMv2'], 
                  color=['gray','r','b'])
for ll in range(len(legnum)):
    ax.set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.83, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_COD_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
ax.set_xlim(time_magic[0],time_magic[-1])

#%% aerosol size distribution, timeseries and mean
fig,ax = plot.timeseries_size([time_magic, time_magic, time_magic],[dmean_magic, np.arange(1,3001), np.arange(1,3001)], 
                          [dN_dlogDp_o.T,dN_dlogDp_m.T,dN_dlogDp_m2.T], figsize=(10,6),
                  xlabel='Time', ylabel=None, title='dN/dlog10Dp (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'])
for ll in range(len(legnum)):
    for ii in range(len(ax)):
        ax[ii].set_xlim(startdate[ll],enddate[ll])
    txt = fig.text(0.1, 0.9, 'trip # '+format(legnum[ll]))
    figname = figpath + 'timeseries_CNsize_'+site+'_leg'+format(legnum[ll])+'.png'
    fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    txt.remove()
for ii in range(len(ax)):
    ax[ii].set_xlim(time_magic[0],time_magic[-1])

# fig,ax = plot.mean_size([dmean_magic, np.arange(1,3001), np.arange(1,3001)], [pdf_obs, pdf_m, pdf_m2],
#                         figsize=(7,5), xlimit=None, ylimit=None, xscale='log',yscale='log',
#                   xlabel='Size (nm)', ylabel=None, title='dN/dlog10Dp (cm$^{-3}$) '+site, legend=['Ship','E3SMv1','E3SMv2'], 
#                   linestyles=['none','-','-'], marker=['.',None,None], color=['k','r','b'])
fig,ax = plot.mean_size_witherror([dmean_magic, np.arange(1,3001), np.arange(1,3001)], 
                                  [dN_dlogDp_o,dN_dlogDp_m,dN_dlogDp_m2],
                        figsize=(7,5), xlimit=(10, 3e3), ylimit=(1e-2,1e4),
                  xlabel='Size (nm)', ylabel=None, title='dN/dlog10Dp (cm$^{-3}$) '+site, 
                  legend=['Ship','E3SMv1','E3SMv2'], 
                  linestyles=['none','-','-'], marker=['o',None,None], color=['k','r','b'])
figname = figpath + 'CNsize_mean_'+site+'.png'
fig.savefig(figname,dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% histogram
w0 = np.ones_like(ccn2_magic)/sum(~np.isnan(ccn2_magic.data))
w1 = np.ones_like(ccn2_m_magic)/sum(~np.isnan(ccn2_m_magic.data))
w2 = np.ones_like(ccn2_m2_magic)/sum(~np.isnan(ccn2_m2_magic.data))
fig,ax = plot.hist([ccn2_magic,ccn2_m_magic,ccn2_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,610,30), 
                    title = 'CCN (SS=0.2%) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_CCN2_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cpc10_magic)/sum(~np.isnan(cpc10_magic.data))
w1 = np.ones_like(ncn10_m_magic)/sum(~np.isnan(ncn10_m_magic.data))
w2 = np.ones_like(ncn10_m2_magic)/sum(~np.isnan(ncn10_m2_magic.data))
fig,ax = plot.hist([cpc10_magic,ncn10_m_magic,ncn10_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,2100,100), 
                    title = 'CN (>10nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_CN10_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(uhsas100_magic)/sum(~np.isnan(uhsas100_magic.data))
w1 = np.ones_like(ncn100_m_magic)/sum(~np.isnan(ncn100_m_magic.data))
w2 = np.ones_like(ncn100_m2_magic)/sum(~np.isnan(ncn100_m2_magic.data))
fig,ax = plot.hist([uhsas100_magic,ncn100_m_magic,ncn100_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,610,30), 
                    title = 'CN (>100nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_CN100_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(lwp_magic)/sum(~np.isnan(lwp_magic.data))
w00 = np.ones_like(lwp_sat_magic)/sum(~np.isnan(lwp_sat_magic.data))
w1 = np.ones_like(lwp_m_magic)/sum(~np.isnan(lwp_m_magic.data))
w2 = np.ones_like(lwp_m2_magic)/sum(~np.isnan(lwp_m2_magic.data))
fig,ax = plot.hist([lwp_magic,lwp_sat_magic,lwp_m_magic,lwp_m2_magic],  weights=[w0,w00,w1,w2], 
                    legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], bins=np.arange(10,410,20), 
                    title = 'LWP '+site, ylabel='Fraction', xlabel='g/m$^2$')
fig.savefig(figpath + 'histogram_LWP_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cld_magic)/sum(~np.isnan(cld_magic.data))
w00 = np.ones_like(cld_sat_magic)/sum(~np.isnan(cld_sat_magic.data))
w1 = np.ones_like(cld_m_magic)/sum(~np.isnan(cld_m_magic.data))
w2 = np.ones_like(cld_m2_magic)/sum(~np.isnan(cld_m2_magic.data))
fig,ax = plot.hist([cld_magic,cld_sat_magic, cld_m_magic,cld_m2_magic],  weights=[w0,w00,w1,w2], 
                    legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],  bins=np.arange(0,101,5), 
                    title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')
fig.savefig(figpath + 'histogram_CF_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(nd_wu_magic)/sum(~np.isnan(nd_wu_magic.data))
w00 = np.ones_like(nd_sat_magic)/sum(~np.isnan(nd_sat_magic.data))
w1 = np.ones_like(nd_m_magic)/sum(~np.isnan(nd_m_magic.data))
w2 = np.ones_like(nd_m2_magic)/sum(~np.isnan(nd_m2_magic.data))
fig,ax = plot.hist([nd_wu_magic, nd_sat_magic, nd_m_magic,nd_m2_magic],  weights=[w0,w00,w1,w2], 
                    legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], bins=np.arange(10,210,5), 
                    title = 'Nd '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath + 'histogram_Nd_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(reff_wu_magic)/sum(~np.isnan(reff_wu_magic.data))
w00 = np.ones_like(reff_sat_magic)/sum(~np.isnan(reff_sat_magic.data))
w1 = np.ones_like(reff_m_magic)/sum(~np.isnan(reff_m_magic.data))
w2 = np.ones_like(reff_m2_magic)/sum(~np.isnan(reff_m2_magic.data))
fig,ax = plot.hist([reff_wu_magic, reff_sat_magic,reff_m_magic,reff_m2_magic], weights=[w0,w00,w1,w2], 
                    legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], bins=np.arange(4,28,1), 
                    title = 'Reff '+site, ylabel='Fraction', xlabel='$\mu$m')
fig.savefig(figpath + 'histogram_Reff_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cod_sat_magic)/sum(~np.isnan(cod_sat_magic.data))
w1 = np.ones_like(cod_m_magic)/sum(~np.isnan(cod_m_magic.data))
w2 = np.ones_like(cod_m2_magic)/sum(~np.isnan(cod_m2_magic.data))
fig,ax = plot.hist([cod_sat_magic,cod_m_magic,cod_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'], bins=np.arange(0,61,3), 
                    title = 'Cloud Optical Depth '+site, ylabel='Fraction', xlabel='N/A')
fig.savefig(figpath + 'histogram_COD_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% sensitivity to different LWP threshold
w0 = np.ones_like(cld_5)/sum(~np.isnan(cld_5.data))
w1 = np.ones_like(cld_m_magic)/sum(~np.isnan(cld_m_magic.data))
w2 = np.ones_like(cld_m2_magic)/sum(~np.isnan(cld_m2_magic.data))
fig,ax = plot.hist([cld_5,cld_m_magic,cld_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship (thres_lwp=5)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
                    title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

w0 = np.ones_like(cld_10)/sum(~np.isnan(cld_10.data))
w1 = np.ones_like(cld_m_magic)/sum(~np.isnan(cld_m_magic.data))
w2 = np.ones_like(cld_m2_magic)/sum(~np.isnan(cld_m2_magic.data))
fig,ax = plot.hist([cld_10,cld_m_magic,cld_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship (thres_lwp=10)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
                    title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

w0 = np.ones_like(cld_20)/sum(~np.isnan(cld_20.data))
w1 = np.ones_like(cld_m_magic)/sum(~np.isnan(cld_m_magic.data))
w2 = np.ones_like(cld_m2_magic)/sum(~np.isnan(cld_m2_magic.data))
fig,ax = plot.hist([cld_20,cld_m_magic,cld_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship (thres_lwp=20)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
                    title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

w0 = np.ones_like(cld_30)/sum(~np.isnan(cld_30.data))
w1 = np.ones_like(cld_m_magic)/sum(~np.isnan(cld_m_magic.data))
w2 = np.ones_like(cld_m2_magic)/sum(~np.isnan(cld_m2_magic.data))
fig,ax = plot.hist([cld_30,cld_m_magic,cld_m2_magic],  weights=[w0,w1,w2], 
                    legend = ['Ship (thres_lwp=30)','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,101,5), 
                    title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel='%')

#%% precentile with latitude
latbin = np.arange(21.5,34,1)

plot.percentile_lat([T_magic,T_m_magic,T_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, #ylimit=(0,3000), 
                  xlabel='Latitude', ylabel=None, title = 'Temperature (C) '+site,
                  legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
fig.savefig(figpath + 'percentile_lat_T_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([RH_magic,RH_m_magic,RH_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, #ylimit=(0,3000), 
                  xlabel='Latitude', ylabel=None, title = 'RH (%) '+site,
                  legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
fig.savefig(figpath + 'percentile_lat_RH_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([Ps_magic,Ps_m_magic,Ps_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=(1000,1030), 
                  xlabel='Latitude', ylabel=None, title = 'Ps (hPa) '+site,
                  legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
fig.savefig(figpath + 'percentile_lat_Ps_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([cpc10_magic,ncn10_m_magic,ncn10_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=(0,3000), 
                  xlabel='Latitude', ylabel=None, title = 'CN (>10nm) (cm$^{-3}$) '+site,
                  legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
fig.savefig(figpath + 'percentile_lat_CN10_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([ccn2_magic,ccn2_m_magic,ccn2_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'CCN (SS=0.2%) (cm$^{-3}$) '+site,
                  legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
fig.savefig(figpath + 'percentile_lat_CCN2_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([uhsas100_magic,ncn100_m_magic,ncn100_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'CN (>100nm) (cm$^{-3}$) '+site,
                  legend = ['Ship','E3SMv1','E3SMv2'], color=['k','r','b'])
fig.savefig(figpath + 'percentile_lat_CN100_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([lwp_magic,lwp_sat_magic,lwp_m_magic,lwp_m2_magic], [lat_magic,lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'LWP (g/m$^2$) '+site,
                  legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'])
fig.savefig(figpath + 'percentile_lat_LWP_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([cld_magic,cld_sat_magic,cld_m_magic,cld_m2_magic], [lat_magic,lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'Cloud Fraction (%) '+site,
                  legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'])
fig.savefig(figpath + 'percentile_lat_CF_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([nd_wu_magic,nd_sat_magic,nd_m_magic,nd_m2_magic], [lat_magic,lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'Nd (cm$^{-3}$) '+site,
                  legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'])
fig.savefig(figpath + 'percentile_lat_Nd_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([reff_wu_magic,reff_sat_magic,reff_m_magic,reff_m2_magic], [lat_magic,lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'Cloud Effective Radius ($\mu$m) '+site,
                  legend = ['Ship','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'])
fig.savefig(figpath + 'percentile_lat_Reff_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

plot.percentile_lat([cod_sat_magic,cod_m_magic,cod_m2_magic], [lat_magic,lat_magic,lat_magic], latbin, 
                    figsize=(10,3), xlimit=None, ylimit=None, 
                  xlabel='Latitude', ylabel=None, title = 'Cloud Optical Depth (N/A) '+site,
                  legend = ['Satellite','E3SMv1','E3SMv2'], color=['gray','r','b'])
fig.savefig(figpath + 'percentile_lat_COD_'+site+'.png', dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% mean statistics

# calc.bias_corrcoef_RMSE(T_m_magic, T_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_T_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(RH_m_magic, RH_magic,label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_RH_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(Ps_m_magic,Ps_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_Ps_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn10_m_magic,cpc10_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CN10_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn100_m_magic,uhsas100_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CN100_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ccn2_m_magic,ccn2_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m_magic,lwp_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_LWP_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m_magic,cld_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_CF_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_m_magic,nd_wu_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_Nd_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_m_magic,reff_wu_magic, label1='E3SMv1', label2='Ship',
#                         outfile=figpath+'statistics_Reff_E3SMv1vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cod_m_magic,cod_sat_magic, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_COD_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m_magic,lwp_sat_magic, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_LWP_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m_magic,cld_sat_magic, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_CF_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_m_magic,nd_sat_magic, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_Nd_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_m_magic,reff_sat_magic, label1='E3SMv1', label2='Satllite',
#                         outfile=figpath+'statistics_Reff_E3SMv1vsSat_'+site+'.txt')

# calc.bias_corrcoef_RMSE(T_m2_magic, T_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_T_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(RH_m2_magic, RH_magic,label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_RH_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(Ps_m2_magic,Ps_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_Ps_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn10_m2_magic,cpc10_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CN10_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ncn100_m2_magic,uhsas100_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CN100_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ccn2_m2_magic,ccn2_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m2_magic,lwp_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_LWP_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m2_magic,cld_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_CF_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_m2_magic,nd_wu_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_Nd_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_m2_magic,reff_wu_magic, label1='E3SMv2', label2='Ship',
#                         outfile=figpath+'statistics_Reff_E3SMv2vsShip_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cod_m2_magic,cod_sat_magic, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_COD_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_m2_magic,lwp_sat_magic, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_LWP_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cld_m2_magic,cld_sat_magic, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_CF_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_m2_magic,nd_sat_magic, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_Nd_E3SMv2vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_m2_magic,reff_sat_magic, label1='E3SMv2', label2='Satllite',
#                         outfile=figpath+'statistics_Reff_E3SMv2vsSat_'+site+'.txt')

# calc.mean_std_percentiles([T_magic,T_m_magic,T_m2_magic],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_T_'+site+'.txt')
# calc.mean_std_percentiles([RH_magic,RH_m_magic,RH_m2_magic],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_RH_'+site+'.txt')
# calc.mean_std_percentiles([Ps_magic,Ps_m_magic,Ps_m2_magic],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Ps_'+site+'.txt')
# calc.mean_std_percentiles([cpc10_magic,ncn10_m_magic,ncn10_m2_magic],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CN10_'+site+'.txt')
# calc.mean_std_percentiles([uhsas100_magic,ncn100_m_magic,ncn100_m2_magic],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CN100_'+site+'.txt')
# calc.mean_std_percentiles([ccn2_magic,ccn2_m_magic,ccn2_m2_magic],['Ship','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CCN2_'+site+'.txt')
# calc.mean_std_percentiles([lwp_magic,lwp_sat_magic,lwp_m_magic,lwp_m2_magic],['Ship','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_LWP_'+site+'.txt')
# calc.mean_std_percentiles([cld_magic,cld_sat_magic,cld_m_magic,cld_m2_magic],['Ship','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_CF_'+site+'.txt')
# calc.mean_std_percentiles([nd_wu_magic,nd_sat_magic,nd_m_magic,nd_m2_magic],['Ship','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Nd_'+site+'.txt')
# calc.mean_std_percentiles([reff_wu_magic,reff_sat_magic,reff_m_magic,reff_m2_magic],['Ship','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Reff_'+site+'.txt')
# calc.mean_std_percentiles([cod_sat_magic,cod_m_magic,cod_m2_magic],['Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_COD_'+site+'.txt')


#%% joint histogram
fig,ax = plot.jointhist([uhsas100_magic,ncn100_m_magic,ncn100_m2_magic], [ccn2_magic,ccn2_m_magic,ccn2_m2_magic], 
                    xedges=np.arange(0,600,30),yedges=np.arange(0,600,30), normalize_x=True,
                    xlabel='CN (>100nm) (cm$^{-3}$)', ylabel='CCN (SS=0.2%) (cm$^{-3}$)', vmax=0.5,
                   title=['Ship','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CN100_CCN2_ship_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ccn2_magic,ccn2_magic,ccn2_m_magic,ccn2_m2_magic], [nd_wu_magic,nd_sat_magic,nd_m_magic,nd_m2_magic],
                    xedges=np.arange(0,600,30),yedges=np.arange(0,400,20), normalize_x=True,
                    xlabel='CCN (SS=0.2%) (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', vmax=0.4,
                   title=['Ship','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CCN2_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([nd_wu_magic,nd_sat_magic,nd_m_magic,nd_m2_magic],[lwp_magic,lwp_sat_magic,lwp_m_magic,lwp_m2_magic], 
                    xedges=np.arange(0,400,20),yedges=np.arange(0,300,10), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', vmax=0.4,
                   title=['Ship','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_LWP_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([nd_wu_magic,nd_sat_magic,nd_m_magic,nd_m2_magic],[reff_wu_magic,reff_sat_magic,reff_m_magic,reff_m2_magic],
                    xedges=np.arange(0,400,20),yedges=np.arange(2,20,1), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', vmax=0.25,
                   title=['Ship','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_Reff_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([cod_sat_magic,cod_sat_magic,cod_m_magic,cod_m2_magic],[lwp_magic,lwp_sat_magic,lwp_m_magic,lwp_m2_magic], 
                    xedges=np.arange(0,40,2),yedges=np.arange(0,300,10), normalize_x=True,
                    xlabel='Cloud Optical Depth (N/A)', ylabel='LWP (g/m$^2$)', vmax=0.25,
                   title=['Ship','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_COD_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter plot

fig,ax = plot.scatter([nd_wu_magic.data,nd_sat_magic.data,nd_m_magic.data,nd_m2_magic.data], 
                      [ccn2_magic.data,ccn2_magic.data,ccn2_m_magic.data,ccn2_m2_magic.data],
                     xlimit=(0,600), ylimit=(0,600),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', title=['Ship','Satellite','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_Nd_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.scatter([uhsas100_magic.data,ncn100_m_magic.data,ncn100_m2_magic.data], 
                      [ccn2_magic.data,ccn2_m_magic.data,ccn2_m2_magic.data],
                     xlimit=(0,600), ylimit=(0,600),
                    xlabel='Surface CN (>100nm) (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', title=['Ship','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_CN100_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% heatmaps

# xedges=np.exp(np.arange(np.log(10),6.5,0.5))
# yedges=np.exp(np.arange(np.log(10),6.5,0.5))
fig,ax = plot.heatmap([nd_wu_magic.data,nd_sat_magic.data,nd_m_magic.data,nd_m2_magic.data],
                      [lwp_magic.data,lwp_sat_magic.data,lwp_m_magic.data,lwp_m2_magic.data],
                      [albedo_magic.data,albedo_magic.data,albedo_m_magic.data,albedo_m2_magic.data],vmax=60,
                    xedges=np.arange(0,500,40), yedges=np.arange(10,300,20),
                    # xedges=xedges, yedges=yedges, 
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', zlabel='TOA Albedo (%)',
                    title=['Ground','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'heatmap_CCN2_LWP_Nd_sfc_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)