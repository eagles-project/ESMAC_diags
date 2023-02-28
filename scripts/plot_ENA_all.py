"""
script to generate all plots for ENA surface data

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
site = 'ENA'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'

# path of output figures
figpath= '../figures/'+site+'/sfc_toa/'
if not os.path.exists(figpath):
    os.makedirs(figpath)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = sorted(glob.glob(prep_sfc_path + 'sfc_ACSM_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_acsm = obsdata['time'].load()
org = obsdata['org'].load()
so4 = obsdata['so4'].load()
nh4 = obsdata['nh4'].load()
no3 = obsdata['no3'].load()
chl = obsdata['chl'].load()
obsdata.close()

lst = sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn2 = obsdata['ccn2_fit'].load()
ccn2_2 = obsdata['ccn2_m'].load()
obsdata.close()

lst = sorted(glob.glob(prep_sfc_path + 'sfc_CPC_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_cpc = obsdata['time'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'sfc_CPC_'+site+'_withmask_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_mask = obsdata['time'].load()
cpc10_withmask = obsdata['cpc_masked'].load()
mask_flag = obsdata['f_valid'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'cod_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_cod = obsdata['time'].load()
cod = obsdata['cod'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'LWP_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_lwp = obsdata['time'].load()
lwp_armbe = obsdata['lwp_armbe'].load()
lwp_mfrsr = obsdata['lwp_mfrsr'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'Ndrop_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ndrop = obsdata['time'].load()
ndrop = obsdata['nd'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'reff_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_reff = obsdata['time'].load()
reff = obsdata['reff'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'Nd_Reff_Wu_etal_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_wu = obsdata['time'].load()
reff_wu = obsdata['reff'].load()
nd_wu = obsdata['nd'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'sfc_radiation_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_rad = obsdata['time'].load()
lwdnsfc = obsdata['lwdn'].load()
swdnsfc = obsdata['swdn'].load()
lwupsfc = obsdata['lwup'].load()
swupsfc = obsdata['swup'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'precip_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_pr = obsdata['time'].load()
precip = obsdata['precip'].load()
obsdata.close()

lst = sorted(glob.glob(prep_sfc_path + 'sfc_UHSAS_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_uhsas = obsdata['time'].load()
size_uhsas = obsdata['size'].load()
uhsas100 = obsdata['uhsas100'].load()
uhsas_all = obsdata['uhsas_all'].load()
dlogDp_uhsas = obsdata['dlogDp'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'totcld_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_totcld = obsdata['time'].load()
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_visst = obsdata['tot_cld_visst'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sfc_path + 'cloud_2d_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_cf = obsdata['time'].load()
height_o = obsdata['height'].load()
cloud_2d = obsdata['cloud'].load()
obsdata.close()

# satellite data
lst = sorted(glob.glob(prep_sat_path + 'cod_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
cod_sat = obsdata['cod'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'Reff_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
reff_sat = obsdata['reff'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'LWP_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
lwp_sat = obsdata['lwp'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'Nd_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
nd_sat = obsdata['Nd'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'lwflx_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
lwnettoa = obsdata['lwnettoa'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'swflx_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
swnettoa = obsdata['swnettoa'].load()
obsdata.close()
lst = sorted(glob.glob(prep_sat_path + 'albedo_VISSTgrid_'+site+'*.nc'))
obsdata = xr.open_mfdataset(lst)
time_visst = obsdata['time'].load()
albedo_sat = obsdata['albedo'].load()
obsdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv1_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
bc_m = modeldata['bc'].load()
dst_m = modeldata['dst'].load()
mom_m = modeldata['mom'].load()
pom_m = modeldata['pom'].load()
ncl_m = modeldata['ncl'].load()
so4_m = modeldata['so4'].load()
soa_m = modeldata['soa'].load()
ccn2_m = modeldata['CCN4'].load()
ncn10_m = modeldata['NCN10'].load()
ncn100_m = modeldata['NCN100'].load()
CNsize_m = modeldata['NCNall'].load()
cod_m = modeldata['cod'].load()
reff_m = modeldata['reff'].load()
lwp_m = modeldata['TGCLDLWP'].load()
nd_m = modeldata['Nd_mean'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
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
org_m = pom_m + mom_m + soa_m

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
bc_m2 = modeldata['bc'].load()
dst_m2 = modeldata['dst'].load()
mom_m2 = modeldata['mom'].load()
pom_m2 = modeldata['pom'].load()
ncl_m2 = modeldata['ncl'].load()
so4_m2 = modeldata['so4'].load()
soa_m2 = modeldata['soa'].load()
ccn2_m2 = modeldata['CCN4'].load()
ncn10_m2 = modeldata['NCN10'].load()
ncn100_m2 = modeldata['NCN100'].load()
CNsize_m2 = modeldata['NCNall'].load()
cod_m2 = modeldata['cod'].load()
reff_m2 = modeldata['reff'].load()
lwp_m2 = modeldata['TGCLDLWP'].load()
nd_m2 = modeldata['Nd_mean'].load()
precip_m2 = modeldata['PRECT'].load()
cld_m2 = modeldata['CLDTOT'].load()
lwdnsfc_m2 = modeldata['FLDS'].load()
lwnetsfc_m2 = modeldata['FLNS'].load()
lwnettoa_m2 = modeldata['FLNT'].load()
lwuptoa_m2 = modeldata['FLUT'].load()
swdnsfc_m2 = modeldata['FSDS'].load()
swnetsfc_m2 = modeldata['FSNS'].load()
swdntoa_m2 = modeldata['SOLIN'].load()
swnettoa_m2 = modeldata['FSNT'].load()
swuptoa_m2 = modeldata['FSUTOA'].load()
modeldata.close()
lwupsfc_m2 = lwnetsfc_m2 + lwdnsfc_m2
swupsfc_m2 = swdnsfc_m2 - swnetsfc_m2
albedo_m2 = swuptoa_m2/swdntoa_m2*100
org_m2 = pom_m2 + mom_m2 + soa_m2

filename = prep_model_path + 'E3SMv1_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
height_m = modeldata['height'].load()
cloud_m = modeldata['cloud_z'].load()
modeldata.close()


filename = prep_model_path + 'E3SMv2_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
height_m2 = modeldata['height'].load()
cloud_m2 = modeldata['cloud_z'].load()
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# divided by dlogDp in size distribution
dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
CNsize_m = CNsize_m.T/dlogDp_e3sm
CNsize_m2 = CNsize_m2.T/dlogDp_e3sm
uhsas_all = uhsas_all / dlogDp_uhsas

pdf_uhsas = np.nanmean(uhsas_all,axis=0)
pdf_m = np.nanmean(CNsize_m,axis=0)
pdf_m2 = np.nanmean(CNsize_m2,axis=0)

ndrop[ndrop<10] = np.nan
nd_wu[nd_wu<10] = np.nan
nd_sat[nd_sat<10] = np.nan
nd_m[nd_m<10] = np.nan
nd_m2[nd_m2<10] = np.nan
ndrop[ndrop>800] = np.nan
nd_sat[nd_sat>800] = np.nan
nd_wu[nd_wu>800] = np.nan
nd_m[nd_m>800] = np.nan
nd_m2[nd_m2>800] = np.nan

lwp_armbe[lwp_armbe<10] = np.nan
lwp_mfrsr[lwp_mfrsr<10] = np.nan
lwp_sat[lwp_sat<10] = np.nan
lwp_m[lwp_m<10] = np.nan
lwp_m2[lwp_m2<10] = np.nan

cod[cod<2] = np.nan
cod_sat[cod_sat<2] = np.nan
cod_m[cod_m<2] = np.nan
cod_m2[cod_m2<2] = np.nan
cod[cod>100] = np.nan
cod_sat[cod_sat>100] = np.nan
cod_m[cod_m>100] = np.nan
cod_m2[cod_m2>100] = np.nan

# unit change:
precip_m = precip_m*3600*1000   # m/s to mm/hr
precip_m2 = precip_m2*3600*1000   # m/s to mm/hr
cloud_m = cloud_m*100  # fraction to %
cloud_m2 = cloud_m2*100  # fraction to %
height_o = height_o.data*0.001   # m to km
height_m = height_m.data*0.001   # m to km
height_m2 = height_m2.data*0.001   # m to km

# set a small threshold of E3SM precipitation
precip[precip<0.02] = 0
precip_m[precip_m<0.02] = 0
precip_m2[precip_m2<0.02] = 0

ccn2[ccn2<0] = np.nan
# # remove aerosol mask data
ccn2[mask_flag<0.5] = np.nan
uhsas100[mask_flag<0.5] = np.nan
org[mask_flag<0.5] = np.nan
so4[mask_flag<0.5] = np.nan

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%% bar plot
datagroup0 = [org,so4,nh4,no3,chl, [], []]
datagroup1 = [org_m, so4_m, [], [], [], bc_m, dst_m]
datagroup2 = [org_m2, so4_m2, [], [], [], bc_m2, dst_m2]
dataall=[datagroup0,datagroup1, datagroup2,]
labelall = ['Organic', 'SO$_4$', 'NH$_4$', 'NO$_3$', 'Chl', 'BC', 'Dust']
colorall = ['limegreen', 'red', 'lightblue', 'orange', 'cyan', 'k', 'silver']
fig,ax = plot.bar(dataall, datalabel=['Obs','E3SMv1','E3SMv2',], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
                  title='Aerosol Composition  '+site, varlabel= labelall, colorall=colorall)
fig.savefig(figpath+'bar_composition_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


#%% diurnal cycle
fig,ax = plot.diurnalcycle([org,org_m,org_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Organic '+site, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'diurnalcycle_org_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([so4,so4_m,so4_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Sulfate '+site, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'diurnalcycle_so4_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([ccn2,ccn2_m,ccn2_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='0.2%CCN '+site, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([cpc10_withmask,ncn10_m,ncn10_m2], legend = ['CPC_masked','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='CN(>10nm) '+site, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_CPC10_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([uhsas100,ncn100_m,ncn100_m2], legend = ['UHSAS100','E3SMv1','E3SMv2'], 
                        title='CN(>100nm) '+site, color=['k','r','b'], xlabel='Time (UTC)',ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_CN100_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle( [cod, cod_sat, cod_m, cod_m2], 
                            legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], 
                        title='Cloud optical depth '+site, xlabel='Time (UTC)', ylabel=None)
fig.savefig(figpath+'diurnalcycle_cod_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([lwp_armbe,lwp_sat, lwp_m, lwp_m2], 
                            legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                        title='LWP '+site,  xlabel='Time (UTC)',ylabel="g/m$^2$")
fig.savefig(figpath+'diurnalcycle_LWP_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([ndrop, nd_wu, nd_sat, nd_m,nd_m2], 
                           legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                          title='Nd '+site, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([reff, reff_wu, reff_sat, reff_m, reff_m2], 
                            legend = ['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                        title='droplet effective radius '+site, xlabel='Time (UTC)', ylabel='$\mu$m')
fig.savefig(figpath+'diurnalcycle_reff_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle( [precip,precip_m,precip_m2], 
                            legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        nozero_percentile=True, title='Precipitation '+site, xlabel='Time (UTC)',ylabel='mm/hr')
fig.savefig(figpath+'diurnalcycle_precip_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([cld_arscl,cld_visst,cld_m,cld_m2],
                            legend = ['ARSCL','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                            title='Total cloud fraction '+site, xlabel='Time (UTC)', ylabel="%")
fig.savefig(figpath+'diurnalcycle_totcld_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#%%
fig,ax = plot.diurnalcycle_2d([uhsas_all.T, CNsize_m.T, CNsize_m2.T], 
                              y=[size_uhsas,np.arange(1,3001), np.arange(1,3001)], 
                              title= ['UHSAS','E3SMv1','E3SMv2'],
                              levellist=np.arange(0,3100,100), xlabel='Time (UTC)', ylabel='Diameter (nm)', 
                              ylimit=(3,1000),cmap='jet')
for ax_i in ax:
    ax_i.set_yscale('log')
fig.savefig(figpath+'diurnalcycle_aerosol_size_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle_2d([cloud_2d, cloud_m, cloud_m2], 
                              y = [height_o, height_m, height_m2],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)',  cmap='jet',
                        levellist=np.arange(0,26,1),
                          title= ['Obs', 'E3SMv1', 'E3SMv2',])
fig.savefig(figpath+'diurnalcycle_cloud2d_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% seasonal cycle
fig,ax = plot.seasonalcycle([org,org_m,org_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Organic '+site, ylabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'seasonalcycle_org_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.seasonalcycle([so4,so4_m,so4_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Sulfate '+site, ylabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'seasonalcycle_so4_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle([ccn2_2,ccn2_m,ccn2_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='0.2%CCN '+site, ylabel='cm$^{-3}$')
fig.savefig(figpath+'seasonalcycle_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle([cpc10_withmask,ncn10_m,ncn10_m2], legend = ['CPC_masked','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='CN(>10nm) '+site, ylabel='cm$^{-3}$')
fig.savefig(figpath+'seasonalcycle_CPC10_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.seasonalcycle([uhsas100,ncn100_m,ncn100_m2], legend = ['UHSAS100','E3SMv1','E3SMv2'], 
                        title='CN(>100nm) '+site, color=['k','r','b'], ylabel='cm$^{-3}$')
fig.savefig(figpath+'seasonalcycle_CN100_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle( [cod, cod_sat, cod_m, cod_m2], 
                            legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], 
                        title='Cloud optical depth '+site, ylabel=None)
fig.savefig(figpath+'seasonalcycle_cod_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle([lwp_armbe,lwp_sat, lwp_m, lwp_m2], 
                            legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                        title='LWP '+site,  ylabel="g/m$^2$")
fig.savefig(figpath+'seasonalcycle_LWP_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle([ndrop, nd_wu, nd_sat, nd_m,nd_m2], 
                           legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                          title='Nd '+site, ylabel='cm$^{-3}$')
fig.savefig(figpath+'seasonalcycle_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle([reff, reff_wu, reff_sat, reff_m, reff_m2], 
                            legend = ['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                        title='droplet effective radius '+site, ylabel='$\mu$m')
fig.savefig(figpath+'seasonalcycle_reff_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle( [precip,precip_m,precip_m2], 
                            legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        nozero_percentile=True, title='Precipitation '+site, ylabel='mm/hr')
fig.savefig(figpath+'seasonalcycle_precip_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle([cld_arscl,cld_visst,cld_m,cld_m2],
                            legend = ['ARSCL','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                            title='Total cloud fraction '+site, ylabel="%")
fig.savefig(figpath+'seasonalcycle_totcld_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#%%
fig,ax = plot.seasonalcycle_2d([uhsas_all.T, CNsize_m.T, CNsize_m2.T], 
                              y=[size_uhsas,np.arange(1,3001), np.arange(1,3001)], 
                              title= ['UHSAS','E3SMv1','E3SMv2'],
                              ylabel='Diameter (nm)',  levellist=np.arange(0,5500,200), 
                              ylimit=(3,1000),cmap='jet')
for ax_i in ax:
    ax_i.set_yscale('log')
fig.savefig(figpath+'seasonalcycle_aerosol_size_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle_2d([cloud_2d, cloud_m, cloud_m2], 
                              y = [height_o, height_m, height_m2],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)',  cmap='jet',
                        levellist=np.arange(0,41,1),
                          title= ['Obs', 'E3SMv1', 'E3SMv2',])
fig.savefig(figpath+'seasonalcycle_cloud2d_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% 1d histogram

w1 = np.ones_like(org)/sum(~np.isnan(org.data))
w2 = np.ones_like(org_m)/sum(~np.isnan(org_m.data))
w3 = np.ones_like(org_m2)/sum(~np.isnan(org_m2.data))
fig,ax = plot.hist([org,org_m,org_m2], weights=[w1,w2,w3], bins=np.arange(0,3,0.1),
                    legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
                    title = 'Total Organic '+site, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_org_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w1 = np.ones_like(so4)/sum(~np.isnan(so4.data))
w2 = np.ones_like(so4_m)/sum(~np.isnan(so4_m.data))
w3 = np.ones_like(so4_m2)/sum(~np.isnan(so4_m2.data))
fig,ax = plot.hist([so4,so4_m,so4_m2], weights=[w1,w2,w3], bins=np.arange(0,3,0.1),
                    legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
                    title = 'Sulfate '+site, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_SO4_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w1 = np.ones_like(ccn2)/sum(~np.isnan(ccn2.data))
w2 = np.ones_like(ccn2_m)/sum(~np.isnan(ccn2_m.data))
w3 = np.ones_like(ccn2_m2)/sum(~np.isnan(ccn2_m2.data))
fig,ax = plot.hist([ccn2,ccn2_m,ccn2_m2], weights=[w1,w2,w3], bins=np.arange(0,500,20),
                    legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
                    title = 'CCN (SS=0.2%) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


w0 = np.ones_like(cpc10_withmask)/sum(~np.isnan(cpc10_withmask.data))
w1 = np.ones_like(ncn10_m)/sum(~np.isnan(ncn10_m.data))
w2 = np.ones_like(ncn10_m2)/sum(~np.isnan(ncn10_m2.data))
fig,ax = plot.hist([cpc10_withmask,ncn10_m,ncn10_m2], weights=[w0,w1,w2], bins=np.arange(0,6000,200),
                    legend = ['CPC_masked','E3SMv1','E3SMv2'], color=['k','r','b'],
                    title='Aerosol number (>10nm) '+site,ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CPC10_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(uhsas100)/sum(~np.isnan(uhsas100.data))
w1 = np.ones_like(ncn100_m)/sum(~np.isnan(ncn100_m.data))
w2 = np.ones_like(ncn100_m2)/sum(~np.isnan(ncn100_m2.data))
fig,ax = plot.hist([uhsas100,ncn100_m,ncn100_m2], weights=[w0,w1,w2], bins=np.arange(0,500,20),
                    legend = ['UHSAS100','E3SMv1','E3SMv2'], color=['k','r','b'],
                    title='Aerosol number (>100nm) '+site,  ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CN100_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cod)/sum(~np.isnan(cod.data))
w00 = np.ones_like(cod_sat)/sum(~np.isnan(cod_sat.data))
w1 = np.ones_like(cod_m)/sum(~np.isnan(cod_m.data))
w2 = np.ones_like(cod_m2)/sum(~np.isnan(cod_m2.data))
fig,ax = plot.hist( [cod, cod_sat, cod_m, cod_m2], weights=[w0,w00,w1,w2], 
                    legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                    title='Cloud Optical Depth '+site, bins=np.arange(2,51,2), ylabel='Fraction', xlabel='N/A')
fig.savefig(figpath+'hist_cod_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(lwp_armbe)/sum(~np.isnan(lwp_armbe.data))
# w0 = np.ones_like(lwp_mfrsr)/sum(~np.isnan(lwp_mfrsr.data))
w00 = np.ones_like(lwp_sat)/sum(~np.isnan(lwp_sat.data))
w1 = np.ones_like(lwp_m)/sum(~np.isnan(lwp_m.data))
w2 = np.ones_like(lwp_m2)/sum(~np.isnan(lwp_m2.data))
fig,ax = plot.hist([lwp_mfrsr, lwp_sat, lwp_m, lwp_m2], weights=[w0,w00,w1,w2], 
                    legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                    title='LWP '+site, bins=np.arange(10,410,20), ylabel='Fraction', xlabel="g/m$^2$")
fig.savefig(figpath+'hist_LWP_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(ndrop)/sum(~np.isnan(ndrop.data))
w00 = np.ones_like(nd_wu)/sum(~np.isnan(nd_wu.data))
w000 = np.ones_like(nd_sat)/sum(~np.isnan(nd_sat.data))
w1 = np.ones_like(nd_m)/sum(~np.isnan(nd_m.data))
w2 = np.ones_like(nd_m2)/sum(~np.isnan(nd_m2.data))
fig,ax = plot.hist([ndrop,nd_wu,nd_sat,nd_m,nd_m2],  weights=[w0,w00,w000,w1,w2], 
                   legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                    title = 'Nd '+site, bins=np.arange(0,330,20), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(reff)/sum(~np.isnan(reff.data))
w00 = np.ones_like(reff_wu)/sum(~np.isnan(reff_wu.data))
w000 = np.ones_like(reff_sat)/sum(~np.isnan(reff_sat.data))
w1 = np.ones_like(reff_m)/sum(~np.isnan(reff_m.data))
w2 = np.ones_like(reff_m2)/sum(~np.isnan(reff_m2.data))
fig,ax = plot.hist([reff,reff_wu,reff_sat,reff_m,reff_m2], weights=[w0,w00,w000,w1,w2], 
                    legend = ['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                    title = 'Cloud Effective Radius '+site, bins=np.arange(4,31,1), ylabel='Fraction', xlabel='$\mu$m')
fig.savefig(figpath+'hist_reff_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

pr0 = precip[precip!=0]
prm = precip_m[precip_m!=0]
prm2 = precip_m[precip_m2!=0]
w0 = np.ones_like(pr0)/sum(~np.isnan(pr0.data))
w1 = np.ones_like(prm)/sum(~np.isnan(prm.data))
w2 = np.ones_like(prm2)/sum(~np.isnan(prm2.data))
fig,ax = plot.hist( [pr0,prm,prm2], weights=[w0,w1,w2], legend = ['Obs','E3SMv1','E3SMv2'], 
                    color=['k','r','b'],  bins=np.arange(0,.5,.02), 
                    title = 'Precipitation '+site, ylabel='Fraction', xlabel='mm/hr')
fig.savefig(figpath+'hist_precip_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cld_arscl)/sum(~np.isnan(cld_arscl.data))
w00 = np.ones_like(cld_visst)/sum(~np.isnan(cld_visst.data))
w1 = np.ones_like(cld_m)/sum(~np.isnan(cld_m.data))
w2 = np.ones_like(cld_m2)/sum(~np.isnan(cld_m2.data))
fig,ax = plot.hist([cld_arscl,cld_visst,cld_m,cld_m2], 
                    weights=[w0,w00,w1,w2],  bins=np.arange(0,101,5), 
                    legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                      title = 'Cloud Fraction '+site, ylabel='Fraction', xlabel="%")
fig.savefig(figpath+'hist_totcld_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% mean size distribution
# fig,ax = plot.mean_size([size_uhsas,np.arange(1,3001),np.arange(1,3001)], 
#             [pdf_uhsas,  pdf_m, pdf_m2], 
#             legend = ['UHSAS','E3SMv1','E3SMv2'],color=['k','r','b'], 
#             marker=['o',None,None], linestyles=['none','-','-'],
#             xlimit=(2, 2e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
#             title = 'Mean Aerosol Size Distribution '+site)
# fig.savefig(figpath+'mean_aerosol_size_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.mean_size_witherror([size_uhsas,np.arange(1,3001),np.arange(1,3001)], 
            [uhsas_all,  CNsize_m, CNsize_m2], 
            legend = ['UHSAS','E3SMv1','E3SMv2'],color=['k','r','b'], 
            marker=['o',None,None], linestyles=['none','-','-'],
            xlimit=(10, 3e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
            title = 'Mean Aerosol Size Distribution '+site)
fig.savefig(figpath+'mean_aerosol_size_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statistics
# calc.mean_std_percentiles([org,org_m,org_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_ORG_'+site+'.txt')
# calc.mean_std_percentiles([so4, so4_m, so4_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_SO4_'+site+'.txt')
# calc.mean_std_percentiles([ccn2,ccn2_m,ccn2_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CCN2_'+site+'.txt')
# calc.mean_std_percentiles([cpc10_withmask,ncn10_m,ncn10_m2],legend=['CPC_masked','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CPC10_'+site+'.txt')
# calc.mean_std_percentiles([uhsas100, ncn100_m, ncn100_m2],legend=['UHSAS','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CN100_'+site+'.txt')
# calc.mean_std_percentiles([cod,cod_sat, cod_m, cod_m2],legend=['MFRSR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_COD_'+site+'.txt')
# calc.mean_std_percentiles([reff,reff_wu,reff_sat,reff_m,reff_m2],
#                           legend=['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Reff_'+site+'.txt')
# calc.mean_std_percentiles([lwp_mfrsr,lwp_armbe,lwp_sat,lwp_m,lwp_m2],
#                           legend=['MFRSR','ARMBE','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_LWP_'+site+'.txt')
# calc.mean_std_percentiles([ndrop,nd_wu,nd_sat,nd_m,nd_m2],
#                           legend=['Ndrop','Wu_etal','Nd_satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Nd_'+site+'.txt')
# calc.mean_std_percentiles([precip,precip_m,precip_m2],legend=['Obs','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Precip_'+site+'.txt')
# calc.mean_std_percentiles([cld_arscl,cld_visst,cld_tsi,cld_m,cld_m2],
#                           legend=['ARSCL','Satellite','TSI','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_totcld_'+site+'.txt')


# calc.bias_corrcoef_RMSE(org,org_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_ORG_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(org,org_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_ORG_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(so4, so4_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_SO4_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(so4, so4_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_SO4_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(ccn2,ccn2_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ccn2,ccn2_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(cpc10_withmask,ncn10_m,label1='CPC_masked',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(cpc10_withmask,ncn10_m2,label1='CPC_masked',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(uhsas100, ncn100_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN100_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(uhsas100, ncn100_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN100_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(lwp_armbe, lwp_m,label1='ARMBE',label2='E3SMv1', 
#                         outfile=figpath+'statistics_lwp_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(lwp_armbe, lwp_m2,label1='ARMBE',label2='E3SMv2', 
#                         outfile=figpath+'statistics_lwp_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(ndrop, nd_m,label1='Ndrop',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(ndrop, nd_m2,label1='Ndrop',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(nd_sat, nd_m,label1='Satellite',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(nd_sat, nd_m2,label1='Satellite',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsSat_'+site+'.txt')

# calc.bias_corrcoef_RMSE(reff, reff_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsOBS_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff, reff_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsOBS_'+site+'.txt')

# calc.bias_corrcoef_RMSE(reff_sat, reff_m,label1='Satellite',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsSat_'+site+'.txt')
# calc.bias_corrcoef_RMSE(reff_sat, reff_m2,label1='Satellite',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsSat_'+site+'.txt')

#%% joint histogram
fig,ax = plot.jointhist([uhsas100,ncn100_m,ncn100_m2], [ccn2,ccn2_m,ccn2_m2], 
                    xedges=np.arange(0,600,40),yedges=np.arange(0,600,40), normalize_x=True,
                    xlabel='CN (>100nm) (cm$^{-3}$)', ylabel='CCN (SS=0.2%) (cm$^{-3}$)', vmax=0.5,
                    title=['Ground','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CN100_CCN2_ship_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ccn2,ccn2,ccn2,ccn2_m,ccn2_m2], 
                        [ndrop,nd_wu,nd_sat,nd_m,nd_m2],
                    xedges=np.arange(0,500,30),yedges=np.arange(0,300,20), normalize_x=True,
                    xlabel='CCN (SS=0.2%) (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', vmax=0.4,
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CCN2_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ndrop,nd_wu,nd_sat,nd_m,nd_m2],
                        [lwp_armbe,lwp_armbe,lwp_sat,lwp_m,lwp_m2], 
                    xedges=np.arange(0,300,20),yedges=np.arange(0,300,20), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', vmax=0.4,
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_LWP_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ndrop,nd_wu,nd_sat,nd_m,nd_m2],
                        [reff,reff_wu,reff_sat,reff_m,reff_m2],
                    xedges=np.arange(0,300,20),yedges=np.arange(4,25,2), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', vmax=0.5,
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_Reff_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([cod_sat,cod_sat,cod_m,cod_m2],[lwp_armbe,lwp_sat,lwp_m,lwp_m2], 
                    xedges=np.arange(0,40,3),yedges=np.arange(0,300,20), normalize_x=True,
                    xlabel='Cloud Optical Depth (N/A)', ylabel='LWP (g/m$^2$)', vmax=0.25,
                    title=['Ground','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_COD_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


#%% scatter plot

fig,ax = plot.scatter([ndrop.data, nd_wu.data,nd_sat.data,nd_m.data,nd_m2.data], 
                      [ccn2.data,ccn2.data,ccn2.data,ccn2_m.data,ccn2_m2.data],
                      xlimit=(0,300), ylimit=(0,600),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=False)
fig.savefig(figpath+'scatter_Nd_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.scatter([uhsas100.data,ncn100_m.data,ncn100_m2.data], 
                      [ccn2.data,ccn2_m.data,ccn2_m2.data],
                      xlimit=(0,800), ylimit=(0,800),
                    xlabel='Surface CN (>100nm) (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                    title=['Obs','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_CN100_CCN2_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% heatmaps

# xedges=np.exp(np.arange(np.log(10),6.5,0.5))
# yedges=np.exp(np.arange(np.log(10),6.5,0.5))
fig,ax = plot.heatmap([ndrop.data, nd_wu.data,nd_sat.data,nd_m.data,nd_m2.data],
                      [lwp_armbe.data,lwp_armbe.data,lwp_sat.data,lwp_m.data,lwp_m2.data],
                      [albedo_sat,albedo_sat,albedo_sat,albedo_m,albedo_m2],vmax=60,
                    xedges=np.arange(0,300,20), yedges=np.arange(10,300,20),
                    # xedges=xedges, yedges=yedges, 
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', zlabel='TOA Albedo (%)',
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'heatmap_Albedo_vs_Nd_LWP_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)