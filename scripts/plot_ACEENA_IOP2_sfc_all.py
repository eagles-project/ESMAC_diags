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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings
# set site name and datapath

# set site name.
site = 'ACEENA'

prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'

time_aceena = pd.date_range(start='2018-01-15', end='2018-02-19', freq="3600s")
IOP = 'IOP2'
            
# path of output figures
figpath= '../figures/'+site+'/sfc_toa/'
if not os.path.exists(figpath):
    os.makedirs(figpath)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
obsdata = xr.open_dataset(prep_sfc_path + 'sfc_ACSM_'+site+'.nc')
org = obsdata['org'].load()
so4 = obsdata['so4'].load()
nh4 = obsdata['nh4'].load()
no3 = obsdata['no3'].load()
chl = obsdata['chl'].load()
obsdata.close()
org_aceena = org.sel(time=time_aceena)
so4_aceena = so4.sel(time=time_aceena)
nh4_aceena = nh4.sel(time=time_aceena)
no3_aceena = no3.sel(time=time_aceena)
chl_aceena = chl.sel(time=time_aceena)

obsdata = xr.open_dataset(prep_sfc_path + 'totcld_'+site+'.nc')
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_visst = obsdata['tot_cld_visst'].load()
obsdata.close()
cld_arscl_aceena = cld_arscl.sel(time=time_aceena)
cld_tsi_aceena = cld_tsi.sel(time=time_aceena)
cld_visst_aceena = cld_visst.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'cloud_2d_'+site+'.nc')
height_o = obsdata['height'].load()
cloud_2d = obsdata['cloud'].load()
obsdata.close()
cloud_2d_aceena = cloud_2d.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'sfc_CCN_'+site+'.nc')
ccn2 = obsdata['CCN2'].load()
obsdata.close()
ccn2_aceena = ccn2.sel(time=time_aceena)

obsdata = xr.open_dataset(prep_sfc_path + 'sfc_CPC_'+site+'.nc')
cpc10 = obsdata['cpc10'].load()
obsdata.close()
cpc10_aceena = cpc10.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'sfc_CPC_'+site+'_withmask.nc')
time_mask = obsdata['time'].load()
cpc_mask = obsdata['cpc_masked'].load()
mask_flag = obsdata['f_valid'].load()
obsdata.close()
cpc10_withmask_aceena = cpc_mask.sel(time=time_aceena)
mask_flag = mask_flag.sel(time=time_aceena)

obsdata = xr.open_dataset(prep_sfc_path + 'cod_'+site+'.nc')
cod = obsdata['cod'].load()
obsdata.close()
cod_aceena = cod.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'LWP_'+site+'.nc')
lwp_armbe = obsdata['lwp_armbe'].load()
lwp_mfrsr = obsdata['lwp_mfrsr'].load()
obsdata.close()
lwp_armbe_aceena = lwp_armbe.sel(time=time_aceena)
lwp_mfrsr_aceena = lwp_mfrsr.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'Ndrop_'+site+'.nc')
ndrop = obsdata['cdnc'].load()
obsdata.close()
ndrop_aceena = ndrop.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'reff_'+site+'.nc')
reff = obsdata['reff'].load()
obsdata.close()
reff_aceena = reff.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sfc_path + 'Nd_Reff_Wu_etal_'+site+'.nc')
time_wu = obsdata['time'].load()
reff_wu = obsdata['reff'].load()
nd_wu = obsdata['cdnc'].load()
obsdata.close()
reff_wu_aceena = reff_wu.sel(time=time_aceena)
nd_wu_aceena = nd_wu.sel(time=time_aceena)

obsdata = xr.open_dataset(prep_sfc_path + 'precip_'+site+'.nc')
precip = obsdata['precip'].load()
obsdata.close()
precip_aceena = precip.sel(time=time_aceena)

obsdata = xr.open_dataset(prep_sfc_path + 'sfc_radiation_'+site+'.nc')
lwdnsfc = obsdata['lwdn'].load()
swdnsfc = obsdata['swdn'].load()
lwupsfc = obsdata['lwup'].load()
swupsfc = obsdata['swup'].load()
obsdata.close()
lwnetsfc = lwupsfc - lwdnsfc
swnetsfc = swdnsfc - swupsfc
lwnetsfc_aceena = lwnetsfc.sel(time=time_aceena)
swnetsfc_aceena = swnetsfc.sel(time=time_aceena)

obsdata = xr.open_dataset(prep_sfc_path + 'sfc_UHSAS_'+site+'.nc')
size_uhsas = obsdata['size'].load()
dmin_aceena = obsdata['size_low'].load()
dmax_aceena = obsdata['size_high'].load()
uhsas100 = obsdata['uhsas100'].load()
uhsas_all = obsdata['uhsas_all'].load()
obsdata.close()
uhsas100_aceena = uhsas100.sel(time=time_aceena)
uhsasall_aceena = uhsas_all.sel(time=time_aceena)
dlogDp_uhsas = np.mean(np.log10(dmax_aceena/dmin_aceena))

obsdata = xr.open_mfdataset(prep_sfc_path + 'cloudheight_ARSCL_'+site+'.nc')
cth = obsdata['cth'].load()
cbh = obsdata['cbh'].load()
cths = obsdata['cths'].load()
obsdata.close()
cth_aceena = cth.sel(time=time_aceena)
cbh_aceena = cbh.sel(time=time_aceena)

# satellite data
obsdata = xr.open_dataset(prep_sat_path + 'albedo_VISSTgrid_'+site+'.nc')
albedo_sat = obsdata['albedo'].load()
obsdata.close()
albedo_aceena = albedo_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'cloudfraction_VISSTgrid_'+site+'.nc')
cfall_sat = obsdata['cldtot'].load()
cflow_sat = obsdata['cldlow'].load()
obsdata.close()
cfall_sat_aceena = cfall_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc')
ctt_sat = obsdata['ctt'].load()
cth_sat = obsdata['cth'].load()
obsdata.close()
ctt_sat_aceena = ctt_sat.sel(time=time_aceena)
cth_sat_aceena = cth_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'cod_VISSTgrid_'+site+'.nc')
cod_sat = obsdata['cod'].load()
obsdata.close()
cod_sat_aceena = cod_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc')
reff_sat = obsdata['reff'].load()
obsdata.close()
reff_sat_aceena = reff_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc')
lwp_sat = obsdata['lwp'].load()
obsdata.close()
lwp_sat_aceena = lwp_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc')
nd_sat = obsdata['Nd'].load()
obsdata.close()
nd_sat_aceena = nd_sat.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc')
lwnettoa = obsdata['lwnettoa'].load()
obsdata.close()
lwnettoa_aceena = lwnettoa.sel(time=time_aceena)
obsdata = xr.open_dataset(prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc')
swnettoa = obsdata['swnettoa'].load()
obsdata.close()
swnettoa_aceena = swnettoa.sel(time=time_aceena)

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
ncn3_m = modeldata['NCN3'].load()
ncn10_m = modeldata['NCN10'].load()
ncn100_m = modeldata['NCN100'].load()
CNsize_m = modeldata['NCNall'].load()
cod_m = modeldata['cod'].load()
reff_m = modeldata['reff'].load()
lwp_m = modeldata['TGCLDLWP'].load()
nd_m = modeldata['Nd_mean'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
cbh_m = modeldata['cbh'].load()
cth_m = modeldata['cth'].load()
Hcld_m = modeldata['clddepth'].load()
cldlow_m = modeldata['CLDLOW'].load()
cldmid_m = modeldata['CLDMED'].load()
cldhgh_m = modeldata['CLDHGH'].load()
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
bc_m_aceena = bc_m.sel(time=time_aceena)
dst_m_aceena = dst_m.sel(time=time_aceena)
org_m_aceena = org_m.sel(time=time_aceena)
so4_m_aceena = so4_m.sel(time=time_aceena)
ncl_m_aceena = ncl_m.sel(time=time_aceena)
ccn2_m_aceena = ccn2_m.sel(time=time_aceena)
ncn3_m_aceena = ncn3_m.sel(time=time_aceena)
ncn10_m_aceena = ncn10_m.sel(time=time_aceena)
ncn100_m_aceena = ncn100_m.sel(time=time_aceena)
CNsize_m_aceena = CNsize_m.sel(time=time_aceena)
cod_m_aceena = cod_m.sel(time=time_aceena)
reff_m_aceena = reff_m.sel(time=time_aceena)
lwp_m_aceena = lwp_m.sel(time=time_aceena)
nd_m_aceena = nd_m.sel(time=time_aceena)
precip_m_aceena = precip_m.sel(time=time_aceena)
cld_m_aceena = cld_m.sel(time=time_aceena)
cbh_m_aceena = cbh_m.sel(time=time_aceena)
cth_m_aceena = cth_m.sel(time=time_aceena)
Hcld_m_aceena = Hcld_m.sel(time=time_aceena)
cldlow_m_aceena = cldlow_m.sel(time=time_aceena)
cldmid_m_aceena = cldmid_m.sel(time=time_aceena)
cldhgh_m_aceena = cldhgh_m.sel(time=time_aceena)
lwnetsfc_m_aceena = lwnetsfc_m.sel(time=time_aceena)
lwnettoa_m_aceena = lwnettoa_m.sel(time=time_aceena)
swnetsfc_m_aceena = swnetsfc_m.sel(time=time_aceena)
swnettoa_m_aceena = swnettoa_m.sel(time=time_aceena)
albedo_m_aceena = albedo_m.sel(time=time_aceena)

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
ncn3_m2 = modeldata['NCN3'].load()
ncn10_m2 = modeldata['NCN10'].load()
ncn100_m2 = modeldata['NCN100'].load()
CNsize_m2 = modeldata['NCNall'].load()
cod_m2 = modeldata['cod'].load()
reff_m2 = modeldata['reff'].load()
lwp_m2 = modeldata['TGCLDLWP'].load()
nd_m2 = modeldata['Nd_mean'].load()
precip_m2 = modeldata['PRECT'].load()
cld_m2 = modeldata['CLDTOT'].load()
cbh_m2 = modeldata['cbh'].load()
cth_m2 = modeldata['cth'].load()
Hcld_m2 = modeldata['clddepth'].load()
cldlow_m2 = modeldata['CLDLOW'].load()
cldmid_m2 = modeldata['CLDMED'].load()
cldhgh_m2 = modeldata['CLDHGH'].load()
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
bc_m2_aceena = bc_m2.sel(time=time_aceena)
dst_m2_aceena = dst_m2.sel(time=time_aceena)
org_m2_aceena = org_m2.sel(time=time_aceena)
so4_m2_aceena = so4_m2.sel(time=time_aceena)
ncl_m2_aceena = ncl_m2.sel(time=time_aceena)
ccn2_m2_aceena = ccn2_m2.sel(time=time_aceena)
ncn3_m2_aceena = ncn3_m2.sel(time=time_aceena)
ncn10_m2_aceena = ncn10_m2.sel(time=time_aceena)
ncn100_m2_aceena = ncn100_m2.sel(time=time_aceena)
CNsize_m2_aceena = CNsize_m2.sel(time=time_aceena)
cod_m2_aceena = cod_m2.sel(time=time_aceena)
reff_m2_aceena = reff_m2.sel(time=time_aceena)
lwp_m2_aceena = lwp_m2.sel(time=time_aceena)
nd_m2_aceena = nd_m2.sel(time=time_aceena)
precip_m2_aceena = precip_m2.sel(time=time_aceena)
cld_m2_aceena = cld_m2.sel(time=time_aceena)
cbh_m2_aceena = cbh_m2.sel(time=time_aceena)
cth_m2_aceena = cth_m2.sel(time=time_aceena)
Hcld_m2_aceena = Hcld_m2.sel(time=time_aceena)
cldlow_m2_aceena = cldlow_m2.sel(time=time_aceena)
cldmid_m2_aceena = cldmid_m2.sel(time=time_aceena)
cldhgh_m2_aceena = cldhgh_m2.sel(time=time_aceena)
lwnetsfc_m2_aceena = lwnetsfc_m2.sel(time=time_aceena)
lwnettoa_m2_aceena = lwnettoa_m2.sel(time=time_aceena)
swnetsfc_m2_aceena = swnetsfc_m2.sel(time=time_aceena)
swnettoa_m2_aceena = swnettoa_m2.sel(time=time_aceena)
albedo_m2_aceena = albedo_m2.sel(time=time_aceena)

filename = prep_model_path + 'E3SMv1_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
height_m = modeldata['height'].load()
cf_e3sm = modeldata['cloud_z'].load()
modeldata.close()
cloud_m_aceena = cf_e3sm.sel(time=time_aceena)

filename = prep_model_path + 'E3SMv2_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
height_m2 = modeldata['height'].load()
cf_e3sm2 = modeldata['cloud_z'].load()
modeldata.close()
cloud_m2_aceena = cf_e3sm2.sel(time=time_aceena)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# divided by dlogDp in size distribution
dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
CNsize_m_aceena = CNsize_m_aceena.T/dlogDp_e3sm
CNsize_m2_aceena = CNsize_m2_aceena.T/dlogDp_e3sm
uhsasall_aceena = uhsasall_aceena / dlogDp_uhsas

pdf_uhsas_aceena = np.nanmean(uhsasall_aceena,axis=0)
pdf_m_aceena = np.nanmean(CNsize_m_aceena,axis=0)
pdf_m2_aceena = np.nanmean(CNsize_m2_aceena,axis=0)

ndrop_aceena[ndrop_aceena<10] = np.nan
nd_wu_aceena[nd_wu_aceena<10] = np.nan
nd_sat_aceena[nd_sat_aceena<10] = np.nan
nd_m_aceena[nd_m_aceena<10] = np.nan
nd_m2_aceena[nd_m2_aceena<10] = np.nan
ndrop_aceena[ndrop_aceena>500] = np.nan
nd_sat_aceena[nd_sat_aceena>500] = np.nan
nd_wu_aceena[nd_wu_aceena>500] = np.nan
nd_m_aceena[nd_m_aceena>500] = np.nan
nd_m2_aceena[nd_m2_aceena>500] = np.nan

lwp_armbe_aceena[lwp_armbe_aceena<20] = np.nan
lwp_mfrsr_aceena[lwp_mfrsr_aceena<20] = np.nan
lwp_sat_aceena[lwp_sat_aceena<20] = np.nan
lwp_m_aceena[lwp_m_aceena<20] = np.nan
lwp_m2_aceena[lwp_m2_aceena<20] = np.nan

cod_aceena[cod_aceena<2] = np.nan
cod_sat_aceena[cod_sat_aceena<2] = np.nan
cod_m_aceena[cod_m_aceena<2] = np.nan
cod_m2_aceena[cod_m2_aceena<2] = np.nan
cod_aceena[cod_aceena>100] = np.nan
cod_sat_aceena[cod_sat_aceena>100] = np.nan
cod_m_aceena[cod_m_aceena>100] = np.nan
cod_m2_aceena[cod_m2_aceena>100] = np.nan

# unit change:
precip_m_aceena = precip_m_aceena*3600*1000   # m/s to mm/hr
precip_m2_aceena = precip_m2_aceena*3600*1000   # m/s to mm/hr
cloud_m_aceena = cloud_m_aceena*100  # fraction to %
cloud_m2_aceena = cloud_m2_aceena*100  # fraction to %
height_o = height_o.data*0.001   # m to km
height_m = height_m.data*0.001   # m to km
height_m2 = height_m2.data*0.001   # m to km

# set a small threshold of E3SM precipitation
precip_m_aceena[precip_m_aceena<0.02] = 0
precip_m2_aceena[precip_m2_aceena<0.02] = 0

# remove aerosol mask data
ccn2_aceena[mask_flag<0.5] = np.nan
uhsas100_aceena[mask_flag<0.5] = np.nan
org_aceena[mask_flag<0.5] = np.nan
so4_aceena[mask_flag<0.5] = np.nan

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%% bar plot
datagroup0 = [org_aceena,so4_aceena,nh4_aceena,no3_aceena,chl_aceena, [], []]
datagroup1 = [org_m_aceena, so4_m_aceena, [], [], [], bc_m_aceena, dst_m_aceena]
datagroup2 = [org_m2_aceena, so4_m2_aceena, [], [], [], bc_m2_aceena, dst_m2_aceena]
dataall=[datagroup0,datagroup1, datagroup2,]
labelall = ['Organic', 'SO$_4$', 'NH$_4$', 'NO$_3$', 'Chl', 'BC', 'Dust']
colorall = ['limegreen', 'red', 'lightblue', 'orange', 'cyan', 'k', 'silver']
fig,ax = plot.bar(dataall, datalabel=['Obs','E3SMv1','E3SMv2',], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
                  title='Aerosol Composition  '+site+' '+IOP, varlabel= labelall, colorall=colorall)
fig.savefig(figpath+'bar_composition_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% timeseries
fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [org_aceena,org_m_aceena,org_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Total Organic '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [so4_aceena,so4_m_aceena,so4_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'],  
                          title='Sulfate '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [ccn2_aceena,ccn2_m_aceena,ccn2_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='0.2%CCN '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [cpc10_withmask_aceena,ncn10_m_aceena,ncn10_m2_aceena], 
                          legend = ['CPC_masked','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='CN(>10nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [uhsas100_aceena,ncn100_m_aceena,ncn100_m2_aceena], 
                        legend = ['UHSAS','E3SMv1','E3SMv2'], color=['k','r','b'],
                        title='CN(>100nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena,time_aceena], [cod_aceena,cod_sat_aceena,cod_m_aceena,cod_m2_aceena], 
                          legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], #marker='.',
                        title='cloud optical depth '+site+' '+IOP, xlabel=None, ylabel=None)
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena,time_aceena], [lwp_armbe_aceena, lwp_sat_aceena, lwp_m_aceena, lwp_m2_aceena], 
                        legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                        title='LWP '+site+' '+IOP,xlabel=None, ylabel="g/m$^2$")
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena,time_aceena,time_aceena], 
                         [ndrop_aceena,nd_wu_aceena,nd_sat_aceena, nd_m_aceena, nd_m2_aceena], 
                          legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'], #marker='.',
                          title='Nd '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena,time_aceena,time_aceena], 
                         [reff_aceena,reff_wu_aceena,reff_sat_aceena,reff_m_aceena,reff_m2_aceena], 
                          legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'], #marker='.',
                        title='Reff '+site+' '+IOP,xlabel=None, ylabel='$\mu$m')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [precip_aceena,precip_m_aceena,precip_m2_aceena],  
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='Precip '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='mm/hr')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [lwnetsfc_aceena,lwnetsfc_m_aceena,lwnetsfc_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='Sfc. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_LWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [swnetsfc_aceena,swnetsfc_m_aceena,swnetsfc_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='Sfc. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_SWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [lwnettoa_aceena,lwnettoa_m_aceena,lwnettoa_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='TOA. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_LWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena], [swnettoa_aceena,swnettoa_m_aceena,swnettoa_m2_aceena], 
                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='TOA. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_SWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time_aceena,time_aceena,time_aceena,time_aceena], [cld_arscl_aceena,cld_visst_aceena,cld_m_aceena,cld_m2_aceena], 
                        legend = ['ARSCL','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                        title='Cloud fraction '+site+' '+IOP,xlabel=None, ylabel="%")
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%%
fig,ax = plot.timeseries_size([time_aceena,time_aceena,time_aceena], 
                              [size_uhsas, np.arange(1,3001), np.arange(1,3001)], 
                              [uhsasall_aceena.T.data,CNsize_m_aceena.T.data, CNsize_m2_aceena.T.data], 
                              legend = ['UHSAS','E3SMv1','E3SMv2'],
                          ylabel='Diameter (nm)', ylimit=(3,1000),
                          title = 'Aerosol Size Distribution (dN/dlogDp, cm$^{-3}$)')
for ax_i in ax:
    ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries_2d([time_aceena,time_aceena,time_aceena], 
                            [height_o, height_m, height_m2], 
                            [cloud_2d_aceena.T.data, cloud_m_aceena.T.data, cloud_m2_aceena.T.data], 
                              yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (m)',cmap='jet', #ylimit=(3,1000),
                              legend = ['Obs','E3SMv1','E3SMv2'], title = 'Cloud Fraction (%)')
for ax_i in ax:
    ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'cloud_2d_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


#%% diurnal cycle
fig,ax = plot.diurnalcycle([org_aceena,org_m_aceena,org_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Organic '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'diurnalcycle_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([so4_aceena,so4_m_aceena,so4_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                          title='Sulfate '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'diurnalcycle_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([ccn2_aceena,ccn2_m_aceena,ccn2_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='0.2%CCN '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([cpc10_withmask_aceena,ncn10_m_aceena,ncn10_m2_aceena], legend = ['CPC_masked','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='CN(>10nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([uhsas100_aceena,ncn100_m_aceena,ncn100_m2_aceena], legend = ['UHSAS100','E3SMv1','E3SMv2'], 
                        title='CN(>100nm) '+site+' '+IOP, color=['k','r','b'], xlabel='Time (UTC)',ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle( [cod_aceena, cod_sat_aceena, cod_m_aceena, cod_m2_aceena], 
                            legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], 
                        title='Cloud optical depth '+site+' '+IOP, xlabel='Time (UTC)', ylabel=None)
fig.savefig(figpath+'diurnalcycle_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([lwp_armbe_aceena,lwp_sat_aceena, lwp_m_aceena, lwp_m2_aceena], 
                            legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                        title='LWP '+site+' '+IOP,  xlabel='Time (UTC)',ylabel="g/m$^2$")
fig.savefig(figpath+'diurnalcycle_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([ndrop_aceena, nd_wu_aceena, nd_sat_aceena, nd_m_aceena,nd_m2_aceena], 
                           legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                          title='Nd '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([reff_aceena, reff_wu_aceena, reff_sat_aceena, reff_m_aceena, reff_m2_aceena], 
                            legend = ['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                        title='droplet effective radius '+site+' '+IOP, xlabel='Time (UTC)', ylabel='$\mu$m')
fig.savefig(figpath+'diurnalcycle_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle( [precip_aceena,precip_m_aceena,precip_m2_aceena], 
                            legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        nozero_percentile=True, title='Precipitation '+site+' '+IOP, xlabel='Time (UTC)',ylabel='mm/hr')
fig.savefig(figpath+'diurnalcycle_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([lwnetsfc_aceena,lwnetsfc_m_aceena,lwnetsfc_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='Sfc. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)',ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_LWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([swnetsfc_aceena,swnetsfc_m_aceena,swnetsfc_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='Sfc. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_SWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([lwnettoa_aceena,lwnettoa_m_aceena,lwnettoa_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='TOA. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_LWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.diurnalcycle([swnettoa_aceena,swnettoa_m_aceena,swnettoa_m2_aceena], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
                        title='TOA. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_SWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle([cld_arscl_aceena,cld_visst_aceena,cld_m_aceena,cld_m2_aceena],
                            legend = ['ARSCL','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                            title='Total cloud fraction '+site+' '+IOP, xlabel='Time (UTC)', ylabel="%")
fig.savefig(figpath+'diurnalcycle_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle_2d([uhsasall_aceena.T, CNsize_m_aceena.T, CNsize_m2_aceena.T], 
                              y=[size_uhsas,np.arange(1,3001), np.arange(1,3001)], 
                              title= ['UHSAS','E3SMv1','E3SMv2'],
                              levellist=np.arange(0,7500,200), xlabel='Time (UTC)', ylabel='Diameter (nm)', 
                              ylimit=(3,1000),cmap='jet')
for ax_i in ax:
    ax_i.set_yscale('log')
fig.savefig(figpath+'diurnalcycle_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.diurnalcycle_2d([cloud_2d_aceena, cloud_m_aceena, cloud_m2_aceena], 
                              y = [height_o, height_m, height_m2],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)',  cmap='jet',
                        levellist=np.arange(0,41,1),
                          title= ['Obs', 'E3SMv1', 'E3SMv2',])
fig.savefig(figpath+'diurnalcycle_cloud2d_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% 1d histogram

w1 = np.ones_like(org_aceena)/sum(~np.isnan(org_aceena.data))
w2 = np.ones_like(org_m_aceena)/sum(~np.isnan(org_m_aceena.data))
w3 = np.ones_like(org_m2_aceena)/sum(~np.isnan(org_m2_aceena.data))
fig,ax = plot.hist([org_aceena,org_m_aceena,org_m2_aceena], weights=[w1,w2,w3], bins=np.arange(0,1,0.05),
                    legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
                    title = 'Total Organic '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w1 = np.ones_like(so4_aceena)/sum(~np.isnan(so4_aceena.data))
w2 = np.ones_like(so4_m_aceena)/sum(~np.isnan(so4_m_aceena.data))
w3 = np.ones_like(so4_m2_aceena)/sum(~np.isnan(so4_m2_aceena.data))
fig,ax = plot.hist([so4_aceena,so4_m_aceena,so4_m2_aceena], weights=[w1,w2,w3], bins=np.arange(0,1,0.05),
                    legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
                    title = 'Sulfate '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
fig.savefig(figpath+'hist_SO4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w1 = np.ones_like(ccn2_aceena)/sum(~np.isnan(ccn2_aceena.data))
w2 = np.ones_like(ccn2_m_aceena)/sum(~np.isnan(ccn2_m_aceena.data))
w3 = np.ones_like(ccn2_m2_aceena)/sum(~np.isnan(ccn2_m2_aceena.data))
fig,ax = plot.hist([ccn2_aceena,ccn2_m_aceena,ccn2_m2_aceena], weights=[w1,w2,w3], bins=np.arange(0,300,15),
                    legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
                    title = 'CCN (SS=0.2%) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


w0 = np.ones_like(cpc10_withmask_aceena)/sum(~np.isnan(cpc10_withmask_aceena.data))
w1 = np.ones_like(ncn10_m_aceena)/sum(~np.isnan(ncn10_m_aceena.data))
w2 = np.ones_like(ncn10_m2_aceena)/sum(~np.isnan(ncn10_m2_aceena.data))
fig,ax = plot.hist([cpc10_withmask_aceena,ncn10_m_aceena,ncn10_m2_aceena], weights=[w0,w1,w2], bins=np.arange(0,3000,100),
                    legend = ['CPC_masked','E3SMv1','E3SMv2'], color=['k','r','b'],
                    title='Aerosol number (>10nm) '+site+' '+IOP,ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(uhsas100_aceena)/sum(~np.isnan(uhsas100_aceena.data))
w1 = np.ones_like(ncn100_m_aceena)/sum(~np.isnan(ncn100_m_aceena.data))
w2 = np.ones_like(ncn100_m2_aceena)/sum(~np.isnan(ncn100_m2_aceena.data))
fig,ax = plot.hist([uhsas100_aceena,ncn100_m_aceena,ncn100_m2_aceena], weights=[w0,w1,w2], bins=np.arange(0,200,10),
                    legend = ['UHSAS100','E3SMv1','E3SMv2'], color=['k','r','b'],
                    title='Aerosol number (>100nm) '+site+' '+IOP,  ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cod_aceena)/sum(~np.isnan(cod_aceena.data))
w00 = np.ones_like(cod_sat_aceena)/sum(~np.isnan(cod_sat_aceena.data))
w1 = np.ones_like(cod_m_aceena)/sum(~np.isnan(cod_m_aceena.data))
w2 = np.ones_like(cod_m2_aceena)/sum(~np.isnan(cod_m2_aceena.data))
fig,ax = plot.hist( [cod_aceena, cod_sat_aceena, cod_m_aceena, cod_m2_aceena], weights=[w0,w00,w1,w2], 
                    legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                    title='Cloud Optical Depth '+site+' '+IOP, bins=np.arange(0,41,2), ylabel='Fraction', xlabel='N/A')
fig.savefig(figpath+'hist_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(lwp_armbe_aceena)/sum(~np.isnan(lwp_armbe_aceena.data))
# w0 = np.ones_like(lwp_mfrsr)/sum(~np.isnan(lwp_mfrsr.data))
w00 = np.ones_like(lwp_sat_aceena)/sum(~np.isnan(lwp_sat_aceena.data))
w1 = np.ones_like(lwp_m_aceena)/sum(~np.isnan(lwp_m_aceena.data))
w2 = np.ones_like(lwp_m2_aceena)/sum(~np.isnan(lwp_m2_aceena.data))
fig,ax = plot.hist([lwp_mfrsr_aceena, lwp_sat_aceena, lwp_m_aceena, lwp_m2_aceena], weights=[w0,w00,w1,w2], 
                    legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                    title='LWP '+site+' '+IOP, bins=np.arange(10,310,20), ylabel='Fraction', xlabel="g/m$^2$")
fig.savefig(figpath+'hist_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(ndrop_aceena)/sum(~np.isnan(ndrop_aceena.data))
w00 = np.ones_like(nd_wu_aceena)/sum(~np.isnan(nd_wu_aceena.data))
w000 = np.ones_like(nd_sat_aceena)/sum(~np.isnan(nd_sat_aceena.data))
w1 = np.ones_like(nd_m_aceena)/sum(~np.isnan(nd_m_aceena.data))
w2 = np.ones_like(nd_m2_aceena)/sum(~np.isnan(nd_m2_aceena.data))
fig,ax = plot.hist([ndrop_aceena,nd_wu_aceena,nd_sat_aceena,nd_m_aceena,nd_m2_aceena],  weights=[w0,w00,w000,w1,w2], 
                   legend = ['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                    title = 'Nd '+site+' '+IOP, bins=np.arange(0,230,10), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(reff_aceena)/sum(~np.isnan(reff_aceena.data))
w00 = np.ones_like(reff_wu_aceena)/sum(~np.isnan(reff_wu_aceena.data))
w000 = np.ones_like(reff_sat_aceena)/sum(~np.isnan(reff_sat_aceena.data))
w1 = np.ones_like(reff_m_aceena)/sum(~np.isnan(reff_m_aceena.data))
w2 = np.ones_like(reff_m2_aceena)/sum(~np.isnan(reff_m2_aceena.data))
fig,ax = plot.hist([reff_aceena,reff_wu_aceena,reff_sat_aceena,reff_m_aceena,reff_m2_aceena], weights=[w0,w00,w000,w1,w2], 
                    legend = ['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'], color=['k','m','gray','r','b'],
                    title = 'Cloud Effective Radius '+site+' '+IOP, bins=np.arange(5,28,1), ylabel='Fraction', xlabel='$\mu$m')
fig.savefig(figpath+'hist_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

pr0 = precip_aceena[precip_aceena!=0]
prm = precip_m_aceena[precip_m_aceena!=0]
prm2 = precip_m_aceena[precip_m2_aceena!=0]
w0 = np.ones_like(pr0)/sum(~np.isnan(pr0.data))
w1 = np.ones_like(prm)/sum(~np.isnan(prm.data))
w2 = np.ones_like(prm2)/sum(~np.isnan(prm2.data))
fig,ax = plot.hist( [pr0,prm,prm2], weights=[w0,w1,w2], legend = ['Obs','E3SMv1','E3SMv2'], 
                    color=['k','r','b'],  bins=np.arange(0,.5,.02), 
                    title = 'Precipitation '+site+' '+IOP, ylabel='Fraction', xlabel='mm/hr')
fig.savefig(figpath+'hist_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cld_arscl_aceena)/sum(~np.isnan(cld_arscl_aceena.data))
w00 = np.ones_like(cld_visst_aceena)/sum(~np.isnan(cld_visst_aceena.data))
w1 = np.ones_like(cld_m_aceena)/sum(~np.isnan(cld_m_aceena.data))
w2 = np.ones_like(cld_m2_aceena)/sum(~np.isnan(cld_m2_aceena.data))
fig,ax = plot.hist([cld_arscl_aceena,cld_visst_aceena,cld_m_aceena,cld_m2_aceena], 
                    weights=[w0,w00,w1,w2],  bins=np.arange(0,101,5), 
                    legend = ['ARMBE','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                      title = 'Cloud Fraction '+site+' '+IOP, ylabel='Fraction', xlabel="%")
fig.savefig(figpath+'hist_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% mean size distribution
fig,ax = plot.mean_size([size_uhsas,np.arange(1,3001),np.arange(1,3001)], 
            [pdf_uhsas_aceena,  pdf_m_aceena, pdf_m2_aceena], 
            legend = ['UHSAS','E3SMv1','E3SMv2'],color=['k','r','b'], 
            marker=['o',None,None], linestyles=['none','-','-'],
            xlimit=(2, 2e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
            title = 'Mean Aerosol Size Distribution '+site)
fig.savefig(figpath+'mean_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statistics
# calc.mean_std_percentiles([org_aceena,org_m_aceena,org_m2_aceena],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_ORG_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([so4_aceena, so4_m_aceena, so4_m2_aceena],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_SO4_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn2_aceena,ccn2_m_aceena,ccn2_m2_aceena],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CCN2_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cpc10_withmask_aceena,ncn10_m_aceena,ncn10_m2_aceena],legend=['CPC_masked','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([uhsas100_aceena, ncn100_m_aceena, ncn100_m2_aceena],legend=['UHSAS','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CN100_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cod_aceena,cod_sat_aceena, cod_m_aceena, cod_m2_aceena],legend=['MFRSR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_COD_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([reff_aceena,reff_wu_aceena,reff_sat_aceena,reff_m_aceena,reff_m2_aceena],
#                           legend=['MFRSR','Wu_etal','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Reff_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([lwp_mfrsr_aceena,lwp_armbe_aceena,lwp_sat_aceena,lwp_m_aceena,lwp_m2_aceena],
#                           legend=['MFRSR','ARMBE','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_LWP_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ndrop_aceena,nd_wu_aceena,nd_sat_aceena,nd_m_aceena,nd_m2_aceena],
#                           legend=['Ndrop','Wu_etal','Nd_satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Nd_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([precip_aceena,precip_m_aceena,precip_m2_aceena],legend=['Obs','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Precip_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cld_arscl_aceena,cld_visst_aceena,cld_tsi_aceena,cld_m_aceena,cld_m2_aceena],
#                           legend=['ARSCL','Satellite','TSI','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_totcld_'+site+'_'+IOP+'.txt')


# calc.bias_corrcoef_RMSE(org_aceena,org_m_aceena,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_ORG_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(org_aceena,org_m2_aceena,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_ORG_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(so4_aceena, so4_m_aceena,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_SO4_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(so4_aceena, so4_m2_aceena,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_SO4_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn2_aceena,ccn2_m_aceena,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn2_aceena,ccn2_m2_aceena,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(cpc10_withmask_aceena,ncn10_m_aceena,label1='CPC_masked',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(cpc10_withmask_aceena,ncn10_m2_aceena,label1='CPC_masked',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(uhsas100_aceena, ncn100_m_aceena,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN100_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(uhsas100_aceena, ncn100_m2_aceena,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN100_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(lwp_armbe_aceena, lwp_m_aceena,label1='ARMBE',label2='E3SMv1', 
#                         outfile=figpath+'statistics_lwp_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(lwp_armbe_aceena, lwp_m2_aceena,label1='ARMBE',label2='E3SMv2', 
#                         outfile=figpath+'statistics_lwp_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ndrop_aceena, nd_m_aceena,label1='Ndrop',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ndrop_aceena, nd_m2_aceena,label1='Ndrop',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(nd_sat_aceena, nd_m_aceena,label1='Satellite',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsSat_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(nd_sat_aceena, nd_m2_aceena,label1='Satellite',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsSat_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(reff_aceena, reff_m_aceena,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(reff_aceena, reff_m2_aceena,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(reff_sat_aceena, reff_m_aceena,label1='Satellite',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsSat_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(reff_sat_aceena, reff_m2_aceena,label1='Satellite',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsSat_'+site+'_'+IOP+'.txt')

#%% joint histogram
fig,ax = plot.jointhist([uhsas100_aceena,ncn100_m_aceena,ncn100_m2_aceena], [ccn2_aceena,ccn2_m_aceena,ccn2_m2_aceena], 
                    xedges=np.arange(0,200,10),yedges=np.arange(0,200,10), normalize_x=True,
                    xlabel='CN (>100nm) (cm$^{-3}$)', ylabel='CCN (SS=0.2%) (cm$^{-3}$)', vmax=0.5,
                    title=['Ground','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CN100_CCN2_ship_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ccn2_aceena,ccn2_aceena,ccn2_aceena,ccn2_m_aceena,ccn2_m2_aceena], 
                        [ndrop_aceena,nd_wu_aceena,nd_sat_aceena,nd_m_aceena,nd_m2_aceena],
                    xedges=np.arange(0,200,20),yedges=np.arange(0,200,20), normalize_x=True,
                    xlabel='CCN (SS=0.2%) (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', vmax=0.4,
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_CCN2_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ndrop_aceena,nd_wu_aceena,nd_sat_aceena,nd_m_aceena,nd_m2_aceena],
                        [lwp_armbe_aceena,lwp_armbe_aceena,lwp_sat_aceena,lwp_m_aceena,lwp_m2_aceena], 
                    xedges=np.arange(0,200,20),yedges=np.arange(0,300,20), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', vmax=0.4,
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_LWP_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([ndrop_aceena,nd_wu_aceena,nd_sat_aceena,nd_m_aceena,nd_m2_aceena],
                        [reff_aceena,reff_wu_aceena,reff_sat_aceena,reff_m_aceena,reff_m2_aceena],
                    xedges=np.arange(0,200,20),yedges=np.arange(4,25,2), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', vmax=0.3,
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_Reff_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.jointhist([cod_sat_aceena,cod_sat_aceena,cod_m_aceena,cod_m2_aceena],[lwp_armbe_aceena,lwp_sat_aceena,lwp_m_aceena,lwp_m2_aceena], 
                    xedges=np.arange(0,40,3),yedges=np.arange(0,300,20), normalize_x=True,
                    xlabel='Cloud Optical Depth (N/A)', ylabel='LWP (g/m$^2$)', vmax=0.25,
                    title=['Ground','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'jointhist_COD_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter plot

fig,ax = plot.scatter([ndrop_aceena.data, nd_wu_aceena.data,nd_sat_aceena.data,nd_m_aceena.data,nd_m2_aceena.data], 
                      [ccn2_aceena.data,ccn2_aceena.data,ccn2_aceena.data,ccn2_m_aceena.data,ccn2_m2_aceena.data],
                      xlimit=(0,250), ylimit=(0,350),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=False)
fig.savefig(figpath+'scatter_Nd_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
fig,ax = plot.scatter([uhsas100_aceena.data,ncn100_m_aceena.data,ncn100_m2_aceena.data], 
                      [ccn2_aceena.data,ccn2_m_aceena.data,ccn2_m2_aceena.data],
                      xlimit=(0,300), ylimit=(0,300),
                    xlabel='Surface CN (>100nm) (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                    title=['Obs','E3SMv1','E3SMv2'],
                linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_CN100_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% heatmaps

# xedges=np.exp(np.arange(np.log(10),6.5,0.5))
# yedges=np.exp(np.arange(np.log(10),6.5,0.5))
fig,ax = plot.heatmap([ndrop_aceena.data, nd_wu_aceena.data,nd_sat_aceena.data,nd_m_aceena.data,nd_m2_aceena.data],
                      [lwp_armbe_aceena.data,lwp_armbe_aceena.data,lwp_sat_aceena.data,lwp_m_aceena.data,lwp_m2_aceena.data],
                      [albedo_aceena,albedo_aceena,albedo_aceena,albedo_m_aceena,albedo_m2_aceena],vmax=60,
                    xedges=np.arange(0,250,20), yedges=np.arange(0,250,20),
                    # xedges=xedges, yedges=yedges, 
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', zlabel='TOA Albedo (%)',
                    title=['Ndrop','Wu_etal','Satellite','E3SMv1','E3SMv2'])
fig.savefig(figpath+'heatmap_Albedo_vs_Nd_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
