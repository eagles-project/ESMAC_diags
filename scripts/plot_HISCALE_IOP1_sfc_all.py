"""
script to generate all plots for HISCALE surface data

"""
import os
import glob
import yaml
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

config_file = '../config/config.yml'
stream = open(config_file, "r")
config = yaml.full_load(stream)

prep_model_path = '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/prep_data/'+site+'/model/sfc_prof/3600s/'
prep_sfc_path = '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/prep_data/'+site+'/surface/300s/'
prep_sat_path = '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/prep_data/'+site+'/satellite/3600s/'

time_hiscale = pd.date_range(start='2016-04-25', end='2016-05-21', freq="3600s")
IOP = 'IOP1'
            
# path of output figures
figpath= '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/figures/'+site+'/sfc_toa/'
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
org_hiscale = org.sel(time=time_hiscale)
so4_hiscale = so4.sel(time=time_hiscale)
nh4_hiscale = nh4.sel(time=time_hiscale)
no3_hiscale = no3.sel(time=time_hiscale)
chl_hiscale = chl.sel(time=time_hiscale)

obsdata = xr.open_dataset(prep_sfc_path + 'totcld_'+site+'.nc')
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
# cld_visst = obsdata['tot_cld_visst'].load()
obsdata.close()
cld_arscl_hiscale = cld_arscl.sel(time=time_hiscale)
cld_tsi_hiscale = cld_tsi.sel(time=time_hiscale)
# cld_visst_hiscale = cld_visst.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'cloud_2d_'+site+'.nc')
height_o = obsdata['height'].load()
cloud_2d = obsdata['cloud'].load()
obsdata.close()
cloud_2d_hiscale = cloud_2d.sel(time=time_hiscale)
obsdata = xr.open_mfdataset(sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'_*.nc')))
ccn2 = obsdata['CCN2'].load()
obsdata.close()
ccn2_hiscale = ccn2.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'sfc_CPC_'+site+'.nc')
cpc10 = obsdata['cpc10'].load()
cpc3 = obsdata['cpc3'].load()
obsdata.close()
cpc3_hiscale = cpc3.sel(time3=time_hiscale)
cpc10_hiscale = cpc10.sel(time10=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'cod_'+site+'.nc')
cod = obsdata['cod'].load()
obsdata.close()
cod_hiscale = cod.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'LWP_'+site+'.nc')
lwp = obsdata['lwp'].load()
obsdata.close()
lwp_hiscale = lwp.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'Ndrop_'+site+'.nc')
ndrop = obsdata['cdnc'].load()
obsdata.close()
ndrop_hiscale = ndrop.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'reff_'+site+'.nc')
reff = obsdata['reff'].load()
obsdata.close()
reff_hiscale = reff.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sfc_path + 'precip_'+site+'.nc')
precip_tbrg = obsdata['precip_tbrg'].load()
precip_pwd = obsdata['precip_pwd'].load()
precip_pars = obsdata['precip_pars'].load()
obsdata.close()
precip_tbrg_hiscale = precip_tbrg.sel(time=time_hiscale)
precip_pwd_hiscale = precip_pwd.sel(time=time_hiscale)
precip_pars_hiscale = precip_pars.sel(time=time_hiscale)

obsdata = xr.open_dataset(prep_sfc_path + 'sfc_radiation_'+site+'.nc')
lwdnsfc = obsdata['lwdn'].load()
swdnsfc = obsdata['swdn'].load()
lwupsfc = obsdata['lwup'].load()
swupsfc = obsdata['swup'].load()
obsdata.close()
lwdnsfc_hiscale = lwdnsfc.sel(time=time_hiscale)
swdnsfc_hiscale = swdnsfc.sel(time=time_hiscale)
lwupsfc_hiscale = lwupsfc.sel(time=time_hiscale)
swupsfc_hiscale = swupsfc.sel(time=time_hiscale)
lwnetsfc_hiscale = lwupsfc_hiscale - lwdnsfc_hiscale
swnetsfc_hiscale = swdnsfc_hiscale - swupsfc_hiscale

obsdata = xr.open_dataset(prep_sfc_path + 'sfc_UHSAS_'+site+'.nc')
size_uhsas = obsdata['size'].load()
dmin_hiscale = obsdata['size_low'].load()
dmax_hiscale = obsdata['size_high'].load()
uhsas100 = obsdata['uhsas100'].load()
uhsas_all = obsdata['uhsas_all'].load()
obsdata.close()
uhsas100_hiscale = uhsas100.sel(time=time_hiscale)
uhsasall_hiscale = uhsas_all.sel(time=time_hiscale)
dlogDp_uhsas = np.mean(np.log10(dmax_hiscale/dmin_hiscale))

obsdata = xr.open_dataset(prep_sfc_path + 'sfc_SMPS_'+site+'_IOP1.nc')
time1 = obsdata['time'].load()
smps100_1 = obsdata['smps100_dlogDp'].load()
smpsall_1 = obsdata['dN_dlogDp'].load()
size1 = obsdata['size'].load()
obsdata.close()
obsdata = xr.open_dataset(prep_sfc_path + 'sfc_SMPS_'+site+'_IOP2.nc')
time2 = obsdata['time'].load()
smps100_2 = obsdata['smps100_dlogDp'].load()
smpsall_2 = obsdata['dN_dlogDp'].load()
size2 = obsdata['size'].load()
obsdata.close()
time_smps = xr.concat((time1,time2),dim='time')
size_smps = size1
smps100 = xr.concat((smps100_1,smps100_2),dim='time')
smps_all = xr.concat((smpsall_1,smpsall_2.interp(size=smpsall_1['size'])),dim='time')
# SMPS data is already dN/dlogDp, total number concentration must multiply by dlogDp
dlogDp_smps = np.mean(np.log10(size_smps[1:].data/size_smps[0:-1].data))
smps100 = smps100 * dlogDp_smps
smps_all = smps_all * dlogDp_smps
smps100_hiscale = smps100.sel(time=time_hiscale)
smpsall_hiscale = smps_all.sel(time=time_hiscale)

obsdata = xr.open_mfdataset(prep_sfc_path + 'cloudheight_ARSCL_'+site+'.nc')
cth = obsdata['cth'].load()
cbh = obsdata['cbh'].load()
cths = obsdata['cths'].load()
obsdata.close()
cth_hiscale = cth.sel(time=time_hiscale)
cbh_hiscale = cbh.sel(time=time_hiscale)

# satellite data
obsdata = xr.open_dataset(prep_sat_path + 'albedo_VISSTgrid_'+site+'.nc')
albedo_sat = obsdata['albedo'].load()
obsdata.close()
albedo_hiscale = albedo_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'cloudfraction_VISSTgrid_'+site+'.nc')
cfall_sat = obsdata['cldtot'].load()
cflow_sat = obsdata['cldlow'].load()
obsdata.close()
cfall_sat_hiscale = cfall_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc')
ctt_sat = obsdata['ctt'].load()
cth_sat = obsdata['cth'].load()
obsdata.close()
ctt_sat_hiscale = ctt_sat.sel(time=time_hiscale)
cth_sat_hiscale = cth_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'cod_VISSTgrid_'+site+'.nc')
cod_sat = obsdata['cod'].load()
obsdata.close()
cod_sat_hiscale = cod_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc')
reff_sat = obsdata['reff'].load()
obsdata.close()
reff_sat_hiscale = reff_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc')
lwp_sat = obsdata['lwp'].load()
obsdata.close()
lwp_sat_hiscale = lwp_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc')
nd_sat = obsdata['Nd'].load()
obsdata.close()
nd_sat_hiscale = nd_sat.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc')
lwnettoa = obsdata['lwnettoa'].load()
obsdata.close()
lwnettoa_hiscale = lwnettoa.sel(time=time_hiscale)
obsdata = xr.open_dataset(prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc')
swnettoa = obsdata['swnettoa'].load()
obsdata.close()
swnettoa_hiscale = swnettoa.sel(time=time_hiscale)

# E3SM data
filename = prep_model_path +site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
if config['aerosol_output'] == True:
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
if config['tau3d_output'] == True:
            cod_m = modeldata['cod'].load()
if config['reff_output'] == True:
            reff_m = modeldata['reff'].load()
lwp_m = modeldata[config['LWP']].load()
nd_m = modeldata['Nd_mean'].load()
if config['convectiveparam'] == True:
            precip_m = modeldata[config['PRECIPSFCTOT']].load()
else:
            precip_m = modeldata[config['PRECIPSFCLIQ']].load()
cld_m = modeldata[config['CLDTOT']].load()
cbh_m = modeldata['cbh'].load()
cth_m = modeldata['cth'].load()
Hcld_m = modeldata['clddepth'].load()
cldlow_m = modeldata[config['CLDLOW']].load()
cldmid_m = modeldata[config['CLDMED']].load()
cldhgh_m = modeldata[config['CLDHGH']].load()
if config['netradiation_output'] == True:
            lwnetsfc_m = modeldata['FLNS'].load()
            lwnettoa_m = modeldata['FLNT'].load()
            swnetsfc_m = modeldata['FSNS'].load()
            swnettoa_m = modeldata['FSNT'].load()
else:
            lwupsfc_m = modeldata['LWUPSFC'].load()
            swupsfc_m = modeldata['SWUPSFC'].load()
            # lwdntoa_m = 
lwdnsfc_m = modeldata[config['LWDOWNSFC']].load()
lwuptoa_m = modeldata[config['LWUPTOA']].load()
swdnsfc_m = modeldata[config['SWDOWNSFC']].load()            
swdntoa_m = modeldata[config['SWDOWNTOA']].load()
swuptoa_m = modeldata[config['SWUPTOA']].load()
modeldata.close()

if config['netradiation_output'] == True:
            lwupsfc_m = lwnetsfc_m + lwdnsfc_m
            swupsfc_m = swdnsfc_m - swnetsfc_m
else:
            lwnetsfc_m = lwupsfc_m - lwdnsfc_m
            # lwnettoa_m = lwdntoa_m - lwuptoa_m
            swnetsfc_m = swdnsfc_m - swupsfc_m
            swnettoa_m = swdntoa_m - swuptoa_m
albedo_m = swuptoa_m/swdntoa_m*100
if config['aerosol_output'] == True:
            org_m = pom_m + mom_m + soa_m
            bc_m_hiscale = bc_m.sel(time=time_hiscale)
            dst_m_hiscale = dst_m.sel(time=time_hiscale)
            org_m_hiscale = org_m.sel(time=time_hiscale)
            so4_m_hiscale = so4_m.sel(time=time_hiscale)
            ncl_m_hiscale = ncl_m.sel(time=time_hiscale)
            ccn2_m_hiscale = ccn2_m.sel(time=time_hiscale)
            ncn3_m_hiscale = ncn3_m.sel(time=time_hiscale)
            ncn10_m_hiscale = ncn10_m.sel(time=time_hiscale)
            ncn100_m_hiscale = ncn100_m.sel(time=time_hiscale)
            CNsize_m_hiscale = CNsize_m.sel(time=time_hiscale)
if config['tau3d_output'] == True:
            cod_m_hiscale = cod_m.sel(time=time_hiscale)
if config['reff_output'] == True:
            reff_m_hiscale = reff_m.sel(time=time_hiscale)
lwp_m_hiscale = lwp_m.sel(time=time_hiscale)
nd_m_hiscale = nd_m.sel(time=time_hiscale)
precip_m_hiscale = precip_m.sel(time=time_hiscale)
cld_m_hiscale = cld_m.sel(time=time_hiscale)
cbh_m_hiscale = cbh_m.sel(time=time_hiscale)
cth_m_hiscale = cth_m.sel(time=time_hiscale)
Hcld_m_hiscale = Hcld_m.sel(time=time_hiscale)
cldlow_m_hiscale = cldlow_m.sel(time=time_hiscale)
cldmid_m_hiscale = cldmid_m.sel(time=time_hiscale)
cldhgh_m_hiscale = cldhgh_m.sel(time=time_hiscale)
lwnetsfc_m_hiscale = lwnetsfc_m.sel(time=time_hiscale)
# lwnettoa_m_hiscale = lwnettoa_m.sel(time=time_hiscale)
swnetsfc_m_hiscale = swnetsfc_m.sel(time=time_hiscale)
swnettoa_m_hiscale = swnettoa_m.sel(time=time_hiscale)
albedo_m_hiscale = albedo_m.sel(time=time_hiscale)

# filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
# modeldata = xr.open_dataset(filename)
# time_m2 = modeldata['time'].load()
# bc_m2 = modeldata['bc'].load()
# dst_m2 = modeldata['dst'].load()
# mom_m2 = modeldata['mom'].load()
# pom_m2 = modeldata['pom'].load()
# ncl_m2 = modeldata['ncl'].load()
# so4_m2 = modeldata['so4'].load()
# soa_m2 = modeldata['soa'].load()
# ccn2_m2 = modeldata['CCN4'].load()
# ncn3_m2 = modeldata['NCN3'].load()
# ncn10_m2 = modeldata['NCN10'].load()
# ncn100_m2 = modeldata['NCN100'].load()
# CNsize_m2 = modeldata['NCNall'].load()
# cod_m2 = modeldata['cod'].load()
# reff_m2 = modeldata['reff'].load()
# lwp_m2 = modeldata['TGCLDLWP'].load()
# nd_m2 = modeldata['Nd_mean'].load()
# precip_m2 = modeldata['PRECT'].load()
# cld_m2 = modeldata['CLDTOT'].load()
# cbh_m2 = modeldata['cbh'].load()
# cth_m2 = modeldata['cth'].load()
# Hcld_m2 = modeldata['clddepth'].load()
# cldlow_m2 = modeldata['CLDLOW'].load()
# cldmid_m2 = modeldata['CLDMED'].load()
# cldhgh_m2 = modeldata['CLDHGH'].load()
# lwdnsfc_m2 = modeldata['FLDS'].load()
# lwnetsfc_m2 = modeldata['FLNS'].load()
# lwnettoa_m2 = modeldata['FLNT'].load()
# lwuptoa_m2 = modeldata['FLUT'].load()
# swdnsfc_m2 = modeldata['FSDS'].load()
# swnetsfc_m2 = modeldata['FSNS'].load()
# swdntoa_m2 = modeldata['SOLIN'].load()
# swnettoa_m2 = modeldata['FSNT'].load()
# swuptoa_m2 = modeldata['FSUTOA'].load()
# modeldata.close()
# lwupsfc_m2 = lwnetsfc_m2 + lwdnsfc_m2
# swupsfc_m2 = swdnsfc_m2 - swnetsfc_m2
# albedo_m2 = swuptoa_m2/swdntoa_m2*100
# org_m2 = pom_m2 + mom_m2 + soa_m2
# bc_m2_hiscale = bc_m2.sel(time=time_hiscale)
# dst_m2_hiscale = dst_m2.sel(time=time_hiscale)
# org_m2_hiscale = org_m2.sel(time=time_hiscale)
# so4_m2_hiscale = so4_m2.sel(time=time_hiscale)
# ncl_m2_hiscale = ncl_m2.sel(time=time_hiscale)
# ccn2_m2_hiscale = ccn2_m2.sel(time=time_hiscale)
# ncn3_m2_hiscale = ncn3_m2.sel(time=time_hiscale)
# ncn10_m2_hiscale = ncn10_m2.sel(time=time_hiscale)
# ncn100_m2_hiscale = ncn100_m2.sel(time=time_hiscale)
# CNsize_m2_hiscale = CNsize_m2.sel(time=time_hiscale)
# cod_m2_hiscale = cod_m2.sel(time=time_hiscale)
# reff_m2_hiscale = reff_m2.sel(time=time_hiscale)
# lwp_m2_hiscale = lwp_m2.sel(time=time_hiscale)
# nd_m2_hiscale = nd_m2.sel(time=time_hiscale)
# precip_m2_hiscale = precip_m2.sel(time=time_hiscale)
# cld_m2_hiscale = cld_m2.sel(time=time_hiscale)
# cbh_m2_hiscale = cbh_m2.sel(time=time_hiscale)
# cth_m2_hiscale = cth_m2.sel(time=time_hiscale)
# Hcld_m2_hiscale = Hcld_m2.sel(time=time_hiscale)
# cldlow_m2_hiscale = cldlow_m2.sel(time=time_hiscale)
# cldmid_m2_hiscale = cldmid_m2.sel(time=time_hiscale)
# cldhgh_m2_hiscale = cldhgh_m2.sel(time=time_hiscale)
# lwnetsfc_m2_hiscale = lwnetsfc_m2.sel(time=time_hiscale)
# lwnettoa_m2_hiscale = lwnettoa_m2.sel(time=time_hiscale)
# swnetsfc_m2_hiscale = swnetsfc_m2.sel(time=time_hiscale)
# swnettoa_m2_hiscale = swnettoa_m2.sel(time=time_hiscale)
# albedo_m2_hiscale = albedo_m2.sel(time=time_hiscale)

filename = prep_model_path +site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
height_m = modeldata['height'].load()
cf_e3sm = modeldata['cloud_z'].load()
modeldata.close()
cloud_m_hiscale = cf_e3sm.sel(time=time_hiscale)

# filename = prep_model_path + 'E3SMv2_'+site+'_profiles.nc'
# modeldata = xr.open_dataset(filename)
# height_m2 = modeldata['height'].load()
# cf_e3sm2 = modeldata['cloud_z'].load()
# modeldata.close()
# cloud_m2_hiscale = cf_e3sm2.sel(time=time_hiscale)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# divided by dlogDp in size distribution
if config['aerosol_output'] == True:
            dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
            CNsize_m_hiscale = CNsize_m_hiscale.T/dlogDp_e3sm
            # CNsize_m2_hiscale = CNsize_m2_hiscale.T/dlogDp_e3sm
            pdf_m_hiscale = np.nanmean(CNsize_m_hiscale,axis=0)
            # pdf_m2_hiscale = np.nanmean(CNsize_m2_hiscale,axis=0)
smpsall_hiscale = smpsall_hiscale / dlogDp_smps
uhsasall_hiscale = uhsasall_hiscale / dlogDp_uhsas
pdf_uhsas_hiscale = np.nanmean(uhsasall_hiscale,axis=0)
pdf_smps_hiscale = np.nanmean(smpsall_hiscale,axis=0)

ndrop_hiscale[ndrop_hiscale<10] = np.nan
nd_sat_hiscale[nd_sat_hiscale<10] = np.nan
nd_m_hiscale[nd_m_hiscale<10] = np.nan
# nd_m2_hiscale[nd_m2_hiscale<10] = np.nan
ndrop_hiscale[ndrop_hiscale>500] = np.nan
nd_sat_hiscale[nd_sat_hiscale>500] = np.nan
nd_m_hiscale[nd_m_hiscale>500] = np.nan
# nd_m2_hiscale[nd_m2_hiscale>500] = np.nan

lwp_hiscale[lwp_hiscale<20] = np.nan
lwp_sat_hiscale[lwp_sat_hiscale<20] = np.nan
lwp_m_hiscale[lwp_m_hiscale<20] = np.nan
# lwp_m2_hiscale[lwp_m2_hiscale<20] = np.nan

cod_hiscale[cod_hiscale<2] = np.nan
cod_sat_hiscale[cod_sat_hiscale<2] = np.nan
cod_hiscale[cod_hiscale>100] = np.nan
cod_sat_hiscale[cod_sat_hiscale>100] = np.nan
if config['tau3d_output'] == True:
            cod_m_hiscale[cod_m_hiscale<2] = np.nan
            # cod_m2_hiscale[cod_m2_hiscale<2] = np.nan
            cod_m_hiscale[cod_m_hiscale>100] = np.nan
            # cod_m2_hiscale[cod_m2_hiscale>100] = np.nan

# unit change:
precip_m_hiscale = precip_m_hiscale*3600*1000   # m/s to mm/hr
# precip_m2_hiscale = precip_m2_hiscale*3600*1000   # m/s to mm/hr
cloud_m_hiscale = cloud_m_hiscale*100  # fraction to %
# cloud_m2_hiscale = cloud_m2_hiscale*100  # fraction to %
height_o = height_o.data*0.001   # m to km
height_m = height_m.data*0.001   # m to km
# height_m2 = height_m2.data*0.001   # m to km

# set a small threshold of E3SM precipitation
precip_m_hiscale[precip_m_hiscale<0.02] = 0
# precip_m2_hiscale[precip_m2_hiscale<0.02] = 0


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if config['aerosol_output'] == True:
            #%% bar plot
            datagroup0 = [org_hiscale,so4_hiscale,nh4_hiscale,no3_hiscale,chl_hiscale, [], []]
            datagroup1 = [org_m_hiscale, so4_m_hiscale, [], [], [], bc_m_hiscale, dst_m_hiscale]
            # datagroup2 = [org_m2_hiscale, so4_m2_hiscale, [], [], [], bc_m2_hiscale, dst_m2_hiscale]
            # dataall=[datagroup0,datagroup1, datagroup2,]
            dataall=[datagroup0,datagroup1,]
            labelall = ['Organic', 'SO$_4$', 'NH$_4$', 'NO$_3$', 'Chl', 'BC', 'Dust']
            colorall = ['limegreen', 'red', 'lightblue', 'orange', 'cyan', 'k', 'silver']
            # fig,ax = plot.bar(dataall, datalabel=['Obs','E3SMv1','E3SMv2',], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
            fig,ax = plot.bar(dataall, datalabel=['Obs','Model',], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
                              title='Aerosol Composition  '+site+' '+IOP, varlabel= labelall, colorall=colorall)
            fig.savefig(figpath+'bar_composition_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            #%% timeseries
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [org_hiscale,org_m_hiscale,org_m2_hiscale], 
            #                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'],
            fig,ax = plot.timeseries([time_hiscale,time_hiscale], [org_hiscale,org_m_hiscale], 
                                      legend = ['Obs','Model'], color=['k','r'], 
                                      title='Total Organic '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [so4_hiscale,so4_m_hiscale,so4_m2_hiscale], 
            #                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'],  
            fig,ax = plot.timeseries([time_hiscale,time_hiscale], [so4_hiscale,so4_m_hiscale], 
                                      legend = ['Obs','Model'], color=['k','r','b'],  
                                      title='Sulfate '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], 
            #                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.timeseries([time_hiscale,time_hiscale], [ccn2_hiscale,ccn2_m_hiscale], 
                                      legend = ['Obs','Model'], color=['k','r'], 
                                    title='0.2%CCN '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [cpc3_hiscale,ncn3_m_hiscale,ncn3_m2_hiscale], 
            #                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.timeseries([time_hiscale,time_hiscale], [cpc3_hiscale,ncn3_m_hiscale], 
                                      legend = ['Obs','Model'], color=['k','r'], 
                                      title='CN(>3nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_CPC3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale], 
            #                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.timeseries([time_hiscale,time_hiscale], [cpc10_hiscale,ncn10_m_hiscale], 
                                      legend = ['Obs','Model'], color=['k','r'], 
                                      title='CN(>10nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [smps100_hiscale,uhsas100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], 
            #                         legend = ['SMPS','UHSAS','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
            fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [smps100_hiscale,uhsas100_hiscale,ncn100_m_hiscale], 
                                    legend = ['SMPS','UHSAS','Model'], color=['k','gray','r'],
                                    title='CN(>100nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['tau3d_output'] == True:
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [cod_hiscale,cod_sat_hiscale,cod_m_hiscale,cod_m2_hiscale], 
            #                           legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], #marker='.',
            fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [cod_hiscale,cod_sat_hiscale,cod_m_hiscale], 
                                      legend = ['MFRSR','Satellite','Model'], color=['k','gray','r'], #marker='.',
                                    title='cloud optical depth '+site+' '+IOP, xlabel=None, ylabel=None)
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [lwp_hiscale, lwp_sat_hiscale, lwp_m_hiscale, lwp_m2_hiscale], 
#                         legend = ['MWR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [lwp_hiscale, lwp_sat_hiscale, lwp_m_hiscale], 
                        legend = ['MWR','Satellite','Model'], color=['k','gray','r'],
                        title='LWP '+site+' '+IOP,xlabel=None, ylabel="g/m$^2$")
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [ndrop_hiscale,nd_sat_hiscale, nd_m_hiscale, nd_m2_hiscale], 
#                           legend = ['Ndrop','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], #marker='.',
fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [ndrop_hiscale,nd_sat_hiscale, nd_m_hiscale], 
                          legend = ['Ndrop','Satellite','Model'], color=['k','gray','r'], #marker='.',
                          title='Nd '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['reff_output'] == True:
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale],  
            #                         legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],marker='.',
            fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [reff_hiscale,reff_sat_hiscale,reff_m_hiscale],  
                                    legend = ['MFRSR','Satellite','Model'], color=['k','gray','r'],marker='.',
                                    title='Reff '+site+' '+IOP,xlabel=None, ylabel='$\mu$m')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale,time_hiscale], [precip_tbrg_hiscale,precip_pwd_hiscale,precip_pars_hiscale,precip_m_hiscale,precip_m2_hiscale],  
#                           legend = ['Tipping Bucket', 'PWD', 'Disdrometer', 'E3SMv1','E3SMv2'], color=['k','gray','silver','r','b'], 
fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [precip_tbrg_hiscale,precip_pwd_hiscale,precip_pars_hiscale,precip_m_hiscale],  
                          legend = ['Tipping Bucket', 'PWD', 'Disdrometer', 'Model'], color=['k','gray','silver','r'], 
                        title='Precip '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='mm/hr')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [lwnetsfc_hiscale,lwnetsfc_m_hiscale,lwnetsfc_m2_hiscale], 
#                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
fig,ax = plot.timeseries([time_hiscale,time_hiscale], [lwnetsfc_hiscale,lwnetsfc_m_hiscale], 
                          legend = ['Obs','Model'], color=['k','r'], 
                        title='Sfc. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_LWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [swnetsfc_hiscale,swnetsfc_m_hiscale,swnetsfc_m2_hiscale], 
#                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
fig,ax = plot.timeseries([time_hiscale,time_hiscale], [swnetsfc_hiscale,swnetsfc_m_hiscale], 
                          legend = ['Obs','Model'], color=['k','r'], 
                        title='Sfc. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_SWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['netradiation_output'] == True:
            # fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [lwnettoa_hiscale,lwnettoa_m_hiscale,lwnettoa_m2_hiscale], 
            #                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.timeseries([time_hiscale,time_hiscale], [lwnettoa_hiscale,lwnettoa_m_hiscale], 
                                      legend = ['Obs','Model'], color=['k','r'], 
                                    title='TOA. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'timeseries_LWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [swnettoa_hiscale,swnettoa_m_hiscale,swnettoa_m2_hiscale], 
#                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
fig,ax = plot.timeseries([time_hiscale,time_hiscale], [swnettoa_hiscale,swnettoa_m_hiscale], 
                          legend = ['Obs','Model'], color=['k','r'], 
                        title='TOA. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_SWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale,time_hiscale], [cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale,cld_m2_hiscale], 
#                         legend = ['ARSCL','TSI','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.timeseries([time_hiscale,time_hiscale,time_hiscale], [cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale], 
                        legend = ['ARSCL','TSI','Model'], color=['k','gray','r'],
                        title='Cloud fraction '+site+' '+IOP,xlabel=None, ylabel="%")
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%%
if config['aerosol_output'] == True:
            # fig,ax = plot.timeseries_size([time_hiscale,time_hiscale,time_hiscale,time_hiscale], 
            #                               [size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
            #                               [uhsasall_hiscale.T.data, smpsall_hiscale.T.data, CNsize_m_hiscale.T.data, CNsize_m2_hiscale.T.data], 
            #                               legend = ['UHSAS','SMPS','E3SMv1','E3SMv2'],
            fig,ax = plot.timeseries_size([time_hiscale,time_hiscale,time_hiscale], 
                                          [size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
                                          [uhsasall_hiscale.T.data, smpsall_hiscale.T.data, CNsize_m_hiscale.T.data], 
                                          legend = ['UHSAS','SMPS','Model'],
                                      ylabel='Diameter (nm)', ylimit=(3,1000),
                                      title = 'Aerosol Size Distribution (dN/dlogDp, cm$^{-3}$)')
            # for ax_i in ax:
            #     ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            #     ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            fig.savefig(figpath+'aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries_2d([time_hiscale,time_hiscale,time_hiscale], 
#                             [height_o, height_m, height_m2], 
#                             [cloud_2d_hiscale.T.data, cloud_m_hiscale.T.data, cloud_m2_hiscale.T.data], 
#                               legend = ['Obs','E3SMv1','E3SMv2'], title = 'Cloud Fraction (%)'),
fig,ax = plot.timeseries_2d([time_hiscale,time_hiscale], [height_o, height_m], [cloud_2d_hiscale.T.data, cloud_m_hiscale.T.data],
                            legend = ['Obs','Model'], title = 'Cloud Fraction (%)', yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (m)',cmap='jet') #, ylimit=(3,1000)
# for ax_i in ax:
#     ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#     ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'cloud_2d_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


#%% diurnal cycle
if config['aerosol_output'] == True:
            # fig,ax = plot.diurnalcycle([org_hiscale,org_m_hiscale,org_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.diurnalcycle([org_hiscale,org_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                                      title='Organic '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
            fig.savefig(figpath+'diurnalcycle_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.diurnalcycle([so4_hiscale,so4_m_hiscale,so4_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.diurnalcycle([so4_hiscale,so4_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                                      title='Sulfate '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
            fig.savefig(figpath+'diurnalcycle_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.diurnalcycle([ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.diurnalcycle([ccn2_hiscale,ccn2_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                                    title='0.2%CCN '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
            fig.savefig(figpath+'diurnalcycle_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.diurnalcycle([cpc3_hiscale,ncn3_m_hiscale,ncn3_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.diurnalcycle([cpc3_hiscale,ncn3_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                                    title='CN(>3nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
            fig.savefig(figpath+'diurnalcycle_CPC3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.diurnalcycle([cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.diurnalcycle([cpc10_hiscale,ncn10_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                                    title='CN(>10nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
            fig.savefig(figpath+'diurnalcycle_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.diurnalcycle([smps100_hiscale,uhsas100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], legend = ['SMPS100','UHSAS100','E3SMv1','E3SMv2'], 
                                    # title='CN(>100nm) '+site+' '+IOP, color=['k','gray','r','b'], xlabel='Time (UTC)',ylabel='cm$^{-3}$')
            fig,ax = plot.diurnalcycle([smps100_hiscale,uhsas100_hiscale,ncn100_m_hiscale], legend = ['SMPS100','UHSAS100','Model'], 
                                    title='CN(>100nm) '+site+' '+IOP, color=['k','gray','r'], xlabel='Time (UTC)',ylabel='cm$^{-3}$')
            fig.savefig(figpath+'diurnalcycle_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['tau3d_output'] == True:
            # fig,ax = plot.diurnalcycle( [cod_hiscale, cod_sat_hiscale, cod_m_hiscale, cod_m2_hiscale], 
            #                             legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'], 
            fig,ax = plot.diurnalcycle( [cod_hiscale, cod_sat_hiscale, cod_m_hiscale], 
                                        legend = ['MFRSR','Satellite','Model'], color=['k','gray','r'], 
                                    title='Cloud optical depth '+site+' '+IOP, xlabel='Time (UTC)', ylabel=None)
            fig.savefig(figpath+'diurnalcycle_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([lwp_hiscale,lwp_sat_hiscale, lwp_m_hiscale, lwp_m2_hiscale], 
#                             legend = ['MWR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.diurnalcycle([lwp_hiscale,lwp_sat_hiscale, lwp_m_hiscale], 
                            legend = ['MWR','Satellite','Model'], color=['k','gray','r'],
                        title='LWP '+site+' '+IOP,  xlabel='Time (UTC)',ylabel="g/m$^2$")
fig.savefig(figpath+'diurnalcycle_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([ndrop_hiscale, nd_sat_hiscale, nd_m_hiscale,nd_m2_hiscale], 
#                             legend = ['Ndrop', 'Satellite', 'E3SMv1', 'E3SMv2'], color=['k','gray','r','b'], 
fig,ax = plot.diurnalcycle([ndrop_hiscale, nd_sat_hiscale, nd_m_hiscale], 
                            legend = ['Ndrop', 'Satellite', 'Model'], color=['k','gray','r'], 
                          title='Nd '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
fig.savefig(figpath+'diurnalcycle_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['reff_output'] == True:
            # fig,ax = plot.diurnalcycle([reff_hiscale, reff_sat_hiscale, reff_m_hiscale, reff_m2_hiscale], 
            #                             legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
            fig,ax = plot.diurnalcycle([reff_hiscale, reff_sat_hiscale, reff_m_hiscale], 
                                        legend = ['MFRSR','Satellite','Model'], color=['k','gray','r'],
                                    title='droplet effective radius '+site+' '+IOP, xlabel='Time (UTC)', ylabel='$\mu$m')
            fig.savefig(figpath+'diurnalcycle_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle( [precip_tbrg_hiscale,precip_pwd_hiscale,precip_pars_hiscale,precip_m_hiscale,precip_m2_hiscale], 
#                             legend = ['Tipping Bucket','PWD','Disdrometer','E3SMv1','E3SMv2'], color=['k','gray','silver','r','b'], 
fig,ax = plot.diurnalcycle( [precip_tbrg_hiscale,precip_pwd_hiscale,precip_pars_hiscale,precip_m_hiscale], 
                            legend = ['Tipping Bucket','PWD','Disdrometer','Model'], color=['k','gray','silver','r'], 
                        nozero_percentile=True, title='Precipitation '+site+' '+IOP, xlabel='Time (UTC)',ylabel='mm/hr')
fig.savefig(figpath+'diurnalcycle_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([lwnetsfc_hiscale,lwnetsfc_m_hiscale,lwnetsfc_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
fig,ax = plot.diurnalcycle([lwnetsfc_hiscale,lwnetsfc_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                        title='Sfc. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)',ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_LWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([swnetsfc_hiscale,swnetsfc_m_hiscale,swnetsfc_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
fig,ax = plot.diurnalcycle([swnetsfc_hiscale,swnetsfc_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                        title='Sfc. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_SWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['netradiation_output'] == True:
            # fig,ax = plot.diurnalcycle([lwnettoa_hiscale,lwnettoa_m_hiscale,lwnettoa_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
            fig,ax = plot.diurnalcycle([lwnettoa_hiscale,lwnettoa_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                                    title='TOA. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
            fig.savefig(figpath+'diurnalcycle_LWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([swnettoa_hiscale,swnettoa_m_hiscale,swnettoa_m2_hiscale], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
fig,ax = plot.diurnalcycle([swnettoa_hiscale,swnettoa_m_hiscale], legend = ['Obs','Model'], color=['k','r'], 
                        title='TOA. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
fig.savefig(figpath+'diurnalcycle_SWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale,cld_m2_hiscale],
#                             legend = ['ARSCL','TSI','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.diurnalcycle([cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale],
                            legend = ['ARSCL','TSI','Model'], color=['k','gray','r'],
                            title='Total cloud fraction '+site+' '+IOP, xlabel='Time (UTC)', ylabel="%")
fig.savefig(figpath+'diurnalcycle_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['aerosol_output'] == True:
            # fig,ax = plot.diurnalcycle_2d([uhsasall_hiscale.T, smpsall_hiscale.T, CNsize_m_hiscale.T, CNsize_m2_hiscale.T], 
            #                               y=[size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
            #                               title= ['UHSAS','SMPS','E3SMv1','E3SMv2'],
            fig,ax = plot.diurnalcycle_2d([uhsasall_hiscale.T, smpsall_hiscale.T, CNsize_m_hiscale.T], 
                                          y=[size_uhsas,size_smps, np.arange(1,3001)], 
                                          title= ['UHSAS','SMPS','Model'],
                                          levellist=np.arange(0,11500,200), xlabel='Time (UTC)', ylabel='Diameter (nm)', 
                                          ylimit=(3,1000),cmap='jet')
            for ax_i in ax:
                ax_i.set_yscale('log')
            fig.savefig(figpath+'diurnalcycle_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle_2d([cloud_2d_hiscale, cloud_m_hiscale, cloud_m2_hiscale], 
#                               y = [height_o, height_m, height_m2],
#                           title= ['Obs', 'E3SMv1', 'E3SMv2',]),
fig,ax = plot.diurnalcycle_2d([cloud_2d_hiscale, cloud_m_hiscale], y = [height_o, height_m], title= ['Obs', 'Model',],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)',  cmap='jet', levellist=np.arange(0,45,1))
fig.savefig(figpath+'diurnalcycle_cloud2d_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% 1d histogram
if config['aerosol_output'] == True:
            w1 = np.ones_like(org_hiscale)/sum(~np.isnan(org_hiscale.data))
            w2 = np.ones_like(org_m_hiscale)/sum(~np.isnan(org_m_hiscale.data))
            # w3 = np.ones_like(org_m2_hiscale)/sum(~np.isnan(org_m2_hiscale.data))
            # fig,ax = plot.hist([org_hiscale,org_m_hiscale,org_m2_hiscale], weights=[w1,w2,w3], bins=np.arange(0,10,0.3),
            #                     legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
            fig,ax = plot.hist([org_hiscale,org_m_hiscale], weights=[w1,w2], bins=np.arange(0,10,0.3),
                                legend =['Obs','Model',], color=['k','r'],
                                title = 'Total Organic '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
            fig.savefig(figpath+'hist_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            w1 = np.ones_like(so4_hiscale)/sum(~np.isnan(so4_hiscale.data))
            w2 = np.ones_like(so4_m_hiscale)/sum(~np.isnan(so4_m_hiscale.data))
            # w3 = np.ones_like(so4_m2_hiscale)/sum(~np.isnan(so4_m2_hiscale.data))
            # fig,ax = plot.hist([so4_hiscale,so4_m_hiscale,so4_m2_hiscale], weights=[w1,w2,w3], bins=np.arange(0,6,0.2),
            #                     legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
            fig,ax = plot.hist([so4_hiscale,so4_m_hiscale], weights=[w1,w2], bins=np.arange(0,6,0.2),
                                legend =['Obs','Model',], color=['k','r'],
                                title = 'Sulfate '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
            fig.savefig(figpath+'hist_SO4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            w1 = np.ones_like(ccn2_hiscale)/sum(~np.isnan(ccn2_hiscale.data))
            w2 = np.ones_like(ccn2_m_hiscale)/sum(~np.isnan(ccn2_m_hiscale.data))
            # w3 = np.ones_like(ccn2_m2_hiscale)/sum(~np.isnan(ccn2_m2_hiscale.data))
            # fig,ax = plot.hist([ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], weights=[w1,w2,w3], bins=np.arange(0,1800,50),
            #                     legend =['Obs','E3SMv1','E3SMv2',], color=['k','r','b'],
            fig,ax = plot.hist([ccn2_hiscale,ccn2_m_hiscale], weights=[w1,w2], bins=np.arange(0,1800,50),
                                legend =['Obs','Model',], color=['k','r'],
                                title = 'CCN (SS=0.2%) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
            fig.savefig(figpath+'hist_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            
            w0 = np.ones_like(cpc10_hiscale)/sum(~np.isnan(cpc10_hiscale.data))
            w1 = np.ones_like(ncn10_m_hiscale)/sum(~np.isnan(ncn10_m_hiscale.data))
            # w2 = np.ones_like(ncn10_m2_hiscale)/sum(~np.isnan(ncn10_m2_hiscale.data))
            # fig,ax = plot.hist([cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale], weights=[w0,w1,w2], bins=np.arange(0,22000,1000),
            #                     legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'],
            fig,ax = plot.hist([cpc10_hiscale,ncn10_m_hiscale], weights=[w0,w1], bins=np.arange(0,22000,1000),
                                legend = ['Obs','Model'], color=['k','r'],
                                title='Aerosol number (>10nm) '+site+' '+IOP,ylabel='Fraction', xlabel='cm$^{-3}$')
            fig.savefig(figpath+'hist_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            w0 = np.ones_like(smps100_hiscale)/sum(~np.isnan(smps100_hiscale.data))
            w1 = np.ones_like(ncn100_m_hiscale)/sum(~np.isnan(ncn100_m_hiscale.data))
            # w2 = np.ones_like(ncn100_m2_hiscale)/sum(~np.isnan(ncn100_m2_hiscale.data))
            # fig,ax = plot.hist([smps100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], weights=[w0,w1,w2], bins=np.arange(0,2100,100),
            #                     legend = ['SMPS100','E3SMv1','E3SMv2'], color=['k','r','b'],
            fig,ax = plot.hist([smps100_hiscale,ncn100_m_hiscale], weights=[w0,w1], bins=np.arange(0,2100,100),
                                legend = ['SMPS100','Model'], color=['k','r'],
                                title='Aerosol number (>100nm) '+site+' '+IOP,  ylabel='Fraction', xlabel='cm$^{-3}$')
            fig.savefig(figpath+'hist_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['tau3d_output'] == True:
            w0 = np.ones_like(cod_hiscale)/sum(~np.isnan(cod_hiscale.data))
            w00 = np.ones_like(cod_sat_hiscale)/sum(~np.isnan(cod_sat_hiscale.data))
            w1 = np.ones_like(cod_m_hiscale)/sum(~np.isnan(cod_m_hiscale.data))
            # w2 = np.ones_like(cod_m2_hiscale)/sum(~np.isnan(cod_m2_hiscale.data))
            # fig,ax = plot.hist( [cod_hiscale, cod_sat_hiscale, cod_m_hiscale, cod_m2_hiscale], weights=[w0,w00,w1,w2], 
            #                     legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
            fig,ax = plot.hist( [cod_hiscale, cod_sat_hiscale, cod_m_hiscale], weights=[w0,w00,w1], 
                                legend = ['MFRSR','Satellite','Model'], color=['k','gray','r'],
                                title='Cloud Optical Depth '+site+' '+IOP, bins=np.arange(0,61,3), ylabel='Fraction', xlabel='N/A')
            fig.savefig(figpath+'hist_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(lwp_hiscale)/sum(~np.isnan(lwp_hiscale.data))
w00 = np.ones_like(lwp_sat_hiscale)/sum(~np.isnan(lwp_sat_hiscale.data))
w1 = np.ones_like(lwp_m_hiscale)/sum(~np.isnan(lwp_m_hiscale.data))
# w2 = np.ones_like(lwp_m2_hiscale)/sum(~np.isnan(lwp_m2_hiscale.data))
# fig,ax = plot.hist([lwp_hiscale, lwp_sat_hiscale, lwp_m_hiscale, lwp_m2_hiscale], weights=[w0,w00,w1,w2], 
#                     legend = ['MWR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.hist([lwp_hiscale, lwp_sat_hiscale, lwp_m_hiscale], weights=[w0,w00,w1], 
                    legend = ['MWR','Satellite','Model'], color=['k','gray','r'],
                    title='LWP '+site+' '+IOP, bins=np.arange(10,410,20), ylabel='Fraction', xlabel="g/m$^2$")
fig.savefig(figpath+'hist_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(ndrop_hiscale)/sum(~np.isnan(ndrop_hiscale.data))
w00 = np.ones_like(nd_sat_hiscale)/sum(~np.isnan(nd_sat_hiscale.data))
w1 = np.ones_like(nd_m_hiscale)/sum(~np.isnan(nd_m_hiscale.data))
# w2 = np.ones_like(nd_m2_hiscale)/sum(~np.isnan(nd_m2_hiscale.data))
# fig,ax = plot.hist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],  weights=[w0,w00,w1,w2], 
#                     legend = ['Ndrop','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.hist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale],  weights=[w0,w00,w1], 
                    legend = ['Ndrop','Satellite','Model'], color=['k','gray','r'],
                    title = 'Nd '+site+' '+IOP, bins=np.arange(0,410,20), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['reff_output'] == True:
            w0 = np.ones_like(reff_hiscale)/sum(~np.isnan(reff_hiscale.data))
            w00 = np.ones_like(reff_sat_hiscale)/sum(~np.isnan(reff_sat_hiscale.data))
            w1 = np.ones_like(reff_m_hiscale)/sum(~np.isnan(reff_m_hiscale.data))
            # w2 = np.ones_like(reff_m2_hiscale)/sum(~np.isnan(reff_m2_hiscale.data))
            # fig,ax = plot.hist([reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale], weights=[w0,w00,w1,w2], 
            #                     legend = ['MFRSR','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
            fig,ax = plot.hist([reff_hiscale,reff_sat_hiscale,reff_m_hiscale], weights=[w0,w00,w1], 
                                legend = ['MFRSR','Satellite','Model'], color=['k','gray','r'],
                                title = 'Cloud Effective Radius '+site+' '+IOP, bins=np.arange(4,28,1), ylabel='Fraction', xlabel='$\mu$m')
            fig.savefig(figpath+'hist_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

pr0 = precip_tbrg_hiscale[precip_tbrg_hiscale!=0]
pr1 = precip_pwd_hiscale[precip_pwd_hiscale!=0]
pr2 = precip_pars_hiscale[precip_pars_hiscale!=0]
prm = precip_m_hiscale[precip_m_hiscale!=0]
# prm2 = precip_m_hiscale[precip_m2_hiscale!=0]
w0 = np.ones_like(pr0)/sum(~np.isnan(pr0.data))
w1 = np.ones_like(pr1)/sum(~np.isnan(pr1.data))
w2 = np.ones_like(pr2)/sum(~np.isnan(pr2.data))
wm1 = np.ones_like(prm)/sum(~np.isnan(prm.data))
# wm2 = np.ones_like(prm2)/sum(~np.isnan(prm2.data))
# fig,ax = plot.hist( [pr0,pr1,pr2,prm,prm2], weights=[w0,w1,w2,wm1,wm2], legend = ['Tipping Bucket','PWD','Disdrometer','E3SMv1','E3SMv2'], 
#                     color=['k','gray','silver','r','b'],  bins=np.arange(0,5,.1), 
fig,ax = plot.hist( [pr0,prm], weights=[w0,w1,w2,wm1], legend = ['Tipping Bucket','PWD','Disdrometer','Model'], 
                    color=['k','gray','silver','r'],  bins=np.arange(0,5,.1), 
                    title = 'Precipitation '+site+' '+IOP, ylabel='Fraction', xlabel='mm/hr')
fig.savefig(figpath+'hist_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(cld_arscl_hiscale)/sum(~np.isnan(cld_arscl_hiscale.data))
w00 = np.ones_like(cld_tsi_hiscale)/sum(~np.isnan(cld_tsi_hiscale.data))
w1 = np.ones_like(cld_m_hiscale)/sum(~np.isnan(cld_m_hiscale.data))
# w2 = np.ones_like(cld_m2_hiscale)/sum(~np.isnan(cld_m2_hiscale.data))
# fig,ax = plot.hist([cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale,cld_m2_hiscale], 
#                     weights=[w0,w00,w1,w2],  bins=np.arange(0,101,5), 
#                     legend = ['ARMBE','TSI','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
fig,ax = plot.hist([cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale], 
                    weights=[w0,w00,w1],  bins=np.arange(0,101,5), 
                    legend = ['ARMBE','TSI','Model'], color=['k','gray','r'],
                      title = 'Cloud Fraction '+site+' '+IOP, ylabel='Fraction', xlabel="%")
fig.savefig(figpath+'hist_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['aerosol_output'] == True:
            #%% mean size distribution
            # fig,ax = plot.mean_size([size_uhsas,size_smps,np.arange(1,3001),np.arange(1,3001)], 
            #             [pdf_uhsas_hiscale, pdf_smps_hiscale, pdf_m_hiscale, pdf_m2_hiscale], 
            #             legend = ['UHSAS','SMPS','E3SMv1','E3SMv2'],color=['k','gray','r','b'], 
            #             marker=['o','+',None,None], linestyles=['none','none','-','-'],
            fig,ax = plot.mean_size([size_uhsas,size_smps,np.arange(1,3001)], 
                        [pdf_uhsas_hiscale, pdf_smps_hiscale, pdf_m_hiscale], 
                        legend = ['UHSAS','SMPS','Model'],color=['k','gray','r'], 
                        marker=['o','+',None], linestyles=['none','none','-'],
                        xlimit=(2, 2e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                        title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
            fig.savefig(figpath+'mean_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statistics
# calc.mean_std_percentiles([org_hiscale,org_m_hiscale,org_m2_hiscale],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_ORG_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([so4_hiscale, so4_m_hiscale, so4_m2_hiscale],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_SO4_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CCN2_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cpc3_hiscale,ncn3_m_hiscale,ncn3_m2_hiscale],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CPC3_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cpc10_hiscale,ncn10_m_hiscale,ncn10_m2_hiscale],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([uhsas100_hiscale, smps100_hiscale, ncn100_m_hiscale, ncn100_m2_hiscale],legend=['UHSAS','SMPS','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CN100_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cod_hiscale,cod_sat_hiscale, cod_m_hiscale, cod_m2_hiscale],legend=['MFRSR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_COD_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale],legend=['MFRSR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Reff_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale],
#                           legend=['MWR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_LWP_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],legend=['Ndrop','Nd_satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Nd_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([precip_tbrg_hiscale,precip_pwd_hiscale,precip_pars_hiscale,precip_m_hiscale,precip_m2_hiscale],legend=['Tipping Bucket','PWD','Disdrometer','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Precip_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cld_arscl_hiscale,cld_tsi_hiscale,cld_m_hiscale,cld_m2_hiscale],
#                           legend=['ARSCL','TSI','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_totcld_'+site+'_'+IOP+'.txt')


# calc.bias_corrcoef_RMSE(org_hiscale,org_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_ORG_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(org_hiscale,org_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_ORG_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(so4_hiscale, so4_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_SO4_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(so4_hiscale, so4_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_SO4_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn2_hiscale,ccn2_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn2_hiscale,ccn2_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(cpc3_hiscale,ncn3_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN3nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(cpc3_hiscale,ncn3_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN3nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(cpc10_hiscale,ncn10_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(cpc10_hiscale,ncn10_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN10nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(smps100_hiscale, ncn100_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN100_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(smps100_hiscale, ncn100_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN100_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(lwp_hiscale, lwp_m_hiscale,label1='ARMBE',label2='E3SMv1', 
#                         outfile=figpath+'statistics_lwp_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(lwp_hiscale, lwp_m2_hiscale,label1='ARMBE',label2='E3SMv2', 
#                         outfile=figpath+'statistics_lwp_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ndrop_hiscale, nd_m_hiscale,label1='Ndrop',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ndrop_hiscale, nd_m2_hiscale,label1='Ndrop',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(nd_sat_hiscale, nd_m_hiscale,label1='Satellite',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsSat_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(nd_sat_hiscale, nd_m2_hiscale,label1='Satellite',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsSat_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(reff_hiscale, reff_m_hiscale,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(reff_hiscale, reff_m2_hiscale,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(reff_sat_hiscale, reff_m_hiscale,label1='Satellite',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsSat_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(reff_sat_hiscale, reff_m2_hiscale,label1='Satellite',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsSat_'+site+'_'+IOP+'.txt')

#%% joint histogram
if config['aerosol_output'] == True:
            # fig,ax = plot.jointhist([uhsas100_hiscale,ncn100_m_hiscale,ncn100_m2_hiscale], [ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], 
            #                     title=['Ground','E3SMv1','E3SMv2']),
            fig,ax = plot.jointhist([uhsas100_hiscale,ncn100_m_hiscale], [ccn2_hiscale,ccn2_m_hiscale], title=['Ground','Model'],
                                xedges=np.arange(0,800,40),yedges=np.arange(0,800,40), normalize_x=True,
                                xlabel='CN (>100nm) (cm$^{-3}$)', ylabel='CCN (SS=0.2%) (cm$^{-3}$)', vmax=0.5)
            fig.savefig(figpath+'jointhist_CN100_CCN2_ship_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.jointhist([ccn2_hiscale,ccn2_hiscale,ccn2_m_hiscale,ccn2_m2_hiscale], 
            #                         [ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],
            #                     title=['Ground','Satellite','E3SMv1','E3SMv2']),
            fig,ax = plot.jointhist([ccn2_hiscale,ccn2_hiscale,ccn2_m_hiscale], [ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale], title=['Ground','Satellite','Model'],
                                xedges=np.arange(0,500,30),yedges=np.arange(0,300,20), normalize_x=True,
                                xlabel='CCN (SS=0.2%) (cm$^{-3}$)', ylabel='Nd (cm$^{-3}$)', vmax=0.4)
            fig.savefig(figpath+'jointhist_CCN2_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.jointhist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],
#                         [lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale], 
#                     title=['Ground','Satellite','E3SMv1','E3SMv2']),
fig,ax = plot.jointhist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale], [lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale], title=['Ground','Satellite','Model'],
                    xedges=np.arange(0,300,20),yedges=np.arange(0,300,20), normalize_x=True,
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', vmax=0.4)
fig.savefig(figpath+'jointhist_LWP_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['reff_output'] == True:
            # fig,ax = plot.jointhist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale,nd_m2_hiscale],
            #                         [reff_hiscale,reff_sat_hiscale,reff_m_hiscale,reff_m2_hiscale],
            #                     title=['Ground','Satellite','E3SMv1','E3SMv2']),
            fig,ax = plot.jointhist([ndrop_hiscale,nd_sat_hiscale,nd_m_hiscale], [reff_hiscale,reff_sat_hiscale,reff_m_hiscale], title=['Ground','Satellite','Model'],
                                xedges=np.arange(0,300,20),yedges=np.arange(4,25,1), normalize_x=True,
                                xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m)', vmax=0.25)
            fig.savefig(figpath+'jointhist_Reff_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if config['tau3d_output'] == True:
            # fig,ax = plot.jointhist([cod_sat_hiscale,cod_sat_hiscale,cod_m_hiscale,cod_m2_hiscale],[lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale], 
            #                     title=['Ground','Satellite','E3SMv1','E3SMv2']),
            fig,ax = plot.jointhist([cod_sat_hiscale,cod_sat_hiscale,cod_m_hiscale],[lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale], title=['Ground','Satellite','Model'],
                                xedges=np.arange(0,40,3),yedges=np.arange(0,300,20), normalize_x=True,
                                xlabel='Cloud Optical Depth (N/A)', ylabel='LWP (g/m$^2$)', vmax=0.25)
            fig.savefig(figpath+'jointhist_COD_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter plot
if config['aerosol_output'] == True:
            # fig,ax = plot.scatter([ndrop_hiscale.data, nd_sat_hiscale.data,nd_m_hiscale.data,nd_m2_hiscale.data], 
            #                       [ccn2_hiscale.data,ccn2_hiscale.data,ccn2_m_hiscale.data,ccn2_m2_hiscale.data],
            #                     title=['Ground','Satellite','E3SMv1','E3SMv2'],
            fig,ax = plot.scatter([ndrop_hiscale.data, nd_sat_hiscale.data,nd_m_hiscale.data], [ccn2_hiscale.data,ccn2_hiscale.data,ccn2_m_hiscale.data],
                                title=['Ground','Satellite','Model'], xlimit=(0,300), ylimit=(0,600),
                                xlabel='Nd (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                            linear_fit=True, intercept=False)
            fig.savefig(figpath+'scatter_Nd_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
            
            # fig,ax = plot.scatter([smps100_hiscale.data,ncn100_m_hiscale.data,ncn100_m2_hiscale.data], 
            #                       [ccn2_hiscale.data,ccn2_m_hiscale.data,ccn2_m2_hiscale.data],
            #                     title=['Ground','E3SMv1','E3SMv2'],
            fig,ax = plot.scatter([smps100_hiscale.data,ncn100_m_hiscale.data], [ccn2_hiscale.data,ccn2_m_hiscale.data],
                                title=['Ground','Model'], xlimit=(0,800), ylimit=(0,800),
                                xlabel='Surface CN (>100nm) (cm$^{-3}$)', ylabel='Surface CCN (SS=0.2%) (cm$^{-3}$)', 
                            linear_fit=True, intercept=True)
            fig.savefig(figpath+'scatter_CN100_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% heatmaps
# xedges=np.exp(np.arange(np.log(10),6.5,0.5))
# yedges=np.exp(np.arange(np.log(10),6.5,0.5))
# fig,ax = plot.heatmap([ndrop_hiscale.data, nd_sat_hiscale.data,nd_m_hiscale.data,nd_m2_hiscale.data],
#                       [lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale,lwp_m2_hiscale],
#                       [albedo_hiscale,albedo_hiscale,albedo_m_hiscale,albedo_m2_hiscale],vmax=60,
fig,ax = plot.heatmap([ndrop_hiscale.data, nd_sat_hiscale.data,nd_m_hiscale.data],
                      [lwp_hiscale,lwp_sat_hiscale,lwp_m_hiscale],
                      [albedo_hiscale,albedo_hiscale,albedo_m_hiscale],vmax=60,
                    xedges=np.arange(0,300,20), yedges=np.arange(10,300,20),
                    # xedges=xedges, yedges=yedges, 
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', zlabel='TOA Albedo (%)',
                    # title=['Ground','Satellite','E3SMv1','E3SMv2'])
                    title=['Ground','Satellite','Model'])
fig.savefig(figpath+'heatmap_Albedo_vs_Nd_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
