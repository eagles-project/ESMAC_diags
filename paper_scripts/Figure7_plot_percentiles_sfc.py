
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
from CBcolors import CB_color_cycle

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set output path for plots
figpath= '../figures/'

#%% read variables in different field campaigns
# set site name.
site = 'HISCALE'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
time_hiscale = np.hstack((  \
    pd.date_range(start='2016-04-25', end='2016-05-21', freq="3600s"), \
    pd.date_range(start='2016-08-28', end='2016-09-23', freq="3600s")  \
    ))
            
lst = sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn2 = obsdata['CCN2'].load()
obsdata.close()
ccn2_hiscale = ccn2.sel(time=time_hiscale)
filename = prep_sfc_path + 'sfc_CPC_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cpc = obsdata['time'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()
cpc10_hiscale = cpc10.sel(time=time_hiscale)
filename = prep_sfc_path + 'cod_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cod = obsdata['time'].load()
cod = obsdata['cod'].load()
obsdata.close()
cod_hiscale = cod.sel(time=time_hiscale)
filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_lwp = obsdata['time'].load()
lwp_armbe = obsdata['lwp_armbe'].load()
lwp_mfrsr = obsdata['lwp_mfrsr'].load()
obsdata.close()
lwp_armbe_hiscale = lwp_armbe.sel(time=time_hiscale)
lwp_mfrsr_hiscale = lwp_mfrsr.sel(time=time_hiscale)
filename = prep_sfc_path + 'Ndrop_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_ndrop = obsdata['time'].load()
ndrop = obsdata['cdnc'].load()
obsdata.close()
ndrop_hiscale = ndrop.sel(time=time_hiscale)
filename = prep_sfc_path + 'reff_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_reff = obsdata['time'].load()
reff = obsdata['reff'].load()
obsdata.close()
reff_hiscale = reff.sel(time=time_hiscale)
filename = prep_sfc_path + 'precip_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_pr = obsdata['time'].load()
precip = obsdata['precip'].load()
obsdata.close()
precip_hiscale = precip.sel(time=time_hiscale)
filename = prep_sfc_path + 'sfc_radiation_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_rad = obsdata['time'].load()
lwdnsfc = obsdata['lwdn'].load()
swdnsfc = obsdata['swdn'].load()
lwupsfc = obsdata['lwup'].load()
swupsfc = obsdata['swup'].load()
obsdata.close()
lwdnsfc_hiscale = lwdnsfc.sel(time=time_hiscale)
swdnsfc_hiscale = swdnsfc.sel(time=time_hiscale)
lwupsfc_hiscale = lwupsfc.sel(time=time_hiscale)
swupsfc_hiscale = swupsfc.sel(time=time_hiscale)

filename = prep_sfc_path + 'sfc_SMPS_'+site+'_IOP1.nc'
obsdata = xr.open_dataset(filename)
time1 = obsdata['time'].load()
smps100_1 = obsdata['smps100_dlogDp'].load()
size1 = obsdata['size'].load()
obsdata.close()
filename = prep_sfc_path + 'sfc_SMPS_'+site+'_IOP2.nc'
obsdata = xr.open_dataset(filename)
time2 = obsdata['time'].load()
smps100_2 = obsdata['smps100_dlogDp'].load()
size2 = obsdata['size'].load()
obsdata.close()
time_smps = xr.concat((time1,time2),dim='time')
size_smps = size1
smps100 = xr.concat((smps100_1,smps100_2),dim='time')
# SMPS data is already dN/dlogDp, total number concentration must multiply by dlogDp
dlogDp_smps = np.log10(size_smps[1:].data/size_smps[0:-1].data)
smps100 = smps100 * np.mean(dlogDp_smps)
smps100_hiscale = smps100.sel(time=time_hiscale)

filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_totcld = obsdata['time'].load()
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_visst = obsdata['tot_cld_visst'].load()
obsdata.close()
cld_arscl_hiscale = cld_arscl.sel(time=time_hiscale)
cld_tsi_hiscale = cld_tsi.sel(time=time_hiscale)
cld_visst_hiscale = cld_visst.sel(time=time_hiscale)

obsdata = xr.open_mfdataset(prep_sfc_path + 'cloudheight_ARSCL_'+site+'*.nc')
time_arscl = obsdata['time']
cth = obsdata['cth'].load()
cths = obsdata['cths'].load()
obsdata.close()
cth_hiscale = cth.sel(time=time_hiscale)

# satellite data
filename = prep_sat_path + 'cod_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cod_sat = obsdata['cod'].load()
obsdata.close()
cod_sat_hiscale = cod_sat.sel(time=time_hiscale)
filename = prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
reff_sat = obsdata['reff'].load()
obsdata.close()
reff_sat_hiscale = reff_sat.sel(time=time_hiscale)
filename = prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwp_sat = obsdata['lwp'].load()
obsdata.close()
lwp_sat_hiscale = lwp_sat.sel(time=time_hiscale)
filename = prep_sat_path + 'IWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
iwp_sat = obsdata['iwp'].load()
obsdata.close()
iwp_sat_hiscale = iwp_sat.sel(time=time_hiscale)
filename = prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
nd_sat = obsdata['Nd'].load()
obsdata.close()
nd_sat_hiscale = nd_sat.sel(time=time_hiscale)
filename = prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwnettoa = obsdata['lwnettoa'].load()
obsdata.close()
lwnettoa_hiscale = lwnettoa.sel(time=time_hiscale)
filename = prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
swnettoa = obsdata['swnettoa'].load()
obsdata.close()
swnettoa_hiscale = swnettoa.sel(time=time_hiscale)
filename = prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cth_sat = obsdata['cth'].load()
ctt_sat = obsdata['ctt'].load()
obsdata.close()
cth_sat_hiscale = cth_sat.sel(time=time_hiscale)
ctt_sat_hiscale = ctt_sat.sel(time=time_hiscale)

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
iwp_m = modeldata['TGCLDIWP'].load()
nd_m = modeldata['Nd_mean'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
cth_m = modeldata['cth'].load()
ctt_m = modeldata['ctt'].load()
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
ccn2_m_hiscale = ccn2_m.sel(time=time_hiscale)
ncn10_m_hiscale = ncn10_m.sel(time=time_hiscale)
ncn100_m_hiscale = ncn100_m.sel(time=time_hiscale)
cod_m_hiscale = cod_m.sel(time=time_hiscale)
reff_m_hiscale = reff_m.sel(time=time_hiscale)
lwp_m_hiscale = lwp_m.sel(time=time_hiscale)
iwp_m_hiscale = iwp_m.sel(time=time_hiscale)
nd_m_hiscale = nd_m.sel(time=time_hiscale)
precip_m_hiscale = precip_m.sel(time=time_hiscale)
cld_m_hiscale = cld_m.sel(time=time_hiscale)
cth_m_hiscale = cth_m.sel(time=time_hiscale)
ctt_m_hiscale = ctt_m.sel(time=time_hiscale)
lwdnsfc_m_hiscale = lwdnsfc_m.sel(time=time_hiscale)
lwupsfc_m_hiscale = lwupsfc_m.sel(time=time_hiscale)
lwnettoa_m_hiscale = lwnettoa_m.sel(time=time_hiscale)
swdnsfc_m_hiscale = swdnsfc_m.sel(time=time_hiscale)
swupsfc_m_hiscale = swupsfc_m.sel(time=time_hiscale)
swnettoa_m_hiscale = swnettoa_m.sel(time=time_hiscale)


# set site name.
site = 'ACEENA'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
time_aceena = np.hstack((  \
    pd.date_range(start='2017-06-21', end='2017-07-21', freq="3600s"), \
    pd.date_range(start='2018-01-15', end='2018-02-19', freq="3600s")  \
    ))
            
filename = prep_sfc_path + 'sfc_CCN_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_ccn = obsdata['time'].load()
ccn2 = obsdata['CCN2'].load()
obsdata.close()
ccn2_aceena = ccn2.sel(time=time_aceena)
filename = prep_sfc_path + 'sfc_CPC_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cpc = obsdata['time'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()
cpc10_aceena = cpc10.sel(time=time_aceena)
filename = prep_sfc_path + 'sfc_CPC_'+site+'_withmask.nc'
obsdata = xr.open_dataset(filename)
time_mask = obsdata['time'].load()
cpc_mask = obsdata['cpc_masked'].load()
obsdata.close()
cpc10_withmask_aceena = cpc_mask.sel(time=time_aceena)
filename = prep_sfc_path + 'cod_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cod = obsdata['time'].load()
cod = obsdata['cod'].load()
obsdata.close()
cod_aceena = cod.sel(time=time_aceena)
filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_lwp = obsdata['time'].load()
lwp_armbe = obsdata['lwp_armbe'].load()
lwp_mfrsr = obsdata['lwp_mfrsr'].load()
obsdata.close()
lwp_armbe_aceena = lwp_armbe.sel(time=time_aceena)
lwp_mfrsr_aceena = lwp_mfrsr.sel(time=time_aceena)
filename = prep_sfc_path + 'Ndrop_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_ndrop = obsdata['time'].load()
ndrop = obsdata['cdnc'].load()
obsdata.close()
ndrop_aceena = ndrop.sel(time=time_aceena)
filename = prep_sfc_path + 'reff_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_reff = obsdata['time'].load()
reff = obsdata['reff'].load()
obsdata.close()
reff_aceena = reff.sel(time=time_aceena)
filename = prep_sfc_path + 'Nd_Reff_Wu_etal_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_wu = obsdata['time'].load()
reff_wu = obsdata['reff'].load()
nd_wu = obsdata['cdnc'].load()
obsdata.close()
reff_wu_aceena = reff_wu.sel(time=time_aceena)
nd_wu_aceena = nd_wu.sel(time=time_aceena)
filename = prep_sfc_path + 'precip_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_pr = obsdata['time'].load()
precip = obsdata['precip'].load()
obsdata.close()
precip_aceena = precip.sel(time=time_aceena)
filename = prep_sfc_path + 'sfc_radiation_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_rad = obsdata['time'].load()
lwdnsfc = obsdata['lwdn'].load()
swdnsfc = obsdata['swdn'].load()
lwupsfc = obsdata['lwup'].load()
swupsfc = obsdata['swup'].load()
obsdata.close()
lwdnsfc_aceena = lwdnsfc.sel(time=time_aceena)
swdnsfc_aceena = swdnsfc.sel(time=time_aceena)
lwupsfc_aceena = lwupsfc.sel(time=time_aceena)
swupsfc_aceena = swupsfc.sel(time=time_aceena)

filename = prep_sfc_path + 'sfc_UHSAS_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_pr = obsdata['time'].load()
uhsas100 = obsdata['uhsas100'].load()
obsdata.close()
uhsas100_aceena = uhsas100.sel(time=time_aceena)
filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_totcld = obsdata['time'].load()
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_visst = obsdata['tot_cld_visst'].load()
obsdata.close()
cld_arscl_aceena = cld_arscl.sel(time=time_aceena)
cld_tsi_aceena = cld_tsi.sel(time=time_aceena)
cld_visst_aceena = cld_visst.sel(time=time_aceena)
obsdata = xr.open_mfdataset(prep_sfc_path + 'cloudheight_ARSCL_'+site+'*.nc')
time_arscl = obsdata['time']
cth = obsdata['cth'].load()
cths = obsdata['cths'].load()
obsdata.close()
cth_aceena = cth.sel(time=time_aceena)

# satellite data
filename = prep_sat_path + 'cod_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cod_sat = obsdata['cod'].load()
obsdata.close()
cod_sat_aceena = cod_sat.sel(time=time_aceena)
filename = prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
reff_sat = obsdata['reff'].load()
obsdata.close()
reff_sat_aceena = reff_sat.sel(time=time_aceena)
filename = prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwp_sat = obsdata['lwp'].load()
obsdata.close()
lwp_sat_aceena = lwp_sat.sel(time=time_aceena)
filename = prep_sat_path + 'IWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
iwp_sat = obsdata['iwp'].load()
obsdata.close()
iwp_sat_aceena = iwp_sat.sel(time=time_aceena)
filename = prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
nd_sat = obsdata['Nd'].load()
obsdata.close()
nd_sat_aceena = nd_sat.sel(time=time_aceena)
filename = prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwnettoa = obsdata['lwnettoa'].load()
obsdata.close()
lwnettoa_aceena = lwnettoa.sel(time=time_aceena)
filename = prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
swnettoa = obsdata['swnettoa'].load()
obsdata.close()
swnettoa_aceena = swnettoa.sel(time=time_aceena)
filename = prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cth_sat = obsdata['cth'].load()
ctt_sat = obsdata['ctt'].load()
obsdata.close()
cth_sat_aceena = cth_sat.sel(time=time_aceena)
ctt_sat_aceena = ctt_sat.sel(time=time_aceena)

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
iwp_m = modeldata['TGCLDIWP'].load()
nd_m = modeldata['Nd_mean'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
cth_m = modeldata['cth'].load()
ctt_m = modeldata['ctt'].load()
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
ccn2_m_aceena = ccn2_m.sel(time=time_aceena)
ncn10_m_aceena = ncn10_m.sel(time=time_aceena)
ncn100_m_aceena = ncn100_m.sel(time=time_aceena)
cod_m_aceena = cod_m.sel(time=time_aceena)
reff_m_aceena = reff_m.sel(time=time_aceena)
lwp_m_aceena = lwp_m.sel(time=time_aceena)
iwp_m_aceena = iwp_m.sel(time=time_aceena)
nd_m_aceena = nd_m.sel(time=time_aceena)
precip_m_aceena = precip_m.sel(time=time_aceena)
cld_m_aceena = cld_m.sel(time=time_aceena)
cth_m_aceena = cth_m.sel(time=time_aceena)
ctt_m_aceena = ctt_m.sel(time=time_aceena)
lwdnsfc_m_aceena = lwdnsfc_m.sel(time=time_aceena)
lwupsfc_m_aceena = lwupsfc_m.sel(time=time_aceena)
lwnettoa_m_aceena = lwnettoa_m.sel(time=time_aceena)
swdnsfc_m_aceena = swdnsfc_m.sel(time=time_aceena)
swupsfc_m_aceena = swupsfc_m.sel(time=time_aceena)
swnettoa_m_aceena = swnettoa_m.sel(time=time_aceena)


# set site name.
site = 'MAGIC'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/ship/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
            
filename = prep_sfc_path + 'CCN_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
ccn2_magic = obsdata['CCN2'].load()
lon_magic = obsdata['lon']
lat_magic = obsdata['lat']
obsdata.close()
filename = prep_sfc_path + 'CN_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
cpc10_magic = obsdata['CPC10'].load()
uhsas100_magic = obsdata['UHSAS100'].load()
obsdata.close()
filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
lwp_magic = obsdata['lwp'].load()
obsdata.close()
filename = prep_sfc_path + 'Nd_Reff_Wu_etal_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
reff_wu_magic = obsdata['reff'].load()
nd_wu_magic = obsdata['cdnc'].load()
obsdata.close()

filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_magic = obsdata['time'].load()
cld_magic = obsdata['cldfrac'].load()
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
filename = prep_sat_path + 'IWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
iwp_sat_magic = obsdata['iwp'].load()
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
filename = prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cth_sat_magic = obsdata['cth'].load()
ctt_sat_magic = obsdata['ctt'].load()
obsdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv2_'+site+'_ship.nc'
modeldata = xr.open_dataset(filename)
time_m_magic = modeldata['time'].load()
ccn2_m_magic = modeldata['CCN4'].load()
ncn10_m_magic = modeldata['NCN10'].load()
ncn100_m_magic = modeldata['NCN100'].load()
cod_m_magic = modeldata['cod'].load()
reff_m_magic = modeldata['reff'].load()
lwp_m_magic = modeldata['TGCLDLWP'].load()
iwp_m_magic = modeldata['TGCLDIWP'].load()
nd_m_magic = modeldata['Nd_mean'].load()
precip_m_magic = modeldata['PRECT'].load()
cld_m_magic = modeldata['CLDTOT'].load()
cth_m_magic = modeldata['cth'].load()
ctt_m_magic = modeldata['ctt'].load()
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


# set site name.
site = 'MARCUS'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/ship/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
            
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
filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
lwp_marcus = obsdata['lwp'].load()
obsdata.close()

filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_marcus = obsdata['time'].load()
cld_marcus = obsdata['cldfrac'].load()
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
filename = prep_sat_path + 'IWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
iwp_sat_marcus = obsdata['iwp'].load()
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
filename = prep_sat_path + 'cloudtop_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cth_sat_marcus = obsdata['cth'].load()
ctt_sat_marcus = obsdata['ctt'].load()
obsdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv2_'+site+'_ship.nc'
modeldata = xr.open_dataset(filename)
time_m_marcus = modeldata['time'].load()
ccn2_m_marcus = modeldata['CCN4'].load()
ncn10_m_marcus = modeldata['NCN10'].load()
ncn100_m_marcus = modeldata['NCN100'].load()
cod_m_marcus = modeldata['cod'].load()
reff_m_marcus = modeldata['reff'].load()
lwp_m_marcus = modeldata['TGCLDLWP'].load()
iwp_m_marcus = modeldata['TGCLDIWP'].load()
nd_m_marcus = modeldata['Nd_mean'].load()
precip_m_marcus = modeldata['PRECT'].load()
cld_m_marcus = modeldata['CLDTOT'].load()
cth_m_marcus = modeldata['cth'].load()
ctt_m_marcus = modeldata['ctt'].load()
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

ndrop_hiscale[ndrop_hiscale<10] = np.nan
nd_sat_hiscale[nd_sat_hiscale<10] = np.nan
nd_m_hiscale[nd_m_hiscale<10] = np.nan
ndrop_aceena[ndrop_aceena<10] = np.nan
nd_wu_aceena[nd_wu_aceena<10] = np.nan
nd_sat_aceena[nd_sat_aceena<10] = np.nan
nd_m_aceena[nd_m_aceena<10] = np.nan
nd_wu_magic[nd_wu_magic<10] = np.nan
nd_sat_magic[nd_sat_magic<10] = np.nan
nd_m_magic[nd_m_magic<10] = np.nan
nd_sat_marcus[nd_sat_marcus<10] = np.nan
nd_m_marcus[nd_m_marcus<10] = np.nan
ndrop_hiscale[ndrop_hiscale>800] = np.nan
nd_sat_hiscale[nd_sat_hiscale>800] = np.nan
nd_m_hiscale[nd_m_hiscale>800] = np.nan
ndrop_aceena[ndrop_aceena>800] = np.nan
nd_wu_aceena[nd_wu_aceena>800] = np.nan
nd_sat_aceena[nd_sat_aceena>800] = np.nan
nd_m_aceena[nd_m_aceena>800] = np.nan
nd_wu_magic[nd_wu_magic>800] = np.nan
nd_sat_magic[nd_sat_magic>800] = np.nan
nd_m_magic[nd_m_magic>800] = np.nan
nd_sat_marcus[nd_sat_marcus>800] = np.nan
nd_m_marcus[nd_m_marcus>800] = np.nan

lwp_armbe_hiscale[lwp_armbe_hiscale<20] = np.nan
lwp_mfrsr_hiscale[lwp_mfrsr_hiscale<20] = np.nan
lwp_sat_hiscale[lwp_sat_hiscale<20] = np.nan
lwp_m_hiscale[lwp_m_hiscale<20] = np.nan
lwp_armbe_aceena[lwp_armbe_aceena<20] = np.nan
lwp_mfrsr_aceena[lwp_mfrsr_aceena<20] = np.nan
lwp_sat_aceena[lwp_sat_aceena<20] = np.nan
lwp_m_aceena[lwp_m_aceena<20] = np.nan
lwp_magic[lwp_magic<20] = np.nan
lwp_sat_magic[lwp_sat_magic<20] = np.nan
lwp_m_magic[lwp_m_magic<20] = np.nan
lwp_marcus[lwp_marcus<20] = np.nan
lwp_sat_marcus[lwp_sat_marcus<20] = np.nan
lwp_m_marcus[lwp_m_marcus<20] = np.nan

cod_hiscale[cod_hiscale<2] = np.nan
cod_sat_hiscale[cod_sat_hiscale<2] = np.nan
cod_m_hiscale[cod_m_hiscale<2] = np.nan
cod_aceena[cod_aceena<2] = np.nan
cod_sat_aceena[cod_sat_aceena<2] = np.nan
cod_m_aceena[cod_m_aceena<2] = np.nan
cod_sat_magic[cod_sat_magic<2] = np.nan
cod_m_magic[cod_m_magic<2] = np.nan
cod_sat_marcus[cod_sat_marcus<2] = np.nan
cod_m_marcus[cod_m_marcus<2] = np.nan
cod_hiscale[cod_hiscale>100] = np.nan
cod_sat_hiscale[cod_sat_hiscale>100] = np.nan
cod_m_hiscale[cod_m_hiscale>100] = np.nan
cod_aceena[cod_aceena>100] = np.nan
cod_sat_aceena[cod_sat_aceena>100] = np.nan
cod_m_aceena[cod_m_aceena>100] = np.nan
cod_sat_magic[cod_sat_magic>100] = np.nan
cod_m_magic[cod_m_magic>100] = np.nan
cod_sat_marcus[cod_sat_marcus>100] = np.nan
cod_m_marcus[cod_m_marcus>100] = np.nan
# remove MAGIC data near continent
cpc10_magic[lon_magic>-122]=np.nan
ncn10_m_magic[lon_magic>-122]=np.nan
uhsas100_magic[lon_magic>-122]=np.nan
ncn100_m_magic[lon_magic>-122]=np.nan
nd_wu_magic[lon_magic>-122]=np.nan
nd_sat_magic[lon_magic>-122]=np.nan
nd_m_magic[lon_magic>-122]=np.nan
reff_wu_magic[lon_magic>-122]=np.nan
reff_sat_magic[lon_magic>-122]=np.nan
reff_m_magic[lon_magic>-122]=np.nan
lwp_m_magic[lon_magic>-122]=np.nan
lwp_magic[lon_magic>-122]=np.nan
lwp_sat_magic[lon_magic>-122]=np.nan
cod_sat_magic[lon_magic>-122]=np.nan


# index of overcasting low clouds
cth_sat_hiscale = cth_sat_hiscale*1000  # km to m
cth_sat_aceena = cth_sat_aceena*1000  # km to m
cth_sat_magic = cth_sat_magic*1000  # km to m
cth_sat_marcus = cth_sat_marcus*1000  # km to m

i_hiscale_sfc = np.full(cld_arscl_hiscale.shape, True)
i_hiscale_sat = np.full(cld_visst_hiscale.shape, True)
i_hiscale_m = np.full(cld_m_hiscale.shape, True)
i_aceena_sfc = np.full(cld_arscl_aceena.shape, True)
i_aceena_sat = np.full( cld_visst_aceena.shape, True)
i_aceena_m = np.full(cld_m_aceena.shape, True)
i_magic_sfc = np.full( cld_magic.shape, True)
i_magic_sat = np.full( cld_sat_magic.shape, True)
i_magic_m = np.full( cld_m_magic.shape, True)
i_marcus_sfc = np.full(cld_marcus.shape, True)
i_marcus_sat = np.full( cld_sat_marcus.shape, True)
i_marcus_m = np.full(cld_m_marcus.shape, True)

i_hiscale_sfc[cld_arscl_hiscale<90] = False
i_hiscale_sat[cld_visst_hiscale<90] = False
i_hiscale_m[cld_m_hiscale<90] = False
i_aceena_sfc[cld_arscl_aceena<90] = False 
i_aceena_sat[cld_visst_aceena<90] = False
i_aceena_m[cld_m_aceena<90] = False
i_magic_sfc[cld_magic<90] = False
i_magic_sat[cld_sat_magic<90] = False 
i_magic_m[cld_m_magic<90] = False
i_marcus_sfc[cld_marcus<90] = False
i_marcus_sat[cld_sat_marcus<90] = False
i_marcus_m[cld_m_marcus<90] = False

i_hiscale_sfc[cth_hiscale>4000] = False
i_hiscale_sat[cth_sat_hiscale>4000] = False
i_hiscale_m[cth_m_hiscale>4000] = False
i_aceena_sfc[cth_aceena>4000] = False 
i_aceena_sat[cth_sat_aceena>4000] = False
i_aceena_m[cth_m_aceena>4000] = False
i_magic_sfc[cth_sat_magic>4000] = False
i_magic_sat[cth_sat_magic>4000] = False 
i_magic_m[cth_m_magic>4000] = False
i_marcus_sfc[cth_sat_marcus>4000] = False
i_marcus_sat[cth_sat_marcus>4000] = False
i_marcus_m[cth_m_marcus>4000] = False


i_hiscale_sfc[iwp_sat_hiscale>10] = False
i_hiscale_sat[iwp_sat_hiscale>10] = False
i_hiscale_m[iwp_m_hiscale>10] = False
i_aceena_sfc[iwp_sat_aceena>10] = False 
i_aceena_sat[iwp_sat_aceena>10] = False
i_aceena_m[iwp_m_aceena>10] = False
i_magic_sfc[iwp_sat_magic>10] = False
i_magic_sat[iwp_sat_magic>10] = False 
i_magic_m[iwp_m_magic>10] = False
i_marcus_sfc[iwp_sat_marcus>10] = False
i_marcus_sat[iwp_sat_marcus>10] = False
i_marcus_m[iwp_m_marcus>10] = False

i_hiscale_sfc[ctt_sat_hiscale<273] = False
i_hiscale_sat[ctt_sat_hiscale<273] = False
i_hiscale_m[ctt_m_hiscale<273] = False
i_aceena_sfc[ctt_sat_aceena<273] = False 
i_aceena_sat[ctt_sat_aceena<273] = False
i_aceena_m[ctt_m_aceena<273] = False
i_magic_sfc[ctt_sat_magic<273] = False
i_magic_sat[ctt_sat_magic<273] = False 
i_magic_m[ctt_m_magic<273] = False
i_marcus_sfc[ctt_sat_marcus<273] = False
i_marcus_sat[ctt_sat_marcus<273] = False
i_marcus_m[ctt_m_marcus<273] = False


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot, all data
if not os.path.exists(figpath):
    os.makedirs(figpath)
    

#%% overcasting low cloud only
fig,ax = plot.percentiles([ccn2_hiscale[i_hiscale_sfc], ccn2_m_hiscale[i_hiscale_m], ], 
                          [ccn2_aceena[i_aceena_sfc], ccn2_m_aceena[i_aceena_m], ], 
                          [ccn2_magic[i_magic_sfc], ccn2_m_magic[i_magic_m], ], 
                          [ccn2_marcus[i_marcus_sfc], ccn2_m_marcus[i_marcus_m], ],
                    title='0.2% CCN (cm$^{-3}$)', figsize=(10,2.5),ylimit=(0,1200), 
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1], 'k'], legend=None)
# fig.savefig(figpath+'percentiles_CCN_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([cpc10_hiscale[i_hiscale_sfc], ncn10_m_hiscale[i_hiscale_m], ], 
                          [cpc10_withmask_aceena[i_aceena_sfc], ncn10_m_aceena[i_aceena_m], ],
                          [cpc10_magic[i_magic_sfc], ncn10_m_magic[i_magic_m], ],
                          [cpc10_marcus[i_marcus_sfc], ncn10_m_marcus[i_marcus_m],],
                    title='(a) CN (>10nm) (cm$^{-3}$)', figsize=(10,2.5), #ylimit=(0,8000),
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1], 'k'], legend=None)
# fig.savefig(figpath+'percentiles_CN10_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([smps100_hiscale[i_hiscale_sfc], ncn100_m_hiscale[i_hiscale_m],], 
                          [uhsas100_aceena[i_aceena_sfc], ncn100_m_aceena[i_aceena_m],],
                          [uhsas100_magic[i_magic_sfc], ncn100_m_magic[i_magic_m],],
                          [uhsas100_marcus[i_marcus_sfc], ncn100_m_marcus[i_marcus_m], ],
                    title='(b) CN (>100nm) (cm$^{-3}$)',figsize=(10,2.5),ylimit=(0,1200),
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1], 'k'], legend=None)
# fig.savefig(figpath+'percentiles_CN100_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([ndrop_hiscale[i_hiscale_sfc], nd_sat_hiscale[i_hiscale_sat], nd_m_hiscale[i_hiscale_m], ], 
                          [ndrop_aceena[i_aceena_sfc], nd_sat_aceena[i_aceena_sat], nd_m_aceena[i_aceena_m],],
                          [nd_wu_magic[i_magic_sfc], nd_sat_magic[i_magic_sat], nd_m_magic[i_magic_m],],
                          [nd_m_marcus*np.nan, nd_sat_marcus[i_marcus_sat], nd_m_marcus[i_marcus_m],],
                    title='(c) Nd (cm$^{-3}$)',figsize=(10,2.5),ylimit=(0,500),
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1],CB_color_cycle[0], 'k'], legend=['Ground/Ship','Satellite','E3SMv2',])
# fig.savefig(figpath+'percentiles_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([lwp_armbe_hiscale[i_hiscale_sfc], lwp_sat_hiscale[i_hiscale_sat], lwp_m_hiscale[i_hiscale_m],], 
                          [lwp_armbe_aceena[i_aceena_sfc], lwp_sat_aceena[i_aceena_sat], lwp_m_aceena[i_aceena_m],],
                          [lwp_magic[i_magic_sfc], lwp_sat_magic[i_magic_sat], lwp_m_magic[i_magic_m],],
                          [lwp_marcus[i_marcus_sfc], lwp_sat_marcus[i_marcus_sat], lwp_m_marcus[i_marcus_m],],
                    title='(d) LWP (g/m$^2$)',figsize=(10,2.5),
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1],CB_color_cycle[0], 'k'], legend=None)
# fig.savefig(figpath+'percentiles_lwp_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([reff_hiscale[i_hiscale_sfc], reff_sat_hiscale[i_hiscale_sat], reff_m_hiscale[i_hiscale_m],], 
                          [reff_aceena[i_aceena_sfc], reff_sat_aceena[i_aceena_sat], reff_m_aceena[i_aceena_m],],
                          [reff_wu_magic[i_magic_sfc], reff_sat_magic[i_magic_sat], reff_m_magic[i_magic_m], ],
                          [reff_m_marcus*np.nan, reff_sat_marcus[i_marcus_sat], reff_m_marcus[i_marcus_m], ],
                    title='(e) Reff ($\mu$m)',figsize=(10,2.5),#ylimit=(0,500),
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1],CB_color_cycle[0], 'k'], legend=None)
# fig.savefig(figpath+'percentiles_Reff_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([cod_hiscale[i_hiscale_sfc], cod_sat_hiscale[i_hiscale_sat], cod_m_hiscale[i_hiscale_m],], 
                          [cod_aceena[i_aceena_sfc], cod_sat_aceena[i_aceena_sat], cod_m_aceena[i_aceena_m],],
                          [cod_m_magic*np.nan, cod_sat_magic[i_magic_sat], cod_m_magic[i_magic_m], ],
                          [cod_m_marcus*np.nan, cod_sat_marcus[i_marcus_sat], cod_m_marcus[i_marcus_m], ],
                    title='(f) Cloud Optical Depth (unitless)',figsize=(10,2.5),
                    xlabel=['HI-SCALE','ACE-ENA','MAGIC','MARCUS'], 
                    color=[CB_color_cycle[1],CB_color_cycle[0], 'k'], legend=None)
# fig.savefig(figpath+'percentiles_cod_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
