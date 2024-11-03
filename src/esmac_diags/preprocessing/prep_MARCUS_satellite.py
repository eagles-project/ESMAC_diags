"""
prepare satellite data from MARCUS
options of output data into coarser resolution
"""

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag
from esmac_diags.subroutines.time_format_change import datetime2cday
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_VISST, calc_clouddepth_VISST, insolation

#%% test settings
# shipmetpath = '../../../data/MARCUS/obs/ship/maraadmetX1.b1/'
# visstgridpath = '../../../data/MARCUS/obs/visst/grid/'
# predatapath = 'C:/Users/tang357/Downloads/prep_data/MARCUS/satellite/'
shipmetpath = '../../../raw_data/obs/MARCUS/ship/maraadmetX1.b1/'
visstgridpath = '../../../raw_data/obs/MARCUS/visst/grid/'
predatapath = '../../../prep_data/MARCUS/satellite/'
dt=3600

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# def prep_VISST_grid(shipmetpath, visstgridpath, predatapath, dt=3600):
#     """
#     prepare VISST-satellite data in grid level (0.5x0.5 degrees)

#     Parameters
#     ----------
#     shipmetpath : str
#         input path for ship location data
#     visstgridpath : str
#         input datapath
#     predatapath : str
#         output datapath
#     dt : float
#         time resolution (unit: sec) of output

#     Returns
#     -------
#     None.

#     """
                                   
#     if not os.path.exists(predatapath):
#         os.makedirs(predatapath)

#%% read in location data
print('read in ship location data:')
lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
if len(lst)==0:
    raise ValueError('cannot find any data')

shipdata = xr.open_mfdataset(lst, combine='by_coords')
shiptime = shipdata['time'].load()
lon0 = shipdata['lon'].load()
qc_lon = shipdata['qc_lon'].load()
lat0 = shipdata['lat'].load()
qc_lat = shipdata['qc_lat'].load()
shipdata.close()

lat0 = qc_mask_qcflag(lat0, qc_lat)
lon0 = qc_mask_qcflag(lon0, qc_lon)

idx = np.logical_and(~np.isnan(lat0), ~np.isnan(lon0))
shiptime = shiptime[idx]
lon0 = lon0[idx]
lat0 = lat0[idx]

#%% read in satellite data
print('read in satellite data:')
lst = glob.glob(os.path.join(visstgridpath, '*visstgrid*.cdf'))
lst.sort()

# first data
visstdata = xr.open_dataset(lst[0])
vissttime = visstdata['time']
lat = visstdata['latitude']
lon = visstdata['longitude']
        
solar_zenith_xyt = visstdata['solar_zenith_angle']
clearsky_vis_reflectance_xyt = visstdata['clearsky_vis_reflectance']
vis_reflectance_all_xyt = visstdata['visible_reflectance'][:,:,:,0]
vis_reflectance_clr_xyt = visstdata['visible_reflectance'][:,:,:,1]
lwp_xyt = visstdata['water_path'][:,:,:,1]
iwp_xyt = visstdata['water_path'][:,:,:,0]
sfc_net_sw_xyt = visstdata['surface_net_shortwave_flux']
sfc_net_lw_xyt = visstdata['surface_net_longwave_flux']
sfc_down_sw_xyt = visstdata['surface_down_shortwave_flux']
sfc_down_lw_xyt = visstdata['surface_down_longwave_flux']
reff_liq_xyt = visstdata['particle_size'][:,:,:,1]
cod_liq_linavg_xyt = visstdata['optical_depth_linear'][:,:,:,2]
cod_liq_logavg_xyt = visstdata['optical_depth_log'][:,:,:,2]
cod_linavg_xyt = visstdata['optical_depth_linear'][:,:,:,0]
cod_logavg_xyt = visstdata['optical_depth_log'][:,:,:,0]
ctt_liq_xyt = visstdata['cloud_temperature'][:,:,:,2]
ctp_liq_xyt = visstdata['cloud_pressure_top'][:,:,:,2]
cth_liq_xyt = visstdata['cloud_height_top'][:,:,:,2]
ctt_all_xyt = visstdata['cloud_temperature'][:,:,:,0]
ctp_all_xyt = visstdata['cloud_pressure_top'][:,:,:,0]
cth_all_xyt = visstdata['cloud_height_top'][:,:,:,0]
cf_all_xyt = visstdata['cloud_percentage'][:,:,:,0]
cf_liq_xyt = visstdata['cloud_percentage'][:,:,:,2]
cf_allz_xyt = visstdata['cloud_percentage_level'][:,:,:,0]
cf_low_xyt = visstdata['cloud_percentage_level'][:,:,:,1]
cf_mid_xyt = visstdata['cloud_percentage_level'][:,:,:,2]
cf_high_xyt = visstdata['cloud_percentage_level'][:,:,:,3]
bb_lw_all_xyt = visstdata['broadband_longwave_flux'][:,:,:,0]
bb_sw_albedo_all_xyt = visstdata['broadband_shortwave_albedo'][:,:,:,0]
bb_lw_clr_xyt = visstdata['broadband_longwave_flux'][:,:,:,1]
bb_sw_albedo_clr_xyt = visstdata['broadband_shortwave_albedo'][:,:,:,1]
visstdata.close()

# get location of ship for VISST time
lon_ship = np.interp(vissttime, shiptime, lon0)
lat_ship = np.interp(vissttime, shiptime, lat0)
lon_ship_all = np.array(lon_ship)
lat_ship_all = np.array(lat_ship)

# choose the data at the ship location
lon.load()
lat.load()
x_idx = abs(lon-lon_ship[0]).argmin()
y_idx = abs(lat-lat_ship[0]).argmin()
solar_zenith = solar_zenith_xyt[0, y_idx, x_idx]
clearsky_vis_reflectance = clearsky_vis_reflectance_xyt[0, y_idx, x_idx]
vis_reflectance_all = vis_reflectance_all_xyt[0, y_idx, x_idx]
vis_reflectance_clr = vis_reflectance_clr_xyt[0, y_idx, x_idx]
lwp = lwp_xyt[0, y_idx, x_idx]
iwp = iwp_xyt[0, y_idx, x_idx]
sfc_net_sw = sfc_net_sw_xyt[0, y_idx, x_idx]
sfc_net_lw = sfc_net_lw_xyt[0, y_idx, x_idx]
sfc_down_sw = sfc_down_sw_xyt[0, y_idx, x_idx]
sfc_down_lw = sfc_down_lw_xyt[0, y_idx, x_idx]
reff_liq = reff_liq_xyt[0, y_idx, x_idx]
cod_liq_linavg = cod_liq_linavg_xyt[0, y_idx, x_idx]
cod_liq_logavg = cod_liq_logavg_xyt[0, y_idx, x_idx]
cod_linavg = cod_linavg_xyt[0, y_idx, x_idx]
cod_logavg = cod_logavg_xyt[0, y_idx, x_idx]
ctt_liq = ctt_liq_xyt[0, y_idx, x_idx]
ctp_liq = ctp_liq_xyt[0, y_idx, x_idx]
cth_liq = cth_liq_xyt[0, y_idx, x_idx]
ctt_all = ctt_all_xyt[0, y_idx, x_idx]
ctp_all = ctp_all_xyt[0, y_idx, x_idx]
cth_all = cth_all_xyt[0, y_idx, x_idx]
cf_all = cf_all_xyt[0, y_idx, x_idx]
cf_liq = cf_liq_xyt[0, y_idx, x_idx]
cf_allz = cf_allz_xyt[0, y_idx, x_idx]
cf_low = cf_low_xyt[0, y_idx, x_idx]
cf_mid = cf_mid_xyt[0, y_idx, x_idx]
cf_high = cf_high_xyt[0, y_idx, x_idx]
bb_lw_all = bb_lw_all_xyt[0, y_idx, x_idx]
bb_lw_clr = bb_lw_clr_xyt[0, y_idx, x_idx]
bb_sw_albedo_all = bb_sw_albedo_all_xyt[0, y_idx, x_idx]
bb_sw_albedo_clr = bb_sw_albedo_clr_xyt[0, y_idx, x_idx]
for tt in range(1,len(vissttime)):
    x_idx = abs(lon-lon_ship[tt]).argmin()
    y_idx = abs(lat-lat_ship[tt]).argmin()
    solar_zenith = xr.concat([solar_zenith, solar_zenith_xyt[tt, y_idx, x_idx]], dim="time")
    clearsky_vis_reflectance = xr.concat([clearsky_vis_reflectance, clearsky_vis_reflectance_xyt[0, y_idx, x_idx]], dim="time")
    vis_reflectance_all = xr.concat([vis_reflectance_all, vis_reflectance_all_xyt[0, y_idx, x_idx]], dim="time")
    vis_reflectance_clr = xr.concat([vis_reflectance_clr, vis_reflectance_clr_xyt[0, y_idx, x_idx]], dim="time")
    lwp = xr.concat([lwp, lwp_xyt[0, y_idx, x_idx]], dim="time")
    iwp = xr.concat([iwp, iwp_xyt[0, y_idx, x_idx]], dim="time")
    sfc_net_sw = xr.concat([sfc_net_sw, sfc_net_sw_xyt[0, y_idx, x_idx]], dim="time")
    sfc_net_lw = xr.concat([sfc_net_lw, sfc_net_lw_xyt[0, y_idx, x_idx]], dim="time")
    sfc_down_sw = xr.concat([sfc_down_sw, sfc_down_sw_xyt[0, y_idx, x_idx]], dim="time")
    sfc_down_lw = xr.concat([sfc_down_lw, sfc_down_lw_xyt[0, y_idx, x_idx]], dim="time")
    reff_liq = xr.concat([reff_liq, reff_liq_xyt[0, y_idx, x_idx]], dim="time")
    cod_liq_linavg = xr.concat([cod_liq_linavg, cod_liq_linavg_xyt[0, y_idx, x_idx]], dim="time")
    cod_liq_logavg = xr.concat([cod_liq_logavg, cod_liq_logavg_xyt[0, y_idx, x_idx]], dim="time")
    cod_linavg = xr.concat([cod_linavg,cod_linavg_xyt[0, y_idx, x_idx]], dim="time")
    cod_logavg = xr.concat([cod_logavg, cod_logavg_xyt[0, y_idx, x_idx]], dim="time")
    ctt_liq = xr.concat([ctt_liq, ctt_liq_xyt[0, y_idx, x_idx]], dim="time")
    ctp_liq = xr.concat([ctp_liq, ctp_liq_xyt[0, y_idx, x_idx]], dim="time")
    cth_liq = xr.concat([cth_liq, cth_liq_xyt[0, y_idx, x_idx]], dim="time")
    ctt_all = xr.concat([ctt_all, ctt_all_xyt[0, y_idx, x_idx]], dim="time")
    ctp_all = xr.concat([ctp_all, ctp_all_xyt[0, y_idx, x_idx]], dim="time")
    cth_all = xr.concat([cth_all, cth_all_xyt[0, y_idx, x_idx]], dim="time")
    cf_all = xr.concat([cf_all, cf_all_xyt[0, y_idx, x_idx]], dim="time")
    cf_liq = xr.concat([cf_liq, cf_liq_xyt[0, y_idx, x_idx]], dim="time")
    cf_allz = xr.concat([cf_allz, cf_allz_xyt[0, y_idx, x_idx]], dim="time")
    cf_low = xr.concat([cf_low, cf_low_xyt[0, y_idx, x_idx]], dim="time")
    cf_mid = xr.concat([cf_mid, cf_mid_xyt[0, y_idx, x_idx]], dim="time")
    cf_high = xr.concat([cf_high, cf_high_xyt[0, y_idx, x_idx]], dim="time")
    bb_lw_all = xr.concat([bb_lw_all, bb_lw_all_xyt[0, y_idx, x_idx]], dim="time")
    bb_lw_clr = xr.concat([bb_lw_clr, bb_lw_clr_xyt[0, y_idx, x_idx]], dim="time")
    bb_sw_albedo_all = xr.concat([bb_sw_albedo_all, bb_sw_albedo_all_xyt[0, y_idx, x_idx]], dim="time")
    bb_sw_albedo_clr = xr.concat([bb_sw_albedo_clr, bb_sw_albedo_clr_xyt[0, y_idx, x_idx]], dim="time")

for filename in lst[1:]:
    print(filename)
    visstdata = xr.open_dataset(filename)
    vissttime2 = visstdata['time']
    lat = visstdata['latitude']
    lon = visstdata['longitude']
            
    solar_zenith_xyt = visstdata['solar_zenith_angle']
    clearsky_vis_reflectance_xyt = visstdata['clearsky_vis_reflectance']
    vis_reflectance_all_xyt = visstdata['visible_reflectance'][:,:,:,0]
    vis_reflectance_clr_xyt = visstdata['visible_reflectance'][:,:,:,1]
    lwp_xyt = visstdata['water_path'][:,:,:,1]
    iwp_xyt = visstdata['water_path'][:,:,:,0]
    sfc_net_sw_xyt = visstdata['surface_net_shortwave_flux']
    sfc_net_lw_xyt = visstdata['surface_net_longwave_flux']
    sfc_down_sw_xyt = visstdata['surface_down_shortwave_flux']
    sfc_down_lw_xyt = visstdata['surface_down_longwave_flux']
    reff_liq_xyt = visstdata['particle_size'][:,:,:,1]
    cod_liq_linavg_xyt = visstdata['optical_depth_linear'][:,:,:,2]
    cod_liq_logavg_xyt = visstdata['optical_depth_log'][:,:,:,2]
    cod_linavg_xyt = visstdata['optical_depth_linear'][:,:,:,0]
    cod_logavg_xyt = visstdata['optical_depth_log'][:,:,:,0]
    ctt_liq_xyt = visstdata['cloud_temperature'][:,:,:,2]
    ctp_liq_xyt = visstdata['cloud_pressure_top'][:,:,:,2]
    cth_liq_xyt = visstdata['cloud_height_top'][:,:,:,2]
    ctt_all_xyt = visstdata['cloud_temperature'][:,:,:,0]
    ctp_all_xyt = visstdata['cloud_pressure_top'][:,:,:,0]
    cth_all_xyt = visstdata['cloud_height_top'][:,:,:,0]
    cf_all_xyt = visstdata['cloud_percentage'][:,:,:,0]
    cf_liq_xyt = visstdata['cloud_percentage'][:,:,:,2]
    cf_allz_xyt = visstdata['cloud_percentage_level'][:,:,:,0]
    cf_low_xyt = visstdata['cloud_percentage_level'][:,:,:,1]
    cf_mid_xyt = visstdata['cloud_percentage_level'][:,:,:,2]
    cf_high_xyt = visstdata['cloud_percentage_level'][:,:,:,3]
    bb_lw_all_xyt = visstdata['broadband_longwave_flux'][:,:,:,0]
    bb_sw_albedo_all_xyt = visstdata['broadband_shortwave_albedo'][:,:,:,0]
    bb_lw_clr_xyt = visstdata['broadband_longwave_flux'][:,:,:,1]
    bb_sw_albedo_clr_xyt = visstdata['broadband_shortwave_albedo'][:,:,:,1]
    visstdata.close()
    
    # get location of ship for VISST time
    lon_ship = np.interp(vissttime2, shiptime, lon0)
    lat_ship = np.interp(vissttime2, shiptime, lat0)
    lon_ship_all = np.hstack((lon_ship_all, lon_ship))
    lat_ship_all = np.hstack((lat_ship_all, lat_ship))
    
    # choose the data at the ship location
    lon.load()
    lat.load()
    for tt in range(len(vissttime2)):
        x_idx = abs(lon-lon_ship[tt]).argmin()
        y_idx = abs(lat-lat_ship[tt]).argmin()
        vissttime = xr.concat([vissttime, vissttime2[tt]], dim="time")
        solar_zenith = xr.concat([solar_zenith, solar_zenith_xyt[tt, y_idx, x_idx]], dim="time")
        clearsky_vis_reflectance = xr.concat([clearsky_vis_reflectance, clearsky_vis_reflectance_xyt[0, y_idx, x_idx]], dim="time")
        vis_reflectance_all = xr.concat([vis_reflectance_all, vis_reflectance_all_xyt[0, y_idx, x_idx]], dim="time")
        vis_reflectance_clr = xr.concat([vis_reflectance_clr, vis_reflectance_clr_xyt[0, y_idx, x_idx]], dim="time")
        lwp = xr.concat([lwp, lwp_xyt[0, y_idx, x_idx]], dim="time")
        iwp = xr.concat([iwp, iwp_xyt[0, y_idx, x_idx]], dim="time")
        sfc_net_sw = xr.concat([sfc_net_sw, sfc_net_sw_xyt[0, y_idx, x_idx]], dim="time")
        sfc_net_lw = xr.concat([sfc_net_lw, sfc_net_lw_xyt[0, y_idx, x_idx]], dim="time")
        sfc_down_sw = xr.concat([sfc_down_sw, sfc_down_sw_xyt[0, y_idx, x_idx]], dim="time")
        sfc_down_lw = xr.concat([sfc_down_lw, sfc_down_lw_xyt[0, y_idx, x_idx]], dim="time")
        reff_liq = xr.concat([reff_liq, reff_liq_xyt[0, y_idx, x_idx]], dim="time")
        cod_liq_linavg = xr.concat([cod_liq_linavg, cod_liq_linavg_xyt[0, y_idx, x_idx]], dim="time")
        cod_liq_logavg = xr.concat([cod_liq_logavg, cod_liq_logavg_xyt[0, y_idx, x_idx]], dim="time")
        cod_linavg = xr.concat([cod_linavg,cod_linavg_xyt[0, y_idx, x_idx]], dim="time")
        cod_logavg = xr.concat([cod_logavg, cod_logavg_xyt[0, y_idx, x_idx]], dim="time")
        ctt_liq = xr.concat([ctt_liq, ctt_liq_xyt[0, y_idx, x_idx]], dim="time")
        ctp_liq = xr.concat([ctp_liq, ctp_liq_xyt[0, y_idx, x_idx]], dim="time")
        cth_liq = xr.concat([cth_liq, cth_liq_xyt[0, y_idx, x_idx]], dim="time")
        ctt_all = xr.concat([ctt_all, ctt_all_xyt[0, y_idx, x_idx]], dim="time")
        ctp_all = xr.concat([ctp_all, ctp_all_xyt[0, y_idx, x_idx]], dim="time")
        cth_all = xr.concat([cth_all, cth_all_xyt[0, y_idx, x_idx]], dim="time")
        cf_all = xr.concat([cf_all, cf_all_xyt[0, y_idx, x_idx]], dim="time")
        cf_liq = xr.concat([cf_liq, cf_liq_xyt[0, y_idx, x_idx]], dim="time")
        cf_allz = xr.concat([cf_allz, cf_allz_xyt[0, y_idx, x_idx]], dim="time")
        cf_low = xr.concat([cf_low, cf_low_xyt[0, y_idx, x_idx]], dim="time")
        cf_mid = xr.concat([cf_mid, cf_mid_xyt[0, y_idx, x_idx]], dim="time")
        cf_high = xr.concat([cf_high, cf_high_xyt[0, y_idx, x_idx]], dim="time")
        bb_lw_all = xr.concat([bb_lw_all, bb_lw_all_xyt[0, y_idx, x_idx]], dim="time")
        bb_lw_clr = xr.concat([bb_lw_clr, bb_lw_clr_xyt[0, y_idx, x_idx]], dim="time")
        bb_sw_albedo_all = xr.concat([bb_sw_albedo_all, bb_sw_albedo_all_xyt[0, y_idx, x_idx]], dim="time")
        bb_sw_albedo_clr = xr.concat([bb_sw_albedo_clr, bb_sw_albedo_clr_xyt[0, y_idx, x_idx]], dim="time")


print('load all data')
solar_zenith.load()
clearsky_vis_reflectance.load()
vis_reflectance_all.load()
vis_reflectance_clr.load()
lwp.load()
iwp.load()
sfc_net_sw.load()
sfc_net_lw.load()
sfc_down_sw.load()
sfc_down_lw.load()
reff_liq.load()
cod_liq_linavg.load()
cod_liq_logavg.load()
cod_linavg.load()
cod_logavg.load()
ctt_liq.load()
ctp_liq.load()
cth_liq.load()
ctt_all.load()
ctp_all.load()
cth_all.load()
cf_all.load()
cf_liq.load()
cf_allz.load()
cf_low.load()
cf_mid.load()
cf_high.load()
bb_lw_all.load()
bb_lw_clr.load()
bb_sw_albedo_all.load()
bb_sw_albedo_clr.load()

#%% calculate TOA SW flux from albedo
print('calculate some variables:')
# change time to calendar day
calday = datetime2cday(vissttime.data)
# calculate insolation
ins = insolation(calday, lon_ship_all, lat_ship_all, leap_year='leap')

# calculate net SW flux
bb_sw_all = ins * (1 - bb_sw_albedo_all*0.01)

#%% retrieve CDNC
# lwp = lwp.data
# ctt = ctt_liq.data
# cod = cod_liq_linavg.data
H = calc_clouddepth_VISST(lwp.data*0.001, ctt_liq.data, adiabaticity=0.8)
H_ad = calc_clouddepth_VISST(lwp.data*0.001, ctt_liq.data, adiabaticity=1.0)
Nd = calc_cdnc_VISST(lwp.data*0.001, ctt_liq.data, cod_liq_linavg.data, adiabaticity=0.8)
Nd_ad = calc_cdnc_VISST(lwp.data*0.001, ctt_liq.data, cod_liq_linavg.data, adiabaticity=1.0)

#filter out columns with ice and bad retrievals
H_array = np.array(H)
H_ad_array = np.array(H_ad)
Nd_array = np.array(Nd)
Nd_ad_array = np.array(Nd_ad)

ind = np.array(iwp > 0)
H_array[ind] = np.nan
H_ad_array[ind] = np.nan
Nd_array[ind] = np.nan
Nd_ad_array[ind] = np.nan

ind = np.isinf(Nd_array)
H_array[ind] = np.nan
H_ad_array[ind] = np.nan
Nd_array[ind] = np.nan
Nd_ad_array[ind] = np.nan

#%% re-shape the data into coarser resolution
time_new = pd.date_range(start='2017-10-21', end='2018-03-23 23:59:00', freq=str(int(dt))+"s")  # MARCUS time period

Nd_new = avg_time_1d(vissttime, Nd_array, time_new, arraytype='numpy')
H_new = avg_time_1d(vissttime, H, time_new, arraytype='numpy')
lwp_new = avg_time_1d(vissttime, lwp, time_new, arraytype='xarray')
iwp_new = avg_time_1d(vissttime, iwp.data, time_new, arraytype='xarray')
swnetsfc_new = avg_time_1d(vissttime, sfc_net_sw.data, time_new, arraytype='xarray')
lwnetsfc_new = avg_time_1d(vissttime, sfc_net_lw.data, time_new, arraytype='xarray')
swdnsfc_new = avg_time_1d(vissttime, sfc_down_sw.data, time_new, arraytype='xarray')
lwdnsfc_new = avg_time_1d(vissttime, sfc_down_lw.data, time_new, arraytype='xarray')
reff_new = avg_time_1d(vissttime, reff_liq.data, time_new, arraytype='xarray')
cod_new = avg_time_1d(vissttime, cod_linavg.data, time_new, arraytype='xarray')
codlog_new = avg_time_1d(vissttime, cod_logavg.data, time_new, arraytype='xarray')
cf_all_new = avg_time_1d(vissttime, cf_allz.data, time_new, arraytype='xarray')
cf_low_new = avg_time_1d(vissttime, cf_low.data, time_new, arraytype='xarray')
cf_mid_new = avg_time_1d(vissttime, cf_mid.data, time_new, arraytype='xarray')
cf_high_new = avg_time_1d(vissttime, cf_high.data, time_new, arraytype='xarray')
ctt_new = avg_time_1d(vissttime, ctt_all.data, time_new, arraytype='xarray')
ctp_new = avg_time_1d(vissttime, ctp_all.data, time_new, arraytype='xarray')
cth_new = avg_time_1d(vissttime, cth_all.data, time_new, arraytype='xarray')
ctt_liq_new = avg_time_1d(vissttime, ctt_liq.data, time_new, arraytype='xarray')
ctp_liq_new = avg_time_1d(vissttime, ctp_liq.data, time_new, arraytype='xarray')
cth_liq_new = avg_time_1d(vissttime, cth_liq.data, time_new, arraytype='xarray')
lw_new = avg_time_1d(vissttime, bb_lw_all.data, time_new, arraytype='xarray')
sw_new = avg_time_1d(vissttime, bb_sw_all.data, time_new, arraytype='xarray')
albedo_new = avg_time_1d(vissttime, bb_sw_albedo_all.data, time_new, arraytype='xarray')

#%% output file
outfile = predatapath + 'Nd_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'Nd': (['time'], np.float32(Nd_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['Nd'].attrs["long_name"] = 'cloud droplet number concentration'
ds['Nd'].attrs["units"] = '#/cm3'

ds.attrs["title"] = 'cloud droplet number concentration retrieved from VISST 0.5x0.5 data'
ds.attrs["description"] = 'retrieved following Bennartz 2007, assuming adiabaticity = 0.8'
ds.attrs["reference"] = 'https://doi.org/10.1029/2006JD007547'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'Hcld_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'Hcld': (['time'], np.float32(H_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['Hcld'].attrs["long_name"] = 'cloud depth for liquid cloud only'
ds['Hcld'].attrs["units"] = 'm'

ds.attrs["title"] = 'liquid cloud depth retrieved from VISST 0.5x0.5 data'
ds.attrs["description"] = 'retrieved from LWP and cloud top temperature assuming adiabaticity = 0.8'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'LWP_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'lwp': (['time'], np.float32(lwp_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['lwp'].attrs["long_name"] = 'liquid water path'
ds['lwp'].attrs["units"] = 'g/m2'

ds.attrs["title"] = 'liquid water path from VISST 0.5x0.5 data'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'IWP_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'iwp': (['time'], np.float32(iwp_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['iwp'].attrs["long_name"] = 'ice water path'
ds['iwp'].attrs["units"] = 'g/m2'

ds.attrs["title"] = 'ice water path from VISST 0.5x0.5 data'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'Reff_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'reff': (['time'], np.float32(reff_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['reff'].attrs["long_name"] = 'effective radius for liquid clouds'
ds['reff'].attrs["units"] = 'um'

ds.attrs["title"] = 'liquid clouds effective radius from VISST 0.5x0.5 data'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'cod_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'cod': (['time'], np.float32(cod_new)),
                'cod_log': (['time'], np.float32(codlog_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['cod'].attrs["long_name"] = 'cloud optical depth for liquid clouds, linear average'
ds['cod'].attrs["units"] = 'N/A'
ds['cod_log'].attrs["long_name"] = 'cloud optical depth for liquid clouds, log average'
ds['cod_log'].attrs["units"] = 'N/A'

ds.attrs["title"] = 'liquid clouds opical depth from VISST 0.5x0.5 data'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'cloudfraction_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'cldtot': (['time'], np.float32(cf_all_new)),
                'cldhigh': (['time'], np.float32(cf_high_new)),
                'cldmid': (['time'], np.float32(cf_mid_new)),
                'cldlow': (['time'], np.float32(cf_low_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['cldtot'].attrs["long_name"] = 'cloud fraction for all heights'
ds['cldtot'].attrs["units"] = '%'
ds['cldhigh'].attrs["long_name"] = 'cloud fraction for high clouds (>6km)'
ds['cldhigh'].attrs["units"] = '%'
ds['cldmid'].attrs["long_name"] = 'cloud fraction for middle clouds (2-6km)'
ds['cldmid'].attrs["units"] = '%'
ds['cldlow'].attrs["long_name"] = 'cloud fraction for low clouds (0-2km)'
ds['cldlow'].attrs["units"] = '%'

ds.attrs["title"] = 'cloud fraction from VISST 0.5x0.5 data'
ds.attrs["date"] = ttt.ctime(ttt.time())
ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'cloudtop_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'ctt': (['time'], np.float32(ctt_new)),
                'cth': (['time'], np.float32(cth_new)),
                'ctp': (['time'], np.float32(ctp_new)),
                'ctt_liq': (['time'], np.float32(ctt_liq_new)),
                'cth_liq': (['time'], np.float32(cth_liq_new)),
                'ctp_liq': (['time'], np.float32(ctp_liq_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['ctt'].attrs["long_name"] = 'cloud top temperature'
ds['ctt'].attrs["units"] = 'K'
ds['ctp'].attrs["long_name"] = 'cloud top pressure'
ds['ctp'].attrs["units"] = 'hPa'
ds['cth'].attrs["long_name"] = 'cloud top height'
ds['cth'].attrs["units"] = 'km'
ds['ctt_liq'].attrs["long_name"] = 'cloud top temperature for liquid clouds'
ds['ctt_liq'].attrs["units"] = 'K'
ds['ctp_liq'].attrs["long_name"] = 'cloud top pressure for liquid clouds'
ds['ctp_liq'].attrs["units"] = 'hPa'
ds['cth_liq'].attrs["long_name"] = 'cloud top height for liquid clouds'
ds['cth_liq'].attrs["units"] = 'km'

ds.attrs["title"] = 'cloud top temperature, pressure and height from VISST 0.5x0.5 data'
ds.attrs["description"] = 'for liquid clouds only'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'lwflx_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'lwnettoa': (['time'], np.float32(lw_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['lwnettoa'].attrs["long_name"] = 'net LW flux at TOA'
ds['lwnettoa'].attrs["units"] = 'W/m2'

ds.attrs["title"] = 'net longwave flux at TOA from VISST 0.5x0.5 data'
ds.attrs["description"] = 'upward positive'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

#
outfile = predatapath + 'swflx_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'swnettoa': (['time'], np.float32(sw_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['swnettoa'].attrs["long_name"] = 'net SW flux at TOA'
ds['swnettoa'].attrs["units"] = 'W/m2'

ds.attrs["title"] = 'net shortwave flux at TOA from VISST 0.5x0.5 data'
ds.attrs["description"] = 'calculate from insolation and SW albedo, downward positive'
ds.attrs["date"] = ttt.ctime(ttt.time())

ds.to_netcdf(outfile, mode='w')

outfile = predatapath + 'albedo_VISSTgrid_MARCUS.nc'
print('output file '+outfile)
ds = xr.Dataset({
                'albedo': (['time'], np.float32(albedo_new)),
                },
                 coords={'time': ('time', time_new)})
#assign attributes
ds['time'].attrs["long_name"] = "Time"
ds['time'].attrs["standard_name"] = "time"
ds['albedo'].attrs["long_name"] = 'broadband shortwave albedo at TOA'
ds['albedo'].attrs["units"] = '%'
ds.attrs["title"] = 'broadband_shortwave_albedo at TOA from VISST 0.5x0.5 data'
ds.attrs["date"] = ttt.ctime(ttt.time())
ds.to_netcdf(outfile, mode='w')
