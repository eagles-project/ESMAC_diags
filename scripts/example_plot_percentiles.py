"""
example to generate percentiles for multiple variables

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
# set output path for plots
figpath= '../figures/'+site+'/surface/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
filename = prep_obs_path + 'Ndrop_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_nd = obsdata['time'].load()
ndrop = obsdata['cdnc'].load()
obsdata.close()

filename = prep_obs_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_lwp = obsdata['time'].load()
lwp_sfc = obsdata['lwp_armbe'].load()
lwp2_sfc = obsdata['lwp_mfrsr'].load()
obsdata.close()

filename = prep_obs_path + 'reff_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_reff = obsdata['time'].load()
reff_sfc = obsdata['reff'].load()
obsdata.close()

filename = prep_obs_path + 'cod_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cod = obsdata['time'].load()
cod_sfc = obsdata['cod'].load()
obsdata.close()

obsdata = xr.open_dataset(prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc')
time_sat = obsdata['time'].load()
nd_sat = obsdata['Nd'].load()
obsdata.close()

obsdata = xr.open_dataset(prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc')
time_sat = obsdata['time'].load()
reff_sat = obsdata['reff'].load()
obsdata.close()

obsdata = xr.open_dataset(prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc')
time_sat = obsdata['time'].load()
lwp_sat = obsdata['lwp'].load()
obsdata.close()

obsdata = xr.open_dataset(prep_sat_path + 'cod_VISSTgrid_'+site+'.nc')
time_sat = obsdata['time'].load()
cod_sat = obsdata['cod'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv1_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
nd_m2 = modeldata['Nd_mean'].load()
lwp_m2 = modeldata['TGCLDLWP'].load()
iwp_m2 = modeldata['TGCLDIWP'].load()
cod_m2 = modeldata['cod'].load()
reff_m2 = modeldata['reff'].load()
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

ndrop[ndrop<10] = np.nan
nd_sat[nd_sat<10] = np.nan
nd_m2[nd_m2<10] = np.nan

ndrop[ndrop>2000] = np.nan
nd_sat[nd_sat>2000] = np.nan
nd_m2[nd_m2>2000] = np.nan

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
fig,ax = plot.percentiles([reff_sfc, reff_sat, reff_m2], 
                          [lwp_sfc*0.1, lwp_sat*0.1,lwp_m2*0.1], 
                       [cod_sfc, cod_sat,cod_m2],
                       [ndrop*0.1, nd_sat*0.1,nd_m2*0.1],
                    xlabel=['Reff ($\mu$m)', 'LWP (x10 g/m$^2$)', 'COD (NA)', 'Nd (x10 cm$^{-3}$)'], 
                    color=['k','b','r'], legend=['SFC','VISST','E3SMv1'])

fig.savefig(figpath+'percentiles_Reff_LWP_COD_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   
