"""
example to generate seasonal cycle of 2d variable (cloud)

"""
import os
import glob
import numpy as np
import xarray as xr
import esmac_diags.plotting.plot_esmac_diags as plot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'SGP'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/surface/'
# set output path for plots
figpath= '../figures/'+site+'/surface/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = glob.glob(prep_obs_path + 'cloud_2d_'+site+'*.nc')
obsdata = xr.open_mfdataset(lst)
time_cf = obsdata['time'].load()
height1 = obsdata['height'].load()
cf_obs = obsdata['cloud'].load()
obsdata.close()


filename = prep_model_path + 'E3SMv1_'+site+'_profiles.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
heightm = modeldata['height'].load()
cf_e3sm = modeldata['cloud_z'].load()
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment
cf_e3sm = cf_e3sm*100  # fraction to %
height1 = height1.data*0.001   # m to km
heightm = heightm.data*0.001   # m to km


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot 
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.seasonalcycle_2d([cf_obs, cf_e3sm], y = [height1, heightm],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)',  cmap='jet',
                        levellist=np.arange(0,41,1),
                         title= ['Cloud Fraction (%) at '+site+' - Obs', 'Cloud Fraction (%) at '+site+' - E3SMv1'])
fig.savefig(figpath+'seasonalcycle_cloud2d_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   
