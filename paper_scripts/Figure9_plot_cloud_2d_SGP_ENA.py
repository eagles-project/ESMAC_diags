"""
example to generate 2-d histogram plot for cloud fraction

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
site = 'ENA'


# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/surface/'
# set output path for plots
figpath= '../figures/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = glob.glob(prep_obs_path + 'cloud_2d_'+site+'*.nc')
obsdata = xr.open_mfdataset(lst)
time_cf = obsdata['time'].load()
height1 = obsdata['height'].load()
cf_obs = obsdata['cloud'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv2_'+site+'_profiles.nc'
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
#%% diurnal and seasonal cycles

fig,ax = plot.diurnalcycle_2d([cf_obs, cf_e3sm, ], y = [height1, heightm, ],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)', cmap='hot_r',
                        levellist=np.arange(0,29,1),
                         title= ['Obs', 'E3SMv2', ])
# fig.savefig(figpath+'diurnalcycle_cloud2d_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.seasonalcycle_2d([cf_obs, cf_e3sm, ], y = [height1, heightm, ],
                        yticks=[0,3,6,9,12], ylimit=(0,12), ylabel='Height (km)', cmap='hot_r',
                        levellist=np.arange(0,41,1),
                        # title= ['Obs', 'E3SMv2',]
                         )
# fig.savefig(figpath+'seasonalcycle_cloud2d_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
    
    
#%% histogram    
xdata = [cf_obs.data.flatten(), cf_e3sm.data.flatten()]
ndata = len(xdata)
h = np.repeat(np.reshape(height1.data,(1,-1)),len(cf_obs.time),axis=0)  # height
ydata = [h.flatten() for i in range(ndata)]
xedges = np.arange(0.,101.,5.)
xedges[0]=0.1
yedges = np.array([0,.2,.4,.6,.8,1,1.5,2,2.5,3,4,5,6,7,8,10,12])
# yedges = np.arange(18)
weights = []
for mm in range(ndata):
        hnum,x,y=np.histogram2d((xdata[mm]*0.0),ydata[mm],bins=[np.array([-1,1]), yedges])
        weights.append(1./hnum.T)

if site=='SGP':
    fig,ax = plot.jointhist(xdata, ydata, xedges=xedges, yedges=yedges, weight=weights,
                vmin=0.0, vmax=0.11, cmap='viridis',# normalize_x=False,
                xlabel='Cloud Fraction (%)', ylabel='Height (km)', xlimit=(0,100),
               # title=['Obs', 'E3SMv2',]
                )
    # fig.savefig(figpath+'histogram_cloud_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

if site=='ENA':
    fig,ax = plot.jointhist(xdata, ydata, xedges=xedges, yedges=yedges, weight=weights,
                vmin=0.0, vmax=0.22, cmap='viridis',# normalize_x=False,
                xlabel='Cloud Fraction (%)', ylabel='Height (km)', xlimit=(0,100),
              #  title=['Obs','E3SMv2',]
                )
    # fig.savefig(figpath+'histogram_cloud_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
