
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set output path for plots
figpath= '../figures/'

#%% read variables in different field campaigns
# set site name.
site = 'HISCALE'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/flight/'

lst = glob.glob(prep_obs_path + 'CPC_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_hiscale = obsdata['time'].load()
cpc10_hiscale = obsdata['cpc10'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'CCN_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_hiscale = obsdata['time'].load()
ccn2_hiscale = obsdata['CCN2'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'PCASP100_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_hiscale = obsdata['time'].load()
pcasp100_hiscale = obsdata['pcasp100'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'WCM_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_hiscale = obsdata['time'].load()
lwc_hiscale = obsdata['LWC'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'mergedSD_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_hiscale = obsdata['time'].load()
nd_hiscale = obsdata['Nd'].load()
obsdata.close()

lst2 = glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc')
modeldata = xr.open_mfdataset(lst2)
timem_hiscale = modeldata['time'].load()
ccn2_m_hiscale = modeldata['CCN4'].load()
lwc_m_hiscale = modeldata['LWC'].load()
nd_m_hiscale = modeldata['ICWNC'].load()
ncn10_m_hiscale = modeldata['NCN10'].load()
ncn100_m_hiscale = modeldata['NCN100'].load()
modeldata.close()

# set site name.
site = 'ACEENA'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/flight/'

lst = glob.glob(prep_obs_path + 'CPC_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_aceena = obsdata['time'].load()
cpc10_aceena = obsdata['cpc10'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'CCN_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_aceena = obsdata['time'].load()
ccn1_aceena = obsdata['CCN1'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'PCASP100_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_aceena = obsdata['time'].load()
pcasp100_aceena = obsdata['pcasp100'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'WCM_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_aceena = obsdata['time'].load()
lwc_aceena = obsdata['LWC'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'mergedSD_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_aceena = obsdata['time'].load()
nd_aceena = obsdata['Nd'].load()
obsdata.close()

lst2 = glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc')
modeldata = xr.open_mfdataset(lst2)
timem_aceena = modeldata['time'].load()
ccn1_m_aceena = modeldata['CCN3'].load()
lwc_m_aceena = modeldata['LWC'].load()
nd_m_aceena = modeldata['ICWNC'].load()
ncn10_m_aceena = modeldata['NCN10'].load()
ncn100_m_aceena = modeldata['NCN100'].load()
modeldata.close()

# set site name.
site = 'CSET'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/flight/'

lst = glob.glob(prep_obs_path + 'CN_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_cset = obsdata['time'].load()
cpc10_cset = obsdata['CPC10'].load()
uhsas100_cset = obsdata['UHSAS100'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'LWC_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_cset = obsdata['time'].load()
lwc_cset = obsdata['LWC'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'Nd_size_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_cset = obsdata['time'].load()
nd_cset = obsdata['Nd'].load()
obsdata.close()

lst1 = glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc')
modeldata = xr.open_mfdataset(lst1)
timem_cset = modeldata['time'].load()
ccn1_m_cset = modeldata['CCN3'].load()
lwc_m_cset = modeldata['LWC'].load()
nd_m_cset = modeldata['ICWNC'].load()
ncn10_m_cset = modeldata['NCN10'].load()
ncn100_m_cset = modeldata['NCN100'].load()
modeldata.close()

# set site name.
site = 'SOCRATES'
# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/flight/'

lst = glob.glob(prep_obs_path + 'CN_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
cpc10_socrates = obsdata['CPC10'].load()
uhsas100_socrates = obsdata['UHSAS100'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'CCN_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
ccn1_socrates = obsdata['CCN1_spec'].load()
ccn1_socrates_2 = obsdata['CCN1_scan'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'LWC_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_socrates = obsdata['time'].load()
lwc_socrates = obsdata['LWC'].load()
obsdata.close()
lst = glob.glob(prep_obs_path + 'Nd_size_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_socrates = obsdata['time'].load()
nd_socrates = obsdata['Nd'].load()
obsdata.close()

lst1 = glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc')
modeldata = xr.open_mfdataset(lst1)
timem_socrates = modeldata['time'].load()
ccn1_m_socrates = modeldata['CCN3'].load()
lwc_m_socrates = modeldata['LWC'].load()
nd_m_socrates = modeldata['ICWNC'].load()
ncn10_m_socrates = modeldata['NCN10'].load()
ncn100_m_socrates = modeldata['NCN100'].load()
modeldata.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment
nd_hiscale = nd_hiscale*0.001   # #/L to #/cm3
nd_aceena = nd_aceena*0.001   # #/L to #/cm3
nd_cset = nd_cset*0.001   # #/L to #/cm3
nd_socrates = nd_socrates*0.001   # #/L to #/cm3

nd_hiscale[nd_hiscale<10] = np.nan
nd_m_hiscale[nd_m_hiscale<10] = np.nan
nd_aceena[nd_aceena<10] = np.nan
nd_m_aceena[nd_m_aceena<10] = np.nan
nd_cset[nd_cset<10] = np.nan
nd_m_cset[nd_m_cset<10] = np.nan
nd_socrates[nd_socrates<10] = np.nan
nd_m_socrates[nd_m_socrates<10] = np.nan

lwc_hiscale[lwc_hiscale<0.02] = np.nan
lwc_m_hiscale[lwc_m_hiscale<0.02] = np.nan
lwc_aceena[lwc_aceena<0.02] = np.nan
lwc_m_aceena[lwc_m_aceena<0.02] = np.nan
lwc_cset[lwc_cset<0.02] = np.nan
lwc_m_cset[lwc_m_cset<0.02] = np.nan
lwc_socrates[lwc_socrates<0.02] = np.nan
lwc_m_socrates[lwc_m_socrates<0.02] = np.nan

# cpc10_hiscale[cpc10_hiscale>2e4] = np.nan
# cpc10_aceena[cpc10_aceena>1e4] = np.nan
# cpc10_cset[cpc10_cset>1e4] = np.nan
# cpc10_socrates[cpc10_socrates>1e4] = np.nan

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
fig,ax = plot.percentiles([ccn2_hiscale, ccn2_m_hiscale], 
                          [ccn1_aceena, ccn1_m_aceena, ],
                          [ccn1_m_cset*np.nan, ccn1_m_cset, ],
                          [ccn1_socrates, ccn1_m_socrates,],
                    title='CCN (cm$^{-3}$)', ylimit=(0,1000), figsize=(10,2.5),
                    xlabel=['HI-SCALE (0.2%)','ACE-ENA (0.1%)','CSET (0.1%)','SOCRATES (0.1%)'], 
                    color=['k','r',], legend=None)
# fig.savefig(figpath+'percentiles_CCN_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([cpc10_hiscale, ncn10_m_hiscale, ], 
                          [cpc10_aceena, ncn10_m_aceena, ],
                          [cpc10_cset, ncn10_m_cset, ],
                          [cpc10_socrates, ncn10_m_socrates, ],
                    title='(a) CN (>10nm) (cm$^{-3}$)', ylimit=(0,8000),figsize=(10,2.5),
                    xlabel=['HI-SCALE','ACE-ENA','CSET','SOCRATES'], 
                    color=['k','r',], legend=['Obs','E3SMv2'])
# fig.savefig(figpath+'percentiles_CN10_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([pcasp100_hiscale, ncn100_m_hiscale, ], 
                          [pcasp100_aceena, ncn100_m_aceena, ],
                          [uhsas100_cset, ncn100_m_cset, ],
                          [uhsas100_socrates, ncn100_m_socrates, ],
                    title='(b) CN (>100nm) (cm$^{-3}$)',ylimit=(0,1200),figsize=(10,2.5),
                    xlabel=['HI-SCALE','ACE-ENA','CSET','SOCRATES'], 
                    color=['k','r',], legend=None)
# fig.savefig(figpath+'percentiles_CN100_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([nd_hiscale, nd_m_hiscale, ], 
                          [nd_aceena, nd_m_aceena, ],
                          [nd_cset, nd_m_cset, ],
                          [nd_socrates, nd_m_socrates, ],
                    title='(c) Nd (cm$^{-3}$)',figsize=(10,2.5),
                    xlabel=['HI-SCALE','ACE-ENA','CSET','SOCRATES'], 
                    color=['k','r',], legend=None)
# fig.savefig(figpath+'percentiles_Nd_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.percentiles([lwc_hiscale, lwc_m_hiscale, ], 
                          [lwc_aceena, lwc_m_aceena, ],
                          [lwc_cset, lwc_m_cset, ],
                          [lwc_socrates, lwc_m_socrates, ],
                    title='(d) LWC (g/m$^3$)',figsize=(10,2.5),
                    xlabel=['HI-SCALE','ACE-ENA','CSET','SOCRATES'], 
                    color=['k','r'], legend=None)
# fig.savefig(figpath+'percentiles_LWC_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
