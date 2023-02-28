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
import matplotlib.pyplot as plt

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'
# path of prepared files
prep_model_path = '../prep_data/' +site+'/model/'
prep_flight_path = '../prep_data/'+site+'/flight/'
# set output path for plots
figpath= '../figures/'

obsdata = xr.open_mfdataset(prep_flight_path + 'CCN_'+site+'*.nc')
time_ccn = obsdata['time']
ccn2 = obsdata['CCN2'].load()
ccn5 = obsdata['CCN5'].load()
# ss1 = obsdata['SS2'].load()
# ss5 = obsdata['SS5'].load()
obsdata.close()

obsdata = xr.open_mfdataset(prep_flight_path + 'PCASP100_'+site+'*.nc')
time_pcasp = obsdata['time']
cn100 = obsdata['pcasp100'].load()
obsdata.close()

lst = glob.glob(prep_flight_path + 'mergedSD_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_nd = obsdata['time'].load()
nd = obsdata['Nd'].load()
obsdata.close()
nd = nd*0.001  # #/L to #/cm3

# E3SM data
lst = prep_model_path + 'E3SMv2_'+site+'_flight_*.nc'
modeldata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_m = modeldata['time'].load()
ccn1_m = modeldata['CCN3'].load()
ccn2_m = modeldata['CCN4'].load()
cn100_m = modeldata['NCN100'].load()
nd_m = modeldata['ICWNC'].load()
nd_bin_m = modeldata['Nd_bin'].load()
modeldata.close()
nd_all_m = np.sum(nd_bin_m,axis=0)

# lst = prep_model_path + 'E3SMv2_'+site+'_flight_*.nc'
# modeldata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
# time_m2 = modeldata['time'].load()
# ccn1_m2 = modeldata['CCN3'].load()
# ccn2_m2 = modeldata['CCN4'].load()
# cn100_m2 = modeldata['NCN100'].load()
# nd_m2 = modeldata['ICWNC'].load()
# nd_bin_m2 = modeldata['Nd_bin'].load()
# modeldata.close()
# nd_all_m2 = np.sum(nd_bin_m2,axis=0)


#%% specific data treatments

ccn2[ccn2==0] = np.nan
cn100 = cn100.interp(time=time_ccn)
nd = nd.interp(time=time_ccn)
ccn1_m = ccn1_m.interp(time=time_ccn)
ccn2_m = ccn2_m.interp(time=time_ccn)
cn100_m = cn100_m.interp(time=time_ccn)
nd_m = nd_m.interp(time=time_ccn)
# ccn1_m2 = ccn1_m2.interp(time=time_ccn)
# ccn2_m2 = ccn2_m2.interp(time=time_ccn)
# cn100_m2 = cn100_m2.interp(time=time_ccn)
# nd_m2 = nd_m2.interp(time=time_ccn)


nd[nd<10] = np.nan
nd[nd>800] = np.nan
nd_m[nd_m<10] = np.nan
nd_m[nd_m>800] = np.nan
# nd_all_m[nd_all_m<10] = np.nan
# nd_all_m[nd_all_m>800] = np.nan
# nd_m2[nd_m2<10] = np.nan
# nd_m2[nd_m2>800] = np.nan
# nd_all_m2[nd_all_m2<10] = np.nan
# nd_all_m2[nd_all_m2>800] = np.nan

print(np.sum(~np.isnan(nd.data)))
print(np.sum(~np.isnan(nd_m.data)))
# print(np.sum(~np.isnan(nd_m2.data)))


#%%
# w0 = np.ones_like(nd[:])/sum(~np.isnan(nd[:].data))
# w1 = np.ones_like(nd_m[:])/sum(~np.isnan(nd_m[:].data))
# w2 = np.ones_like(nd_m2[:])/sum(~np.isnan(nd_m2[:].data))
# fig,ax = plot.hist([nd[:],nd_m[:],nd_m2[:]], 
#                     weights=[w0,w1,w2], bins=np.arange(0,810,20), 
#                     legend = ['Ndrop','E3SMv2','E3SMv2'], color=['k','r','b'],
#                     title = 'Nd '+site, ylabel='Fraction', xlabel='cm$^{-3}$')

# w0 = np.ones_like(ccn2[:])/sum(~np.isnan(ccn2[:].data))
# w1 = np.ones_like(ccn2_m[:])/sum(~np.isnan(ccn2_m[:].data))
# w2 = np.ones_like(ccn2_m2[:])/sum(~np.isnan(ccn2_m2[:].data))
# fig,ax = plot.hist([ccn2[:],ccn2_m[:],ccn2_m2[:]], 
#                     weights=[w0,w1,w2], bins=np.arange(0,1110,30), 
#                     legend = ['obs','E3SMv2','E3SMv2'], color=['k','r','b'],
#                     title = '0.2% CCN '+site, ylabel='Fraction', xlabel='cm$^{-3}$')

# w0 = np.ones_like(cn100[:])/sum(~np.isnan(cn100[:].data))
# w1 = np.ones_like(cn100_m[:])/sum(~np.isnan(cn100_m[:].data))
# w2 = np.ones_like(cn100_m2[:])/sum(~np.isnan(cn100_m2[:].data))
# fig,ax = plot.hist([cn100[:],cn100_m[:],cn100_m2[:]], 
#                     weights=[w0,w1,w2], bins=np.arange(0,1110,30), 
#                     legend = ['obs','E3SMv2','E3SMv2'], color=['k','r','b'],
#                     title = 'CN (>100nm) '+site, ylabel='Fraction', xlabel='cm$^{-3}$')

#%%
# Nd vs surface CCN
fig,ax = plot.scatter([np.log(ccn2[:].data), np.log(ccn2_m[:].data), ], 
                        [np.log(nd[:].data), np.log(nd_m[:].data),], 
                    linear_fit=True, intercept=True, xlimit=(1,7), ylimit=(1,7), 
                    xlabel='0.2%CCN (cm$^{-3}$)', 
                    ylabel='Nd (cm$^{-3}$)', title=['Aircraft','E3SMv2',], )
for ax0 in ax:
    ax0.set_xlim(np.log(10), np.log(800))
    ax0.set_ylim(np.log(10), np.log(800))
    ax0.set_xticks(np.log([10,30,100,300]))
    ax0.set_xticklabels((10, 30,100, 300))
    ax0.set_yticks(np.log([10,30,100,300]))
    ax0.set_yticklabels((10, 30,100, 300))
# fig,ax = plot.scatter([ccn2[:].data, ccn2_m[:].data, ccn2_m2[:].data], 
#                         [nd[:].data, nd_m[:].data, nd_m2[:].data], figsize=(18,5),
#                     linear_fit=True, intercept=False,
#                     xlimit=(0,800), ylimit=(0,800), xlabel='0.2%CCN (cm$^{-3}$)', 
#                     ylabel='Nd (cm$^{-3}$)', title=['Obs','E3SMv2','E3SMv2'], )
# for ax0 in ax:
#     ax0.set_xscale('log')
#     ax0.set_yscale('log')
#     ax0.set_xlim((10, 300))
#     ax0.set_ylim((10, 300))
regress,sample = calc.linear_regress([np.log(ccn2[:].data), np.log(ccn2_m[:].data), ], 
                        [np.log(nd[:].data), np.log(nd_m[:].data), ],
                        figpath+'linearfit_Nd_CCN_aircraft_'+site+'.txt',legend=['Surface','E3SMv2',],
                        labelx='ln(CCN2)',labely='ln(Nd)')
for nn in range(len(ax)):
    ax[nn].text(np.log(25),np.log(400), 'R = '+format(regress[nn][2],'3.2f'))
    print(regress[nn][3])


# fig,ax = plot.jointhist([ccn2[:].data, ccn2_m[:].data, ccn2_m2[:].data], 
#                         [nd[:].data, nd_m[:].data, nd_m2[:].data],
#                         xedges=np.exp(np.arange(np.log(10),7,0.2)), yedges=np.exp(np.arange(np.log(10),7,0.2)), 
#                         normalize_x=True, vmin=0, vmax=0.3,
#                     xlabel=['0.2%CCN (cm$^{-3}$)','0.2%CCN (cm$^{-3}$)','0.2%CCN (cm$^{-3}$)'], 
#                     ylabel='Nd (cm$^{-3}$)', title=['Obs','E3SMv2','E3SMv2'], )
# for ax0 in ax[0,:]:
#     ax0.set_xscale('log')
#     ax0.set_yscale('log')
#     ax0.set_xticks([10,30,100,300,1000])
#     ax0.set_yticks([10,30,100,300,1000])
#     ax0.set_yticklabels([10,30,100,300,1000])
# for ax0 in ax[1,:]:
#     ax0.set_xscale('log')
#     ax0.set_xticks([10,30,100,300,1000])
#     ax0.set_xticklabels([10,30,100,300,1000])



# regress,sample = calc.linear_regress([np.log(ccn2[:].data), np.log(ccn2_m[:].data), np.log(ccn2_m2[:].data)], 
#                         [np.log(nd[:].data), np.log(nd_m[:].data), np.log(nd_m2[:].data)],
#                         figpath+'linearfit_Nd_CCN_aircraft_'+site+'.txt',legend=['Surface','E3SMv2','E3SMv2'],
#                         labelx='ln(CCN2)',labely='ln(Nd)')
# x = np.arange(20,700,50)
# for nn in range(3):
#     y = np.exp(regress[nn][0] * np.log(x) + regress[nn][1])
#     ax[0,nn].plot(x, y, color='r')
