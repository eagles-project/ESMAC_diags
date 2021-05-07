# calculate statistics (mean, bias, correlation, RMSE) of Aerosol number concentration
# for aircraft measurements
# compare models and CPC measurements


import sys
sys.path.insert(1,'../subroutines/')

# import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import glob
from read_aircraft import read_cpc
from read_netcdf import read_merged_size,read_extractflight

#%% settings

from settings import campaign, cpcpath,merged_size_path, Model_List, IOP, \
    E3SM_aircraft_path, figpath_aircraft_statistics

import os
if not os.path.exists(figpath_aircraft_statistics):
    os.makedirs(figpath_aircraft_statistics)
   
missing_value = -999999.

#%% find files for flight information

lst = glob.glob(merged_size_path+'merged_bin_*'+campaign+'*.nc')
lst.sort()

if len(lst)==0:
    print('ERROR: cannot find any file at '+merged_size_path)
    error

# choose files for specific IOP
if campaign=='HiScale':
    if IOP=='IOP1':
        lst=lst[0:17]
    elif IOP=='IOP2':
        lst=lst[17:]
    elif IOP[0:4]=='2016':
        a=lst[0].split('_'+campaign+'_')
        lst = glob.glob(a[0]+'*'+IOP+'*')
        lst.sort()
elif campaign=='ACEENA':
    if IOP=='IOP1':
        lst=lst[0:20]
    elif IOP=='IOP2':
        lst=lst[20:]
    elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
        a=lst[0].split('_'+campaign+'_')
        lst = glob.glob(a[0]+'*'+IOP+'*')
        lst.sort()
else:
    print('ERROR: campaign name is not recognized: '+campaign)
    error

if len(lst)==0:
    print('ERROR: cannot find any file for '+IOP)
    error
    
#%% read all data

cpc10_o = np.empty(0)
cpc3_o = np.empty(0)
cpc10_m = []
cpc3_m = []
nmodels=len(Model_List)
for mm in range(nmodels):
    cpc10_m.append(np.empty(0))
    cpc3_m.append(np.empty(0))
    
print('reading '+format(len(lst))+' files to calculate the statistics: ')

for filename in lst:
    
    # get date info:        
    date=filename[-12:-3]
    if date[-1]=='a':
        flightidx=1
    else:
        flightidx=2
    print(date)
    
    #% read in flight information
    (time,size,cvi,timeunit,cunit,long_name)=read_merged_size(filename,'CVI_inlet')
    (time,size,cflag,timeunit,cunit,long_name)=read_merged_size(filename,'cld_flag')
    (time,size,height,timeunit,zunit,long_name)=read_merged_size(filename,'height')
    time=np.ma.compressed(time)
    
    
    #%% read in CPC measurements
    
    if campaign=='HiScale':
        filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
    elif campaign=='ACEENA':
        filename_c=glob.glob(cpcpath+'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')    
    filename_c.sort()
    # read in data
    if len(filename_c)==1 or len(filename_c)==2: # some days have two flights
        (cpc,cpclist)=read_cpc(filename_c[flightidx-1])
        if np.logical_and(campaign=='ACEENA', date=='20180216a'):
            cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1)
        elif np.logical_and(campaign=='HiScale', date=='20160425a'):
            cpc=np.insert(cpc,0,cpc[:,0],axis=1)
            cpc[0,0]=cpc[0,0]-1
        time_cpc = cpc[0,:]
        cpc10 = cpc[1,:]
        cpc3 = cpc[2,:]
    elif len(filename_c)==0:
        time_cpc=time
        cpc10=np.nan*np.empty([len(time)])
        cpc3=np.nan*np.empty([len(time)])
    else:
        print('find too many files, check: ')
        print(filename_c)
        error
    
    # some quality checks
    cpc3[cpc3<100]=np.nan
    cpc10[cpc10<10]=np.nan
    
    cpc10_o=np.hstack((cpc10_o, cpc10))
    cpc3_o=np.hstack((cpc3_o, cpc3))
    
    #%% read in Models
    for mm in range(nmodels):
        filename_m = E3SM_aircraft_path+'Aircraft_vars_'+campaign+'_'+Model_List[mm]+'_'+date+'.nc'
    
        (timem,heightm,cpc_m,timeunitm,ncn_unit,ncn_longname)=read_extractflight(filename_m,'NCN')
        (timem,heightm,cpcu_m,timeunitm,ncnu_unit,ncnu_longname)=read_extractflight(filename_m,'NUCN')
        if len(cpc10)!=len(cpc_m):
            print('CPC and MAM have different dimensions! check')
            print(cpc10.shape,cpc_m.shape)
            errors
        
        cpc10_m[mm] = np.hstack((cpc10_m[mm], cpc_m*1e-6))    # #/m3 to #/cm3
        cpc3_m[mm] = np.hstack((cpc3_m[mm], cpcu_m*1e-6))    # #/m3 to #/cm3
    
#%% calculate statistics

# select only valid data in obs and the corresponding data in models
idx10 = ~np.isnan(cpc10_o)
idx3 = ~np.isnan(cpc3_o)

mean10 = [None]*(nmodels+1)
mean3 = [None]*(nmodels+1)
bias10 = [None]*(nmodels)
bias3 = [None]*(nmodels)
corr10 = [None]*(nmodels)
corr3 = [None]*(nmodels)
rmse10 = [None]*(nmodels)
rmse3 = [None]*(nmodels)
    
if sum(idx10)/len(idx10)<0.1:   # two few observation available
    # for obs
    mean10[nmodels] = missing_value
    mean3[nmodels] = missing_value
    # for models
    for mm in range(nmodels):
        mean10[mm] = np.nanmean(cpc10_m[mm][idx10])
        mean3[mm] = np.nanmean(cpc3_m[mm][idx3])
        bias10[mm] =  missing_value
        bias3[mm] =  missing_value
        corr10[mm] = [missing_value, missing_value]
        corr3[mm] = [missing_value, missing_value]
        rmse10[mm] = missing_value
        rmse3[mm] = missing_value
else:
    # for obs
    mean10[nmodels] = np.nanmean(cpc10_o[idx10])
    mean3[nmodels] = np.nanmean(cpc3_o[idx3])
    # for models
    for mm in range(nmodels):
        mean10[mm] = np.nanmean(cpc10_m[mm][idx10])
        mean3[mm] = np.nanmean(cpc3_m[mm][idx3])
        bias10[mm] = mean10[mm] - mean10[nmodels]
        bias3[mm] = mean3[mm] - mean3[nmodels]
        c10 = scipy.stats.pearsonr(cpc10_m[mm][idx10],cpc10_o[idx10])
        c3 = scipy.stats.pearsonr(cpc3_m[mm][idx3],cpc3_o[idx3])
        corr10[mm] = [c10[0],c10[1]]
        corr3[mm] = [c3[0],c3[1]]
        rmse10[mm] = np.sqrt(((cpc10_m[mm][idx10]-cpc10_o[idx10])**2).mean())
        rmse3[mm] = np.sqrt(((cpc3_m[mm][idx3]-cpc3_o[idx3])**2).mean())

#%% write out files

outfile = figpath_aircraft_statistics+'statistics_CN10nm_'+campaign+'_'+IOP+'.txt'
print('write statistics to file '+outfile)

with open(outfile, 'w') as f:
    f.write('statistics of Aerosol Number Concentration comparing with CPC(>10nm). sample size '+format(sum(idx10))+'\n')
    line1 = list(Model_List)
    line1.insert(0,' --- ')
    line1.append('OBS')
    for ii in range(len(line1)):
        f.write(format(line1[ii],'10s')+', ')
    # write mean
    f.write('\n mean,\t')
    for ii in range(len(mean10)):
        f.write(format(mean10[ii],'10.2f')+', ')
    # write bias
    f.write('\n bias,\t')
    for ii in range(len(bias10)):
        f.write(format(bias10[ii],'10.2f')+', ')
    # write rmse
    f.write('\n RMSE,\t')
    for ii in range(len(rmse10)):
        f.write(format(rmse10[ii],'10.2f')+', ')
    # write correlation
    f.write('\n corrcoef,\t')
    for ii in range(len(rmse10)):
        f.write(format(corr10[ii][0],'10.4f')+', ')
    # write p value of correlation
    f.write('\n P_corr,\t')
    for ii in range(len(rmse10)):
        f.write(format(corr10[ii][1],'10.2f')+', ')
        
        
outfile = figpath_aircraft_statistics+'statistics_CN3nm_'+campaign+'_'+IOP+'.txt'
print('write statistics to file '+outfile)

with open(outfile, 'w') as f:
    f.write('statistics of Aerosol Number Concentration comparing with CPC(>3nm). sample size '+format(sum(idx3))+'\n')
    line1 = list(Model_List)
    line1.insert(0,' --- ')
    line1.append('OBS')
    for ii in range(len(line1)):
        f.write(format(line1[ii],'10s')+', ')
    # write mean
    f.write('\n mean,\t')
    for ii in range(len(mean3)):
        f.write(format(mean3[ii],'10.2f')+', ')
    # write bias
    f.write('\n bias,\t')
    for ii in range(len(bias3)):
        f.write(format(bias3[ii],'10.2f')+', ')
    # write rmse
    f.write('\n RMSE,\t')
    for ii in range(len(rmse3)):
        f.write(format(rmse3[ii],'10.2f')+', ')
    # write correlation
    f.write('\n corrcoef,\t')
    for ii in range(len(rmse3)):
        f.write(format(corr3[ii][0],'10.4f')+', ')
    # write p value of correlation
    f.write('\n P_corr,\t')
    for ii in range(len(rmse3)):
        f.write(format(corr3[ii][1],'10.2f')+', ')