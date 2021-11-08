
import sys
sys.path.insert(1,'../subroutines/')

import numpy as np
from time_format_change import  yyyymmdd2cday, cday2mmdd
from read_netcdf import read_E3SM
from netCDF4 import Dataset


#%% settings

dateall = ['2015-10','2015-11','2015-12','2016-01','2016-02','2016-03','2016-04','2016-05','2016-06','2016-07',\
           '2016-08','2016-09','2016-10','2016-11','2016-12','2017-01','2017-02','2017-03','2017-04','2017-05',\
           '2017-06','2017-07','2017-08','2017-09','2017-10','2017-11','2017-12','2018-01','2018-02','2018-03']

E3SM_input_path =  '/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_h0/'   
E3SM_output_path = '../../figures/'
model='E3SMv1'

#%%  process data for each day
for date in dateall:
    
    print(date)
    filename_input = E3SM_input_path+'E3SMv1_2014-2018.cam.h0.'+date+'.nc'
    
    if date=='2015-10':
        (timem,lonm,timeunitm,lonmunit,lonmname)=read_E3SM(filename_input,'lon')
        (timem,latm,timeunitm,latmunit,latmname)=read_E3SM(filename_input,'lat')
        aodall = np.empty((0,len(lonm)))
    (timem,aod,timeunitm,aodunit,aodname)=read_E3SM(filename_input,'AODVIS')
    aod[aod>1e10]=np.nan
    aodall = np.vstack((aodall,np.ma.filled(aod,np.nan)))
    
aodmean = np.nanmean(aodall,0)
if len(aodmean)!=len(lonm):
    print(aodall.shape)
    error
    
    
# %% output extacted file
outputname = 'AOD_mean_global_'+model+'.nc'
print('output to this file: '+E3SM_output_path+outputname)

# define filename
f = Dataset(E3SM_output_path+outputname, 'w', format='NETCDF4')

# define dimensions
t = f.createDimension('ncol', None)  # unlimited

# create variable list
lat_o = f.createVariable("lat","f8",("ncol",))
lon_o = f.createVariable("lon","f8",("ncol",))
var_o = f.createVariable("AODmean",'f8',("ncol",))

# write data
lat_o[:] = latm
lon_o[:] = lonm
var_o[:] = aodmean

# attributes
lat_o.units = latmunit
lon_o.units = lonmunit
var_o.units = aodunit
var_o.long_name = aodname

# global attributes
import time as ttt
f.description = model+" extact global mean AOD"
f.history = "Created by Shuaiqi at " + ttt.ctime(ttt.time())

f.close()
    

        
        
        
