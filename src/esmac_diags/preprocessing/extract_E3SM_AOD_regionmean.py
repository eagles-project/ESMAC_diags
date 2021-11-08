
import sys

import numpy as np
from ..subroutines.time_format_change import  yyyymmdd2cday, cday2mmdd
from read_netcdf import read_E3SM
from netCDF4 import Dataset


#%% settings

from settings import campaign, E3SM_h3_path, E3SM_h3_filehead, Model_List

E3SM_region_path = '../../figures/'

    
if campaign=='HISCALE':
    E3SMdomain_range='260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    start_date='2016-04-25'
    end_date='2016-09-22'
elif campaign=='ACEENA':
    E3SMdomain_range='330e_to_335e_37n_to_42n'   
    start_date='2017-06-20'
    end_date='2018-02-19'
elif campaign=='CSET':
    E3SMdomain_range='202e_to_240e_19n_to_40n' 
    start_date='2015-07-01'
    end_date='2015-08-15'   
elif campaign=='SOCRATES':
    E3SMdomain_range='133e_to_164e_42s_to_63s'    
    start_date='2018-01-15'
    end_date='2018-02-24'
elif campaign=='MAGIC':
    E3SMdomain_range='202e_to_243e_20n_to_35n'    # domain range in E3SM regional output 
    start_date='2012-10-01'
    end_date='2013-09-30'
elif campaign=='MARCUS':
    E3SMdomain_range='60e_to_160e_42s_to_70s'    
    start_date='2017-10-01'
    end_date='2018-04-20'
else:
    print('ERROR: please specify domain info for '+campaign)
    error

# change start date into calendar day
cday1 = yyyymmdd2cday(start_date,'noleap')
cday2 = yyyymmdd2cday(end_date,'noleap')

year0 = start_date[0:4]
if start_date[0:4]!=end_date[0:4]:
    cday2=cday2+365

for mm in range(len(Model_List)):
    model=Model_List[mm]
    #%%  process data for each day
    for cday in range(cday1,cday2+1):
        if cday>365:
            mmdd=cday2mmdd(cday-365)
            date=end_date[0:4]+'-'+mmdd[0:2]+'-'+mmdd[2:4]
        else:
            mmdd=cday2mmdd(cday)
            date=year0+'-'+mmdd[0:2]+'-'+mmdd[2:4]
        
        print(date)
        filename_input = E3SM_h3_path[mm]+E3SM_h3_filehead[mm]+'.cam.h3.'+date+'-00000.nc'
        
        if cday==cday1:
            (timem,lonm,timeunitm,lonmunit,lonmname)=read_E3SM(filename_input,'lon_'+E3SMdomain_range)
            (timem,latm,timeunitm,latmunit,latmname)=read_E3SM(filename_input,'lat_'+E3SMdomain_range)
            aodall = np.empty((0,len(lonm)))
        (timem,aod,timeunitm,aodunit,aodname)=read_E3SM(filename_input,'AODVIS_'+E3SMdomain_range)
        aod[aod>1e10]=np.nan
        aodall = np.vstack((aodall,np.ma.filled(aod,np.nan)))
        
    aodmean = np.nanmean(aodall,0)
    if len(aodmean)!=len(lonm):
        print(aodall.shape)
        error
        
        
    # %% output extacted file
    outputname = 'AOD_mean_'+campaign+'_'+model+'.nc'
    print('output to this file: '+E3SM_region_path+outputname)
    
    # define filename
    f = Dataset(E3SM_region_path+outputname, 'w', format='NETCDF4')
    
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
    f.description = model+" extact mean AOD during "+campaign
    f.history = "Created by Shuaiqi at " + ttt.ctime(ttt.time())
    
    f.close()
    

        
        
        
