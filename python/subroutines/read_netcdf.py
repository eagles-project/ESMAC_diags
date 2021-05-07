# -*- coding: utf-8 -*-

#%% read variables from E3SM output related to height
# filename='../data/MAM/cutoff_number_MAM5.nc'
# varname='T'
def read_E3SM_z(filename,varname):
    
    from netCDF4 import Dataset
    
    f = Dataset(filename,'r')
    
    # read in variables
    height = f.variables['height'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    
    f.close()
    
    return(time,height,data,timeunit,dataunit,long_name)


#%% E3SM data without Z3
def read_E3SM(filename,varname):
    
    from netCDF4 import Dataset
    
    f = Dataset(filename,'r')
    
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    if type(varname) is str:
        d_id = f.variables[varname]
        data = d_id[:]
        dataunit = d_id.units
        long_name = d_id.long_name
    elif type(varname) is list:
        data=list()
        dataunit=list()
        long_name=list()
        for vv in range(len(varname)):
            d_id=f.variables[varname[vv]]
            data.append(d_id[:])
            dataunit.append(d_id.units)
            long_name.append(d_id.long_name)
    
    f.close()
    
    return(time,data,timeunit,dataunit,long_name)

#%%
# filename='../data/MAM/NCN_MAM5_20160907_R1_L2.nc'
# varname='NCN'
def read_extractflight(filename,varname):
    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    height = f.variables['height'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    
    f.close()
    return(time,height,data,timeunit,dataunit,long_name)


#%% read merged size distribution data
def read_merged_size(filename,varname):
    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    size = f.variables['size'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    
    f.close()
    return(time,size,data,timeunit,dataunit,long_name)