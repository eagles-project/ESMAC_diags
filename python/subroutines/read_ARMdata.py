# read ARM standard NETCDF files

#%% composition in acsm
# filename='../data/arm-cpcu/iop1/sgpaoscpcuS01.b1.20160422.000000.nc'
def read_acsm(filename,varname):
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    # flag = f.variables['qc_'+varname][:]
    # data[flag!=0]=np.nan
    f.close()
    return(time,data,timeunit,dataunit)

#%% armbecldrad data
# filename='../../data/ACEENA/obs/profile/enaarmbecldrad/enaarmbecldradC1.c1.20170101.003000.nc'
# varname='cld_frac'
def read_armbe(filename,varname):
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    height = f.variables['height'][:]
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    # flag = f.variables['qc_'+varname][:]
    # data[flag!=0]=np.nan
    f.close()
    return(time,height,data,timeunit,dataunit)

#%%  CCN
# filename='../data/ccn_aaf/enaaafccn2colaF1.b1.20180121.094335.nc'

def read_ccn(filename):    

    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    # lon = f.variables['lon'][:]
    # lat = f.variables['lat'][:]
    # alt = f.variables['alt'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    SS = f.variables['supersaturation_calculated'][:]
    d_id = f.variables['N_CCN']
    qc_ccn = f.variables['qc_N_CCN'][:]
    ccn = d_id[:]
    dataunit = d_id.units
    f.close()
    ccn[qc_ccn!=0] = -9999.
    return(time,timeunit,ccn,dataunit,SS)

# MAGIC
def read_ccn_magic(filename):    

    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    SS = f.variables['CCN_ss_set'][:]
    d_id = f.variables['N_CCN']
    ccn = d_id[:]
    dataunit = d_id.units
    f.close()
    return(time,timeunit,ccn,dataunit,SS)

#%% CPC and CPCu
# filename='../data/arm-cpcu/iop1/sgpaoscpcuS01.b1.20160422.000000.nc'
def read_cpc(filename):
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables['concentration']
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    try:
        flag = f.variables['qc_concentration'][:]
        data[np.logical_and(flag>0, flag<=20)]=np.nan
    except:
        data[data<0]=np.nan
    f.close()
    return(time,data,timeunit,dataunit)

#%%  CVI
# filename='../data/inletcvi/enaaafinletcviF1.c1.20170718.083145.nc'
# varname='size_distribution'

def read_cvi_aceena(filename):    

    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    alt = f.variables['alt'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables['cvi_mode']
    cvimode = d_id[:]
    d_id = f.variables['inlet_selector']
    cviinlet = d_id[:]
    d_id = f.variables['enhancement_factor']
    enhance_factor = d_id[:]
    d_id = f.variables['inlet_dilution_factor']
    dilution_factor = d_id[:]
    f.close()
    return(time,lon,lat,alt,timeunit,cvimode,cviinlet,enhance_factor,dilution_factor)

#%% met
# filename='../data/arm-met/iop1/sgpmetE13.b1.20160502.000000.cdf'
# varname='wdir_vec_mean'
def read_met(filename,varname):
    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    f.close()
    return(time,data,timeunit,dataunit,long_name)

#%% mwrret 
def read_mwr(filename,varname):
    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    try:
        flag = f.variables['qc_'+varname][:]
    except:
        flag = data*0
    f.close()
    return(time,data,timeunit,dataunit,flag)

#%% pblh from sonde
# filename='../data/arm-pblh/sgppblhtsonde1mcfarlC1.c1.20160501.052900.cdf'
def read_pblhtsonde1(filename):
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    height = f.variables['height_ss'][:]
    T = f.variables['air_temp'][:]
    p = f.variables['atm_pres'][:]
    p2 = f.variables['pressure_gridded'][:]
    rh = f.variables['rh'][:]
    wspd = f.variables['wspd'][:]
    pbl_type = ['heffter','liu_liang','bulk_richardson_pt25','bulk_richardson_pt5']
    pblh = np.nan*np.zeros(len(pbl_type))
    long_name=list()
    for i in range(len(pbl_type)):
        d_id = f.variables['pbl_height_'+pbl_type[i]]
        qc = f.variables['qc_pbl_height_'+pbl_type[i]][:]
        if qc==0:
            pblh[i] = d_id[:]
        long_name.append(d_id.long_name)
    f.close()
    
    height2=np.interp(p,np.flip(p2),np.flip(height))
    return(time,timeunit,height2,p,T,rh,wspd,pblh,pbl_type,long_name)

#%% pblh from mpl
# filename='../data/arm-pblh/sgppblhtmpl1sawyerliC1.c1.20160830.000002.nc'
def read_pblhtmpl1(filename):
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename,'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    height = f.variables['height'][:]
    d_id = f.variables['annealing_pbl_height_sawyer_li']
    qc = f.variables['qc_annealing_pbl_height_sawyer_li'][:]
    pblh = d_id[:]
    pblh[qc!=0] = np.nan
    f.close()
    return(time,timeunit,pblh)

#%%  BNL SMPS data
# filename='../data/bnl-smps/sgpaossmpsS01.a1.20160513.000000.nc'
# varname='total_concentration'

def read_smps_bnl(filename,varname):    

    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    diameter = f.variables['diameter_midpoint'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    f.close()
    return(time,diameter,data,timeunit,dataunit,long_name)



#%%  UHSAS
# filename='../data/arm-uhsas/sgpaosuhsasS01.a1.20160501.000003.nc'
# varname='size_distribution'

def read_uhsas(filename):    

    from netCDF4 import Dataset
    f = Dataset(filename,'r')
    # read in variables
    dmin = f.variables['lower_size_limit'][:]
    dmax = f.variables['upper_size_limit'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    try:        # direct read concentration
        d_id = f.variables['concentration']
        data = d_id[:]
        dataunit = d_id.units
        long_name = d_id.long_name
    except:     # calculate from raw count
        import numpy as np
        raw_count = f.variables['size_distribution'][:]
        flow_rate = f.variables['sampling_volume'][:]/60.    # cc/min to cc/s
        sample_time = 10    # sample interval is 10s
        data=np.full(raw_count.shape,np.nan)
        for bb in range(data.shape[1]):
            data[:,bb] = raw_count[:,bb] /flow_rate /sample_time
        dataunit='1/cc'
        long_name='size distribution'
    f.close()
    return(time,dmin,dmax,data,timeunit,dataunit,long_name)
