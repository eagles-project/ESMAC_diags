
#%% READ AMS composition data
# filename='../../data/HiScale/obs/aircraft/shilling-ams\\HiScaleAMS_G1_20160425_R0.ict'
def read_ams(filename):    
    import numpy as np
        
    f=open(filename,'r')
    
    # read in data:
    
    h='aaa'
    for ii in range(40):
        h=f.readline()
    while(h[0:7]!='dat_ams'):
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2,varlist)


#%% read in CPC data (mei-cpc)

# filename='../data/mei-cpc/CPC_G1_20160830143515_R2_HiScale001s.ict.txt'
def read_cpc(filename):    
    import numpy as np
        
    f=open(filename,'r')
    
    # read in data:
    
    h='aaa'
    while h[0:19]!='Start_UTC,CPC_Conc_':
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    
    f.close()
    data2[data2<-9990]=np.nan
    return(data2,varlist)

#%% read in CCN data (mei-ccn)

# filename='../data/mei-ccn/CCN_G1_20160518170015_R2_HiScale001s.ict'
def read_ccn_hiscale(filename):    
    import numpy as np
        
    f=open(filename,'r')
    
    # read in data:
    
    h='aaa'
    while h[0:14]!='Start_UTC,DT_A':
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2,varlist)


#%% read CVI
# filename='../data/pekour-cvi/CVI_G1_20160518170015_R4_HISCALE_001s.ict.txt'
def read_cvi_hiscale(filename):    
    import numpy as np
    f=open(filename,'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2,varlist)

#%% read in FIMS data (wang-fims)
def read_fims(filename):
    
    import numpy as np
    
    f=open(filename,'r')
    
    # read in data:
    
    h='aaa'
    while h[0:10]!='UTC,n_Dp_1':
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    
    f.close()
    
    # data2[data2<-9990]=np.nan
    
    # multiply dlnDp since FIMS data are dN/dlnDp 
    # data2[1:31,:]=data2[1:31,:]*0.1272
    return(data2,varlist)


#%% read FIMS bins information (wang-fims)
# filename='../data/wang-fims/HISCALE_FIMS_bins_R1.dat'
def read_fims_bin(filename):
    import numpy as np
    f=open(filename,'r')
    h=f.readline()
    h=h.strip()
    dmean=[]
    dmin=[]
    dmax=[]
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        dmin.append(source[0])
        dmean.append(source[1])
        dmax.append(source[2])
    f.close()

    return(dmean,dmin,dmax)

#%% read in IWG1 data (mei-iwg1)

# filename='../data/mei-iwg1/aaf.iwg1001s.g1.hiscale.20160511a.a2.txt'
# try:
def read_iwg1(filename):
    
    import numpy as np
    
    f=open(filename,'r')
    varname=f.readline()
    varunit=f.readline()
    varname=varname.strip()
    varname=varname.split(',')
    varunit=varunit.strip()
    varunit=varunit.split(',')
    a = np.column_stack((np.asarray(varname),np.asarray(varunit)))
    data=[]
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(columns[i])
        data.append(source)
    
    f.close()
    # data[data<-9990]=np.nan
    return(data,a)
    
# except:
#     print('error')

#%% read kappa
# filename='../data/Kappa/FIMS_kappa_IOP2_part2/20160830aAir_data_table_FIMS_kappa_col_A.dat'
def read_kappa(filename):    
    import numpy as np
    f=open(filename,'r')
    varname=f.readline()
    # varunit=f.readline()
    varname=varname.strip()
    varname=varname.split()
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split()
        if columns==[]:
            continue
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    
    f.close()
    return(data2,varname)

#%% read OPC
# filename='../data/opciso/OPCISO_G1_20170707103233_R3_ACEENA_001s.ict'
def read_opc(filename):    
    import numpy as np
    f=open(filename,'r')
    h='aaa'
    while h[0:26]!='UPPER_BIN_SIZE_micrometer:':
        h=f.readline()
    h=h.strip()
    d_max = h[26:].split(',')
    d_max=[float(x) for x in d_max]
    h=f.readline()
    h=h.strip()
    d_min = h[26:].split(',')
    d_min=[float(x) for x in d_min]
    h=f.readline()
    h=h.strip()
    d_center = h[24:].split(',')
    d_center=[float(x) for x in d_center]
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2,np.array(d_min),np.array(d_max),np.array(d_center),varlist)

#%% check with processed data
def read_processed(filename):
    # filename='../data/merged-bin/avesize_0830a_leg03.dat'
    import numpy as np
    f=open(filename)
    h=f.readline()
    ii=0
    for line in f:
        line=line.strip()
        columns = line.split()
        source=[]
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        if ii==0:
            fims=np.asarray(source)
        elif ii<30:
            fims=np.column_stack((fims,source))
        elif ii==30:
            pcasp=np.asarray(source)
        elif ii>30 and ii<60:
            pcasp=np.column_stack((pcasp,source))
        elif ii==60:
            merge=np.asarray(source)
        elif ii>60:
            merge=np.column_stack((merge,source))
        ii=ii+1
    f.close()
    return (merge,fims,pcasp)



#%% read PCASP
# filename='../data/tomlinson-pcasp/pcasp_g1_20160511163038_R2_L1_hiscale001s.ict.txt'
def read_pcasp(filename):    
    import numpy as np
    f=open(filename,'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    f.close()
    data2[data2<-9990]=np.nan
    return(data2,varlist)
    
#%% read research flight data from NCAR
def read_RF_NCAR(filename,varname):
    
    from netCDF4 import Dataset
    
    f = Dataset(filename,'r')
    
    # read in variables
    t_id = f.variables['Time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    try:
        cellsize = d_id.CellSizes
    except:
        cellsize = 'N/A'
    try:
        cellunit = d_id.CellSizeUnits
    except:
        cellunit = 'N/A'
    
    f.close()
    
    return(time,data,timeunit,dataunit,long_name,cellsize,cellunit)

#%% SO2
# filename='../data/springston-tracegases/SO2/SO2_G1_20160510_R2.ict'
def read_so2(filename):
    import numpy as np
    f=open(filename,'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    f.close()
    data2[data2<-9990]=np.nan
    return(data2,varlist)
    
#%% read UHSAS
# filename='../data/tomlinson-uhsas/uhsasa_g1_20160425155810_R1_L1_hiscale001s.ict.txt'
def read_uhsas(filename):    
    import numpy as np
    f=open(filename,'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    f.close()
    data2[data2<-9990]=np.nan
    return(data2,varlist)

#%% read LWC from WCM
# filename='../data/matthews-wcm/WCM_G1_20160909150045_R2_HISCALE001s.ict.txt'
def read_wcm(filename):    
    import numpy as np
    f=open(filename,'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if 'data2' in locals():
            data2=np.column_stack((data2,source))
        else:
            data2=np.asarray(source)
    f.close()
    data2[data2<-9990]=np.nan
    return(data2,varlist)


# In[test read in data]
# filename = '../data/wang-fims/FIMS_G1_20160830_R1_L1_HISCALE_001s.ict'

# (x,y)=read_fims(filename)

# filename='../data/mei-iwg1/aaf.iwg1001s.g1.hiscale.20160511a.a2.txt'
# (x,y)=read_iwg1(filename)