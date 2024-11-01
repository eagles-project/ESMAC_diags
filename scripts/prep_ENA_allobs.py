"""
prepare all observational data from ENA
inlcude:
    prep_ENA_sfc
    prep_ENA_satellite
"""

import os
import sys
import yaml
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_ENA_sfc as sfc
import esmac_diags.preprocessing.prep_ENA_satellite as sat


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

# data path for obs
obs_input_path = config['obs_input_path']

# surface data path
acsmpath = obs_input_path+'ENA/enaaosacsm/'
armbepath = obs_input_path+'ENA/enaarmbe/'
arsclpath = obs_input_path+'ENA/enaarsclkazrbnl1kolliasC1.c0/'
ccnsfcpath = obs_input_path+'ENA/enaaosccn1colspectraC1.b1/'
cpcpath = obs_input_path+'ENA/enaaoscpcf/'
aerosolmaskpath = obs_input_path+'ENA/ENA_AerosolMask/'
uhsaspath = obs_input_path+'ENA/enaaosuhsasC1.b1/'
mfrsrpath = obs_input_path+'ENA/enamfrsrcldod1minC1.c1/'
ndroppath = obs_input_path+'ENA/enandrop/'
WUpath = obs_input_path+'ENA/Wu_etal_retrieval/'
# satellite data path
visstgridpath = obs_input_path+'ENA/visst/grid/'
visstpixpath = obs_input_path+'ENA/visst/pix_3x3/'

# output data path
prep_data_path = '../prep_data/ENA/'

# time frequencies
surface_dt = config['obs_surface_dt']
satellite_dt = config['obs_satellite_dt']

# other settings
height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare surface data. output time in 1hr (dt=3600s) resolution
for year in range(2018,2018):
    # print('prepare surface data:')
    # sfc.prep_ACSM(acsmpath, prep_data_path+'surface/',year, dt=surface_dt)            # aerosol composition
    # sfc.prep_ccn(ccnsfcpath, prep_data_path+'surface/',year, dt=surface_dt)              # CCN number concentration
    # sfc.prep_CPC(cpcpath, prep_data_path+'surface/',year, dt=surface_dt)    # aerosol number concentration (>10nm)
    # sfc.prep_CPC_withENAmask(aerosolmaskpath, prep_data_path+'surface/',year, dt=surface_dt)# aerosol number concentration with mask
    # sfc.prep_CNsize_UHSAS(uhsaspath, prep_data_path+'surface/',year, dt=surface_dt)   # aerosol size distribution from UHSAS
    # sfc.prep_cloud_2d(armbepath, prep_data_path+'surface/', height_out,year, dt=surface_dt)   # 2D cloud fraction
    # sfc.prep_cloudheight_ARSCL(arsclpath, prep_data_path+'surface/',year, dt=surface_dt)   # cloud height 
    # sfc.prep_mfrsr_cod(mfrsrpath,  prep_data_path+'surface/',year, dt=surface_dt)     # cloud optical depth from MFRSR
    # sfc.prep_mfrsr_Reff(mfrsrpath,  prep_data_path+'surface/',year, dt=surface_dt)    # cloud effective radius from MFRSR
    # sfc.prep_LWP(armbepath, mfrsrpath, prep_data_path+'surface/',year, dt=surface_dt) # cloud liquid water path
    # sfc.prep_LTS(armbepath, arsclpath, prep_data_path+'surface/',year, dt=surface_dt)            # lower tropospheric stability
    # sfc.prep_precip(armbepath, prep_data_path+'surface/',year, dt=surface_dt)         # surface precipitation
    # sfc.prep_radiation(armbepath, prep_data_path+'surface/',year, dt=surface_dt)      # surface radiation
    # sfc.prep_totcld(armbepath, prep_data_path+'surface/',year, dt=surface_dt)         # cloud fraction. from ARSCL, TSI and satellite sources
    # sfc.prep_Ndrop(ndroppath, prep_data_path+'surface/',year, dt=surface_dt)          # cloud droplet number retrieval from ARM Ndrop VAP
    # sfc.prep_Nd_WU(WUpath, prep_data_path+'surface/',year, dt=surface_dt)          # cloud droplet number retrieval from Wu's retrieval
    # # prepare satellite data. output time in 1hr (dt=3600s) resolution
    print('prepare satellite data:')
    sat.prep_VISST_grid(visstgridpath, prep_data_path+'satellite/',year, dt=satellite_dt)     # VISST 0.5x0.5 degree gridded data
    sat.prep_VISST_pixel(visstpixpath, prep_data_path+'satellite/',year, dt=satellite_dt)     # VISST 4km pixel-level data
