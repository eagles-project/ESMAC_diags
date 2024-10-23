"""
prepare all observational data from ACEENA
inlcude:
    prep_ACEENA_flight
    prep_ACEENA_sfc
    prep_ACEENA_satellite
"""

import os
import sys
import yaml
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_ACEENA_flight as air
import esmac_diags.preprocessing.prep_ACEENA_sfc as sfc
import esmac_diags.preprocessing.prep_ACEENA_satellite as sat

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# flight data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

# data path for obs
obs_input_path = config['obs_input_path']

# aircraft data path
iwgpath = obs_input_path+'ACEENA/aircraft/IWG/'
amspath = obs_input_path+'ACEENA/aircraft/shilling-hrfams/'
beasdpath = obs_input_path+'ACEENA/aircraft/beasd/'
ccnairpath = obs_input_path+'ACEENA/aircraft/ccn_aaf/'
wcmpath = obs_input_path+'ACEENA/aircraft/wcm_ACEENA/'
fimspath = obs_input_path+'ACEENA/aircraft/FIMS/'
pcasppath = obs_input_path+'ACEENA/aircraft/pcasp_g1/'
cvipath = obs_input_path+'ACEENA/aircraft/inletcvi/'
opcpath = obs_input_path+'ACEENA/aircraft/opciso/'
cpcairpath = obs_input_path+'ACEENA/aircraft/cpc_aaf/'
mergeSDpath = obs_input_path+'ACEENA/aircraft/mergedSD/'
merged_CNsize_path = obs_input_path+'ACEENA/aircraft/merged_bin/'

# surface data path
acsmpath = obs_input_path+'ACEENA/surface/arm_acsm/'
armbepath = obs_input_path+'ACEENA/profile/armbe/'
arsclpath = obs_input_path+'ACEENA/profile/arscl/'
arsclbndpath = obs_input_path+'ACEENA/profile/arsclbnd/'
ccnsfcpath = obs_input_path+'ACEENA/surface/arm_aosccn1/'
cpcpath = obs_input_path+'ACEENA/surface/arm_cpcf/'
aerosolmaskpath = obs_input_path+'ENA/ENA_AerosolMask/'
uhsaspath = obs_input_path+'ACEENA/surface/arm_uhsas/'
mfrsrpath = obs_input_path+'ACEENA/surface/arm_mfrsr/'
metpath = obs_input_path+'ACEENA/obs/surface/arm_met/'
parspath =  obs_input_path+'ACEENA/obs/surface/arm_pars/'
radfluxpath =  obs_input_path+'ACEENA/obs/surface/arm_radflux/'
ndroppath = obs_input_path+'ACEENA/surface/enandrop/'
Wuetalpath = obs_input_path+'ACEENA/surface/Wu_etal_retrieval/'

# satellite data path
visstgridpath = obs_input_path+'ACEENA/satellite/visst/grid/'
visstpixpath = obs_input_path+'ACEENA/satellite/visst/pix_3x3/'

# output data path
prep_data_path = '../prep_data/ACEENA/'

# time frequencies
aircraft_dt = config['obs_aircraft_dt']
surface_dt = config['obs_surface_dt']
satellite_dt = config['obs_satellite_dt']

# other settings
height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare flight data; output time in aircraft_dt resolution
# print('prepare aircraft data:')
# air.prep_AMS(amspath, iwgpath, prep_data_path+'flight/', dt=aircraft_dt)       # aerosol composition
# air.prep_beasd(beasdpath,iwgpath, prep_data_path+'flight/', dt=aircraft_dt)               # merged aerosol size distribution
# air.prep_CCN(ccnairpath, iwgpath, prep_data_path+'flight/', dt=aircraft_dt)       # CCN number concentration
# air.prep_CPC(cpcairpath, iwgpath, prep_data_path+'flight/', dt=aircraft_dt)       # aerosol number concentration (>3 or 10nm)
# air.prep_mergeSD(mergeSDpath, iwgpath, prep_data_path+'flight/', dt=aircraft_dt) # merged cloud size distribution
# air.prep_mergesize_withCPC_ACEENA(cpcairpath, fimspath, pcasppath, opcpath, 
#                             iwgpath, cvipath, prep_data_path+'flight/', dt=aircraft_dt) # merged aerosol size distribution
# air.prep_PCASP100(pcasppath, iwgpath, prep_data_path+'flight/', dt=aircraft_dt)# aerosol number concentration (>100nm)
# air.prep_WCM(wcmpath, iwgpath, prep_data_path+'flight/', dt=aircraft_dt)       # cloud liquid water content

# prepare surface data. output time in surface_dt resolution
print('prepare surface data:')
sfc.prep_ACSM(acsmpath, prep_data_path+'surface/', dt=surface_dt)            # aerosol composition
# sfc.prep_ccn(ccnsfcpath, prep_data_path+'surface/', dt=surface_dt)              # CCN number concentration
# sfc.prep_CPC(cpcpath, prep_data_path+'surface/', dt=surface_dt)    # aerosol number concentration (>3 or 10nm)
# sfc.prep_CPC_withENAmask(aerosolmaskpath, prep_data_path+'surface/', dt=surface_dt)    # aerosol mask data
# sfc.prep_CNsize_UHSAS(uhsaspath, prep_data_path+'surface/', dt=surface_dt)   # aerosol size distribution from UHSAS
# sfc.prep_cloud_2d(armbepath, prep_data_path+'surface/', height_out, dt=surface_dt)   # 2D cloud fraction
# sfc.prep_cloudheight_ARSCL(arsclpath, prep_data_path+'surface/', dt=surface_dt)   # cloud height 
# sfc.prep_totcld(armbepath, prep_data_path+'surface/', dt=surface_dt)         # cloud fraction. from ARSCL, TSI and satellite sources
# sfc.prep_LWP(armbepath, mfrsrpath, prep_data_path+'surface/', dt=surface_dt) # cloud liquid water path
# sfc.prep_Ndrop(ndroppath, prep_data_path+'surface/', dt=surface_dt)          # cloud droplet number retrieval from ARM Ndrop VAP
# sfc.prep_Nd_WU(Wuetalpath, prep_data_path+'surface/', dt=surface_dt)           # cloud droplet number retrieval from Wu et al. algorithm
# sfc.prep_mfrsr_cod(mfrsrpath,  prep_data_path+'surface/', dt=surface_dt)     # cloud optical depth from MFRSR
# sfc.prep_mfrsr_Reff(mfrsrpath,  prep_data_path+'surface/', dt=surface_dt)    # cloud effective radius from MFRSR
# sfc.prep_radiation(armbepath, prep_data_path+'surface/', dt=surface_dt)      # surface radiation
# sfc.prep_LTS(armbepath, prep_data_path+'surface/', dt=surface_dt)            # lower tropospheric stability
# sfc.prep_precip(armbepath, prep_data_path+'surface/', dt=surface_dt)         # surface precipitation

# # prepare satellite data. output time in satellite_dt resolution
# print('prepare satellite data:')
# sat.prep_VISST_grid(visstgridpath, prep_data_path+'satellite/', dt=satellite_dt)     # VISST 0.5x0.5 degree gridded data
# sat.prep_VISST_pixel(visstpixpath, prep_data_path+'satellite/', dt=satellite_dt)     # VISST 4km pixel-level data
