# -*- coding: utf-8 -*-

"""
prepare all observational data from CSET
inlcude:
    prep_CSET_flight
"""

import os
import sys
import yaml
import esmac_diags
import esmac_diags.preprocessing.prep_CSET_flight as air

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

# data path for obs
obs_input_path = config['obs_input_path']

# flight data path
RFpath = obs_input_path+'CSET/aircraft/aircraft_lowrate/'

# output data path
prep_data_path = '/global/cscratch1/sd/sqtang/EAGLES/ESMAC_Diags_v2/prep_data/CSET/'

# time frequencies
aircraft_dt = config['obs_aircraft_dt']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare flight data; output time in 1min (dt=60s) resolution
print('prepare aircraft data:')
air.prep_CNsize(RFpath, prep_data_path+'flight/', dt=aircraft_dt)    # aerosol size distribution
air.prep_CN(RFpath, prep_data_path+'flight/', dt=aircraft_dt)        # CN number concentration
air.prep_LWC(RFpath, prep_data_path+'flight/', dt=aircraft_dt)       # cloud liquid water content
air.prep_Nd(RFpath, prep_data_path+'flight/', dt=aircraft_dt)        # cloud size distribution
