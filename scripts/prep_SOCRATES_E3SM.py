# -*- coding: utf-8 -*-

"""
prepare E3SM data from SOCRATES
inlcude:
    prep_SOCRATES_E3SM
"""

import os
import sys
import yaml
import esmac_diags
import esmac_diags.preprocessing.prep_SOCRATES_E3SM as prep

import warnings
warnings.filterwarnings("ignore")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

input_path = config['model_input_path']
input_filehead = config['model_input_filehead']
# input_path = '../raw_data/model/'
# input_filehead = 'E3SMv2_SO'
output_path = '../prep_data/SOCRATES/model/'
output_filehead = 'E3SMv2_SOCRATES'

# flight data path
obs_input_path = config['obs_input_path']
RFpath = obs_input_path + 'SOCRATES/aircraft/aircraft_lowrate/'

# time frequencies
aircraft_dt = config['model_aircraft_dt']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
prep.prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, RFpath, dt=aircraft_dt)
