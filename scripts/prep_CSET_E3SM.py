# -*- coding: utf-8 -*-

"""
prepare E3SM data from CSET
inlcude:
    prep_CSET_E3SM
"""

import esmac_diags
import esmac_diags.preprocessing.prep_CSET_E3SM as prep

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
# input_filehead = 'E3SMv2_CSET'
output_path = '../prep_data/CSET/model/'
output_filehead = 'E3SMv2_CSET'

# flight data path
obs_input_path = config['obs_input_path']
RFpath = obs_input_path + 'CSET/aircraft/aircraft_lowrate/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
prep.prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, RFpath, dt=60)
