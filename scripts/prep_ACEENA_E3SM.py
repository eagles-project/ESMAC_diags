"""
prepare all observational data from ACEENA
inlcude:
    prep_ACEENA_E3SM in src/esmac_diags/preprocessing/
"""

import os
import sys
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_ACEENA_E3SM as prep

import warnings
warnings.filterwarnings("ignore")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settings
# Load configuration file
config_file = sys.argv[1]
config = load_config(config_file)
input_path = config['model_input_path']
input_filehead = config['model_input_filehead']
# input_path = '../raw_data/rrm/ena_rrm/'
# input_filehead = 'ena_ne32x32pg2'
output_path = '../prep_data/ACEENA/model/'
output_filehead = 'E3SMv2_ACEENA'

# iwg data path for aircraft information
obs_input_path = config['obs_input_path']
iwgpath = obs_input_path + 'aircraft/IWG/'

# time frequencies
aircraft_dt = config['model_aircraft_dt']
surface_dt = config['model_surface_dt']
profile_dt = config['model_profile_dt']

# vertical coordinates for output
lev_out=np.arange(25.,1001,25.)    #pressure
height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                    1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                    3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                    10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output time in 1min (dt=60s) resolution for flight track and 1hr (dt=3600s) for other data
# prep.prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, iwgpath, dt=aircraft_dt, config)
prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, dt=surface_dt, config)
# prep.prep_E3SM_profiles(input_path, input_filehead, output_path, output_filehead, height_out, lev_out=lev_out, dt=profile_dt, config)
