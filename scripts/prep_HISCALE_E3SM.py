"""
prepare all observational data from HISCALE
inlcude:
    prep_HISCALE_E3SM in src/esmac_diags/preprocessing/
"""

import os
import sys
import yaml
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_HISCALE_E3SM as prep

import warnings
warnings.filterwarnings("ignore")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settings
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

input_path = config['model_input_path']
input2d_filehead = config['model_2d_input_filehead']
input3d_filehead = config['model_3d_input_filehead']
# input_path = '../raw_data/model/'
# input_filehead = 'E3SMv1_SGP_ENA_2011_2020'

output_path = '/pscratch/sd/a/avarble/eagles/ESMAC_DIAG/prep_data/HISCALE/model/'
output_filehead = 'HISCALE'

# iwg data path for aircraft information
obs_input_path = config['obs_input_path']
iwgpath = obs_input_path + 'HISCALE/aircraft/mei-iwg1/'

# time frequencies
aircraft_dt = config['model_aircraft_dt']
surface_dt = config['model_surface_dt']
profile_dt = config['model_profile_dt']

# vertical coordinates for output
lev_out=np.arange(25.,1001,25.)
height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                    1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                    3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                    10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output time in 1min (dt=60s) resolution for flight track and 1hr (dt=3600s) for other data
# prep.prep_E3SM_flight(input_path, input2d_filehead, input3d_filehead, output_path, output_filehead, iwgpath, dt=aircraft_dt, config=config)
# prep.prep_E3SM_sfc(input_path, input2d_filehead, input3d_filehead, output_path, output_filehead, dt=surface_dt, config=config)
prep.prep_E3SM_profiles(input_path, input2d_filehead, input3d_filehead, output_path, output_filehead, height_out, lev_out=lev_out, dt=profile_dt, config=config)
