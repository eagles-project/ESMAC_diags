# -*- coding: utf-8 -*-

"""
prepare all observational data from MAGIC
inlcude:
    prep_MAGIC_ship
"""

import os
import sys
import yaml
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_MAGIC_ship as ship
import esmac_diags.preprocessing.prep_MAGIC_satellite as sat

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

# data path for obs
obs_input_path = config['obs_input_path']

# ship data path
shipmetpath = obs_input_path+'MAGIC/ship/magmarinemetM1.b1/'
mwrpath = obs_input_path+'MAGIC/ship/magmwrret1liljclouM1.s2/'
cpcpath = obs_input_path+'MAGIC/ship/magaoscpcfM1.a1/'
ccnpath = obs_input_path+'MAGIC/ship/magaosccn100M1.a1/'
uhsaspath = obs_input_path+'MAGIC/ship/magaosuhsasM1.a1/'
Ndpath = obs_input_path+'MAGIC/ship/Cloud_Micro_Retrieval/'
visstpixpath = obs_input_path+'MAGIC/visst/pix/'

# output data path
prep_data_path = '../prep_data/MAGIC/'

# time frequencies
surface_dt = config['obs_surface_dt']
satellite_dt = config['obs_satellite_dt']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare ship data. output time in 1hr (dt=3600s) resolution
# print('prepare ship data:')
# ship.prep_CCN(shipmetpath, ccnpath, prep_data_path+'ship/', dt=surface_dt)              # CCN number concentration
# ship.prep_CN(shipmetpath, cpcpath, uhsaspath, prep_data_path+'ship/', dt=surface_dt)    # aerosol number concentration (>3 or 10nm)
# ship.prep_CNsize(shipmetpath, uhsaspath, prep_data_path+'ship/', dt=surface_dt)   # aerosol size distribution from UHSAS
# ship.prep_MWR(shipmetpath, mwrpath, prep_data_path+'ship/', dt=surface_dt) # cloud liquid water path
# ship.prep_Nd_Wu_etal(Ndpath, prep_data_path+'ship/', dt=surface_dt)          # cloud droplet number retrieval from Wu et al.

print('prepare satellite data:')
sat.prep_VISST_pixel(shipmetpath, visstpixpath, prep_data_path+'satellite/', dt=satellite_dt)              # CCN number concentration
