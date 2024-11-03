
"""
prepare all observational data from MARCUS
inlcude:
    prep_MARCUS_ship
"""

import os
import sys
import yaml
import esmac_diags.preprocessing.prep_MARCUS_ship as ship
import esmac_diags.preprocessing.prep_MARCUS_satellite as sat

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

# data path for obs
obs_input_path = config['obs_input_path']

# ship data path
shipmetpath = obs_input_path+'MARCUS/ship/maraadmetX1.b1/'
mwrpath = obs_input_path+'MARCUS/ship/marmwrret1liljclouM1.s2/'
cpcpath = obs_input_path+'MARCUS/ship/maraoscpcf1mM1.b1/'
ccnpath = obs_input_path+'MARCUS/ship/maraosccn1colavgM1.b1/'
uhsaspath = obs_input_path+'MARCUS/ship/maraosuhsasM1.a1/'
exhaustfreepath = obs_input_path+'MARCUS/ship/ship_exhaustfree/'
visstgridpath = obs_input_path+'MARCUS/visst/grid/'

# output data path
prep_data_path = '../prep_data/MARCUS/'

# time frequencies
surface_dt = config['obs_surface_dt']
satellite_dt = config['obs_satellite_dt']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare ship data. output time in 1hr (dt=3600s) resolution
print('prepare ship data:')
# ship.prep_CCN(shipmetpath, ccnpath, prep_data_path+'ship/', dt=surface_dt)              # CCN number concentration
# ship.prep_CN(shipmetpath, cpcpath, uhsaspath, prep_data_path+'ship/', dt=surface_dt)    # aerosol number concentration
# ship.prep_CNsize(shipmetpath, uhsaspath, prep_data_path+'ship/', dt=surface_dt)   # aerosol size distribution from UHSAS
# ship.prep_CCN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path+'ship/', dt=surface_dt)              # CCN number concentration
# ship.prep_CN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path+'ship/', dt=surface_dt)    # aerosol number concentration
# ship.prep_CNsize_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path+'ship/', dt=surface_dt)   # aerosol size distribution from UHSAS
# ship.prep_MWR(shipmetpath, mwrpath, prep_data_path+'ship/', dt=surface_dt) # cloud liquid water path

print('prepare satellite data:')
sat.prep_VISST_grid(shipmetpath, visstgridpath, prep_data_path+'satellite/', dt=satellite_dt)              # CCN number concentration
