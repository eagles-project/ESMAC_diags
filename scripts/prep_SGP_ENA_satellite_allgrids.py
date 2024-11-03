
import os
import sys
import yaml
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_ENA_satellite as enasat
import esmac_diags.preprocessing.prep_SGP_satellite as sgpsat


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# Load configuration file
config_file = sys.argv[1]
stream = open(config_file, "r")
config = yaml.full_load(stream)

# data path for obs
obs_input_path = config['obs_input_path']

# time frequencies
satellite_dt = config['obs_satellite_dt']

# ENA
# satellite data path
visstgridpath = obs_input_path+'ENA/visst/grid/'

# output data path
prep_data_path = '../prep_data/ENA/'

#for year in range(2016,2019):
#   enasat.prep_VISST_grid_allgrids(visstgridpath, prep_data_path+'satellite/',year, dt=satellite_dt)     # VISST 0.5x0.5 degree gridded data

# SGP
# satellite data path
visstgridpath = obs_input_path+'SGP/visst/grid/'

# output data path
prep_data_path = '../prep_data/SGP/'

for year in range(2020,2021):
    sgpsat.prep_VISST_grid_allgrids(visstgridpath, prep_data_path+'satellite/',year, dt=satellite_dt)     # VISST 0.5x0.5 degree gridded data
