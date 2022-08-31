
import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_ENA_satellite as enasat
import esmac_diags.preprocessing.prep_SGP_satellite as sgpsat


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# ENA

# satellite data path
visstgridpath = '../raw_data/obs/ENA/visst/grid/'
# output data path
prep_data_path = '../prep_data/ENA/'

#for year in range(2016,2019):
#   enasat.prep_VISST_grid_allgrids(visstgridpath, prep_data_path+'satellite/',year, dt=3600)     # VISST 0.5x0.5 degree gridded data


# SGP
# satellite data path
visstgridpath = '../raw_data/obs/SGP/visst/grid/'
# output data path
prep_data_path = '../prep_data/SGP/'

for year in range(2020,2021):
    sgpsat.prep_VISST_grid_allgrids(visstgridpath, prep_data_path+'satellite/',year, dt=3600)     # VISST 0.5x0.5 degree gridded data
