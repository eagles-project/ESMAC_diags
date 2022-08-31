
"""
prepare all observational data from MARCUS
inlcude:
    prep_MARCUS_ship
"""

import esmac_diags.preprocessing.prep_MARCUS_ship as ship
import esmac_diags.preprocessing.prep_MARCUS_satellite as sat

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# ship data path
shipmetpath = '../raw_data/obs/MARCUS/ship/maraadmetX1.b1/'
mwrpath = '../raw_data/obs/MARCUS/ship/marmwrret1liljclouM1.s2/'
cpcpath = '../raw_data/obs/MARCUS/ship/maraoscpcf1mM1.b1/'
ccnpath = '../raw_data/obs/MARCUS/ship/maraosccn1colavgM1.b1/'
uhsaspath = '../raw_data/obs/MARCUS/ship/maraosuhsasM1.a1/'
exhaustfreepath = '../raw_data/obs/MARCUS/ship/ship_exhaustfree/'
visstgridpath = '../raw_data/obs/MARCUS/visst/grid/'

# output data path
prep_data_path = '../prep_data/MARCUS/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare ship data. output time in 1hr (dt=3600s) resolution
print('prepare ship data:')
# ship.prep_CCN(shipmetpath, ccnpath, prep_data_path+'ship/', dt=3600)              # CCN number concentration
# ship.prep_CN(shipmetpath, cpcpath, uhsaspath, prep_data_path+'ship/', dt=3600)    # aerosol number concentration
# ship.prep_CNsize(shipmetpath, uhsaspath, prep_data_path+'ship/', dt=3600)   # aerosol size distribution from UHSAS
# ship.prep_CCN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path+'ship/', dt=3600)              # CCN number concentration
# ship.prep_CN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path+'ship/', dt=3600)    # aerosol number concentration
# ship.prep_CNsize_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path+'ship/', dt=3600)   # aerosol size distribution from UHSAS
# ship.prep_MWR(shipmetpath, mwrpath, prep_data_path+'ship/', dt=3600) # cloud liquid water path

print('prepare satellite data:')
sat.prep_VISST_grid(shipmetpath, visstgridpath, prep_data_path+'satellite/', dt=3600)              # CCN number concentration