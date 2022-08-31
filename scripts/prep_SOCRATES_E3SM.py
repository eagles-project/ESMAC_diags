# -*- coding: utf-8 -*-

"""
prepare E3SM data from SOCRATES
inlcude:
    prep_SOCRATES_E3SM
"""

import esmac_diags
import esmac_diags.preprocessing.prep_SOCRATES_E3SM as prep

import warnings
warnings.filterwarnings("ignore")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
input_path = '../raw_data/model/'
output_path = '../prep_data/SOCRATES/model/'
input_filehead = 'E3SMv2_SO'
output_filehead = 'E3SMv2_SOCRATES'

# flight data path
RFpath = '../raw_data/obs/SOCRATES/aircraft/aircraft_lowrate/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
prep.prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, RFpath, dt=60)
