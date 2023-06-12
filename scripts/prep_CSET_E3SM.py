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
input_path = '../raw_data/model/'
output_path = '../prep_data/CSET/model/'
input_filehead = 'E3SMv2_CSET'
output_filehead = 'E3SMv2_CSET'

# flight data path
RFpath = '../raw_data/obs/CSET/aircraft/aircraft_lowrate/'

E3SMdomain_range = '202e_to_240e_19n_to_40n'    # domain range in E3SM regional output

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
prep.prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, RFpath, E3SMdomain_range='_'+E3SMdomain_range, dt=60)
