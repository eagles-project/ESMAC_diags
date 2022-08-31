"""
prepare all observational data for ENA
inlcude:
    prep_ENA_E3SM in src/esmac_diags/preprocessing/
"""

import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_ENA_E3SM as prep

import warnings
warnings.filterwarnings("ignore")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settings
input_path = '../raw_data/model/'
output_path = '../prep_data/ENA/model/'
input_filehead = 'E3SMv2_SGP_ENA_2011_2020'
output_filehead = 'E3SMv2_ENA'

# vertical coordinates for output
lev_out=np.arange(25.,1001,25.)
height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                    1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                    3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                    10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output time in 1hr (dt=3600s) for other data
prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, dt=3600)
prep.prep_E3SM_profiles(input_path, input_filehead, output_path, output_filehead, height_out, lev_out=lev_out, dt=3600)

input_filehead = 'E3SMv1_SGP_ENA_2011_2020'
output_filehead = 'E3SMv1_ENA'
prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, dt=3600)
prep.prep_E3SM_profiles(input_path, input_filehead, output_path, output_filehead, height_out, lev_out=lev_out, dt=3600)
