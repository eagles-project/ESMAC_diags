"""
prepare all observational data from HISCALE
inlcude:
    prep_HISCALE_E3SM in src/esmac_diags/preprocessing/
"""

import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_HISCALE_E3SM as prep

import warnings
warnings.filterwarnings("ignore")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settings
for iop in ['1', '2']:
    for grid in ['3km', 'ne30']:
        input_path = '../raw_data/model/cus'+iop+'_'+grid+'/'
        output_path = '../prep_data/HISCALE/model/'
        if grid=='3km':
            input_filehead = 'CUS'+iop+'.F2010-CICE.cus_ne32x32pg2'
        elif grid =='ne30':
            input_filehead = 'CUS'+iop+'.F2010-CICE.ne30pg2_EC30to60E2r2'
        output_filehead = 'hiscale_RRM'+grid+'_IOP'+iop

        # iwg data path for aircraft information
        iwgpath = '../raw_data/obs/HISCALE/aircraft/mei-iwg1/'

        # vertical coordinates for output
        lev_out=np.arange(25.,1001,25.)
        height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                    1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                    3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                    10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # output time in 1min (dt=60s) resolution for flight track and 1hr (dt=3600s) for other data
        prep.prep_E3SM_flight(input_path, input_filehead, output_path, output_filehead, iwgpath, dt=60)
        prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, dt=3600)
        prep.prep_E3SM_profiles(input_path, input_filehead, output_path, output_filehead, height_out, lev_out=lev_out, dt=3600)

