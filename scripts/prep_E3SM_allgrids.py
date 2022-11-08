
import esmac_diags
import esmac_diags.preprocessing.prep_E3SM_allgrids_SGP_ENA as prep

import warnings
warnings.filterwarnings("ignore")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% settings
site = 'ENA'
input_path = '../raw_data/model/'
output_path = '../prep_data/'+site+'/model/'

# output time in 1hr (dt=3600s) for other data
input_filehead = 'E3SMv1_SGP_ENA_2011_2020'
output_filehead = 'E3SMv1_'+site
prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)

# input_filehead = 'E3SMv2_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv2_'+site
# prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)

site = 'SGP'
input_path = '../raw_data/model/'
output_path = '../prep_data/'+site+'/model/'

# output time in 1hr (dt=3600s) for other data
input_filehead = 'E3SMv1_SGP_ENA_2011_2020'
output_filehead = 'E3SMv1_'+site
prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)

# input_filehead = 'E3SMv2_SGP_ENA_2011_2020'
# output_filehead = 'E3SMv2_'+site
# prep.prep_E3SM_sfc(input_path, input_filehead, output_path, output_filehead, site, dt=3600)

