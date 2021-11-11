"""
script to generate all plots

Instruction:
    edit the first section "user-specified settings" for 
"""
import numpy as np
from esmac_diags.preprocessing import *


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# user-specified settings
settings = {}

# set field campaign name. More settings on specific field campaigns are in next section
# HISCALE, ACEENA, CSET, SOCRATES, MAGIC, MARCUS
settings['campaign'] = 'MAGIC'

# set model names. 
settings['Model_List'] = ['E3SMv1']
# settings['Model_List'] = ['E3SMv1','EAMv1_CONUS_RRM']

# set line colors for each model. corresponding to the Model_List
settings['color_model'] = ['r','b','g']

# set field campaign IOPs. Only used for HISCALE and ACEENA. 
# IOP1/IOP2 
settings['IOP'] = 'IOP1'

########## set filepath for preprocessing. If you don't run preprocessing, ignore this part
# path of E3SM model data (h3) for preprocessing. same length of Model_List
# settings['E3SM_hourly_path'] = ['/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_hourly/']
settings['E3SM_hourly_path'] = \
    ['/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_hourly/', \
    '/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_hourly/']
settings['E3SM_hourly_filehead'] = ['E3SMv1','EAMv1_CONUS_RRM']
#############

# Please specify the path of your data and for figure output
settings['datapath'] = '/global/homes/s/sqtang/EAGLES/ESMAC_diags/data/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def add_other_setting(settings):
    """
    add other settings for different field campaigns

    Parameters
    ----------
    settings : dictionary
        all setting variables

    Returns
    -------
    settings

    """
    datapath = settings['datapath']+settings['campaign']
    
    
    # other settings for different field campaigns
    if settings['campaign']=='HISCALE':
        settings['site']='SGP'
        # lat/lon at SGP
        settings['lat0'] = 36.6059
        settings['lon0'] = 360-97.48792     # 0-360
        # time periods for IOPs. needed in preprocessing of surface data
        if settings['IOP']=='IOP1':
            settings['start_date']='2016-04-25'
            settings['end_date']='2016-05-29'
        elif settings['IOP']=='IOP2':
            settings['start_date']='2016-08-27'
            settings['end_date']='2016-09-22'
        #### observational data path. ######
        settings['merged_size_path'] = datapath+'/obs/aircraft/merged_bin/'
        settings['iwgpath'] = datapath+'/obs/aircraft/mei-iwg1/'
        settings['fimspath'] = datapath+'/obs/aircraft/wang-fims/'
        settings['pcasppath'] = datapath+'/obs/aircraft/tomlinson-pcasp/'
        settings['cvipath'] = datapath+'/obs/aircraft/pekour-cvi/'
        #### pre-processed model data path ######
        settings['E3SM_sfc_path'] = datapath+'/model/surface/'
        settings['E3SM_aircraft_path'] = datapath+'/model/flighttrack/'
        settings['E3SM_profile_path'] = datapath+'/model/profile/'
        ########################
        
    elif settings['campaign']=='ACEENA':
        settings['site']='ENA'
        # lat/lon for ENA
        settings['lat0'] = 39.09527
        settings['lon0'] = 360-28.0339
        # time periods for IOPs. needed in preprocessing of surface data
        if settings['IOP']=='IOP1':
            settings['start_date']='2017-06-20'
            settings['end_date']='2017-07-20'
        elif settings['IOP']=='IOP2':
            settings['start_date']='2018-01-21'
            settings['end_date']='2018-02-19'
        #### observational data path. ######
        settings['merged_size_path']=datapath+'/obs/aircraft/merged_bin/'
        settings['iwgpath'] = datapath+'/obs/aircraft/IWG/'
        settings['cvipath'] = datapath+'/obs/aircraft/inletcvi/'
        settings['fimspath'] = datapath+'/obs/aircraft/FIMS/'
        settings['pcasppath'] = datapath+'/obs/aircraft/pcasp_g1/'
        settings['opcpath'] = datapath+'/obs/aircraft/opciso/'
        #### pre-processed model data path ######
        settings['E3SM_sfc_path'] = datapath+'/model/surface/'
        settings['E3SM_aircraft_path'] = datapath+'/model/flighttrack/'
        settings['E3SM_profile_path'] = datapath+'/model/profile/'
        #########
    elif settings['campaign']=='MAGIC':
        settings['site']='MAG'
        # reference lat/lon
        settings['lat0']=30.
        settings['lon0']=230.
        #### observational data path. ######
        settings['shipmetpath'] = datapath+'/obs/ship/raynolds-marmet/'
        #### pre-processed model data path ######    
        settings['E3SM_ship_path'] = datapath+'/model/shiptrack/'
        settings['E3SM_profile_path'] = datapath+'/model/profile/'
        ################
    elif settings['campaign']=='MARCUS':
        settings['site'] = 'MAR'
        # reference lat/lon
        settings['lat0'] = -40.
        settings['lon0'] = 120.
        #### observational data path. ######
        settings['shipmetpath'] = datapath+'/obs/ship/maraadmetX1.b1/'
        #### pre-processed model data path ######    
        settings['E3SM_ship_path'] = datapath+'/model/shiptrack/'
        settings['E3SM_profile_path'] = datapath+'/model/profile/'
        ############
    elif settings['campaign']=='CSET':
        # lat/lon at the airport
        settings['lat0'] = 38.5564
        settings['lon0'] = 360-121.3120
        #### observational data path. ######
        settings['RFpath'] = datapath+'/obs/aircraft/aircraft_lowrate/'
        #### pre-processed model data path ###### 
        settings['E3SM_aircraft_path'] = datapath+'/model/flighttrack/'
        settings['E3SM_profile_path'] = datapath+'/model/profile/'
        ###############
    elif settings['campaign']=='SOCRATES':
        # lat/lon at the airport
        settings['lat0'] = -42.8371
        settings['lon0'] = 147.5054
        #### observational data path. ######
        settings['RFpath'] = datapath+'/obs/aircraft/aircraft_lowrate/'
        #### pre-processed model data path ###### 
        settings['E3SM_aircraft_path'] = datapath+'/model/flighttrack/'
        settings['E3SM_profile_path'] = datapath+'/model/profile/'
        ################
    else:
        raise ValueError("does not recognize this campaign: "+settings['campaign'])
    return(settings)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# running command
all_settings = add_other_setting(settings)
prep_obs_mergesize_HISCALE.run_prep(all_settings)
prep_E3SM_flighttrack_bins.run_prep(all_settings)
    