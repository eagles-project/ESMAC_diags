# settings of the aerosol diagnostic package

import numpy as np

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# global settings

############ these settings will be replaced by the settings in scripts_*.csh #############
# set field campaign name. More settings on specific field campaigns are in next section
campaign = 'HISCALE'
# set model names. up to three
Model_List = ['EAMv1_CONUS_RRM']
# set line colors for each model. corresponding to the Model_List
color_model = ['b','g']
# set IOP that the statistics, pdf and percentiles are averaged for. Only available for HISCALE and ACEENA
# IOP1/IOP2 
IOP = 'IOP2'
############ these settings will be replaced by the settings in scripts_*.csh #############


# path of the diagnostic package
package_path = '../../'

# path of E3SM model data (h3) for preprocessing. list with the same length of Model_List
E3SM_h3_path=[]
E3SM_h3_filehead=[]     # filename before .cam.h3.yyyy-mm-dd.00000.nc
for mm in Model_List:
    E3SM_h3_path.append('/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_h3/')
    if campaign=='MAGIC':
        E3SM_h3_filehead.append(mm+'_2012-2013')
    else:
#        E3SM_h3_filehead.append(mm+'_2014-2018')
        E3SM_h3_filehead.append(mm)
    #E3SM_h3_path.append('/qfs/projects/eagles/zhan524/simulations/compy_F20TRC5-CMIP6_ne30_EG1_R2_'+mm+'/h3/')
    #E3SM_h3_filehead.append('compy_F20TRC5-CMIP6_ne30_EG1_R2_'+mm)

# path of output figures
figpath_aircraft_timeseries = package_path+'figures/'+campaign+'/aircraft/timeseries/'
figpath_aircraft_statistics = package_path+'figures/'+campaign+'/aircraft/statistics/'
figpath_ship_timeseries = package_path+'figures/'+campaign+'/ship/timeseries/'
figpath_ship_statistics = package_path+'figures/'+campaign+'/ship/statistics/'
figpath_sfc_timeseries = package_path+'figures/'+campaign+'/surface/timeseries/'
figpath_sfc_statistics = package_path+'figures/'+campaign+'/surface/statistics/'
figpath_profile_timeseries = package_path+'figures/'+campaign+'/profile/timeseries/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings for different field campaigns

# set location and time information
if campaign=='HISCALE':
    site='SGP'
    # lat/lon at SGP
    lat0 = 36.6059
    lon0 = 360-97.48792     # 0-360
    # bin of flight heights to calculate percentiles
    height_bin = np.arange(300,4300,300)
    # height_bin = np.arange(400,4300,200)
    
    # time periods for IOPs. needed in preprocessing of surface data
    if IOP=='IOP1':
        start_date='2016-04-25'
        end_date='2016-05-29'
    elif IOP=='IOP2':
        start_date='2016-08-27'
        end_date='2016-09-22'
        
    # observational data path. 
    # aircraf measurements merged_bin data are used for all plot_flight_*.py to provide flight/cloud/CVI info
    merged_size_path=package_path+'data/'+campaign+'/obs/aircraft/merged_bin/'
    iwgpath = package_path+'data/'+campaign+'/obs/aircraft/mei-iwg1/'
    fimspath = package_path+'data/'+campaign+'/obs/aircraft/wang-fims/'
    pcasppath = package_path+'data/'+campaign+'/obs/aircraft/tomlinson-pcasp/'
    cvipath = package_path+'data/'+campaign+'/obs/aircraft/pekour-cvi/'
    cpcpath = package_path+'data/'+campaign+'/obs/aircraft/mei-cpc/'
    ccnpath = package_path+'data/'+campaign+'/obs/aircraft/mei-ccn/'
    amspath = package_path+'data/'+campaign+'/obs/aircraft/shilling-ams/'
    wcmpath = package_path+'data/'+campaign+'/obs/aircraft/matthews-wcm/'
    # surface measurements
    smps_pnnl_path = package_path+'data/'+campaign+'/obs/surface/pnnl-smps/'
    smps_bnl_path = package_path+'data/'+campaign+'/obs/surface/bnl-smps/'
    nanosmps_bnl_path = package_path+'data/'+campaign+'/obs/surface/bnl-nanosmps/'
    uhsassfcpath = package_path+'data/'+campaign+'/obs/surface/arm-uhsas/'
    cpcsfcpath = package_path+'data/'+campaign+'/obs/surface/arm-cpc/'
    cpcusfcpath = package_path+'data/'+campaign+'/obs/surface/arm-cpcu/'
    ccnsfcpath = package_path+'data/'+campaign+'/obs/surface/arm-ccn/'
    metpath = package_path+'data/'+campaign+'/obs/surface/arm-met/'
    acsmpath = package_path+'data/'+campaign+'/obs/surface/arm_acsm/'
    # vertical profile measurements
    armbepath = package_path+'data/'+campaign+'/obs/profile/sgparmbecldrad/'
    
    # PBLH data needed for plot_flight_pdf_percentile_SeparatePBLH_hiscale.py only
    pblhpath = package_path+'data/'+campaign+'/obs/profile/arm-pblh/'
    dlpath = package_path+'data/'+campaign+'/obs/profile/dl-pblh/'
    
    # model path
    # pre-processed model path
    E3SM_sfc_path = package_path+'data/'+campaign+'/model/surface/'
    E3SM_aircraft_path = package_path+'data/'+campaign+'/model/flighttrack/'
    E3SM_profile_path = package_path+'data/'+campaign+'/model/profile/'
    
    
elif campaign=='ACEENA':
    site='ENA'
    # lat/lon for ENA
    lat0 = 39.09527
    lon0 = 360-28.0339
    # bin of flight heights to calculate percentiles
    height_bin = np.arange(100,4300,300)
    
    # time periods for IOPs. needed in preprocessing of surface data
    if IOP=='IOP1':
        start_date='2017-06-20'
        end_date='2017-07-20'
    elif IOP=='IOP2':
        start_date='2018-01-21'
        end_date='2018-02-19'
    
    # observational data path. 
    # aircraf measurements merged_bin data are used for all plot_flight_*.py to provide flight/cloud/CVI info
    merged_size_path=package_path+'data/'+campaign+'/obs/aircraft/merged_bin/'
    iwgpath = package_path+'data/'+campaign+'/obs/aircraft/IWG/'
    fimspath = package_path+'data/'+campaign+'/obs/aircraft/FIMS/'
    pcasppath = package_path+'data/'+campaign+'/obs/aircraft/pcasp_g1/'
    cvipath = package_path+'data/'+campaign+'/obs/aircraft/inletcvi/'
    opcpath = package_path+'data/'+campaign+'/obs/aircraft/opciso/'
    cpcpath = package_path+'data/'+campaign+'/obs/aircraft/cpc_aaf/'
    ccnpath = package_path+'data/'+campaign+'/obs/aircraft/ccn_aaf/'
    amspath = package_path+'data/'+campaign+'/obs/aircraft/shilling-hrfams/'
    wcmpath = package_path+'data/'+campaign+'/obs/aircraft/wcm_ACEENA/'
    # surface measurements
    uhsassfcpath = package_path+'data/'+campaign+'/obs/surface/arm_uhsas/'
    cpcsfcpath = package_path+'data/'+campaign+'/obs/surface/arm_cpcf/'
    cpcusfcpath = 'N/A'
    ccnsfcpath = package_path+'data/'+campaign+'/obs/surface/arm_aosccn1/'
    metpath = package_path+'data/'+campaign+'/obs/surface/arm_met/'
    acsmpath = package_path+'data/'+campaign+'/obs/surface/arm_acsm/'
    # vertical profile measurements
    armbepath = package_path+'data/'+campaign+'/obs/profile/enaarmbecldrad/'
    
    # model path
    # pre-processed model path    
    E3SM_sfc_path = package_path+'data/'+campaign+'/model/surface/'
    E3SM_aircraft_path = package_path+'data/'+campaign+'/model/flighttrack/'
    E3SM_profile_path = package_path+'data/'+campaign+'/model/profile/'
    
elif campaign=='MAGIC':
    site='MAG'
    
    # bin of latitude to calculate ship track composite
    latbin = np.arange(21.5,34,1)
    
    # reference lat/lon
    lat0=30.
    lon0=230.
    
    # observational data path. 
    # ship measurements
    shipmetpath=package_path+'data/'+campaign+'/obs/ship/raynolds-marmet/'
    shipccnpath=package_path+'data/'+campaign+'/obs/ship/magaosccn100M1.a1/'
    shipcpcpath=package_path+'data/'+campaign+'/obs/ship/magaoscpcfM1.a1/'
    shipmwrpath=package_path+'data/'+campaign+'/obs/ship/magmwrret1liljclouM1.s2/'
    shipuhsaspath=package_path+'data/'+campaign+'/obs/ship/magaosuhsasM1.a1/'
    
    # model path
    # pre-processed model path    
    E3SM_ship_path = package_path+'data/'+campaign+'/model/shiptrack/'
    E3SM_profile_path = package_path+'data/'+campaign+'/model/profile/'
    
elif campaign=='MARCUS':
    site='MAR'
    
    # bin of latitude to calculate ship track composite
    latbin = np.arange(-68.5,-42,1)
    
    # reference lat/lon
    lat0=-40.
    lon0=120.
    
    
    # observational data path. 
    # ship measurements
    shipmetpath=package_path+'data/'+campaign+'/obs/ship/maraadmetX1.b1/'
    shipccnpath=package_path+'data/'+campaign+'/obs/ship/maraosccn1colavgM1.b1/'
    shipcpcpath=package_path+'data/'+campaign+'/obs/ship/maraoscpcf1mM1.b1/'
    shipmwrpath=package_path+'data/'+campaign+'/obs/ship/marmwrret1liljclouM1.s2/'
    shipuhsaspath=package_path+'data/'+campaign+'/obs/ship/maraosuhsasM1.a1/'
    
    # model path
    # pre-processed model path    
    E3SM_ship_path = package_path+'data/'+campaign+'/model/shiptrack/'
    E3SM_profile_path = package_path+'data/'+campaign+'/model/profile/'
    
elif campaign=='CSET':
    # bin of flight heights to calculate percentiles
    height_bin = np.arange(200,8000,400)
    # bin of latitude to calculate composite percentiles, same as MAGIC
    latbin = np.arange(22.5,39,1)
    
    # lat/lon at the airport
    lat0 = 38.5564
    lon0 = 360-121.3120
    
    # observational data path. 
    # aircraft measurements
    RFpath=package_path+'data/'+campaign+'/obs/aircraft/aircraft_lowrate/'
    ccnpath='N/A'
    
    # model path
    # pre-processed model path    
    E3SM_aircraft_path = package_path+'data/'+campaign+'/model/flighttrack/'
    E3SM_profile_path = package_path+'data/'+campaign+'/model/profile/'
    
elif campaign=='SOCRATES':
    # bin of flight heights to calculate percentiles
    height_bin = np.arange(200,8000,400)
    # bin of latitude to calculate composite percentiles
    latbin = np.arange(-63.5,-42,1)
    # height_bin = np.arange(200,7000,400)
    # lat/lon at the airport
    lat0 = -42.8371
    lon0 = 147.5054
    
    # observational data path. 
    # aircraft measurements
    RFpath=package_path+'data/'+campaign+'/obs/aircraft/aircraft_lowrate/'
    ccnpath=package_path+'data/'+campaign+'/obs/aircraft/CCN/'
    
    # model path
    # pre-processed model path    
    E3SM_aircraft_path = package_path+'data/'+campaign+'/model/flighttrack/'
    E3SM_profile_path = package_path+'data/'+campaign+'/model/profile/'
    
else:
    raise ValueError("does not recognize this campaign: "+campaign)

