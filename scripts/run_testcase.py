# settings of the aerosol diagnostic package

from matplotlib.pyplot import plot
import numpy as np
import esmac_diags

from esmac_diags.plotting import plot_flight_pdf_AerosolSize, plot_flight_track_height


settings = {}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# global settings
############ these settings will be replaced by the settings in scripts_*.csh #############
# set field campaign name. More settings on specific field campaigns are in next section
settings['campaign'] = 'ACEENA'
# set model names. up to three
settings['Model_List'] = ['E3SMv1']
# set line colors for each model. corresponding to the Model_List
settings['color_model'] = ['r','b','g']
# set IOP that the statistics, pdf and percentiles are averaged for. Only available for HISCALE and ACEENA
# IOP1/IOP2 
settings['IOP'] = 'IOP1'
############ these settings will be replaced by the settings in scripts_*.csh #############


# path of the diagnostic package
package_path = '../../'

# # path of E3SM model data (h3) for preprocessing. list with the same length of Model_List
# settings['E3SM_h3_path']=[]
# settings['E3SM_h3_filehead']=[]     # filename before .cam.h3.yyyy-mm-dd.00000.nc
# for mm in settings['Model_List']:
#     E3SM_h3_path.append('/global/cscratch1/sd/sqtang/EAGLES/E3SM_output/E3SMv1_h3/')
#     if campaign=='MAGIC':
#         E3SM_h3_filehead.append(mm+'_2012-2013')
#     else:
# #        E3SM_h3_filehead.append(mm+'_2014-2018')
#         E3SM_h3_filehead.append(mm)
#     #E3SM_h3_path.append('/qfs/projects/eagles/zhan524/simulations/compy_F20TRC5-CMIP6_ne30_EG1_R2_'+mm+'/h3/')
#     #E3SM_h3_filehead.append('compy_F20TRC5-CMIP6_ne30_EG1_R2_'+mm)

# path of output figures
settings['figpath_aircraft_timeseries'] = package_path+'testcase/figures/'
settings['figpath_aircraft_statistics'] = package_path+'testcase/figures/'
settings['figpath_ship_timeseries'] = package_path+'testcase/figures/'
settings['figpath_ship_statistics'] = package_path+'testcase/figures/'
settings['figpath_sfc_timeseries'] = package_path+'testcase/figures/'
settings['figpath_sfc_statistics'] = package_path+'testcase/figures/'
settings['figpath_profile_timeseries'] = package_path+'testcase/figures/'
settings['RFpath']=package_path+'data/'+settings['campaign']+'/obs/aircraft/aircraft_lowrate/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings for different field campaigns

# set location and time information
if  settings['campaign']=='ACEENA':
    settings['site']='ENA'
    # lat/lon for ENA
    settings['lat0'] = 39.09527
    settings['lon0'] = 360-28.0339
    # bin of flight heights to calculate percentiles
    settings['height_bin'] = np.arange(100,4300,300)
    
    # time periods for IOPs. needed in preprocessing of surface data
    if settings['IOP']=='IOP1':
        settings['start_date']='2017-06-30'
        settings['end_date']='2017-06-31'
    elif settings['IOP']=='IOP2':
        settings['start_date']='2018-01-21'
        settings['end_date']='2018-02-19'
    
    # observational data path. 
    # aircraf measurements merged_bin data are used for all plot_flight_*.py to provide flight/cloud/CVI info
    settings['merged_size_path']=package_path+'testcase/data/obs/'
    settings['iwgpath'] = package_path+'testcase/data/obs/'
    settings['cvipath'] = package_path+'testcase/data/obs/'
    settings['amspath'] = package_path+'testcase/data/obs/AMS/'
    
    # model path
    # pre-processed model path    
    settings['E3SM_aircraft_path'] = package_path+'testcase/data/model/'
   
    
else:
    raise ValueError("Test case should only for ACEENA. Current campaign is: "+campaign)

# Long list of blah.run_plot()

### Create our actual plots
plot_flight_pdf_AerosolSize.run_plot(settings)
plot_flight_track_height.run_plot(settings)


