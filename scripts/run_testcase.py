"""
script to run a test case
compare the figures generated at testcase/figures/ with testcase/figures_verify
to makesure testcase works as expected
"""

from esmac_diags.plotting import plot_flight_timeseries_AerosolComposition, plot_flight_track_height


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# user-specified settings
settings = {}

# Please change the path of the diagnostic package to your own
package_path = '/global/homes/s/sqtang/EAGLES/ESMAC_diags/'

# set field campaign name. use ACEENA for the test case
settings['campaign'] = 'ACEENA'
# set model names. 
settings['Model_List'] = ['E3SMv1']
# set line colors for each model. corresponding to the Model_List
settings['color_model'] = ['r','b','g']
# set field campaign IOPs. Only used for HISCALE and ACEENA. 
# IOP1/IOP2 
settings['IOP'] = 'IOP1'


# path of output figures
settings['figpath_aircraft_timeseries'] = package_path+'testcase/figures/'
settings['figpath_aircraft_statistics'] = package_path+'testcase/figures/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# other settings for different field campaigns

# set location and time information
if  settings['campaign']=='ACEENA':
    settings['site']='ENA'
    # lat/lon for ENA
    settings['lat0'] = 39.09527
    settings['lon0'] = 360-28.0339
    
    # observational data path. 
    settings['merged_size_path']=package_path+'testcase/data/obs/'
    settings['iwgpath'] = package_path+'testcase/data/obs/'
    settings['cvipath'] = package_path+'testcase/data/obs/'
    settings['amspath'] = package_path+'testcase/data/obs/AMS/'
    
    # model path
    # pre-processed model path    
    settings['E3SM_aircraft_path'] = package_path+'testcase/data/model/'

else:
    raise ValueError("Test case should only for ACEENA. Current campaign is: "+settings['campaign'])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Long list of blah.run_plot()

### Create our actual plots
plot_flight_track_height.run_plot(settings)
plot_flight_timeseries_AerosolComposition.run_plot(settings)

