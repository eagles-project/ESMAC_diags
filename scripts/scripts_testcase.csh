#!/bin/csh


# This script makes user-specified plots comparing model simulations with ARM measurements.
#
#
# ############################################################
# # Step 1: change settings in settings.py                   #
# #    such as campaign name and model names                 #
# ############################################################
#
# # load modules. Tested version is Python 3.6.7 (Constance) and Python 3.8.5 (NERSC)
 module load python
 
# # this should be consistent with settings.py
# set field campaign name. More settings on specific field campaigns are in next section
set campaign = 'ACEENA'   # HISCALE, ACEENA, CSET, SOCRATES, MAGIC, MARCUS
#foreach campaign ('MAGIC' 'MARCUS' 'CSET' 'SOCRATES')
 
# set model names. up to three
# set Model_List = "['CTRL','Nuc','NucSoaCond']"
set Model_List = "['E3SMv1']"

# set plotting line colors for each model. corresponding to the Model_List
# set color_model = "['b','r','g']"
set color_model = "['r','b','g']"

# set IOP (or flight date) that the statistics, pdf and percentiles are averaged for.
# options: IOP1, IOP2, ALL, 20160830b
set IOP = 'IOP1'
# foreach IOP ('IOP1' 'IOP2')

# ############################################################
# # Step 2: update settings.py with the above settings       #
# ############################################################

sed -i "s/^campaign = .*/campaign = '$campaign'/" settings_testcase.py
sed -i "s/^Model_List = .*/Model_List = $Model_List/" settings_testcase.py
sed -i "s/^color_model = .*/color_model = $color_model/" settings_testcase.py
sed -i "s/^IOP = .*/IOP = '$IOP'/" settings_testcase.py

# remove ^M in the file
sed -i "s/\r//g" settings_testcase.py
#

# ############################################################
# # Step 3: start plotting                                   #
# ############################################################

cp settings_testcase.py ../python/plotting/settings.py

echo '***** start testcase ********'
echo 'enter the plotting directory: ../python/plotting/'
cd ../python/plotting/

############# evaluate with flight measurements ################
echo '**********************************************'
echo 'plotting flight infomation'
python plot_flight_track_height.py
# timeseries comparison for each flight
echo '**********************************************'
echo 'plotting aerosol composition timeseries for flight'
python plot_flight_timeseries_AerosolComposition.py


###########################################################
#   end 
###########################################################
echo '*********** finishing testcase **************'
cd ../../scripts/

# end

exit

