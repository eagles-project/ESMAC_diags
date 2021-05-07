#!/bin/csh

#SBATCH -A sooty2
#SBATCH -t 02:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J job_tang357
#SBATCH -p short
#SBATCH -o a.out
#SBATCH -e a.err


#
# ############################################################
# # Step 1: change settings in settings.py                   #
# #    such as campaign name and model names                 #
# ############################################################
#
# # load modules. Tested version is Python 3.6.7
 module load python/anaconda3.6
 
# # this should be consistent with settings.py
# set field campaign name. More settings on specific field campaigns are in next section
set campaign = 'HiScale'   # HiScale, ACEENA
 
# set model names. up to three
# set Model_List = "['CTRL','Nuc','NucSoaCond']"
set Model_List = "['CTRL','NucSoaCond']"

# set plotting line colors for each model. corresponding to the Model_List
# set color_model = "['b','r','g']"
set color_model = "['b','g']"

# set IOP (or flight date) that the statistics, pdf and percentiles are averaged for.
# options: IOP1, IOP2, ALL, 20160830b
set IOP = 'IOP1'

# ############################################################
# # Step 2: update settings.py with the above settings       #
# ############################################################

sed -i "s/^campaign = .*/campaign = '$campaign'/" settings.py
sed -i "s/^Model_List = .*/Model_List = $Model_List/" settings.py
sed -i "s/^color_model = .*/color_model = $color_model/" settings.py
sed -i "s/^IOP = .*/IOP = '$IOP'/" settings.py

# remove ^M in the file
sed -i "s/\r//g" settings.py
#

# ############################################################
# # Step 3: start plotting                                   #
# ############################################################

cp settings.py ../python/plotting/settings.py
echo '***** start plotting ********'
echo 'enter the plotting directory: ../python/plotting/'
cd ../python/plotting/

############# evaluate with flight measurements ################
echo 'plotting flight infomation'
python plot_flight_track_height.py
# timeseries comparison for each flight
echo 'plotting CN number timeseries for flight'
python plot_flight_timeseries_CN.py
echo 'plotting timeseries of aerosol PDF'
python contour_flight_timeseries_AerosolSize.py
echo 'plotting CCN number timeseries for flight'
python plot_flight_timeseries_CCN.py
echo 'plotting aerosol composition timeseries for flight'
python plot_flight_timeseries_AerosolComposition.py
# mean statistics for the entire IOP
echo 'calculate statistics of CN number for flight'
python calc_statistic_flight_CN.py
echo 'plotting mean aerosol PDF'
python plot_flight_pdf_AerosolSize.py
# vertical profiles or percentiles
echo 'plotting percentiles of CN number with height'
python plot_flight_percentile_z_CN.py
echo 'plotting percentiles of CCN number with height'
python plot_flight_percentile_z_CCN.py
echo 'plotting percentiles of Aerosol COmposition with height'
python plot_flight_percentile_z_AerosolComposition.py
# vertical profiles or percentiles of cloud
echo 'plotting vertical profile of cloud'
python plot_profile_cloud.py
echo 'plotting vertical profile of cloud frequency'
python plot_flight_profile_z_CldFreq.py
echo 'plotting vertical profile of cloud LWC'
python plot_flight_profile_z_LWC.py
# specific plotting separated by PBLH or clouds
if ($campaign == 'ACEENA') then
    echo 'plotting aerosol PDF and percentile separated by near surface, near cloud, above cloud'
    python plot_flight_pdf_percentile_SeparateCloud_aceena.py
endif
if ($campaign == 'HiScale') then
    echo 'plotting aerosol PDF and percentile separated by below/above PBLH'
    python plot_flight_pdf_percentile_SeparatePBLH_hiscale.py
endif

############# evaluate with surface measurements ################


# timeseries
echo 'plotting CN number  timeseries at surface'
python plot_sfc_timeseries_CN.py
echo 'plotting CCN number timeseries at surface'
python plot_sfc_timeseries_CCN.py
echo 'plotting aerosol composition timeseries at surface'
python plot_sfc_timeseries_AerosolComposition.py
echo 'plotting timeseries of aerosol PDF'
python contour_sfc_timeseries_AerosolSize.py
# diurnal cycle
echo 'plotting diurnalcycle of CN number at surface'
python plot_sfc_diurnalcycle_CN.py
echo 'plotting diurnalcycle of CCN number at surface'
python plot_sfc_diurnalcycle_CCN.py
echo 'plotting diurnalcycle of Aerosol COmposition at surface'
python plot_sfc_diurnalcycle_AerosolComposition.py
echo 'plotting diurnal cycle of aerosol PDF'
python contour_sfc_diurnalcycle_AerosolSize.py
# mean statistics
echo 'calculate statistics of CN number at surface'
python calc_statistic_sfc_CN.py
echo 'plotting mean aerosol PDF'
python plot_sfc_pdf_AerosolSize.py

###########################################################
#   end 
###########################################################
echo '*********** end plotting **************'
cd ../../scripts/
exit
