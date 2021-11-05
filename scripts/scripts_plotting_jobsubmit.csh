#!/bin/csh

#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --qos=regular
#SBATCH --constraint=knl
#SBATCH --account=m3525
#SBATCH --output=a.out
#SBATCH --error=a.err

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
#set campaign = 'ACEENA'   # HISCALE, ACEENA, CSET, SOCRATES, MAGIC, MARCUS
foreach campaign ('MAGIC' 'MARCUS' 'CSET' 'SOCRATES' 'HISCALE' 'ACEENA')
 
# set model names. up to three
# set Model_List = "['CTRL','Nuc','NucSoaCond']"
set Model_List = "['E3SMv1']"

# set plotting line colors for each model. corresponding to the Model_List
# set color_model = "['b','r','g']"
set color_model = "['r','b','g']"

# set IOP (or flight date) that the statistics, pdf and percentiles are averaged for.
# options: IOP1, IOP2, ALL, 20160830b
set IOP = 'IOP1'
#foreach IOP ('IOP1' 'IOP2')

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
if (($campaign == 'HISCALE') || ($campaign == 'ACEENA') || ($campaign == 'CSET') ||($campaign == 'SOCRATES')) then
echo '**********************************************'
echo 'plotting flight infomation'
python plot_flight_track_height.py
# timeseries comparison for each flight
echo '**********************************************'
echo 'plotting CN number timeseries for flight'
python plot_flight_timeseries_CN.py
echo '**********************************************'
echo 'plotting timeseries of aerosol PDF'
python contour_flight_timeseries_AerosolSize.py
echo '**********************************************'
echo 'plotting CCN number timeseries for flight'
python plot_flight_timeseries_CCN.py
# mean statistics for the entire IOP
echo '**********************************************'
echo 'calculate statistics of CN number for flight'
python calc_statistic_flight_CN.py
echo '**********************************************'
echo 'plotting mean aerosol PDF'
python plot_flight_pdf_AerosolSize.py
# vertical profiles or percentiles
echo '**********************************************'
echo 'plotting percentiles of CN number with height'
python plot_flight_percentile_z_CN.py
echo '**********************************************'
echo 'plotting percentiles of CCN number with height'
python plot_flight_percentile_z_CCN.py
# vertical profiles or percentiles of cloud
echo '**********************************************'
echo 'plotting vertical profile of cloud frequency'
python plot_flight_profile_z_CldFreq.py
echo '**********************************************'
echo 'plotting vertical profile of cloud LWC'
python plot_flight_profile_z_LWC.py
if (($campaign == 'CSET') ||($campaign == 'SOCRATES')) then
    echo '**********************************************'
    echo 'plotting flight height percentile in latitude bins'
    python plot_flight_percentile_lat_cldfreq.py
    echo '**********************************************'
    echo 'plotting flight CN percentile in latitude bins'
    python plot_flight_percentile_lat_CN.py
    echo '**********************************************'
    echo 'plotting flight CCN percentile in latitude bins'
    python plot_flight_percentile_lat_CCN.py
endif
if (($campaign == 'HISCALE') ||($campaign == 'ACEENA')) then
    echo '**********************************************'
    echo 'plotting vertical profile of cloud'
    python plot_profile_cloud.py
    echo '**********************************************'
    echo 'plotting aerosol composition timeseries for flight'
    python plot_flight_timeseries_AerosolComposition.py
    echo '**********************************************'
    echo 'plotting percentiles of Aerosol COmposition with height'
    python plot_flight_percentile_z_AerosolComposition.py
endif
# specific plotting separated by PBLH or clouds
if ($campaign == 'ACEENA') then
    echo '**********************************************'
    echo 'plotting aerosol PDF and percentile separated by near surface, near cloud, above cloud'
    python plot_flight_pdf_percentile_SeparateCloud_aceena.py
endif
if ($campaign == 'HISCALE') then
    echo '**********************************************'
    echo 'plotting aerosol PDF and percentile separated by below/above PBLH'
    python plot_flight_pdf_percentile_SeparatePBLH_hiscale.py
endif
endif   # end evaluate with flight measurements


############# evaluate with surface measurements ################
if (($campaign == 'HISCALE') || ($campaign == 'ACEENA')) then
# timeseries
echo '**********************************************'
echo 'plotting CN number  timeseries at surface'
python plot_sfc_timeseries_CN.py
echo '**********************************************'
echo 'plotting CCN number timeseries at surface'
python plot_sfc_timeseries_CCN.py
echo '**********************************************'
echo 'plotting aerosol composition timeseries at surface'
python plot_sfc_timeseries_AerosolComposition.py
echo '**********************************************'
echo 'plotting timeseries of aerosol PDF'
python contour_sfc_timeseries_AerosolSize.py
# diurnal cycle
echo '**********************************************'
echo 'plotting diurnalcycle of CN number at surface'
python plot_sfc_diurnalcycle_CN.py
echo '**********************************************'
echo 'plotting diurnalcycle of CCN number at surface'
python plot_sfc_diurnalcycle_CCN.py
echo '**********************************************'
echo 'plotting diurnalcycle of Aerosol COmposition at surface'
python plot_sfc_diurnalcycle_AerosolComposition.py
echo '**********************************************'
echo 'plotting diurnal cycle of aerosol PDF'
python contour_sfc_diurnalcycle_AerosolSize.py
# mean statistics
echo '**********************************************'
echo 'calculate statistics of CN number at surface'
python calc_statistic_sfc_CN.py
echo '**********************************************'
echo 'plotting mean aerosol PDF'
python plot_sfc_pdf_AerosolSize.py
echo '**********************************************'
echo 'plotting fraction of surface aerosol composition'
python plot_sfc_pie_AerosolComposition.py
endif # end evaluate with surface measurements


############# evaluate with ship measurements ################
if (($campaign == 'MAGIC') ||($campaign == 'MARCUS')) then
# timeseries
echo '**********************************************'
echo 'plotting meterological fields timeseries for ship measurements'
python plot_ship_timeseries_met.py
echo '**********************************************'
echo 'plotting CN number timeseries for ship measurements'
python plot_ship_timeseries_CN.py
echo '**********************************************'
echo 'plotting CCN number timeseries for ship measurements'
python plot_ship_timeseries_CCN.py
# statistics
echo '**********************************************'
echo 'plotting meterological fields percentiles in latitude for ship measurements'
python plot_ship_percentile_lat_met.py
echo '**********************************************'
echo 'plotting CN number percentiles in latitude for ship measurements'
python plot_ship_percentile_lat_CN.py
echo '**********************************************'
echo 'plotting CCN number percentiles in latitude for ship measurements'
python plot_ship_percentile_lat_CCN.py
echo '**********************************************'
echo 'plotting LWP composition in latitude for ship measurements'
python plot_ship_percentile_lat_LWP.py
echo '**********************************************'
echo 'plotting mean aerosol size distribution for ship measurements'
python plot_ship_pdf_AerosolSize.py
echo '**********************************************'
echo 'calculate mean statistics of CN for ship measurements'
python calc_statistic_ship_CN.py
echo '**********************************************'
echo 'plotting timeseries of aerosol size distribution for ship measurements'
python contour_ship_timeseries_AerosolSize.py
endif # end evaluate with ship measurements

###########################################################
#   end 
###########################################################
echo '*********** end plotting **************'
cd ../../scripts/

end

exit

