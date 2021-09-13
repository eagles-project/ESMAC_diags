#!/bin/csh


# Kai Zhang provides aerosol data from E3SM simulations that are extracted for aircraft tracks and for selected ARM sites. However, those data are column variables and cover ~3 degrees around the ARM sites. This script extract the variables along the flight track, or the surface variables at the ARM site for direct comparison with ARM measurements.
#
#
# ############################################################
# # Step 1: change settings in settings.py                   #
# #    such as campaign name and model names                 #
# #    campaign name also needs to be set here               #
# ############################################################
#
# # load modules. Tested version is Python 3.6.7 (Constance) and Python 3.8.5 (NERSC)
 module load python
 
# # this should be consistent with settings.py
# set field campaign name. More settings on specific field campaigns are in next section
set campaign = 'HISCALE'   # HISCALE, ACEENA, CSET, SOCRATES, MAGIC, MARCUS
 
# set model names. up to three
# set Model_List = "['CTRL','Nuc','NucSoaCond']"
set Model_List = "['EAMv1_CONUS_RRM']"

# set plotting line colors for each model. corresponding to the Model_List
# set color_model = "['b','r','g']"
set color_model = "['b','g']"

# set IOP (or flight date) that the statistics, pdf and percentiles are averaged for.
# options: IOP1, IOP2, ALL, 20160830b
# set IOP = 'IOP1'
foreach IOP ('IOP1' 'IOP2')


# ############################################################
# # Step 2: update settings.py with the above settings       #
# ############################################################

sed -i "s/^campaign = .*/campaign = '$campaign'/" settings.py
sed -i "s/^Model_List = .*/Model_List = $Model_List/" settings.py
sed -i "s/^color_model = .*/color_model = $color_model/" settings.py
sed -i "s/^IOP = .*/IOP = '$IOP'/" settings.py

# remove ^M in the file
sed -i "s/\r//g" settings.py

# ############################################################
# # Step 3: preprocessing obs and/or model data              #
# ############################################################

cp settings.py ../python/preprocessing/settings.py
echo '***** start preprocessing ********'
echo 'enter the preprocess directory: ../python/preprocessing/'
cd ../python/preprocessing/

# # for observation
# merge observed aerosol sizes from several aircraft instruments
#echo '**** merge aerosol size distribution: ****'
#if ($campaign == 'HiScale') then
#    python prep_obs_mergesize_HiScale.py
#else if ($campaign == 'ACEENA')
#    python prep_obs_mergesize_ACEENA.py
#else
#    echo 'ERROR: not recognize campaign name'
#endif
    
# for models
if (($campaign == 'HISCALE') || ($campaign == 'ACEENA')) then
echo '**********************************************'
echo '**** extract aerosol size distribution at Surface ****'
python prep_E3SM_sfc_bins.py
echo '**********************************************'
echo '**** extract all other variables at Surface ****'
python prep_E3SM_sfc_allvars.py
echo '**********************************************'
echo '**** extract vertical profiles at ARM site ****'
python prep_E3SM_profile_allvars.py
endif

if (($campaign == 'HISCALE') || ($campaign == 'ACEENA') || ($campaign == 'CSET') ||($campaign == 'SOCRATES')) then
echo '**********************************************'
echo '**** extract aerosol size distribution for aircraft tracks ****'
python prep_E3SM_flighttrack_bins.py
echo '**********************************************'
echo '**** extract all other variables for aircraft tracks ****'
python prep_E3SM_flighttrack_allvars.py
endif
        
if (($campaign == 'MAGIC') ||($campaign == 'MARCUS')) then
echo '**********************************************'
echo '**** extract all other variables along ship tracks ****'
python prep_E3SM_shiptrack_allvars.py
echo '**********************************************'
echo '**** extract aerosol size distribution along ship tracks ****'
python prep_E3SM_shiptrack_bins.py
echo '**********************************************'
echo '**** extract vertical profiles along ship tracks ****'
python prep_E3SM_shiptrack_profiles.py
endif

# ############################################################
# # Step 4: end of preprocessing                             #
# ############################################################

echo '***** finished ******'

cd ../../scripts/

end

exit

