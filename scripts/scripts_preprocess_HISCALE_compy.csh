#!/bin/csh

#SBATCH -A sooty2
#SBATCH -t 02:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J job_tang357
#SBATCH -p short
#SBATCH -o a.out
#SBATCH -e a.err



# This script extract E3SM variables along the flight track, ship track or surface variables at the ARM site for direct comparison with ARM measurements.
#
#
# ############################################################
# # Step 1: change settings in settings.py                   #
# #    such as campaign name and model names                 #
# #    campaign name also needs to be set here               #
# ############################################################
#
# # load modules. Tested version is Python 3.6.7
#  module load python   # for NERSC
 module load python/anaconda3.6  # for constance 
 
# # this should be consistent with settings.py
# set field campaign name. More settings on specific field campaigns are in next section
set campaign = 'HISCALE'   # HISCALE, ACEENA,MAGIC,MARCUS,CSET,SOCRATES
 
# set model names. 
set Model_List = "['CTRL','nonucl','nonucl_pbl','NucSoaCond','NucSoaCond_noDpMin','NucSoaCond_DpMin5']"
# set Model_List = "['E3SMv1']"

# set plotting line colors for each model. corresponding to the Model_List
# set color_model = ['b','r','g']
set color_model = "['b','g']"

# set IOP (or flight date) that the statistics, pdf and percentiles are averaged for.
# options: IOP1, IOP2, ALL, 20160830b
#set IOP = 'IOP2'
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
echo '**** extract aerosol size distribution at Surface ****'
python prep_E3SM_sfc_bins.py
echo '**** extract all other variables at Surface ****'
python prep_E3SM_sfc_allvars.py
echo '**** extract vertical profiles at ARM site ****'
python prep_E3SM_profile_allvars.py
echo '**** extract aerosol size distribution for aircraft tracks ****'
python prep_E3SM_flighttrack_bins.py
echo '**** extract all other variables for aircraft tracks ****'
python prep_E3SM_flighttrack_allvars.py
        
#echo '**** extract all other variables along ship tracks ****'
#python prep_E3SM_shiptrack_allvars.py
#echo '**** extract aerosol size distribution along ship tracks ****'
#python prep_E3SM_shiptrack_bins.py
#echo '**** extract vertical profiles along ship tracks ****'
#python prep_E3SM_shiptrack_profiles.py

# ############################################################
# # Step 4: end of preprocessing                             #
# ############################################################

echo '***** finished ******'

cd ../../scripts/

end

exit

