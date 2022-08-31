#!/bin/csh

#SBATCH --nodes=1
#SBATCH --time=31:25:00
#SBATCH --qos=regular
#SBATCH --constraint=knl
#SBATCH --account=m3525
#SBATCH --output=a.out
#SBATCH --error=a.err

#module load python
#conda activate esmac_diags

#python prep_HISCALE_allobs.py
#python prep_ACEENA_allobs.py
#python prep_MAGIC_allobs.py
#python prep_MARCUS_allobs.py
#python prep_CSET_allobs.py
#python prep_SOCRATES_allobs.py

#python prep_ENA_allobs.py
#python prep_SGP_allobs.py

python prep_HISCALE_E3SM.py
python prep_ACEENA_E3SM.py
python prep_MAGIC_E3SM.py
python prep_MARCUS_E3SM.py
#python prep_CSET_E3SM.py
#python prep_SOCRATES_E3SM.py

#python prep_ENA_E3SM.py
#python prep_SGP_E3SM.py

exit
