#!/bin/csh

#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --qos=regular
#SBATCH --constraint=knl
#SBATCH --account=m3525
#SBATCH --output=a.out
#SBATCH --error=a.err

#module load python
#conda activate esmac_diags

python run_plotting.py

exit
