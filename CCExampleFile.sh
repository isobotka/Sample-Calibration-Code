#!/bin/bash
#SBATCH --account=def-cipriano   # replace this with your own account
#SBATCH --time=00-01:00          # time (DD-HH:MM)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

##Email Notifications
##SBATCH --mail-user=
##SBATCH --mail-type=ALL

module load gcc/9.3.0 r/4.1.2 gsl/2.6            # Adjust version and add the gcc module used for installing packages.
export R_LIBS=~/local/R_libs/

R CMD BATCH --no-save --no-restore Ontario-Diabetes-Microsim-Model.R