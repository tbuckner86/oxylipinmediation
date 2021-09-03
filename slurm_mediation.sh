#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=teresa_mediation
#SBATCH --output=teresa_mediation.log
time Rscript $HOME/mediation_TB.R