#!/bin/bash

module load R

Rscript clustering.R
#Rscript annotation.R

#sbatch --cpus-per-task=50 --mem=100g --time=88:00:00 SC_data_sbatch.sh
#sbatch --cpus-per-task=20 --mem=50g --time=24:00:00 sbatch.sh