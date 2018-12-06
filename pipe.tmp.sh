#!/bin/bash
#SBATCH -c 2                               # 1 core
#SBATCH -t 1-12:05                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH -p priopark                           # Run in short partition
#SBATCH --mem=20000M                           # Run in short partition
#SBATCH --account=park_contrib
#SBATCH -o ../hostname_brain_50x.%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH -o ../hostname_brain_50x.%j.err                # File to which STDOUT + STDERR will be written, including job ID in filename

module load gcc
module load R/3.4.1

##50-200x data:
#python feature_extraction.python3.6.1.py input.bed.2 200x_features_forR 200x_train/ /home/yd65/tools/MosaicHunter/resources/human_g1k_v37_decoy.fasta
Rscript Train_RFmodel.R demo/phasable_trainset demo/Phase_model.rds Phase
Rscript Train_RFmodel.R demo/phasable_trainset demo/Refine_model.rds Refine
