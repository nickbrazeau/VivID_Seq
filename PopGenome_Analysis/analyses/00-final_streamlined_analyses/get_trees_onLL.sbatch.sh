#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=24g
#SBTACH --job-name=gettrees
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/PopGenome_Analysis"); source("analyses/00-final_streamlined_analyses/fit_tree_backend.R")'
