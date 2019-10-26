#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=11-00:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/PopGenome_Analysis"); source("analyses/04-ape_trees_all_backend.R")'
Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/PopGenome_Analysis"); source("analyses/05-ape_trees_clonal_backend.R")'
