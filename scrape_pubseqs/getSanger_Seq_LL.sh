#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/scrape_pubseqs"); source("03-scrape_ncbi_sanger.R")'
