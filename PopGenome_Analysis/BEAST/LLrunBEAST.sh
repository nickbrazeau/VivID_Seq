#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=11-00:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

# note, .xml file was generated using BEAUti GUI
# as recommended per BEAST tutorials

BEASTfile=/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/PopGenome_Analysis/BEAST/vivid_beast.xml

/usr/local/Cellar/beast/1.10.4/bin/beast -threads 1 $BEASTfile
