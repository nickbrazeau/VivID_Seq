#!/usr/bin/env bash
rsync -avr /Users/nickbrazeau/Documents/GitHub/VivID_Seq/PopGenome_Analysis/data/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/PopGenome_Analysis/data/
# from longleaf
rsync -av nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/PopGenome_Analysis/all_smpl_tree_results ./
