#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/labstrains_variants
WD=/pine/scr/n/f/nfb/Projects/VivID_Seq/vcfs_gatk_lab_strains/
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/get_labstrains.snake \
	--configfile $ROOT/config.yaml \
	--printshellcmds \
	--directory $WD \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
#	--dryrun -p
