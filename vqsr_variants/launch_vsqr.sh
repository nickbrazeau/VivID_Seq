#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/vqsr_variants
WD=/pine/scr/n/f/nfb/Projects/VivID_Seq/vcfs_vqsr_variants/
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/run_vqsr.snake \
	--configfile $ROOT/config_vqsr.yaml \
	--printshellcmds \
	--directory $WD \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
#	--dryrun -p --unlock
