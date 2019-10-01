#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/vqsr_variants
WD=/pine/scr/n/f/nfb/Projects/VivID_Seq/vqsr_variants
NODES=18
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/run_combine_vsqr.snake \
	--configfile config_combine_vqsr.yaml \
	--printshellcmds \
	--directory $WD \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
#	--dryrun -p
