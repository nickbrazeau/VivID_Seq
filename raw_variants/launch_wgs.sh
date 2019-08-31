#! /bin/bash

ROOT=/proj/ideel/julianog/users/apm/tz_zb # root directory for project (non-scratch)
WD=$SCRATCH/tz_zb # working directory for alignments (scratch)
NODES=1028 # max number of cluster nodes
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/align_wgs.snake \
	--configfile config_wgs.yaml \
	--printshellcmds \
	--directory $WD \
	--cluster launch.py \
	-j $NODES \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	#--dryrun
