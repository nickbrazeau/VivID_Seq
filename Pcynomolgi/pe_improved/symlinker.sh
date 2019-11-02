#!/bin/bash

##----------------------------------------------------------------------------##
## Authors: Nick Brazeau
## This script is used take the output from APM's rename_htsf_fastq.py
## and produce the file architecture that corresponds to the expectations for
## APM's align_wgs.snake
##----------------------------------------------------------------------------##

if [ $# -eq 0 ]; then
  echo "Usage: symlinker [OPTIONS]"
  echo "Try symlinker --help"
  exit 1
fi

##----------------------------------------------------------------------------##
## Collect inputs
##----------------------------------------------------------------------------#

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -h|--help)
    HELP=true
    ;;
    -I|--input)
    INPUT="$2"
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

##----------------------------------------------------------------------------##
## Help file
##----------------------------------------------------------------------------##


if [ "$HELP" = true ]; then

	echo "This is a symlink wrapper. Will create symlinks based on the input file from'"
	echo "APM's rename_htsf_fastq.py (syd.out)."
	echo "Usage: symlinker [OPTIONS]"
	echo ""
	echo "  -I, --input		Expects a tab delimited file with 3 columns. First column contains sample name, second column contains original fastq path (R1/R2 seperate lines), and third column is the path for the symlink."
	exit 0

fi



##----------------------------------------------------------------------------##
## write out symlinks
##----------------------------------------------------------------------------##
mkdir -p fastq
len=`wc -l $INPUT |awk '{print $1}'`
for i in `seq 1 ${len}`;
  do
    smpl=`cat $INPUT | cut -f1 -d " " | sed -n ${i}p`
    path=`cat $INPUT | cut -f2 -d " " | sed -n ${i}p`
    link=`cat $INPUT | cut -f3 -d " " | sed -n ${i}p`
    echo "make sampledir=${smpl} from=${path} to=fastq/${link}"
    mkdir -p fastq/${smpl}
    ln -s ${path} fastq/${link}

done
