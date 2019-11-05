#!/usr/bin/env bash

outdir=/pine/scr/n/f/nfb/Projects/VivID_Seq/Erbo1944
smpl=/pine/scr/n/f/nfb/Projects/VivID_Seq/wgs_se_improved_global/aln/merged/SRS1607662.bam
region=/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/regions/mtdna.bed
ref=/proj/ideel/resources/genomes/Pvivax/genomes/PvP01.fasta

mkdir -p $outdir

bcftools mpileup $smpl -R $region  --fasta-ref $ref -Ou | \
bcftools call -c  --variants-only > $outdir/Erbo1944.vcf
