#!/bin/bash
outfold=$3 
base_output=$4
cell_names_folder=$5
bamfile=$6
samtools_path=$7
gparallel_path=$8
$gparallel_path $samtools_path view -bo $outfold/{2}/{1}.bam -D CB:$cell_names_folder/{2} $bamfile {1} :::: ${1} ${2}