#!/bin/bash
unset module;
outfold=$3;
base_output=$4;
cell_names_folder=$5;
bamfile=$6;
samtools_path=$7;
echo $samtools_path 2
$samtools_path view -bo $outfold/${2%.txt}/${1}.bam -D CB:$cell_names_folder/${2} $bamfile $1;