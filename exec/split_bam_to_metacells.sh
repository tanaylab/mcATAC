#!/bin/bash

samtools_path="/home/feshap/src/samtools-1.15.1/samtools"

bamfile=$1
cell_names_folder=$2
base_output=$3
gnu_parallel_path=$4
cell_name_files=$(ls $cell_names_folder)
outfold=${base_output}
cwd=$(pwd)
bam_chroms=$($samtools_path view -H $bamfile | \
                egrep @SQ | \
                cut -f2 | \
                sed -e "s/^SN://" | \
                egrep chr)
printf "%s\n" "${bam_chroms[@]}" > bam_chroms.txt
printf "%s\n" "${cell_name_files[@]}" > cell_name_files.txt

if [ ! -d $outfold ] 
then 
    mkdir $outfold
fi
for cnf in $cell_name_files; do
    cnf_fold=$outfold/$cnf
    if [ ! -d $cnf_fold ] 
    then 
        mkdir $cnf_fold
    fi
done

bash $gnu_parallel_path bam_chroms.txt \
                                    cell_name_files.txt \
                                    $outfold \
                                    $base_output \
                                    $cell_names_folder \
                                    $bamfile \
                                    $samtools_path
wait
for cnf in $cell_name_files; do
    cnf_fold=$outfold/$cnf
    cd $cnf_fold
    rm -f ${cnf%.txt}.bam
    samtools merge -f $cwd/$outfold/${cnf%.txt}.bam ./*.bam &
    cd $cwd
done
wait
for cnf in $cell_name_files; do
    cnf_fold=$outfold/$cnf
    rm -rf $cnf_fold
done
rm -f bam_chroms.txt cell_name_files.txt