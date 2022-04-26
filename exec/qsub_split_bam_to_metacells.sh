#!/bin/bash
unset module;
samtools_path="/home/feshap/src/samtools-1.15.1/samtools";

bamfile=$1
cell_names_folder=$2
base_output=$3
cell_name_files=($(ls $cell_names_folder))
cwd=$(pwd)
outfold=$cwd/${base_output}
bam_chroms=($($samtools_path view -H $bamfile | \
                egrep @SQ | \
                cut -f2 | \
                sed -e "s/^SN://" | \
                egrep chr))
printf "%s\n" "${bam_chroms[@]}" > bam_chroms.txt
printf "%s\n" "${cell_name_files[@]}" > cell_name_files.txt

if [ ! -d $outfold ] 
then 
    mkdir $outfold
fi
for cnf in ${cell_name_files[@]}; do
    cnf_fold=$outfold/${cnf%.txt};
    # echo $cnf_fold
    if [ ! -d $cnf_fold ] 
    then 
        mkdir $cnf_fold
    fi
done

for chrom in ${bam_chroms[@]}; do
    for cnf in ${cell_name_files[@]}; do
        qsub -terse -cwd -S /bin/bash -N ${base_output}_${cnf%.txt} -V ~/src/mcATAC/exec/qsub_command.sh $chrom $cnf $outfold $base_output $cell_names_folder $bamfile $samtools_path
    done
done

# bash $gnu_parallel_path bam_chroms.txt \
#                                     cell_name_files.txt \
#                                     $outfold \
#                                     $base_output \
#                                     $cell_names_folder \
#                                     $bamfile \
#                                     $samtools_path
wait
# for cnf in $cell_name_files; do
#     cnf_fold=$outfold/${cnf%.txt}
#     cd $cnf_fold
#     rm -f ${cnf%.txt}.bam
#     samtools merge -f $cwd/$outfold/${cnf%.txt}.bam ./*.bam &
#     cd $cwd
# done
# wait
# for cnf in $cell_name_files; do
#     cnf_fold=$outfold/$cnf
#     rm -rf $cnf_fold
# done
# rm -f bam_chroms.txt cell_name_files.txt