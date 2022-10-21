#!/bin/bash

#PBS -N priori
#PBS -l nodes=1:ppn=14
#PBS -l vmem=32g
#PBS -l mem=32g
#PBS -l walltime=24:00:00

tool=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/pipeline_new.R
funcs=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/pipeline_new_funcs.R

# If all files are in a single folder:
# : <<'END'
infile_dir=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/emedgene_final/
output_dir=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/tests/emedgene/

mkdir "$output_dir"
files=$(find $infile_dir -name '*.tsv.gz' -exec basename {} \;)

for file in $files
do  
    infile="$infile_dir$file"
    genome=$(echo $file | awk -F '_fx.' '{print $1}')
    outpath="$output_dir$genome"
    mkdir "$outpath"
    echo "The genome is $genome"
    echo "The infile is $infile"
    echo "The outpath is $outpath"
    qsub -v tool=$tool,funcs=$funcs,infile=$infile,genome=$genome,outpath=$outpath /hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/get_prioritization_results.sh
done
# END

# If files are in their separate folder and the folders are named after their sample name:
: <<'END'
infile_dir=/hpf/largeprojects/tcagstor/projects/wgs_rae_yeung/dragen_analysis/mis-c/ # path to folder
output_dir=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/results/mis-c/

mkdir "$output_dir"

for folder in "$infile_dir"*
do 
    genome=$(basename $folder)
    cur_dir="$infile_dir$genome"
    infile=$(find $cur_dir -name '*SUBSET*')
    outpath="$output_dir$genome"
    mkdir "$outpath"
    echo "The genome is $genome"
    echo "The cur_dir is $cur_dir"
    echo "The infile is $infile"
    echo "The outpath is $outpath"
    qsub -v tool=$tool,funcs=$funcs,infile=$infile,genome=$genome,outpath=$outpath /hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/get_prioritization_results.sh
done
END