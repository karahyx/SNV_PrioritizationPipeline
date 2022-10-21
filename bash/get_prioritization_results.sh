#!/bin/bash

#PBS -N priori
#PBS -l nodes=1:ppn=2
#PBS -l vmem=100g
#PBS -l mem=100g
#PBS -l file=100g
#PBS -l walltime=24:00:00
#PBS -d /hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/tests/emedgene/

module load R/4.2.1

tool=$tool
funcs=$funcs
infile=$infile
genome=$genome
outpath=$outpath

echo "The tool is $tool"
echo "The funcs is $funcs"
echo "The infile is $infile"
echo "The genome is $genome"
echo "The outpath is $outpath"

Rscript $tool $funcs $infile $genome $outpath