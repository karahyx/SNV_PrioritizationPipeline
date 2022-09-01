#!/bin/bash

#PBS -N priori
#PBS -l nodes=1:ppn=14
#PBS -l vmem=100g
#PBS -l mem=36g
#PBS -l file=100g
#PBS -l walltime=24:00:00
module load R/4.2.1

tool=$tool
infile=$infile
genome=$genome
outpath=$outpath

echo "The tool is $tool"
echo "The infile is $infile"
echo "The genome is $genome"
echo "The outpath is $outpath"

Rscript $tool $infile $genome $outpath
