#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=50G
#SBATCH --tmp=50G
#SBATCH -t 1:00:00

module load R/4.2.1

tool=$tool
root=$root
funcs=$funcs
infile=$infile
genome=$genome
stats=$stats
outpath=$outpath

echo "The tool is $tool"
echo "The tool is $root"
echo "The funcs is $funcs"
echo "The infile is $infile"
echo "The genome is $genome"
echo "The tool is $stats"
echo "The outpath is $outpath"

Rscript $tool $root $funcs $infile $genome $stats $outpath