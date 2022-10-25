#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=100G
#SBATCH --tmp=100G
#SBATCH -t 1:00:00

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
