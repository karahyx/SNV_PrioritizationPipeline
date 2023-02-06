#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=40G
#SBATCH --tmp=40G
#SBATCH -t 24:00:00

module load R/4.2.1

tool=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/R/pipeline/v230116/Pipeline_v17_main_fam_230116.R

root_path=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/ # args[1]
funcs=script/R/pipeline/v230116/Pipeline_v17_funcs_fam_230116.R # args[2]
infile=/hpf/largeprojects/tcagstor/tools/dragen_dnaseq_prod/test_family_AZ/outbase/family_annotations/AZ/AZ.hard-filtered.vcf.gz.annovar.out_FINAL_rev27.7_hg38.tsv.gz # args[3]
child=NA24385_b # args[4]
pedigree=AZ # args[5]
stats=script/R/stats/Pipeline_stats.Rmd # args[6]
outpath=tests/pipeline_v230116/FAM # args[7]

mkdir "$root_path$outpath"

Rscript $tool $root_path $funcs $infile $child $pedigree $stats $outpath
