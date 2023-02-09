#!/bin/bash

#SBATCH -N 1 -c 2
#SBATCH --mem=32G
#SBATCH --tmp=32G
#SBATCH -t 24:00:00

root_path=/hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/ # args[1]
main_path=new_script/script/R/pipeline/v230116/Pipeline_v17_main_230116.R 
funcs_path=new_script/script/R/pipeline/v230116/Pipeline_v17_funcs_230116.R # args[2]
stats_path=new_script/script/R/stats/Pipeline_stats.Rmd # args[5]

tool="$root_path$main_path" # args[0]

# If all files are in a single folder:
# : <<'END'
infile_dir=/hpf/largeprojects/tcagstor/projects/wgs_rae_yeung/ICA_V2/MIS-C_269samples_20221121/annotated/ 
output_dir=results/marla_230208/ 

mkdir "$root_path$output_dir"
files=$(find $infile_dir -name '*tsv.gz' -exec basename {} \;)

for file in $files
do  
    infile_path="$infile_dir$file" # args[3]
    genome=$(echo $file | awk -F '.' '{print $1}') # args[4]
    output_path="$output_dir$genome" # args[6]
    mkdir "$output_path"
    echo "The genome is $genome"
    echo "The infile_path is $infile_path"
    echo "The output_path is $output_path"
    sbatch --export=tool=$tool,root=$root_path,funcs=$funcs_path,infile=$infile_path,genome=$genome,stats=$stats_path,outpath=$output_path /hpf/largeprojects/tcagstor/scratch/kara.han/PrioritizationPipeline/new_script/script/slurm/v230116/get_prioritization_results_230116.sh
done
# END