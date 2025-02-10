#!/bin/bash

sample_names=("A006850389_232503" "A006850389_232504" "A006850389_232505" "A006850389_232506")
slurm_log_folder="/data/public/jkoubele/ercc1/ercc1-mouse/slurm_logs"

for sample in "${sample_names[@]}"; do
    echo "$sample"
    sbatch --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
    --nodelist=beyer-n02,beyer-n03,beyer-n04,beyer-n05 \
    /data/public/jkoubele/ercc1/ercc1-mouse/scripts/align_cellranger.sh \
    -i "/data/public/jkoubele/ercc1/data/CCG_data/" \
    -o "/data/public/jkoubele/ercc1/data/cellranger_output/" \
    -d "/data/public/jkoubele/cellranger_docker/cellranger_docker_image.tar" \
    -g "/data/public/jkoubele/refdata-gex-mm10-2020-A/" \
    -s "$sample"
done