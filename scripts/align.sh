#!/bin/bash

#SBATCH --job-name=align
#SBATCH --ntasks=15

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path> -g <genome_folder> -b <barcode_folder>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
genome_folder=""
barcode_folder=""

# Parse command line arguments
while getopts ":i:o:d:g:b:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        b )
            barcode_folder=$OPTARG
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$input_folder" ] || [ -z "$output_folder" ] \
|| [ -z "$docker_image_path" ] || [ -z "$genome_folder" ] || [ -b "$barcode_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run STAR aligner
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder \
-v "$barcode_folder":/barcode_folder \
--security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "STAR \
--runThreadN 14 \
--soloType CB_UMI_Simple \
--genomeDir /genome_folder/STAR_index \
--readFilesIn /input_folder/R2.fastq.gz /input_folder/R1.fastq.gz \
--readFilesCommand zcat \
--soloCBwhitelist /barcode_folder/3M-3pgex-may-2023.txt \
--soloUMIlen 12 \
--soloBarcodeReadLength 0 \
--outFileNamePrefix /output_folder/ \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 50000000000 \
--outSAMattributes All CR CY UR UY GX GN gx gn CB UB sM sS sQ sF \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloFeatures GeneFull; \
samtools index /output_folder/Aligned.sortedByCoord.out.bam; \
chmod 777 -R /output_folder"