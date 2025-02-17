#!/bin/bash

#SBATCH --job-name=align_cellranger
#SBATCH --ntasks=32

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path> -g <genome_folder> -s <sample_name>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
genome_folder=""
sample_name=""

# Parse command line arguments
while getopts ":i:o:d:g:s:" opt; do
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
        s )
            sample_name=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ] \
|| [ -z "$genome_folder" ] || [ -z "$sample_name" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

docker load -i "$docker_image_path"
# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^cellranger$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run Cellranger
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder \
--security-opt seccomp=unconfined \
cellranger /bin/sh -c "/cellranger-9.0.0/cellranger count \
--id $sample_name \
--create-bam true \
--transcriptome /genome_folder \
--fastqs /input_folder \
--sample $sample_name \
--output-dir /output_folder/$sample_name;
chmod 777 -R /output_folder"