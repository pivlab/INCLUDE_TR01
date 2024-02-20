#!/bin/env bash
# Download GTEx data
set -u

# Ensure the script is executed with two arguments.
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <URL> <output_path_without_gz_extension>"
    exit 1
fi

# Assign command line arguments to variables.
gtex_link=$1
output=$2

# Validate the URL to ensure it's not empty.
if [[ -z "$gtex_link" ]]; then
    echo "Error: The URL is empty."
    exit 1
fi

# Check if wget is installed.
if ! command -v wget &> /dev/null; then
    echo "Error: wget could not be found. Please install wget and try again."
    exit 1
fi

# Append .gz to the output filename.
output=${output}.gz

# Extract the directory path from the provided output path.
output_dir=$(dirname ${output})

# Create the output directory if it doesn't already exist.
mkdir -p $output_dir

# Check if the target file already exists to avoid unnecessary downloads.
if [ -f "$output" ]; then
    echo "The file $output already exists. Skipping download."
else
    # Download the file from the specified URL and save it to the output path with the .gz extension.
    wget $gtex_link -O $output
    echo "Download completed: $output"
fi

# Decompress the downloaded file.
if [ -f "$output" ]; then
    gzip -d $output
    echo "Decompression completed: ${output%.gz}"
else
    echo "Error: The file $output does not exist and cannot be decompressed."
fi