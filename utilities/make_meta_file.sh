#!/bin/bash

set -e

input_folders=()
output_file=""

while getopts "i:o:" opt; do
  case $opt in
    i)
      IFS=',' read -r -a input_folders <<< "$OPTARG"
      ;;
    o)
      output_file="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z "$output_file" ] || [ ${#input_folders[@]} -eq 0 ]; then
  echo "Usage: $0 -o output_file -i input_folder1[,input_folder2,...]"
  echo "Example: $0 -o /path/to/output_file.txt -i /path/to/input_folder1,/path/to/input_folder2"
  exit 1
fi

> "$output_file"

for folder in "${input_folders[@]}"; do
  for file in "$folder"/*; do
    if [ -e "$file" ]; then
      path=$file
      label=$(basename "$folder")
      echo -e "$path\t$label" >> "$output_file"
    fi
  done
done