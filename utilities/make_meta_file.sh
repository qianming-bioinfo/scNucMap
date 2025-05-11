#!/bin/bash

set -e

input_folders=()
output_file=""

show_help() {
  echo ""
  echo "Usage: sh $0 -o output_file -i input_folder1[,input_folder2,...]"
  echo "Example: sh $0 -o /path/to/output_file.txt -i /path/to/input_folder1,/path/to/input_folder2"
  echo ""
  echo "Options:"
  echo "  -i    Comma-separated list of input folders"
  echo "  -o    Path to the output file"
  echo "  -h    Show this help message and exit"
}

while getopts "i:o:h" opt; do
  case $opt in
    i)
      IFS=',' read -r -a input_folders <<< "$OPTARG"
      ;;
    o)
      output_file="$OPTARG"
      ;;
    h)
      show_help
      exit 0
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
  show_help
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