#!/bin/bash

# Check if the directory path is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <directory_path>"
  exit 1
fi

DIRECTORY=$1

# Check if the directory exists
if [ ! -d "$DIRECTORY" ]; then
  echo "Error: Directory does not exist."
  exit 1
fi

# Finding all IPYNB files recursively within the directory and clearing their outputs
find "$DIRECTORY" -name "*.ipynb" -print0 | while IFS= read -r -d $'\0' notebook; do
  if [ -f "$notebook" ]; then
    # Clear the outputs using nbconvert
    jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace "$notebook"
    echo "Cleared outputs for $notebook"
  fi
done

echo "All notebook outputs in $DIRECTORY have been cleared."
