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

# Reverting directories to read, write, and execute for owner and read, execute for others
find "$DIRECTORY" -type d -exec chmod 755 {} \;

# Reverting files to read and write for owner and read for others
find "$DIRECTORY" -type f -exec chmod 644 {} \;

echo "Permissions have been reverted for $DIRECTORY."
