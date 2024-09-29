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

# Making directories read and execute only for all users
find "$DIRECTORY" -type d -exec chmod 555 {} \;

# Making files read-only for all users
find "$DIRECTORY" -type f -exec chmod 444 {} \;

echo "All files and directories in $DIRECTORY are now read-only."
