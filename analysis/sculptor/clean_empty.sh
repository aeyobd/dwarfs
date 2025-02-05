#!/bin/bash

# Function to check if a directory contains only directories and no files
contains_only_directories() {
    local dir="$1"
    # Check if there are any files in the directory
    if [ -n "$(find "$dir" -type f)" ]; then
        return 1  # Contains files, return false
    fi

    return 0  # Empty directory, return false
}

remove_empty_dir() {
    local dir="$1"
    echo "Removing empty directory: $dir"
    for subdir in "$dir"/*/; do
        remove_empty_dir "$subdir"
    done

    rmdir "$dir"
}

# Function to recursively remove directories that only contain directories
remove_empty_dirs() {
    local dir="$1"
    
    for subdir in "$dir"/*/; do
        echo $subdir
        # Check if the current directory contains only directories
        if contains_only_directories "$subdir"; then
            # Prompt for confirmation before deleting
            read -p "Directory '$subdir' contains only directories. Delete it? [y/N] " confirm
            if [[ "$confirm" =~ ^[Yy]$ ]]; then
                echo "Deleting directory: $subdir"
                remove_empty_dir "$subdir"
            else
                echo "Skipping directory: $subdir"
            fi
        fi
    done
}
