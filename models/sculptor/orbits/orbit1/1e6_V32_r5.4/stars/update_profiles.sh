#!/bin/bash

directory="./"
extension=".fits"
new_ext="_profile.toml"
scriptname="calc_2d_density.jl"

for file in "$directory"/*"$extension"; do
    filename=$(basename "$file")
    filename_noext="${filename%.*}"
    new_file="${directory}${filename_noext}${new_ext}"

    if [[ ! -f "$new_file" || "$file" -nt "$new_file" ]]; then
        echo "Processing $filename"
        
        "$scriptname" -s "$file" "$new_file"
    else
        echo "Skipping $filename"
    fi
done
