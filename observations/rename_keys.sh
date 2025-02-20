#!/bin/bash

# Define the key substitutions (add more as needed)
substitutions=(
    # Rename 'rh' to 'r_h' while preserving whitespace
    's/^rh\([[:space:]]*\)=/r_h\1=/'
    's/^rh_em\([[:space:]]*\)=/r_h_em\1=/'
    's/^rh_ep\([[:space:]]*\)=/r_h_ep\1=/'
    's/^rh_err\([[:space:]]*\)=/r_h_err\1=/'
    's/^rh_ref\([[:space:]]*\)=/r_h_ref\1=/'
)

# Find and process all observed_properties.toml files
find . -type f -path "*/observed_properties.toml" -print0 | while IFS= read -r -d $'\0' file; do
    cp "$file" "$file.bak"
    echo "$file"

    # Build sed command with all substitutions
    sed_args=("-i")
    for sub in "${substitutions[@]}"; do
        sed_args+=("-e" "$sub")
    done
    sed_args+=("$file")

    # Execute sed command
    if sed --version &>/dev/null; then  # GNU sed
        sed "${sed_args[@]}"
    else                                # BSD sed (macOS)
        sed -i '' "${sed_args[@]:1}"    # Remove empty argument for -i
    fi
done

