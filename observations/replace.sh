shopt -s globstar

#
# Replace these with your actual key names
OLD_KEY="LLR_cut"
NEW_KEY="LLR_min"

# Process all TOML files in density_profiles directories
for file in **/density_profiles/*.toml; do
    if [[ -f "$file" ]]; then
         if grep -q -E "^[[:space:]]*${OLD_KEY}[[:space:]]*=" "$file"; then
            # Create backup only if key exists, then perform substitution
            sed -i.bak -E "s/^([[:space:]]*)${OLD_KEY}([[:space:]]*=)/\\1${NEW_KEY}\\2/" "$file"
            echo "Updated: $file (backup created at $file.bak)"
        else
            echo "Skipped: $file (no matching key found)"
        fi
    fi
done
