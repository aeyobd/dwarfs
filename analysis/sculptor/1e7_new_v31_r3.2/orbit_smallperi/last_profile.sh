set -xe 

out_path="combined.hdf5"

idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' orbital_properties.toml)

mass_profile.jl $out_path/$idx_f mass_profile_unbound.toml -u

