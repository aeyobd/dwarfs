#for d in "sculptor" "ursa_minor" "fornax" "sextans1" "leo1" "leo2" "carina" "draco"; do 
#    echo "$d"
#    (cd "$d/density_profiles" && julia ../../make_profile.jl jax.toml;)
#done

# for d in */; do 
#     echo "$d"
#     if [ -f "$d/density_profiles/jax_2c.toml" ]; then
#         (cd "$d/density_profiles" && julia ../../make_profile.jl jax_2c.toml;)
#     fi
# done
