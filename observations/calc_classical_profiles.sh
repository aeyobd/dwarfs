for d in "sculptor" "ursa_minor" "fornax" "sextans1" "leo1" "leo2" "carina" "draco"; do 
    echo "$d"
    cd "$d/density_profiles" 
    mv jax_LL_0.toml jax_LLR_0.toml

    for prof in "jax" "jax_LLR_0" "jax_nospace" "simple"; do
        echo "calculating $prof"
        rm $prof.log
        rm ${prof}_profile.toml
        julia ../../make_profile.jl $prof.toml > >(tee -a $prof.log) 2> >(tee -a $prof.log >&2)
    done
    #julia ../../make_profile.jl jax_LL_0.toml;
    cd -
done

# for d in */; do 
#     echo "$d"
#     if [ -f "$d/density_profiles/jax_2c.toml" ]; then
#         (cd "$d/density_profiles" && julia ../../make_profile.jl jax_2c.toml;)
#     fi
# done
