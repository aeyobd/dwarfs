using LilGuys
import TOML
using Printf


include("table_utils.jl")



function (@main)(ARGS)
    print_all(["ra", "dec", "r_h", "ellipticity", "position_angle"])
end

function print_all(keys)
    for galaxy in ["sculptor", "ursa_minor", "fornax", "leo1", "leo2", "sextans1", "carina", "draco", "canes_venatici1", "crater2"]
        df = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxy/observed_properties.toml"))

        println(galaxy)

        for key in keys
            @printf "%16s %8.4f\n" key get(df, key, NaN) 
        end

    end
end

