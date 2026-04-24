using LilGuys
import TOML
using Printf


include("table_utils.jl")



function (@main)(ARGS)
    print_all(["distance", "distance_em", "distance_ref", "pmra", "pmdec", "radial_velocity", "radial_velocity_em", "radial_velocity_ref", "sigma_v", "sigma_v_ref", "sigma_v_em"])
end

function print_all(keys)
    for galaxy in ["sculptor", "ursa_minor", "fornax", "leo1", "leo2", "sextans1", "carina", "draco", "canes_venatici1", "crater2", "antlia2", "sagittarius"]
        df = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxy/observed_properties.toml"))
        df_pace = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/all/pace/properties_$galaxy.toml"))
        if galaxy == "sagittarius"
            df_mv = Dict()
        else
            df_mv = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/all/MV2020/properties_$galaxy.toml"))
        end

        println(galaxy)

        for key in keys
            if endswith(key, "_ref")
                @printf "%16s %8s %8s %8s\n" key get(df, key, NaN) get(df_mv, key, NaN) get(df_pace, key, NaN)
            else
                @printf "%16s %8.2f %8.2f %8.2f\n" key get(df, key, NaN) get(df_mv, key, NaN) get(df_pace, key, NaN)
            end
        end

    end
end

