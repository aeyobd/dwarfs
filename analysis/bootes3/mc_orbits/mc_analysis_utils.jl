using LilGuys


function get_coord_err(obs_props)
    df_out = Dict{String, Float64}()

    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
        df_out[key] = LilGuys.get_uncertainty(obs_props, key)
    end

    return df_out
end

