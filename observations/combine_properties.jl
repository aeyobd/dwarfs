import LilGuys as lguys
using CSV, DataFrames
import TOML

required_keys = ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity", "sigma_v", "Mv"]

function get_err(df, key)
    if haskey(df,  key * "_err")
        return df[key * "_err"]
    elseif haskey(df, key * "_em")
        return max(df[key * "_em"], df[key * "_ep"])
    else
        return 0
    end
end


function combine_properties()
    df_all = DataFrame()

    for galaxyname in readdir(".")
        dir = galaxyname
        file = dir * "/observed_properties.toml"

        if !isdir(dir)
            continue
        end

        if !isfile(file)
            println("Skipping $dir")
            continue
        end

        df = TOML.parsefile(file)

        if all([haskey(df, k) for k in required_keys])

            row = DataFrame(
                galaxyname=galaxyname, 
                ra=df["ra"], 
                ra_err=get_err(df, "ra"),
                dec=df["dec"], 
                dec_err=get_err(df, "dec"),
                distance=df["distance"],
                distance_modulus_err = get_err(df, "distance_modulus"),
                pmra=df["pmra"], 
                pmra_err = get_err(df, "pmra"),
                pmdec=df["pmdec"], 
                pmdec_err = get_err(df, "pmdec"),
                radial_velocity=df["radial_velocity"], 
                radial_velocity_err = get_err(df, "radial_velocity"),
                sigma_v = df["sigma_v"],
                sigma_v_err = get_err(df, "sigma_v"),
                Mv=df["Mv"],
                ellipticity=df["ellipticity"], 
            )

            append!(df_all, row)

        else
            println("incomplete info for $dir")
        end
    end


    df_all[!, :index] = 1:length(df_all.galaxyname)
    
    CSV.write("observed_properties.csv", df_all)
    
    filt = @. !isnan(df_all.distance) & !isnan(df_all.pmra) & !isnan(df_all.pmdec) & !isnan(df_all.radial_velocity)  & !isnan(df_all.Mv)
    df_complete = df_all[filt, :]
    CSV.write("observed_properties_complete.csv", df_complete)


    return
end


function (@main)(ARGS)
    combine_properties()
end
