using LilGuys
import TOML
using Printf


include("table_utils.jl")



function (@main)(ARGS)
    print_all(["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity", "r_h", "sigma_v"])
end

short_refs = Dict(
    "McConnachie2012AJ....144....4M" => 1,
    "An2024MNRAS.532.3713A" => 2,
    "MV2020" => 3,
    "Vivas2022ApJ...926...78V" => 4,
    "Oakes2022ApJ...929..116O" => 5,
    "Stetson2014PASP..126..616S" => 6,
    "Karczmarek2015AJ....150...90K" => 7,
    "Bhardwaj2024AJ....167..247B" => 8,
    "Kuehn2008ApJ...674L..81K" => 9,
    "Ji2021ApJ...921...32J" => 10,
    "munoz+2018" => 11,
    "Walker2009AJ....137.3100W" => 12,
    "Mateo2008ApJ...675..201M" => 13,
    "Spencer2017ApJ...836..202S" => 14,
   )

function to_name(s)
    s = replace(s, "_" => " ")
    s = titlecase(s)
    s = replace(s, "1" => " I")
    s = replace(s, "2" => " II")
end

function print_all(keys)
    for galaxy in ["sagittarius", "antlia2", "fornax", "leo1", "leo2", "sextans1", "carina", "draco", "canes_venatici1", "crater2"]
        df = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxy/observed_properties.toml"))

        refs = ""
        @printf "%16s" to_name(galaxy)
        for key in keys
            val = df[key]
            ref = df[key * "_ref"]

            if key == "ra"
                sval = @sprintf "%0.4f" val
            elseif key == "dec"
                sval = @sprintf " %+0.4f" val
            else
                sval = @sprintf "%0.2f" val
            end

            sval = "\$$sval\$"
            @printf " & %12s" sval
            refs *= "$(short_refs[ref]);"
        end
        refs = refs[1:end-1] # remove trailing colon

        println(" & $refs \\\\")
    end

    println()
    println()

end

