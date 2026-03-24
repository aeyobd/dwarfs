using LilGuys
import TOML
using Printf
using OrderedCollections


include("table_utils.jl")



function (@main)(ARGS)
    print_all(["ra", "dec", "distance_modulus", "pmra", "pmdec", "radial_velocity", "r_h", "sigma_v"])
    println()
    println()
    for (ref, idx) in short_refs
        println("($idx) \\citet{$(short_cite_keys[ref])},")
    end
end

short_refs = OrderedDict(
    "MV2020" => 1,
    "McConnachie2012AJ....144....4M" => 2,
    "An2024MNRAS.532.3713A" => 3,
    "Vivas2022ApJ...926...78V" => 4,
    "Ji2021ApJ...921...32J" => 5,
    "Oakes2022ApJ...929..116O" => 6,
    "munoz+2018" => 7,
    "Walker2009AJ....137.3100W" => 8,
    "Stetson2014PASP..126..616S" => 9,
    "Mateo2008ApJ...675..201M" => 10,
    "Spencer2017ApJ...836..202S" => 11,
    "Karczmarek2015AJ....150...90K" => 12,
    "Bhardwaj2024AJ....167..247B" => 13,
    "Kuehn2008ApJ...674L..81K" => 14,
   )
short_cite_keys = Dict(
    "McConnachie2012AJ....144....4M" => "mcconnachie2012",
    "An2024MNRAS.532.3713A" => "an+walker+pace2024",
    "MV2020" => "MV2020a",
    "Vivas2022ApJ...926...78V" => "vivas+2022",
    "Oakes2022ApJ...929..116O" => "oakes+2022",
    "Stetson2014PASP..126..616S" => "stetson+2014",
    "Karczmarek2015AJ....150...90K" => "karczmarek+2015",
    "Bhardwaj2024AJ....167..247B" => "bhardwaj+2024",
    "Kuehn2008ApJ...674L..81K" => "kuehn+2008",
    "Ji2021ApJ...921...32J" => "ji+2021",
    "munoz+2018" => "munoz+2018",
    "Walker2009AJ....137.3100W" => "WMO2009",
    "Mateo2008ApJ...675..201M" => "mateo+olszewski+walker2008",
    "Spencer2017ApJ...836..202S" => "spencer+2017a",
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

        refs = []
        @printf "%16s" to_name(galaxy)
        for key in keys
            val = df[key]
            if key == "distance_modulus"
                ref = df["distance" * "_ref"]
            else
                ref = df[key * "_ref"]
            end

            if key == "ra"
                sval = @sprintf "%0.4f" val
            elseif key == "dec"
                sval = @sprintf " %+0.4f" val
            elseif key == "r_h"
                err = get_err(df, key)
                sval = @sprintf " %0.2f\\pm%0.2f" val err
            elseif key == "sigma_v"
                err = get_err(df, key)
                sval = @sprintf " %0.2f\\pm%0.2f" val err
            elseif key == "distance_modulus"
                err = get_err(df, key)
                sval = @sprintf " %0.2f\\pm%0.2f" val err
            else
                sval = @sprintf "%0.2f" val
            end

            sval = "\$$sval\$"
            @printf " & %12s" sval
            append!(refs, short_refs[ref])
        end
        refs = unique(refs)
        ref_str = join(refs, ";")

        println(" & $ref_str \\\\")
    end

    println()
    println()

end

function get_err(df, key)
    errkey = key * "_err"
    if errkey in keys(df)
        return df[errkey]
    end

    return (df[key * "_em"] + df[key * "_ep"])/2
end

