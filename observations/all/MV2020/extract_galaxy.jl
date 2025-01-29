import TOML
using LilGuys
import CSV
using DataFrames
using OrderedCollections

"""
    dm_to_distance(dm) 

Converts a distance modulus (magnitudes) to a distance in kpc.
"""
function dm_to_distance(dm)
    return 10^(dm / 5 + 1 - 3)
end


function main(ARGS)
    df = LilGuys.read_fits("NearbyGalaxies_Jan2021_PUBLIC.fits")

    for galaxyname in readdir("../../")
        try
            extract_galaxy(galaxyname, df)
        catch e
            println("Error processing $galaxyname: $e")
        end
    end
end

function select_galaxy(galaxyname, df)

    galaxyname = to_mv_name(galaxyname)

    if galaxyname ∉ df.GalaxyName
        if "*" * galaxyname ∈ df.GalaxyName
            galaxyname = "*" * galaxyname
        elseif galaxyname * "(1)" ∈ df.GalaxyName
            galaxyname = galaxyname * "(1)"
        elseif match(r"[a-zA-Z]1$", galaxyname) !== nothing && (galaxyname[1:end-1] * "(1)" ∈ df.GalaxyName)
            galaxyname = galaxyname[1:end-1] * "(1)"
        else
            @warn("Could not find $galaxyname")
            return
        end
    end

    matches = df[df.GalaxyName .== galaxyname, :]
    if size(matches, 1) != 1
        @warn ("Found $(size(matches, 1)) matches for $galaxyname")
        return 
    end
    galaxy = matches[1, :]

    println(galaxyname)
    return galaxy
end

function extract_galaxy(galaxyname, df)
    galaxy = select_galaxy(galaxyname, df)
    if galaxy === nothing
        return
    end

    ref = "MV2020"
    dist = dm_to_distance(galaxy.dmod)

    properties = OrderedDict(
        "ra" => parse_ra_to_degrees(galaxy.RA),
        "ra_em" => 0.0,
        "ra_ep" => 0.0,
        "ra_ref" => ref,
        "dec" => parse_dec_to_degrees(galaxy.Dec),
        "dec_em" => 0.0,
        "dec_ep" => 0.0,
        "dec_ref" => ref,
        "distance" => dist,
        "distance_em" => dist - dm_to_distance(galaxy.dmod - galaxy["dmod-"]),
        "distance_ep" => dm_to_distance(galaxy["dmod+"] + galaxy.dmod) - dist,
        "distance_modulus" => galaxy.dmod,
        "distance_modulus_em" => galaxy["dmod-"],
        "distance_modulus_ep" => galaxy["dmod+"],
        "distance_ref" => ref,
        "pmra" => galaxy.pmra,
        "pmra_em" => galaxy["epmra-"],
        "pmra_ep" => galaxy["epmra+"],
        "pmra_ref" => ref,
        "pmdec" => galaxy.pmdec,
        "pmdec_em" => galaxy["epmdec-"],
        "pmdec_ep" => galaxy["epmdec+"],
        "pmdec_ref" => ref,
        "radial_velocity" => galaxy["vh"],
        "radial_velocity_em" => galaxy["vh-"],
        "radial_velocity_ep" => galaxy["vh+"],
        "radial_velocity_ref" => ref,
        "sigma_v" => galaxy["sigma_s"],
        "sigma_v_em" => galaxy["sigma_s-"],
        "sigma_v_ep" => galaxy["sigma_s+"],
        "sigma_v_ref" => ref,
        "position_angle" => galaxy["PA"],
        "position_angle_em" => galaxy["PA-"],
        "position_angle_ep" => galaxy["PA+"],
        "postiion_angle_ref" => ref,
        "ellipticity" => galaxy["e=1-b/a"],
        "ellipticity_em" => galaxy["e-"],
        "ellipticity_ep" => galaxy["e+"],
        "ellipticity_ref" => ref,
        "Mv" => galaxy["Vmag"] - galaxy.dmod,
        "Mv_em" => galaxy["Vmag-"] + galaxy["dmod+"],
        "Mv_ep" => galaxy["Vmag+"] + galaxy["dmod-"],
        "mv" => galaxy["Vmag"],
        "mv_em" => galaxy["Vmag-"],
        "mv_ep" => galaxy["Vmag+"],
        "mv_ref" => ref,
        "rh" => galaxy["rh"],
        "rh_em" => galaxy["rh-"],
        "rh_ep" => galaxy["rh+"],
        "rh_ref" => ref,
        "metallicity" => galaxy["[Fe/H]"],
        "metallicity_em" => galaxy["feh-"],
        "metallicity_ep" => galaxy["feh+"],
        "metallicity_ref" => ref,
        "MV2020_references" => galaxy.References,
       )

    for key in keys(properties)
        if properties[key] == 999.0
            properties[key] = NaN
        end
    end

    open("properties_$(galaxyname).toml", "w") do io
        println("Writing properties for $galaxyname")
        TOML.print(io, properties)
    end
end


function to_mv_name(galaxyname)
    if startswith(galaxyname, "DES")
        return galaxyname
    end

    newname = titlecase(galaxyname)
    newname =  replace(newname, "_" => "")


    return newname

end

function parse_ra_to_degrees(ra_str::String)
    parts = split(ra_str, ':')
    h = parse(Float64, parts[1])
    m = length(parts) >= 2 ? parse(Float64, parts[2]) : 0.0
    s = length(parts) >= 3 ? parse(Float64, parts[3]) : 0.0
    total_hours = h + m/60 + s/3600
    return total_hours * 15.0  # Convert hours to degrees
end


function parse_dec_to_degrees(dec_str::String)
    parts = split(dec_str, ':')
    d = parse(Float64, parts[1])
    sign = d < 0 ? -1 : 1
    d_abs = abs(d)
    m = length(parts) >= 2 ? parse(Float64, parts[2]) : 0.0
    s = length(parts) >= 3 ? parse(Float64, parts[3]) : 0.0
    return sign * (d_abs + m/60 + s/3600)
end


@main
