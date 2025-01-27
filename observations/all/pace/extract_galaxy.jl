import TOML
using LilGuys
import CSV
using DataFrames
using OrderedCollections

function main(ARGS)
    df = CSV.read("pace_all.csv", DataFrame)

    for galaxyname in readdir("../../")
        try
            extract_galaxy(galaxyname, df)
        catch e
            println("Error processing $galaxyname: $e")
        end
    end
end


function extract_galaxy(galaxyname, df)
    matches = df[df.key .== to_pace_name(galaxyname), :]
    if size(matches, 1) != 1
        @warn ("Found $(size(matches, 1)) matches for $galaxyname")
        return 
    end
    galaxy = matches[1, :]

    properties = OrderedDict(
        "ra" => galaxy.ra,
        "ra_err" => 0.0,
        "ra_ref" => galaxy.ref_structure,
        "dec" => galaxy.dec,
        "dec_err" => 0.0,
        "dec_ref" => galaxy.ref_structure,
        "distance" => galaxy.distance,
        "distance_err" => max(galaxy.distance_em, galaxy.distance_ep),
        "distance_ref" => galaxy.ref_distance,
        "pmra" => galaxy.pmra,
        "pmra_err" => max(galaxy.pmra_em, galaxy.pmra_ep),
        "pmra_ref" => galaxy.ref_proper_motion,
        "pmdec" => galaxy.pmdec,
        "pmdec_err" => max(galaxy.pmdec_em, galaxy.pmdec_ep),
        "pmdec_ref" => galaxy.ref_proper_motion,
        "radial_velocity" => galaxy.vlos_systemic,
        "radial_velocity_err" => max(galaxy.vlos_systemic_em, galaxy.vlos_systemic_ep),
        "radial_velocity_ref" => galaxy.ref_vlos,
        "sigma_v" => galaxy.vlos_sigma,
        "sigma_v_err" => max(galaxy.vlos_sigma_em, galaxy.vlos_sigma_ep),
        "sigma_v_ref" => galaxy.ref_vlos,
        "PA" => galaxy.position_angle,
        "PA_err" => max(galaxy.position_angle_em, galaxy.position_angle_ep),
        "PA_ref" => galaxy.ref_structure,
        "ellipticity" => galaxy.ellipticity,
        "ellipticity_err" => max(galaxy.ellipticity_em, galaxy.ellipticity_ep),
        "ellipticity_ref" => galaxy.ref_structure,
        "Mv" => galaxy.M_V,
        "Mv_err" => max(galaxy.M_V_em, galaxy.M_V_ep),
       )
    for key in keys(properties)
        if ismissing(properties[key])
            properties[key] = NaN
        end
    end

    icrs = LilGuys.ICRS(
        ra = properties["ra"],
        dec = properties["dec"],
        distance = properties["distance"],
        pmra = properties["pmra"],
        pmdec = properties["pmdec"],
        radial_velocity = properties["radial_velocity"]
       )

    gsr = LilGuys.transform(LilGuys.GSR, icrs)

    properties["pmra_gsr"] = gsr.pmra
    properties["pmdec_gsr"] = gsr.pmdec
    properties["radial_velocity_gsr"] = gsr.radial_velocity
    properties["theta_pm_gsr"] = -39.64 # atand(pmra_gsr, pmdec_gsr)


    open("properties_$(galaxyname).toml", "w") do io
        TOML.print(io, properties)
    end
end


function to_pace_name(galaxyname)
    if startswith(galaxyname, "canesvenatici")
        galaxyname = replace(galaxyname, "canesvenatici" => "canes_venatici")
    end

    if galaxyname == "leo_t"
        return galaxyname
    end

    if match(r"\d+", galaxyname) === nothing
        return galaxyname * "_1"
    else
        return replace(galaxyname, r"(\D)(\d+)" => s"\1_\2")
    end
end



@main
