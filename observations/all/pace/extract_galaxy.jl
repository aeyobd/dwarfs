import TOML
using LilGuys
import CSV
using DataFrames
using OrderedCollections

function main(ARGS)
    df = CSV.read("pace_all.csv", DataFrame)

    df2 = CSV.read("gc_mw_new.csv", DataFrame)

    df = vcat(df, df2)
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
        "ra_em" => 0.0,
        "ra_ep" => 0.0,
        "ra_ref" => galaxy.ref_structure,
        "dec" => galaxy.dec,
        "dec_em" => 0.0,
        "dec_ep" => 0.0,
        "dec_ref" => galaxy.ref_structure,
        "distance" => galaxy.distance,
        "distance_em" => galaxy.distance_em,
        "distance_ep" => galaxy.distance_ep,
        "distance_modulus" => galaxy.distance_modulus,
        "distance_modulus_em" => galaxy.distance_modulus_em,
        "distance_modulus_ep" => galaxy.distance_modulus_ep,
        "distance_ref" => galaxy.ref_distance,
        "pmra" => galaxy.pmra,
        "pmra_em" => galaxy.pmra_em,
        "pmra_ep" => galaxy.pmra_ep,
        "pmra_ref" => galaxy.ref_proper_motion,
        "pmdec" => galaxy.pmdec,
        "pmdec_em" => galaxy.pmdec_em,
        "pmdec_ep" => galaxy.pmdec_ep,
        "pmdec_ref" => galaxy.ref_proper_motion,
        "radial_velocity" => galaxy.vlos_systemic,
        "radial_velocity_em" => galaxy.vlos_systemic_em,
        "radial_velocity_ep" => galaxy.vlos_systemic_ep,
        "radial_velocity_ref" => galaxy.ref_vlos,
        "sigma_v" => galaxy.vlos_sigma,
        "sigma_v_em" => galaxy.vlos_sigma_em,
        "sigma_v_ep" => galaxy.vlos_sigma_ep,
        "sigma_v_ref" => galaxy.ref_vlos,
        "position_angle" => galaxy.position_angle,
        "position_angle_em" => galaxy.position_angle_em,
        "position_angle_ep" => galaxy.position_angle_ep,
        "postiion_angle_ref" => galaxy.ref_structure,
        "ellipticity" => galaxy.ellipticity,
        "ellipticity_em" => galaxy.ellipticity_em,
        "ellipticity_ep" => galaxy.ellipticity_ep,
        "ellipticity_ref" => galaxy.ref_structure,
        "Mv" => galaxy.M_V,
        "Mv_em" => galaxy.M_V_em,
        "Mv_ep" => galaxy.M_V_ep,
        "mv" => galaxy.apparent_magnitude_v,
        "mv_em" => galaxy.apparent_magnitude_v_em,
        "mv_ep" => galaxy.apparent_magnitude_v_ep,
        "mv_ref" => galaxy.ref_m_v,
        "rh" => galaxy.rhalf,
        "rh_em" => galaxy.rhalf_em,
        "rh_ep" => galaxy.rhalf_ep,
        "rh_ref" => galaxy.ref_structure,
        "metallicity" => galaxy.metallicity_spectroscopic,
        "metallicity_em" => galaxy.metallicity_spectroscopic_em,
        "metallicity_ep" => galaxy.metallicity_spectroscopic_ep,
        "metallicity_ref" => galaxy.ref_metallicity_spectroscopic,
       )

    for key in keys(properties)
        if ismissing(properties[key])
            properties[key] = NaN
        end
    end

    open("properties_$(galaxyname).toml", "w") do io
        TOML.print(io, properties)
    end
end


function to_pace_name(galaxyname)
    if startswith(galaxyname, "canesvenatici")
        galaxyname = replace(galaxyname, "canesvenatici" => "canes_venatici")
    end

    if galaxyname ∈ ["leo_t", "smc", "lmc"]
        return galaxyname
    end

    if match(r"\d+", galaxyname) === nothing
        return galaxyname * "_1"
    else
        return replace(galaxyname, r"(\D)(\d+)" => s"\1_\2")
    end
end



@main
