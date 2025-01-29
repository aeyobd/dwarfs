using TOML
import LilGuys
using OrderedCollections
using Printf


function main(ARGS)
    for galaxyname in readdir("../..")
        MV20 = safe_load_toml("../MV2020/properties_$galaxyname.toml")
        pace = safe_load_toml("../pace/properties_$galaxyname.toml")

        newfile = "properties_$galaxyname.toml"

        properties = nothing
        println("trying $galaxyname")
        if MV20 !== nothing 
            properties = MV20
            println("\tloaded MV2020")
        elseif pace !== nothing
            properties = pace
            println("\tloaded pace")
        else
            println("\tskipping $galaxyname")
            continue
        end

        properties = combine_properties(properties, pace)
        update_derived!(properties)
        properties = reorder_keys(properties)

        open(newfile, "w") do f
            TOML.print(f, properties)
        end
    end
end

function safe_load_toml(filename)
    if isfile(filename)
        return TOML.parsefile(filename)
    else
        return nothing
    end
end

function is_better_measurement(MV20, pace, key; tol=1e-2)
    ep = MV20["$(key)_ep"]
    em = MV20["$(key)_em"]
    ep1 = pace["$(key)_ep"]
    em1 = pace["$(key)_em"]

    if (isnan(ep1) || isnan(em1)) && (!isnan(ep) && !isnan(em))
        return false
    end
    if (isnan(ep) || isnan(em)) && (!isnan(ep1) && !isnan(em1))
        return true
    end
    return (ep > ep1 * (1+tol)) && (em > em1 * (1+tol))
end


function replace_measurement!(best, pace, key; ref_key="$(key)_ref")
    for suffix in ["", "_em", "_ep"]
        best[key * suffix] = pace[key * suffix]
    end
    best[ref_key] = pace[ref_key]
end

function combine_properties(MV20, pace)
    if pace === nothing
        return MV20
    end

    best = copy(MV20)

    for key in ["sigma_v", "mv", "metallicity", "distance_modulus"]
        if is_better_measurement(best, pace, key)
            println("\tupdating $key")
            @printf("\t\tMV20: %8.4f + %8.4f - %8.4f\n", best[key], best["$(key)_ep"], best["$(key)_em"])
            @printf("\t\tpace: %8.4f + %8.4f - %8.4f\n", pace[key], pace["$(key)_ep"], pace["$(key)_em"])
            if key == "distance_modulus"
                replace_measurement!(best, pace, key, ref_key="distance_ref")
            else
                replace_measurement!(best, pace, key)
            end
        end
    end

    return best
end

function update_derived!(properties)
    dist = dm_to_distance(properties["distance_modulus"])
    properties["distance"] = dist
    properties["distance_em"] = dist - dm_to_distance(properties["distance_modulus"] - properties["distance_modulus_em"])
    properties["distance_ep"] = dm_to_distance(properties["distance_modulus"] + properties["distance_modulus_ep"]) - dist

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

end

"""
    dm_to_distance(dm) 

Converts a distance modulus (magnitudes) to a distance in kpc.
"""
function dm_to_distance(dm)
    return 10^(dm / 5 + 1 - 3)
end


function reorder_keys(best)
    new_properties = OrderedDict()
    
    for key in keys_order
        if haskey(best, key)
            new_properties[key] = best[key]
        else
            println("\tmissing key $key")
        end
    end

    if any(keys(best) .∉ [keys_order])
        println("\textra keys: ", setdiff(keys(best), keys_order))
    end

    return new_properties
end



const keys_order = ["ra", "ra_em", "ra_ep", "ra_ref", "dec", "dec_em", "dec_ep", "dec_ref", "distance", "distance_em", "distance_ep", "distance_modulus", "distance_modulus_em", "distance_modulus_ep", "distance_ref", "pmra", "pmra_em", "pmra_ep", "pmra_ref", "pmdec", "pmdec_em", "pmdec_ep", "pmdec_ref", "radial_velocity", "radial_velocity_em", "radial_velocity_ep", "radial_velocity_ref", "sigma_v", "sigma_v_em", "sigma_v_ep", "sigma_v_ref", "position_angle", "position_angle_em", "position_angle_ep", "postiion_angle_ref", "ellipticity", "ellipticity_em", "ellipticity_ep", "ellipticity_ref", "Mv", "Mv_em", "Mv_ep", "mv", "mv_em", "mv_ep", "mv_ref", "rh", "rh_em", "rh_ep", "rh_ref", "metallicity", "metallicity_em", "metallicity_ep", "metallicity_ref", "pmra_gsr", "pmdec_gsr", "radial_velocity_gsr", "theta_pm_gsr", "MV2020_references"]
@main
