import TOML
using Measurements
import LilGuys

function print_props(galaxyname)
    println(galaxyname)
    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors

    for prop in ["ra", "dec", "distance_modulus", "distance", "pmra", "pmdec", "radial_velocity", "sigma_v", "R_h", "ellipticity", "position_angle", "Mv"]
        ref = get(obs_props, prop*"_ref", "")
        println(prop, "\t", obs_props[prop], "\t", ref)
    end
end



print_props("sculptor")
println()
print_props("ursa_minor")
