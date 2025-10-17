using LilGuys
import TOML
include("utils.jl")


function get_orbit(galaxyname, orbitname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, orbitname * ".toml"))
end

function get_orbit(galaxyname, modelname, starsname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "orbital_properties.toml"))
end

function get_orbit_i(galaxyname, modelname, starsname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "simulation/orbit.toml"))
end


function print_orbit(args...)
    orbit_props = get_orbit(args...)
    orbit_props_i = get_orbit_i(args...)

    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
        println(key, "\t", round(orbit_props_i[key], sigdigits=4))
    end

    println("t_i (nbody)", "\t", orbit_props["t_f_gyr"])
    t_i = get(orbit_props_i, "t_i", NaN)
    println("t_i (point)", "\t", t_i * T2GYR)

    for key in ["x_i", "y_i", "z_i"]
        println(key, "\t", round(orbit_props_i[key], sigdigits=4))
    end


    if radii([orbit_props_i["v_x_i"], orbit_props_i["v_y_i"], orbit_props_i["v_z_i"]]) > 10
        v_scale = 1
    else
        v_scale = V2KMS
    end

    for key in ["v_x_i", "v_y_i", "v_z_i"]
        println(key, "\t", round(orbit_props_i[key] * v_scale, digits=2))
    end

    println()
end

function print_point_orbit(args...)
    orbit_props_i = get_orbit(args...)

    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity", "x_i", "y_i", "z_i"]
        println(key, "\t", round(orbit_props_i[key], sigdigits=4))
    end

    for key in ["v_x_i", "v_y_i", "v_z_i"]
        println(key, "\t", round(orbit_props_i[key] * V2KMS, digits=2))
    end
    println()
end


println("sculptor --- smallperi")
print_orbit(modelnames["scl_smallperi"]...)

println("sculptor --- lmc flyby")
print_orbit(modelnames["scl_lmc"]...)

println("ursa minor --- smallperi")
print_orbit(modelnames["umi_smallperi"]...)


println("sculptor -- mw impact")
print_orbit(modelnames["mw_impact"]...)
#print_point_orbit("sculptor", "vasiliev24_L3M11_9Gyr_special_cases/orbit_smallperilmc")

