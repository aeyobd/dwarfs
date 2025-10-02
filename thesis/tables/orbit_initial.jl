using LilGuys
import TOML


function get_orbit(galaxyname, orbitname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, orbitname * ".toml"))
end

function get_orbit(galaxyname, haloname, orbitname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, orbitname, "orbital_properties.toml"))
end

function get_orbit_i(galaxyname, haloname, orbitname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, orbitname, "simulation/orbit.toml"))
end


function print_orbit(args...)
    orbit_props = get_orbit(args...)
    orbit_props_i = get_orbit_i(args...)

    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity", "x_i", "y_i", "z_i"]
        println(key, "\t", round(orbit_props_i[key], sigdigits=4))
    end

    for key in ["v_x_i", "v_y_i", "v_z_i"]
        println(key, "\t", round(orbit_props_i[key] * V2KMS, digits=2))
    end

    println("t_i (nbody)", "\t", orbit_props["t_f_gyr"])
    println("t_i (point)", "\t", orbit_props_i["t_i"] * T2GYR)
end

function print_point_orbit(args...)
    orbit_props_i = get_orbit(args...)

    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity", "x_i", "y_i", "z_i"]
        println(key, "\t", round(orbit_props_i[key], sigdigits=4))
    end

    for key in ["v_x_i", "v_y_i", "v_z_i"]
        println(key, "\t", round(orbit_props_i[key] * V2KMS, digits=2))
    end
end


println("sculptor --- smallperi")
print_orbit("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi")

println("sculptor --- lmc flyby")
print_orbit("sculptor", "1e7_new_v25_r2.5", "smallperilmc")

println("ursa minor --- smallperi")
print_orbit("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5")


println("sculptor point")
#print_point_orbit("sculptor", "vasiliev24_L3M11_2x_special_cases/orbit_smallperilmc")
#print_point_orbit("sculptor", "EP2020_special_cases/orbit_smallperi")
#print_point_orbit("ursa_minor", "EP2020_special_cases/orbit_smallperi")
