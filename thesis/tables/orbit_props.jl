using LilGuys
import TOML
using PyFITS
import StatsBase: quantile
import DataFrames: disallowmissing

function calc_r_J(halo, rho_peri)
    return 10 ^ LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10^lr) - 3*rho_peri, log10(halo.r_s))
end


function mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end


function get_orbit(galaxyname, haloname, orbitname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, orbitname, "orbital_properties.toml"))
    pot = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, orbitname, "simulation/agama_potential.ini"))
end


function get_samples(galaxyname, potname)
    props = read_fits( joinpath(ENV["DWARFS_ROOT"], "orbits", 
                                galaxyname, potname, "orbital_properties.fits"))
end


function print_orbit(args...)
    orbit_props = get_orbit(args...)

    for key in ["pericentre", "apocentre", "t_last_peri", "idx_peris", "distance_f"]
        val = orbit_props[key]
        if val isa Real
            val = round(orbit_props[key], sigdigits=4)
        end
        println(key, "\t", val)
    end
end

function print_orbit_samples(args...)
    samples = get_samples(args...)


    for key in ["pericentre", "apocentre", "distance"]
        print_quantity(key, samples[!, key])
    end
    print_quantity("time last peri", samples[!, "time_last_peri"] * T2GYR)
end

function print_quantity(key, x)
    x = disallowmissing(replace(x, missing => NaN))
    l, m, h = quantile(x, [0.16, 0.5, 0.84])
    println(key, "\t", m, "\t", l-m, "\t", h-m)
end

println("sculptor --- smallperi")
println("random samples")
print_orbit_samples("sculptor", "EP2020")
println("smallperi")
print_orbit("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi")

println("sculptor --- lmc flyby")
print_orbit_samples("sculptor", "vasiliev24_L3M11")
println("example")
print_orbit("sculptor", "1e7_new_v25_r2.5", "smallperilmc")

println("ursa minor --- smallperi")
print_orbit_samples("ursa_minor", "EP2020")
print_orbit("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5")

