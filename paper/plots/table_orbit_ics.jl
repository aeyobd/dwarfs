using LilGuys
import TOML
using PyFITS
import StatsBase: quantile
import DataFrames: disallowmissing
using Printf
import Agama

include("table_utils.jl")



function print_radius(label, r, distance)
    r_arcmin = LilGuys.kpc2arcmin.(r, distance)
    print_quantity(label * " / kpc", r)
    print_quantity(label * " / arcmin", r_arcmin)
end

function mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end

function get_orbit(galaxyname, modelname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "simulations", galaxyname, modelname, "orbit.toml"))
    orbit_props
end



function print_orbit(galaxyname, modelname, starsname)
    orbit_props = get_orbit(galaxyname, modelname)

    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
        val = orbit_props[key]
        if val isa Real
            val = round(orbit_props[key], sigdigits=5)
        end
        println(key, "\t", val)
    end


    println()
end




println("sculptor smallperi ")
print_orbit(modelnames["scl_smallperi"]...,)


println("scl lmc ")
print_orbit(modelnames["scl_lmc"]...,)

println("umi.5 ")
print_orbit(modelnames["umi_smallperi"]...,)

println("umi.1")
print_orbit("ursa_minor", "1e5_v38_r4.0/orbit_smallperi.1", "")
