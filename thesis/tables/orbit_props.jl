using LilGuys
import TOML
using PyFITS
import StatsBase: quantile
import DataFrames: disallowmissing
using Printf
import Agama

include("utils.jl")


function calc_r_J(halo, rho_peri)
    return 10 ^ LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10^lr) - 3*rho_peri, log10(halo.r_s))
end


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
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "orbital_properties.toml"))
    orbit_props
end


function get_samples(galaxyname, potname)
    props = read_fits( joinpath(ENV["DWARFS_ROOT"], "orbits", 
                                galaxyname, potname, "orbital_properties.fits"))
end


function print_orbit(galaxyname, modelname, starsname, halo)
    orbit_props = get_orbit(galaxyname, modelname)

    for key in ["pericentre", "apocentre", "t_last_peri", "idx_peris", "distance_f"]
        val = orbit_props[key]
        if val isa Real
            val = round(orbit_props[key], sigdigits=4)
        end
        println(key, "\t", val)
    end


	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)
    pot = Agama.Potential(file = joinpath(modeldir, "simulation/agama_potential.ini"))
    rho_peri = mean_density(pot, orbit_props["pericentre"], Agama.AgamaUnits())
    r_J = calc_r_J(halo, rho_peri)
    println("r_J", "\t", r_J)
    println("r_J_arcmin\t", LilGuys.kpc2arcmin(r_J, orbit_props["distance_f"]))
end


function print_orbit_samples(galaxyname, modelname, halo)
    samples = get_samples(galaxyname, modelname)

    for key in ["pericentre", "apocentre", "n_peris", "period_apo", "time_last_peri", "distance"]
        if key ∈ ["time_last_peri", "period_apo"]
            print_quantity(key, samples[!, key] * T2GYR)
        else
            print_quantity(key, samples[!, key])
        end
    end


	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)
    pot = Agama.Potential(file = joinpath(modeldir, "agama_potential.ini"))
    units = Agama.AgamaUnits()
    rho_peri = mean_density(pot, samples.pericentre, units)
    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J", r_J, samples.distance)
end



println("sculptor --- samples")
println("random samples")
scl_halo = LilGuys.NFW(v_circ_max=31/V2KMS, r_circ_max=3.2)
print_orbit_samples(modelnames["scl_orbits"]..., scl_halo)
println("sculptor --- smallperi")
print_orbit(modelnames["scl_smallperi"]..., scl_halo)

umi_halo = LilGuys.NFW(v_circ_max=38/V2KMS, r_circ_max=4)
println("ursa minor --- smallperi")
print_orbit_samples(modelnames["umi_orbits"]..., umi_halo)
println("ursa minor --- smallperi")
print_orbit(modelnames["umi_smallperi"]..., umi_halo)

