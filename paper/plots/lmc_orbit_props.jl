using LilGuys
import TOML
using PyFITS
import StatsBase: quantile, median
import DataFrames: disallowmissing
using Printf
import Agama

include("table_utils.jl")

units = Agama.VASILIEV_UNITS

function calc_r_J(halo, rho_peri)
    return 10 ^ LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10^lr) - 3*rho_peri, log10(halo.r_s))
end


function print_radius(label, r, distance)
    r_arcmin = LilGuys.kpc2arcmin.(r, distance)
    print_quantity(label * " / kpc", r)
    print_quantity(label * " / arcmin", r_arcmin)
end

function LilGuys.mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end

function get_orbit(galaxyname, modelname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "orbital_properties.toml"))
    orbit_props
end

function get_orbit_lmc(galaxyname, modelname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "orbital_properties_lmc.toml"))
    orbit_props
end

function get_samples(galaxyname, potname)
    props = read_fits( joinpath(ENV["DWARFS_ROOT"], "orbits", 
                                galaxyname, potname, "orbital_properties.fits"))
end


function print_orbit(galaxyname, modelname, starsname, halo)
    orbit_props = get_orbit(galaxyname, modelname)
    orbit_props_lmc = get_orbit_lmc(galaxyname, modelname)

    for key in ["pericentre", "apocentre", "t_last_peri", "idx_peris", "distance_f"]
        val = orbit_props[key]
        if val isa Real
            val = round(orbit_props[key], sigdigits=4)
        end
        println(key, "\t", val)
    end

    for key in ["pericentre", "t_last_peri"]
        println(key, "\t", orbit_props_lmc[key])
    end


	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)
    pot = Agama.Potential(file = joinpath(modeldir, "potential_mw_init.ini"))
    rho_peri = LilGuys.mean_density(pot, orbit_props["pericentre"], units)
    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J", r_J, orbit_props["distance_f"])


    pot = Agama.Potential(file = joinpath(modeldir, "potential_lmc_init.ini"))
    rho_peri = LilGuys.mean_density(pot, orbit_props_lmc["pericentre"], units)
    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J", r_J, orbit_props["distance_f"])


    obs_props = get_obs_props(galaxyname)
    σv = obs_props["sigma_v"] 
    r_break = LilGuys.break_radius.(orbit_props["t_last_peri"]/T2GYR, σv/V2KMS)
    print_radius("r_break", r_break, orbit_props["distance_f"])
    r_break = LilGuys.break_radius.(orbit_props_lmc["t_last_peri"]/T2GYR, σv/V2KMS)
    print_radius("r_break", r_break, orbit_props["distance_f"])
end


function print_orbit_samples(galaxyname, modelname, halo)
    samples = get_samples(galaxyname, modelname)

    for key in ["pericentre", "apocentre", "time_last_peri", "distance", "n_peris", "period_apo", "pericentre_lmc", "time_last_peri_lmc"]
        if key ∈ ["time_last_peri", "time_last_peri_lmc"]
            print_quantity(key, samples[!, key] * T2GYR)
        else
            print_quantity(key, samples[!, key] * 1.0)
        end
    end


	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)
    pot = Agama.Potential(file = joinpath(modeldir, "potential_mw_init.ini"))
    rho_peri = LilGuys.mean_density(pot, samples.pericentre, units)
    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J", r_J, samples.distance)

    pot = Agama.Potential(file = joinpath(modeldir, "potential_lmc_init.ini"))
    rho_peri = LilGuys.mean_density(pot, samples.pericentre_lmc, units)
    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J_lmc", r_J, samples.distance)

    obs_props = get_obs_props(galaxyname)
    σv = obs_props["sigma_v"] .+ obs_props["sigma_v_err"] * randn(size(samples, 1))
    r_break = LilGuys.break_radius.(samples.time_last_peri, σv/V2KMS)
    print_radius("r_break", r_break, samples.distance)
    r_break = LilGuys.break_radius.(samples.time_last_peri_lmc, σv/V2KMS)
    print_radius("r_break_lmc", r_break, samples.distance)
end



println("sculptor --- lmc orbits")
halo = NFW(r_circ_max=2.5, v_circ_max=25/V2KMS)
print_orbit_samples(modelnames["scl_lmc_orbits"]..., halo)
println("sculptor --- lmc flyby")
print_orbit(modelnames["scl_lmc"]..., halo)


println("umi --- lmc orbits")
halo = NFW(r_circ_max=4, v_circ_max=38/V2KMS)
print_orbit_samples(modelnames["umi_lmc_orbits"]..., halo)
