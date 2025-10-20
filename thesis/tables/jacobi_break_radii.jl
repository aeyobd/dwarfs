using LilGuys

using OrderedCollections
using DataFrames
using Printf
using PyFITS
import TOML
import Agama

quantile = LilGuys.quantile

include("utils.jl")

function calc_r_J(halo, rho_peri)
    return 10 ^ LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10^lr) - 3*rho_peri, log10(halo.r_s))
end


function print_radius(label, r, df)
    r_arcmin = LilGuys.kpc2arcmin.(r, df.distance)
    @printf "%s kpc   \t%0.2f ± %0.2f\n" label LilGuys.mean(r) LilGuys.std(r)
    @printf "%s arcmin\t%0.2f ± %0.2f\n" label LilGuys.mean(r_arcmin) LilGuys.std(r_arcmin)
end

function mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end


replace_missings(x) = ifelse.(ismissing.(x), NaN, x)

function analyze_lmc(galaxyname, modelname, halo; units=Agama.VASILIEV_UNITS)

    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))


	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)

    df = read_fits(joinpath(modeldir, "orbital_properties.fits"))

    pot = Agama.Potential(file = joinpath(modeldir, "potential_mw_init.ini"))
    pot_lmc = Agama.Potential(file = joinpath(modeldir, "potential_lmc_init.ini"))

    rho_peri = mean_density(pot, df.pericentre, units)
    rho_peri_lmc = mean_density(pot_lmc, df.pericentre_lmc, units)

    r_J = calc_r_J.(halo, rho_peri)
    r_J_lmc = calc_r_J.(halo, rho_peri_lmc)
    print_radius("r_J", r_J, df)
    print_radius("r_J_lmc", r_J_lmc, df)


    σv = obs_props["sigma_v"] .+ randn(size(df, 1)) * obs_props["sigma_v_err"]

	
    r_break = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri)
    r_break_lmc = LilGuys.break_radius.(σv / V2KMS, replace_missings(df.time_last_peri_lmc))
	
    print_radius("r_break", r_break, df)
    print_radius("r_break_lmc", r_break_lmc, df)
    print_radius("r_break_lmc", r_break_lmc, df)

	print_quantity("pericentre", df.pericentre)
	print_quantity("pericentre lmc", df.pericentre_lmc)
	print_quantity("apocentre", df.apocentre)
	print_quantity("apocentre_lmc", df.apocentre_lmc)
	print_quantity("time last", df.time_last_peri * T2GYR)
	print_quantity("time last lmc", replace_missings(df.time_last_peri_lmc) * T2GYR)
end

function analyze_mw(galaxyname, modelname, halo; units=Agama.AgamaUnits())
    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))


	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)

    df = read_fits(joinpath(modeldir, "orbital_properties.fits"))

    pot = Agama.Potential(file = joinpath(modeldir, "agama_potential.ini"))

    rho_peri = mean_density(pot, df.pericentre, units)

    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J", r_J, df)


    σv = obs_props["sigma_v"] .+ randn(size(df, 1)) * obs_props["sigma_v_err"]

	t_last_peri =  df.time_last_peri
	t_last_peri[DataFrames.ismissing.(t_last_peri)] .== -1
	
    r_break = LilGuys.break_radius.(σv / V2KMS, t_last_peri)
    print_radius("r_break", r_break, df)
	print_quantity("pericentre", df.pericentre)
	print_quantity("apocentre", df.apocentre)
	print_quantity("time last", df.time_last_peri * T2GYR)
	#print_quantity("periods", floor.(10 ./ (df.period * T2GYR)))
end

halo_scl = NFW(r_circ_max=3.2, v_circ_max=31/V2KMS)
analyze_mw(modelnames["scl_orbits"]..., halo_scl)
halo_scl = NFW(r_circ_max=2.5, v_circ_max=25/V2KMS)
analyze_lmc(modelnames["scl_lmc_orbits"]..., halo_scl)


halo_umi = NFW(r_circ_max=4, v_circ_max=38/V2KMS)
analyze_mw(modelnames["umi_orbits"]..., halo_scl)
