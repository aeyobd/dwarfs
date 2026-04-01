### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 791d6268-2c81-11f1-9366-e7a4214e567f
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie

	import CSV
	import DataFrames: DataFrame
	import TOML
end

# ╔═╡ ba2b61b8-ca92-49d0-9ea6-922a5e450c6a
phot_system = :delve

# ╔═╡ 5cdee12d-4899-41d7-aac0-d7538e5ea6de
samplesname = Dict(
	:delve => "samples.delve_matched_filter_6deg.mcmc_plummer.csv",
	:gaia => "samples.G_20.0.mcmc_sersic.csv",
)[phot_system]

# ╔═╡ 39cd58d0-6d8e-40a4-8492-c30fe1782531
mag_limit = Dict(
	:gaia => 20,
	:delve => 23
)[phot_system]

# ╔═╡ bebc5ba6-870f-473a-a05f-471f854b1ef8
import StatsBase: median

# ╔═╡ f6e51d53-506e-4c17-8c6d-b044a9bef0cf
obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")

# ╔═╡ e4fdfb9b-73a8-410b-b86f-cdb1d71019a1
obs_props = TOML.parsefile(joinpath(obs_dir, "observed_properties.toml"))

# ╔═╡ 0f977fa3-24ff-4f40-9407-0f783ea9ef26
dm_0 = obs_props["distance_modulus"]

# ╔═╡ 4de8788f-1c4b-4fa5-ab70-02c1fd1c5c4a
dm_err = obs_props["distance_modulus_em"]

# ╔═╡ d26b9e4b-d68f-4b72-a080-1d5c66248e19
module Utils
	obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")
	include(joinpath(obs_dir, "delve_utils.jl"))
end

# ╔═╡ a4b3e72d-607f-43b0-869e-096646046eb2
md"""
# Utiliies
"""

# ╔═╡ 15e5f243-1f14-479f-978d-bb09b9e3a72b
function get_all_isochrones_gaia(M_H::Real, age::Real=12)
    iso_columns = string.(split("Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	mbolmag  Gmag    G_BPbrmag  G_BPftmag  G_RPmag", 
	r"\s+"))

    all_isochrones = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/padova/isochrone.gaiadr2.weiler2018.$(age)Gyrs.dat"),DataFrame,
					 comment="#", ignorerepeated=true, delim = ' ', header=iso_columns)

	all_isochrones
end


# ╔═╡ c725772e-e4c3-4fad-95ac-d3ed601012f1
function get_all_isochrones_decam( age::Real=12;)
	ISO_HEADER = "Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label   McoreTP C_O  period0  period1  period2  period3  period4  pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	mbolmag  umag    gmag    rmag    imag    zmag    Ymag"
	
    iso_columns = string.(split(ISO_HEADER, r"\s+"))

    all_isochrones = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/padova/isocahrone.decam.$(age)Gyrs.dat"),DataFrame,
					 comment="#", ignorerepeated=true, delim = ' ', header=iso_columns)

	

	all_isochrones
end

# ╔═╡ df906ecc-29a2-4c4b-a93d-265e1e381061
all_isochrones_gaia = Dict(
	age => get_all_isochrones_gaia(age) for age in [8, 9, 10, 11, 12, 13]
)

# ╔═╡ fbb7df70-3f62-48f0-a3b3-c64b356cc90d
function get_isochrone(M_H, age, system=:gaia; stage_max=5)
	if system == :gaia
		isos = all_isochrones_gaia[age]
	end
	M_Hs = unique(isos.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]

	filt = isapprox.(isos.MH, M_H_adopted)
	filt .&= isos.label .<= stage_max # only keep through HB

	return isos[filt, :]
end

# ╔═╡ 310f64f4-11ba-4400-ad55-e23f17954161
function get_isochrone_interp(iso; cols=["gmag", "rma", "Gmag", "G_BPftmag", "G_RPmag"])
	interps = Dict(col => LilGuys.lerp(iso.Mini, iso[!, col]) for col in cols if col in names(iso))
end

# ╔═╡ e20663bd-028a-47ac-8efb-5f893321daa5
md"""
## Samplers
"""

# ╔═╡ 121b92b6-4bb7-4800-94c7-f57684b9f184
function sample_kroupa(N; mass_max, mass_min=0.08)
	masses = LinRange(mass_min, mass_max, 10_000)
	imf =  Utils.kroupa_imf.(masses)
	N_cdf = cumsum(imf)
	x = N_cdf .- N_cdf[1]
	x ./= x[end]
	mass_func =  LilGuys.lerp(x, masses)

	return mass_func.(rand(N))

end

# ╔═╡ 41aa41c8-8bab-42c9-b548-39aab368e8aa
function create_kroupa_sampler(mass_max; mass_min=0.08)
	masses = LinRange(mass_min, mass_max, 10_000)
	imf =  Utils.kroupa_imf.(masses)
	N_cdf = cumsum(imf)
	x = N_cdf .- N_cdf[1]
	x ./= x[end]
	mass_func =  LilGuys.lerp(x, masses)

	function sampler()
		u = rand()
		return mass_func(u)
	end

end

# ╔═╡ 5b2975d7-5d15-4ffe-b328-ba7f64549a25
function find_tot_number_stars(N_obs; mag_limit, M_H=-2.1, age=12, N_max=300_000, distance_modulus=18.34)
	iso = get_isochrone(M_H, age)
	iso_interp = get_isochrone_interp(iso)
	M_max = maximum(iso.Mini)
	sampler = create_kroupa_sampler(M_max, mass_min=0.08)

	count = 0
	mass_tot = 0
	L_tot = 0
	N_tot = 0
	gcol = "Gmag"
	
	for i in 1:N_max
		mass = sampler()
		mag = iso_interp[gcol](mass) + distance_modulus

		N_tot += 1
		L_tot += 10^(-0.4*mag)
		mass_tot += mass
		count += mag < mag_limit
		if count >= N_obs
			break
		end
	end

	if N_tot == N_max
		@warn "failed to converge"
		@info count / N_obs
		return NaN, NaN, NaN
	end

	return N_tot, mass_tot, -2.5*log10(L_tot) - distance_modulus
end

# ╔═╡ 4574c42f-e5c4-466a-a209-6fb30eab8cb4
function sample_tot_stars(obs_samples, N_samples; kwargs...)
	N_tots = Vector{Int}(undef, N_samples)
	MVs = Vector{Float64}(undef, N_samples)
	M_tots = Vector{Float64}(undef, N_samples)
	
	for i in 1:N_samples
		distance_modulus = dm_0 + dm_err*randn()
		M_H = rand(-1.8:-0.01:-2.1)
		age = rand(9:12)

		N_tots[i], M_tots[i], MVs[i] = find_tot_number_stars(obs_samples[i]; 
			distance_modulus=distance_modulus,
			M_H = M_H,
			age = age,
			kwargs...)
		
		if (i % 100) == 0
			@info "completed $i / $N_samples"
		end
	end

	return N_tots, M_tots, MVs
end

# ╔═╡ ce015a4e-3c9c-405f-afa6-59324216c591
md"""
# The main drag
"""

# ╔═╡ 37f3a337-05a3-401d-91ed-c5c94eba899a
mcmc_samples = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/mcmc", samplesname), DataFrame)

# ╔═╡ 4a6a1c7a-a4b2-468f-9db8-a814cdeee666
N_obs_samples = mcmc_samples.N_memb

# ╔═╡ 0cbecf30-0449-4d03-9874-105c75bd7830
N_tots, M_tots, M_Vs = sample_tot_stars(N_obs_samples, 300, mag_limit=mag_limit)

# ╔═╡ f1af266c-86d9-4649-bd6d-cd5c2b695892
md"""
# Analysis and saving
"""

# ╔═╡ 6ac311e3-268c-4e2d-95ca-77a643ce1bbb
hist(N_obs_samples)

# ╔═╡ f66db8b0-1aef-4419-bbcf-84c9ef11d477
hist(N_tots / 1e3)

# ╔═╡ bbd1541c-4fdd-40df-b03e-a1a60730336b
hist(M_tots / 1e3)

# ╔═╡ 443e2da9-5620-4139-bf58-c217d9708ef3
hist(M_Vs )

# ╔═╡ ad6ce2cc-292a-45c1-b885-9d6927971c35
find_tot_number_stars(100, mag_limit=21, N_max=1e5)

# ╔═╡ 87653c99-3492-42b2-87c6-ca693d307493
median(M_Vs), LilGuys.quantile(M_Vs, [0.16, 0.84]) .- median(M_Vs)

# ╔═╡ e602c637-f6a0-4e90-aeb5-79ef55c6dc39
median(M_tots), LilGuys.quantile(M_tots, [0.16, 0.84]) .- median(M_tots)

# ╔═╡ 055f1d47-0436-4267-a2fd-139c358ecd84
median(N_tots), LilGuys.quantile(N_tots, [0.16, 0.84]) .- median(N_tots)

# ╔═╡ 06a65372-9145-49e8-a4e2-43b552ab96cf
median(N_obs_samples), LilGuys.quantile(N_obs_samples, [0.16, 0.84]) .- median(N_obs_samples)

# ╔═╡ 5e2c55c0-7ac2-430b-84fb-cac861132fde
md"""
# Sanity checks
"""

# ╔═╡ ca8329d8-66f0-432b-a0a2-6d856c2e5fac
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale=log10,
			 yscale=log10, 
			 xticks=[0.01, 0.1, 1.0, 2.0],
			 yticks=[0.01, 0.1, 1, 10])
	
	masses = sample_kroupa(10_000, mass_max=2, mass_min=0.01)
		
	stephist!(masses, bins=logrange(0.01, 2, 100), normalization=:pdf)

	x = LinRange(0.01, 2, 1000)
	A = LilGuys.integrate(Utils.kroupa_imf, 0.01, 2)
	
	lines!(x, Utils.kroupa_imf.(x) / A )

	fig

end

# ╔═╡ 01cb561b-0e1f-4443-bc2c-8943d4c944cc
iso = get_isochrone(-2.1, 12)

# ╔═╡ 1b003181-ecc0-4cb7-923d-eccab5ddac10
iso[!, :Gmag]

# ╔═╡ e24c999a-ca90-45af-818a-d936ed293cb0
let
	iso = get_isochrone(-2.1, 12)

	f =	lines(iso.G_BPftmag - iso.G_RPmag, iso.Gmag .+ obs_props["distance_modulus"])
	iso_interp = get_isochrone_interp(iso)
	@info maximum(iso.Mass)

	x = LinRange(0.08, 0.92, 10000)
	lines!(iso_interp["G_BPftmag"].(x) - iso_interp["G_RPmag"].(x), iso_interp["Gmag"].(x) .+ obs_props["distance_modulus"])

	f
end

# ╔═╡ Cell order:
# ╠═ba2b61b8-ca92-49d0-9ea6-922a5e450c6a
# ╠═5cdee12d-4899-41d7-aac0-d7538e5ea6de
# ╠═39cd58d0-6d8e-40a4-8492-c30fe1782531
# ╠═791d6268-2c81-11f1-9366-e7a4214e567f
# ╠═bebc5ba6-870f-473a-a05f-471f854b1ef8
# ╠═f6e51d53-506e-4c17-8c6d-b044a9bef0cf
# ╠═e4fdfb9b-73a8-410b-b86f-cdb1d71019a1
# ╠═0f977fa3-24ff-4f40-9407-0f783ea9ef26
# ╠═4de8788f-1c4b-4fa5-ab70-02c1fd1c5c4a
# ╠═d26b9e4b-d68f-4b72-a080-1d5c66248e19
# ╠═a4b3e72d-607f-43b0-869e-096646046eb2
# ╠═15e5f243-1f14-479f-978d-bb09b9e3a72b
# ╠═c725772e-e4c3-4fad-95ac-d3ed601012f1
# ╠═df906ecc-29a2-4c4b-a93d-265e1e381061
# ╠═fbb7df70-3f62-48f0-a3b3-c64b356cc90d
# ╠═310f64f4-11ba-4400-ad55-e23f17954161
# ╠═e20663bd-028a-47ac-8efb-5f893321daa5
# ╠═121b92b6-4bb7-4800-94c7-f57684b9f184
# ╠═41aa41c8-8bab-42c9-b548-39aab368e8aa
# ╠═5b2975d7-5d15-4ffe-b328-ba7f64549a25
# ╠═4574c42f-e5c4-466a-a209-6fb30eab8cb4
# ╠═ce015a4e-3c9c-405f-afa6-59324216c591
# ╠═37f3a337-05a3-401d-91ed-c5c94eba899a
# ╠═4a6a1c7a-a4b2-468f-9db8-a814cdeee666
# ╠═1b003181-ecc0-4cb7-923d-eccab5ddac10
# ╠═0cbecf30-0449-4d03-9874-105c75bd7830
# ╠═f1af266c-86d9-4649-bd6d-cd5c2b695892
# ╠═6ac311e3-268c-4e2d-95ca-77a643ce1bbb
# ╠═f66db8b0-1aef-4419-bbcf-84c9ef11d477
# ╠═bbd1541c-4fdd-40df-b03e-a1a60730336b
# ╠═443e2da9-5620-4139-bf58-c217d9708ef3
# ╠═ad6ce2cc-292a-45c1-b885-9d6927971c35
# ╠═87653c99-3492-42b2-87c6-ca693d307493
# ╠═e602c637-f6a0-4e90-aeb5-79ef55c6dc39
# ╠═055f1d47-0436-4267-a2fd-139c358ecd84
# ╠═06a65372-9145-49e8-a4e2-43b552ab96cf
# ╟─5e2c55c0-7ac2-430b-84fb-cac861132fde
# ╠═ca8329d8-66f0-432b-a0a2-6d856c2e5fac
# ╠═01cb561b-0e1f-4443-bc2c-8943d4c944cc
# ╠═e24c999a-ca90-45af-818a-d936ed293cb0
