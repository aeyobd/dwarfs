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

	import StatsBase: median
end

# ╔═╡ ba2b61b8-ca92-49d0-9ea6-922a5e450c6a
sample = :gaia_sersic

# ╔═╡ 6eb7f393-3098-4e86-bab8-f17ae6f5c8a5
phot_system = if sample == :delve
	:delve
else
	:gaia
end

# ╔═╡ 5cdee12d-4899-41d7-aac0-d7538e5ea6de
samplesname = Dict(
	:delve => "samples.delve_matched_filter_6deg.g22.mcmc_ell.csv",
	:gaia => "../final_derivations/mcmc/samples.j24_1c.mcmc_ell.csv",
	:gaia_exp => "../final_derivations/mcmc/samples.j24_1c.mcmc_exp_ell.csv",
	:gaia_sersic => "../final_derivations/mcmc/samples.j24_1c.mcmc_sersic_ell.csv",
	:gaia_bright => "samples.best_sample.G_20.simple.L_0.01.mcmc_sersic_ell.csv",
)[sample]

# ╔═╡ 03c38e6b-f404-4337-99be-e100ec4f6d07
if phot_system == :gaia
	gcol = "Gmag"
else
	gcol = "gmag"
end

# ╔═╡ 39cd58d0-6d8e-40a4-8492-c30fe1782531
mag_limit = Dict(
	:gaia => 21,
	:gaia_exp => 21,
	:gaia_sersic => 21,
	:delve => 22,
	:gaia_bright => 20,
)[sample]

# ╔═╡ d26b9e4b-d68f-4b72-a080-1d5c66248e19
module Utils # only used for obs_dir and kroupa_imf
	obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")
	include(joinpath(obs_dir, "delve_utils.jl"))
end

# ╔═╡ 68adcf15-2683-4456-bf1e-0a1cd6907676
CairoMakie.activate!(type=:png)

# ╔═╡ 78e8db76-054e-46ae-9ef0-5837b4811a9b
read_mist_file = Utils.read_mist_file

# ╔═╡ f6e51d53-506e-4c17-8c6d-b044a9bef0cf
obs_dir = Utils.obs_dir

# ╔═╡ f0894db7-2a6c-49c3-a97e-6a1fc04856be
kroupa_imf = Utils.kroupa_imf

# ╔═╡ e4fdfb9b-73a8-410b-b86f-cdb1d71019a1
obs_props = TOML.parsefile(joinpath(obs_dir, "observed_properties.toml"))

# ╔═╡ 0f977fa3-24ff-4f40-9407-0f783ea9ef26
dm_0 = obs_props["distance_modulus"]

# ╔═╡ 4de8788f-1c4b-4fa5-ab70-02c1fd1c5c4a
dm_err = obs_props["distance_modulus_em"]

# ╔═╡ 20a73a2f-7326-4cee-8489-8b0b3a217828
@assert obs_props["distance_modulus_em"] == obs_props["distance_modulus_ep"]

# ╔═╡ a4b3e72d-607f-43b0-869e-096646046eb2
md"""
# Utiliies
"""

# ╔═╡ 64f34c1f-6dbd-4cdf-a575-183cf0ede2ea
read_mist_file(joinpath(obs_dir, "../padova", "isochrone.ubvri.12Gyrs.dat"))

# ╔═╡ e51ca29f-470a-4be1-a6d9-2f96567d2bbf
all_isochrones_gaia = Dict(
	age => read_mist_file(joinpath(obs_dir, "../padova", "isochrone.gaiadr3.$(age)Gyrs.dat")) for age in [8, 9, 10, 12]
);

# ╔═╡ 3d0c8e48-a331-4768-bc28-08e28807af25
all_isochrones_delve = Dict(
	age => read_mist_file(joinpath(obs_dir, "../padova", "isochrone.decam.$(age)Gyrs.dat")) for age in [10, 12]
);

# ╔═╡ 6b9c74ee-e337-4f77-9f46-d2daff6d06fa
isochrones = phot_system == :gaia ? all_isochrones_gaia : all_isochrones_delve

# ╔═╡ d29a927f-51b2-4156-afac-945b179c3e9f
all_isochrones_ubvri = Dict(
	age => read_mist_file(joinpath(obs_dir, "../padova", "isochrone.ubvri.$(age)Gyrs.dat")) for age in [8, 9, 10, 12]
);

# ╔═╡ fbb7df70-3f62-48f0-a3b3-c64b356cc90d
function get_isochrone(isochrones, age, M_H; stage_max=5)

	if age ∉ keys(isochrones)
		throw(KeyError("$age not in available ages"))
	end
	
	isos = isochrones[age]
	@assert isapprox(isos.logAge[1],  log10(age) + 9, atol=0.01)
	if (M_H < minimum(isos.MH)) || (M_H > maximum(isos.MH))
		throw(DomainError(M_H, "metallicity out of isochrone range: $(extrema(isos.MH))"))
	end

	M_Hs = unique(isos.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]

	filt = isapprox.(isos.MH, M_H_adopted)
	filt .&= isos.label .<= stage_max # only keep through HB

	return isos[filt, :]
end

# ╔═╡ a2d4ab77-ff9b-45c7-9cad-4f34986f1296
function interpolate_magnitude(iso, mag_col)
	return LilGuys.lerp(iso.Mini, iso[!, mag_col])
end

# ╔═╡ e20663bd-028a-47ac-8efb-5f893321daa5
md"""
## Samplers
"""

# ╔═╡ 41aa41c8-8bab-42c9-b548-39aab368e8aa
function find_tot_number_stars(N_obs; mag_limit, N_max=300_000, distance_modulus=dm_0, mag_interp, mag_v_interp, mass_sampler)

	count = 0
	mass_tot = 0
	L_tot = 0
	N_tot = 0
	
	for i in 1:N_max
		mass = mass_sampler()
		mag = mag_interp(mass) + distance_modulus
		count += mag < mag_limit

		mag_v = mag_v_interp(mass) + distance_modulus 

		N_tot += 1
		L_tot += 10^(-0.4*mag_v)
		mass_tot += mass

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

# ╔═╡ 5b2975d7-5d15-4ffe-b328-ba7f64549a25
function find_tot_number_stars(N_obs, age, M_H; all_isochrones=isochrones,	gcol = gcol, kwargs...)

	iso = get_isochrone(all_isochrones, age, M_H, )
	mag_interp = interpolate_magnitude(iso, gcol)
	iso_v = get_isochrone(all_isochrones_ubvri, age, M_H)
	mag_v_interp = interpolate_magnitude(iso_v, "Vmag")

	
	M_max = maximum(iso.Mini)
	mass_sampler = Utils.create_kroupa_sampler(M_max, mass_min=0.08)


	return find_tot_number_stars(N_obs; mass_sampler=mass_sampler, mag_interp=mag_interp, mag_v_interp=mag_v_interp, kwargs...)
end

# ╔═╡ bc06d4d0-30b2-423e-bb49-f632db0f3a8d
get_isochrone(all_isochrones_gaia, 10, -2.177)

# ╔═╡ 5ee6fd79-3320-4993-ab3b-4501279645e4
function sample_tot_stars(obs_samples, N_samples; 
						  M_H, age, all_isochrones=isochrones,
						  gcol = gcol,
						  kwargs...)
	N_tots = Vector{Int}(undef, N_samples)
	MVs = Vector{Float64}(undef, N_samples)
	M_tots = Vector{Float64}(undef, N_samples)


	iso = get_isochrone(all_isochrones, age, M_H, )
	mag_interp = interpolate_magnitude(iso, gcol)
	iso_v = get_isochrone(all_isochrones_ubvri, age, M_H)
	mag_v_interp = interpolate_magnitude(iso_v, "Vmag")

	
	M_max = maximum(iso.Mini)
	mass_sampler = Utils.create_kroupa_sampler(M_max, mass_min=0.08)

	
	for i in 1:N_samples
		distance_modulus = dm_0 + dm_err*randn()

		N_tots[i], M_tots[i], MVs[i] = find_tot_number_stars(obs_samples[i]; 
			distance_modulus=distance_modulus,
			mass_sampler=mass_sampler, 
			mag_interp=mag_interp, mag_v_interp=mag_v_interp,
			kwargs...)
		
		if (i % 100) == 0
			@info "completed $i / $N_samples"
		end
	end

	return N_tots, M_tots, MVs
end

# ╔═╡ 4574c42f-e5c4-466a-a209-6fb30eab8cb4
function sample_tot_stars_age_MH(obs_samples, N_samples; 
						  M_H_sampler = ()->-2.19, age_sampler = ()->12, 
						  kwargs...)
	N_tots = Vector{Int}(undef, N_samples)
	MVs = Vector{Float64}(undef, N_samples)
	M_tots = Vector{Float64}(undef, N_samples)
	
	for i in 1:N_samples
		distance_modulus = dm_0 + dm_err*randn()
		M_H = M_H_sampler()
		age = age_sampler

		M_H = -2.19
		age = 12

		N_tots[i], M_tots[i], MVs[i] = find_tot_number_stars(
			obs_samples[i], age, M_H; 
			distance_modulus=distance_modulus,
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
N_obs_samples = mcmc_samples.N_sat

# ╔═╡ f4602567-9b79-4e01-8ffe-160d898cb9a9
Nsamples = min(max(1_000, length(N_obs_samples)), 300)

# ╔═╡ 0cbecf30-0449-4d03-9874-105c75bd7830
N_tots, M_tots, M_Vs = sample_tot_stars(N_obs_samples, Nsamples, M_H=-2.19, age=12, mag_limit=mag_limit)

# ╔═╡ dd5d0ecf-01c7-4c1b-a2d3-f65ac3ba32ab
M_Vs_free_age_MH = sample_tot_stars_age_MH(N_obs_samples, Nsamples, mag_limit=mag_limit, 
								   age_sampler = () -> rand([8, 9, 10, 12]), 
								   M_H_sampler= () -> rand(-1.5:-0.01:-2.19))[3]

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

# ╔═╡ a35b0151-d3c9-49fe-ab28-73eacc9c43f1
isochrones

# ╔═╡ ad6ce2cc-292a-45c1-b885-9d6927971c35
find_tot_number_stars(100, 12, -2.19, mag_limit=21, N_max=1e5)

# ╔═╡ 6ceba1bf-c390-4944-82d9-82edfd9a71c8
module MCMCUtils
	obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")

	include(joinpath(obs_dir, "mcmc_utils.jl"))
end

# ╔═╡ d50a8cbf-a6ca-4048-a593-1c645bff9b63
MCMCUtils.to_measurement(N_obs_samples)

# ╔═╡ 15369c73-7a50-443e-99a5-f03d692c86da
MCMCUtils.to_measurement(N_tots) / 1e3

# ╔═╡ 7195e648-af49-434b-920d-f26ccb2c741e
MCMCUtils.to_measurement(M_tots) / 1e3

# ╔═╡ 08ed0979-a3a1-42c6-b15f-206b25914d59
MCMCUtils.to_measurement(M_Vs)

# ╔═╡ 40628a8d-5757-420f-b02a-7d970e390131
MCMCUtils.to_measurement(M_Vs_free_age_MH)

# ╔═╡ 5e2c55c0-7ac2-430b-84fb-cac861132fde
md"""
# Sanity checks
"""

# ╔═╡ df906ecc-29a2-4c4b-a93d-265e1e381061
all_isochrones_gaia_old = Dict(
	age => read_mist_file(joinpath(obs_dir, "../padova", "isochrone.gaiadr2.weiler2018.$(age)Gyrs.dat")) for age in [8, 9, 10, 12]
)

# ╔═╡ 7e6bdccb-85d1-423c-b36f-e7fb49197556
all_isochrones_gaia[10]

# ╔═╡ ead5f264-f9b9-4f0c-9ab2-424207ef9f38
ssp_padova = read_mist_file(joinpath(ENV["DWARFS_ROOT"], "observations/padova/ssp1e4.gaiadr3.12Gyrs.fe_h_m2.dat"))

# ╔═╡ 0dbb147f-2db0-4c3a-91b5-fb38ee540be8
median(N_obs_samples) / LilGuys.mean(ssp_padova.Gmag .+ dm_0 .< mag_limit)

# ╔═╡ 046f04e2-e382-401a-9325-777b6fe4cfea
let
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10)

	stephist!(ssp_padova.Gmag .+ dm_0, bins=100)


	iso = get_isochrone(all_isochrones_gaia, 12, -2)
	iso_interp = interpolate_magnitude(iso, "Gmag")
	M_max = maximum(iso.Mini)
	sampler = Utils.create_kroupa_sampler(M_max, mass_min=0.09)

	mass = [sampler() for _ in 1:size(ssp_padova, 1)]
	mag = iso_interp.(mass) .+ dm_0

	stephist!(mag, bins=100)
		
	fig
end

# ╔═╡ ca8329d8-66f0-432b-a0a2-6d856c2e5fac
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale=log10,
			 yscale=log10, 
			 xticks=[0.01, 0.1, 1.0, 2.0],
			 yticks=[0.01, 0.1, 1, 10])

	sampler = Utils.create_kroupa_sampler(2, mass_min=0.01)
	masses = [sampler() for _ in 1:10_000]
		
	stephist!(masses, bins=logrange(0.01, 2, 100), normalization=:pdf)

	x = LinRange(0.01, 2, 1000)
	A = LilGuys.integrate(Utils.kroupa_imf, 0.01, 2)
	
	lines!(x, Utils.kroupa_imf.(x) / A )

	fig

end

# ╔═╡ 01cb561b-0e1f-4443-bc2c-8943d4c944cc
iso = get_isochrone(all_isochrones_gaia, 12, -2.19)

# ╔═╡ 1b003181-ecc0-4cb7-923d-eccab5ddac10
iso[!, :Gmag]

# ╔═╡ b8659aa3-17f9-48a2-b8ab-803056c9f600
let
	age = 8
	fe_h = -0.3
	iso1 = get_isochrone(all_isochrones_gaia, age, fe_h)
	iso2 = get_isochrone(all_isochrones_gaia_old, age, fe_h)
	iso3 = get_isochrone(all_isochrones_ubvri, age, fe_h)

	f =	lines(iso1.Mini, iso1.Gmag)

	lines!(iso2.Mini, iso2.Gmag, linestyle=:dot)
	lines!(iso3.Mini, iso3.Vmag)

	f
end

# ╔═╡ e24c999a-ca90-45af-818a-d936ed293cb0
let
	fe_h = -2.18
	age = 12
	
	iso = get_isochrone(all_isochrones_gaia, age, fe_h)

	iso2 = get_isochrone(all_isochrones_gaia_old, age, fe_h)

	f =	lines(iso.G_BPmag - iso.G_RPmag, iso.Gmag .+ obs_props["distance_modulus"])


	lines!(iso2.G_BPbrmag - iso2.G_RPmag, iso2.Gmag .+ obs_props["distance_modulus"])
		
	BP = interpolate_magnitude(iso, "G_BPmag")
	RP = interpolate_magnitude(iso, "G_RPmag")
	G = interpolate_magnitude(iso, "Gmag")

	x = LinRange(0.08, 0.92, 10000)
	lines!(BP.(x) - RP.(x), G.(x) .+ obs_props["distance_modulus"])

	f
end

# ╔═╡ e37f2a49-a2b2-433a-8ce6-c5f99e202f02


# ╔═╡ c7c17bbf-c117-4d43-87f7-3dc721076596
M_L_samples = let
	Nsamples = 1_000
	M_L = Vector{Float64}(undef, Nsamples)

	iso = get_isochrone(all_isochrones_ubvri, 12, -2.19)
	iso_interp = interpolate_magnitude(iso, "Vmag")
	mass_sampler = Utils.create_kroupa_sampler(iso.Mini[end])

	for i in 1:Nsamples
		masses = [mass_sampler() for _ in 1:50_000]
		mags = iso_interp.(masses)
		if i % 100 == 0
		end
		L = sum(LilGuys.mag_to_L.(mags))

		M_L[i] = sum(masses) / L
	end

	M_L
end

# ╔═╡ 2b861d38-d705-47d8-99a8-bf9065583fc0
hist(M_L_samples)

# ╔═╡ Cell order:
# ╠═ba2b61b8-ca92-49d0-9ea6-922a5e450c6a
# ╠═6eb7f393-3098-4e86-bab8-f17ae6f5c8a5
# ╠═5cdee12d-4899-41d7-aac0-d7538e5ea6de
# ╠═03c38e6b-f404-4337-99be-e100ec4f6d07
# ╠═f4602567-9b79-4e01-8ffe-160d898cb9a9
# ╠═39cd58d0-6d8e-40a4-8492-c30fe1782531
# ╠═791d6268-2c81-11f1-9366-e7a4214e567f
# ╠═d26b9e4b-d68f-4b72-a080-1d5c66248e19
# ╠═68adcf15-2683-4456-bf1e-0a1cd6907676
# ╠═78e8db76-054e-46ae-9ef0-5837b4811a9b
# ╠═f6e51d53-506e-4c17-8c6d-b044a9bef0cf
# ╠═f0894db7-2a6c-49c3-a97e-6a1fc04856be
# ╠═e4fdfb9b-73a8-410b-b86f-cdb1d71019a1
# ╠═0f977fa3-24ff-4f40-9407-0f783ea9ef26
# ╠═4de8788f-1c4b-4fa5-ab70-02c1fd1c5c4a
# ╠═20a73a2f-7326-4cee-8489-8b0b3a217828
# ╟─a4b3e72d-607f-43b0-869e-096646046eb2
# ╠═64f34c1f-6dbd-4cdf-a575-183cf0ede2ea
# ╠═e51ca29f-470a-4be1-a6d9-2f96567d2bbf
# ╠═3d0c8e48-a331-4768-bc28-08e28807af25
# ╠═6b9c74ee-e337-4f77-9f46-d2daff6d06fa
# ╠═d29a927f-51b2-4156-afac-945b179c3e9f
# ╠═fbb7df70-3f62-48f0-a3b3-c64b356cc90d
# ╠═a2d4ab77-ff9b-45c7-9cad-4f34986f1296
# ╟─e20663bd-028a-47ac-8efb-5f893321daa5
# ╠═41aa41c8-8bab-42c9-b548-39aab368e8aa
# ╠═5b2975d7-5d15-4ffe-b328-ba7f64549a25
# ╠═bc06d4d0-30b2-423e-bb49-f632db0f3a8d
# ╠═5ee6fd79-3320-4993-ab3b-4501279645e4
# ╠═4574c42f-e5c4-466a-a209-6fb30eab8cb4
# ╟─ce015a4e-3c9c-405f-afa6-59324216c591
# ╠═37f3a337-05a3-401d-91ed-c5c94eba899a
# ╠═4a6a1c7a-a4b2-468f-9db8-a814cdeee666
# ╠═1b003181-ecc0-4cb7-923d-eccab5ddac10
# ╠═0cbecf30-0449-4d03-9874-105c75bd7830
# ╠═dd5d0ecf-01c7-4c1b-a2d3-f65ac3ba32ab
# ╠═f1af266c-86d9-4649-bd6d-cd5c2b695892
# ╠═6ac311e3-268c-4e2d-95ca-77a643ce1bbb
# ╠═f66db8b0-1aef-4419-bbcf-84c9ef11d477
# ╠═bbd1541c-4fdd-40df-b03e-a1a60730336b
# ╠═443e2da9-5620-4139-bf58-c217d9708ef3
# ╠═a35b0151-d3c9-49fe-ab28-73eacc9c43f1
# ╠═ad6ce2cc-292a-45c1-b885-9d6927971c35
# ╠═0dbb147f-2db0-4c3a-91b5-fb38ee540be8
# ╠═6ceba1bf-c390-4944-82d9-82edfd9a71c8
# ╠═d50a8cbf-a6ca-4048-a593-1c645bff9b63
# ╠═15369c73-7a50-443e-99a5-f03d692c86da
# ╠═7195e648-af49-434b-920d-f26ccb2c741e
# ╠═08ed0979-a3a1-42c6-b15f-206b25914d59
# ╠═40628a8d-5757-420f-b02a-7d970e390131
# ╟─5e2c55c0-7ac2-430b-84fb-cac861132fde
# ╠═df906ecc-29a2-4c4b-a93d-265e1e381061
# ╠═7e6bdccb-85d1-423c-b36f-e7fb49197556
# ╠═ead5f264-f9b9-4f0c-9ab2-424207ef9f38
# ╠═046f04e2-e382-401a-9325-777b6fe4cfea
# ╠═ca8329d8-66f0-432b-a0a2-6d856c2e5fac
# ╠═01cb561b-0e1f-4443-bc2c-8943d4c944cc
# ╠═b8659aa3-17f9-48a2-b8ab-803056c9f600
# ╠═e24c999a-ca90-45af-818a-d936ed293cb0
# ╠═e37f2a49-a2b2-433a-8ce6-c5f99e202f02
# ╠═c7c17bbf-c117-4d43-87f7-3dc721076596
# ╠═2b861d38-d705-47d8-99a8-bf9065583fc0
