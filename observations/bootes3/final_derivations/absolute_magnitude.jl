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
	using OrderedCollections
	import TOML

	import StatsBase: median
end

# ╔═╡ 5cdee12d-4899-41d7-aac0-d7538e5ea6de
samplesname = "../final_derivations/mcmc/samples.j24_1c.mcmc_ell.csv"

# ╔═╡ f4602567-9b79-4e01-8ffe-160d898cb9a9
Nsamples = 10_000

# ╔═╡ 39cd58d0-6d8e-40a4-8492-c30fe1782531
mag_limit = 21

# ╔═╡ 10e1a046-ccf3-4428-a754-eff0220ad9e5
FIGDIR = "figures"

# ╔═╡ d26b9e4b-d68f-4b72-a080-1d5c66248e19
module Utils
	include("../delve_utils.jl")
end

# ╔═╡ 68adcf15-2683-4456-bf1e-0a1cd6907676
CairoMakie.activate!(type=:png)

# ╔═╡ 78e8db76-054e-46ae-9ef0-5837b4811a9b
read_mist_file = Utils.read_mist_file

# ╔═╡ f6e51d53-506e-4c17-8c6d-b044a9bef0cf
obs_dir = ".."

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

# ╔═╡ 6ceba1bf-c390-4944-82d9-82edfd9a71c8
module MCMCUtils
	obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")

	include(joinpath(obs_dir, "mcmc_utils.jl"))
end

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

# ╔═╡ 5ee6fd79-3320-4993-ab3b-4501279645e4
function sample_tot_stars(obs_samples, N_samples; 
						  M_H, age, all_isochrones=all_isochrones_gaia,
						  gcol = "Gmag",
						  kwargs...)
	N_tots = Vector{Int}(undef, N_samples)
	MVs = Vector{Float64}(undef, N_samples)
	M_tots = Vector{Float64}(undef, N_samples)
	dm = Vector{Float64}(undef, N_samples)


	iso = get_isochrone(all_isochrones, age, M_H, )
	mag_interp = interpolate_magnitude(iso, gcol)
	iso_v = get_isochrone(all_isochrones_ubvri, age, M_H)
	mag_v_interp = interpolate_magnitude(iso_v, "Vmag")

	
	M_max = maximum(iso.Mini)
	mass_sampler = Utils.create_kroupa_sampler(M_max, mass_min=0.08)

	
	for i in 1:N_samples
		distance_modulus = dm_0 + dm_err*randn()
		dm[i] = distance_modulus

		N_tots[i], M_tots[i], MVs[i] = find_tot_number_stars(obs_samples[i]; 
			distance_modulus=distance_modulus,
			mass_sampler=mass_sampler, 
			mag_interp=mag_interp, mag_v_interp=mag_v_interp,
			kwargs...)
		
		if (i % 100) == 0
			@info "completed $i / $N_samples"
		end
	end

	return DataFrame(
		:N_tot => N_tots, 
		:M_tot => M_tots, 
		:M_V => MVs,
		:distance_modulus => dm,
		:N_memb => obs_samples[1:N_samples],
		:M_L => M_tots ./ LilGuys.mag_to_L.(MVs),
		:M_ave => M_tots ./ N_tots,
		:m_v => MVs .+ dm,
	)
end

# ╔═╡ ce015a4e-3c9c-405f-afa6-59324216c591
md"""
# The main drag
"""

# ╔═╡ 37f3a337-05a3-401d-91ed-c5c94eba899a
mcmc_samples = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/mcmc", samplesname), DataFrame)

# ╔═╡ 4a6a1c7a-a4b2-468f-9db8-a814cdeee666
N_obs_samples = mcmc_samples.N_sat

# ╔═╡ 0cbecf30-0449-4d03-9874-105c75bd7830
df_luminosity = sample_tot_stars(N_obs_samples, Nsamples, M_H=-2.19, age=12, mag_limit=mag_limit)

# ╔═╡ 722fffd6-f538-41ae-8f1e-31637f8456ca
md"""
## Outputs
"""

# ╔═╡ ecd06b6d-6871-41f0-8e7f-b98ada858c46
summary = let
	df = OrderedDict()
	for col in names(df_luminosity)
		meas = MCMCUtils.to_measurement(df_luminosity[!, col])

		df[col] = middle(meas)
		df[col * "_el"] = lower_error(meas)
		df[col * "_ep"] = upper_error(meas)

	end

	df
end

# ╔═╡ 52becf7f-ec7d-46ef-bcab-3c6a1c6fa409
open("mcmc/luminosity_summary.toml", "w") do f
	TOML.print(f, summary)
end

# ╔═╡ f5e029fd-43e0-4bde-b8d9-0f0ad752e7d0
CSV.write("mcmc/luminosity_samples.csv", df_luminosity)

# ╔═╡ f1af266c-86d9-4649-bd6d-cd5c2b695892
md"""
# Analysis and saving
"""

# ╔═╡ 8e818c94-8692-4dd5-b913-9fed16b38015
import PairPlots: pairplot

# ╔═╡ 17caa9d0-c00b-444a-ae89-b74c00334de2
@savefig "corner.luminosity" pairplot(df_luminosity[:, [:N_memb, :distance_modulus, :M_V, :M_tot, :M_L ]])

# ╔═╡ 89ccacbd-a568-4985-ab41-93613ce60830
M_obs = df_luminosity.M_V + df_luminosity.distance_modulus

# ╔═╡ 6fd9fd56-1619-41d0-882e-763aa2f891ca
pairplot(df_luminosity[:, [:m_v]])

# ╔═╡ Cell order:
# ╠═5cdee12d-4899-41d7-aac0-d7538e5ea6de
# ╠═f4602567-9b79-4e01-8ffe-160d898cb9a9
# ╠═39cd58d0-6d8e-40a4-8492-c30fe1782531
# ╠═791d6268-2c81-11f1-9366-e7a4214e567f
# ╠═10e1a046-ccf3-4428-a754-eff0220ad9e5
# ╠═d26b9e4b-d68f-4b72-a080-1d5c66248e19
# ╠═68adcf15-2683-4456-bf1e-0a1cd6907676
# ╠═78e8db76-054e-46ae-9ef0-5837b4811a9b
# ╠═f6e51d53-506e-4c17-8c6d-b044a9bef0cf
# ╠═f0894db7-2a6c-49c3-a97e-6a1fc04856be
# ╠═e4fdfb9b-73a8-410b-b86f-cdb1d71019a1
# ╠═0f977fa3-24ff-4f40-9407-0f783ea9ef26
# ╠═4de8788f-1c4b-4fa5-ab70-02c1fd1c5c4a
# ╠═20a73a2f-7326-4cee-8489-8b0b3a217828
# ╠═6ceba1bf-c390-4944-82d9-82edfd9a71c8
# ╟─a4b3e72d-607f-43b0-869e-096646046eb2
# ╠═64f34c1f-6dbd-4cdf-a575-183cf0ede2ea
# ╠═e51ca29f-470a-4be1-a6d9-2f96567d2bbf
# ╠═d29a927f-51b2-4156-afac-945b179c3e9f
# ╠═fbb7df70-3f62-48f0-a3b3-c64b356cc90d
# ╠═a2d4ab77-ff9b-45c7-9cad-4f34986f1296
# ╟─e20663bd-028a-47ac-8efb-5f893321daa5
# ╠═41aa41c8-8bab-42c9-b548-39aab368e8aa
# ╠═5ee6fd79-3320-4993-ab3b-4501279645e4
# ╟─ce015a4e-3c9c-405f-afa6-59324216c591
# ╠═37f3a337-05a3-401d-91ed-c5c94eba899a
# ╠═4a6a1c7a-a4b2-468f-9db8-a814cdeee666
# ╠═0cbecf30-0449-4d03-9874-105c75bd7830
# ╟─722fffd6-f538-41ae-8f1e-31637f8456ca
# ╠═ecd06b6d-6871-41f0-8e7f-b98ada858c46
# ╠═52becf7f-ec7d-46ef-bcab-3c6a1c6fa409
# ╠═f5e029fd-43e0-4bde-b8d9-0f0ad752e7d0
# ╟─f1af266c-86d9-4649-bd6d-cd5c2b695892
# ╠═8e818c94-8692-4dd5-b913-9fed16b38015
# ╠═17caa9d0-c00b-444a-ae89-b74c00334de2
# ╠═89ccacbd-a568-4985-ab41-93613ce60830
# ╠═6fd9fd56-1619-41d0-882e-763aa2f891ca
