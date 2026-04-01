### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ fac66140-2d45-11f1-8f44-a5fe29dff650
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using PyFITS
	import PairPlots
	
	using Turing
end

# ╔═╡ c2fe9c8c-092c-45de-be8d-5d25dd445663
using DataFrames: DataFrame

# ╔═╡ 6bcac06e-ca2e-4f50-acc8-beecb327e295
using LilGuys

# ╔═╡ a21feb25-2b7a-4261-aef7-0099a0ca69eb
mock_number = 2

# ╔═╡ 196959fc-12e2-4ffb-bdd2-0164cf99711f
if mock_number == 1
	ra0 = 205
	dec0 = 22
	PA = 24
	R_h = 35
	ell = 0.35
	fe_h = -2.1
	age = 12
	N_sampled = 100_000
	prof = LilGuys.Plummer(r_s=R_h)
	
elseif mock_number == 2
	ra0 = 202.5
	dec0 = 27
	PA = 24
	R_h = 35
	ell = 0.0
	fe_h = -1.8
	age = 12
	N_sampled = 30_000
	dm = 18.2
	prof = LilGuys.Exp2D(R_s=R_h)
end

# ╔═╡ 90ad7849-7094-400d-afe2-a244f0369b84
R_max = 6*60

# ╔═╡ b8f236d3-4a35-4b47-b6d5-c7b85a96c116
md"""
# Setup
"""

# ╔═╡ 82942e69-7bdb-42ed-876c-76d0564ada97
module Utils
	include("delve_utils.jl")
end

# ╔═╡ 3ba43fc5-c713-47ea-b4a9-7fc1e14bf476
CairoMakie.activate!(type=:png)

# ╔═╡ 0a314494-6ed5-4cb1-b846-97dc02f048eb
allstars = read_fits("data/delve_dr2_good.fits")

# ╔═╡ 580a0f4d-850e-4191-92b8-31aa19a2043f
md"""
### Selecting a subset
We use the 12 degree field and select a reasonable circular subset to inject the mocks into
"""

# ╔═╡ 4db0e285-8ac4-4b2a-89b2-a2996d60f882
stars = let
	xi, eta = LilGuys.to_tangent(allstars.ra, allstars.dec, ra0, dec0)
	filt = @. sqrt(xi^2 + eta^2) < R_max / 60

	df = allstars[filt, :]
	df.xi, df.eta = 60 .* LilGuys.to_tangent(df.ra, df.dec, ra0, dec0)
	df
end

# ╔═╡ 6ae142ee-2d7f-436b-be8d-75bae3ed2849
let
	f = scatter(allstars.ra, allstars.dec, markersize=0.3, alpha=0.3)
	scatter!(stars.ra, stars.dec, markersize=0.5, alpha=0.05)
	f
end

# ╔═╡ 179b3161-3510-40ab-bfee-e906115e0fdb
md"""
## Creating the mock population
Creating the mock is as simple as
1. Sampling masses from the IMF
2. Convergin masses into mags using an isochrone
3. Adding typical observed errors
4. Sampling radii from the specified profile
5. converting radii back into tangent plane coordinates, given the system structural properties
"""

# ╔═╡ 49fca556-cdcf-406a-9c52-2a356e217b02
iso = Utils.get_isochrone(fe_h, age, stage_max=6)

# ╔═╡ 84c9e14a-e1d4-4b3b-bdbd-9065f74acf76
function sample_kroupa(N; mass_max, mass_min=0.08)
	masses = LinRange(mass_min, mass_max, 10_000)
	imf =  Utils.kroupa_imf.(masses)
	N_cdf = cumsum(imf)
	x = N_cdf .- N_cdf[1]
	x ./= x[end]
	mass_func =  LilGuys.lerp(x, masses)

	return mass_func.(rand(N))

end

# ╔═╡ 6e68eda4-d18f-439c-b1bf-be23142e6986
function get_isochrone_interp(iso; cols=["gmag", "rmag", "Gmag", "G_BPftmag", "G_RPmag"])
	interps = Dict(col => LilGuys.lerp(iso.Mini, iso[!, col]) for col in cols if col in names(iso))
end 

# ╔═╡ c6862981-6155-41ab-9b0c-092608e8de30
iso_interp = get_isochrone_interp(iso)

# ╔═╡ 58b0e3b0-2cb6-4bb5-83c3-8c02fd1656d2
function sample_mags(N, mag_min=23.0)
	masses = sample_kroupa(N, mass_max=iso.Mini[end])
	gmag = iso_interp["gmag"].(masses) .+ dm
	rmag = iso_interp["rmag"].(masses) .+ dm

	g_err = Utils.delve_g_err.(gmag)
	gr_err = Utils.delve_gr_err.(gmag) 
	r_err = @. sqrt(gr_err^2 - g_err^2)
	d_g = g_err .* randn(N)
	d_r = r_err .* randn(N)
	gmag_obs = gmag .+ d_g 
	rmag_obs = rmag .+ d_r


	filt = gmag .< mag_min
	return DataFrame(Dict(
		:Mini => masses[filt], 
		:gmag => gmag_obs[filt], 
		:rmag => rmag_obs[filt]
	)
	)
end

# ╔═╡ f13275da-fc59-4141-b523-2d527a3b0baa
df_mags = sample_mags(N_sampled)

# ╔═╡ 18e8a064-acb5-40d8-8cbf-9cf4829ffdb1
md"""
The input N_samples controls the number of mass draws. However, the number of observed stars is much less (since most stars are below the magnitude limit)
"""

# ╔═╡ c45e15bc-511c-45a9-80c0-21df5e2a7edf
N_star = size(df_mags, 1)

# ╔═╡ 96288a8c-8ffb-423b-8b0b-2eee4af08cc0
sum(df_mags.gmag .< 20)

# ╔═╡ a84dd15f-a19d-4743-a7c2-43313feee860
sum(df_mags.gmag .< 21)

# ╔═╡ 3954b4fc-6633-4bce-a289-10f065bac88a
md"""
### Sampling density profile & converting to tangent plane coordinates
"""

# ╔═╡ 3746aed6-692d-46b6-bcfa-402dd16fe093
R_ell = LilGuys.sample_surface_density(x->surface_density(prof, x), N_star)

# ╔═╡ e446a8fb-eb33-45aa-b769-f26bfc71f766
θs = 360 * rand(N_star)

# ╔═╡ 9338f12a-9bc7-4a18-be99-f7cac7f247ef
aspect = sqrt(1 - ell)

# ╔═╡ 707880b1-f6b7-4eff-b0eb-44712c8b8628
xi_p = R_ell .* cosd.(θs) *aspect 

# ╔═╡ 709d1b81-50f4-484a-bd2f-31d052b97a77
eta_p = R_ell .* sind.(θs) / aspect

# ╔═╡ c2d9c678-df68-40a0-9603-0b7d49a4cb49
xi = xi_p * cosd(PA) + eta_p * sind(PA)

# ╔═╡ 56375606-f5bb-4b4c-927e-239d9b0dd69c
eta = -xi_p * sind(PA) + eta_p * cosd(PA)

# ╔═╡ 6a25ed78-d747-4536-be73-091641fb9513
dm

# ╔═╡ d4b19fe0-1b18-4232-a8ea-c1a900d1e94a
R_ell_2 = LilGuys.calc_R_ell(xi, eta, ell, PA)

# ╔═╡ 50aeb85a-e414-4b86-a809-fbb3dd2a66fa
@assert LilGuys.std(R_ell_2 ./ R_ell) < 1e-8

# ╔═╡ ee4b2cde-d1f3-4d22-9147-6cbbd9649bc4
df_new = let
	df = vcat(
			stars,
			DataFrame(Dict(
				:xi => xi,
				:eta => eta,
				:gmag => df_mags.gmag,
				:rmag => df_mags.rmag,
				:Mini => df_mags.Mini,
			)),
			cols=:union
		)

	radec = LilGuys.from_tangent.(df.xi/60, df.eta/60, ra0, dec0)
	
	df[:, :ra], df[:, :dec] = first.(radec), last.(radec)

	df[:, :R] = @. sqrt(df.xi^2 + df.eta^2)

	df[df.R .< R_max, :]
end

# ╔═╡ 95493f25-7039-424e-80e9-a8e27394ee7a
f_sat = LilGuys.mean(.!ismissing.(df_new.Mini))

# ╔═╡ ba84df00-6f0c-4335-9d54-f9c6a17601e0
N_star

# ╔═╡ d8c122ff-4d85-453e-8b56-c0db61a182f7
N_star / size(df_new, 1)

# ╔═╡ 7c95c1f8-60e8-45dc-adfa-3b6331c96983
LilGuys.mean(.!ismissing.(df_new.Mini)[df_new.gmag .< 21])

# ╔═╡ f2602eb9-2ce0-4f0a-93c1-a9f921b69643
md"""
# What does the sample look like?
"""

# ╔═╡ 19078c3d-dbe7-4ba6-9461-8b4bff8c6d78
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "xi", 
			 ylabel = "eta", 
			 aspect = DataAspect(),
			 limits=(-R_max, R_max, -R_max, R_max))

	scatter!(df_new.xi, df_new.eta, markersize=0.3, color=:black, alpha=0.3)


	fig
end

# ╔═╡ f1aa86de-7d10-4f17-bf79-3debde4f1750
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "xi", 
			 ylabel = "eta", 
			 aspect = DataAspect(),
			 limits=(-R_max, R_max, -R_max, R_max))

	scatter!(stars.xi, stars.eta, markersize=0.3, color=:black, alpha=0.3)

	scatter!(xi, eta, markersize=1, color=COLORS[2])

	fig
end

# ╔═╡ dddebbd1-110d-4bf7-b973-e308d5b8a646
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 yreversed=true)

	scatter!(df_new.gmag .- df_new.rmag, df_new.gmag, markersize=0.3, alpha=0.3, color=:black)

	fig
end

# ╔═╡ 8e6911ea-9f07-49b5-9109-f41f7f9a3773
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 yreversed=true)

	scatter!(df_new.gmag .- df_new.rmag, df_new.gmag, markersize=0.3, alpha=0.3, color=:black)

	scatter!(df_mags.gmag .- df_mags.rmag, df_mags.gmag, color=COLORS[2], markersize=1)
	# lines!(iso.gmag .- iso.rmag, iso.gmag .+ dm)

	fig
end

# ╔═╡ cf114609-3cad-47fa-90c5-ae705aabe6db
write_fits("data/delve_dr2_good_mock.fits", df_new, overwrite=true)

# ╔═╡ Cell order:
# ╠═a21feb25-2b7a-4261-aef7-0099a0ca69eb
# ╠═196959fc-12e2-4ffb-bdd2-0164cf99711f
# ╠═90ad7849-7094-400d-afe2-a244f0369b84
# ╠═b8f236d3-4a35-4b47-b6d5-c7b85a96c116
# ╠═fac66140-2d45-11f1-8f44-a5fe29dff650
# ╠═82942e69-7bdb-42ed-876c-76d0564ada97
# ╠═c2fe9c8c-092c-45de-be8d-5d25dd445663
# ╠═6bcac06e-ca2e-4f50-acc8-beecb327e295
# ╠═3ba43fc5-c713-47ea-b4a9-7fc1e14bf476
# ╠═0a314494-6ed5-4cb1-b846-97dc02f048eb
# ╠═580a0f4d-850e-4191-92b8-31aa19a2043f
# ╠═6ae142ee-2d7f-436b-be8d-75bae3ed2849
# ╠═4db0e285-8ac4-4b2a-89b2-a2996d60f882
# ╠═179b3161-3510-40ab-bfee-e906115e0fdb
# ╠═49fca556-cdcf-406a-9c52-2a356e217b02
# ╠═c6862981-6155-41ab-9b0c-092608e8de30
# ╠═84c9e14a-e1d4-4b3b-bdbd-9065f74acf76
# ╠═6e68eda4-d18f-439c-b1bf-be23142e6986
# ╠═58b0e3b0-2cb6-4bb5-83c3-8c02fd1656d2
# ╠═f13275da-fc59-4141-b523-2d527a3b0baa
# ╠═18e8a064-acb5-40d8-8cbf-9cf4829ffdb1
# ╠═c45e15bc-511c-45a9-80c0-21df5e2a7edf
# ╠═96288a8c-8ffb-423b-8b0b-2eee4af08cc0
# ╠═a84dd15f-a19d-4743-a7c2-43313feee860
# ╠═3954b4fc-6633-4bce-a289-10f065bac88a
# ╠═3746aed6-692d-46b6-bcfa-402dd16fe093
# ╠═e446a8fb-eb33-45aa-b769-f26bfc71f766
# ╠═9338f12a-9bc7-4a18-be99-f7cac7f247ef
# ╠═707880b1-f6b7-4eff-b0eb-44712c8b8628
# ╠═709d1b81-50f4-484a-bd2f-31d052b97a77
# ╠═c2d9c678-df68-40a0-9603-0b7d49a4cb49
# ╠═56375606-f5bb-4b4c-927e-239d9b0dd69c
# ╠═6a25ed78-d747-4536-be73-091641fb9513
# ╠═d4b19fe0-1b18-4232-a8ea-c1a900d1e94a
# ╠═50aeb85a-e414-4b86-a809-fbb3dd2a66fa
# ╠═ee4b2cde-d1f3-4d22-9147-6cbbd9649bc4
# ╠═95493f25-7039-424e-80e9-a8e27394ee7a
# ╠═ba84df00-6f0c-4335-9d54-f9c6a17601e0
# ╠═d8c122ff-4d85-453e-8b56-c0db61a182f7
# ╠═7c95c1f8-60e8-45dc-adfa-3b6331c96983
# ╠═f2602eb9-2ce0-4f0a-93c1-a9f921b69643
# ╠═19078c3d-dbe7-4ba6-9461-8b4bff8c6d78
# ╠═f1aa86de-7d10-4f17-bf79-3debde4f1750
# ╠═dddebbd1-110d-4bf7-b973-e308d5b8a646
# ╠═8e6911ea-9f07-49b5-9109-f41f7f9a3773
# ╠═cf114609-3cad-47fa-90c5-ae705aabe6db
