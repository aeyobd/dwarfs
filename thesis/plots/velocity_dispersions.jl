### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ 51f50a8d-1684-47f8-af6d-4cea827a4961
using DataFrames, CSV

# ╔═╡ 9309c10c-6ba3-436c-b975-36d26dafb821
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 2bacd818-4985-4922-85a3-716bdfda5146
import DensityEstimators: histogram2d

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ e084a02a-f445-422f-b0b3-700a44bf204c
function get_R_h(galaxyname)

	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))
	R_h = obs_props["R_h"]
end

# ╔═╡ b9e03109-9f49-4897-b26c-e31698a5fe49
function get_r_b(galaxyname, haloname, orbitname, starsname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/stars/$starsname/")

	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")

	σv = prof_f.annotations["sigma_v"]
	if lmc
		props = TOML.parsefile(model_dir * "../../orbital_properties_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "../../orbital_properties.toml")
	end

	dist_f =  TOML.parsefile(model_dir * "../../orbital_properties.toml")["distance_f"]

	
	dt = props["t_last_peri"]
	r_b = LilGuys.break_radius(σv / V2KMS, dt / T2GYR)
	R_h = get_R_h(galaxyname)

	return LilGuys.kpc2arcmin(r_b, dist_f)	/ R_h
end

# ╔═╡ fe3bc6ee-14ed-4006-b3ec-f068d2492da4
function get_r_j(galaxyname, haloname, orbitname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/")

	if lmc
		props = TOML.parsefile(model_dir * "jacobi_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "jacobi.toml")
	end
	R_h = get_R_h(galaxyname)

	return props["r_J"] / R_h
end

# ╔═╡ 5abb6cec-e947-4fc1-9848-760e50bd5628
function get_stars_final(galaxyname, haloname, orbitname, starsname, filename="final.fits")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/stars/$starsname/")

	return read_fits(joinpath(modeldir, filename))
end

# ╔═╡ 3ed68a46-04aa-4fea-b484-d9f3cc158738
smallfontsize=0.8*theme(:fontsize)[]

# ╔═╡ cb898eeb-0803-42a7-a2c9-d6e2b95f8945
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ 629e397d-a4ae-469d-b81c-f92069d0b28a
import DensityEstimators

# ╔═╡ 8ae7dc9a-9f56-4540-bbe4-d98f8ea2bdc9
import StatsBase: weights

# ╔═╡ 5acf2849-d5bf-49df-94fe-7c3fe738f921
"""
Given a set of radii and velocity
"""
function calc_σv(r_ell, rv, mass;  r_max = 30)
	filt = r_ell .< r_max
	vel = rv[filt]

	return LilGuys.std(vel, weights(mass)[filt])
end

# ╔═╡ 0cc4c7b4-d707-4f3d-9a82-1b6eb4119143
function plot_σv_prof!(galaxyname, args...; kwargs...)
	stars = get_stars_final(galaxyname, args...)


	x = log10.(stars.r_ell)
	w = stars.weights
	v = stars.radial_velocity

	filt = @. !isnan(x)

	x = x[filt]
	w = w[filt]
	v = v[filt]

	bins = 80
	
	r_bins = DensityEstimators.make_bins(x, (-1, 3), DensityEstimators.bins_equal_number, n=bins)

	
	σs = Vector{Float64}(undef, bins)
	err = Vector{Float64}(undef, bins)
	Ns = Vector{Float64}(undef, bins)
	
	for i in 1:bins
		filt = r_bins[i] .<= x .< r_bins[i+1]
		σs[i] = calc_σv(x[filt], v[filt], w[filt], r_max=Inf)
		N = sum(filt)
		if N > 1
			err[i] = σs[i] * sqrt(2/(N - 1))
		else
			err[i] = NaN
		end
		Ns[i] = N
	end

	lines!(LilGuys.midpoints(r_bins), log10.(σs); kwargs...)


	
end

# ╔═╡ 908d84f6-18b2-4f33-abb4-537cd2fdea52
df_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/velocities/processed/vz_r_ell_binned.rv_combined_x_wide_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ be9267a9-eebc-4edf-bde4-c224341f0ff5
df_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/velocities/processed/vz_r_ell_binned.rv_combined_x_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ 6368512b-30ed-4457-b419-69dfef2721c7
function plot_σv_prof(args...; kwargs...)
	
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "log R / arcmin",
		ylabel=L"$\log\,\sigma_\textrm{v, los}$ / km s$^{-1}$",
		limits=((0.0, 2.5), (0, 2))
	)

	plot_σv_prof!(args...; kwargs...)

	if args[1] == "sculptor"
		df = df_scl
	else
		df = df_umi
	end
	
	errorscatter!(log10.(df_scl.x), log10.(df_scl.σ), yerror=df_scl.σ_em ./ df_scl.σ ./ log(10), color=:black, markersize=6)

	fig
	
end

# ╔═╡ affe66b1-dc16-4c07-b945-88513cf554e8
let 
	fig = plot_σv_prof("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi", "exp2d_rs0.10", label="smallperi")

	plot_σv_prof!("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi", "plummer_rs0.20", label="smallperi-Plummer")
	
	plot_σv_prof!("sculptor", "1e7_new_v25_r2.5", "smallperilmc", "exp2d_rs0.11", label="LMC flyby")

	plot_σv_prof!("sculptor", "1e6_new_v31_r3.2", "L3M11_9Gyr_smallperi.a4", "exp2d_rs0.10", label="MW impact")

	plot_σv_prof!("sculptor", "1e6_v43_r5_beta0.2_a4", "orbit_smallperi", "exp2d_rs0.13", label="anisotropy")


	axislegend(position=:lt)

	ylims!(0.6, 1.6)

	@savefig "scl_sigma_v_profiles"

	fig
end

# ╔═╡ 108e5761-fdc8-4729-96b1-26b9fc886606
let 
	fig = plot_σv_prof("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5", "exp2d_rs0.10", label="smallperi")

	plot_σv_prof!("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5", "plummer_rs0.20", label="smallperi-Plummer")


	plot_σv_prof!("ursa_minor", "1e6_v37_r5.0", "orbit_mean.2", "exp2d_rs0.10", label="mean")

	axislegend(position=:lt)

	ylims!(0.6, 1.6)

	@savefig "umi_sigma_v_profiles"


	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═51f50a8d-1684-47f8-af6d-4cea827a4961
# ╠═2bacd818-4985-4922-85a3-716bdfda5146
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═9309c10c-6ba3-436c-b975-36d26dafb821
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e084a02a-f445-422f-b0b3-700a44bf204c
# ╠═b9e03109-9f49-4897-b26c-e31698a5fe49
# ╠═fe3bc6ee-14ed-4006-b3ec-f068d2492da4
# ╠═5abb6cec-e947-4fc1-9848-760e50bd5628
# ╠═3ed68a46-04aa-4fea-b484-d9f3cc158738
# ╠═cb898eeb-0803-42a7-a2c9-d6e2b95f8945
# ╠═5acf2849-d5bf-49df-94fe-7c3fe738f921
# ╠═629e397d-a4ae-469d-b81c-f92069d0b28a
# ╠═8ae7dc9a-9f56-4540-bbe4-d98f8ea2bdc9
# ╠═6368512b-30ed-4457-b419-69dfef2721c7
# ╠═0cc4c7b4-d707-4f3d-9a82-1b6eb4119143
# ╠═affe66b1-dc16-4c07-b945-88513cf554e8
# ╠═108e5761-fdc8-4729-96b1-26b9fc886606
# ╠═908d84f6-18b2-4f33-abb4-537cd2fdea52
# ╠═be9267a9-eebc-4edf-bde4-c224341f0ff5
