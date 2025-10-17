### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 5bf5bd8d-65db-43e0-a999-5289d5cd1acf
using OrderedCollections

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

# ╔═╡ 629e397d-a4ae-469d-b81c-f92069d0b28a
import DensityEstimators

# ╔═╡ 8ae7dc9a-9f56-4540-bbe4-d98f8ea2bdc9
import StatsBase: weights

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ 5abb6cec-e947-4fc1-9848-760e50bd5628
function get_stars_final(galaxyname, modelname, starsname, filename="final.fits")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/stars/$starsname/")

	return read_fits(joinpath(modeldir, filename))
end

# ╔═╡ 3ed68a46-04aa-4fea-b484-d9f3cc158738
smallfontsize=0.8*theme(:fontsize)[]

# ╔═╡ cb898eeb-0803-42a7-a2c9-d6e2b95f8945
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ 5acf2849-d5bf-49df-94fe-7c3fe738f921
"""
Given a set of radii and velocity
"""
function calc_σv(r_ell, rv, mass;  r_max = 30)
	filt = r_ell .< r_max
	vel = rv[filt]

	return LilGuys.std(vel, weights(mass)[filt])
end

# ╔═╡ dd10823e-cf9a-4c5f-9cab-8d8f30effa33
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ cd0ea281-28b5-46fe-8ca2-23e44c623da6
function plot_density_f!(galaxyname, modelname, starsname; R_shift=0, norm_shift=0, kwargs...)
	prof_i, prof_f, norm = ModelUtils.load_stellar_profiles(galaxyname, modelname, starsname, norm_shift=norm_shift)

	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/stars/$starsname/")
	dist_f =  TOML.parsefile(model_dir * "../../orbital_properties.toml")["distance_f"]
	
	lines!(prof_f.log_R .+ R_shift, prof_f.log_Sigma; kwargs...)

end

# ╔═╡ 0cc4c7b4-d707-4f3d-9a82-1b6eb4119143
function plot_σv_prof!(galaxyname, args...; R_shift=0, kwargs...)
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

	lines!(LilGuys.midpoints(r_bins), log10.(σs) .- 1/2 * R_shift; kwargs...)


	
end

# ╔═╡ 908d84f6-18b2-4f33-abb4-537cd2fdea52
df_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/velocities/processed/vz_r_ell_binned.rv_combined_x_wide_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ be9267a9-eebc-4edf-bde4-c224341f0ff5
df_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/velocities/processed/vz_r_ell_binned.rv_combined_x_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ ed028d05-86b4-4752-9462-2a837fe11e55
function plot_σv_obs!(galaxyname)
	if galaxyname == "sculptor"
		df = df_scl
	else
		df = df_umi
	end
	
	errorscatter!(log10.(df.x), log10.(df.σ), yerror=df.σ_em ./ df.σ ./ log(10), color=:black, markersize=6)
end

# ╔═╡ 6368512b-30ed-4457-b419-69dfef2721c7
function plot_σv_prof(args...; kwargs...)
	
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "log R / arcmin",
		ylabel=L"$\log\,\sigma_\textrm{v, los}$ / km s$^{-1}$",
		limits=((0.0, 2.5), (0, 2))
	)

	plot_σv_prof!(args...; kwargs...)


	plot_σv_obs!(args[1])

	fig
	
end

# ╔═╡ 767186a1-090f-426c-862c-892e5e728ba7
function compare_both(galaxyname, modelnames; density_kwargs = Dict()) 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel="log radius / arcmin", ylabel = "log surface density", limits=(-0.5, 3, -3, 2))

	ModelUtils.plot_expected_profile!(galaxyname)

	for (i, (label, model)) in enumerate(modelnames)
		kwargs = get(density_kwargs, label, (;))
		plot_density_f!(model..., label=label; color=COLORS[i], kwargs...)
	end
	
	axislegend(position=:lb)
	hidexdecorations!(ax, ticks=false, minorticks=false)


	ax_v = Axis(fig[2, 1], xlabel = "log radius / arcmin", ylabel = "log velocity dispersion")
	
	plot_σv_obs!(galaxyname)


	for (i, (label, model)) in enumerate(modelnames)
		kwargs = get(density_kwargs, label, Dict{Symbol, Any}())
		plot_σv_prof!(model..., color=COLORS[i], label=label)
	end



	ylims!(0.6, 1.2)

	xlims!(0, 2.5)
	linkxaxes!(ax, ax_v)

	fig
end

# ╔═╡ 5c56c73e-d8db-4894-b0cf-f7f1b0e01a72
md"""
# Plot
"""

# ╔═╡ 9c3caa8f-cb12-4b40-b091-9e570af0f074
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 2c5eca31-bea1-4805-b869-877ec54828ca
plot_σv_prof("sculptor", "1e4_exp_M3e-4_r0.1/orbit_smallperi", "stars")

# ╔═╡ 4a7222a7-4739-4059-9b91-1e317e2e2662
plot_σv_prof("sculptor", "1e4_exp_M3e-4_r0.1/smallperilmc", "stars")

# ╔═╡ 8b791ca5-3747-48f4-b410-a5b3523777a5
compare_both("ursa_minor", OrderedDict(
	"dm_free" => ("1e4_exp_M4e-4_r0.13/orbit_smallperi", "stars")
)
			)

# ╔═╡ 23837d49-6bec-4a86-ad39-0cb66c655f02
ModelUtils.compare_both("sculptor", "1e4_exp_M3e-4_r0.1/orbit_smallperi", "stars", r_j=true)

# ╔═╡ c61b48ac-62a4-466e-8238-a7c8ad5638f5
ModelUtils.compare_both("sculptor", "1e4_exp_M3e-4_r0.1/smallperilmc", "stars", r_j=true)

# ╔═╡ 1d917d95-bcb3-41b7-a00c-b7b6959d69bc
ModelUtils.compare_both("ursa_minor", "1e4_exp_M4e-4_r0.13/orbit_smallperi", "stars", norm_shift=-2)

# ╔═╡ 57b7ae43-1ab1-49ab-a915-e0064fade974
ModelUtils.compare_both(modelnames["oblate"]..., norm_shift=-0.8)

# ╔═╡ 8888a121-61e0-46ba-81f6-5b4047b29620
@savefig "scl_extra_densities_sigma_vs" compare_both("sculptor", OrderedDict(
	"fiducial" => modelnames["scl_smallperi"],
	"MW impact" => modelnames["mw_impact"],
	"anisotropic" =>  modelnames["anisotropic"],
	"oblate" => modelnames["oblate"],
	# "dm free" => ("1e4_exp_M3e-4_r0.1/orbit_smallperi", "stars")
),
	density_kwargs=Dict(
		"anisotropic" => (; norm_shift=0, R_shift=-0.25),
		"oblate" => (; norm_shift=-0.5, R_shift=-0.1),
	))

# ╔═╡ 108e5761-fdc8-4729-96b1-26b9fc886606


# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═5bf5bd8d-65db-43e0-a999-5289d5cd1acf
# ╠═51f50a8d-1684-47f8-af6d-4cea827a4961
# ╠═2bacd818-4985-4922-85a3-716bdfda5146
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═9309c10c-6ba3-436c-b975-36d26dafb821
# ╠═629e397d-a4ae-469d-b81c-f92069d0b28a
# ╠═8ae7dc9a-9f56-4540-bbe4-d98f8ea2bdc9
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═5abb6cec-e947-4fc1-9848-760e50bd5628
# ╠═3ed68a46-04aa-4fea-b484-d9f3cc158738
# ╠═cb898eeb-0803-42a7-a2c9-d6e2b95f8945
# ╠═5acf2849-d5bf-49df-94fe-7c3fe738f921
# ╠═dd10823e-cf9a-4c5f-9cab-8d8f30effa33
# ╠═cd0ea281-28b5-46fe-8ca2-23e44c623da6
# ╠═ed028d05-86b4-4752-9462-2a837fe11e55
# ╠═6368512b-30ed-4457-b419-69dfef2721c7
# ╠═0cc4c7b4-d707-4f3d-9a82-1b6eb4119143
# ╠═767186a1-090f-426c-862c-892e5e728ba7
# ╠═908d84f6-18b2-4f33-abb4-537cd2fdea52
# ╠═be9267a9-eebc-4edf-bde4-c224341f0ff5
# ╟─5c56c73e-d8db-4894-b0cf-f7f1b0e01a72
# ╠═9c3caa8f-cb12-4b40-b091-9e570af0f074
# ╠═2c5eca31-bea1-4805-b869-877ec54828ca
# ╠═4a7222a7-4739-4059-9b91-1e317e2e2662
# ╠═8b791ca5-3747-48f4-b410-a5b3523777a5
# ╠═23837d49-6bec-4a86-ad39-0cb66c655f02
# ╠═c61b48ac-62a4-466e-8238-a7c8ad5638f5
# ╠═1d917d95-bcb3-41b7-a00c-b7b6959d69bc
# ╠═57b7ae43-1ab1-49ab-a915-e0064fade974
# ╠═8888a121-61e0-46ba-81f6-5b4047b29620
# ╠═108e5761-fdc8-4729-96b1-26b9fc886606
