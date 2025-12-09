### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ dd10823e-cf9a-4c5f-9cab-8d8f30effa33
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 5acf2849-d5bf-49df-94fe-7c3fe738f921
"""
Given a set of radii and velocity
"""
function calc_σv(r_ell, rv, mass;  r_max = 30)
	filt = r_ell .< r_max
	vel = rv[filt]

	return LilGuys.std(vel, weights(mass)[filt])
end

# ╔═╡ be52a89e-0709-4c93-83f2-c2a3f177274d
function calc_sigma_v_prof(stars, w = stars.weights)
	x = log10.(stars.r_ell)
	
	v = stars.radial_velocity

	filt = @. !isnan(x)

	x = x[filt]
	w = w[filt]
	v = v[filt]

	bins = 80
	
	r_bins = DensityEstimators.make_bins(x, (-1, 3), DensityEstimators.bins_equal_number, n=bins)

	
	σs = Vector{Float64}(undef, bins)
	Ns = Vector{Float64}(undef, bins)
	Ms = Vector{Float64}(undef, bins)
	
	for i in 1:bins
		filt = r_bins[i] .<= x .< r_bins[i+1]
		σs[i] = LilGuys.std(v[filt], weights(w[filt]))
		N = sum(filt)
		Ns[i] = N
		Ms[i] = sum(w[filt])
	end

	return r_bins, σs, Ns, Ms
end

# ╔═╡ 0cc4c7b4-d707-4f3d-9a82-1b6eb4119143
function plot_σv_prof!(galaxyname, args...; filename="final.fits", R_shift=0, kwargs...)
	stars = get_stars_final(galaxyname, args..., filename)


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


	lines!(LilGuys.midpoints(r_bins), σs ./ sqrt(10^R_shift); kwargs...)


	
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
	
	errorscatter!(log10.(df.x), (df.σ), yerror=df.σ_em , color=:black, markersize=6)
end

# ╔═╡ 9792d244-c295-479f-bbc5-157c5a07c300
function plot_σv_prof(gs::Makie.GridPosition, args...; kwargs...)
	ax = Axis(gs,
		xlabel = "log R / arcmin",
		ylabel=L"$\sigma_\textrm{v, los}$ / km s$^{-1}$",
		limits=((0.0, 2.5), (4, 15))
	)

	plot_σv_prof!(args...; kwargs...)
	# plot_σv_prof!(args...; filename="initial.fits", kwargs...)


	plot_σv_obs!(args[1])

	ax
	
end

# ╔═╡ 6368512b-30ed-4457-b419-69dfef2721c7
function plot_σv_prof(args...; limits=((0.0, 2.5), (4, 15)), kwargs...)
	
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "log R / arcmin",
		ylabel=L"$\Sigma_\textrm{v, los}$ / km s$^{-1}$",
		limits=limits
	)

	plot_σv_prof!(args...; kwargs...)
	plot_σv_prof!(args...; filename="initial.fits", linestyle=:dot, kwargs...)


	plot_σv_obs!(args[1])

	fig
	
end

# ╔═╡ 5c56c73e-d8db-4894-b0cf-f7f1b0e01a72
md"""
# Plot
"""

# ╔═╡ 9c3caa8f-cb12-4b40-b091-9e570af0f074
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 2c5eca31-bea1-4805-b869-877ec54828ca
plot_σv_prof(modelnames["scl_smallperi"]..., limits=((0.0, 5), (4, 50)))

# ╔═╡ 14d85541-d6aa-4d00-8c37-f84126cb89b3
plot_σv_prof(modelnames["scl_smallperi_2exp"]..., limits=((0.0, 5), (4, 50)))

# ╔═╡ a212281f-2789-441b-9d0b-5712fef3c8ec
plot_σv_prof(modelnames["umi_smallperi"]...)

# ╔═╡ 405c74a8-075d-48ab-9ae8-827779952a0b
plot_σv_prof(modelnames["umi_smallperi_2exp"]...)

# ╔═╡ 26b99d47-5024-4619-b62f-e69ab24361fd


# ╔═╡ b9ba3889-b77a-406d-b488-2c7bd7db00e0
function plot_prof!(prof; kwargs...)
	lines!(LilGuys.log_radii(prof), LilGuys.log_surface_density(prof); kwargs...)
end

# ╔═╡ 4abe36c0-d10e-4d9b-8786-8a7f7aea593a
function compare_density(gs, modelname_exp; modelname_plummer=nothing, limits=(-0.5, 2.5, -6, 2.5), lmc=false, y_low_R_h=nothing, y_low=nothing)

	galaxy = modelname_exp[1]
	ax = Axis(gs,
		xlabel = "log Radii / arcmin",
		ylabel = L"log $\Sigma$ / stars arcmin$^{-2}$",
		limits = limits,
			  title=Dict("sculptor"=>"Sculptor", "ursa_minor" => "Ursa Minor")[galaxy]
			 )
	ModelUtils.plot_expected_profile!(galaxy)

	prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_exp...)

	plot_prof!(prof_i, color=COLORS[1], linewidth=smalllinewidth, linestyle=:dot, label="exp initial")
	plot_prof!(prof_f, color=COLORS[1], linewidth=smalllinewidth, linestyle=:solid, label="exp final")
	lines!([NaN], [NaN], linewidth=0, label=" ")


	if !isnothing(modelname_plummer)
		prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_plummer...)
	
		plot_prof!(prof_i, color=COLORS[5], linestyle=:dot, label="2exp initial")
		plot_prof!(prof_f, color=COLORS[5], linestyle=:solid, label="2exp final")
	end


	limits!(limits...)
	if isnothing(y_low)
		y_low = limits[3]
	end
	r_b = ModelUtils.get_r_b(modelname_exp..., lmc=lmc)
	ModelUtils.plot_r_break_arrow!(r_b, y_low)

	r_j = ModelUtils.get_r_j(modelname_exp[1:2]..., lmc=lmc)
	ModelUtils.plot_r_jacobi_arrow!(r_j, y_low)

	R_h = ModelUtils.get_R_h(galaxy)
	if isnothing(y_low_R_h)
		y_low_R_h = limits[3]
	end
	ModelUtils.plot_R_h_arrow!(R_h, y_low_R_h)
	
	ax
end

# ╔═╡ 2be156ef-6a65-473a-b89e-b9fa8d0e83e6
@savefig "density_sigma_i_f_2exp" let
	fig = Figure(size=(3.5, 4.5) .* 72)


	ax_scl = compare_density(fig[1,1], modelnames["scl_smallperi"], modelname_plummer=modelnames["scl_smallperi_2exp"], y_low=-4)	
	hidexdecorations!(ax_scl, ticks=false, minorticks=false)

	ax_umi = compare_density(fig[1, 2], modelnames["umi_smallperi"], modelname_plummer=modelnames["umi_smallperi_2exp"], y_low=-4)
	hidedecorations!(ax_umi, ticks=false, minorticks=false)

	ax_scl_sigma = plot_σv_prof(fig[2,1], modelnames["scl_smallperi_2exp"]..., color=COLORS[5])

	plot_σv_prof!(modelnames["scl_smallperi_2exp"]..., color=COLORS[5], filename="initial.fits", linestyle=:dot)

	plot_σv_prof!(modelnames["scl_smallperi"]..., color=COLORS[1], linewidth=smalllinewidth)

	plot_σv_prof!(modelnames["scl_smallperi"]..., color=COLORS[1], linewidth=smalllinewidth, filename="initial.fits", linestyle=:dot)



	ax_umi_sigma = plot_σv_prof(fig[2,2], modelnames["umi_smallperi_2exp"]..., color=COLORS[5])

	plot_σv_prof!(modelnames["umi_smallperi_2exp"]..., color=COLORS[5], filename="initial.fits", linestyle=:dot)

	plot_σv_prof!(modelnames["umi_smallperi"]..., color=COLORS[1], linewidth=smalllinewidth)

	plot_σv_prof!(modelnames["umi_smallperi"]..., color=COLORS[1], linewidth=smalllinewidth, filename="initial.fits", linestyle=:dot)


	hideydecorations!(ax_umi_sigma, ticks=false, minorticks=false)


	linkaxes!(ax_scl, ax_umi)
	linkxaxes!(ax_scl, ax_scl_sigma)
	linkxaxes!(ax_umi, ax_umi_sigma)

	Legend(fig[3, :], ax_scl, position=:lb, tellwidth=false, nbanks=3)

	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)


	rowsize!(fig.layout, 1, Aspect(1, 1))
	rowsize!(fig.layout, 2, Aspect(1, 1))

	resize_to_layout!()
	fig

end

# ╔═╡ e8edd587-dcd4-4943-9149-a94b567329be
scl_stars = get_stars_final(modelnames["scl_smallperi"]..., "initial.fits")

# ╔═╡ 8814b79f-d296-4750-9d46-98971b654b84
outer_stars = scl_stars[scl_stars.r_ell .> 120, :]

# ╔═╡ 5be6fd0a-ad4e-44ef-8532-6e371fcd3550
import StatsBase: sample

# ╔═╡ a9299df3-5bbd-4387-802d-779c8e1ca67f
stars_sampled = outer_stars[sample(1:size(outer_stars, 1), weights(outer_stars.weights), 100000), :]

# ╔═╡ e3e5c02d-eca3-4a57-97a8-c64848a6272c
let
	fig = Figure()
	ax = Axis(fig[1,1])

	# hist2d!(log10.(stars_sampled.r_ell), stars_sampled.radial_velocity, bins=100)

	N_bins = 15
	r_bins = LinRange(-0.5, 2.5, N_bins+1)


	for i in 1:N_bins
		filt = scl_stars.weights .> 0
		filt .&= scl_stars.r_ell .> 10 .^ r_bins[i]
		filt .&= scl_stars.r_ell .< 10 .^ r_bins[i+1]

		if sum(filt) .< 2
			continue

		end
		h = Arya.DensityEstimators.histogram(scl_stars.radial_velocity[filt], weights=scl_stars.weights[filt], errors=:weighted, normalization=:pdf)

		lines!(midpoints(h.bins), (h.values) ./ maximum(h.values) * diff(r_bins)[1] .+ (r_bins[i]), color=:black)
		hlines!(r_bins[i], color=:black)

	end

	fig
end

# ╔═╡ 6732e1e2-e5af-4fea-acfd-b2d47a6967e4
minimum(scl_stars.weights)

# ╔═╡ 3576d271-3807-4162-a90f-26584090ebad
prof = calc_sigma_v_prof(scl_stars)

# ╔═╡ 61241dc9-06a2-4830-94ec-332395c064b5
scatter(midpoints(prof[1]), prof[2])

# ╔═╡ 33018f0c-9776-4030-b1c1-f2e446ed88bc
scatter(midpoints(prof[1]), log10.(prof[4]))

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
# ╠═dd10823e-cf9a-4c5f-9cab-8d8f30effa33
# ╠═5acf2849-d5bf-49df-94fe-7c3fe738f921
# ╠═ed028d05-86b4-4752-9462-2a837fe11e55
# ╠═9792d244-c295-479f-bbc5-157c5a07c300
# ╠═6368512b-30ed-4457-b419-69dfef2721c7
# ╠═be52a89e-0709-4c93-83f2-c2a3f177274d
# ╠═0cc4c7b4-d707-4f3d-9a82-1b6eb4119143
# ╠═908d84f6-18b2-4f33-abb4-537cd2fdea52
# ╠═be9267a9-eebc-4edf-bde4-c224341f0ff5
# ╟─5c56c73e-d8db-4894-b0cf-f7f1b0e01a72
# ╠═9c3caa8f-cb12-4b40-b091-9e570af0f074
# ╠═2c5eca31-bea1-4805-b869-877ec54828ca
# ╠═14d85541-d6aa-4d00-8c37-f84126cb89b3
# ╠═a212281f-2789-441b-9d0b-5712fef3c8ec
# ╠═405c74a8-075d-48ab-9ae8-827779952a0b
# ╠═26b99d47-5024-4619-b62f-e69ab24361fd
# ╠═4abe36c0-d10e-4d9b-8786-8a7f7aea593a
# ╠═b9ba3889-b77a-406d-b488-2c7bd7db00e0
# ╠═2be156ef-6a65-473a-b89e-b9fa8d0e83e6
# ╠═e8edd587-dcd4-4943-9149-a94b567329be
# ╠═8814b79f-d296-4750-9d46-98971b654b84
# ╠═5be6fd0a-ad4e-44ef-8532-6e371fcd3550
# ╠═a9299df3-5bbd-4387-802d-779c8e1ca67f
# ╠═e3e5c02d-eca3-4a57-97a8-c64848a6272c
# ╠═6732e1e2-e5af-4fea-acfd-b2d47a6967e4
# ╠═3576d271-3807-4162-a90f-26584090ebad
# ╠═61241dc9-06a2-4830-94ec-332395c064b5
# ╠═33018f0c-9776-4030-b1c1-f2e446ed88bc
