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

# ╔═╡ 314da249-8a3f-47e9-ac40-f86406ec9955
using CSV, DataFrames

# ╔═╡ 71335ca6-b58b-4c68-8467-81cd4f481b1f
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ ad8f5ae4-5540-4f1a-8ca0-7f2ddd96fe03
import TOML

# ╔═╡ 7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
import StatsBase: median

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ f27c2500-b92c-4de0-8482-4a8acfa2667f
module RVUtils
	include(joinpath(ENV["DWARFS_ROOT"], "observations/rv_utils.jl"))
end

# ╔═╡ 62d634a7-d1f2-4e49-b478-1521a8e172e3
obs_props_scl = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ 0555866b-6bce-44d2-b860-60a3318474ca
obs_props_umi = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/observed_properties.toml"))

# ╔═╡ 3b5e60d4-1ac4-4b8d-8284-fb1c53e966c6
stars = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_combined_x_wide_2c_psat_0.2.fits")

# ╔═╡ f9f00f0c-0924-4abe-8c16-a58d46c76faa
stars_umi = read_fits(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/rv_combined_x_2c_psat_0.2.fits")

# ╔═╡ 91f9f767-7d35-4a99-9493-d9bb3ef6e28c
function make_label(galaxy, samples; digits=1, digits_sigma=digits)
	μ = median(samples.μ)
	μ_err = round.(LilGuys.quantile(samples.μ, [0.16, 0.84]) .- μ, digits=digits)
	μ = round(μ, digits=digits)

	

	σ = median(samples.σ)
	σ_err = round.(LilGuys.quantile(samples.σ, [0.16, 0.84]) .- σ, digits=digits_sigma)

	σ = round(σ, digits=digits_sigma)

	s = L"""
%$galaxy\\
$\mu = %$(μ)_{%$(μ_err[1])}^{+%$(μ_err[2])}$\\
$\sigma = %$(σ)_{%$(σ_err[1])}^{+%$(σ_err[2])}$
"""

	return s
end

# ╔═╡ d98d9fcf-feee-4274-b90c-fa22b6b171da
dv_scl = RVUtils.rv_gsr_shift(obs_props_scl["ra"], obs_props_scl["dec"])

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
samples = let
	df = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_both.rv_combined_x_wide_2c_psat_0.2.csv", DataFrame)
	df.μ .+= dv_scl
	df
end

# ╔═╡ 0a22b941-0651-405e-9624-e9c344b7d93c
v0 = median(samples.μ)

# ╔═╡ 06746433-dea4-485a-b5cd-e0ef2055c314
σ = median(samples.σ)

# ╔═╡ 43829ff1-e0be-476f-8ce8-2eff6eaf51b4
v0-3σ,v0+3σ

# ╔═╡ 63a2a5f3-133b-4bde-b6ea-7c3e76be761b
(extrema(stars.vz[stars.vz_err .< 2]) .- v0) ./ σ

# ╔═╡ e10a7697-58f8-4a65-bd7c-47b0771b811d
make_label("Sculptor", samples, digits_sigma=1)

# ╔═╡ 4f2d821d-c4ff-44fd-b494-cc74050049fe
dv_umi = RVUtils.rv_gsr_shift(obs_props_umi["ra"], obs_props_umi["dec"])

# ╔═╡ 21bdc277-d3b0-48c6-984b-538bda4348a4
samples_umi = let
	df = CSV.read(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/mcmc_samples_both.rv_combined_x_2c_psat_0.2.csv", DataFrame)
	df.μ .+= dv_umi
	df

end

# ╔═╡ 0ff4cce0-601b-4532-8c75-1e8961631fa1
v0_umi = median(samples_umi.μ)

# ╔═╡ f7999d84-cfe7-4803-b135-c6c7c7050066
σ_umi = median(samples_umi.σ)

# ╔═╡ cb1ff71e-48df-41b3-9da1-8b29c9118e05
(extrema(stars_umi.vz[stars_umi.vz_err .< 2])  .- v0_umi )./ σ_umi

# ╔═╡ 1acafc19-e86a-4f18-89da-833ae16664e9
smallfontsize = theme(:fontsize)[]*0.7

# ╔═╡ 9e5c76bd-ca64-4824-afb3-b1ea898f3b13
df_gradient_umi = CSV.read(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/vz_r_ell_binned.rv_combined_x_2c_psat_0.2.csv", DataFrame)

# ╔═╡ 26b3a353-eaf4-43c3-b55e-b3c80dbc2488
df_gradient_scl = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/vz_r_ell_binned.rv_combined_x_wide_2c_psat_0.2.csv", DataFrame)

# ╔═╡ 21a633a6-4e53-4ff7-a115-537742f8bdb2
function plot_obs_sigma!(df)
	x_mid = log10.(df.x)
	log_bin_errs = [ (x_mid[i] .- log10.(df.x_low[i]), log10(df.x_high[i]) .- x_mid[i]) for i in eachindex(x_mid)]
	y_err = collect(zip(df.σ_em ./ df.σ ./ log(10), df.σ_ep ./ df.σ ./ log(10)))

	errorscatter!(x_mid, log10.(df.σ), xerror=log_bin_errs, yerror=y_err, color=:black)
end

# ╔═╡ a7ef05e7-52ca-45b6-a3f2-7cf451d532bb
function plot_samples_sigma!(df_Rell; skip=800)
	x = LinRange(-1, 3, 100)

	
	for i in 1:skip:size(df_Rell, 1)
		y = df_Rell.σ[i] .* 10 .^ ((x .- 1.0) .* df_Rell.dlσ_dlR[i])
		lines!(x, log10.(y), alpha=0.03, color=COLORS[3], linestyle=:solid)
	end

	hlines!(log10(median(df_Rell.σ)), color=:black, linewidth=theme(:linewidth)[]/2)
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "scl_umi_rv_fits" let
	fig = Figure()
	v0 = median(samples.μ)
	σ = median(samples.σ)
	ax = Axis(fig[1,1], 
		ylabel = "density",
			  limits=(v0-4σ, v0+4σ, 0, nothing)
	)



	bins = LinRange(v0-4σ, v0+4σ, 40)
	_, values, errs = LilGuys.histogram(stars.vz .+ dv_scl, bins, normalization=:pdf)

	RVUtils.plot_samples!(samples, LinRange(v0-4σ, v0+4σ, 100), thin=800, color=COLORS[3], linestyle=:solid)
	errorscatter!(midpoints(bins), values, yerror=errs, color=:black)

	text!(0.05, 0.9, text=make_label("Sculptor", samples, digits_sigma=2), space=:relative, align=(:left, :top), justification=:left, fontsize=smallfontsize)



	v0 = median(samples_umi.μ)
	σ = median(samples_umi.σ)
	ax = Axis(fig[2,1], 
		xlabel = L"$\textrm{v}_\textrm{gsr}'$ / km\,s$^{-1}$",
		ylabel = "density",
		limits=(v0-4σ, v0+4σ, 0, nothing)

	)


	bins = LinRange(v0-5σ, v0+5σ, 40)
	_, values, errs = LilGuys.histogram(stars_umi.vz .+ dv_umi, bins, normalization=:pdf)

	RVUtils.plot_samples!(samples_umi, LinRange(v0-4σ, v0+4σ, 100), thin=800, color=COLORS[3], linestyle=:solid)
	errorscatter!(midpoints(bins), values, yerror=errs, color=:black)

	text!(0.05, 0.9, text=make_label("Ursa Minor", samples_umi), space=:relative, align=(:left, :top), justification=:left, fontsize=smallfontsize)




	# velocity gradients
	ax_sigma = Axis(fig[1, 2],
		#xlabel = L"$\log\ R$ / arcmin",
		ylabel = L"$\log\ \sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	text!(0, 1, offset=(6, -6), text="Sculptor", space=:relative, align=(:left, :top), fontsize=smallfontsize)

	plot_samples_sigma!(samples)

	plot_obs_sigma!(df_gradient_scl)	

	xlims!(-1, 2.5)
	ylims!(0.8, 1.1)
	annotation!(0, 32, log10(obs_props_scl["R_h"]), 0.8)
	text!(log10(obs_props_scl["R_h"]), 0.8, text=L"R_h", offset=(6, 18), fontsize=smallfontsize)


	ax2 = Axis(fig[2,2],
		xlabel = L"$\log\ R$ / arcmin",
		ylabel = L"$\log\ \sigma_{v, \textrm{los}}$ / km s$^{-1}$"

			  )

	plot_samples_sigma!(samples_umi)

	plot_obs_sigma!(df_gradient_umi)	
	text!(0, 1, offset=(6, -6), text="Ursa Minor", space=:relative, align=(:left, :top), fontsize=smallfontsize)

	annotation!(0, 32, log10(obs_props_umi["R_h"]), 0.6)
	text!(log10(obs_props_umi["R_h"]), 0.6, text=L"R_h", offset=(6, 18), fontsize=smallfontsize)
	ylims!(0.6, 1.2)
	xlims!(-1, 2.5)

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═ad8f5ae4-5540-4f1a-8ca0-7f2ddd96fe03
# ╠═314da249-8a3f-47e9-ac40-f86406ec9955
# ╠═71335ca6-b58b-4c68-8467-81cd4f481b1f
# ╠═7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f27c2500-b92c-4de0-8482-4a8acfa2667f
# ╠═62d634a7-d1f2-4e49-b478-1521a8e172e3
# ╠═0555866b-6bce-44d2-b860-60a3318474ca
# ╠═3b5e60d4-1ac4-4b8d-8284-fb1c53e966c6
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═f9f00f0c-0924-4abe-8c16-a58d46c76faa
# ╠═21bdc277-d3b0-48c6-984b-538bda4348a4
# ╠═0a22b941-0651-405e-9624-e9c344b7d93c
# ╠═06746433-dea4-485a-b5cd-e0ef2055c314
# ╠═43829ff1-e0be-476f-8ce8-2eff6eaf51b4
# ╠═63a2a5f3-133b-4bde-b6ea-7c3e76be761b
# ╠═0ff4cce0-601b-4532-8c75-1e8961631fa1
# ╠═f7999d84-cfe7-4803-b135-c6c7c7050066
# ╠═cb1ff71e-48df-41b3-9da1-8b29c9118e05
# ╠═91f9f767-7d35-4a99-9493-d9bb3ef6e28c
# ╠═e10a7697-58f8-4a65-bd7c-47b0771b811d
# ╠═d98d9fcf-feee-4274-b90c-fa22b6b171da
# ╠═4f2d821d-c4ff-44fd-b494-cc74050049fe
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═1acafc19-e86a-4f18-89da-833ae16664e9
# ╠═9e5c76bd-ca64-4824-afb3-b1ea898f3b13
# ╠═26b3a353-eaf4-43c3-b55e-b3c80dbc2488
# ╠═21a633a6-4e53-4ff7-a115-537742f8bdb2
# ╠═a7ef05e7-52ca-45b6-a3f2-7cf451d532bb
