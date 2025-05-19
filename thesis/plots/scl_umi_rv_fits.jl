### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
import StatsBase: median

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ f27c2500-b92c-4de0-8482-4a8acfa2667f
module RVUtils
	include(joinpath(ENV["DWARFS_ROOT"], "observations/rv_utils.jl"))
end

# ╔═╡ 3b5e60d4-1ac4-4b8d-8284-fb1c53e966c6
stars = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_combined_x_wide_2c_psat_0.2.fits")

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
samples = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_vz.rv_combined_x_wide_2c_psat_0.2.csv", DataFrame)

# ╔═╡ 6d558d6c-60b6-4361-b875-c821078c0675


# ╔═╡ f9f00f0c-0924-4abe-8c16-a58d46c76faa
stars_umi = read_fits(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/rv_combined_x_2c_psat_0.2.fits")

# ╔═╡ 21bdc277-d3b0-48c6-984b-538bda4348a4
samples_umi = CSV.read(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/mcmc_samples_vz.rv_combined_x_2c_psat_0.2.csv", DataFrame)

# ╔═╡ 0a22b941-0651-405e-9624-e9c344b7d93c
v0 = median(samples.μ)

# ╔═╡ 06746433-dea4-485a-b5cd-e0ef2055c314
σ = median(samples.σ)

# ╔═╡ 43829ff1-e0be-476f-8ce8-2eff6eaf51b4
v0-3σ,v0+3σ

# ╔═╡ 63a2a5f3-133b-4bde-b6ea-7c3e76be761b
(extrema(stars.vz[stars.vz_err .< 2]) .- v0) ./ σ

# ╔═╡ 0ff4cce0-601b-4532-8c75-1e8961631fa1
v0_umi = median(samples_umi.μ)

# ╔═╡ f7999d84-cfe7-4803-b135-c6c7c7050066
σ_umi = median(samples_umi.σ)

# ╔═╡ cb1ff71e-48df-41b3-9da1-8b29c9118e05
(extrema(stars_umi.vz[stars_umi.vz_err .< 2])  .- v0_umi )./ σ_umi

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

# ╔═╡ e10a7697-58f8-4a65-bd7c-47b0771b811d
make_label("Sculptor", samples, digits_sigma=1)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "scl_umi_rv_fits" let
	fig = Figure()
	v0 = median(samples.μ)
	σ = median(samples.σ)
	ax = Axis(fig[1,1], 
		xlabel = L"$v_\textrm{gsr}'$ / km\,s$^{-1}$",
		ylabel = "density",
			  limits=(v0-4σ, v0+4σ, 0, nothing)
	)



	bins = LinRange(v0-4σ, v0+4σ, 40)
	_, values, errs = LilGuys.histogram(stars.vz, bins, normalization=:pdf)

	RVUtils.plot_samples!(samples, LinRange(v0-4σ, v0+4σ, 100), thin=800, color=COLORS[3])
	errorscatter!(midpoints(bins), values, yerror=errs, color=:black)

	text!(0.05, 0.9, text=make_label("Sculptor", samples, digits_sigma=2), space=:relative, align=(:left, :top), justification=:left)



	v0 = median(samples_umi.μ)
	σ = median(samples_umi.σ)
	ax = Axis(fig[2,1], 
		xlabel = L"$v_\textrm{gsr}'$ / km\,s$^{-1}$",
		ylabel = "density",
		limits=(v0-4σ, v0+4σ, 0, nothing)

	)


	bins = LinRange(v0-5σ, v0+5σ, 40)
	_, values, errs = LilGuys.histogram(stars_umi.vz, bins, normalization=:pdf)

	RVUtils.plot_samples!(samples_umi, LinRange(v0-4σ, v0+4σ, 100), thin=800, color=COLORS[3])
	errorscatter!(midpoints(bins), values, yerror=errs, color=:black)

	text!(0.05, 0.9, text=make_label("Ursa Minor", samples_umi), space=:relative, align=(:left, :top), justification=:left)

	
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═314da249-8a3f-47e9-ac40-f86406ec9955
# ╠═71335ca6-b58b-4c68-8467-81cd4f481b1f
# ╠═7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f27c2500-b92c-4de0-8482-4a8acfa2667f
# ╠═3b5e60d4-1ac4-4b8d-8284-fb1c53e966c6
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═6d558d6c-60b6-4361-b875-c821078c0675
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
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
