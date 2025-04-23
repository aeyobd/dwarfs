### A Pluto.jl notebook ###
# v0.20.5

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
stars = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_members_all.fits")

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
samples = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_vz.csv", DataFrame)

# ╔═╡ 0a22b941-0651-405e-9624-e9c344b7d93c
v0 = median(samples.μ)

# ╔═╡ 06746433-dea4-485a-b5cd-e0ef2055c314
σ = median(samples.σ)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "scl_vz_hist" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$v_z$ / km\,s$^{-1}$",
		ylabel = "density"
	)

	bins = LinRange(v0-5σ, v0+5σ, 40)
	_, values, errs = LilGuys.histogram(stars.vz, bins, normalization=:pdf)

	RVUtils.plot_samples!(samples, LinRange(v0-5σ, v0+5σ, 100), thin=100, )
	errorscatter!(midpoints(bins), values, yerror=errs, color=COLORS[2], size=5, linewidth=2)
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
# ╠═0a22b941-0651-405e-9624-e9c344b7d93c
# ╠═06746433-dea4-485a-b5cd-e0ef2055c314
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
