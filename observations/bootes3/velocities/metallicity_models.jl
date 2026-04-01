### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 93838644-cad6-4df3-b554-208b7afeb3b8
using PyFITS

# ╔═╡ 72f1febc-c6ea-449a-8cec-cd0e49c4e20c
using DataFrames

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
using PairPlots

# ╔═╡ bd6dfd17-02ee-4855-be37-fecfdab6776f
using LilGuys; FIGDIR = "figures"

# ╔═╡ 3ed8c28f-5908-42dc-a56b-24a9b2685a07
md"""
## Inputs
"""

# ╔═╡ 428c3d92-e49c-426e-b328-2d2a8d4c4159
rv_file = "rv_deimos_x_2c_psat_0.2.fits"

# ╔═╡ 8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
n_samples = 1_000

# ╔═╡ 0b8f0193-f1e5-471e-8c34-498ac827be63
n_threads = 11

# ╔═╡ ad79710a-0392-4b75-94ac-23881376a70b
sampler = NUTS(0.65)

# ╔═╡ 680e7f76-cb4d-40d6-9a9f-d4672427a633
md"""
## derived
"""

# ╔═╡ 86fe351f-ef12-474a-85cc-c10c22a65e77
FIGSUFFIX  = "." * splitext(basename(rv_file))[1]

# ╔═╡ 7330c75e-1bf9-476a-8274-ebc86d555e6f
md"""
# RV sample models
"""

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ 34e43f4a-bcff-41cb-92c4-0c8d600fd053
import CSV

# ╔═╡ 2a422e88-fc0d-4a89-a841-42f3c5c8dace
import KernelDensity

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 2ab0018f-1628-4f5e-b7da-370eb20c00d0
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ 5ec475a1-14bb-40f6-856a-69fa9efe087a
⊕ = RVUtils.:⊕

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/bootes3/observed_properties.toml")

# ╔═╡ 8b0d5ad4-2858-47dd-93c4-224fa9e016c8
md"""
# Data processing
"""

# ╔═╡ eac815ec-680c-42f9-aa40-fa8e87daf2b4
rv_meas = read_fits("processed/$rv_file")

# ╔═╡ cb511bb7-ee81-4c74-b07e-1d1b268305d8
memb_stars = rv_meas[(rv_meas.PSAT_RV .> 0.5) .& (rv_meas.FE_H .> -999), :]

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership
"""

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.FE_H)

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f


# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ f22d0148-d804-4611-a30b-baee561b20b7
@model function metallicity_model(fe_h, fe_h_err)
	μ ~ Uniform(-4, 0)
	σ ~ Uniform(0, 2)

	m = fill(μ, length(fe_h))
	s = @. sqrt(σ^2 + fe_h_err^2)
	fe_h ~ MvNormal(m, s)
end

# ╔═╡ 10a90125-a1e4-44c9-afa0-843d46455ea2
model = metallicity_model(memb_stars.FE_H, memb_stars.FE_H_err)

# ╔═╡ 9ffcc0b0-11a1-4962-b64e-59a731f22bd8
samples = sample(model, sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ f8904d41-c448-417c-bbe9-6626a6283b67
df_out = DataFrame(samples)

# ╔═╡ 9069d73b-2086-448d-8387-f8a8b8d8dcf9
median(rv_meas.FE_H)

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
pairplot(samples)

# ╔═╡ 1f251e67-7fed-4375-8bed-697704302e43
summary = RVUtils.summarize(samples)

# ╔═╡ 6dd003db-bd29-4872-94f2-445f5c4536ba

let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(-3, -2, -1.1, 5))
	RVUtils.plot_samples!(df_out, LinRange(-3, -2, 100), thin=100)

	scatter!(memb_stars.FE_H, fill(-1, length(memb_stars.FE_H)))
	
	fig
end

# ╔═╡ c1cf3a80-1c26-4417-9ddf-1046f3f5f608
CSV.write("processed/mcmc_samples_feh$FIGSUFFIX.csv", df_out)

# ╔═╡ 1acdc61b-fb5f-449d-ba86-46525c881a39
md"""
# Writing Information
"""

# ╔═╡ Cell order:
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═428c3d92-e49c-426e-b328-2d2a8d4c4159
# ╠═8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
# ╠═0b8f0193-f1e5-471e-8c34-498ac827be63
# ╠═ad79710a-0392-4b75-94ac-23881376a70b
# ╠═680e7f76-cb4d-40d6-9a9f-d4672427a633
# ╠═86fe351f-ef12-474a-85cc-c10c22a65e77
# ╟─7330c75e-1bf9-476a-8274-ebc86d555e6f
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═34e43f4a-bcff-41cb-92c4-0c8d600fd053
# ╠═2a422e88-fc0d-4a89-a841-42f3c5c8dace
# ╠═93838644-cad6-4df3-b554-208b7afeb3b8
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═bd6dfd17-02ee-4855-be37-fecfdab6776f
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╠═5ec475a1-14bb-40f6-856a-69fa9efe087a
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╟─8b0d5ad4-2858-47dd-93c4-224fa9e016c8
# ╠═eac815ec-680c-42f9-aa40-fa8e87daf2b4
# ╠═cb511bb7-ee81-4c74-b07e-1d1b268305d8
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╟─6734991c-16c0-4424-a2bb-84bfa811121f
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═f22d0148-d804-4611-a30b-baee561b20b7
# ╠═10a90125-a1e4-44c9-afa0-843d46455ea2
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═f8904d41-c448-417c-bbe9-6626a6283b67
# ╠═9069d73b-2086-448d-8387-f8a8b8d8dcf9
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═1f251e67-7fed-4375-8bed-697704302e43
# ╠═6dd003db-bd29-4872-94f2-445f5c4536ba
# ╠═c1cf3a80-1c26-4417-9ddf-1046f3f5f608
# ╟─1acdc61b-fb5f-449d-ba86-46525c881a39
