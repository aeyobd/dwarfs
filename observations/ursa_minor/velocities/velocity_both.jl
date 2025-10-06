### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 09d570a4-56f1-4ff9-990d-c02534f7351e
using OrderedCollections

# ╔═╡ aea6dd0e-930c-4ccf-9fd4-276c50c2d861
rv_file = "rv_combined_x_2c_psat_0.2_bin.fits"

# ╔═╡ ff6df679-7ec8-4822-a6ad-3d5faf590765
n_samples = 10000

# ╔═╡ f282dbef-9276-4c1c-af3f-366d8b51aa38
n_threads = 16

# ╔═╡ 729206f5-72fb-4592-bc6d-92f39d9ca305
sampler = NUTS(0.65)

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ b5336da2-2d1b-4640-a674-8a3f4d40b78b
md"""
# Pacckages
"""

# ╔═╡ 964ab896-5d9d-42ed-82c6-d7790ce9c871
import CSV

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

# ╔═╡ ed55f077-6e5b-439a-afa4-83946ff5e401
FIGSUFFIX  = "." * splitext(basename(rv_file))[1]

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
memb_stars = let
	rv_meas = read_fits(joinpath(data_dir, rv_file))

	# filter!(row->!ismissing(row.pmra), rv_meas)
	RVUtils.add_gsr!(rv_meas, 
					 distance=obs_properties["distance"],
					pmra=obs_properties["pmra"], pmdec=obs_properties["pmdec"])
	
	rv_meas
end

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ b3ba1ba1-be90-4f13-b9f5-5ee51c55c299
md"""
# Both model
"""

# ╔═╡ df91e538-b283-470a-9a54-167968eac287
R_h = obs_properties["R_h"]

# ╔═╡ 7d958afa-df2f-4c06-b537-2a4121a343c6
model_both = RVUtils.model_vel_gradient_both(memb_stars.vz, memb_stars.vz_err, memb_stars.xi, memb_stars.eta, memb_stars.R_ell, R_h=R_h)

# ╔═╡ fa545320-7119-4fe3-82a9-32a6b7f8bc5a
samples_both = sample(model_both, sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ 22455687-817b-41bb-8a89-ae1ebb145e88
df_both = let
	df = DataFrame(samples_both)
	df[:, :A]
	df[:, :B]
	df[!, :r_grad] = @. 60 * (df.A ⊕ df.B )
	df[!, :Θ_grad] = @. atand(df.A, df.B) 

	df
end

# ╔═╡ 34844c3e-fae5-4163-9f13-d41738116858
@savefig "both_corner" pairplot(samples_both)

# ╔═╡ 9ffd5532-3ce2-4bfc-87f4-bf9501690e44
summary_both = RVUtils.summarize(samples_both)

# ╔═╡ 363081fa-20aa-43dd-a737-356504858e0d
bf_gradient_both = RVUtils.bayes_evidence(model_both, df_both, ["A", "B"])

# ╔═╡ abf72545-93d3-4e9e-b8cb-c4fc2aa1e5d8
bf_Rell_both = RVUtils.bayes_evidence(model_both, df_both, "dlσ_dlR")

# ╔═╡ d717656b-9b53-4bed-b433-477dead6ec09
CSV.write("processed/mcmc_samples_both$FIGSUFFIX.csv", df_both)

# ╔═╡ bfd8425e-5cec-4f2a-a18e-16be3faadc07
θs_both = mod1.(df_both.Θ_grad, 360.) .- 360

# ╔═╡ 4bccb365-8019-4659-a580-374c21112959
θ_m = median(θs_both)

# ╔═╡ a394f7b0-d443-4a2f-ba41-864ca559ac66
θ_err = quantile(θs_both, [0.16, 0.5, 0.84]) .- θ_m

# ╔═╡ 0ea60392-33f4-480f-9002-7dad289b96b5
r_grad_m = median(df_both.r_grad)

# ╔═╡ bb11d2b2-72a7-4e37-89bf-092b20dc0805
r_grad_err = quantile(df_both.r_grad, [0.16, 0.5, 0.84]) .- r_grad_m

# ╔═╡ d5615552-caf8-4c0c-a17c-502c0f8198dc
md"""
# Summaries
"""

# ╔═╡ 0dbb27cf-8a8c-4521-bda4-5768d8a02176
function OrderedCollections.OrderedDict(summary_vz::DataFrame)
	
	df =  OrderedDict(string(col) => summary_vz[!, col] for col in names(summary_vz))

	for key in keys(df)
		if eltype(df[key]) == Symbol
			df[key] = string.(df[key])
		end
	end

	df
end

# ╔═╡ 010b6aa7-e3d0-4441-aac5-6ab87c053e33
df_summaries = OrderedDict(
	"summary" => summary_both |> OrderedDict, 
	"bf_gradient" => bf_gradient_both,
	"bf_rell" => bf_Rell_both,
	"Nmemb" => length(memb_stars.RV), 
	"median_err" => median(memb_stars.RV_err),
	"R_grad_median" => r_grad_m,
	"R_grad_el" => r_grad_err[1],
	"R_grad_ep" => r_grad_err[end],
	"theta_grad_median" => θ_m,
	"theta_grad_el" => -θ_err[1],
	"theta_grad_ep" => θ_err[3]
)

# ╔═╡ c6edc3b5-accc-44c6-afa5-d4b7ed17fd65
open("processed/mcmc_properties_both$FIGSUFFIX.toml", "w") do f
	TOML.print(f, df_summaries)
end

# ╔═╡ Cell order:
# ╠═aea6dd0e-930c-4ccf-9fd4-276c50c2d861
# ╠═ff6df679-7ec8-4822-a6ad-3d5faf590765
# ╠═f282dbef-9276-4c1c-af3f-366d8b51aa38
# ╠═729206f5-72fb-4592-bc6d-92f39d9ca305
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╟─b5336da2-2d1b-4640-a674-8a3f4d40b78b
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═964ab896-5d9d-42ed-82c6-d7790ce9c871
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
# ╠═ed55f077-6e5b-439a-afa4-83946ff5e401
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╟─b3ba1ba1-be90-4f13-b9f5-5ee51c55c299
# ╠═df91e538-b283-470a-9a54-167968eac287
# ╠═7d958afa-df2f-4c06-b537-2a4121a343c6
# ╠═fa545320-7119-4fe3-82a9-32a6b7f8bc5a
# ╠═22455687-817b-41bb-8a89-ae1ebb145e88
# ╠═34844c3e-fae5-4163-9f13-d41738116858
# ╠═9ffd5532-3ce2-4bfc-87f4-bf9501690e44
# ╠═363081fa-20aa-43dd-a737-356504858e0d
# ╠═abf72545-93d3-4e9e-b8cb-c4fc2aa1e5d8
# ╠═d717656b-9b53-4bed-b433-477dead6ec09
# ╠═bfd8425e-5cec-4f2a-a18e-16be3faadc07
# ╠═4bccb365-8019-4659-a580-374c21112959
# ╠═a394f7b0-d443-4a2f-ba41-864ca559ac66
# ╠═0ea60392-33f4-480f-9002-7dad289b96b5
# ╠═bb11d2b2-72a7-4e37-89bf-092b20dc0805
# ╟─d5615552-caf8-4c0c-a17c-502c0f8198dc
# ╠═09d570a4-56f1-4ff9-990d-c02534f7351e
# ╠═010b6aa7-e3d0-4441-aac5-6ab87c053e33
# ╠═0dbb27cf-8a8c-4521-bda4-5768d8a02176
# ╠═c6edc3b5-accc-44c6-afa5-d4b7ed17fd65
