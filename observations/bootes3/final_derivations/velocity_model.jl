### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	using LilGuys; FIGDIR = "figures"

	using OrderedCollections, DataFrames
	import CSV
	using PyFITS
	import TOML

	import StatsBase: quantile, mean, std, median, sem
	import DensityEstimators: histogram, bins_equal_number

	using Turing, PairPlots
end

# ╔═╡ 3ed8c28f-5908-42dc-a56b-24a9b2685a07
md"""
## Inputs
"""

# ╔═╡ 428c3d92-e49c-426e-b328-2d2a8d4c4159
rv_file = "rv_deimos_geha_x_2c.fits"
	

# ╔═╡ 36a5bc6b-0e87-4464-9c81-d29246b855ab
j24_sample = "2c"

# ╔═╡ 8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
n_samples = 10_000

# ╔═╡ 0b8f0193-f1e5-471e-8c34-498ac827be63
n_threads = 16

# ╔═╡ ad79710a-0392-4b75-94ac-23881376a70b
sampler = NUTS(0.65)

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

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 2ab0018f-1628-4f5e-b7da-370eb20c00d0
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../velocities/processed"

# ╔═╡ a02dfeb7-1289-4f23-970b-35b064274515
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 70e77c50-38e4-4ac7-b2e8-9ad4ba0cb2b6
rv0 = obs_properties["radial_velocity"]

# ╔═╡ 65da9292-661a-4602-ab24-11e2bcbc441b
σ0 = obs_properties["sigma_v"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 6eef05fd-6247-4655-b3cf-f8429b38077c
j24 = RVUtils.read_gaia_stars("../data/j24_$j24_sample.fits", obs_properties)

# ╔═╡ fb8251ae-7689-4f23-acd2-f162e45b7bf2
rv_all = read_fits(joinpath(data_dir, rv_file))

# ╔═╡ e80ee6db-5e67-4e25-b4e6-8e0525c3ad51
rv_meas = RVUtils.xmatch_and_clean(rv_all, j24, obs_properties, require_match=false)

# ╔═╡ 9ad10d28-3811-449d-854d-b8159f0051a0
md"""
## Membership
"""

# ╔═╡ d10065d4-60f1-456d-99d9-73d87a159bd1
memb_filt =  (rv_meas.P_SAT_deimos .> 0.5) .& (rv_meas.Var .<= 0)

# ╔═╡ 198e206c-9358-458f-87d6-dc5c2b8ca6ca
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9ffcc0b0-11a1-4962-b64e-59a731f22bd8
samples_gsr = sample(RVUtils.model_vel_1c(memb_stars.vz, memb_stars.vz_err), sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ f8904d41-c448-417c-bbe9-6626a6283b67
df_gsr = let
	df = DataFrame(samples_gsr)
	df[!, :μ_icrs] = df.μ .+ Δv_gsr
	df

end

# ╔═╡ d551b237-d645-4adf-808d-c56022de2cb7
@savefig "rv_sigma_corner" pairplot(df_gsr[!, [:μ_icrs, :σ]])

# ╔═╡ 1f251e67-7fed-4375-8bed-697704302e43
summary_vz = RVUtils.summarize(samples_gsr)

# ╔═╡ 815d6d49-1f90-422a-9f46-aa2208404615
df_icrs = let
	df = copy(df_gsr)
	df[!, :μ] = df.μ_icrs
	df
end

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	xlims = rv0 .+ (-σ0, σ0) .* 4

	errorscatter!((memb_stars.vz .+ Δv_gsr), zeros(length(memb_stars.vz)) .- 0.01*rand(length(memb_stars.vz)), xerror=memb_stars.vz_err,  color=COLORS[2], label="observations")

	RVUtils.plot_samples!(df_icrs, LinRange(xlims..., 100), thin=160)
	xlims!(xlims...)

	lines!([NaN], [NaN], color=:black, alpha=0.3, label="MCMC samples")

	axislegend()
	@savefig "RV_hist_mcmc"
	fig
end

# ╔═╡ 5501db91-1c97-4b3a-b508-b7c14948acaa
scatter(memb_stars.vz .+ Δv_gsr, memb_stars.RV .- (memb_stars.vz .+ Δv_gsr))

# ╔═╡ 6dd003db-bd29-4872-94f2-445f5c4536ba
stephist(Float64.(memb_stars.vz), bins=9, normalization=:pdf)

# ╔═╡ c1cf3a80-1c26-4417-9ddf-1046f3f5f608
CSV.write("mcmc/mcmc_samples_vz.csv", samples_gsr)

# ╔═╡ 1acdc61b-fb5f-449d-ba86-46525c881a39
md"""
# Writing Information
"""

# ╔═╡ ba6c51c9-6d8f-4eca-af33-c75d0a5a5b37
df_summaries = summary_vz 

# ╔═╡ a393eeca-c9a4-412f-ae86-35ee1aca4d51
function OrderedCollections.OrderedDict(summary_vz::DataFrame)
	
	df =  OrderedDict(string(col) => summary_vz[!, col] for col in names(summary_vz))

	for key in keys(df)
		if eltype(df[key]) == Symbol
			df[key] = string.(df[key])
		end
	end

	df
end

# ╔═╡ d71f6eba-7d64-4212-91c3-707a664c6b0b
open("mcmc/mcmc_properties_vz.toml", "w") do f
	TOML.print(f, OrderedDict(df_gsr))
end

# ╔═╡ Cell order:
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═428c3d92-e49c-426e-b328-2d2a8d4c4159
# ╠═36a5bc6b-0e87-4464-9c81-d29246b855ab
# ╠═8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
# ╠═0b8f0193-f1e5-471e-8c34-498ac827be63
# ╠═ad79710a-0392-4b75-94ac-23881376a70b
# ╠═86fe351f-ef12-474a-85cc-c10c22a65e77
# ╟─7330c75e-1bf9-476a-8274-ebc86d555e6f
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═a02dfeb7-1289-4f23-970b-35b064274515
# ╠═70e77c50-38e4-4ac7-b2e8-9ad4ba0cb2b6
# ╠═65da9292-661a-4602-ab24-11e2bcbc441b
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═6eef05fd-6247-4655-b3cf-f8429b38077c
# ╠═fb8251ae-7689-4f23-acd2-f162e45b7bf2
# ╠═e80ee6db-5e67-4e25-b4e6-8e0525c3ad51
# ╟─9ad10d28-3811-449d-854d-b8159f0051a0
# ╠═d10065d4-60f1-456d-99d9-73d87a159bd1
# ╠═198e206c-9358-458f-87d6-dc5c2b8ca6ca
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═f8904d41-c448-417c-bbe9-6626a6283b67
# ╠═d551b237-d645-4adf-808d-c56022de2cb7
# ╠═1f251e67-7fed-4375-8bed-697704302e43
# ╠═815d6d49-1f90-422a-9f46-aa2208404615
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═5501db91-1c97-4b3a-b508-b7c14948acaa
# ╠═6dd003db-bd29-4872-94f2-445f5c4536ba
# ╠═c1cf3a80-1c26-4417-9ddf-1046f3f5f608
# ╟─1acdc61b-fb5f-449d-ba86-46525c881a39
# ╠═ba6c51c9-6d8f-4eca-af33-c75d0a5a5b37
# ╠═a393eeca-c9a4-412f-ae86-35ee1aca4d51
# ╠═d71f6eba-7d64-4212-91c3-707a664c6b0b
