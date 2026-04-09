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

# ╔═╡ 0532010e-7832-44d4-a7ee-5a6f6ee9d7da
using OrderedCollections

# ╔═╡ 3ed8c28f-5908-42dc-a56b-24a9b2685a07
md"""
## Inputs
"""

# ╔═╡ 428c3d92-e49c-426e-b328-2d2a8d4c4159
rv_file = let
	"rv_deimos_geha_x_2c.fits"
	# "rv_deimos_x_2c_sigma_3.fits"
	# "rv_deimos_x_2c_psat_0.2.fits"
end

# ╔═╡ 8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
n_samples = 10_000

# ╔═╡ dead238a-f726-4c26-9402-c42797f08f06
n_samples_extra = 300

# ╔═╡ 0b8f0193-f1e5-471e-8c34-498ac827be63
n_threads = 16

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

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
σv = obs_properties["sigma_v"]

# ╔═╡ 00cd5ad7-ab0e-428f-ba00-05ef6fc86806
R_h = obs_properties["R_h"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ 8b0d5ad4-2858-47dd-93c4-224fa9e016c8
md"""
# Data processing
"""

# ╔═╡ eac815ec-680c-42f9-aa40-fa8e87daf2b4
rv_meas = read_fits("processed/$rv_file")

# ╔═╡ 236e3900-4c55-49c0-ad5d-eaa06d292c5c
if :PSAT_RV ∈ names(rv_meas)
	@assert all([rv_meas.PSAT_RV .> 0.5, :])
end

# ╔═╡ cb511bb7-ee81-4c74-b07e-1d1b268305d8
memb_stars = rv_meas

# ╔═╡ 82d6afc3-50cc-4e85-924e-dc91ae14f4de
sort(rv_meas, "ra")

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership
"""

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ 55ce0f69-8a96-4bbb-a59f-ee6503624ea6
md"""
# Numbers
"""

# ╔═╡ 3377f632-713d-4fec-84a9-b0211b02cb43
median(memb_stars.RV_err)

# ╔═╡ 31a6c2e4-538c-4adc-bbda-5043680b17f7
extrema(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f
md"""
## MCMC Priors
"""

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err), Prior(), n_samples)
)

# ╔═╡ ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
pairplot(prior_samples[:, [:μ, :σ]])

# ╔═╡ 74ad07df-f15c-436f-b390-ce95b27f7fab
let
	fig, ax = FigAxis(
		xgridvisible=false,
		ygridvisible=false,
		xlabel=xlabel,
		ylabel="density",
		title="priors"
	)
	
	RVUtils.plot_samples!(prior_samples, LinRange(-130, 130, 100))

	fig
end

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 318d29b9-4c84-4d38-af6d-048518952970
samples = DataFrame(sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err, μ_0_prior=Δv_gsr), sampler, MCMCThreads(), n_samples, n_threads))

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)

	errorscatter!(memb_stars.RV, zeros(length(memb_stars.RV)) .+ 0.003*randn(length(memb_stars.vz)), xerror=memb_stars.RV_err,  color=COLORS[2])
	RVUtils.plot_samples!(samples, LinRange(160, 220, 100), thin=150)

	fig
end

# ╔═╡ 61e15c47-c454-48db-95f9-02abe052676e
mean(memb_stars.RV)

# ╔═╡ d5938fc3-9c8a-4e33-8401-500b4201df14
sem(memb_stars.RV)

# ╔═╡ d0ae48e2-8389-4641-b311-cfb4944b0851
std(memb_stars.RV)

# ╔═╡ 49e7d483-1dd1-4406-9bce-58d6e4412e7b
mean(memb_stars.radial_velocity_gsr)

# ╔═╡ 9ffcc0b0-11a1-4962-b64e-59a731f22bd8
samples_gsr = sample(RVUtils.model_vel_1c(memb_stars.vz, memb_stars.vz_err), sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ f8904d41-c448-417c-bbe9-6626a6283b67
df_gsr = DataFrame(samples_gsr)

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
@savefig "rv_sigma_corner" pairplot(samples_gsr)

# ╔═╡ 1f251e67-7fed-4375-8bed-697704302e43
summary_vz = RVUtils.summarize(samples_gsr)

# ╔═╡ a5f1a339-6093-49ea-bd10-b513688d668c
median(df_gsr.μ) + Δv_gsr

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)

	errorscatter!((memb_stars.vz), zeros(length(memb_stars.vz)) .+ 0.003*randn(length(memb_stars.vz)), xerror=memb_stars.vz_err,  color=COLORS[2])

	RVUtils.plot_samples!(df_gsr, LinRange(210, 260, 100), thin=160)

	@savefig "RV_hist_mcmc"
	fig
end

# ╔═╡ 6dd003db-bd29-4872-94f2-445f5c4536ba
stephist(Float64.(memb_stars.vz), bins=9, normalization=:pdf)


# ╔═╡ c1cf3a80-1c26-4417-9ddf-1046f3f5f608
CSV.write("processed/mcmc_samples_vz$FIGSUFFIX.csv", samples_gsr)

# ╔═╡ 1acdc61b-fb5f-449d-ba86-46525c881a39
md"""
# Writing Information
"""

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

# ╔═╡ ba6c51c9-6d8f-4eca-af33-c75d0a5a5b37
df_summaries = OrderedDict(
	"vz" => summary_vz |> OrderedDict, 
	# "gradient" => summary_gradient |> OrderedDict,
	# "rell" => summary_Rell |> OrderedDict,
	# "both" => summary_both |> OrderedDict,
	# "bf_gradient" => bf_gradient,
	# "bf_rell" => bf_sigma_Rell,
	# "bf_gradient_both" => bf_gradient_both,
	# "bf_rell_both" => bf_Rell_both,
	# "R_grad_median" => r_grad_m_both,
	# "R_grad_el" => r_grad_both_err[1],
	# "R_grad_ep" => r_grad_both_err[end],
	# "theta_grad_median" => θ_m_both,
	# "theta_grad_el" => -θ_both_err[1],
	# "theta_grad_ep" => θ_both_err[3]
)

# ╔═╡ d71f6eba-7d64-4212-91c3-707a664c6b0b
open("processed/mcmc_properties_$FIGSUFFIX.toml", "w") do f
	TOML.print(f, df_summaries)
end

# ╔═╡ Cell order:
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═428c3d92-e49c-426e-b328-2d2a8d4c4159
# ╠═8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
# ╠═dead238a-f726-4c26-9402-c42797f08f06
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
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═00cd5ad7-ab0e-428f-ba00-05ef6fc86806
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╟─8b0d5ad4-2858-47dd-93c4-224fa9e016c8
# ╠═eac815ec-680c-42f9-aa40-fa8e87daf2b4
# ╠═236e3900-4c55-49c0-ad5d-eaa06d292c5c
# ╠═cb511bb7-ee81-4c74-b07e-1d1b268305d8
# ╠═82d6afc3-50cc-4e85-924e-dc91ae14f4de
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╟─55ce0f69-8a96-4bbb-a59f-ee6503624ea6
# ╠═3377f632-713d-4fec-84a9-b0211b02cb43
# ╠═31a6c2e4-538c-4adc-bbda-5043680b17f7
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╟─6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═d0ae48e2-8389-4641-b311-cfb4944b0851
# ╠═49e7d483-1dd1-4406-9bce-58d6e4412e7b
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═f8904d41-c448-417c-bbe9-6626a6283b67
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═1f251e67-7fed-4375-8bed-697704302e43
# ╠═a5f1a339-6093-49ea-bd10-b513688d668c
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═6dd003db-bd29-4872-94f2-445f5c4536ba
# ╠═c1cf3a80-1c26-4417-9ddf-1046f3f5f608
# ╟─1acdc61b-fb5f-449d-ba86-46525c881a39
# ╠═0532010e-7832-44d4-a7ee-5a6f6ee9d7da
# ╠═ba6c51c9-6d8f-4eca-af33-c75d0a5a5b37
# ╠═a393eeca-c9a4-412f-ae86-35ee1aca4d51
# ╠═d71f6eba-7d64-4212-91c3-707a664c6b0b
