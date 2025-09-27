### A Pluto.jl notebook ###
# v0.20.18

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
rv_file = "rv_combined_x_2c_psat_0.2.fits"

# ╔═╡ ff6df679-7ec8-4822-a6ad-3d5faf590765
n_samples = 1000

# ╔═╡ f282dbef-9276-4c1c-af3f-366d8b51aa38
n_threads = 4

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

# ╔═╡ 32674b48-0118-414f-9d48-2b1c44c02885
md"""
# Numbers
"""

# ╔═╡ 29425bbc-05b6-434f-b4b3-17ca39ebf830
median(memb_stars.RV_err)

# ╔═╡ d8b97f8b-f9c3-4f3e-942f-e25a4f27fff6
length(memb_stars.RV)

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 318d29b9-4c84-4d38-af6d-048518952970
samples = DataFrame(sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err, μ_0_prior=Δv_gsr), sampler, MCMCThreads(), n_samples, n_threads))

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ e93c66f7-394a-46bf-96c9-475494979548
memb_stars

# ╔═╡ a162219f-df5a-41d8-bf54-927a355f6431
write_fits(joinpath(data_dir, "rv_members_all$FIGSUFFIX.fits"), memb_stars, overwrite=true)

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(samples, LinRange(-280, -180, 100), thin=400)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

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

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
@savefig "rv_sigma_corner" pairplot(samples_gsr)

# ╔═╡ 5f290753-2e2e-4d8b-b1e0-f4beed7fba25
summary_vz = RVUtils.summarize(samples_gsr)

# ╔═╡ bcc6a64f-aba6-4a82-b8df-1597c2c978e8
df_gsr = DataFrame(samples_gsr)

# ╔═╡ a5f1a339-6093-49ea-bd10-b513688d668c
median(df_gsr.μ) + Δv_gsr

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(df_gsr, LinRange(-120, -40, 100), thin=400)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 65fa819a-f73f-4954-b089-a569e7c9b113
CSV.write("processed/mcmc_samples_vz$FIGSUFFIX.csv", df_gsr)

# ╔═╡ 97dea677-17bb-434f-8174-c7b1fc09b329
md"""
# Sigma with Rell
"""

# ╔═╡ f08b0fc3-12c2-4d85-83f3-3fcd9af6641b
model_Rell = RVUtils.model_vel_sigma_R(memb_stars.vz, memb_stars.vz_err, memb_stars.R_ell)

# ╔═╡ b1cd8c53-9b3f-4f07-907d-7f8dd7207902
samples_Rell = sample(model_Rell, sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ a28ee611-3971-4641-b9ec-f338ea3b76ef
summary_Rell = RVUtils.summarize(samples_Rell)

# ╔═╡ b92d373a-657e-4242-9740-213462046331
df_Rell = DataFrame(samples_Rell)

# ╔═╡ 2ab87fa6-92ee-43ab-ad3b-b906dc01654f
@savefig "sigma_Rell_corner" pairplot(samples_Rell)

# ╔═╡ 74a99227-6c33-489b-8fe0-43145048acf1
median(df_Rell.μ) + Δv_gsr

# ╔═╡ f8b877ad-7c7e-4960-a7f3-2f97476d5573
samples_prior_Rell = sample(model_Rell, Prior(), 10000)

# ╔═╡ ba874a4f-ddfe-4a64-9e77-35faa94e6993
bf_sigma_Rell = RVUtils.bayes_evidence(model_Rell, df_Rell, "dlσ_dlR")

# ╔═╡ f11486fa-a88c-4790-a55e-1a6fa5033140
md"""
In the plot below, we just want to make sure that the KDE density estimate looks reasonable at zero
"""

# ╔═╡ 97b2d72c-26aa-429b-90f6-e3ff7d49facf
CSV.write("processed/mcmc_samples_Rell$FIGSUFFIX.csv", df_Rell)

# ╔═╡ 062994bc-fb0c-4f08-b7f7-7cc6714bad1e
md"""
## Gradient
"""

# ╔═╡ 5ff274e6-d57f-441c-9f33-9d622c536c6a
model_gradient = RVUtils.model_vel_gradient(memb_stars.vz, memb_stars.vz_err, memb_stars.xi, memb_stars.eta)

# ╔═╡ cb440be7-1f0e-4a20-b15a-82b1820d1ced
samples_gradient = sample(model_gradient, sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
@savefig "gradient_corner" pairplot(samples_gradient)

# ╔═╡ 184b4a5d-cbab-44b2-9620-bf928ad81d0e
df_gradient = let
	df = DataFrame(samples_gradient)
	df[:, :A] 
	df[:, :B] 
	df[!, :r_grad] = @. 60 * ( df.A ⊕ df.B )
	df[!, :Θ_grad] = @. atand(df.A, df.B) 

	df
end

# ╔═╡ 0a4a2fcd-e0aa-4ab3-b07e-5f57734b2c6b
@savefig "gradient_cyl_corner" pairplot(df_gradient[:, [:μ, :σ, :r_grad, :Θ_grad]])

# ╔═╡ 99362018-8762-40df-b77d-f768286041a6
BF_gradient = RVUtils.bayes_evidence(model_gradient, df_gradient, ["A", "B"])

# ╔═╡ 88f2918e-e126-420a-96a2-5746a8010f73
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ 4a473039-79f0-4d77-aa0c-681e2fba4f4c
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ ed35eb68-74f7-4009-9b68-dfca2ea547af
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

# ╔═╡ 0ca7dc1b-3b41-4089-9c89-20c6e48213ea
@savefig "v_gradient_derived" let
	fig = Figure()

	ax=Axis(fig[1,1];
		  limits=(-12, 12, -12, 12), 
		  aspect=DataAspect(),
		  xlabel = L"$\partial \,v_z / \partial \xi$ / km\,s$^{-1}$\,degree$^{-1}$",
		  ylabel = L"$\partial \,v_z / \partial \eta$ / km\,s$^{-1}$\,degree$^{-1}$",
		xreversed=true
		 )
	
	scatter!(60df_gradient.A, 60df_gradient.B, alpha=0.1, markersize=1, 
		rasterize=4,
	   )

	scatter!(0, 0, color=:black)
	arrows!([0], [0], [5gsr0.pmra], [5gsr0.pmdec])
	arrows!([0], [0], [-5pm_gsr_induced.pmra], [-5pm_gsr_induced.pmdec])

	fig
end

# ╔═╡ 2a422e88-fc0d-4a89-a841-42f3c5c8dace
import KernelDensity

# ╔═╡ 4d0dbbb0-05bc-4f26-be42-b3226f972e28
kde_Rell = KernelDensity.kde(df_Rell.dlσ_dlR)

# ╔═╡ 935cdb6d-56d9-4e95-9983-3ed70fbca11f
let
	fig = Figure()
	ax = Axis(fig[1,1],
		yscale=log10, 
		yticks = Makie.automatic,
		ylabel = "density",
		xlabel = "dlσ/dlR"
	)
		
	
	lines!(kde_Rell)
	vlines!(0, color=:black)

	x = df_Rell.dlσ_dlR[df_Rell.dlσ_dlR .< quantile(df_Rell.dlσ_dlR, 0.001)]
	scatter!(x, 10^-6 .* (1 .+ rand(length(x))), color=COLORS[2], markersize=1)
	fig
end

# ╔═╡ f523cf48-82bf-4b20-9d7c-215bbe10a193
kde = KernelDensity.kde((df_gradient.A, df_gradient.B))

# ╔═╡ 2b88915c-7222-44be-a483-b967ea131b80
log(pdf(kde, 0, 0) ./ lguys.gaussian(0, 0., 0.1)^2)

# ╔═╡ b827e765-646c-4928-9f66-c64e7a20539f
sum(df_gradient.A .> 0)

# ╔═╡ 6371d804-cc73-4ce1-9b36-79fa61780d75
median(atand.(df_gradient.B ./ df_gradient.A))

# ╔═╡ 70a22ef4-2eb1-4094-94e3-5a13fb51b9e6
median(sqrt.(df_gradient.B .^ 2 .+ df_gradient.A .^ 2)) * 60

# ╔═╡ 183f572a-bc0f-435b-a656-2ee2a3057559
quantile(sqrt.(df_gradient.B .^ 2 .+ df_gradient.A .^ 2), [0.16, 0.5, 0.84]) * 60

# ╔═╡ d3fb7136-7600-4782-ba97-f2f785fb3c0a
quantile(atand.(df_gradient.B ./ df_gradient.A), [0.16, 0.5, 0.84]) 

# ╔═╡ 3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 60, normalization=:pdf)
	
	RVUtils.plot_samples!(DataFrame(samples_gradient), LinRange(-110, -30, 100), thin=400)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 128c8ab8-6033-40fe-8b68-c633816df9a2
CSV.write("processed/mcmc_samples_gradient$FIGSUFFIX.csv", df_gradient)

# ╔═╡ 82a0e58a-30a4-4e42-b9c1-cb184eb551aa
md"""
# Misc plots
"""

# ╔═╡ 3a69f395-3c2d-4357-89af-5963d5fa79b8
let
	fig, ax = FigAxis(
		xlabel=L"$\log r_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)
	
	scatter!(log10.(memb_stars.R_ell), memb_stars.RV)


	fig
end

# ╔═╡ 9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
not = !

# ╔═╡ f2313731-5b83-42d6-b624-c618bfb0bb5c
rv_mean_gsr = median(df_gsr.μ)

# ╔═╡ 766146e5-d91a-47fe-969c-a4dc48c457b6
hist(log10.(memb_stars.vz_err))

# ╔═╡ eb16cb5b-562f-4ef5-8a3b-247664d47f52
r_max_tangent = 30

# ╔═╡ 36c37e0b-924d-4a7f-b3ab-81e79e27a059
tangent_bins = LinRange(-r_max_tangent, r_max_tangent, 20)

# ╔═╡ 8f555e41-58d4-4e16-8f4a-1c0f58589b08
outside_bins = abs.(memb_stars.xi) .> r_max_tangent .|| abs.(memb_stars.eta) .> r_max_tangent

# ╔═╡ 8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
rv_mean = icrs0.radial_velocity

# ╔═╡ 53da82d5-2e69-4f88-8043-6694d57cdd91
import StatsBase: weights

# ╔═╡ b5533db0-a734-4d37-9d75-24471634f855
memb_stars

# ╔═╡ 4eea0a17-257e-4d0e-88df-9ff4858771b1
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{arcmin}",
		ylabel = L"\eta / \textrm{arcmin}",
		aspect=DataAspect(),
		xreversed=true,
	)

	w = 1 ./ memb_stars.vz_err .^2

	bins = (33, 25)
	k1 = Arya.histogram2d(memb_stars.xi,memb_stars.eta, bins, weights= w .* memb_stars.vz)
	k2 = Arya.histogram2d(memb_stars.xi,memb_stars.eta, bins, weights=w)

	k1.values ./= k2.values

	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		)
	Colorbar(fig[1, 2], p, label="GSR radial velocity / km/s",
)
	@savefig "vlos_xi_eta_hist"

	fig
end

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
	"vz" => summary_vz |> OrderedDict, 
	"bf_gradient" => BF_gradient,
	"bf_rell" => bf_sigma_Rell,
	"Nmemb" => length(memb_stars.RV), 
	"median_err" => median(memb_stars.RV_err)
)

# ╔═╡ c6edc3b5-accc-44c6-afa5-d4b7ed17fd65
open("processed/mcmc_properties$FIGSUFFIX.toml", "w") do f
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
# ╟─32674b48-0118-414f-9d48-2b1c44c02885
# ╠═29425bbc-05b6-434f-b4b3-17ca39ebf830
# ╠═d8b97f8b-f9c3-4f3e-942f-e25a4f27fff6
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═e93c66f7-394a-46bf-96c9-475494979548
# ╠═a162219f-df5a-41d8-bf54-927a355f6431
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═d0ae48e2-8389-4641-b311-cfb4944b0851
# ╠═49e7d483-1dd1-4406-9bce-58d6e4412e7b
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═a5f1a339-6093-49ea-bd10-b513688d668c
# ╠═5f290753-2e2e-4d8b-b1e0-f4beed7fba25
# ╠═bcc6a64f-aba6-4a82-b8df-1597c2c978e8
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═65fa819a-f73f-4954-b089-a569e7c9b113
# ╟─97dea677-17bb-434f-8174-c7b1fc09b329
# ╠═f08b0fc3-12c2-4d85-83f3-3fcd9af6641b
# ╠═b1cd8c53-9b3f-4f07-907d-7f8dd7207902
# ╠═a28ee611-3971-4641-b9ec-f338ea3b76ef
# ╠═b92d373a-657e-4242-9740-213462046331
# ╠═2ab87fa6-92ee-43ab-ad3b-b906dc01654f
# ╠═74a99227-6c33-489b-8fe0-43145048acf1
# ╠═f8b877ad-7c7e-4960-a7f3-2f97476d5573
# ╠═ba874a4f-ddfe-4a64-9e77-35faa94e6993
# ╠═4d0dbbb0-05bc-4f26-be42-b3226f972e28
# ╟─f11486fa-a88c-4790-a55e-1a6fa5033140
# ╠═935cdb6d-56d9-4e95-9983-3ed70fbca11f
# ╠═97b2d72c-26aa-429b-90f6-e3ff7d49facf
# ╠═062994bc-fb0c-4f08-b7f7-7cc6714bad1e
# ╠═5ff274e6-d57f-441c-9f33-9d622c536c6a
# ╠═cb440be7-1f0e-4a20-b15a-82b1820d1ced
# ╠═e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
# ╠═0a4a2fcd-e0aa-4ab3-b07e-5f57734b2c6b
# ╠═99362018-8762-40df-b77d-f768286041a6
# ╠═184b4a5d-cbab-44b2-9620-bf928ad81d0e
# ╠═0ca7dc1b-3b41-4089-9c89-20c6e48213ea
# ╠═4a473039-79f0-4d77-aa0c-681e2fba4f4c
# ╠═88f2918e-e126-420a-96a2-5746a8010f73
# ╠═ed35eb68-74f7-4009-9b68-dfca2ea547af
# ╠═2a422e88-fc0d-4a89-a841-42f3c5c8dace
# ╠═f523cf48-82bf-4b20-9d7c-215bbe10a193
# ╠═2b88915c-7222-44be-a483-b967ea131b80
# ╠═b827e765-646c-4928-9f66-c64e7a20539f
# ╠═6371d804-cc73-4ce1-9b36-79fa61780d75
# ╠═70a22ef4-2eb1-4094-94e3-5a13fb51b9e6
# ╠═183f572a-bc0f-435b-a656-2ee2a3057559
# ╠═d3fb7136-7600-4782-ba97-f2f785fb3c0a
# ╠═3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
# ╠═128c8ab8-6033-40fe-8b68-c633816df9a2
# ╟─82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═f2313731-5b83-42d6-b624-c618bfb0bb5c
# ╠═766146e5-d91a-47fe-969c-a4dc48c457b6
# ╠═eb16cb5b-562f-4ef5-8a3b-247664d47f52
# ╠═36c37e0b-924d-4a7f-b3ab-81e79e27a059
# ╠═8f555e41-58d4-4e16-8f4a-1c0f58589b08
# ╠═8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
# ╠═53da82d5-2e69-4f88-8043-6694d57cdd91
# ╠═b5533db0-a734-4d37-9d75-24471634f855
# ╠═4eea0a17-257e-4d0e-88df-9ff4858771b1
# ╟─d5615552-caf8-4c0c-a17c-502c0f8198dc
# ╠═09d570a4-56f1-4ff9-990d-c02534f7351e
# ╠═010b6aa7-e3d0-4441-aac5-6ab87c053e33
# ╠═0dbb27cf-8a8c-4521-bda4-5768d8a02176
# ╠═c6edc3b5-accc-44c6-afa5-d4b7ed17fd65
