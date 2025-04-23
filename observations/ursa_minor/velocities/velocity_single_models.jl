### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 35c87efe-6899-444d-952f-94e5e5342846
using Measurements

# ╔═╡ c05a4e8e-c6d2-43e6-9b87-679f4ecee84c
study = "s18"

# ╔═╡ b3ad64bc-7dbb-4c27-a91d-a9f2342f98da
begin 
	FIGDIR = "./figures/"
	FIGSUFFIX = ".$study"
	using LilGuys
end

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ 9039df68-748a-41be-9a32-e2bf45bee85a
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 5fd838b6-e7b1-4449-8f7a-b26a77042f0f
⊕ = RVUtils.:⊕

# ╔═╡ 48bb25f2-19a7-4145-87a3-b7083c10d186
import KernelDensity

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 8863e6f7-6aaa-4370-919f-fe6afd3847bf
θ_orbit = obs_properties["theta_pm_gsr"]

# ╔═╡ d476513e-6211-44a4-ade5-0ab184a9829d
j24 = read_fits("../processed/best_sample.fits")

# ╔═╡ 8188ae3b-2064-4073-b9d0-c04f50b78707


# ╔═╡ 1bc7adb7-fe85-4878-9527-c5d15dc761b1
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ e391b59d-b8d9-4726-8e4d-8476f6d62800
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ cd71ec9b-4b8c-47ca-a0c8-3c6cbc257927
σv = obs_properties["sigma_v"]

# ╔═╡ 03636f8d-0fc9-4f20-8a84-b92f9bd27ffb
f_sat = mean(j24.PSAT)

# ╔═╡ 0930d5db-2b03-4cf1-ac54-80a592d959ed
rv_meas_all = read_fits("processed/rv_combined.fits")

# ╔═╡ 8dd53e5a-279b-4356-b7f8-dcfdbf2cd532
sum(.!ismissing.(rv_meas_all.RV_s18))

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
rv_meas = let 
	rv_meas = filter(x->!ismissing(x["RV_$study"]), rv_meas_all)

	# rename columns
	rv_meas[:, :RV] = rv_meas[!, "RV_$study"]
	rv_meas[:, :RV_err] .= rv_meas[!, "RV_err_$study"]
	rv_meas[:, :RV_sigma] .= replace(rv_meas[!, "RV_sigma_$study"], missing=>NaN)
	rv_meas[:, :RV_count] .= rv_meas[!, "RV_count_$study"]
	rv_meas[:, :F_scatter] .= rv_meas[:, "F_scatter_$study"]
	rv_meas[:, :F_match] .= rv_meas[:, "F_match_$study"]

	
	RVUtils.add_gsr!(rv_meas, 
					 distance=obs_properties["distance"],
					pmra=obs_properties["pmra"], pmdec=obs_properties["pmdec"])
	
	rv_meas[:, :vz] .= rv_meas.radial_velocity_gsr .+ rv_meas.delta_rv
	rv_meas[:, :vz_err] .= rv_meas.RV_err .⊕ rv_meas.delta_rv_err
	
	RVUtils.add_PSAT_RV!(rv_meas; sigma_v=σv, radial_velocity_gsr=rv0, f_sat=f_sat)
	
	rv_meas
end

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership criteria
"""

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = rv_meas.PSAT_RV .> 0.2# ignore spatial 
	memb_filt .&= rv_meas.F_scatter

	memb_filt
end

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ 29405c62-01f3-4f6e-8d92-3b6c36551080
hist(rv_meas.RV, bins=60)

# ╔═╡ 372adffe-ece1-4c26-a7ad-e897ba6b829b
md"""
# Numbers
"""

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 130c9e30-d900-436e-9749-c915edc43b3d
median(memb_stars.RV_err)

# ╔═╡ 398dbf3c-21a2-4785-864e-0f0a04186361
hist(memb_stars.RV_err)

# ╔═╡ 4c1af687-cb83-4441-99b3-6de928edc656
sum(rv_meas.F_scatter)

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; μ_0_prior=0, N=3_000, p=0.16, burn=0.2)
	samples = DataFrame(sample(RVUtils.model_vel_1c(rv, rv_err, μ_0_prior=μ_0_prior), NUTS(0.65), MCMCThreads(), N, 16))
	burn = round(Int, burn * N)
	samples=samples[burn:end, :]

	μ = median(samples.μ)
	μ_p = quantile(samples.μ, [p, 1-p]) 
	μ_err = (μ - μ_p[1], μ_p[2] - μ)

	σ = median(samples.σ)
	σ_p = quantile(samples.σ, [p, 1-p])
	σ_err = (σ - σ_p[1], σ_p[2] - σ)

	return μ, σ, μ_err, σ_err
end

# ╔═╡ 9380d9d1-58b2-432d-8528-d247cf5724e9
function plot_sample_normal_fit(sample, props; kwargs...)
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density";
		kwargs...
	)
	h = Arya.histogram(Float64.(sample.RV), normalization=:pdf)
	
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=:black)

	μ, σ, _, _ = props
	x_model = LinRange(80, 140, 1000)
	y_model = lguys.gaussian.(x_model, μ, σ)
	lines!(x_model, y_model)

	fig
end

# ╔═╡ 318d29b9-4c84-4d38-af6d-048518952970
samples = DataFrame(sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err), NUTS(0.65), MCMCThreads(), 10000, 16))

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ e93c66f7-394a-46bf-96c9-475494979548
memb_stars

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(samples, LinRange(-280, -200, 100), thin=15)
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
mean(memb_stars.RV)

# ╔═╡ 9ffcc0b0-11a1-4962-b64e-59a731f22bd8
samples_gsr = DataFrame(sample(RVUtils.model_vel_1c(memb_stars.vz, memb_stars.vz_err), NUTS(0.65), MCMCThreads(), 10000, 16))

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
pairplot(samples_gsr[:, [:μ, :σ]])

# ╔═╡ 212cee31-a04d-47c1-b56e-4a7b06322b99
median(samples_gsr.μ) + Δv_gsr

# ╔═╡ 0eafe553-4434-4b2b-9646-912ac99beccf


# ╔═╡ d321f8ac-1044-45ec-8e1c-a2d8395b6917
# ╠═╡ disabled = true
#=╠═╡
icrs = lguys.ICRS(ra=obs_properties["ra"], dec=obs_properties["dec"], pmdec=obs_properties["pmdec"], pmra=obs_properties["pmra"], distance=obs_properties["distance"], radial_velocity=obs_properties["radial_velocity"])
  ╠═╡ =#

# ╔═╡ 931ed52e-5e7a-4692-b568-ae26ea44b638
#=╠═╡
lguys.transform(lguys.GSR, icrs)
  ╠═╡ =#

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.vz), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(samples_gsr, LinRange(-110, -50, 100), thin=15)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 18765306-a840-4173-99e1-ee7c291d2652
md"""
# Sigma with Rell
"""

# ╔═╡ 833c2cb2-fc23-4bc8-8f01-0b5de4fb9c5e
model_Rell = RVUtils.model_vel_sigma_R(memb_stars.vz, memb_stars.vz_err, memb_stars.R_ell)

# ╔═╡ 399867cf-5b89-4ed6-8686-d71c6a9d8996
samples_Rell = sample(model_Rell, NUTS(0.65), MCMCThreads(), 1_000, 16)

# ╔═╡ 8f5d079b-19f4-4791-9f23-1c09350d5de6
summary_Rell = RVUtils.summarize(samples_Rell)

# ╔═╡ 132116fa-03f1-4768-ad46-062839b5d2b1
df_Rell = DataFrame(samples_Rell)

# ╔═╡ f93b0032-3000-4d04-8f08-d25ad9448e07
@savefig "sigma_Rell_corner" pairplot(samples_Rell)

# ╔═╡ 8df22f07-4319-4b60-9aae-a0c31b07564d
median(df_Rell.μ) + Δv_gsr

# ╔═╡ e2aa3d66-d72d-4510-a81c-8b46d6ed90e7
samples_prior_Rell = sample(model_Rell, Prior(), 10000)

# ╔═╡ 68250c58-8b67-42ce-b447-70f87aabecbc
RVUtils.bayes_evidence(model_Rell, df_Rell, "dlσ_dlR")

# ╔═╡ 6f734c62-3545-48fb-ae95-9877a7de151e
kde_Rell = KernelDensity.kde(df_Rell.dlσ_dlR)

# ╔═╡ 1bb3f1ff-5800-4a95-a09d-5de2d1136163
md"""
In the plot below, we just want to make sure that the KDE density estimate looks reasonable at zero
"""

# ╔═╡ d3a8e602-85a5-441f-aa23-db0bc0b5dcc6
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

# ╔═╡ 79b2bc05-793a-48ce-91c9-774c23ff5e25
md"""
## Gradient
"""

# ╔═╡ 87937ec6-459d-4eb6-824d-68b41f9ed486
model_gradient = RVUtils.model_vel_gradient(memb_stars.vz, memb_stars.vz_err, memb_stars.xi, memb_stars.eta)

# ╔═╡ 3c68dc7c-939e-4f13-8d1c-61afa884ec4b
samples_gradient = sample(model_gradient, NUTS(0.65), MCMCThreads(), 1000, 16)

# ╔═╡ af5c3c0b-4cf6-41f7-be18-d61ff52d625d
@savefig "gradient_corner" pairplot(samples_gradient)

# ╔═╡ a1a8bd58-034a-42bb-ac1c-1dafb4507d5c
df_gradient = let
	df = DataFrame(samples_gradient)
	df[:, :A] 
	df[:, :B] 
	df[!, :r_grad] = @. 60 * ( df.A ⊕ df.B )
	df[!, :Θ_grad] = @. atand(df.A, df.B) 

	df
end

# ╔═╡ 9d2ccb07-548b-43fe-ae66-b3c44fb8a98b
@savefig "gradient_cyl_corner" pairplot(df_gradient[:, [:μ, :σ, :r_grad, :Θ_grad]])

# ╔═╡ 9724df49-645c-43c3-8385-346f78ac82af
RVUtils.bayes_evidence(model_gradient, df_gradient, ["A", "B"])

# ╔═╡ 993526b6-bde8-4011-bbda-84d044841980
RVUtils.bayes_evidence(model_gradient, df_gradient, ["A", "B"])

# ╔═╡ d5994630-7cc6-4430-b67d-ad090ccccad9
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ 84db6a39-ea47-4b1d-8adc-333a45b626dd
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ 1b4f5ff9-328f-4cc4-a07c-ba3a3862290c
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

# ╔═╡ 7791974e-d917-4925-9630-2faa62b993c0
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

	   )

	scatter!(0, 0, color=:black)
	arrows!([0], [0], [5gsr0.pmra], [5gsr0.pmdec])
	arrows!([0], [0], [-5pm_gsr_induced.pmra], [-5pm_gsr_induced.pmdec])

	fig
end

# ╔═╡ 80b38274-cb7b-4aec-9005-3204cee1b268
# ╠═╡ disabled = true
#=╠═╡
import KernelDensity
  ╠═╡ =#

# ╔═╡ fb13e9e3-4746-4d33-8036-306efc4295a9
kde = KernelDensity.kde((df_gradient.A, df_gradient.B))

# ╔═╡ 1fb67165-f697-4ff9-bee8-d528551dd45d


# ╔═╡ f7fcb16a-4298-435f-aa9f-80f000e99876
log(pdf(kde, 0, 0) ./ lguys.gaussian(0, 0., 0.1)^2)

# ╔═╡ 773d2361-33ee-4b8e-8664-8cb63efe2f68
sum(df_gradient.A .> 0)

# ╔═╡ 969a3aab-2394-4ab1-bfa1-b11811214e05
median(atand.(df_gradient.B ./ df_gradient.A))

# ╔═╡ a9a9dfad-2b6b-4da6-8d66-1eac90ea456c
median(sqrt.(df_gradient.B .^ 2 .+ df_gradient.A .^ 2)) * 60

# ╔═╡ 430f54ed-619d-4d72-a9a9-69e14af9252f
quantile(sqrt.(df_gradient.B .^ 2 .+ df_gradient.A .^ 2), [0.16, 0.5, 0.84]) * 60

# ╔═╡ 652c230d-957a-4795-ab52-9cd892962d04
quantile(atand.(df_gradient.B ./ df_gradient.A), [0.16, 0.5, 0.84])

# ╔═╡ dc3e206b-01a5-4fc5-ac01-4a04f0261947
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 60, normalization=:pdf)
	
	RVUtils.plot_samples!(DataFrame(samples_gradient), LinRange(-110, -30, 100), thin=15)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties along orbit
"""

# ╔═╡ c2735c49-2892-46ac-bcf8-7cdcef409f44
function calc_binned_mu_sigma(x, y, yerr, bins; kwargs...)

	if !issorted(bins)
		error("bins must be sorted")
	end

	N = length(bins) - 1
	μs = Vector{Float64}(undef, N)
	σs = Vector{Float64}(undef, N)
	μ_errs = Vector{Tuple{Float64, Float64}}(undef, N)
	σ_errs = Vector{Tuple{Float64, Float64}}(undef, N)
	xs = Vector{Float64}(undef, N)

	
	for i in 1:N
		filt = x .>= bins[i]
		filt .&= x .< bins[i+1]
		@info "calculating bin $i"

		μs[i], σs[i], μ_errs[i], σ_errs[i] = fit_rv_sigma(y[filt], yerr[filt]; kwargs...)
		xs[i] = median(x[filt])
	end

	return DataFrame(
		x=xs,
		x_low = bins[1:end-1],
		x_high = bins[2:end],
		μ=μs, 
		σ = σs, 
		μ_err = μ_errs, 
		σ_err = σ_errs,
		
	)
end	

# ╔═╡ 8b21cc49-ca17-4844-8238-e27e9752bee7
bins = bins_equal_number(memb_stars.R_ell, n=10)

# ╔═╡ f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
df_r_ell_z = calc_binned_mu_sigma(memb_stars.R_ell, memb_stars.vz, memb_stars.vz_err, bins)

# ╔═╡ 38da4da1-74f5-4661-89f2-4b25562a1faf
function scatter_range!(df_r_ell)
	errorscatter!(df_r_ell.x, df_r_ell.μ, yerror=df_r_ell.μ_err, color=:black)
	
	errorbars!(df_r_ell.x, df_r_ell.μ .+ df_r_ell.σ, df_r_ell.x .- df_r_ell.x_low, df_r_ell.x_high .- df_r_ell.x,  direction = :x, color=:black)
	errorbars!(df_r_ell.x, df_r_ell.μ .- df_r_ell.σ, df_r_ell.x .- df_r_ell.x_low, df_r_ell.x_high .- df_r_ell.x, direction = :x, color=:black)
end

# ╔═╡ e05ec20a-3165-4360-866e-3e8cae8665e5
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.R_ell, memb_stars.vz, color=COLORS[3], alpha=0.1)

	scatter_range!(df_r_ell_z)

	fig
end

# ╔═╡ 6bd9dd7d-7185-4969-a13c-7ee41ac5c8fa
let
	fig = Figure()
	ax = Axis(fig[1,1],
			  xlabel = "Rell / arcmin",
			  ylabel = "sigma"
			 )


	errorscatter!(df_r_ell_z.x, df_r_ell_z.σ, yerror=df_r_ell_z.σ_err)

	fig
end

# ╔═╡ eb16cb5b-562f-4ef5-8a3b-247664d47f52
r_max_tangent = 0.4

# ╔═╡ 36c37e0b-924d-4a7f-b3ab-81e79e27a059
tangent_bins = LinRange(-r_max_tangent, r_max_tangent, 20)

# ╔═╡ 8f555e41-58d4-4e16-8f4a-1c0f58589b08
outside_bins = abs.(memb_stars.xi) .> r_max_tangent .|| abs.(memb_stars.eta) .> r_max_tangent

# ╔═╡ 53da82d5-2e69-4f88-8043-6694d57cdd91
import StatsBase: weights

# ╔═╡ b5533db0-a734-4d37-9d75-24471634f855
memb_stars

# ╔═╡ 4eea0a17-257e-4d0e-88df-9ff4858771b1
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
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
		colorrange=(rv0 - 3σv, rv0 + 3σv)
		)
	Colorbar(fig[1, 2], p, label=L"$v_z$ / km/s",
)

		@savefig "vlos_xi_eta_hist"

	fig
end

# ╔═╡ Cell order:
# ╠═c05a4e8e-c6d2-43e6-9b87-679f4ecee84c
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═93838644-cad6-4df3-b554-208b7afeb3b8
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╠═9039df68-748a-41be-9a32-e2bf45bee85a
# ╠═b3ad64bc-7dbb-4c27-a91d-a9f2342f98da
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═35c87efe-6899-444d-952f-94e5e5342846
# ╠═5fd838b6-e7b1-4449-8f7a-b26a77042f0f
# ╠═48bb25f2-19a7-4145-87a3-b7083c10d186
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═8863e6f7-6aaa-4370-919f-fe6afd3847bf
# ╠═d476513e-6211-44a4-ade5-0ab184a9829d
# ╠═8188ae3b-2064-4073-b9d0-c04f50b78707
# ╠═e391b59d-b8d9-4726-8e4d-8476f6d62800
# ╠═1bc7adb7-fe85-4878-9527-c5d15dc761b1
# ╠═cd71ec9b-4b8c-47ca-a0c8-3c6cbc257927
# ╠═03636f8d-0fc9-4f20-8a84-b92f9bd27ffb
# ╠═0930d5db-2b03-4cf1-ac54-80a592d959ed
# ╠═8dd53e5a-279b-4356-b7f8-dcfdbf2cd532
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═29405c62-01f3-4f6e-8d92-3b6c36551080
# ╟─372adffe-ece1-4c26-a7ad-e897ba6b829b
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═130c9e30-d900-436e-9749-c915edc43b3d
# ╠═398dbf3c-21a2-4785-864e-0f0a04186361
# ╠═4c1af687-cb83-4441-99b3-6de928edc656
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═9380d9d1-58b2-432d-8528-d247cf5724e9
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═e93c66f7-394a-46bf-96c9-475494979548
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═d0ae48e2-8389-4641-b311-cfb4944b0851
# ╠═49e7d483-1dd1-4406-9bce-58d6e4412e7b
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═212cee31-a04d-47c1-b56e-4a7b06322b99
# ╠═0eafe553-4434-4b2b-9646-912ac99beccf
# ╠═d321f8ac-1044-45ec-8e1c-a2d8395b6917
# ╠═931ed52e-5e7a-4692-b568-ae26ea44b638
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═18765306-a840-4173-99e1-ee7c291d2652
# ╠═833c2cb2-fc23-4bc8-8f01-0b5de4fb9c5e
# ╠═399867cf-5b89-4ed6-8686-d71c6a9d8996
# ╠═8f5d079b-19f4-4791-9f23-1c09350d5de6
# ╠═132116fa-03f1-4768-ad46-062839b5d2b1
# ╠═f93b0032-3000-4d04-8f08-d25ad9448e07
# ╠═8df22f07-4319-4b60-9aae-a0c31b07564d
# ╠═e2aa3d66-d72d-4510-a81c-8b46d6ed90e7
# ╠═68250c58-8b67-42ce-b447-70f87aabecbc
# ╠═6f734c62-3545-48fb-ae95-9877a7de151e
# ╠═1bb3f1ff-5800-4a95-a09d-5de2d1136163
# ╠═d3a8e602-85a5-441f-aa23-db0bc0b5dcc6
# ╠═79b2bc05-793a-48ce-91c9-774c23ff5e25
# ╠═87937ec6-459d-4eb6-824d-68b41f9ed486
# ╠═3c68dc7c-939e-4f13-8d1c-61afa884ec4b
# ╠═af5c3c0b-4cf6-41f7-be18-d61ff52d625d
# ╠═9d2ccb07-548b-43fe-ae66-b3c44fb8a98b
# ╠═9724df49-645c-43c3-8385-346f78ac82af
# ╠═a1a8bd58-034a-42bb-ac1c-1dafb4507d5c
# ╠═7791974e-d917-4925-9630-2faa62b993c0
# ╠═993526b6-bde8-4011-bbda-84d044841980
# ╠═84db6a39-ea47-4b1d-8adc-333a45b626dd
# ╠═d5994630-7cc6-4430-b67d-ad090ccccad9
# ╠═1b4f5ff9-328f-4cc4-a07c-ba3a3862290c
# ╠═80b38274-cb7b-4aec-9005-3204cee1b268
# ╠═fb13e9e3-4746-4d33-8036-306efc4295a9
# ╠═1fb67165-f697-4ff9-bee8-d528551dd45d
# ╠═f7fcb16a-4298-435f-aa9f-80f000e99876
# ╠═773d2361-33ee-4b8e-8664-8cb63efe2f68
# ╠═969a3aab-2394-4ab1-bfa1-b11811214e05
# ╠═a9a9dfad-2b6b-4da6-8d66-1eac90ea456c
# ╠═430f54ed-619d-4d72-a9a9-69e14af9252f
# ╠═652c230d-957a-4795-ab52-9cd892962d04
# ╠═dc3e206b-01a5-4fc5-ac01-4a04f0261947
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═e05ec20a-3165-4360-866e-3e8cae8665e5
# ╠═6bd9dd7d-7185-4969-a13c-7ee41ac5c8fa
# ╠═eb16cb5b-562f-4ef5-8a3b-247664d47f52
# ╠═36c37e0b-924d-4a7f-b3ab-81e79e27a059
# ╠═8f555e41-58d4-4e16-8f4a-1c0f58589b08
# ╠═53da82d5-2e69-4f88-8043-6694d57cdd91
# ╠═b5533db0-a734-4d37-9d75-24471634f855
# ╠═4eea0a17-257e-4d0e-88df-9ff4858771b1
