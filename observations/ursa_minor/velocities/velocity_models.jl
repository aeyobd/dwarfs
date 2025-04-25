### A Pluto.jl notebook ###
# v0.20.6

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

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

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

# ╔═╡ 188f3bd3-5a99-400f-bb65-80c8657c81c6
j24 = read_fits("../processed/best_sample.fits")

# ╔═╡ 1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
f_sat = mean(j24.PSAT)

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
σv = obs_properties["sigma_v"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ 9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
θ_orbit = obs_properties["theta_pm_gsr"]

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
rv_meas = let
	rv_meas = read_fits(joinpath(data_dir, "rv_combined.fits"))

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
## Membership
"""

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = rv_meas.PSAT_RV .> 0.2# ignore spatial 

	memb_filt .&= rv_meas.F_scatter

	memb_filt
end

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 31a6c2e4-538c-4adc-bbda-5043680b17f7
extrema(memb_stars.RV)

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ 7e0304d2-2483-4cec-9dfc-ffb1896acd99
hist(rv_meas.PSAT_RV)

# ╔═╡ 32674b48-0118-414f-9d48-2b1c44c02885
md"""
# Numbers
"""

# ╔═╡ f9d4eade-c648-4f20-8403-07be993fb8c1
sum(.!ismissing.(memb_stars.RV_graces))

# ╔═╡ 29425bbc-05b6-434f-b4b3-17ca39ebf830
median(memb_stars.RV_err)

# ╔═╡ 34ad3ce7-9c48-4c10-9d02-32e6ae99d4fa
length(rv_meas.F_scatter)

# ╔═╡ 4c89d38b-552a-400d-9a8d-2bce4fb8da98
sum(rv_meas.F_scatter)

# ╔═╡ d8b97f8b-f9c3-4f3e-942f-e25a4f27fff6
length(memb_stars.RV)

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f
md"""
## MCMC Priors
"""

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err), Prior(), 10000)
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

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; μ_0_prior=0, N=3_000, p=0.16)
	samples = DataFrame(sample(RVUtils.model_vel_1c(rv, rv_err, μ_0_prior=μ_0_prior), NUTS(0.65), MCMCThreads(), N, 16))

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
samples = DataFrame(sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err, μ_0_prior=Δv_gsr), NUTS(0.65), MCMCThreads(), 10000, 16))

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ e93c66f7-394a-46bf-96c9-475494979548
memb_stars

# ╔═╡ a162219f-df5a-41d8-bf54-927a355f6431
write_fits(joinpath(data_dir, "rv_members_all.fits"), memb_stars, overwrite=true)

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(samples, LinRange(-280, -180, 100), thin=15)
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
samples_gsr = sample(RVUtils.model_vel_1c(memb_stars.vz, memb_stars.vz_err), NUTS(0.65), MCMCThreads(), 10000, 16)

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
@savefig "rv_sigma_corner" pairplot(samples_gsr)

# ╔═╡ 5f290753-2e2e-4d8b-b1e0-f4beed7fba25
summary_vz = RVUtils.summarize(samples_gsr)

# ╔═╡ bcc6a64f-aba6-4a82-b8df-1597c2c978e8
df_gsr = DataFrame(samples_gsr)

# ╔═╡ a5f1a339-6093-49ea-bd10-b513688d668c
median(df_gsr.μ) + Δv_gsr

# ╔═╡ d321f8ac-1044-45ec-8e1c-a2d8395b6917
icrs = lguys.ICRS(ra=obs_properties["ra"], dec=obs_properties["dec"], pmdec=obs_properties["pmdec"], pmra=obs_properties["pmra"], distance=obs_properties["distance"], radial_velocity=obs_properties["radial_velocity"])

# ╔═╡ 1bc7adb7-fe85-4878-9527-c5d15dc761b1
Δv2 = lguys.transform(lguys.GSR, lguys.ICRS(ra=obs_properties["ra"], dec=obs_properties["dec"], radial_velocity=0), ).radial_velocity

# ╔═╡ 931ed52e-5e7a-4692-b568-ae26ea44b638
lguys.transform(lguys.GSR, icrs)

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(df_gsr, LinRange(-120, -40, 100), thin=160)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 97dea677-17bb-434f-8174-c7b1fc09b329
md"""
# Sigma with Rell
"""

# ╔═╡ f08b0fc3-12c2-4d85-83f3-3fcd9af6641b
model_Rell = RVUtils.model_vel_sigma_R(memb_stars.vz, memb_stars.vz_err, memb_stars.R_ell)

# ╔═╡ b1cd8c53-9b3f-4f07-907d-7f8dd7207902
samples_Rell = sample(model_Rell, NUTS(0.65), MCMCThreads(), 1_000, 16)

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

# ╔═╡ 062994bc-fb0c-4f08-b7f7-7cc6714bad1e
md"""
## Gradient
"""

# ╔═╡ 5ff274e6-d57f-441c-9f33-9d622c536c6a
model_gradient = RVUtils.model_vel_gradient(memb_stars.vz, memb_stars.vz_err, memb_stars.xi, memb_stars.eta)

# ╔═╡ cb440be7-1f0e-4a20-b15a-82b1820d1ced
samples_gradient = sample(model_gradient, NUTS(0.65), MCMCThreads(), 1000, 16)

# ╔═╡ e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
@savefig "gradient_corner" pairplot(samples_gradient)

# ╔═╡ 0a4a2fcd-e0aa-4ab3-b07e-5f57734b2c6b
@savefig "gradient_cyl_corner" pairplot(df_gradient[:, [:μ, :σ, :r_grad, :Θ_grad]])

# ╔═╡ 184b4a5d-cbab-44b2-9620-bf928ad81d0e
df_gradient = let
	df = DataFrame(samples_gradient)
	df[:, :A] 
	df[:, :B] 
	df[!, :r_grad] = @. 60 * ( df.A ⊕ df.B )
	df[!, :Θ_grad] = @. atand(df.A, df.B) 

	df
end

# ╔═╡ 99362018-8762-40df-b77d-f768286041a6
BF_gradient = RVUtils.bayes_evidence(model_gradient, df_gradient, ["A", "B"])

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

	   )

	scatter!(0, 0, color=:black)
	arrows!([0], [0], [5gsr0.pmra], [5gsr0.pmdec])
	arrows!([0], [0], [-5pm_gsr_induced.pmra], [-5pm_gsr_induced.pmdec])

	fig
end

# ╔═╡ 88f2918e-e126-420a-96a2-5746a8010f73
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ 4a473039-79f0-4d77-aa0c-681e2fba4f4c
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ ed35eb68-74f7-4009-9b68-dfca2ea547af
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

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

# ╔═╡ d688d2e5-faca-4b14-801a-d58b08fd6654
let
	dr = 10
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{arcmin}",
		ylabel = L"\eta / \textrm{arcmin}",
		aspect=DataAspect()
	)

	p = scatter!(memb_stars.xi, memb_stars.eta, color=memb_stars.radial_velocity_gsr,
		colorrange=(rv_mean_gsr - dr, rv_mean_gsr + dr),
		colormap=:bluesreds,
		markersize = 1 ./ memb_stars.RV_err
	)

	Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

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

# ╔═╡ 6b59d6e9-833b-4582-b5fb-f0a1f69a16c1
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		aspect=DataAspect(),
		xreversed = true,
	)


	k1 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins)


	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=Reverse(:greys)
		)

	scatter!(memb_stars.xi[outside_bins], memb_stars.eta[outside_bins],
		color=:grey
	)
	Colorbar(fig[1, 2], p, label="density",
)
	fig
end

# ╔═╡ 89bf44ef-1ff5-443e-b8be-a3d1571b72e3
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		aspect=DataAspect(),
		xreversed = true
	)

	w = 1 ./ memb_stars.RV_err

	k1 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins, weights= w .* memb_stars.RV)
	k2 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins, weights=w)

	k1.values ./= k2.values

	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		colorrange=(rv_mean - 20, rv_mean + 20)
		)

	scatter!(memb_stars.xi[outside_bins], memb_stars.eta[outside_bins],
		color = memb_stars.RV[outside_bins],
		colormap=:bluesreds,
		colorrange=(rv_mean - 20, rv_mean + 20),
			 markersize=1
	)
	Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s",
)
	fig
end

# ╔═╡ 53da82d5-2e69-4f88-8043-6694d57cdd91
import StatsBase: weights

# ╔═╡ 106482c9-a9f9-4b6e-95a9-614ab7991e23
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		xreversed=true,
		aspect=DataAspect()
	)


	w = 1 ./ memb_stars.vz_err .^ 2
	rv_mean = mean(memb_stars.vz, weights(w))

	k1 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins, weights= w .* memb_stars.vz)
	k2 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins, weights=w)

	k1.values ./= k2.values

	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		colorrange=(rv_mean - 20, rv_mean + 20)
		)

	scatter!(memb_stars.xi[outside_bins], memb_stars.eta[outside_bins],
		color = memb_stars.vz[outside_bins],
		colormap=:bluesreds,
		colorrange=(rv_mean - 20, rv_mean + 20),
			markersize=2
	)
	Colorbar(fig[1, 2], p, label=L"$v_{z}$ / km/s",
)
	fig
end

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

# ╔═╡ 66329d60-ce6a-44ea-9642-2d8b85dfad32


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
	"Nqual" => sum(rv_meas.F_scatter)
)

# ╔═╡ c6edc3b5-accc-44c6-afa5-d4b7ed17fd65
open("processed/mcmc_properties.toml", "w") do f
	TOML.print(f, df_summaries)
end

# ╔═╡ Cell order:
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
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
# ╠═188f3bd3-5a99-400f-bb65-80c8657c81c6
# ╠═1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═31a6c2e4-538c-4adc-bbda-5043680b17f7
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═7e0304d2-2483-4cec-9dfc-ffb1896acd99
# ╠═32674b48-0118-414f-9d48-2b1c44c02885
# ╠═f9d4eade-c648-4f20-8403-07be993fb8c1
# ╠═29425bbc-05b6-434f-b4b3-17ca39ebf830
# ╠═34ad3ce7-9c48-4c10-9d02-32e6ae99d4fa
# ╠═4c89d38b-552a-400d-9a8d-2bce4fb8da98
# ╠═d8b97f8b-f9c3-4f3e-942f-e25a4f27fff6
# ╟─6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═9380d9d1-58b2-432d-8528-d247cf5724e9
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
# ╠═d321f8ac-1044-45ec-8e1c-a2d8395b6917
# ╠═1bc7adb7-fe85-4878-9527-c5d15dc761b1
# ╠═931ed52e-5e7a-4692-b568-ae26ea44b638
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═97dea677-17bb-434f-8174-c7b1fc09b329
# ╠═f08b0fc3-12c2-4d85-83f3-3fcd9af6641b
# ╠═b1cd8c53-9b3f-4f07-907d-7f8dd7207902
# ╠═a28ee611-3971-4641-b9ec-f338ea3b76ef
# ╠═b92d373a-657e-4242-9740-213462046331
# ╠═2ab87fa6-92ee-43ab-ad3b-b906dc01654f
# ╠═74a99227-6c33-489b-8fe0-43145048acf1
# ╠═f8b877ad-7c7e-4960-a7f3-2f97476d5573
# ╠═ba874a4f-ddfe-4a64-9e77-35faa94e6993
# ╠═4d0dbbb0-05bc-4f26-be42-b3226f972e28
# ╠═f11486fa-a88c-4790-a55e-1a6fa5033140
# ╠═935cdb6d-56d9-4e95-9983-3ed70fbca11f
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
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═e05ec20a-3165-4360-866e-3e8cae8665e5
# ╟─82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═f2313731-5b83-42d6-b624-c618bfb0bb5c
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═766146e5-d91a-47fe-969c-a4dc48c457b6
# ╠═eb16cb5b-562f-4ef5-8a3b-247664d47f52
# ╠═36c37e0b-924d-4a7f-b3ab-81e79e27a059
# ╠═8f555e41-58d4-4e16-8f4a-1c0f58589b08
# ╠═8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
# ╠═6b59d6e9-833b-4582-b5fb-f0a1f69a16c1
# ╠═89bf44ef-1ff5-443e-b8be-a3d1571b72e3
# ╠═53da82d5-2e69-4f88-8043-6694d57cdd91
# ╠═106482c9-a9f9-4b6e-95a9-614ab7991e23
# ╠═b5533db0-a734-4d37-9d75-24471634f855
# ╠═4eea0a17-257e-4d0e-88df-9ff4858771b1
# ╠═66329d60-ce6a-44ea-9642-2d8b85dfad32
# ╠═d5615552-caf8-4c0c-a17c-502c0f8198dc
# ╠═09d570a4-56f1-4ff9-990d-c02534f7351e
# ╠═010b6aa7-e3d0-4441-aac5-6ab87c053e33
# ╠═0dbb27cf-8a8c-4521-bda4-5768d8a02176
# ╠═c6edc3b5-accc-44c6-afa5-d4b7ed17fd65
