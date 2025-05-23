### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 50488b8f-6886-4191-8778-af66929f1445
rv_file = "rv_combined_x_wide_2c_psat_0.2.fits"

# ╔═╡ 4827c99d-6326-42fb-a329-d56992f6d20d
j24_sample = "wide_2c"

# ╔═╡ 8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
n_samples = 10000

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

# ╔═╡ 1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
f_sat = TOML.parsefile("../j+24_fsat.toml")[j24_sample]["f_sat"]

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
σv = obs_properties["sigma_v"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ 9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
θ_orbit = obs_properties["theta_pm_gsr"]

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
memb_stars = let
	rv_meas = read_fits(joinpath(data_dir, rv_file))

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


# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═╡ disabled = true
#=╠═╡
memb_stars = rv_meas
  ╠═╡ =#

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

# ╔═╡ 45d748e7-42a0-471e-b290-e3469c0aa283
memb_stars[.!ismissing.(memb_stars.RV_gmos), [:source_id, :ra, :dec, :PSAT]]

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

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; μ_0_prior=0, N=3_000, p=0.16)
	samples = DataFrame(sample(RVUtils.model_vel_1c(rv, rv_err, μ_0_prior=μ_0_prior), sampler, MCMCThreads(), N, n_threads))

	μ = median(samples.μ)
	μ_p = quantile(samples.μ, [p, 1-p]) 
	μ_err = (μ - μ_p[1], μ_p[2] - μ)

	σ = median(samples.σ)
	σ_p = quantile(samples.σ, [p, 1-p])
	σ_err = (σ - σ_p[1], σ_p[2] - σ)

	return μ, σ, μ_err, σ_err
end

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
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	RVUtils.plot_samples!(samples, LinRange(70, 150, 100), thin=15)
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
	h = histogram(Float64.(memb_stars.vz), 60, normalization=:pdf)
	
	RVUtils.plot_samples!(df_gsr, LinRange(0, 130, 100), thin=160)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	@savefig "RV_hist_mcmc"
	fig
end

# ╔═╡ c1cf3a80-1c26-4417-9ddf-1046f3f5f608
CSV.write("processed/mcmc_samples_vz$FIGSUFFIX.csv", samples_gsr)

# ╔═╡ 433d1d51-4b2f-4bc0-807a-813eb215231a
md"""
# Rell sigma
"""

# ╔═╡ b32a6521-dae5-4d2a-89ab-6e8b3bc43be3
model_Rell = RVUtils.model_vel_sigma_R(memb_stars.vz, memb_stars.vz_err, memb_stars.R_ell)

# ╔═╡ 93ee52ad-2c6d-4eac-85cd-ad2a88a70428
samples_Rell = sample(model_Rell, sampler, MCMCThreads(), n_samples, n_threads)

# ╔═╡ b983cbd3-a90f-4ed1-b0a7-fd21e903f81d
summary_Rell = RVUtils.summarize(samples_Rell)

# ╔═╡ 527ccde2-64db-40f5-9b77-8898a92b3b6a
df_Rell = DataFrame(samples_Rell)

# ╔═╡ 820b82a1-63a4-4a83-8d6e-5659ef75ce2d
@savefig "sigma_Rell_corner" pairplot(samples_Rell)

# ╔═╡ 924369b0-aff6-46e4-a0bf-c6766ff93cbf
median(df_Rell.μ) + Δv_gsr

# ╔═╡ 4f00f714-6b16-41c7-a3c8-f2743d6023bf


# ╔═╡ 137bbc9c-1e02-4f78-acd1-0e30cf982615
bf_sigma_Rell = RVUtils.bayes_evidence(model_Rell, df_Rell, "dlσ_dlR")

# ╔═╡ 3b7cb3ea-3684-4df9-b9e6-919fa11e9ffa
kde_Rell = KernelDensity.kde(df_Rell.dlσ_dlR)

# ╔═╡ f8cc1ebb-617a-46a0-ab22-2f84d75da2ea
md"""
In the plot below, we just want to make sure that the KDE density estimate looks reasonable at zero
"""

# ╔═╡ 9158e64f-5a8d-48bf-8043-8e93d41bd626
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

# ╔═╡ 4a83224a-1f2a-4f7c-b670-6a694fc8826a
KernelDensity.default_bandwidth(df_Rell.dlσ_dlR)

# ╔═╡ d051635b-123c-45e2-94da-6511f179a2bb
CSV.write("processed/mcmc_samples_Rell_sigma$FIGSUFFIX.csv", df_Rell)

# ╔═╡ 062994bc-fb0c-4f08-b7f7-7cc6714bad1e
md"""
## Gradient
"""

# ╔═╡ 25d9799f-5e89-4ec3-995c-a7309efb798f
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
	df[!, :r_grad] = @. 60 * (df.A ⊕ df.B )
	df[!, :Θ_grad] = @. atand(df.A, df.B) 

	df
end

# ╔═╡ e79d9d82-84cd-4bc5-9ea8-be6b07cacf6d
@savefig "gradient_cyl_corner" pairplot(df_gradient[:, [:μ, :σ, :r_grad, :Θ_grad]])

# ╔═╡ f2aa15e2-824b-4f0e-8ad0-42611abe06da
median(df_gradient.μ) + Δv_gsr

# ╔═╡ ad607c95-0aea-4d84-8f18-e5b929b9bfca
bf_gradient = RVUtils.bayes_evidence(model_gradient, df_gradient, ["A", "B"])

# ╔═╡ 65f7b7e0-2439-477e-8036-a28a6903955f
θs = mod1.(df_gradient.Θ_grad, 360.) .- 360

# ╔═╡ 96b41e0b-b9f3-41ba-bc56-f3febd8833d3
quantile(θs, [0.16, 0.5, 0.84]) .- median(θs)

# ╔═╡ cbf50d46-219b-489d-b359-6b473b85735d
median(θs)

# ╔═╡ 4d7e7b96-74a0-4066-88cf-739c043c7f47
summary_gradient = RVUtils.summarize(samples_gradient)

# ╔═╡ 8555e608-53c1-40d3-b21e-413af8953c30
md"""
code below validates induced PM gradient (should be approx 2).
"""

# ╔═╡ 7a9d2e0f-4dcb-4b69-9f68-d43d6dde8bf2
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ 3195286d-d85f-43a3-aa25-dae2134f570b
xi_rot, eta_rot = lguys.to_orbit_coords(memb_stars.ra, memb_stars.dec, obs_properties["ra"], obs_properties["dec"], θ_m) .* 60


# ╔═╡ 88f2918e-e126-420a-96a2-5746a8010f73
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ c48bb30d-1186-4940-b061-91f53e8335e1
vec_pm = LilGuys.pm2kms.([icrs0.pmra, icrs0.pmdec], icrs0.distance) /(180/π)

# ╔═╡ 0ca7dc1b-3b41-4089-9c89-20c6e48213ea
@savefig "v_gradient_derived" let
	fig = Figure()

	ax=Axis(fig[1,1];
		  limits=(-10, 10, -10, 10), 
		  aspect=DataAspect(),
		  xlabel = L"$\partial \,v_z / \partial\, \xi$ / km\,s$^{-1}$\,degree$^{-1}$",
		  ylabel = L"$\partial \,v_z / \partial\, \eta$ / km\,s$^{-1}$\,degree$^{-1}$",
		xreversed=true
		 )
	
	scatter!(60df_gradient.A, 60df_gradient.B, alpha=0.1, markersize=1, 

	   )

	scatter!(0, 0, color=:black)
	arrows!([0], [0], [vec_pm[1]], [vec_pm[2]])

	fig
end

# ╔═╡ af0d2050-b42e-4a7f-aabb-5b08d23381e9
lguys.transform(ICRS, lguys.GSR(ra=icrs0.ra + 2/vec_pm[1] *cos(icrs0.dec), dec=icrs0.dec + 2/vec_pm[2], distance=icrs0.distance, radial_velocity=0)).radial_velocity  .- Δv_gsr

# ╔═╡ 4a473039-79f0-4d77-aa0c-681e2fba4f4c
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ ed35eb68-74f7-4009-9b68-dfca2ea547af
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

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

# ╔═╡ 51c7c22e-405b-4a47-beb1-b0fcf3efa391
KernelDensity.default_bandwidth((df_gradient.A, df_gradient.B))

# ╔═╡ 39bb4bb9-8445-4ee1-a3af-639d8fa96f65
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	
	contour!(kde.x, kde.y, asinh.(kde.density ./ 1e-5), levels=100,)

	scatter!(0, 0, color=:black)


	fig
end

# ╔═╡ 3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 60, normalization=:pdf)
	
	RVUtils.plot_samples!(DataFrame(samples_gradient), LinRange(40, 110, 100), thin=15)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 7fce47a1-16be-4dad-afb5-0fb04bd91355
CSV.write("processed/mcmc_samples_gradient$FIGSUFFIX.csv", df_gradient)

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

# ╔═╡ 03928d42-ee61-4480-82d3-a482842f8521
CSV.write("processed/vz_r_ell_binned$FIGSUFFIX.csv", df_r_ell_z)

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

# ╔═╡ 31aa8fc5-1415-4c44-9b92-a7d097181639
md"""
# Binned properties with radius
"""

# ╔═╡ 2f5db664-fc99-41ac-a726-f21dd5d88ad4
df_xi_p = calc_binned_mu_sigma(xi_rot, memb_stars.vz, memb_stars.vz_err, bins_equal_number(xi_rot, n=10))

# ╔═╡ 301ad1d5-83a1-418d-aa06-583f7212541d
CSV.write("processed/vz_xi_p_binned$FIGSUFFIX.csv", df_xi_p)

# ╔═╡ 090acae4-1209-49e6-882c-20ac2c972dd5
df_eta_p = calc_binned_mu_sigma(eta_rot, memb_stars.vz, memb_stars.vz_err, bins_equal_number(eta_rot, n=10))

# ╔═╡ fd0a74a1-6513-4612-8181-745d5b7c3f4c
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.xi, memb_stars.radial_velocity_gsr, color=COLORS[3], alpha=0.1)

	scatter_range!(df_xi_p)

	fig
end

# ╔═╡ 3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
@savefig "vel_gradient_binned" let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degree",
		ylabel = L"$\mu_{v,z}$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x ./ 60, df_xi_p.μ, yerror=df_xi_p.μ_err, color=:black)

	for i in 1:400:size(df_gradient, 1)
		M = 60df_gradient.B[i] .* sind(θ_m) .+ 60df_gradient.A[i] .* cosd(θ_m)
		x=LinRange(-0.4, 0.4, 100)
		y = M*x  .+ df_gradient.μ[i]
		lines!(x, y, alpha=0.03, color=COLORS[1])
	end
			
	fig
end

# ╔═╡ bf8a9940-1c8b-4fdd-a231-1861721ea8e9
sind(θ_m), cosd(θ_m)

# ╔═╡ 14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
let
	fig, ax = FigAxis(
		xlabel = L"$\eta'$ / arcmin",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(df_eta_p.x, df_eta_p.μ, yerror=df_eta_p.μ_err, color=:black)

	fig
end

# ╔═╡ d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / arcmin",
		ylabel = L"$\sigma_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x, df_xi_p.σ, yerror=df_xi_p.σ_err, color=:black)


		
	@savefig "sigma_v_vs_eta_p"
	fig
end

# ╔═╡ dbea26bd-539f-4eda-a4de-ed00efdc93a9
df_gradient

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ 0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
σ_m = median(samples.σ)

# ╔═╡ f82d2ff7-7a7f-4520-811e-126f3f4f5349
let
	fig, ax = FigAxis(
		xlabel = L"$R_\textrm{ell}$ / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errorscatter!(midpoints(bins), df_r_ell_z.σ, yerror=df_r_ell_z.σ_err, xerror=bin_errs, color=:black)
	hlines!(σ_m)

	@savefig "sigma_v_vs_rell"
	fig
end

# ╔═╡ e50a1968-9f91-4751-a7a9-0a4768c9ab7e
let
	fig, ax = FigAxis(
		xlabel = L"$\log\ R$ / arcmin",
		ylabel = L"$\log\ \sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	x_mid = midpoints(log10.(bins))
	log_bin_errs = [(x_mid[i] - log10(bins[i]), log10(bins[i+1]) - x_mid[i]) for i in eachindex(x_mid)]
	@info log_bin_errs
	
	errorscatter!(x_mid, log10.(df_r_ell_z.σ), yerror=maximum.(df_r_ell_z.σ_err) ./ df_r_ell_z.σ ./ log(10), xerror=log_bin_errs, color=:black)
	hlines!(log10.(σ_m), color=:black)

	for i in 1:400:size(df_Rell, 1)
		x=LinRange(-0.5, 2, 100)
		y = df_Rell.σ[i] .* 10 .^ ((x .- 1.0) .* df_Rell.dlσ_dlR[i])
		lines!(x, log10.(y), alpha=0.03, color=COLORS[1])
	end
	@savefig "log_sigma_log_R"
	fig
end

# ╔═╡ 319bd778-7e17-4bd7-856f-d6785b287219
quantile(samples.σ, [0.16, 0.84]) .- σ_m

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

	p = scatter!(memb_stars.xi, memb_stars.eta, color=memb_stars.vz,
		colorrange=(rv_mean_gsr - dr, rv_mean_gsr + dr),
		colormap=:bluesreds,
		markersize = 1 ./ memb_stars.vz_err
	)

	Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ 687b79fb-b0ce-4a93-8c0e-83b77aec081e
let
	fig, ax = FigAxis(
		xlabel = L"\xi' / \textrm{arcmin}",
		ylabel = L"RV",
		#aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false
	)

	p = scatter!(xi_rot, memb_stars.RV,
		alpha=0.1
	)

	#Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ 9c66756b-af33-48fb-9947-201d53cdeba3
let
	fig, ax = FigAxis(
		xlabel = L"\xi' / \textrm{arcmin}",
		ylabel = L"absolute radial velocity / km s$^{-1}$",
		#aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false
	)

	p = scatter!(xi_rot,  memb_stars.radial_velocity_gsr,
		alpha=0.1
	)

	#Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ e834933b-ef6c-4347-a24d-787ee65bc9b5
let
	fig, ax = FigAxis(
		xlabel = L"\eta' / \textrm{degree}",
		ylabel = L"RV",
		#aspect=DataAspect()
	)

	p = scatter!(eta_rot, memb_stars.RV
	)

	#Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ eb16cb5b-562f-4ef5-8a3b-247664d47f52
r_max_tangent = 30

# ╔═╡ 36c37e0b-924d-4a7f-b3ab-81e79e27a059
tangent_bins = LinRange(-r_max_tangent, r_max_tangent, 20)

# ╔═╡ 8f555e41-58d4-4e16-8f4a-1c0f58589b08
outside_bins = abs.(memb_stars.xi) .> r_max_tangent .|| abs.(memb_stars.eta) .> r_max_tangent

# ╔═╡ 8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
rv_mean = 111.0

# ╔═╡ 6b59d6e9-833b-4582-b5fb-f0a1f69a16c1
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		aspect=DataAspect()
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
		aspect=DataAspect()
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
		colorrange=(50, 100)
		)
	Colorbar(fig[1, 2], p, label="GSR radial velocity / km/s",
)
	@savefig "vlos_xi_eta_hist"

	fig
end

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
	"gradient" => summary_gradient |> OrderedDict,
	"rell" => summary_Rell |> OrderedDict,
	"bf_gradient" => bf_gradient,
	"bf_rell" => bf_sigma_Rell,
)

# ╔═╡ d71f6eba-7d64-4212-91c3-707a664c6b0b
open("processed/mcmc_properties_$FIGSUFFIX.toml", "w") do f
	TOML.print(f, df_summaries)
end

# ╔═╡ Cell order:
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═50488b8f-6886-4191-8778-af66929f1445
# ╠═4827c99d-6326-42fb-a329-d56992f6d20d
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
# ╠═1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╟─55ce0f69-8a96-4bbb-a59f-ee6503624ea6
# ╠═3377f632-713d-4fec-84a9-b0211b02cb43
# ╠═31a6c2e4-538c-4adc-bbda-5043680b17f7
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═45d748e7-42a0-471e-b290-e3469c0aa283
# ╟─6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
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
# ╠═c1cf3a80-1c26-4417-9ddf-1046f3f5f608
# ╠═433d1d51-4b2f-4bc0-807a-813eb215231a
# ╠═b32a6521-dae5-4d2a-89ab-6e8b3bc43be3
# ╠═93ee52ad-2c6d-4eac-85cd-ad2a88a70428
# ╠═b983cbd3-a90f-4ed1-b0a7-fd21e903f81d
# ╠═527ccde2-64db-40f5-9b77-8898a92b3b6a
# ╠═820b82a1-63a4-4a83-8d6e-5659ef75ce2d
# ╠═924369b0-aff6-46e4-a0bf-c6766ff93cbf
# ╠═4f00f714-6b16-41c7-a3c8-f2743d6023bf
# ╠═137bbc9c-1e02-4f78-acd1-0e30cf982615
# ╠═3b7cb3ea-3684-4df9-b9e6-919fa11e9ffa
# ╠═f8cc1ebb-617a-46a0-ab22-2f84d75da2ea
# ╠═9158e64f-5a8d-48bf-8043-8e93d41bd626
# ╠═4a83224a-1f2a-4f7c-b670-6a694fc8826a
# ╠═d051635b-123c-45e2-94da-6511f179a2bb
# ╟─062994bc-fb0c-4f08-b7f7-7cc6714bad1e
# ╠═25d9799f-5e89-4ec3-995c-a7309efb798f
# ╠═cb440be7-1f0e-4a20-b15a-82b1820d1ced
# ╠═e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
# ╠═184b4a5d-cbab-44b2-9620-bf928ad81d0e
# ╠═e79d9d82-84cd-4bc5-9ea8-be6b07cacf6d
# ╠═f2aa15e2-824b-4f0e-8ad0-42611abe06da
# ╠═ad607c95-0aea-4d84-8f18-e5b929b9bfca
# ╠═65f7b7e0-2439-477e-8036-a28a6903955f
# ╠═96b41e0b-b9f3-41ba-bc56-f3febd8833d3
# ╠═cbf50d46-219b-489d-b359-6b473b85735d
# ╠═2b88915c-7222-44be-a483-b967ea131b80
# ╠═4d7e7b96-74a0-4066-88cf-739c043c7f47
# ╠═0ca7dc1b-3b41-4089-9c89-20c6e48213ea
# ╠═c48bb30d-1186-4940-b061-91f53e8335e1
# ╟─8555e608-53c1-40d3-b21e-413af8953c30
# ╠═af0d2050-b42e-4a7f-aabb-5b08d23381e9
# ╠═7a9d2e0f-4dcb-4b69-9f68-d43d6dde8bf2
# ╠═3195286d-d85f-43a3-aa25-dae2134f570b
# ╠═4a473039-79f0-4d77-aa0c-681e2fba4f4c
# ╠═88f2918e-e126-420a-96a2-5746a8010f73
# ╠═ed35eb68-74f7-4009-9b68-dfca2ea547af
# ╠═f523cf48-82bf-4b20-9d7c-215bbe10a193
# ╠═b827e765-646c-4928-9f66-c64e7a20539f
# ╠═6371d804-cc73-4ce1-9b36-79fa61780d75
# ╠═70a22ef4-2eb1-4094-94e3-5a13fb51b9e6
# ╠═183f572a-bc0f-435b-a656-2ee2a3057559
# ╠═d3fb7136-7600-4782-ba97-f2f785fb3c0a
# ╠═51c7c22e-405b-4a47-beb1-b0fcf3efa391
# ╠═39bb4bb9-8445-4ee1-a3af-639d8fa96f65
# ╠═3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
# ╠═7fce47a1-16be-4dad-afb5-0fb04bd91355
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
# ╠═03928d42-ee61-4480-82d3-a482842f8521
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═e05ec20a-3165-4360-866e-3e8cae8665e5
# ╠═f82d2ff7-7a7f-4520-811e-126f3f4f5349
# ╠═e50a1968-9f91-4751-a7a9-0a4768c9ab7e
# ╟─31aa8fc5-1415-4c44-9b92-a7d097181639
# ╠═2f5db664-fc99-41ac-a726-f21dd5d88ad4
# ╠═301ad1d5-83a1-418d-aa06-583f7212541d
# ╠═090acae4-1209-49e6-882c-20ac2c972dd5
# ╠═fd0a74a1-6513-4612-8181-745d5b7c3f4c
# ╠═3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
# ╠═bf8a9940-1c8b-4fdd-a231-1861721ea8e9
# ╠═14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
# ╠═d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
# ╠═dbea26bd-539f-4eda-a4de-ed00efdc93a9
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
# ╠═0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
# ╠═319bd778-7e17-4bd7-856f-d6785b287219
# ╟─82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═f2313731-5b83-42d6-b624-c618bfb0bb5c
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═687b79fb-b0ce-4a93-8c0e-83b77aec081e
# ╠═9c66756b-af33-48fb-9947-201d53cdeba3
# ╠═e834933b-ef6c-4347-a24d-787ee65bc9b5
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
# ╟─1acdc61b-fb5f-449d-ba86-46525c881a39
# ╠═0532010e-7832-44d4-a7ee-5a6f6ee9d7da
# ╠═ba6c51c9-6d8f-4eca-af33-c75d0a5a5b37
# ╠═a393eeca-c9a4-412f-ae86-35ee1aca4d51
# ╠═d71f6eba-7d64-4212-91c3-707a664c6b0b
