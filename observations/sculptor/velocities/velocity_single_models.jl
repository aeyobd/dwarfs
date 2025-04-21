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
study = "w09"

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

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 5fd838b6-e7b1-4449-8f7a-b26a77042f0f
⊕ = RVUtils.:⊕

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 8863e6f7-6aaa-4370-919f-fe6afd3847bf
θ_orbit = obs_properties["theta_pm_gsr"]

# ╔═╡ d476513e-6211-44a4-ade5-0ab184a9829d
j24 = read_fits("../processed/best_sample.fits")

# ╔═╡ e391b59d-b8d9-4726-8e4d-8476f6d62800
rv0 = obs_properties["radial_velocity"] .- RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ cd71ec9b-4b8c-47ca-a0c8-3c6cbc257927
σv = obs_properties["sigma_v"]

# ╔═╡ 03636f8d-0fc9-4f20-8a84-b92f9bd27ffb
f_sat = mean(j24.PSAT)

# ╔═╡ 5b752b84-39a7-483d-96dd-a9107e330308
rv_meas_all = read_fits("processed/rv_combined.fits")


# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
rv_meas = let 
	rv_meas = filter(x->!ismissing(x["RV_$study"]), rv_meas_all)

	
	rv_meas[:, :RV] = rv_meas[!, "RV_$study"]
	rv_meas[:, :RV_err] .= rv_meas[!, "RV_err_$study"]

	obs = [lguys.ICRS(ra=r.ra, dec=r.dec, 
		distance=obs_properties["distance"], 
						pmra=r.pmra, pmdec=r.pmdec, radial_velocity=r.RV)
		for r in eachrow(rv_meas)
	]
		
	obs_gsr = lguys.to_frame(lguys.transform.(lguys.GSR, obs))

	
	rv_meas[!, :pmra_gsr] = obs_gsr.pmra
	rv_meas[!, :pmdec_gsr] = obs_gsr.pmdec
	rv_meas[!, :radial_velocity_gsr] = obs_gsr.radial_velocity

	Δ_gsr = RVUtils.rv_gsr_shift.(rv_meas.ra, rv_meas.dec)
	
	rv_meas[:, :vz] .= rv_meas.radial_velocity_gsr .+ rv_meas.delta_rv

	rv_meas[:, :vz_err] .= rv_meas.RV_err .⊕ rv_meas.delta_rv_err
	
	rv_meas[:, :RV_sigma] .= rv_meas[!, "RV_sigma_$study"]
	rv_meas[:, :RV_count] .= rv_meas[!, "RV_count_$study"]
	rv_meas[:, :F_RV] .= rv_meas[:, "F_scatter_$study"]
	rv_meas[:, :F_match] .= rv_meas[:, "F_match_$study"]

	rv_meas[:, :L_RV_SAT] = RVUtils.L_RV_SAT.(rv_meas.vz, rv_meas.RV_err, rv0, σv)
	rv_meas[:, :L_RV_BKD] = RVUtils.L_RV_BKD.(rv_meas.vz, rv_meas.ra, rv_meas.dec)

	rv_meas[:, :PSAT_RV] = RVUtils.PSAT_RV(rv_meas, f_sat)
	rv_meas
end

# ╔═╡ 811088e4-1ae4-467d-907d-242db314e08c
md"""
## counts
"""

# ╔═╡ 196845bf-5959-4223-9e39-75afa8cb5f0b
sum(rv_meas.RV_count)

# ╔═╡ e93d0bc6-ba43-465b-bec2-e1b1b1e24d10
sum(rv_meas.F_RV)

# ╔═╡ 1166808e-972f-4867-8a96-f84315aba7f6
sum(rv_meas.F_match)

# ╔═╡ c0fb98ca-bf95-451f-bb96-29a5c0b30730
sum(rv_meas.F_match .& rv_meas.F_RV)

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership criteria
"""

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = rv_meas.PSAT_RV .> 0.2# ignore spatial 
	memb_filt .&= rv_meas.F_RV
	memb_filt .&= rv_meas.F_match

	memb_filt
end

# ╔═╡ 616f378a-3063-41fc-8fdd-513c0bae9803
sum(memb_filt)

# ╔═╡ 59395e95-e04b-4782-b7d6-ef1945249b3f


# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f
md"""
## MCMC Priors
"""

# ╔═╡ 1b97d0e5-7a77-44a5-b609-ed8945cd959c


# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err), Prior(), 10000)
)

# ╔═╡ 55504311-fbbd-4351-b9e2-ecab053ba207
function plot_samples!(samples, x;
		thin=10, color=:black, alpha=nothing, kwargs...)

	alpha = 1 / (size(samples, 1))^(1/3)
	for sample in eachrow(samples)[1:thin:end]
		y = lguys.gaussian.(x, sample.μ, sample.σ)
		lines!(x, y, color=color, alpha=alpha)
	end
end

# ╔═╡ ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
pairplot(prior_samples[:, [:μ, :σ]])

# ╔═╡ 46d4bc0b-e6f7-455a-b491-eb85f0cd0e64


# ╔═╡ 74ad07df-f15c-436f-b390-ce95b27f7fab
let
	fig, ax = FigAxis(
		xgridvisible=false,
		ygridvisible=false,
		xlabel=xlabel,
		ylabel="density",
		title="priors"
	)
	
	plot_samples!(prior_samples, LinRange(-200, 200, 100))

	fig
end

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

# ╔═╡ 149b31b4-c451-42ff-88ec-6a01b4695e55
scatter(memb_stars.RV, memb_stars.vz .- memb_stars.radial_velocity_gsr, color=memb_stars.xi)

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ e93c66f7-394a-46bf-96c9-475494979548
memb_stars

# ╔═╡ 21a71cf7-efba-4b48-b280-92d6f1ea5d6d
begin 
	df_out = copy(memb_stars)
	for colname in names(df_out)
		col = df_out[!, colname]
		if eltype(col) === Any
			col = map(col) do x
				if x isa String
					x
				elseif isnan(x)
					""
				else
					error()
				end
			end
			df_out[!, colname] = col
		end
	end
	df_out
end

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	plot_samples!(samples, LinRange(70, 150, 100), thin=15)
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

# ╔═╡ d321f8ac-1044-45ec-8e1c-a2d8395b6917
# ╠═╡ disabled = true
#=╠═╡
icrs = lguys.ICRS(ra=obs_properties["ra"], dec=obs_properties["dec"], pmdec=obs_properties["pmdec"], pmra=obs_properties["pmra"], distance=obs_properties["distance"], radial_velocity=obs_properties["radial_velocity"])
  ╠═╡ =#

# ╔═╡ 1bc7adb7-fe85-4878-9527-c5d15dc761b1
Δv = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 212cee31-a04d-47c1-b56e-4a7b06322b99
median(samples_gsr.μ) + Δv

# ╔═╡ 0eafe553-4434-4b2b-9646-912ac99beccf
median(samples_gsr.μ .+ Δv)

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
	
	plot_samples!(samples_gsr, LinRange(40, 110, 100), thin=15)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 062994bc-fb0c-4f08-b7f7-7cc6714bad1e
md"""
## Gradient
"""

# ╔═╡ 42c24e51-77d2-4c3e-a367-0f52c9e52d0d
"""
Fits a normal (gaussian) distribution to 3d data with errors (to include spatial gradient)
"""
@model function normal_gradient(x, xerr, ξ, η; μ_min=90, μ_max=120)
	μ ~ Uniform(μ_min, μ_max)
	A ~ Normal(0, 0.1)
	B ~ Normal(0, 0.1)
	m = @. μ + A*ξ + B*η
	σ ~ LogNormal(2.5, 1) # approx 1 - 100 km / s, very broad but should cover all
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(m, s)
end

# ╔═╡ cb440be7-1f0e-4a20-b15a-82b1820d1ced
samples_gradient = sample(normal_gradient(memb_stars.vz, memb_stars.RV_err, memb_stars.xi, memb_stars.eta, μ_min=30, μ_max=100), 
						  NUTS(0.65), MCMCThreads(), 2000, 16)

# ╔═╡ e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
@savefig "gradient_corner" pairplot(samples_gradient)

# ╔═╡ d2143947-1cd2-4ec6-b88e-ef6e0caa2c78
sqrt((0.073/0.017)^2 + (0.043 / 0.026)^2)

# ╔═╡ 184b4a5d-cbab-44b2-9620-bf928ad81d0e
df_gradient = DataFrame(samples_gradient)

# ╔═╡ f506e011-6b13-4529-b6d2-305e1b00205d
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ 2e904b39-5b9b-4b9f-9497-7aabd648f7b2
vec_pm = LilGuys.pm2kms.([icrs0.pmra, icrs0.pmdec], icrs0.distance) /(180/π)

# ╔═╡ 0ca7dc1b-3b41-4089-9c89-20c6e48213ea
@savefig "v_gradient_derived" let
	fig = Figure()

	ax=Axis(fig[1,1];
		  limits=(-12, 12, -12, 12), 
		  aspect=DataAspect(),
		  xlabel = L"$\partial \,v_z / \partial\,\xi$ / km\,s$^{-1}$\,degree$^{-1}$",
		  ylabel = L"$\partial \,v_z / \partial\,\eta$ / km\,s$^{-1}$\,degree$^{-1}$",
		xreversed=true
		 )
	
	scatter!(60df_gradient.A, 60df_gradient.B, alpha=0.1, markersize=1, 

	   )

	scatter!(0, 0, color=:black)
	arrows!([0], [0], [vec_pm[1]], [vec_pm[2]])

	fig
end

# ╔═╡ d6a6ee30-ae1f-461c-a67b-c4b599246120
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ ceb1d7cb-2430-415d-98ea-a0b5ef1a4a42
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

# ╔═╡ 0bab042f-ecc5-46be-8a68-d510e650862c
# ╠═╡ disabled = true
#=╠═╡
pm_gsr_induced = lguys.transform(lguys.ICRS, lguys.GSR(ra=icrs0.ra, dec=icrs0.dec, pmra=0, pmdec=0, radial_velocity=0, distance=icrs0.distance, ))
  ╠═╡ =#

# ╔═╡ 2a422e88-fc0d-4a89-a841-42f3c5c8dace
import KernelDensity

# ╔═╡ f523cf48-82bf-4b20-9d7c-215bbe10a193
kde = KernelDensity.kde((df_gradient.A, df_gradient.B))

# ╔═╡ 2b88915c-7222-44be-a483-b967ea131b80
log(pdf(kde, 0, 0) ./ pdf(MvNormal([0,0], 0.1), [0,0]))

# ╔═╡ 333d9066-3943-4a7c-9a49-efb8af24ee0f
pdf(MvNormal([0,0], [0.1, 0.1]), [0,0])

# ╔═╡ b41d9be3-5da3-4a64-82c8-60b536c2f166
 lguys.gaussian(0, 0., 0.1)^2

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
	h = histogram(Float64.(memb_stars.vz), 30, normalization=:pdf)
	
	plot_samples!(DataFrame(samples_gradient), LinRange(40, 110, 100), thin=15)
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

# ╔═╡ c50f68d7-74c3-4c36-90c5-a5262982ed9f


# ╔═╡ f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
df_r_ell_z = calc_binned_mu_sigma(memb_stars.R_ell, memb_stars.vz, memb_stars.RV_err, bins)

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
df_xi_p = calc_binned_mu_sigma(memb_stars.xi_p, memb_stars.vz, memb_stars.vz_err, bins_equal_number(memb_stars.xi_p, n=10))

# ╔═╡ 090acae4-1209-49e6-882c-20ac2c972dd5
df_eta_p = calc_binned_mu_sigma(memb_stars.eta_p, memb_stars.vz, memb_stars.vz_err, bins_equal_number(memb_stars.eta_p, n=10))

# ╔═╡ fd0a74a1-6513-4612-8181-745d5b7c3f4c
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degrees",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(xi_p, memb_stars.radial_velocity_gsr, color=COLORS[3], alpha=0.1)

	scatter_range!(df_xi_p)

	fig
end

# ╔═╡ 3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
let
	fig, ax = FigAxis(
		xlabel = L"$\xi$ / arcmin",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x, df_xi_p.μ, yerror=df_xi_p.μ_err, color=:black)

	fig
end

# ╔═╡ 14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
let
	fig, ax = FigAxis(
		xlabel = L"$\eta'$ / arcmin",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(60df_eta_p.x, df_eta_p.μ, yerror=df_eta_p.μ_err, color=:black)

	fig
end

# ╔═╡ d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$",
		ylabel = L"$\sigma_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x, df_xi_p.σ, yerror=df_xi_p.σ_err, color=:black)

	fig
end

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ 0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
σ_m = median(samples.σ)

# ╔═╡ f82d2ff7-7a7f-4520-811e-126f3f4f5349
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errorscatter!(midpoints(bins), df_r_ell_z.σ, yerror=df_r_ell_z.σ_err, xerror=bin_errs, color=:black)
	hlines!(σ_m)

	fig
end

# ╔═╡ 319bd778-7e17-4bd7-856f-d6785b287219
quantile(samples.σ, [0.16, 0.84]) .- σ_m

# ╔═╡ 24ae8277-9644-40e5-b2ab-f4fc9584823c


# ╔═╡ 82a0e58a-30a4-4e42-b9c1-cb184eb551aa
md"""
# Misc plots
"""

# ╔═╡ 9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
not = !

# ╔═╡ f2313731-5b83-42d6-b624-c618bfb0bb5c
rv_mean_gsr = median(samples_gsr.μ)

# ╔═╡ d688d2e5-faca-4b14-801a-d58b08fd6654
let
	dr = 10
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
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

# ╔═╡ eb16cb5b-562f-4ef5-8a3b-247664d47f52
r_max_tangent = 0.4

# ╔═╡ 36c37e0b-924d-4a7f-b3ab-81e79e27a059
tangent_bins = LinRange(-r_max_tangent, r_max_tangent, 20)

# ╔═╡ 8f555e41-58d4-4e16-8f4a-1c0f58589b08
outside_bins = abs.(memb_stars.xi) .> r_max_tangent .|| abs.(memb_stars.eta) .> r_max_tangent

# ╔═╡ 8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
rv_mean = rv_mean_gsr

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
		colorrange=(50, 100)
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
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═35c87efe-6899-444d-952f-94e5e5342846
# ╠═5fd838b6-e7b1-4449-8f7a-b26a77042f0f
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═8863e6f7-6aaa-4370-919f-fe6afd3847bf
# ╠═d476513e-6211-44a4-ade5-0ab184a9829d
# ╠═e391b59d-b8d9-4726-8e4d-8476f6d62800
# ╠═cd71ec9b-4b8c-47ca-a0c8-3c6cbc257927
# ╠═03636f8d-0fc9-4f20-8a84-b92f9bd27ffb
# ╠═5b752b84-39a7-483d-96dd-a9107e330308
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╠═811088e4-1ae4-467d-907d-242db314e08c
# ╠═196845bf-5959-4223-9e39-75afa8cb5f0b
# ╠═e93d0bc6-ba43-465b-bec2-e1b1b1e24d10
# ╠═1166808e-972f-4867-8a96-f84315aba7f6
# ╠═c0fb98ca-bf95-451f-bb96-29a5c0b30730
# ╠═616f378a-3063-41fc-8fdd-513c0bae9803
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═59395e95-e04b-4782-b7d6-ef1945249b3f
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═1b97d0e5-7a77-44a5-b609-ed8945cd959c
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═55504311-fbbd-4351-b9e2-ecab053ba207
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═46d4bc0b-e6f7-455a-b491-eb85f0cd0e64
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═9380d9d1-58b2-432d-8528-d247cf5724e9
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═149b31b4-c451-42ff-88ec-6a01b4695e55
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═e93c66f7-394a-46bf-96c9-475494979548
# ╠═21a71cf7-efba-4b48-b280-92d6f1ea5d6d
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═d0ae48e2-8389-4641-b311-cfb4944b0851
# ╠═49e7d483-1dd1-4406-9bce-58d6e4412e7b
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═212cee31-a04d-47c1-b56e-4a7b06322b99
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═0eafe553-4434-4b2b-9646-912ac99beccf
# ╠═d321f8ac-1044-45ec-8e1c-a2d8395b6917
# ╠═1bc7adb7-fe85-4878-9527-c5d15dc761b1
# ╠═931ed52e-5e7a-4692-b568-ae26ea44b638
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═062994bc-fb0c-4f08-b7f7-7cc6714bad1e
# ╠═42c24e51-77d2-4c3e-a367-0f52c9e52d0d
# ╠═cb440be7-1f0e-4a20-b15a-82b1820d1ced
# ╠═e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
# ╠═d2143947-1cd2-4ec6-b88e-ef6e0caa2c78
# ╠═184b4a5d-cbab-44b2-9620-bf928ad81d0e
# ╠═2e904b39-5b9b-4b9f-9497-7aabd648f7b2
# ╠═0ca7dc1b-3b41-4089-9c89-20c6e48213ea
# ╠═f506e011-6b13-4529-b6d2-305e1b00205d
# ╠═d6a6ee30-ae1f-461c-a67b-c4b599246120
# ╠═ceb1d7cb-2430-415d-98ea-a0b5ef1a4a42
# ╠═0bab042f-ecc5-46be-8a68-d510e650862c
# ╠═2a422e88-fc0d-4a89-a841-42f3c5c8dace
# ╠═f523cf48-82bf-4b20-9d7c-215bbe10a193
# ╠═2b88915c-7222-44be-a483-b967ea131b80
# ╠═333d9066-3943-4a7c-9a49-efb8af24ee0f
# ╠═b41d9be3-5da3-4a64-82c8-60b536c2f166
# ╠═b827e765-646c-4928-9f66-c64e7a20539f
# ╠═6371d804-cc73-4ce1-9b36-79fa61780d75
# ╠═70a22ef4-2eb1-4094-94e3-5a13fb51b9e6
# ╠═183f572a-bc0f-435b-a656-2ee2a3057559
# ╠═d3fb7136-7600-4782-ba97-f2f785fb3c0a
# ╠═3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═c50f68d7-74c3-4c36-90c5-a5262982ed9f
# ╠═f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═e05ec20a-3165-4360-866e-3e8cae8665e5
# ╠═f82d2ff7-7a7f-4520-811e-126f3f4f5349
# ╟─31aa8fc5-1415-4c44-9b92-a7d097181639
# ╠═2f5db664-fc99-41ac-a726-f21dd5d88ad4
# ╠═090acae4-1209-49e6-882c-20ac2c972dd5
# ╠═fd0a74a1-6513-4612-8181-745d5b7c3f4c
# ╠═3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
# ╠═14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
# ╠═d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
# ╠═0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
# ╠═319bd778-7e17-4bd7-856f-d6785b287219
# ╠═24ae8277-9644-40e5-b2ab-f4fc9584823c
# ╟─82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═f2313731-5b83-42d6-b624-c618bfb0bb5c
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═eb16cb5b-562f-4ef5-8a3b-247664d47f52
# ╠═36c37e0b-924d-4a7f-b3ab-81e79e27a059
# ╠═8f555e41-58d4-4e16-8f4a-1c0f58589b08
# ╠═8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
# ╠═53da82d5-2e69-4f88-8043-6694d57cdd91
# ╠═b5533db0-a734-4d37-9d75-24471634f855
# ╠═4eea0a17-257e-4d0e-88df-9ff4858771b1
