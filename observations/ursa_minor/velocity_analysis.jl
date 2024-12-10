### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 72f1febc-c6ea-449a-8cec-cd0e49c4e20c
using DataFrames

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
using PairPlots

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ fe953ab7-1c13-4bd9-acfd-ec033a481e3d
rv_high = -210

# ╔═╡ 49e9aefb-a1cb-4211-8dc3-2344ef3427e5
rv_low = -280

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("observed_properties.toml")

# ╔═╡ 9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
θ_orbit = obs_properties["theta_pm_gsr"]

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
begin
	rv_meas = lguys.read_fits(joinpath(data_dir, "umi_averaged_rv.fits"))

	rv_meas[!, :xi_p], rv_meas[!, :eta_p] = lguys.to_orbit_coords(rv_meas.ra, rv_meas.dec, obs_properties["ra"], obs_properties["dec"], θ_orbit)
	rv_meas

end

# ╔═╡ 4c1a5aac-4bad-4eba-aa61-ccd317113633
fig_dir = "./figures/"

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = rv_meas.PSAT .> 0.2
	memb_filt .&= rv_high .> rv_meas.RV .> rv_low

	quality_score = rv_meas.RV_std ./ rv_meas.RV_err ./ sqrt.(rv_meas.RV_count)
	quality_score[isnan.(rv_meas.RV_std)] .= 0.

	memb_filt .&= quality_score .< 5
end

# ╔═╡ 4233a458-db79-418c-ae73-f2090be632d4
quality_score

# ╔═╡ 74b10a3e-1342-454f-8eed-77b371f81edf
hist(quality_score[isfinite.(quality_score)], bins=100, axis=(;yscale=log10, xlabel="quality score", ylabel="count"))

# ╔═╡ d8800a31-1ed3-422f-ac51-90f18cf61c29
sum(quality_score .< Inf)

# ╔═╡ 015f2345-74d4-4296-b0dc-35d6e13ad8cc
sum(quality_score .< 5)

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
"""
Fits a normal (gaussian) distribution to 1d data with errors.
"""
@model function normal_dist(x, xerr; μ_min=rv_low, μ_max=rv_high)
	μ ~ Uniform(μ_min, μ_max)
	σ ~ LogNormal(2.5, 1) # approx 1 - 100 km / s, very broad but should cover all
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(fill(μ, length(x)), s)
end

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(normal_dist(memb_stars.RV, memb_stars.RV_err), Prior(), 10000)
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

# ╔═╡ 74ad07df-f15c-436f-b390-ce95b27f7fab
let
	fig, ax = FigAxis(
		xgridvisible=false,
		ygridvisible=false,
		xlabel=xlabel,
		ylabel="density",
		title="priors"
	)
	
	plot_samples!(prior_samples, LinRange(rv_low, rv_high, 100))

	fig
end

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; μ_min=rv_low, μ_max=rv_high, N=3_000, p=0.16, burn=0.2)
	samples = DataFrame(sample(normal_dist(rv, rv_err, μ_min=μ_min, μ_max=μ_max), NUTS(0.65), N))
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
samples = DataFrame(sample(normal_dist(memb_stars.RV, memb_stars.RV_err), NUTS(0.65), 10000))

# ╔═╡ 149b31b4-c451-42ff-88ec-6a01b4695e55
scatter(memb_stars.RV_gsr, memb_stars.RV .- memb_stars.RV_gsr, color=memb_stars.xi_p)

# ╔═╡ b18e4622-41e0-4700-9e4b-3dbebeefea53
describe(samples)

# ╔═╡ b0a99549-0183-4893-9ba9-0435a5b499f4


# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ e93c66f7-394a-46bf-96c9-475494979548
memb_stars

# ╔═╡ a162219f-df5a-41d8-bf54-927a355f6431
lguys.write_fits(joinpath(data_dir, "sculptor_memb_rv.fits"), memb_stars)

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	plot_samples!(samples, LinRange(rv_low, rv_high, 100), thin=15)
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 61e15c47-c454-48db-95f9-02abe052676e
mean(memb_stars.RV)

# ╔═╡ d5938fc3-9c8a-4e33-8401-500b4201df14
sem(memb_stars.RV)

# ╔═╡ d0ae48e2-8389-4641-b311-cfb4944b0851
std(memb_stars.RV)

# ╔═╡ 49e7d483-1dd1-4406-9bce-58d6e4412e7b
mean(memb_stars.RV_gsr)

# ╔═╡ f56b1046-e63c-41dd-abcd-44f9e19b5cad
hist(memb_stars.RV_gsr[isfinite.(memb_stars.RV_gsr)])

# ╔═╡ 2f12d15b-20a0-411a-97e2-4aa38338d1b9
memb_best = memb_stars[isfinite.(memb_stars.VZ), :]

# ╔═╡ 081a4ece-a4f7-467d-a763-17668f7c1b82
samples_best = DataFrame(sample(normal_dist(memb_best.RV, memb_best.RV_err), NUTS(0.65), 10000))

# ╔═╡ 7107bc46-b93d-44e0-92fd-cf5ddf978b5e
pairplot(samples_best[:, [:μ, :σ]])

# ╔═╡ 9ffcc0b0-11a1-4962-b64e-59a731f22bd8
samples_gsr = DataFrame(sample(normal_dist(memb_best.RV_gsr, memb_best.RV_err, μ_min=-110, μ_max=-40), NUTS(0.65), 10000))

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
pairplot(samples_gsr[:, [:μ, :σ]])

# ╔═╡ 50abc283-02f6-41d6-a6ec-57ab5a19657a
std(memb_stars.RV)

# ╔═╡ 9129c4bb-8d0c-4680-8cdb-0e5abfffbfb6
std(memb_best.RV)

# ╔═╡ 43222bee-c3ec-4004-9745-29f5194a2ae4
std(memb_best.VZ)

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
bins = bins_equal_number(memb_stars.r_ell, n=10)

# ╔═╡ c50f68d7-74c3-4c36-90c5-a5262982ed9f
df_r_ell = calc_binned_mu_sigma(memb_stars.r_ell, memb_stars.RV, memb_stars.RV_err, bins)

# ╔═╡ f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
df_r_ell_z = calc_binned_mu_sigma(memb_best.r_ell, memb_best.VZ, memb_best.VZ_err, bins, μ_min=-110, μ_max=-40)

# ╔═╡ 45d217cc-bc7e-463a-a528-229b16e9d112
hist(memb_best.VZ)

# ╔═╡ 38da4da1-74f5-4661-89f2-4b25562a1faf
function scatter_range!(df_r_ell)
	errscatter!(df_r_ell.x, df_r_ell.μ, yerr=df_r_ell.μ_err, color=:black)
	
	errorbars!(df_r_ell.x, df_r_ell.μ .+ df_r_ell.σ, df_r_ell.x .- df_r_ell.x_low, df_r_ell.x_high .- df_r_ell.x,  direction = :x, color=:black)
	errorbars!(df_r_ell.x, df_r_ell.μ .- df_r_ell.σ, df_r_ell.x .- df_r_ell.x_low, df_r_ell.x_high .- df_r_ell.x, direction = :x, color=:black)
end

# ╔═╡ 86776e68-d47f-43ed-b37f-432c864050bb
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.r_ell, memb_stars.RV, color=COLORS[3], alpha=0.1)

	scatter_range!(df_r_ell)

	fig
end

# ╔═╡ e05ec20a-3165-4360-866e-3e8cae8665e5
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.r_ell, memb_stars.VZ, color=COLORS[3], alpha=0.1)

	scatter_range!(df_r_ell_z)

	fig
end

# ╔═╡ 31aa8fc5-1415-4c44-9b92-a7d097181639
md"""
# Binned properties with radius
"""

# ╔═╡ 4825643b-d904-4231-8f44-8cd934f23795
vgsr_low = -110

# ╔═╡ 51f25214-45bc-44e8-82c3-447c3d8dba09
vgsr_high = -40

# ╔═╡ 6a778219-f65a-4d81-9fc8-2f98a54d05c2
samples_vz = DataFrame(sample(normal_dist(memb_best.VZ, memb_best.VZ_err, μ_min=vgsr_low, μ_max=vgsr_high), NUTS(0.65), 10000))

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV_gsr), 30, normalization=:pdf)
	
	plot_samples!(samples_vz, LinRange(-110, -40, 100), thin=15)
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 327d3bc0-67ee-42bd-9eb1-7702431cde2a
pairplot(samples_vz[:, [:μ, :σ]])

# ╔═╡ 2f5db664-fc99-41ac-a726-f21dd5d88ad4
df_xi_p = calc_binned_mu_sigma(memb_best.xi_p, memb_best.VZ, memb_best.VZ_err, bins_equal_number(memb_best.xi_p, n=10), μ_min=vgsr_low, μ_max=vgsr_high)

# ╔═╡ 090acae4-1209-49e6-882c-20ac2c972dd5
df_eta_p = calc_binned_mu_sigma(memb_best.eta_p, memb_best.VZ, memb_best.VZ_err, bins_equal_number(memb_best.eta_p, n=10), μ_min=vgsr_low, μ_max=vgsr_high)

# ╔═╡ fd0a74a1-6513-4612-8181-745d5b7c3f4c
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degrees",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.xi_p, memb_stars.VZ, color=COLORS[3], alpha=0.1)

	scatter_range!(df_xi_p)

	fig
end

# ╔═╡ 3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / arcmin",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errscatter!(60df_xi_p.x, df_xi_p.μ, yerr=df_xi_p.μ_err, color=:black)

	fig
end

# ╔═╡ 14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
let
	fig, ax = FigAxis(
		xlabel = L"$\eta'$ / arcmin",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errscatter!(60df_eta_p.x, df_eta_p.μ, yerr=df_eta_p.μ_err, color=:black)

	fig
end

# ╔═╡ d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$",
		ylabel = L"$\sigma_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errscatter!(df_xi_p.x, df_xi_p.σ, yerr=df_xi_p.σ_err, color=:black)

	fig
end

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ 0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
σ_m = median(samples.σ)

# ╔═╡ 614f3f09-1880-491d-b41e-4e229330d66f
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errscatter!(df_r_ell.x, df_r_ell.σ, yerr=df_r_ell.σ_err, color=:black)
	hlines!(σ_m)

	fig
end

# ╔═╡ f82d2ff7-7a7f-4520-811e-126f3f4f5349
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errscatter!(midpoints(bins), df_r_ell_z.σ, yerr=df_r_ell_z.σ_err, xerr=bin_errs, color=:black)
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

# ╔═╡ 3a69f395-3c2d-4357-89af-5963d5fa79b8
let
	fig, ax = FigAxis(
		xlabel=L"$\log r_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)
	
	scatter!(log10.(memb_stars.r_ell), memb_stars.RV)


	fig
end

# ╔═╡ 9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
not = !

# ╔═╡ 30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
sum(not.(ismissing.(memb_stars.RV_s24)))

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

	p = scatter!(memb_stars.xi, memb_stars.eta, color=memb_stars.RV_gsr,
		colorrange=(rv_mean_gsr - dr, rv_mean_gsr + dr),
		colormap=:bluesreds,
		markersize = 600 ./ memb_stars.RV
	)

	Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ f72b9d3c-b165-4a2b-b008-d77ffd7d59d1
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		aspect=DataAspect()
	)

	p = scatter!(memb_stars.xi, memb_stars.eta, color=memb_stars.eta_p,
		colormap=:bluesreds,
		markersize = 600 ./ memb_stars.RV
	)

	Colorbar(fig[1, 2], p, label=L"\eta'")
	fig
end

# ╔═╡ 687b79fb-b0ce-4a93-8c0e-83b77aec081e
let
	fig, ax = FigAxis(
		xlabel = L"\xi' / \textrm{degree}",
		ylabel = L"RV",
		#aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false
	)

	p = scatter!(memb_stars.xi_p, memb_stars.RV,
		alpha=0.1
	)

	#Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ 9c66756b-af33-48fb-9947-201d53cdeba3
let
	fig, ax = FigAxis(
		xlabel = L"\xi' / \textrm{degree}",
		ylabel = L"absolute radial velocity / km s$^{-1}$",
		#aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false
	)

	p = scatter!(memb_stars.xi_p,  memb_stars.RV_gsr,
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

	p = scatter!(memb_stars.eta_p, memb_stars.RV
	)

	#Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s")
	fig
end

# ╔═╡ eb16cb5b-562f-4ef5-8a3b-247664d47f52
r_max_tangent = 0.4

# ╔═╡ 36c37e0b-924d-4a7f-b3ab-81e79e27a059
tangent_bins = LinRange(-r_max_tangent, r_max_tangent, 20)

# ╔═╡ 8f555e41-58d4-4e16-8f4a-1c0f58589b08
outside_bins = abs.(memb_stars.xi) .> r_max_tangent .|| abs.(memb_stars.eta) .> r_max_tangent

# ╔═╡ b94fce57-e85b-4991-9026-235cc0cb7b9b
outside_bins_best = abs.(memb_best.xi) .> r_max_tangent .|| abs.(memb_best.eta) .> r_max_tangent

# ╔═╡ 8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
rv_mean = mean(samples.μ)

# ╔═╡ f22576e4-2ea0-4ade-b87d-50c811bb604f
let
	fig, ax = FigAxis(
		xlabel=L"$\Delta$ radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV) .- rv_mean, 30, normalization=:pdf)
	scatter!(midpoints(h.bins), h.values,  label="ICRS")
	
	h = histogram(Float64.(memb_best.RV) .- rv_mean, 30, normalization=:pdf)
	scatter!(midpoints(h.bins), h.values,  label="ICRS")

	
	h = histogram(Float64.(memb_stars.RV_gsr) .- rv_mean_gsr, 30, normalization=:pdf)
	scatter!(midpoints(h.bins), h.values, label="GSR")

	h = histogram(Float64.(memb_stars.VZ) .- rv_mean_gsr, 30, normalization=:pdf)
	scatter!(midpoints(h.bins), h.values, label="VZ")



	axislegend()
	fig
end

# ╔═╡ be8712b3-8b4f-42cd-8293-2a4c6cf370d9
obs_properties

# ╔═╡ 319b4678-e7fa-4ec4-a1b9-8e6ccc786e94
gc_cen = lguys.transform(lguys.GSR, lguys.ICRS(ra=obs_properties["ra"], dec=obs_properties["dec"], distance=obs_properties["distance"], pmra=obs_properties["pmra"], pmdec=obs_properties["pmdec"], radial_velocity=rv_mean))

# ╔═╡ e479a40b-17e7-422b-970c-1876393475af
gc_cen.radial_velocity

# ╔═╡ 1ae62158-7b52-4761-a67a-10c754ea18da
rv_mean_gsr

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
		colorrange=(rv_mean - 20, rv_mean + 20)
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
	memb_stars = memb_best

	w = 1 ./ memb_stars.VZ_err .^ 2
	rv_mean = mean(memb_stars.VZ, weights(w))

	k1 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins, weights= w .* memb_stars.VZ)
	k2 = Arya.histogram2d(memb_stars.xi, memb_stars.eta, tangent_bins, weights=w)

	k1.values ./= k2.values

	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		colorrange=(rv_mean - 20, rv_mean + 20)
		)

	scatter!(memb_stars.xi[outside_bins_best], memb_stars.eta[outside_bins_best],
		color = memb_stars.VZ[outside_bins_best],
		colormap=:bluesreds,
		colorrange=(rv_mean - 20, rv_mean + 20)
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
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		aspect=DataAspect(),
		xreversed=true,
	)

	w = 1 ./ memb_stars.VZ_err

	bins = (33, 25)
	k1 = Arya.histogram2d(memb_stars.xi,memb_stars.eta, bins, weights= w .* memb_stars.VZ)
	k2 = Arya.histogram2d(memb_stars.xi,memb_stars.eta, bins, weights=w)

	k1.values ./= k2.values

	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		colorrange=(vgsr_low, vgsr_high)
		)
	Colorbar(fig[1, 2], p, label="GSR radial velocity / km/s",
)

	fig
end

# ╔═╡ 7178e5b9-cc42-4933-970a-4707ba69dbe9
let
	fig, ax = FigAxis()

	memb_stars = memb_best
	
	p = arrows!(memb_stars.ra, memb_stars.dec, memb_stars.pmra_gsr, memb_stars.pmdec_gsr, 				
		color=memb_stars.RV_gsr,
		colorrange=(vgsr_low, vgsr_high),
		colormap=:bluesreds,
		lengthscale=0.1
	)

	#Colorbar(fig[1, 2], p)
	fig
end


# ╔═╡ Cell order:
# ╠═6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═fe953ab7-1c13-4bd9-acfd-ec033a481e3d
# ╠═49e9aefb-a1cb-4211-8dc3-2344ef3427e5
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╠═9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═4c1a5aac-4bad-4eba-aa61-ccd317113633
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═4233a458-db79-418c-ae73-f2090be632d4
# ╠═74b10a3e-1342-454f-8eed-77b371f81edf
# ╠═d8800a31-1ed3-422f-ac51-90f18cf61c29
# ╠═015f2345-74d4-4296-b0dc-35d6e13ad8cc
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═1b97d0e5-7a77-44a5-b609-ed8945cd959c
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═55504311-fbbd-4351-b9e2-ecab053ba207
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═9380d9d1-58b2-432d-8528-d247cf5724e9
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═149b31b4-c451-42ff-88ec-6a01b4695e55
# ╠═b18e4622-41e0-4700-9e4b-3dbebeefea53
# ╠═b0a99549-0183-4893-9ba9-0435a5b499f4
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═e93c66f7-394a-46bf-96c9-475494979548
# ╠═a162219f-df5a-41d8-bf54-927a355f6431
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═081a4ece-a4f7-467d-a763-17668f7c1b82
# ╠═7107bc46-b93d-44e0-92fd-cf5ddf978b5e
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═d0ae48e2-8389-4641-b311-cfb4944b0851
# ╠═49e7d483-1dd1-4406-9bce-58d6e4412e7b
# ╠═f56b1046-e63c-41dd-abcd-44f9e19b5cad
# ╠═2f12d15b-20a0-411a-97e2-4aa38338d1b9
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═6a778219-f65a-4d81-9fc8-2f98a54d05c2
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═327d3bc0-67ee-42bd-9eb1-7702431cde2a
# ╠═f22576e4-2ea0-4ade-b87d-50c811bb604f
# ╠═50abc283-02f6-41d6-a6ec-57ab5a19657a
# ╠═9129c4bb-8d0c-4680-8cdb-0e5abfffbfb6
# ╠═43222bee-c3ec-4004-9745-29f5194a2ae4
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═c50f68d7-74c3-4c36-90c5-a5262982ed9f
# ╠═f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
# ╠═45d217cc-bc7e-463a-a528-229b16e9d112
# ╠═30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═86776e68-d47f-43ed-b37f-432c864050bb
# ╠═e05ec20a-3165-4360-866e-3e8cae8665e5
# ╠═614f3f09-1880-491d-b41e-4e229330d66f
# ╠═f82d2ff7-7a7f-4520-811e-126f3f4f5349
# ╟─31aa8fc5-1415-4c44-9b92-a7d097181639
# ╠═4825643b-d904-4231-8f44-8cd934f23795
# ╠═51f25214-45bc-44e8-82c3-447c3d8dba09
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
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═f2313731-5b83-42d6-b624-c618bfb0bb5c
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═f72b9d3c-b165-4a2b-b008-d77ffd7d59d1
# ╠═687b79fb-b0ce-4a93-8c0e-83b77aec081e
# ╠═9c66756b-af33-48fb-9947-201d53cdeba3
# ╠═e834933b-ef6c-4347-a24d-787ee65bc9b5
# ╠═eb16cb5b-562f-4ef5-8a3b-247664d47f52
# ╠═36c37e0b-924d-4a7f-b3ab-81e79e27a059
# ╠═8f555e41-58d4-4e16-8f4a-1c0f58589b08
# ╠═b94fce57-e85b-4991-9026-235cc0cb7b9b
# ╠═8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
# ╠═be8712b3-8b4f-42cd-8293-2a4c6cf370d9
# ╠═319b4678-e7fa-4ec4-a1b9-8e6ccc786e94
# ╠═e479a40b-17e7-422b-970c-1876393475af
# ╠═1ae62158-7b52-4761-a67a-10c754ea18da
# ╠═6b59d6e9-833b-4582-b5fb-f0a1f69a16c1
# ╠═89bf44ef-1ff5-443e-b8be-a3d1571b72e3
# ╠═53da82d5-2e69-4f88-8043-6694d57cdd91
# ╠═106482c9-a9f9-4b6e-95a9-614ab7991e23
# ╠═b5533db0-a734-4d37-9d75-24471634f855
# ╠═4eea0a17-257e-4d0e-88df-9ff4858771b1
# ╠═7178e5b9-cc42-4933-970a-4707ba69dbe9
