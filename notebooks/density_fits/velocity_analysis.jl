### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../../data/"

# ╔═╡ 9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
θ_orbit = -40.095

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("/astro/dboyea/dwarfs/sculptor_obs_properties.toml")

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
begin
	rv_meas = lguys.load_fits(joinpath(data_dir, "sculptor_all_rv.fits"))

	rv_meas[!, :xi_p], rv_meas[!, :eta_p] = lguys.to_orbit_coords(rv_meas.ra, rv_meas.dec, obs_properties["ra"], obs_properties["dec"], θ_orbit)


	obs = [lguys.ICRS(r.ra, r.dec, obs_properties["distance"], 
						r.pmra, r.pmdec, r.RV)
	for r in eachrow(rv_meas)]
		
	obs_gsr = lguys.to_frame(lguys.transform.(lguys.GSR, obs))

	
	rv_meas[!, :pmra_gsr] = obs_gsr.pmra
	rv_meas[!, :pmdec_gsr] = obs_gsr.pmdec
	rv_meas[!, :radial_velocity_gsr] = obs_gsr.radial_velocity

	rv_meas

end

# ╔═╡ 4c1a5aac-4bad-4eba-aa61-ccd317113633
fig_dir = "./figures/"

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = rv_meas.PSAT .> 0.2
	memb_filt .&= 150 .> rv_meas.RV .> 60

	quality_score = rv_meas.RV_std ./ rv_meas.RV_err
	quality_score[isnan.(rv_meas.RV_std)] .= 0

	memb_filt .&= quality_score .< 5
end

# ╔═╡ 74b10a3e-1342-454f-8eed-77b371f81edf
lines(histogram(quality_score, normalization=:none))

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
@model function normal_dist(x, xerr; μ_min=90, μ_max=120)
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
	
	plot_samples!(prior_samples, LinRange(80, 130, 100))

	fig
end

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; μ_min=90, μ_max=120, N=3_000, p=0.16, burn=0.2)
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
scatter(memb_stars.radial_velocity_gsr, memb_stars.RV .- memb_stars.radial_velocity_gsr, color=memb_stars.xi_p)

# ╔═╡ b18e4622-41e0-4700-9e4b-3dbebeefea53
describe(samples)

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
	
	plot_samples!(samples, LinRange(70, 150, 100), thin=15)
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
mean(memb_stars.radial_velocity_gsr)

# ╔═╡ 9ffcc0b0-11a1-4962-b64e-59a731f22bd8
samples_gsr = DataFrame(sample(normal_dist(memb_stars.radial_velocity_gsr, memb_stars.RV_err, μ_min=30, μ_max=100), NUTS(0.65), 10000))

# ╔═╡ 3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
pairplot(samples_gsr[:, [:μ, :σ]])

# ╔═╡ d321f8ac-1044-45ec-8e1c-a2d8395b6917
icrs = lguys.ICRS(ra=obs_properties["ra"], dec=obs_properties["dec"], pmdec=obs_properties["pmdec"], pmra=obs_properties["pmra"], distance=obs_properties["distance"], radial_velocity=obs_properties["radial_velocity"])

# ╔═╡ 931ed52e-5e7a-4692-b568-ae26ea44b638
lguys.transform(lguys.GSR, icrs)

# ╔═╡ 7514e306-ddc1-44a3-9242-5b12cf2a1536
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 30, normalization=:pdf)
	
	plot_samples!(samples_gsr, LinRange(40, 110, 100), thin=15)
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=COLORS[6])

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

	
	for i in 1:N
		filt = x .>= bins[i]
		filt .&= x .< bins[i+1]
		@info "calculating bin $i"

		μs[i], σs[i], μ_errs[i], σ_errs[i] = fit_rv_sigma(y[filt], yerr[filt]; kwargs...)
	end

	return DataFrame(
		x=midpoints(bins),
		x_err = diff(bins)/2,
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

# ╔═╡ 38da4da1-74f5-4661-89f2-4b25562a1faf
function scatter_range!(df_r_ell)
	errscatter!(df_r_ell.x, df_r_ell.μ, yerr=df_r_ell.μ_err, color=:black)
	
	errorbars!(df_r_ell.x, df_r_ell.μ .+ df_r_ell.σ, df_r_ell.x_err, direction = :x, color=:black)
	errorbars!(df_r_ell.x, df_r_ell.μ .- df_r_ell.σ, df_r_ell.x_err, direction = :x, color=:black)
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

# ╔═╡ 2f5db664-fc99-41ac-a726-f21dd5d88ad4
df_xi_p = calc_binned_mu_sigma(memb_stars.xi_p, memb_stars.radial_velocity_gsr, memb_stars.RV_err, bins_equal_number(memb_stars.xi_p, n=10), μ_min=50, μ_max=90)

# ╔═╡ fd0a74a1-6513-4612-8181-745d5b7c3f4c
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degrees",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.xi_p, memb_stars.radial_velocity_gsr, color=COLORS[3], alpha=0.1)

	scatter_range!(df_xi_p)

	fig
end

# ╔═╡ 3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errscatter!(df_xi_p.x, df_xi_p.μ, yerr=df_xi_p.μ_err, xerr=df_xi_p.x_err, color=:black)

	fig
end

# ╔═╡ d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$",
		ylabel = L"$\sigma_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errscatter!(df_xi_p.x, df_xi_p.σ, yerr=df_xi_p.σ_err, xerr=df_xi_p.x_err, color=:black)

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

	errscatter!(midpoints(bins), df_r_ell.σ, yerr=df_r_ell.σ_err, xerr=bin_errs, color=:black)
	hlines!(σ_m)

	fig
end

# ╔═╡ 319bd778-7e17-4bd7-856f-d6785b287219
quantile(samples.σ, [0.16, 0.84]) .- σ_m

# ╔═╡ 31aa8fc5-1415-4c44-9b92-a7d097181639
md"""
# Binned properties with radius
"""

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
sum(not.(ismissing.(memb_stars.RV_gmos)))

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

	p = scatter!(memb_stars.xi_p,  memb_stars.radial_velocity_gsr,
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
		colorrange=(rv_mean - 20, rv_mean + 20)
	)
	Colorbar(fig[1, 2], p, label="heliocentric radial velocity / km/s",
)
	fig
end

# ╔═╡ b5533db0-a734-4d37-9d75-24471634f855


# ╔═╡ 4eea0a17-257e-4d0e-88df-9ff4858771b1
let
	fig, ax = FigAxis(
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		aspect=DataAspect()
	)

	w = 1 ./ memb_stars.RV_err

	bins = (33, 25)
	k1 = Arya.histogram2d(memb_stars.xi,memb_stars.eta, bins, weights= w .* memb_stars.radial_velocity_gsr)
	k2 = Arya.histogram2d(memb_stars.xi,memb_stars.eta, bins, weights=w)

	k1.values ./= k2.values

	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		# colorrange=(91, 131)
		)
	Colorbar(fig[1, 2], p, label="absolute radial velocity / km/s",
)
	fig
end

# ╔═╡ 7178e5b9-cc42-4933-970a-4707ba69dbe9
let
	fig, ax = FigAxis()

	p = arrows!(memb_stars.ra, memb_stars.dec, memb_stars.pmra_gsr, memb_stars.pmdec_gsr, 				
		color=memb_stars.radial_velocity_gsr,
		colorrange=(50, 90),
		colormap=:bluesreds,
		lengthscale=0.1
	)

	#Colorbar(fig[1, 2], p)
	fig
end


# ╔═╡ Cell order:
# ╠═6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╠═9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═4c1a5aac-4bad-4eba-aa61-ccd317113633
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
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
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═e93c66f7-394a-46bf-96c9-475494979548
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═d0ae48e2-8389-4641-b311-cfb4944b0851
# ╠═49e7d483-1dd1-4406-9bce-58d6e4412e7b
# ╠═9ffcc0b0-11a1-4962-b64e-59a731f22bd8
# ╠═3f9dfb6e-e4d5-4895-b8ac-f5304dd15088
# ╠═d321f8ac-1044-45ec-8e1c-a2d8395b6917
# ╠═931ed52e-5e7a-4692-b568-ae26ea44b638
# ╠═7514e306-ddc1-44a3-9242-5b12cf2a1536
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═c50f68d7-74c3-4c36-90c5-a5262982ed9f
# ╠═30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═86776e68-d47f-43ed-b37f-432c864050bb
# ╠═614f3f09-1880-491d-b41e-4e229330d66f
# ╠═2f5db664-fc99-41ac-a726-f21dd5d88ad4
# ╠═fd0a74a1-6513-4612-8181-745d5b7c3f4c
# ╠═3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
# ╠═d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
# ╠═0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
# ╠═319bd778-7e17-4bd7-856f-d6785b287219
# ╟─31aa8fc5-1415-4c44-9b92-a7d097181639
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
# ╠═8bee2f6a-c65d-4f89-8a4e-5f5789a7b03d
# ╠═6b59d6e9-833b-4582-b5fb-f0a1f69a16c1
# ╟─89bf44ef-1ff5-443e-b8be-a3d1571b72e3
# ╠═b5533db0-a734-4d37-9d75-24471634f855
# ╠═4eea0a17-257e-4d0e-88df-9ff4858771b1
# ╠═7178e5b9-cc42-4933-970a-4707ba69dbe9
