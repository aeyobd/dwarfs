### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using LilGuys
end

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ e74e6a96-a18e-40bf-ade8-07ce0e30e5c4
using Distributions

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
include("../utils/gaia_filters.jl")

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	galaxy = "antlia2"
end

# ╔═╡ 40755684-2881-44f8-9fe6-09c211bbfe45
all_profiles = false

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 109b91c2-c0ad-4683-849b-a56da5766ee4
import StatsBase: sem

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(galaxy, "figures"); FIGSUFFIX=".mcmc_hist_fast"

# ╔═╡ 8bbd8c9f-2d35-47d9-8c30-0c86c14d6ce4


# ╔═╡ 57a19d65-59d9-46cf-8916-d9ac3a4dc92b
datafile = let
	dir = joinpath(galaxy, "data")

	filenames = ["jensen+24_1c.fits", "jensen+24_2c.fits", "jensen+24_wide.fits", "j24_1c.fits", "j24_2c.fits"]
	filename = ""
	for file in filenames
		if isfile(joinpath(dir, file))
			filename =  joinpath(dir, file)
		end
	end

	filename
end

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type="svg", pt_per_unit=2)

# ╔═╡ 2ab34a36-542c-41da-a83d-6eb67cae253c


# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,r_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ f078d929-62be-48f2-acdb-808357801a7b
md"""
note, we do change the first and last bins to be something more workable (just double the last bin)
"""

# ╔═╡ 4b77612f-e124-4d77-974c-d40c2f5a37ff
begin 
	bins = TOML.parsefile(joinpath(galaxy, "processed/info.mcmc_hist_fast.toml"))["bins"]
	bins[1] = bins[2] / 2
	bins[end] = bins[end-1] * 2

	bins
end

# ╔═╡ b9cf2c23-6c9a-44f5-9740-22caf4959831
obs_props = TOML.parsefile(joinpath(galaxy, "observed_properties.toml"))

# ╔═╡ c58ed939-07c1-4d34-90a4-3c5d9cc9910b
all_stars = read_gaia_stars(GaiaFilterParams(filename = datafile, obs_props))

# ╔═╡ 91867748-9b36-4f62-9310-8b778935776b
best_stars = filter(r->r.F_BEST == 1, all_stars)

# ╔═╡ 43d159c4-1599-4138-96b9-c9447a8d6315
md"""
The idea with this model is to use the likelihood

``
{\cal L} \propto {\rm Prior}(\theta) \prod_i (f(\theta, x_i) P_{\rm memb,\ \it i} + (1-f(\theta, x_i)) P_{\rm bg,\it\ i})
``

where $f$ is a model for the fraction of stars in a small region belonging to the satellite.
"""

# ╔═╡ e60ea082-818d-46f8-bf50-30a932f97ef2
prof_auto = LilGuys.StellarProfile(best_stars.R_ell[best_stars.PSAT .> 0.2], errors=:weighted, normalization=:none)

# ╔═╡ 1d9b3718-45d7-4765-ac7d-017dbf939ec8
md"""
# Histogram model
"""

# ╔═╡ e7bf7a13-afca-49b8-8771-f7914adb347b
function hist_fractions(bins, params, x)
	i = DE.bin_indices(x, bins)
	i = max.(i, 1)
	i = min.(i, length(bins)-1)
	return params[i]
end

# ╔═╡ 27c5a02f-0be6-4d0b-9c7f-99be12912732
pos_err = (sem(filter(r->r.PSAT .> 0.2, best_stars).xi) + sem(filter(r->r.PSAT .> 0.2, best_stars).eta))/2

# ╔═╡ 6f016a8e-38ae-4f05-a7ee-c292ac0e5741
df_chains = CSV.read(joinpath(galaxy, "processed", "samples.mcmc_hist_fast.csv"), DataFrame)

# ╔═╡ 33e02945-2ab5-4f55-b3d5-e8855d86f120
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		yscale = log10,
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = log_r_label,
		ylabel = "counts",
		limits=(nothing, nothing, 1, 3e4),
	)
	
	stephist!(log10.(best_stars.R_ell[best_stars.PSAT .> 0.2]), bins=log10.(bins[2:end-1]), label="PSAT > 0.2")
	stephist!(log10.(best_stars.R_ell), bins=log10.(bins), label="best")

	axislegend(position=:lt)
	fig
end

# ╔═╡ 64c31de5-41fa-47e8-9209-61901e05f1be
bins

# ╔═╡ 9ea0c184-4a09-4293-a6d7-a963cf1d3d58
maximum(best_stars.xi .⊕ best_stars.eta) * 60

# ╔═╡ 984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
r_max = maximum(best_stars.xi .⊕ best_stars.eta) * 60 * sqrt(1 - obs_props["ellipticity"])

# ╔═╡ 5feaa9b0-ca44-48d6-80d1-16303a34005f
prof_simple = LilGuys.StellarProfile(best_stars.R_ell[best_stars.PSAT .> 0.2], bins=log10.(bins)[bins .< r_max], errors=:weighted, normalization=:none)

# ╔═╡ 9bee7d00-d281-4d2d-ab9a-cbd3333f0e8e
sum(prof_simple.counts)

# ╔═╡ 82f20e46-90c7-4925-83c9-3f49a909664d
sum(prof_simple.counts)

# ╔═╡ ab124afe-c94a-44d9-8869-bd4d4cee3fbd
prof_weighted = LilGuys.StellarProfile(best_stars.R_ell, weights=best_stars.PSAT, bins=log10.(bins)[bins .< r_max], errors=:weighted,  normalization=:none)

# ╔═╡ babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
prof_bernoulli = LilGuys.StellarProfile(best_stars.R_ell, weights=best_stars.PSAT, bins=log10.(bins)[bins .< r_max], errors=:bernoulli,  normalization=:none)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ e87b5ee8-e281-43e3-8172-4848018e2cd2
function read_error(key)
	if key * "_em" ∈ keys(obs_props)
		return max(obs_props[key * "_em"], obs_props[key * "_ep"])
	else
		return obs_props[key * "_err"]
	end
end

# ╔═╡ c2e842eb-06dd-48da-8948-366983f240a9
quantile(LogitNormal(-6.0, 5.0), [0.001, 0.01, 0.5, 0.84, 0.99, 0.999])

# ╔═╡ d254b67b-8017-4bbf-ba04-50b5b76f48f0
prof_counts = LilGuys.StellarProfile(best_stars.R_ell, bins=log10.(bins),  errors=:weighted)

# ╔═╡ 939aa449-26ea-4d2d-8773-42a73c3d0410
if all_profiles
	skip = 1
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	Ns = size(best_stars, 1)
	
	global fsats = Matrix{Float64}(undef, Ns, Nc)
	global psats = Matrix{Float64}(undef, Ns, Nc)
	global profiles = Vector{LilGuys.StellarProfile}(undef, Nc)

	for i in 1:skip:Nc
		if i % 200 == 0
			@info "calculated $i / $Nc"
		end
		
		params = [df_chains[i, "params[$j]"] for j in 1:Nbins]
		radii = best_stars.R_ell

		f = hist_fractions(bins, params, radii)
		fsats[:, i] .= f

		Ls = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT
		Lb = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD
		psat =  @. f*Ls / (f*Ls + (1-f) * Lb)
		psats[:, i] .= psat

		prof = LilGuys.StellarProfile(radii, bins=log10.(bins), weights=psat, errors=:weighted)

		profiles[i] = prof
	end

end

# ╔═╡ f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
let
	skip = 1
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	Ns = size(best_stars, 1)
	
	global log_Sigmas_fast = Matrix{Float64}(undef, Nbins, Nc)
	areas = diff(π * bins .^2)
	counts = prof_counts.counts

	for i in 1:skip:Nc
		if i % 200 == 0
			@info "calculated $i / $Nc"
		end
		
		params = [df_chains[i, "params[$j]"] for j in 1:Nbins]
		log_Sigmas_fast[:, i] .= log10.(counts .* params ./ areas)
	end

end

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range on either end to allow shifts in parameter to matter less.
"""

# ╔═╡ 500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
logit(x) = 1/(1 + exp(-x))

# ╔═╡ 35fa912d-953f-4d71-b4de-32cc30e307ba
psat_err = dropdims(std(psats, dims=2), dims=2)

# ╔═╡ da737b3e-f28a-4917-91d6-eea97458ddb0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "PSAT J+24", ylabel = "PSAT MCMC")

	skip = 100

	for i in 1:skip:size(psats, 2)
		scatter!(best_stars.PSAT, psats[:, i], color=:black, alpha=0.1, markersize=2, rasterize=true)
	end
	@savefig "j+24_vs_mcmc"
	
	fig
end

# ╔═╡ f2c43c30-b069-4d17-8d61-be801d85b245
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel =log_Sigma_label,
		limits = (nothing, nothing, -8, nothing)
	)
	skip = 10
	
	for h in profiles[1:skip:end]
		scatter!(h.log_R[2:end-1] .+ 0.01*randn(length(h.log_R) - 2), h.log_Sigma[2:end-1], color=:black, alpha=0.1, markersize=1)
	end

	@savefig "scatter_profiles"

	fig
end

# ╔═╡ 08424d6a-6a21-45d2-9b1a-d28e129c8158
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = L"Gamma",
		limits = (nothing, nothing, -2, 2)
	)
	skip = 10
	
	for h in profiles[1:skip:end]
		scatter!(h.log_R[2:end-1] .+ 0.01*randn(length(h.log_R) - 2), h.Gamma[2:end-1], color=:black, alpha=0.1, markersize=1)
	end


	fig
end

# ╔═╡ e52cac33-9611-437b-98b3-bdf20fc771e0
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = "Gamma max",
		limits = (nothing, nothing, 0, nothing)
	)
	skip = 10
	
	for h in profiles[1:skip:end]
		scatter!(h.log_R[2:end-1] .+ 0.01*randn(length(h.log_R) - 2), h.Gamma_max[2:end-1], color=:black, alpha=0.1, markersize=1)
	end


	fig
end

# ╔═╡ af3db174-950d-46be-8d9c-61df70beb40e
truncated(Normal(), lower=0)

# ╔═╡ 36bb1876-5937-44fe-bb1c-add10470371d
pvalue = 0.16

# ╔═╡ 6470edfd-bd06-43b7-aa5e-a80645af78fb
psat_r_mean = dropdims(mean(psats, dims=2), dims=2)

# ╔═╡ 7cb6c373-99c6-4d89-8401-14f2fcede4d5
psat_em = psat_r_mean .- quantile.(eachrow(psats), pvalue)

# ╔═╡ 6f0fe4b8-37a8-4c62-9600-8ddfd63c8830
psat_ep = -psat_r_mean .+ quantile.(eachrow(psats), 1-pvalue)

# ╔═╡ ad85b789-8f43-489b-a41c-a966e08d78ad
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 10

	scatter!(log10.(best_stars.r_ell), psat_r_mean.- best_stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ c48f8ad0-a8c8-46ef-9b56-6f4ff5ac7a83
all_Sigmas = hcat([prof.log_Sigma for prof in profiles]...)

# ╔═╡ 084a491c-e53a-4371-9c5b-f9450ec5deec
sum(isnan.(all_Sigmas))

# ╔═╡ 428d53fd-39cd-40f3-ac6b-389f637e351e
all_Sigmas[isnan.(all_Sigmas)] .= -Inf

# ╔═╡ 5a1c8268-7ca9-45d8-8c0e-ec6bdf41dece
begin
	log_Sigmas = dropdims(median(log_Sigmas_fast, dims=2), dims=2)
	log_Sigma_em = log_Sigmas .- quantile.(eachrow(log_Sigmas_fast),pvalue)
	log_Sigma_ep =  quantile.(eachrow(log_Sigmas_fast), 1-pvalue) .- log_Sigmas
end

# ╔═╡ 9451cfcd-582e-40d6-af78-939f4b9ac161
begin
	log_Sigmas_slow = dropdims(median(all_Sigmas, dims=2), dims=2)
	log_Sigma_em_int_slow = log_Sigmas .- quantile.(eachrow(all_Sigmas),pvalue)
	log_Sigma_ep_int_slow =  quantile.(eachrow(all_Sigmas), 1-pvalue) .- log_Sigmas

	log_Sigma_em_stat_slow = dropdims(median(hcat([prof.log_Sigma_em for prof in profiles]...), dims=2), dims=2)
	log_Sigma_ep_stat_slow = dropdims(median(hcat([prof.log_Sigma_ep for prof in profiles]...), dims=2), dims=2)

	log_Sigma_em_slow = log_Sigma_em_int_slow .+ log_Sigma_em_stat_slow
	log_Sigma_ep_slow = log_Sigma_ep_int_slow .+ log_Sigma_ep_stat_slow
end

# ╔═╡ 81e5a128-c1ed-4a8a-8bb6-e4efd0593c80
sum(best_stars.PSAT)

# ╔═╡ 80bcca59-4721-45b6-bd99-cdf2be43c8d7
bins

# ╔═╡ d98b5b58-42a9-4119-bec2-b8435b35eb03


# ╔═╡ 8abe3d92-f9e8-4e89-96ec-14eaf9482e64
Nsteps = size(df_chains, 1)

# ╔═╡ 94d6a13f-1165-4d8f-a40a-a4cfa9224093
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel =log_Sigma_label,
		limits = (nothing, nothing, -8, nothing)
	)
	skip = 10
	
	for i in 1:skip:Nsteps
		y = log_Sigmas_fast[2:end-1, i]
		x = prof_counts.log_R[2:end-1] .+ 0.01*randn(length(prof_counts.log_R) - 2)
		
		filt = isfinite.(y)
		scatter!(x[filt], y[filt], color=:black, alpha=0.1, markersize=1)
	end

	@savefig "scatter_profiles"

	fig
end

# ╔═╡ 75930ff5-1add-4b1c-a495-64c50835ec4b
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		yscale=log10,
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = L"$\log\, r_\textrm{ell}$ / arcmin",
		ylabel = L"f_\textrm{sat}(r_\textrm{ell})",
		#limits=(log10(bins_robust[2]), log10(bins_robust[end-1]), nothing, nothing,)
	)

	xs = []
	ys = []
	yerrs = []
	
	for i in 1:length(bins)-2 # avoid last bin
		x = midpoints(bins)[i]
		yy = df_chains[:, "params[$i]"]
		y = median(yy)
		push!(xs, x)
		push!(ys, y)
		ye = quantile(yy, [0.16, 0.84])
		ye = (y-ye[1], ye[2]-y)
		push!(yerrs, ye)
	end

	
	errorscatter!(log10.(xs), (ys), yerror=yerrs)

	@savefig "f_sat_vs_r_ell"
	fig
end

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ 8c565a28-84bc-4bc7-8a0e-b2e0dff76665
if !isdir(joinpath(galaxy, "processed"))
	mkdir(joinpath(galaxy, "processed"))
end

# ╔═╡ d1009491-62de-42d7-89ad-b41b8385aaaf
samplesout = joinpath(galaxy, "processed", "samples$FIGSUFFIX.csv")

# ╔═╡ a93579f2-d37d-492c-95e6-02bd7491dc8c
profout = joinpath(galaxy, "processed", "profile$FIGSUFFIX.toml")

# ╔═╡ 639c7d31-8693-45dd-88be-492b124804e9
infoout = joinpath(galaxy, "processed", "info$FIGSUFFIX.log")

# ╔═╡ 28b67692-45a9-435d-8862-882b0f06f947
starsout = joinpath(galaxy, "processed", "stars$FIGSUFFIX.fits")

# ╔═╡ bf8d5936-1a4e-47c2-bb22-531ab344b8ad
filt_max = (bins[2:end] .< r_max)

# ╔═╡ 4482bead-75a5-413f-b442-e60ab3d3ec88
prof_mc = LilGuys.StellarProfile(
	R_units = "arcmin",
	log_Sigma=log_Sigmas[filt_max],
	log_Sigma_em=log_Sigma_em[filt_max],
	log_Sigma_ep=log_Sigma_ep[filt_max],
	log_R_bins=log10.(bins[bins.< r_max]), # use normal bins for these
	log_R=midpoints(log10.(bins))[filt_max],
)

# ╔═╡ 7c02d24f-7d59-410b-99bc-93513afa4c5d
let
	fig = Figure(size=(5*72, 3*72))
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		limits = (nothing, nothing, -12, nothing)
	)
	jitter=0.01

	filt = isfinite.(prof_mc.log_Sigma)
	errorscatter!(prof_mc.log_R[filt], prof_mc.log_Sigma[filt] ,
		yerror = collect(zip(prof_mc.log_Sigma_em, prof_mc.log_Sigma_ep))[filt],
		label="mc hist")
	
	LilGuys.plot_density_prof!(ax, prof_weighted, label="weighted")
	LilGuys.plot_density_prof!(ax, prof_simple, label="cut")
	LilGuys.plot_density_prof!(ax, prof_auto, label="auto")

	Legend(fig[1,2], ax)


	ax_res = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = "residual",
		limits=(nothing, nothing, -1, 1)
	)


	filt = isfinite.(prof_mc.log_Sigma)
	errorscatter!(prof_mc.log_R[filt], prof_mc.log_Sigma[filt] .- prof_simple.log_Sigma[filt], 
		yerror = collect(zip(prof_mc.log_Sigma_em, prof_mc.log_Sigma_ep))[filt],
		label="mc hist")
	
	errorscatter!(prof_weighted.log_R .+ jitter, prof_weighted.log_Sigma .- prof_simple.log_Sigma, 
		yerror = collect(zip(prof_weighted.log_Sigma_em, prof_weighted.log_Sigma_ep)),
		label="weighted")
	
	errorscatter!(prof_simple.log_R .- jitter, prof_simple.log_Sigma .- prof_simple.log_Sigma, 
		yerror = collect(zip(prof_simple.log_Sigma_em, prof_simple.log_Sigma_ep)),
		label="simple.")

	rowsize!(fig.layout, 2, Relative(0.3))
	rowgap!(fig.layout, 0)
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
	linkxaxes!(ax_res, ax)

	@savefig "density_versus_simple"
end

# ╔═╡ a562c141-8626-4242-9f6a-c381a4da619b
function median_of(profiles, key)
	return median.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	))[filt_max]
end

# ╔═╡ ed16c378-2786-4f83-8152-144838dedd19
function err_of(profiles, key)
	m = median_of(profiles, key)
	A = hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)

	@info "$key = $(sum(isnan.(A))) nans"
	if key == "Sigma"
		A[isnan.(A)] .= 0
	end
	if key == "Gamma_max"
		A[isnan.(A)] .= Inf
	end

	h = quantile.(eachrow(A), 1-pvalue)

	l = quantile.(eachrow(A), pvalue)

	e = median.(eachrow(
		hcat([getproperty(prof, Symbol(key * "_err")) for prof in profiles]...)
	))

	return (@. max(h[filt_max]-m, m-l[filt_max]) + e[filt_max])
end

# ╔═╡ 6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
open(profout, "w") do f
	print(f, prof_mc)
end

# ╔═╡ 8b05449b-6d03-40cc-b4d8-9230303c15c1
let
	stars_new = copy(best_stars)
	stars_new[!, :PSAT] = psat_r_mean
	stars_new[!, :PSAT_em] = psat_em
	stars_new[!, :PSAT_ep] = psat_ep
	stars_new[!, :f_sat] = median.(eachrow(fsats))

	LilGuys.write_fits(starsout, stars_new)
end	

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═40755684-2881-44f8-9fe6-09c211bbfe45
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═109b91c2-c0ad-4683-849b-a56da5766ee4
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═28550d54-ad88-4df5-ac04-f46a15588ff8
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═8bbd8c9f-2d35-47d9-8c30-0c86c14d6ce4
# ╠═57a19d65-59d9-46cf-8916-d9ac3a4dc92b
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═2ab34a36-542c-41da-a83d-6eb67cae253c
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╠═04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╠═f078d929-62be-48f2-acdb-808357801a7b
# ╠═4b77612f-e124-4d77-974c-d40c2f5a37ff
# ╠═91867748-9b36-4f62-9310-8b778935776b
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╠═5feaa9b0-ca44-48d6-80d1-16303a34005f
# ╠═9bee7d00-d281-4d2d-ab9a-cbd3333f0e8e
# ╠═ab124afe-c94a-44d9-8869-bd4d4cee3fbd
# ╠═babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
# ╠═e60ea082-818d-46f8-bf50-30a932f97ef2
# ╠═82f20e46-90c7-4925-83c9-3f49a909664d
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═e7bf7a13-afca-49b8-8771-f7914adb347b
# ╠═27c5a02f-0be6-4d0b-9c7f-99be12912732
# ╠═e74e6a96-a18e-40bf-ade8-07ce0e30e5c4
# ╠═6f016a8e-38ae-4f05-a7ee-c292ac0e5741
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═64c31de5-41fa-47e8-9209-61901e05f1be
# ╠═9ea0c184-4a09-4293-a6d7-a963cf1d3d58
# ╠═984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═e87b5ee8-e281-43e3-8172-4848018e2cd2
# ╠═c2e842eb-06dd-48da-8948-366983f240a9
# ╠═d254b67b-8017-4bbf-ba04-50b5b76f48f0
# ╠═939aa449-26ea-4d2d-8773-42a73c3d0410
# ╠═f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
# ╠═1d5de2ad-db30-4d18-845c-687bdf7ff2d1
# ╠═500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
# ╠═35fa912d-953f-4d71-b4de-32cc30e307ba
# ╠═7cb6c373-99c6-4d89-8401-14f2fcede4d5
# ╠═6f0fe4b8-37a8-4c62-9600-8ddfd63c8830
# ╠═da737b3e-f28a-4917-91d6-eea97458ddb0
# ╠═ad85b789-8f43-489b-a41c-a966e08d78ad
# ╠═f2c43c30-b069-4d17-8d61-be801d85b245
# ╠═94d6a13f-1165-4d8f-a40a-a4cfa9224093
# ╠═08424d6a-6a21-45d2-9b1a-d28e129c8158
# ╠═e52cac33-9611-437b-98b3-bdf20fc771e0
# ╠═af3db174-950d-46be-8d9c-61df70beb40e
# ╠═36bb1876-5937-44fe-bb1c-add10470371d
# ╠═6470edfd-bd06-43b7-aa5e-a80645af78fb
# ╠═c48f8ad0-a8c8-46ef-9b56-6f4ff5ac7a83
# ╠═084a491c-e53a-4371-9c5b-f9450ec5deec
# ╠═428d53fd-39cd-40f3-ac6b-389f637e351e
# ╠═9451cfcd-582e-40d6-af78-939f4b9ac161
# ╠═5a1c8268-7ca9-45d8-8c0e-ec6bdf41dece
# ╠═81e5a128-c1ed-4a8a-8bb6-e4efd0593c80
# ╠═80bcca59-4721-45b6-bd99-cdf2be43c8d7
# ╠═7c02d24f-7d59-410b-99bc-93513afa4c5d
# ╠═d98b5b58-42a9-4119-bec2-b8435b35eb03
# ╠═8abe3d92-f9e8-4e89-96ec-14eaf9482e64
# ╠═75930ff5-1add-4b1c-a495-64c50835ec4b
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═639c7d31-8693-45dd-88be-492b124804e9
# ╠═28b67692-45a9-435d-8862-882b0f06f947
# ╠═4482bead-75a5-413f-b442-e60ab3d3ec88
# ╠═bf8d5936-1a4e-47c2-bb22-531ab344b8ad
# ╠═a562c141-8626-4242-9f6a-c381a4da619b
# ╠═ed16c378-2786-4f83-8152-144838dedd19
# ╠═6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
# ╠═8b05449b-6d03-40cc-b4d8-9230303c15c1
