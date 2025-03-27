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
	using Turing
end

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 2f62c5c2-e397-463b-9e73-f554c31a7b85
using PairPlots

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
include("../utils/gaia_filters.jl")

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	galaxy = "crater2"
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 109b91c2-c0ad-4683-849b-a56da5766ee4
import StatsBase: sem

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(galaxy, "figures"); FIGSUFFIX=".mcmc_hist_nostruct"

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

# ╔═╡ c58ed939-07c1-4d34-90a4-3c5d9cc9910b


# ╔═╡ b9cf2c23-6c9a-44f5-9740-22caf4959831
obs_props = TOML.parsefile(joinpath(galaxy, "observed_properties.toml"))

# ╔═╡ f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
all_stars = read_gaia_stars(GaiaFilterParams(filename = datafile,
	ra=obs_props["ra"],
	dec = obs_props["dec"],
	ellipticity = obs_props["ellipticity"],
	PA = obs_props["position_angle"]
))

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

# ╔═╡ 1d9b3718-45d7-4765-ac7d-017dbf939ec8
md"""
# Histogram model
"""

# ╔═╡ e7bf7a13-afca-49b8-8771-f7914adb347b
function hist_fractions(bins, params, x)
	i = DE.bin_indices(x, bins)
	return params[i]
end

# ╔═╡ 27c5a02f-0be6-4d0b-9c7f-99be12912732
pos_err = (sem(filter(r->r.PSAT .> 0.2, best_stars).xi) + sem(filter(r->r.PSAT .> 0.2, best_stars).eta))/2

# ╔═╡ ccd115f0-5a6e-4a65-a45f-d97cdbbbe02c
bin_width_min = 0.05

# ╔═╡ 8d3bbb76-47e1-4a36-981f-487b657bc74a
N_per_bin_min = max(round(Int, LilGuys.Interface.default_n_per_bin(best_stars.r_ell[best_stars.PSAT .> 0.2], nothing)), 2)

# ╔═╡ 1a5e8de2-0d74-48a9-aead-0855602734f3
bins = 10 .^ LilGuys.Interface.bins_both(log10.(best_stars.r_ell), nothing, bin_width=bin_width_min, num_per_bin=N_per_bin_min)


# ╔═╡ 8e665feb-1a41-440c-a813-163cbf3be4f8
Nmemb = sum(best_stars.PSAT)

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
	
	stephist!(log10.(best_stars.r_ell[best_stars.PSAT .> 0.2]), bins=log10.(bins), label="PSAT > 0.2")
	stephist!(log10.(best_stars.r_ell), bins=log10.(bins), label="best")

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
prof_simple = LilGuys.StellarProfile(best_stars.r_ell[best_stars.PSAT .> 0.2], bins=log10.(bins)[bins .< r_max], errors=:weighted)

# ╔═╡ 82f20e46-90c7-4925-83c9-3f49a909664d
sum(prof_simple.counts)

# ╔═╡ ab124afe-c94a-44d9-8869-bd4d4cee3fbd
prof_weighted = LilGuys.StellarProfile(best_stars.r_ell, weights=best_stars.PSAT, bins=log10.(bins)[bins .< r_max], errors=:weighted)

# ╔═╡ babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
prof_bernoulli = LilGuys.StellarProfile(best_stars.r_ell, weights=best_stars.PSAT, bins=log10.(bins)[bins .< r_max], errors=:bernoulli)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function hist_model(r_b::Vector{Int}, Lsat::Vector{Float64}, Lbkd::Vector{Float64}, Nb)
	params ~ filldist(LogitNormal(-6.0, 5.0), Nb)

	f = params[r_b]

	LL = sum(@. log10.(f*Lsat + (1-f) * Lbkd) )
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 562fe2c7-3ed4-4adb-af7f-1eca2aa35d8b
Nbins = length(bins) - 1	

# ╔═╡ ba5ceef7-ac10-46f7-a89b-30f38b6ddec1
r_b = DE.bin_indices(best_stars.r_ell, bins)

# ╔═╡ df84d941-7ddf-4c1c-96bb-1858b49bb710
Lsat = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT

# ╔═╡ 97cbf4b6-025d-461b-b82e-044f86b713c1
Lbg = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD

# ╔═╡ cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
model = hist_model(r_b, Lsat, Lbg, Nbins)

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range on either end to allow shifts in parameter to matter less.
"""

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
chains = sample(model, NUTS(0.25), MCMCThreads(), 1000, 4) 

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
summarize(chains)

# ╔═╡ 500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
logit(x) = 1/(1 + exp(-x))

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fig = Figure(size=(3*72, 10*72))

	Nvar = length(bins)-1
	for i in 1:Nvar
		ax = Axis(fig[i, 1])
		for c in 1:size(chains, 3)
			y = chains[:, i, c]
			lines!(logit.(y), linewidth=0.1)
		end

		if i < Nvar
			hidexdecorations!(ax)
		end
	end
	rowgap!(fig.layout, 0)

	fig

end

# ╔═╡ af3db174-950d-46be-8d9c-61df70beb40e
truncated(Normal(), lower=0)

# ╔═╡ 36bb1876-5937-44fe-bb1c-add10470371d
pvalue = 0.16

# ╔═╡ d98b5b58-42a9-4119-bec2-b8435b35eb03


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
	
	for i in 1:length(bins)-1
		x = midpoints(bins)[i]
		y = median(chains[:, i, :])
		push!(xs, x)
		push!(ys, y)
		ye = quantile(chains[:, i, :], [0.16, 0.84])
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
filt_max = bins[2:end] .< r_max

# ╔═╡ a562c141-8626-4242-9f6a-c381a4da619b
function median_of(profiles, key)
	return median.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	))[filt_max]
end

# ╔═╡ ed16c378-2786-4f83-8152-144838dedd19
function err_of(profiles, key)
	m = median_of(profiles, key)
	h = quantile.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	), 1-pvalue)

	l = quantile.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	), pvalue)

	e = median.(eachrow(
		hcat([getproperty(prof, Symbol(key * "_err")) for prof in profiles]...)
	))

	return (@. max(h[filt_max]-m, m-l[filt_max]) + e[filt_max])
end

# ╔═╡ 36db7e1d-4b48-4510-99d5-7d567ac70d5d
df_out = DataFrame(chains)

# ╔═╡ ae9b79f2-c8ab-476c-a1f4-9323d64b20a4
df_out[:, "params[31]"]

# ╔═╡ 349a1568-5d9d-4f44-a334-cb89857d5c4b
let
	global psat_r = []
	global profiles = LilGuys.StellarProfile[]
	global fsats = []

	skip = 1
	
	for i in 1:skip:size(df_out, 1)
		params = [df_out[i, "params[$j]"] for j in 1:Nbins]
		
		f = hist_fractions(bins, params, best_stars.r_ell)
		push!(fsats, f)

		Ls = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT
		Lb = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD
		psat = @. f*Ls / (f*Ls + (1-f) * Lb)
		push!(psat_r, psat)

		prof = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat, errors=:weighted)

		push!(profiles, prof)
	end

	psat_r = hcat(psat_r...)
end

# ╔═╡ 35fa912d-953f-4d71-b4de-32cc30e307ba
psat_err = dropdims(std(psat_r, dims=2), dims=2)

# ╔═╡ da737b3e-f28a-4917-91d6-eea97458ddb0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "PSAT J+24", ylabel = "PSAT MCMC")

	skip = 100

	for i in 1:skip:size(psat_r, 2)
		scatter!(best_stars.PSAT, psat_r[:, i], color=:black, alpha=0.1, markersize=2, rasterize=true)
	end
	@savefig "j+24_vs_mcmc"
	
	fig
end

# ╔═╡ 08424d6a-6a21-45d2-9b1a-d28e129c8158
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		limits = (nothing, nothing, -8, nothing)
	)
	skip = 10
	
	for h in profiles[1:skip:end]
		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.log_Sigma, color=:black, alpha=0.1, markersize=1)
	end

	@savefig "scatter_profiles"

	fig
end

# ╔═╡ e52cac33-9611-437b-98b3-bdf20fc771e0
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = L"\Gamma",
		limits = (nothing, nothing, -10, 10)
	)
	skip = 10
	
	for h in profiles[1:skip:end]
		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.Gamma, color=:black, alpha=0.1, markersize=1)
	end
	@savefig "scatter_profiles"


	fig
end

# ╔═╡ bc7ab4bd-5dff-4b08-b2d8-017836c45819
psat_r

# ╔═╡ 6470edfd-bd06-43b7-aa5e-a80645af78fb
psat_r_mean = dropdims(mean(psat_r, dims=2), dims=2)

# ╔═╡ 7cb6c373-99c6-4d89-8401-14f2fcede4d5
psat_em = psat_r_mean .- quantile.(eachrow(psat_r), pvalue)

# ╔═╡ 6f0fe4b8-37a8-4c62-9600-8ddfd63c8830
psat_ep = -psat_r_mean .+ quantile.(eachrow(psat_r), 1-pvalue)

# ╔═╡ ad85b789-8f43-489b-a41c-a966e08d78ad
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 10

	scatter!(log10.(best_stars.r_ell), psat_r_mean.- best_stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ 972e743a-cb4c-40b9-9a00-7a23cd7520e7
scatter(psat_r_mean, psat_r_mean .- best_stars.PSAT)

# ╔═╡ 9451cfcd-582e-40d6-af78-939f4b9ac161
begin
	all_Sigmas = hcat([prof.log_Sigma for prof in profiles]...)
	log_Sigmas = dropdims(median(all_Sigmas, dims=2), dims=2)
	log_Sigma_em_int = log_Sigmas .- quantile.(eachrow(all_Sigmas),pvalue)
	log_Sigma_ep_int =  quantile.(eachrow(all_Sigmas), 1-pvalue) .- log_Sigmas

	log_Sigma_em_stat = dropdims(median(hcat([prof.log_Sigma_em for prof in profiles]...), dims=2), dims=2)
	log_Sigma_ep_stat = dropdims(median(hcat([prof.log_Sigma_ep for prof in profiles]...), dims=2), dims=2)

	log_Sigma_em = log_Sigma_em_int .+ log_Sigma_em_stat
	log_Sigma_ep = log_Sigma_ep_int .+ log_Sigma_ep_stat
end

# ╔═╡ 4482bead-75a5-413f-b442-e60ab3d3ec88
prof_mc = LilGuys.StellarProfile(
	r_units = "arcmin",
	log_Sigma=log_Sigmas[filt_max],
	log_Sigma_em=log_Sigma_em[filt_max],
	log_Sigma_ep=log_Sigma_ep[filt_max],
	log_r_bins=log10.(bins[bins .< r_max]), # use normal bins for these
	log_r=midpoints(log10.(bins))[filt_max],
	counts=median_of(profiles, "counts"),
	mass_in_annulus=median_of(profiles, "mass_in_annulus"),
	mass_in_annulus_err=err_of(profiles, "mass_in_annulus"),
	M_in=median_of(profiles, "M_in"),
	M_in_err=err_of(profiles, "M_in"),
	Sigma=median_of(profiles, "Sigma"),
	Sigma_err=err_of(profiles, "Sigma"),
	Sigma_m=median_of(profiles, "Sigma_m"),
	Sigma_m_err=err_of(profiles, "Sigma_m"),
	Gamma=median_of(profiles, "Gamma"),
	Gamma_err=err_of(profiles, "Gamma"),
	Gamma_max=median_of(profiles, "Gamma_max"),
	Gamma_max_err=err_of(profiles, "Gamma_max"),
)

# ╔═╡ 7c02d24f-7d59-410b-99bc-93513afa4c5d
let
	fig = Figure(size=(5*72, 3*72))
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
	)
	jitter=0.01

	LilGuys.plot_density_prof!(ax, prof_mc, label="mc hist")
	LilGuys.plot_density_prof!(ax, prof_weighted, label="weighted")
	LilGuys.plot_density_prof!(ax, prof_simple, label="cut")
	
	Legend(fig[1,2], ax)


	ax_res = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = "residual",
		limits=(nothing, nothing, -1, 1)
	)


	errorscatter!(prof_mc.log_r, prof_mc.log_Sigma .- prof_simple.log_Sigma, 
		yerror = collect(zip(prof_mc.log_Sigma_em, prof_mc.log_Sigma_ep)),
		label="mc hist")
	
	errorscatter!(prof_weighted.log_r .+ jitter, prof_weighted.log_Sigma .- prof_simple.log_Sigma, 
		yerror = collect(zip(prof_weighted.log_Sigma_em, prof_weighted.log_Sigma_ep)),
		label="weighted")
	
	errorscatter!(prof_simple.log_r .- jitter, prof_simple.log_Sigma .- prof_simple.log_Sigma, 
		yerror = collect(zip(prof_simple.log_Sigma_em, prof_simple.log_Sigma_ep)),
		label="simple.")

	rowsize!(fig.layout, 2, Relative(0.3))
	rowgap!(fig.layout, 0)
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)

	@savefig "density_versus_simple"
end

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 52d66c81-64cc-4eee-bad9-b9728c38c97b
psat_err

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
	stars_new[!, :f_sat] = median.(eachrow(hcat(fsats...)))

	LilGuys.write_fits(starsout, stars_new)
end	

# ╔═╡ f9e8d5d8-dd6c-4d7c-bedb-02a48d241ca8
median.(eachrow(hcat(fsats...)))

# ╔═╡ 77f564ee-f461-4107-b042-ae6d14ba9867
open(infoout, "w") do f
	println(f)
	println(f, "N per bin: $N_per_bin_min")
	println(f, "min binwidth: $bin_width_min")
	println(f)
	
	println(f, "bins")
	println(f, bins)


	println(f)
	println(f, "chains")
	println(f, DataFrame(summarize(chains)))
end

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═109b91c2-c0ad-4683-849b-a56da5766ee4
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═2f62c5c2-e397-463b-9e73-f554c31a7b85
# ╠═28550d54-ad88-4df5-ac04-f46a15588ff8
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═57a19d65-59d9-46cf-8916-d9ac3a4dc92b
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═2ab34a36-542c-41da-a83d-6eb67cae253c
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╠═04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╠═91867748-9b36-4f62-9310-8b778935776b
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╠═5feaa9b0-ca44-48d6-80d1-16303a34005f
# ╠═ab124afe-c94a-44d9-8869-bd4d4cee3fbd
# ╠═babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
# ╠═82f20e46-90c7-4925-83c9-3f49a909664d
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═e7bf7a13-afca-49b8-8771-f7914adb347b
# ╠═27c5a02f-0be6-4d0b-9c7f-99be12912732
# ╠═1a5e8de2-0d74-48a9-aead-0855602734f3
# ╠═ccd115f0-5a6e-4a65-a45f-d97cdbbbe02c
# ╠═8d3bbb76-47e1-4a36-981f-487b657bc74a
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═64c31de5-41fa-47e8-9209-61901e05f1be
# ╠═9ea0c184-4a09-4293-a6d7-a963cf1d3d58
# ╠═984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═562fe2c7-3ed4-4adb-af7f-1eca2aa35d8b
# ╠═cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
# ╠═ba5ceef7-ac10-46f7-a89b-30f38b6ddec1
# ╠═df84d941-7ddf-4c1c-96bb-1858b49bb710
# ╠═97cbf4b6-025d-461b-b82e-044f86b713c1
# ╠═1d5de2ad-db30-4d18-845c-687bdf7ff2d1
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
# ╠═ae9b79f2-c8ab-476c-a1f4-9323d64b20a4
# ╠═349a1568-5d9d-4f44-a334-cb89857d5c4b
# ╠═35fa912d-953f-4d71-b4de-32cc30e307ba
# ╠═7cb6c373-99c6-4d89-8401-14f2fcede4d5
# ╠═6f0fe4b8-37a8-4c62-9600-8ddfd63c8830
# ╠═da737b3e-f28a-4917-91d6-eea97458ddb0
# ╠═ad85b789-8f43-489b-a41c-a966e08d78ad
# ╠═08424d6a-6a21-45d2-9b1a-d28e129c8158
# ╠═e52cac33-9611-437b-98b3-bdf20fc771e0
# ╠═af3db174-950d-46be-8d9c-61df70beb40e
# ╠═bc7ab4bd-5dff-4b08-b2d8-017836c45819
# ╠═36bb1876-5937-44fe-bb1c-add10470371d
# ╠═6470edfd-bd06-43b7-aa5e-a80645af78fb
# ╠═972e743a-cb4c-40b9-9a00-7a23cd7520e7
# ╠═9451cfcd-582e-40d6-af78-939f4b9ac161
# ╠═7c02d24f-7d59-410b-99bc-93513afa4c5d
# ╠═d98b5b58-42a9-4119-bec2-b8435b35eb03
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
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═52d66c81-64cc-4eee-bad9-b9728c38c97b
# ╠═6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
# ╠═8b05449b-6d03-40cc-b4d8-9230303c15c1
# ╠═f9e8d5d8-dd6c-4d7c-bedb-02a48d241ca8
# ╠═77f564ee-f461-4107-b042-ae6d14ba9867
