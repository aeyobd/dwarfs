### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
	using Turing
end

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 2f62c5c2-e397-463b-9e73-f554c31a7b85
using PairPlots

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
include("../utils/gaia_filters.jl")

# ╔═╡ cb7f2108-173a-4882-a95a-9111d1330176
md"""
This notebook compares different MCMC methods to check how much they might differ.
"""

# ╔═╡ 109b91c2-c0ad-4683-849b-a56da5766ee4
import StatsBase: sem

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
@bind galaxy TextField(default="draco")

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(galaxy, "figures")

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

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

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


# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

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

# ╔═╡ dd1b6b08-8821-4c3d-b937-c0d0705d2c33
md"""
# Histogram model
"""

# ╔═╡ 05333927-2605-4e8c-afd4-a2bfb120b8ce
function hist_fractions(bins, params, x)
	i = DE.bin_indices(x, bins)
	return params[i]
end

# ╔═╡ 69b19f98-7a40-4874-821e-5dcf58fa70a6
@model function hist_model(data, bins;)
	Nb = length(bins) - 1
	
	θ ~ filldist(LogitNormal(-0.0, 4.0), Nb)
	
	radii = data.r_ell
	r_b = DE.bin_indices(radii, bins)
	f = θ[r_b]

	LL = sum(@. log(
		f*(data.L_PM_SAT * data.L_CMD_SAT) 
		+ (1-f) * data.L_PM_BKD * data.L_CMD_BKD
	))
	Turing.@addlogprob!(LL)
end

# ╔═╡ b14e1b89-6cb2-4c26-b0fe-1e2335b07d05
pos_err = (sem(filter(r->r.PSAT .> 0.2, best_stars).xi) + sem(filter(r->r.PSAT .> 0.2, best_stars).eta))/2

# ╔═╡ c5c505d3-ed61-41bc-aaa5-b6c7037273ab
bin_width_min = 0.05

# ╔═╡ b393ac76-4961-41ed-b847-900ee81364b1
Nmemb = sum(best_stars.PSAT)

# ╔═╡ 56b1b79d-7990-4a88-8ce0-273b211cbf75
N_per_bin_min = max(round(Int, LilGuys.Interface.default_n_per_bin(best_stars.r_ell[best_stars.PSAT .> 0.2], nothing)), 2)

# ╔═╡ 66dd7757-7f10-4202-ad5f-a193b57cdf41
bins = 10 .^ LilGuys.Interface.bins_both(log10.(best_stars.r_ell), nothing, bin_width=bin_width_min, num_per_bin=N_per_bin_min)

# ╔═╡ cc7b4837-1974-4098-86d6-95b1dace8a72
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

# ╔═╡ c17a4353-564c-45d2-a561-ac71fc340baf
10 .^ LilGuys.Interface.bins_both(log10.(best_stars.r_ell), nothing, bin_width=0.1, num_per_bin=10)

# ╔═╡ df7cb725-380e-409f-bad6-1f1fca77a5b5
model = hist_model(best_stars, bins)

# ╔═╡ 85f68616-28da-4d2d-83cc-3d59e6c433f1
chain = sample(model, NUTS(0.65), 100)

# ╔═╡ 6136af16-85ed-4f30-b995-eb6c1f4b6af6
summarize(chain)

# ╔═╡ 0793f4fd-0e8c-4afb-a4cc-e9ded55bebe3
@time LilGuys.calc_r_ell(best_stars.xi, best_stars.eta, obs_props["ellipticity"], obs_props["position_angle"])

# ╔═╡ 1ff82516-8b0a-4a95-bbca-f0230926f136
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		yscale=log10,
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = L"$\log\, r_\textrm{ell}$ / arcmin",
		ylabel = L"f_\textrm{sat}(r_\textrm{ell})"
	)

	xs = []
	ys = []
	yerrs = []
	
	for i in 1:length(bins)-1
		x = midpoints(bins)[i]
		y = median(chain[:, i, :])
		push!(xs, x)
		push!(ys, y)
		ye = quantile(chain[:, i, :], [0.16, 0.84])
		ye = (y-ye[1], ye[2]-y)
		push!(yerrs, ye)
	end

	
	errorscatter!(log10.(xs), (ys), yerror=yerrs)

	fig
end

# ╔═╡ 5a272c79-0f27-40d2-a197-2b6db3f937be
bins

# ╔═╡ 18e0757d-c4eb-469d-b601-119d3f03fb7c
maximum(best_stars.xi .⊕ best_stars.eta) * 60

# ╔═╡ 3040d802-1c29-4d7d-874e-c6f0b198baee
maximum(best_stars.r_ell)

# ╔═╡ 79f36951-8b07-4569-a792-aba6db9c3d41
r_max = 240 * sqrt(1 - obs_props["ellipticity"])

# ╔═╡ 5feaa9b0-ca44-48d6-80d1-16303a34005f
prof_simple = LilGuys.StellarProfile(best_stars.r_ell[best_stars.PSAT .> 0.2], bins=log10.(bins)[bins .< r_max], errors=:weighted)

# ╔═╡ 82f20e46-90c7-4925-83c9-3f49a909664d
sum(prof_simple.counts)

# ╔═╡ ab124afe-c94a-44d9-8869-bd4d4cee3fbd
prof_weighted = LilGuys.StellarProfile(best_stars.r_ell, weights=best_stars.PSAT, bins=log10.(bins)[bins .< r_max], errors=:weighted)

# ╔═╡ babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
prof_bernoulli = LilGuys.StellarProfile(best_stars.r_ell, weights=best_stars.PSAT, bins=log10.(bins)[bins .< r_max], errors=:bernoulli)

# ╔═╡ 455471ee-bf3e-441f-bd9d-b53252c5639f
let
	fig = Figure()
	ax = Axis(fig[1,1])
	skip = 1
	global psats = []
	
	for i in 1:skip:length(chain)
		params = chain.value[i, :, 1]
		f = hist_fractions(bins, params, best_stars.r_ell)

		Ls = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT
		Lb = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD
		psat = @. f*Ls / (f*Ls + (1-f) * Lb)

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat, errors=:bernoulli)

		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.log_Sigma, color=:black, alpha=0.1, markersize=1)

		push!(psats, psat)
	end

	LilGuys.plot_density_prof!(ax, prof_simple)

	psats = hcat(psats...)
	fig
end

# ╔═╡ 546ac410-5695-408c-81e0-d041fa12750d
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label
	)
	skip = 1

	ys = []
	for i in 1:skip:length(chain)
		psat = psats[:, i]

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat, errors=:weighted)

		push!(ys, h.log_Sigma)
	end

	ys = hcat(ys...)
	ys_m = dropdims(mean(ys, dims=2), dims=2)
	ys_e = dropdims(std(ys, dims=2), dims=2)

	errorscatter!(prof_simple.log_r, prof_simple.log_Sigma, yerror=prof_simple.log_Sigma_err)
	filt = bins[2:end] .< r_max
	errorscatter!(log10.(midpoints(bins))[filt], ys_m[filt], yerror=ys_e[filt])

	fig
end

# ╔═╡ 9846b747-2d15-408b-858d-b7e7f8229a04
mean(hcat(psats...), dims=2)

# ╔═╡ 8179ee4d-ece5-41ee-9814-bda858e0a3a8
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "psat jax", ylabel = "psat me")

	skip = 10

	for i in 1:skip:length(chain)
		scatter!(best_stars.PSAT, psats[:, i], color=:black, alpha=0.1, markersize=5, rasterize=true)
	end
	
	fig
end

# ╔═╡ 381bb765-9c64-4a7e-bac6-7b69de22b1db
psat_mean = dropdims(mean(psats, dims=2), dims=2)

# ╔═╡ efeb208c-afb8-4165-8897-7d356cf1a5ce
psat_err = dropdims(std(psats, dims=2), dims=2)

# ╔═╡ 3870d70e-5c09-491e-be99-cf0aadca6ed6
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 10

	scatter!(log10.(best_stars.r_ell), psat_mean.- best_stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ 1b189fd5-e9c4-47d0-b389-04d2a66a307c
md"""
## Spline model
"""

# ╔═╡ 28fe8c1b-de0b-4158-a883-bde96ff92825
pkgversion(DE)

# ╔═╡ 5d3b661b-ee45-42e5-89b1-3283e9b66c14
import BSplineKit as BSK

# ╔═╡ 9aa3ffb7-17ce-47ea-9bac-a27d656fa4a2
spline_basis = BSK.BSplineBasis(BSK.BSplineOrder(3), copy(bins))

# ╔═╡ 4ff739d1-0515-48ba-94e2-2012fd0bf865
logit(x) = 1/(1 + exp(-x))

# ╔═╡ b77d9b07-a50f-40d3-b50d-8e2935ea6476
function spline_fsat(spline_basis, coeffs, radii)
	spline =  BSK.Spline(spline_basis, coeffs)
	f = @. logit.(spline.(radii))
	return f
end

# ╔═╡ 12b74feb-1b58-4d31-bfd5-a9e0bb3768ef
@model function spline_model(data, spline_basis;)
	Nb = length(spline_basis) 
	
	θ ~ MvNormal(zeros(Nb), 2.0)
	
	radii = data.r_ell
	spline = BSK.Spline(spline_basis, θ)

	f = @. logit.(spline.(radii))

	LL = sum(@. log(
		f*(data.L_PM_SAT * data.L_CMD_SAT) 
		+ (1-f) * data.L_PM_BKD * data.L_CMD_BKD
	))
	Turing.@addlogprob!(LL)
end

# ╔═╡ 4632607a-ead9-4e3c-9ea4-e7959c7affc0
model_s = spline_model(best_stars, spline_basis)

# ╔═╡ a070ea1d-ecd6-41fc-afc8-21430e800b77
chain_s = sample(model_s, NUTS(), 100)

# ╔═╡ 6bf9eb27-7654-4196-a6f6-5b232db348c7
spline_basis.t

# ╔═╡ 7d80653c-83c6-4e28-bc88-48cceeda297a
length(spline_basis)

# ╔═╡ 1765d2dd-9e59-480d-a636-7a89147fc52b
BSK.Spline(spline_basis, collect(chain_s.value[1, 1:length(spline_basis), 1]))

# ╔═╡ e225d7ec-29fa-4b33-85b7-c453c1c5a639
chain_s[1, :, 1].value

# ╔═╡ 12a55136-bb11-4fd4-b76f-b71c360dca67
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = L"$\log\, r_\textrm{ell}$ / arcmin",
		ylabel = L"f_\textrm{sat}"
	)

	xs = 10 .^ LinRange(log10(spline_basis.t[begin]), log10(spline_basis.t[end]), 300)
	ys_a = []
	
	for i in 1:length(chain_s)
		ys = spline_fsat(spline_basis, collect(chain_s.value[i, 1:length(spline_basis), 1]), xs)
		push!(ys_a, ys)

	end
	ys_a = hcat(ys_a...)
	ys = dropdims(median(ys_a, dims=2), dims=2)
	ys_l = quantile.(eachrow(ys_a), 0.16)
	ys_h = quantile.(eachrow(ys_a), 0.84)
	ys_err = collect(zip(ys_l, ys_h))
	
	band!(log10.(xs), ys_l, ys_h, alpha=0.5)
	lines!(log10.(xs), ys)

	fig
end

# ╔═╡ 199a2765-c2f8-46dc-826b-d60875df9bb2
let
	fig = Figure()
	ax = Axis(fig[1,1])
	skip = 1
	global psats_s = []
	
	for i in 1:skip:length(chain_s)
		params = collect(chain_s.value[i, 1:length(spline_basis), 1])
		f = spline_fsat(spline_basis, params, best_stars.r_ell)

		Ls = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT
		Lb = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD
		psat = @. f*Ls / (f*Ls + (1-f) * Lb)

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat, errors=:weighted)

		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.log_Sigma, color=:black, alpha=0.1, markersize=1)

		push!(psats_s, psat)
	end

	psats_s = hcat(psats_s...)
	fig
end

# ╔═╡ c3827324-b69c-4ddf-83cc-3568ef9b5229
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label
	)
	skip = 1

	ys = []
	ys_e_stat = []
	for i in 1:skip:length(chain)
		psat = psats_s[:, i]

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat, errors=:weighted)

		push!(ys, h.log_Sigma)
		push!(ys_e_stat, h.log_Sigma_err)
	end

	ys = hcat(ys...)
	ys_m = dropdims(mean(ys, dims=2), dims=2)
	ys_e = dropdims(std(ys, dims=2), dims=2)
	ys_e_s = dropdims(mean(hcat(ys_e_stat...), dims=2), dims=2)

	errorscatter!(prof_simple.log_r, prof_simple.log_Sigma, yerror=prof_simple.log_Sigma_err, label="2component PSAT")
	filt = bins[2:end] .< r_max
	errorscatter!(log10.(midpoints(bins))[filt], ys_m[filt], yerror=(ys_e .+ ys_e_s)[filt], label="spline spatial")

	axislegend(position=:lb)
	fig
end

# ╔═╡ 5fdeab29-1424-47ae-b95a-53fc3d6daf86
psat_mean_s = dropdims(mean(psats_s, dims=2), dims=2)

# ╔═╡ f83d1311-af43-4b03-82aa-04b6e69f498b
psat_mean_err_s = dropdims(std(psats_s, dims=2), dims=2)

# ╔═╡ e4ab0464-af03-4e05-86b2-1b542d2c9334
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 10

	scatter!(log10.(best_stars.r_ell), psat_mean_s.- best_stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ eca07277-1540-4102-bef0-214e40123a98
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "psat me", ylabel = "psat me - jax")

	skip = 10

	scatter!(psat_mean_s, psat_mean_s.- best_stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ 63dbc12e-fb6d-4ef6-9645-b217298f6cae
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "psat me", ylabel = "psat me err")

	skip = 10

	scatter!(psat_mean_s, psat_mean_err_s, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function hist_model_robust(data, bins; 
		pos_err=pos_err, 
		position_angle=obs_props["position_angle"], 
		position_angle_err =obs_props["position_angle_err"], 
		ellipticity = obs_props["ellipticity"],
		ellipticity_err = obs_props["ellipticity_err"],
	)
	Nb = length(bins) - 1

	d_xi = rand(Normal(0.0, pos_err))
	d_eta = rand(Normal(0.0, pos_err))
	pos_ang = rand(Normal(position_angle, position_angle_err))
	ell = rand(truncated(Normal(ellipticity, ellipticity_err), lower=0, upper=0.99))
	
	θ ~ filldist(LogitNormal(-2.0, 2.0), Nb)

	xi = data.xi .+ d_xi
	eta = data.eta .+ d_eta

	radii = 60LilGuys.calc_r_ell(xi, eta, ell, pos_ang)
	
	r_b = DE.bin_indices(radii, bins)
	r_b = min.(r_b, Nb)
	r_b = max.(r_b, 1)

	f = θ[r_b]

	LL = sum(@. log10.(f*data.L_PM_SAT*data.L_CMD_SAT + (1-f) * data.L_PM_BKD*data.L_CMD_BKD))
	Turing.@addlogprob!(LL)
end

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range since stars can shift in r ell a bit more now.
"""

# ╔═╡ 98c06e5a-9789-4bd8-9e5b-0b39b5230186
bins_robust = [0; bins]

# ╔═╡ cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
model_robust = hist_model_robust(best_stars, bins_robust)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
chain_robust = sample(model_robust, HMC(0.05, 10), 100)

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
summarize(chain_robust)

# ╔═╡ 08424d6a-6a21-45d2-9b1a-d28e129c8158
let
	fig = Figure()
	ax = Axis(fig[1,1])
	skip = 1
	global psat_r = []
	
	for i in 1:skip:length(chain_robust)
		params = chain_robust.value[i, :, 1]
		f = hist_fractions(bins_robust, params, best_stars.r_ell)

		Ls = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT
		Lb = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD
		psat = @. f*Ls / (f*Ls + (1-f) * Lb)

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins_robust[2:end-1]), weights=psat, errors=:bernoulli)

		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.log_Sigma, color=:black, alpha=0.1, markersize=1)

		push!(psat_r, psat)
	end

	errorscatter!(prof_simple.log_r, prof_simple.log_Sigma, yerror=prof_simple.log_Sigma_err)

	psat_r = hcat(psat_r...)
	fig
end

# ╔═╡ 42844a8c-c8fc-4c6d-9826-02b88f6eff62
log10.(bins_robust)

# ╔═╡ af3db174-950d-46be-8d9c-61df70beb40e
truncated(Normal(), lower=0)

# ╔═╡ 7c02d24f-7d59-410b-99bc-93513afa4c5d
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		limits=(nothing, (-8, 0))
	)
	skip = 1

	ys = []
	for i in 1:skip:length(chain_robust)
		psat = psat_r[:, i]

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat, errors=:weighted)

		push!(ys, h.log_Sigma)
	end

	ys = hcat(ys...)
	ys_m = dropdims(mean(ys, dims=2), dims=2)
	ys_e = dropdims(std(ys, dims=2), dims=2)

	errorscatter!(prof_simple.log_r, prof_simple.log_Sigma, yerror=prof_simple.log_Sigma_err)
	filt = bins[2:end] .< r_max
	errorscatter!(log10.(midpoints(bins))[filt], ys_m[filt], yerror=ys_e[filt])

	fig
end

# ╔═╡ 75930ff5-1add-4b1c-a495-64c50835ec4b
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		yscale=log10,
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = L"$\log\, r_\textrm{ell}$ / arcmin",
		ylabel = L"f_\textrm{sat}(r_\textrm{ell})"
	)

	xs = []
	ys = []
	yerrs = []
	
	for i in 1:length(bins)-1
		x = midpoints(bins)[i]
		y = median(chain_robust[:, i, :])
		push!(xs, x)
		push!(ys, y)
		ye = quantile(chain_robust[:, i, :], [0.16, 0.84])
		ye = (y-ye[1], ye[2]-y)
		push!(yerrs, ye)
	end

	
	errorscatter!(log10.(xs), (ys), yerror=yerrs)

	fig
end

# ╔═╡ bb97c5cc-35f4-42db-b5fa-c087e4527e0e
md"""
# Adding extra fit parameters
"""

# ╔═╡ 01aa60ce-3849-4f3e-81d9-090c6aa8d104
@model function hist_model_extra(data, bins; 
		pos_err=pos_err, 
		position_angle=obs_props["position_angle"], 
		position_angle_err = read_error("position_angle"), 
		ellipticity = obs_props["ellipticity"],
		ellipticity_err = read_error("ellipticity"),
	)
	Nb = length(bins) - 1

	d_xi = rand(Normal(0.0, pos_err))
	d_eta = rand(Normal(0.0, pos_err))
	pos_ang ~ Normal(position_angle, position_angle_err)
	ell ~ truncated(Normal(ellipticity, ellipticity_err), lower=0, upper=0.99)
	
	params ~ filldist(LogitNormal(-6.0, 5.0), Nb)

	xi = data.xi .+ d_xi
	eta = data.eta .+ d_eta

	radii = 60LilGuys.calc_r_ell(xi, eta, ell, pos_ang)
	
	r_b = DE.bin_indices(radii, bins)

	f = params[r_b]

	LL = sum(@. log10.(f*data.L_PM_SAT*data.L_CMD_SAT + (1-f) * data.L_PM_BKD*data.L_CMD_BKD))
	Turing.@addlogprob!(LL)
end

# ╔═╡ 51b06c5f-47e2-4734-a819-91393900f51f
model_extra = hist_model_extra(best_stars, bins_robust)

# ╔═╡ 2c58b55d-e0c6-4998-9b7b-614593ac8e00
chain_extra = sample(model_extra, NUTS(0.25), 100)

# ╔═╡ 969ccf42-1d72-4a21-98f6-041b5cac8666
summarize(chain_extra)

# ╔═╡ Cell order:
# ╠═cb7f2108-173a-4882-a95a-9111d1330176
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═109b91c2-c0ad-4683-849b-a56da5766ee4
# ╠═2f62c5c2-e397-463b-9e73-f554c31a7b85
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═57a19d65-59d9-46cf-8916-d9ac3a4dc92b
# ╠═28550d54-ad88-4df5-ac04-f46a15588ff8
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═2ab34a36-542c-41da-a83d-6eb67cae253c
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╠═04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╠═91867748-9b36-4f62-9310-8b778935776b
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╠═5feaa9b0-ca44-48d6-80d1-16303a34005f
# ╠═ab124afe-c94a-44d9-8869-bd4d4cee3fbd
# ╠═babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
# ╠═82f20e46-90c7-4925-83c9-3f49a909664d
# ╠═dd1b6b08-8821-4c3d-b937-c0d0705d2c33
# ╠═05333927-2605-4e8c-afd4-a2bfb120b8ce
# ╠═69b19f98-7a40-4874-821e-5dcf58fa70a6
# ╠═b14e1b89-6cb2-4c26-b0fe-1e2335b07d05
# ╠═66dd7757-7f10-4202-ad5f-a193b57cdf41
# ╠═c5c505d3-ed61-41bc-aaa5-b6c7037273ab
# ╠═b393ac76-4961-41ed-b847-900ee81364b1
# ╠═56b1b79d-7990-4a88-8ce0-273b211cbf75
# ╠═cc7b4837-1974-4098-86d6-95b1dace8a72
# ╠═c17a4353-564c-45d2-a561-ac71fc340baf
# ╠═df7cb725-380e-409f-bad6-1f1fca77a5b5
# ╠═85f68616-28da-4d2d-83cc-3d59e6c433f1
# ╠═6136af16-85ed-4f30-b995-eb6c1f4b6af6
# ╠═0793f4fd-0e8c-4afb-a4cc-e9ded55bebe3
# ╠═1ff82516-8b0a-4a95-bbca-f0230926f136
# ╠═5a272c79-0f27-40d2-a197-2b6db3f937be
# ╠═18e0757d-c4eb-469d-b601-119d3f03fb7c
# ╠═3040d802-1c29-4d7d-874e-c6f0b198baee
# ╠═79f36951-8b07-4569-a792-aba6db9c3d41
# ╠═455471ee-bf3e-441f-bd9d-b53252c5639f
# ╠═546ac410-5695-408c-81e0-d041fa12750d
# ╠═9846b747-2d15-408b-858d-b7e7f8229a04
# ╠═8179ee4d-ece5-41ee-9814-bda858e0a3a8
# ╠═381bb765-9c64-4a7e-bac6-7b69de22b1db
# ╠═efeb208c-afb8-4165-8897-7d356cf1a5ce
# ╠═3870d70e-5c09-491e-be99-cf0aadca6ed6
# ╠═1b189fd5-e9c4-47d0-b389-04d2a66a307c
# ╠═28fe8c1b-de0b-4158-a883-bde96ff92825
# ╠═5d3b661b-ee45-42e5-89b1-3283e9b66c14
# ╠═9aa3ffb7-17ce-47ea-9bac-a27d656fa4a2
# ╠═4ff739d1-0515-48ba-94e2-2012fd0bf865
# ╠═b77d9b07-a50f-40d3-b50d-8e2935ea6476
# ╠═12b74feb-1b58-4d31-bfd5-a9e0bb3768ef
# ╠═4632607a-ead9-4e3c-9ea4-e7959c7affc0
# ╠═a070ea1d-ecd6-41fc-afc8-21430e800b77
# ╠═6bf9eb27-7654-4196-a6f6-5b232db348c7
# ╠═7d80653c-83c6-4e28-bc88-48cceeda297a
# ╠═1765d2dd-9e59-480d-a636-7a89147fc52b
# ╠═e225d7ec-29fa-4b33-85b7-c453c1c5a639
# ╠═12a55136-bb11-4fd4-b76f-b71c360dca67
# ╠═199a2765-c2f8-46dc-826b-d60875df9bb2
# ╠═c3827324-b69c-4ddf-83cc-3568ef9b5229
# ╠═5fdeab29-1424-47ae-b95a-53fc3d6daf86
# ╠═f83d1311-af43-4b03-82aa-04b6e69f498b
# ╠═e4ab0464-af03-4e05-86b2-1b542d2c9334
# ╠═eca07277-1540-4102-bef0-214e40123a98
# ╠═63dbc12e-fb6d-4ef6-9645-b217298f6cae
# ╠═24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
# ╠═1d5de2ad-db30-4d18-845c-687bdf7ff2d1
# ╠═98c06e5a-9789-4bd8-9e5b-0b39b5230186
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═08424d6a-6a21-45d2-9b1a-d28e129c8158
# ╠═42844a8c-c8fc-4c6d-9826-02b88f6eff62
# ╠═af3db174-950d-46be-8d9c-61df70beb40e
# ╠═7c02d24f-7d59-410b-99bc-93513afa4c5d
# ╠═75930ff5-1add-4b1c-a495-64c50835ec4b
# ╠═bb97c5cc-35f4-42db-b5fa-c087e4527e0e
# ╠═01aa60ce-3849-4f3e-81d9-090c6aa8d104
# ╠═51b06c5f-47e2-4734-a819-91393900f51f
# ╠═2c58b55d-e0c6-4998-9b7b-614593ac8e00
# ╠═969ccf42-1d72-4a21-98f6-041b5cac8666
