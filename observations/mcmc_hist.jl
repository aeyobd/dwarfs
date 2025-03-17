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

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
include("../utils/gaia_filters.jl")

# ╔═╡ 109b91c2-c0ad-4683-849b-a56da5766ee4
import StatsBase: sem

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
@bind galaxy TextField(default="ursa_minor")

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

# ╔═╡ 1d9b3718-45d7-4765-ac7d-017dbf939ec8
md"""
# Histogram model
"""

# ╔═╡ e7bf7a13-afca-49b8-8771-f7914adb347b
function hist_fractions(bins, params, x)
	i = DE.bin_indices(x, bins)
	return params[i]
end

# ╔═╡ 05b5489b-3e2f-4641-b079-9ced3677effe
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

# ╔═╡ 27c5a02f-0be6-4d0b-9c7f-99be12912732
pos_err = sem(filter(r->r.PSAT .> 0.2, best_stars).xi)

# ╔═╡ 1a5e8de2-0d74-48a9-aead-0855602734f3
bins = 10 .^ LilGuys.Interface.bins_equal_width(log10.(best_stars.r_ell), nothing, bin_width=0.1)

# ╔═╡ 10671310-6c4b-4c89-84ce-f7713dc778e9
model = hist_model(best_stars, bins)

# ╔═╡ 46bc408f-61f6-406c-aa18-061cfe8c65e2
chain = sample(model, NUTS(0.65), 100)

# ╔═╡ b3aa2259-645c-45d5-9b87-a876ae34ac91
summarize(chain)

# ╔═╡ 84ebf8d9-70a9-48c4-abb7-a058e4a0133c
@time LilGuys.calc_r_ell(best_stars.xi, best_stars.eta, obs_props["ellipticity"], obs_props["position_angle"])

# ╔═╡ 313e2b97-5d18-4c47-8cbb-68a018768976
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

# ╔═╡ 64c31de5-41fa-47e8-9209-61901e05f1be
bins

# ╔═╡ 9ea0c184-4a09-4293-a6d7-a963cf1d3d58
maximum(best_stars.xi .⊕ best_stars.eta) * 60

# ╔═╡ 93a1a541-5b01-4a14-8445-3acb67f37e8c
maximum(best_stars.r_ell)

# ╔═╡ 220b9de2-c87b-4b8f-a1ae-e25166638f7e
size(chain)

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

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat)

		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.log_Sigma, color=:black, alpha=0.1, markersize=1)

		push!(psats_s, psat)
	end

	psats_s = hcat(psats_s...)
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

	d_xi ~ Normal(0.0, pos_err)
	d_eta ~ Normal(0.0, pos_err)
	pos_ang ~ Normal(position_angle, position_angle_err)
	ell ~ truncated(Normal(ellipticity, ellipticity_err), lower=0)
	
	θ ~ filldist(LogitNormal(-2.0, 2.0), Nb)
	
	xi = data.xi .+ d_xi
	eta = data.eta .+ d_eta
	radii = LilGuys.calc_r_ell(xi, eta, ell, pos_ang)
	
	r_b = DE.bin_indices(radii, bins)
	f = θ[r_b]

	LL = sum(@. log10.(f*data.L_PM_SAT*data.L_CMD_SAT + (1-f) * data.L_PM_BKD*datafile.L_CMD_BKD))
	Turing.@addlogprob!(LL)
end

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range since stars can shift in r ell a bit more now.
"""

# ╔═╡ 5feaa9b0-ca44-48d6-80d1-16303a34005f
prof_simple = LilGuys.StellarProfile(best_stars.r_ell[best_stars.PSAT .> 0.2], bins=log10.(bins)[bins .< r_max])

# ╔═╡ 82f20e46-90c7-4925-83c9-3f49a909664d
sum(prof_simple.counts)

# ╔═╡ 84b6fdeb-3d30-4a3c-aab3-01e5e3fd2cde
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

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat)

		scatter!(h.log_r .+ 0.01*randn(length(h.log_r)), h.log_Sigma, color=:black, alpha=0.1, markersize=1)

		push!(psats, psat)
	end

	errorscatter!(prof_simple.log_r, prof_simple.log_Sigma, yerror=prof_simple.log_Sigma_err)

	psats = hcat(psats...)
	fig
end

# ╔═╡ 7478087f-a618-4e18-9597-58e22872f811
mean(hcat(psats...), dims=2)

# ╔═╡ da737b3e-f28a-4917-91d6-eea97458ddb0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "psat jax", ylabel = "psat me")

	skip = 10

	for i in 1:skip:length(chain)
		scatter!(best_stars.PSAT, psats[:, i], color=:black, alpha=0.1, markersize=5, rasterize=true)
	end
	
	fig
end

# ╔═╡ 10633dd3-857f-40a3-aa46-eabbd838b8a8
psat_mean = dropdims(mean(psats, dims=2), dims=2)

# ╔═╡ ad85b789-8f43-489b-a41c-a966e08d78ad
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 10

	scatter!(log10.(best_stars.r_ell), psat_mean.- best_stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
	fig
end

# ╔═╡ 7a2b2b3f-e3ab-4c49-a6d3-baec6688ec80
scatter(psat_mean, psat_mean .- best_stars.PSAT)

# ╔═╡ 35fa912d-953f-4d71-b4de-32cc30e307ba
psat_err = dropdims(std(psats, dims=2), dims=2)

# ╔═╡ 4a58490d-bc6f-4e5a-86f9-fed2d2842f5b
scatter(psat_mean, psat_err)

# ╔═╡ 124ec9a7-2b03-4c4d-8855-a3fbbd3d36c3
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

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat)

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

		h = LilGuys.StellarProfile(best_stars.r_ell, bins=log10.(bins), weights=psat)

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

# ╔═╡ 98c06e5a-9789-4bd8-9e5b-0b39b5230186
bins_robust = [0; bins; r_max]

# ╔═╡ cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
model_robust = hist_model_robust(best_stars, bins_robust)

# ╔═╡ 60d6f866-a05f-4233-b207-335dda4a82c4


# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
chain_robust = sample(model_robust, NUTS(0.65), 100)

# ╔═╡ af3db174-950d-46be-8d9c-61df70beb40e
truncated(Normal(), lower=0)

# ╔═╡ 984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
r_max = 240 * sqrt(1 - obs_props["ellipticity"])

# ╔═╡ 60fefea3-642d-4302-aace-afac44184426
# ╠═╡ disabled = true
#=╠═╡
r_max = maximum(all_stars.r_ell) * (1 + obs_props["ellipticity"])
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═109b91c2-c0ad-4683-849b-a56da5766ee4
# ╠═2f62c5c2-e397-463b-9e73-f554c31a7b85
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
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
# ╠═82f20e46-90c7-4925-83c9-3f49a909664d
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═e7bf7a13-afca-49b8-8771-f7914adb347b
# ╠═05b5489b-3e2f-4641-b079-9ced3677effe
# ╠═27c5a02f-0be6-4d0b-9c7f-99be12912732
# ╠═1a5e8de2-0d74-48a9-aead-0855602734f3
# ╠═10671310-6c4b-4c89-84ce-f7713dc778e9
# ╠═46bc408f-61f6-406c-aa18-061cfe8c65e2
# ╠═b3aa2259-645c-45d5-9b87-a876ae34ac91
# ╠═84ebf8d9-70a9-48c4-abb7-a058e4a0133c
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═313e2b97-5d18-4c47-8cbb-68a018768976
# ╠═64c31de5-41fa-47e8-9209-61901e05f1be
# ╠═9ea0c184-4a09-4293-a6d7-a963cf1d3d58
# ╠═93a1a541-5b01-4a14-8445-3acb67f37e8c
# ╠═984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
# ╠═84b6fdeb-3d30-4a3c-aab3-01e5e3fd2cde
# ╠═124ec9a7-2b03-4c4d-8855-a3fbbd3d36c3
# ╠═7478087f-a618-4e18-9597-58e22872f811
# ╠═da737b3e-f28a-4917-91d6-eea97458ddb0
# ╠═10633dd3-857f-40a3-aa46-eabbd838b8a8
# ╠═35fa912d-953f-4d71-b4de-32cc30e307ba
# ╠═ad85b789-8f43-489b-a41c-a966e08d78ad
# ╠═4a58490d-bc6f-4e5a-86f9-fed2d2842f5b
# ╠═7a2b2b3f-e3ab-4c49-a6d3-baec6688ec80
# ╠═220b9de2-c87b-4b8f-a1ae-e25166638f7e
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
# ╠═60fefea3-642d-4302-aace-afac44184426
# ╠═60d6f866-a05f-4233-b207-335dda4a82c4
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═af3db174-950d-46be-8d9c-61df70beb40e
