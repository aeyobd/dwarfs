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

# ╔═╡ 48c2def6-574e-436b-8b9e-a608e68aec5c
begin
	import Pkg; Pkg.activate()

	using Arya, CairoMakie
	CairoMakie.activate!(pt_per_unit=5, type="svg")
end

# ╔═╡ 5d232d36-a6cf-41a5-8f18-c676e9079ce0
using OrderedCollections

# ╔═╡ b078d7d2-71fb-4f7b-ba25-6b488976f0a6
using PlutoUI

# ╔═╡ bd452643-7d45-4a85-9a49-e259a9eec32e
using Distributions

# ╔═╡ 1a18648e-881e-479c-baac-5abf7e7766ec
md"""
# Robust density profile

The goal of this notebook is to calculate a density profile while sampling over observational uncertanties to better represent the errors
"""

# ╔═╡ 8252b891-f46d-448d-99ae-93d2892af3c5
 md"""
 Galaxyname = $(@bind galaxyname confirm(PlutoUI.TextField(default="draco")))
 """

# ╔═╡ f5953fef-726a-45d0-b3e9-a23ddfe642f0
begin 
	import LilGuys as lguys
	using LilGuys
	FIGDIR = joinpath(galaxyname, "figures")
	FIGSUFFIX = ".profile"
end

# ╔═╡ f0a80ece-1889-4894-89b6-110ed593348f
Arya.update_figsize!(3.25)

# ╔═╡ a4e53ac6-c989-4fd8-ad4e-8b39cdaf8ead
import TOML

# ╔═╡ d8485cff-9df4-4325-8771-68a34b9da9bc
if isfile("$galaxyname/data/jensen+24_2c.fits")
	datafile = "$galaxyname/data/jensen+24_2c.fits"
elseif isfile("$galaxyname/data/j24_2c.fits")
	datafile = "$galaxyname/data/j24_2c.fits"
elseif isfile("$galaxyname/data/j24_1c.fits")
	datafile = "$galaxyname/data/j24_1c.fits"
else
	datafile = "$galaxyname/data/jensen+24_1c.fits"
end

# ╔═╡ e676c5b3-ab8c-46c5-a44f-6dc56aa9b51e
md"""
## Utility functions
"""

# ╔═╡ 03f3c61f-4d5a-4e13-a3c1-43821ecf80e1
function get_mean_err(d::Dict, key)
	mean = d[key]
	if "$(key)_em" ∈ keys(d)
		err = max(d["$(key)_em"], d["$(key)_ep"])
	else
		err = d[key * "_err"]
	end

	return mean, err
end

# ╔═╡ 75e4bca8-5664-4620-8bb4-a14ee7b6da3b
PSAT_dist = Beta(4, 5)

# ╔═╡ fdc148e5-4a39-45d2-9e08-046a3f75f118
mode(PSAT_dist)

# ╔═╡ d057d08d-7736-44f6-905f-2e784d2da509
PSAT_0 = median(PSAT_dist)

# ╔═╡ 03df7329-7d29-476a-a63d-7c6e344d6a31
quantile(PSAT_dist, [0.001, 0.999])

# ╔═╡ 720ebad9-d01c-4773-8b74-182541811138
plot(PSAT_dist)

# ╔═╡ 6bfe7b25-c82f-4153-b0ee-cd3ef23cb516


# ╔═╡ 3f7c4462-906f-46d2-9431-fb2ab12ddc4f
let
	fig = Figure()
	ax = Axis(fig[1,1])
	fig
end

# ╔═╡ 29bea80f-a782-4052-9025-fc6f83f44a21
observed_properties = TOML.parsefile(galaxyname * "/observed_properties.toml")

# ╔═╡ 5ef66a21-eee1-4d43-b57b-90770aa55840
log_r_ell_label = L"\log\ r_\textrm{ell}\,/\,\textrm{arcmin}"

# ╔═╡ df077acf-93d6-4277-9c4a-ac967b55c54a
log_Σ_label = L"$\log \Sigma$\,/\,stars\ arcmin$^{-2}$"

# ╔═╡ 8d13d560-6620-479b-98ca-f07554158f35
begin 
	all_stars = lguys.read_fits(datafile)

	all_stars[:, :r_ell_old] = all_stars[:, :r_ell]
	all_stars.xi, all_stars.eta = LilGuys.to_tangent(all_stars.ra, all_stars.dec, observed_properties["ra"], observed_properties["dec"])

	all_stars.r_ell = 60LilGuys.calc_r_ell(all_stars.xi, all_stars.eta, observed_properties["ellipticity"], observed_properties["position_angle"])
	all_stars
end

# ╔═╡ edbc5116-73c8-4fe8-b87f-0f0c635244e7
r_ell_max = (1-observed_properties["ellipticity"]) * maximum(all_stars.r_ell)

# ╔═╡ cc134060-b7dd-49dc-90db-097bfcc91ed2
members = all_stars[all_stars.PSAT .> PSAT_0, :]

# ╔═╡ e7314a16-e397-4f54-8279-b514d269e5f8
members_nospace = all_stars[all_stars.PSAT_NOSPACE .> PSAT_0, :]

# ╔═╡ c418a5ff-2b46-4842-831b-fd244302544a
cen = lguys.to_tangent(lguys.calc_centre2D(members.ra, members.dec, "mean")..., observed_properties["ra"], observed_properties["dec"])

# ╔═╡ cdb7dbee-b395-4991-9b64-26632f568b31
md"""
each component of the cen_errs tuple below should be similar in magnitude
"""

# ╔═╡ 289e86f4-312c-4f6e-9a1f-4e0ccd3d30be
cen_errs =  (cen..., 
lguys.std(members.eta) / sqrt(length(members.eta)), 
lguys.std(members.xi) / sqrt(length(members.xi))
)

# ╔═╡ a10a7b16-1c40-4e86-9ae6-b01aa76a223f
cen_err = maximum(abs.(cen_errs))

# ╔═╡ 5fa9d63d-1ce8-4c1e-8871-85e5017abadd
function get_num_per_bin(x)
	num_per_bin = round(Int64, LilGuys.Interface.default_n_per_bin(x))
	if length(x) < 2
		num_per_bin = 1
	end
	return num_per_bin
end


# ╔═╡ 72d6571b-9bb6-4eed-bd21-c2307dd38b8c
function get_bins(r_ell; bin_width = 0.05, num_per_bin=nothing)

	if num_per_bin === nothing
		num_per_bin = get_num_per_bin(log10.(r_ell))
	end

	num_per_bin = round(Int, num_per_bin)
	bins = LilGuys.Interface.bins_both(log10.(r_ell), nothing, bin_width=bin_width, num_per_bin=num_per_bin)
	bins[end] = min(bins[end], log10(r_ell_max))

	return bins
end

# ╔═╡ 0c483a09-1acb-4c66-9937-8a3a79c62ed3
bins = get_bins(members.r_ell)

# ╔═╡ 444b9717-e774-4b60-b6b4-13bb222dbce1
bins_nospace = get_bins(members_nospace.r_ell)

# ╔═╡ 893c53db-4f6b-4aee-a472-ccc3df73ea9e
function perturbed_profile(; d_xi=0, d_eta=0, psat_min=PSAT_0, d_θ=0, d_ell=0.0,
		bins=bins, psat_col=:PSAT, 
	)
	
	members = all_stars[all_stars[!, psat_col] .> psat_min, :]
	xi, eta = members.xi, members.eta

	ell =  observed_properties["ellipticity"] + d_ell
	θ = observed_properties["position_angle"] + d_θ

	xi_s = xi .+ d_xi
	eta_s =  eta .+ d_eta
	
	x = 60 * lguys.calc_r_ell(xi_s, eta_s, ell, θ) 
	
	x = x[r_ell_max .> x .> 0]
	
	prof = LilGuys.StellarProfile(x, normalization=:none, bins=bins)

	return prof
end

# ╔═╡ 478eae36-58ff-41db-b5e2-906806a0e082
ell_0, ell_err = get_mean_err(observed_properties, "ellipticity")

# ╔═╡ 6bb428cc-dc9c-4825-9e67-4c4092e480e4
theta_0, theta_err = get_mean_err(observed_properties, "position_angle")

# ╔═╡ faa38734-8998-40df-a170-ac6162d146a6
prof_fiducial = LilGuys.StellarProfile(members.r_ell, normalization=:none, bins=(bins))

# ╔═╡ d92ab14c-81b8-48d7-874f-152d2225e3cc
prof_nospace = LilGuys.StellarProfile(all_stars[all_stars.PSAT_NOSPACE .> PSAT_0, :r_ell], normalization=:none, bins=(bins))

# ╔═╡ c238d456-287b-4d3a-86f1-0b790a85a145
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	sum_sq = zeros(length(prof_fiducial.log_r))
	N = 1000
	for _ in 1:N
		
		prof = perturbed_profile(
			d_ell=randn() * ell_err,
			d_θ = randn() * theta_err,
			d_xi = randn() * cen_err, 
			d_eta = randn() * cen_err, 
			psat_min = rand(PSAT_dist)
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.01))
		sres = @.  (prof.log_Sigma - prof_fiducial.log_Sigma)^2
		sres[isnan.(sres)] .= prof_fiducial.log_Sigma[isnan.(sres)] .^ 2
		sum_sq .+= sres
	end

	y_err_tot = @. sqrt(sum_sq / N)
	
	prof = prof_fiducial
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=y_err_tot, color=COLORS[2])
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[3])

	LilGuys.hide_grid!(ax)
	
	@savefig "density_profiles_mc"
end

# ╔═╡ 92b0d9d0-ddc2-4726-b771-9fa9ae486547
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	bins = bins_nospace
	sum_sq = zeros(length(prof_fiducial.log_r))
	N = 1000
	for _ in 1:N
		
		prof = perturbed_profile(
			d_ell=randn() * ell_err,
			d_θ = randn() * theta_err,
			d_xi = randn() * cen_err, 
			d_eta = randn() * cen_err, 
			psat_min = rand(PSAT_dist),
			psat_col = :PSAT_NOSPACE,
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.01))
		sres = @.  (prof.log_Sigma - prof_nospace.log_Sigma)^2
		sres[isnan.(sres)] .= prof_nospace.log_Sigma[isnan.(sres)] .^ 2
		sum_sq .+= sres
	end

	prof = prof_nospace

	y_err_tot = @. sqrt(sum_sq / N) .+ prof.log_Sigma_err
	
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=y_err_tot, color=COLORS[2])
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[3])

	LilGuys.hide_grid!(ax)
	
	@savefig "density_profiles_mc_nospace"
end

# ╔═╡ 55df8e1d-32ce-43ba-96de-fcfd4b8cbb7e
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	N = 1000

	all_log_Σ = Matrix{Float64}(undef, N, length(bins) - 1)
	all_log_Σ_err = Matrix{Float64}(undef, N, length(bins) - 1)

	for i in 1:N
		
		prof = perturbed_profile(
			d_ell=randn() * ell_err,
			d_θ = randn() * theta_err,
			d_xi = randn() * cen_err, 
			d_eta = randn() * cen_err, 
			psat_min = rand(PSAT_dist),
			psat_col = :PSAT_NOSPACE,
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.01))

		all_log_Σ[i, :] .= prof.log_Sigma
		all_log_Σ_err[i, :] .= prof.log_Sigma_err
		all_log_Σ_err[i, isnan.(prof.log_Sigma_err)] .= prof_nospace.log_Sigma[ isnan.(prof.log_Sigma_err)]
	end

	all_log_Σ[isnan.(all_log_Σ)] .= -Inf

	Σ_median = median(all_log_Σ, dims=1)
	Σ_ep = quantile.(eachcol(all_log_Σ), 0.84) .- Σ_median'
	Σ_em = Σ_median' .- quantile.(eachcol(all_log_Σ), 0.16)
	Σ_median = dropdims(Σ_median, dims=1)
	Σ_estat = dropdims(median(all_log_Σ_err, dims=1), dims=1)
	Σ_estat[isnan.(Σ_estat)] .= Σ_median[isnan.(Σ_estat)]
	println(size(Σ_ep))
	Σ_ep = dropdims(Σ_ep, dims=2)
	Σ_em = dropdims(Σ_em, dims=2)

	errorscatter!(prof.log_r, Σ_median, yerror=collect(zip(Σ_em .+ Σ_estat, Σ_ep.+ Σ_estat)), 
		color=COLORS[3],
		label="bootstrap", 
		size=3
	)

	prof = prof_nospace
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[2], label = "simple PSAT_NOSPACE", size=3,
	)

	axislegend(position=:lb)


	LilGuys.hide_grid!(ax)
	
	@savefig "density_profiles_mc_nospace"
end

# ╔═╡ ff72787b-75e9-4d3b-ac8d-a028bd5bc49c
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	N = 1000

	all_log_Σ = Matrix{Float64}(undef, N, length(bins) - 1)
	all_log_Σ_err = Matrix{Float64}(undef, N, length(bins) - 1)

	for i in 1:N
		
		prof = perturbed_profile(
			d_ell=randn() * ell_err,
			d_θ = randn() * theta_err,
			d_xi = randn() * cen_err, 
			d_eta = randn() * cen_err, 
			psat_min = rand(PSAT_dist),
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.01))

		all_log_Σ[i, :] .= prof.log_Sigma
		all_log_Σ_err[i, :] .= prof.log_Sigma_err
		all_log_Σ_err[i, isnan.(prof.log_Sigma_err)] .= prof_fiducial.log_Sigma[ isnan.(prof.log_Sigma_err)]
	end

	all_log_Σ[isnan.(all_log_Σ)] .= -Inf

	Σ_median = median(all_log_Σ, dims=1)
	Σ_ep = quantile.(eachcol(all_log_Σ), 0.84) .- Σ_median'
	Σ_em = Σ_median' .- quantile.(eachcol(all_log_Σ), 0.16)
	Σ_median = dropdims(Σ_median, dims=1)
	Σ_estat = dropdims(median(all_log_Σ_err, dims=1), dims=1)
	Σ_estat[isnan.(Σ_estat)] .= Σ_median[isnan.(Σ_estat)]
	println(size(Σ_ep))
	Σ_ep = dropdims(Σ_ep, dims=2)
	Σ_em = dropdims(Σ_em, dims=2)

	errorscatter!(prof.log_r, Σ_median, yerror=collect(zip(Σ_em .+ Σ_estat, Σ_ep.+ Σ_estat)), 
		color=COLORS[3],
		label = "bootstrap"
	)

	prof = prof_fiducial
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[2], label="simple PSAT")


	axislegend(position=:lb)
	LilGuys.hide_grid!(ax)
	
	@savefig "density_profiles_mc"
end

# ╔═╡ f7e0591e-fb61-4a7b-9b31-c84e9f677507
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
		title = "PSAT_only"
	)

	for _ in 1:1000
		
		prof = perturbed_profile(
			psat_min = rand(PSAT_dist)
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.03))
	end

	prof = prof_fiducial
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[2])

	LilGuys.hide_grid!(ax)

	fig
end

# ╔═╡ 90828e61-8be0-414c-b1b5-69575f0591ad
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
		title = "cen only"
	)

	sum_sq = zeros(length(prof_fiducial.log_r))
	N = 1000
	for _ in 1:N
		
		prof = perturbed_profile(
			d_xi = randn() * cen_err, 
			d_eta = randn() * cen_err, 
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.03))
		sum_sq .+= @. (prof.log_Sigma - prof_fiducial.log_Sigma)^2
	end

	y_err_tot = @. sqrt(sum_sq / N)
	
	prof = prof_fiducial
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=y_err_tot, color=COLORS[2])
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[3])

	LilGuys.hide_grid!(ax)
	
	fig
end

# ╔═╡ 422e9f8d-56ff-4198-be0e-296d47501a5d
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
		title = "structural only"
	)

	sum_sq = zeros(length(prof_fiducial.log_r))
	N = 1000
	for _ in 1:N
		
		prof = perturbed_profile(
			d_ell=randn() * ell_err,
			d_θ = randn() * theta_err,
		)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.03))
		sum_sq .+= @. (prof.log_Sigma - prof_fiducial.log_Sigma)^2
	end

	y_err_tot = @. sqrt(sum_sq / N)
	
	prof = prof_fiducial
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=y_err_tot, color=COLORS[2])
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[3])

	LilGuys.hide_grid!(ax)
	
	fig
end

# ╔═╡ 01e81c57-2466-4084-bd61-f094c25bc620
md"""
# Systematics
"""

# ╔═╡ 8e1c3bb1-1521-4419-b11b-194dc92ca5f1
best_stars = all_stars[all_stars.F_BEST .== 1., :]

# ╔═╡ 44d238cb-2ac9-4dc2-9a68-58037a32dbe9
η_bins_bg = 0.3

# ╔═╡ c5aa6063-9a61-4a7e-bc2a-5f3b047f1850
r_ell_min_bg = 3 * observed_properties["r_h"]

# ╔═╡ 1b9c9581-1f88-4ba0-979b-4e0efb4b794f
let
	global log_Σ_bg
	global log_Σ_bg_err

	filt_bg = best_stars.r_ell .< r_ell_max
	filt_bg .&= best_stars.r_ell .> r_ell_min_bg

	N = sum(filt_bg)

	A = π * ( r_ell_max^2 - r_ell_min_bg^2)
	log_Σ_bg = log10.(N / A)

	x = best_stars.r_ell[filt_bg]
	bins = get_bins(x)
	prof_bg = LilGuys.StellarProfile(x, normalization=:none, bins=(bins))

	
	log_Σ_bg_err = lguys.std(prof_bg.log_Sigma) / sqrt(length(prof_bg.log_Sigma)) .+ 1/sqrt(N) / log(10)
end

# ╔═╡ 21247186-a26f-4cce-97cb-32d0c8689dea
function plot_density_for_filter(filters; sequential=false, normalization=:none, axislabel="")
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	names = collect(keys(filters))
	for i in 1:length(filters)
		name = names[i]
		x = filters[name]
		
		x = x[r_ell_max .> x .> 0]
		if length(x) < 2
			println("skipping with < 2 members: ", "\t", name)
			continue
		end
		
		println(name, "\t", length(x))
		bins = get_bins(x)

		if length(bins) < 2
			println("skipping with < 2 bins: ", "\t", name)
			continue
		end
		
		prof = LilGuys.StellarProfile(x, normalization=normalization, bins=(bins))

		if sequential
			 kwargs = (; color=i, colorrange=(1, length(names)))
			else
			 kwargs = (; )
		end

		
		scatterlines!(prof.log_r, prof.log_Sigma, label=name; kwargs...)
	end

	axislegend(ax, axislabel, position=:lb)

	return fig
end

# ╔═╡ 5d55a3c9-5d9a-4a8a-a21e-c23616a0450e
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,

	)


	
	x = (best_stars.r_ell)
	x = x[r_ell_max .> x .> 0]
	prof = LilGuys.StellarProfile(x, normalization=:none, bins=get_bins(x, num_per_bin=η_bins_bg*get_num_per_bin(log10.(x))))
	
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, label="all", )

	y = 10 .^ prof.log_Sigma
	y_err = prof.log_Sigma_err .* y .* log(10)

	y_bg = 10 .^ log_Σ_bg
	y_bg_err = y_bg * log_Σ_bg_err * log(10)
	y_new = y .- y_bg

	filt = 1:(findfirst(y_new .<= 0)-1)
	if sum(filt) > 1
		@info filt
		@info y_new
	
		y_new = y_new[filt]
		y_err = y_err[filt]
	
		y_new_em = min.(log10.(y_new) .- log10.(max.(0, y_new .- y_err .- y_bg_err)), 10.0)
		y_new_ep = log10.(y_new .+ y_err .+ y_bg_err) .- log10.(y_new)
		y_err = collect(zip(y_new_em, y_new_ep))[filt]
	
		@info y_err
		
		errorscatter!(prof.log_r[filt], log10.(y_new), yerror=y_err, label="bg-subtracted")

	end

	hlines!(log_Σ_bg, color=:black, label="background")
	hspan!(log_Σ_bg - log_Σ_bg_err, log_Σ_bg + log_Σ_bg_err, color=(:black, 0.1), )

	filt = best_stars.PSAT .> 0.2
	x = (best_stars.r_ell[filt])
	x = x[r_ell_max .> x .> 0]
	prof = LilGuys.StellarProfile(x, normalization=:none, bins=get_bins(x))
	errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, label="PSAT > 0.2")

	axislegend(position=:lb)
	ylims!(ax, minimum(prof.log_Sigma .- 2prof.log_Sigma_err), nothing)
	@savefig "density_profiles_bg_subtraction"

	fig
end

# ╔═╡ 297fc423-8d87-432d-9f29-39ebf99f2047
let 
	fig = plot_density_for_filter(OrderedDict(
		"nonmembers" => best_stars.r_ell[best_stars.PSAT .< 0.1],
		"members" => best_stars.r_ell[best_stars.PSAT .> 0.2],
		"cmd" => best_stars.r_ell[best_stars.PSAT_CMD .> 0.2],
		"space" => best_stars.r_ell[best_stars.PSAT_S .> 0.2],
		"PM" => best_stars.r_ell[best_stars.PSAT_PM .> 0.2],
		"PM+CMD" => best_stars.r_ell[best_stars.PSAT_NOSPACE .> 0.2],
		
	),
)
	@savefig "density_profiles_by_component"

	fig
end

# ╔═╡ 971707bc-78b9-4d2f-8063-ee549497bac6
md"""
How much does PSAT affect the resulting density profile?
"""

# ╔═╡ 3b2ef255-9b36-4dae-88e2-4ea60555e197
let

	threshholds = [0.01, 0.1, 0.2, 0.5, 0.9, 0.99]

	fig = plot_density_for_filter(
		OrderedDict([string(pmin) => best_stars.r_ell[best_stars.PSAT .>= pmin] for pmin in threshholds]),
		sequential=true,
		axislabel="PSAT"
	)

	@savefig "density_profiles_by_PSAT"
	fig
end

# ╔═╡ 86a832d0-0771-472a-9cb5-d8dbbe058d3a
let

	
	threshholds = LilGuys.quantile(members.phot_g_mean_mag, LinRange(0, 1, 5))
	threshholds = round.(threshholds, digits=1)

	Nc = 4
	fig = plot_density_for_filter(
		OrderedDict(
			["$(threshholds[i]) - $(threshholds[i+1])" => members.r_ell[threshholds[i] .<= members.phot_g_mean_mag .< threshholds[i+1]] for i in 1:Nc]),
		sequential=true,
		axislabel="G mag"
	)

	@savefig "density_profiles_by_gmag"
	fig
end

# ╔═╡ bc7a8d01-3611-4ddb-bf9d-af9029f51d08
cen

# ╔═╡ 264f58ce-97cb-4e7b-a8cb-4934ada85e6b
let
	xi, eta = members.xi .- cen[1], members.eta .- cen[2]
	filters = OrderedDict(
		"r_ell" => members.r_ell_old .* observed_properties["r_h"] * sqrt(1 - observed_properties["ellipticity"]),
		"r_circ" => @.(60 * sqrt(members.xi^2 + members.eta^2)),
		"r_ell_me" => 60*lguys.calc_r_ell(xi, eta, observed_properties["ellipticity"], observed_properties["position_angle"]) ,
		"recentred" => 60*lguys.calc_r_ell(xi, eta, observed_properties["ellipticity"], observed_properties["position_angle"]) 
	)

	fig = plot_density_for_filter(filters)

	@savefig "density_profiles_r_method"

	fig
end

# ╔═╡ 6a33edf0-d555-49ad-8e3a-51da44ce8689
Arya.update_fontsize!(10)

# ╔═╡ a00d540b-d7ba-4ecc-878d-2737e7066e87
let

	fig = plot_density_for_filter(OrderedDict(
		"fiducial" => members.r_ell,
		"circular" => @.(60 * sqrt(members.xi^2 + members.eta^2)),
		"nospace" => members_nospace.r_ell,
		"bright" => members.r_ell[members.phot_g_mean_mag .< 20],
		"PSAT > 0.8" => members.r_ell[members.PSAT .> 0.8],
		),
		#normalization=:mass,
	)

	@savefig "density_profiles_medley"
	fig
end

# ╔═╡ Cell order:
# ╟─1a18648e-881e-479c-baac-5abf7e7766ec
# ╟─8252b891-f46d-448d-99ae-93d2892af3c5
# ╠═48c2def6-574e-436b-8b9e-a608e68aec5c
# ╠═f0a80ece-1889-4894-89b6-110ed593348f
# ╠═a4e53ac6-c989-4fd8-ad4e-8b39cdaf8ead
# ╠═5d232d36-a6cf-41a5-8f18-c676e9079ce0
# ╠═f5953fef-726a-45d0-b3e9-a23ddfe642f0
# ╠═b078d7d2-71fb-4f7b-ba25-6b488976f0a6
# ╠═bd452643-7d45-4a85-9a49-e259a9eec32e
# ╠═d8485cff-9df4-4325-8771-68a34b9da9bc
# ╟─e676c5b3-ab8c-46c5-a44f-6dc56aa9b51e
# ╠═03f3c61f-4d5a-4e13-a3c1-43821ecf80e1
# ╠═75e4bca8-5664-4620-8bb4-a14ee7b6da3b
# ╠═fdc148e5-4a39-45d2-9e08-046a3f75f118
# ╠═d057d08d-7736-44f6-905f-2e784d2da509
# ╠═03df7329-7d29-476a-a63d-7c6e344d6a31
# ╠═720ebad9-d01c-4773-8b74-182541811138
# ╠═6bfe7b25-c82f-4153-b0ee-cd3ef23cb516
# ╠═3f7c4462-906f-46d2-9431-fb2ab12ddc4f
# ╠═29bea80f-a782-4052-9025-fc6f83f44a21
# ╠═5ef66a21-eee1-4d43-b57b-90770aa55840
# ╠═df077acf-93d6-4277-9c4a-ac967b55c54a
# ╠═8d13d560-6620-479b-98ca-f07554158f35
# ╠═edbc5116-73c8-4fe8-b87f-0f0c635244e7
# ╠═cc134060-b7dd-49dc-90db-097bfcc91ed2
# ╠═e7314a16-e397-4f54-8279-b514d269e5f8
# ╠═c418a5ff-2b46-4842-831b-fd244302544a
# ╠═cdb7dbee-b395-4991-9b64-26632f568b31
# ╠═289e86f4-312c-4f6e-9a1f-4e0ccd3d30be
# ╠═a10a7b16-1c40-4e86-9ae6-b01aa76a223f
# ╠═72d6571b-9bb6-4eed-bd21-c2307dd38b8c
# ╠═5fa9d63d-1ce8-4c1e-8871-85e5017abadd
# ╠═0c483a09-1acb-4c66-9937-8a3a79c62ed3
# ╠═444b9717-e774-4b60-b6b4-13bb222dbce1
# ╠═893c53db-4f6b-4aee-a472-ccc3df73ea9e
# ╠═478eae36-58ff-41db-b5e2-906806a0e082
# ╠═6bb428cc-dc9c-4825-9e67-4c4092e480e4
# ╠═faa38734-8998-40df-a170-ac6162d146a6
# ╠═d92ab14c-81b8-48d7-874f-152d2225e3cc
# ╠═c238d456-287b-4d3a-86f1-0b790a85a145
# ╠═92b0d9d0-ddc2-4726-b771-9fa9ae486547
# ╠═55df8e1d-32ce-43ba-96de-fcfd4b8cbb7e
# ╠═ff72787b-75e9-4d3b-ac8d-a028bd5bc49c
# ╠═f7e0591e-fb61-4a7b-9b31-c84e9f677507
# ╠═90828e61-8be0-414c-b1b5-69575f0591ad
# ╠═422e9f8d-56ff-4198-be0e-296d47501a5d
# ╠═01e81c57-2466-4084-bd61-f094c25bc620
# ╠═8e1c3bb1-1521-4419-b11b-194dc92ca5f1
# ╠═1b9c9581-1f88-4ba0-979b-4e0efb4b794f
# ╠═44d238cb-2ac9-4dc2-9a68-58037a32dbe9
# ╠═c5aa6063-9a61-4a7e-bc2a-5f3b047f1850
# ╠═21247186-a26f-4cce-97cb-32d0c8689dea
# ╠═5d55a3c9-5d9a-4a8a-a21e-c23616a0450e
# ╠═297fc423-8d87-432d-9f29-39ebf99f2047
# ╠═971707bc-78b9-4d2f-8063-ee549497bac6
# ╠═3b2ef255-9b36-4dae-88e2-4ea60555e197
# ╠═86a832d0-0771-472a-9cb5-d8dbbe058d3a
# ╠═bc7a8d01-3611-4ddb-bf9d-af9029f51d08
# ╠═264f58ce-97cb-4e7b-a8cb-4934ada85e6b
# ╠═6a33edf0-d555-49ad-8e3a-51da44ce8689
# ╠═a00d540b-d7ba-4ecc-878d-2737e7066e87
