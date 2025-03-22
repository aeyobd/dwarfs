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
	galaxy = "leo1"
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 109b91c2-c0ad-4683-849b-a56da5766ee4
import StatsBase: sem

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(galaxy, "figures"); FIGSUFFIX=".mcmc_hist"

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

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

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

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type="svg", pt_per_unit=2)

# ╔═╡ cd338225-85b9-45db-86cc-979d5b77a532
Arya.update_figsize!(3.25)

# ╔═╡ 33e02945-2ab5-4f55-b3d5-e8855d86f120
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		yscale = log10,
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = log_r_label,
		ylabel = "counts",
		limits=(nothing, nothing, 1, 3e4),
		xlabelsize=12
	)
	
	stephist!(log10.(best_stars.r_ell[best_stars.PSAT .> 0.2]), bins=log10.(bins), label="PSAT > 0.2")
	stephist!(log10.(best_stars.r_ell), bins=log10.(bins), label="best")

	axislegend(position=:lt)
	@savefig "hist_bin_choices"
	fig
end

# ╔═╡ 64c31de5-41fa-47e8-9209-61901e05f1be
bins

# ╔═╡ 9ea0c184-4a09-4293-a6d7-a963cf1d3d58
maximum(best_stars.xi .⊕ best_stars.eta) * 60

# ╔═╡ 984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
r_max = maximum(best_stars.xi .⊕ best_stars.eta) * 60 * sqrt(1 - obs_props["ellipticity"])

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

# ╔═╡ 40c85a03-c607-433c-b902-13a4948f2a15
function perturbed_radii(data; 
		pos_err=pos_err, 
		position_angle=obs_props["position_angle"], 
		position_angle_err = read_error("position_angle"), 
		ellipticity = obs_props["ellipticity"],
		ellipticity_err = read_error("ellipticity"),
	)

	d_xi = rand(Normal(0.0, pos_err))
	d_eta = rand(Normal(0.0, pos_err))
	pos_ang = rand(Normal(position_angle, position_angle_err))
	ell = rand(truncated(Normal(ellipticity, ellipticity_err), lower=0, upper=0.99))
	
	xi = data.xi .+ d_xi
	eta = data.eta .+ d_eta

	radii = 60LilGuys.calc_r_ell(xi, eta, ell, pos_ang)

	return radii
end

# ╔═╡ 82b3b170-729f-463e-9063-756a34840555
perturbed_radii(best_stars)

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function hist_model_robust(data, bins; kwargs...)
	Nb = length(bins) - 1

	radii = perturbed_radii(data; kwargs...)
	
	params ~ filldist(LogitNormal(-6.0, 5.0), Nb)
	r_b = DE.bin_indices(radii, bins)

	f = params[r_b]

	LL = sum(@. log10.(f*data.L_PM_SAT*data.L_CMD_SAT 
		+ (1-f) * data.L_PM_BKD*data.L_CMD_BKD
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range on either end to allow shifts in parameter to matter less.
"""

# ╔═╡ 98c06e5a-9789-4bd8-9e5b-0b39b5230186
begin 
	bins_robust = copy(bins)
	bins_robust[1] = 0
	bins_robust[end]=Inf
	bins_robust
end

# ╔═╡ cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
model_robust = hist_model_robust(best_stars, bins_robust)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
chain_robust = sample(model_robust, NUTS(0.25), MCMCThreads(), 1000, 4) 

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
summarize(chain_robust)

# ╔═╡ 500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
logit(x) = 1 / (1 + exp(-x))

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fig = Figure(size=(3*72, 10*72))

	Nvar = length(bins_robust)-1
	for i in 1:Nvar
		ax = Axis(fig[i, 1])
		for c in 1:size(chain_robust, 3)
			y = chain_robust[:, i, c]
			lines!(logit.(y), linewidth=0.1)
		end

		if i < Nvar
			hidexdecorations!(ax)
		end
	end
	rowgap!(fig.layout, 0)

	@savefig "chains"
	
	fig

end

# ╔═╡ 42844a8c-c8fc-4c6d-9826-02b88f6eff62
log10.(bins_robust)

# ╔═╡ af3db174-950d-46be-8d9c-61df70beb40e
truncated(Normal(), lower=0)

# ╔═╡ 36bb1876-5937-44fe-bb1c-add10470371d
pvalue = 0.16

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
		y = median(chain_robust[:, i, :])
		push!(xs, x)
		push!(ys, y)
		ye = quantile(chain_robust[:, i, :], [pvalue, 1-pvalue])
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

# ╔═╡ b1487f12-76b2-4a2f-8d50-4e80f90b1d9d
binsout = joinpath(galaxy, "processed", "info$FIGSUFFIX.toml")

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
df_out = DataFrame(chain_robust)

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 77f564ee-f461-4107-b042-ae6d14ba9867
open(infoout, "w") do f
	println(f)
	println(f, "N per bin: $N_per_bin_min")
	println(f, "min binwidth: $bin_width_min")
	println(f)
	
	println(f, "bins robust")
	println(f, bins_robust)


	println(f)
	println(f, "chains")
	println(f, DataFrame(summarize(chain_robust)))
end

# ╔═╡ ebdfb5dc-b996-4ead-a61f-ebb9f93adf0e
open(binsout, "w") do f
	TOML.print(f, OrderedDict(
		"N per bin" => N_per_bin_min,
		"min binwidth" => bin_width_min,
		"bins" => bins_robust,
		"pos_err" => pos_err, 
		"position_angle" => obs_props["position_angle"], 
		"position_angle_err" => read_error("position_angle"), 
		"ellipticity" => obs_props["ellipticity"],
		"ellipticity_err" => read_error("ellipticity"),
	)
	)
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
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╠═04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
# ╠═91867748-9b36-4f62-9310-8b778935776b
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═e7bf7a13-afca-49b8-8771-f7914adb347b
# ╠═27c5a02f-0be6-4d0b-9c7f-99be12912732
# ╠═1a5e8de2-0d74-48a9-aead-0855602734f3
# ╠═ccd115f0-5a6e-4a65-a45f-d97cdbbbe02c
# ╠═8d3bbb76-47e1-4a36-981f-487b657bc74a
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═cd338225-85b9-45db-86cc-979d5b77a532
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═64c31de5-41fa-47e8-9209-61901e05f1be
# ╠═9ea0c184-4a09-4293-a6d7-a963cf1d3d58
# ╠═984af5ed-2c77-4f85-9c9e-3f0ed4f36fcf
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═e87b5ee8-e281-43e3-8172-4848018e2cd2
# ╠═c2e842eb-06dd-48da-8948-366983f240a9
# ╠═40c85a03-c607-433c-b902-13a4948f2a15
# ╠═82b3b170-729f-463e-9063-756a34840555
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═cee7c244-42e1-4aae-93d4-5a4e8e5f4b09
# ╠═1d5de2ad-db30-4d18-845c-687bdf7ff2d1
# ╠═98c06e5a-9789-4bd8-9e5b-0b39b5230186
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
# ╠═42844a8c-c8fc-4c6d-9826-02b88f6eff62
# ╠═af3db174-950d-46be-8d9c-61df70beb40e
# ╠═36bb1876-5937-44fe-bb1c-add10470371d
# ╠═75930ff5-1add-4b1c-a495-64c50835ec4b
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═639c7d31-8693-45dd-88be-492b124804e9
# ╠═b1487f12-76b2-4a2f-8d50-4e80f90b1d9d
# ╠═bf8d5936-1a4e-47c2-bb22-531ab344b8ad
# ╠═a562c141-8626-4242-9f6a-c381a4da619b
# ╠═ed16c378-2786-4f83-8152-144838dedd19
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═77f564ee-f461-4107-b042-ae6d14ba9867
# ╠═ebdfb5dc-b996-4ead-a61f-ebb9f93adf0e
