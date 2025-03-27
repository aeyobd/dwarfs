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
FIGDIR = joinpath(galaxy, "figures"); FIGSUFFIX=".mcmc_hist_fast"

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
all_stars = read_gaia_stars(GaiaFilterParams(obs_props, filename = datafile,))

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
N_per_bin_min = max(round(Int, LilGuys.Interface.default_n_per_bin(best_stars.R_ell[best_stars.PSAT .> 0.2], nothing)), 2)

# ╔═╡ 1a5e8de2-0d74-48a9-aead-0855602734f3
bins = 10 .^ LilGuys.Interface.bins_both(log10.(best_stars.R_ell), nothing, bin_width=bin_width_min, num_per_bin=N_per_bin_min)


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
	
	stephist!(log10.(best_stars.R_ell[best_stars.PSAT .> 0.2]), bins=log10.(bins), label="PSAT > 0.2")
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

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
#  model
"""

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function hist_model(Lsat::Vector{Float64}, Lbkd::Vector{Float64})
	f ~ LogitNormal(-6.0, 5.0)


	LL = sum(@. log10(f*Lsat + (1-f) * Lbkd) )
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 562fe2c7-3ed4-4adb-af7f-1eca2aa35d8b
Nbins = length(bins) - 1	

# ╔═╡ ba5ceef7-ac10-46f7-a89b-30f38b6ddec1
r_b = DE.bin_indices(best_stars.R_ell, bins)

# ╔═╡ df84d941-7ddf-4c1c-96bb-1858b49bb710
Lsat = best_stars.L_CMD_SAT .* best_stars.L_PM_SAT

# ╔═╡ 97cbf4b6-025d-461b-b82e-044f86b713c1
Lbg = best_stars.L_CMD_BKD .* best_stars.L_PM_BKD

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range on either end to allow shifts in parameter to matter less.
"""

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
begin 
	chains = []
	for i in 1:Nbins
		filt = r_b .== i
		model = hist_model(Lsat[filt], Lbg[filt])
		chain = sample(model, NUTS(0.65), MCMCThreads(), 1000, 4) 
		push!(chains, chain)
	end
end

# ╔═╡ 500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
logit(x) = 1/(1 + exp(-x))

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fig = Figure(size=(3*72, 10*72))

	Nvar = length(bins)-1
	for i in 1:Nvar
		ax = Axis(fig[i, 1])
		for c in 1:size(chains[1], 3)
			y = chains[i][:, 1, c]
			lines!(logit.(y), linewidth=0.1)
		end

		if i < Nvar
			hidexdecorations!(ax)
		end
	end
	rowgap!(fig.layout, 0)

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

# ╔═╡ 639c7d31-8693-45dd-88be-492b124804e9
infoout = joinpath(galaxy, "processed", "info$FIGSUFFIX.toml")

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
begin 
	df_out = DataFrame()
	df_out[!, "iteration"] = DataFrame(chains[1])[!, "iteration"]
	df_out[!, "chain"] = DataFrame(chains[1])[!, "chain"]
	df_out[!, "lp"] .= 0

	for i in eachindex(chains)
		df = DataFrame(chains[i])
		df_out[!, "params[$i]"] = df[!, "f"]
		df_out[!, "lp"] .+= df.lp
		@assert df_out.chain == df.chain
		@assert df_out.iteration == df.iteration
	end

	df_out
end

# ╔═╡ 5e084776-a5f7-444a-9516-178cf610f584
vcat([DataFrame(summarize(chain)) for chain in chains]...)

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 77f564ee-f461-4107-b042-ae6d14ba9867
open(infoout, "w") do f
	TOML.print(f, OrderedDict(
		"N per bin" => N_per_bin_min,
		"min binwidth" => bin_width_min,
		"bins" => bins,
		"pos_err" => pos_err, 
		"position_angle" => obs_props["position_angle"], 
		"ellipticity" => obs_props["ellipticity"],
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
# ╠═24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═562fe2c7-3ed4-4adb-af7f-1eca2aa35d8b
# ╠═ba5ceef7-ac10-46f7-a89b-30f38b6ddec1
# ╠═df84d941-7ddf-4c1c-96bb-1858b49bb710
# ╠═97cbf4b6-025d-461b-b82e-044f86b713c1
# ╠═1d5de2ad-db30-4d18-845c-687bdf7ff2d1
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═500f6f5f-f1dc-4739-a40e-fa0b9fe78a51
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═639c7d31-8693-45dd-88be-492b124804e9
# ╠═a562c141-8626-4242-9f6a-c381a4da619b
# ╠═ed16c378-2786-4f83-8152-144838dedd19
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═5e084776-a5f7-444a-9516-178cf610f584
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═77f564ee-f461-4107-b042-ae6d14ba9867
