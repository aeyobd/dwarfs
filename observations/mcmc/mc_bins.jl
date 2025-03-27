### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
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
	using CairoMakie
	using Arya
end

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 2f62c5c2-e397-463b-9e73-f554c31a7b85
using PairPlots

# ╔═╡ 287cea0a-2745-4720-88b5-96d92742b382
md"""
This notebook checks the variance in assigned bin indicies due to uncertainty in structural parameters. Typically, this can vary ~ 
"""

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	@bind galaxy confirm(TextField(default="leo2"))
end

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath("..", galaxy, "figures"); FIGSUFFIX=".mcmc_hist"
end

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
bin_width = 0.05

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 109b91c2-c0ad-4683-849b-a56da5766ee4
import StatsBase: sem

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

# ╔═╡ 066c7b30-f818-4e4a-8db8-c8bac469f558
module MCMCUtils
	include("mcmc_utils.jl")
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
obs_props = MCMCUtils.get_obs_props(galaxy)

# ╔═╡ e06cd77a-cd95-40f1-997f-d2e90adcba67
stars = MCMCUtils.get_fits(galaxy, obs_props)

# ╔═╡ 08ab794f-02a9-4835-98ca-b0c01a19c60a
num_per_bin_default = max(round(Int, LilGuys.Interface.default_n_per_bin(stars.R_ell[stars.PSAT .> 0.2], nothing)), 2)

# ╔═╡ 11756aee-a4be-46d8-ad0e-3e2fc0089d2d
@bind num_per_bin NumberField(1:1000, default=num_per_bin_default)

# ╔═╡ 43d159c4-1599-4138-96b9-c9447a8d6315
md"""
The idea with this model is to use the likelihood

``
{\cal L} \propto {\rm Prior}(\theta) \prod_i (f(\theta, x_i) P_{\rm memb,\ \it i} + (1-f(\theta, x_i)) P_{\rm bg,\it\ i})
``

where $f$ is a model for the fraction of stars in a small region belonging to the satellite.
"""

# ╔═╡ e609b6e9-343e-487f-8b0b-fd79c7783e8c
stars.R_ell

# ╔═╡ 8e665feb-1a41-440c-a813-163cbf3be4f8
Nmemb = sum(stars.PSAT)

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ cd338225-85b9-45db-86cc-979d5b77a532
Arya.update_figsize!(3.25)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Analysis
"""

# ╔═╡ 028fb0dd-0076-4061-92cc-be08ec1b92d8
mc_params = MCMCUtils.StructuralParams(stars, obs_props, num_per_bin=num_per_bin)

# ╔═╡ 91867748-9b36-4f62-9310-8b778935776b
mc_params

# ╔═╡ 69327dc7-a82a-441e-901f-7474f4520575
bins = mc_params.bins

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
	
	stephist!(log10.(stars.R_ell[stars.PSAT .> 0.2]), bins=log10.(bins), label="PSAT > 0.2")
	stephist!(log10.(stars.R_ell), bins=log10.(bins), label="best")

	axislegend(position=:lt)
	@savefig "hist_bin_choices"
	fig
end

# ╔═╡ b515beab-8624-4ef1-ad0c-c6f548c74c0b
mc_params.bins

# ╔═╡ 82b3b170-729f-463e-9063-756a34840555
MCMCUtils.perturbed_radii(stars, mc_params)

# ╔═╡ 5d75c6ae-7e74-4ba0-8e4d-9e1a77e26748
function sample_bins(data=stars, params=mc_params, bins=bins)
	radii = MCMCUtils.perturbed_radii(data, mc_params)
	bin_idx = DE.bin_indices(radii, bins)
end

# ╔═╡ a06cc670-efb0-4ad3-b8fd-5b58ac6cb003
begin
	Nsamples = 1000
	bin_idxs = Matrix{Float64}(undef, length(stars.xi), Nsamples)
	for i in 1:Nsamples
		bin_idxs[:, i] = sample_bins()
	end
	bin_idxs
end

# ╔═╡ 79cf9232-1abe-4466-832b-67207fbc2e11
best_idx = round.(LilGuys.mean(eachcol(bin_idxs)))

# ╔═╡ 39f5e329-c32e-44bc-829e-f62a616cc88c
d_idx = bin_idxs .- best_idx

# ╔═╡ 9b303e04-62f2-4e72-a4fc-2019645764dd
scatter(log10.(stars.R_ell), LilGuys.std(eachcol(d_idx)), rasterize=true)

# ╔═╡ 9d12630d-7180-4d72-aa83-182b61d6fe79
scatter(log10.(stars.R_ell), best_idx, rasterize=true)

# ╔═╡ 0b859f82-3c1a-4be7-b2be-e742ce3f1f9c
hist(LilGuys.std(eachcol(bin_idxs)))

# ╔═╡ 067aaf2b-e6d2-4466-bd4d-d1b7f1aa17a6
let
	fig = Figure()
	ax = Axis(fig[1,1],
			  yscale=log10,
			  yticks = Makie.automatic,
			 )
		
	hist!(vec(d_idx))
	fig
end

# ╔═╡ 586c11c9-5933-4c5a-9e4f-3184e4a8c84f
let
	fig = Figure()
	ax = Axis(fig[1,1],
			  yscale=log10,
			  yticks = Makie.automatic,
			 )

	filt = stars.PSAT .> 0.2
	hist!(vec(d_idx[filt, :]))
	fig
end

# ╔═╡ 30ab89ef-6ce5-4a02-b2ae-9040b441c58d
md"""
Statistical summary of how much these variances matter
"""

# ╔═╡ cea39c5b-5e80-46e2-b41f-560a05a8bcac
LilGuys.mean(vec(d_idx[stars.PSAT .> 0.2, :]) .== 0)

# ╔═╡ 0cafe9c4-fb47-4acc-a514-686caa693a69
LilGuys.mean(vec(d_idx[stars.PSAT .> 0.2, :]) .< 2)

# ╔═╡ 34d1e570-3745-4508-9995-74750f434253
LilGuys.mean(vec(d_idx[stars.PSAT .> 0.2, :]) .< 4)

# ╔═╡ 176ffc99-bdc1-4f50-836f-5cf504610ce9
LilGuys.mean(vec(d_idx[stars.PSAT .> 0.2, :]) .< 10)

# ╔═╡ 4f2787f0-5f6b-4a31-a7b0-60c8ebe9ef8c
LilGuys.mean(vec(d_idx[stars.PSAT .> 0.2, :]) .< 20)

# ╔═╡ c162fbb4-21c6-4b91-833b-b7911fa83699
md"""
# Data writing
"""

# ╔═╡ a352dddb-c3a6-49ef-8238-43fe808c675d
gal_dir = joinpath("..", galaxy, "mcmc")

# ╔═╡ 3a556e21-bf2e-4d89-b97c-e09d428bfa9c
if !isdir(gal_dir)
	mkdir(gal_dir)
end

# ╔═╡ 45c7b5f1-b118-4ef6-afbd-2ee755713b75
outfile_bins = joinpath(gal_dir, "default_bins.toml")

# ╔═╡ be1c03e9-c152-4e5d-8106-38516f4dcf0f
open(outfile_bins, "w") do f
	d = LilGuys.struct_to_dict(mc_params)

	TOML.print(f, d)
	@info "wrote bins to $outfile_bins"

	d
end

# ╔═╡ Cell order:
# ╠═287cea0a-2745-4720-88b5-96d92742b382
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═08ab794f-02a9-4835-98ca-b0c01a19c60a
# ╠═11756aee-a4be-46d8-ad0e-3e2fc0089d2d
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═109b91c2-c0ad-4683-849b-a56da5766ee4
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═2f62c5c2-e397-463b-9e73-f554c31a7b85
# ╠═28550d54-ad88-4df5-ac04-f46a15588ff8
# ╠═066c7b30-f818-4e4a-8db8-c8bac469f558
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═e06cd77a-cd95-40f1-997f-d2e90adcba67
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╠═91867748-9b36-4f62-9310-8b778935776b
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╠═e609b6e9-343e-487f-8b0b-fd79c7783e8c
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═cd338225-85b9-45db-86cc-979d5b77a532
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═028fb0dd-0076-4061-92cc-be08ec1b92d8
# ╠═69327dc7-a82a-441e-901f-7474f4520575
# ╠═b515beab-8624-4ef1-ad0c-c6f548c74c0b
# ╠═82b3b170-729f-463e-9063-756a34840555
# ╠═5d75c6ae-7e74-4ba0-8e4d-9e1a77e26748
# ╠═a06cc670-efb0-4ad3-b8fd-5b58ac6cb003
# ╠═39f5e329-c32e-44bc-829e-f62a616cc88c
# ╠═9d12630d-7180-4d72-aa83-182b61d6fe79
# ╠═9b303e04-62f2-4e72-a4fc-2019645764dd
# ╠═79cf9232-1abe-4466-832b-67207fbc2e11
# ╠═0b859f82-3c1a-4be7-b2be-e742ce3f1f9c
# ╠═067aaf2b-e6d2-4466-bd4d-d1b7f1aa17a6
# ╠═586c11c9-5933-4c5a-9e4f-3184e4a8c84f
# ╠═30ab89ef-6ce5-4a02-b2ae-9040b441c58d
# ╠═cea39c5b-5e80-46e2-b41f-560a05a8bcac
# ╠═0cafe9c4-fb47-4acc-a514-686caa693a69
# ╠═34d1e570-3745-4508-9995-74750f434253
# ╠═176ffc99-bdc1-4f50-836f-5cf504610ce9
# ╠═4f2787f0-5f6b-4a31-a7b0-60c8ebe9ef8c
# ╠═c162fbb4-21c6-4b91-833b-b7911fa83699
# ╠═a352dddb-c3a6-49ef-8238-43fe808c675d
# ╠═3a556e21-bf2e-4d89-b97c-e09d428bfa9c
# ╠═45c7b5f1-b118-4ef6-afbd-2ee755713b75
# ╠═be1c03e9-c152-4e5d-8106-38516f4dcf0f
