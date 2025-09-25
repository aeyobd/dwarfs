### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 3d31db99-0d07-42fc-a437-215e29913021
using OrderedCollections: OrderedDict

# ╔═╡ 287cea0a-2745-4720-88b5-96d92742b382
md"""
This notebook checks the variance in assigned bin indicies due to uncertainty in structural parameters. Typically, this can vary ~ 
"""

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	md"""
	$(@bind galaxy confirm(TextField(default="galaxy")))

	$(@bind bin_width confirm(NumberField(0:0.1:1, 0.05)))

	"""
end

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath("..", galaxy, "figures"); FIGSUFFIX=".mcmc_hist"
end

# ╔═╡ 11756aee-a4be-46d8-ad0e-3e2fc0089d2d
num_per_bin = 20

# ╔═╡ e53cdf92-f17d-40ae-9785-c00b7e33aeba
PSAT_min = 0.2

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

# ╔═╡ a955c503-c834-4719-89d0-a0004c33ea15
@info sum(stars.PSAT .> 0.2)

# ╔═╡ f9ea040c-0d80-4ccc-a990-010d72a5647e
sum(stars.PSAT)

# ╔═╡ 5c58e059-7e5c-41b6-accd-7020beda42f8
f_cen = LilGuys.mean(stars.PSAT[stars.R_ell .< obs_props["r_h"]/3])

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
mc_params = MCMCUtils.StructuralParams(stars, obs_props, num_per_bin=num_per_bin, bin_width=bin_width)

# ╔═╡ da4218ef-34ed-4464-81ad-fb8fa216fb9a
@info length(mc_params.bins)

# ╔═╡ ff4e7ea9-fbd9-4e07-87ca-934ff8961393
mc_params.bins

# ╔═╡ a6060154-db75-456c-903b-959c508dd7c0
@assert minimum(diff(log10.(mc_params.bins))) >= bin_width * (1-1e-10)

# ╔═╡ f6a14386-2da3-4dc7-9fbe-ccd632b7ae53
diff(log10.(mc_params.bins))

# ╔═╡ 99e00997-c816-4cda-9643-378a5eff7972
counts_in_bin = LilGuys.histogram(stars.R_ell, mc_params.bins)[2]

# ╔═╡ 9ea81473-2bde-43e2-b092-d2ce1b7ff5e6
@assert minimum(counts_in_bin) >= num_per_bin

# ╔═╡ 81876bbe-a200-49ce-9155-5e3c98bc1db5
MCMCUtils.centring_stats(stars, PSAT_min=0.2)

# ╔═╡ dbd3a8ce-bd01-492c-a741-00f310091e7e
MCMCUtils.centring_stats(stars, PSAT_min=0.8)

# ╔═╡ e7af37c8-801a-41ba-8f67-50011a144b5d
MCMCUtils.centring_stats(stars, PSAT_min=0.99)

# ╔═╡ 6e10ff2d-6621-4f82-8bdf-42fa606bbb8f
MCMCUtils.centring_stats(stars, PSAT_min=0.999)

# ╔═╡ a8289faa-8f9a-463d-b0ea-f20cd900cee0
MCMCUtils.centring_stats(stars, PSAT_min=0.9999)

# ╔═╡ 3cc1745a-8d6a-4560-8815-26f049538c9e
members = stars[stars.PSAT .> PSAT_min, :]

# ╔═╡ 374439b8-767c-4872-991a-90298c9a6dd6
N_min_clean = max(round(Int, LilGuys.Interface.default_n_per_bin(members.R_ell, nothing)), 2)

# ╔═╡ 08ab794f-02a9-4835-98ca-b0c01a19c60a
num_per_bin_default = max(round(Int, LilGuys.Interface.default_n_per_bin(members.R_ell, nothing)), 2) / f_cen

# ╔═╡ 62d39b10-f078-4cd3-b4c6-0945386ac507
num_per_bin_default

# ╔═╡ 33e02945-2ab5-4f55-b3d5-e8855d86f120
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		yscale = log10,
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = log_r_label,
		ylabel = "counts",
		limits=(nothing, nothing, 1, 1e6),
		xlabelsize=12
	)
	
	stephist!(log10.(members.R_ell), bins=log10.(mc_params.bins), label="PSAT > 0.2")
	stephist!(log10.(stars.R_ell), bins=log10.(mc_params.bins), label="best")
	scatter!(log10.(mc_params.bins), ones(length(mc_params.bins)))
	axislegend(position=:lt)
	@savefig "hist_bin_choices"
	fig
end

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

# ╔═╡ 360f23d8-934f-4a40-9e3a-a0aa75e4d495
outfile_params = joinpath(gal_dir, "default_bins_params.toml")

# ╔═╡ be1c03e9-c152-4e5d-8106-38516f4dcf0f
open(outfile_bins, "w") do f
	d = LilGuys.struct_to_dict(mc_params)

	TOML.print(f, d)
	@info "wrote bins to $outfile_bins"

	d
end

# ╔═╡ 03bc0ea8-7f1e-4f1e-988e-45fbf3daa52c
open(outfile_params, "w") do f
	d = OrderedDict(
		"f_sat_cen" => f_cen,
		"PSAT" => PSAT_min,
		"num_per_bin" => num_per_bin,
		"bin_width" => bin_width,
		"centring_error" => LilGuys.struct_to_dict(
			MCMCUtils.centring_stats(stars, PSAT_min=PSAT_min))
	)

	TOML.print(f, d)
	@info "wrote info to $outfile_params"

	d
end

# ╔═╡ e2162d45-66b7-42d0-acf5-e3c0fa6b5db2
LilGuys.struct_to_dict(
			MCMCUtils.centring_stats(stars, PSAT_min=PSAT_min))

# ╔═╡ Cell order:
# ╠═287cea0a-2745-4720-88b5-96d92742b382
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═11756aee-a4be-46d8-ad0e-3e2fc0089d2d
# ╠═62d39b10-f078-4cd3-b4c6-0945386ac507
# ╠═e53cdf92-f17d-40ae-9785-c00b7e33aeba
# ╠═a955c503-c834-4719-89d0-a0004c33ea15
# ╠═f9ea040c-0d80-4ccc-a990-010d72a5647e
# ╠═da4218ef-34ed-4464-81ad-fb8fa216fb9a
# ╠═ff4e7ea9-fbd9-4e07-87ca-934ff8961393
# ╠═5c58e059-7e5c-41b6-accd-7020beda42f8
# ╠═374439b8-767c-4872-991a-90298c9a6dd6
# ╠═08ab794f-02a9-4835-98ca-b0c01a19c60a
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
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╠═e609b6e9-343e-487f-8b0b-fd79c7783e8c
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═cd338225-85b9-45db-86cc-979d5b77a532
# ╠═a6060154-db75-456c-903b-959c508dd7c0
# ╠═f6a14386-2da3-4dc7-9fbe-ccd632b7ae53
# ╠═9ea81473-2bde-43e2-b092-d2ce1b7ff5e6
# ╠═99e00997-c816-4cda-9643-378a5eff7972
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═028fb0dd-0076-4061-92cc-be08ec1b92d8
# ╠═81876bbe-a200-49ce-9155-5e3c98bc1db5
# ╠═dbd3a8ce-bd01-492c-a741-00f310091e7e
# ╠═e7af37c8-801a-41ba-8f67-50011a144b5d
# ╠═6e10ff2d-6621-4f82-8bdf-42fa606bbb8f
# ╠═a8289faa-8f9a-463d-b0ea-f20cd900cee0
# ╠═3cc1745a-8d6a-4560-8815-26f049538c9e
# ╠═c162fbb4-21c6-4b91-833b-b7911fa83699
# ╠═a352dddb-c3a6-49ef-8238-43fe808c675d
# ╠═3a556e21-bf2e-4d89-b97c-e09d428bfa9c
# ╠═45c7b5f1-b118-4ef6-afbd-2ee755713b75
# ╠═360f23d8-934f-4a40-9e3a-a0aa75e4d495
# ╠═be1c03e9-c152-4e5d-8106-38516f4dcf0f
# ╠═3d31db99-0d07-42fc-a437-215e29913021
# ╠═03bc0ea8-7f1e-4f1e-988e-45fbf3daa52c
# ╠═e2162d45-66b7-42d0-acf5-e3c0fa6b5db2
