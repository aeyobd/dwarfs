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

	FIGDIR = joinpath(galaxy, "figures"); FIGSUFFIX=".mcmc_hist"
end

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
bin_width = 0.05

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

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

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
module GaiaFilters
	include("../utils/gaia_filters.jl")
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
all_stars = GaiaFilters.read_gaia_stars(
	GaiaFilters.GaiaFilterParams(obs_props, filename = datafile,)
)

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

# ╔═╡ 1a5e8de2-0d74-48a9-aead-0855602734f3
# ╠═╡ disabled = true
#=╠═╡
bins = 10 .^ LilGuys.Interface.bins_both(log10.(best_stars.r_ell), nothing, bin_width=bin_width_min, num_per_bin=N_per_bin_min)

  ╠═╡ =#

# ╔═╡ 41afbd13-c22f-46f3-8d1f-43285ccf867e
bins = 10 .^ LilGuys.bins_equal_width(log10.(best_stars.R_ell), nothing, bin_width=bin_width)

# ╔═╡ 8e665feb-1a41-440c-a813-163cbf3be4f8
Nmemb = sum(best_stars.PSAT)

# ╔═╡ c8ca665d-0c63-442f-8b67-8ec175d38b06
hist_example = DE.histogram(best_stars.R_ell, bins)

# ╔═╡ fda2c541-04ef-44bc-a0d0-30854cfe4d71


# ╔═╡ 9eb25479-df26-4a86-bb6b-9b14f404ee93
sum(best_stars.PSAT .> 0.2)

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
	
	stephist!(log10.(best_stars.R_ell[best_stars.PSAT .> 0.2]), bins=log10.(bins), label="PSAT > 0.2")
	stephist!(log10.(best_stars.R_ell), bins=log10.(bins), label="best")

	axislegend(position=:lt)
	@savefig "hist_bin_choices"
	fig
end

# ╔═╡ dac9bbd7-8fac-4fab-8e74-7e1fea166616
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,

	)

	prof = LilGuys.StellarDensityProfile(best_stars.R_ell, bins=log10.(bins))
	prof.log_Sigma[.!isfinite.(prof.log_Sigma)] .= minimum(prof.log_Sigma[isfinite.(prof.log_Sigma)]) .- 1
	scatter!(prof.log_R, prof.log_Sigma)
	vlines!(log10(obs_props["r_h"]))
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
#=╠═╡
filt_max = bins[2:end] .< R_max
  ╠═╡ =#

# ╔═╡ ebdfb5dc-b996-4ead-a61f-ebb9f93adf0e
open(binsout, "w") do f
	TOML.print(f, OrderedDict(
		"binwidth" => bin_width,
		"bins" => bins,
	)
	)
end

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═28550d54-ad88-4df5-ac04-f46a15588ff8
# ╠═57a19d65-59d9-46cf-8916-d9ac3a4dc92b
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╠═f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
# ╠═91867748-9b36-4f62-9310-8b778935776b
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═e7bf7a13-afca-49b8-8771-f7914adb347b
# ╠═1a5e8de2-0d74-48a9-aead-0855602734f3
# ╠═41afbd13-c22f-46f3-8d1f-43285ccf867e
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═c8ca665d-0c63-442f-8b67-8ec175d38b06
# ╠═fda2c541-04ef-44bc-a0d0-30854cfe4d71
# ╠═9eb25479-df26-4a86-bb6b-9b14f404ee93
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═dac9bbd7-8fac-4fab-8e74-7e1fea166616
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═639c7d31-8693-45dd-88be-492b124804e9
# ╠═b1487f12-76b2-4a2f-8d50-4e80f90b1d9d
# ╠═bf8d5936-1a4e-47c2-bb22-531ab344b8ad
# ╠═ebdfb5dc-b996-4ead-a61f-ebb9f93adf0e
