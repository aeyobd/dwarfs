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

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	@bind galaxy confirm(TextField(default="galaxy"))
end

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath(galaxy, "figures"); FIGSUFFIX=".bins_test"
end

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
@bind bin_width confirm(NumberField(0.01:0.01:0.5, default=0.05))

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators as DE

# ╔═╡ 57a19d65-59d9-46cf-8916-d9ac3a4dc92b
datafile = let
	dir = joinpath(galaxy, "data")

	filenames = ["jensen+24_1c.fits", "jensen+24_2c.fits", "j24_1c.fits", "j24_2c.fits"]
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

# ╔═╡ a823c38b-e843-4896-bb3f-448d2400356f
filt_params = GaiaFilters.GaiaFilterParams(obs_props, filename = datafile,)

# ╔═╡ f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
all_stars = GaiaFilters.read_gaia_stars(
	filt_params
)

# ╔═╡ 2e3681c1-015c-40d5-a933-b5bc4f6b3f85
best_stars = GaiaFilters.select_members(all_stars, filt_params)

# ╔═╡ 9ca9652b-64e8-4105-a9f6-6b499a31e6b1
R_ell = maximum(best_stars.R_ell)

# ╔═╡ d897930c-f087-4be2-b9a6-a89d364b3cd6
@bind bin_max confirm(NumberField(0:0.01:300, default=R_ell))

# ╔═╡ 7c46769e-0fd8-4001-a00d-1aed26552801
members = filter(r->r.PSAT > 0.2, best_stars)

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

# ╔═╡ 41afbd13-c22f-46f3-8d1f-43285ccf867e
bins = reverse(10 .^ (log10(bin_max):-bin_width:log10(minimum(best_stars.R_ell))))

# ╔═╡ 4446b9c4-e309-43df-a06d-ad9effc03a0c
@assert diff(log10.(bins)) ≈ fill(0.05, length(bins)-1)

# ╔═╡ 86613400-9276-4d94-9230-da996c4373de
@assert bins[end] ≈ bin_max

# ╔═╡ 8e665feb-1a41-440c-a813-163cbf3be4f8
Nmemb = sum(best_stars.PSAT)

# ╔═╡ c8ca665d-0c63-442f-8b67-8ec175d38b06
hist_example = DE.histogram(best_stars.R_ell, bins)

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
	
	stephist!(log10.(members.R_ell), bins=log10.(bins), label="PSAT > 0.2")
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

	prof = LilGuys.SurfaceDensityProfile(members.R_ell, bins=log10.(bins))
	prof.log_Sigma[.!isfinite.(prof.log_Sigma)] .= minimum(prof.log_Sigma[isfinite.(prof.log_Sigma)]) .- 1
	LilGuys.plot_log_Σ!(ax, prof)

	prof = LilGuys.SurfaceDensityProfile(best_stars.R_ell, bins=log10.(bins))
	prof.log_Sigma[.!isfinite.(prof.log_Sigma)] .= minimum(prof.log_Sigma[isfinite.(prof.log_Sigma)]) .- 1
	LilGuys.plot_log_Σ!(ax, prof)
	vlines!(log10(obs_props["r_h"]))
	fig
end

# ╔═╡ 8a23338b-f994-4710-b27c-2d12f608f3c1


# ╔═╡ 21ad2129-b2e2-459c-ab44-629ae16b027c
log10.(bins)

# ╔═╡ fef42ea7-6465-4493-956c-7848eaba9200
minimum(log10.(best_stars.R_ell))

# ╔═╡ 52a087d1-2a47-4b58-b131-5c09061cae76
log10(R_ell)

# ╔═╡ b2d864b4-e1a9-499d-9f44-4e93614d41f4
print(log10.(filter(x->x<=R_ell, bins)))

# ╔═╡ 5056894a-7582-4c1c-b881-888d46e1e2c3
md"""
```
inherits = "bins_eqw.toml"
```
"""

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═d897930c-f087-4be2-b9a6-a89d364b3cd6
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
# ╠═a823c38b-e843-4896-bb3f-448d2400356f
# ╠═f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
# ╠═2e3681c1-015c-40d5-a933-b5bc4f6b3f85
# ╠═9ca9652b-64e8-4105-a9f6-6b499a31e6b1
# ╠═7c46769e-0fd8-4001-a00d-1aed26552801
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═41afbd13-c22f-46f3-8d1f-43285ccf867e
# ╠═4446b9c4-e309-43df-a06d-ad9effc03a0c
# ╠═86613400-9276-4d94-9230-da996c4373de
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═c8ca665d-0c63-442f-8b67-8ec175d38b06
# ╠═9eb25479-df26-4a86-bb6b-9b14f404ee93
# ╠═33e02945-2ab5-4f55-b3d5-e8855d86f120
# ╠═dac9bbd7-8fac-4fab-8e74-7e1fea166616
# ╠═8a23338b-f994-4710-b27c-2d12f608f3c1
# ╠═21ad2129-b2e2-459c-ab44-629ae16b027c
# ╠═fef42ea7-6465-4493-956c-7848eaba9200
# ╠═52a087d1-2a47-4b58-b131-5c09061cae76
# ╠═b2d864b4-e1a9-499d-9f44-4e93614d41f4
# ╠═5056894a-7582-4c1c-b881-888d46e1e2c3
