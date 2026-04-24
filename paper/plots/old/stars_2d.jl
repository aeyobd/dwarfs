### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 1f1780d7-c6ed-4a92-b263-f8e43e9fab68
md"""
Monolithic figure creation haha

"""

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ e93918f6-0c72-4dec-9cdf-97f185c0bceb
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 56365cfd-16e0-4567-a85c-f59372a611c7
function plot_stars(modelname)
	fig = Figure()
	Σ0 = - ModelUtils.get_normalization(ModelUtils.load_expected_density_profile(modelname[1]))

	_, _, norm = ModelUtils.load_stellar_profiles(modelname...)

	colorrange = [Σ0-9, Σ0+1]

	r_b = ModelUtils.get_r_b(modelname...)
	
	p = ModelUtils.plot_stars_2d(fig[1,1], modelname..., colorrange=colorrange, norm=norm, r_b=r_b)

	Colorbar(fig[1, 2], p, label="log stellar density")

	rowsize!(fig.layout, 1, Aspect(1, 1))
	resize_to_layout!()
	fig
end

# ╔═╡ 7948bbd7-2dda-46c8-a850-d80944ba9096
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ dbed5ec4-6030-4667-a87b-e8100deebe4a
@savefig "umi_stars_2d" plot_stars(modelnames["umi_smallperi_plummer"],)

# ╔═╡ 7132c6c8-8b28-4d3a-9cac-b0bfcfa09cec
@savefig "scl_stars_2d" plot_stars(modelnames["scl_lmc_plummer"],)

# ╔═╡ Cell order:
# ╠═1f1780d7-c6ed-4a92-b263-f8e43e9fab68
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e93918f6-0c72-4dec-9cdf-97f185c0bceb
# ╠═56365cfd-16e0-4567-a85c-f59372a611c7
# ╠═7948bbd7-2dda-46c8-a850-d80944ba9096
# ╠═dbed5ec4-6030-4667-a87b-e8100deebe4a
# ╠═7132c6c8-8b28-4d3a-9cac-b0bfcfa09cec
