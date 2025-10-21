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
- `scl_smallperi_i_f`
- `scl_plummer_i_f`
- `scl_lmc_i_f`
- `scl_lmc_plummer_i_f`
- `umi_smallperi_i_f`
- `umi_plummer_i_f`

Appendix
- `scl_impact_i_f`
"""

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ e93918f6-0c72-4dec-9cdf-97f185c0bceb
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 7948bbd7-2dda-46c8-a850-d80944ba9096
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ c2eed577-ff2c-4a28-b043-3ba7f4a68d3a
compare_both = ModelUtils.compare_both

# ╔═╡ 3317cfec-be70-418c-944a-0a6741f03cf3
@savefig "scl_smallperi_i_f" compare_both(modelnames["scl_smallperi"]...,  r_j=true, norm_shift=0.2, 	break_height=-3, title="Sculptor: smallperi-exponential")

# ╔═╡ fc45829c-bf28-4d7b-9bea-e223e4b29e7d
@savefig "scl_plummer_i_f" compare_both(modelnames["scl_smallperi_plummer"]..., norm_shift=0,  r_j=true, break_height=-2, title="Sculptor: smallperi-Plummer")

# ╔═╡ 5f989e91-8870-462e-a363-767be64dac5c
@savefig "scl_lmc_i_f" compare_both(modelnames["scl_lmc"]..., norm_shift=0.1, lmc=true, r_j=true, break_height=-7.590274001669055, title="Sculptor: LMC–exponential")

# ╔═╡ e9a0ce9f-27a5-47dc-8c4d-22887f6a7fc1
@savefig "scl_lmc_plummer_i_f" compare_both(modelnames["scl_lmc_plummer"]..., norm_shift=0.1, lmc=true, r_j=true, break_height=-7.590274001669055, title="Sculptor: LMC–Plummer")

# ╔═╡ 77c057e4-78da-4201-b9f8-cb3135508a93
md"""
# Ursa Minor
"""

# ╔═╡ 3263677b-e0e3-4d30-8e28-309de99a8603
@savefig "umi_smallperi_i_f" compare_both(modelnames["umi_smallperi"]...,  r_j=true, norm_shift=0.0, 	break_height=-5.5, title="Ursa Minor: smallperi-exponential")

# ╔═╡ d0ce116f-aa51-43d8-9748-806eb44b8409
@savefig "umi_plummer_i_f" compare_both(modelnames["umi_smallperi_plummer"]...,  r_j=true, norm_shift=0.0,break_height=-3, title="Ursa Minor: smallperi-Plummer")

# ╔═╡ 8979aeee-5c1d-45bb-ac5b-0d0655cae3a8
md"""
## Ursa Minor extra
"""

# ╔═╡ 18f2e3da-b749-4134-9604-4e2318ad82ee
ModelUtils.compare_both("ursa_minor", "1e6_v37_r5.0/orbit_mean.2", "exp2d_rs0.10",  r_j=true, norm_shift=-0.6, 	break_height=-4)

# ╔═╡ f11acd53-a04d-4e1a-b7f5-d66b6011e3a2
@savefig "scl_impact_i_f" compare_both(modelnames["mw_impact"]...,  title="Sculptor: MW impact", break_height=-4, r_j=true)

# ╔═╡ d8c7cb86-1712-4547-8552-3bf485c86cca
md"""
## Sculptor tests
"""

# ╔═╡ c7b50c9b-1bfd-4abf-89cd-7c8cd1dc45c0
compare_both(modelnames["oblate"]...,  title="Sculptor: oblate", norm_shift=-0.5)

# ╔═╡ 2a3b1df5-7a44-420c-b303-856346df40cf
compare_both(modelnames["anisotropic"]...,  title="Sculptor: anisotropic")

# ╔═╡ 91d8ebb5-47b8-44ed-a334-96dbc6514c44
ModelUtils.compare_both("sculptor", "1e6_new_v43_r7/orbit_smallperi.3", "exp2d_rs0.10", norm_shift=-0.7, 	break_height=-4, title="Sculptor heavy")

# ╔═╡ fafd07cb-2378-421b-b6e5-8525c71f9c8d
compare_both("sculptor", "1e6_new_v43_r7/orbit_smallperi.3", "exp2d_rs0.10", norm_shift=-0.7, 	break_height=-4, title="Sculptor heavy")

# ╔═╡ e231dce2-1969-4e86-a519-e9e2eb91d9f6
ModelUtils.compare_both("sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10", norm_shift=0.0, 	break_height=-4, title="Sculptor low-resolution")

# ╔═╡ 8e4c5116-7de6-4f8d-bc09-ff3aa4ffc50f
ModelUtils.compare_both("sculptor", "1e4_exp_M3e-4_r0.1/orbit_smallperi", "stars", title="Sculptor: dm free", norm_shift=0)

# ╔═╡ Cell order:
# ╠═1f1780d7-c6ed-4a92-b263-f8e43e9fab68
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e93918f6-0c72-4dec-9cdf-97f185c0bceb
# ╠═7948bbd7-2dda-46c8-a850-d80944ba9096
# ╠═c2eed577-ff2c-4a28-b043-3ba7f4a68d3a
# ╠═3317cfec-be70-418c-944a-0a6741f03cf3
# ╠═fc45829c-bf28-4d7b-9bea-e223e4b29e7d
# ╠═5f989e91-8870-462e-a363-767be64dac5c
# ╠═e9a0ce9f-27a5-47dc-8c4d-22887f6a7fc1
# ╟─77c057e4-78da-4201-b9f8-cb3135508a93
# ╠═3263677b-e0e3-4d30-8e28-309de99a8603
# ╠═d0ce116f-aa51-43d8-9748-806eb44b8409
# ╟─8979aeee-5c1d-45bb-ac5b-0d0655cae3a8
# ╠═18f2e3da-b749-4134-9604-4e2318ad82ee
# ╠═f11acd53-a04d-4e1a-b7f5-d66b6011e3a2
# ╟─d8c7cb86-1712-4547-8552-3bf485c86cca
# ╠═c7b50c9b-1bfd-4abf-89cd-7c8cd1dc45c0
# ╠═2a3b1df5-7a44-420c-b303-856346df40cf
# ╠═91d8ebb5-47b8-44ed-a334-96dbc6514c44
# ╠═fafd07cb-2378-421b-b6e5-8525c71f9c8d
# ╠═e231dce2-1969-4e86-a519-e9e2eb91d9f6
# ╠═8e4c5116-7de6-4f8d-bc09-ff3aa4ffc50f
