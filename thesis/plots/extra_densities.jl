### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 9309c10c-6ba3-436c-b975-36d26dafb821
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 1f1780d7-c6ed-4a92-b263-f8e43e9fab68
md"""
Monolithic figure creation haha
"""

# ╔═╡ 2bacd818-4985-4922-85a3-716bdfda5146
import DensityEstimators: histogram2d

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ e93918f6-0c72-4dec-9cdf-97f185c0bceb
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ c2eed577-ff2c-4a28-b043-3ba7f4a68d3a
function plot_density_f!(galaxyname, modelname, starsname; norm_shift=0, kwargs...)
	prof_i, prof_f, norm = ModelUtils.load_stellar_profiles(galaxyname, modelname, starsname, norm_shift=norm_shift)

	lines!(prof_f.log_R, prof_f.log_Sigma; kwargs...)

end

# ╔═╡ bfec5b65-905a-4a4a-858c-67c3091e2aaa
let
	fig = Figure()

	ax = Axis(fig[1,1])

	plot_density_f!("sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10")

	plot_density_f!("sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10")


	plot_density_f!("sculptor", "1e6_new_v31_r3.2/L3M11_9Gyr_smallperi.a4", "exp2d_rs0.10", label="MW impact")

	plot_density_f!("sculptor", "1e6_v43_r5_beta0.2_a4/orbit_smallperi", "exp2d_rs0.10", label="anisotropy")
	plot_density_f!("sculptor", "1e6_v48_r7_oblate_0.5/orbit_smallperi", "exp2d_rs0.10", label="oblate")
	fig

end

# ╔═╡ Cell order:
# ╠═1f1780d7-c6ed-4a92-b263-f8e43e9fab68
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═2bacd818-4985-4922-85a3-716bdfda5146
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═9309c10c-6ba3-436c-b975-36d26dafb821
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e93918f6-0c72-4dec-9cdf-97f185c0bceb
# ╠═c2eed577-ff2c-4a28-b043-3ba7f4a68d3a
# ╠═bfec5b65-905a-4a4a-858c-67c3091e2aaa
