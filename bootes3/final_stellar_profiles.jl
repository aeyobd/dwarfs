### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 5183c43a-3db0-11f1-8a22-fff0545dbb1d
begin
	import Pkg; Pkg.activate()

	using LilGuys
	import TOML

	using CairoMakie, Arya
end

# ╔═╡ 9fcbcfa0-2a04-41c8-9aa1-4033d294b63f
function load_i_f_stars(modelname)
	return LilGuys.SurfaceDensityProfile(modelname)

# ╔═╡ Cell order:
# ╠═5183c43a-3db0-11f1-8a22-fff0545dbb1d
# ╠═9fcbcfa0-2a04-41c8-9aa1-4033d294b63f
