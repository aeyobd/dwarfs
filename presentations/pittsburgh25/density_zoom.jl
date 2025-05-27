### A Pluto.jl notebook ###
# v0.20.8

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
include("./style.jl")

# ╔═╡ c638b2a8-faf0-4eee-b436-9ee22304e12a
import TOML

# ╔═╡ eda9cbe2-bafa-45e1-9494-202be90c6321
module Projection
	include(joinpath(ENV["HOME"], "LilGuys.jl/scripts/project_2d.jl"))
end

# ╔═╡ 10ef498c-5f47-4691-83f0-3b104ac29e00
module Animations
	include(joinpath(ENV["HOME"], "LilGuys.jl/scripts/animate_dm.jl"))
end

# ╔═╡ 683ae5c4-6f81-4150-bac0-4f8907d39fd3
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/1e7_V31_r3.2/orbit_smallperi")

# ╔═╡ f87386e4-7c1e-4125-a007-599d5099607f
idx_f = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))["idx_f"]

# ╔═╡ 57f6f815-8ee3-4263-9714-55eaccca645b
stars = LilGuys.read_hdf5_table(joinpath(modeldir, "../stars/exp2d_rs0.10/probabilities_stars.hdf5"))

# ╔═╡ e696f042-0c10-4b66-b169-9323a16d3c86
out = Output(modeldir, weights=stars.probability)

# ╔═╡ 0460ba75-559a-476a-b520-1a75d4971b8a
snap = out[idx_f]

# ╔═╡ fafd671d-2c4b-4087-945f-0913d2fb881b
x_vec = [sind(5), cosd(5), 0]

# ╔═╡ 47337350-e506-4727-a005-16ab355c4834
y_vec = [-sind(5), 0, cosd(5)]

# ╔═╡ a9cd55b6-b046-403d-aeeb-28e8a8179c9a
import LinearAlgebra: ⋅

# ╔═╡ 89cec6c2-67d3-4b68-932d-02154b0c9579
xycen = [x_vec, y_vec] .⋅ [snap.x_cen]

# ╔═╡ 5e9318e4-2537-4418-9308-9cc8f71f1ac1
plotrange = 10

# ╔═╡ 9fa746f9-e29f-4c99-ad96-7963838efc15
norm = 10^3.72

# ╔═╡ f307feed-b309-4b01-9b39-4cb7f7e3f993
fg_color = :grey

# ╔═╡ b2dbb5b2-2ef1-4c5a-9313-cce8c5d03bbc
let
	x, y, w = Projection.get_xy(out, idx_f, x_vec=x_vec, y_vec=y_vec)

	N = 500
	bins = (LinRange(xycen[1] - plotrange, xycen[1] + plotrange, N),
			LinRange(xycen[2] - plotrange, xycen[2] + plotrange, N)
		   )
	hist = Projection.project_points(x, y, w, bins)


	fig = Figure(backgroundcolor=:transparent)
	ax = Axis(fig[1,1], topspinecolor=fg_color, bottomspinecolor=fg_color, leftspinecolor=fg_color, rightspinecolor=fg_color)

	hidedecorations!()
	cmax = maximum(log10.(hist))
	p = image!(extrema.(bins)..., log10.(hist), colorrange=(cmax-5, cmax))

	#Colorbar(fig[1,2], p, width=theme(:fontsize))


	@info plotrange * 2 * 0.25
	lines!([0.1, 0.35], [0.1, 0.1], space=:relative, color=fg_color, linewidth=theme(:linewidth)[] )
	text!(0.1, 0.1, text="5 kpc", space=:relative, color=fg_color, font="hero", fontsize=theme(:fontsize)[])


	colsize!(fig.layout, 1, Aspect(1, 1.))

	fig

end

# ╔═╡ d742c9b3-cc81-4f4c-84b1-5140abbae162
let
	x, y, w = Projection.get_xy(out, idx_f, x_vec=x_vec, y_vec=y_vec)

	N = 500
	bins = (LinRange(xycen[1] - plotrange, xycen[1] + plotrange, N),
			LinRange(xycen[2] - plotrange, xycen[2] + plotrange, N)
		   )
	hist = Projection.project_points(x, y, ones(size(x)), bins)

	fig = Figure(backgroundcolor=:transparent)
	ax = Axis(fig[1,1], topspinecolor=fg_color, bottomspinecolor=fg_color, leftspinecolor=fg_color, rightspinecolor=fg_color)
	hidedecorations!()

	cmax = maximum(log10.(hist))

	p = image!(extrema.(bins)..., log10.(hist), colormap=:devon, colorrange=(cmax-5, cmax))

	#Colorbar(fig[1,2], p)

	colsize!(fig.layout, 1, Aspect(1, 1.))
	fig
end

# ╔═╡ 1c58a4de-106d-4ee6-a222-a786fcd1063f
md"""
This plot just plots the region for the density zoom ins
"""

# ╔═╡ b943187d-8fa3-472a-8de2-80312e487207
let
	fig = Figure(backgroundcolor=:transparent, figure_padding=0)
	ax = Axis(fig[1,1], #topspinecolor=:white, bottomspinecolor=:white, leftspinecolor=:white, rightspinecolor=:white,
			  limits=(-125,125,-125,125), aspect=DataAspect(), backgroundcolor=:transparent)
	hidedecorations!()

	x1 = xycen[1] - plotrange
	x2 = xycen[1] + plotrange
	y1, y2 = xycen[2] - plotrange, xycen[2] + plotrange

	lines!([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], color=fg_color, linewidth=theme(:linewidth)[]/2)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═c638b2a8-faf0-4eee-b436-9ee22304e12a
# ╠═eda9cbe2-bafa-45e1-9494-202be90c6321
# ╠═10ef498c-5f47-4691-83f0-3b104ac29e00
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═683ae5c4-6f81-4150-bac0-4f8907d39fd3
# ╠═e696f042-0c10-4b66-b169-9323a16d3c86
# ╠═f87386e4-7c1e-4125-a007-599d5099607f
# ╠═0460ba75-559a-476a-b520-1a75d4971b8a
# ╠═57f6f815-8ee3-4263-9714-55eaccca645b
# ╠═fafd671d-2c4b-4087-945f-0913d2fb881b
# ╠═47337350-e506-4727-a005-16ab355c4834
# ╠═a9cd55b6-b046-403d-aeeb-28e8a8179c9a
# ╠═89cec6c2-67d3-4b68-932d-02154b0c9579
# ╠═5e9318e4-2537-4418-9308-9cc8f71f1ac1
# ╠═9fa746f9-e29f-4c99-ad96-7963838efc15
# ╠═f307feed-b309-4b01-9b39-4cb7f7e3f993
# ╠═b2dbb5b2-2ef1-4c5a-9313-cce8c5d03bbc
# ╠═d742c9b3-cc81-4f4c-84b1-5140abbae162
# ╠═1c58a4de-106d-4ee6-a222-a786fcd1063f
# ╠═b943187d-8fa3-472a-8de2-80312e487207
