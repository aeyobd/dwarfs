### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ b0a0dddc-fb5a-4bb2-b048-a54859d0b703
using DataFrames, CSV

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 913e0316-a04b-4270-ba31-0ba0f7fdd705
galaxyname = "ursa_minor"

# ╔═╡ 34c0f4f8-ef34-4b58-bcdd-7de69b58db2d
import DensityEstimators

# ╔═╡ bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
import TOML

# ╔═╡ 280e0cef-896a-4c14-a1ac-a611418be57e
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
xy_iso = CSV.read("resources/EP2020_iso_xy.csv", DataFrame, tasks=1)

# ╔═╡ 53cdcc20-96d4-4bd5-8028-df9a38af71ae
xz_iso = CSV.read("resources/EP2020_iso_xz.csv", DataFrame, tasks=1)

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b9600d93-5946-4380-a2ae-9b5f673bbaf5
modelname = if galaxyname == "sculptor"
	joinpath(modelnames["scl_smallperi"][1:2]...)
elseif galaxyname == "ursa_minor"
	joinpath(modelnames["umi_smallperi"][1:2]...)
end

# ╔═╡ 0f71807d-d698-4164-9f30-49af8dd8ba55
point_orbit = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "simulations", modelname, "orbit.toml"))

# ╔═╡ b593188f-1618-44d4-b38c-d5be448e927c
pos_final = LilGuys.position(LilGuys.transform(Galactocentric, ICRS(point_orbit)))

# ╔═╡ 4adea59d-a467-4b64-a7f3-e1444c2b35c0
starsname = if galaxyname == "sculptor"
	modelnames["scl_smallperi"][3]
elseif galaxyname == "ursa_minor"
	modelnames["umi_smallperi"][3]
end

# ╔═╡ 32db23d9-7959-41ac-aff4-b63df5e4b94a
figname = Dict(
	"sculptor" => "scl",
	"ursa_minor" => "umi"
)[galaxyname]

# ╔═╡ 21d0e542-08ea-4f56-a07d-f1a8e7e69019
figtitle = Dict(
	"sculptor" => "Sculptor",
	"ursa_minor" => "Ursa Minor"
)[galaxyname]

# ╔═╡ df61eda6-6d70-4ce4-8f6b-6936dc12283d
stars = LilGuys.read_hdf5_table(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "../stars", starsname, "probabilities_stars.hdf5"))

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
out = Output(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname), weights=stars.probability);

# ╔═╡ a3be2d61-98eb-4037-afb4-4155ba24cc21
orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "orbital_properties.toml"))

# ╔═╡ d50c120c-d72a-4fdc-b18a-9ec185c25c04
orbit = Orbit(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "centres.hdf5"))

# ╔═╡ b2cdaa3a-1c20-4609-b6f7-5b0a5f572102
import StatsBase

# ╔═╡ 5a40b893-021b-46e5-a115-0284e13ae7ae
bins = LinRange(-150, 150, 512)

# ╔═╡ 04c8c564-f768-4845-ac22-898173a10582
function get_histogram(snap, bins=bins; weights=nothing)
	x = snap.positions[2, :]
	y = snap.positions[3, :]

	if eltype(bins) <: Real
		bins = (bins, bins)
	end
	
	if weights === nothing                                                              
        h1 = StatsBase.fit(StatsBase.Histogram, (x, y), bins)                                      
    else                                                                          
        h1 = StatsBase.fit(StatsBase.Histogram, (x, y), StatsBase.weights(weights), bins) 
    end         
	
	return bins, StatsBase.normalize(h1, mode=:density).weights
end

# ╔═╡ 4ae004da-7e17-4f7c-a11d-6ad97dcbe6dd
orbit_props["idx_apos"], orbit_props["idx_peris"]

# ╔═╡ d57c2310-706f-4a0f-af6e-f4dae0cd8376
orbit_props["idx_f"]

# ╔═╡ 962d5024-6afa-452c-bb04-85046971ee2a
idxs = if galaxyname == "sculptor"
	[orbit_props["idx_f"]]
elseif galaxyname == "ursa_minor"
	[orbit_props["idx_f"]]
end

# ╔═╡ 971bb9b8-6807-40f7-ad3e-c9e62ab4f3f1
function plot_today!()
	scatter!(pos_final[2], pos_final[3], color=COLORS[3], markersize=theme(:markersize)[]/2)
end

# ╔═╡ 800caee7-deee-4884-92f8-54ce6bb439f2
idx_apos = [1; orbit_props["idx_apos"]]

# ╔═╡ 2e0655e6-1351-4495-b56c-bc70bc478d01
@assert idx_apos |> issorted

# ╔═╡ d5591e1b-ec3d-46fc-856e-ce8d1cf6bb03
function plot_orbit_trace!(orbit, idx; d_idx=10)
	if idx-d_idx < idx_apos[1]
		idx_last = 1
	else
		idx_last = filter(x->x .<= idx - d_idx, idx_apos)[end]
	end
	if 	idx + d_idx > idx_apos[end]
		idx_next = length(orbit)
	else
		idx_next = filter(x->x .>= idx .+ d_idx, idx_apos)[1]
	end

	x = orbit.positions[2, idx_last:idx-1]
	y = orbit.positions[3, idx_last:idx-1]
	lines!(x, y, linewidth=theme(:linewidth)[]/2, color=(:white, 0.3), linestyle=:dot)

	# x = orbit.positions[2, idx:idx_next-2]
	# y = orbit.positions[3, idx:idx_next-2]
	# lines!(x, y, linestyle=:dot, linewidth=theme(:linewidth)[]/2, color=(:white, 0.3))

end

# ╔═╡ c689f2b1-eded-4284-9b40-67c72016de8c
function plot_scalebar!(scale_length=50, plotrange=bins[end]-bins[1])
	length_relative = scale_length / plotrange

	x0, y0 = 0.05, 0.05
	lines!([x0, x0 + length_relative], [y0, y0], color=:white, space=:relative, linewidth=theme(:linewidth)[] / 2)
	text!(x0, y0, text="$scale_length kpc", color=:white, space=:relative, fontsize=0.8 * theme(:fontsize)[], )
end

# ╔═╡ c98b7028-73ca-4fda-8ac9-15ace0af0b6c
plotrange = 5

# ╔═╡ 59c11f0e-798f-4102-a338-e1db7029d59c
function get_zoom_histogram(snap, bins=nothing; weights=nothing, plotrange=plotrange, N = 128)

	xycen = snap.x_cen[2:3]
	
	if isnothing(bins)
		
		bins = (LinRange(xycen[1] - plotrange, xycen[1] + plotrange, N),
				LinRange(xycen[2] - plotrange, xycen[2] + plotrange, N)
			   )
	end

	
	hist = get_histogram(snap, bins, weights=weights)

	return hist
end

# ╔═╡ 3178e7d1-151b-4c67-9fe9-57742ca01802
md"""
# Zoom 
"""

# ╔═╡ 65e70122-16ed-4242-b407-043e0d8276f2
x_vec = [sind(0), cosd(0), 0]

# ╔═╡ 93ae35f0-b6f6-4249-8355-e83731601c1e
y_vec = [-sind(0), 0, cosd(0)]

# ╔═╡ fa0a5bf8-7f96-4a53-82dc-a368c12c5241
idx_f = orbit_props["idx_f"]

# ╔═╡ 903528cb-4e12-4659-9ccc-8377aa7250f9
snap_f = out[idx_f]

# ╔═╡ 35635766-eca4-495e-bde1-ba8fb151552d
snap_i = out[1]

# ╔═╡ 3b999809-7bf7-48f0-bd6e-f66be85ba314
colorrange_dm = (-5, 0) .+ log10(maximum(get_zoom_histogram(snap_i)[2]))

# ╔═╡ 9cb67d7e-4c97-4d7f-8705-f737c8f19b88
colorrange_stars = (-5, 0) .+ log10(maximum(get_zoom_histogram(snap_i, weights=snap_i.weights)[2]))

# ╔═╡ 0e9e54b8-0f14-4bb0-8262-83f9b63955e7
import LinearAlgebra: ⋅

# ╔═╡ 385d0cee-6758-49ea-933d-82b654b9e1a0
fg_color = :grey

# ╔═╡ 2c745440-049c-4a15-b4ca-ade29b9e69db
function inset_axis(gs; kwargs...)

	ax = Axis(gs, width=Relative(0.2), height=Relative(0.2), backgroundcolor=:black, 
				leftspinecolor=fg_color, rightspinecolor=fg_color, 
			  bottomspinecolor=fg_color, topspinecolor=fg_color; kwargs...)

	hidedecorations!()
	ax
end

# ╔═╡ 14546575-917a-4799-a6eb-d84d8a890c89
function plot_hist!(binshist; kwargs...)
	bins, hist = binshist
	image!(extrema.(bins)..., log10.(hist);  kwargs...)
end

# ╔═╡ f3da6ca5-e988-4f6b-8efa-18fe32887dab
function plot_both(snap; kwargs...)
	fig = Figure()
	ax = Axis(fig[1,1], backgroundcolor=:black,)

	plot_both(ax, snap; kwargs...)

	hidedecorations!()

	rowsize!(fig.layout, 1, Aspect(1, 1))

	resize_to_layout!(fig)

	fig
end

# ╔═╡ bfa1ee5a-473d-4f21-ac37-4665ca1eade6
function to_transparent_cmap(color)
	return ([Makie.RGBAf(color.r, color.g, color.b, alpha) for alpha in LinRange(0, 1., 100)])
end

# ╔═╡ 5ace10dd-3d95-4a30-bfd0-1c7e846af2c7
colormap_stars = to_transparent_cmap(RGBf(1., 1., 1.))

# ╔═╡ d546247e-7c5c-4452-a1b6-60cbc9fd133f
function to_black_cmap(color)
	return ([Makie.RGBf(color.r*alpha, color.g*alpha, color.b*alpha) for alpha in LinRange(0, 1., 100)])
end


# ╔═╡ 124f1c38-3d15-42c3-a48c-dbdf14aba497
colormap_dm = (to_black_cmap(COLORS[5]))

# ╔═╡ 666eff06-c57a-4219-ac5e-e68e6b860882
function plot_xy_density!(snap)

	h = get_histogram(snap)
	plot_hist!(h, colorrange=colorrange_dm, colormap = colormap_dm)

	h = get_histogram(snap, weights=snap.weights)
	plot_hist!(h, colorrange=colorrange_stars, colormap = colormap_stars)

end


# ╔═╡ e5d118a8-b0dc-4e4d-93b9-02b129c0fe1f
function plot_both(ax, snap)

	hist = get_zoom_histogram(snap)

	
	hidedecorations!()

	p = plot_hist!(hist, colormap=colormap_dm, colorrange=colorrange_dm)

	hist = get_zoom_histogram(snap, weights=snap.weights)

	p = plot_hist!(hist, colormap=colormap_stars,
		colorrange=colorrange_stars)

	ax
end

# ╔═╡ f251ed8f-c908-47b8-96ca-174ecc73f30f
let
	fig = Figure()

	ax = Axis(fig[1, 1], backgroundcolor=:black)

	idx = orbit_props["idx_f"]

	plot_xy_density!(snap_f)
	plot_orbit_trace!(orbit, idx)
	plot_scalebar!(50, bins[end]-bins[1])


	# MW isocontour. same if using x or y
	poly!(xz_iso.x, xz_iso.z, color=COLORS[8])


	# zoom box
	xycen = snap_f.x_cen[2:3]
	x1 = xycen[1] - plotrange
	x2 = xycen[1] + plotrange
	y1, y2 = xycen[2] - plotrange, xycen[2] + plotrange

	lines!([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], color=fg_color, linewidth=theme(:linewidth)[]/2)


	# legend
	text!(0.95, 0.05, text="stars", color=:white, space=:relative, align=(:right, :center), offset=(0, 10))
	text!(0.95, 0.05, text="dark matter", color=COLORS[5], space=:relative, align=(:right, :center))


	hidexdecorations!()
	hideydecorations!()

	ax_initial = inset_axis(fig[1,1], halign=0.05, valign=0.95)
	plot_both(ax_initial, snap_i)
	text!(0.1, 0.9, text="initial", color=:white, space=:relative, fontsize=8, align=(:left, :top))
	plot_scalebar!(1, 2*plotrange)




	ax_final = inset_axis(fig[1,1], halign=0.95, valign=0.95)
	plot_both(ax_final, snap_f)
	text!(0.1, 0.9, text="final", color=:white, space=:relative, fontsize=8, align=(:left, :top))

	rowsize!(fig.layout, 1, Aspect(1, 1))



	resize_to_layout!(fig)

	@savefig "umi_zoom_image"

	fig
end

# ╔═╡ 51d13973-739e-4697-9f18-4b25334213a6
plot_both(snap_i)

# ╔═╡ cffa0546-2adc-42a0-ae0d-beab2809f1f0
plot_both(snap_f)

# ╔═╡ 0bf02f93-581b-4ffd-88a8-6f2236f7a62e
RGBAf(COLORS[1].r, COLORS[1].g, COLORS[1].b, 0.5)

# ╔═╡ 37d40762-62b3-4211-8933-1c7cd1b8bc71
promote(COLORS[1].r, COLORS[1].g, COLORS[1].b, 0.5)

# ╔═╡ Cell order:
# ╠═913e0316-a04b-4270-ba31-0ba0f7fdd705
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═34c0f4f8-ef34-4b58-bcdd-7de69b58db2d
# ╠═b0a0dddc-fb5a-4bb2-b048-a54859d0b703
# ╠═bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═280e0cef-896a-4c14-a1ac-a611418be57e
# ╠═4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
# ╠═53cdcc20-96d4-4bd5-8028-df9a38af71ae
# ╠═0f71807d-d698-4164-9f30-49af8dd8ba55
# ╠═b593188f-1618-44d4-b38c-d5be448e927c
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═b9600d93-5946-4380-a2ae-9b5f673bbaf5
# ╠═4adea59d-a467-4b64-a7f3-e1444c2b35c0
# ╠═32db23d9-7959-41ac-aff4-b63df5e4b94a
# ╠═21d0e542-08ea-4f56-a07d-f1a8e7e69019
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═df61eda6-6d70-4ce4-8f6b-6936dc12283d
# ╠═a3be2d61-98eb-4037-afb4-4155ba24cc21
# ╠═d50c120c-d72a-4fdc-b18a-9ec185c25c04
# ╠═b2cdaa3a-1c20-4609-b6f7-5b0a5f572102
# ╠═5a40b893-021b-46e5-a115-0284e13ae7ae
# ╠═04c8c564-f768-4845-ac22-898173a10582
# ╠═59c11f0e-798f-4102-a338-e1db7029d59c
# ╠═666eff06-c57a-4219-ac5e-e68e6b860882
# ╠═3b999809-7bf7-48f0-bd6e-f66be85ba314
# ╠═9cb67d7e-4c97-4d7f-8705-f737c8f19b88
# ╠═4ae004da-7e17-4f7c-a11d-6ad97dcbe6dd
# ╠═d57c2310-706f-4a0f-af6e-f4dae0cd8376
# ╠═962d5024-6afa-452c-bb04-85046971ee2a
# ╠═971bb9b8-6807-40f7-ad3e-c9e62ab4f3f1
# ╠═800caee7-deee-4884-92f8-54ce6bb439f2
# ╠═2e0655e6-1351-4495-b56c-bc70bc478d01
# ╠═d5591e1b-ec3d-46fc-856e-ce8d1cf6bb03
# ╠═c689f2b1-eded-4284-9b40-67c72016de8c
# ╠═c98b7028-73ca-4fda-8ac9-15ace0af0b6c
# ╠═f251ed8f-c908-47b8-96ca-174ecc73f30f
# ╠═2c745440-049c-4a15-b4ca-ade29b9e69db
# ╠═3178e7d1-151b-4c67-9fe9-57742ca01802
# ╠═65e70122-16ed-4242-b407-043e0d8276f2
# ╠═93ae35f0-b6f6-4249-8355-e83731601c1e
# ╠═fa0a5bf8-7f96-4a53-82dc-a368c12c5241
# ╠═903528cb-4e12-4659-9ccc-8377aa7250f9
# ╠═35635766-eca4-495e-bde1-ba8fb151552d
# ╠═0e9e54b8-0f14-4bb0-8262-83f9b63955e7
# ╠═385d0cee-6758-49ea-933d-82b654b9e1a0
# ╠═14546575-917a-4799-a6eb-d84d8a890c89
# ╠═e5d118a8-b0dc-4e4d-93b9-02b129c0fe1f
# ╠═f3da6ca5-e988-4f6b-8efa-18fe32887dab
# ╠═51d13973-739e-4697-9f18-4b25334213a6
# ╠═cffa0546-2adc-42a0-ae0d-beab2809f1f0
# ╠═124f1c38-3d15-42c3-a48c-dbdf14aba497
# ╠═5ace10dd-3d95-4a30-bfd0-1c7e846af2c7
# ╠═bfa1ee5a-473d-4f21-ac37-4665ca1eade6
# ╠═d546247e-7c5c-4452-a1b6-60cbc9fd133f
# ╠═0bf02f93-581b-4ffd-88a8-6f2236f7a62e
# ╠═37d40762-62b3-4211-8933-1c7cd1b8bc71
