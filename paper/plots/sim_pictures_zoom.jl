### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ 3ad8154b-9c07-4666-8b0c-b123276fae7e
include(joinpath(ENV["DWARFS_ROOT"], "orbits/orbit_utils.jl"))

# ╔═╡ 913e0316-a04b-4270-ba31-0ba0f7fdd705
galaxyname = "sculptor_lmc"

# ╔═╡ bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
import TOML

# ╔═╡ 280e0cef-896a-4c14-a1ac-a611418be57e
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ a9023e51-4e24-4526-9e29-bd6e4bfa2b24
module PictureUtils
	include("picture_utils.jl")
end

# ╔═╡ e62f58df-67fd-4cf6-a1a1-fc4a11ec334b
md"""
## Isodensity plots
"""

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b9600d93-5946-4380-a2ae-9b5f673bbaf5
modelname = if galaxyname == "sculptor"
	joinpath(modelnames["scl_smallperi"][1:2]...)
elseif galaxyname == "ursa_minor"
	joinpath(modelnames["umi_smallperi"][1:2]...)
elseif galaxyname == "sculptor_lmc"
	joinpath(modelnames["scl_lmc"][1:2]...)

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
elseif galaxyname == "sculptor_lmc"
	modelnames["scl_lmc"][3]
end

# ╔═╡ 32db23d9-7959-41ac-aff4-b63df5e4b94a
figname = Dict(
	"sculptor" => "scl",
	"ursa_minor" => "umi",
	"sculptor_lmc" => "scl_lmc"
)[galaxyname]

# ╔═╡ 21d0e542-08ea-4f56-a07d-f1a8e7e69019
figtitle = Dict(
	"sculptor" => "Sculptor",
	"sculptor_lmc" => "Sculptor (MW+LMC)",

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

# ╔═╡ 5a40b893-021b-46e5-a115-0284e13ae7ae
bins = LinRange(-150, 150, 512)

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

end

# ╔═╡ cec9cbd6-9e8f-4375-9105-8ced9089dd5e
smallfontsize = 0.8 * theme(:fontsize)[]

# ╔═╡ 13137339-95c2-4a23-b775-6eecb59c2ce6
function plot_lmc_orbit!(galaxyname)
	if galaxyname == "sculptor_lmc"
		lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/vasiliev24_L3M11")) 

		filt = lmc_orbit.times * T2GYR .> -3
		lines!(lmc_orbit.positions[2, filt], lmc_orbit.positions[3, filt], color=COLORS[3], )

		scatter!(lmc_orbit.positions[2, end], lmc_orbit.positions[3, end], color=COLORS[3])
		text!(lmc_orbit.positions[2, end], lmc_orbit.positions[3, end], text="LMC", color=COLORS[3], align=(:right, :center), offset=(-smallfontsize/2, 0), fontsize=smallfontsize)

	end

end
		

# ╔═╡ 2b4e16b0-059c-4dff-b2a9-55684a94134a
md"""
# Plot utils
"""

# ╔═╡ c98b7028-73ca-4fda-8ac9-15ace0af0b6c
plotrange = PictureUtils.plotrange

# ╔═╡ fa0a5bf8-7f96-4a53-82dc-a368c12c5241
idx_f = orbit_props["idx_f"]

# ╔═╡ 903528cb-4e12-4659-9ccc-8377aa7250f9
snap_f = out[idx_f]

# ╔═╡ 35635766-eca4-495e-bde1-ba8fb151552d
snap_i = out[1]

# ╔═╡ 3b999809-7bf7-48f0-bd6e-f66be85ba314
colorrange_dm = (-5, 0) .+ log10(maximum(PictureUtils.get_zoom_histogram(snap_i)[2]))

# ╔═╡ 9cb67d7e-4c97-4d7f-8705-f737c8f19b88
colorrange_stars = (-5, 0) .+ log10(maximum(PictureUtils.get_zoom_histogram(snap_i, weights=snap_i.weights)[2]))

# ╔═╡ 385d0cee-6758-49ea-933d-82b654b9e1a0
fg_color = :grey

# ╔═╡ e5d118a8-b0dc-4e4d-93b9-02b129c0fe1f
function plot_zoom(ax, snap)
	hidedecorations!()
	
	hist = PictureUtils.get_zoom_histogram(snap)

	p = PictureUtils.plot_hist!(hist, interpolate=true, colormap=PictureUtils.colormap_dm, colorrange=colorrange_dm)

	hist = PictureUtils.get_zoom_histogram(snap, weights=snap.weights)

	p = PictureUtils.plot_hist!(hist, interpolate=true, colormap=PictureUtils.colormap_stars,
		colorrange=colorrange_stars)

	ax
end

# ╔═╡ f3da6ca5-e988-4f6b-8efa-18fe32887dab
function plot_zoom(snap; kwargs...)
	fig = Figure()
	ax = Axis(fig[1,1], backgroundcolor=:black,)

	plot_zoom(ax, snap; kwargs...)

	hidedecorations!()

	rowsize!(fig.layout, 1, Aspect(1, 1))

	resize_to_layout!(fig)

	fig
end

# ╔═╡ 2c745440-049c-4a15-b4ca-ade29b9e69db
function inset_axis(gs; kwargs...)

	ax = Axis(gs, width=Relative(0.2), height=Relative(0.2), backgroundcolor=:black, 
				leftspinecolor=fg_color, rightspinecolor=fg_color, 
			  bottomspinecolor=fg_color, topspinecolor=fg_color; kwargs...)

	hidedecorations!()
	ax
end

# ╔═╡ 72493e30-56eb-4766-a7f5-b089edb2921b
light_grey = RGBAf(0.8, 0.8, 0.8)

# ╔═╡ f251ed8f-c908-47b8-96ca-174ecc73f30f
let

	fig = PictureUtils.plot_frame(snap_f, colorrange_dm=colorrange_dm, colorrange_stars=colorrange_stars, legend=true, scalebar_color=light_grey)
	plot_orbit_trace!(orbit, orbit_props["idx_f"])
	plot_lmc_orbit!(galaxyname)

	hidexdecorations!()
	hideydecorations!()
	
	# zoom box
	xycen = snap_f.x_cen[2:3]
	x1 = xycen[1] - plotrange
	x2 = xycen[1] + plotrange
	y1, y2 = xycen[2] - plotrange, xycen[2] + plotrange

	lines!([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], color=fg_color, linewidth=theme(:linewidth)[]/2)

	

	ax_initial = inset_axis(fig[1,1], halign=0.05, valign=0.95)
	plot_zoom(ax_initial, snap_i)
	text!(0.1, 0.9, text="initial", color=light_grey, space=:relative, fontsize=8, align=(:left, :top))
	PictureUtils.plot_scalebar!(1, 2*plotrange, color=light_grey)




	ax_final = inset_axis(fig[1,1], halign=0.95, valign=0.95)
	plot_zoom(ax_final, snap_f)
	text!(0.1, 0.9, text="final", color=light_grey, space=:relative, fontsize=8, align=(:left, :top))
	PictureUtils.plot_scalebar!(1, 2*plotrange, color=light_grey)



	resize_to_layout!(fig)

	@savefig "$(figname)_zoom_image"

	fig
end

# ╔═╡ Cell order:
# ╠═913e0316-a04b-4270-ba31-0ba0f7fdd705
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═280e0cef-896a-4c14-a1ac-a611418be57e
# ╠═a9023e51-4e24-4526-9e29-bd6e4bfa2b24
# ╟─e62f58df-67fd-4cf6-a1a1-fc4a11ec334b
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
# ╠═5a40b893-021b-46e5-a115-0284e13ae7ae
# ╠═3b999809-7bf7-48f0-bd6e-f66be85ba314
# ╠═9cb67d7e-4c97-4d7f-8705-f737c8f19b88
# ╠═4ae004da-7e17-4f7c-a11d-6ad97dcbe6dd
# ╠═d57c2310-706f-4a0f-af6e-f4dae0cd8376
# ╠═962d5024-6afa-452c-bb04-85046971ee2a
# ╠═971bb9b8-6807-40f7-ad3e-c9e62ab4f3f1
# ╠═800caee7-deee-4884-92f8-54ce6bb439f2
# ╠═2e0655e6-1351-4495-b56c-bc70bc478d01
# ╠═d5591e1b-ec3d-46fc-856e-ce8d1cf6bb03
# ╠═3ad8154b-9c07-4666-8b0c-b123276fae7e
# ╠═13137339-95c2-4a23-b775-6eecb59c2ce6
# ╠═cec9cbd6-9e8f-4375-9105-8ced9089dd5e
# ╟─2b4e16b0-059c-4dff-b2a9-55684a94134a
# ╠═c98b7028-73ca-4fda-8ac9-15ace0af0b6c
# ╠═fa0a5bf8-7f96-4a53-82dc-a368c12c5241
# ╠═903528cb-4e12-4659-9ccc-8377aa7250f9
# ╠═35635766-eca4-495e-bde1-ba8fb151552d
# ╠═385d0cee-6758-49ea-933d-82b654b9e1a0
# ╠═e5d118a8-b0dc-4e4d-93b9-02b129c0fe1f
# ╠═f3da6ca5-e988-4f6b-8efa-18fe32887dab
# ╠═2c745440-049c-4a15-b4ca-ade29b9e69db
# ╠═72493e30-56eb-4766-a7f5-b089edb2921b
# ╠═f251ed8f-c908-47b8-96ca-174ecc73f30f
