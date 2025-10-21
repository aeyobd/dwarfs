### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 5fa9aec4-f876-11ef-1eed-bd8da0982d2b
begin
	using Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	using PlutoUI
end

# ╔═╡ 01240d81-5b6e-41b6-a2be-32b87a6568ef
include("./paper_style.jl")

# ╔═╡ 4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
modelname = "1e6/orbit_20_100"

# ╔═╡ 75c691d6-77ab-4715-a645-88352879313c
starsname = "exp2d_rs0.25"

# ╔═╡ 6feb26d5-8b1a-4913-95a2-d57c237fdadb
idx_f = 79

# ╔═╡ 258525e6-7fa1-46b3-95f6-950c2b82aafc
idx_i = 1

# ╔═╡ 5d32c845-5cce-430a-945a-957cd7400559
Mtot = 5e5

# ╔═╡ ea3d8d00-11ee-4b60-9b6a-58a3e396a572
md"""
# Setup
"""

# ╔═╡ a4ffe2d7-5ade-4167-92e7-e49badb81671
import DataFrames: DataFrame

# ╔═╡ 221247cd-e447-400c-8491-b7248787f7a7
import TOML

# ╔═╡ 0089a9dd-314a-47be-b505-1f0491b248b2
import DensityEstimators: histogram2d

# ╔═╡ 4c1ab53f-69b4-437c-9713-dc22bfce0db3
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ b51a4ddd-90fd-415f-8a48-e4f35aeee801
FIGDIR = "figures"

# ╔═╡ 81684bfc-3316-42c2-a2a9-d9db98858edd
using LilGuys; FIGDIR

# ╔═╡ b8d70c2e-30bc-43f6-965f-c68399ba54d1
md"""
## Data Loading
"""

# ╔═╡ 70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "idealized", modelname)

# ╔═╡ 76fdb3ab-60a6-42ea-bb71-0204d37284ce
stars_dir_in = joinpath(model_dir, "../stars/$starsname")

# ╔═╡ b1edbdfc-298b-40e2-b6c0-687a866fe180
stars_dir_out = joinpath(model_dir, "stars/$starsname")

# ╔═╡ bc9bc9e3-51c5-44ee-b04b-56c66a15246b
colormap = Reverse(:Greys)

# ╔═╡ 529fb229-e633-4c73-a994-86e1130b5c9e
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
orbit_props["idx_peris"]

# ╔═╡ 74c6a109-7b1b-4a84-8744-feef39f6f2a6
orbit_props["idx_apos"]

# ╔═╡ 025a0ed4-90f9-469e-80ca-bb5258100eb9
df_probs = LilGuys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")

# ╔═╡ 6f8aea0a-91e6-40f4-863f-33a8cca3aee4
p_idx = df_probs.index; probabilities = df_probs.probability

# ╔═╡ 932fd989-6856-42cf-a6e8-d5453af07c64
@assert isperm(p_idx)

# ╔═╡ 5ea02b02-1ab2-4ed3-92d6-214a1772c01e
out = Output(model_dir, weights=probabilities)

# ╔═╡ b1452bd1-fdc2-4f8f-9d8b-774971c98d82
snap_i = out[1]

# ╔═╡ 7af97b96-2bf5-47d7-973b-dc7ad8c46840
snap_f = out[idx_f]

# ╔═╡ abe2deee-32e3-43d5-8232-10aeb797bac1
md"""
# Analysis
"""

# ╔═╡ a2e625a6-8a9d-4ce1-b6a4-4b64c6e445b7
function get_delta_t(idx)
	if idx <=orbit_props["idx_peris"][1]
		return 0
	end
	return out.times[idx] - filter(x->x < out.times[idx], orbit_props["t_peris"])[end]
end

# ╔═╡ fedadc2b-2ea7-4bf0-9d47-e8f2875d974d
function get_r_b(idx)
	return LilGuys.break_radius(LilGuys.σv_1d(out[idx]), get_delta_t(idx))
end

# ╔═╡ 7f2cf5d2-aa67-428e-83b9-09b04cb28082
r_break = get_r_b(idx_f)

# ╔═╡ bce9a19a-4261-46c7-b10c-1e1be5e5592e
LilGuys.std(snap_f.velocities[1, :, :], snap_f.weights) * snap_f.time/2 * 0.55 # should be nearly the break radius provided this is first post-peri apo

# ╔═╡ f80faf4b-5740-44d7-9160-a40028922b55
@assert sum(snap_f.weights) ≈ 1

# ╔═╡ 24a53a8d-71d7-4f8d-9a76-52ccbf41fd17
md"""
# Processing
"""

# ╔═╡ 77adef1c-3ba1-4c12-8acd-40024afc55d6
import LinearAlgebra: normalize, ×, ⋅

# ╔═╡ 6979b8f2-7116-4e45-91e6-96e9caa61a97
md"""
coordinate system with galaxy at [x, 0, 0] and angular momentum in [0, 0, z] direction
"""

# ╔═╡ 1bb36506-b428-4336-9d3d-5f5beb5847c3
function get_vectors(snap)
	x_hat = normalize(snap.x_cen)
	y_hat = normalize(snap.v_cen .- (x_hat ⋅ snap.v_cen) * x_hat)
	z_hat = x_hat × y_hat

	return [x_hat y_hat z_hat]'
end

# ╔═╡ 558f4739-5a4a-49bf-929f-b7e1dbbb677c
vec = get_vectors(snap_f)

# ╔═╡ 1739773e-6a3b-4308-ab5d-18253a4d76e8
@assert isapprox(vec * vec', [1 0 0 ; 0 1 0 ; 0 0 1]) "invalid basis"

# ╔═╡ 60f32d04-9554-49c6-8fb0-5a46385617e5
vec * (snap_f.x_cen × snap_f.v_cen) #should be positive z

# ╔═╡ 2e0325ed-08ed-437d-b97e-86a18ebf965e
get_vectors(snap_f) * snap_f.x_cen # should be positive x only

# ╔═╡ 576aa32f-2802-4f46-a034-5fd0acf8d0fb
function make_stars(snap; Mtot=Mtot)


	A = get_vectors(snap)

	xyz = A * snap.positions
	v_xyz = A * snap.velocities * V2KMS
	xyz0 = A*snap.x_cen

	coords = [LilGuys.Cartesian{ICRS}(xyz[:, i], v_xyz[:, i]) for i in 1:length(snap)]
	coords_icrs = LilGuys.transform.(ICRS, coords)

	ra = [c.ra for c in coords_icrs]
	dec = [c.dec for c in coords_icrs]

	icrs0 = LilGuys.transform(ICRS, LilGuys.Cartesian{ICRS}(A*snap.x_cen))
	@info icrs0
	@info "centre", A*snap.x_cen, A*snap.v_cen
	xi, eta = LilGuys.to_tangent(ra, dec, icrs0.ra, icrs0.dec)

	df = DataFrame(
		:index => snap.index,
		:xi => xi,
		:eta => eta,
		:x => xyz[1, :],
		:y => xyz[2, :],
		:z => xyz[3, :],		
		:weights => snap.weights * Mtot,
	)

	df[!, :R] = @. sqrt(df.y^2 + df.z^2)
	df[!, :xi] = (df.y .- xyz0[2]) 
	df[!, :eta] = (df.z .- xyz0[3])
	df[!, :R_ell] = @. sqrt(df.xi .^2 .+df.eta .^2)

	df
end
	

# ╔═╡ 511f581f-8199-46ad-b8d3-6bbc7367d125
stars_i = make_stars(snap_i)

# ╔═╡ 9764853a-c413-4668-a439-a3bb5d10ffe6
stars_f = make_stars(snap_f)

# ╔═╡ 37bf0c8d-9ef3-4cf8-b818-ff37f2d856eb
profiles = [LilGuys.SurfaceDensityProfile(stars_i.R, weights=stars_i.weights), LilGuys.SurfaceDensityProfile(stars_f.R, weights=stars_f.weights)]


# ╔═╡ 9f7b9f37-077c-4682-9754-1eb6dd1182c8
@assert sum(stars_f.weights) ≈ Mtot

# ╔═╡ bef7da54-c316-4c86-b858-b6b00758022d
md"""
# Plots
"""

# ╔═╡ eb1df5f9-46b7-409c-a052-9d3ae8a01db8
labels = ["initial", "final"]

# ╔═╡ c649bd38-960b-450d-b236-c17991181fbb
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ 2b697ba5-19de-45f9-a31a-695f52455c21
smallfontsize = theme(:fontsize)[] * 0.8

# ╔═╡ 2d2b181b-5dca-4321-945b-05fe116a08d3
Σ_max = maximum(middle.(profiles[1].log_Sigma))

# ╔═╡ fab04bb1-7ab6-404a-8813-23b530562348
densityrange = (Σ_max-10, Σ_max)

# ╔═╡ 961101e4-2662-4cb3-9138-c10a035d029d
orbit_direction = (get_vectors(snap_f) * snap_f.v_cen)

# ╔═╡ 64754d5f-af87-4834-a55c-69a665ed2808
function plot_stars_2d(gs, stars; R_h=NaN, R_h_label=false, 
        bins=100, r_max=4*60, colormap=Reverse(:Greys), r_b=nothing, orbit_direction=nothing, norm=0, colorrange=densityrange)

	ax = Axis(gs, 
			 xlabel = L"$y'$ / kpc ", 
			 ylabel = L"$z'$ / kpc ", 
			)

	
    bins = LinRange(-r_max, r_max, bins)
	h = histogram2d(stars.y, stars.z, bins, weights=stars.weights * 10^norm, normalization=:density)

	@info "density maximum: $(maximum(h.values))"
	p = heatmap!(h.xbins, h.ybins, log10.(h.values), 
                 colorrange=colorrange, colormap=colormap)

	ModelUtils.plot_r_b_circle!(r_b)

    if !isnothing(orbit_direction)
        dx, dy = 5 .* orbit_direction[[2, 3]]
		@info dx, dy
        ModelUtils.plot_orbit_arrow!(ax, dx, dy, r_max)
    end

    return p
end

# ╔═╡ 93852850-104e-4359-bf2f-96b2ce4f69ef
function plot_density_i_f!(ax, profiles)
	for (i, prof) in enumerate(profiles)
		lines!(ax, LilGuys.log_radii(prof), LilGuys.log_surface_density(prof),  color=COLORS[1], linestyle=[:dot, :solid][i], label = labels[i])

	end
end

# ╔═╡ de8bbe10-0969-492e-937c-61288ea7bc5d
function v_rad_hist(snap, bins=80)

	mass = snap.weights

	v_rad = LilGuys.radial_velocities(snap)
	logr = log10.(radii(snap))
	x_bins, v_bins, err = LilGuys.histogram(logr, bins, weights=v_rad .* mass, errors=:weighted)
	_, counts, _ = LilGuys.histogram(logr, x_bins, weights=mass, errors=:weighted)


	return x_bins, v_bins ./ counts, err ./ counts
end

# ╔═╡ dfab9f31-9a62-4740-b3a0-e5611b06ad6b
function plot_v_rad_hist!(ax, snap_i, snap_f)
	for i in 1:2
		snap = [snap_i, snap_f][i]
		x, y, ye = v_rad_hist(snap)
		lines!(ax, midpoints(x), y * V2KMS, color=COLORS[1], linestyle=[:dot, :solid][i])
	end

	hlines!(0, color=:black, linewidth=theme(:linewidth)[]/2)

end

# ╔═╡ 5dbbd046-8764-484e-89b4-37dee5fb0555
@savefig "idealized_break_radius" let
	fig = Figure(figure_padding=(12, 12, 24, 12))
	r_max = 10

	p = plot_stars_2d(fig[1,1], stars_i, r_max=r_max,)
	
	plot_stars_2d(fig[1,2], stars_f, r_max=r_max, r_b=-r_break, orbit_direction=orbit_direction)

	hideydecorations!(ticks=false, minorticks=false)


	Colorbar(fig[1, 3], p, label=L"$\log\,\Sigma_\star$ / M$_\odot$\,kpc$^{-2}$", ticks=Makie.automatic)
	rowsize!(fig.layout, 1, Aspect(1, 1))


	gs_lower = GridLayout(fig[2, :])
	ax = Axis(gs_lower[1, 1],
		limits = ((-1.5, log10(r_max)), densityrange),
		xlabel = L"log $R$ / kpc",
		ylabel = L"log $\Sigma_\star$ / M$_\odot$ kpc$^{-2}$",
	)

	plot_density_i_f!(ax, profiles)
	ModelUtils.plot_r_break_arrow!(r_break, 2.7)
	axislegend(position=:lb)



	ax_v = Axis(gs_lower[1, 2],
		xlabel=L"log $R$ / kpc",
		ylabel=L"\langle v_\textrm{rad}\rangle\ / \ \textrm{km\ s}^{-1}",
		limits=((-1.5, log10(r_max)), (-10, 10)),
		ylabelpadding=-4
	)
	plot_v_rad_hist!(ax_v, snap_i, snap_f)
	ModelUtils.plot_r_break_arrow!(r_break, 0)

	rowsize!(fig.layout, 2, Aspect(1, 1))

	resize_to_layout!()
	fig
end

# ╔═╡ b40c15f6-f8f9-4c41-aaeb-ac43728ce33e


# ╔═╡ Cell order:
# ╠═4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
# ╠═75c691d6-77ab-4715-a645-88352879313c
# ╠═6feb26d5-8b1a-4913-95a2-d57c237fdadb
# ╠═258525e6-7fa1-46b3-95f6-950c2b82aafc
# ╠═c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
# ╠═74c6a109-7b1b-4a84-8744-feef39f6f2a6
# ╠═5d32c845-5cce-430a-945a-957cd7400559
# ╟─ea3d8d00-11ee-4b60-9b6a-58a3e396a572
# ╠═5fa9aec4-f876-11ef-1eed-bd8da0982d2b
# ╠═a4ffe2d7-5ade-4167-92e7-e49badb81671
# ╠═221247cd-e447-400c-8491-b7248787f7a7
# ╠═0089a9dd-314a-47be-b505-1f0491b248b2
# ╠═81684bfc-3316-42c2-a2a9-d9db98858edd
# ╠═4c1ab53f-69b4-437c-9713-dc22bfce0db3
# ╠═b51a4ddd-90fd-415f-8a48-e4f35aeee801
# ╠═01240d81-5b6e-41b6-a2be-32b87a6568ef
# ╟─b8d70c2e-30bc-43f6-965f-c68399ba54d1
# ╠═70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
# ╠═76fdb3ab-60a6-42ea-bb71-0204d37284ce
# ╠═b1edbdfc-298b-40e2-b6c0-687a866fe180
# ╠═bc9bc9e3-51c5-44ee-b04b-56c66a15246b
# ╠═529fb229-e633-4c73-a994-86e1130b5c9e
# ╠═025a0ed4-90f9-469e-80ca-bb5258100eb9
# ╠═6f8aea0a-91e6-40f4-863f-33a8cca3aee4
# ╠═932fd989-6856-42cf-a6e8-d5453af07c64
# ╠═5ea02b02-1ab2-4ed3-92d6-214a1772c01e
# ╠═b1452bd1-fdc2-4f8f-9d8b-774971c98d82
# ╠═7af97b96-2bf5-47d7-973b-dc7ad8c46840
# ╟─abe2deee-32e3-43d5-8232-10aeb797bac1
# ╠═a2e625a6-8a9d-4ce1-b6a4-4b64c6e445b7
# ╠═fedadc2b-2ea7-4bf0-9d47-e8f2875d974d
# ╠═37bf0c8d-9ef3-4cf8-b818-ff37f2d856eb
# ╠═7f2cf5d2-aa67-428e-83b9-09b04cb28082
# ╠═bce9a19a-4261-46c7-b10c-1e1be5e5592e
# ╠═f80faf4b-5740-44d7-9160-a40028922b55
# ╟─24a53a8d-71d7-4f8d-9a76-52ccbf41fd17
# ╠═77adef1c-3ba1-4c12-8acd-40024afc55d6
# ╟─6979b8f2-7116-4e45-91e6-96e9caa61a97
# ╠═1bb36506-b428-4336-9d3d-5f5beb5847c3
# ╠═558f4739-5a4a-49bf-929f-b7e1dbbb677c
# ╠═1739773e-6a3b-4308-ab5d-18253a4d76e8
# ╠═60f32d04-9554-49c6-8fb0-5a46385617e5
# ╠═2e0325ed-08ed-437d-b97e-86a18ebf965e
# ╠═576aa32f-2802-4f46-a034-5fd0acf8d0fb
# ╠═511f581f-8199-46ad-b8d3-6bbc7367d125
# ╠═9764853a-c413-4668-a439-a3bb5d10ffe6
# ╠═9f7b9f37-077c-4682-9754-1eb6dd1182c8
# ╟─bef7da54-c316-4c86-b858-b6b00758022d
# ╠═eb1df5f9-46b7-409c-a052-9d3ae8a01db8
# ╠═c649bd38-960b-450d-b236-c17991181fbb
# ╠═2b697ba5-19de-45f9-a31a-695f52455c21
# ╠═fab04bb1-7ab6-404a-8813-23b530562348
# ╠═2d2b181b-5dca-4321-945b-05fe116a08d3
# ╠═961101e4-2662-4cb3-9138-c10a035d029d
# ╠═5dbbd046-8764-484e-89b4-37dee5fb0555
# ╠═64754d5f-af87-4834-a55c-69a665ed2808
# ╠═93852850-104e-4359-bf2f-96b2ce4f69ef
# ╠═dfab9f31-9a62-4740-b3a0-e5611b06ad6b
# ╠═de8bbe10-0969-492e-937c-61288ea7bc5d
# ╠═b40c15f6-f8f9-4c41-aaeb-ac43728ce33e
