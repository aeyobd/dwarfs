### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 5fa9aec4-f876-11ef-1eed-bd8da0982d2b
begin
	using Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	using PlutoUI
end

# ╔═╡ f9bc723e-5c29-4a21-bf72-55fd8396a807
include("./paper_style.jl")

# ╔═╡ 01240d81-5b6e-41b6-a2be-32b87a6568ef
include("./paper_style.jl")

# ╔═╡ 4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
modelname = "1e6/orbit_10_100"

# ╔═╡ 75c691d6-77ab-4715-a645-88352879313c
starsname = "exp2d_rs0.25"

# ╔═╡ 6feb26d5-8b1a-4913-95a2-d57c237fdadb
idx_f = 48

# ╔═╡ 258525e6-7fa1-46b3-95f6-950c2b82aafc
idx_i = 1

# ╔═╡ ea3d8d00-11ee-4b60-9b6a-58a3e396a572
md"""
# Setup
"""

# ╔═╡ 221247cd-e447-400c-8491-b7248787f7a7
import TOML

# ╔═╡ 0089a9dd-314a-47be-b505-1f0491b248b2
import DensityEstimators: histogram2d

# ╔═╡ b51a4ddd-90fd-415f-8a48-e4f35aeee801
FIGDIR = "figures"

# ╔═╡ 81684bfc-3316-42c2-a2a9-d9db98858edd
using LilGuys; FIGDIR

# ╔═╡ 70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "idealized", modelname)

# ╔═╡ 76fdb3ab-60a6-42ea-bb71-0204d37284ce
stars_dir_in = joinpath(model_dir, "../stars/$starsname")

# ╔═╡ b1edbdfc-298b-40e2-b6c0-687a866fe180
stars_dir_out = joinpath(model_dir, "stars/$starsname")

# ╔═╡ 01a46bcb-7d35-446e-8de8-656e1fec9a54
idxs = [1, idx_f]

# ╔═╡ bc9bc9e3-51c5-44ee-b04b-56c66a15246b
colormap = Reverse(:Greys)

# ╔═╡ b8d70c2e-30bc-43f6-965f-c68399ba54d1
md"""
## Data Loading
"""

# ╔═╡ 529fb229-e633-4c73-a994-86e1130b5c9e
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
orbit_props["idx_peris"]

# ╔═╡ 74c6a109-7b1b-4a84-8744-feef39f6f2a6
orbit_props["idx_apos"]

# ╔═╡ 98fe6fe6-7db2-460a-8eec-65cdf316bc6e
idx_peri = orbit_props["idx_peris"][4]

# ╔═╡ 0aa28c43-ebda-4261-8292-408039394a84
r_J_kpc = TOML.parsefile(joinpath(model_dir, "jacobi.toml"))["r_J_kpc"]

# ╔═╡ 025a0ed4-90f9-469e-80ca-bb5258100eb9
df_probs = LilGuys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")

# ╔═╡ 6f8aea0a-91e6-40f4-863f-33a8cca3aee4
p_idx = df_probs.index; probabilities = df_probs.probability

# ╔═╡ 932fd989-6856-42cf-a6e8-d5453af07c64
@assert isperm(p_idx)

# ╔═╡ 5ea02b02-1ab2-4ed3-92d6-214a1772c01e
out = Output(model_dir, weights=probabilities)

# ╔═╡ a11b127d-f565-43b2-8417-55c2814c367c
times = out.times[idxs] .- out.times[idx_peri]

# ╔═╡ abe2deee-32e3-43d5-8232-10aeb797bac1
md"""
# Analysis
"""

# ╔═╡ de8bbe10-0969-492e-937c-61288ea7bc5d
function v_rad_hist(snap, bins=80)

	mass = snap.weights

	v_rad = LilGuys.radial_velocities(snap)
	logr = log10.(radii(snap))
	x_bins, v_bins, err = LilGuys.histogram(logr, bins, weights=v_rad .* mass, errors=:weighted)
	_, counts, _ = LilGuys.histogram(logr, x_bins, weights=mass, errors=:weighted)


	return x_bins, v_bins ./ counts, err ./ counts
end

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

# ╔═╡ 37bf0c8d-9ef3-4cf8-b818-ff37f2d856eb
profiles = [LilGuys.SurfaceDensityProfile(out[i]) for i in idxs]


# ╔═╡ eb1df5f9-46b7-409c-a052-9d3ae8a01db8
labels = ["initial", "final"]

# ╔═╡ 7f2cf5d2-aa67-428e-83b9-09b04cb28082
r_bs = [get_r_b(i) for i in idxs]

# ╔═╡ 7af97b96-2bf5-47d7-973b-dc7ad8c46840
snap_f = out[idx_f]

# ╔═╡ 5d32c845-5cce-430a-945a-957cd7400559
Mtot = 3e6

# ╔═╡ f80faf4b-5740-44d7-9160-a40028922b55
sum(snap_f.weights)

# ╔═╡ 90ff0dc7-2169-4af5-9ea6-d9fbde4cce45
sum(LilGuys.surface_density(profiles[1]) .* diff(π .* LilGuys.radii_bins(profiles[1]) .^ 2 ))

# ╔═╡ c649bd38-960b-450d-b236-c17991181fbb
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ 2b697ba5-19de-45f9-a31a-695f52455c21
smallfontsize = theme(:fontsize)[] * 0.8

# ╔═╡ 486d933f-7ee5-4388-9f42-7285dfc59df5
recent_orbit = out.x_cen[:, idx_f-2:idx_f+2] .- out.x_cen[:, idx_f]

# ╔═╡ 6714b625-0dfe-488a-bc7b-fd8100be7b76
v_cen = out.v_cen[:, idx_f] / radii(out.v_cen[:, idx_f])

# ╔═╡ edd72799-a2fe-4528-b1d7-8ae416096049
r_b = get_r_b(idx_f)

# ╔═╡ 3a55922b-a042-4870-9c2e-111c2f6e0870
past_orbit = [recent_orbit[1, 1] -v_cen[1] * r_b
			  recent_orbit[2, 1] -v_cen[2] * r_b
			 ]

# ╔═╡ 62ccf5f7-f891-4714-8e87-5b66c44c7019
future_orbit = [recent_orbit[1, end] v_cen[1] * r_b
			  recent_orbit[2, end] v_cen[2] * r_b
			 ]

# ╔═╡ bef7da54-c316-4c86-b858-b6b00758022d
md"""
# Plots
"""

# ╔═╡ fab04bb1-7ab6-404a-8813-23b530562348
densityrange = (0, 8)

# ╔═╡ 93852850-104e-4359-bf2f-96b2ce4f69ef
function plot_density_i_f!(ax, profiles)
	
	for (i, prof) in enumerate(profiles)
		lines!(ax, LilGuys.log_radii(prof), LilGuys.log_surface_density(prof) .+ log10(Mtot),  color=COLORS[1], linestyle=[:dot, :solid][i], label = labels[i])


		if i == 2
			x = log10(r_bs[i])
			y = LilGuys.lerp(LilGuys.log_radii(prof), middle.(LilGuys.log_surface_density(prof)))(x) .+ log10(Mtot)
			annotation!(ax, 0, 36, x, y,  color=COLORS[3], )
			text!(x, y, offset=(5, 28), text="break", color=COLORS[3],rotation=0π/2, align=(:left, :top), fontsize=smallfontsize)
		end

	end
end

# ╔═╡ dfab9f31-9a62-4740-b3a0-e5611b06ad6b
function plot_v_rad_hist!(ax, snap_i, snap_f)
	for i in 1:2
		snap = [snap_i, snap_f][i]
		x, y, ye = v_rad_hist(snap)
		lines!(ax, midpoints(x), y * V2KMS, color=COLORS[1], linestyle=[:dot, :solid][i])
		if i == 2
			h = 0.2 * 20 

			annotation!(0, -40, log10(r_bs[i]), 0,  color=COLORS[3])
		end
	end

	hlines!(0, color=:black, linewidth=theme(:linewidth)[]/2)

end

# ╔═╡ 8fd72897-bcea-4f57-a00d-5b5fa7a995ba
function plot_stars_2d(gs, snap_f; r_b=nothing, r_max=4, title="")
	ax = Axis(gs, xlabel=L"\Delta x\,/\,\textrm{kpc}", ylabel = L"\Delta y\,/\,\textrm{kpc}", title=title)
	bins = (LinRange(-r_max, r_max, 100), LinRange(-r_max, r_max, 100))

	dx = snap_f.positions[1, :] .- snap_f.x_cen[1]
	dy = snap_f.positions[2, :] .- snap_f.x_cen[2]
	w = snap_f.weights
	v = LilGuys.radial_velocities(snap_f)

	h = histogram2d(dx, dy, bins, weights=w.* Mtot, normalization=:density)

	p = heatmap!(ax, h.xbins, h.ybins, log10.(h.values), colorrange=densityrange, colormap=colormap)

	mw_hat = -snap_f.x_cen[1:2] ./ radii(snap_f.x_cen)


	if !isnothing(r_b)
		arc!((0,0), r_b, 0, 2π, color=COLORS[3], linewidth=smalllinewidth, linestyle=:dash)
		text!(-r_b, 0, color=COLORS[3], text="break", align=(:center, :bottom), offset=(-3, 0), rotation=π/2)
	end

	limits!(-r_max, r_max, -r_max, r_max)
	p
end

# ╔═╡ 5dbbd046-8764-484e-89b4-37dee5fb0555
@savefig "idealized_break_radius" let
	fig = Figure(figure_padding=(12, 12, 24, 12))

	p = plot_stars_2d(fig[1, 1], out[1], title="initial")
	p = plot_stars_2d(fig[1, 2], snap_f, r_b=get_r_b(idx_f), title="final")
	hideydecorations!(ticks=false, minorticks=false)
	lines!(past_orbit[1, :], past_orbit[2, :], color=COLORS[1], linestyle=:solid, linewidth=smalllinewidth)
	lines!(future_orbit[1, :], future_orbit[2, :], color=COLORS[1], linestyle=:dot, linewidth=smalllinewidth)

	text!(-1.5*r_b*v_cen[1], -1.5*r_b*v_cen[2], text="orbit", color=COLORS[1], rotation=π + atan(v_cen[2], v_cen[1]))
	
	Colorbar(fig[1, 3], p, label=L"$\log\,\Sigma_\star$ / M$_\odot$\,kpc$^{-2}$", ticks=Makie.automatic)
	rowsize!(fig.layout, 1, Aspect(1, 1))


	gs_lower = GridLayout(fig[2, :])
	ax = Axis(gs_lower[1, 1],
		limits = ((-1.5, 1.0), densityrange),
		xlabel = L"log $R$ / kpc",
		ylabel = L"log $\Sigma_\star$ / M$_\odot$ kpc$^{-2}$",
	)

	plot_density_i_f!(ax, profiles)
	axislegend(position=:lb)



	ax_v = Axis(gs_lower[1, 2],
		xlabel=L"log $R$ / kpc",
		ylabel=L"\langle v_\textrm{rad}\rangle\ / \ \textrm{km\ s}^{-1}",
		limits=((-1.5, 0.8), (-10, 10)),
		ylabelpadding=-4
	)
	plot_v_rad_hist!(ax_v, out[1], snap_f)

	rowsize!(fig.layout, 2, Aspect(1, 1))

	resize_to_layout!()
	fig
end

# ╔═╡ 9944c205-ace3-4ff6-b00c-79be8f98054f


# ╔═╡ Cell order:
# ╠═4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
# ╠═75c691d6-77ab-4715-a645-88352879313c
# ╠═6feb26d5-8b1a-4913-95a2-d57c237fdadb
# ╠═258525e6-7fa1-46b3-95f6-950c2b82aafc
# ╠═c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
# ╠═74c6a109-7b1b-4a84-8744-feef39f6f2a6
# ╠═98fe6fe6-7db2-460a-8eec-65cdf316bc6e
# ╠═f9bc723e-5c29-4a21-bf72-55fd8396a807
# ╟─ea3d8d00-11ee-4b60-9b6a-58a3e396a572
# ╠═5fa9aec4-f876-11ef-1eed-bd8da0982d2b
# ╠═221247cd-e447-400c-8491-b7248787f7a7
# ╠═0089a9dd-314a-47be-b505-1f0491b248b2
# ╠═81684bfc-3316-42c2-a2a9-d9db98858edd
# ╠═b51a4ddd-90fd-415f-8a48-e4f35aeee801
# ╠═01240d81-5b6e-41b6-a2be-32b87a6568ef
# ╠═70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
# ╠═76fdb3ab-60a6-42ea-bb71-0204d37284ce
# ╠═b1edbdfc-298b-40e2-b6c0-687a866fe180
# ╠═01a46bcb-7d35-446e-8de8-656e1fec9a54
# ╠═a11b127d-f565-43b2-8417-55c2814c367c
# ╠═bc9bc9e3-51c5-44ee-b04b-56c66a15246b
# ╟─b8d70c2e-30bc-43f6-965f-c68399ba54d1
# ╠═529fb229-e633-4c73-a994-86e1130b5c9e
# ╠═0aa28c43-ebda-4261-8292-408039394a84
# ╠═025a0ed4-90f9-469e-80ca-bb5258100eb9
# ╠═932fd989-6856-42cf-a6e8-d5453af07c64
# ╠═6f8aea0a-91e6-40f4-863f-33a8cca3aee4
# ╠═5ea02b02-1ab2-4ed3-92d6-214a1772c01e
# ╟─abe2deee-32e3-43d5-8232-10aeb797bac1
# ╠═de8bbe10-0969-492e-937c-61288ea7bc5d
# ╠═a2e625a6-8a9d-4ce1-b6a4-4b64c6e445b7
# ╠═fedadc2b-2ea7-4bf0-9d47-e8f2875d974d
# ╠═37bf0c8d-9ef3-4cf8-b818-ff37f2d856eb
# ╠═eb1df5f9-46b7-409c-a052-9d3ae8a01db8
# ╠═7f2cf5d2-aa67-428e-83b9-09b04cb28082
# ╠═7af97b96-2bf5-47d7-973b-dc7ad8c46840
# ╠═5d32c845-5cce-430a-945a-957cd7400559
# ╠═f80faf4b-5740-44d7-9160-a40028922b55
# ╠═90ff0dc7-2169-4af5-9ea6-d9fbde4cce45
# ╠═c649bd38-960b-450d-b236-c17991181fbb
# ╠═2b697ba5-19de-45f9-a31a-695f52455c21
# ╠═486d933f-7ee5-4388-9f42-7285dfc59df5
# ╠═6714b625-0dfe-488a-bc7b-fd8100be7b76
# ╠═edd72799-a2fe-4528-b1d7-8ae416096049
# ╠═3a55922b-a042-4870-9c2e-111c2f6e0870
# ╠═62ccf5f7-f891-4714-8e87-5b66c44c7019
# ╟─bef7da54-c316-4c86-b858-b6b00758022d
# ╠═fab04bb1-7ab6-404a-8813-23b530562348
# ╠═5dbbd046-8764-484e-89b4-37dee5fb0555
# ╠═93852850-104e-4359-bf2f-96b2ce4f69ef
# ╠═dfab9f31-9a62-4740-b3a0-e5611b06ad6b
# ╠═8fd72897-bcea-4f57-a00d-5b5fa7a995ba
# ╠═9944c205-ace3-4ff6-b00c-79be8f98054f
