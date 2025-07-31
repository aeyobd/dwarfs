### A Pluto.jl notebook ###
# v0.20.13

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
modelname = "1e6_v31_r3.2/orbit_15_100"

# ╔═╡ 75c691d6-77ab-4715-a645-88352879313c
starsname = "exp2d_rs0.10"

# ╔═╡ 6feb26d5-8b1a-4913-95a2-d57c237fdadb
idx_f = 98

# ╔═╡ 258525e6-7fa1-46b3-95f6-950c2b82aafc
idx_i = 1

# ╔═╡ ea3d8d00-11ee-4b60-9b6a-58a3e396a572
md"""
# Setup
"""

# ╔═╡ 221247cd-e447-400c-8491-b7248787f7a7
import TOML

# ╔═╡ 70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "idealized", modelname)

# ╔═╡ 76fdb3ab-60a6-42ea-bb71-0204d37284ce
stars_dir_in = joinpath(model_dir, "../stars/$starsname")

# ╔═╡ b1edbdfc-298b-40e2-b6c0-687a866fe180
stars_dir_out = joinpath(model_dir, "stars/$starsname")

# ╔═╡ 01a46bcb-7d35-446e-8de8-656e1fec9a54
	idxs = [1, 18]

# ╔═╡ b51a4ddd-90fd-415f-8a48-e4f35aeee801
FIGDIR = "figures"

# ╔═╡ 81684bfc-3316-42c2-a2a9-d9db98858edd
using LilGuys; FIGDIR

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


# ╔═╡ 1cc6563a-5945-4cdf-9cf9-ae54bfb24499
(profiles[1])

# ╔═╡ eb1df5f9-46b7-409c-a052-9d3ae8a01db8
labels = ["initial", "final"]

# ╔═╡ 7f2cf5d2-aa67-428e-83b9-09b04cb28082
r_bs = [get_r_b(i) for i in idxs]

# ╔═╡ 325c0885-9f69-4ed1-b787-303a3fcdda11
colors = [COLORS[1], COLORS[2], COLORS[3]]

# ╔═╡ 7af97b96-2bf5-47d7-973b-dc7ad8c46840
snap_f = out[idx_f]

# ╔═╡ 5d32c845-5cce-430a-945a-957cd7400559
Mtot = 3e6

# ╔═╡ f80faf4b-5740-44d7-9160-a40028922b55
sum(snap_f.weights)

# ╔═╡ 90ff0dc7-2169-4af5-9ea6-d9fbde4cce45
sum(LilGuys.surface_density(profiles[1]) .* diff(π .* LilGuys.radii_bins(profiles[1]) .^ 2 ))

# ╔═╡ 697e034a-bf10-4afa-b95c-37dcdc3b6cb5
arrowcolor = :grey

# ╔═╡ 957b8d3f-03ca-4efa-b774-555b1ef203a6
import DensityEstimators

# ╔═╡ 57e7cd2b-e936-4a04-b4ed-0032a5c6b4df
let
	dx = snap_f.positions[1, :] .- snap_f.x_cen[1]
	dy = snap_f.positions[2, :] .- snap_f.x_cen[2]
	w = snap_f.weights

	Rmax = 0.1

	filt = @. dx^2 + dy^2 < Rmax^2

	M = sum(w[filt])

	@info log10(M / (π*Rmax^2) * Mtot)

	h = DensityEstimators.histogram2d(dx, dy, [-Rmax, 0, Rmax], weights=w * Mtot, normalization=:density)
	h.values
end

# ╔═╡ 5dbbd046-8764-484e-89b4-37dee5fb0555
@savefig "idealized_break_radius" let
	fig = Figure()
	ax = Axis(fig[1,1],
		limits = (-1.5, 1.0, 0, 8),
		xlabel = L"log $R$ / kpc",
		ylabel = L"log $\Sigma_\star$ / M$_\odot$ kpc$^{-2}$",
		# yticks = -4:2:2
	)

	t_f = 0
	
	for (i, prof) in enumerate(profiles)
		lines!(LilGuys.log_radii(prof), LilGuys.log_surface_density(prof) .+ log10(Mtot),  color=colors[1], linestyle=[:dot, :solid][i], label = labels[i])


		if i == 2
			x = log10(r_bs[i])
			y = LilGuys.lerp(LilGuys.log_radii(prof), middle.(LilGuys.log_surface_density(prof)))(x) .+ log10(Mtot)
			dy = 1
			h = 0.2 * 8
			t = round(times[i] * T2GYR * 1e3, digits=0)
			arrows2d!([x], [y+h], [0], [-h],  color=colors[3], align=0.5, minshaftlength = 0)
			text!([x], [y+h * 1.5],color=colors[3], text=L"r_\textrm{break}", align=(:left, :bottom), fontsize=14)
		end

	end

	axislegend(position=:lb)


	
	ax_v = Axis(fig[2,1],
		xlabel=L"log $R$ / kpc",
		ylabel=L"\langle v_\textrm{rad}\rangle\ / \ \textrm{km\ s}^{-1}",
		limits=((-1.5, 0.8), (-10, 10))
	)

	for i in eachindex(idxs)
		snap = out[idxs[i]]
		x, y, ye = v_rad_hist(snap)
		lines!(midpoints(x), y * V2KMS, color=colors[1], linestyle=[:dot, :solid][i])
		if i == 2
			

			h = 0.2 * 20 

			arrows2d!([log10(r_bs[i])], [-h], [0], [h],  color=colors[3], align=0.5, minshaftlength = 0)
		end
	end

	hlines!(0, color=:black, linewidth=theme(:linewidth)[]/2)


	linkxaxes!(ax, ax_v)


	ax3 = Axis(fig[1, 2],
			   xlabel = L"$\Delta x$ / kpc",
			   ylabel = L"$\Delta y$ / kpc",

			  )
	bins = (LinRange(-2.5, 2.5, 100), LinRange(-2.5, 2.5, 100))

	dx = snap_f.positions[1, :] .- snap_f.x_cen[1]
	dy = snap_f.positions[2, :] .- snap_f.x_cen[2]
	w = snap_f.weights
	v = LilGuys.radial_velocities(snap_f)

	h = DensityEstimators.histogram2d(dx, dy, bins, weights=w.*Mtot, normalization=:density)
	#p = hist2d!(dx, dy, weights=w .* Mtot, bins=bins, colorscale=log10, colorrange=(3, nothing), normalization=:density)
	p = heatmap!(h, colorscale=log10, colorrange=(1e3, 1e8),)

	mw_hat = -snap_f.x_cen[1:2] ./ radii(snap_f.x_cen)
	
	arrows2d!([0], [0], mw_hat[1:1], mw_hat[2:2], align=-1.5, minshaftlength=0, color=arrowcolor)
	text!(mw_hat[1]*2.5, mw_hat[2]*2.5, text="MW", color=arrowcolor, align=(:left, :top))


	arc!((0,0), r_bs[2], 0, 2π, color=COLORS[3])


	Colorbar(fig[1, 3], p, label=L"$\Sigma_\star$ / M$_\odot$\,kpc$^{-2}$", ticks=Makie.automatic,)


	


	ax_rv2d = Axis(fig[2, 2],
			   xlabel = L"$\Delta x$ / kpc",
			   ylabel = L"$\Delta y$ / kpc",)


	R = 2.5
	h_sum = DensityEstimators.histogram2d(dx, dy, bins, weights=w .* v .* V2KMS,)
	h_norm = DensityEstimators.histogram2d(dx, dy, bins, weights=w )

	h_ave = h_sum 
	h_ave.values ./= h_norm.values

	p = heatmap!(h_ave, colorrange=(-10, 10), colormap=:bluesreds)
	arc!((0,0), r_bs[2], 0, 2π, color=COLORS[3])



	v_hat = snap_f.v_cen[1:2] ./ radii(snap_f.v_cen)

	arrows2d!([0], [0], v_hat[1:1], v_hat[2:2], align=-1.5, minshaftlength=0, color=arrowcolor)
	text!(v_hat[1]*2.5, v_hat[2]*2.5, text="orbit", color=arrowcolor)


	colsize!(fig.layout, 1, Aspect(1,1))
	colsize!(fig.layout, 2, Aspect(1,1))
	Colorbar(fig[2, 3], p,  label=L"\langle v_\textrm{rad}\rangle\ / \ \textrm{km\ s}^{-1}",
)
	
	fig
end

# ╔═╡ d3c68651-4f68-4f59-b410-91c877275e92
snap_f.v_cen

# ╔═╡ bccf2ecd-325b-4ac7-a930-d33c7b07b0d1
snap_f.x_cen

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
# ╠═81684bfc-3316-42c2-a2a9-d9db98858edd
# ╠═01240d81-5b6e-41b6-a2be-32b87a6568ef
# ╠═70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
# ╠═76fdb3ab-60a6-42ea-bb71-0204d37284ce
# ╠═b1edbdfc-298b-40e2-b6c0-687a866fe180
# ╠═01a46bcb-7d35-446e-8de8-656e1fec9a54
# ╠═a11b127d-f565-43b2-8417-55c2814c367c
# ╠═b51a4ddd-90fd-415f-8a48-e4f35aeee801
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
# ╠═1cc6563a-5945-4cdf-9cf9-ae54bfb24499
# ╠═eb1df5f9-46b7-409c-a052-9d3ae8a01db8
# ╠═7f2cf5d2-aa67-428e-83b9-09b04cb28082
# ╠═325c0885-9f69-4ed1-b787-303a3fcdda11
# ╠═7af97b96-2bf5-47d7-973b-dc7ad8c46840
# ╠═5d32c845-5cce-430a-945a-957cd7400559
# ╠═f80faf4b-5740-44d7-9160-a40028922b55
# ╠═90ff0dc7-2169-4af5-9ea6-d9fbde4cce45
# ╠═57e7cd2b-e936-4a04-b4ed-0032a5c6b4df
# ╠═697e034a-bf10-4afa-b95c-37dcdc3b6cb5
# ╠═5dbbd046-8764-484e-89b4-37dee5fb0555
# ╠═957b8d3f-03ca-4efa-b774-555b1ef203a6
# ╠═d3c68651-4f68-4f59-b410-91c877275e92
# ╠═bccf2ecd-325b-4ac7-a930-d33c7b07b0d1
