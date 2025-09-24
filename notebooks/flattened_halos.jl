### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 88d39cc6-95a4-11f0-1b8a-e519835822e7
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 83fea91d-83fe-436a-8e4e-f21a61468562
snap = Snapshot("../analysis/isolation/1e6_oblate_0.5/fiducial/combined.hdf5/1")

# ╔═╡ 8a457a1f-05c9-4cf5-8244-18a364782b74
snap_2 = Snapshot("../analysis/isolation/1e6_oblate_0.5/fiducial/combined.hdf5/-1")

# ╔═╡ 6278077a-b027-4dc9-8304-197b049346e6
CairoMakie.activate!(type=:png)

# ╔═╡ da071214-ff07-438b-a677-b54687beb8ae
LilGuys.plot_xyz(snap.positions, plot=:scatter, markersize=1, alpha=0.1)

# ╔═╡ 67d1b39c-c314-4313-8e5a-5073d9b6a11d
LilGuys.plot_xyz(snap_2.positions, plot=:scatter, markersize=1, alpha=0.1)

# ╔═╡ 75db6089-d5a2-419a-8354-d2a996fc19d4
import DensityEstimators: histogram2d

# ╔═╡ 280caba5-b73d-43f6-9792-68166000bb8d
function plot_isocontours(snap)
	fig = Figure()

	bins = -10:0.5:10
	
	for i in 1:3
		x = snap.positions[i, :]
		y = snap.positions[1 + ((i+1) % 3), :]

		h = histogram2d(x, y, (bins, bins))
		h_min = 1e-10

		h_max = log10(maximum(h.values))
		ax = Axis(fig[1,i])
		contour!(midpoints(h.xbins), midpoints(h.ybins), log10.(h.values .+ h_min), levels=LinRange(h_max-5, h_max, 40), linestyle=:solid)

		hidedecorations!(ax)
	end

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	fig
end


# ╔═╡ b642d4d5-7291-4c1f-a1eb-5616b4869229
plot_isocontours(snap)

# ╔═╡ 3f32a9dc-8937-404d-9dbf-aa4f0d50b9a5
plot_isocontours(snap_2)

# ╔═╡ 565176e9-cf17-4752-958c-d3c5b5efadc5
prof = LilGuys.DensityProfile(snap)

# ╔═╡ 1a263d8a-7b22-4779-bdb0-29bee53290f0
prof_f = LilGuys.DensityProfile(snap_2)

# ╔═╡ 5f703522-29b8-490f-a5f8-d875988d5fca
let 
	f = lines(prof.log_r, log10.(prof.rho))

	lines!(prof_f.log_r, log10.(prof_f.rho))
	ylims!(-10, 2)

	f
end

# ╔═╡ Cell order:
# ╠═88d39cc6-95a4-11f0-1b8a-e519835822e7
# ╠═83fea91d-83fe-436a-8e4e-f21a61468562
# ╠═8a457a1f-05c9-4cf5-8244-18a364782b74
# ╠═6278077a-b027-4dc9-8304-197b049346e6
# ╠═da071214-ff07-438b-a677-b54687beb8ae
# ╠═67d1b39c-c314-4313-8e5a-5073d9b6a11d
# ╠═75db6089-d5a2-419a-8354-d2a996fc19d4
# ╠═280caba5-b73d-43f6-9792-68166000bb8d
# ╠═b642d4d5-7291-4c1f-a1eb-5616b4869229
# ╠═3f32a9dc-8937-404d-9dbf-aa4f0d50b9a5
# ╠═565176e9-cf17-4752-958c-d3c5b5efadc5
# ╠═1a263d8a-7b22-4779-bdb0-29bee53290f0
# ╠═5f703522-29b8-490f-a5f8-d875988d5fca
