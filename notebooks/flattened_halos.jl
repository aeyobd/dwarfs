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
snap = Snapshot("../agama/halos_spheroidal/sculptor_1e6.hdf5")

# ╔═╡ 6278077a-b027-4dc9-8304-197b049346e6
CairoMakie.activate!(type=:png)

# ╔═╡ da071214-ff07-438b-a677-b54687beb8ae
LilGuys.plot_xyz(snap.positions, plot=:scatter, markersize=1, alpha=0.1)

# ╔═╡ 75db6089-d5a2-419a-8354-d2a996fc19d4
import DensityEstimators: histogram2d

# ╔═╡ 280caba5-b73d-43f6-9792-68166000bb8d
let
	fig = Figure()

	bins = -10:0.5:10
	
	for i in 1:3
		x = snap.positions[i, :]
		y = snap.positions[1 + ((i+1) % 3), :]

		h = histogram2d(x, y, (bins, bins))
		h_min = 1e-10

		h_max = log10(maximum(h.values))
		ax = Axis(fig[1,i])
		contour!(midpoints(h.xbins), midpoints(h.ybins), log10.(h.values .+ h_min), levels=LinRange(h_max-5, h_max, 20), linestyle=:solid)

		hidedecorations!(ax)
	end

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	fig
end


# ╔═╡ 565176e9-cf17-4752-958c-d3c5b5efadc5
prof = LilGuys.DensityProfile(snap)

# ╔═╡ 5f703522-29b8-490f-a5f8-d875988d5fca
let 
	f = lines(prof.log_r, log10.(prof.rho))

	ylims!(-10, 2)

	f
end

# ╔═╡ Cell order:
# ╠═88d39cc6-95a4-11f0-1b8a-e519835822e7
# ╠═83fea91d-83fe-436a-8e4e-f21a61468562
# ╠═6278077a-b027-4dc9-8304-197b049346e6
# ╠═da071214-ff07-438b-a677-b54687beb8ae
# ╠═75db6089-d5a2-419a-8354-d2a996fc19d4
# ╠═280caba5-b73d-43f6-9792-68166000bb8d
# ╠═565176e9-cf17-4752-958c-d3c5b5efadc5
# ╠═5f703522-29b8-490f-a5f8-d875988d5fca
