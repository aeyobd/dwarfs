### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ e644b842-38b5-11ef-2c88-6bbbba8f7e2f
begin
	using Pkg; Pkg.activate()

	import LilGuys as lguys
	using CairoMakie
	using Arya
end

# ╔═╡ 09defccf-1bba-49a9-83e2-a752ab12a118
md"""
# Test Profile

Tests the creation of an N-body model from either zeno or Agama
"""

# ╔═╡ 8bb36e00-a92b-4389-abfe-32866090fffe
snap = lguys.Snapshot("test.hdf5")

# ╔═╡ 7176f499-3ac5-4cce-8451-e6d383bcaa71
r, rho = lguys.calc_ρ_hist(lguys.calc_r(snap), 100, weights=snap.masses)

# ╔═╡ 5ed8fe7b-77da-48f2-9157-5a4c14b011fb
sum(snap.masses[lguys.calc_r(snap) .< 1])

# ╔═╡ 1254ea0a-7c80-4123-9eed-de05d588e8e4
halo = lguys.NFW()

# ╔═╡ 3bf23b50-fce4-4da3-b33d-204428f20cfc
let
	fig, ax = FigAxis(
		xscale=log10,
		yscale=log10
	)

	scatter!(midpoints(r), rho)

	r_model = 10 .^ LinRange(-1, 3, 1000)
	ρ_model = lguys.calc_ρ.(halo, r_model)

	lines!(r_model, ρ_model )
	fig
end

# ╔═╡ 61b7f562-f361-4e36-bf38-ce69a582cd46
snap.velocities

# ╔═╡ 3b18f5b3-c1e6-4719-b489-b9004eb1d1a8
[ [1,2]j [3, 4] [5,6]]

# ╔═╡ Cell order:
# ╠═09defccf-1bba-49a9-83e2-a752ab12a118
# ╠═e644b842-38b5-11ef-2c88-6bbbba8f7e2f
# ╠═8bb36e00-a92b-4389-abfe-32866090fffe
# ╠═7176f499-3ac5-4cce-8451-e6d383bcaa71
# ╠═5ed8fe7b-77da-48f2-9157-5a4c14b011fb
# ╠═1254ea0a-7c80-4123-9eed-de05d588e8e4
# ╠═3bf23b50-fce4-4da3-b33d-204428f20cfc
# ╠═61b7f562-f361-4e36-bf38-ce69a582cd46
# ╠═3b18f5b3-c1e6-4719-b489-b9004eb1d1a8
