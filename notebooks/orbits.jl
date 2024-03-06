### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 1f568548-84e8-11ee-0abf-cd4651e92589
begin
	import Pkg
	Pkg.activate()
	import LilGuys as lguys
	using Plots; #plotly()
	using PlutoUI

	using CSV 
	using DataFrames
end

# ╔═╡ 553b9e95-2748-4651-b108-0e177b286322
include("/cosma/home/durham/dc-boye1/dwarfs/scripts/orbit.jl")

# ╔═╡ 3a360b8d-e830-47ab-bfbc-ce147b15d09d
md"""
This notebook uses Galpy to compute some example orbits in a milky way potential. This is useful for rapidly development and exploration of orbits.
"""

# ╔═╡ 6607d695-788a-4667-91ec-b1d290fb9ea5
begin
	Φ_mw = MWPotential()
	r_c = LinRange(0, 50, 100)
	vs = V_circ(Φ_mw, r_c)
	plot(r_c, vs)
	xlabel!("r/kpc")
	ylabel!("Vcirc")
	savefig("v_circ.pdf")
end

# ╔═╡ 08fde7d1-b529-4a37-a77f-c27ea3bdc30d
begin 
	ddist = 0 #6.151974
	dpm_ra = 0# -0.0123
	dpm_dec = 0 #-0.04785
	drv = 0# 1.8
end

# ╔═╡ aeaae20c-4ee7-4f66-8b1f-52d3afdda205
begin
	distance = 86 +  ddist
	pm_ra = 0.099 + dpm_ra
	pm_dec = -0.160 + dpm_dec
	radial_velocity = 111.4 + drv
	ra = 15.03917
	dec = -33.70917
	obs = lguys.Observation(ra, dec, distance, pm_ra, pm_dec, radial_velocity)

end

# ╔═╡ 74eb54e2-77b3-4fed-9408-65547b7bb295
begin
	phase = lguys.transform(lguys.Galactocentric, obs)
	phase2 = lguys.Galactocentric(phase.position, -phase.velocity)
end

# ╔═╡ 0cbbba9a-8009-4e88-a19d-41fb121ca994
begin
	t, positions, velocities, phi = calc_orbit(phase2)
	t = t.value
	positions = reverse(positions, dims=2)
	velocities = -reverse(velocities, dims=2)
end

# ╔═╡ 3aee96d3-e86f-4b1a-b228-68cfc2c48729
begin 
	r = lguys.calc_r(positions)
	v = lguys.calc_r(velocities)
end

# ╔═╡ ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
md"""
The current pericenter is $(minimum(r)) and apocenter $(maximum(r))
"""

# ╔═╡ 9e8b5a2f-8033-4123-a0fc-78e43c78603b
begin
	plot(t, r)
	hline!([minimum(r), maximum(r)])
end

# ╔═╡ 50a2c268-d886-4d8b-b79e-6a204a2f0d98
phase_i = lguys.Galactocentric(positions[:, 1], velocities[:, 1])

# ╔═╡ f742c296-b2c6-4fd8-8d87-a0ad84b3b525
begin
	df = DataFrame(t = t , x=positions[1, :], y=positions[2, :], z=positions[3, :], 
		vx=velocities[1, :], vy=velocities[2, :], vz=velocities[3, :], phi=phi)
	CSV.write("sculptor_orbits.csv", df)
end

# ╔═╡ f1d9d39b-4d93-4359-9742-59ef67a531e3
lguys.plot_xyz(positions)

# ╔═╡ 8d23a56a-f90e-4b62-94d4-b4d3609542da
lguys.plot_xyz(velocities)

# ╔═╡ 0de9032e-e2d3-4cca-925f-c7686f3487e3

begin
plot(r, v)
end

# ╔═╡ 2b4696b8-05f6-489b-b8fb-49f5529f4abc
positions[:, 1]

# ╔═╡ dc63b947-1d17-44b5-bec4-3a0b6218e176
positions[:, end]

# ╔═╡ 899c01ec-4d13-4b32-9bcf-390e1dcd72dd
velocities[:, 1]  ./ lguys.V0

# ╔═╡ f804f503-afcb-49b9-a1c2-aad84089290f
velocities[:, end]  ./ lguys.V0

# ╔═╡ Cell order:
# ╟─3a360b8d-e830-47ab-bfbc-ce147b15d09d
# ╠═1f568548-84e8-11ee-0abf-cd4651e92589
# ╠═553b9e95-2748-4651-b108-0e177b286322
# ╟─6607d695-788a-4667-91ec-b1d290fb9ea5
# ╠═08fde7d1-b529-4a37-a77f-c27ea3bdc30d
# ╠═aeaae20c-4ee7-4f66-8b1f-52d3afdda205
# ╠═74eb54e2-77b3-4fed-9408-65547b7bb295
# ╟─ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
# ╠═0cbbba9a-8009-4e88-a19d-41fb121ca994
# ╠═3aee96d3-e86f-4b1a-b228-68cfc2c48729
# ╠═9e8b5a2f-8033-4123-a0fc-78e43c78603b
# ╠═50a2c268-d886-4d8b-b79e-6a204a2f0d98
# ╠═f742c296-b2c6-4fd8-8d87-a0ad84b3b525
# ╠═f1d9d39b-4d93-4359-9742-59ef67a531e3
# ╠═8d23a56a-f90e-4b62-94d4-b4d3609542da
# ╠═0de9032e-e2d3-4cca-925f-c7686f3487e3
# ╠═2b4696b8-05f6-489b-b8fb-49f5529f4abc
# ╠═dc63b947-1d17-44b5-bec4-3a0b6218e176
# ╠═899c01ec-4d13-4b32-9bcf-390e1dcd72dd
# ╠═f804f503-afcb-49b9-a1c2-aad84089290f
