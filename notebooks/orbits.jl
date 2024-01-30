### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1f568548-84e8-11ee-0abf-cd4651e92589
begin
	import Pkg
	Pkg.activate()
	import LilGuys as lguys
	using Plots; plotly()
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
	r_c = LinRange(0, 20, 100)
	vs = V_circ(Φ_mw, r_c)
	plot(r_c, vs)
	xlabel!("r/kpc")
	ylabel!("Vcirc")
end

# ╔═╡ ffa7eb1e-b515-4bf0-9f68-c8ea1c5191c8
begin
c = 9.545
Mhalo = 115 * 1/(log(1+c) - (c)/(1+c))
end

# ╔═╡ fc6cc98e-58c9-46e1-8a77-874cb078614d
σ_slider = Slider(-3:0.1:3, default=0) # a 3 σ slider for standard normal distributions...

# ╔═╡ 27c9d926-06ff-45b1-895e-a93e12749452
@bind dpm_dec σ_slider

# ╔═╡ a8af2878-70df-4eca-a535-e0b3dd720a81
@bind dpm_ra σ_slider

# ╔═╡ ea6a411f-428e-4aed-98a8-c2fdf783cfca
@bind ddistance σ_slider

# ╔═╡ a1e19b57-4cfc-4081-993d-8ad6ca7426de
@bind dradial_velocity σ_slider

# ╔═╡ aeaae20c-4ee7-4f66-8b1f-52d3afdda205
begin
	distance = 86 + 6* ddistance
	pm_ra =(0.099 + 0.002*dpm_ra)
	pm_dec = (-0.160 + 0.002*dpm_dec)
	radial_velocity = (111.4 +0.37 * dradial_velocity)
	ra=15.03917
	dec=-33.70917
	obs = lguys.Observation(ra, dec, distance, pm_ra, pm_dec, radial_velocity)

end

# ╔═╡ 74eb54e2-77b3-4fed-9408-65547b7bb295
begin
	phase = lguys.transform(lguys.Galactocentric, obs)
	phase2 = lguys.Galactocentric(phase.position, -phase.velocity)
end

# ╔═╡ 98bb6eb8-ead1-456d-b906-612763a2ddd6
md"""
- pm-ra = $(pm_ra)
- pm-dec = $(pm_dec)
- rv = $(radial_velocity)
- distance = $(distance)
"""

# ╔═╡ 05ad80bc-7892-4c0c-b843-f6b84da2fcc3
import LinearAlgebra: transpose!

# ╔═╡ 0cbbba9a-8009-4e88-a19d-41fb121ca994
begin
	t, positions, velocities, phi = calc_orbit(phase2)
	t2, positions2, velocities2, phi2 = calc_orbit(phase)
	t = [a[1] for a in t]
end

# ╔═╡ 3aee96d3-e86f-4b1a-b228-68cfc2c48729
begin 
	r = lguys.calc_r(positions)
	v = lguys.calc_r(velocities)

	r2 = lguys.calc_r(positions2)
	v2 = lguys.calc_r(velocities2)
end

# ╔═╡ ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
md"""
The current pericenter is $(minimum(r)) and apocenter $(maximum(r))
"""

# ╔═╡ 9e8b5a2f-8033-4123-a0fc-78e43c78603b
begin
	plot(t, r)

	plot!(-t2.value, r2)
	hline!([minimum(r), maximum(r)])
end

# ╔═╡ 0de9032e-e2d3-4cca-925f-c7686f3487e3

begin
plot(r, v)
plot!(r2, v2)
end

# ╔═╡ 73b03786-b61c-4a13-8260-b85f9ad1857b
lguys.scatter_xyz(positions[:, 1:2], reshape(phase.position, 3, 1))

# ╔═╡ cbc81993-5d73-4831-91e0-ca7039293f4a
lguys.scatter_xyz(velocities[:, 1:2], reshape(phase.velocity, 3, 1))

# ╔═╡ f742c296-b2c6-4fd8-8d87-a0ad84b3b525
begin
	df = DataFrame(t = t , x=positions[1, :], y=positions[2, :], z=positions[3, :], 
		vx=velocities[1, :], vy=velocities[2, :], vz=velocities[3, :], phi=phi)
	CSV.write("sculptor_orbits.csv", df)
end

# ╔═╡ f1d9d39b-4d93-4359-9742-59ef67a531e3
lguys.plot_xyz(positions)

# ╔═╡ Cell order:
# ╟─3a360b8d-e830-47ab-bfbc-ce147b15d09d
# ╠═1f568548-84e8-11ee-0abf-cd4651e92589
# ╠═553b9e95-2748-4651-b108-0e177b286322
# ╠═6607d695-788a-4667-91ec-b1d290fb9ea5
# ╠═ffa7eb1e-b515-4bf0-9f68-c8ea1c5191c8
# ╠═fc6cc98e-58c9-46e1-8a77-874cb078614d
# ╠═aeaae20c-4ee7-4f66-8b1f-52d3afdda205
# ╠═74eb54e2-77b3-4fed-9408-65547b7bb295
# ╟─98bb6eb8-ead1-456d-b906-612763a2ddd6
# ╠═27c9d926-06ff-45b1-895e-a93e12749452
# ╠═a8af2878-70df-4eca-a535-e0b3dd720a81
# ╠═ea6a411f-428e-4aed-98a8-c2fdf783cfca
# ╠═a1e19b57-4cfc-4081-993d-8ad6ca7426de
# ╟─ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
# ╠═05ad80bc-7892-4c0c-b843-f6b84da2fcc3
# ╠═0cbbba9a-8009-4e88-a19d-41fb121ca994
# ╠═3aee96d3-e86f-4b1a-b228-68cfc2c48729
# ╠═9e8b5a2f-8033-4123-a0fc-78e43c78603b
# ╠═0de9032e-e2d3-4cca-925f-c7686f3487e3
# ╠═73b03786-b61c-4a13-8260-b85f9ad1857b
# ╠═cbc81993-5d73-4831-91e0-ca7039293f4a
# ╠═f742c296-b2c6-4fd8-8d87-a0ad84b3b525
# ╠═f1d9d39b-4d93-4359-9742-59ef67a531e3
