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

	u = lguys.u
	gp = lguys.galpy.potential

	using CSV 
	using DataFrames
end

# ╔═╡ 3a360b8d-e830-47ab-bfbc-ce147b15d09d
md"""
This notebook uses Galpy to compute some example orbits in a milky way potential. This is useful for rapidly development and exploration of orbits.
"""

# ╔═╡ 6607d695-788a-4667-91ec-b1d290fb9ea5
begin
	Φ_mw = lguys.MWPotential()
	r_c = LinRange(0, 20, 100)
	vs = lguys.V_circ(Φ_mw, r_c)
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

# ╔═╡ f6e46008-b849-4f2b-9296-ec28120e8866
lguys.PhasePoint(zeros(3), ones(3)).z

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
phase = lguys.to_galcen(obs)
phase2 = lguys.PhasePoint(phase.position, -phase.velocity)

obs_r = lguys.to_sky(phase2)
end

# ╔═╡ 0c1a14da-5d01-4937-9df4-91b4bb27f78f
phase2

# ╔═╡ ca0d411b-5dbc-4fa1-b4da-c6107a8088fa
lguys.to_galcen(obs_r)

# ╔═╡ 98bb6eb8-ead1-456d-b906-612763a2ddd6
md"""
- pm-ra = $(pm_ra)
- pm-dec = $(pm_dec)
- rv = $(radial_velocity)
- distance = $(distance)
"""

# ╔═╡ cb724c3b-8ad9-466c-be2e-610ad5c6251a
phase.velocity

# ╔═╡ 0cbbba9a-8009-4e88-a19d-41fb121ca994
begin
	t, x, y, z, vx, vy, vz, phi = lguys.calc_orbit(phase)
	t1, x1, y1, z1, vx1, vy1, vz1, phi1 = lguys.calc_orbit(phase2)

	r = @. sqrt(x^2 + y^2 + z^2)
	v = @. sqrt(vx^2 + vy^2 + vz^2)

	r1 = @. sqrt(x1^2 + y1^2 + z1^2)
	v1 = @. sqrt(vx1^2 + vy1^2 + vz1^2)

end

# ╔═╡ ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
md"""
The current pericenter is $(minimum(r)) and apocenter $(maximum(r))
"""

# ╔═╡ 9e8b5a2f-8033-4123-a0fc-78e43c78603b
begin
	plot(t.value, r)

	plot!(-t1.value, r1)
	hline!([minimum(r), maximum(r)])
end

# ╔═╡ 41d84183-6a2e-4617-864c-48a8cc49965a
[plot(x, y), plot(x, z), plot(y, z)]

# ╔═╡ 84102b3f-b7cf-445f-a7a7-6cc71638ac90
[plot(vx, vy), plot(vx, vz), plot(vy, vz)]

# ╔═╡ 0de9032e-e2d3-4cca-925f-c7686f3487e3

begin
plot(r, v)
plot!(r1, v1)
end

# ╔═╡ 8024cb70-2830-4b38-b7d1-ced4a0e49742
v

# ╔═╡ a4fcb1e7-012c-43df-a23a-086b72aed39d
v1

# ╔═╡ 52fce3f3-ff1b-4908-a676-71747cc64cf2
begin
plot(x, z)
plot!(x1, z1)
end

# ╔═╡ f742c296-b2c6-4fd8-8d87-a0ad84b3b525
begin
	df = DataFrame(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, phi=phi)
	CSV.write("sculptor_orbits.csv", df)
end

# ╔═╡ Cell order:
# ╟─3a360b8d-e830-47ab-bfbc-ce147b15d09d
# ╠═1f568548-84e8-11ee-0abf-cd4651e92589
# ╠═6607d695-788a-4667-91ec-b1d290fb9ea5
# ╠═ffa7eb1e-b515-4bf0-9f68-c8ea1c5191c8
# ╠═fc6cc98e-58c9-46e1-8a77-874cb078614d
# ╠═0c1a14da-5d01-4937-9df4-91b4bb27f78f
# ╠═ca0d411b-5dbc-4fa1-b4da-c6107a8088fa
# ╠═f6e46008-b849-4f2b-9296-ec28120e8866
# ╠═aeaae20c-4ee7-4f66-8b1f-52d3afdda205
# ╠═74eb54e2-77b3-4fed-9408-65547b7bb295
# ╟─98bb6eb8-ead1-456d-b906-612763a2ddd6
# ╠═27c9d926-06ff-45b1-895e-a93e12749452
# ╠═a8af2878-70df-4eca-a535-e0b3dd720a81
# ╠═ea6a411f-428e-4aed-98a8-c2fdf783cfca
# ╠═a1e19b57-4cfc-4081-993d-8ad6ca7426de
# ╟─ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
# ╠═cb724c3b-8ad9-466c-be2e-610ad5c6251a
# ╠═0cbbba9a-8009-4e88-a19d-41fb121ca994
# ╠═9e8b5a2f-8033-4123-a0fc-78e43c78603b
# ╠═41d84183-6a2e-4617-864c-48a8cc49965a
# ╠═84102b3f-b7cf-445f-a7a7-6cc71638ac90
# ╠═0de9032e-e2d3-4cca-925f-c7686f3487e3
# ╠═8024cb70-2830-4b38-b7d1-ced4a0e49742
# ╠═a4fcb1e7-012c-43df-a23a-086b72aed39d
# ╠═52fce3f3-ff1b-4908-a676-71747cc64cf2
# ╠═f742c296-b2c6-4fd8-8d87-a0ad84b3b525
