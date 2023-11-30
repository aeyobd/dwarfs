### A Pluto.jl notebook ###
# v0.19.32

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
	using LilGuys
	using PyCall
	using Plots
	using PlutoUI
end

# ╔═╡ 7625daf9-5e45-41a1-aebd-96656bf0e86f
begin
	galpy = pyimport("galpy")
	pot = pyimport("galpy.potential")
	conversion = pyimport("galpy.util.conversion")
	u = pyimport("astropy.units")
	orbit = pyimport("galpy.orbit")
	coords = pyimport("astropy.coordinates")
end

# ╔═╡ d30862dc-00ca-4ccf-a602-0b67147ffd1f
begin
	M0 = 1e10*u.M_sun
	R0 = 1*u.kpc
	V0 = 207.4*u.km/u.s
end

# ╔═╡ 95071915-95ec-4d56-9fbc-9ab6ef9f4ed9
md"""
| Component  | Values                    |      |
| ---------- | ------------------------- | ---- |
| thin disk  | M=5.9, a=3.9, b=0.31      |      |
| thick disk | M=2, a=4.4, b=0.92        |      |
| bulge      | M = 2.1, a = 1.3          |      |
| halo       | Mvir=115, r=20.2, c=9.545 |      |
"""

# ╔═╡ ffa7eb1e-b515-4bf0-9f68-c8ea1c5191c8
begin
c = 9.545
Mhalo = 115 * 1/(log(1+c) - (c)/(1+c))
end

# ╔═╡ 7507f841-51a7-495b-a19a-d3da6751b140
begin
	hp = pot.HernquistPotential(amp=2.1M0, a=1.3R0)
	np = pot.NFWPotential(amp=Mhalo*M0, a=20R0)
	thin_p = pot.MiyamotoNagaiPotential(amp=5.9M0, a=3.9R0, b=0.31R0)
	thick_p = pot.MiyamotoNagaiPotential(amp=2M0, a=4.4R0, b=0.92R0)
	tp = hp + np + thick_p + thin_p
end

# ╔═╡ 6607d695-788a-4667-91ec-b1d290fb9ea5
begin
	r_c = LinRange(0, 20, 100)
	v = pot.vcirc(tp, r_c * u.kpc)
	plot(r_c, v)
	xlabel!("r/kpc")
	ylabel!("Vcirc")
end

# ╔═╡ fc6cc98e-58c9-46e1-8a77-874cb078614d
σ_slider = Slider(-3:0.1:3, default=0)

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
	obs = coords.SkyCoord(ra=15.03917*u.deg, dec=-33.70917*u.deg,
		distance= distance*u.kpc,
		pm_ra_cosdec =pm_ra*u.mas/u.yr,
		pm_dec = pm_dec*u.mas/u.yr,
		radial_velocity=radial_velocity* u.km/u.s)

	pp = obs.transform_to(LilGuys.galcen_frame)

end

# ╔═╡ d4616277-e1d8-40ab-8d3d-c5a493ab5615
begin
	o = orbit.Orbit(obs)
end

# ╔═╡ 98bb6eb8-ead1-456d-b906-612763a2ddd6
md"""
- pm-ra = $(pm_ra)
- pm-dec = $(pm_dec)
- rv = $(radial_velocity)
- distance = $(distance)
"""

# ╔═╡ 8fb855e6-246f-44dc-a862-b28409e5883e
md"""
- x = $(pp.x[])
- y = $(pp.y[])
- z = $(pp.z[])
- vx = $(pp.v_x[])
- vy = $(pp.v_y[])
- vz = $(pp.v_z[])

"""

# ╔═╡ 1fbe302e-4f2f-47ae-963c-d94baed9e30c
1

# ╔═╡ 0cbbba9a-8009-4e88-a19d-41fb121ca994
begin
	ts = LinRange(0, 10, 1000)
o.integrate(ts * u.Gyr, tp, method="odeint")
x = o.x(ts * u.Gyr)
y = o.y(ts * u.Gyr)
z = o.z(ts * u.Gyr)
end

# ╔═╡ 9e8b5a2f-8033-4123-a0fc-78e43c78603b
begin
	r = @. sqrt(x^2 + y^2 + z^2)
	plot(ts, r)
	hline!([minimum(r), maximum(r)])
end

# ╔═╡ ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
md"""
The current pericenter is $(minimum(r)) and apocenter $(maximum(r))
"""

# ╔═╡ 41d84183-6a2e-4617-864c-48a8cc49965a
plot(x, y)

# ╔═╡ 56d10644-402c-4cb2-bf04-eee362820fbb
plot(x, z)

# ╔═╡ Cell order:
# ╠═1f568548-84e8-11ee-0abf-cd4651e92589
# ╠═7625daf9-5e45-41a1-aebd-96656bf0e86f
# ╠═d30862dc-00ca-4ccf-a602-0b67147ffd1f
# ╠═7507f841-51a7-495b-a19a-d3da6751b140
# ╟─95071915-95ec-4d56-9fbc-9ab6ef9f4ed9
# ╠═6607d695-788a-4667-91ec-b1d290fb9ea5
# ╠═ffa7eb1e-b515-4bf0-9f68-c8ea1c5191c8
# ╠═fc6cc98e-58c9-46e1-8a77-874cb078614d
# ╠═d4616277-e1d8-40ab-8d3d-c5a493ab5615
# ╠═27c9d926-06ff-45b1-895e-a93e12749452
# ╠═a8af2878-70df-4eca-a535-e0b3dd720a81
# ╠═ea6a411f-428e-4aed-98a8-c2fdf783cfca
# ╠═a1e19b57-4cfc-4081-993d-8ad6ca7426de
# ╟─98bb6eb8-ead1-456d-b906-612763a2ddd6
# ╟─8fb855e6-246f-44dc-a862-b28409e5883e
# ╠═aeaae20c-4ee7-4f66-8b1f-52d3afdda205
# ╠═1fbe302e-4f2f-47ae-963c-d94baed9e30c
# ╟─ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
# ╠═0cbbba9a-8009-4e88-a19d-41fb121ca994
# ╠═9e8b5a2f-8033-4123-a0fc-78e43c78603b
# ╠═41d84183-6a2e-4617-864c-48a8cc49965a
# ╠═56d10644-402c-4cb2-bf04-eee362820fbb
