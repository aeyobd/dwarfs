### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 1f568548-84e8-11ee-0abf-cd4651e92589
begin
	import Pkg
	Pkg.activate()

	using CairoMakie

	import LilGuys as lguys
	import LilGuys.Plots as LP

	using PlutoUI

	using Arya
	using PythonCall
end

# ╔═╡ 3a360b8d-e830-47ab-bfbc-ce147b15d09d
md"""
# Orbits jl

This notebook uses Agama to compute some example orbits in a milky way potential. The notebook is somewhat interactive which is nice for exploring orbit families and such

"""

# ╔═╡ a2467e1c-8976-4812-9c23-6940c5ebb9c7
agama = pyimport("agama")

# ╔═╡ b1502d98-23df-4e9b-bbbd-077951bb1b78
pot_file = "/astro/dboyea/dwarfs/agama/potentials/EP2020.ini"

# ╔═╡ 9cbe93e1-984f-4454-956f-c8968eb7846b
Φ_mw = agama.Potential(pot_file)

# ╔═╡ 70344aa2-4c70-4643-a131-e39a10eb0540
function v_circ(Φ, r)
	M = pyconvert(Vector{Float64}, Φ.enclosedMass(r))
	
	return sqrt.(M ./ r)
end

# ╔═╡ 17cc9a81-9a32-489d-8393-b83b4a97f0b2
Φ_mw.enclosedMass(1:10)

# ╔═╡ 6607d695-788a-4667-91ec-b1d290fb9ea5
let
	fig, ax = FigAxis(
		xlabel = "r / kpc",
		ylabel = L"circular velocity / km\,s$^{-1}$",
	)
	
	r_c = LinRange(0, 100, 1000)
	vs = v_circ(Φ_mw, r_c)
	lines!(r_c, vs * lguys.V2KMS)

	fig

end

# ╔═╡ 7c8ab08f-a910-48f2-a1aa-16c9d244b677
md"""
## Inputs
"""

# ╔═╡ fb541fbf-e030-455f-be7c-6ae36a6af906
import TOML

# ╔═╡ f6597677-fef5-45fe-b26f-cf1ea7d6df90
obs_props = TOML.parsefile("/astro/dboyea/dwarfs/observations/sculptor/observed_properties.toml")

# ╔═╡ 9814b558-0fd6-43ef-b145-4d94136bada7
ddist = 0 * obs_props["distance_err"]

# ╔═╡ af1bb0be-59e4-40b2-af88-f7428a125c24
dpm_ra = 0 * obs_props["pmra_err"]

# ╔═╡ eb2e8a83-f2c2-4ac2-b4e4-3c78e8ed21e8
dpm_dec = 0 * obs_props["pmdec_err"]

# ╔═╡ 741ee80a-fac4-4494-a546-26534ba0719e
drv = 0 * obs_props["radial_velocity_err"]

# ╔═╡ aeaae20c-4ee7-4f66-8b1f-52d3afdda205
obs = lguys.ICRS(
	distance = obs_props["distance"] +  ddist,
	pmra = obs_props["pmra"] + dpm_ra,
	pmdec = obs_props["pmdec"] + dpm_dec,
	radial_velocity = obs_props["radial_velocity"] + drv,
	ra = obs_props["ra"],
	dec = obs_props["dec"],

)

# ╔═╡ 15c9c09a-a628-4b50-a8c1-e55e7d4aa1af
md"""
## Plots
"""

# ╔═╡ 7e12d914-d858-47f6-9bc1-415d216178af
md"""
## Calculations
"""

# ╔═╡ f2764f05-a333-4947-94b7-134073f392ec
phase = lguys.transform(lguys.Galactocentric, obs)

# ╔═╡ f7d0cb1a-05e9-45ca-bac3-cd7e418bf68b
posvel = [phase.x, phase.y, phase.z, 
	((phase.v_x, phase.v_y, phase.v_z) ./ lguys.V2KMS )...]

# ╔═╡ 59ccd356-b519-4d85-a15a-c805d834f821
orbit = agama.orbit(potential=Φ_mw, ic=posvel, time=10/lguys.T2GYR, trajsize=1_000)

# ╔═╡ 61296a89-2518-4f2a-95ec-6c969aebaa2a
t = pyconvert(Vector{Float64}, orbit[0])

# ╔═╡ 0cbbba9a-8009-4e88-a19d-41fb121ca994
posvel_orbit = pyconvert(Matrix{Float64}, orbit[1])'

# ╔═╡ 29859b57-fce5-4517-9be4-b2b8d92b1858
positions = posvel_orbit[1:3, :]

# ╔═╡ 78a1cf0c-b0c5-45b7-ad69-75c7d49c1d6a
r = lguys.calc_r(positions)

# ╔═╡ e53d69e3-58f8-4fba-a419-fb850aa1dc10
peri = minimum(r)

# ╔═╡ fbf20220-c502-41cf-8f92-592e33a023d6
apo = maximum(r)

# ╔═╡ 9e8b5a2f-8033-4123-a0fc-78e43c78603b
let
	fig, ax = FigAxis(
		xlabel="time / Gyr",
		ylabel = "radius / kpc"
	)
	
	lines!(t * lguys.T2GYR, r)
	
	hlines!([minimum(r), maximum(r)], color=COLORS[2])
	fig
end

# ╔═╡ ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
md"""
The current pericenter is $(minimum(r)) and apocenter $(maximum(r))
"""

# ╔═╡ f1d9d39b-4d93-4359-9742-59ef67a531e3
LP.plot_xyz(positions)

# ╔═╡ 9e08147d-42cc-4556-9360-f639ba00e28b
velocities = posvel_orbit[4:6, :]

# ╔═╡ d6267c97-2f4d-4bff-baa5-3784ad3416ee
v = lguys.calc_r(velocities)

# ╔═╡ 0de9032e-e2d3-4cca-925f-c7686f3487e3

let
	fig, ax = FigAxis(
		xlabel = "r / kpc",
		ylabel = L"v / km\,s$^{-1}$"
	)
	scatter!(r, v)
	fig
end

# ╔═╡ 50a2c268-d886-4d8b-b79e-6a204a2f0d98
phase_i = lguys.Galactocentric(positions[:, 1], velocities[:, 1])

# ╔═╡ 8d23a56a-f90e-4b62-94d4-b4d3609542da
LP.plot_xyz(velocities)

# ╔═╡ f742c296-b2c6-4fd8-8d87-a0ad84b3b525
#df = DataFrame(t = t , x=positions[1, :], y=positions[2, :], z=positions[3, :], vx=velocities[1, :], vy=velocities[2, :], vz=velocities[3, :])

# ╔═╡ Cell order:
# ╟─3a360b8d-e830-47ab-bfbc-ce147b15d09d
# ╠═1f568548-84e8-11ee-0abf-cd4651e92589
# ╠═a2467e1c-8976-4812-9c23-6940c5ebb9c7
# ╠═b1502d98-23df-4e9b-bbbd-077951bb1b78
# ╠═9cbe93e1-984f-4454-956f-c8968eb7846b
# ╠═70344aa2-4c70-4643-a131-e39a10eb0540
# ╠═17cc9a81-9a32-489d-8393-b83b4a97f0b2
# ╟─6607d695-788a-4667-91ec-b1d290fb9ea5
# ╟─7c8ab08f-a910-48f2-a1aa-16c9d244b677
# ╠═fb541fbf-e030-455f-be7c-6ae36a6af906
# ╠═f6597677-fef5-45fe-b26f-cf1ea7d6df90
# ╠═aeaae20c-4ee7-4f66-8b1f-52d3afdda205
# ╠═9814b558-0fd6-43ef-b145-4d94136bada7
# ╠═af1bb0be-59e4-40b2-af88-f7428a125c24
# ╠═eb2e8a83-f2c2-4ac2-b4e4-3c78e8ed21e8
# ╠═741ee80a-fac4-4494-a546-26534ba0719e
# ╟─15c9c09a-a628-4b50-a8c1-e55e7d4aa1af
# ╠═e53d69e3-58f8-4fba-a419-fb850aa1dc10
# ╠═fbf20220-c502-41cf-8f92-592e33a023d6
# ╠═78a1cf0c-b0c5-45b7-ad69-75c7d49c1d6a
# ╠═d6267c97-2f4d-4bff-baa5-3784ad3416ee
# ╟─9e8b5a2f-8033-4123-a0fc-78e43c78603b
# ╠═50a2c268-d886-4d8b-b79e-6a204a2f0d98
# ╠═f1d9d39b-4d93-4359-9742-59ef67a531e3
# ╠═8d23a56a-f90e-4b62-94d4-b4d3609542da
# ╠═0de9032e-e2d3-4cca-925f-c7686f3487e3
# ╟─7e12d914-d858-47f6-9bc1-415d216178af
# ╠═f2764f05-a333-4947-94b7-134073f392ec
# ╠═f7d0cb1a-05e9-45ca-bac3-cd7e418bf68b
# ╠═59ccd356-b519-4d85-a15a-c805d834f821
# ╟─ccfc0c3c-008e-4cc4-8f90-1f319b8559d3
# ╠═61296a89-2518-4f2a-95ec-6c969aebaa2a
# ╠═0cbbba9a-8009-4e88-a19d-41fb121ca994
# ╠═29859b57-fce5-4517-9be4-b2b8d92b1858
# ╠═9e08147d-42cc-4556-9360-f639ba00e28b
# ╠═f742c296-b2c6-4fd8-8d87-a0ad84b3b525
