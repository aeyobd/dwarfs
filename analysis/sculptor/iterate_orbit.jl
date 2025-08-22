### A Pluto.jl notebook ###
# v0.20.15

using Markdown
using InteractiveUtils

# ╔═╡ bb303f02-7f11-11f0-241c-a503f150d8aa
begin 
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ a7b71a4d-0779-49e0-9781-b4c03877e2dc
using Measurements

# ╔═╡ 6ac60f86-633d-4429-b3bc-3030b0190447
using OrderedCollections

# ╔═╡ d592d2bc-b3a6-45be-ab14-cc883e76b93e
md"""
# Gradient based orbit adjustments
The idea with this notebook is to adjust an orbit using an emperical gradient to match a specified target.

Specifically, by taking a central orbit and 6 perturbations (in orthogonal directions ideally), we can compute the Jacobian relating the initial and final positions to eachother. Since our goal is to match the final position and time, we take the timestep with the closest to final time as our final position.

Improvements would allow more flexibility in the naming of orbits (currently must all have the same suffixes) and overspecified jacobian so we can solve for a more stable solution.

"""

# ╔═╡ 2fdd2888-9deb-440c-920f-847cee3723b4
import TOML

# ╔═╡ 3585c1bb-3735-45d9-a85b-d7ad392a1fdc
import Agama

# ╔═╡ 2da73c6c-1ba0-4a99-b782-c3daf8c75d46
import LilGuys: positions, velocities

# ╔═╡ 20e97233-e9ad-4172-a6b9-328ae8074969
md"""
# Inputs

The `basename` specifies the base n-body model path which all the companion orbits are only suffixes of. 

The `expected_orbit_file` specifies the orbit which at `time_match` we want the final position to be about the same.


"""

# ╔═╡ 1ad1053e-d674-4ada-8086-017cb77f6af5
basename = "1e5_new_v31_r3.2/L3M11_9Gyr_smallperi.1"

# ╔═╡ 777b4e2e-a9a5-4e92-9f7e-33c0dc07d3fb
expected_orbit_file = "1e6_new_v31_r3.2/L3M11_9Gyr_smallperi/simulation/orbit.csv"

# ╔═╡ 75922729-c1b7-4aa3-91af-5d9515f92711
time_match = -4 / T2GYR

# ╔═╡ 714ab700-5aa8-45ff-94c6-0d762b72504a
md"""
# Data Loading

Each model only requires a `centres.hdf5` file satesfying our centres file format. The central orbit also need the `simulation/agama_potential.ini` to load in the potential and the `simulation/orbit.toml` to load in the orbit properties.

"""

# ╔═╡ dab3ec1d-5c45-4081-b778-7ae6e61bc77f
function load_orbit(extension)
	filename = joinpath(basename * extension, "centres.hdf5")
	return Orbit(filename)
end

# ╔═╡ a7b4ab24-c4a2-47c1-b0a4-4799c222e626
orbit_expected = Orbit(expected_orbit_file)

# ╔═╡ 98bf26ce-a086-44bc-8a94-2968d0ee34aa
pot = Agama.Potential(file = joinpath(basename, "simulation/agama_potential.ini"))

# ╔═╡ 788a12f3-4078-4e06-8027-7fb259d6e7fb
orbits = [
	load_orbit(""),
	load_orbit("_R"),
	load_orbit("_R_m"),
	load_orbit("_z"),
	load_orbit("_z_m"),
	load_orbit("_phi"),
	load_orbit("_phi_m"),
	load_orbit("_v_R"),
	load_orbit("_v_R_m"),
	load_orbit("_v_z"),
	load_orbit("_v_z_m"),
	load_orbit("_v_phi"),
	load_orbit("_v_phi_m"),
]

# ╔═╡ 4b773131-2161-4beb-bf01-fdb63229cadc
md"""
# Plots of the orbits

TODO: label these better
"""

# ╔═╡ 4fd2d8ce-4727-450f-a185-994fb2b184b6
CairoMakie.activate!(px_per_unit=2, type=:png)

# ╔═╡ 0444d257-aa50-4115-b290-af5e0272d752
md"""
# Analysis

We define coordinate frames based on the initial position of the central orbit, and the expected final position. In all cases, we define the frame in terms of a rotated coordinate frame such that `R` points in the same vector as the galaxy's position, `z` is the direction of (instantaneous) angular momentum, and `phi` is the orbital direction. The `U` coordinates are the initial positions transformed into the relative initial coordinate frame, and the `W` are similar except for the final coordinate frame.


"""

# ╔═╡ 50ed57c0-b48a-4ced-a1a4-97c316e3bcfd
import LinearAlgebra: normalize, ×

# ╔═╡ 6a9a0584-a949-49b5-ad74-8417ed935eef
function orbit_unit_vectors(pos_old, vel_old)
	z_hat = normalize(pos_old × vel_old)
	R_hat = normalize(pos_old)
	phi_hat = R_hat × z_hat
	
	return [R_hat z_hat phi_hat]'
end

# ╔═╡ e17a406f-3776-4eec-87bc-0ce275dee59f
idx_f = argmin(abs.(orbit_expected.times .- time_match))

# ╔═╡ 642b4b4e-85eb-4358-8c04-6f9a542fc5df
coord_sys_i = orbit_unit_vectors(orbits[1].positions[:, 1], orbits[1].velocities[:, 1])

# ╔═╡ 1d31c563-e394-4fb9-8251-ed7a9cf8bade
coord_sys_f = orbit_unit_vectors(orbit_expected.positions[:, idx_f], orbit_expected.velocities[:, idx_f])

# ╔═╡ 519d7e7b-377c-4b26-aa19-58088c3cd958
pos_i = coord_sys_i * orbits[1].positions[:, 1]

# ╔═╡ 3f58e0a1-8d81-44bb-be22-c6b3f325cc62
vel_i = coord_sys_i * orbits[1].velocities[:, 1] 

# ╔═╡ 8a05c190-8967-441a-8b20-ef74a0d3a84f
pos_f = coord_sys_f * orbit_expected.positions[:, idx_f]

# ╔═╡ 8855cd39-fbea-4e81-b9db-ac066322d76e
vel_f = coord_sys_f * orbit_expected.velocities[:, idx_f]

# ╔═╡ 91c68e1a-73f7-45a2-b9af-eae7527f405f
v_scale = radii(pos_f) / radii(vel_f)

# ╔═╡ 08e4a7f0-5aee-4326-83f5-082f47e9016a
function initial_projected(orbit)
	return [coord_sys_i * orbit.positions[:, 1] .- pos_i; v_scale * coord_sys_i * orbit.velocities[:, 1] .- v_scale * vel_i]
end

# ╔═╡ 0cbb32f7-fbdb-47f7-a12d-bfdc89d8e37f
function initial_unprojected(u)
	return inv(coord_sys_i) * (u[1:3] + pos_i),   inv(coord_sys_i) * (u[4:6]/v_scale + vel_i)
end

# ╔═╡ c5d3b354-f4b6-4f59-b3c2-c2946efb8b5f
function final_projected(orbit)
	idx = argmin(abs.(orbit.times .- time_match))
	return [coord_sys_f * orbit.positions[:, idx] .- pos_f; v_scale * coord_sys_f * orbit.velocities[:, idx] .- v_scale * vel_f]
end

# ╔═╡ 7812c9d9-92ea-47f8-b836-59249aa6820d
md"""
## Calculations

Given now our initial and final respective coordinate frames, the goal (following vasiliev2024) is to solve for the Jacobian so we can write

``
w \approx J u + w_0
``

where `w_0` is the current `w` vector for the central orbit and ``J`` is the Jacobian.
"""

# ╔═╡ 98ef0760-126a-425a-bcda-fe42079d8b76
O_Ut = hcat(ones(length(orbits)), hcat(initial_projected.(orbits)...)')

# ╔═╡ b484641d-5e48-453a-a140-64844296e588
Wt = hcat(final_projected.(orbits)...)'

# ╔═╡ 7a3ea09a-d1e0-4f32-adde-786f4144cd7a
xi_J = O_Ut \ Wt

# ╔═╡ d2251d57-0c11-45b4-864d-b49b371ddd35
xi = xi_J[1, :]

# ╔═╡ 0f442405-224d-4476-9bec-14127d22dcf5
final_projected(orbits[1])

# ╔═╡ 8b08ce23-bfd9-403b-936d-0f4d867d6731
J = xi_J[2:end, :]'

# ╔═╡ ccc085e6-1571-4f87-9772-791a7750455f
J * initial_projected(orbits[4]) + xi

# ╔═╡ 7378c4c4-2f27-4e3b-9ab0-f8a2bdc492cc
final_projected(orbits[4])

# ╔═╡ d6f7504a-e067-4772-a538-a88efa190e0f
u_next = -inv(J) * xi

# ╔═╡ f4405647-c79d-4a0d-a16f-1572bafa507a
initial_unprojected(inv(J) * (final_projected(orbits[2]) .- xi))

# ╔═╡ 474b504a-b4fa-4ff2-a0d2-70ee3129db05
orbits[2].positions[:, 1], orbits[2].velocities[:, 1]

# ╔═╡ 3939f0d2-545a-46ef-895e-00039eb092cd
pos_next, vel_next = initial_unprojected(u_next)

# ╔═╡ 3bd61ff5-4c9b-49e9-8f55-fcf45c793164
@assert all(initial_unprojected(initial_projected(orbits[2])) .≈ (orbits[2].positions[:, 1], orbits[2].velocities[:, 1]))

# ╔═╡ 83ada54e-e3f2-4f2d-840e-ec7d3d3237b3
orbit_next = LilGuys.agama_orbit(pot, Galactocentric(pos_next, vel_next * V2KMS), timerange=(orbit_expected.times[1], orbit_expected.times[end]), agama_units=Agama.VASILIEV_UNITS)

# ╔═╡ 6e055fa8-47a1-4ba3-9924-e211a13ad63d
LilGuys.plot_xyz(positions(orbit_expected), positions(orbit_next), positions.(orbits)...)

# ╔═╡ df5c9933-51ea-4cc8-a4fe-4a0aeffcba74
let
	fig = Figure()

	ax = Axis(fig[1,1], xlabel="time", ylabel="radius")
	for orbit in orbits
		lines!(orbit.times, radii(orbit))
	end

	lines!(orbit_expected.times, radii(orbit_expected))

	lines!(orbit_next.times, radii(orbit_next), color=:black)

	vlines!(time_match)
	fig

end

# ╔═╡ bc479a76-7d79-4046-80ae-c553fe3af892
md"""
# Data Writing
"""

# ╔═╡ 5fcb209f-a748-41c6-b1bd-90cd7b63a5dc
write(joinpath(basename, "orbit_next.csv"), orbit_next)

# ╔═╡ f7144212-8f81-4396-862e-af255afb2a9b
obs_props = let
	df = TOML.parsefile(joinpath(basename, "simulation/orbit.toml"))
	dm = LilGuys.kpc2dm(df["distance"] ± df["distance_err"])
	df["distance_modulus"] = Measurements.value(dm)
	df["distance_modulus_err"] = Measurements.uncertainty(dm)

	df["ra_err"] = 0
	df["dec_err"] = 0
	df
end

# ╔═╡ 13a5ac47-c607-4e4a-ae53-e558c74a8944
md"""
Orbit

"x_i" 131.857

"y_i" -271.604

"z_i" 74.7348

"v_x_i"-1.19936

"v_y_i" -7.18343

"v_z_i" -8.5247
"""

# ╔═╡ 672c3758-5690-4688-8d2f-56621cc2d985
open(joinpath(basename, "orbit_next.toml"), "w") do f
	d = OrderedDict(
		"lastmodel" => basename,
		"x_i" => orbit_next.positions[1, 1],
		"y_i" => orbit_next.positions[2, 1],
		"z_i" => orbit_next.positions[3, 1],
		"v_x_i" => orbit_next.velocities[1, 1] * V2KMS,
		"v_y_i" => orbit_next.velocities[2, 1] * V2KMS,
		"v_z_i" => orbit_next.velocities[3, 1] * V2KMS,
	)

	for k in ["ra", "dec", "pmra", "pmdec", "distance", "radial_velocity"]
		d[k] = obs_props[k]
		if k * "_err" ∈ keys(obs_props)
			d[k*"_err"] = obs_props[k * "_err"]
		elseif k * "_em" ∈ keys(obs_props)
			d[k*"_em"] = obs_props[k * "_em"]
			d[k*"_ep"] = obs_props[k * "_ep"]
		else
			@warn("error associated with $k missing")
		end
	end
	
	TOML.print(f, d)
	d
end

# ╔═╡ Cell order:
# ╟─d592d2bc-b3a6-45be-ab14-cc883e76b93e
# ╠═bb303f02-7f11-11f0-241c-a503f150d8aa
# ╠═2fdd2888-9deb-440c-920f-847cee3723b4
# ╠═3585c1bb-3735-45d9-a85b-d7ad392a1fdc
# ╠═2da73c6c-1ba0-4a99-b782-c3daf8c75d46
# ╟─20e97233-e9ad-4172-a6b9-328ae8074969
# ╠═1ad1053e-d674-4ada-8086-017cb77f6af5
# ╠═777b4e2e-a9a5-4e92-9f7e-33c0dc07d3fb
# ╠═75922729-c1b7-4aa3-91af-5d9515f92711
# ╟─714ab700-5aa8-45ff-94c6-0d762b72504a
# ╠═dab3ec1d-5c45-4081-b778-7ae6e61bc77f
# ╠═a7b4ab24-c4a2-47c1-b0a4-4799c222e626
# ╠═98bf26ce-a086-44bc-8a94-2968d0ee34aa
# ╠═788a12f3-4078-4e06-8027-7fb259d6e7fb
# ╟─4b773131-2161-4beb-bf01-fdb63229cadc
# ╠═4fd2d8ce-4727-450f-a185-994fb2b184b6
# ╠═6e055fa8-47a1-4ba3-9924-e211a13ad63d
# ╠═df5c9933-51ea-4cc8-a4fe-4a0aeffcba74
# ╠═0444d257-aa50-4115-b290-af5e0272d752
# ╠═50ed57c0-b48a-4ced-a1a4-97c316e3bcfd
# ╠═6a9a0584-a949-49b5-ad74-8417ed935eef
# ╠═e17a406f-3776-4eec-87bc-0ce275dee59f
# ╠═642b4b4e-85eb-4358-8c04-6f9a542fc5df
# ╠═1d31c563-e394-4fb9-8251-ed7a9cf8bade
# ╠═519d7e7b-377c-4b26-aa19-58088c3cd958
# ╠═3f58e0a1-8d81-44bb-be22-c6b3f325cc62
# ╠═8a05c190-8967-441a-8b20-ef74a0d3a84f
# ╠═8855cd39-fbea-4e81-b9db-ac066322d76e
# ╠═08e4a7f0-5aee-4326-83f5-082f47e9016a
# ╠═91c68e1a-73f7-45a2-b9af-eae7527f405f
# ╠═0cbb32f7-fbdb-47f7-a12d-bfdc89d8e37f
# ╠═c5d3b354-f4b6-4f59-b3c2-c2946efb8b5f
# ╠═7812c9d9-92ea-47f8-b836-59249aa6820d
# ╠═98ef0760-126a-425a-bcda-fe42079d8b76
# ╠═b484641d-5e48-453a-a140-64844296e588
# ╠═7a3ea09a-d1e0-4f32-adde-786f4144cd7a
# ╠═d2251d57-0c11-45b4-864d-b49b371ddd35
# ╠═0f442405-224d-4476-9bec-14127d22dcf5
# ╠═8b08ce23-bfd9-403b-936d-0f4d867d6731
# ╠═ccc085e6-1571-4f87-9772-791a7750455f
# ╠═7378c4c4-2f27-4e3b-9ab0-f8a2bdc492cc
# ╠═d6f7504a-e067-4772-a538-a88efa190e0f
# ╠═f4405647-c79d-4a0d-a16f-1572bafa507a
# ╠═474b504a-b4fa-4ff2-a0d2-70ee3129db05
# ╠═3939f0d2-545a-46ef-895e-00039eb092cd
# ╠═3bd61ff5-4c9b-49e9-8f55-fcf45c793164
# ╠═83ada54e-e3f2-4f2d-840e-ec7d3d3237b3
# ╟─bc479a76-7d79-4046-80ae-c553fe3af892
# ╠═a7b71a4d-0779-49e0-9781-b4c03877e2dc
# ╠═6ac60f86-633d-4429-b3bc-3030b0190447
# ╠═5fcb209f-a748-41c6-b1bd-90cd7b63a5dc
# ╠═f7144212-8f81-4396-862e-af255afb2a9b
# ╟─13a5ac47-c607-4e4a-ae53-e558c74a8944
# ╠═672c3758-5690-4688-8d2f-56621cc2d985
