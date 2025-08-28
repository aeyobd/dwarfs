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
basename = "1e5_new_v31_r3.2/L3M11_9Gyr_smallperi.3b/"

# ╔═╡ 777b4e2e-a9a5-4e92-9f7e-33c0dc07d3fb
expected_orbit_file = "1e6_new_v31_r3.2/L3M11_9Gyr_smallperi/simulation/orbit.csv"

# ╔═╡ 75922729-c1b7-4aa3-91af-5d9515f92711
time_match = -0 / T2GYR

# ╔═╡ fe5c1723-d829-4679-b4af-c7eb24c6297c
stepsize = 1

# ╔═╡ 8601fbd4-1636-449c-b625-f85b82f666a7
secondary_method = true

# ╔═╡ 714ab700-5aa8-45ff-94c6-0d762b72504a
md"""
# Data Loading

Each model only requires a `centres.hdf5` file satesfying our centres file format. The central orbit also need the `simulation/agama_potential.ini` to load in the potential and the `simulation/orbit.toml` to load in the orbit properties.

"""

# ╔═╡ dab3ec1d-5c45-4081-b778-7ae6e61bc77f
function load_orbit(extension)
	filename = joinpath(basename, "orbit" * extension, "centres.hdf5")
	return Orbit(filename)
end

# ╔═╡ a7b4ab24-c4a2-47c1-b0a4-4799c222e626
orbit_expected = Orbit(expected_orbit_file)

# ╔═╡ 98bf26ce-a086-44bc-8a94-2968d0ee34aa
pot = Agama.Potential(file = joinpath(basename, "orbit_0", "simulation/agama_potential.ini"))

# ╔═╡ 788a12f3-4078-4e06-8027-7fb259d6e7fb
orbits = [
	load_orbit("_0"),
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

# ╔═╡ 6fb2b2dd-949e-431e-b3e8-f9f573452353
@assert all(orbit.times[1] ≈ orbits[1].times[1] for orbit in orbits)

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
import LinearAlgebra: normalize, ×, ⋅,  norm

# ╔═╡ 6a9a0584-a949-49b5-ad74-8417ed935eef
function orbit_unit_vectors(pos_old, vel_old)
	z_hat = normalize(pos_old × vel_old)
	R_hat = normalize(pos_old)
	phi_hat = R_hat × z_hat
	
	return [R_hat z_hat phi_hat]'
end

# ╔═╡ e17a406f-3776-4eec-87bc-0ce275dee59f
idx_f = argmin(abs.(orbit_expected.times .- time_match))

# ╔═╡ a70ed721-3c14-40b1-98b8-3f0103b4e237
LilGuys.plot_xyz(positions(orbit_expected)[:, idx_f:idx_f], positions(orbits[1])[:, argmin(abs.(orbits[1].times .- time_match))], plot=:scatter)

# ╔═╡ 5dd21ae3-34ea-44a1-ba7a-7f24f4544970
LilGuys.plot_xyz(velocities(orbit_expected)[:, idx_f:idx_f], velocities(orbits[1])[:, argmin(abs.(orbits[1].times .- time_match))], plot=:scatter)

# ╔═╡ 642b4b4e-85eb-4358-8c04-6f9a542fc5df
coord_sys_i = orbit_unit_vectors(orbits[1].positions[:, 1], orbits[1].velocities[:, 1])

# ╔═╡ 1d31c563-e394-4fb9-8251-ed7a9cf8bade
coord_sys_f = orbit_unit_vectors(orbit_expected.positions[:, idx_f], orbit_expected.velocities[:, idx_f])

# ╔═╡ 519d7e7b-377c-4b26-aa19-58088c3cd958
pos_i = coord_sys_i * orbits[1].positions[:, 1]

# ╔═╡ 2eab8160-5ba2-45fc-930d-8b666516f5f9
[coord_sys_i * orbit.positions[:, 1] for orbit in orbits]

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
	pos_p = coord_sys_f * orbit.positions[:, idx] .- pos_f
	vel_p = v_scale * coord_sys_f * orbit.velocities[:, idx] .- v_scale * vel_f
	return [pos_p; vel_p]
end

# ╔═╡ 7812c9d9-92ea-47f8-b836-59249aa6820d
md"""
# Calculations

Given now our initial and final respective coordinate frames, the goal (following vasiliev2024) is to solve for the Jacobian so we can write

``
w \approx J u + w_0
``

where `w_0` is the current `w` vector for the central orbit and ``J`` is the Jacobian.
"""

# ╔═╡ 6d18bf9a-c46c-4bb8-99e8-53244496405a
md"""
## Forward method
"""

# ╔═╡ f67a9e33-8a74-4ee8-a396-4eedbe402ccb
Ut = hcat(initial_projected.(orbits)...)'

# ╔═╡ b484641d-5e48-453a-a140-64844296e588
Wt = hcat(final_projected.(orbits)...)'

# ╔═╡ 98ef0760-126a-425a-bcda-fe42079d8b76
O_Ut = hcat(ones(length(orbits)), Ut)

# ╔═╡ 7a3ea09a-d1e0-4f32-adde-786f4144cd7a
xi_J = O_Ut \ Wt

# ╔═╡ d2251d57-0c11-45b4-864d-b49b371ddd35
xi = xi_J[1, :]

# ╔═╡ 8b08ce23-bfd9-403b-936d-0f4d867d6731
J = xi_J[2:end, :]'

# ╔═╡ d6f7504a-e067-4772-a538-a88efa190e0f
u_next_1 = -inv(J) * xi

# ╔═╡ 9bc17881-c7af-49c2-972d-4cc01f471944
md"""
## Secondary method
In this version, the equation is rearranged to solve for the inverse-Jacobian and next initial coordinates given the final coordinates
"""

# ╔═╡ 02939c82-f800-49d7-a0b4-05b135fa24b1
O_Wt = hcat(ones(length(orbits)), Wt)

# ╔═╡ fb1df444-cf2e-4cc3-89fe-91aa00fc3151
u_Jinv = O_Wt \ Ut

# ╔═╡ 393e8a37-cf04-47f1-aeb4-3014323cba5b
u_next_2 = u_Jinv[1, :]

# ╔═╡ 4c44d3e8-f29c-45c0-ab85-23159d889d9f
J_2 = inv(u_Jinv[2:end, :]')

# ╔═╡ 506b20fe-88f5-4768-81b8-0da2731cb023
md"""
## New orbit
"""

# ╔═╡ 3216748d-d787-4470-9a6c-4220c84a9979
u_next = secondary_method ? u_next_2 : u_next_1

# ╔═╡ 3939f0d2-545a-46ef-895e-00039eb092cd
pos_next, vel_next = initial_unprojected(stepsize * u_next)

# ╔═╡ b10f8873-19d8-4e8a-8166-d0f92e8fa920
orbits[1].positions[:, 1], V2KMS*orbits[1].velocities[:, 1]

# ╔═╡ e63b7dc3-98cf-4c22-b562-c0b28a5012c3
vel_next * V2KMS .- orbits[1].velocities[:, 1] *V2KMS

# ╔═╡ 54d04008-0c98-4986-aba6-65e4869c42ae
pos_next .- orbits[1].positions[:, 1]

# ╔═╡ 83ada54e-e3f2-4f2d-840e-ec7d3d3237b3
orbit_next = LilGuys.agama_orbit(pot, Galactocentric(pos_next, vel_next * V2KMS), timerange=(orbits[1].times[1], orbit_expected.times[end]), agama_units=Agama.VASILIEV_UNITS)

# ╔═╡ dee31bd0-7fa0-4bbd-9621-39702c5a81bc
LilGuys.plot_xyz(positions(orbit_expected), positions(orbit_next), positions(orbits[1]))

# ╔═╡ 6e055fa8-47a1-4ba3-9924-e211a13ad63d
LilGuys.plot_xyz(positions(orbit_expected), positions(orbit_next), positions.(orbits)...)

# ╔═╡ df5c9933-51ea-4cc8-a4fe-4a0aeffcba74
let
	fig = Figure()

	ax = Axis(fig[1,1], xlabel="time", ylabel="radius")
	for orbit in orbits
		lines!(orbit.times, radii(orbit))
	end

	lines!(orbit_expected.times, radii(orbit_expected), color=:red)

	lines!(orbit_next.times, radii(orbit_next), color=:black)

	vlines!(time_match)
	fig

end

# ╔═╡ 929ebd0e-394a-4fdf-9788-785b9319595a
md"""
## Comparison
"""

# ╔═╡ 23a18c59-0bae-43e4-8b4f-3022b23eba14
u_next_1

# ╔═╡ a909ae85-728f-4d92-82c2-d4c4d0bb61a8
u_next_2

# ╔═╡ dc0345cc-1ffd-47bf-8f6e-2549ebb8c064
J

# ╔═╡ c00695e0-72d1-4d55-9d9f-0d89b872aa63
J_2

# ╔═╡ 09aee1e6-8a7a-4e21-90de-9ec650cf613a
md"""
# Sanity Checks
"""

# ╔═╡ 8f347586-dbc3-4eae-b6f7-ebbbb471a237
round.(O_Ut[:, 2:end], digits=3)

# ╔═╡ db741309-cdc6-4518-be61-a755e16b63d7
u_next

# ╔═╡ 6a8b06af-a460-4a48-936a-ef657b403e7e
Wt .- xi'

# ╔═╡ c1c322d9-d3c5-4ce5-8293-cbd68bad1310
@assert isapprox(J * u_next_1 + xi, zeros(6), atol=1e-10)

# ╔═╡ b222b35b-c549-49a5-9d60-d6950815c5f3
O_Wt * u_Jinv .- Ut

# ╔═╡ bb6d3ef7-3bd5-4996-a131-a1cbf962d715
Ut

# ╔═╡ b8e18450-8fd2-4e02-96c6-1666f0a4f9d1
@assert isapprox(J * u_next_1 .+ xi, zeros(6), atol=1e-10)

# ╔═╡ f4405647-c79d-4a0d-a16f-1572bafa507a
for orbit in orbits
	println(initial_unprojected(inv(J) * (final_projected(orbits[2]) .- xi)))
	println(orbits[2].positions[:, 1], orbits[2].velocities[:, 1])
end

# ╔═╡ 3bd61ff5-4c9b-49e9-8f55-fcf45c793164
@assert all(initial_unprojected(initial_projected(orbits[2])) .≈ (orbits[2].positions[:, 1], orbits[2].velocities[:, 1]))

# ╔═╡ d6820e63-441c-44d7-b9fc-a2f2a1535f05
for orbit in orbits
	u_act =final_projected(orbit)
	u_pred = J * initial_projected(orbit) + xi
	println("actual: ", u_act)
	println("predicted: ", u_pred)
	println("Δ\t", norm(u_act - u_pred))
end

# ╔═╡ bc479a76-7d79-4046-80ae-c553fe3af892
md"""
# Data Writing
"""

# ╔═╡ f7144212-8f81-4396-862e-af255afb2a9b
obs_props = let
	df = TOML.parsefile(joinpath(basename, "orbit_0", "simulation/orbit.toml"))
	dm = LilGuys.kpc2dm(df["distance"] ± df["distance_err"])
	df["distance_modulus"] = Measurements.value(dm)
	df["distance_modulus_err"] = Measurements.uncertainty(dm)

	df["ra_err"] = 0
	df["dec_err"] = 0
	df
end

# ╔═╡ 7a1c3dfc-037f-486b-ba91-2b3f5d36cf61
outfile = let
	filename = "orbit_next"
	if stepsize != 1
		filename *= "_$stepsize"
	end

	if secondary_method
		filename *= "_b"
	end

end

# ╔═╡ 5fcb209f-a748-41c6-b1bd-90cd7b63a5dc
write(joinpath(basename, "$outfile.csv"), orbit_next)

# ╔═╡ 672c3758-5690-4688-8d2f-56621cc2d985
open(joinpath(basename, "$outfile.toml"), "w") do f
	d = OrderedDict(
		"lastmodel" => basename,
		"x_i" => orbit_next.positions[1, 1],
		"y_i" => orbit_next.positions[2, 1],
		"z_i" => orbit_next.positions[3, 1],
		"v_x_i" => orbit_next.velocities[1, 1] * V2KMS,
		"v_y_i" => orbit_next.velocities[2, 1] * V2KMS,
		"v_z_i" => orbit_next.velocities[3, 1] * V2KMS,
		"t_match" => time_match,
		"method" => secondary_method ? "backwards" : "forwards"
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

# ╔═╡ d80a4e14-a901-4095-92a2-09fc478a5d34
joinpath(basename, "$outfile.toml")

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
# ╠═fe5c1723-d829-4679-b4af-c7eb24c6297c
# ╠═8601fbd4-1636-449c-b625-f85b82f666a7
# ╟─714ab700-5aa8-45ff-94c6-0d762b72504a
# ╠═dab3ec1d-5c45-4081-b778-7ae6e61bc77f
# ╠═a7b4ab24-c4a2-47c1-b0a4-4799c222e626
# ╠═98bf26ce-a086-44bc-8a94-2968d0ee34aa
# ╠═788a12f3-4078-4e06-8027-7fb259d6e7fb
# ╠═6fb2b2dd-949e-431e-b3e8-f9f573452353
# ╟─4b773131-2161-4beb-bf01-fdb63229cadc
# ╠═4fd2d8ce-4727-450f-a185-994fb2b184b6
# ╠═a70ed721-3c14-40b1-98b8-3f0103b4e237
# ╠═5dd21ae3-34ea-44a1-ba7a-7f24f4544970
# ╠═dee31bd0-7fa0-4bbd-9621-39702c5a81bc
# ╠═6e055fa8-47a1-4ba3-9924-e211a13ad63d
# ╠═df5c9933-51ea-4cc8-a4fe-4a0aeffcba74
# ╠═0444d257-aa50-4115-b290-af5e0272d752
# ╠═50ed57c0-b48a-4ced-a1a4-97c316e3bcfd
# ╠═6a9a0584-a949-49b5-ad74-8417ed935eef
# ╠═e17a406f-3776-4eec-87bc-0ce275dee59f
# ╠═642b4b4e-85eb-4358-8c04-6f9a542fc5df
# ╠═1d31c563-e394-4fb9-8251-ed7a9cf8bade
# ╠═519d7e7b-377c-4b26-aa19-58088c3cd958
# ╠═2eab8160-5ba2-45fc-930d-8b666516f5f9
# ╠═3f58e0a1-8d81-44bb-be22-c6b3f325cc62
# ╠═8a05c190-8967-441a-8b20-ef74a0d3a84f
# ╠═8855cd39-fbea-4e81-b9db-ac066322d76e
# ╠═08e4a7f0-5aee-4326-83f5-082f47e9016a
# ╠═91c68e1a-73f7-45a2-b9af-eae7527f405f
# ╠═0cbb32f7-fbdb-47f7-a12d-bfdc89d8e37f
# ╠═c5d3b354-f4b6-4f59-b3c2-c2946efb8b5f
# ╠═7812c9d9-92ea-47f8-b836-59249aa6820d
# ╟─6d18bf9a-c46c-4bb8-99e8-53244496405a
# ╠═f67a9e33-8a74-4ee8-a396-4eedbe402ccb
# ╠═b484641d-5e48-453a-a140-64844296e588
# ╠═98ef0760-126a-425a-bcda-fe42079d8b76
# ╠═7a3ea09a-d1e0-4f32-adde-786f4144cd7a
# ╠═d2251d57-0c11-45b4-864d-b49b371ddd35
# ╠═8b08ce23-bfd9-403b-936d-0f4d867d6731
# ╠═d6f7504a-e067-4772-a538-a88efa190e0f
# ╟─9bc17881-c7af-49c2-972d-4cc01f471944
# ╠═02939c82-f800-49d7-a0b4-05b135fa24b1
# ╠═fb1df444-cf2e-4cc3-89fe-91aa00fc3151
# ╠═393e8a37-cf04-47f1-aeb4-3014323cba5b
# ╠═4c44d3e8-f29c-45c0-ab85-23159d889d9f
# ╟─506b20fe-88f5-4768-81b8-0da2731cb023
# ╠═3216748d-d787-4470-9a6c-4220c84a9979
# ╠═3939f0d2-545a-46ef-895e-00039eb092cd
# ╠═b10f8873-19d8-4e8a-8166-d0f92e8fa920
# ╠═e63b7dc3-98cf-4c22-b562-c0b28a5012c3
# ╠═54d04008-0c98-4986-aba6-65e4869c42ae
# ╠═83ada54e-e3f2-4f2d-840e-ec7d3d3237b3
# ╟─929ebd0e-394a-4fdf-9788-785b9319595a
# ╠═23a18c59-0bae-43e4-8b4f-3022b23eba14
# ╠═a909ae85-728f-4d92-82c2-d4c4d0bb61a8
# ╠═dc0345cc-1ffd-47bf-8f6e-2549ebb8c064
# ╠═c00695e0-72d1-4d55-9d9f-0d89b872aa63
# ╟─09aee1e6-8a7a-4e21-90de-9ec650cf613a
# ╠═8f347586-dbc3-4eae-b6f7-ebbbb471a237
# ╠═db741309-cdc6-4518-be61-a755e16b63d7
# ╠═6a8b06af-a460-4a48-936a-ef657b403e7e
# ╠═c1c322d9-d3c5-4ce5-8293-cbd68bad1310
# ╠═b222b35b-c549-49a5-9d60-d6950815c5f3
# ╠═bb6d3ef7-3bd5-4996-a131-a1cbf962d715
# ╠═b8e18450-8fd2-4e02-96c6-1666f0a4f9d1
# ╠═f4405647-c79d-4a0d-a16f-1572bafa507a
# ╠═3bd61ff5-4c9b-49e9-8f55-fcf45c793164
# ╠═d6820e63-441c-44d7-b9fc-a2f2a1535f05
# ╟─bc479a76-7d79-4046-80ae-c553fe3af892
# ╠═a7b71a4d-0779-49e0-9781-b4c03877e2dc
# ╠═6ac60f86-633d-4429-b3bc-3030b0190447
# ╠═5fcb209f-a748-41c6-b1bd-90cd7b63a5dc
# ╠═f7144212-8f81-4396-862e-af255afb2a9b
# ╠═7a1c3dfc-037f-486b-ba91-2b3f5d36cf61
# ╠═672c3758-5690-4688-8d2f-56621cc2d985
# ╠═d80a4e14-a901-4095-92a2-09fc478a5d34
