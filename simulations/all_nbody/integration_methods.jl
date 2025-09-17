### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 020abc48-8f4a-11f0-0fc3-a5343e45722c
begin
	import Pkg; Pkg.activate()

	import Agama

	using LilGuys
	using Arya
	using CairoMakie

	using CSV, DataFrames
end

# ╔═╡ af408ef4-04bc-40ba-80c8-3ba98faecac5
module NBody 
	include("nbody_utils.jl")
end

# ╔═╡ 256d6c67-6115-4b91-87e3-57a8fd9d1122
md"""
## Simple tests
"""

# ╔═╡ 95479b75-7efb-46b7-b180-755b14c6b10b
h1 = NBody.Keplerian(1)

# ╔═╡ cf8526ad-5203-4e5b-bcd0-a419d3f91d4e
h2 = NBody.Keplerian(1)

# ╔═╡ 50f8b06b-66bf-4f11-b47f-bbf306d4803b
orbits_test = NBody.integrate_particles(
	[NBody.Point(-50, 0, 0), NBody.Point(50, 0, 0)],
	[NBody.Point(0, 1/sqrt(2 * 100), 0), NBody.Point(0, -1/sqrt(2*100), 0)],
	[h1, h2],
	NBody.force_nbody(),
	timestep=1.5, 
	timerange=(0, 10_000)
)

# ╔═╡ 5b65d7e3-588f-48ef-b4a6-b32ef21c4a35
let
	pos = LilGuys.positions.(orbits_test)
	
	fig = Figure()
	ax = Axis(fig[1,1])
	lines!(pos[1][1, :], pos[1][2, :])
	lines!(pos[2][1, :], pos[2][2, :])

	fig
end

# ╔═╡ e13f15b9-0552-459a-ac92-122e29fd9334
md"""
# Setup
"""

# ╔═╡ fd73df12-81ef-47c7-b05b-fe6a208fa8aa
pot = NBody.get_potential("vasiliev24/L3M11/potential")

# ╔═╡ cd6b0591-f1b9-4bd4-9def-18b12ca99e04
pot_c = NBody.get_potential_c("vasiliev24/L3M11/potential")

# ╔═╡ 25b3d59d-f6c4-4556-b341-1659c56760eb
pot_static = NBody.get_potential("vasiliev24/L3M11/potential_mw_init")

# ╔═╡ 8c9c0c16-097b-4db4-b81d-09cff946de1d
units = Agama.VASILIEV_UNITS

# ╔═╡ d2dede17-1014-40db-919b-14754a12a85d


# ╔═╡ 633417cd-475d-4935-a8a0-0847aef02745
md"""
# Nbody
"""

# ╔═╡ 13dbab77-5a87-4afb-af93-c7e77e0bb338


# ╔═╡ 0d14d82d-251a-44a3-bcf4-3187de663cd2
Agama.acceleration(pot, Point(100,1,1), units)

# ╔═╡ c396d0d0-92b0-49e0-9910-b77b31289509
pot_simple = NFW(r_circ_max = 20, v_circ_max = 1.0 )

# ╔═╡ 2bf58188-d575-4fd6-9391-4eebc1386288
function force_potential_nbody_friction(accelerations, positions, velocities, halos, time)

	for i in 1:length(positions)
		dyn_fric = LilGuys.ChandrashakarDynamicalFriction(r_s=halos[i].r_s, σv=x->σv(radii(x)), M=LilGuys.M200(halos[i]), ρ = x->Agama.density(pot, x, units), Λ = exp(4))
		
		accelerations[i] += (LilGuys.acceleration(dyn_fric, positions[i], velocities[i]) .+
			Agama.acceleration(pot, positions[i], units, t=time))
	end

	return accelerations
end

# ╔═╡ 3cf365f6-0938-4808-9d7b-fbbbac16f8b2


# ╔═╡ 6241fb59-d0ed-4716-a57a-5e55694a741f
md"""
## Initial Conditions
"""

# ╔═╡ f3ac66c9-6eb0-4c6e-9e7b-b4675829b0de
obs_props = let
	df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv"), DataFrame)

	df[!, :distance_modulus] = LilGuys.kpc2dm.(df.distance)
	df[!, :ra_err] .= 0
	df[!, :dec_err] .= 0
	df
end

# ╔═╡ 3b14d264-7939-4228-9753-9be27b33fccd
ic_example = NBody.sample_initial_conditions(obs_props)

# ╔═╡ ae8be670-ff7b-4b09-8e61-c596a50638f9
md"""
# Orbit integration
"""

# ╔═╡ bb2db152-1b1a-4e9f-aad3-ec9392f665c8
@time NBody.integrate_particles(
	ic_example..., NBody.force_potential(pot, units),
	timestep=0.1, 
)

# ╔═╡ 14d684ad-fd24-46b6-9f7a-5c647675a356
@time NBody.integrate_particles(
	ic_example..., NBody.force_potential(pot_c, units),
	timestep=0.1,
)

# ╔═╡ 2b389ab9-d83f-4b24-9a5d-6c8980bafe8a
@time NBody.integrate_particles(
	ic_example..., NBody.force_potential_nbody(pot_c, units),
	timestep=0.1,
)

# ╔═╡ a82d85de-4e52-4404-a63a-73c4ef6d9651
pot_simple

# ╔═╡ 1b3f1906-caf1-4d60-9dff-1403c754ce1b
pot_simple_agama = Agama.Potential(type="NFW", mass=pot_simple.M_s, scaleRadius=pot_simple.r_s)

# ╔═╡ ffdcc788-a53f-461c-98e4-70cf75b65c5f
pot_simple_agama_c = NBody.CAgama.AgamaPotential(type="NFW", mass=pot_simple.M_s, scaleRadius=pot_simple.r_s)

# ╔═╡ 45e68e58-f1a4-46ab-a677-fec704f91298
@time NBody.integrate_particles(
	ic_example..., 
	NBody.force_potential(pot_simple_agama_c, units),
	timestep=0.1,
)

# ╔═╡ 3330842b-3513-4901-8f3f-5f46724df1cb
@time NBody.integrate_particles(
	ic_example..., 
	NBody.force_potential(pot_simple, units),
	timestep=1.0,
)

# ╔═╡ fa6bf729-3a2a-4e47-8550-09bb71efa368
@time  integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.1,
	skip = 100, 
	force_extra! = force_extra_potential_c(pot_simple_agama_c, Agama.AgamaUnits())
)

# ╔═╡ d57ae910-4039-4ead-b134-cee1b3499bc4
@time  integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.1,
	skip = 100, 
	force_extra! = force_extra_potential(pot_simple)
)

# ╔═╡ aee858c1-8234-4c9e-8964-5e3be36b4e54
@time orbits_fast = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.1,
	skip = 100, 
	force_extra! = force_extra_potential(pot_simple)
)

# ╔═╡ 9294def8-7950-4104-8f6b-32d7ea27149a
@time traj_no_interact = integrate_particles(
	ic_example..., 
	NBody.force_potential(pot, units),
	timestep=0.2,

)

# ╔═╡ decc6840-84ab-451e-a629-4fbefe3c2f8a
traj_evolving = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.2,
	force_extra! = force_extra_potential(pot, units)
)

# ╔═╡ 2375b143-4886-467f-ac89-bc681ded65a7
# ╠═╡ disabled = true
#=╠═╡
traj_both = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.2,
	force_extra! =force_extra
)
  ╠═╡ =#

# ╔═╡ b5230a6e-bf55-4485-9a12-b66b3d75f511


# ╔═╡ 66abdeca-d843-4d3a-9ca3-eb5880451209


# ╔═╡ ad1f032e-23f1-4f60-96e3-e1ea85087489
function plot_traj(traj_both; kwargs...)
	fig = Figure()
	
	Axis(fig[1,1], aspect=DataAspect(); kwargs...)

	for i in 1:length(traj_both)
		pos = LilGuys.positions(traj_both[i])
		lines!(pos[2, :], pos[3, :])
	end
	
	fig
end

# ╔═╡ ad0354b2-5218-48a1-a992-d807625ee7f6
plot_traj(traj_evolving, limits=(-200, 200, -200, 200))

# ╔═╡ 6305cae3-5ee4-45c8-ae8b-2f1eb699382c
plot_traj(traj_evolving .- LilGuys.resample.(traj_evolving_hr, [traj_evolving[1].times]))

# ╔═╡ 36055224-3e48-4236-82d9-8787fb39ddbb
let
	fig = Figure()
	ax = Axis(fig[1,1])


	orbits = traj_evolving .- LilGuys.resample.(traj_evolving_hr, [traj_evolving[1].times])

	for orbit in orbits
		lines!(orbit.times, log10.(LilGuys.radii(orbit)))
	end

	ylims!(-1, 3)
	fig
end

# ╔═╡ a71b2c02-53f2-4e22-8a5e-bc4605c0d600
let
	fig = Figure()
	ax = Axis(fig[1,1])


	orbits = traj_evolving .- LilGuys.resample.(traj_no_interact, [traj_evolving[1].times])

	for orbit in orbits
		lines!(orbit.times, log10.(LilGuys.radii(orbit)))
	end

	ylims!(-1, 3)
	fig
end

# ╔═╡ 8d8416cb-5b76-48eb-9f17-0ebbb98cc6da
plot_traj(traj_evolving_hr, limits=(-200, 200, -200, 200))

# ╔═╡ f4d5b775-90ff-47c2-93fa-98207b550777
plot_traj(traj_evolving_lr, limits=(-200, 200, -200, 200))

# ╔═╡ ef83b7f7-b190-4575-ad7d-87cce6070f5b
plot_traj(traj_no_interact, limits=(-200, 200, -200, 200))

# ╔═╡ 15223439-9fd7-4600-ba50-c934b7e8ba86
md"""
# Galaxy by galaxy comparison
"""

# ╔═╡ 7055d232-3e2f-4e74-a99c-a021e8fd15bd
galaxynames = obs_props.galaxyname

# ╔═╡ d080a8a4-1822-4835-ab45-377e9d77360f
function get_galaxy(orbits, galaxyname)
	idx = findall(galaxynames .== [galaxyname]) |> only

	return orbits[idx]
end

# ╔═╡ 769b22b6-6f63-4bd9-9668-1486c827bb80
function compare_orbits(orbits...; galaxyname)

	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())

	for orbit in orbits
		o = get_galaxy(orbit, galaxyname)

		lines!(o.positions[2, :], o.positions[3, :])

	end


	fig
end

# ╔═╡ de0bf40e-8064-4c64-b7ff-e9ea97411c84
get_galaxy(traj_evolving, "sculptor")

# ╔═╡ 4be6dfb7-dc52-4205-a46e-7289a77bfd9a
compare_orbits(traj_static, orbits_fast, galaxyname="bootes3")

# ╔═╡ fa3295dc-2bd1-415b-b406-8690048207a3
compare_orbits(traj_no_interact, traj_evolving, traj_evolving_hr, galaxyname="sculptor")

# ╔═╡ 53ba5094-2999-46cc-be41-2243a047246f
compare_orbits(traj_no_interact, traj_evolving, traj_evolving_hr, galaxyname="ursa_minor")

# ╔═╡ 3c372f29-afba-49b9-bf2f-8e0e7fd337e0
compare_orbits(traj_simple, traj_evolving, traj_evolving_hr, galaxyname="sculptor")

# ╔═╡ 594b7506-4d43-41f2-967e-420af104628a


# ╔═╡ Cell order:
# ╠═020abc48-8f4a-11f0-0fc3-a5343e45722c
# ╠═af408ef4-04bc-40ba-80c8-3ba98faecac5
# ╠═256d6c67-6115-4b91-87e3-57a8fd9d1122
# ╠═95479b75-7efb-46b7-b180-755b14c6b10b
# ╠═cf8526ad-5203-4e5b-bcd0-a419d3f91d4e
# ╠═50f8b06b-66bf-4f11-b47f-bbf306d4803b
# ╠═5b65d7e3-588f-48ef-b4a6-b32ef21c4a35
# ╠═e13f15b9-0552-459a-ac92-122e29fd9334
# ╠═fd73df12-81ef-47c7-b05b-fe6a208fa8aa
# ╠═cd6b0591-f1b9-4bd4-9def-18b12ca99e04
# ╠═25b3d59d-f6c4-4556-b341-1659c56760eb
# ╠═8c9c0c16-097b-4db4-b81d-09cff946de1d
# ╠═d2dede17-1014-40db-919b-14754a12a85d
# ╟─633417cd-475d-4935-a8a0-0847aef02745
# ╠═13dbab77-5a87-4afb-af93-c7e77e0bb338
# ╠═0d14d82d-251a-44a3-bcf4-3187de663cd2
# ╠═c396d0d0-92b0-49e0-9910-b77b31289509
# ╠═2bf58188-d575-4fd6-9391-4eebc1386288
# ╠═3cf365f6-0938-4808-9d7b-fbbbac16f8b2
# ╟─6241fb59-d0ed-4716-a57a-5e55694a741f
# ╠═f3ac66c9-6eb0-4c6e-9e7b-b4675829b0de
# ╠═3b14d264-7939-4228-9753-9be27b33fccd
# ╟─ae8be670-ff7b-4b09-8e61-c596a50638f9
# ╠═bb2db152-1b1a-4e9f-aad3-ec9392f665c8
# ╠═14d684ad-fd24-46b6-9f7a-5c647675a356
# ╠═2b389ab9-d83f-4b24-9a5d-6c8980bafe8a
# ╠═45e68e58-f1a4-46ab-a677-fec704f91298
# ╠═a82d85de-4e52-4404-a63a-73c4ef6d9651
# ╠═1b3f1906-caf1-4d60-9dff-1403c754ce1b
# ╠═ffdcc788-a53f-461c-98e4-70cf75b65c5f
# ╠═3330842b-3513-4901-8f3f-5f46724df1cb
# ╠═fa6bf729-3a2a-4e47-8550-09bb71efa368
# ╠═d57ae910-4039-4ead-b134-cee1b3499bc4
# ╠═aee858c1-8234-4c9e-8964-5e3be36b4e54
# ╠═9294def8-7950-4104-8f6b-32d7ea27149a
# ╠═decc6840-84ab-451e-a629-4fbefe3c2f8a
# ╠═2375b143-4886-467f-ac89-bc681ded65a7
# ╠═b5230a6e-bf55-4485-9a12-b66b3d75f511
# ╠═66abdeca-d843-4d3a-9ca3-eb5880451209
# ╠═ad1f032e-23f1-4f60-96e3-e1ea85087489
# ╠═ad0354b2-5218-48a1-a992-d807625ee7f6
# ╠═6305cae3-5ee4-45c8-ae8b-2f1eb699382c
# ╠═36055224-3e48-4236-82d9-8787fb39ddbb
# ╠═a71b2c02-53f2-4e22-8a5e-bc4605c0d600
# ╠═8d8416cb-5b76-48eb-9f17-0ebbb98cc6da
# ╠═f4d5b775-90ff-47c2-93fa-98207b550777
# ╠═ef83b7f7-b190-4575-ad7d-87cce6070f5b
# ╠═15223439-9fd7-4600-ba50-c934b7e8ba86
# ╠═7055d232-3e2f-4e74-a99c-a021e8fd15bd
# ╠═d080a8a4-1822-4835-ab45-377e9d77360f
# ╠═769b22b6-6f63-4bd9-9668-1486c827bb80
# ╠═de0bf40e-8064-4c64-b7ff-e9ea97411c84
# ╠═4be6dfb7-dc52-4205-a46e-7289a77bfd9a
# ╠═fa3295dc-2bd1-415b-b406-8690048207a3
# ╠═53ba5094-2999-46cc-be41-2243a047246f
# ╠═3c372f29-afba-49b9-bf2f-8e0e7fd337e0
# ╠═594b7506-4d43-41f2-967e-420af104628a
