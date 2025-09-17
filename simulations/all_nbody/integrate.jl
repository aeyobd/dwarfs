### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 020abc48-8f4a-11f0-0fc3-a5343e45722c
begin
	import Pkg; Pkg.activate()

	import Agama

	using LilGuys
	#using Arya
	using CairoMakie

	using CSV, DataFrames
end

# ╔═╡ 9d18bd8f-33fe-40af-807c-7185f97552ff
md"""
This notebook runs a convergence test on the orbits
"""

# ╔═╡ ab6e67aa-6bfc-4391-8ce1-7fb2302f43fa
module NBody 
	include("nbody_utils.jl")
end

# ╔═╡ 8a40a212-951d-4a47-806d-935e88d7b4db
CairoMakie.activate!(type=:png)

# ╔═╡ 3a168281-c517-4281-a47e-9b1612588d0c
import Arya

# ╔═╡ fd73df12-81ef-47c7-b05b-fe6a208fa8aa
pot = NBody.get_potential_c("vasiliev24/L3M11/potential")

# ╔═╡ 8c9c0c16-097b-4db4-b81d-09cff946de1d
units = Agama.VASILIEV_UNITS

# ╔═╡ 633417cd-475d-4935-a8a0-0847aef02745
md"""
# Nbody
"""

# ╔═╡ 6241fb59-d0ed-4716-a57a-5e55694a741f
md"""
## Initial Conditions
"""

# ╔═╡ 17e00154-85f5-48bc-af21-c5a65816638b
obs_props = let
	df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv"), DataFrame)

	df[!, :distance_modulus] = LilGuys.kpc2dm.(df.distance)
	df[!, :ra_err] .= 0
	df[!, :dec_err] .= 0
	df
end

# ╔═╡ ae8be670-ff7b-4b09-8e61-c596a50638f9
md"""
# Orbit integration
"""

# ╔═╡ 9294def8-7950-4104-8f6b-32d7ea27149a
ic_example = NBody.sample_initial_conditions(obs_props)

# ╔═╡ 308f0ce4-9db8-4030-850b-d181231f8fb1
masses_initial = [mass(h, 10*h.r_s) for h in ic_example[3]]

# ╔═╡ d1042dfa-b0b3-436a-ae05-df911b0c6091
maximum(masses_initial)

# ╔═╡ e7fa630b-60fb-43bd-811e-983a5a6d89b5
sort(LilGuys.M200.(ic_example[3]))

# ╔═╡ decc6840-84ab-451e-a629-4fbefe3c2f8a
traj_evolving = NBody.integrate_particles(
	ic_example...,
	NBody.force_potential_nbody(pot, units),
	timestep=0.1,
)

# ╔═╡ 0c7dbdce-a5de-41cf-affe-92f90cc33aea
traj_evolving_lr = NBody.integrate_particles(
	ic_example...,
	NBody.force_potential_nbody(pot, units),
	timestep=0.2,
)

# ╔═╡ 05cbefe2-1a19-4ed6-8a5b-b1b166c0d9e7
pot_py = NBody.get_potential("vasiliev24/L3M11/potential")

# ╔═╡ ddd37893-ae77-4b01-83f6-cc040ec60748
pot_static_py = NBody.get_potential("vasiliev24/L3M11/potential_mw_init")

# ╔═╡ 12d37752-507c-4a88-bd72-3b2f3192eab7
f_fric_nbody = 	NBody.force_potential_nbody_friction(pot_py, units, ic_example[3], pot_sigma=pot_static_py)

# ╔═╡ a7ae84c1-438d-443b-a650-46f3dd792946
traj_fric_nbody = NBody.integrate_particles(
	ic_example[1], ic_example[2], ic_example[3],
	f_fric_nbody,
	timestep=0.1,
	verbose=true
)

# ╔═╡ 02bf6dce-a580-4a2d-8acc-ae496422da53
f_fric = NBody.force_potential_friction(pot_py, units, ic_example[3], pot_sigma=pot_static_py)

# ╔═╡ 5f58bda1-2f01-484f-82d7-c9d30135b382
traj_fric = NBody.integrate_particles(
	ic_example[1], ic_example[2], ic_example[3],
	f_fric
	,
	timestep=0.1,
	verbose=true
)

# ╔═╡ fa5621bd-d2b2-4904-9e73-dbcd8ebf3547
traj_evolving_hr = NBody.integrate_particles(
	ic_example...,
	NBody.force_potential_nbody(pot, units),
	timestep=0.04,
)

# ╔═╡ 34c78afb-3719-4caa-99ea-df07229e60ef
traj_evolving_hhr = NBody.integrate_particles(
	ic_example...,
	NBody.force_potential_nbody(pot, units),
	timestep=0.02,
)

# ╔═╡ 0c9056a3-5ce5-45d1-84f6-732d9b46c7f6
traj_static = NBody.integrate_particles(
	ic_example...,
	NBody.force_potential(pot, units),
	timestep = 0.1,
)

# ╔═╡ 8218e54b-47d7-422a-a90e-c875b5927af1
md"""
# Plots
"""

# ╔═╡ 2e566943-e4af-42c3-a9a9-711d64ecd9fb
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

# ╔═╡ c51eb0a7-971b-4ec8-884a-ac8a545f5a2e
plot_traj(traj_evolving_hr, limits=(-200, 200, -200, 200))

# ╔═╡ cf5d50ac-bcde-4c29-a77c-c74203d50493
plot_traj(traj_fric, limits=(-200, 200, -200, 200))

# ╔═╡ 6305cae3-5ee4-45c8-ae8b-2f1eb699382c
plot_traj(traj_evolving .- LilGuys.resample.(traj_evolving_hr, [traj_evolving[1].times]))

# ╔═╡ d4565cef-6c88-4fec-8815-318ae1af10f2
plot_traj(traj_evolving_hr .- LilGuys.resample.(traj_evolving_hhr, [traj_evolving_hr[1].times]))

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

# ╔═╡ e023cfc7-5702-437f-8963-beb0e649b777
let
	fig = Figure()
	ax = Axis(fig[1,1])


	orbits = traj_evolving_hr .- LilGuys.resample.(traj_evolving_hhr, [traj_evolving_hr[1].times])

	for orbit in orbits
		lines!(orbit.times, log10.(LilGuys.radii(orbit)))
	end

	ylims!(-1, 3)
	fig
end

# ╔═╡ 15223439-9fd7-4600-ba50-c934b7e8ba86
md"""
# Galaxy by galaxy comparison
"""

# ╔═╡ 7055d232-3e2f-4e74-a99c-a021e8fd15bd
galaxynames = obs_props.galaxyname

# ╔═╡ ae7d8d0d-7ab4-4484-8224-06a5eadf347d
galaxynames

# ╔═╡ c059a0c9-5e92-40d5-8aaa-af7cdf958a80
galaxynames[sortperm(LilGuys.M200.(ic_example[3]))]

# ╔═╡ d080a8a4-1822-4835-ab45-377e9d77360f
function get_galaxy(orbits, galaxyname)
	idx = findall(galaxynames .== [galaxyname]) |> only

	return orbits[idx]
end

# ╔═╡ 769b22b6-6f63-4bd9-9668-1486c827bb80
function compare_orbits(orbits...; galaxyname)

	fig = Figure(size=(6*72, 4*72))
	
	ax_xy = Axis(fig[1,1], )
	ax_yz = Axis(fig[1,2], )
	ax_xz = Axis(fig[1,3], )

	for orbit in orbits
		o = get_galaxy(orbit, galaxyname)
		lines!(ax_xy, o.positions[1, :], o.positions[2, :])
		lines!(ax_yz, o.positions[2, :], o.positions[3, :])
		lines!(ax_xz, o.positions[1, :], o.positions[3, :])
	end


	ax_rt = Axis(fig[2, 1:2])

	for orbit in orbits
		o = get_galaxy(orbit, galaxyname)
		lines!(ax_rt, o.times, radii(o))
	end


	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 2, Aspect(1, 1.0))

	linkaxes!(ax_xy, ax_yz, ax_xz)

	
	fig
end

# ╔═╡ 3c372f29-afba-49b9-bf2f-8e0e7fd337e0
compare_orbits(traj_evolving_lr, traj_evolving, traj_evolving_hr, traj_evolving_hhr, galaxyname="sculptor")

# ╔═╡ d41660c5-f49c-4a2d-8bfb-bd4e55eea0f6
compare_orbits(traj_evolving, traj_static, traj_fric, traj_fric_nbody, galaxyname="sculptor")

# ╔═╡ 594b7506-4d43-41f2-967e-420af104628a
compare_orbits(traj_evolving_lr, traj_evolving, traj_evolving_hr, traj_evolving_hhr, galaxyname="ursa_minor")

# ╔═╡ 5f277b73-01ec-4f8b-9bcb-eba1e5b31c2e
compare_orbits(traj_evolving, traj_static, traj_fric, traj_fric_nbody, galaxyname="ursa_minor")

# ╔═╡ 2874bb7a-5487-4ef7-9ee4-f0b7c5367b2a
compare_orbits(traj_evolving, traj_evolving_hr, traj_evolving_hhr, galaxyname="bootes3")

# ╔═╡ Cell order:
# ╠═9d18bd8f-33fe-40af-807c-7185f97552ff
# ╠═020abc48-8f4a-11f0-0fc3-a5343e45722c
# ╠═ab6e67aa-6bfc-4391-8ce1-7fb2302f43fa
# ╠═8a40a212-951d-4a47-806d-935e88d7b4db
# ╠═3a168281-c517-4281-a47e-9b1612588d0c
# ╠═308f0ce4-9db8-4030-850b-d181231f8fb1
# ╠═d1042dfa-b0b3-436a-ae05-df911b0c6091
# ╠═ae7d8d0d-7ab4-4484-8224-06a5eadf347d
# ╠═fd73df12-81ef-47c7-b05b-fe6a208fa8aa
# ╠═8c9c0c16-097b-4db4-b81d-09cff946de1d
# ╟─633417cd-475d-4935-a8a0-0847aef02745
# ╟─6241fb59-d0ed-4716-a57a-5e55694a741f
# ╠═17e00154-85f5-48bc-af21-c5a65816638b
# ╟─ae8be670-ff7b-4b09-8e61-c596a50638f9
# ╠═9294def8-7950-4104-8f6b-32d7ea27149a
# ╠═e7fa630b-60fb-43bd-811e-983a5a6d89b5
# ╠═c059a0c9-5e92-40d5-8aaa-af7cdf958a80
# ╠═decc6840-84ab-451e-a629-4fbefe3c2f8a
# ╠═0c7dbdce-a5de-41cf-affe-92f90cc33aea
# ╠═05cbefe2-1a19-4ed6-8a5b-b1b166c0d9e7
# ╠═ddd37893-ae77-4b01-83f6-cc040ec60748
# ╠═12d37752-507c-4a88-bd72-3b2f3192eab7
# ╠═a7ae84c1-438d-443b-a650-46f3dd792946
# ╠═02bf6dce-a580-4a2d-8acc-ae496422da53
# ╠═5f58bda1-2f01-484f-82d7-c9d30135b382
# ╠═fa5621bd-d2b2-4904-9e73-dbcd8ebf3547
# ╠═34c78afb-3719-4caa-99ea-df07229e60ef
# ╠═0c9056a3-5ce5-45d1-84f6-732d9b46c7f6
# ╟─8218e54b-47d7-422a-a90e-c875b5927af1
# ╠═2e566943-e4af-42c3-a9a9-711d64ecd9fb
# ╠═ad0354b2-5218-48a1-a992-d807625ee7f6
# ╠═c51eb0a7-971b-4ec8-884a-ac8a545f5a2e
# ╠═cf5d50ac-bcde-4c29-a77c-c74203d50493
# ╠═6305cae3-5ee4-45c8-ae8b-2f1eb699382c
# ╠═d4565cef-6c88-4fec-8815-318ae1af10f2
# ╠═36055224-3e48-4236-82d9-8787fb39ddbb
# ╠═e023cfc7-5702-437f-8963-beb0e649b777
# ╠═15223439-9fd7-4600-ba50-c934b7e8ba86
# ╠═7055d232-3e2f-4e74-a99c-a021e8fd15bd
# ╠═d080a8a4-1822-4835-ab45-377e9d77360f
# ╠═769b22b6-6f63-4bd9-9668-1486c827bb80
# ╠═3c372f29-afba-49b9-bf2f-8e0e7fd337e0
# ╠═d41660c5-f49c-4a2d-8bfb-bd4e55eea0f6
# ╠═594b7506-4d43-41f2-967e-420af104628a
# ╠═5f277b73-01ec-4f8b-9bcb-eba1e5b31c2e
# ╠═2874bb7a-5487-4ef7-9ee4-f0b7c5367b2a
