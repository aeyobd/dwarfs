### A Pluto.jl notebook ###
# v0.20.15

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

# ╔═╡ d59771ec-effe-4112-a7b3-57e69ce72437
using StaticVectors

# ╔═╡ 3a168281-c517-4281-a47e-9b1612588d0c
import Arya

# ╔═╡ 97437249-b50e-4176-a738-ab38d2ba2e51
Point = Values{3, Float64}

# ╔═╡ 2d4cd4c6-81f1-43b8-b1b9-cfc200c042f1
function radii(p::Point)
	return sqrt(p[1]^2 + p[2] ^ 2 + p[3]^2)
end

# ╔═╡ c2919a3e-11de-4aa5-876a-f7c5340fde4a
struct Keplerian
	mass::Float64
end

# ╔═╡ 4351a2c0-5a6c-448a-86a2-a45598e26c65
function acceleration(h::Keplerian, pos::Point)
	r = radii(pos)
	return -pos /r * h.mass / r^2
end

# ╔═╡ c408e4eb-8301-45b2-a685-eadb6e4310c2
function acceleration(h::NFW, pos::Point)
	M = mass(h, radii(pos))
	r = radii(pos)
	return -pos / r * M / r^2
end

# ╔═╡ 819d3930-4054-4d34-a258-ac179da10587
function check_args(initial_positions, initial_velocities, halos)
	N = length(initial_positions)
	@assert N == length(halos) == length(initial_velocities)
end

# ╔═╡ f5e2bfa4-2507-430c-93af-3f55a889eb0b
[1, undef] .=== undef

# ╔═╡ 19cc98c1-23c0-49c8-b1b1-bc5cd8fae2ad
function collect_into_orbits(all_positions, all_velocities, all_times::Vector)
	Nt = length(all_times)
	Nh = length(all_positions[1])

	orbits = Vector{Orbit}(undef, Nh)

	for i in 1:Nh
		positions = hcat((all_positions[j][i] for j in 1:Nt)...)
		velocities = hcat((all_velocities[j][i] for j in 1:Nt)...)

		orbits[i] = Orbit(positions=positions, velocities=velocities, times=all_times)
	end

	return orbits
end

# ╔═╡ 183c6686-0eb7-4619-9c4e-0f58c3ac4cf8
function nbody_acceleration!(accelerations::Vector{Point}, positions::Vector{Point}, halos)
	N = length(positions)

	for i in 1:N
		for j in 1:N 
			if j != i
				accelerations[i] += acceleration(halos[j], positions[i] .- positions[j])
			end
		end
	end


	return accelerations
end

# ╔═╡ 260e2627-50e7-40d1-bb63-0fb08c25f4ef
function integrate_nbody(initial_positions::AbstractVector{<:Point}, initial_velocities::AbstractVector{<:Point}, halos; 
						 force_extra! = (a, p, v, h, t) -> 0.,
						 skip = 100, 
						 timestep = 0.1, timerange=(0, -10/T2GYR), 
						 method = LilGuys.step_kdk, reuse_acceleration = true,
	)
	
	check_args(initial_positions, initial_velocities, halos)

	N = length(halos)
	time_i, time_f = timerange

	time_total = time_f - time_i
	itermax = ceil(Int, abs(time_total) / timestep)
	num_to_save = ceil(Int, itermax / skip) + 2
	
	all_positions = [Vector{Point}(undef, N) for _ in 1:num_to_save]
	all_velocities = [Vector{Point}(undef, N) for _ in 1:num_to_save]
	all_times = Vector{Float64}(undef, num_to_save)

	dt = timestep * sign(time_total)

	positions = deepcopy(initial_positions)
	velocities = deepcopy(initial_velocities)
	time = time_i

	accelerations = fill(Point(0, 0, 0), N)
	nbody_acceleration!(accelerations, positions, halos)


	all_positions[1] = initial_positions
	all_velocities[1] = initial_velocities
	all_times[1] = time_i
	
	is_done = false

	num_saved = 1
	for i in 1:itermax
		# Kick
		velocities .+= accelerations * dt / 2
		
		# Drift
		positions .+= velocities * dt

		# update acceleration
		fill!(accelerations, Point(0, 0, 0))
		nbody_acceleration!(accelerations, positions, halos)
		force_extra!(accelerations, positions, velocities, halos, time)

		# kick
		velocities .+= accelerations * dt / 2

		# tick
		time += dt 

		# check
		if (time > time_f > time_i) || (time .< time_i < time_f)
			is_done = true
		end

		# save
		if (i % skip == 0) || is_done
			@info "completed $i / $itermax steps"

			if is_done 
				idx = num_saved + 1
			else
				idx = (i ÷ skip) + 1
			end
			
			all_positions[idx] = copy(positions)
			all_velocities[idx] = copy(velocities)
			all_times[idx] = copy(time)
			num_saved += 1
		end


		if is_done
			break
		end
	end


	return collect_into_orbits(all_positions[1:num_saved], all_velocities[1:num_saved], all_times[1:num_saved])

end

# ╔═╡ 97cdb179-79f4-4a11-9d2d-b8287bbaaef0
function integrate_particles(initial_positions::AbstractVector{<:Point}, initial_velocities::AbstractVector{<:Point}, halos; 
						 force! = (a, p, v, h, t) -> 0.,
						 skip = 100, 
						 timestep = 0.1, timerange=(0, -10/T2GYR), 
						 method = LilGuys.step_kdk, reuse_acceleration = true,
						 kwargs...
	)
	
	check_args(initial_positions, initial_velocities, halos)

	N = length(halos)
	time_i, time_f = timerange

	time_total = time_f - time_i
	itermax = ceil(Int, abs(time_total) / timestep)
	num_to_save = ceil(Int, itermax / skip) + 2
	
	all_positions = [Vector{Point}(undef, N) for _ in 1:num_to_save]
	all_velocities = [Vector{Point}(undef, N) for _ in 1:num_to_save]
	all_times = Vector{Float64}(undef, num_to_save)

	dt = timestep * sign(time_total)

	positions = deepcopy(initial_positions)
	velocities = deepcopy(initial_velocities)
	time = time_i

	accelerations = fill(Point(0, 0, 0), N)
	nbody_acceleration!(accelerations, positions, halos)


	all_positions[1] = initial_positions
	all_velocities[1] = initial_velocities
	all_times[1] = time_i
	
	is_done = false

	num_saved = 1
	for i in 1:itermax
		# Kick
		velocities .+= accelerations * dt / 2
		
		# Drift
		positions .+= velocities * dt

		# update acceleration
		fill!(accelerations, Point(0, 0, 0))
		force!(accelerations, positions, velocities, halos, time)

		# kick
		velocities .+= accelerations * dt / 2

		# tick
		time += dt 

		# check
		if (time > time_f > time_i) || (time .< time_i < time_f)
			is_done = true
		end

		# save
		if (i % skip == 0) || is_done	
			@info "completed $i / $itermax steps"
			if is_done 
				idx = num_saved + 1
			else
				idx = (i ÷ skip) + 1
			end
			
			all_positions[idx] = copy(positions)
			all_velocities[idx] = copy(velocities)
			all_times[idx] = copy(time)
			num_saved += 1
		end


		if is_done
			break
		end
	end


	return collect_into_orbits(all_positions[1:num_saved], all_velocities[1:num_saved], all_times[1:num_saved])

end

# ╔═╡ 95479b75-7efb-46b7-b180-755b14c6b10b
h1 = Keplerian(1)

# ╔═╡ 1b4ee556-9506-4b1a-9a5c-70ff850a0e32
acceleration(h1, Point(3, 4, 0))

# ╔═╡ cf8526ad-5203-4e5b-bcd0-a419d3f91d4e
h2 = Keplerian(1)

# ╔═╡ 50f8b06b-66bf-4f11-b47f-bbf306d4803b
orbits_test = integrate_nbody(
	[Point(-50, 0, 0), Point(50, 0, 0)],
	[Point(0, 1/sqrt(2 * 100), 0), Point(0, -1/sqrt(2*100), 0)],
	[h1, h2],
	timestep=1.5, 
	timerange=(0, 10_000)
)

# ╔═╡ a200c354-ce4f-4a92-85ea-95941fd89ff0
ceil(10_001 / 1.5)

# ╔═╡ e13f15b9-0552-459a-ac92-122e29fd9334
function get_potential(potname; kwargs...)
	Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/", potname * ".ini"); kwargs...)
end

# ╔═╡ fd73df12-81ef-47c7-b05b-fe6a208fa8aa
pot = get_potential("vasiliev24/L3M11/potential")

# ╔═╡ 25b3d59d-f6c4-4556-b341-1659c56760eb
pot_static = get_potential("vasiliev24/L3M11/potential_mw_init")

# ╔═╡ 8c9c0c16-097b-4db4-b81d-09cff946de1d
units = Agama.VASILIEV_UNITS

# ╔═╡ d2dede17-1014-40db-919b-14754a12a85d
function get_sigma_v(pot)
	gm = Agama.GalaxyModel(pot, Agama.DistributionFunction(pot, type="QuasiSpherical"))


	N = 100
	Rs = 10 .^ LinRange(-1, 3, N)

	sigmas = gm._py.moments(Agama.mat2py([zeros(N) Rs zeros(N)]'), dens=false) |> Agama.py2mat

	sigma3 = sigmas[1, :] .+ sigmas[2, :] + sigmas[3, :]

	sigma1 = sigma3 ./ 3
	return LilGuys.lerp(Rs, sqrt.(sigma1) .* Agama.velocity_scale(units))
end

# ╔═╡ 5b65d7e3-588f-48ef-b4a6-b32ef21c4a35
let
	pos = LilGuys.positions.(orbits_test)
	
	fig = Figure()
	ax = Axis(fig[1,1])
	lines!(pos[1][1, :], pos[1][2, :])
	lines!(pos[2][1, :], pos[2][2, :])

	fig
end

# ╔═╡ 633417cd-475d-4935-a8a0-0847aef02745
md"""
# Nbody
"""

# ╔═╡ 13dbab77-5a87-4afb-af93-c7e77e0bb338
function force_extra(accelerations, positions, velocities, halos, time)

	for i in 1:length(positions)
		dyn_fric = LilGuys.ChandrashakarDynamicalFriction(r_s=halos[i].r_s, σv=x->σv(radii(x)), M=LilGuys.M200(halos[i]), ρ = x->Agama.density(pot, x, units), Λ = exp(4))
		
		accelerations[i] += Point(LilGuys.acceleration(dyn_fric, positions[i], velocities[i]) .+
			Agama.acceleration(pot, hcat(positions[i], units, t=time))
	end

	return accelerations
end

# ╔═╡ 0d14d82d-251a-44a3-bcf4-3187de663cd2
Agama.acceleration(pot, Point(100,1,1), units)

# ╔═╡ c396d0d0-92b0-49e0-9910-b77b31289509
pot_simple = NFW(r_circ_max = 20, v_circ_max = 1.0 )

# ╔═╡ d91fae74-2665-4325-b23a-6e3e930f3e7d
function force_extra_potential(pot, units)
	function f(accelerations, positions, velocities, halos, time)
		accelerations .+= Point.(eachcol(Agama.acceleration(pot, hcat(positions...), units, t=time)))
	end
end

# ╔═╡ d82ab1c9-829b-43bb-b464-5051dccfcab1
function force_extra_potential(pot::NFW)
	function f(accelerations, positions, velocities, halos, time)
		for i in 1:length(positions)
			accelerations[i] += Point(acceleration(pot, positions[i]))
		end
	end
end

# ╔═╡ cef2eb86-d6dd-4339-9f03-78a3819a52d7
function force_potential_nbody(pot, units)
	function f(accelerations, positions, velocities, halos, time)
		nbody_acceleration!(accelerations, positions, halos)

		for i in 1:length(positions)
			accelerations[i] += Point(Agama.acceleration(pot, positions[i], units, t=time))
		end
	end

end

# ╔═╡ ff81b9f1-1d09-4fb7-bcce-dd472490842c
function force_nbody()
	function f(accelerations, positions, velocities, halos, time)
		nbody_acceleration!(accelerations, positions, halos)
	end

end

# ╔═╡ 3cf365f6-0938-4808-9d7b-fbbbac16f8b2
σv = get_sigma_v(pot_static)

# ╔═╡ 2bf58188-d575-4fd6-9391-4eebc1386288
function force_potential_nbody_friction(accelerations, positions, velocities, halos, time)

	for i in 1:length(positions)
		dyn_fric = LilGuys.ChandrashakarDynamicalFriction(r_s=halos[i].r_s, σv=x->σv(radii(x)), M=LilGuys.M200(halos[i]), ρ = x->Agama.density(pot, x, units), Λ = exp(4))
		
		accelerations[i] += Point(LilGuys.acceleration(dyn_fric, positions[i], velocities[i]) .+
			Agama.acceleration(pot, positions[i], units, t=time))
	end

	return accelerations
end

# ╔═╡ 6241fb59-d0ed-4716-a57a-5e55694a741f
md"""
## Initial Conditions
"""

# ╔═╡ 17e00154-85f5-48bc-af21-c5a65816638b
obs_props = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv"), DataFrame)

# ╔═╡ 321ed0a2-b98e-49b0-8742-3e4c40f35778
icrs_i = LilGuys.coords_from_df(obs_props)

# ╔═╡ 95a20ae7-3a1c-4ccd-ac96-06a844390524
 M_L_star = 2.0

# ╔═╡ 39cec73b-f1fe-4928-b9d6-ba90279b0dfa
L = @. LilGuys.mag_to_L(obs_props.Mv)

# ╔═╡ d59df207-5415-487a-8571-7c5eb87d329e
 Mstar = @. L / M2MSUN * M_L_star

# ╔═╡ 9c322ae0-9c97-48ad-a648-72f728fb57c3
gc_i = LilGuys.transform.(Galactocentric, icrs_i)

# ╔═╡ 11eea25e-a9c8-438e-9462-8d8f4ee487b1
pos_i = [Point(LilGuys.position(p)) for p in gc_i]

# ╔═╡ bed9d886-1b77-44eb-a765-c212a3fb8a5e
vel_i = [Point(LilGuys.velocity(p) / V2KMS) for p in gc_i]

# ╔═╡ e276efad-dbf0-451b-b6f1-4d61d13c73a9
vcircmax = LilGuys.vel_from_M_s_fattahi.(Mstar)
   

# ╔═╡ eda4039b-e07b-4cc7-8b33-6f2e06e41937
 rcircmax = LilGuys.Ludlow.solve_rmax.(vcircmax)
    

# ╔═╡ b5ccff3a-e8b3-4738-9024-fe72882b585c
halos = [LilGuys.NFW(v_circ_max=v, r_circ_max=r) for (v, r) in zip(vcircmax, rcircmax)]

# ╔═╡ ae8be670-ff7b-4b09-8e61-c596a50638f9
md"""
# Orbit integration
"""

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

# ╔═╡ 15223439-9fd7-4600-ba50-c934b7e8ba86
md"""
# Galaxy by galaxy comparison
"""

# ╔═╡ 7055d232-3e2f-4e74-a99c-a021e8fd15bd
galaxynames = obs_props.galaxyname

# ╔═╡ f3ac66c9-6eb0-4c6e-9e7b-b4675829b0de
Nbody = length(galaxynames)

# ╔═╡ bb2db152-1b1a-4e9f-aad3-ec9392f665c8
@time integrate_particles(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=1.0,
	force! = force_extra_potential(pot, units)
)

# ╔═╡ aee858c1-8234-4c9e-8964-5e3be36b4e54
@time orbits_fast = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.1,
	skip = 100, 
	force_extra! = force_extra_potential(pot_simple)
)

# ╔═╡ 9294def8-7950-4104-8f6b-32d7ea27149a
traj_no_interact = integrate_particles(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.2,
	force! = force_extra_potential(pot, units)
)

# ╔═╡ ef83b7f7-b190-4575-ad7d-87cce6070f5b
plot_traj(traj_no_interact, limits=(-200, 200, -200, 200))

# ╔═╡ decc6840-84ab-451e-a629-4fbefe3c2f8a
traj_evolving = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.2,
	force_extra! = force_extra_potential(pot, units)
)

# ╔═╡ ad0354b2-5218-48a1-a992-d807625ee7f6
plot_traj(traj_evolving, limits=(-200, 200, -200, 200))

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
traj_evolving_lr = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.5,
	force_extra! = force_extra_potential(pot, units)
)

# ╔═╡ f4d5b775-90ff-47c2-93fa-98207b550777
plot_traj(traj_evolving_lr, limits=(-200, 200, -200, 200))

# ╔═╡ 66abdeca-d843-4d3a-9ca3-eb5880451209
traj_evolving_hr = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep=0.1,
	force_extra! = force_extra_potential(pot, units)
)

# ╔═╡ 6305cae3-5ee4-45c8-ae8b-2f1eb699382c
plot_traj(traj_evolving .- LilGuys.resample.(traj_evolving_hr, [traj_evolving[1].times]))

# ╔═╡ 8d8416cb-5b76-48eb-9f17-0ebbb98cc6da
plot_traj(traj_evolving_hr, limits=(-200, 200, -200, 200))

# ╔═╡ b0b53c44-553f-423b-966a-525cc1845f4b
traj_static = integrate_nbody(
	pos_i[1:Nbody], vel_i[1:Nbody], halos[1:Nbody],
	timestep = 0.2,
	force_extra! = force_extra_potential(pot_static, units)
)

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

# ╔═╡ Cell order:
# ╠═020abc48-8f4a-11f0-0fc3-a5343e45722c
# ╠═3a168281-c517-4281-a47e-9b1612588d0c
# ╠═d59771ec-effe-4112-a7b3-57e69ce72437
# ╠═97437249-b50e-4176-a738-ab38d2ba2e51
# ╠═2d4cd4c6-81f1-43b8-b1b9-cfc200c042f1
# ╠═c2919a3e-11de-4aa5-876a-f7c5340fde4a
# ╠═4351a2c0-5a6c-448a-86a2-a45598e26c65
# ╠═c408e4eb-8301-45b2-a685-eadb6e4310c2
# ╠═819d3930-4054-4d34-a258-ac179da10587
# ╠═260e2627-50e7-40d1-bb63-0fb08c25f4ef
# ╠═97cdb179-79f4-4a11-9d2d-b8287bbaaef0
# ╠═f5e2bfa4-2507-430c-93af-3f55a889eb0b
# ╠═19cc98c1-23c0-49c8-b1b1-bc5cd8fae2ad
# ╠═183c6686-0eb7-4619-9c4e-0f58c3ac4cf8
# ╠═95479b75-7efb-46b7-b180-755b14c6b10b
# ╠═1b4ee556-9506-4b1a-9a5c-70ff850a0e32
# ╠═cf8526ad-5203-4e5b-bcd0-a419d3f91d4e
# ╠═50f8b06b-66bf-4f11-b47f-bbf306d4803b
# ╠═a200c354-ce4f-4a92-85ea-95941fd89ff0
# ╠═e13f15b9-0552-459a-ac92-122e29fd9334
# ╠═fd73df12-81ef-47c7-b05b-fe6a208fa8aa
# ╠═25b3d59d-f6c4-4556-b341-1659c56760eb
# ╠═8c9c0c16-097b-4db4-b81d-09cff946de1d
# ╠═d2dede17-1014-40db-919b-14754a12a85d
# ╠═5b65d7e3-588f-48ef-b4a6-b32ef21c4a35
# ╟─633417cd-475d-4935-a8a0-0847aef02745
# ╠═13dbab77-5a87-4afb-af93-c7e77e0bb338
# ╠═0d14d82d-251a-44a3-bcf4-3187de663cd2
# ╠═c396d0d0-92b0-49e0-9910-b77b31289509
# ╠═d91fae74-2665-4325-b23a-6e3e930f3e7d
# ╠═d82ab1c9-829b-43bb-b464-5051dccfcab1
# ╠═cef2eb86-d6dd-4339-9f03-78a3819a52d7
# ╠═ff81b9f1-1d09-4fb7-bcce-dd472490842c
# ╠═2bf58188-d575-4fd6-9391-4eebc1386288
# ╠═3cf365f6-0938-4808-9d7b-fbbbac16f8b2
# ╟─6241fb59-d0ed-4716-a57a-5e55694a741f
# ╠═17e00154-85f5-48bc-af21-c5a65816638b
# ╠═321ed0a2-b98e-49b0-8742-3e4c40f35778
# ╠═95a20ae7-3a1c-4ccd-ac96-06a844390524
# ╠═39cec73b-f1fe-4928-b9d6-ba90279b0dfa
# ╠═d59df207-5415-487a-8571-7c5eb87d329e
# ╠═9c322ae0-9c97-48ad-a648-72f728fb57c3
# ╠═11eea25e-a9c8-438e-9462-8d8f4ee487b1
# ╠═bed9d886-1b77-44eb-a765-c212a3fb8a5e
# ╠═e276efad-dbf0-451b-b6f1-4d61d13c73a9
# ╠═eda4039b-e07b-4cc7-8b33-6f2e06e41937
# ╠═b5ccff3a-e8b3-4738-9024-fe72882b585c
# ╠═f3ac66c9-6eb0-4c6e-9e7b-b4675829b0de
# ╠═ae8be670-ff7b-4b09-8e61-c596a50638f9
# ╠═bb2db152-1b1a-4e9f-aad3-ec9392f665c8
# ╠═aee858c1-8234-4c9e-8964-5e3be36b4e54
# ╠═9294def8-7950-4104-8f6b-32d7ea27149a
# ╠═decc6840-84ab-451e-a629-4fbefe3c2f8a
# ╠═2375b143-4886-467f-ac89-bc681ded65a7
# ╠═b5230a6e-bf55-4485-9a12-b66b3d75f511
# ╠═66abdeca-d843-4d3a-9ca3-eb5880451209
# ╠═b0b53c44-553f-423b-966a-525cc1845f4b
# ╠═ad1f032e-23f1-4f60-96e3-e1ea85087489
# ╠═ad0354b2-5218-48a1-a992-d807625ee7f6
# ╠═6305cae3-5ee4-45c8-ae8b-2f1eb699382c
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
