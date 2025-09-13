### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 020abc48-8f4a-11f0-0fc3-a5343e45722c
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie

	import Agama
	using CSV, DataFrames
end

# ╔═╡ c2919a3e-11de-4aa5-876a-f7c5340fde4a
struct Keplerian
	mass::Float64
end

# ╔═╡ 4351a2c0-5a6c-448a-86a2-a45598e26c65
function acceleration(h::Keplerian, pos::AbstractVector)
	r = radii(pos)
	return -pos /r * h.mass / r^3
end

# ╔═╡ c408e4eb-8301-45b2-a685-eadb6e4310c2
function acceleration(h::NFW, pos::AbstractVector)
	M = mass(h, radii(pos))
	r = radii(pos)
	return -pos /r * M / r^3
end

# ╔═╡ 819d3930-4054-4d34-a258-ac179da10587
function check_args(initial_positions, initial_velocities, halos)
	N = size(initial_positions, 2)
	@assert N == length(halos) == size(initial_velocities, 2)
	LilGuys.@assert_3vector initial_positions
	LilGuys.@assert_3vector initial_velocities
end

# ╔═╡ 183c6686-0eb7-4619-9c4e-0f58c3ac4cf8
function nbody_acceleration(positions, velocities, halos)
	accelerations = zeros(size(positions)...)
	N = size(positions, 2)

	for i in 1:N
		for j in 1:N 
			if j != i
				accelerations[:, i] .+= acceleration(halos[j], positions[:, i] .- positions[:, j])
			end
		end
	end


	return accelerations
end

# ╔═╡ 260e2627-50e7-40d1-bb63-0fb08c25f4ef
function integrate_nbody(initial_positions::AbstractMatrix, initial_velocities::AbstractMatrix, halos; 
						 force_extra = (p, v, h, t) -> 0.,
						 skip = 100, 
						 dt_max = 0.1, dt_min=0.1, timerange=(0, -10/T2GYR), η = 0.003,
						 method = LilGuys.step_kdk, reuse_acceleration = true,
						 kwargs...
	)
	
	check_args(initial_positions, initial_velocities, halos)

	N = length(halos)
	time_i, time_f = timerange

	time_total = time_f - time_i
	itermax = ceil(Int, abs(time_total) / dt_min)

	all_positions = [initial_positions]
	all_velocities = [initial_velocities]
	all_times = Float64[time_i]

	dt = dt_min * sign(time_total)
	is_done = false

	accelerations = nbody_acceleration(initial_positions, initial_velocities, halos)
	accelerations .+= force_extra(initial_positions, initial_velocities, halos, time_i)
	positions = copy(all_positions[1])
	velocities = copy(all_velocities[1])
	time = copy(time_i)
	
	for i in 1:itermax
		# leapfrog
		velocities .+= accelerations * dt / 2
		positions .+= velocities * dt
		accelerations = nbody_acceleration(positions, velocities, halos)
		accelerations .+= force_extra(positions, velocities, halos, time)

		velocities .+= accelerations * dt / 2

		time += dt 

		if i % skip == 0
			push!(all_positions, copy(positions))
			push!(all_velocities, copy(velocities))
			push!(all_times, copy(time))
		end

		if (time > time_f > time_i) || (time .< time_i < time_f)
			is_done = true
		end

		if is_done
			break
		end
	end

	return all_positions, all_velocities, all_times
end

# ╔═╡ 95479b75-7efb-46b7-b180-755b14c6b10b
h1 = Keplerian(1)

# ╔═╡ 1b4ee556-9506-4b1a-9a5c-70ff850a0e32
acceleration(h1, [1, 0, 1])

# ╔═╡ cf8526ad-5203-4e5b-bcd0-a419d3f91d4e
h2 = Keplerian(1)

# ╔═╡ 50f8b06b-66bf-4f11-b47f-bbf306d4803b
pos, vel, t = integrate_nbody(
	[-15. 15
	0 0
	0 0],
	[0. 0
	0.4/1.4/8 -0.4/1.4/8
	0 0],
	[h1, h2],
	dt_min=0.05,
	dt_max = 5
)

# ╔═╡ 661e7da3-5477-49eb-a3a7-25ff05c82889
LilGuys.v_circ(h1, 30)

# ╔═╡ 224f3f2c-d134-4299-b64c-7d97e8dcdc07
nbody_acceleration([0. 100
	0 0
	0 0],
	[0. 0
	0.05 -0.05
	0 0],
	[h1, h2],)

# ╔═╡ 3b5ff468-b3c5-4ea2-b6c6-8e0cebf744a1
acceleration(h2, [100, 0, 0])

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
	ps = cat(pos..., dims=3)
	
	f = lines(ps[1, 1,:], ps[2, 1, :])
	lines!(ps[1, 2, :], ps[2, 2, :])

	f
end

# ╔═╡ d51a04fe-e845-4c7d-929d-4a1e59726018
cat(pos..., dims=3)

# ╔═╡ 633417cd-475d-4935-a8a0-0847aef02745
md"""
# Nbody
"""

# ╔═╡ d91fae74-2665-4325-b23a-6e3e930f3e7d
function force_extra_potential(pot)
	function f(positions, velocities, halos, time)
		accelerations = zeros(size(positions)...)
		for i in 1:size(positions, 2)
			accelerations[:, i] = Agama.acceleration(pot, positions[:, i], units, t=time)
		end
		return accelerations
	end
end

# ╔═╡ 3cf365f6-0938-4808-9d7b-fbbbac16f8b2
σv = get_sigma_v(pot_static)

# ╔═╡ 13dbab77-5a87-4afb-af93-c7e77e0bb338
function force_extra(positions, velocities, halos, time)

	accelerations = zeros(size(positions)...)
	for i in 1:size(positions, 2)
		dyn_fric = LilGuys.ChandrashakarDynamicalFriction(r_s=halos[i].r_s, σv=x->σv(radii(x)), M=LilGuys.M200(halos[i]), ρ = x->Agama.density(pot, x, units), Λ = exp(4))
		
		accelerations[:, i] = LilGuys.acceleration(dyn_fric, positions[:, i], velocities[:, i]) .+
			Agama.acceleration(pot, positions[:, i], units, t=time)
	end

	return accelerations
end

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
pos_i = hcat([LilGuys.position(p) for p in gc_i]...)

# ╔═╡ bed9d886-1b77-44eb-a765-c212a3fb8a5e
vel_i = hcat([LilGuys.velocity(p) for p in gc_i]...) ./ V2KMS

# ╔═╡ e276efad-dbf0-451b-b6f1-4d61d13c73a9
vcircmax = LilGuys.vel_from_M_s_fattahi.(Mstar)
   

# ╔═╡ eda4039b-e07b-4cc7-8b33-6f2e06e41937
 rcircmax = LilGuys.Ludlow.solve_rmax.(vcircmax)
    

# ╔═╡ b5ccff3a-e8b3-4738-9024-fe72882b585c
halos = [LilGuys.NFW(v_circ_max=v, r_circ_max=r) for (v, r) in zip(vcircmax, rcircmax)]

# ╔═╡ f3ac66c9-6eb0-4c6e-9e7b-b4675829b0de
Nbody = 57

# ╔═╡ 2375b143-4886-467f-ac89-bc681ded65a7
traj_both = integrate_nbody(
	pos_i[:, 1:Nbody], vel_i[:, 1:Nbody], halos[1:Nbody],
	dt_min=0.2,
	dt_max = 5,
	force_extra=force_extra
)

# ╔═╡ decc6840-84ab-451e-a629-4fbefe3c2f8a
traj_evolving = integrate_nbody(
	pos_i[:, 1:Nbody], vel_i[:, 1:Nbody], halos[1:Nbody],
	dt_min=0.2,
	dt_max = 5,
	force_extra=force_extra_potential(pot)
)

# ╔═╡ b0b53c44-553f-423b-966a-525cc1845f4b
traj_static = integrate_nbody(
	pos_i[:, 1:Nbody], vel_i[:, 1:Nbody], halos[1:Nbody],
	dt_min=0.2,
	dt_max = 5,
	force_extra=force_extra_potential(pot_static)
)

# ╔═╡ ad1f032e-23f1-4f60-96e3-e1ea85087489
function plot_traj(traj_both)
	pos_all = traj_both[1]

	ps = cat(pos_all..., dims=3)
	fig = Figure()
	Axis(fig[1,1], aspect=DataAspect())

	for i in 1:size(pos_all[1], 2)
		lines!(ps[1, i,:], ps[2, i, :])
	end
	
	fig
end

# ╔═╡ 19b61f50-c343-4890-b46f-efed03f935b4
plot_traj(traj_both)

# ╔═╡ 5102235b-6129-472a-be05-a9a5042a4546
plot_traj(traj_static)

# ╔═╡ ad0354b2-5218-48a1-a992-d807625ee7f6
plot_traj(traj_evolving)

# ╔═╡ 8d8416cb-5b76-48eb-9f17-0ebbb98cc6da


# ╔═╡ Cell order:
# ╠═020abc48-8f4a-11f0-0fc3-a5343e45722c
# ╠═c2919a3e-11de-4aa5-876a-f7c5340fde4a
# ╠═4351a2c0-5a6c-448a-86a2-a45598e26c65
# ╠═c408e4eb-8301-45b2-a685-eadb6e4310c2
# ╠═819d3930-4054-4d34-a258-ac179da10587
# ╠═260e2627-50e7-40d1-bb63-0fb08c25f4ef
# ╠═183c6686-0eb7-4619-9c4e-0f58c3ac4cf8
# ╠═95479b75-7efb-46b7-b180-755b14c6b10b
# ╠═1b4ee556-9506-4b1a-9a5c-70ff850a0e32
# ╠═cf8526ad-5203-4e5b-bcd0-a419d3f91d4e
# ╠═50f8b06b-66bf-4f11-b47f-bbf306d4803b
# ╠═661e7da3-5477-49eb-a3a7-25ff05c82889
# ╠═224f3f2c-d134-4299-b64c-7d97e8dcdc07
# ╠═3b5ff468-b3c5-4ea2-b6c6-8e0cebf744a1
# ╠═e13f15b9-0552-459a-ac92-122e29fd9334
# ╠═fd73df12-81ef-47c7-b05b-fe6a208fa8aa
# ╠═25b3d59d-f6c4-4556-b341-1659c56760eb
# ╠═8c9c0c16-097b-4db4-b81d-09cff946de1d
# ╠═d2dede17-1014-40db-919b-14754a12a85d
# ╠═5b65d7e3-588f-48ef-b4a6-b32ef21c4a35
# ╠═d51a04fe-e845-4c7d-929d-4a1e59726018
# ╟─633417cd-475d-4935-a8a0-0847aef02745
# ╠═13dbab77-5a87-4afb-af93-c7e77e0bb338
# ╠═d91fae74-2665-4325-b23a-6e3e930f3e7d
# ╠═3cf365f6-0938-4808-9d7b-fbbbac16f8b2
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
# ╠═2375b143-4886-467f-ac89-bc681ded65a7
# ╠═decc6840-84ab-451e-a629-4fbefe3c2f8a
# ╠═b0b53c44-553f-423b-966a-525cc1845f4b
# ╠═ad1f032e-23f1-4f60-96e3-e1ea85087489
# ╠═19b61f50-c343-4890-b46f-efed03f935b4
# ╠═5102235b-6129-472a-be05-a9a5042a4546
# ╠═ad0354b2-5218-48a1-a992-d807625ee7f6
# ╠═8d8416cb-5b76-48eb-9f17-0ebbb98cc6da
