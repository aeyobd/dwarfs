import Agama
import CAgama

using LilGuys
using StaticArrays
using CSV, DataFrames

import LilGuys: acceleration

Point = SVector{3, Float64}

function LilGuys.radii(p::Point)
	return sqrt(p[1]^2 + p[2]^2 + p[3]^2)
end


struct Keplerian
	mass::Float64
end

struct Massless
end

function acceleration(h::Keplerian, pos::Point)
	r = radii(pos)
	return -pos /r * h.mass / r^2
end

function acceleration(h::Massless, pos::Point, vel=nothing)
    return Point(0, 0, 0)
end

function acceleration(h::LilGuys.NFW, pos::Point; trunc=10)
	r = radii(pos)
	M = mass(h, min(r, trunc * h.r_s))
	return -pos / r * M / r^2
end


function acceleration(pot::Agama.Potential, x, units=Agama.AgamaUnits(); time=0.)
    return Agama.acceleration(pot, x, units; t=time)
end


function acceleration(pot::CAgama.AgamaPotential, x::Point, units; time=0.)
    pot, deriv, deriv2 = CAgama.eval_potential(pot, x / Agama.length_scale(units), time=time / Agama.time_scale(units), calc_deriv=true, calc_deriv2=false)

	return -deriv * Agama.acceleration_scale(units)
end



function check_args(initial_positions, initial_velocities, halos)
	N = length(initial_positions)
	@assert N == length(halos) == length(initial_velocities)
end



function integrate_particles(initial_positions::AbstractVector{<:Point}, initial_velocities::AbstractVector{<:Point}, halos, force!; 
         save_step = 5,
         timestep = 0.1, timerange=(0, -10/T2GYR), 
         reuse_acceleration = true,
         verbose = false
	)
	
	check_args(initial_positions, initial_velocities, halos)

	N = length(halos)
	time_i, time_f = timerange

	time_total = time_f - time_i
	itermax = ceil(Int, abs(time_total) / timestep)
    skip = ceil(Int, save_step / timestep)
	num_to_save = ceil(Int, itermax / skip) + 2
	
	all_positions = [Vector{Point}(undef, N) for _ in 1:num_to_save]
	all_velocities = [Vector{Point}(undef, N) for _ in 1:num_to_save]
	all_times = Vector{Float64}(undef, num_to_save)

	dt = timestep * sign(time_total)

	positions = deepcopy(initial_positions)
	velocities = deepcopy(initial_velocities)
	time = time_i

	accelerations = fill(Point(0, 0, 0), N)
	force!(accelerations, positions, velocities, halos, time)

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
            time = time_f
			is_done = true
		end

		# save
		if (i % skip == 0) || is_done	
            if verbose
                @info "completed $i / $itermax steps"
            end

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


function collect_into_orbits(all_positions::Vector{Vector{Point}}, all_velocities::Vector{Vector{Point}}, all_times::Vector)
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


function get_potential(potname; kwargs...)
	Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/", potname * ".ini"); kwargs...)
end

function get_potential_c(potname; kwargs...)
	CAgama.AgamaPotential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/", potname * ".ini"); kwargs...)
end



function get_sigma_v(pot, units)
	gm = Agama.GalaxyModel(pot, Agama.DistributionFunction(pot, type="QuasiSpherical"))

	N = 100
	Rs = 10 .^ LinRange(-1, 3, N)

	sigmas = gm._py.moments(Agama.mat2py([zeros(N) Rs zeros(N)]'), dens=false) |> Agama.py2mat

	sigma3 = sigmas[1, :] .+ sigmas[2, :] + sigmas[3, :]

	sigma1 = sigma3 ./ 3
	interp = LilGuys.lerp(Rs, sqrt.(sigma1) .* Agama.velocity_scale(units))
    return x-> interp(radii(x))
end

function get_rho(pot, units)
	N = 100
	Rs = 10 .^ LinRange(-1, 3, N)

    pos = [zeros(N) Rs zeros(N)]'
    ρs = Agama.density(pot, pos, units)
    interp = LilGuys.lerp(Rs, ρs)
    return x-> interp(radii(x))
end


# Force methods

function potential_acceleration!(accelerations::Vector{Point}, pot::Agama.Potential, positions, time; units)
    accelerations .+= (eachcol(acceleration(pot, hcat(positions...), units, time=time)))
    accelerations
end

function potential_acceleration!(accelerations::Vector{Point}, pot::CAgama.AgamaPotential, positions, time; units)
    for i in eachindex(accelerations)
        accelerations[i] += acceleration(pot, positions[i], units, time=time)
    end
    accelerations
end

function potential_acceleration!(accelerations::Vector{Point}, pot, positions, time)
    for i in eachindex(accelerations)
        accelerations[i] += acceleration(pot, positions[i], time=time)
    end
    accelerations
end


function force_potential_nbody(pot, units)
	function f(accelerations, positions, velocities, halos, time)
		nbody_acceleration!(accelerations, positions, halos)
        potential_acceleration!(accelerations, pot, positions, time; units=units)
	end
end

function force_nbody()
	function f(accelerations, positions, velocities, halos, time)
		nbody_acceleration!(accelerations, positions, halos)
	end
end

function force_potential(pot, units)
	function f(accelerations, positions, velocities, halos, time)
        potential_acceleration!(accelerations, pot, positions, time; units=units)
	end
end


function force_potential(pot)
	function f(accelerations, positions, velocities, halos, time)
        potential_acceleration!(accelerations, pot, positions, time)
	end
end

function force_dyn_friction!(accelerations, positions, velocities, halos; rev=false)

	for i in 1:length(positions)
        dyn_fric = halos[i]
        factor = rev ? -1 : 1
        accelerations[i] += factor * acceleration(dyn_fric, positions[i], velocities[i])
	end

	return accelerations
end

function make_dyn_fric_models(pot_sigma, units, halos; f_σ=nothing, f_ρ=nothing)
    if isnothing(f_σ)
        f_σ = get_sigma_v(pot_sigma, units)
    end
    if isnothing(f_ρ)
        f_ρ = get_rho(pot_sigma, units)
    end

    f_fric = []
    for i in eachindex(halos)
        M = LilGuys.M200(halos[i])
        r_s = halos[i].r_s
        dyn_fric = LilGuys.ChandrashakarDynamicalFriction(r_s=r_s, σv=f_σ, M=M, ρ=f_ρ)
        push!(f_fric, dyn_fric)
    end

    return f_fric
end


function force_potential_nbody_friction(pot, units, halos; pot_sigma=pot, f_σ=nothing, f_ρ=nothing)

    f_fric = make_dyn_fric_models(pot_sigma, units, halos, f_σ=f_σ, f_ρ=f_ρ)
	function f(accelerations, positions, velocities, halos, time)
		nbody_acceleration!(accelerations, positions, halos)
        potential_acceleration!(accelerations, pot, positions, time; units=units)
        force_dyn_friction!(accelerations, positions, velocities, f_fric)
	end
end

function force_potential_friction(pot, units, halos; pot_sigma=pot, f_σ=nothing, f_ρ=nothing)

    f_fric = make_dyn_fric_models(pot_sigma, units, halos, f_σ=f_σ, f_ρ=f_ρ)
	function f(accelerations, positions, velocities, halos, time)
        potential_acceleration!(accelerations, pot, positions, time; units=units)
        force_dyn_friction!(accelerations, positions, velocities, f_fric)
	end
end


"""
    write_orbits(output, orbits; N_max)

Write the first `N_max` orbits to a file "orbits.hdf5" in `output`.
"""
function write_orbits(filename, orbits; N_max=1000)
    N_max = min(length(orbits), N_max)
    structs = [(string(i) => orbit) for (i, orbit) in enumerate(orbits[1:N_max])]

    LilGuys.write_structs_to_hdf5(filename, structs)
end
