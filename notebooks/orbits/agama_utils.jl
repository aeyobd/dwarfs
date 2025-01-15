using SpecialFunctions: erf
using PythonCall
agama = pyimport("agama")
u = pyimport("astropy.units")

using LilGuys
using CairoMakie

import Base: -, reverse
potential_dir = ENV["DWARFS_ROOT"] * "/agama/potentials/"


# vasiliev units
V_V2KMS = 1
V_M2MSUN = 232_500
V_R2KPC = 1
V_T2GYR = 0.97779

F = Float64

"""
Represents an orbit.
"""
Base.@kwdef struct Orbit
    time::Vector{F}
    position::Matrix{F}
    velocity::Matrix{F}
    acceleration::Union{Matrix{F}, Nothing} = nothing
    del_a::Union{Matrix{F}, Nothing} = nothing
    pericenter::F = minimum(calc_r(position))
    apocenter::F = maximum(calc_r(position))
end


function (-)(a::Orbit, b::Orbit)
    @assert isapprox(a.time, b.time, rtol=1e-6)

    return Orbit(time=a.time, position=a.position - b.position, velocity=a.velocity - b.velocity)
end


function reverse(a::Orbit)
    return Orbit(time=reverse(a.time), position=reverse(a.position, dims=2), velocity=reverse(a.velocity, dims=2))
end

function resample(a::Orbit, time::AbstractVector)
    if a.time[2] < a.time[1]
        a = reverse(a)
    end
    x = LilGuys.lerp(a.time, a.position[1, :]).(time)
    y = LilGuys.lerp(a.time, a.position[2, :]).(time)
    z = LilGuys.lerp(a.time, a.position[3, :]).(time)
    v_x = LilGuys.lerp(a.time, a.velocity[1, :]).(time)
    v_y = LilGuys.lerp(a.time, a.velocity[2, :]).(time)
    v_z = LilGuys.lerp(a.time, a.velocity[3, :]).(time)

    return Orbit(time=time, position=[x y z]', velocity=[v_x v_y v_z]')
end

function load_agama_potential(filename)
    return agama.Potential(joinpath(potential_dir, filename))
end


"""
    calc_orbit(phase, potential; N=10_000)

Given an inital phase position and agama potential, computes the orbit.
Time is in specified units.
"""
function calc_orbit(coords, pot; N=10_001, time=10, units = :code)
    ic = make_agama_init(coords, units=units)

    if time isa Real
        o = agama.orbit(ic=ic, time=time, potential=pot, trajsize=N)
        pytime = o[0]
        pyposvel = o[1]
        time = pyconvert(Vector{Float64}, pytime)
    elseif time isa AbstractVector
        time0 = time[1]
        tottime = time[end] - time[1]
        o = agama.orbit(ic=ic, time=tottime, timestart=time0, potential=pot, dtype="object")
        pyposvel = o(time)
    end
    return from_agama_orbit(time, pyposvel, units=units)
end


function calc_orbits(coords::AbstractVector, pot; N=10_000, time=10, units=:code)
    ic = make_agama_init(coords, units=units)

    Np = length(coords)

    if time isa Real
        o = agama.orbit(ic=ic, time=time, potential=pot, trajsize=N)
        time = pyconvert(Vector{Float64}, o[0][0])
        pyposvels = [o[i][1] for i in 0:Np-1]
    elseif time isa AbstractVector
        time0 = time[1]
        tottime = time[end] - time[1]
        o = agama.orbit(ic=ic, time=tottime, timestart=time0, potential=pot, dtype="object")
        pyposvels = [oo(time) for oo in o]
    end
    return ["$i" => from_agama_orbit(time, pyposvels[i], units=units) for i in eachindex(pyposvels)]
end



function from_agama_orbit(time, posvel; units=:code)
    posvel = pyconvert(Matrix{Float64}, posvel)'

    position = posvel[1:3, :]
    velocity = posvel[4:6, :]

    if units == :physical
        velocity = velocity ./ V2KMS
        time = time ./ T2GYR
    elseif units == :vasiliev
        velocity = velocity ./ V2KMS
        time = time ./ T2GYR
    end
    return Orbit(time=time, position=position, velocity=velocity)
end

function make_agama_init(coord::LilGuys.Galactocentric; units=:code)
    x = coord.x
    y = coord.y
    z = coord.z
    v_x = coord.v_x
    v_y = coord.v_y
    v_z = coord.v_z
    if units == :code
        v_x /= V2KMS
        v_y /= V2KMS
        v_z /= V2KMS
    elseif units == :vasiliev
        v_x /= V_V2KMS
        v_y /= V_V2KMS
        v_z /= V_V2KMS
    elseif units == :physical
        # pass
    end

    return [x, y, z, v_x, v_y, v_z]
end


function make_agama_init(coords::AbstractVector{<:LilGuys.Galactocentric}; units=:code)
    x = [coord.x for coord in coords]
    y = [coord.y for coord in coords]
    z = [coord.z for coord in coords]
    v_x = [coord.v_x for coord in coords]
    v_y = [coord.v_y for coord in coords]
    v_z = [coord.v_z for coord in coords]

    if units == :code
        v_x ./= V2KMS
        v_y ./= V2KMS
        v_z ./= V2KMS

    elseif units == :vasiliev
        v_x ./= V_V2KMS
        v_y ./= V_V2KMS
        v_z ./= V_V2KMS
    elseif units == :physical
        # pass
    end
    return [x y z v_x v_y v_z]
end

function calc_v_circ_pot(pot, r; vasiliev_units = false)
    r = vec(r)
    M = pot.enclosedMass(r) 
    M = pyconvert(Vector{F}, M)

    v_circ = @. sqrt(M / r)

    if vasiliev_units
        v_circ *= V_V2KMS / V2KMS
    end


    return v_circ
end




@doc raw"""
    a_dyn_friction(pos, vel; r_s, σv, ρ, M, Λ)

Computes the dynamical friction based on 

``
\frac{d{\bf v}}{dt} = -\frac{4\pi\,G^2\,M\,\rho\,\ln\Lambda}{v^2} \left({\rm erf}(X) - \frac{2X}{\sqrt\pi} \exp(-X^2)\right) \frac{{\bf v}}{v}
``
"""
function a_dyn_friction(pos, vel; r_s=nothing, σv, ρ, M, Λ=nothing  )
    G = 1
    v = calc_r(vel)
    r = calc_r(pos)
    
    if Λ === nothing
        if r_s < 8
            ϵ = 0.45r_s
        else
            ϵ = 2.2r_s - 14
        end

        Λ = r / ϵ
    end

    X = v / (√2 * σv(r))
    fX = erf(X) - 2X/√π * exp(-X^2)

    return -4π*G^2 * M * ρ(pos) * max(log(Λ), 0) * fX * vel ./ v^3
end


"""
    calc_σv_interp(pot; log_r = LinRange(-3, 5, 100000))

Computes the velocity dispersion as a function of radius by integrating the Jeans equation.

"""
function calc_σv_interp(pot; log_r = LinRange(-3, 5, 100000))
    x0 = [1/√2, 0., 1/√2] # direction

    if !issorted(log_r)
        @warn "log_r is not sorted. Sorting."
        log_r = sort(log_r)
    end

    log_r = reverse(log_r)
    radii = 10 .^ log_r

    positions = x0' .* radii
    acc = pyconvert(Matrix{Float64}, pot.force(positions))
    rho = pyconvert(Array{Float64}, pot.density(positions))

    a = calc_r(acc')
    
    dr = abs.(LilGuys.gradient(radii))

    σ2 = 1 ./ rho .* cumsum(rho .* a .* dr)

    return LilGuys.lerp(radii, sqrt.(σ2))
end



"""
    leap_frog(gc, acceleration; params...)

Computes the orbit of a Galactocentric object using the leap frog method.

Parameters:
- `gc::Galactocentric`: the initial conditions
- `acceleration::Function`: the acceleration function
- `dt_max::Real=0.1`: the maximum timestep
- `dt_min::Real=0.001`: the minimum timestep, will exit if timestep is below this
- `time::Real=-10/T2GYR`: the time to integrate to
- `timebegin::Real=0`: the time to start integrating from
- `timestep::Symbol=:adaptive`: the timestep to use, either `:adaptive` or a real number
- `η::Real=0.01`: the adaptive timestep parameter

"""
function leap_frog(gc, acceleration; 
        dt_max=0.1, dt_min=0.001, timebegin=0, time=-10/T2GYR, 
        timestep=:adaptive, η=0.01
    )

    if timestep isa Real
        dt_min = timestep
    end

    t = timebegin

    Nt = round(Int, abs(time / dt_min))
    positions = Vector{Vector{Float64}}()
    velocities = Vector{Vector{Float64}}()
    accelerations = Vector{Vector{Float64}}()
    times = Float64[]

    push!(positions, [gc.x, gc.y, gc.z])
    push!(velocities, [gc.v_x, gc.v_y, gc.v_z] / V2KMS)
    push!(accelerations, acceleration(positions[1], velocities[1], t))
    push!(times, 0.)
    is_done = false
    backwards = time < 0

    for i in 1:Nt
        pos = positions[i]
        vel = velocities[i]
        acc = acceleration(pos, vel, t)
        
        if timestep == :adaptive
            dt = min(sqrt(η / calc_r(acc)), dt_max)
        elseif timestep isa Real
            dt = timestep
        end
        
        if backwards
            dt *= -1
        end
        if abs(dt) < dt_min
            @warn "timestep below minimum timestep"
            break
        end

        if (backwards && t + dt <= time ) || (!backwards && t + dt >= time)
            dt = time - t
            is_done = true
        end

        vel_h = vel + 1/2 * dt * acc
        pos_new = pos + dt*vel_h
        acc = acceleration(pos_new, vel_h, t + dt/2)
        vel_new = vel_h + 1/2 * dt * acc

        push!(positions, pos_new)
        push!(velocities, vel_new)
        push!(accelerations, acc)
        t = times[i] + dt
        push!(times, t)

        if is_done
            break
        end
    end

    positions_matrix = hcat(positions...)
    velocities_matrix = hcat(velocities...)
    accelerations_matrix = hcat(accelerations...)
    return Orbit(time=times, position=positions_matrix, velocity=velocities_matrix, acceleration=accelerations_matrix)
end
