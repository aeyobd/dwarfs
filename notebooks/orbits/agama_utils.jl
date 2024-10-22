using PythonCall
agama = pyimport("agama")
u = pyimport("astropy.units")

using LilGuys
using CairoMakie

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
    pericenter::F = NaN
    apocenter::F = NaN
end


function load_agama_potential(filename)
    return agama.Potential(joinpath(potential_dir, filename))
end


"""
    calc_orbit(phase, potential; N=10_000)

Given an inital phase position and agama potential, computes the orbit.
Returns a ve
"""
function calc_orbit(phase::LilGuys.Galactocentric, pot; N=10_000, time=10, units = :code)

    if units == :code
        ic = [phase.x, phase.y, phase.z, 
            ((phase.v_x, phase.v_y, phase.v_z) ./ V2KMS)...]
        time = time ./ T2GYR
    elseif units == :vasiliev
        ic = [phase.x, phase.y, phase.z, 
            ((phase.v_x, phase.v_y, phase.v_z) ./ V_V2KMS)...]
        time  = time ./ V_T2GYR
    elseif units == :physical
        ic = [phase.x,
              phase.y,
              phase.z,
              phase.v_x,
              phase.v_y,
              phase.v_z]
    end

    if time isa Real
        o = agama.orbit(ic=ic, time=time, potential=pot, trajsize=N)
        pytime = o[0]
        pyposvel = o[1]
        time = pyconvert(Vector{Float64}, pytime)
    elseif time isa AbstractVector
        time0 = time[1]
        tottime = time[end] - time[1]
        println("start time: ", time0)
        println("integration time: ", tottime)
        o = agama.orbit(ic=ic, time=tottime, timestart=time0, potential=pot, dtype="object")
        pyposvel = o(time)
    end
    posvel = pyconvert(Matrix{Float64}, pyposvel)'

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



function axis_r_t(gp; kwargs...)
    Axis(gp;
         xlabel = "time / Gyr",
         ylabel = "r / kpc",
         kwargs...
        )
end


function axis_R_z(gp; kwargs...)
    Axis(gp;
         xlabel = "R / kpc",
         ylabel = "z / kpc",
         aspect = DataAspect(),
         kwargs...
        )
end


function axis_y_z(gp; kwargs...)
    Axis(gp;
         xlabel = "y / kpc",
         ylabel = "z / kpc",
         aspect = DataAspect(),
         kwargs...
        )
end

function plot_R_z!(ax, orbit::Orbit; kwargs...)
    R = @. sqrt(orbit.position[1, :]^2 + orbit.position[2, :]^2)
    z = orbit.position[3, :]

    lines!(ax, R, z; kwargs...)
end


function plot_y_z!(ax, orbit::Orbit; kwargs...)
    y = orbit.position[2, :]
    z = orbit.position[3, :]
    x = orbit.position[1, :]

    x ./= maximum(abs.(x))
    x .+= 1
    x .*= 5

    lines!(ax, y, z; kwargs...)
end


function plot_r_t!(ax, orbit::Orbit; kwargs...)
    r = calc_r(orbit.position)
    t = orbit.time * T2GYR

    lines!(ax, t, r; kwargs...)
end



function plot_y_z(orbit::Orbit; kwargs...)
    fig = Figure()
    ax = axis_y_z(fig[1, 1]; kwargs...)
    p = plot_y_z!(ax, orbit, color=orbit.time)
    return Makie.FigureAxisPlot(fig, ax, p)
end


function plot_R_z(orbit::Orbit; kwargs...)
    fig = Figure()
    ax = axis_R_z(fig[1, 1]; kwargs...)
    p = plot_R_z!(ax, orbit, color=orbit.time)
    return Makie.FigureAxisPlot(fig, ax, p)
end

function plot_r_t(orbit::Orbit; kwargs...)
    fig = Figure()
    ax = axis_r_t(fig[1, 1]; kwargs...)
    p = plot_r_t!(ax, orbit)
    return Makie.FigureAxisPlot(fig, ax, p)
end



function plot_r_t(orbits; kwargs...)
    fig = Figure()
    ax = axis_r_t(fig[1, 1]; kwargs...)

    for (label, orbit) in orbits
        plot_r_t!(ax, orbit, label=label)
    end

    axislegend()
    return fig
end


function plot_y_z(orbits; kwargs...)
    fig = Figure()
    ax = axis_y_z(fig[1, 1]; kwargs...)

    for (label, orbit) in orbits
        plot_y_z!(ax, orbit, label=label)
    end

    axislegend()
    return fig
end


function ax_v_circ(gp; kwargs...)
    Axis(gp;
         xlabel = "log r / kpc",
         ylabel = "v_circ / km/s",
         kwargs...
        )
end

function plot_v_circ(pot; vasiliev_units = false, kwargs...)
    fig = Figure()
    ax = ax_v_circ(fig[1, 1]; kwargs...)

    log_r = LinRange(-2, 2, 100)
    r = 10 .^ log_r
    v_circ = calc_v_circ_pot(pot, r, vasiliev_units=vasiliev_units)

    p = lines!(ax, log_r, v_circ * V2KMS)
    return Makie.FigureAxisPlot(fig, ax, p)
end


# TODO 
# add function to plot Kz (1.1 kpc)
# function to plot Vcirc of potential
# maybe plot isopotential contours
#
#
