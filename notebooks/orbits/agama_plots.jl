using LilGuys
using CairoMakie



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



function plot_r_t(orbits; legend=true, kwargs...)
    fig = Figure()
    ax = axis_r_t(fig[1, 1])

    for (label, orbit) in orbits
        plot_r_t!(ax, orbit, label=label; kwargs...)
    end

    if legend
        axislegend()
    end

    return fig
end


function plot_y_z(orbits; legend=true, kwargs...)
    fig = Figure()
    ax = axis_y_z(fig[1, 1]; kwargs...)

    for (label, orbit) in orbits
        plot_y_z!(ax, orbit, label=label)
    end

    if legend
        Legend(fig[1,2], ax)
    end
    return fig
end

function plot_xyz(orbits::Vector{<:Pair{<:Any, Orbit}}; kwargs...)
    pos = getproperty.(last.(orbits), :position)
    labels = first.(orbits)

    return LilGuys.Plots.plot_xyz(pos...; labels=labels, kwargs...)
end


function ax_v_circ(gp; kwargs...)
    Axis(gp;
         xlabel = "log r / kpc",
         ylabel = "v_circ / km/s",
         kwargs...
        )
end

function plot_v_circ!(pot; vasiliev_units = false, log_r=LinRange(-1, 2.5, 100), log=true, kwargs...)
    r = 10 .^ log_r
    v_circ = calc_v_circ_pot(pot, r, vasiliev_units=vasiliev_units)

    if log
        x = log_r
    else
        x = r
    end
    p = lines!(x, v_circ * V2KMS; kwargs...)
    return p
end
