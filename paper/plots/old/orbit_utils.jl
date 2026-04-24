using CairoMakie
using LilGuys: mean
using LilGuys
using Arya
import TOML

strokewidth = @lift $(theme(:linewidth)) / 2

"""
    plot_rt!(ax, times, rs; plot, kwargs...)

Plots the radii at each given time onto the specified axes.
"""
function plot_rt!(axes, times, rs; plot=lines!, kwargs...)
    plot(axes, times, rs; kwargs...)
end


"""
    plot_rt!(ax, orbit(s); plot, kwargs...)

Plots the radii for each timestep in the orbit(s). 
"""
function plot_rt!(axes, orbit::Orbit; time_min=-Inf, kwargs...)
    filt = orbit.times  .> time_min
    plot_rt!(axes, orbit.times[filt] * T2GYR, radii(orbit)[filt]; kwargs...)
end


function plot_rt!(axes, orbits::AbstractVector{<:Orbit}; rasterize=true, time_min=-Inf, kwargs...)
    rs = Float64[]
    ts = Float64[]
    for orbit in orbits
        filt = orbit.times .> time_min
        rs = vcat(rs, radii(orbit)[filt], NaN)
        ts = vcat(ts, orbit.times[filt]*T2GYR, NaN)
    end

    plot_rt!(axes, ts, rs; rasterize=rasterize, kwargs...)
end



"""
    plot_present_position!(axes, positions, velocities)

Plot the present position as a dot point.
"""
function plot_xyz_today!(axes, orbit::Orbit, index=1; kwargs...)
    x0 = orbit.positions[:, index]
    v0 = orbit.velocities[:, index]
    plot_xyz_today!(axes, x0; kwargs...)
end

function plot_xyz_today!(axes, orbits, index=1; kwargs...)
    xyz_mean = LilGuys.mean([o.positions[:, index] for o in orbits])

    return plot_xyz_today!(axes, xyz_mean; kwargs...)
end


"""
    plot_xyz_today!(axes, x0; kwargs...)

Plots the xy, yz, and zx position of the provided 3-vecotr onto a 3-tuple of xyz axes slices.
"""
function plot_xyz_today!(axes, x0::AbstractVector{<:Real}; strokewidth=strokewidth, kwargs...)
    ax_xy, ax_yz, ax_xz = axes
	i = 1
	j = 2
    scatter!(ax_xy, x0[i], x0[j]; strokewidth=strokewidth, kwargs...)

	i = 2
	j = 3
    scatter!(ax_yz, x0[i], x0[j]; strokewidth=strokewidth, kwargs...)

	i = 1
	j = 3
    scatter!(ax_xz, x0[i], x0[j]; strokewidth=strokewidth, kwargs...)
	axes
end


function plot_rt_today!(ax_rt, orbit::Orbit, index=1; kwargs...)
    x0 = orbit.positions[:, index]
    v0 = orbit.velocities[:, index]
    r0 = radii(x0)
    t0 = orbit.times[index] * T2GYR

    plot_rt_today!(ax_rt, t0, r0; kwargs...)
end


function plot_rt_today!(ax_rt, orbits, index=1; kwargs...)
    rs = [radii(o.positions[:, index]) for o in orbits]
    r0 = mean(rs)
    t0 = orbits[1].times[index] * T2GYR
    return plot_rt_today!(ax_rt, t0, r0; kwargs...)
end


function plot_rt_today!(ax_rt, t0::Real, r0::Real; strokewidth=strokewidth, kwargs...)
    scatter!(ax_rt, t0, r0; strokewidth=strokewidth, kwargs...)
end



function plot_xyz_sun!(axes; kwargs...)
    X_SUN = [-8.1219733661223, 0.0, 0.0208]

    plot_xyz_today!(axes, X_SUN; marker=:star5, color=COLORS[9], kwargs...)
end

function load_best_orbit(galaxyname, modelname)
    modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)

	best_orbit = LilGuys.Orbit(joinpath(modeldir, "centres.hdf5"))
		
	idx_f = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))["idx_f"]
	best_orbit.times .-= best_orbit.times[idx_f]

	return best_orbit[best_orbit.times .<= 0]
end
