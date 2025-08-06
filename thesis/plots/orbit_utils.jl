using CairoMakie
using LilGuys: mean

import LinearAlgebra: ⋅

"""
    plot_present_position!(axes, positions, velocities)

Plot the present position as an arrow pointing at the present velocity
"""
function plot_present_position!(axes, x0::AbstractVector, v0::AbstractVector, )

	ax = axes[1]
	i = 1
	j = 2
	arrowhead!(ax, x0[i], x0[j], v0[i], v0[j]; )

	ax = axes[2]
	i = 2
	j = 3
	arrowhead!(ax, x0[i], x0[j], v0[i], v0[j]; )

	ax = axes[3]
	i = 1
	j = 3
	arrowhead!(ax, x0[i], x0[j], v0[i], v0[j]; )

	axes
end


"""
    arrowhead!(ax, x, y, v_x, v_y; length=0.15)

Plot an arrowhead pointing in dataspace v_x and v_y
"""
function arrowhead!(ax, x, y, v_x, v_y; markersize=1.5 * theme(:markersize)[], notch=0.2, kwargs...)

	(xlims, ylims) = ax.limits[]
	ax_scale = sqrt((xlims[2] - xlims[1])^2 + (ylims[2] - ylims[1]))

    # this is just so arrowhead shows up due to bug
    length = 1
	scale = length*ax_scale
	v_x_norm, v_y_norm = (v_x, v_y) ./ sqrt(v_x^2 + v_y^2)
	
    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=markersize, notch=notch))
        
	annotation!(ax, x - scale*v_x_norm, y - scale*v_y_norm, x, y;
                style=style, labelspace=:data, linewidth=0, shrink=(0., 0.))
end




function axes_xyz_flat(fig, limits, units="kpc")
    ax_xy = Axis(fig[1,1],
        xlabel = "x / $units",
        ylabel = "y / $units",
        limits = limits[[1,2]],
        aspect = DataAspect(),
    )

    ax_yz = Axis(fig[1,2],
        xlabel = "y / $units",
        ylabel = "z / $units",
        limits = limits[[2,3]],
        aspect = DataAspect(),
    )

    ax_xz = Axis(fig[1,3],
        xlabel = "x / $units",
        ylabel = "z / $units",
        limits = limits[[1,3]],
        aspect = DataAspect(),
    )

    rowsize!(fig.layout, 1, Aspect(1,1))

    return [ax_xy, ax_yz, ax_xz]
end


function plot_xyz!(axes, positions::AbstractMatrix;
        plot=lines!, kwargs...)
    ax_xy, ax_yz, ax_xz = axes

    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]

    plot(ax_xy, x, y; kwargs...)
    plot(ax_yz, y, z; kwargs...)
    plot(ax_xz, x, z; kwargs...)
end


function plot_xyz!(axes, orbit::Orbit; time_min=-Inf, kwargs...)
    time_filt = orbit.times .> time_min
    plot_xyz!(axes, orbit.positions[:, time_filt]; kwargs...)
end

function plot_xyz!(axes, orbits::AbstractVector{<:Orbit}; kwargs...)
    for orbit in orbits
        plot_xyz!(axes, orbit; kwargs...)
    end
end



function plot_rt!(axis, times, rs; plot=lines!, kwargs...)
    plot(axis, times, rs; kwargs...)
end


function plot_rt!(axes, orbit::Orbit; kwargs...)
    plot_rt!(axes, orbit.times * T2GYR, radii(orbit); kwargs...)
end

function plot_rt!(axes, orbits::AbstractVector{<:Orbit}; kwargs...)
    for orbit in orbits
        plot_rt!(axes, orbit; kwargs...)
    end
end



"""
    plot_present_position!(axes, positions, velocities)

Plot the present position as an arrow pointing at the present velocity
"""
function plot_xyz_today!(axes, orbit::Orbit, index=1; kwargs...)
    ax_xy, ax_yz, ax_xz = axes
    x0 = orbit.positions[:, index]
    v0 = orbit.velocities[:, index]

	i = 1
	j = 2
	arrowhead!(ax_xy, x0[i], x0[j], v0[i], v0[j]; )

	i = 2
	j = 3
	arrowhead!(ax_yz, x0[i], x0[j], v0[i], v0[j]; )

	i = 1
	j = 3
	arrowhead!(ax_xz, x0[i], x0[j], v0[i], v0[j]; )
	axes
end

function plot_rt_today!(ax_rt, orbit::Orbit, index=1; kwargs...)
    x0 = orbit.positions[:, index]
    v0 = orbit.velocities[:, index]
    r0 = radii(x0)
    vr = (x0 ⋅ v0) / r0 / T2GYR
    t0 = orbit.times[index] * T2GYR
    
    arrowhead!(ax_rt, t0, r0, 1, vr)
end
