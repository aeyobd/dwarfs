using CairoMakie
using LilGuys: mean


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
function arrowhead!(ax, x, y, v_x, v_y; markersize=theme(:markersize)[], notch=0.2, kwargs...)

	(xlims, ylims) = ax.limits[]
	ax_scale = sqrt((xlims[2] - xlims[1])^2 + (ylims[2] - ylims[1]))

    # this is just so arrowhead shows up due to bug
    length = 1
	scale = length*ax_scale
	v_x_norm, v_y_norm = (v_x, v_y) ./ sqrt(v_x^2 + v_y^2)
	
    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(length=markersize, notch=notch))
        
	annotation!(ax, x - scale*v_x_norm, y - scale*v_y_norm, x, y;
                style=style, labelspace=:data, linewidth=0, shrink=(0., -markersize/2), kwargs...)
end
