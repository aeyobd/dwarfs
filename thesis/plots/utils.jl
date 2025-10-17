import StatsBase: median
import Random
using Makie
using OrderedCollections
using CSV, DataFrames
using LilGuys
using Arya
import TOML

function get_obs_props(galaxyname) 
    return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))
end

function get_M_star(galaxyname)
    obs_props = get_obs_props(galaxyname)
	M_star = LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"]
	return M_star / M2MSUN
end

function ellipse!(radius, ellipticity, position_angle; x0=0, y0=0, kwargs...)
	t = LinRange(0, 2π, 1000) 

	h = 1/sqrt(1 - ellipticity)
	a = radius * h
	b = radius / h

	x = @. a * cos(t) 
	y = @. b * sin(t)

	θ = deg2rad(position_angle)
	x1 = @. x * sin(θ) - y*cos(θ) + x0
	y1 = @. x * cos(θ) + y*sin(θ) + y0

	rs = LilGuys.calc_R_ell(x1 .- x0, y1 .- y0, ellipticity, position_angle)

	@assert all(rs .≈ radius)
	lines!(x1, y1; kwargs...)
end


"""
    rotation_factor(ax[, log=false)

Given an axis (with pre-specified limits?), computes
the facto to wich to multiply a slope to convert 
a geometric slope into a screen slope
"""
function rotation_factor(ax, log=false)
	if log
		limits = ax.finallimits
		x1 = @lift $limits.origin[1]
		x2 = @lift $x1 + $limits.widths[1]
		y1 = @lift $limits.origin[2]
		y2 = @lift $y1 + $limits.widths[2]

		Δx_data = @lift log10($x2) - log10($x1)
		Δy_data = @lift log10($y2) - log10($y1)
	else
		Δx_data = @lift $(ax.finallimits).widths[2]
		Δy_data = @lift $(ax.finallimits).widths[2]
	end
	Δx_screen = @lift $(ax.scene.viewport).widths[1]
	Δy_screen = @lift $(ax.scene.viewport).widths[2]
	correction_factor = @lift $Δx_data / $Δx_screen * $Δy_screen / $Δy_data
	return correction_factor
end


function log_derivative(f, x0; h=0.001)
	y0 = f(x0)
	return (log10(f(x0 * 10^h)) - log10(y0)) / h
end



function text_along_line_log!(x, y, x_0; text, h=0.03, kwargs...)
	f = LilGuys.lerp(x, y)
	y_0 = f(x_0)
	dy = log_derivative(f, x_0, h=h)
	rf = rotation_factor(Makie.current_axis(), true)
	θ = @lift atan($rf * dy)

	text!(x_0, y_0, text=text, rotation=θ; kwargs...)
end
