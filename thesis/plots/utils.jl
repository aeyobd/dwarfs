import StatsBase: median
import Random
using Makie
using OrderedCollections
using CSV, DataFrames
using LilGuys
using Arya
import TOML


"""
    get_obs_props(galaxyname::String)

Retrieve the observed properties of the galaxy (from the observations/ directory) as a Dictionary
"""
function get_obs_props(galaxyname::String) 
    return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))
end


"""
    get_M_star(galaxyname)

Calculate the stellar mass (in code units) of the galaxy name based
on the Mv magnitude and M_L_s ratio in the observed properties
"""
function get_M_star(galaxyname::String)
    obs_props = get_obs_props(galaxyname)
	M_star = LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"]
	return M_star / M2MSUN
end


"""
    ellipse!(radius::Real, ellipticity::Real, position_angle::Real; x0, y0, kwargs...)

Plot an ellipse of the specified radius, ellipticity, and position angle, 
centred at x0 and y0. kwargs... passed to `lines!`
"""
function ellipse!(radius::Real, ellipticity::Real, position_angle::Real; 
        x0::Real=0, y0::Real=0, kwargs...)
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

Given an axis (with pre-specified limits?), compute
the factor to wich to multiply a slope to convert 
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
		Δx_data = @lift $(ax.finallimits).widths[1]
		Δy_data = @lift $(ax.finallimits).widths[2]
	end
	Δx_screen = @lift $(ax.scene.viewport).widths[1]
	Δy_screen = @lift $(ax.scene.viewport).widths[2]
	correction_factor = @lift $Δx_data / $Δx_screen * $Δy_screen / $Δy_data
	return correction_factor
end


"""
    log_derivative(f, x0::Real; h::Real)

Compute the derivative of f (a function) at x0 in log-log space. 
"""
function log_derivative(f, x0::Real; h::Real=0.001)
	y0 = f(x0)
	return (log10(f(x0 * 10^h)) - log10(y0)) / h
end



"""
    text_along_line_log!(x, y, x_0; text, h, kwargs...)

Given vectores of x and y positions, plot text following the x-y curve (in log-log space) at a position x_0.
"""
function text_along_line_log!(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, x_0::Real; text::AbstractString, h::Real=0.03, kwargs...)
	f = LilGuys.lerp(x, y)
	y_0 = f(x_0)
	dy = log_derivative(f, x_0, h=h)
	rf = rotation_factor(Makie.current_axis(), true)
	θ = @lift atan($rf * dy)

	text!(x_0, y_0, text=text, rotation=θ; kwargs...)
end
