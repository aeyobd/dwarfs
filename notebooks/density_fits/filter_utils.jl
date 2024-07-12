import TOML
using Arya
using Makie

OptF = Union{Float64, Nothing}
F = Float64

"""
The parameters used to calculate density from a dataframe
of Gaia observations.
"""
Base.@kwdef struct DensityParams
    filename::String
    ra::Float64
    dec::Float64
    ellipticity::F
    rh::F
    PA::F
    dist::F
    dist_err::F
    PSAT_min::OptF = nothing
    ruwe_max::OptF = nothing
    g_min::OptF = nothing
    g_max::OptF = nothing
    max_ang_dist::OptF = nothing
    n_sigma_dist::OptF = nothing
    pmra::OptF = nothing
    pmdec::OptF = nothing
    filt_r_max::Bool = true

    dpm::OptF = nothing
    n_sigma_pm::OptF = nothing

    cmd_cut::Union{Array,Nothing} = nothing
end



function DensityParams(dict::Dict; kwargs...)
    d2 =  NamedTuple{Tuple(Symbol.(keys(dict)))}(values(dict))

    return DensityParams(; d2..., kwargs...)
end

"""
    load_stars(filename, params)

Given a fits file, loads the stars and calculates the tangent plane
"""
function load_stars(filename, params)
	all_stars_unfiltered = DataFrame(FITS(filename)[2])
	xi, eta = lguys.to_tangent(all_stars_unfiltered.ra, all_stars_unfiltered.dec, params.ra, params.dec)
	
	r_ell = 60lguys.calc_r_ell(xi, eta, params.ellipticity, params.PA)

    for col in ["r_ell", "xi", "eta"]
        if col ∈ names(all_stars_unfiltered)
            all_stars_unfiltered[!, Symbol(col * "_original")] = all_stars_unfiltered[!, col]
        end
    end

	all_stars_unfiltered[!, :r_ell] = r_ell
	all_stars_unfiltered[!, :xi] = xi
	all_stars_unfiltered[!, :eta] = eta

	return all_stars_unfiltered
end



"""
Reads a TOML param file and includes inhereitance
"""
function read_file(filename)
	f = TOML.parsefile(filename)

	if "inherits" ∈ keys(f)
		f1 = read_file(dirname(filename) * "/" * f["inherits"])
		delete!(f, "inherits")
        f = merge(f1, f)
	end

    for (k, v) in f
        if v === NaN
            f[k] = nothing
        end
    end

	return f
end



"""
given an all_stars dataframe with xi and eta in degrees, 
plots the tangent plane with proper labels
"""
function plot_all_tangent(all_stars; kwargs...)
    fig = Figure()

    plot_all_tangent!(fig[1,1], all_stars; kwargs...)

    return fig
end

function plot_all_tangent!(grid::GridPosition, all_stars; scale=1, units="degrees", r_max=nothing, kwargs...)

    x = scale*all_stars.xi
    y = scale*all_stars.eta

    if r_max === nothing
        r_max = max(maximum(abs.(x)), maximum(abs.(y)))
    end

    ax = Axis(grid,
        xlabel=L"\xi / \textrm{%$units}", ylabel=L"\eta / \textrm{%$units}",
        aspect=1,
        limits=(-r_max, r_max, -r_max, r_max),
        xgridvisible=false, ygridvisible=false,
        xreversed=true,
    )

    plot_all_tangent!(ax, all_stars; scale=scale, kwargs...)

    return ax
end


function plot_all_tangent!(ax::Axis, all_stars; scale=1, kwargs...)

    x = scale*all_stars.xi
    y = scale*all_stars.eta
    

    p = scatter!(ax, x, y; kwargs...)

    return p
end



function load_fits(filename)
	f = FITS(filename)
	all_stars = DataFrame(f[2])
	close(f)
	return all_stars
end



"""
Given ra0, dec0, 
calculates and adds xi and eta to a `stars` DataFrame
"""
function add_xi_eta!(stars, ra0, dec0)
	xi, eta = lguys.to_tangent(stars.ra, stars.dec, ra0, dec0)
	
	stars[:, "xi"] = xi
	stars[:, "eta"] = eta
	stars
end



"""
Given a `stars` DataFrame, applies the filter function func
passing params... along (provided no params are non) 
and pretty prints some notes
"""
function apply_filter(df, func, params...)
	if any(params .== nothing)
		filt = trues(size(df, 1))
	else
		filt = func(df, params...)
	end
	println("filter cuts $func \t", sum(map(!, filt)))
	return filt
end


function min_filter(x, attr, cut)
	return x[:, attr] .> cut
end

max_filter(x, attr, cut) = x[:, attr] .< cut


function select_members(all_stars, params)

    filters = [
        psat_filter,
        g_filter,
        ruwe_filter,
        cmd_filter,
        ang_dist_filter,
        parallax_filter,
        pm_filter,
        r_ell_filter
    ]


    filt = trues(size(all_stars, 1))

    for f in filters
        filt .&= f(all_stars, params)
    end

	println(sum(filt), " stars remaining")
	return all_stars[filt, :]
	
end


function ruwe_filter(all_stars, ruwe_max)
    return all_stars.ruwe .< ruwe_max
end

function ruwe_filter(all_stars, params::DensityParams)
    apply_filter(all_stars, ruwe_filter, params.ruwe_max)
end


function g_filter(all_stars, g_min, g_max)
    return (all_stars.phot_g_mean_mag .> g_min) .& (all_stars.phot_g_mean_mag .< g_max)
end

function g_filter(all_stars, params::DensityParams)
    return apply_filter(all_stars, g_filter, params.g_min, params.g_max)
end


function cmd_filter(all_stars, cmd_cut)
	cmd_cut_m = reshape(cmd_cut, 2, :)
	filt_cmd = is_point_in_polygon.(zip(all_stars.bp_rp, all_stars.phot_g_mean_mag), [cmd_cut_m])
end


function cmd_filter(all_stars, params::DensityParams)
    return apply_filter(all_stars, cmd_filter, params.cmd_cut)
end


function r_ell_filter(all_stars, ra0, dec0, ellipticity, PA)
    r_ell_max = 60*lguys.calc_r_max(all_stars.ra, all_stars.dec, ellipticity, PA, centre=(ra0, dec0))

    println("max r_ell = ", r_ell_max)
    return all_stars.r_ell .< r_ell_max
end


function r_ell_filter(all_stars, params::DensityParams)
    if params.filt_r_max
        return r_ell_filter(all_stars, params.ra, params.dec, params.ellipticity, params.PA)
    else
        return trues(size(all_stars, 1))
    end
end






function ang_dist_filter(all_stars, ra0, dec0, max_ang_dist)
    filt_ang_dist = @. (
        max_ang_dist ^2
        > (all_stars.ra - ra0)^2 * cosd(dec0)^2 
        + (all_stars.dec - dec0)^2
        )
end

function ang_dist_filter(all_stars, params::DensityParams)
    return apply_filter(all_stars, ang_dist_filter, params.ra, params.dec, params.max_ang_dist)
end



function pm_filter(all_stars, pmra, pmdec, dpm, n_sigma)
    σx = @. sqrt(all_stars.pmra_error^2 + dpm^2)
    σy = @. sqrt(all_stars.pmdec_error^2 + dpm^2)

    Δx = pmra .- all_stars.pmra
    Δy = pmdec .- all_stars.pmdec

    dist = @. sqrt(Δx^2/σx^2 + Δy^2/σy^2)
    return dist .< n_sigma
end


function pm_filter(all_stars, params::DensityParams)
    pm_filter(all_stars, params.pmra, params.pmdec, params.dpm, params.n_sigma_pm)
end


function psat_filter(all_stars, psat_min)
    println("number nan PSAT         \t", sum(isnan.(all_stars.PSAT)))
    println("number exactly zero PSAT\t", sum(all_stars.PSAT .== 0))
    println("number > zero           \t", sum(all_stars.PSAT .> 0))
    println("number == 1             \t", sum(all_stars.PSAT .== 1))

    println("total                   \t", length(all_stars.PSAT))
	return all_stars.PSAT .> psat_min
end

function psat_filter(all_stars, params::DensityParams)
    return apply_filter(all_stars, psat_filter, params.PSAT_min)
end


function parallax_filter(all_stars, dist, dist_err, n_sigma_dist)
	parallax = 1/dist
	parallax_err = 1/dist * dist_err / dist_err

	sigma = @. sqrt(all_stars.parallax_error^2 + parallax_err^2)
	
	filt_parallax = @. (
    abs(all_stars.parallax - parallax) <  sigma * n_sigma_dist
    )
end


function parallax_filter(all_stars, params::DensityParams)
    return apply_filter(all_stars, parallax_filter, params.dist, params.dist_err, params.n_sigma_dist)
end


function is_point_in_polygon(point, polygon)
    x, y = point
    inside = false
    N = size(polygon, 2)  # Number of vertices in the polygon
    j = N  # Start with the last vertex
    for i in 1:N
        xi, yi = polygon[1, i], polygon[2, i]
        xj, yj = polygon[1, j], polygon[2, j]
        
        # Check if point intersects with polygon edge
        intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
        if intersect
            inside = !inside
        end
        j = i
    end
    return inside
end
