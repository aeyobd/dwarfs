import TOML
using Arya
using Makie
using FITSIO
import LilGuys as lguys
import Polyhedra
using DataFrames: DataFrame, names, rename!
using LinearAlgebra: norm, dot

# types used in this code.
OptF = Union{Float64, Nothing}
F = Float64


"""
The parameters used to calculate density from a dataframe
of Gaia observations.
"""
Base.@kwdef struct GaiaFilterParams
    "The name of the datafile to load"
    filename::String

    "The right ascension of the centre of the galaxy"
    ra::F
    "The declination of the centre of the galaxy"
    dec::F
    "The ellipticity of the galaxy"
    ellipticity::F
    "The position angle of the galaxy in degrees (North to East or + in ra)"
    PA::F
    "The maximum angular distance from the centre of the galaxy in degrees"
    max_ang_dist::OptF = nothing
    "Whether to filter stars by the maximum complete r_ell bin"
    filt_r_max::Bool = true

    "The distance to the galaxy in kpc. If the distance is not set, then `n_sigma_dist` sets a simple cut on parallax / error"
    dist::OptF = nothing
    "The error in the distance to the galaxy in kpc"
    dist_err::OptF = nothing
    "The number of sigma which the parallax must be consistent with the distance"
    n_sigma_dist::OptF = nothing

    "The minimum value of PSAT to include in the density calculation"
    PSAT_min::OptF = nothing

    PSAT_col::String = "PSAT"

    "The maximum value of Gaia's renormalised unit weight error"
    ruwe_max::OptF = nothing
    "The minimum value of G magnitude"
    g_min::OptF = nothing
    "The maximum value of G magnitude"
    g_max::OptF = nothing
    "A flattened list of tuples of the vertices of the CMD cut polygon in (bp_rp, g) space"
    cmd_cut::Union{Array,Nothing} = nothing

    "x dimension of CMD"
    cmd_x = "bp_rp"
    "CMD magnitude column"
    cmd_y = "phot_g_mean_mag"
    "The first column to subtract from the second to get the CMD x"
    cmd_x1 = nothing
    "The second column to subtract from the first to get the CMD x"
    cmd_x2 = nothing

    "The pmra mean of the galaxy"
    pmra::OptF = nothing
    "The pmdec mean of the galaxy"
    pmdec::OptF = nothing
    "The assumed intrinsic proper motion error"
    dpm::OptF = nothing
    """The number of sigma to include stars in the pm filter.
    If dpm is set but not this, then the pm filter is a simple distance cutoff
    """
    n_sigma_pm::OptF = nothing

    remove_qso::Bool = false
    remove_galaxy::Bool = false
end



function GaiaFilterParams(dict::Dict; kwargs...)
    d2 =  NamedTuple{Tuple(Symbol.(keys(dict)))}(values(dict))

    return GaiaFilterParams(; d2..., kwargs...)
end

function GaiaFilterParams(filename::String; kwargs...)
    return GaiaFilterParams(read_paramfile(filename); kwargs...)
end


"""
    read_gaia_stars(filename, params; θ=nothing)

Given a fits file, reads the data andd adds columns for the tangent plane
(`xi`, `eta`, in degrees) and the elliptical radius `r_ell` (in arcmin).
If these columns already exist, they are renamed to `xi_original`, `eta_original`, and `r_ell_original`.

If θ is specified, than the orbital coordinate frame is also calculated as xi_p, eta_p.
"""
function read_gaia_stars(filename, params; θ=nothing)
    df = lguys.read_fits(filename)

    for col in ["r_ell", "xi", "eta"]
        if col ∈ names(df)
            df[!, Symbol(col * "_original")] = df[!, col]
        end
    end

    add_xi_eta!(df, params.ra, params.dec)
    r_ell = 60lguys.calc_r_ell(df.xi, df.eta, params.ellipticity, params.PA)
    df[!, :r_ell] = r_ell

    if params.dist !== nothing
        add_pm_gsr!(df, params.dist)
    end

    if θ !== nothing
    	df[:, :xi_p], df[:, :eta_p] = lguys.to_orbit_coords(df.ra, df.dec, 
                                                            params.ra, params.dec, θ)
    end

    if params.PSAT_col !== "PSAT"
        df[!, :PSAT] = df[!, params.PSAT_col]
    end

    if params.cmd_y !== "phot_g_mean_mag"
        df[!, :G] = df[!, params.cmd_y]
    else
        df[!, :G] = df.phot_g_mean_mag
    end

    if params.cmd_x1 !== nothing
        df[!, :bp_rp] = df[!, params.cmd_x1] - df[!, params.cmd_x2]
    elseif params.cmd_x !== "bp_rp"
        df[!, :bp_rp] = df[!, params.cmd_x]
    end


    return df
end



"""
    read_paramfile(filename)

Reads a TOML param file and includes inhereitance.
If the filename includes a key `inherits`, then the file specified by that key is also read and merged, with any duplicate keys being overwritten by the primary `filename` file.
"""
function read_paramfile(filename::String)
    f = TOML.parsefile(filename)

    if "inherits" ∈ keys(f)
        if dirname(filename) == ""
            f1 = read_paramfile(f["inherits"])
        else
            f1 = read_paramfile(dirname(filename) * "/" * f["inherits"])
        end
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
    add_xi_eta!(stars, ra0, dec0)

Given ra0, dec0, 
calculates and adds xi and eta to a `stars` DataFrame
"""
function add_xi_eta!(stars, ra0::Real, dec0::Real)
	xi, eta = lguys.to_tangent(stars.ra, stars.dec, ra0, dec0)
	
	stars[:, "xi"] = xi
	stars[:, "eta"] = eta
	stars
end


"""
    add_pm_gsr!(all_stars, distance=1)

Adds the GSR proper motions to a DataFrame of stars
In theory, distance should not matter ?
"""
function add_pm_gsr!(all_stars, distance=1)
    icrs = [lguys.ICRS(ra=all_stars.ra[i], dec=all_stars.dec[i], distance=distance, 
        pmra=all_stars.pmra[i], pmdec=all_stars.pmdec[i], radial_velocity=0.)
        for i in 1:size(all_stars, 1)]

    gsr = lguys.transform.(lguys.GSR, icrs)

    all_stars[:, :pmra_gsr] = getproperty.(gsr, :pmra)
    all_stars[:, :pmdec_gsr] = getproperty.(gsr, :pmdec)
end



"""
    apply_filter(stars, func, params...)

Given a `stars` DataFrame, applies the filter function func
passing params... along (provided no params are non) 
and pretty prints some notes
"""
function apply_filter(df::DataFrame, func::Function, params...)
	if any(params .== nothing)
		filt = trues(size(df, 1))
	else
		filt = func(df, params...)
	end
    println("filter $func cuts \t", sum(.!filt))
	return filt
end


"""
    select_members(all_stars, params)

Given a DataFrame of all stars and a set of parameters,
applies a series of filters to select the members of the galaxy
"""
function select_members(all_stars, params::GaiaFilterParams)

    filters = [
        psat_filter,
        g_filter,
        ruwe_filter,
        cmd_filter,
        ang_dist_filter,
        parallax_filter,
        pm_filter,
        r_ell_filter,
        qso_filter,
        galaxy_filter,
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

function ruwe_filter(all_stars, params::GaiaFilterParams)
    apply_filter(all_stars, ruwe_filter, params.ruwe_max)
end


function g_filter(all_stars, g_min, g_max)
    return (all_stars.G .> g_min) .& (all_stars.G .< g_max)
end

function g_filter(all_stars, params::GaiaFilterParams)
    return apply_filter(all_stars, g_filter, params.g_min, params.g_max)
end


function qso_filter(all_stars)
    filt = .!all_stars.in_qso_candidates
    return filt
end


function qso_filter(all_stars, params::GaiaFilterParams)
    if params.remove_qso
        return apply_filter(all_stars, qso_filter)
    else
        return trues(size(all_stars, 1))
    end
end

function galaxy_filter(all_stars)
    filt = .!all_stars.in_galaxy_candidates
    return filt
end

function galaxy_filter(all_stars, params::GaiaFilterParams)
    if params.remove_galaxy
        return apply_filter(all_stars, galaxy_filter)
    else
        return trues(size(all_stars, 1))
    end
end


function cmd_filter(all_stars, cmd_cut)
	cmd_cut_m = reshape(cmd_cut, 2, :)
	filt_cmd = is_point_in_polygon.(zip(all_stars.bp_rp, all_stars.G), [cmd_cut_m])
end


function cmd_filter(all_stars, params::GaiaFilterParams)
    return apply_filter(all_stars, cmd_filter, params.cmd_cut)
end


function r_ell_filter(all_stars, ra0, dec0, ellipticity, PA)
    r_ell_max = 60*calc_r_max(all_stars.ra, all_stars.dec, ellipticity, PA, centre=(ra0, dec0))

    println("max r_ell = ", r_ell_max)
    return all_stars.r_ell .< r_ell_max
end


function r_ell_filter(all_stars, params::GaiaFilterParams)
    if params.filt_r_max
        return apply_filter(all_stars, r_ell_filter, params.ra, params.dec, params.ellipticity, params.PA)
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

function ang_dist_filter(all_stars, params::GaiaFilterParams)
    return apply_filter(all_stars, ang_dist_filter, params.ra, params.dec, params.max_ang_dist)
end



function pm_filter(all_stars, pmra, pmdec, dpm, n_sigma)
    σx = all_stars.pmra_error .⊕ dpm
    σy = all_stars.pmdec_error .⊕ dpm

    Δx = pmra .- all_stars.pmra
    Δy = pmdec .- all_stars.pmdec

    dist = @. ⊕(Δx/σx, Δy/σy)
    return dist .< n_sigma
end


function pm_simple_filter(all_stars, pmra, pmdec, dpm)
    δx = pmra .- all_stars.pmra
    δy = pmdec .- all_stars.pmdec
    dist = @. δx ⊕ δy
    return dist .< dpm
end


function pm_filter(all_stars, params::GaiaFilterParams)
    if params.n_sigma_pm === nothing
        return apply_filter(all_stars, pm_simple_filter, params.pmra, params.pmdec, params.dpm)
    else
        return apply_filter(all_stars, pm_filter, params.pmra, params.pmdec, params.dpm, params.n_sigma_pm)
    end
end



function psat_filter(all_stars, psat_min)
    println("number nan PSAT         \t", sum(isnan.(all_stars.PSAT)))
    println("number exactly zero PSAT\t", sum(all_stars.PSAT .== 0))
    println("number > zero           \t", sum(all_stars.PSAT .> 0))
    println("number == 1             \t", sum(all_stars.PSAT .== 1))

    println("total                   \t", length(all_stars.PSAT))
	return all_stars.PSAT .> psat_min
end

function psat_filter(all_stars, params::GaiaFilterParams)
    return apply_filter(all_stars, psat_filter, params.PSAT_min)
end


function parallax_filter(all_stars, dist, dist_err, n_sigma_dist)
	parallax = 1/dist
	parallax_err = 1/dist * dist_err / dist_err

	sigma = all_stars.parallax_error .⊕ parallax_err
	
	filt_parallax = @. (
    abs(all_stars.parallax - parallax) <  sigma * n_sigma_dist
    )
end

function parallax_simple_filter(all_stars, n_sigma_dist)
    return @. abs(all_stars.parallax) < n_sigma_dist * all_stars.parallax_error
end


function parallax_filter(all_stars, params::GaiaFilterParams)
    if params.dist === nothing
        return apply_filter(all_stars, parallax_simple_filter, params.n_sigma_dist)
    else
        return apply_filter(all_stars, parallax_filter, params.dist, params.dist_err, params.n_sigma_dist)
    end
end



"""
    is_point_in_polygon(point, polygon)

returns if the given point is inside the polygon specified as a 2xN matrix of points in 2D.
"""
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



@doc raw"""
	bivariate_normal(x, y, μx, μy, σx, σy, ρ)

A bivariate normal distribution on input vectors x, y assuming a mean (μx, μy), standard deviations of (σx, σy), and a correlation ρ ∈ [-1, 1].

``
\frac{1}{2\pi \sigma_x \sigma_y \sqrt{1 - \rho^2}} \exp\left[-\frac{1}{2(1-\rho^2)}\left(\frac{(x-\mu_x)^2}{\sigma_x^2} + \frac{(y-\mu_y)^2}{\sigma_y^2} - 2\rho \frac{(x-\mu_x)(y-\mu_y)}{\sigma_x\sigma_y}\right)\right]
``

"""
function bivariate_normal(x, y, μx, μy, σx, σy, ρ)
	A = 1 / (2π * σx * σy * √(1 - ρ^2))

	zx = (x - μx) / σx
	zy = (y - μy) / σy

	return A * exp(-1/(2*(1 - ρ^2)) * (
		zx^2 + zy^2 - 2ρ * zx * zy
	))
end


@doc raw"""
    bivariate_z(x, y, μx, μy, σx, σy, ρ)

The z-score of a bivariate normal distribution on input vectors x, y assuming a mean (μx, μy), standard deviations of (σx, σy), and a correlation ρ ∈ [-1, 1].

``
z = \frac{1}{1 - \rho^2} \left(z_x^2 + z_y^2 - 2\rho z_x z_y\right)
``

where `z_x = (x - μ_x) / σ_x` and `z_y = (y - μ_y) / σ_y`.

"""
function bivariate_z(x, y, μx, μy, σx, σy, ρ)
	zx = (x - μx) / σx
	zy = (y - μy) / σy

	return 1/(1 - ρ^2) * (
		zx^2 + zy^2 - 2ρ * zx * zy
	)
end



@doc raw"""
	⊕(x, y)

Add x and y in quadrature. (i.e. standard propogation of errors by addition)

``
x \oplus y \equiv \sqrt{x^2 + y^2}
``
"""
function ⊕(x::Real, y::Real)
	return sqrt(x^2 + y^2)
end


"""
    calc_r_max(ra, dec, args...; centre="mean", weights=nothing)

    Calculates the approximate maximum elliptical radius which is complete
in a set of points in ra, dec space.
"""
function calc_r_max(ra, dec, args...; 
        centre="mean",
        weights=nothing
    )

    ra0, dec0 = lguys.calc_centre2D(ra, dec, centre, weights)

    x, y = lguys.to_tangent(ra, dec, ra0, dec0)
    x_p, y_p = lguys.shear_points_to_ellipse(x, y, args...)

    if length(args) == 3
        a, b, _ = args
        aspect = b/a
    else
        aspect = lguys.ellipticity_to_aspect(args[1])
    end
    if aspect < 1
        aspect = 1/aspect
    end
    
    hull = convex_hull(x_p, y_p)
    r_max = min_distance_to_polygon(hull...)
        # @warn "r_ell: Convex hull not defined. Using max radius. Load Polyhedra to enable convex hull."
    r_max_simple = maximum(@. sqrt(x_p^2 + y_p^2)) ./ sqrt(aspect)
    @info "r_ell: r_max = $r_max, r_max_simple = $r_max_simple"

    return r_max
end

function min_distance_to_polygon(x, y)
    min_dist = Inf
    N = length(x)
    for i in 1:N
        a = [x[i], y[i]]
        j = mod1(i + 1, N)
        b = [x[j], y[j]]

        dist = distance_to_segment(a, b)

        min_dist = min(min_dist, dist)
    end

    return min_dist
end

    
"""
    distance_to_segment(a, b, p)

Distance from point `p` to the line segment defined by `a` and `b`.
all points are 2D vectors.
"""
function distance_to_segment(a, b, p=zeros(2))
    a = vec(a)
    b = vec(b)

    # work in origin at p
    a -= p
    b -= p

    # is the segment a point?
    l = norm(a - b)
    if l == 0
        return norm(a)  
    end

    # line unit vector
    n = (a - b) / l
    # projection along line
    t = dot(a, n) 

    if t < 0
        closest_point = a
    elseif t > l
        closest_point = b
    else
        closest_point = a - t * n
    end

    dist = norm(closest_point)

    return dist
end



"""
    convex_hull(x, y)
Given a vector of x and y coordinates, returns
the convex hull bounding the points.
"""
function convex_hull(xi, eta)
    @assert length(xi) == length(eta)
    filt = .!isnan.(xi) .& .!isnan.(eta)
    @assert sum(filt) > 2

    ps = [[x, e] for (x, e) in zip(xi[filt], eta[filt])]
    p = Polyhedra.convexhull(ps...)
    b = Polyhedra.planar_hull(p).points.points
    return first.(b), last.(b)
end
