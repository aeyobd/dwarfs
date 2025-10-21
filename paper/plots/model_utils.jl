import TOML
using LilGuys
using CairoMakie
using PyFITS

using Arya
import DensityEstimators: histogram2d
import LinearAlgebra: normalize

smallfontsize = @lift 0.8*$(theme(:fontsize))
smalllinewidth = @lift 0.5 * $(theme(:linewidth))


"""
    compare_both(galaxyname, modelname, starsname, <keyword arguments)

Makes a comparison plot of the 2D densities of the model 
before and after tidal evolution,
and the initial and final density as compared with the observed density profile

# Arguments
- `norm_shift`: Artificially shift the calculated normalization of the density profile in order to better agree with data
- `lmc`: calculate break and jacobi radii with respect to the lmc
- `title`: A plot title
- `kwargs...`: passed to compare profiles.

"""
function compare_both(galaxyname::String, modelname::String, star::String; 
        norm_shift::Real=0, lmc::Bool=false, title::String="", 
        kwargs...)
    # setup
	r_b = get_r_b(galaxyname, modelname, star, lmc=lmc)
    @info "break radius $r_b"

	prof_i, prof_f, norm = load_stellar_profiles(galaxyname, modelname, star, 
                                                 norm_shift=norm_shift)

    Σ0 = -get_normalization(load_expected_density_profile(galaxyname))
    logdensityrange = [Σ0 - 9, Σ0 + 1]

    @info "density range $logdensityrange"

    # plots
    fig = Figure(size=(3.5*72, 4.5*72))

	p = plot_stars_2d(fig[1,1], galaxyname, modelname, star, 
                        r_b=r_b, norm=norm, colorrange=logdensityrange)

    fig.content[1].aspect = DataAspect()
    hidedecorations!()
    Colorbar(fig[1,2], p, label="log surface density", ticks=Makie.automatic)

    # Density profiles
	compare_profiles(fig[2,:], prof_i, prof_f, r_b; 
                     galaxyname=galaxyname, 
                     modelname=modelname,
                     lmc=lmc,
                     logdensityrange=logdensityrange, 
                     kwargs...)

	Makie.Label(fig[0, :], title)
	#rowsize!(fig.layout, 1, Aspect(1, 1.0))
	fig
end


"""
    compare_profiles(gs, prof_i, prof, r_b; 
    galaxyname, time_i, break_height, r_j, plot_final, logdensityrange)

"""
function compare_profiles(gs, prof_i, prof, r_b; 
        logdensityrange, galaxyname=nothing, break_height=0, 
        r_j=false, plot_final=true, modelname, lmc=false,
    )
	ax = Axis(gs, 
		xlabel = "log Radius / arcmin",
		ylabel = "log surface density",
        limits = ((-0.5 , log10(240)), tuple(logdensityrange...)),
	)

    t_i = get_time_ini(galaxyname, modelname)
	
	lines!(prof_i.log_R, prof_i.log_Sigma, 
           color=COLORS[3], linestyle=:dot,
			label=L"initial ($t = %$(round(t_i, digits=1))$\,Gyr)")

	if plot_final
		lines!(prof.log_R, prof.log_Sigma, 
               color=COLORS[3], linestyle=:solid, 
               label=L"final ($t = 0.0\,$Gyr)")
	end

    plot_expected_profile!(galaxyname)

    # Add annotations
    if isnothing(break_height)
        break_height = logdensityrange[1]
    end
    plot_r_break_arrow!(r_b, break_height)

	if r_j
		r_j = get_r_j(galaxyname, modelname, lmc=lmc)
        @info "r_j = $r_j"
    else
        r_j = nothing
    end

    plot_r_jacobi_arrow!(r_j, break_height)
    if !isnothing(galaxyname)
        R_h = get_R_h(galaxyname)
        plot_R_h_arrow!(R_h, logdensityrange[1])
    end

	axislegend(position=:lb)
end


function compare_profiles(prof_i, prof, r_b; kwargs...)
	fig = Figure()
	compare_profiles(fig[1,1], prof_i, prof, r_b; kwargs...)
	return fig
end


function plot_expected_profile!(galaxyname)
    if !isnothing(galaxyname)
        prof_expected = load_expected_density_profile(galaxyname)

        errorscatter!(prof_expected.log_R, prof_expected.log_Sigma,
                      yerror = error_interval.(prof_expected.log_Sigma),
            label="observed",
            color=:black
        )
    end
end


function plot_r_break_arrow!(r_b, break_height)
	if !isnothing(break_height)
		annotation!(0, 30, log10(r_b), break_height, 
					color=COLORS[3], 
				   )
		
		text!(log10(r_b), break_height, 
			  text="break", 
			  rotation=π/2, 
			  offset=(0., 30.), 
			  color=COLORS[3], 
			  align=(:left, :center), 
			  fontsize=0.8 * theme(:fontsize)[]
			)
	end
end

function plot_R_h_arrow!(R_h, height)
	if isnothing(R_h)
        return
    end
		
    annotation!(0, 30, log10(R_h), height, color=:grey,)

    text!(log10(R_h), height, color=:grey, text=L"R_h",
          fontsize=smallfontsize, offset=(0, 30), align=(:center, :bottom))
end


function plot_r_jacobi_arrow!(r_j, break_height)
	if isnothing(r_j)
        return
    end

    annotation!(0, 30, log10(r_j), break_height)

    text!(log10(r_j), break_height, 
          text="Jacobi", 
          rotation=π/2, 
          offset=(0., 30.), 
          align=(:left, :center),
          fontsize=0.8 * theme(:fontsize)[]
         )
end


function plot_stars_2d(gs, stars; 
        R_h=NaN, R_h_label=false, 
        bins=100, r_max=4*60, colormap=Reverse(:Greys), r_b=nothing, 
        orbit_direction=nothing, norm=0, colorrange, position_multiplier=60)

	ax = Axis(gs, 
        xlabel = L"\xi \, / \, \textrm{arcmin}", 
        ylabel = L"\eta \, / \, \textrm{arcmin}",
        xreversed = true,
    )

    bins = LinRange(-r_max, r_max, bins)
    xi_am = stars.xi * position_multiplier
    eta_am = stars.eta * position_multiplier

	h = histogram2d(xi_am, eta_am, bins, 
        weights=stars.weights * 10^norm, normalization=:density)


    @info "density maximum: $(log10(maximum(h.values)))"

	p = heatmap!(h.xbins, h.ybins, log10.(h.values), 
                 colorrange=colorrange, colormap=colormap)

    plot_R_h_circle!(R_h, label=R_h_label)
    plot_r_b_circle!(r_b)

    if !isnothing(orbit_direction)
        dx, dy = orbit_direction
        plot_orbit_arrow!(ax, dx, dy, r_max)
    end

    return p
end


function plot_stars_2d(gs, galaxyname, modelname, starsname; initial=false, kwargs...)
    filename = initial ? "initial.fits" : "final.fits"
	stars = get_stars_final(galaxyname, modelname, starsname, filename)
    R_h = get_R_h(galaxyname)
    orbit_direction = initial ? nothing : get_orbit_direction(galaxyname, modelname) 

    return plot_stars_2d(gs, stars; R_h=R_h, R_h_label=initial, orbit_direction=orbit_direction, kwargs...)
end


function plot_R_h_circle!(R_h::Real; label=false)
	arc!((0,0), 6R_h, 0, 2π, color=:white, linewidth=smalllinewidth)
	if label
        text!(-6R_h, 0, offset=(smallfontsize[]/2, 0), text=L"6R_h", align=(:left, :center), color=:white, fontsize=smallfontsize)
	end
end


function plot_R_h_circle!(galaxyname::String; label=false)
	R_h = get_R_h(galaxyname)
    plot_R_h_circle!(R_h, label=label)
end


function plot_r_b_circle!(r_b)
	if isnothing(r_b)
        return
    end

    arc!((0, 0), r_b, 0, 2π, 
         color=COLORS[3], linestyle=:dash, linewidth=smalllinewidth)

    text!(r_b, 0, text="break", 
          color=COLORS[3], fontsize=smallfontsize, offset=(-smallfontsize[]/4, 0),
          align=(:center, :bottom), rotation=π/2)
end


function plot_orbit_arrow!(ax, dx, dy, r_max; r_label=30, kwargs...)
    dx, dy = normalize([dx, dy])
    dx_label = @lift $(ax.xreversed) ? dx * r_label : -dx * r_label
    dy_label = @lift $(ax.yreversed) ? dy * r_label : -dy * r_label

	annotation!(dx_label, dy_label, r_max*dx, r_max*dy, color=COLORS[1])
end


# Data Loading functions


function load_expected_density_profile(galaxyname)
	prof = SurfaceDensityProfile(joinpath(ENV["DWARFS_ROOT"],
        "observations", galaxyname, "density_profiles/fiducial_profile.toml"))
		
	prof = LilGuys.filter_empty_bins(prof)

	prof
end


function load_stellar_profiles(galaxyname, modelname, starsname; norm_shift=0)
    modeldir = get_starsdir_out(galaxyname, modelname, starsname)

    prof_i = SurfaceDensityProfile(joinpath(modeldir, "initial_profile.toml")) |> LilGuys.filter_empty_bins
    prof_f = SurfaceDensityProfile(joinpath(modeldir, "final_profile.toml")) |> LilGuys.filter_empty_bins

	prof_expected = load_expected_density_profile(galaxyname)

	dy = get_normalization(prof_f) + norm_shift - get_normalization(prof_expected)

	prof_i = LilGuys.scale(prof_i, 1, 10^dy)
	prof_f = LilGuys.scale(prof_f, 1, 10^dy)
	
	return prof_i, prof_f, dy
end


function get_normalization(prof)
	return -LilGuys.mean(prof.log_Sigma[1:5]) |> middle
end


"""
    get_R_h(galaxyname::String)

get the value of the (circularized) half-light radius of the galaxy.
"""
function get_R_h(galaxyname::String)
    filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
                        "observed_properties.toml")
    obs_props = TOML.parsefile(filename)
	R_h = obs_props["R_h"]
end


"""
    get_r_b(galaxyname, modelname, starsname; lmc=false)


Get the break radius of the specified model
"""
function get_r_b(galaxyname::String, modelname::String, starsname::String; lmc=false)
    modeldir = get_starsdir_out(galaxyname, modelname, starsname)
    prof_f = SurfaceDensityProfile(joinpath(modeldir, "final_profile.toml"))
	σv = prof_f.annotations["sigma_v"]

	if lmc
        props = TOML.parsefile(joinpath(modeldir, "../../orbital_properties_lmc.toml"))
	else
        props = TOML.parsefile(joinpath(modeldir, "../../orbital_properties.toml"))
	end

	dt = props["t_last_peri"]
	r_b = LilGuys.break_radius(σv / V2KMS, dt / T2GYR)

    # same with or without lmc
    dist_f =  TOML.parsefile(joinpath(modeldir, "../../orbital_properties.toml"))["distance_f"]
	
	return LilGuys.kpc2arcmin(r_b, dist_f)
end


"""
    get_r_j(galaxyname, modelname; lmc=false)

Get the jacobi radius of the specified model
"""
function get_r_j(galaxyname::String, modelname::String; lmc=false)
    modeldir = get_modeldir(galaxyname, modelname)

	if lmc
        props = TOML.parsefile(joinpath(modeldir, "jacobi_lmc.toml"))
	else
        props = TOML.parsefile(joinpath(modeldir, "jacobi.toml"))
	end

	return props["r_J"] 
end


function get_modeldir(galaxyname::String, modelname::String)
    modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)
end


function get_starsdir_out(galaxyname::String, modelname::String, starsname::String)
    modeldir = joinpath(get_modeldir(galaxyname, modelname), "stars", starsname)
end


function get_time_ini(galaxyname, modelname)
    modeldir = get_modeldir(galaxyname, modelname)
	TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))["t_f_gyr"] 
end

function get_stars_final(galaxyname, modelname, starsname, filename="final.fits")
    starsdir = get_starsdir_out(galaxyname, modelname, starsname)
	return read_fits(joinpath(starsdir, filename))
end



function get_orbit_direction(galaxyname, modelname)
    modeldir = get_modeldir(galaxyname, modelname)
	orbit = LilGuys.Orbit(joinpath(modeldir, "centres.hdf5"))
	idx_f = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))["idx_f"]

    gsr = get_gsr_position(orbit, idx_f)
    sanity_check_orbit_pm(orbit, idx_f)
	return (gsr.pmra, gsr.pmdec)
end


function sanity_check_orbit_pm(orbit, idx_f)
    # sanity checks
    gsr = get_gsr_position(orbit, idx_f)
	dx, dy = (gsr.pmra, gsr.pmdec) ./ sqrt(gsr.pmra^2 + gsr.pmdec^2)
    @info "proper motion: $dx, $dy"
    gsr_old = get_gsr_position(orbit, idx_f-1)
    dx_test = cosd(gsr.dec) * (gsr.ra - gsr_old.ra)
    dy_test = gsr.dec - gsr_old.dec
    dx_test, dy_test = (dx_test, dy_test) ./ sqrt(dx_test^2 + dy_test^2)

    @info "last motion: $dx_test, $dy_test"
end

function get_gsr_position(orbit, idx_f)
    pos_f = orbit.positions[:, idx_f]
    vel_f = orbit.velocities[:, idx_f] 
    gc = Galactocentric(pos_f, vel_f * V2KMS)
    gsr = LilGuys.transform(GSR, gc)
end
