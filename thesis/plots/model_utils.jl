import TOML
using LilGuys
using CairoMakie
using PyFITS

using Arya
import DensityEstimators: histogram2d
import LinearAlgebra: normalize

smallfontsize = @lift 0.8*$(theme(:fontsize))
smalllinewidth = @lift 0.5 * $(theme(:linewidth))
logdensityrange = (-8, 2)

function get_R_h(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))
	R_h = obs_props["R_h"]
end

function get_r_b(galaxyname, modelname, starsname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/stars/$starsname/")

	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")

	σv = prof_f.annotations["sigma_v"]
	if lmc
		props = TOML.parsefile(model_dir * "../../orbital_properties_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "../../orbital_properties.toml")
	end

	dist_f =  TOML.parsefile(model_dir * "../../orbital_properties.toml")["distance_f"]

	
	dt = props["t_last_peri"]
	r_b = LilGuys.break_radius(σv / V2KMS, dt / T2GYR)
	R_h = get_R_h(galaxyname)

	return LilGuys.kpc2arcmin(r_b, dist_f)
end


function get_r_j(galaxyname, modelname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/")

	if lmc
		props = TOML.parsefile(model_dir * "jacobi_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "jacobi.toml")
	end

	return props["r_J"] 
end


function load_expected_density_profile(galaxyname; scale_by_R_h=false)
	prof = SurfaceDensityProfile(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/density_profiles/fiducial_profile.toml")
		
	prof =  LilGuys.filter_empty_bins(prof)

    if scale_by_R_h
        R_h = get_R_h(galaxyname)
        prof = LilGuys.scale(prof, 1/R_h, 1)
    end
	prof
end


function load_stellar_profiles(galaxyname, modelname, starsname; norm_shift=0)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/stars/$starsname/")

	prof_i = SurfaceDensityProfile(model_dir * "initial_profile.toml") |> LilGuys.filter_empty_bins
	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml") |> LilGuys.filter_empty_bins

	prof_expected = load_expected_density_profile(galaxyname)

	dy = get_normalization(prof_f) + norm_shift .- get_normalization(prof_expected)

	dy = middle(dy)
	prof_i = LilGuys.scale(prof_i, 1, 10^dy)
	prof_f = LilGuys.scale(prof_f, 1, 10^dy)
	
	return prof_i,  prof_f, dy
end


function get_normalization(prof)
	return -LilGuys.mean(prof.log_Sigma[1:5])
end


function compare_both(galaxyname, modelname::String, star; norm_shift=0, lmc=false, r_j=nothing, title="", r_max=4*60, break_height=nothing, kwargs...)
    # setup
	r_b = get_r_b(galaxyname, modelname, star, lmc=lmc)
    @info "break radius $r_b"
	R_h = get_R_h(galaxyname)
	if !isnothing(r_j)
		r_j = get_r_j(galaxyname, modelname, lmc=lmc)
	end

	prof_i, prof_f, norm = load_stellar_profiles(galaxyname, modelname, star, norm_shift=norm_shift)
    Σ0 = -get_normalization(load_expected_density_profile(galaxyname)).middle
    logdensityrange = [Σ0 - 9, Σ0 + 1]
    if isnothing(break_height)
        break_height = logdensityrange[1]
    end

    @info "density range $logdensityrange"

	fig = Figure()

	p = plot_stars_2d(fig[1,2], galaxyname, modelname, star, r_b=r_b, norm=norm, colorrange=logdensityrange)
	hideydecorations!()

    #xlims!(r_max, -r_max)
    #ylims!(r_max, -r_max)
	p = plot_stars_2d(fig[1,1], galaxyname, modelname, star, initial=true,norm=norm, colorrange=logdensityrange)

    t_i = get_time_ini(galaxyname, modelname)
	Colorbar(fig[1,3], p, label="log surface density", ticks=Makie.automatic)
	compare_profiles(fig[2,1:3], prof_i, prof_f, r_b; galaxyname=galaxyname, r_j=r_j, t_i=t_i, logdensityrange=logdensityrange, break_height=break_height, kwargs...)


	Makie.Label(fig[0, :], title)
	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 2, Aspect(1, 1.5))
	resize_to_layout!(fig)
	fig
end



function compare_profiles(halo::String, orbit::String, star; lmc=false, norm_shift=0, r_j=nothing, kwargs...)
	prof_i, prof_f, dy= load_profile(halo, orbit, star; norm_shift=norm_shift)
	r_b = get_r_b(halo, orbit, star, lmc=lmc)
	if !isnothing(r_j)
		r_j = get_r_j(halo, orbit, lmc=lmc)
	end
	compare_profiles(prof_i, prof_f, r_b; galaxyname=galaxyname, r_j=r_j, kwargs...)
end


function compare_profiles(gs, prof_i, prof, r_b; galaxyname=nothing, t_i, break_height=0, r_j=nothing, plot_final=true, logdensityrange=logdensityrange) 
	ax = Axis(gs, 
		xlabel="log Radius / arcmin",
		ylabel = "log surface density",
        limits=((-0.5 , log10(240)), tuple(logdensityrange...)),
	)
	
	lines!(prof_i.log_R, prof_i.log_Sigma, color=COLORS[3], linestyle=:dot,
			label=L"initial ($t = %$(round(t_i, digits=1))$\,Gyr)")
	if plot_final
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[3], linestyle=:solid,
				label=L"final ($t = 0.0\,$Gyr)")
	end

    plot_expected_profile!(galaxyname)
    plot_r_break_arrow!(r_b, break_height)
    plot_r_jacobi_arrow!(r_j, break_height)
    if !isnothing(galaxyname)
        R_h = get_R_h(galaxyname)
        plot_R_h_arrow!(R_h, logdensityrange[1])
    end

	axislegend(position=:lb)
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
					linewidth=theme(:linewidth)[]/2,  
					style=Ann.Styles.LineArrow(head = Ann.Arrows.Head()),
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
		
    annotation!(0, 36, log10(R_h), height, color=:grey, linewidth=theme(:linewidth)[]/2, text=L"R_h")
end


function plot_r_jacobi_arrow!(r_j, break_height)
	if !isnothing(r_j)
		annotation!(0, 30, log10(r_j), break_height,
					linewidth=theme(:linewidth)[]/2,
					style=Ann.Styles.LineArrow(head = Ann.Arrows.Head()),
				   )

		text!(log10(r_j), break_height, 
			  text="Jacobi", 
			  rotation=π/2, 
			  offset=(0., 30.), 
			  align=(:left, :center),
			  fontsize=0.8 * theme(:fontsize)[]
			 )

	end
end

function compare_profiles(prof_i, prof, r_b; kwargs...)
	fig = Figure()
	compare_profiles(fig[1,1], prof_i, prof, r_b; kwargs...)
	return fig
end

function get_time_ini(galaxyname, modelname)
    model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname")
	
	TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))["t_f_gyr"] 
end

function get_stars_final(galaxyname, modelname, starsname, filename="final.fits")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/stars/$starsname/")

	return read_fits(joinpath(modeldir, filename))
end


function plot_stars_2d(gs, stars; R_h=NaN, R_h_label=false, 
        bins=100, r_max=4*60, colormap=Reverse(:Greys), r_b=nothing, orbit_direction=nothing, norm=0, colorrange=logdensityrange, position_multiplier=60)

	ax = Axis(gs, 
			  xlabel = L"\xi \, / \, \textrm{arcmin}", 
			  ylabel = L"\eta \, / \, \textrm{arcmin}",
			  xreversed = true,
			)

	
    bins = LinRange(-r_max, r_max, bins)
	h = histogram2d(stars.xi*position_multiplier, stars.eta*position_multiplier, bins, weights=stars.weights * 10^norm, normalization=:density)

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
    orbit_direction = initial ? nothing :  get_orbit_direction(galaxyname, modelname) 

    return plot_stars_2d(gs, stars; R_h=R_h, R_h_label=initial, orbit_direction=orbit_direction, kwargs...)
end

function plot_R_h_circle!(R_h::Real; label=false)
	arc!((0,0), 6R_h, 0, 2π, color=:white, linewidth=theme(:linewidth)[]/2)
	if label
        text!(-6R_h, 0, offset=(smallfontsize[]/2, 0), text=L"6R_h", align=(:left, :center), color=:white, smallfontsize)
	end
end

function plot_R_h_circle!(galaxyname; initial=false)
	R_h = get_R_h(galaxyname)
    plot_R_h_circle!(galaxyname, label=initial)
end

function plot_r_b_circle!(r_b)
	if isnothing(r_b)
        return
    end

    arc!((0, 0), r_b, 0, 2π, color=COLORS[3], linestyle=:dash, linewidth=smalllinewidth)
    text!(r_b, 0, text="break", color=COLORS[3], fontsize=smallfontsize, offset=(-smallfontsize[]/4, 0), align=(:center, :bottom), rotation=π/2)
end


function get_orbit_direction(galaxyname, modelname)
    model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname")
	orbit = LilGuys.Orbit(joinpath(model_dir, "centres.hdf5"))

	idx_f = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))["idx_f"]
	gcs = Galactocentric.(orbit.positions[1, :], orbit.positions[2, :], orbit.positions[3, :], V2KMS*orbit.velocities[1, :], V2KMS*orbit.velocities[2, :], V2KMS*orbit.velocities[3, :])
	icrs = LilGuys.transform.(ICRS, gcs)


	gsr = LilGuys.transform(GSR, icrs[idx_f])
	dx, dy = (gsr.pmra, gsr.pmdec) ./ sqrt(gsr.pmra^2 + gsr.pmdec^2)
    @info "proper motion: $(gsr.pmra), $(gsr.pmdec)"
    @info "last motion: $((icrs[idx_f].ra - icrs[idx_f-1].ra)*cosd(icrs[idx_f].dec)), $(icrs[idx_f].dec - icrs[idx_f-1].dec)"
    return dx, dy
end



function plot_orbit_arrow!(ax, dx, dy, r_max; r_label=36, kwargs...)
    dx, dy = normalize([dx, dy])
    dx_label = @lift $(ax.xreversed) ? dx * r_label : -dx * r_label
    dy_label = @lift $(ax.yreversed) ? dy * r_label : -dy * r_label

	annotation!(dx_label, dy_label, r_max*dx, r_max*dy, color=COLORS[1])
end
