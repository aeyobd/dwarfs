using Arya
using Makie
using Measurements
using KernelDensity

red = COLORS[6]


"""
    @savefig name [fig]

Saves the curent figure (assumed to be `fig` if not passed to the directory `fig_dir` with a pdf extension
"""
macro savefig(name, fig=nothing)
	if fig === nothing
		fig = esc(:fig)
	end
	
	return quote
		filename = joinpath(fig_dir, $name) * ".pdf"
		save(filename, $fig)
		@info "saved figure to $filename"
	end
end


"""
    xieta_axis(gridspec)

Creates a tangent plane axis to plot xi/eta coordinates into
assuming the units are in degrees.
"""
function xieta_axis(gs; kwargs...)
	return Axis(gs;
		xlabel=L"$\xi$ / degrees",
		ylabel=L"$\eta$ / degrees",
		aspect=DataAspect(),
		xgridvisible=false, 
		ygridvisible=false,
        xreversed=true,
        kwargs...
	)
end



"""
    pm_axis(gridspec; dpm=11, kwargs...)

Creates a proper motion space axis (mura versus mudec) 
assuming units will be mas/yr.
Extra kwargs are passed to Axis(), and dpm controls the 
limits of the axis (assumed to be centred at zero).
"""
function pm_axis(gp; dpm=11, kwargs...)
	return Axis(gp;
		xlabel=L"$\mu_{\alpha*}$ / mas\,yr$^{-1}$",
		ylabel=L"$\mu_\delta$ / mas\,yr$^{-1}$",
		aspect=DataAspect(),
		limits= dpm .* (-1, 1, -1, 1),
		xgridvisible=false,
		ygridvisible=false,
        kwargs...
	)
end


"""
    cmd_axis(gridspec)

Creates a Gaia CMD axis at the given gridspec.
"""
function cmd_axis(gs)
	return Axis(gs,
		xlabel = "Bp-Rp",
		ylabel = "G",
		yreversed=true,
		limits = (-0.2, 2, 15, 21),
		xgridvisible=false,
		ygridvisible=false,
	)
end

"""
    plot_tangent!(grid, all_stars, members=nothing; markersize=2, kwargs...)

Plots the tangent plane coordinates of all stars in `all_stars` on the grid `grid`.
If members is not nothing, it also plots the members in red
"""
function plot_tangent!(grid, all_stars, members=nothing; markersize=2, kwargs...)
    ax = Axis(grid, 
        xlabel = "xi / degrees",
        ylabel = "eta / degrees",
    )
        
    scatter!(ax, all_stars.xi, all_stars.eta, markersize=2,
        color=(:black, 0.2))

    if !isnothing(members )
        scatter!(ax, members.xi, members.eta; 
        markersize=markersize, color=red, kwargs...)
    end
    
    ax.xgridvisible = false
    ax.ygridvisible = false
    
    return ax
end



function plot_pms!(grid, all_stars, members=nothing; da=10, markersize=2, kwargs...)
    ax = Axis(grid, limits=(-da, da, -da, da), aspect=1,
        xlabel=L"\mu_{\alpha *} / \mathrm{ mas\, yr^{-1}}", 
        ylabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}")
    scatter!(ax, all_stars.pmra, all_stars.pmdec, 
        color=(:black, 0.2), markersize=1)

    if !isnothing(members)
        scatter!(ax, members.pmra, members.pmdec, 
            color=red, markersize=markersize; kwargs...)
    end

    ax
end


function plot_cmd!(grid, all_stars, members=nothing; markersize=2, kwargs...)
    ax = Axis(grid, aspect=1,
        limits=(-0.5, 3, 15, 22), yreversed=true,
        xlabel="bp - rp", 
        ylabel="G",)
    scatter!(ax, all_stars.bp_rp, all_stars.phot_g_mean_mag, 
        color=(:black, 0.2), markersize=1)

    if !isnothing(members)
        scatter!(ax, members.bp_rp, members.phot_g_mean_mag;
            color=red, markersize=markersize, kwargs...)
    end
    
    ax.xgridvisible = false
    ax.ygridvisible = false

    ax
end

function plot_parallax!(grid, all_stars, members=nothing; markersize=2, kwargs...)
    da = 15
    ax = Axis(grid, aspect=1, limits=(-10, 10, 0, 4),
        xlabel=L"\varpi / \mathrm{ mas }", 
        ylabel=L"\delta\varpi / \mathrm{ mas }")
    scatter!(ax, all_stars.parallax, all_stars.parallax_error, 
        color=(:grey, 0.2), markersize=1)

    if !isnothing(members)
        scatter!(ax, members.parallax, members.parallax_error; 
            color=red, markersize=markersize, kwargs...)
    end
    ax
end


function plot_density!(grid::GridPosition, all_stars, sample=nothing)
    ax = Axis(grid,
        xlabel="log radius / arcmin",
        ylabel=L"\log\;\Sigma \, / \, \textrm{stars arcmin^{-2}}",
        limits=((-1., 2.3), (-4, 2.2))
    )
    
    plot_density!(ax, all_stars, yerr=false, color=:black)

    if sample !== nothing
        y_end = plot_density!(ax, sample, color=red)
        N = length(sample.r_ell)

        hlines!(value.(y_end))    
        
        text!(ax, 0.1, 0.1, space=:relative, 
            text="$N stars ")
    
        text!(ax, -1, value.(y_end), text=        
            L"\log\Sigma_\textrm{bg} = %$y_end", fontsize=14)
        
    end

    return ax
end


function plot_density!(ax::Axis, sample; bins=nothing, yerr=true, kwargs...)
    x, y, y_err = get_density(sample, bins=bins)
    
    if !yerr
        y_err = nothing
    end
    errscatter!(ax, x, y, yerr=y_err; kwargs...)

    y_end = (y .± yerr)[end-3:end]
    y_end = minimum(y_end)
end


function get_density(df; bins=nothing)
    r = df.r_ell
    props = lguys.StellarProfile(r, bins=bins, normalization=:none)

    println("stars left ", length(r))
    println("counts in last bin ", props.counts[end-2: end])
    println("densities ", (props.log_Sigma .± props.log_Sigma_err)[end-5:end])
    
    return props.log_r, props.log_Sigma, props.log_Sigma_err
end


function plot_all(all_stars, members=nothing; markersize=2, kwargs...)
    fig = Figure()
    plot_tangent!(fig[1, 1], all_stars, members; markersize=markersize, kwargs...)
    plot_pms!(fig[1, 2], all_stars, members; markersize=markersize, kwargs...)
    plot_cmd!(fig[2, 1], all_stars, members; markersize=markersize, kwargs...)
    plot_parallax!(fig[2, 2], all_stars, members; markersize=markersize, kwargs...)
    return fig
end


function isocontours!(x, y; x_bw, y_bw, scale=:log, density_scale=1, kwargs...)
    k = kde((x, y), bandwidth = (x_bw, y_bw))

    if scale == :log
        density = log10.(k.density .+ density_scale)
    elseif scale == :linear
        density = k.density ./ density_scale
    elseif scale == :asinh
        density = asinh.(k.density ./ density_scale)
    else
        throw(ArgumentError("scale must be :log, :linear or :asinh"))
    end

    return contour!(k.x, k.y, density; kwargs...)
end
