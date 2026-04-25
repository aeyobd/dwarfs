import StatsBase
using DataFrames, CSV

using CairoMakie
using Arya
using Printf

plotrange = 5

px_per_pt = 4.166666666666666
update_theme!(px_per_unit=1)
CairoMakie.activate!(px_per_unit=px_per_pt, pt_per_unit=1) # 300 DPI

xz_iso = CSV.read(joinpath( @__DIR__,  "resources/EP2020_iso_xz.csv"), DataFrame, ntasks=1)


function get_histogram(snap, bins; weights=nothing)
	x = snap.positions[2, :]
	y = snap.positions[3, :]

	if eltype(bins) <: Real
		bins = (bins, bins)
	end
	
	if weights === nothing                                                              
        h1 = StatsBase.fit(StatsBase.Histogram, (x, y), bins)                                      
    else                                                                          
        h1 = StatsBase.fit(StatsBase.Histogram, (x, y), StatsBase.weights(weights), bins) 
    end         
	
	return bins, StatsBase.normalize(h1, mode=:density).weights
end


function plot_xy_density!(snap, bins; colorrange_dm, colorrange_stars)

	h = get_histogram(snap, bins)
	plot_hist!(h, colorrange=colorrange_dm, colormap = colormap_dm)

	h = get_histogram(snap, bins, weights=snap.weights)
	plot_hist!(h, colorrange=colorrange_stars, colormap = colormap_stars)

end


function plot_hist!(binshist; interpolate=false, kwargs...)
	bins, hist = binshist
	image!(extrema.(bins)..., log10.(hist);  interpolate=interpolate, kwargs...)
end


function to_transparent_cmap(color)
	return ([Makie.RGBAf(color.r, color.g, color.b, alpha) for alpha in LinRange(0, 1., 100)])
end

function to_black_cmap(color)
	return ([Makie.RGBf(color.r*alpha, color.g*alpha, color.b*alpha) for alpha in LinRange(0, 1., 100)])
end

function plot_scalebar!(scale_length=50, plotrange=bins[end]-bins[1]; color=:grey)
	length_relative = scale_length / plotrange

	x0, y0 = 0.05, 0.05
	lines!([x0, x0 + length_relative], [y0, y0], color=color, space=:relative, linewidth=theme(:linewidth)[] / 2)
	text!(x0, y0, text="$scale_length kpc", color=color, space=:relative, fontsize=0.8 * theme(:fontsize)[], )
end


function add_time!(ax, time; fontsize=0.8*theme(:fontsize)[])
    label = @sprintf("today %+2.1f Gyr", time)
    label = replace(label, "-" => "– ") # nicer minus sign
    label = replace(label, "+" => "+ ") # correct spacing

    text!(ax, 0.05, 0.95, text=label, color=:grey, align=(:left, :top),
        fontsize=fontsize, space=:relative)
end


purple = Makie.colorant"#c274ff"

colormap_dm = (to_black_cmap(purple))

colormap_stars = to_transparent_cmap(RGBf(1., 1., 1.))



function plot_frame(snap; colorrange_dm, colorrange_stars,
        scalebar_color=:grey,
        bins = LinRange(-150, 150, 1001),
        legend = false
    )
    N = length(bins) - 1
    if N != 1000
        @warn "using a different number of bins than 1000, textsizes may be wrong"
    end

    fig = Figure(figure_padding=(0,0,0,0), backgroundcolor=:black, px_per_unit=1, pt_per_unit=1/3, size=(N/px_per_pt, N/px_per_pt))


	ax = Axis(fig[1, 1], backgroundcolor=:black, 
			 limits = (extrema(bins), extrema(bins)),
            )


	plot_xy_density!(snap, bins; colorrange_dm=colorrange_dm, colorrange_stars=colorrange_stars)
	plot_scalebar!(50, bins[end]-bins[1], color=scalebar_color)

	# MW isocontour. same if using x or y
	poly!(xz_iso.x, xz_iso.z, color=COLORS[8])
	
	# legend
    if legend
        text!(0.95, 0.05, text="stars", color=:white, space=:relative, align=(:right, :center), offset=(0, 10))
        text!(0.95, 0.05, text="dark matter", color=(purple, 1), space=:relative, align=(:right, :center))
    end


	hidexdecorations!()
	hideydecorations!()

    resize_to_layout!()
	fig
end


function get_zoom_histogram(snap, bins=nothing; weights=nothing, plotrange=plotrange, N = 128)

	xycen = snap.x_cen[2:3]
	
	if isnothing(bins)
		
		bins = (LinRange(xycen[1] - plotrange, xycen[1] + plotrange, N),
				LinRange(xycen[2] - plotrange, xycen[2] + plotrange, N)
			   )
	end

	
	hist = get_histogram(snap, bins, weights=weights)

	return hist
end
