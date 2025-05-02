import StatsBase: median
import Random

plot_labels = OrderedDict(
	:xi => L"$\xi$\,/\,degrees",
	:eta => L"$\eta$\,/\,degrees",
	:xi_am => L"$\xi$\,/\,arcmin",
	:eta_am => L"$\eta$\,/\,arcmin",
	:G => "G (mag)",
	:bp_rp => "BP – RP (mag)",
	:pmra => L"$\mu_{\alpha*}$ / mas yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas yr$^{-1}$",
)


function ellipse!(radius, ellipticity, position_angle; kwargs...)
	t = LinRange(0, 2π, 1000) 

	h = 1/sqrt(1 - ellipticity)
	a = radius * h
	b = radius / h

	x = @. a * cos(t) 
	y = @. b * sin(t)

	θ = deg2rad(position_angle)
	x1 = @. x * sin(θ) - y*cos(θ)
	y1 = @. x * cos(θ) + y*sin(θ)

	rs = LilGuys.calc_R_ell(x1, y1, ellipticity, position_angle)

	@assert all(rs .≈ radius)
	lines!(x1, y1; kwargs...)
end




SELECTION_COLORS = [colorant"#cecece"; to_colormap(:YlGnBu_5)[end-2:end]]

function compare_j24_samples(datasets, scatter_kwargs, observed_properties) 
	fig = Figure(
		size = (4*72, 6*72),
	)

	# tangent
    dθ = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))
	ax = Axis(fig[1, 1], 
		xlabel=plot_labels[:xi_am], 
		ylabel=plot_labels[:eta_am], 
		#aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false,
		xreversed = true
	)


	for (label, df) in datasets
		scatter!(df.xi, df.eta; scatter_kwargs[label]...)
	end

	ellipse!(3observed_properties["r_h"], observed_properties["ellipticity"], observed_properties["position_angle"], color=:black, linewidth=1)
	text!(-4observed_properties["r_h"], 0, text=L"3R_h", color=:black)

	
	axislegend(position=:lt)


	grid_low = fig[2, 1] = GridLayout()
	# cmd
	ax =  Axis(grid_low[1,1], 
		yreversed=true,
		xlabel=plot_labels[:bp_rp],
		ylabel=plot_labels[:G],
		limits=(-0.5, 2.5, 15, 21),
		xticks = 0:2,
		xgridvisible=false,
		ygridvisible=false,
		#aspect = 1,
	)

	for (label, df) in datasets
		scatter!(df.bp_rp, df.phot_g_mean_mag; scatter_kwargs[label]...)
	end

	df = datasets[:members]
	errorscatter!([-0.125], [15.75], xerror=[median(df.dBP .+ df.dRP)], yerror=[median(df.dG)], color=:black, size=0, linewidth=1)

	# proper motions
	ax = Axis(grid_low[1,2],
		xlabel = plot_labels[:pmra],
		ylabel = plot_labels[:pmdec],
		#aspect=DataAspect(),
		limits=(-10, 10, -10, 10),
		xgridvisible=false,
		ygridvisible=false,
		)
	
	for (label, df) in datasets
		scatter!(df.pmra, df.pmdec; scatter_kwargs[label]...)
	end

	df = datasets[:members]
	errorscatter!([-7.5], [7.5], xerror=[median(df.pmra_error)], yerror=[median(df.pmdec_error)], color=:black, size=0, linewidth=1)


	rowsize!(grid_low, 1, Aspect(1, 1))

	rowsize!(fig.layout, 1, Aspect(1, 1))

	#resize_to_layout!(fig)
	fig
end
