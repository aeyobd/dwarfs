import StatsBase: median
import Random
using Makie
using OrderedCollections
using CSV, DataFrames
using LilGuys
using Arya


SELECTION_COLORS = [colorant"#cecece"; to_colormap(:YlGnBu_5)[end-2:end]]
plot_labels = OrderedDict(
	:xi => L"$\xi$\,/\,degree",
	:eta => L"$\eta$\,/\,degree",
	:xi_am => L"$\xi$\,/\,arcmin",
	:eta_am => L"$\eta$\,/\,arcmin",
	:G => L"$G$ (mag)",
    :bp_rp => L"$G_\textrm{BP} – G_\textrm{RP}$ (mag)",
	:pmra => L"$\mu_{\alpha*}$\,/\,mas\,yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$\,/\,mas\,yr$^{-1}$",
)

using PythonCall

Dustmaps = pyimport("dustmaps.sfd")

DUSTMAP = Dustmaps.SFDQuery()

SkyCoord = pyimport("astropy.coordinates").SkyCoord


"""
    get_extinction(ra, dec, bp_rp)

"""
function get_extinction(ra, dec, bp_rp)
    icrss = SkyCoord(ra, dec, unit="degree", frame="icrs")

    ebv = pyconvert(Vector, DUSTMAP(icrss))

    A0=3.1*ebv

    Ag, Ab, Ar = let
        cg=[0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]
        cb=[1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
        cr=[0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]
        terms = hcat(ones(length(ra)), bp_rp, (bp_rp) .^2.,(bp_rp) .^3., A0, A0 .^2., (bp_rp) .*A0)

        kg = terms * (cg)
        kb =  terms * cb
        kr = terms * cr

        A0 .*kg, A0 .*kb, A0 .*kr
    end

    return Ag, Ab, Ar
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





function compare_j24_samples(datasets, scatter_kwargs, observed_properties; 
        age=12, legend_position=:lt) 
	fig = Figure(
		size = (4*72, 6*72),
	)

    all_stars = first(datasets)[2]
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

    plot_labeled_R_h_ellipse!(observed_properties, 3)
    plot_labeled_R_h_ellipse!(observed_properties, 6)

	
    axislegend(position=legend_position)


	grid_low = fig[2, 1] = GridLayout()
	# cmd
    df_best = datasets[(collect∘keys)(datasets)[1]]
    Gmin = minimum(df_best.G) - 0.2
    Gmax = 21
	ax =  Axis(grid_low[1,1], 
		yreversed=true,
		xlabel=plot_labels[:bp_rp],
		ylabel=plot_labels[:G],
		limits=(-0.5, 2.5, Gmin, Gmax),
		xticks = 0:2,
		xgridvisible=false,
		ygridvisible=false,
		#aspect = 1,
	)

	for (label, df) in datasets
        Ag, Ab, Ar = get_extinction(df.ra, df.dec, df.bp_rp)
		scatter!(df.bp_rp .- Ab .+ Ar, df.phot_g_mean_mag .- Ag; scatter_kwargs[label]...)
	end

    G_errorbar = Gmin + (Gmax - Gmin) / 8
	df = datasets[:members]
	errorscatter!([-0.125], [G_errorbar], xerror=[median(df.dBP .+ df.dRP)], yerror=[median(df.dG)], color=:black, markersize=0, linewidth=1)

    isochrone = get_isochrone(observed_properties["metallicity"], age)
    G_iso = isochrone.Gmag .+ observed_properties["distance_modulus"]
    BP_RP_iso = isochrone.G_BPftmag .- isochrone.G_RPmag
    lines!(BP_RP_iso, G_iso, color=COLORS[2], linewidth=1)


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
	errorscatter!([-7.5], [7.5], xerror=[median(df.pmra_error)], yerror=[median(df.pmdec_error)], color=:black, markersize=0, linewidth=1)

    obs = LilGuys.collapse_errors(observed_properties)
    xe = maximum(error_interval(obs["pmra"]))
    ye = maximum(error_interval(obs["pmra"]))
    errorscatter!([observed_properties["pmra"]], [observed_properties["pmdec"]],
                  xerror = [xe],
                  yerror = [ye],
        linewidth=1,
        color = COLORS[2],
        markersize=3
       )

	rowsize!(grid_low, 1, Aspect(1, 1))

	rowsize!(fig.layout, 1, Aspect(1, 1))

	#resize_to_layout!(fig)
	fig
end


function plot_labeled_R_h_ellipse!(observed_properties, n=3)
    ellipse!(n*observed_properties["R_h"], observed_properties["ellipticity"], observed_properties["position_angle"], color=COLORS[2], linewidth=1)
    b = n * observed_properties["R_h"] * sqrt(1 - observed_properties["ellipticity"])
    θ = observed_properties["position_angle"]
    text!(b*cosd(θ), -b*sind(θ), text=L"%$(n)R_h", rotation=deg2rad(θ-90), 
          align = (:center, :bottom), color=COLORS[2], fontsize=10)
end

function get_isochrone(M_H, age=12)
    iso_columns = string.(split("Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	mbolmag  Gmag    G_BPbrmag  G_BPftmag  G_RPmag", 
	r"\s+"))

    all_isochrones = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/padova/isochrone.gaiadr2.weiler2018.$(age)Gyrs.dat"),DataFrame,
					  comment="#", ignorerepeated=true, delim = ' ', header=iso_columns)

	
	M_Hs = unique(all_isochrones.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]
	@info "using metallicity $M_H_adopted"

	filt = isapprox.(all_isochrones.MH, M_H_adopted)
	filt .&= all_isochrones.label .< 4

	isochrone =  all_isochrones[filt, :]
    # add nans between stages
    #insert!(isochrone, findlast(isochrone.label .== 3) + 1, fill(NaN, size(isochrone, 2)), promote=true)
    isochrone
end
