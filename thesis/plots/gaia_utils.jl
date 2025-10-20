
using PythonCall

Dustmaps = pyimport("dustmaps.sfd")
DUSTMAP = Dustmaps.SFDQuery()
SkyCoord = pyimport("astropy.coordinates").SkyCoord

include("utils.jl")


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


"""
    compare_j24_samples(
        datasets::AbstractDict, 
        scatter_kwargs::AbstractDict, 
        observed_properties::AbstractDict; 
        age::Real=12, 
        legend_position=:lt, 
        title::AbstractString=""
    )

Compare J+24 samples by plotting tangent-plane, CMD, and PM diagrams side by side.
"""
function compare_j24_samples(datasets::AbstractDict, scatter_kwargs::AbstractDict, observed_properties::AbstractDict; 
        age=12, legend_position=:lt, title="") 

	fig = Figure(
        size = (6 * 72, 5*72)
	)

    if title != ""
        Label(fig[0, 1:2], title, 
              fontsize=theme(:fontsize)[]*1.2, tellwidth=false, font=:bold)
    end

    gs_left = GridLayout(fig[1,1])
    gs_right = GridLayout(fig[1, 2])
    colsize!(fig.layout, 2, Relative(2/5))

	# tangent
    ax = plot_tangent(gs_left[1,1], datasets, scatter_kwargs, observed_properties)
    rowsize!(gs_left, 1, Aspect(1, 1.0))

    Legend(gs_left[2, 1], ax, tellwidth=false, tellheight=true, valign=:top)

    plot_cmd(gs_right[1,1], datasets, scatter_kwargs, observed_properties; age=age)
    plot_pm(gs_right[2,1], datasets, scatter_kwargs, observed_properties)

    rowsize!(gs_right, 1, Aspect(1, 1.0))
	rowsize!(gs_right, 2, Aspect(1, 1.0))

	resize_to_layout!(fig)
	fig
end



"""
    plot_tangent(
        gs::GridLayout, 
        datasets::AbstractDict, 
        scatter_kwargs::AbstractDict, 
        observed_properties::AbstractDict
    )

Plots the tangent-plane projection (ξ, η) for the provided datasets, 
including labeled half-light ellipses from observed properties.
"""
function plot_tangent(gs, datasets::AbstractDict, scatter_kwargs::AbstractDict, observed_properties::AbstractDict)
    all_stars = datasets[:best]

    dθ = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))

	ax = Axis(gs,
		xlabel=plot_labels[:xi_am], 
		ylabel=plot_labels[:eta_am], 
		aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false,
		xreversed = true
	)

	for (label, df) in datasets
		scatter!(df.xi, df.eta; scatter_kwargs[label]...)
	end

    plot_labeled_R_h_ellipse!(observed_properties, 3)
    plot_labeled_R_h_ellipse!(observed_properties, 6)

    return ax
end


"""
    plot_cmd(
        gs::GridLayout, 
        datasets::AbstractDict, 
        scatter_kwargs::AbstractDict, 
        observed_properties::AbstractDict; 
        age::Real
    )

Plots the color–magnitude diagram (CMD) with extinction correction and isochrone overlay.
"""
function plot_cmd(gs, datasets::AbstractDict, scatter_kwargs::AbstractDict, observed_properties::AbstractDict; age)
    df_best = datasets[:best]
    Gmin = minimum(df_best.G) - 0.2
    Gmax = 21

	ax =  Axis(gs,
		yreversed=true,
		xlabel=plot_labels[:bp_rp],
		ylabel=plot_labels[:G],
		limits=(-0.5, 2.5, Gmin, Gmax),
		xticks = 0:2,
	)

	for (label, df) in datasets
        Ag, Ab, Ar = get_extinction(df.ra, df.dec, df.bp_rp)
		scatter!(df.bp_rp .- Ab .+ Ar, df.phot_g_mean_mag .- Ag; scatter_kwargs[label]...)
	end


    # plot member errorbar
    G_errorbar = Gmin + (Gmax - Gmin) / 8
	df = datasets[:members]
	errorscatter!([-0.125], [G_errorbar], 
                  xerror=[median(df.dBP .+ df.dRP)], yerror=[median(df.dG)], 
                  color=:black, markersize=0, linewidth=1)

    # plot isochrone
    isochrone = get_isochrone(observed_properties["metallicity"], age)
    G_iso = isochrone.Gmag .+ observed_properties["distance_modulus"]
    BP_RP_iso = isochrone.G_BPftmag .- isochrone.G_RPmag
    lines!(BP_RP_iso, G_iso, color=COLORS[2], linewidth=1)

    return ax
end



"""
    plot_pm(
        gs::GridLayout, 
        datasets::AbstractDict, 
        scatter_kwargs::AbstractDict, 
        observed_properties::AbstractDict
    )

Plots proper-motion (μₐ*, μ_δ) distributions for the given datasets and overplots
the observed mean proper motion with uncertainties.
"""
function plot_pm(gs, datasets::AbstractDict, scatter_kwargs::AbstractDict, observed_properties::AbstractDict)
	ax = Axis(gs,
		xlabel = plot_labels[:pmra],
		ylabel = plot_labels[:pmdec],
		limits=(-10, 10, -10, 10),
		xgridvisible=false,
		ygridvisible=false,
		)
	
	for (label, df) in datasets
		scatter!(df.pmra, df.pmdec; scatter_kwargs[label]...)
	end

    # scatter median membership error
	df = datasets[:members]
	errorscatter!([-7.5], [7.5], xerror=[median(df.pmra_error)], yerror=[median(df.pmdec_error)], color=:black, markersize=0, linewidth=1)

    # plot observed proper motion
    obs = LilGuys.collapse_errors(observed_properties)
    xe = maximum(error_interval(obs["pmra"]))
    ye = maximum(error_interval(obs["pmdec"]))

    errorscatter!([observed_properties["pmra"]], [observed_properties["pmdec"]],
        xerror = [xe], yerror = [ye],
        linewidth = 1,
        color = COLORS[2],
        markersize = 3
       )

    return ax
end


"""
    plot_labeled_R_h_ellipse!(
        observed_properties::AbstractDict, 
        n::Integer=3
    )

Plots an ellipse at radius n×Rₕ with the observed ellipticity and position angle,
and labels it accordingly.
"""
function plot_labeled_R_h_ellipse!(observed_properties::AbstractDict, n::Real=3)
    R_h = observed_properties["R_h"]
    ell = observed_properties["ellipticity"]
    θ = observed_properties["position_angle"]

    ellipse!(n*R_h, ell, θ, color=COLORS[2], linewidth=1)

    b = n * R_h * sqrt(1 - ell)
    text!(b*cosd(θ), -b*sind(θ), text=L"%$(n)R_h", rotation=deg2rad(θ-90), 
          align = (:center, :bottom), color=COLORS[2], fontsize=10)
end



"""
    get_isochrone(M_H::Real, age::Real=12)

Loads and filters a Padova isochrone of a given metallicity and age.
Returns a DataFrame containing stellar evolution quantities and magnitudes.
"""

function get_isochrone(M_H::Real, age::Real=12)
    iso_columns = string.(split("Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	mbolmag  Gmag    G_BPbrmag  G_BPftmag  G_RPmag", 
	r"\s+"))

    all_isochrones = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/padova/isochrone.gaiadr2.weiler2018.$(age)Gyrs.dat"),DataFrame,
					  comment="#", ignorerepeated=true, delim = ' ', header=iso_columns)

	
	M_Hs = unique(all_isochrones.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]
	@info "using isochrone with metallicity $M_H_adopted, age $age"

	filt = isapprox.(all_isochrones.MH, M_H_adopted)
	filt .&= all_isochrones.label .< 4 # only keep through RGB

	all_isochrones[filt, :]
end



"""
    get_extinction(ra, dec, bp_rp)

Computes Gaia DR2 G, BP, and RP extinction values using the SFD dust map
and color-dependent extinction coefficients.
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



const default_styles = Dict(
	:best => (;	alpha=1, markersize=0.5, color=SELECTION_COLORS[1], 
		label="all" => (;markersize=2),
		rasterize=2,
	),
	:members_nospace => (;
		alpha=1, markersize=1,
		label = "CMD + PM" =>(alpha=1, markersize=2),
		color=SELECTION_COLORS[2],
		strokecolor=SELECTION_COLORS[2],
		strokewidth=0.3
	),
	:members => (;
		markersize=3,
		label = L"fiducial ($P_\textrm{sat} > 0.2$)" =>(alpha=1, markersize=2*2),
		marker=:rect,
		color = SELECTION_COLORS[3],
		alpha=1,
		strokecolor = SELECTION_COLORS[2],
		strokewidth=0.0
	),
	:rv => (;
		markersize=4,
		marker=:diamond,
		label = "RV members" =>(alpha=1, markersize=2.5*2),
		color = SELECTION_COLORS[4],
		strokecolor = :black,
		strokewidth=0.0
	),
	:rv_distant => (;
		markersize=8,
		marker=:star5,
		color = SELECTION_COLORS[4],
		strokewidth=1,
		strokecolor=COLORS[4]
	),
)
