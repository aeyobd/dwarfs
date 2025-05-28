import StatsBase: median
import Random
using Makie
using OrderedCollections
using CSV, DataFrames
using LilGuys
using Arya


SELECTION_COLORS = [colorant"#808080"; to_colormap(:YlGnBu_5)[end-2:end]]
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



function plot_tangent!(gs, datasets, scatter_kwargs, observed_properties)
    all_stars = first(datasets)[2]
	# tangent
    dθ = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))
	ax = Axis(gs,
		xlabel=plot_labels[:xi_am], 
		ylabel=plot_labels[:eta_am], 
		#aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false,
		xreversed = true
	)


	for (label, df) in datasets
		scatter!(df.xi, df.eta; scatter_kwargs[label]...)
        if label == :members
            ellipse!(6observed_properties["R_h"], observed_properties["ellipticity"], observed_properties["position_angle"], color=:black)
        end
	end

    b = 6observed_properties["R_h"] * sqrt(1 - observed_properties["ellipticity"])
    θ = observed_properties["position_angle"]
    text!(b*cosd(θ), -b*sind(θ), text=L"6R_h", rotation=deg2rad(θ-90), 
          align = (:center, :bottom), color=:black, fontsize=0.8*theme(:fontsize)[])


    return ax
end


function compare_j24_samples(datasets, scatter_kwargs, observed_properties; 
        age=12)
	fig = Figure(
		size = (1720, 800),
	)

    ax = plot_tangent!(fig[1,1], datasets, scatter_kwargs, observed_properties)

	# cmd
    df_best = datasets[(collect∘keys)(datasets)[1]]
    Gmin = minimum(df_best.G) - 0.2
    Gmax = 21

    G_errorbar = Gmin + (Gmax - Gmin) / 8
	df = datasets[:members]
    #errorscatter!([-0.125], [G_errorbar], xerror=[median(df.dBP .+ df.dRP)], yerror=[median(df.dG)], color=:black, markersize=0)

    isochrone = get_isochrone(observed_properties["metallicity"], age)
    G_iso = isochrone.Gmag .+ observed_properties["distance_modulus"]
    BP_RP_iso = isochrone.G_BPftmag .- isochrone.G_RPmag


	ax_cmd =  Axis(fig[1,2], 
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

        if label == :members
            lines!(BP_RP_iso, G_iso, color=:black)
        end
	end



	# proper motions
	ax_pm = Axis(fig[1,3],
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
    #errorscatter!([-7.5], [7.5], xerror=[median(df.pmra_error)], yerror=[median(df.pmdec_error)], color=:black, markersize=0)

    obs = LilGuys.collapse_errors(observed_properties)
    xe = maximum(error_interval(obs["pmra"]))
    ye = maximum(error_interval(obs["pmra"]))
    scatter!([observed_properties["pmra"]], [observed_properties["pmdec"]],
             color=:black
       )


	rowsize!(fig.layout, 1, Aspect(1, 1))

    Legend(fig[2, 2], ax, tellwidth=false, tellheight=true, nbanks=3)
	#resize_to_layout!(fig)
	fig
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



function integrate_isodensity(pot, x0=[8., 0.]; x_direction=2, y_direction=3, s_factor=0.01, h_factor=0.01, h0=0.0001)
	x = x0[1]
	y = x0[2]

	xs = [x]
	ys = [y]
	h = h0
    s = h0

	θ = atan(y, x)
	x_vec = zeros(3)
	x_vec[x_direction] = 1

	y_vec = zeros(3)
	y_vec[y_direction] = 1
	ρ(x) = Agama.density(pot, x)
	ρ_0 = ρ(x_vec * x0[1] .+ y_vec * x0[2])
    dlρ_max = 0
	
    dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
    dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h

    s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
    h = h_scale * s

	for i in 1:100000
		x0 = zeros(3)
		x0[x_direction] = x
		x0[y_direction] = y
		dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
		dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h
		x += ds .* dy
		y += -ds .* dx

		push!(xs, x)
		push!(ys, y)

        s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
        h = h_scale * s

		θ_new = atan(y, x)
        dlρ_max = max(dlρ_max, abs(log10(ρ(x0) - ρ_0)))

		if θ_new > 0 && (θ < 0)
			break
		end

		θ = θ_new
	end
    @info "max log rel error = $dlρ_max"

    return xs, ys
end


X_SUN = [-LilGuys.GalactocentricFrame().d, 0, 0]
function plot_sun!(; x_direction=2, y_direction=3)
	scatter!(X_SUN[x_direction], X_SUN[y_direction], marker=:star5, color=COLORS[9])
end
		
