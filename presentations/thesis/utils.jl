import StatsBase: median
import Random
using Makie
using OrderedCollections
using CSV, DataFrames
using LilGuys
using Arya
import Agama


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



function plot_tangent!(gs, datasets, scatter_kwargs, observed_properties; R_ell=[6])
    all_stars = first(datasets)[2]
	ax = TangentAxis(gs, all_stars)


	for (label, df) in datasets
		scatter!(df.xi, df.eta; scatter_kwargs[label]...)
        if label == :members
        end
	end

    for R in R_ell
        plot_R_h_ell!(observed_properties, R)
    end

    return ax
end


function plot_R_h_ell!(observed_properties, R=6)
    b = R* observed_properties["R_h"] * sqrt(1 - observed_properties["ellipticity"])
    θ = observed_properties["position_angle"]
    text!(b*cosd(θ), -b*sind(θ), text=L"%$(R)R_h", rotation=deg2rad(θ-90), 
          align = (:center, :top), color=:black, fontsize=0.8*theme(:fontsize)[])

    ellipse!(R * observed_properties["R_h"], observed_properties["ellipticity"], observed_properties["position_angle"], color=:black)
end


function compare_j24_samples(datasets, scatter_kwargs, observed_properties; 
        age=12, Gmin=nothing, references=true)
	fig = Figure(
		size = (1720, 800),
	)

    if references
        R_ell = [6]
    else
        R_ell  = []
    end
    ax = plot_tangent!(fig[1,1], datasets, scatter_kwargs, observed_properties, R_ell=R_ell)

    ax_cmd =  CMDAxis(fig[1,2], first(datasets)[2], Gmin=Gmin)

	for (label, df) in datasets
        scatter!(df.bp_rp_corrected, df.G_corrected; scatter_kwargs[label]...)
	end

    if references
        plot_isochrone!(observed_properties, age)
    end


	# proper motions
    ax_pm = PMAxis(fig[1,3])
	
	for (label, df) in datasets
		scatter!(df.pmra, df.pmdec; scatter_kwargs[label]...)
	end

    if references
        plot_obs_pm!(observed_properties)
    end


    # styling
	rowsize!(fig.layout, 1, Aspect(1, 1))
    Legend(fig[2, 2], ax, tellwidth=false, tellheight=true, nbanks=3)
	#resize_to_layout!(fig)
	fig
end


function plot_isochrone!(observed_properties, age)
	isochrone = Utils.get_isochrone(observed_properties["metallicity"], age)
    G_iso = isochrone.Gmag .+ observed_properties["distance_modulus"]
    BP_RP_iso = isochrone.G_BPftmag .- isochrone.G_RPmag
	lines!(BP_RP_iso, G_iso, color=:black)
end


function TangentAxis(gs, all_stars)
	dθ = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))

	ax = Axis(gs,
		xlabel=Utils.plot_labels[:xi_am], 
		ylabel=Utils.plot_labels[:eta_am], 
		#aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false,
		xreversed = true
	)
end

function CMDAxis(gs, all_stars; Gmin=nothing)
    if isnothing(Gmin)
        Gmin = minimum(all_stars.G) - 0.2
    end
    Gmax = 21
	ax_cmd =  Axis(gs, 
		yreversed=true,
		xlabel=Utils.plot_labels[:bp_rp],
		ylabel=Utils.plot_labels[:G],
		limits=(-0.5, 2.5, Gmin, Gmax),
		xticks = 0:2,
		xgridvisible=false,
		ygridvisible=false,
		#aspect = 1,
	)

end

function PMAxis(gs)
	ax_pm = Axis(gs,
		xlabel = plot_labels[:pmra],
		ylabel = plot_labels[:pmdec],
		#aspect=DataAspect(),
		limits=(-10, 10, -10, 10),
		xgridvisible=false,
		ygridvisible=false,
		)
end



function correct_extinction!(df)
    Ag, Ab, Ar = get_extinction(df.ra, df.dec, df.bp_rp)
    df[!, :bp_rp_corrected] = df.bp_rp .- Ab .+ Ar 
    df[!, :G_corrected] = df.phot_g_mean_mag .- Ag
    return df
end


function plot_obs_pm!(observed_properties)
    obs = LilGuys.collapse_errors(observed_properties)
    xe = maximum(error_interval(obs["pmra"]))
    ye = maximum(error_interval(obs["pmra"]))
    scatter!([observed_properties["pmra"]], [observed_properties["pmdec"]],
             color=:black
       )
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

X_SUN = [-LilGuys.GalactocentricFrame().d, 0, 0]


function integrate_isodensity(pot, initial=[-X_SUN[1], 0.]; x_direction=2, y_direction=3, s_scale=0.01, h_scale=0.01, h0=0.0001)
	x = initial[1]
	y = initial[2]

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
    ρ_0 = ρ(x_vec * x .+ y_vec * y)
    dlρ_max = 0

    x0 = zeros(3)
    x0[x_direction] = x
    x0[y_direction] = y
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
		x += s .* dy
		y += -s .* dx

		push!(xs, x)
		push!(ys, y)

        s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
        h = h_scale * s

		θ_new = atan(y, x)
        dlρ_max = max(dlρ_max, abs(log10(ρ(x0)) - log10(ρ_0)))

		if θ_new > 0 && (θ < 0)
			break
		end

		θ = θ_new
	end
    @info "max log rel error = $dlρ_max"

    return xs, ys
end


function plot_sun!(; x_direction=2, y_direction=3)
    scatter!(X_SUN[x_direction], X_SUN[y_direction], marker=:star5, color=COLORS[9], strokewidth=theme(:linewidth)[]/4, strokecolor=:black)
end
		
