### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ d5bec398-03e3-11ef-0930-f3bd4f3c64fd
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using GLMakie
	using Measurements
	#using KernelDensity
	
	import SciPy
	using QuadGK
	
	import LinearAlgebra: diag
	
	import LilGuys as lguys
	using Arya
	
	using JSON
	import YAML
end

# ╔═╡ 542dd574-194a-49a2-adfc-a56db2ebea31
using Tables

# ╔═╡ 48caecb2-180c-4ce4-a57b-6fed82328b01
md"""
Given a set of yaml parameters,
this notebook simply filters the J23 observations according to the specifications.

Use the calc_density... and the fit_profile notebooks to then analyze the sample
"""

# ╔═╡ 47b0d3e6-a79b-4f49-b847-e708e5d6aabf
md"""
 # setup
"""

# ╔═╡ acb9ae92-924b-4723-8bd7-d775595b24c3
COLORS = Arya.COLORS;

# ╔═╡ 95a463d3-fdec-4fa6-9f4f-6d4a485aebf1
F = Float64

# ╔═╡ d6cbd9cb-bfa7-46c1-970a-ab3fb3740c48
OptF = Union{F, Nothing}

# ╔═╡ e3f52a32-9a21-4c8f-8e5c-9ca3e7dbc331
begin
	Base.@kwdef struct DensityParams
		filename::String
		ra::Float64
		dec::Float64
		ecc::F # can likely remove ecc, rh, PA
		rh::F
		PA::F
		dist::F
		dist_err::F
		PSAT_min::OptF = nothing
		ruwe_max::OptF = nothing
		g_min::OptF = nothing
		g_max::OptF = nothing
		max_ang_dist::OptF = nothing
		n_sigma_dist::OptF = nothing
		pmra::OptF = nothing
		pmdec::OptF = nothing
		dpm::OptF = nothing

		cmd_cut::Union{Array,Nothing} = nothing
	end

	function DensityParams(dict::Dict; kwargs...)
		d2 =  NamedTuple{Tuple(Symbol.(keys(dict)))}(values(dict))

		return DensityParams(; d2..., kwargs...)
	end
end

# ╔═╡ ff92927e-b078-45fd-9c13-1ce5a009d0bb
red = COLORS[6]

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "sculptor/fiducial"

# ╔═╡ 4093a7d6-2f74-4c37-a4a8-270934ede924
md"""
# functions
"""

# ╔═╡ fc4efd97-a140-4997-b668-904aa7ff5d30
function read_file(filename)
	f = YAML.load_file(filename)

	if "inherits" ∈ keys(f)
		f1 = read_file(dirname(filename) * "/" * f["inherits"])
		delete!(f, "inherits")
		merge!(f, f1)
	end

	return f
end

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
begin 
	param_file = "$name.yml"

	params_json =read_file(param_file)
	params = DensityParams(params_json)
end

# ╔═╡ 695d532c-86d1-4b24-b7af-600a8ca29687
function plot_all_tangent!(ax, all_stars; scale=1, kwargs...)
    x = scale*all_stars.xi
    y = scale*all_stars.eta 
    return scatter!(ax, x, y; kwargs...)
end

# ╔═╡ 07235d51-10e1-4408-a4d1-cd2079fadb75
function plot_all_tangent(all_stars; scale=1, units="degrees", r_max=nothing, kwargs...)
    fig = Figure()

    x = scale*all_stars.xi
    y = scale*all_stars.eta
    
    if r_max === nothing
        r_max = max(maximum(abs.(x)), maximum(abs.(y)))
    end
    
    ax = Axis(fig[1,1], 
        xlabel=L"\xi / \textrm{%$units}", ylabel=L"\eta / \textrm{%$units}",
        aspect=1,
        limits=(-r_max, r_max, -r_max, r_max)
    )

    p = plot_all_tangent!(ax, all_stars; scale=scale, kwargs...) 

    return Makie.FigureAxisPlot(fig, ax, p)
end

# ╔═╡ 32fd9b79-1a8b-4a69-9115-9065dd61ced2
function load_fits(filename)
	f = FITS(filename)
	all_stars = DataFrame(f[2])
	close(f)
	return all_stars
end

# ╔═╡ 2dbdc00f-7f78-4cc7-b505-99115d8a1042
function add_xi_eta!(stars, ra0, dec0)
	xi, eta = lguys.to_tangent(stars.ra, stars.dec, ra0, dec0)
	
	stars[:, "xi"] = xi
	stars[:, "eta"] = eta
	stars
end

# ╔═╡ 2d8b2740-f9f2-4cca-823f-4b6491c31fe4
function apply_filter(df, func, params...)
	if any(params .== nothing)
		filt = trues(size(df, 1))
	else
		filt = func(df, params...)
	end
	println("filter cuts $func \t", sum(map(!, filt)))
	return filt
end

# ╔═╡ 2ce9bedf-1ecb-4773-af65-1feff0f42f76
function min_filter(x, attr, cut)
	return x[:, attr] .> cut
end

# ╔═╡ ba162d65-2ee8-4390-9450-3975c649a05b
max_filter(x, attr, cut) = x[:, attr] .< cut

# ╔═╡ 60444f7f-71ee-4886-b02b-66bdc5324f99
md"""
# Filtering
"""

# ╔═╡ a0683e1e-5210-40fd-8841-1a6a315d3efe
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

# ╔═╡ 029bb274-df79-4f8a-bdd2-13f293e38279
function cmd_filter(all_stars, cmd_cut)
	cmd_cut_m = reshape(cmd_cut, 2, :)
	filt_cmd = is_point_in_polygon.(zip(all_stars.bp_rp, all_stars.phot_g_mean_mag), [cmd_cut_m])
end

# ╔═╡ 4467fa31-7171-404e-a225-e3c27afa0f7d
function ang_dist_filter(all_stars, ra0, dec0, max_ang_dist)
filt_ang_dist = @. (
   	max_ang_dist ^2
    > (all_stars.ra - ra0)^2 * cosd(dec0)^2 
    + (all_stars.dec - dec0)^2
    )
end

# ╔═╡ d36ca6c6-f9bc-47c4-b728-e6badaf866b3
function pm_filter(all_stars, pmra, pmdec, dpm)
	filt_pm = @. (
	    dpm^2 
	    > (pmra - all_stars.pmra)^2 
	    + (pmdec - all_stars.pmdec)^2 
	    )
end

# ╔═╡ d86bb8bb-8d7a-4a22-964c-f3c8d90cc9f0
function psat_filter(all_stars, psat_min)
    println("number nan PSAT         \t", sum(isnan.(all_stars.PSAT)))
    println("number exactly zero PSAT\t", sum(all_stars.PSAT .== 0))
    println("number > zero           \t", sum(all_stars.PSAT .> 0))
    println("number == 1             \t", sum(all_stars.PSAT .== 1))

    println("total                   \t", length(all_stars.PSAT))
	return all_stars.PSAT .> psat_min
end

# ╔═╡ 9ee9ad05-41d2-4e26-8252-1ae322947bb1
function parallax_filter(all_stars, dist, dist_err, n_sigma_dist)
	parallax = 1/dist
	parallax_err = 1/dist * dist_err / dist_err

	sigma = @. sqrt(all_stars.parallax_error^2 + parallax_err^2)
	
	filt_parallax = @. (
    abs(all_stars.parallax - parallax) <  sigma * n_sigma_dist
    )
end

# ╔═╡ c933ca4e-993c-488d-9f1a-c735d411c559
function radius_filter(all_stars, ecc, r_max=maximum(all_stars.r_ell))
	r_max = r_max * (1 - ecc)
	return all_stars.r_ell .< r_max
end

# ╔═╡ 914bdd71-da5c-4bf8-894d-12f64df5ca02
function select_members(all_stars, params)
	filt = apply_filter(all_stars, psat_filter, params.PSAT_min)

	filt .&= apply_filter(all_stars, min_filter, :phot_g_mean_mag, params.g_min)
	filt .&= apply_filter(all_stars, max_filter, :phot_g_mean_mag, params.g_max)
	filt .&= apply_filter(all_stars, max_filter, :ruwe, params.ruwe_max)
	filt .&= apply_filter(all_stars, cmd_filter, params.cmd_cut)
	filt .&= apply_filter(all_stars, ang_dist_filter, params.ra, params.dec, params.max_ang_dist)
	filt .&= apply_filter(all_stars, parallax_filter, params.dist, params.dist_err, params.n_sigma_dist)

	filt .&= apply_filter(all_stars, pm_filter, params.pmra, params.pmdec, params.dpm)
	filt .&= apply_filter(all_stars, radius_filter, params.ecc)


	println(sum(filt), " stars remaining")
	return all_stars[filt, :]
	
end

# ╔═╡ 753a7958-e12a-486e-b65b-7b6a8002a400
function load_and_filter(params)
	all_stars = load_fits(params.filename)

	add_xi_eta!(all_stars, params.ra, params.dec)
	members = select_members(all_stars, params)
	return all_stars, members
end

# ╔═╡ 44a44f97-9115-4610-9706-33acf065d0e7
all_stars, members = load_and_filter(params)

# ╔═╡ 0c498087-0184-4da2-a079-e972dd987712
md"""
The next three plots compare how different the r_ell and xi and eta calculated here are from what is (presumably) in the given catalogue.
"""

# ╔═╡ 52bb6b36-736a-45a8-b1e1-7f174b366ec8
let
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel="probability", 
		ylabel="count",
		yscale=log10)

	hist!(ax, all_stars.PSAT[all_stars.PSAT .>= 0], 
		bins=20, label="j24")
	#stephist!(b22.Pmemb, label="b22")
	f
end

# ╔═╡ f890216a-2e4e-4f44-92ee-ded0eaa17a68
params.rh

# ╔═╡ d1cc0201-f5da-4747-ae20-b14a03f1abd6
members

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
let 
	fig, ax, p = plot_all_tangent(all_stars, markersize=2,
        color=(:grey, 0.2))
	plot_all_tangent!(ax, members, markersize=2, color=red)
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	fig
end

# ╔═╡ 0010dffc-9717-4747-b7c2-2e396097399b
let 
	fig, ax, p = scatter(members.ra, members.dec, markersize=3,
        color=(:grey, 0.2))
	dra = 1
	ax.limits = (params.ra .+ dra * [-1, 1] ./ cosd(params.dec), params.dec .+ dra*[-1, 1])
	ax.aspect =  1 #cosd(params.dec)
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	ax.xlabel="ra / degrees"
	ax.ylabel = "dec / degrees"
	fig
end

# ╔═╡ d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
let 
	fig, ax, p = plot_all_tangent(all_stars, markersize=5,
        color=(:grey, 0.2), r_max=30, scale=60, units="arcmin")
	plot_all_tangent!(ax, members, scale=60, markersize=5, color=red)
	fig
end


# ╔═╡ b6424e6f-9b0d-4f29-b53d-0bd814a67139
let	
	fig = Figure()
	da = 60
	ax = Axis(fig[1, 1], 
	    xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}",
	    aspect=1,
	    limits = (-da, da, -da, da))
	
	scatter!(ax, 
	    60*members.xi, 60*members.eta, 
	    color=members.phot_g_mean_mag, 
	    colormap=:greys,
	    markersize=5
	)
	
	
	fig
end

# ╔═╡ bffe11bd-4233-4a7c-9411-0dfb1ac79077
let
	fig = Figure()
	da = 15
	ax = Axis(fig[1, 1], limits=(-da, da, -da, da), aspect=1,
	    xlabel=L"\mu_{\alpha *} / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}")
	scatter!(ax, all_stars.pmra, all_stars.pmdec, 
	    color=(:black, 0.2), markersize=3)
	
	scatter!(ax, members.pmra, members.pmdec, 
	    color=(red, 1), markersize=3)
	fig
end

# ╔═╡ 0f002b56-8b8f-4025-8d7b-fb51423e8da0
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=1,
	    limits=(-0.5, 3, 10, 22), yreversed=true,
	    xlabel="bp - rp", 
	    ylabel="G",)
	scatter!(ax, all_stars.bp_rp, all_stars.phot_g_mean_mag, 
	    color=(:black, 0.2), markersize=1)
	
	scatter!(ax, members.bp_rp, members.phot_g_mean_mag, 
	    color=(red, 1), markersize=1)
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	
	fig
end

# ╔═╡ 049ff11e-c04c-41d9-abf1-ec040b799649
let
	fig = Figure()
	da = 15
	ax = Axis(fig[1, 1], aspect=1, limits=(-10, 10, 0, 4),
	    xlabel=L"\varpi / \mathrm{ mas }", 
	    ylabel=L"\delta\varpi / \mathrm{ mas }")
	scatter!(ax, all_stars.parallax, all_stars.parallax_error, 
	    color=(:grey, 0.2), markersize=1)
	
	scatter!(ax, members.parallax, members.parallax_error, 
	    color=(red, 1), markersize=1)
	fig
end

# ╔═╡ 4873af32-c387-4d42-909d-d39a25f56e24
md"""
# Save sample
"""

# ╔═╡ 5a172f47-589f-4b2c-8180-108b293cebf7
out_name = "$(name)_sample.fits"

# ╔═╡ 00a51936-8a24-4d27-a8fd-ff1d11c0aa91
columnindex

# ╔═╡ 49bd2648-9ed8-49d0-becc-ea58a83d9d22
FITSHeader

# ╔═╡ fafb4f70-7f7c-40de-a34f-69c9a8fe0920
FITSIO.Libcfitsio.CFITSIO

# ╔═╡ 28ff0827-dd3e-43ff-b210-9a45687dd1f8
let
	f1 = FITS(params.filename)
	header = read_header(f1[2])

	cols = FITSIO.colnames(f1[2])
	close(f1)
	
	f = FITS(out_name, "w")

	df = Dict(String(name) => members[:, name] for name in Tables.columnnames(members))
	write(f, df,  header=header)

	println("wrote to $out_name")
	close(f)
	df
end

# ╔═╡ 39f615b4-b62c-493e-8d11-be45910d79a8
md"""
# density calculation
"""

# ╔═╡ 9b3288b4-3c17-4325-9e2c-94f96328f3c3
function plot_rh!()
    vline!([log10.(r_h)], color="grey", s=:dash, z_order=1, label=L"r_h")
end

# ╔═╡ 62388099-b80b-4272-9bde-0c7315b15c19
r = 60 * members.r_ell # arcminutes

# ╔═╡ a254083b-61f2-45ad-b678-c4df6f16964b
let
	fig = Figure()
	ax = PolarAxis(fig[1,1])
	ϕ = atan.(members.eta, members.xi)
	scatter!(ax, ϕ, members.r_ell, markersize=5, alpha=0.1)

	fig
end

# ╔═╡ 5c117b2f-a32c-4afd-9c66-943ab4634e71
dist = params.dist ± params.dist_err

# ╔═╡ b20058a8-a9b6-49ff-b8ff-a2d45c76f645
R_s_over_R_h = 1.6783

# ╔═╡ 39cb37d9-f1d4-419e-9e19-c033bfba8556
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=(-0.8, 1.7, -8, 5))
	scatter!(value.(obs.log_r), value.(obs.Γ),)
	errorbars!(value.(obs.log_r), value.(obs.Γ), err.(obs.Γ) )
	
	lines!(pred.log_r, pred.Γ, label="exponential", color=COLORS[2])
	
	ax.xlabel = log_r_label
	ax.ylabel = L"\Gamma = d\,\log \Sigma / d\,\log r"
	
	fig
end

# ╔═╡ 82d90218-f32e-4b72-a99a-bc2a264d7dce
theme(:colorcycle)

# ╔═╡ a7209445-84e9-435d-9144-90f8eb5e70cb
md"""
# Membership selection effects
"""

# ╔═╡ 45422d53-317c-4824-a41a-4a80b1fbd102
let 
	fig = Figure(size=(900, 500))
	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		xlabel=log_r_label,
	)

	g_cuts = [21.5, 20.5, 20, 19.5, 19, 10]
	labels = []

	ls = []
	for i in 1:length(g_cuts) - 1
		g_l = g_cuts[i+1]
		g_h = g_cuts[i]
		filt = g_cuts[i+1] .< members.phot_g_mean_mag .<= g_cuts[i]
		push!(labels, "$g_l, $g_h")
		memb = members[filt, :]
		println(size(memb, 1))
		println(maximum(memb.r_ell))
		r = memb.r_ell * 60
		obs = calc_properties(r)

		y = log10.(obs.Σ)
		f2 = isfinite.(y)
		y = y[f2]
		x = obs.log_r[f2]
		l = lines!(ax, x, value.(y), color=i, colorrange=(1, length(g_cuts) - 1))
		push!(ls, l)
	end

	Legend(fig[1,2], ls, labels, "G magnitude")

	fig
end

# ╔═╡ 13fb3ebc-50c0-43aa-88e9-1a7543e4e202
hist(members.phot_g_mean_mag)

# ╔═╡ 80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
size(members)

# ╔═╡ c0b3c3f6-0450-4242-9e13-41f9af17e562
let 
	fig = Figure(size=(900,500))
	p_cuts = [nothing, 0.01, 0.05, 0.1, 0.2, 0.5, 0.99]

	Nc = length(p_cuts) 

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label,
		limits=(nothing, (-4, 3))
		
	)

	labels = string.(p_cuts)

	ls = []
	for i in 1:Nc
		params = DensityParams(params_json, PSAT_min=p_cuts[i])
		_, memb = load_and_filter(params)
		r = memb.r_ell * 60
		println(maximum(memb.r_ell))
		obs = calc_properties(r)
		l = lines!(ax, obs.log_r, log10.(value.(obs.Σ * obs.N)), color=i, colorrange=(1, Nc))
		push!(ls, l)
	end

	Legend(fig[1,2], ls, labels, "probability cut")

	fig
end

# ╔═╡ 49ae0572-5d6b-4935-bc95-0a845bb3df2f
md"""
# Background density
"""

# ╔═╡ d7984df8-84b1-41ff-b19b-dd17b1772d4a
r_max = maximum(sqrt.(all_stars.xi .^ 2 + all_stars.eta .^ 2))

# ╔═╡ e033e344-737e-46e8-ab85-5fe33d191f41
"""
A simple density calculation 
"""
function calc_offset_density(dra, ddec, r_cut; n_sigma_dist=3, dpm=1, cmd_cut=cmd_cut_umi)

	
	params = DensityParams(params_json, 
		ra=params_json["ra"] + dra, dec=params_json["dec"] + ddec, PSAT_min=nothing, max_ang_dist=r_cut, ecc=0,
		dpm=dpm,
		cmd_cut=cmd_cut, n_sigma_dist=n_sigma_dist
	)

	_, memb = load_and_filter(params)
	r = memb.r_ell * 60
	
	obs = calc_properties(r)
	return obs
end

# ╔═╡ b6eaa6be-4a23-4357-9ce8-40aa9f16d7f6
let 
	fig = Figure(size=(900,500))

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label
	)

	r_shift = 0.75 * r_max
	r_cut = r_max - r_shift

	
	dras = r_shift*[1, 0, -1, 0, 0]
	ddecs = r_shift*[0, 1, 0, -1, 0]

	Nc = length(dras) 

	Σs = Measurement[]
	
	for i in 1:Nc
		dra = dras[i]
		ddec =  ddecs[i]
		label = "$(round(dra, digits=1)), $(round(ddec, digits=1))"
		obs = calc_offset_density(dra, ddec, r_cut)
		scatter_dens!(obs, label=label)
		if abs(dra^2 + ddec^2) > 0 
			append!(Σs, obs.Σ * obs.N)
		end
	end

	Σ_m = sum(Σs) / length(Σs)
	println("log Sigma background = $(log10.(Σ_m))")
	global log_Σ_bg
	log_Σ_bg = log10(Σ_m)
	hlines!(value.(log_Σ_bg))
	hspan!(value.(log_Σ_bg) .- err.(log_Σ_bg), value.(log_Σ_bg) .+ err.(log_Σ_bg), alpha=0.1)
	
	Legend(fig[1,2], ax, "position offset \n(degrees)", merge=true)

	fig
end

# ╔═╡ f832459e-edcb-48b4-ba3c-1d75a23f51e0
let 
	fig = Figure(size=(900,500))

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label
	)

	r_shift = 0.75 * r_max
	r_cut = r_max - r_shift

	
	dras = r_shift*[1, 0, -1, 0, 0]
	ddecs = r_shift*[0, 1, 0, -1, 0]

	Nc = length(dras) 

	Σs = Measurement[]
	
	for i in 1:Nc
		dra = dras[i]
		ddec =  ddecs[i]
		label = "$(round(dra, digits=1)), $(round(ddec, digits=1))"
		obs = calc_offset_density(dra, ddec, r_cut, dpm=nothing, cmd_cut=nothing, n_sigma_dist=nothing)
		scatter_dens!(obs, label=label)
		if abs(dra^2 + ddec^2) > 0 
			append!(Σs, obs.Σ * obs.N)
		end
	end

	Σ_m = sum(Σs) / length(Σs)
	println("log Sigma background = $(log10.(Σ_m))")
	global log_Σ_bg2
	log_Σ_bg2 = log10(Σ_m)
	hlines!(value.(log_Σ_bg2))
	hspan!(value.(log_Σ_bg2) .- err.(log_Σ_bg2), value.(log_Σ_bg2) .+ err.(log_Σ_bg2), alpha=0.1)
	
	Legend(fig[1,2], ax, "position offset \n(degrees)", merge=true)

	fig
end

# ╔═╡ 48e41a6f-775d-4d9f-850d-df9bd20dcf09
let 
	fig = Figure(size=(900,500))


	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label,
		
	)
	
	params = DensityParams(params_json)
	_, memb = load_and_filter(params)
	r = memb.r_ell * 60
	obs = calc_properties(r)
	scatter_dens!(obs, label="fiducial")

	obs = calc_offset_density(0, 0, 2)
	scatter_dens!(obs, label="simple")

	
	params = DensityParams(params_json, PSAT_min=nothing, ecc=0, )
	_, memb = load_and_filter(params)
	r = memb.r_ell * 60
	obs = calc_properties(r)
	scatter_dens!(obs, label="all")


	hlines!(value.(log_Σ_bg), label="background")
	hspan!(value.(log_Σ_bg) .- err.(log_Σ_bg), value.(log_Σ_bg) .+ err.(log_Σ_bg), alpha=0.1)

	hlines!(value.(log_Σ_bg2),  label="all backbround")
	hspan!(value.(log_Σ_bg2) .- err.(log_Σ_bg2), value.(log_Σ_bg2) .+ err.(log_Σ_bg2), alpha=0.1)

	vlines!(log10.(60r_max * (1-params.ecc^2)))
	Legend(fig[1, 2], ax, merge=true)
	fig
end

# ╔═╡ 44803049-4cc5-4a23-990a-322934ccb076
params_json

# ╔═╡ 45ce7a5d-75d4-4c7f-8233-5b2f7dde3a95
params

# ╔═╡ Cell order:
# ╟─48caecb2-180c-4ce4-a57b-6fed82328b01
# ╟─47b0d3e6-a79b-4f49-b847-e708e5d6aabf
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═acb9ae92-924b-4723-8bd7-d775595b24c3
# ╠═95a463d3-fdec-4fa6-9f4f-6d4a485aebf1
# ╠═d6cbd9cb-bfa7-46c1-970a-ab3fb3740c48
# ╠═e3f52a32-9a21-4c8f-8e5c-9ca3e7dbc331
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╟─4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═fc4efd97-a140-4997-b668-904aa7ff5d30
# ╠═07235d51-10e1-4408-a4d1-cd2079fadb75
# ╠═695d532c-86d1-4b24-b7af-600a8ca29687
# ╠═32fd9b79-1a8b-4a69-9115-9065dd61ced2
# ╠═2dbdc00f-7f78-4cc7-b505-99115d8a1042
# ╠═44a44f97-9115-4610-9706-33acf065d0e7
# ╠═2d8b2740-f9f2-4cca-823f-4b6491c31fe4
# ╠═2ce9bedf-1ecb-4773-af65-1feff0f42f76
# ╠═ba162d65-2ee8-4390-9450-3975c649a05b
# ╠═753a7958-e12a-486e-b65b-7b6a8002a400
# ╟─60444f7f-71ee-4886-b02b-66bdc5324f99
# ╠═914bdd71-da5c-4bf8-894d-12f64df5ca02
# ╠═029bb274-df79-4f8a-bdd2-13f293e38279
# ╠═a0683e1e-5210-40fd-8841-1a6a315d3efe
# ╠═4467fa31-7171-404e-a225-e3c27afa0f7d
# ╠═d36ca6c6-f9bc-47c4-b728-e6badaf866b3
# ╠═d86bb8bb-8d7a-4a22-964c-f3c8d90cc9f0
# ╠═9ee9ad05-41d2-4e26-8252-1ae322947bb1
# ╠═c933ca4e-993c-488d-9f1a-c735d411c559
# ╟─0c498087-0184-4da2-a079-e972dd987712
# ╠═52bb6b36-736a-45a8-b1e1-7f174b366ec8
# ╠═f890216a-2e4e-4f44-92ee-ded0eaa17a68
# ╠═d1cc0201-f5da-4747-ae20-b14a03f1abd6
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
# ╠═0010dffc-9717-4747-b7c2-2e396097399b
# ╠═d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
# ╠═b6424e6f-9b0d-4f29-b53d-0bd814a67139
# ╠═bffe11bd-4233-4a7c-9411-0dfb1ac79077
# ╠═0f002b56-8b8f-4025-8d7b-fb51423e8da0
# ╠═049ff11e-c04c-41d9-abf1-ec040b799649
# ╟─4873af32-c387-4d42-909d-d39a25f56e24
# ╠═5a172f47-589f-4b2c-8180-108b293cebf7
# ╠═00a51936-8a24-4d27-a8fd-ff1d11c0aa91
# ╠═49bd2648-9ed8-49d0-becc-ea58a83d9d22
# ╠═fafb4f70-7f7c-40de-a34f-69c9a8fe0920
# ╠═28ff0827-dd3e-43ff-b210-9a45687dd1f8
# ╠═542dd574-194a-49a2-adfc-a56db2ebea31
# ╟─39f615b4-b62c-493e-8d11-be45910d79a8
# ╠═9b3288b4-3c17-4325-9e2c-94f96328f3c3
# ╠═62388099-b80b-4272-9bde-0c7315b15c19
# ╠═a254083b-61f2-45ad-b678-c4df6f16964b
# ╠═5c117b2f-a32c-4afd-9c66-943ab4634e71
# ╠═b20058a8-a9b6-49ff-b8ff-a2d45c76f645
# ╠═39cb37d9-f1d4-419e-9e19-c033bfba8556
# ╠═82d90218-f32e-4b72-a99a-bc2a264d7dce
# ╟─a7209445-84e9-435d-9144-90f8eb5e70cb
# ╠═45422d53-317c-4824-a41a-4a80b1fbd102
# ╠═13fb3ebc-50c0-43aa-88e9-1a7543e4e202
# ╠═80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
# ╠═c0b3c3f6-0450-4242-9e13-41f9af17e562
# ╟─49ae0572-5d6b-4935-bc95-0a845bb3df2f
# ╠═d7984df8-84b1-41ff-b19b-dd17b1772d4a
# ╠═e033e344-737e-46e8-ab85-5fe33d191f41
# ╠═48e41a6f-775d-4d9f-850d-df9bd20dcf09
# ╠═b6eaa6be-4a23-4357-9ce8-40aa9f16d7f6
# ╠═f832459e-edcb-48b4-ba3c-1d75a23f51e0
# ╠═44803049-4cc5-4a23-990a-322934ccb076
# ╠═45ce7a5d-75d4-4c7f-8233-5b2f7dde3a95
