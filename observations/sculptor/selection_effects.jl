### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ d5bec398-03e3-11ef-0930-f3bd4f3c64fd
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using CairoMakie

	using Polyhedra
	import LilGuys as lguys
	using Arya
end

# ╔═╡ 4cae5cc6-f270-42bf-97a1-067b7f57a7da
include("../../utils/gaia_filters.jl")

# ╔═╡ 48caecb2-180c-4ce4-a57b-6fed82328b01
md"""
# Introduction 

Are there confounding selection effects in the field (e.g. magnitude limits, PSAT, etc.)


"""

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ 3aa7e2b8-d51c-46d0-b961-f8e5d015aa57
cd("processed")

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "fiducial"

# ╔═╡ 9edb0794-5d2f-4070-9801-5b10a3a2cbc8
param_file = "$name.toml"

# ╔═╡ f8779f92-3ae0-474e-907e-1067170b1531
params = GaiaFilterParams(read_paramfile(param_file))

# ╔═╡ 4093a7d6-2f74-4c37-a4a8-270934ede924
md"""
# functions
"""

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
        limits=(-r_max, r_max, -r_max, r_max),
		xgridvisible=false,
		ygridvisible=false
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

# ╔═╡ 753a7958-e12a-486e-b65b-7b6a8002a400
function load_and_filter(params)
	all_stars = load_fits(params.filename)

	members = select_members(all_stars, params)
	return all_stars, members
end

# ╔═╡ 1c3d187a-cbc0-4e1d-8c2d-4a4e76d0402d
md"""
## loading inputs
"""

# ╔═╡ 0de39cdf-5723-4464-9498-a9908877aa42
all_stars_unfiltered = read_gaia_stars(params.filename, params)

# ╔═╡ 9d9d7778-68ba-405b-abb1-9c8a30125947
r_ell_max = 60calc_r_max(all_stars_unfiltered.ra, all_stars_unfiltered.dec, params.ellipticity, params.PA)

# ╔═╡ b824ab20-b42e-4f59-8c87-673796d10606
all_stars = all_stars_unfiltered[all_stars_unfiltered.r_ell .< r_ell_max, :]

# ╔═╡ 44a44f97-9115-4610-9706-33acf065d0e7
members = select_members(all_stars, params)

# ╔═╡ d1cc0201-f5da-4747-ae20-b14a03f1abd6
members

# ╔═╡ bf5ff417-8b1e-4209-86c1-df47483b1d26
md"""
# Validating r_ell / etc
"""

# ╔═╡ 933963cf-48e2-49d7-8c21-25c123bad456
df = all_stars_unfiltered;

# ╔═╡ cf741ca1-ad71-4ff6-925b-ad47c0cbde17
let
	fig, ax = FigAxis(
		limits=(nothing, (0.99, 1.01)),
		ylabel="relative error (jax / me)"
	)

	scatter!(df.xi, df.xi_original ./ df.xi)
	scatter!(df.eta, df.eta_original ./ df.eta)

	fig
end

# ╔═╡ cad5e46d-0088-4a36-bab7-9cc3d6922876
let
	fig, ax = FigAxis(
		ylabel="relative error (jax / me)"
	)

	# jax is in terms of a
	a = 12.33 * sqrt(1 - params.ellipticity)
	scatter!(df.r_ell, a * df.r_ell_original ./ df.r_ell)

	fig
end

# ╔═╡ e64257c1-15bc-4693-9f90-f3952d3559f5
let
	fig, ax = FigAxis(
		ylabel="relative error (jax / me transposed)"
	)

	r_ell_2 = 60lguys.calc_r_ell(df.xi, df.eta, params.ellipticity, -params.PA)

	# jax is in terms of a
	a = 12.33 .* sqrt(1 - params.ellipticity)
	scatter!(df.r_ell, a * df.r_ell_original ./ r_ell_2)

	fig
end

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

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

# ╔═╡ 52bb6b36-736a-45a8-b1e1-7f174b366ec8
let
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel="probability", 
		ylabel="count",
		yscale=log10)

	hist!(ax, all_stars.PSAT[all_stars.PSAT .>= 0], 
		bins=20, label="j24")

	f
end

# ╔═╡ 0aae7fec-f003-4ebc-be9b-8028408c121e
let	
	fig = Figure()
	da = 60
	ax = Axis(fig[1, 1], 
	    xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}",
	    aspect=1,
	    limits = (-da, da, -da, da))
	
	scatter!(ax, 
	    60*members.xi, 60*members.eta, 
	    color=members.PSAT, 
	    markersize=5
	)
	
	fig
end

# ╔═╡ a7209445-84e9-435d-9144-90f8eb5e70cb
md"""
# Membership selection effects
"""

# ╔═╡ 1610647f-7dcf-4ffc-8906-bf28d36ae529
log_r_label = "log r / arcmin"

# ╔═╡ 13fb3ebc-50c0-43aa-88e9-1a7543e4e202
let
	fig, ax = Arya.FigAxis(xlabel="G mag", ylabel="count")
	
	hist!(members.phot_g_mean_mag)

	fig
end

# ╔═╡ 45422d53-317c-4824-a41a-4a80b1fbd102
let 
	fig = Figure()
	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		xlabel=log_r_label,
		limits=((-1, 2), (-7, -2))
	)

	g_cuts = [21.5, 20.5, 20, 19.5, 19, 10]

	for i in 1:length(g_cuts) - 1
		# filter by G mag
		g_l = g_cuts[i+1]
		g_h = g_cuts[i]
		filt = g_cuts[i+1] .< members.phot_g_mean_mag .<= g_cuts[i]
		
		memb = members[filt, :]

		obs = lguys.StellarProfile(memb.r_ell)

		# plot
		lines!(ax, obs.log_r, obs.log_Sigma, 
			color=i, colorrange=(1, length(g_cuts)),
			label="$g_l, $g_h")

		# diagnostics
		println("$g_l, $g_h")
		println("number of members ", size(memb, 1))
		println("r_ell_max ", maximum(memb.r_ell))
		println(maximum(memb.phot_g_mean_mag), " > G > ",  minimum(memb.phot_g_mean_mag))
		println()
	end

	axislegend("G magnitude", position=:lb)

	fig
end

# ╔═╡ 858a2f19-393f-45cd-8aeb-2602c4cbae05
let 
	fig = Figure()
	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		xlabel=log_r_label,
		limits=((-1, 2.3), (-7, -2))
	)

	memb = all_stars[all_stars.PSAT_1C .> 0.2, :]

	obs = lguys.StellarProfile(memb.r_ell)
	lines!(ax, obs.log_r, obs.log_Sigma, 
		label="1 component")

	
	memb = all_stars[all_stars.PSAT_CIRC .> 0.2, :]
	r = @. sqrt(memb.xi^2 + memb.eta^2) * 60
	obs = lguys.StellarProfile(r)
	lines!(ax, obs.log_r, obs.log_Sigma, 
		label="2 component, circ")


	memb = all_stars[all_stars.PSAT_ELL .> 0.2, :]
	obs = lguys.StellarProfile(memb.r_ell)
	lines!(ax, obs.log_r, obs.log_Sigma, 
		label="2 component, ell")


	memb = all_stars[all_stars.PSAT_NOSPACE .> 0.2, :]
	obs = lguys.StellarProfile(memb.r_ell)
	lines!(ax, obs.log_r, obs.log_Sigma, 
		label="no spatial prior")
	
	axislegend(position=:lb)
	
	fig
end

# ╔═╡ 80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
size(members)

# ╔═╡ edb7b9ff-4902-4004-b5b4-50a0f0b0ee52
sort(members.PSAT)

# ╔═╡ c0b3c3f6-0450-4242-9e13-41f9af17e562
let 
	fig = Figure()
	p_cuts = [nothing, 0.0, 0.0001, 0.1,0.9, 0.99, 0.9999]

	Nc = length(p_cuts) 

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label,
		limits=(nothing, (-4, 3))
		
	)

	for i in 1:Nc
		p_cut =  p_cuts[i] 
		filt = apply_filter(all_stars, psat_filter,  p_cut)
		
		r = all_stars[filt, :r_ell]

		
		obs = lguys.StellarProfile(r, normalization=:none)
		
		lines!(ax, obs.log_r, log10.(obs.Sigma), color=i, colorrange=(1, Nc+1),
			label = string(p_cut)
		)
	end

	axislegend("probability cut", position=:lb)

	fig
end

# ╔═╡ 82881a33-8302-472f-80b7-c4040f8fc1e2
md"""
## Convex hull test
"""

# ╔═╡ 50f5a1f7-7537-4166-a582-e9c3453de0e0
x_p, y_p = 60 .* lguys.shear_points_to_ellipse(all_stars_unfiltered.xi, all_stars_unfiltered.eta, params.ellipticity, params.PA)

# ╔═╡ 556d5c09-d531-4ce0-ab01-1810c6db4deb
pol = convex_hull(x_p, y_p)

# ╔═╡ 8c882887-43df-472f-91db-4d9d8e201639
r_max = lguys.min_distance_to_polygon(pol...)

# ╔═╡ 0723118a-5685-4b81-9098-4e702c9f4d8b
poly = convex_hull(x_p, y_p)

# ╔═╡ f9606772-d51d-4701-80dc-5637b3510627
let
	fig, ax = FigAxis(aspect=DataAspect(),
		xlabel="a",
		ylabel="b"
	)


	scatter!(x_p, y_p, alpha=0.2, markersize=5)

	filt = all_stars_unfiltered.r_ell .< r_ell_max
	scatter!(x_p[filt], y_p[filt], alpha=0.2, markersize=5)

	poly!(pol..., color=:transparent, strokecolor=COLORS[2], strokewidth=2)

	poly!(Circle(Point2f(0,0,), r_max), color=:transparent, strokecolor=COLORS[3], strokewidth=2)

	fig
end

# ╔═╡ Cell order:
# ╠═48caecb2-180c-4ce4-a57b-6fed82328b01
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═4cae5cc6-f270-42bf-97a1-067b7f57a7da
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═3aa7e2b8-d51c-46d0-b961-f8e5d015aa57
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═9edb0794-5d2f-4070-9801-5b10a3a2cbc8
# ╠═f8779f92-3ae0-474e-907e-1067170b1531
# ╟─4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═07235d51-10e1-4408-a4d1-cd2079fadb75
# ╠═695d532c-86d1-4b24-b7af-600a8ca29687
# ╠═32fd9b79-1a8b-4a69-9115-9065dd61ced2
# ╠═753a7958-e12a-486e-b65b-7b6a8002a400
# ╟─1c3d187a-cbc0-4e1d-8c2d-4a4e76d0402d
# ╠═0de39cdf-5723-4464-9498-a9908877aa42
# ╠═9d9d7778-68ba-405b-abb1-9c8a30125947
# ╠═b824ab20-b42e-4f59-8c87-673796d10606
# ╠═44a44f97-9115-4610-9706-33acf065d0e7
# ╠═d1cc0201-f5da-4747-ae20-b14a03f1abd6
# ╠═bf5ff417-8b1e-4209-86c1-df47483b1d26
# ╠═933963cf-48e2-49d7-8c21-25c123bad456
# ╠═cf741ca1-ad71-4ff6-925b-ad47c0cbde17
# ╠═cad5e46d-0088-4a36-bab7-9cc3d6922876
# ╠═e64257c1-15bc-4693-9f90-f3952d3559f5
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═b6424e6f-9b0d-4f29-b53d-0bd814a67139
# ╠═52bb6b36-736a-45a8-b1e1-7f174b366ec8
# ╠═0aae7fec-f003-4ebc-be9b-8028408c121e
# ╟─a7209445-84e9-435d-9144-90f8eb5e70cb
# ╠═1610647f-7dcf-4ffc-8906-bf28d36ae529
# ╠═13fb3ebc-50c0-43aa-88e9-1a7543e4e202
# ╠═45422d53-317c-4824-a41a-4a80b1fbd102
# ╠═858a2f19-393f-45cd-8aeb-2602c4cbae05
# ╠═80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
# ╠═edb7b9ff-4902-4004-b5b4-50a0f0b0ee52
# ╠═c0b3c3f6-0450-4242-9e13-41f9af17e562
# ╟─82881a33-8302-472f-80b7-c4040f8fc1e2
# ╠═556d5c09-d531-4ce0-ab01-1810c6db4deb
# ╠═0723118a-5685-4b81-9098-4e702c9f4d8b
# ╠═8c882887-43df-472f-91db-4d9d8e201639
# ╠═50f5a1f7-7537-4166-a582-e9c3453de0e0
# ╠═f9606772-d51d-4701-80dc-5637b3510627
