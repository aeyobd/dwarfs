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
	
	import LinearAlgebra: diag
	
	import LilGuys as lguys
	using Arya
	
end

# ╔═╡ 4cae5cc6-f270-42bf-97a1-067b7f57a7da
include("filter_utils.jl")

# ╔═╡ 48caecb2-180c-4ce4-a57b-6fed82328b01
md"""
# Introduction 

Given a set of TOML parameters,
this notebook simply filters the J+24 or Gaia-like observations according to the specifications.

Use the `calc_density`... and the `fit_profile` notebooks to then analyze the sample
"""

# ╔═╡ 47b0d3e6-a79b-4f49-b847-e708e5d6aabf
md"""
 # setup
"""

# ╔═╡ acb9ae92-924b-4723-8bd7-d775595b24c3
COLORS = Arya.COLORS;

# ╔═╡ ff92927e-b078-45fd-9c13-1ce5a009d0bb
red = COLORS[6]

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "sculptor/fiducial"

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
begin 
	param_file = "$name.toml"

	params_raw = read_file(param_file)
end

# ╔═╡ f8779f92-3ae0-474e-907e-1067170b1531
params = DensityParams(params_raw)

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

# ╔═╡ 0c8fd415-01e6-4295-ab40-b8beb9c82e4c
params

# ╔═╡ d1cc0201-f5da-4747-ae20-b14a03f1abd6
members

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
let 

	fig = plot_all_tangent(all_stars, markersize=2,
        color=(:grey, 0.2))
	
	plot_all_tangent!(current_axis(), members, markersize=2, color=red)

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
	fig = plot_all_tangent(all_stars, markersize=5,
        color=(:grey, 0.2), r_max=30, scale=60, units="arcmin")
	plot_all_tangent!(current_axis(), members, scale=60, markersize=5, color=red)
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

# ╔═╡ 28ff0827-dd3e-43ff-b210-9a45687dd1f8
let
	f1 = FITS(params.filename)
	header = read_header(f1[2])

	cols = FITSIO.colnames(f1[2])
	close(f1)
	
	f = FITS(out_name, "w")

	df = Dict(String(name) => members[:, name] for name in names(members))
	write(f, df,  header=header)

	println("wrote to $out_name")
	close(f)
	df
end

# ╔═╡ 62388099-b80b-4272-9bde-0c7315b15c19
r = 60 * members.r_ell # arcminutes

# ╔═╡ a254083b-61f2-45ad-b678-c4df6f16964b
let
	fig = Figure()
	ax = PolarAxis(fig[1,1])
	ϕ = atan.(members.eta, members.xi)
	scatter!(ax, ϕ, log10.(members.r_ell), markersize=5, alpha=0.2)

	fig
end

# ╔═╡ a7209445-84e9-435d-9144-90f8eb5e70cb
md"""
# Membership selection effects
"""

# ╔═╡ 1610647f-7dcf-4ffc-8906-bf28d36ae529
log_r_label = "log r / arcmin"

# ╔═╡ 648439e3-a771-4c17-a8df-2374a5369109
rs = lguys.calc_r_ell(members.xi, members.eta, params.ecc, params.PA)

# ╔═╡ 2a377969-3658-4b05-87f7-7de33859a588
xi = members.xi; eta = members.eta

# ╔═╡ 50f5a1f7-7537-4166-a582-e9c3453de0e0
x_p, y_p = lguys.shear_points_to_ellipse(xi, eta, params.ecc, params.PA)

# ╔═╡ 0723118a-5685-4b81-9098-4e702c9f4d8b
poly = lguys.convex_hull(x_p, y_p)

# ╔═╡ 13fb3ebc-50c0-43aa-88e9-1a7543e4e202
let
	fig, ax = Arya.FigAxis(xlabel="G mag", ylabel="count")
	
	hist!(members.phot_g_mean_mag)

	fig
end

# ╔═╡ 45422d53-317c-4824-a41a-4a80b1fbd102
let 
	fig = Figure(size=(900, 500))
	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		xlabel=log_r_label,
	)

	g_cuts = [21.5, 20.5, 20, 19.5, 19, 10]

	for i in 1:length(g_cuts) - 1
		# filter by G mag
		g_l = g_cuts[i+1]
		g_h = g_cuts[i]
		filt = g_cuts[i+1] .< members.phot_g_mean_mag .<= g_cuts[i]
		
		memb = members[filt, :]

		# calc properties
		r, r_max = lguys.calc_r_ell_sky(memb.ra, memb.dec, params.ecc, params.PA)
		println(r_max)
		
		obs = lguys.calc_properties(r[r .< r_max])

		# plot
		lines!(ax, obs.log_r, obs.log_Sigma, 
			color=i, colorrange=(1, length(g_cuts) - 1),
			label="$g_l, $g_h")

		# diagnostics
		println("$g_l, $g_h")
		println("number of members ", size(memb, 1))
		println("r_ell_max ", maximum(r))
		println(maximum(memb.phot_g_mean_mag), " > G > ",  minimum(memb.phot_g_mean_mag))
		println()
	end

	Legend(fig[1,2], ax, "G magnitude")

	fig
end

# ╔═╡ 80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
size(members)

# ╔═╡ edb7b9ff-4902-4004-b5b4-50a0f0b0ee52
sort(members.PSAT)

# ╔═╡ c0b3c3f6-0450-4242-9e13-41f9af17e562
let 
	fig = Figure(size=(900,500))
	p_cuts = [nothing, 0.0, 0.0001, 0.1,0.9, 0.99, 0.9999]

	Nc = length(p_cuts) 

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		#xlabel=log_r_label,
		limits=(nothing, (-4, 3))
		
	)

	labels = string.(p_cuts)

	ls = []
	for i in 1:Nc
		params = DensityParams(params_raw, PSAT_min=p_cuts[i])
		_, memb = load_and_filter(params)
		r, r_max = lguys.calc_r_ell_sky(memb.ra, memb.dec, params.ecc, params.PA)
		println(r_max)
		r = r[r .< r_max]
		
		obs = lguys.calc_properties(r, normalization=false)
		l = lines!(ax, obs.log_r, log10.(obs.Sigma), color=i, colorrange=(1, Nc+1))
		push!(ls, l)
	end

	Legend(fig[1,2], ls, labels, "probability cut")

	fig
end

# ╔═╡ 556d5c09-d531-4ce0-ab01-1810c6db4deb
pol = lguys.convex_hull(x_p, y_p)

# ╔═╡ 8c882887-43df-472f-91db-4d9d8e201639
r_max = lguys.min_distance_to_polygon(pol...)

# ╔═╡ f9606772-d51d-4701-80dc-5637b3510627
let
	fig, ax = FigAxis(aspect=DataAspect(),
		xlabel="a",
		ylabel="b"
	)


	scatter!(x_p, y_p, alpha=0.2, markersize=5)
	poly!(pol..., color=nothing, strokecolor=COLORS[2], strokewidth=2)

	poly!(Circle(Point2f(0,0,), r_max), color=nothing, strokecolor=COLORS[3], strokewidth=2)

	fig
end

# ╔═╡ a348b2db-b6c6-49fb-8207-a9e2f9cd7e87
# ╠═╡ disabled = true
#=╠═╡
poly = Polyhedra.planar_hull(Polyhedra.convexhull(ps...))
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═48caecb2-180c-4ce4-a57b-6fed82328b01
# ╠═47b0d3e6-a79b-4f49-b847-e708e5d6aabf
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═acb9ae92-924b-4723-8bd7-d775595b24c3
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═4cae5cc6-f270-42bf-97a1-067b7f57a7da
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╠═f8779f92-3ae0-474e-907e-1067170b1531
# ╟─4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═07235d51-10e1-4408-a4d1-cd2079fadb75
# ╠═695d532c-86d1-4b24-b7af-600a8ca29687
# ╠═32fd9b79-1a8b-4a69-9115-9065dd61ced2
# ╠═44a44f97-9115-4610-9706-33acf065d0e7
# ╠═753a7958-e12a-486e-b65b-7b6a8002a400
# ╟─0c498087-0184-4da2-a079-e972dd987712
# ╠═52bb6b36-736a-45a8-b1e1-7f174b366ec8
# ╠═f890216a-2e4e-4f44-92ee-ded0eaa17a68
# ╠═0c8fd415-01e6-4295-ab40-b8beb9c82e4c
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
# ╠═28ff0827-dd3e-43ff-b210-9a45687dd1f8
# ╠═62388099-b80b-4272-9bde-0c7315b15c19
# ╠═a254083b-61f2-45ad-b678-c4df6f16964b
# ╟─a7209445-84e9-435d-9144-90f8eb5e70cb
# ╠═1610647f-7dcf-4ffc-8906-bf28d36ae529
# ╠═648439e3-a771-4c17-a8df-2374a5369109
# ╠═2a377969-3658-4b05-87f7-7de33859a588
# ╠═50f5a1f7-7537-4166-a582-e9c3453de0e0
# ╠═0723118a-5685-4b81-9098-4e702c9f4d8b
# ╠═13fb3ebc-50c0-43aa-88e9-1a7543e4e202
# ╠═45422d53-317c-4824-a41a-4a80b1fbd102
# ╠═80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
# ╠═edb7b9ff-4902-4004-b5b4-50a0f0b0ee52
# ╠═c0b3c3f6-0450-4242-9e13-41f9af17e562
# ╠═556d5c09-d531-4ce0-ab01-1810c6db4deb
# ╠═8c882887-43df-472f-91db-4d9d8e201639
# ╠═f9606772-d51d-4701-80dc-5637b3510627
# ╠═a348b2db-b6c6-49fb-8207-a9e2f9cd7e87
