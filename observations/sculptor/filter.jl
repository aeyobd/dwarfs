### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d5bec398-03e3-11ef-0930-f3bd4f3c64fd
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using CairoMakie

	
	import LilGuys as lguys
	using Arya
	
end

# ╔═╡ 26e401fa-aa9f-41ae-88fa-815deb6427c7
using KernelDensity

# ╔═╡ 4cae5cc6-f270-42bf-97a1-067b7f57a7da
include("../../utils/gaia_filters.jl")

# ╔═╡ 48caecb2-180c-4ce4-a57b-6fed82328b01
md"""
# Introduction 

Given a set of TOML parameters,
this notebook simply filters the J+24 or Gaia-like observations according to the specifications.

Use the `calc_density`... and the `fit_profile` notebooks to then analyze the sample
"""

# ╔═╡ acb9ae92-924b-4723-8bd7-d775595b24c3
COLORS = Arya.COLORS;

# ╔═╡ ff92927e-b078-45fd-9c13-1ce5a009d0bb
red = COLORS[6]

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ daac0dde-fe80-44e3-9e6c-d00384769710
md"""
Code snips
Gaia
```
SELECT * FROM gaiadr3.gaia_source WHERE 1 = CONTAINS( POINT(15.03917, -33.70917), CIRCLE(ra, dec, 3))
```

Aladin
```
15.03917, -33.70917	
zoom 20 arcmin
```
"""

# ╔═╡ bdef32be-2748-48c4-85ed-1d0d68c38b28
galaxy_dir = "processed"

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "fiducial"

# ╔═╡ f8779f92-3ae0-474e-907e-1067170b1531
params = GaiaFilterParams("$galaxy_dir/$name.toml")

# ╔═╡ 6f4ee4fd-0fe9-4314-a2ec-8b0caa08e8af
md"""
# Loading files
"""

# ╔═╡ 44a44f97-9115-4610-9706-33acf065d0e7
all_stars_unfiltered = load_stars("$galaxy_dir/" *params.filename, params)

# ╔═╡ ce3cead5-b49e-41ff-ae3f-bdcaab68e858
r_ell_max = 60lguys.calc_r_max(all_stars_unfiltered.ra, all_stars_unfiltered.dec, params.ellipticity, params.PA)

# ╔═╡ c3641354-58d7-42e1-97d8-98db1c0bf0ac
all_stars = all_stars_unfiltered

# ╔═╡ bfca40e9-7889-408b-a2ab-d13c22f6891b
sum(all_stars.PSAT .> 0.2)

# ╔═╡ 103b8a58-0d23-42df-a2d8-4ec508c0246a
members = select_members(all_stars, params)

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
		ygridvisible=false,
		xreversed=true,
    )

    p = plot_all_tangent!(ax, all_stars; scale=scale, kwargs...) 

    return Makie.FigureAxisPlot(fig, ax, p)
end

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ 7c7c360b-cc11-43a5-b6ee-f6347aa9ae32
bw = (0.05, 0.05)

# ╔═╡ c4b92086-db8c-421d-8ce2-f55863b1df18
kd = kde([members.xi members.eta]; bandwidth=bw)

# ╔═╡ 3a555897-1d4d-4391-bf9a-a92c317ca34f
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect(),
		xlabel = L"$\xi$ / degrees",
		ylabel = L"$\eta$ / degrees",
		xreversed = true,
		#limits=(-r_ell_max, r_ell_max, -r_ell_max, r_ell_max) ./ 60
	)
	scale = 0.05


	h = heatmap!(kd.x, kd.y, asinh.(kd.density ./ scale))

	Colorbar(fig[1, 2], h, label="asinh density / $scale")
	fig
end

# ╔═╡ 52bb6b36-736a-45a8-b1e1-7f174b366ec8
let
	if :PSAT in names(all_stars)
		f = Figure()
		ax = Axis(f[1, 1],
			xlabel="probability", 
			ylabel="count",
			yscale=log10)
	
		hist!(ax, all_stars.PSAT[all_stars.PSAT .>= 0], 
			bins=20, label="j24")
		vlines!(params.PSAT_min)
		#stephist!(b22.Pmemb, label="b22")
		f
	end
end

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
	dra = 6
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

# ╔═╡ 390ee4b1-9da8-4ac1-b9f4-905462115c38
x_poly = params.cmd_cut[1:2:end]

# ╔═╡ f1e9e618-ecc3-4d32-a23f-9f8e9ff388ba
y_poly = params.cmd_cut[2:2:end]

# ╔═╡ 6fad2fb0-6ba1-4359-83dd-7a342d1c0db8
[x_poly; x_poly[1]]

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

	lines!(ax, [x_poly; x_poly[1]], [y_poly; y_poly[1]])
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
out_name = "processed/$(name)_sample.fits"

# ╔═╡ 28ff0827-dd3e-43ff-b210-9a45687dd1f8
let
	f1 = FITS(joinpath(galaxy_dir, params.filename))
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

# ╔═╡ Cell order:
# ╟─48caecb2-180c-4ce4-a57b-6fed82328b01
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═acb9ae92-924b-4723-8bd7-d775595b24c3
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╠═4cae5cc6-f270-42bf-97a1-067b7f57a7da
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╟─daac0dde-fe80-44e3-9e6c-d00384769710
# ╠═bdef32be-2748-48c4-85ed-1d0d68c38b28
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═f8779f92-3ae0-474e-907e-1067170b1531
# ╟─6f4ee4fd-0fe9-4314-a2ec-8b0caa08e8af
# ╠═44a44f97-9115-4610-9706-33acf065d0e7
# ╠═ce3cead5-b49e-41ff-ae3f-bdcaab68e858
# ╠═c3641354-58d7-42e1-97d8-98db1c0bf0ac
# ╠═bfca40e9-7889-408b-a2ab-d13c22f6891b
# ╠═103b8a58-0d23-42df-a2d8-4ec508c0246a
# ╟─4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═07235d51-10e1-4408-a4d1-cd2079fadb75
# ╠═695d532c-86d1-4b24-b7af-600a8ca29687
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═7c7c360b-cc11-43a5-b6ee-f6347aa9ae32
# ╠═c4b92086-db8c-421d-8ce2-f55863b1df18
# ╠═3a555897-1d4d-4391-bf9a-a92c317ca34f
# ╠═52bb6b36-736a-45a8-b1e1-7f174b366ec8
# ╠═d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
# ╠═0010dffc-9717-4747-b7c2-2e396097399b
# ╠═d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
# ╠═b6424e6f-9b0d-4f29-b53d-0bd814a67139
# ╠═bffe11bd-4233-4a7c-9411-0dfb1ac79077
# ╠═390ee4b1-9da8-4ac1-b9f4-905462115c38
# ╠═f1e9e618-ecc3-4d32-a23f-9f8e9ff388ba
# ╠═6fad2fb0-6ba1-4359-83dd-7a342d1c0db8
# ╠═0f002b56-8b8f-4025-8d7b-fb51423e8da0
# ╠═049ff11e-c04c-41d9-abf1-ec040b799649
# ╠═26e401fa-aa9f-41ae-88fa-815deb6427c7
# ╟─4873af32-c387-4d42-909d-d39a25f56e24
# ╠═5a172f47-589f-4b2c-8180-108b293cebf7
# ╠═28ff0827-dd3e-43ff-b210-9a45687dd1f8
