### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames 
	using CSV
	using CairoMakie
	update_theme!(Axis=(xgridvisible=false, ygridvisible=false))

	import NaNMath as nm
	using KernelDensity
	using Measurements

end

# ╔═╡ 8b3955e5-ab08-4d5a-abe5-8f40120ff8db
using PlutoUI

# ╔═╡ 69c98029-165c-407b-9a63-a27e06e30e45
using Arya

# ╔═╡ ae29bed0-6700-47f1-8952-35e867ce126b
using OrderedCollections

# ╔═╡ 65822789-b261-45be-817f-cba0e37b82ce
include("../utils/gaia_filters.jl")

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sammple

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
md"""
galaxyname: $(@bind galaxyname confirm(TextField(default="galaxy")))
"""

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using LilGuys
	FIGDIR = joinpath(galaxyname, "figures")
	FIGSUFFIX = ".jax_observations"
end

# ╔═╡ da9ca7a7-18b8-49cb-a041-ab1c667920ff
import DensityEstimators: histogram2d

# ╔═╡ 489f6d21-9e9a-4b1e-b05c-c63a44ba1951
import StatsBase: percentile, mean

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ db1264b7-02c3-4b55-ae2b-9ce78aa1304a
import DensityEstimators: histogram

# ╔═╡ 81c677bf-8372-4fd3-9d60-af8377485f1c
import Statistics: cor

# ╔═╡ 53c0567a-c845-486d-9510-1590cc0ee9ff
CairoMakie.activate(type="svg", pt_per_unit=2)

# ╔═╡ cf7aeb0a-d453-462a-b43d-a832567440fd
diverging_cmap = Reverse(:bluesreds)

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ 9ecf79a8-2ed3-40c6-b555-a102250ecbd4
observed_properties = TOML.parsefile(galaxyname * "/observed_properties.toml")

# ╔═╡ 4c8b9c30-f873-4756-a729-aa023b3a804e
readdir("$galaxyname/data")

# ╔═╡ 26cf1867-02be-4d36-8c35-6c58a1feca27
if isfile("$galaxyname/data/jensen+24_2c.fits")
	datafile = "$galaxyname/data/jensen+24_2c.fits"
else
	datafile = "$galaxyname/data/jensen+24_1c.fits"
end

# ╔═╡ d87f28e5-cc94-4afc-9e70-ece1686db22c
filt_params = GaiaFilterParams(filename=datafile, observed_properties)

# ╔═╡ 338be0b8-e732-4202-a20f-9ca841516075
all_stars = read_gaia_stars(filt_params)

# ╔═╡ 7ed70894-48a8-46af-8c78-8a75cb39d957
md"""
estimate of centring uncertanty in this dataset (standard error on mean and present difference between mean centre & assumed centre)

"""

# ╔═╡ 3abea230-2b95-40d3-9851-a91236f75c4a
md"""
each component of the cen_errs tuple below should be similar in magnitude
"""

# ╔═╡ 082a06dd-eeb5-4761-a233-1ee89e8cb819
best_stars = all_stars[all_stars.F_BEST .== 1.0, :]

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = best_stars[best_stars.PSAT .> 0.2, :]

# ╔═╡ 223abc41-5158-49c2-96bf-df55b7be1114
cen = lguys.to_tangent(lguys.calc_centre2D(members.ra, members.dec, "mean")..., observed_properties["ra"], observed_properties["dec"])

# ╔═╡ 32e8e100-1fad-425f-b05c-2310e7ba559c
cen_errs =  (cen..., 
lguys.std(members.eta) / sqrt(length(members.eta)), 
lguys.std(members.xi) / sqrt(length(members.xi))
)

# ╔═╡ b1be5f97-0e61-4dda-87fd-7bb0f43147b6
cen_err = maximum(abs.(cen_errs))

# ╔═╡ 731ea468-5003-44e9-95b8-7fa7ef4b475b
Nmemb = size(members, 1)

# ╔═╡ 60d0e593-88fd-4b4c-9009-cc24a597c6d5
members_nospace = best_stars[best_stars.PSAT_NOSPACE .> 0.2, :]

# ╔═╡ a9d94121-ea6e-416a-bae8-aa93c16bde72
md"""
# Utilities
"""

# ╔═╡ f3d3f7b2-30c2-4f89-bef5-ac1a8720963f
logit(p) = log(p / (1-p))

# ╔═╡ 7cc17173-160c-45f1-a685-c79271679414
xi_label = L"\xi\,/\,degrees"

# ╔═╡ 965217b8-b2a5-485b-83de-cac887065b19
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

# ╔═╡ 2d4da56b-5d0a-49d3-83ae-a90f85192101
θmax = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))

# ╔═╡ 5a9b5529-50f1-4cb5-b716-42180ea77d5f
md"""
# Simple selection plots
"""

# ╔═╡ 77f69d97-f71a-48e9-a048-1bb520222855
md"""
The plots below show various subsamples of the dataset in different planes to get a handle on how generally members are selected using the J+24 algorithm.
"""

# ╔═╡ 12aa5708-dfee-4b48-8835-69055ad82680
md"member markersize: $(@bind member_markersize confirm(PlutoUI.NumberField(1:10, default=5)))"

# ╔═╡ 57f558c1-ca31-44bb-a681-a2ed083a0b70
md"""background: 
$(@bind bg_markersize confirm(PlutoUI.NumberField(0.5:0.1:3; default=8/θmax^2)))
"""

# ╔═╡ e4d42cbd-0166-43d8-838c-27a14956623b
datasets = OrderedDict(
	:best => best_stars,
	:members => members,
)

# ╔═╡ 6d343a21-5bf9-4329-a65d-f6ad30c00c5e
scatter_kwargs = Dict(
	:best => (;	alpha=0.1, markersize=bg_markersize, color=:black, 
		label="all" => (alpha=1, markersize=3),
		rasterize=true,
	),
	:members => (;
		alpha=1, markersize=member_markersize,
		label = "members",
		color=COLORS[2]
	),
)

# ╔═╡ 0c12931c-5fa2-42d3-a129-47a3f55650fd


# ╔═╡ 0e6847cb-e246-49c1-9cb4-0f00f74446cf
sum((all_stars.PSAT .> 0.2) .& (all_stars.F_BEST .== 0))

# ╔═╡ 6c4abf21-34b8-4cec-bc4e-4d221175c624
sum(.!isnan.(all_stars.PSAT) .& (all_stars.F_BEST .== 0))

# ╔═╡ 1473e2d9-6153-4bd3-aa31-c40cf31f300b
sum(isnan.(all_stars.PSAT) .& (all_stars.F_BEST .== 1))

# ╔═╡ 77b7678f-158f-46cb-82b1-2e719ec6895a
function compare_samples(datasets, scatter_kwargs) 
	fig = Figure(
		size = (400, 400),
		fontsize=12,
	)

	# tangent
	dθ = 60θmax
	ax = Axis(fig[1, 1], 
		xlabel=plot_labels[:xi_am], 
		ylabel=plot_labels[:eta_am], 
		aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false
	)


	for (label, df) in datasets
		scatter!(60df.xi, 60df.eta; scatter_kwargs[label]...)
	end

	# cmd
	ax =  Axis(fig[1,2], 
		yreversed=true,
		xlabel=plot_labels[:bp_rp],
		ylabel=plot_labels[:G],
		limits=(-0.5, 2.5, 15, 22),
		xgridvisible=false,
		ygridvisible=false,
		aspect = 1,
	)

	for (label, df) in datasets
		scatter!(df.bp_rp, df.phot_g_mean_mag; scatter_kwargs[label]...)
	end

	# proper motions
	ax = Axis(fig[2,1],
		xlabel = plot_labels[:pmra],
		ylabel = plot_labels[:pmdec],
		aspect=DataAspect(),
		limits=(-10, 10, -10, 10),
		xgridvisible=false,
		ygridvisible=false,
		)
	
	for (label, df) in datasets
		scatter!(df.pmra, df.pmdec; scatter_kwargs[label]...)
	end


	# parallax
	ax = Axis(fig[2,2],
		xlabel=L"\varpi\ /\ \textrm{mas}", 
		ylabel=L"\delta\varpi\ /\ \textrm{mas}",
		limits=(-5, 5, 0, 3),
		aspect = 1,
		xgridvisible=false, ygridvisible=false,
	)

	for (label, df) in datasets
		scatter!(df.parallax, df.parallax_error; scatter_kwargs[label]...)
	end
	
	
	fig
end

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
@savefig "membership_selection_all" compare_samples(
		(
		:all => all_stars,
		:members => members,
	),
	Dict(
		:all => (;	alpha=0.1, markersize=bg_markersize, color=:black, 
			label="all" => (alpha=1, markersize=3),
			rasterize=true,
		),
		:members => (;
			alpha=1, markersize=member_markersize,
			label = "members",
			color=COLORS[2]
		),
	)
)

# ╔═╡ dde9398a-26b4-4a66-a396-f5b6188f9674
@savefig "membership_selection_best" compare_samples(
		(
		:best => best_stars,
		:members => members,
	),
	Dict(
		:best => (;	alpha=0.1, markersize=2bg_markersize, color=:black, 
			label="all" => (alpha=1, markersize=3),
			rasterize=true,
		),
		:members => (;
			alpha=1, markersize=member_markersize,
			label = "members",
			color=COLORS[2]
		),
	)
)

# ╔═╡ a4042ba0-f108-4025-bd99-04e32e687287
compare_samples(
		(
		:best => best_stars,
		:nospace => members_nospace,
		#:members => members,
	),
	Dict(
		:best => (;	alpha=0.1, markersize=2bg_markersize, color=:black, 
			label="all" => (alpha=1, markersize=3),
			rasterize=true,
		),
		:nospace => (;
			alpha=1, markersize=member_markersize * 1.1,
			label = "nospace",
			color=COLORS[1]
		),
		:members => (;
			alpha=1, markersize=member_markersize,
			label = "members",
			color=COLORS[2]
		),
	)
)

# ╔═╡ 722a1f29-be6b-4fca-a9b3-8c304b6c4eb6
let 
	dθ = 60θmax
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel=plot_labels[:xi_am], 
		ylabel=plot_labels[:eta_am], 
		aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false
	)


	scatter!(60best_stars.xi, 60best_stars.eta, 
		alpha=0.2, markersize=1, color=:black, 
		label="best" => (alpha=1, markersize=3)
	)
	
	scatter!(60members_nospace.xi, 60members_nospace.eta, 
		alpha=1, markersize=3, 
		label="members w/o spatial"
	)
	
	scatter!(60members.xi, 60members.eta, 
		alpha=1, markersize=member_markersize,
		label = "members"
	)

	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ 9c1947ef-b9d2-4a48-8e0e-f8d1b8e95124
md"""
## Selected sample checks
"""

# ╔═╡ 63fd1adf-3d9a-4a37-82b5-4f69df85c2a9
md"""
The KDE plot below is mostly another way of visualizing the sample distribution, if there might be anything strange going on. Reminder that we only have $Nmemb stars
"""

# ╔═╡ a5d834ed-6fd2-4f90-96a7-df393aca541e
md"""
Plot of members on sky, should look similar to tangent plane plot / kde above
"""

# ╔═╡ 4a1bf5cb-6083-429f-8747-9ed0478f8f3e
md"""
If there is crowding, the plot below should appear to have deviations in brightness distributions. Otherwise, the distribution should look like a scatterplot
"""

# ╔═╡ 02d19c07-411f-40d6-9ac3-010ebfd4bdfe
let 
	dθ = θmax*60 / 2
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel=L"\xi / \textrm{arcmin}", ylabel=L"\eta / \textrm{arcmin}",
		aspect=1, limits=(-dθ, dθ, -dθ, dθ), xgridvisible=false, ygridvisible=false)
	
	h = scatter!(60*all_stars.xi, 60*all_stars.eta, color=all_stars.phot_g_mean_mag, 
		colormap=:greys, markersize=3)

	Colorbar(fig[1,2], h, label="G mag")
	fig
end

# ╔═╡ ca63e89f-7ee1-49d9-8bc1-afcb755f9ebf
md"""
CMD of members only
"""

# ╔═╡ 83b506bb-f385-464e-8f5c-7adfff45105a
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
	title = "members")

	
	h = scatter!(members.bp_rp, members.phot_g_mean_mag, color=log10.(members.r_ell), 
	)

	Colorbar(f[1, 2], h, label="log r ell")

	LilGuys.hide_grid!(ax)
	f
end

# ╔═╡ 4056fc57-0058-4da5-89af-e372164b2c64
md"""
## Ellipticity evidence
"""

# ╔═╡ c6115ac4-57de-43f1-abfe-adfde1294bb5
function plot_circle!(radius, x=0, y=0; kwargs...)
	t = LinRange(0, 2π, 10000)
	lines!(x .+ radius*cos.(t), y .+ radius*sin.(t); kwargs...)
end

# ╔═╡ 042f8215-b6fe-49db-a166-055567da4bb0
begin 
	θ_bins = LinRange(-π, π, 15)
	θ_bm = lguys.midpoints(θ_bins)
end

# ╔═╡ 64e6f0e6-aca9-46db-bc7c-28cdadb54be2
md"""
If the object was circular, then the total mass in each bin in theta should be about constant. Otherwise, the plot below should be some sine-wave object
"""

# ╔═╡ 117cb774-fe06-4a87-9b81-465790e0b288
md"""
Do the members appear to be distributed anisotropically?
"""

# ╔═╡ 0bea6192-76e4-4d92-a5f6-52e74e4634de
md"""
This plot should be identical to above, just in polar coordinates
"""

# ╔═╡ 3e40a3f8-efbc-4807-b185-22fbb2e99adf
begin 
	xi = members.xi .- cen[1]
	eta = members.eta .- cen[2]
end

# ╔═╡ 1942cac0-2da8-4bec-80d0-a9302ddc9ea1
let 
	global Σ_kde,kde_d, x_kde, y_kde

	kde_d = kde((60xi, 60eta))
	ik = InterpKDE(kde_d)


	r_max = max(abs.(extrema(kde_d.x))..., abs.(extrema(kde_d.y))..., )
	x_kde = LinRange(-r_max, r_max, 1000)
	y_kde = LinRange(-r_max, r_max, 1000)
	Σ_kde = @. pdf([ik], x_kde', y_kde)
	Σ_kde = Σ_kde'
	;
end

# ╔═╡ d49c5f63-e2a4-4ff5-af4a-d4239b50abae
begin
	kde_y = reshape(y_kde, (1, :))
	kde_r = @. sqrt((x_kde')^2 + (y_kde)^2)
	kde_theta = atan.(y_kde, x_kde')
end

# ╔═╡ d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "theta",
		ylabel = "KDE mass / bin"
	)
	hist!(vec(kde_theta), θ_bins, weights=vec(Σ_kde))

	vlines!(observed_properties["position_angle"] * π/180, color=COLORS[2])

	LilGuys.hide_grid!(ax)
	fig
end

# ╔═╡ cdd246ae-e6bb-49a9-9512-085bd37827bd
let 
	density_scale = 1e-2
	contour_levels = 10

	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel = "xi / arcmin",
		ylabel = "eta / arcmin",
		aspect = DataAspect()
	)

	p = heatmap!(x_kde, y_kde, Σ_kde, colormap=:greys, colorscale=x -> asinh(x/density_scale))
	ax.aspect = 1
	contour!(x_kde, y_kde, asinh.(Σ_kde ./ density_scale), levels=contour_levels )


	filt = members.r_ell .> 0
	scatter!(ax, 60xi[filt], 60eta[filt], markersize=3, color=Arya.COLORS[3])

	Colorbar(fig[1, 2], p, ticks=Makie.automatic, minorticks=Makie.automatic,
		label = "density / $density_scale"
	)

	@savefig "kde_density_2d"
	fig
end

# ╔═╡ ec1c1408-5c3c-4bfa-b8bd-f32dbc3379a9
ϕ = atan.(eta, xi)

# ╔═╡ 9282a867-b9c7-4a50-b58a-da3ceb068cad
radii = @. sqrt(xi^2 + eta^2)

# ╔═╡ 4f538de6-8827-454f-a97f-8f6c2cd7ea3f
let
	fig = Figure()

	# limits based on full range of data
	ra = observed_properties["ra"]
	dec = observed_properties["dec"]
	r_max = maximum(radii)
	
	ax = Axis(fig[1,1],
		xlabel="ra / degrees", ylabel="dec / degrees", 
		title="members only", 
		aspect=  1, 
		xgridvisible=false, ygridvisible=false, 
		limits=(ra - r_max/cosd(dec), ra+r_max/cosd(dec), dec-r_max, dec+r_max)
	)
		
	scatter!(members.ra, members.dec)

	fig
end

# ╔═╡ 9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
let
	dθ = maximum(radii)

	f = Figure()
	ax = Axis(f[1, 1], limits=(-dθ,dθ,-dθ,dθ), aspect=1,
		xlabel=L"$\xi$ / degrees",
		ylabel=L"$\eta$ / degrees")
	
	LilGuys.hide_grid!(ax)
	
	for a in LinRange(0.05dθ, dθ, 8)
		plot_circle!(a, color=:black, alpha=0.1)
	end
	
	scatter!(members.xi, members.eta, markersize=3, label="", alpha=1)
	
	f
end

# ╔═╡ 52d3db7f-8bbf-43cf-89f2-16134247d713
let
	fig = Figure()
	ax = PolarAxis(fig[1,1], rlimits=(0, maximum(radii)))

	scatter!(ax,  ϕ, radii)
	fig
end

# ╔═╡ 4ff675ef-8609-49df-bfac-4070f14e3c25
maximum(radii)

# ╔═╡ d2668f44-9cb7-4ab8-af49-1aef2baecdda
md"""
# Sanity checks
"""

# ╔═╡ 12bf64e9-ae49-4f93-8572-660dbc86ffc1
scatter(all_stars.r_ell, all_stars.r_ell_old * observed_properties["r_h"] * sqrt(1 - observed_properties["ellipticity"]), 
	rasterize=true,
	axis = (;
	xlabel = "r ell",
	ylabel = "r ell jax",
	aspect = DataAspect(),
	)
)

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═8b3955e5-ab08-4d5a-abe5-8f40120ff8db
# ╠═da9ca7a7-18b8-49cb-a041-ab1c667920ff
# ╠═489f6d21-9e9a-4b1e-b05c-c63a44ba1951
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═db1264b7-02c3-4b55-ae2b-9ce78aa1304a
# ╠═81c677bf-8372-4fd3-9d60-af8377485f1c
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═53c0567a-c845-486d-9510-1590cc0ee9ff
# ╠═cf7aeb0a-d453-462a-b43d-a832567440fd
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═65822789-b261-45be-817f-cba0e37b82ce
# ╠═4c8b9c30-f873-4756-a729-aa023b3a804e
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═d87f28e5-cc94-4afc-9e70-ece1686db22c
# ╠═338be0b8-e732-4202-a20f-9ca841516075
# ╟─7ed70894-48a8-46af-8c78-8a75cb39d957
# ╠═223abc41-5158-49c2-96bf-df55b7be1114
# ╟─3abea230-2b95-40d3-9851-a91236f75c4a
# ╠═32e8e100-1fad-425f-b05c-2310e7ba559c
# ╠═b1be5f97-0e61-4dda-87fd-7bb0f43147b6
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═60d0e593-88fd-4b4c-9009-cc24a597c6d5
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╟─a9d94121-ea6e-416a-bae8-aa93c16bde72
# ╠═f3d3f7b2-30c2-4f89-bef5-ac1a8720963f
# ╠═7cc17173-160c-45f1-a685-c79271679414
# ╠═965217b8-b2a5-485b-83de-cac887065b19
# ╠═2d4da56b-5d0a-49d3-83ae-a90f85192101
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═12aa5708-dfee-4b48-8835-69055ad82680
# ╠═57f558c1-ca31-44bb-a681-a2ed083a0b70
# ╠═e4d42cbd-0166-43d8-838c-27a14956623b
# ╠═6d343a21-5bf9-4329-a65d-f6ad30c00c5e
# ╠═0c12931c-5fa2-42d3-a129-47a3f55650fd
# ╠═0e6847cb-e246-49c1-9cb4-0f00f74446cf
# ╠═6c4abf21-34b8-4cec-bc4e-4d221175c624
# ╠═1473e2d9-6153-4bd3-aa31-c40cf31f300b
# ╠═77b7678f-158f-46cb-82b1-2e719ec6895a
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
# ╠═dde9398a-26b4-4a66-a396-f5b6188f9674
# ╠═a4042ba0-f108-4025-bd99-04e32e687287
# ╟─722a1f29-be6b-4fca-a9b3-8c304b6c4eb6
# ╟─9c1947ef-b9d2-4a48-8e0e-f8d1b8e95124
# ╠═1942cac0-2da8-4bec-80d0-a9302ddc9ea1
# ╠═d49c5f63-e2a4-4ff5-af4a-d4239b50abae
# ╟─63fd1adf-3d9a-4a37-82b5-4f69df85c2a9
# ╠═cdd246ae-e6bb-49a9-9512-085bd37827bd
# ╠═a5d834ed-6fd2-4f90-96a7-df393aca541e
# ╟─4f538de6-8827-454f-a97f-8f6c2cd7ea3f
# ╟─4a1bf5cb-6083-429f-8747-9ed0478f8f3e
# ╟─02d19c07-411f-40d6-9ac3-010ebfd4bdfe
# ╟─ca63e89f-7ee1-49d9-8bc1-afcb755f9ebf
# ╟─83b506bb-f385-464e-8f5c-7adfff45105a
# ╟─4056fc57-0058-4da5-89af-e372164b2c64
# ╠═c6115ac4-57de-43f1-abfe-adfde1294bb5
# ╠═042f8215-b6fe-49db-a166-055567da4bb0
# ╟─64e6f0e6-aca9-46db-bc7c-28cdadb54be2
# ╠═d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
# ╟─117cb774-fe06-4a87-9b81-465790e0b288
# ╟─9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
# ╟─0bea6192-76e4-4d92-a5f6-52e74e4634de
# ╟─52d3db7f-8bbf-43cf-89f2-16134247d713
# ╠═3e40a3f8-efbc-4807-b185-22fbb2e99adf
# ╠═4ff675ef-8609-49df-bfac-4070f14e3c25
# ╠═ec1c1408-5c3c-4bfa-b8bd-f32dbc3379a9
# ╠═9282a867-b9c7-4a50-b58a-da3ceb068cad
# ╟─d2668f44-9cb7-4ab8-af49-1aef2baecdda
# ╠═12bf64e9-ae49-4f93-8572-660dbc86ffc1
