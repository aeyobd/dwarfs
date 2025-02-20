### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

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

# ╔═╡ 69c98029-165c-407b-9a63-a27e06e30e45
using Arya

# ╔═╡ ae29bed0-6700-47f1-8952-35e867ce126b
using OrderedCollections

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sammple

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
 galaxyname = "draco"

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using LilGuys
	FIGDIR = joinpath(galaxyname, "figures")
	FIGSUFFIX = ".jax_observations"
end

# ╔═╡ cf7aeb0a-d453-462a-b43d-a832567440fd
diverging_cmap = Reverse(:bluesreds)

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

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ b05bb3ea-2492-4da0-bfd5-8aaf5bf936c4
if isdir(galaxyname)
	if ! isdir(FIGDIR)
		mkdir(FIGDIR)
	end
end

# ╔═╡ 9ecf79a8-2ed3-40c6-b555-a102250ecbd4
observed_properties = TOML.parsefile(galaxyname * "/observed_properties.toml")

# ╔═╡ 26cf1867-02be-4d36-8c35-6c58a1feca27
if isfile("$galaxyname/data/j24_2c.fits")
	datafile = "$galaxyname/data/j24_2c.fits"
else
	datafile = "$galaxyname/data/j24_1c.fits"
end

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	all_stars = lguys.read_fits(datafile)

	all_stars[:, :r_ell_old] = all_stars[:, :r_ell]
	all_stars.xi, all_stars.eta = lguys.to_tangent(all_stars.ra, all_stars.dec, observed_properties["ra"], observed_properties["dec"])

	all_stars.r_ell = 60lguys.calc_r_ell(all_stars.xi, all_stars.eta, observed_properties["ellipticity"], observed_properties["position_angle"])
	all_stars[!, :LL_S] = @. log10(all_stars.L_S_SAT / all_stars.L_S_BKD)
	all_stars[!, :LL_PM] = @. log10(all_stars.L_PM_SAT / all_stars.L_PM_BKD)
	all_stars[!, :LL_CMD] = @. log10(all_stars.L_CMD_SAT / all_stars.L_CMD_BKD)
	all_stars[!, :LL] = @. log10(all_stars.L_SAT / all_stars.L_BKD)

	all_stars[!, :LL_norm] = @. abs(all_stars.LL_S) + abs(all_stars.LL_PM) + abs(all_stars.LL_CMD)

	all_stars[!, :f_LL_S] = @. all_stars.LL_S / all_stars.LL
	all_stars[!, :f_LL_PM] = @. all_stars.LL_PM / all_stars.LL
	all_stars[!, :f_LL_CMD] = @. all_stars.LL_CMD / all_stars.LL

	all_stars[!, :f_LL_S_norm] = @. all_stars.LL_S / all_stars.LL_norm
	all_stars[!, :f_LL_PM_norm] = @. all_stars.LL_PM / all_stars.LL_norm
	all_stars[!, :f_LL_CMD_norm] = @. all_stars.LL_CMD / all_stars.LL_norm

	all_stars
end

# ╔═╡ 7ed70894-48a8-46af-8c78-8a75cb39d957
md"""
estimate of centring uncertanty in this dataset (standard error on mean and present difference between mean centre & assumed centre)

"""

# ╔═╡ 3abea230-2b95-40d3-9851-a91236f75c4a
md"""
each component of the cen_errs tuple below should be similar in magnitude
"""

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = all_stars[all_stars.PSAT .> 0.2, :]

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
members_nospace = all_stars[all_stars.PSAT_NOSPACE .> 0.2, :]

# ╔═╡ 082a06dd-eeb5-4761-a233-1ee89e8cb819
best_stars = all_stars[all_stars.F_BEST .== 1.0, :]

# ╔═╡ a9d94121-ea6e-416a-bae8-aa93c16bde72
md"""
# Utilities
"""

# ╔═╡ f3d3f7b2-30c2-4f89-bef5-ac1a8720963f
logit(p) = log(p / (1-p))

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

# ╔═╡ 722a1f29-be6b-4fca-a9b3-8c304b6c4eb6
let 
	dθ = 60θmax
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel=L"\xi / \textrm{arcmin}", ylabel=L"\eta / \textrm{arcmin}", aspect=1, limits=(-dθ, dθ, -dθ, dθ), xgridvisible=false, ygridvisible=false)


	scatter!(60best_stars.xi, 60best_stars.eta, 
		alpha=0.2, markersize=1, color=:black, 
		label="best" => (alpha=1, markersize=3)
	)
	
	scatter!(60members_nospace.xi, 60members_nospace.eta, 
		alpha=1, markersize=3, 
		label="members w/o spatial"
	)
	
	scatter!(60members.xi, 60members.eta, 
		alpha=1, markersize=3,
		label = "members"
	)

	Legend(fig[1, 2], ax)

	@savefig "xi_eta_members"
	fig
end

# ╔═╡ b6b8d7e7-d0a1-45bc-8593-66aea4a65637
let 
	fig = Figure()
	ax =  Axis(fig[1,1], 
		yreversed=true,
		xlabel="bp - rp",
		ylabel = "G",
		limits=(-0.5, 3, 10, 22),
		xgridvisible=false,
		ygridvisible=false
	)

	scatter!(all_stars.bp_rp, all_stars.phot_g_mean_mag, 
		color=:black, markersize=1, alpha=0.1,
		rasterize=true, 
		label="all" => (;markersize=5, alpha=1)
	)

	scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, 
		color=COLORS[1], alpha=0.1, markersize=3, 
		rasterize=true, 
		label="best" => (markersize=5, alpha=1)
	)

	scatter!(members.bp_rp, members.phot_g_mean_mag, 
		color=COLORS[2], 
		label="members"
	)

	axislegend(ax, position=:rt)
	@savefig "cmd_memb"

	fig
end

# ╔═╡ a4e857c0-39d9-4fa0-871f-ccc66cb17c25
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\varpi\ /\ \textrm{mas}", 
		ylabel=L"\delta\varpi\ /\ \textrm{mas}",
		limits=(-4, 5, 0, 3)
	)

	scatter!(all_stars.parallax, all_stars.parallax_error, 
		color=:black, markersize=1, alpha=0.1, 
		label = "all stars" => (;alpha=1),
		rasterize=true
	)

	scatter!(best_stars.parallax, best_stars.parallax_error,
		markersize=3, alpha=0.1,
		rasterize=true, 
		label="best",
	)

	scatter!(members.parallax, members.parallax_error, 
		label="members"
	)


	LilGuys.hide_grid!(ax)
	axislegend()

	@savefig "pi_pierr_memb"
	fig
end

# ╔═╡ d27062d3-9374-4342-bb4e-3461c774da79
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		limits=(-10, 10, -10, 10),
		)
	
	p = scatter!((all_stars.pmra), (all_stars.pmdec), 
		markersize=1, color=:black, alpha=0.1,
		label="all" => (markersize=5, alpha=1),
		rasterize=true
	)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), 
		markersize=2,
		alpha=0.1,
		label="best" => (markersize=5, alpha=1),
		rasterize=true
	)

	p = scatter!((members.pmra), (members.pmdec), 
		markersize=5,
		label="members",
		color=COLORS[2]
	)

	Legend(fig[1,2], ax)
	
	LilGuys.hide_grid!(ax)

	@savefig "pm_members"

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

	Colorbar(fig[1, 2], p, #ticks=[10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0],
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

# ╔═╡ b0aec8b4-410f-4678-ad8c-cb055b743bd3
md"""
# Likelihood components
"""

# ╔═╡ 3c2ee4bb-d594-47e7-a7f3-5b7c9e38b7d4
md"""
## Likelihoods maps
"""

# ╔═╡ b191b7b8-28ed-4aee-b376-b2bdf22abc29
md"""
These maps are just for reference, to see how strongly each component of the likelihood varies in that parameter space.
"""

# ╔═╡ b06445ed-6290-42e8-bdb8-20f0372518c2
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi / degrees",
		ylabel = "eta / dgrees",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.xi), (best_stars.eta), color=(best_stars.L_S_BKD),
		colormap=:bluesreds,
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S_BKD", )
	
	fig
end

# ╔═╡ 2a8e13bd-b8ca-4cf3-b61b-096af218265b
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.xi), (best_stars.eta), color=log10.(best_stars.L_S_SAT),
		markersize=3,
		#colorrange=(-1, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S_SAT", )
	
	fig
end

# ╔═╡ 25dd2e38-1c78-4ef6-95ef-d38a0e80f0e4
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
)
	LilGuys.hide_grid!(ax)

	
	h = scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, color=log10.(best_stars.L_CMD_SAT), 
		markersize = 3
	)

	Colorbar(f[1, 2], h, label="log L_CMD_SAT")
	f
end

# ╔═╡ 1b63c1ee-0a46-49c6-880e-c1742055605e
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
)

	
	h = scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, color=log10.(best_stars.L_CMD_BKD), 
		markersize = 3
	)

	LilGuys.hide_grid!(ax)
	Colorbar(f[1, 2], h, label="log L_CMD_BKD")
	f
end

# ╔═╡ 643773f0-2b03-46f7-8df4-579b3c708909
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
		xlabel="bp - rp",
		ylabel = "G",
	)
	LilGuys.hide_grid!(ax)
	
	h = scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, color=log10.(best_stars.L_CMD_SAT ./ best_stars.L_CMD_BKD), 
		colorrange=(-10, 1),
		markersize=3
	)

	Colorbar(f[1, 2], h, label="log L_CMD_SAT / L_CMD_BKD")
	f
end

# ╔═╡ cad88ab4-1ada-4cc3-8ef0-d4674500d57e
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), color=log10.(best_stars.L_PM_SAT),
		colorrange=(-30, 1),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_PM_SAT", )
	
	fig
end

# ╔═╡ 18252c3c-d4d6-46f2-8277-33a700ab4074
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), color=log10.(best_stars.L_PM_BKD),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_PM_BKD", )
	
	fig
end

# ╔═╡ d47a8d02-6c63-41fd-95b5-85550d1f5a85
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), color=log10.(best_stars.L_PM_SAT ./ best_stars.L_PM_BKD),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_PM_SAT / L_PM_BKD", )
	
	fig
end

# ╔═╡ b7e7e245-3732-453c-a2ec-3fe84c65c012
md"""
## Likelihood ratios
"""

# ╔═╡ 87f29963-2c7f-4d8f-8631-f70197d315c0
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		title = "LLR > -5"
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_S) .- (df.LL_norm),
		markersize=3,
		colorrange=(-10, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S / L_tot_norm", )
	
	fig
end

# ╔═╡ 72aa995b-8dd1-4d08-88c2-5fb527ff7f3e
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		title = "LLR > -10",
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_S) .- (df.LL),
		markersize=3,
		colorrange=(-10, 10),
		colormap=diverging_cmap,
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S / L_tot", )
	
	fig
end

# ╔═╡ 6cdcd105-fa1d-48f5-9861-b722aafa28a3
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_PM) .- (df.LL),
		markersize=3,
		colorrange=(-10, 10),
		colormap=diverging_cmap,
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L pm / L tot", )
	
	fig
end

# ╔═╡ cf8e51a9-194c-4aa2-9766-39256a35ef02
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -5, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_PM) .- (df.LL_norm),
		markersize=3,
		colorrange=(-10, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L pm / L tot norm", )
	
	fig
end

# ╔═╡ e265ec56-d62d-4295-975f-98ba93cf5e16
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_CMD) .- (df.LL),
		markersize=3,
		colorrange=(-10, 10),
		colormap=diverging_cmap,
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L cmd / L tot", )
	
	fig
end

# ╔═╡ e2686165-6499-4558-9da9-91a78624b992
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -5, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_CMD) .- (df.LL_norm),
		markersize=3,
		colorrange=(-10, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L cmd / L tot norm", )
	
	fig
end

# ╔═╡ 1f2aa1ee-b395-4040-84cb-973d57b534d5
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "Log Likelihood",
		ylabel = "counts",
	)

	filt = isfinite.(best_stars.LL)

	stephist!(log10.(best_stars.L_S_SAT[filt]), 
			label="spatial")
	
	stephist!(log10.(best_stars.L_CMD_SAT[filt]),  label="cmd")

	stephist!(log10.(best_stars.L_PM_SAT[filt]),  label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ ca06c336-3950-4831-9f3f-7fe41406bc08
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log likilihood",
		ylabel = "counts",
		title = "members"
	)
	bins = LinRange(-10, 2, 30)
	
	stephist!(log10.(members.L_S_SAT), 
		bins=bins, label="spatial")
	
	stephist!(log10.(members.L_CMD_SAT), bins=bins, label="cmd")
	
	stephist!(log10.(members.L_PM_SAT), bins=bins,  label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ 63a3e02b-0ea0-472d-8c2e-962bd8e770c0
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log likelihood sat / background",
		ylabel = "counts",
		title = "members"
	)
	
	stephist!(log10.(members.L_S_SAT ./ members.L_S_BKD), 
		label="spatial")
	
	stephist!(log10.(members.L_CMD_SAT ./ members.L_CMD_BKD), label="cmd")
	
	stephist!(log10.(members.L_PM_SAT ./ members.L_PM_BKD), label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ 3a52eba9-e07d-4445-a322-1c4b8bb33024
md"""
### Correlations
"""

# ╔═╡ f54c50ae-0e03-4917-899f-81a30f116b3f
function plot_corr(best_stars; markersize=3, title="")
	fig = Figure(
		size = (800, 400)
	)

	kwargs = (;
		markersize=markersize,
	)
	
	local p

	axs = []
	
	for i in 1:3
		col = [:LL_S, :LL_CMD, :LL_PM][i]
		label = ["space", "cmd", "pm"][i]

		ax = Axis(fig[1,i],
			xlabel = L"\log\mathcal{L}_\textrm{%$label}",
			ylabel = L"\log\mathcal{L}",
			aspect=DataAspect(),
			title = i == 2 ? title : ""
		)
		
		
		p = scatter!(best_stars[!, col], best_stars.LL;
			color= best_stars.PSAT ,
			colorrange=(0, 1),
			kwargs...
		)


		# correlation coefficient
		R = cor(best_stars[!, col], best_stars.LL)
		t = "ρ = $(round(R, digits=2))"
		text!(0.1, 0.9, text=t, space=:relative)

		# formatting
		if i > 1
			hideydecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
		LilGuys.hide_grid!(ax)
		push!(axs, ax)

	end

	Colorbar(fig[1,4], p, label = "PSAT")

	linkaxes!(axs...)
	fig
end

# ╔═╡ b4f816f0-0abe-43c2-bf5b-48a2eb2b99c8
plot_corr(members, title="members")

# ╔═╡ dc280eae-4f10-4167-9545-c70d67af2d8f
let
	fig = plot_corr(best_stars[best_stars.LL .> 0, :], markersize=5, title = "LLR > 0")

	@savefig "LL_corr_trunc"

	fig
end

# ╔═╡ fb3a27be-62c7-41ab-ba5d-761efb74c044
let
	fig = plot_corr(best_stars[best_stars.LL .> -10, :], title="LLR > -10")

	fig
end

# ╔═╡ b4363369-4590-4446-9f4a-9f71ff31a870
plot_corr(best_stars[best_stars.LL .> -20, :], title = "LLR > -20")

# ╔═╡ 00d4679d-6c55-4929-bbd3-c796479ba70d
plot_corr(all_stars[all_stars.LL .> -80, :], title = "LLR > -80")

# ╔═╡ 6be3bd9d-c257-444f-87ab-2b3956fc83aa
md"""
## Fractional contributions
"""

# ╔═╡ 4ad0330e-f1d7-4d62-b6f9-f49be1aa0d80
md"""
Still not sure how useful these plots are, but they should represent an approximate metric for how substantially the log-likelihood is influenced by any component.
"""

# ╔═╡ 9a2ca944-dd5c-43f5-b70d-ab3d67b42743
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "LL",
		ylabel = "f spatial",
		limits=(-30, 10, -1.5, 3)
	)
	
	scatter!(best_stars.LL, (best_stars.f_LL_S), markersize=3, label = "best")
	scatter!(members.LL, (members.f_LL_S), label = "members")

	axislegend()
	fig
end

# ╔═╡ 79ac654d-7de9-4a2c-bc5d-2f8c8a7e057f
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "LL",
		ylabel = "f spatial normed",
		limits=(-30, 10, -1, 1)
	)
	
	scatter!(best_stars.LL, (best_stars.f_LL_S_norm), markersize=3)
	scatter!(members.LL, (members.f_LL_S_norm))

	fig
end

# ╔═╡ 479c9cad-8cae-4d9b-87d3-79a816b3af82
md"""
The three plots below verify that the normalized version of likelihood fractions makes sense as a diagnostic.
"""

# ╔═╡ a7dbb8bf-8c82-4ace-9f17-e9aa159a84ff
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f spatial",
		ylabel = "f spatial normed",
		limits=(0, 1.5, 0, 1.5),
		aspect = DataAspect(),
	)
	
	scatter!(abs.(best_stars.f_LL_S), abs.(best_stars.f_LL_S_norm), markersize=3)
	scatter!(abs.(members.f_LL_S), abs.(members.f_LL_S_norm))

	fig
end

# ╔═╡ ec4a8ca8-a4b9-41d7-a369-8049b2a98795
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f CMD",
		ylabel = "f CMD normed",
		limits=(0, 3, 0, 1.5),
		aspect = DataAspect(),
	)
	
	scatter!(abs.(best_stars.f_LL_CMD), abs.(best_stars.f_LL_CMD_norm), markersize=3)
	scatter!(abs.(members.f_LL_CMD), abs.(members.f_LL_CMD_norm))

	fig
end

# ╔═╡ 32fd5ed6-fb41-456c-bd05-d28ad013c2c7
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f PM",
		ylabel = "f PM normed",
		limits=(0, 3, 0, 1.5),
		aspect = DataAspect(),
	)
	
	scatter!(abs.(best_stars.f_LL_PM), abs.(best_stars.f_LL_PM_norm), markersize=3)
	scatter!(abs.(members.f_LL_PM), abs.(members.f_LL_PM_norm))

	fig
end

# ╔═╡ 301f42ce-437f-4d45-a186-42df47303f25
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f PM",
		ylabel = "f PM unnormed",
		limits=(-10, 10, -10, 10),
		aspect = DataAspect(),
	)
	
	scatter!(best_stars.LL_PM .- best_stars.LL_norm,  best_stars.LL_PM .- best_stars.LL, markersize=3)

	fig
end

# ╔═╡ 625795e7-aca7-4bb9-8db6-cfcfb87bddaf
import DensityEstimators as DE

# ╔═╡ 3c05db74-0cc1-4893-a468-b495c0121118
md"""
## Fractional contributions on the sky
"""

# ╔═╡ e3f951ed-3725-4ece-bf59-cb6998fdf170
function plot_f_l_tangent(best_stars; title="", markersize=5)
	fig = Figure(
		size=(800, 500)
	)

	filt = best_stars.PSAT .> 0.2
	df = best_stars[filt, :]


	for i in 1:3
		col = [:f_LL_S_norm, :f_LL_PM_norm, :f_LL_CMD_norm][i]
		label = ["space", "proper motion", "CMD"][i]
		
		ax = Axis(fig[1,i],
			title = label,
			xlabel = "f likelihood",
			ylabel = "counts",
			limits=(-1, 1, 0, nothing),
		)
		LilGuys.hide_grid!(ax)

		h = histogram(df[:, col])
		
		x = [h.bins[1]; midpoints(h.bins); h.bins[end]]
		y = [0; h.values; 0]
		lines!(x, y, color=x,
			colorrange=(-1, 1),
			colormap=diverging_cmap,
		)

		
		ax = Axis(fig[2,i],
			xlabel = "ξ / degree",
			ylabel = "η / degree",
			aspect=DataAspect(),
			#limits = (-2, 2, -2, 2)
			title = i == 2 ? title : ""
			)
		
		p = scatter!(df.xi, df.eta,
			markersize=markersize,
			color = df[:, col],
			colorrange=(-1, 1),
			colormap=diverging_cmap,
		)
		
		LilGuys.hide_grid!(ax)
	end



	rowsize!(fig.layout, 1, Relative(0.2))
	rowgap!(fig.layout, 1, 0)

	return fig
end

# ╔═╡ e4e6239f-9fdd-433b-8a39-ffaeaf96b59d
let
	fig = plot_f_l_tangent(members, title="members")

	@savefig "xi_eta_vs_f"
	fig
end

# ╔═╡ 363dfc70-ec52-44ed-ade4-f88322247871
let
	fig = plot_f_l_tangent(best_stars[best_stars.LL .> 0, :], title="LL > 0")

	@savefig "xi_eta_vs_f"
	fig
end

# ╔═╡ c68ff47f-97de-4bbb-8b0a-42d837da21f7
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = L"\log\,r_\textrm{ell}",
		ylabel = "f spatial",
	)
	
	p = scatter!(log10.(best_stars.r_ell), (best_stars.f_LL_S_norm), 
		color=(best_stars.LL),
		colorrange=(-10, 1),
		markersize=3,
	)
	LilGuys.hide_grid!(ax)
	
	Colorbar(fig[1,2], p, label="log Likelihood sat / background", )

	@savefig "f_space_vs_r_ell"
	fig
end

# ╔═╡ 8d497b75-a47d-42d1-9e54-e15bfb1db805
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "r_ell",
		ylabel = "f spatial",
	)
	
	scatter!((best_stars.r_ell), (best_stars.f_LL_S_norm), markersize=3)

	
	scatter!((members.r_ell), (members.f_LL_S_norm))

	fig
end

# ╔═╡ 21f80519-1215-4aa1-be36-92956903502e
md"""
# Density profiles
"""

# ╔═╡ a4d06da3-0197-4c34-88c9-29099ff9ab46
log_r_ell_label = L"\log\ r_\textrm{ell}\,/\,\textrm{arcmin}"

# ╔═╡ 83b506bb-f385-464e-8f5c-7adfff45105a
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
	title = "members")

	
	h = scatter!(members.bp_rp, members.phot_g_mean_mag, color=log10.(members.r_ell), 
	)

	Colorbar(f[1, 2], h, label=log_r_ell_label)

	LilGuys.hide_grid!(ax)
	f
end

# ╔═╡ 571b98f1-10d0-4012-8bfb-755f26f71b37
log_Σ_label = L"$\log \Sigma$\,/\,stars\ arcmin$^{-2}$"

# ╔═╡ b1724e25-6318-460e-a057-70b4588b3321
md"""
maximum r_ell is a proxy for the maximum semi-major axis of any ellipse in the field. Multiplying by (1-ellipticity) gives the largest semi-minor axis who's ellipse is entirely contained in the field.
"""

# ╔═╡ d1c46aea-e419-4feb-873f-ba035b46cda7
observed_properties["ellipticity"]

# ╔═╡ 642eb4bd-6082-4209-bc92-c4ef5ddaa3fc
r_ell_max = (1-observed_properties["ellipticity"]) * maximum(all_stars.r_ell)

# ╔═╡ 5cf14b86-0a30-410c-955a-82e9e8696e7b
maximum(all_stars.r_ell)

# ╔═╡ c6805532-116c-4e8c-9c40-81161a1334cd
maximum(sqrt.(all_stars.xi .^2 .+ all_stars.eta .^2))  * 60 * sqrt((1 - observed_properties["ellipticity"]))

# ╔═╡ f9edc1ee-2b05-4750-ba26-c12cd95d4041
maximum(sqrt.(all_stars.xi .^2 .+ all_stars.eta .^2)) 

# ╔═╡ 2d96b315-1207-4278-98e4-5b9239773779
function get_num_per_bin(x)
	num_per_bin = round(Int64, LilGuys.Interface.default_n_per_bin(x))
	if length(x) < 2
		num_per_bin = 1
	end
	return num_per_bin
end


# ╔═╡ a507dc20-8326-4f66-aebe-24ea16d4254e
function get_bins(r_ell; bin_width = 0.05, num_per_bin=nothing)

	if num_per_bin === nothing
		num_per_bin = get_num_per_bin(log10.(r_ell))
	end

	num_per_bin = round(Int, num_per_bin)
	bins = LilGuys.Interface.bins_both(log10.(r_ell), nothing, bin_width=bin_width, num_per_bin=num_per_bin)
	bins[end] = min(bins[end], log10(r_ell_max))

	return bins
end

# ╔═╡ 53bc3533-f9fc-4789-9c8e-e5085faacc33
md"""
calculates the background density
"""

# ╔═╡ 3a676557-390a-4042-ad9f-527c292a9da2
η_bins_bg = 0.3

# ╔═╡ d50bc64f-f67b-4498-aea8-81c4b5c258cc
r_ell_min_bg = 3

# ╔═╡ 76fa8518-bc9a-4126-880b-2945d6b8e853
let
	global log_Σ_bg
	global log_Σ_bg_err

	filt_bg = best_stars.r_ell .< r_ell_max
	filt_bg .&= best_stars.r_ell .> r_ell_min_bg

	N = sum(filt_bg)

	A = π * ( r_ell_max^2 - r_ell_min_bg^2)
	log_Σ_bg = log10.(N / A)

	x = best_stars.r_ell[filt_bg]
	bins = get_bins(x)
	prof_bg = LilGuys.StellarProfile(x, normalization=:none, bins=(bins))

	
	log_Σ_bg_err = lguys.std(prof_bg.log_Sigma) / sqrt(length(prof_bg.log_Sigma)) .+ 1/sqrt(N) / log(10)
end

# ╔═╡ 8c253d06-b902-4baa-82f7-6c507f8cf014
function plot_density_for_filter(filters; sequential=false, axislabel="")
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	names = collect(keys(filters))
	for i in 1:length(filters)
		name = names[i]
		x = filters[name]
		
		x = x[r_ell_max .> x .> 0]
		if length(x) < 2
			println("skipping with < 2 members: ", "\t", name)
			continue
		end
		
		println(name, "\t", length(x))
		bins = LilGuys.Interface.bins_both(log10.(x), nothing, bin_width=0.05)

		if length(bins) < 2
			println("skipping with < 2 bins: ", "\t", name)
			continue
		end
		
		prof = LilGuys.StellarProfile(x, normalization=:none, bins=(bins))

		if sequential
			 kwargs = (; color=i, colorrange=(1, length(names)))
			else
			 kwargs = (; )
		end

		
		scatterlines!(prof.log_r, prof.log_Sigma, label=name; kwargs...)
	end

	axislegend(ax, axislabel, position=:lb)

	return fig
end

# ╔═╡ d1fc80e9-55ca-4bc7-b458-7925a9c87d15
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,

	)


	
	x = (best_stars.r_ell)
	x = x[r_ell_max .> x .> 0]
	prof = LilGuys.StellarProfile(x, normalization=:none, bins=get_bins(x, num_per_bin=η_bins_bg*get_num_per_bin(log10.(x))))
	
	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err, label="all", )

	y = 10 .^ prof.log_Sigma
	y_err = prof.log_Sigma_err .* y .* log(10)

	y_bg = 10 .^ log_Σ_bg
	y_bg_err = y_bg * log_Σ_bg_err * log(10)
	y_new = y .- y_bg

	filt = 1:(findfirst(y_new .<= 0)-1)
	if sum(filt) > 1
		@info filt
		@info y_new
	
		y_new = y_new[filt]
		y_err = y_err[filt]
	
		y_new_em = min.(log10.(y_new) .- log10.(max.(0, y_new .- y_err .- y_bg_err)), 10.0)
		y_new_ep = log10.(y_new .+ y_err .+ y_bg_err) .- log10.(y_new)
		y_err = collect(zip(y_new_em, y_new_ep))[filt]
	
		@info y_err
		
		errscatter!(prof.log_r[filt], log10.(y_new), yerr=y_err, label="bg-subtracted")

	end

	hlines!(log_Σ_bg, color=:black, label="background")
	hspan!(log_Σ_bg - log_Σ_bg_err, log_Σ_bg + log_Σ_bg_err, color=(:black, 0.1), )

	filt = best_stars.PSAT .> 0.2
	x = (best_stars.r_ell[filt])
	x = x[r_ell_max .> x .> 0]
	prof = LilGuys.StellarProfile(x, normalization=:none, bins=get_bins(x))
	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err, label="PSAT > 0.2")

	Legend(fig[1,2], ax)

	ylims!(ax, minimum(prof.log_Sigma .- 2prof.log_Sigma_err), nothing)
	@savefig "density_profiles_bg_subtraction"

	fig
end

# ╔═╡ 79ddbc86-69f7-4543-a535-f5563685f301
let 
	fig = plot_density_for_filter(OrderedDict(
		"nonmembers" => best_stars.r_ell[best_stars.PSAT .< 0.1],
		"members" => best_stars.r_ell[best_stars.PSAT .> 0.2],
		"cmd" => best_stars.r_ell[best_stars.PSAT_CMD .> 0.2],
		"space" => best_stars.r_ell[best_stars.PSAT_S .> 0.2],
		"PM" => best_stars.r_ell[best_stars.PSAT_PM .> 0.2],
		"PM+CMD" => best_stars.r_ell[best_stars.PSAT_NOSPACE .> 0.2],
		
	),
)
	@savefig "density_profiles_by_component"

	fig
end
	

# ╔═╡ 0ca3a04c-8dae-4672-8aab-e105ba127ec7
md"""
How much does PSAT affect the resulting density profile?
"""

# ╔═╡ 0fc6a483-72d7-4512-9739-d6149dcc40a9
let

	threshholds = [0.01, 0.1, 0.2, 0.5, 0.9, 0.99]

	fig = plot_density_for_filter(
		OrderedDict([string(pmin) => best_stars.r_ell[best_stars.PSAT .>= pmin] for pmin in threshholds]),
		sequential=true,
		axislabel="PSAT"
	)

	@savefig "density_profiles_by_PSAT"
	fig
end

# ╔═╡ 158c7038-d20d-467b-b10e-3e0a41ad26b6
let

	
	threshholds = LilGuys.quantile(members.phot_g_mean_mag, LinRange(0, 1, 5))
	threshholds = round.(threshholds, digits=1)

	Nc = 4
	fig = plot_density_for_filter(
		OrderedDict(
			["$(threshholds[i]) - $(threshholds[i+1])" => members.r_ell[threshholds[i] .<= members.phot_g_mean_mag .< threshholds[i+1]] for i in 1:Nc]),
		sequential=true,
		axislabel="PSAT"
	)

	@savefig "density_profiles_by_gmag"
	fig
end

# ╔═╡ 766dd278-a87f-44d3-878c-18cbe624a57e
let
	filters = OrderedDict(
		"r_ell" => members.r_ell_old .* observed_properties["r_h"],
		"r_circ" => @.(60 * sqrt(members.xi^2 + members.eta^2)),
		"r_ell_me" => 60*lguys.calc_r_ell(members.xi, members.eta, observed_properties["ellipticity"], observed_properties["position_angle"]) ,
		"recentred" => 60*lguys.calc_r_ell(xi, eta, observed_properties["ellipticity"], observed_properties["position_angle"]) 
	)

	fig = plot_density_for_filter(filters)

	@savefig "density_profiles_r_method"

	fig
end

# ╔═╡ 4ac549ef-994c-4087-a42a-723994bc979e
md"""
The figure below samples over the uncertanties in the centring, position angle, and ellipticity.
"""

# ╔═╡ 724b5ca9-0287-4471-bcb4-64fb459557e6
function get_mean_err(d::Dict, key)
	mean = d[key]
	if "$(key)_em" ∈ keys(d)
		err = max(d["$(key)_em"], d["$(key)_ep"])
	else
		err = d[key * "_err"]
	end

	return mean, err
end

# ╔═╡ 259dbfa9-a2de-494b-b3f4-fc69f37478dc
let
	xi, eta = members.xi, members.eta

	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	ell_0, ell_err = get_mean_err(observed_properties, "ellipticity")
	theta_0, theta_err = get_mean_err(observed_properties, "position_angle")
	for _ in 1:1000
		ell =  ell_0 + randn() * ell_err
		θ = theta_0 + randn() * theta_err

		xi_s = xi .+ randn()*cen_err
		eta_s =  eta .+ randn()*cen_err
		x = 60 * lguys.calc_r_ell(xi_s, eta_s, ell, θ) 
		x = x[r_ell_max .> x .> 0]
		
		bins = get_bins((x))
		prof = LilGuys.StellarProfile(x, normalization=:none, bins=bins)
	
		lines!(prof.log_r, prof.log_Sigma, color=(:black, 0.03))
	end

	x = members.r_ell
	bins = LilGuys.Interface.bins_both(log10.(x), nothing, bin_width=0.05)
	prof = LilGuys.StellarProfile(x, normalization=:none, bins=(bins))

	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err, color=COLORS[2])

	LilGuys.hide_grid!(ax)
	
	@savefig "density_profiles_mc"

	fig
end

# ╔═╡ d2668f44-9cb7-4ab8-af49-1aef2baecdda
md"""
# Sanity checks
"""

# ╔═╡ 12bf64e9-ae49-4f93-8572-660dbc86ffc1
scatter(all_stars.r_ell, all_stars.r_ell_old, 
	axis = (;
	xlabel = "r ell",
	ylabel = "r ell jax",
	)
)

# ╔═╡ 2b6f8f68-61ff-4052-901d-8f839044d969
md"""
The plot below sanity checks the calculation of r_ell_max. This should be the largest elliptical radius bin with complete star counts. The orange region should be a complete circle perfectly inscribed in the oblong, blue region.
"""

# ╔═╡ e34d7c3f-5d98-4b35-8abb-91afc16a2a08
let
	fig = Figure()
	ax = PolarAxis(fig[1,1])
	filt = best_stars.r_ell .< Inf
	df = best_stars[filt, :]
	scatter!(ax, atan.(df.xi, df.eta), df.r_ell, markersize=1)

	
	filt = best_stars.r_ell .< r_ell_max
	df = best_stars[filt, :]
	scatter!(ax, atan.(df.xi, df.eta), df.r_ell, markersize=2)
	fig
end

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═cf7aeb0a-d453-462a-b43d-a832567440fd
# ╠═da9ca7a7-18b8-49cb-a041-ab1c667920ff
# ╠═489f6d21-9e9a-4b1e-b05c-c63a44ba1951
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═db1264b7-02c3-4b55-ae2b-9ce78aa1304a
# ╠═81c677bf-8372-4fd3-9d60-af8377485f1c
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═b05bb3ea-2492-4da0-bfd5-8aaf5bf936c4
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
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
# ╠═2d4da56b-5d0a-49d3-83ae-a90f85192101
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╟─722a1f29-be6b-4fca-a9b3-8c304b6c4eb6
# ╟─b6b8d7e7-d0a1-45bc-8593-66aea4a65637
# ╟─a4e857c0-39d9-4fa0-871f-ccc66cb17c25
# ╟─d27062d3-9374-4342-bb4e-3461c774da79
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
# ╟─b0aec8b4-410f-4678-ad8c-cb055b743bd3
# ╟─3c2ee4bb-d594-47e7-a7f3-5b7c9e38b7d4
# ╟─b191b7b8-28ed-4aee-b376-b2bdf22abc29
# ╟─b06445ed-6290-42e8-bdb8-20f0372518c2
# ╟─2a8e13bd-b8ca-4cf3-b61b-096af218265b
# ╟─25dd2e38-1c78-4ef6-95ef-d38a0e80f0e4
# ╟─1b63c1ee-0a46-49c6-880e-c1742055605e
# ╟─643773f0-2b03-46f7-8df4-579b3c708909
# ╟─cad88ab4-1ada-4cc3-8ef0-d4674500d57e
# ╟─18252c3c-d4d6-46f2-8277-33a700ab4074
# ╟─d47a8d02-6c63-41fd-95b5-85550d1f5a85
# ╟─b7e7e245-3732-453c-a2ec-3fe84c65c012
# ╟─87f29963-2c7f-4d8f-8631-f70197d315c0
# ╟─72aa995b-8dd1-4d08-88c2-5fb527ff7f3e
# ╟─6cdcd105-fa1d-48f5-9861-b722aafa28a3
# ╟─cf8e51a9-194c-4aa2-9766-39256a35ef02
# ╟─e265ec56-d62d-4295-975f-98ba93cf5e16
# ╟─e2686165-6499-4558-9da9-91a78624b992
# ╟─1f2aa1ee-b395-4040-84cb-973d57b534d5
# ╟─ca06c336-3950-4831-9f3f-7fe41406bc08
# ╟─63a3e02b-0ea0-472d-8c2e-962bd8e770c0
# ╟─3a52eba9-e07d-4445-a322-1c4b8bb33024
# ╠═f54c50ae-0e03-4917-899f-81a30f116b3f
# ╠═b4f816f0-0abe-43c2-bf5b-48a2eb2b99c8
# ╟─dc280eae-4f10-4167-9545-c70d67af2d8f
# ╟─fb3a27be-62c7-41ab-ba5d-761efb74c044
# ╟─b4363369-4590-4446-9f4a-9f71ff31a870
# ╟─00d4679d-6c55-4929-bbd3-c796479ba70d
# ╟─6be3bd9d-c257-444f-87ab-2b3956fc83aa
# ╟─4ad0330e-f1d7-4d62-b6f9-f49be1aa0d80
# ╟─9a2ca944-dd5c-43f5-b70d-ab3d67b42743
# ╟─79ac654d-7de9-4a2c-bc5d-2f8c8a7e057f
# ╟─479c9cad-8cae-4d9b-87d3-79a816b3af82
# ╟─a7dbb8bf-8c82-4ace-9f17-e9aa159a84ff
# ╟─ec4a8ca8-a4b9-41d7-a369-8049b2a98795
# ╠═32fd5ed6-fb41-456c-bd05-d28ad013c2c7
# ╠═301f42ce-437f-4d45-a186-42df47303f25
# ╠═625795e7-aca7-4bb9-8db6-cfcfb87bddaf
# ╟─3c05db74-0cc1-4893-a468-b495c0121118
# ╠═e3f951ed-3725-4ece-bf59-cb6998fdf170
# ╠═e4e6239f-9fdd-433b-8a39-ffaeaf96b59d
# ╠═363dfc70-ec52-44ed-ade4-f88322247871
# ╟─c68ff47f-97de-4bbb-8b0a-42d837da21f7
# ╟─8d497b75-a47d-42d1-9e54-e15bfb1db805
# ╟─21f80519-1215-4aa1-be36-92956903502e
# ╠═a4d06da3-0197-4c34-88c9-29099ff9ab46
# ╠═571b98f1-10d0-4012-8bfb-755f26f71b37
# ╟─b1724e25-6318-460e-a057-70b4588b3321
# ╠═d1c46aea-e419-4feb-873f-ba035b46cda7
# ╠═642eb4bd-6082-4209-bc92-c4ef5ddaa3fc
# ╠═5cf14b86-0a30-410c-955a-82e9e8696e7b
# ╠═c6805532-116c-4e8c-9c40-81161a1334cd
# ╠═f9edc1ee-2b05-4750-ba26-c12cd95d4041
# ╠═2d96b315-1207-4278-98e4-5b9239773779
# ╠═a507dc20-8326-4f66-aebe-24ea16d4254e
# ╟─53bc3533-f9fc-4789-9c8e-e5085faacc33
# ╠═76fa8518-bc9a-4126-880b-2945d6b8e853
# ╠═3a676557-390a-4042-ad9f-527c292a9da2
# ╠═d50bc64f-f67b-4498-aea8-81c4b5c258cc
# ╠═8c253d06-b902-4baa-82f7-6c507f8cf014
# ╠═d1fc80e9-55ca-4bc7-b458-7925a9c87d15
# ╠═79ddbc86-69f7-4543-a535-f5563685f301
# ╟─0ca3a04c-8dae-4672-8aab-e105ba127ec7
# ╠═0fc6a483-72d7-4512-9739-d6149dcc40a9
# ╠═158c7038-d20d-467b-b10e-3e0a41ad26b6
# ╟─766dd278-a87f-44d3-878c-18cbe624a57e
# ╟─4ac549ef-994c-4087-a42a-723994bc979e
# ╠═724b5ca9-0287-4471-bcb4-64fb459557e6
# ╠═259dbfa9-a2de-494b-b3f4-fc69f37478dc
# ╟─d2668f44-9cb7-4ab8-af49-1aef2baecdda
# ╠═12bf64e9-ae49-4f93-8572-660dbc86ffc1
# ╟─2b6f8f68-61ff-4052-901d-8f839044d969
# ╟─e34d7c3f-5d98-4b35-8abb-91afc16a2a08
