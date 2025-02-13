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

# ╔═╡ cf7aeb0a-d453-462a-b43d-a832567440fd
diverging_cmap = Reverse(:bluesreds)

# ╔═╡ da9ca7a7-18b8-49cb-a041-ab1c667920ff
import DensityEstimators: histogram2d

# ╔═╡ 489f6d21-9e9a-4b1e-b05c-c63a44ba1951
import StatsBase: percentile, mean

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
galaxyname = "carina"

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using LilGuys
	figdir = joinpath(galaxyname, "figures")
end

# ╔═╡ b05bb3ea-2492-4da0-bfd5-8aaf5bf936c4
if isdir(galaxyname)
	if ! isdir(figdir)
		mkdir(figdir)
	end
end

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	all_stars = lguys.read_fits("$galaxyname/data/j24_1c.fits")

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

# ╔═╡ 223abc41-5158-49c2-96bf-df55b7be1114
cen = lguys.calc_centre2D(all_stars.xi, all_stars.eta, "mean")

# ╔═╡ db1264b7-02c3-4b55-ae2b-9ce78aa1304a
import DensityEstimators: histogram

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = all_stars[all_stars.PSAT .> 0.2, :]

# ╔═╡ 731ea468-5003-44e9-95b8-7fa7ef4b475b
Nmemb = size(members, 1)

# ╔═╡ 3e40a3f8-efbc-4807-b185-22fbb2e99adf
begin 
	xi = members.xi .- cen[1]
	eta = members.eta .- cen[2]
end

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

# ╔═╡ 722a1f29-be6b-4fca-a9b3-8c304b6c4eb6
let 
	dθ = 60θmax
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel=L"\xi / \textrm{arcmin}", ylabel=L"\eta / \textrm{arcmin}", aspect=1, limits=(-dθ, dθ, -dθ, dθ), xgridvisible=false, ygridvisible=false)


	scatter!(60best_stars.xi, 60best_stars.eta, 
		alpha=1, markersize=1, color=(:black, 0.2), label="all")
	
	scatter!(60members_nospace.xi, 60members_nospace.eta, 
		alpha=1, markersize=3, label="members w/o spatial")
	scatter!(60members.xi, 60members.eta, alpha=1, markersize=3,
		label = "members")

	Legend(fig[1, 2], ax)

	@savefig "xi_eta_members.jax_observations"
	fig
end

# ╔═╡ 4a1bf5cb-6083-429f-8747-9ed0478f8f3e
md"""
If there is crowding, the plot below should appear to have deviations in brightness distributions. Otherwise, the distribution should look like a scatterplot
"""

# ╔═╡ 1942cac0-2da8-4bec-80d0-a9302ddc9ea1
let 
	global Σ_kde,kde_d, x_kde, y_kde

	kde_d = kde((60xi, 60eta))
	ik = InterpKDE(kde_d)

	Σ_kde = @. pdf([ik], kde_d.x', kde_d.y)
	Σ_kde = Σ_kde'
	x_kde = kde_d.x
	y_kde = kde_d.y
end

# ╔═╡ d49c5f63-e2a4-4ff5-af4a-d4239b50abae
begin
	kde_y = reshape(kde_d.y, (1, :))
	kde_r = @. sqrt((kde_d.x')^2 + (kde_d.y)^2)
	kde_theta = atan.(kde_d.y, kde_d.x')
end

# ╔═╡ cdd246ae-e6bb-49a9-9512-085bd37827bd
let 
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel = "xi / arcmin",
		ylabel = "eta / arcmin",
		aspect = DataAspect()
	)

	density_scale = 1e-3
	p = heatmap!(x_kde, y_kde, Σ_kde, colormap=:greys, colorscale=x -> asinh(x/density_scale))
	ax.aspect = 1
	contour!(x_kde, y_kde, asinh.(Σ_kde ./ density_scale), levels=10 )


	filt = members.r_ell .> 0
	scatter!(ax, 60xi[filt], 60eta[filt], markersize=3, color=Arya.COLORS[3])

	Colorbar(fig[1, 2], p, #ticks=[10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0],
		label = "density / $density_scale"
	)
	fig
end

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

# ╔═╡ 4f538de6-8827-454f-a97f-8f6c2cd7ea3f
scatter(members.ra, members.dec, 
	axis=(; xlabel="ra / degrees", ylabel="dec / degrees", title="members only")
)

# ╔═╡ 83b506bb-f385-464e-8f5c-7adfff45105a
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
	title = "members")

	
	h = scatter!(members.bp_rp, members.phot_g_mean_mag, color=log10.(members.r_ell), 
	)

	Colorbar(f[1, 2], h, label="log10 ell radius / rh")
	f
end

# ╔═╡ a4e857c0-39d9-4fa0-871f-ccc66cb17c25
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\varpi\ /\ \textrm{mas}", 
		ylabel=L"\delta\varpi\ /\ \textrm{mas}",
		limits=(-4, 5, 0, 3)
	)

	scatter!(all_stars.parallax, all_stars.parallax_error, color=:black, markersize=1, label = "all stars")

	scatter!(best_stars.parallax, best_stars.parallax_error, markersize=3, label="best")

	scatter!(members.parallax, members.parallax_error, label="members")


	LilGuys.hide_grid!(ax)
	axislegend()

	@savefig "pi_pierr_memb"
	fig
end

# ╔═╡ c70fce56-9367-4f6a-92bb-f0d0d7a616e0
function calc_Σ(rs)
	#r_bins = lguys.make_equal_number_bins(rs, 500)
	r_bins = 10 .^ LinRange(log10(minimum(rs)), log10(maximum(rs)), 50)
	
	N_bins = length(r_bins) - 1
	ν = zeros(N_bins)
	r_mid = zeros(N_bins)
	for i in 1:(length(r_bins) - 1)
		fl = r_bins[i] .< rs .< r_bins[i+1]
		ν[i] = sum(fl)
		r_mid[i] = lguys.mean(rs[fl])
	end
	
	Areas = π  * diff(r_bins.^2 )
	Σs = ν ./ Areas
	return r_mid, Σs
end

# ╔═╡ d27062d3-9374-4342-bb4e-3461c774da79
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		limits=(-10, 10, -10, 10)
		)
	p = scatter!((all_stars.pmra), (all_stars.pmdec), 
		markersize=1, color=(:black, 0.3),
		label="all"
	)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), 
		markersize=2,
		label="best"
	)

	p = scatter!((members.pmra), (members.pmdec), 
		markersize=5,
		label="members"
	)

	Legend(fig[1,2], ax)
	
	LilGuys.hide_grid!(ax)

	@savefig "pm_members.jax_observations"

	fig
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

# ╔═╡ d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
histogram(vec(kde_theta), θ_bins , weights=vec(Σ_kde))

# ╔═╡ 9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
let
	dθ = 1

	f = Figure()
	ax = Axis(f[1, 1], limits=(-dθ,dθ,-dθ,dθ), aspect=1,
		xlabel=L"\xi",
		ylabel=L"\eta")
	
	LilGuys.hide_grid!(ax)
	
	for a in LinRange(0.05dθ, dθ, 8)
		t = LinRange(0, 360, 100)
		x = @. a*cosd(t)
		y = @. a*sind(t)
		lines!(x, y, label="", color="black")
	end
	
	scatter!(members.xi, members.eta, markersize=3, label="", alpha=1)
	
	f
end

# ╔═╡ ec1c1408-5c3c-4bfa-b8bd-f32dbc3379a9
ϕ = atan.(eta, xi)

# ╔═╡ 9282a867-b9c7-4a50-b58a-da3ceb068cad
radii = @. sqrt(xi^2 + eta^2)

# ╔═╡ 6bb8fe50-d7c0-4d76-9d05-a403b7e3239f
rh = percentile(radii, 50)

# ╔═╡ 6ddd0b7b-6c51-442e-97c4-ddc321899da1
60rh

# ╔═╡ 4ff675ef-8609-49df-bfac-4070f14e3c25
maximum(radii)

# ╔═╡ 52d3db7f-8bbf-43cf-89f2-16134247d713
let
	fig = Figure()
	ax = PolarAxis(fig[1,1], rlimits=(0, 0.5))

	scatter!(ax, ϕ, radii, markersize=2)
	fig
end

# ╔═╡ 9427aa86-2d46-4890-b2c5-36de3ca83264
r_cuts = [percentile(radii, p) for p in LinRange(0, 100, 8)]

# ╔═╡ 75b8cbb1-519e-4174-a7b9-22e661cbde5c
N_cuts = length(r_cuts) - 1

# ╔═╡ ebb24519-ee2d-4781-a2bd-16ccbc957060
let
	f = Figure(size=(1200, 500))
	ax = Axis(f[1, 1], xlabel=L"\theta / \textrm{radians}", ylabel="normalized density")
	bins = 10
	ps = []
	labels = []
	for i in 1:N_cuts
		filt = r_cuts[i] .<= radii .<= r_cuts[i+1] 
		println(sum(filt))
		if sum(filt) > 0
			h = histogram(ϕ[filt], bins)
			bins = h.bins
			counts = h.values
			counts .*= 1/sum(diff(bins) .* counts) * sum(diff(bins))
 			p = lines!(ax, lguys.midpoints(bins), counts, color=i, colorrange=(1, N_cuts))
			push!(ps, p)

			r_l = round(r_cuts[i], digits=3)
			r_h = round(r_cuts[i+1], digits=3)
			push!(labels, L"r \in [%$r_l, %$r_h ]")
		end
	end

	Legend(f[1,2], ps, labels)
	f
end

# ╔═╡ 7af6fd4f-a7f4-4e7f-a1e3-8899717dd2cd
let
	f = Figure(size=(1200, 500))
	ax = Axis(f[1, 1], xlabel=L"\theta / \textrm{radians}", ylabel="normalized density for rell")
	bins = 10
	ps = []
	labels = []
	for i in 1:N_cuts
		filt = r_cuts[i] .<= members.r_ell * rh .<= r_cuts[i+1] 
		println(sum(filt))
		if sum(filt) > 0
			h = histogram(ϕ[filt], bins)
			h.values .*= 1/sum(diff(h.bins) .* h.values) * sum(diff(h.bins))
 			p = lines!(ax, lguys.midpoints(h.bins), h.values, color=i, colorrange=(1, N_cuts))
			push!(ps, p)

			r_l = round(r_cuts[i], digits=3)
			r_h = round(r_cuts[i+1], digits=3)
			push!(labels, L"rell \in [%$r_l, %$r_h ]")
		end
	end

	Legend(f[1,2], ps, labels)
	f
end

# ╔═╡ 49a5de24-907c-49e1-9747-8be5bb28b64d
let 
	fig, ax, h = hist(ϕ, bins=50, label="")
	ax.xlabel = "angle"
	ax.ylabel = "count"
	fig
end

# ╔═╡ b0aec8b4-410f-4678-ad8c-cb055b743bd3
md"""
# Likelihood components
"""

# ╔═╡ 23e84597-4b11-461b-825d-5203c7b31273
import Statistics: cor

# ╔═╡ 0550a80f-afa1-424a-bfe8-46c426554ab9
extrema(best_stars.LL[.!isnan.(best_stars.LL)])

# ╔═╡ af305fde-e827-4ea6-a102-0e71a112d565
cor(best_stars.LL_S[isfinite.(best_stars.LL)], best_stars.LL[isfinite.(best_stars.LL)])

# ╔═╡ 3c2ee4bb-d594-47e7-a7f3-5b7c9e38b7d4
md"""
## Likelihoods maps
"""

# ╔═╡ b06445ed-6290-42e8-bdb8-20f0372518c2
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.xi), (best_stars.eta), color=(best_stars.L_S_BKD),
		colormap=:bluesreds,
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="L background", )
	
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
		colormap=diverging_cmap,
		markersize=3,
		colorrange=(-1, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L s", )
	
	fig
end

# ╔═╡ 2004a768-d1cf-4641-869c-39796e0f96aa
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.xi), (best_stars.eta), color=log10.(best_stars.L_S_SAT) .- (best_stars.LL),
		colormap=diverging_cmap,
		markersize=3,
		#colorrange=(-1, 1)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L s / L tot", )
	
	fig
end

# ╔═╡ 1f2aa1ee-b395-4040-84cb-973d57b534d5
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log likilihood",
		ylabel = "counts"
	)
	bins = LinRange(-10, 2, 30)
	
	stephist!(log10.(best_stars.L_S_SAT), 
		bins=bins, label="spatial")
	
	stephist!(log10.(best_stars.L_CMD_SAT), bins=bins, label="cmd")
	
	stephist!(log10.(best_stars.L_PM_SAT), bins=bins,  label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ be659eae-528c-4d84-9a18-2a9da57d8ce3
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
	title = "members")

	
	h = scatter!(members.bp_rp, members.phot_g_mean_mag, color=log10.(members.L_CMD_BKD), 
	)

	Colorbar(f[1, 2], h, label="log10 ell radius / rh")
	f
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

	Colorbar(f[1, 2], h, label="log L sat")
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
	Colorbar(f[1, 2], h, label="log L mw")
	f
end

# ╔═╡ 8c048f55-85f0-4c30-b05a-fffeb3856ede
fsat = 0.001

# ╔═╡ 87f29963-2c7f-4d8f-8631-f70197d315c0
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	filt = best_stars.r_ell .< 3
	df = best_stars[filt, :]

	
	p = scatter!((df.xi), (df.eta), color=log10.(df.L_S_SAT) .- log10.(df.L_SAT*fsat .+ (1-fsat)*df.L_BKD),
		colormap=diverging_cmap,
		markersize=3,
		#colorrange=(-1, 1)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L s", )
	
	fig
end

# ╔═╡ 062c6e95-bde2-4489-9be7-1f012408ae32
Ltot = sum(fsat * best_stars.L_SAT .+ best_stars.L_BKD)

# ╔═╡ 359924c4-36d0-4439-8551-5dfc6edf1396
maximum(best_stars.L_S_SAT ./ Ltot)

# ╔═╡ 89850963-a799-4bda-a801-9da181c50d34
sum(fsat * all_stars.L_SAT)

# ╔═╡ 8eb7f0a7-fcf5-46c9-8821-bbce47bb0816
sum(best_stars.L_BKD)

# ╔═╡ c8e13c62-2202-4256-a6f3-9d390513090e
sum(all_stars.L_BKD)

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

	Colorbar(f[1, 2], h, label="log L sat / bg")
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
		colormap=diverging_cmap,
		colorrange=(-30, 1),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L sat", )
	
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
		colormap=diverging_cmap,
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L pm bg", )
	
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

	Colorbar(fig[1,2], p, label="log L pm sat / bg", )
	
	fig
end

# ╔═╡ 3a52eba9-e07d-4445-a322-1c4b8bb33024
md"""
### Correlations
"""

# ╔═╡ f54c50ae-0e03-4917-899f-81a30f116b3f
function plot_corr(best_stars)
	fig = Figure(
		size = (800, 400)
	)

	kwargs = (;
		markersize=3,
	)


	axs = []
	for i in 1:3
		col = [:LL_S, :LL_CMD, :LL_PM][i]
		label = ["space", "cmd", "pm"][i]

		ax = Axis(fig[1,i],
			xlabel = L"\log\mathcal{L}_\textrm{%$label}",
			ylabel = L"\log\mathcal{L}",
			aspect=DataAspect()
		)
		push!(axs, ax)
		
		p = scatter!(best_stars[:, col], (best_stars.LL);
			color=(best_stars.PSAT),
			kwargs...
		)
		LilGuys.hide_grid!(ax)
	
		R = cor(best_stars[:, col], best_stars.LL)
		t = "ρ = $(round(R, digits=2))"
		
		text!(0.1, 0.9, text=t, space=:relative)

		if i > 1
			hideydecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	linkaxes!(axs...)
	fig
end

# ╔═╡ b4f816f0-0abe-43c2-bf5b-48a2eb2b99c8
plot_corr(members)

# ╔═╡ fb3a27be-62c7-41ab-ba5d-761efb74c044
let
	fig = plot_corr(best_stars[best_stars.LL .> -10, :])

	@savefig "LL_corr_trunc.jax_observations"

	fig
end

# ╔═╡ b4363369-4590-4446-9f4a-9f71ff31a870
plot_corr(best_stars[best_stars.LL .> -20, :])

# ╔═╡ 00d4679d-6c55-4929-bbd3-c796479ba70d
plot_corr(all_stars[all_stars.LL .> -80, :])

# ╔═╡ 6be3bd9d-c257-444f-87ab-2b3956fc83aa
md"""
## Fractional contributions
"""

# ╔═╡ 9a2ca944-dd5c-43f5-b70d-ab3d67b42743
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "LL",
		ylabel = "f spatial",
		limits=(-30, 10, -1.5, 3)
	)
	
	scatter!(best_stars.LL, (best_stars.f_LL_S), markersize=3)
	scatter!(members.LL, (members.f_LL_S))

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

# ╔═╡ f952aaac-0fd1-46ae-85fc-02dc294281ff


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

# ╔═╡ 5ec54a1e-fab8-44fc-8fd4-a72d15396c2b
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		limits = (-2, 2, -2, 2)
		)

	filt = best_stars.PSAT .> 0.2
	filt .&= best_stars.f_LL_S_norm .> 0.5
	
	
	p = scatter!((best_stars.xi[filt]), (best_stars.eta[filt]),
		markersize=3,
		color = best_stars.r_ell[filt]
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="rh", )
	
	fig
end

# ╔═╡ 625795e7-aca7-4bb9-8db6-cfcfb87bddaf
import DensityEstimators as DE

# ╔═╡ e4e6239f-9fdd-433b-8a39-ffaeaf96b59d
let
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

		h = DE.histogram(df[:, col])
		
		x = midpoints(h.bins)
		scatter!(x, h.values, color=x,
			colorrange=(-1, 1),
			colormap=:bluesreds,
		)

		
		ax = Axis(fig[2,i],
			xlabel = "ξ / degree",
			ylabel = "η / degree",
			aspect=DataAspect(),
			limits = (-2, 2, -2, 2)
			)
		
		p = scatter!(df.xi, df.eta,
			markersize=3,
			color = df[:, col],
			colorrange=(-1, 1),
			colormap=:bluesreds,
		)
		LilGuys.hide_grid!(ax)
	end



	rowsize!(fig.layout, 1, Relative(0.2))

	@savefig "xi_eta_vs_f.jax_observations"
	fig
end

# ╔═╡ 363dfc70-ec52-44ed-ade4-f88322247871
let
	fig = Figure(
		size=(800, 500)
	)

	filt = best_stars.LL .> -10
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

		h = DE.histogram(df[:, col])
		
		x = midpoints(h.bins)
		scatter!(x, h.values, color=x,
			colorrange=(-1, 1),
			colormap=:bluesreds,
		)

		
		ax = Axis(fig[2,i],
			xlabel = "ξ / degree",
			ylabel = "η / degree",
			aspect=DataAspect(),
			limits = (-2, 2, -2, 2)
			)
		
		p = scatter!(df.xi, df.eta,
			markersize=3,
			color = df[:, col],
			colorrange=(-1, 1),
			colormap=:bluesreds,
		)
		LilGuys.hide_grid!(ax)
	end



	rowsize!(fig.layout, 1, Relative(0.2))

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
	
	Colorbar(fig[1,2], p, label="LL", )
	fig
end

# ╔═╡ 21f80519-1215-4aa1-be36-92956903502e
md"""
# Density profiles
"""

# ╔═╡ 5be4a3eb-3f7a-4228-9126-8bcfe45d85e2
let
	filters = OrderedDict(
		"all" => trues(length(best_stars.PSAT)),
		"nonmembers" => best_stars.PSAT .< 0.1,
		"members" => best_stars.PSAT .> 0.2,
		"PSAT cmd" => best_stars.PSAT_CMD .> 0.2,
		"PSAT S" => best_stars.PSAT_S .> 0.2,
		"PSAT PM" => best_stars.PSAT_PM .> 0.2,
		"PSAT PM+CMD" => best_stars.PSAT_NOSPACE .> 0.2,
		
	)

	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = L"\log\ r_\textrm{ell}\,/\,r_h",
		ylabel = L"$\Sigma_\star$\,/\,stars\,$r_h^{-2}$",
	)

	for (name, filt) in filters
		x = (best_stars.r_ell[filt])
		x = x[x .> 0]
		if length(x) < 2
			println("skipping", "\t", name)
			continue
		end
		println(name, "\t", length(x))
		bins = LilGuys.Interface.bins_both(x, nothing, bin_width=0.05)
		
		prof = LilGuys.StellarProfile(x, normalization=:none, bins=log10.(bins))
	
		scatterlines!(prof.log_r, prof.log_Sigma, label=name)
	end

	axislegend(ax, position=:lb)

	@savefig "density_profiles_by_component.jax_observations"
	fig
end

# ╔═╡ f0ed5063-83f3-4268-88a6-23ebd4d16a4d
let

	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = L"\log\ r_\textrm{ell}\,/\,r_h",
		ylabel = L"$\Sigma_\star$\,/\,stars\,$r_h^{-2}$",
	)

	threshholds = [0.01, 0.1, 0.2, 0.5, 0.9, 0.99]
	for (i, threshhold) in enumerate(threshholds)
		x = (best_stars.r_ell[best_stars.PSAT .>= threshhold])
		name = string(threshhold)
		x = x[x .> 0]
		if length(x) < 2
			println("skipping", "\t", name)
			continue
		end
		
		println(name, "\t", length(x))
		bins = LilGuys.Interface.bins_both(x, nothing, bin_width=0.05)
		
		prof = LilGuys.StellarProfile(x, normalization=:none, bins=log10.(bins))
	
		scatterlines!(prof.log_r, prof.log_Sigma, label=name, color=i, colorrange=(1, length(threshholds)))
	end

	axislegend(ax, position=:lb, title="PSAT")

	@savefig "density_profiles_by_PSAT.jax_observations"
	fig
end

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═cf7aeb0a-d453-462a-b43d-a832567440fd
# ╠═da9ca7a7-18b8-49cb-a041-ab1c667920ff
# ╠═489f6d21-9e9a-4b1e-b05c-c63a44ba1951
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═b05bb3ea-2492-4da0-bfd5-8aaf5bf936c4
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═223abc41-5158-49c2-96bf-df55b7be1114
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═6bb8fe50-d7c0-4d76-9d05-a403b7e3239f
# ╠═6ddd0b7b-6c51-442e-97c4-ddc321899da1
# ╠═db1264b7-02c3-4b55-ae2b-9ce78aa1304a
# ╠═3e40a3f8-efbc-4807-b185-22fbb2e99adf
# ╠═4ff675ef-8609-49df-bfac-4070f14e3c25
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═60d0e593-88fd-4b4c-9009-cc24a597c6d5
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╟─a9d94121-ea6e-416a-bae8-aa93c16bde72
# ╠═f3d3f7b2-30c2-4f89-bef5-ac1a8720963f
# ╠═2d4da56b-5d0a-49d3-83ae-a90f85192101
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─722a1f29-be6b-4fca-a9b3-8c304b6c4eb6
# ╟─4a1bf5cb-6083-429f-8747-9ed0478f8f3e
# ╠═1942cac0-2da8-4bec-80d0-a9302ddc9ea1
# ╠═d49c5f63-e2a4-4ff5-af4a-d4239b50abae
# ╠═cdd246ae-e6bb-49a9-9512-085bd37827bd
# ╟─02d19c07-411f-40d6-9ac3-010ebfd4bdfe
# ╠═4f538de6-8827-454f-a97f-8f6c2cd7ea3f
# ╠═d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
# ╠═83b506bb-f385-464e-8f5c-7adfff45105a
# ╠═a4e857c0-39d9-4fa0-871f-ccc66cb17c25
# ╠═c70fce56-9367-4f6a-92bb-f0d0d7a616e0
# ╠═d27062d3-9374-4342-bb4e-3461c774da79
# ╟─4056fc57-0058-4da5-89af-e372164b2c64
# ╠═c6115ac4-57de-43f1-abfe-adfde1294bb5
# ╠═042f8215-b6fe-49db-a166-055567da4bb0
# ╠═9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
# ╠═ec1c1408-5c3c-4bfa-b8bd-f32dbc3379a9
# ╠═9282a867-b9c7-4a50-b58a-da3ceb068cad
# ╠═52d3db7f-8bbf-43cf-89f2-16134247d713
# ╠═75b8cbb1-519e-4174-a7b9-22e661cbde5c
# ╠═9427aa86-2d46-4890-b2c5-36de3ca83264
# ╠═ebb24519-ee2d-4781-a2bd-16ccbc957060
# ╠═7af6fd4f-a7f4-4e7f-a1e3-8899717dd2cd
# ╠═49a5de24-907c-49e1-9747-8be5bb28b64d
# ╟─b0aec8b4-410f-4678-ad8c-cb055b743bd3
# ╠═23e84597-4b11-461b-825d-5203c7b31273
# ╠═0550a80f-afa1-424a-bfe8-46c426554ab9
# ╠═af305fde-e827-4ea6-a102-0e71a112d565
# ╠═3c2ee4bb-d594-47e7-a7f3-5b7c9e38b7d4
# ╟─b06445ed-6290-42e8-bdb8-20f0372518c2
# ╠═2a8e13bd-b8ca-4cf3-b61b-096af218265b
# ╠═2004a768-d1cf-4641-869c-39796e0f96aa
# ╠═87f29963-2c7f-4d8f-8631-f70197d315c0
# ╠═1f2aa1ee-b395-4040-84cb-973d57b534d5
# ╟─be659eae-528c-4d84-9a18-2a9da57d8ce3
# ╟─25dd2e38-1c78-4ef6-95ef-d38a0e80f0e4
# ╟─1b63c1ee-0a46-49c6-880e-c1742055605e
# ╠═8c048f55-85f0-4c30-b05a-fffeb3856ede
# ╠═062c6e95-bde2-4489-9be7-1f012408ae32
# ╠═359924c4-36d0-4439-8551-5dfc6edf1396
# ╠═89850963-a799-4bda-a801-9da181c50d34
# ╠═8eb7f0a7-fcf5-46c9-8821-bbce47bb0816
# ╠═c8e13c62-2202-4256-a6f3-9d390513090e
# ╠═643773f0-2b03-46f7-8df4-579b3c708909
# ╠═cad88ab4-1ada-4cc3-8ef0-d4674500d57e
# ╠═18252c3c-d4d6-46f2-8277-33a700ab4074
# ╠═d47a8d02-6c63-41fd-95b5-85550d1f5a85
# ╟─3a52eba9-e07d-4445-a322-1c4b8bb33024
# ╠═f54c50ae-0e03-4917-899f-81a30f116b3f
# ╠═b4f816f0-0abe-43c2-bf5b-48a2eb2b99c8
# ╟─fb3a27be-62c7-41ab-ba5d-761efb74c044
# ╠═b4363369-4590-4446-9f4a-9f71ff31a870
# ╠═00d4679d-6c55-4929-bbd3-c796479ba70d
# ╟─6be3bd9d-c257-444f-87ab-2b3956fc83aa
# ╠═9a2ca944-dd5c-43f5-b70d-ab3d67b42743
# ╠═79ac654d-7de9-4a2c-bc5d-2f8c8a7e057f
# ╠═a7dbb8bf-8c82-4ace-9f17-e9aa159a84ff
# ╠═32fd5ed6-fb41-456c-bd05-d28ad013c2c7
# ╠═f952aaac-0fd1-46ae-85fc-02dc294281ff
# ╠═8d497b75-a47d-42d1-9e54-e15bfb1db805
# ╠═5ec54a1e-fab8-44fc-8fd4-a72d15396c2b
# ╠═625795e7-aca7-4bb9-8db6-cfcfb87bddaf
# ╟─e4e6239f-9fdd-433b-8a39-ffaeaf96b59d
# ╟─363dfc70-ec52-44ed-ade4-f88322247871
# ╠═c68ff47f-97de-4bbb-8b0a-42d837da21f7
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╟─21f80519-1215-4aa1-be36-92956903502e
# ╠═5be4a3eb-3f7a-4228-9126-8bcfe45d85e2
# ╠═f0ed5063-83f3-4268-88a6-23ebd4d16a4d
