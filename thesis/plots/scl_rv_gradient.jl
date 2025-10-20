### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ 314da249-8a3f-47e9-ac40-f86406ec9955
using CSV, DataFrames

# ╔═╡ 71335ca6-b58b-4c68-8467-81cd4f481b1f
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 0adc833e-8de6-4ab7-a137-9cb6e5aeba86
include("utils.jl")

# ╔═╡ 7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
import StatsBase: median

# ╔═╡ 19038700-b427-4810-8818-9f789e4aeba1
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ f27c2500-b92c-4de0-8482-4a8acfa2667f
module RVUtils
	include(joinpath(ENV["DWARFS_ROOT"], "observations/rv_utils.jl"))
end

# ╔═╡ 2f067cc4-9793-44e0-885e-2875d36e2a20
filename = "rv_combined_x_wide_2c_psat_0.2"

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
df_gradient = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_both.$filename.csv", DataFrame)

# ╔═╡ 5927a386-fdb1-459f-b042-0f95677e345e
rv_memb = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/$filename.fits")

# ╔═╡ 33ee935f-36f6-424c-9f26-138c871c812a
minimum(rv_memb.PSAT_RV) 

# ╔═╡ 668bea30-c7e4-4f90-a0e8-fa315dc79cd1
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ fc4794cf-e45b-41ee-acac-6d5dcf8c02ab
function rolling_medians(x, y; window=50)
	window -= 1
	filt = sortperm(x)
	xs = x[filt]
	ys = y[filt]

    n = length(x)
    medians = Vector{Float64}(undef, n-window)
	x_medians = Vector{Float64}(undef, n-window)
	n_window =  Vector{Int64}(undef, n-window)

	for i in 1:n-window
		start = i
		fin = start + window 
		medians[i] = median(ys[start:fin])
		x_medians[i] = median(xs[start:fin])
		n_window[i] = length(start:fin)
	end

	@assert all(n_window .== window + 1)
	return x_medians, medians
end

# ╔═╡ cfc844b8-0274-4c38-a676-32dfd9bcb44b
filt = (rv_memb.PSAT_RV .> 0.2)

# ╔═╡ f382453c-51ab-4cb7-84a7-ee20354cce93
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "sculptor", "observed_properties.toml"))

# ╔═╡ 4c05511c-d209-4213-a380-d16c60a834a4
xi_p, eta_p = LilGuys.to_orbit_coords(rv_memb.ra, rv_memb.dec, obs_props["ra"], obs_props["dec"], θ_m)

# ╔═╡ e0764d78-afcf-45fb-b07c-ce0485250f29
σv = obs_props["sigma_v"]

# ╔═╡ c1ad6eab-eb98-470f-ab26-b30316e4b8ce
gsr = LilGuys.transform(GSR, ICRS(obs_props))

# ╔═╡ 29a63484-1998-4ce4-978b-ff26cb264ab5
v0 = gsr.radial_velocity

# ╔═╡ d413132c-d6f3-45d4-8a1e-1306b8ad084b
let
	fig = Figure(size=(5*72, 5*72))

	gs = GridLayout(fig[1,1])
	ax = Axis(gs[1,1],
		xlabel = L"\xi\ /\ \textrm{arcmin}",
		ylabel = L"\eta\ /\ \textrm{arcmin}",
		xreversed=true,
		aspect = DataAspect(), 
			  limits=(-50, 50, -30, 30),
		xticks = -40:20:40,
		yticks = -40:20:40,
	)

	arrowscale = 4
	
	θ_grad = median(df_gradient.Θ_grad)
	r_grad = median(df_gradient.r_grad)
	derived_gradient = arrowscale * r_grad * [sind(θ_grad), cosd(θ_grad)]
	
	pm_gradient = arrowscale * LilGuys.pm2kms.([gsr.pmra, gsr.pmdec], gsr.distance) ./ (180/π)
	θ_pm = atand(pm_gradient[1], pm_gradient[2])
	
	memb_stars = rv_memb

	dv = memb_stars.vz .- v0
	x = memb_stars.xi
	y = memb_stars.eta
	w = 1 ./ memb_stars.vz_err .^ 2

	colorrange=(-3σv, 3σv)


	p = scatter!(x, y,
		color = dv,
		colormap=:bluesreds,
		colorrange=colorrange,
		markersize=2.
	)

	arrows2d!([0], [0], [pm_gradient[1]], [pm_gradient[2]])
	text!([pm_gradient[1]], [pm_gradient[2]], 
		  text="  PM", color=:black, fontsize=10, rotation=π/2+deg2rad(θ_pm), align=(:left, :center),
		 )

	
	arrows2d!([0], [0], [derived_gradient[1]], [derived_gradient[2]], 
			label="gradient", color=COLORS[3])
	text!([derived_gradient[1]], [derived_gradient[2]], 
		  text=L"  rot ($\xi'$)", color=COLORS[3], fontsize=10, rotation=π/2+deg2rad(θ_grad), align=(:left, :center))

	
	ellipse!(1obs_props["r_h"], obs_props["ellipticity"], obs_props["position_angle"], linewidth=1, color=:black)
	
	text!(-obs_props["r_h"], 0, text=L"  $R_h$", align=(:left, :centre))

	Colorbar(gs[1, 2], p, label=L"$(\textrm{v}_\textrm{gsr}' - \bar{\textrm{v}}_\textrm{gsr}')$ / km\,s$^{-1}$", halign=:left)


	gs_lower = GridLayout(fig[2, 1])


	ax = Axis(gs_lower[1, 1],
		xlabel = L"$\xi'$ / arcmin",
		ylabel = L"$\textrm{v}_\textrm{gsr}'$ / km s$^{-1}$",
		limits=(-48, 48, 35, 115)
	)
	scatter!(xi_p[filt] .* 60, rv_memb.vz[filt], markersize=2, label="members", color=:black)

	
	for i in 1:400:size(df_gradient, 1)
		M = df_gradient.B[i] .* sind(df_gradient.Θ_grad[i]) .+ df_gradient.A[i] .* cosd(df_gradient.Θ_grad[i])
		x=LinRange(-24, 24, 100)
		y = M*x  .+ df_gradient.μ[i]
		lines!(x, y, alpha=0.03, color=COLORS[1], linewidth=1, label="MCMC samples" => (; alpha=1.0), linestyle=:solid)
	end

	hlines!(v0, color=:black, linewidth=1, alpha=0.3, label="systemic velocity")

	xs, ys = rolling_medians(xi_p[filt] .* 60 , rv_memb.vz[filt], window=200)
	lines!(xs, ys,  color=COLORS[2], label="rolling median")

	annotation!(0, 36, obs_props["R_h"], 35, text=L"R_h")
	axislegend(unique=true, merge=true, tellwidth=false, position=:lt, backgroundcolor=(:white, 0.5))

	colsize!(gs_lower, 1, Relative(3/3))

	@savefig "scl_rv_scatter_gradient"
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═314da249-8a3f-47e9-ac40-f86406ec9955
# ╠═71335ca6-b58b-4c68-8467-81cd4f481b1f
# ╠═7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
# ╠═19038700-b427-4810-8818-9f789e4aeba1
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f27c2500-b92c-4de0-8482-4a8acfa2667f
# ╠═2f067cc4-9793-44e0-885e-2875d36e2a20
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═5927a386-fdb1-459f-b042-0f95677e345e
# ╠═33ee935f-36f6-424c-9f26-138c871c812a
# ╠═4c05511c-d209-4213-a380-d16c60a834a4
# ╠═668bea30-c7e4-4f90-a0e8-fa315dc79cd1
# ╠═fc4794cf-e45b-41ee-acac-6d5dcf8c02ab
# ╠═0adc833e-8de6-4ab7-a137-9cb6e5aeba86
# ╠═cfc844b8-0274-4c38-a676-32dfd9bcb44b
# ╠═f382453c-51ab-4cb7-84a7-ee20354cce93
# ╠═d413132c-d6f3-45d4-8a1e-1306b8ad084b
# ╠═e0764d78-afcf-45fb-b07c-ce0485250f29
# ╠═c1ad6eab-eb98-470f-ab26-b30316e4b8ce
# ╠═29a63484-1998-4ce4-978b-ff26cb264ab5
