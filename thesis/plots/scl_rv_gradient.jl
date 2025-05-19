### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 0ef64034-7d3e-4270-a36e-f4d6b5f7eaaf
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ 47ddbd9b-6d22-4c11-85d8-9b3f672e5287
gsr = LilGuys.transform(GSR, ICRS(obs_props))

# ╔═╡ 2f067cc4-9793-44e0-885e-2875d36e2a20
filename = "rv_combined_x_wide_2c_psat_0.2"

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
df_gradient = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_gradient.$filename.csv", DataFrame)

# ╔═╡ 5927a386-fdb1-459f-b042-0f95677e345e
rv_memb = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/$filename.fits")

# ╔═╡ 33ee935f-36f6-424c-9f26-138c871c812a
minimum(rv_memb.PSAT_RV) 

# ╔═╡ 668bea30-c7e4-4f90-a0e8-fa315dc79cd1
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ 4c05511c-d209-4213-a380-d16c60a834a4
xi_p, eta_p = LilGuys.to_orbit_coords(rv_memb.ra, rv_memb.dec, obs_props["ra"], obs_props["dec"], θ_m)

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

# ╔═╡ dbebfa62-77df-4e38-8fc6-47e36dc6e7b6
length(50:100)

# ╔═╡ ebea53d4-c421-4c11-9c02-3b8072a3d3f0
LilGuys.mean(rv_memb.vz), gsr.radial_velocity

# ╔═╡ cfc844b8-0274-4c38-a676-32dfd9bcb44b
filt = (rv_memb.PSAT_RV .> 0.2)

# ╔═╡ 51fc1b13-1cf6-4e7d-b156-1a5bd44cf17b
rolling_medians(xi_p[filt], rv_memb.vz[filt])

# ╔═╡ 7c723942-5376-4d34-8acb-3a43c181e60a
sort(xi_p[filt])[[300, end-300]]

# ╔═╡ 75ed46bf-6550-44d4-a15e-fbbd59c1ba12
@savefig "scl_vel_gradient_scatter" let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degree",
		ylabel = L"$v_\textrm{gsr}'$ / km s$^{-1}$",
		limits=(-0.8, 0.8, 35, 115)
	)
	scatter!(xi_p[filt], rv_memb.vz[filt], markersize=2, color=:black, label="members")

	
	for i in 1:400:size(df_gradient, 1)
		M = 60df_gradient.B[i] .* sind(df_gradient.Θ_grad[i]) .+ 60df_gradient.A[i] .* cosd(df_gradient.Θ_grad[i])
		x=LinRange(-0.4, 0.4, 100)
		y = M*x  .+ df_gradient.μ[i]
		lines!(x, y, alpha=0.03, color=COLORS[1], linewidth=1, label="MCMC samples" => (; alpha=1.0))
	end

	hlines!(gsr.radial_velocity, color=:black, linewidth=1, alpha=0.3, label="systematic velocity")

	xs, ys = rolling_medians(xi_p[filt] , rv_memb.vz[filt], window=200)
	lines!(xs, ys,  color=COLORS[2], label="rolling median")

	axislegend(position=:lb, unique=true, merge=true)
	
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
# ╠═0ef64034-7d3e-4270-a36e-f4d6b5f7eaaf
# ╠═47ddbd9b-6d22-4c11-85d8-9b3f672e5287
# ╠═2f067cc4-9793-44e0-885e-2875d36e2a20
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═5927a386-fdb1-459f-b042-0f95677e345e
# ╠═33ee935f-36f6-424c-9f26-138c871c812a
# ╠═4c05511c-d209-4213-a380-d16c60a834a4
# ╠═668bea30-c7e4-4f90-a0e8-fa315dc79cd1
# ╠═fc4794cf-e45b-41ee-acac-6d5dcf8c02ab
# ╠═dbebfa62-77df-4e38-8fc6-47e36dc6e7b6
# ╠═51fc1b13-1cf6-4e7d-b156-1a5bd44cf17b
# ╠═7c723942-5376-4d34-8acb-3a43c181e60a
# ╠═ebea53d4-c421-4c11-9c02-3b8072a3d3f0
# ╠═75ed46bf-6550-44d4-a15e-fbbd59c1ba12
# ╠═cfc844b8-0274-4c38-a676-32dfd9bcb44b
