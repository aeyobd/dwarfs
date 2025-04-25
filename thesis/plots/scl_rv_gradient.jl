### A Pluto.jl notebook ###
# v0.20.6

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

# ╔═╡ 535d6048-d83b-4439-9fbb-40616ad7b35b
using RollingFunctions

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

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
df_gradient = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_gradient.csv", DataFrame)

# ╔═╡ b65c1d81-b681-47e0-ab9a-210dc6a98f6e
df_xi_p = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/vz_xi_p_binned.csv", DataFrame)

# ╔═╡ 5927a386-fdb1-459f-b042-0f95677e345e
rv_memb = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_members_all.fits")

# ╔═╡ 33ee935f-36f6-424c-9f26-138c871c812a
minimum(rv_memb.PSAT_RV) 

# ╔═╡ bc3e07a1-ff18-44ae-8f5e-3abc043c5769
df_xi_p[:, :μ_err2] = [tuple(parse.(Float64, split(s[2:end-1], ", "))...) for s in df_xi_p.μ_err]

# ╔═╡ 668bea30-c7e4-4f90-a0e8-fa315dc79cd1
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ 4c05511c-d209-4213-a380-d16c60a834a4
xi_p, eta_p = LilGuys.to_orbit_coords(rv_memb.ra, rv_memb.dec, obs_props["ra"], obs_props["dec"], θ_m)

# ╔═╡ fc4794cf-e45b-41ee-acac-6d5dcf8c02ab
function rolling_medians(x, y; window=50)
	filt = sortperm(x)
	xs = x[filt]
	ys = y[filt]

    n = length(x)
    medians = Vector{Float64}(undef, n)
	x_medians = Vector{Float64}(undef, n)

	for i in 1:n
		start = max(1, i-window+1)
		fin = min(start + i, n)
		medians[i] = median(ys[start:fin])
		x_medians[i] = median(xs[start:fin])
	end

	return x_medians[window:end-window], medians[window:end-window]
end

# ╔═╡ ebea53d4-c421-4c11-9c02-3b8072a3d3f0
LilGuys.mean(rv_memb.vz), gsr.radial_velocity

# ╔═╡ cfc844b8-0274-4c38-a676-32dfd9bcb44b
filt = (rv_memb.PSAT_RV .> 0.2)

# ╔═╡ 51fc1b13-1cf6-4e7d-b156-1a5bd44cf17b
rolling_medians(xi_p[filt], rv_memb.vz[filt])

# ╔═╡ 75ed46bf-6550-44d4-a15e-fbbd59c1ba12
@savefig "scl_vel_gradient_scatter" let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degree",
		ylabel = L"$v_z$ / km s$^{-1}$",
		limits=(-0.8, 0.8, 35, 115)
	)
	scatter!(xi_p[filt] ./ 60, rv_memb.vz[filt], markersize=2, color=:black, label="members")

	
	for i in 1:400:size(df_gradient, 1)
		M = 60df_gradient.B[i] .* sind(df_gradient.Θ_grad[i]) .+ 60df_gradient.A[i] .* cosd(df_gradient.Θ_grad[i])
		x=LinRange(-0.4, 0.4, 100)
		y = M*x  .+ df_gradient.μ[i]
		lines!(x, y, alpha=0.03, color=COLORS[1], linewidth=1, label="MCMC samples" => (; alpha=1.0))
	end

	hlines!(gsr.radial_velocity, color=:black, linewidth=1, alpha=0.3, label="systematic velocity")

	xs, ys = rolling_medians(xi_p[filt] ./ 60, rv_memb.vz[filt],)
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
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═b65c1d81-b681-47e0-ab9a-210dc6a98f6e
# ╠═5927a386-fdb1-459f-b042-0f95677e345e
# ╠═33ee935f-36f6-424c-9f26-138c871c812a
# ╠═bc3e07a1-ff18-44ae-8f5e-3abc043c5769
# ╠═4c05511c-d209-4213-a380-d16c60a834a4
# ╠═668bea30-c7e4-4f90-a0e8-fa315dc79cd1
# ╠═535d6048-d83b-4439-9fbb-40616ad7b35b
# ╠═fc4794cf-e45b-41ee-acac-6d5dcf8c02ab
# ╠═51fc1b13-1cf6-4e7d-b156-1a5bd44cf17b
# ╠═ebea53d4-c421-4c11-9c02-3b8072a3d3f0
# ╠═75ed46bf-6550-44d4-a15e-fbbd59c1ba12
# ╠═cfc844b8-0274-4c38-a676-32dfd9bcb44b
