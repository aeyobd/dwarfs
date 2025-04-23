### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ f27c2500-b92c-4de0-8482-4a8acfa2667f
module RVUtils
	include(joinpath(ENV["DWARFS_ROOT"], "observations/rv_utils.jl"))
end

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
df_gradient = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_gradient.csv", DataFrame)

# ╔═╡ b65c1d81-b681-47e0-ab9a-210dc6a98f6e
df_xi_p = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/vz_xi_p_binned.csv", DataFrame)

# ╔═╡ bc3e07a1-ff18-44ae-8f5e-3abc043c5769
df_xi_p[:, :μ_err2] = [tuple(parse.(Float64, split(s[2:end-1], ", "))...) for s in df_xi_p.μ_err]

# ╔═╡ 668bea30-c7e4-4f90-a0e8-fa315dc79cd1
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ 75ed46bf-6550-44d4-a15e-fbbd59c1ba12
@savefig "scl_vel_gradient_binned" let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degree",
		ylabel = L"$\bar{v}_z$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x ./ 60, df_xi_p.μ, yerror=df_xi_p.μ_err2, color=:black)

	for i in 1:400:size(df_gradient, 1)
		M = 60df_gradient.B[i] .* sind(df_gradient.Θ_grad[i]) .+ 60df_gradient.A[i] .* cosd(df_gradient.Θ_grad[i])
		x=LinRange(-0.4, 0.4, 100)
		y = M*x  .+ df_gradient.μ[i]
		lines!(x, y, alpha=0.03, color=COLORS[1], linewidth=1)
	end
			
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═314da249-8a3f-47e9-ac40-f86406ec9955
# ╠═71335ca6-b58b-4c68-8467-81cd4f481b1f
# ╠═7fe7b9af-f22c-416e-b9f4-4a4d3cba82df
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f27c2500-b92c-4de0-8482-4a8acfa2667f
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═b65c1d81-b681-47e0-ab9a-210dc6a98f6e
# ╠═bc3e07a1-ff18-44ae-8f5e-3abc043c5769
# ╠═668bea30-c7e4-4f90-a0e8-fa315dc79cd1
# ╠═75ed46bf-6550-44d4-a15e-fbbd59c1ba12
