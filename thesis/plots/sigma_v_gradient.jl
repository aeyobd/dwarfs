### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ d4689051-9ef1-46f6-8438-2c5b1de773bd
obs_props_umi = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/observed_properties.toml"))

# ╔═╡ 799fd372-b4ea-473f-93cb-65eff73b9bbc
df_samples_scl = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/mcmc_samples_Rell_sigma.rv_combined_x_wide_2c_psat_0.2.csv", DataFrame)

# ╔═╡ c714dc7d-64fe-4ce4-a6c9-54fd65923960
df_gradient_scl = CSV.read(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/vz_r_ell_binned.rv_combined_x_wide_2c_psat_0.2.csv", DataFrame)

# ╔═╡ 9235fb6d-2b18-4243-8e98-18ba1f65dfb1
df_samples_umi = CSV.read(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/mcmc_samples_Rell.rv_combined_x_2c_psat_0.2.csv", DataFrame)

# ╔═╡ 272efc45-70ad-4b11-b145-41350c3fa791
df_gradient_umi = CSV.read(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/vz_r_ell_binned.rv_combined_x_2c_psat_0.2.csv", DataFrame)

# ╔═╡ d6236fa7-4f45-49b1-8fc4-82b08f196df5
function plot_obs!(df)
	x_mid = log10.(df.x)
	log_bin_errs = [ (x_mid[i] .- log10.(df.x_low[i]), log10(df.x_high[i]) .- x_mid[i]) for i in eachindex(x_mid)]
	y_err = collect(zip(df.σ_em ./ df.σ ./ log(10), df.σ_ep ./ df.σ ./ log(10)))

	errorscatter!(x_mid, log10.(df.σ), xerror=log_bin_errs, yerror=y_err, color=:black)
end

# ╔═╡ 8e6b1448-6a93-4180-afcb-208e1150a61c
function plot_samples!(df_Rell; skip=800)
	x = LinRange(-1, 3, 100)

	
	for i in 1:skip:size(df_Rell, 1)
		y = df_Rell.σ[i] .* 10 .^ ((x .- 1.0) .* df_Rell.dlσ_dlR[i])
		lines!(x, log10.(y), alpha=0.03, color=COLORS[3], linestyle=:solid)
	end

	hlines!(log10(median(df_Rell.σ)), color=:black, linewidth=theme(:linewidth)[]/2)
end

# ╔═╡ 0433a0dd-6f7b-49dc-998c-90d89a694da3
let
	fig, ax = FigAxis(
		#xlabel = L"$\log\ R$ / arcmin",
		ylabel = L"$\log\ \sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	text!(0, 1, offset=(6, -6), text="Sculptor", space=:relative, align=(:left, :top))

	plot_samples!(df_samples_scl)

	plot_obs!(df_gradient_scl)	

	xlims!(-1, 2.5)
	ylims!(0.8, 1.1)
	annotation!(0, 36, log10(obs_props["R_h"]), 0.8, text=L"R_h")


	ax2 = Axis(fig[2,1],
		xlabel = L"$\log\ R$ / arcmin",
		ylabel = L"$\log\ \sigma_{v, \textrm{los}}$ / km s$^{-1}$"

			  )

	plot_samples!(df_samples_umi)

	plot_obs!(df_gradient_umi)	
	text!(0, 1, offset=(6, -6), text="Ursa Minor", space=:relative, align=(:left, :top))

	annotation!(0, 36, log10(obs_props_umi["R_h"]), 0.6,)
	text!(log10(obs_props_umi["R_h"]), 0.6, text=L"R_h", offset=(6, 18))

	xlims!(-1, 2.5)

	@savefig "sigma_v_gradient" 
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
# ╠═d4689051-9ef1-46f6-8438-2c5b1de773bd
# ╠═799fd372-b4ea-473f-93cb-65eff73b9bbc
# ╠═c714dc7d-64fe-4ce4-a6c9-54fd65923960
# ╠═9235fb6d-2b18-4243-8e98-18ba1f65dfb1
# ╠═272efc45-70ad-4b11-b145-41350c3fa791
# ╠═0433a0dd-6f7b-49dc-998c-90d89a694da3
# ╠═d6236fa7-4f45-49b1-8fc4-82b08f196df5
# ╠═8e6b1448-6a93-4180-afcb-208e1150a61c
