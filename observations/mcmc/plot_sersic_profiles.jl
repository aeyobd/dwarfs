### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using LilGuys
end

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
include("../../utils/gaia_filters.jl")

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	galaxy = "bootes3"
end

# ╔═╡ af15cb9c-c722-4aea-bd2a-155dadf53efd
import StatsBase: median

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath("..", galaxy, "figures"); FIGSUFFIX=".mcmc_hist"

# ╔═╡ 57a19d65-59d9-46cf-8916-d9ac3a4dc92b
datafile = let
	dir = joinpath("..", galaxy, "data")

	filenames = ["jensen+24_1c.fits", "jensen+24_2c.fits", "jensen+24_wide.fits", "j24_1c.fits", "j24_2c.fits"]
	filename = ""
	for file in filenames
		if isfile(joinpath(dir, file))
			filename =  joinpath(dir, file)
		end
	end

	filename
end

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type="svg", pt_per_unit=2)

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ c58ed939-07c1-4d34-90a4-3c5d9cc9910b
prof_fiducial = SurfaceDensityProfile(joinpath("..", galaxy, "density_profiles", "fiducial_profile.toml")) |> LilGuys.filter_empty_bins

# ╔═╡ b9cf2c23-6c9a-44f5-9740-22caf4959831
obs_props = TOML.parsefile(joinpath("..", galaxy, "observed_properties.toml"))

# ╔═╡ 6f016a8e-38ae-4f05-a7ee-c292ac0e5741
df_chains = CSV.read(joinpath("..", galaxy, "mcmc", "samples.mcmc_sersic_position.csv"), DataFrame)

# ╔═╡ 0aa561a6-cd6f-4196-90a2-6aa84a71d45b
f_sat = median(df_chains.f_sat)

# ╔═╡ 87d6917f-d751-4783-8c72-136a0589f998
Ntot = sum(prof_fiducial.counts) / f_sat

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ 86d85142-967b-45a0-b314-fb7a537c1b20
Σ_range = minimum(prof_fiducial.log_Sigma[isfinite.(prof_fiducial.log_Sigma)]).middle - 2, maximum(prof_fiducial.log_Sigma).middle + 1

# ╔═╡ 2fe96489-49df-4362-ac90-b4f10899c534
log_r_range = prof_fiducial.log_R_bins[1] - 0.5, prof_fiducial.log_R_bins[end] + 0.5

# ╔═╡ f2c43c30-b069-4d17-8d61-be801d85b245
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel =log_Sigma_label,
		limits = (log_r_range, Σ_range)
	)
	skip = 100

	for i in 1:skip:size(df_chains, 1)
		row = df_chains[i, :]
		
		x = LinRange(log_r_range..., 100)

		prof = LilGuys.Sersic(n=row.n, R_h=row.R_h, )
		y = LilGuys.surface_density.(prof, 10 .^ x)
		M = LilGuys.mass_2D(prof, 100)

		lines!(x, log10.(y ./ M * Ntot * row.f_sat), color=:black, alpha=0.03)
	end
	LilGuys.plot_log_Σ!(ax, prof_fiducial, color=Arya.COLORS[2])

	@savefig "profiles.mcmc.sersic"
	fig
end

# ╔═╡ 3e9b9577-ee25-4d7f-a8f5-03c7ee467542
import PairPlots: pairplot

# ╔═╡ 1304e2ec-dd1c-4a19-8914-15ba0f93dfe2
pairplot(df_chains[:, [:d_xi, :d_eta]], labels=Dict(:d_xi=>L"\Delta\xi / '", :d_eta=>L"\Delta\eta / '"))

# ╔═╡ bc129d32-4d2c-4ddf-a2f1-f372e13ad64d
pairplot(df_chains[:, [:log_R_h, :n, :position_angle, :ellipticity]])

# ╔═╡ a6af6075-0b8d-41a9-8d13-edb727d978e5
df_chains[!, :log_R_h] = log10.(df_chains.R_h)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═af15cb9c-c722-4aea-bd2a-155dadf53efd
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═57a19d65-59d9-46cf-8916-d9ac3a4dc92b
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╠═b9cf2c23-6c9a-44f5-9740-22caf4959831
# ╠═6f016a8e-38ae-4f05-a7ee-c292ac0e5741
# ╠═87d6917f-d751-4783-8c72-136a0589f998
# ╠═0aa561a6-cd6f-4196-90a2-6aa84a71d45b
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═86d85142-967b-45a0-b314-fb7a537c1b20
# ╠═2fe96489-49df-4362-ac90-b4f10899c534
# ╠═f2c43c30-b069-4d17-8d61-be801d85b245
# ╠═3e9b9577-ee25-4d7f-a8f5-03c7ee467542
# ╠═1304e2ec-dd1c-4a19-8914-15ba0f93dfe2
# ╠═bc129d32-4d2c-4ddf-a2f1-f372e13ad64d
# ╠═a6af6075-0b8d-41a9-8d13-edb727d978e5
