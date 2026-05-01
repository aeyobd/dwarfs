### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 51f6ff42-2350-11f1-9128-7106ae887ea8
begin
	import Pkg
	Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ af7426ac-64a3-42bd-a4dd-918627cc50b0
using CSV, DataFrames

# ╔═╡ 7a9f2106-6ce4-4d29-a0f9-d11090c128d2
using PyFITS

# ╔═╡ 76976bdb-3d81-4986-ba20-8d2629daaf04
CairoMakie.activate!(type=:png)

# ╔═╡ f8bcf091-79a5-4c43-b815-92ac92c499e5
import TOML

# ╔═╡ 647838f8-3523-40cd-b58d-41ec10dd9e0f
module GaiaFilters
	include(joinpath(ENV["DWARFS_ROOT"], "utils/gaia_filters.jl"))
end

# ╔═╡ 002638e8-c36d-4772-b093-a8439f720070
module GaiaPlots
	include("gaia_utils.jl")
end

# ╔═╡ 7a9c7ebc-e692-4e91-9e56-1b83be1154cf
LilGuys.error_interval(x) = (0., 0.)

# ╔═╡ 1b62a47f-8d1b-4271-9226-5abd3a3da299
md"""
# Data loading
"""

# ╔═╡ 126f78fe-7bba-4b8f-ac02-2971957fdb04
vel_memb = read_fits(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3", "velocities/processed/rv_deimos_x_2c_psat_0.2.fits"))

# ╔═╡ d97384ad-5a20-4eb2-b4e5-7aa92cb146b9
obs_props = let
	df = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3", "observed_properties.toml"))
	df
end

# ╔═╡ 1df23863-6068-49b8-891b-2492134a4cd1
params = let
	p = GaiaFilters.read_paramfile(
	joinpath(ENV["DWARFS_ROOT"], "observations/bootes3", "samples/fiducial.toml")
) 
	pop!(p, "profile_kwargs", Dict())

	p["filename"] = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3", "density_profiles", p["filename"])
	p |> LilGuys.dict_to_tuple
end

# ╔═╡ cc0bee7b-2149-4c6d-be50-3d13df448d90
filt_params = GaiaFilters.GaiaFilterParams(obs_props; params... )

# ╔═╡ 471fe182-f422-4744-97a2-f37f24778ec4
all_stars = GaiaFilters.read_gaia_stars(filt_params)

# ╔═╡ 53b258be-c58a-4e2f-bc41-85c2e7eefc10
md"""
# Subsets
"""

# ╔═╡ 387a5f4e-8f7b-4a4f-8502-2ab583b2ffcd
good = GaiaFilters.select_members(all_stars, filt_params)

# ╔═╡ 88a65ca2-f4eb-487c-b5ce-6a9db04ffd70
good_stars = all_stars[all_stars.F_BEST .== 1, :]

# ╔═╡ 2c805d84-55e1-488d-8b71-655d8bb3db2d
sample_pm_space = good_stars[good_stars.LLR_PM .+ good_stars.LLR_S .> 0, :]

# ╔═╡ c10568c6-fe5f-435c-a54a-6b33eab16c80
members = good_stars[good_stars.LLR_nospace .> 1, :]

# ╔═╡ f4d59c5d-b869-4595-a15c-1a172319434b
md"""
# Sanity check plots
"""

# ╔═╡ d9fa38b6-ed7a-432a-9646-5da36a15ef06
style_kwargs = Dict(
	:best=>(; markersize=0.3, color=:grey, label="all" => (; markersize=1)), 
	:members => (; markersize=1, color=COLORS[1], label="PM+CMD" => (;markersize=2)),
	:rv => (; markersize=3, color=COLORS[4], label="velocity members" =>(;markersize=4))
)

# ╔═╡ 93070414-65f1-4cf8-bf7b-6324841158f2
let
	fig = Figure()
	GaiaPlots.plot_tangent(fig[1,1], Dict(:best => all_stars), style_kwargs, nothing)

	fig
end

# ╔═╡ 5a2497e6-e634-4c96-ab96-2ef493676184
let
	fig = Figure()
	GaiaPlots.plot_tangent(fig[1,1], Dict(:best => good_stars), style_kwargs, nothing)

	fig
end

# ╔═╡ 926c3da9-d617-430b-92fc-93da60f89392
let
	fig = Figure()
	GaiaPlots.plot_tangent(fig[1,1], Dict(:best => members), Dict(:best => (; markersize=1, color=:grey)), nothing)

	@savefig "gaia_memb_tangent"
	fig
end

# ╔═╡ 84d406e1-6f62-4e7b-9b51-3a674a6e5556
close_stars = good_stars[good_stars.R_ell .< 5.25 * 30, :]

# ╔═╡ 7b15f1ff-6cb0-4266-9cdf-3ed809b8e262
@savefig "gaia_tangent_cmd_pm" GaiaPlots.compare_j24_samples(Dict(:best => good_stars, :members=>members, :rv => vel_memb), style_kwargs, obs_props)

# ╔═╡ Cell order:
# ╠═51f6ff42-2350-11f1-9128-7106ae887ea8
# ╠═76976bdb-3d81-4986-ba20-8d2629daaf04
# ╠═f8bcf091-79a5-4c43-b815-92ac92c499e5
# ╠═af7426ac-64a3-42bd-a4dd-918627cc50b0
# ╠═647838f8-3523-40cd-b58d-41ec10dd9e0f
# ╠═002638e8-c36d-4772-b093-a8439f720070
# ╠═7a9c7ebc-e692-4e91-9e56-1b83be1154cf
# ╠═7a9f2106-6ce4-4d29-a0f9-d11090c128d2
# ╟─1b62a47f-8d1b-4271-9226-5abd3a3da299
# ╠═126f78fe-7bba-4b8f-ac02-2971957fdb04
# ╠═d97384ad-5a20-4eb2-b4e5-7aa92cb146b9
# ╠═1df23863-6068-49b8-891b-2492134a4cd1
# ╠═cc0bee7b-2149-4c6d-be50-3d13df448d90
# ╠═471fe182-f422-4744-97a2-f37f24778ec4
# ╟─53b258be-c58a-4e2f-bc41-85c2e7eefc10
# ╠═387a5f4e-8f7b-4a4f-8502-2ab583b2ffcd
# ╠═88a65ca2-f4eb-487c-b5ce-6a9db04ffd70
# ╠═2c805d84-55e1-488d-8b71-655d8bb3db2d
# ╠═c10568c6-fe5f-435c-a54a-6b33eab16c80
# ╟─f4d59c5d-b869-4595-a15c-1a172319434b
# ╠═d9fa38b6-ed7a-432a-9646-5da36a15ef06
# ╠═93070414-65f1-4cf8-bf7b-6324841158f2
# ╠═5a2497e6-e634-4c96-ab96-2ef493676184
# ╠═926c3da9-d617-430b-92fc-93da60f89392
# ╠═84d406e1-6f62-4e7b-9b51-3a674a6e5556
# ╠═7b15f1ff-6cb0-4266-9cdf-3ed809b8e262
