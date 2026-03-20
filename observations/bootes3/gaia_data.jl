### A Pluto.jl notebook ###
# v0.20.23

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
	include("../../utils/gaia_filters.jl")
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

# ╔═╡ d9fa38b6-ed7a-432a-9646-5da36a15ef06
style_kwargs = Dict(
	:best=>(; markersize=0.3, color=:grey, label="all"), 
	:members => (; markersize=1, color=COLORS[1], label="PM+CMD"),
	:rv => (; markersize=3, color=COLORS[4], label="velocity members")
)

# ╔═╡ d97384ad-5a20-4eb2-b4e5-7aa92cb146b9
obs_props = let
	df = TOML.parsefile("observed_properties.toml")
	df["position_angle"] += 180
	df
end

# ╔═╡ 1df23863-6068-49b8-891b-2492134a4cd1
params = let
	p = GaiaFilters.read_paramfile(
	joinpath("./density_profiles/fiducial.toml")
) 
	pop!(p, "profile_kwargs", Dict())

	p["filename"] = joinpath("density_profiles", p["filename"])
	p |> LilGuys.dict_to_tuple
end

# ╔═╡ cc0bee7b-2149-4c6d-be50-3d13df448d90
filt_params = GaiaFilters.GaiaFilterParams(obs_props; params... )

# ╔═╡ 471fe182-f422-4744-97a2-f37f24778ec4
all_stars = GaiaFilters.read_gaia_stars(filt_params)

# ╔═╡ 126f78fe-7bba-4b8f-ac02-2971957fdb04
vel_memb = read_fits("velocities/processed/rv_deimos_x_2c_psat_0.2.fits")

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

# ╔═╡ 18dbd5f2-9218-4b5d-806c-61312e8e74e5
sample_pm = good_stars[good_stars.LLR_PM .> 0, :]

# ╔═╡ 6f9e9dad-1e93-4cf1-99a3-17a3f1f11c57
R_cen = 60

# ╔═╡ 1270e9a2-85d2-4d79-968c-5aff48e3c636
R_bkg = 225

# ╔═╡ 5d9621c3-6b9c-4d1d-bfcc-6b9b592a90b2
sample_pm_close = sample_pm[sample_pm.R_ell .< R_cen, :]

# ╔═╡ d3a43028-c191-4b43-8d39-211b84aba4cb
sample_pm_far = sample_pm[R_bkg .< sample_pm.R_ell .< sqrt(R_bkg^2 + R_cen^2), :]

# ╔═╡ c10568c6-fe5f-435c-a54a-6b33eab16c80
members = good_stars[good_stars.LLR_nospace .> 0.0, :]

# ╔═╡ f4d59c5d-b869-4595-a15c-1a172319434b
md"""
# Sanity check plots
"""

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

# ╔═╡ 55fce147-81da-460b-a208-e6e5baedc15a
let
	fig = Figure()
	GaiaPlots.plot_tangent(fig[1,1], Dict(:best => sample_pm), Dict(:best => (; markersize=1, color=:grey)), nothing)

	fig
end

# ╔═╡ 926c3da9-d617-430b-92fc-93da60f89392
let
	fig = Figure()
	GaiaPlots.plot_tangent(fig[1,1], Dict(:best => members), Dict(:best => (; markersize=1, color=:grey)), nothing)

	fig
end

# ╔═╡ f575d602-c082-489e-8520-b8bc8d070369
let
	fig = Figure()
	ax = GaiaPlots.tangent_axis(fig[1,1], 300)

	scatter!(good_stars.xi, good_stars.eta, color=:grey, markersize=0.3)

	scatter!(sample_pm_close.xi, sample_pm_close.eta, alpha=0.5)
	scatter!(sample_pm_far.xi, sample_pm_far.eta, alpha=0.5)
	fig

end

# ╔═╡ c7273018-e832-47ee-971f-3a374216c5ee
let
	fig = Figure()
	ax = GaiaPlots.cmd_axis(fig[1,1])
	
	scatter!(sample_pm_space.bp_rp, sample_pm_space.G, markersize=1, alpha=1, color=:black)
	GaiaPlots.plot_iso!(obs_props)
	fig
end

# ╔═╡ 65ac975e-9add-40d7-ba12-afa1f6acc5db
let
	fig = Figure(size=(4, 3) .* 72)
	ax =  GaiaPlots.plot_cmd(fig[1,1], Dict(:selected => sample_pm_close), Dict(:selected => (;)), obs_props, age=12 )
	ax.title = "centre"


	ax =  GaiaPlots.plot_cmd(fig[1,2], Dict(:selected => sample_pm_far), Dict(:selected => (;)), obs_props, age=12 )
	ax.title = "background"

	hideydecorations!(ticks=false, minorticks=false)

	fig
end

# ╔═╡ 84d406e1-6f62-4e7b-9b51-3a674a6e5556
close_stars = good_stars[good_stars.R_ell .< 5.25 * 30, :]

# ╔═╡ 4dfc1b95-2813-4095-b844-26165a70736f
let
	fig = Figure(size=(4, 3) .* 72)
	ax = Axis(fig[1,1],
		xlabel = "BP - RP",
		ylabel = "G",
		yreversed = true,
			  limits=(-0.5, 2.5, 14, 21.2)
	)

	scatter!(close_stars.bp_rp, close_stars.G, markersize=0.3, alpha=1, color=:grey, label="<5.25R_h" => (;markersize= 2))

	p = scatter!(sample_pm_space.bp_rp, sample_pm_space.G, markersize=2.5, alpha=1,color=sample_pm_space.LLR_CMD, colorrange=(-5, 5), label="PM+Space")


	scatter!(vel_memb.bp_rp, vel_memb.G, markersize=5, alpha=1, color=COLORS[3], label="Vlos memb")


	axislegend(position=:lt)
	Colorbar(fig[1,2], p, label="LLR CMD")

	fig
end

# ╔═╡ 56b67ff8-257b-4e42-b167-f625923f4aa4
GaiaPlots.compare_j24_samples(Dict(:best => good_stars, :members => good_stars), Dict(:best => style_kwargs[:best], :members => style_kwargs[:best]), nothing)

# ╔═╡ 69ca64e6-df98-42d5-bb31-ce6d5831a410
GaiaPlots.compare_j24_samples(Dict(:best => close_stars,), Dict(:best => style_kwargs[:best], :members => style_kwargs[:best]), nothing)

# ╔═╡ 7b15f1ff-6cb0-4266-9cdf-3ed809b8e262
GaiaPlots.compare_j24_samples(Dict(:best => good_stars, :members=>members, :rv => vel_memb), style_kwargs, obs_props)

# ╔═╡ 90eb03de-2afe-4232-ab07-00b607431d6b
md"""
# Background objects
There are two globular clusters in the field. 
Fortunately, the Gaia stars do indeed fall nicely along
the expected (literature) observed properties
"""

# ╔═╡ cd492e8c-8a2e-42a3-866b-f2cc99e559ef
obs_props_ngc5466 = Dict(
	"ra" => 211.36371,
	"dec" => 28.53444,
	"pmra" => -5.41,
	"pmdec" => -0.79,
	"distance_modulus" => LilGuys.kpc2dm(16.9),
	"R_h" => 10.72,
	"ellipticity" => 0,
	"position_angle" => 0,
	"metallicity" => -2.0,
)

# ╔═╡ 1e75de81-9400-488c-8331-9eb30b1f2185
obs_props_ngc5272 = Dict(
	"ra" => 205.54842,
	"dec" => +28.37728,
	"R_h" => 16.2,
	"ellipticity" => 0,
	"pmra" => -0.14,
	"pmdec" => -2.64,
	"position_angle" => 0,
	"distance_modulus" => LilGuys.kpc2dm(10.06),
	"metallicity" => -1.5,
)

# ╔═╡ 3d110341-b33d-4e38-9d85-04ca7d08b605
xi_ngc5466, eta_ngc5466 = 60 .* LilGuys.to_tangent(obs_props_ngc5466["ra"], obs_props_ngc5466["dec"], obs_props["ra"], obs_props["dec"])

# ╔═╡ 06f2e86e-ac5d-4182-8987-d3a3d6a2ff95
sample_ngc5466 = let
	r_ngc = @. sqrt(
		(good_stars.xi - xi_ngc5466)^2 
		+ (good_stars.eta - eta_ngc5466) ^ 2
	)

	df = good_stars[r_ngc .< obs_props_ngc5466["R_h"] * 3, :]

	df.xi .-= xi_ngc5466
	df.eta .-= eta_ngc5466

	df
end

# ╔═╡ dc46ff9c-1ee6-4c38-89af-8fa8d4e5e4f5
xi_ngc5272, eta_ngc5272 = 60 .* LilGuys.to_tangent(obs_props_ngc5272["ra"], obs_props_ngc5272["dec"], obs_props["ra"], obs_props["dec"])

# ╔═╡ f6f67983-1c02-400c-9e1b-1e427dc78755
sample_ngc5272 = let
	r_ngc = @. sqrt(
		(good_stars.xi - xi_ngc5272)^2 
		+ (good_stars.eta - eta_ngc5272) ^ 2
	)

	df = good_stars[r_ngc .< obs_props_ngc5272["R_h"] * 3, :]
	df.xi .-= xi_ngc5272
	df.eta .-= eta_ngc5272

	df
end

# ╔═╡ 709a8dde-8b64-489c-ab42-742032c5d99c
let
	df = copy(good_stars)
	df.xi .-= xi_ngc5466
	df.eta .-= eta_ngc5466

	f = GaiaPlots.compare_j24_samples(Dict(:best => df, :members=>sample_ngc5466, ), style_kwargs, obs_props_ngc5466)

	f.content[1].title = "NGC 5466"
	f.content[1].limits = nothing, nothing

	f
end

# ╔═╡ 850b4a52-b83d-4ab2-a4b0-c137fb20161a
let
	df = copy(good_stars)
	df.xi .-= xi_ngc5272
	df.eta .-= eta_ngc5272

	f = GaiaPlots.compare_j24_samples(Dict(:best => df, :members=>sample_ngc5272, ), style_kwargs, obs_props_ngc5272)

	f.content[1].title = "NGC 5272"
	f.content[1].limits = nothing, nothing
	f
end

# ╔═╡ Cell order:
# ╠═51f6ff42-2350-11f1-9128-7106ae887ea8
# ╠═76976bdb-3d81-4986-ba20-8d2629daaf04
# ╠═f8bcf091-79a5-4c43-b815-92ac92c499e5
# ╠═af7426ac-64a3-42bd-a4dd-918627cc50b0
# ╠═647838f8-3523-40cd-b58d-41ec10dd9e0f
# ╠═002638e8-c36d-4772-b093-a8439f720070
# ╠═7a9c7ebc-e692-4e91-9e56-1b83be1154cf
# ╟─1b62a47f-8d1b-4271-9226-5abd3a3da299
# ╠═d9fa38b6-ed7a-432a-9646-5da36a15ef06
# ╠═d97384ad-5a20-4eb2-b4e5-7aa92cb146b9
# ╠═1df23863-6068-49b8-891b-2492134a4cd1
# ╠═cc0bee7b-2149-4c6d-be50-3d13df448d90
# ╠═471fe182-f422-4744-97a2-f37f24778ec4
# ╠═7a9f2106-6ce4-4d29-a0f9-d11090c128d2
# ╠═126f78fe-7bba-4b8f-ac02-2971957fdb04
# ╠═53b258be-c58a-4e2f-bc41-85c2e7eefc10
# ╠═387a5f4e-8f7b-4a4f-8502-2ab583b2ffcd
# ╠═88a65ca2-f4eb-487c-b5ce-6a9db04ffd70
# ╠═2c805d84-55e1-488d-8b71-655d8bb3db2d
# ╠═18dbd5f2-9218-4b5d-806c-61312e8e74e5
# ╠═6f9e9dad-1e93-4cf1-99a3-17a3f1f11c57
# ╠═1270e9a2-85d2-4d79-968c-5aff48e3c636
# ╠═5d9621c3-6b9c-4d1d-bfcc-6b9b592a90b2
# ╠═d3a43028-c191-4b43-8d39-211b84aba4cb
# ╠═c10568c6-fe5f-435c-a54a-6b33eab16c80
# ╟─f4d59c5d-b869-4595-a15c-1a172319434b
# ╠═93070414-65f1-4cf8-bf7b-6324841158f2
# ╠═5a2497e6-e634-4c96-ab96-2ef493676184
# ╠═55fce147-81da-460b-a208-e6e5baedc15a
# ╠═926c3da9-d617-430b-92fc-93da60f89392
# ╠═f575d602-c082-489e-8520-b8bc8d070369
# ╠═c7273018-e832-47ee-971f-3a374216c5ee
# ╠═4dfc1b95-2813-4095-b844-26165a70736f
# ╠═65ac975e-9add-40d7-ba12-afa1f6acc5db
# ╠═84d406e1-6f62-4e7b-9b51-3a674a6e5556
# ╠═56b67ff8-257b-4e42-b167-f625923f4aa4
# ╠═69ca64e6-df98-42d5-bb31-ce6d5831a410
# ╠═7b15f1ff-6cb0-4266-9cdf-3ed809b8e262
# ╠═90eb03de-2afe-4232-ab07-00b607431d6b
# ╠═cd492e8c-8a2e-42a3-866b-f2cc99e559ef
# ╠═1e75de81-9400-488c-8331-9eb30b1f2185
# ╠═3d110341-b33d-4e38-9d85-04ca7d08b605
# ╠═06f2e86e-ac5d-4182-8987-d3a3d6a2ff95
# ╠═dc46ff9c-1ee6-4c38-89af-8fa8d4e5e4f5
# ╠═f6f67983-1c02-400c-9e1b-1e427dc78755
# ╠═709a8dde-8b64-489c-ab42-742032c5d99c
# ╠═850b4a52-b83d-4ab2-a4b0-c137fb20161a
