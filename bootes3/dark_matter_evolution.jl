### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import TOML

	using LilGuys
	using Arya
	using PyFITS

	using PlutoUI
	import LinearAlgebra: cross
end

# ╔═╡ 3e09298c-86c8-49b6-aaa8-516b9f199478
include(joinpath(ENV["DWARFS_ROOT"], "utils/pluto_utils.jl"))

# ╔═╡ bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
md"""
Analyzes the dark matter particles and profiles for the simulation.
In particular, we plot the initial/final density profiles, 
circular velocity evolution, 
and make some nice projections of the DM in different orientations.
All of the figures are saved to figures directory inside the model analysis directory. 
"""

# ╔═╡ d3bd61d8-1e90-4787-b892-d90717f6be6e
@bind inputs confirm(notebook_inputs(;
	modelname = TextField(60, default="bootes3/1e5_v30_r2.2/orbit_"),
))

# ╔═╡ 99f96d71-b543-4680-a022-2195e6cca897
md"""
# Setup
"""

# ╔═╡ 2b2a1cc7-d005-4bef-b2cf-5d26b8c203a0
CairoMakie.activate!(type=:png)

# ╔═╡ 9c4d9492-64bc-4212-a99d-67cc507e99e0
md"""
## Inputs
"""

# ╔═╡ 3db38875-fe22-4cfd-8c5a-47f4a0fa7f3a
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", inputs.modelname)

# ╔═╡ 1b87d662-da3c-4438-98eb-72dc93e32f6a
FIGDIR = joinpath(model_dir, "figures")

# ╔═╡ c260ee35-7eed-43f4-b07a-df4371397195
readdir(model_dir)

# ╔═╡ d010a230-7331-4afd-86dc-380da0e0f720
halo = LilGuys.load_profile(joinpath(model_dir, "../halo.toml"))

# ╔═╡ 510706ac-ffbd-4996-af9e-67f1b910d51c
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 50bb0dfa-7332-484f-a27d-d8491413ef1e
df_scalars = read_fits(joinpath(model_dir, "profiles_scalars.fits"))

# ╔═╡ 7094bc54-deb4-48a5-bf09-9ee6c684ac3c
out = Output(model_dir)

# ╔═╡ 2470e05f-9215-45e4-88fc-daab0638272f
begin 
	profiles = LilGuys.read_ordered_structs(joinpath(model_dir, "profiles.hdf5"), LilGuys.MassProfile)
	snap_idx = first.(profiles)
	profiles = last.(profiles);
end

# ╔═╡ b0e336df-678a-4406-b294-0c353f3c0c38
# dens_profiles = LilGuys.read_ordered_structs(joinpath(model_dir, "profiles_densities.hdf5"), LilGuys.DensityProfile) .|> last

# ╔═╡ ef37f7b9-15a2-412f-a254-6a0efac05a04
md"""
# Final Position
"""

# ╔═╡ 5115a3fc-e71d-4806-a34a-aaa3b00c2d4e
obs_gc = [Galactocentric(out.x_cen[:, i], out.v_cen[:, i]*V2KMS) for i in eachindex(out)]

# ╔═╡ 0526f791-9aff-4ba1-912c-97f74759e72b
obs_c = LilGuys.to_frame(LilGuys.transform.(ICRS, obs_gc))

# ╔═╡ f835ba02-3e4b-4816-8480-452bfa308f5d
orbit_props

# ╔═╡ a9e79439-16a4-4908-bfe0-f0770cdb26df
md"""
# Mass evolution
"""

# ╔═╡ 53641449-c5b3-45ff-a692-a5cd717c8369
idx_f = min(orbit_props["idx_f"], length(out))

# ╔═╡ 9429b4e0-96d0-4fae-a7d6-44737d568f76
if idx_f < orbit_props["idx_f"]
	@warn "snapshots do not line up:/"
end

# ╔═╡ 7e3df305-9678-447e-a48e-f102cf6ebced
idx_i = 2

# ╔═╡ 9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
snap_i = out[idx_i]

# ╔═╡ 8d127679-401c-439d-913d-e2020df1c600
snap_f = out[idx_f]

# ╔═╡ 4977303f-b958-4d24-9a04-0f2835137d37
times = out.times * T2GYR

# ╔═╡ 485eab53-43ec-4591-a3af-9e4cbfefbbe2
df_scalars.bound_mass[end] ./ df_scalars.bound_mass[1]

# ╔═╡ f3b1fd8e-0591-4d51-94ea-2b4eb65b6a71
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "bound mass",
			 yscale=log10, yticks=Makie.automatic)


	lines!(df_scalars.time * T2GYR .- df_scalars.time[idx_f] * T2GYR, df_scalars.bound_mass ./ df_scalars.bound_mass[1])
	fig
end

# ╔═╡ e14fa4a1-6175-4b9f-ad01-525c1617fe63
md"""
# Evolution of circular Velocity
"""

# ╔═╡ 8dae2e01-652b-4afc-b040-dd2ba1c6eedb
prof_i = MassProfile(snap_i)

# ╔═╡ b64c1caf-9ee0-4633-bd52-0258557b8847
prof_f = MassProfile(snap_f)

# ╔═╡ e6fc3297-c1b7-40a4-b2bb-98490a42604a
v_max = df_scalars.v_circ_max

# ╔═╡ d57501a1-4764-4b23-962f-2d37547d7bcc
r_max = df_scalars.r_circ_max

# ╔═╡ 04ca92d1-f64b-4d2c-a079-30f89866fda9
prof_i

# ╔═╡ 6f7eaf57-4ada-4995-9b88-52c984573b13
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", split(inputs.modelname, "/")[1], "observed_properties.toml"))

# ╔═╡ a098c336-fda6-4932-a0ef-aa0f139a4499
obs_today = let 
	properties_file = joinpath(model_dir, "simulation/orbit.toml")
	
	obs_today = TOML.parsefile(properties_file)

	rh = obs_props["R_h"]

	obs_today["ra_err"] = rh / 60 
	obs_today["dec_err"] = rh / 60
	
	obs_today
end

# ╔═╡ 092460b1-ea96-4609-be91-3625d0e5f459
let

	fig = Figure(size=(4*72, 6*72))
	for i in 1:3
		x, y = [("ra", "dec"), ("pmra", "pmdec"), ("distance", "radial_velocity")][i]
		ax = Axis(fig[i,1],
			xlabel=x,
			ylabel=y
		)
	
		idx = idx_f-3:min(idx_f+3,length(out)) 
	
		scatterlines!(obs_c[idx, x], obs_c[idx, y], color=idx .- idx_f, colormap=:twilight, colorrange=(-3, 3))
		
		errorscatter!([obs_today[x]], [obs_today[y]],
			xerror=[obs_today[x * "_err"]], yerror=[obs_today[y * "_err"]]
		)
	
	end

	@savefig "skyorbit_agreement_today"
	fig
end

# ╔═╡ c73ebaa4-16d9-428c-9967-6f8e67cdc680
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 82621ca0-2e0b-455a-b692-ab9a27c2a589
prof_stars = LilGuys.Exp2D(R_s=LilGuys.arcmin2kpc(obs_props["R_h"] / α_exp, orbit_props["distance_f"]) )

# ╔═╡ db320665-f46d-4aed-a2b2-4b39bcb605c5
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log \; r_\textrm{circ}\ /\ \textrm{kpc}", 
		ylabel=L"$v_\textrm{circ}$ / km s$^{-1}$",
		yscale=log10,
		yticks=[1, 10, 20, 30, 40, 50, 60],
		yminorticks=[1:9; 10:2:60],
		limits=((-2, 3), (1, 120)),
		xgridvisible=false,
		ygridvisible=false
	)

	r_model = 10 .^ LinRange(-2, 3, 1000)
	v_model = v_circ.(halo, r_model)
	lines!(log10.(r_model), v_model * V2KMS, linestyle=:dot, label="analytic")
	
	lines!(log10.(prof_i.radii), LilGuys.circular_velocity(prof_i) * V2KMS, label="initial")

	lines!(log10.(prof_f.radii), LilGuys.circular_velocity(prof_f) * V2KMS, label="final")
	vlines!(log10.(LilGuys.r_h(prof_stars)))


	α = 0.4
	β = 0.65
	x = LinRange(1, 0.01, 100)

	y = @. 2^α * x^β * (1 + x^2)^(-α)
	lines!(log10.(x .* r_max[1]), y .* v_max[1] * V2KMS,  label="EN21",
	color=:black, linestyle=:dash)

	scatter!(log10.(skipmissing(r_max)), skipmissing(v_max) .* V2KMS, color=Arya.COLORS[4], label=L"v_\textrm{circ,\ max}")

		
	axislegend(ax, position=:rt)
	@savefig "v_circ_profiles"
	fig
end

# ╔═╡ 54e28676-c63d-4802-9a6f-4706654f3972
obs_props["distance"]

# ╔═╡ 1630a217-a6c5-439b-a79b-10b9a098ea29
LilGuys.σv_1d(snap_f) * V2KMS

# ╔═╡ cabdac5e-9a2c-47f6-ba7f-c408b3af0f3d
halos = [LilGuys.TruncNFW(r_circ_max=r_max[i], v_circ_max=v_max[i], trunc=20, xi=3) for i in eachindex(r_max)]

# ╔═╡ 90f916db-3390-460e-8639-d5d02603aede
σv_predicted = LilGuys.σv_star_mean.(halos, [prof_stars])

# ╔═╡ 245721a6-01aa-43e7-922d-ed5da02207c1
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel=L"\mathrm{v}\,/\, \text{km\,s^{-1}}",
	limits=(nothing, (0, nothing)))
	x = times .- times[idx_f]
	lines!(x, v_max * V2KMS, label="max")

	lines!(x, σv_predicted * V2KMS, label="stellar dispersion")

	hspan!(obs_props["sigma_v"] - obs_props["sigma_v_em"], obs_props["sigma_v"] + obs_props["sigma_v_ep"], color=:black, alpha=0.1)

	hlines!(obs_props["sigma_v"], color=:black)

	scatter!(0, obs_props["sigma_v"], color=:black)
	@savefig "v_circ_time"

	fig
end

# ╔═╡ 7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-1, 3e2)),
		xlabel=L"\epsilon", ylabel=L"dN/d\epsilon",
			 yticks=Makie.automatic)

	x = LilGuys.specific_energy(snap_i)
	bins, values, errs = LilGuys.histogram(x[x .> 0], 
		normalization=:pdf)
	
	scatterlines!(midpoints(bins), values, label="initial")

	x = LilGuys.specific_energy(snap_f)
	bins, values, errs = LilGuys.histogram(x[x .> 0], 
		normalization=:pdf)
	
	scatterlines!(midpoints(bins), values, label="final")

	@savefig "energy_distribution"
	fig
end

# ╔═╡ c5796d82-013b-4cdc-a625-31249b51197d
md"""
# 2D plots
"""

# ╔═╡ 4801ff80-5761-490a-801a-b263b90d63fd
let
	fig, ax = FigAxis(aspect=1)
	ax.title = "initial"

	bins = 100
	colorrange=(1e-7, 1e-3)
	r_max = 5

	LilGuys.projecteddensity!(snap_i, centre=true, r_max=r_max,
		colorrange=colorrange, colorscale=log10,
		direction1=2, direction2=3,
		bins=bins
	)
	
	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc",
	title="final")

	hm = LilGuys.projecteddensity!(snap_f, centre=true, r_max=r_max,  
		colorrange=colorrange, colorscale=log10,
		direction1=2, direction2=3,
		bins=bins

	)
	
	Colorbar(fig[:, end+1], hm, label="DM density", ticks=Makie.automatic)
	
    rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])

	resize_to_layout!(fig)

	@savefig "xy_cen_projection"

	fig
end

# ╔═╡ 4cd952f3-555d-401b-aa31-8b79a23ca42e
let 
	fig = Figure()
	
	ax =Axis(fig[1, 1], aspect=1, 
		xlabel = "x / kpc", ylabel="z / kpc", title="dark matter",
	)

	colorrange=(1e-6, nothing)

	h = LilGuys.projecteddensity!(snap_f, centre=false, r_max=250, 
		colorrange=colorrange, colorscale=log10,
		xdirection=2, ydirection=3
	)


	# scatter!(snap_f.x_cen[2], snap_f.x_cen[3], markersize=0.3)
	Colorbar(fig[1, 2], h, label="DM density", ticks=Makie.automatic)

	@savefig "xz_fin_projection"

	fig
end

# ╔═╡ Cell order:
# ╠═bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
# ╠═d3bd61d8-1e90-4787-b892-d90717f6be6e
# ╟─99f96d71-b543-4680-a022-2195e6cca897
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═3e09298c-86c8-49b6-aaa8-516b9f199478
# ╠═1b87d662-da3c-4438-98eb-72dc93e32f6a
# ╠═2b2a1cc7-d005-4bef-b2cf-5d26b8c203a0
# ╟─9c4d9492-64bc-4212-a99d-67cc507e99e0
# ╠═82621ca0-2e0b-455a-b692-ab9a27c2a589
# ╠═3db38875-fe22-4cfd-8c5a-47f4a0fa7f3a
# ╠═c260ee35-7eed-43f4-b07a-df4371397195
# ╠═d010a230-7331-4afd-86dc-380da0e0f720
# ╠═510706ac-ffbd-4996-af9e-67f1b910d51c
# ╠═50bb0dfa-7332-484f-a27d-d8491413ef1e
# ╠═7094bc54-deb4-48a5-bf09-9ee6c684ac3c
# ╠═2470e05f-9215-45e4-88fc-daab0638272f
# ╠═b0e336df-678a-4406-b294-0c353f3c0c38
# ╟─ef37f7b9-15a2-412f-a254-6a0efac05a04
# ╠═5115a3fc-e71d-4806-a34a-aaa3b00c2d4e
# ╠═0526f791-9aff-4ba1-912c-97f74759e72b
# ╠═f835ba02-3e4b-4816-8480-452bfa308f5d
# ╠═a098c336-fda6-4932-a0ef-aa0f139a4499
# ╠═092460b1-ea96-4609-be91-3625d0e5f459
# ╟─a9e79439-16a4-4908-bfe0-f0770cdb26df
# ╠═53641449-c5b3-45ff-a692-a5cd717c8369
# ╠═9429b4e0-96d0-4fae-a7d6-44737d568f76
# ╠═7e3df305-9678-447e-a48e-f102cf6ebced
# ╠═9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
# ╠═8d127679-401c-439d-913d-e2020df1c600
# ╠═4977303f-b958-4d24-9a04-0f2835137d37
# ╠═485eab53-43ec-4591-a3af-9e4cbfefbbe2
# ╠═f3b1fd8e-0591-4d51-94ea-2b4eb65b6a71
# ╟─e14fa4a1-6175-4b9f-ad01-525c1617fe63
# ╠═8dae2e01-652b-4afc-b040-dd2ba1c6eedb
# ╠═b64c1caf-9ee0-4633-bd52-0258557b8847
# ╠═e6fc3297-c1b7-40a4-b2bb-98490a42604a
# ╠═d57501a1-4764-4b23-962f-2d37547d7bcc
# ╠═04ca92d1-f64b-4d2c-a079-30f89866fda9
# ╠═db320665-f46d-4aed-a2b2-4b39bcb605c5
# ╠═6f7eaf57-4ada-4995-9b88-52c984573b13
# ╠═c73ebaa4-16d9-428c-9967-6f8e67cdc680
# ╠═54e28676-c63d-4802-9a6f-4706654f3972
# ╠═1630a217-a6c5-439b-a79b-10b9a098ea29
# ╠═cabdac5e-9a2c-47f6-ba7f-c408b3af0f3d
# ╠═90f916db-3390-460e-8639-d5d02603aede
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
# ╟─c5796d82-013b-4cdc-a625-31249b51197d
# ╠═4801ff80-5761-490a-801a-b263b90d63fd
# ╟─4cd952f3-555d-401b-aa31-8b79a23ca42e
