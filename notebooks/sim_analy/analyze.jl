### A Pluto.jl notebook ###
# v0.20.19

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
end

# ╔═╡ ba60baa5-4213-4b24-b7f0-3eb647c5d311
using PyFITS

# ╔═╡ 987c3284-5a8f-463e-9c68-9011b348e076
using PlutoUI

# ╔═╡ bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
md"""
Analyzes the dark matter particles and profiles for the simulation.
In particular, we plot the initial/final density profiles, 
circular velocity evolution, 
and make some nice projections of the DM in different orientations.
All of the figures are saved to figures directory inside the model analysis directory. 
"""

# ╔═╡ 99f96d71-b543-4680-a022-2195e6cca897
md"""
# Setup
"""

# ╔═╡ 2b2a1cc7-d005-4bef-b2cf-5d26b8c203a0
CairoMakie.activate!(type=:png)

# ╔═╡ a16571bd-d231-4d70-a077-3983ea847833
import LinearAlgebra: cross

# ╔═╡ 5f2e646b-a7aa-453a-8afd-30b81ef07ff3
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ d3bd61d8-1e90-4787-b892-d90717f6be6e
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="sculptor"),
	modelname = TextField(60, default="1e6_V31_r3.2/orbit_"),
))

# ╔═╡ 9c4d9492-64bc-4212-a99d-67cc507e99e0
md"""
Inputs
"""

# ╔═╡ 3db38875-fe22-4cfd-8c5a-47f4a0fa7f3a
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", inputs.galaxyname, inputs.modelname)

# ╔═╡ 1b87d662-da3c-4438-98eb-72dc93e32f6a
FIGDIR = joinpath(model_dir, "figures")

# ╔═╡ c260ee35-7eed-43f4-b07a-df4371397195
readdir(model_dir)

# ╔═╡ d010a230-7331-4afd-86dc-380da0e0f720
halo = LilGuys.load_profile(joinpath(model_dir, "../halo.toml"))

# ╔═╡ 7094bc54-deb4-48a5-bf09-9ee6c684ac3c
out =  Output(model_dir)

# ╔═╡ 510706ac-ffbd-4996-af9e-67f1b910d51c
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 2470e05f-9215-45e4-88fc-daab0638272f
begin 
	profiles = LilGuys.read_ordered_structs(joinpath(model_dir, "profiles.hdf5"), LilGuys.MassProfile)
	snap_idx = first.(profiles)
	profiles = last.(profiles);
end

# ╔═╡ b0e336df-678a-4406-b294-0c353f3c0c38
dens_profiles = LilGuys.read_ordered_structs(joinpath(model_dir, "profiles_densities.hdf5"), LilGuys.DensityProfile) .|> last

# ╔═╡ 50bb0dfa-7332-484f-a27d-d8491413ef1e
df_scalars = read_fits(joinpath(model_dir, "profiles_scalars.fits"))

# ╔═╡ a9e79439-16a4-4908-bfe0-f0770cdb26df
md"""
# Mass evolution
"""

# ╔═╡ 53641449-c5b3-45ff-a692-a5cd717c8369
idx_f = min(orbit_props["idx_f"], length(profiles))

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

# ╔═╡ 8dae2e01-652b-4afc-b040-dd2ba1c6eedb
prof_i = profiles[1]

# ╔═╡ b64c1caf-9ee0-4633-bd52-0258557b8847
prof_f = profiles[end]

# ╔═╡ 4977303f-b958-4d24-9a04-0f2835137d37
times = out.times * T2GYR

# ╔═╡ 193e78d7-e902-4b6b-8eae-43bae7e5722b
df_scalars.time

# ╔═╡ f3b1fd8e-0591-4d51-94ea-2b4eb65b6a71
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "bound mass",
			 yscale=log10, yticks=Makie.automatic)


	lines!(df_scalars.time * T2GYR .- df_scalars.time[idx_f] * T2GYR, df_scalars.bound_mass ./ df_scalars.bound_mass[1])
	fig
end

# ╔═╡ 485eab53-43ec-4591-a3af-9e4cbfefbbe2
df_scalars.bound_mass[end] ./ df_scalars.bound_mass[1]

# ╔═╡ 4c042ef0-21a9-4a54-9562-05dc891f1dbf
function enclosed_mass(profile, r)
	f = LilGuys.lerp(profile.radii, middle.(profile.M_in))
	return f(r)
end

# ╔═╡ 04a07b4f-b747-4738-9d0f-18eae2f85baf
Rvir = LilGuys.R200(halo)

# ╔═╡ 0fa11815-3ab0-4b19-9be7-186b7c2c1063
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel=L"Bound mass within $r$ kpc",
		yscale=log10,
		#limits=(nothing, (1e-3, 2)),
		yticks=Makie.automatic
	)

	#ax.yticks = [0.1:0.1:1;]

	for r in [37, 10, 3, 1, 0.3, 0.1, 0.03]
		M_dm_h = enclosed_mass.(profiles, [r])
		@info M_dm_h[end] ./ M_dm_h[1]
		lines!(df_scalars.time, M_dm_h ./ M_dm_h[1], label="$r")
	end
	
	Legend(fig[1, 2], ax, "r/kpc")
	
	@savefig "boundmass"

	fig
end

# ╔═╡ ef2274ba-0220-49f5-8478-d5eac824c3ee


# ╔═╡ e14fa4a1-6175-4b9f-ad01-525c1617fe63
md"""

# Evolution of circular Velocity
"""

# ╔═╡ e6fc3297-c1b7-40a4-b2bb-98490a42604a
v_max = df_scalars.v_circ_max

# ╔═╡ d57501a1-4764-4b23-962f-2d37547d7bcc
r_max = df_scalars.r_circ_max

# ╔═╡ 04ca92d1-f64b-4d2c-a079-30f89866fda9
prof_i

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

# ╔═╡ c068c177-e879-4b8e-b1af-18690af9b334
let 
	i = idx_f - 30
	snap = out[i]
	prof = profiles[argmin(abs.(snap_idx .- i))]

	
	fig = Figure(size=(7*72, 5*72))
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"V_\textrm{circ}",
	title="time = $(out.times[i]  * T2GYR) Gyr")

	prof = MassProfile(snap)
	v = LilGuys.circular_velocity(prof)
	scatter!(log10.(prof.radii), middle.(v) * V2KMS)

	fit = LilGuys.fit_v_r_circ_max(prof.radii, middle.(v))
	r_fit = LinRange(fit.r_min, fit.r_max, 100)
	v_fit = LilGuys._v_circ_max_model(r_fit, [fit.r_circ_max, fit.v_circ_max])

	lines!(log10.(r_fit), v_fit * V2KMS, linewidth=3, color=COLORS[2])

	r_fit = 10 .^ LinRange(log10(prof.radii[1]), log10(prof.radii[end]), 1000)
	v_fit = LilGuys._v_circ_max_model(r_fit, [fit.r_circ_max, fit.v_circ_max])
	lines!(log10.(r_fit), v_fit * V2KMS, color=COLORS[2])
	
	scatter!(log10.(fit.r_circ_max), fit.v_circ_max * V2KMS, color=Arya.COLORS[3])


	ax2 = Axis(fig[2, 1], xlabel="r/kpc", ylabel="residual")

	filt = fit.r_min .<= prof.radii .<= fit.r_max
	dv = v[filt] .- LilGuys._v_circ_max_model(prof.radii[filt], [fit.r_circ_max, fit.v_circ_max])
	scatter!(log10.(prof.radii[filt]), dv * V2KMS)
	vlines!(log10(fit.r_circ_max))
	
    rowsize!(fig.layout, 2, Relative(1/4))

	linkxaxes!(ax, ax2)
	
	fig
end

# ╔═╡ 245721a6-01aa-43e7-922d-ed5da02207c1
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel=L"v_\text{circ} / \text{km\,s^{-1}}",
	limits=(nothing, (0, nothing)))
	x = out.times[snap_idx]
	lines!(x*T2GYR, v_max * V2KMS, label=L"maximum $v_\text{circ}$")
	#scatter!(x, v_h, label=L"r=r_h")
	#axislegend(ax)

	@savefig "v_circ_time"

	fig
end

# ╔═╡ c5796d82-013b-4cdc-a625-31249b51197d
md"""
# Density evolution
"""

# ╔═╡ 3aa62ecf-495a-434b-8008-02783bd5b56e
prof_f_allpart = LilGuys.DensityProfile(snap_f, filt_bound=false)

# ╔═╡ 79d8fbe1-09fb-4ed7-b75d-a9ac3a0fa9c5
dens_prof_i = dens_profiles[idx_i]

# ╔═╡ 71a9e730-68fd-4047-a647-54a95a448c01
dens_prof_f = dens_profiles[idx_f]

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=((-2, 2.5), (-10, 0)))

	lines!(dens_prof_i.log_r, log10.(dens_prof_i.rho), label="initial")
	lines!(dens_prof_f.log_r, log10.(dens_prof_f.rho), label="final")

	
	lines!(prof_f_allpart.log_r, log10.(prof_f_allpart.rho), label="final (all particles)", color=COLORS[2], linestyle=:solid)

	#LP.plot_ρ!(halo, linestyle=:dot, label="NFW", color=:black)

	axislegend(ax, position=:lb)
	#lines!([0, 0] .+ log10.(r_break), [-11.5, -10],  color=:black)
	#scatter!(log10.(r_break), -10, marker=:utriangle, color=:black)

	#text!(L"r_\textrm{break}", position=(log10.(r_break),-11.5), space=:data, rotation=π/2, align=(:left, :baseline))
	@savefig "dm_density"

	fig
	# only include bound points in profile...
end

# ╔═╡ 4a8fb43f-e2b1-4fa3-a99c-aa6fff4b727f
dens_profiles[1].time

# ╔═╡ 871f7679-dbaa-4901-85ab-357b58588d46
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=((-2, 2.5), (-10, 0)))

	colorrange = (prof_i.time, prof_f.time) .* T2GYR

	for (i, prof) in enumerate(dens_profiles[1:1:end])
		lines!(prof.log_r, log10.(middle.(prof.rho)), color=profiles[i].time * T2GYR, colorrange=colorrange)
	end
	
	Colorbar(fig[1,2], colorrange=colorrange, label="time / Gyr")

	@savefig "density_all_snapshots"

	fig
	# only include bound points in profile...
end

# ╔═╡ 15eecc05-c7fb-4424-94ac-819450e4f6ec
df_scalars.bound_mass[idx_f]

# ╔═╡ 48e54b34-4b22-4609-8928-ba6d8d027370
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"$\log\ r_\text{circ}$ / kpc", ylabel=L"\log\ v_\text{circ}\ /\ \text{km\,s^{-1}}",
		limits=(-2, 2, 0.8, 1.6)
	)

	colorrange = (prof_i.time, prof_f.time) .* T2GYR

	for prof in profiles[1:1:end]
		lines!(log10.(prof.radii), log10.(LilGuys.circular_velocity(prof) * V2KMS), color=prof.time * T2GYR, colorrange=colorrange)
	end
	
	Colorbar(fig[1,2], colorrange=colorrange, label="time / Gyr")

	@savefig "vcirc_rcirc_all_snapshots"

	fig
	# only include bound points in profile...
end

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

# ╔═╡ f7f8ed80-c715-43db-bebe-e62b14173ac6


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

# ╔═╡ 7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-1, 1e2)),
		xlabel=L"\epsilon", ylabel=L"dN/d\epsilon",
			 yticks=Makie.automatic)

	x = LilGuys.specific_energy(snap_i)
	bins, values, errs = LilGuys.histogram(x[x .> 0], 
		normalization=:pdf)
	
	scatter!(midpoints(bins), values, label="initial")

	x = LilGuys.specific_energy(snap_f)
	bins, values, errs = LilGuys.histogram(x[x .> 0], 
		normalization=:pdf)
	
	scatter!(midpoints(bins), values, label="final")

	@savefig "energy_distribution"
	fig
end

# ╔═╡ e385f8cf-faa1-49bf-88c9-15c3d2489f90
if length(out[1].index) < 1e7 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\beta",
	limits=((-2, 2.5), (-2, 1))
	)

	ax2 = Axis(fig[2,1], xlabel="log r / kpc", ylabel=L"\sigma",
	limits=((-2, 2.5), (0, 40))
	)

	for i in 1:30:length(out)
		snap = out[i]
		bins = LilGuys.quantile(radii(snap), LinRange(0, 1, 100))
	
		σs, β = LilGuys.β_prof(snap, r_bins=bins)
		x = log10.(midpoints(bins))
		lines!(ax, x, β, color=i, colorrange=(1, length(out)))
		
		lines!(ax2, x, sqrt.(σs) .* V2KMS, color=i, colorrange=(1, length(out)))
	end


	fig
	# only include bound points in profile...
end

# ╔═╡ 5af581bc-613c-4728-9626-dbef0ebaef7d
let 
	snap_test = deepcopy(snap_i)
	v_test = LilGuys.speeds(snap_test)

	bins = LilGuys.quantile(radii(snap_test), LinRange(0, 1, 30))

	println(1 .- LilGuys.β_prof(snap_test, r_bins=bins)[2])

	snap_test.positions .-= snap_test.x_cen
	snap_test.velocities .-= snap_test.v_cen
	snap_test.x_cen .= zeros(3)
	snap_test.v_cen .= zeros(3)

	println(1 .- LilGuys.β_prof(snap_test, r_bins=bins)[2])
end

# ╔═╡ Cell order:
# ╟─bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
# ╠═d3bd61d8-1e90-4787-b892-d90717f6be6e
# ╟─99f96d71-b543-4680-a022-2195e6cca897
# ╠═ba60baa5-4213-4b24-b7f0-3eb647c5d311
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═2b2a1cc7-d005-4bef-b2cf-5d26b8c203a0
# ╠═987c3284-5a8f-463e-9c68-9011b348e076
# ╠═a16571bd-d231-4d70-a077-3983ea847833
# ╠═1b87d662-da3c-4438-98eb-72dc93e32f6a
# ╠═5f2e646b-a7aa-453a-8afd-30b81ef07ff3
# ╟─9c4d9492-64bc-4212-a99d-67cc507e99e0
# ╠═3db38875-fe22-4cfd-8c5a-47f4a0fa7f3a
# ╠═c260ee35-7eed-43f4-b07a-df4371397195
# ╠═d010a230-7331-4afd-86dc-380da0e0f720
# ╠═7094bc54-deb4-48a5-bf09-9ee6c684ac3c
# ╠═510706ac-ffbd-4996-af9e-67f1b910d51c
# ╠═2470e05f-9215-45e4-88fc-daab0638272f
# ╠═b0e336df-678a-4406-b294-0c353f3c0c38
# ╠═50bb0dfa-7332-484f-a27d-d8491413ef1e
# ╟─a9e79439-16a4-4908-bfe0-f0770cdb26df
# ╠═53641449-c5b3-45ff-a692-a5cd717c8369
# ╠═9429b4e0-96d0-4fae-a7d6-44737d568f76
# ╠═7e3df305-9678-447e-a48e-f102cf6ebced
# ╠═9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
# ╠═8d127679-401c-439d-913d-e2020df1c600
# ╠═8dae2e01-652b-4afc-b040-dd2ba1c6eedb
# ╠═b64c1caf-9ee0-4633-bd52-0258557b8847
# ╠═4977303f-b958-4d24-9a04-0f2835137d37
# ╠═193e78d7-e902-4b6b-8eae-43bae7e5722b
# ╠═f3b1fd8e-0591-4d51-94ea-2b4eb65b6a71
# ╠═485eab53-43ec-4591-a3af-9e4cbfefbbe2
# ╠═4c042ef0-21a9-4a54-9562-05dc891f1dbf
# ╠═04a07b4f-b747-4738-9d0f-18eae2f85baf
# ╠═0fa11815-3ab0-4b19-9be7-186b7c2c1063
# ╠═ef2274ba-0220-49f5-8478-d5eac824c3ee
# ╟─e14fa4a1-6175-4b9f-ad01-525c1617fe63
# ╠═e6fc3297-c1b7-40a4-b2bb-98490a42604a
# ╠═d57501a1-4764-4b23-962f-2d37547d7bcc
# ╠═04ca92d1-f64b-4d2c-a079-30f89866fda9
# ╠═db320665-f46d-4aed-a2b2-4b39bcb605c5
# ╠═c068c177-e879-4b8e-b1af-18690af9b334
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╟─c5796d82-013b-4cdc-a625-31249b51197d
# ╠═3aa62ecf-495a-434b-8008-02783bd5b56e
# ╠═79d8fbe1-09fb-4ed7-b75d-a9ac3a0fa9c5
# ╠═71a9e730-68fd-4047-a647-54a95a448c01
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═4a8fb43f-e2b1-4fa3-a99c-aa6fff4b727f
# ╠═871f7679-dbaa-4901-85ab-357b58588d46
# ╠═15eecc05-c7fb-4424-94ac-819450e4f6ec
# ╠═48e54b34-4b22-4609-8928-ba6d8d027370
# ╟─4801ff80-5761-490a-801a-b263b90d63fd
# ╠═f7f8ed80-c715-43db-bebe-e62b14173ac6
# ╠═4cd952f3-555d-401b-aa31-8b79a23ca42e
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
# ╠═e385f8cf-faa1-49bf-88c9-15c3d2489f90
# ╠═5af581bc-613c-4728-9626-dbef0ebaef7d
