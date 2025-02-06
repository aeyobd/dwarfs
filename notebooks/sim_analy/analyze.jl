### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import TOML

	using LilGuys
	using Arya
end

# ╔═╡ bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
md"""
Analyzes the dark matter particles and profiles for the simulation.
In particular, we plot the initial/final density profiles, 
circular velocity evolution, 
and make some nice projections of the DM in different orientations.
All of the figures are saved to figures directory inside the model analysis directory. 
"""

# ╔═╡ 9c4d9492-64bc-4212-a99d-67cc507e99e0
md"""
Inputs
"""

# ╔═╡ 14279a79-bf66-4b34-bf9f-735ff2886ea5
# model_dir = "/astro/dboyea/dwarfs/analysis/sculptor/1e7_V31_r3.2/orbit_smallperi"
model_dir = "/astro/dboyea/dwarfs/analysis/sculptor/1e7_V31_r4.2/vasiliev24_L3M11_2x_smallperilmc"

# ╔═╡ c260ee35-7eed-43f4-b07a-df4371397195
readdir(model_dir)

# ╔═╡ d010a230-7331-4afd-86dc-380da0e0f720
halo = LilGuys.load_profile(joinpath(model_dir, "../halo.toml"))

# ╔═╡ 7094bc54-deb4-48a5-bf09-9ee6c684ac3c
out =  Output(model_dir)

# ╔═╡ 1b87d662-da3c-4438-98eb-72dc93e32f6a
figdir = joinpath(model_dir, "figures")

# ╔═╡ 63b7c3a2-247e-41b3-8a52-b92fd7a3cffe
mkpath(figdir)

# ╔═╡ 510706ac-ffbd-4996-af9e-67f1b910d51c
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ a9e79439-16a4-4908-bfe0-f0770cdb26df
md"""
# Mass evolution
"""

# ╔═╡ 53641449-c5b3-45ff-a692-a5cd717c8369
idx_f = orbit_props["idx_f"]

# ╔═╡ 7e3df305-9678-447e-a48e-f102cf6ebced
idx_i = 2

# ╔═╡ 9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
snap_i = out[idx_i]

# ╔═╡ 8d127679-401c-439d-913d-e2020df1c600
snap_f = out[idx_f]

# ╔═╡ 2470e05f-9215-45e4-88fc-daab0638272f
begin 
	profiles = LilGuys.read_structs_from_hdf5(joinpath(model_dir, "profiles.hdf5"), LilGuys.MassProfile3D)

	snap_idx = parse.(Int, first.(profiles))

	profiles = last.(profiles)

	profiles = profiles[sortperm(snap_idx)]
	snap_idx = sort(snap_idx)
end

# ╔═╡ 8dae2e01-652b-4afc-b040-dd2ba1c6eedb
prof_i = profiles[1]

# ╔═╡ b64c1caf-9ee0-4633-bd52-0258557b8847
prof_f = profiles[end]

# ╔═╡ 4977303f-b958-4d24-9a04-0f2835137d37
times = out.times * T2GYR

# ╔═╡ 0fa11815-3ab0-4b19-9be7-186b7c2c1063
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel=L"Bound mass within $r$ kpc",
		yscale=log10,
		#limits=(nothing, (1e-3, 2)),
	)

	#ax.yticks = [0.1:0.1:1;]

	for r in [Inf, 10, 1]
		M_dm_h = LilGuys.calc_M_in(out, r)
		scatter!(times[1:10:end], M_dm_h ./ M_dm_h[1], label="$r")
	end
	
	Legend(fig[1, 2], ax)
	
	@savefig "boundmass"

	fig
end

# ╔═╡ e14fa4a1-6175-4b9f-ad01-525c1617fe63
md"""

# Evolution of circular Velocity
"""

# ╔═╡ e6fc3297-c1b7-40a4-b2bb-98490a42604a
v_max = [f.v_circ_max for f in profiles]

# ╔═╡ d57501a1-4764-4b23-962f-2d37547d7bcc
r_max = [f.r_circ_max for f in profiles]

# ╔═╡ db320665-f46d-4aed-a2b2-4b39bcb605c5
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log \; r_\textrm{circ}\ /\ \textrm{kpc}", 
		ylabel=L"$v_\textrm{circ}$ / km s$^{-1}$",
		yscale=log10,
		yticks=[1, 10, 20, 30, 40, 50, 60],
		yminorticks=[1:9; 10:2:60],
		limits=((-2, 3), (5, 50)),
		xgridvisible=false,
		ygridvisible=false
	)

	r_model = 10 .^ LinRange(-2, 3, 1000)
	v_model = calc_v_circ.(halo, r_model)
	lines!(log10.(r_model), v_model * V2KMS, linestyle=:dot, label="analytic")
	
	lines!(log10.(prof_i.r_circ), prof_i.v_circ * V2KMS, label="initial")

	lines!(log10.(prof_f.r_circ), prof_f.v_circ * V2KMS, label="final")


	α = 0.4
	β = 0.65
	x = LinRange(1, 0.1, 100)

	y = @. 2^α * x^β * (1 + x^2)^(-α)
	lines!(log10.(x .* r_max[1]), y .* v_max[1] * V2KMS,  label="EN21",
	color=:black, linestyle=:dash)

	scatter!(log10.(r_max), v_max * V2KMS, color=Arya.COLORS[4], label=L"v_\textrm{circ,\ max}")

		
	axislegend(ax, position=:rt)
	@savefig "v_circ_profiles"
	fig
end

# ╔═╡ c068c177-e879-4b8e-b1af-18690af9b334
let 
	i = length(out) - 10
	snap = out[i]
	prof = profiles[argmin(abs.(snap_idx .- i))]

	
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"V_\textrm{circ}",
	title="time = $(out.times[i]  * T2GYR) Gyr")

	r, v = LilGuys.calc_v_circ(snap)
	scatter!(log10.(r), v * V2KMS)

	fit = LilGuys.fit_v_r_circ_max(snap)
	r_fit = LinRange(fit.r_min, fit.r_max, 100)
	v_fit = LilGuys._v_circ_max_model(r_fit, [fit.r_circ_max, fit.v_circ_max])

	lines!(log10.(r_fit), v_fit * V2KMS, linewidth=3, color=COLORS[2])

	r_fit = 10 .^ LinRange(prof.log_r[1], prof.log_r[end], 1000)
	v_fit = LilGuys._v_circ_max_model(r_fit, [fit.r_circ_max, fit.v_circ_max])
	lines!(log10.(r_fit), v_fit * V2KMS, color=COLORS[2])
	
	scatter!(log10.(fit.r_circ_max), fit.v_circ_max * V2KMS, color=Arya.COLORS[3])


	ax2 = Axis(fig[2, 1], xlabel="r/kpc", ylabel="residual")

	filt = fit.r_min .<= r .<= fit.r_max
	dv = v[filt] .- LilGuys._v_circ_max_model(r[filt], [fit.r_circ_max, fit.v_circ_max])
	scatter!(log10.(r[filt]), dv * V2KMS)
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
prof_f_allpart = LilGuys.MassProfile3D(snap_f, filt_bound=false)

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=((-2, 2.5), (-10, 0)))

	lines!(prof_i.log_r, log10.(prof_i.rho), label="initial")
	lines!(prof_f.log_r, log10.(prof_f.rho), label="final")

	
	lines!(prof_f_allpart.log_r, log10.(prof_f_allpart.rho), label="final (all particles)", color=COLORS[2], linestyle=:dash)

	#LP.plot_ρ!(halo, linestyle=:dot, label="NFW", color=:black)

	axislegend(ax)
	#lines!([0, 0] .+ log10.(r_break), [-11.5, -10],  color=:black)
	#scatter!(log10.(r_break), -10, marker=:utriangle, color=:black)

	#text!(L"r_\textrm{break}", position=(log10.(r_break),-11.5), space=:data, rotation=π/2, align=(:left, :baseline))
	@savefig "dm_density"

	fig
	# only include bound points in profile...
end

# ╔═╡ 871f7679-dbaa-4901-85ab-357b58588d46
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=((-2, 2.5), (-10, 0)))

	colorrange = (prof_i.time, prof_f.time) .* T2GYR

	for prof in profiles[1:5:end]
		lines!(prof.log_r, log10.(prof.rho), color=prof.time * T2GYR, colorrange=colorrange)
	end
	
	Colorbar(fig[1,2], colorrange=colorrange, label="time / Gyr")

	@savefig "density_all_snapshots"

	fig
	# only include bound points in profile...
end

# ╔═╡ 48e54b34-4b22-4609-8928-ba6d8d027370
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"$\log\ r_\text{circ}$ / kpc", ylabel=L"\log\ v_\text{circ}\ /\ \text{km\,s^{-1}}",
		limits=(-2, 2, 0.8, 1.6)
	)

	colorrange = (prof_i.time, prof_f.time) .* T2GYR

	for prof in profiles[1:1:end]
		lines!(log10.(prof.r_circ), log10.(prof.v_circ * V2KMS), color=prof.time * T2GYR, colorrange=colorrange)
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
	
	Colorbar(fig[:, end+1], hm, label="DM density")
	
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

	Colorbar(fig[1, 2], h, label="DM density")

	@savefig "xz_fin_projection"

	fig
end

# ╔═╡ 7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-1, 1e2)),
		xlabel=L"\epsilon", ylabel=L"dN/d\epsilon")

	x = LilGuys.calc_ϵ(snap_i)
	bins, values, errs = LilGuys.histogram(x[x .> 0], 
		normalization=:pdf)
	
	scatter!(midpoints(bins), values, label="initial")

	x = LilGuys.calc_ϵ(snap_f)
	bins, values, errs = LilGuys.histogram(x[x .> 0], 
		normalization=:pdf)
	
	scatter!(midpoints(bins), values, label="final")

	@savefig "energy_distribution"
	fig
end

# ╔═╡ Cell order:
# ╠═bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╟─9c4d9492-64bc-4212-a99d-67cc507e99e0
# ╠═14279a79-bf66-4b34-bf9f-735ff2886ea5
# ╠═c260ee35-7eed-43f4-b07a-df4371397195
# ╠═d010a230-7331-4afd-86dc-380da0e0f720
# ╠═7094bc54-deb4-48a5-bf09-9ee6c684ac3c
# ╠═1b87d662-da3c-4438-98eb-72dc93e32f6a
# ╠═63b7c3a2-247e-41b3-8a52-b92fd7a3cffe
# ╠═510706ac-ffbd-4996-af9e-67f1b910d51c
# ╟─a9e79439-16a4-4908-bfe0-f0770cdb26df
# ╠═53641449-c5b3-45ff-a692-a5cd717c8369
# ╠═7e3df305-9678-447e-a48e-f102cf6ebced
# ╠═9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
# ╠═8d127679-401c-439d-913d-e2020df1c600
# ╠═8dae2e01-652b-4afc-b040-dd2ba1c6eedb
# ╠═b64c1caf-9ee0-4633-bd52-0258557b8847
# ╠═2470e05f-9215-45e4-88fc-daab0638272f
# ╠═4977303f-b958-4d24-9a04-0f2835137d37
# ╟─0fa11815-3ab0-4b19-9be7-186b7c2c1063
# ╟─e14fa4a1-6175-4b9f-ad01-525c1617fe63
# ╠═e6fc3297-c1b7-40a4-b2bb-98490a42604a
# ╠═d57501a1-4764-4b23-962f-2d37547d7bcc
# ╠═db320665-f46d-4aed-a2b2-4b39bcb605c5
# ╠═c068c177-e879-4b8e-b1af-18690af9b334
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╟─c5796d82-013b-4cdc-a625-31249b51197d
# ╠═3aa62ecf-495a-434b-8008-02783bd5b56e
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═871f7679-dbaa-4901-85ab-357b58588d46
# ╠═48e54b34-4b22-4609-8928-ba6d8d027370
# ╟─4801ff80-5761-490a-801a-b263b90d63fd
# ╠═f7f8ed80-c715-43db-bebe-e62b14173ac6
# ╟─4cd952f3-555d-401b-aa31-8b79a23ca42e
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
