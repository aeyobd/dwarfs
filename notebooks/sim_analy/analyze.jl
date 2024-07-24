### A Pluto.jl notebook ###
# v0.19.45

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

# ╔═╡ d7aaba0e-1ba9-4349-b2f0-c047bb49bcd7
using KernelDensity

# ╔═╡ bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
md"""
Analyzes the dark matter particles and profiles for the simulation.
"""

# ╔═╡ 589f2e18-d697-4fa1-a1bc-eb3c83ad73bf
save = CairoMakie.save

# ╔═╡ 7f29d14d-0e3f-4e64-b726-6b8fd5bb7548
LP = LilGuys.Plots

# ╔═╡ 9c4d9492-64bc-4212-a99d-67cc507e99e0
md"""
Inputs
"""

# ╔═╡ 14279a79-bf66-4b34-bf9f-735ff2886ea5
model_dir = "/astro/dboyea/sculptor/orbits/1e6/orbit1/V70_r0.4"

# ╔═╡ d010a230-7331-4afd-86dc-380da0e0f720
halo = NFW(; LilGuys.dict_to_tuple(TOML.parsefile(joinpath(model_dir, "halo.toml"))["profile"])...)

# ╔═╡ d971556d-8b66-4b2d-9dc0-31799f94b10a
skip = 10

# ╔═╡ 7094bc54-deb4-48a5-bf09-9ee6c684ac3c
out =  Output(model_dir)

# ╔═╡ 1b87d662-da3c-4438-98eb-72dc93e32f6a
figures_dir = joinpath(model_dir, "figures")

# ╔═╡ 63b7c3a2-247e-41b3-8a52-b92fd7a3cffe
mkpath(figures_dir)

# ╔═╡ 510706ac-ffbd-4996-af9e-67f1b910d51c
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ a9e79439-16a4-4908-bfe0-f0770cdb26df
md"""
# Mass evolution
"""

# ╔═╡ 53641449-c5b3-45ff-a692-a5cd717c8369
idx_f = orbit_props["idx_f"]

# ╔═╡ 7e3df305-9678-447e-a48e-f102cf6ebced
idx_i = 1

# ╔═╡ 9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
snap_i = out[idx_i]

# ╔═╡ 8d127679-401c-439d-913d-e2020df1c600
snap_f = out[idx_f]

# ╔═╡ 8dae2e01-652b-4afc-b040-dd2ba1c6eedb
prof_i = LilGuys.calc_profile(snap_i)

# ╔═╡ b64c1caf-9ee0-4633-bd52-0258557b8847
prof_f = LilGuys.calc_profile(snap_f)

# ╔═╡ 2470e05f-9215-45e4-88fc-daab0638272f
profiles = [LilGuys.calc_profile(snap) for snap in out[1:skip:end]]

# ╔═╡ e8956092-e811-4af8-bcb0-2fc9829ca817
LilGuys.calc_profile(out[140])

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

	for r in [0.3, 1, 10, Inf]
		M_dm_h = LilGuys.get_M_h(out, r)
		scatter!(times[1:10:end], M_dm_h ./ M_dm_h[1], label="$r")
	end
	
	Legend(fig[1, 2], ax)
	
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
		xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"$v_\textrm{circ}$ / km s$^{-1}$",
		yscale=log10,
		yticks=[1, 10, 20, 30, 40, 50, 60],
		yminorticks=[1:9; 10:2:60],
		limits=(nothing, (1, 75)),
		xgridvisible=false,
		ygridvisible=false
	)

	r_model = 10 .^ LinRange(-2, 2, 1000)
	v_model = calc_v_circ.(halo, r_model)
	lines!(log10.(r_model), v_model * V2KMS, linestyle=:dot, label="NFW")
	
	lines!(prof_i.log_r, prof_i.v_circ * V2KMS, label="initial")

	lines!(prof_f.log_r, prof_f.v_circ * V2KMS, label="final")


	α = 0.4
	β = 0.65
	x = LinRange(1, 0.1, 100)

	y = @. 2^α * x^β * (1 + x^2)^(-α)
	lines!(log10.(x .* r_max[1]), y .* v_max[1] * V2KMS,  label="EN21",
	color=:black, linestyle=:dash)

	lines!(log10.(r_max), v_max * V2KMS, color=Arya.COLORS[4], label=L"v_\textrm{circ,\ max}")

		
	axislegend(ax, position=:rb)
	save(joinpath(figures_dir, "v_circ_profiles.pdf"), fig)
	fig
end

# ╔═╡ c068c177-e879-4b8e-b1af-18690af9b334
let 
	i = length(out) - 80
	snap = out[i]

	
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"V_\textrm{circ}",
	title="time = $(out.times[i]  * T2GYR) Gyr")

	r, v = LilGuys.calc_v_circ(snap)
	scatter!(log10.(r), v * V2KMS)

	fit = LilGuys.fit_v_r_circ_max(snap)
	r_fit = LinRange(fit.r_min, fit.r_max, 100)
	v_fit = LilGuys.v_circ_max_model(r_fit, [fit.r_circ_max, fit.v_circ_max])

	lines!(log10.(r_fit), v_fit * V2KMS, linewidth=3)
	
	scatter!(log10.(fit.r_circ_max), fit.v_circ_max * V2KMS, color=Arya.COLORS[3])


	ax2 = Axis(fig[2, 1], xlabel="r/kpc", ylabel="residual")

	filt = fit.r_min .<= r .<= fit.r_max
	dv = v[filt] .- LilGuys.v_circ_max_model(r[filt], [fit.r_circ_max, fit.v_circ_max])
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
	x = out.times[1:skip:end] * T2GYR
	scatter!(x, v_max * V2KMS, label=L"maximum $v_\text{circ}$")
	#scatter!(x, v_h, label=L"r=r_h")
	axislegend(ax)

	save(joinpath(figures_dir, "v_circ_time.pdf"), fig)

	fig
end

# ╔═╡ c5796d82-013b-4cdc-a625-31249b51197d
md"""
# Density evolution
"""

# ╔═╡ 6cc4868a-8ef6-4d29-9d7d-f04504d6b157
snap_i.masses

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=(nothing, (-13, 0)))

	lines!(prof_i.log_r, log10.(prof_i.rho), label="initial")
	lines!(prof_f.log_r, log10.(prof_f.rho), label="final")

	LP.plot_ρ_dm!(snap_f, label="final (all particles)", filt_bound=false, color=COLORS[2], linestyle=:dash)

	LP.plot_ρ!(halo, linestyle=:dot, label="NFW")

	axislegend(ax)
	#lines!([0, 0] .+ log10.(r_break), [-11.5, -10],  color=:black)
	#scatter!(log10.(r_break), -10, marker=:utriangle, color=:black)

	#text!(L"r_\textrm{break}", position=(log10.(r_break),-11.5), space=:data, rotation=π/2, align=(:left, :baseline))
	save(joinpath(figures_dir, "density.pdf"), fig)

	fig
	# only include bound points in profile...
end

# ╔═╡ 4cd952f3-555d-401b-aa31-8b79a23ca42e
let 
	fig, ax = LP.xy_axis()

	bins = LinRange(-10, 10, 100)
	colorrange=(1e-6, nothing)

	LP.projected_density!(snap_f, centre=false, r_max=130, 
		colorrange=colorrange, colorscale=log10,
		direction1=2, direction2=3
	)
	fig
end

# ╔═╡ 8c0c017d-647d-49ed-95ce-1bc85725ca79


# ╔═╡ 4801ff80-5761-490a-801a-b263b90d63fd
let
	fig, ax = LP.xy_axis()
	ax.title = "initial"

	bins = LinRange(-10, 10, 100)
	colorrange=(1e-3, 1e3)

	LP.projected_density!(out[1], centre=true)
	
	bins = (out.x_cen[1, idx_i]  .+ bins,  out.x_cen[2, idx_i]  .+bins)
	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc",
	title="final")

	bins = LinRange(-10, 10, 100)
	bins = (out.x_cen[1, idx_f]  .+ bins,  out.x_cen[2, idx_f]  .+bins)
	hm = Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], bins = bins, colorscale=log10, colorrange=colorrange)
	
	Colorbar(fig[:, end+1], hm, label="DM density")
	
    rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])

	resize_to_layout!(fig)

	save(joinpath(figures_dir, "xy_cen_projection.pdf"), fig)

	fig
end

# ╔═╡ fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
let
	fig = Figure()
	r_max = 130
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc", title="dark matter",
	limits=(-r_max, r_max, -r_max, r_max))
	bins = LinRange(-r_max, r_max, 200)
	
	h = Arya.hist2d!(ax, snap_f.positions[2, :], snap_f.positions[3, :], bins = bins, colorscale=log10, colorrange=(1e-1, nothing))

	#scatter!(snap_f.x_cen[2], snap_f.x_cen[3])
	Colorbar(fig[1, 2], h, label="DM density")

	save(joinpath(figures_dir, "xy_projection.pdf"), fig)
	fig
end

# ╔═╡ 7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e2, 1e6)),
		xlabel=L"\epsilon", ylabel="count")
	stephist!(LilGuys.calc_ϵ(snap_i), )
	es = LilGuys.calc_ϵ(snap_f)
	es = es[es .> 0]
	stephist!(es)

	save(joinpath(figures_dir, "energy distribution.pdf"), fig)
	fig
end

# ╔═╡ Cell order:
# ╟─bafc8bef-6646-4b2f-9ac0-2ac09fbcb8e1
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═589f2e18-d697-4fa1-a1bc-eb3c83ad73bf
# ╠═7f29d14d-0e3f-4e64-b726-6b8fd5bb7548
# ╟─9c4d9492-64bc-4212-a99d-67cc507e99e0
# ╠═14279a79-bf66-4b34-bf9f-735ff2886ea5
# ╠═d010a230-7331-4afd-86dc-380da0e0f720
# ╠═d971556d-8b66-4b2d-9dc0-31799f94b10a
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
# ╠═e8956092-e811-4af8-bcb0-2fc9829ca817
# ╠═4977303f-b958-4d24-9a04-0f2835137d37
# ╠═0fa11815-3ab0-4b19-9be7-186b7c2c1063
# ╟─e14fa4a1-6175-4b9f-ad01-525c1617fe63
# ╠═e6fc3297-c1b7-40a4-b2bb-98490a42604a
# ╠═d57501a1-4764-4b23-962f-2d37547d7bcc
# ╠═db320665-f46d-4aed-a2b2-4b39bcb605c5
# ╟─c068c177-e879-4b8e-b1af-18690af9b334
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╟─c5796d82-013b-4cdc-a625-31249b51197d
# ╠═6cc4868a-8ef6-4d29-9d7d-f04504d6b157
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═4cd952f3-555d-401b-aa31-8b79a23ca42e
# ╠═8c0c017d-647d-49ed-95ce-1bc85725ca79
# ╠═4801ff80-5761-490a-801a-b263b90d63fd
# ╠═fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
# ╠═d7aaba0e-1ba9-4349-b2f0-c047bb49bcd7
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
