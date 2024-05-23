### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 340ffbbe-17bd-11ef-35c6-63505bb128b7
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ cf6a7cbb-1034-4026-a3f3-1e854d2929e2
using FITSIO, Tables

# ╔═╡ 8115808d-b2d5-4207-93e0-e3e4afd64c90
Arya.hist2d(randn(100), randn(100))

# ╔═╡ c17b6fb0-3fdc-437e-84f8-cb0c81764898
import YAML

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = YAML.load_file("obital_properties.yml")

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
starsfile = "../../isolation/1e6/stars/exp2d_small_stars.hdf5"

# ╔═╡ f9fe37ef-de81-4d69-9308-cda968851ed2
begin 
	using HDF5

	f = h5open(starsfile)
	p_idx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
	close(f)
	
end

# ╔═╡ ef3481f8-2505-4c04-a04a-29bdf34e9abb
outfile = splitext(basename(starsfile))[1] * "_today.fits"

# ╔═╡ 172588cc-ae22-440e-8488-f508aaf7ce96
rel_p_cut = 1e-5

# ╔═╡ 973955ad-3214-42cf-831f-a782f1d2434a
idx_i = 1 

# ╔═╡ 8f9ee2a7-de43-43e0-8257-93ecc630044f
idx_f = orbit_props["idx_f"]

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
begin 
	out =  lguys.Output("out/combined.hdf5")
	
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
	out.x_cen .= x_cen
	out.v_cen .= v_cen

	out
end

# ╔═╡ 6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
snap_f = out[idx_f]

# ╔═╡ 396a53a3-de0f-4d97-9693-40f3757d66f9
snap_i = out[idx_i]

# ╔═╡ 5edf4585-a842-4898-b018-e36255f6fee9
supertypes(NTuple)

# ╔═╡ 77479cd4-513c-4603-9aa0-1acd964c403a
let
	fig = Figure()
	r_max = 100
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc", title="stars",
	limits=(-r_max, r_max, -r_max, r_max))
	bins = LinRange(-r_max, r_max, 300)

	hm = Arya.hist2d!(ax, snap_f.positions[2, :], snap_f.positions[3, :], bins = bins, colorscale=log10, colorrange=(1e-10, 0.06), weights=probabilities[snap_f.index])

	Colorbar(fig[:, end+1], hm)
	fig
end

# ╔═╡ e33245ad-7562-47af-94b9-9708af45b2a4
mlog10 = Makie.pseudolog10

# ╔═╡ 7bc2c15c-33f7-43f3-a47f-ca39ffc22071
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bin_range = LinRange(-2, 2, 100)
	colorrange =(1e-10, 1e-2)

	bins = (x_cen[1, idx_i]  .+ bin_range,  x_cen[2, idx_i]  .+ bin_range)
		
	probs = probabilities[snap_i.index]

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=probs, bins = bins, colorscale=mlog10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")

	probs = probabilities[snap_f.index]

	bins = (x_cen[1, idx_f]  .+ bin_range,  x_cen[2, idx_f]  .+ bin_range)

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=probs, bins = bins, colorscale=mlog10, colorrange=colorrange)
	
	fig
end

# ╔═╡ 8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-8, 1)),
		xlabel=L"\epsilon", ylabel="count")
	stephist!(lguys.calc_ϵ(snap_i), weights=probabilities[snap_i.index])
	es = lguys.calc_ϵ(snap_f)
	filt = es .> 0
	es = es[filt]
	probs = probabilities[snap_f.index][filt]
	stephist!(es, weights = probs)

	fig

end

# ╔═╡ c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
md"""
# 3D density profiles
"""

# ╔═╡ a1138e64-5fa5-4a0e-aeef-487ee78a7adc
function plot_ρ_s!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions, p_idx)
	rs = lguys.calc_r(pos, snap.x_cen)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=probabilities)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 1866c280-89c3-4a71-9dbf-50b583360145
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel = L"\log \rho_\star", 
		limits=((-1.9, 1), (-8, 2)))

	plot_ρ_s!(snap_i, label="initial")
	plot_ρ_s!(snap_f, label="final")
	
	axislegend(ax)
	fig
end

# ╔═╡ f5a3ea2f-4f19-4b0e-af55-74902f2c6485
md"""
# Sky projection
"""

# ╔═╡ 1f722acb-f7b9-4d6c-900e-11eae85e0708
begin 
	p_min = rel_p_cut * maximum(probabilities)
	snap_f_stars = snap_f[probabilities[snap_f.index] .> p_min]
	println("number of final stars: ", length(snap_f_stars))
	obs_pred = lguys.to_sky(snap_f_stars)
	m_star_f = probabilities[snap_f_stars.index]
end

# ╔═╡ fe56d883-1e4b-400c-b5f2-03ac058bb0ad
begin
	obs_df = DataFrame(
		Dict(
			String(key) => [getproperty(o, key) for o in obs_pred]
			for key in [:ra, :dec, :pm_ra, :pm_dec, :distance, :radial_velocity]
		)
	)

	obs_df[!, "probability"] = m_star_f

	
end

# ╔═╡ ead40db8-ab7c-4167-b2ca-b00b5cb3cde4
df = Dict(String(name) => obs_df[:, name] for name in names(obs_df))

# ╔═╡ 7e588ae3-89f3-4b91-8963-f6bf4391a859
FITS(outfile, "w") do f
	write(f, df)
end

# ╔═╡ 2612e3a2-6e6e-494e-b140-720dd2db6ec2
obs = YAML.load_file("../../mc_orbits/orbit1.yml")

# ╔═╡ 8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
let 
	fig = Figure()

	dy = 2
	dx = dy * 1/cosd(obs["dec"])
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=(obs["ra"] .+ (-dx, dx), obs["dec"] .+ (-dy, dy)),
		aspect = 1/cosd(obs["dec"])
	)

	x = [o.ra for o in obs_pred]
	y = [o.dec for o in obs_pred]

	
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)
	errscatter!([obs["ra"]], [obs["dec"]], color=:red)
	
	fig
end

# ╔═╡ 76d3b8ee-94f3-4d8e-adeb-b7f7346a3c5b
obs

# ╔═╡ 9d74ffbb-4c31-4062-9434-7755f53e4da0
let
	dr = 0.05
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\mu_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\mu_\delta / \textrm{mas\,yr^{-1}}",
	limits=(obs["pm_ra"] .+ (-dr, dr), obs["pm_dec"] .+ (-dr, dr)),
	aspect=DataAspect())

	x = [o.pm_ra for o in obs_pred]
	y = [o.pm_dec for o in obs_pred]
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)
	errscatter!([obs["pm_ra"]], [obs["pm_dec"]], xerr=[obs["pm_ra_err"]], yerr=[obs["pm_dec_err"]], color=:red)
	
	fig
end

# ╔═╡ 6fe6deb4-ae44-4ca0-9617-95841fdaf791
let
	fig = Figure()
	dx = 2
	dy = 20
	
	ax = Axis(fig[1,1],
		xlabel="distance / kpc",
		ylabel = "radial velocity / km/s",
		limits=(obs["distance"] .+ (-dx, dx), obs["radial_velocity"] .+ (-dy, dy))
	)


	
	x = [o.distance for o in obs_pred]
	y = [o.radial_velocity for o in obs_pred]
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)

	errscatter!([obs["distance"]], [obs["radial_velocity"]], xerr=[obs["distance"]], yerr=[obs["radial_velocity"]], color=:red)

	
	fig
end

# ╔═╡ 59202a69-fb7e-4b10-854d-5f37dbf1b1fa
let
	fig = Figure()
	ax = Axis(fig, limits=((1, 2), nothing))
	ax.limits.val
end

# ╔═╡ Cell order:
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═8115808d-b2d5-4207-93e0-e3e4afd64c90
# ╠═c17b6fb0-3fdc-437e-84f8-cb0c81764898
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═ef3481f8-2505-4c04-a04a-29bdf34e9abb
# ╠═172588cc-ae22-440e-8488-f508aaf7ce96
# ╠═973955ad-3214-42cf-831f-a782f1d2434a
# ╠═6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
# ╠═8f9ee2a7-de43-43e0-8257-93ecc630044f
# ╠═396a53a3-de0f-4d97-9693-40f3757d66f9
# ╠═f9fe37ef-de81-4d69-9308-cda968851ed2
# ╠═6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
# ╠═5edf4585-a842-4898-b018-e36255f6fee9
# ╠═77479cd4-513c-4603-9aa0-1acd964c403a
# ╠═e33245ad-7562-47af-94b9-9708af45b2a4
# ╠═7bc2c15c-33f7-43f3-a47f-ca39ffc22071
# ╠═8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
# ╟─c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
# ╠═a1138e64-5fa5-4a0e-aeef-487ee78a7adc
# ╠═1866c280-89c3-4a71-9dbf-50b583360145
# ╟─f5a3ea2f-4f19-4b0e-af55-74902f2c6485
# ╠═1f722acb-f7b9-4d6c-900e-11eae85e0708
# ╠═fe56d883-1e4b-400c-b5f2-03ac058bb0ad
# ╠═7e588ae3-89f3-4b91-8963-f6bf4391a859
# ╠═cf6a7cbb-1034-4026-a3f3-1e854d2929e2
# ╠═ead40db8-ab7c-4167-b2ca-b00b5cb3cde4
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
# ╠═76d3b8ee-94f3-4d8e-adeb-b7f7346a3c5b
# ╠═9d74ffbb-4c31-4062-9434-7755f53e4da0
# ╠═6fe6deb4-ae44-4ca0-9617-95841fdaf791
# ╠═59202a69-fb7e-4b10-854d-5f37dbf1b1fa
