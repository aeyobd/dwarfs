### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ 8ed2e715-1095-4c99-9ecb-54fecbc27d7f
pwd()

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = TOML.parsefile("orbital_properties.toml")

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
starsfile = "../../isolation/1e6/stars/king_stars.hdf5"

# ╔═╡ f9fe37ef-de81-4d69-9308-cda968851ed2
begin 
	using HDF5

	f = h5open(starsfile)
	p_idx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
	close(f)
	
end

# ╔═╡ 377284f2-dcee-44d3-9a04-728605cea92a
md"""
Given a stellar probability file, calculates initial-final density profiles, 
and projects stars onto the sky
"""

# ╔═╡ 5e12d306-a430-4b15-b3a7-d4806a5856cd
name = splitext(basename(starsfile))[1]

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
	colorrange =(1e-5, nothing)

	bins = (x_cen[1, idx_i]  .+ bin_range,  x_cen[2, idx_i]  .+ bin_range)
		
	probs = probabilities[snap_i.index]

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=probs, bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")

	probs = probabilities[snap_f.index]

	bins = (x_cen[1, idx_f]  .+ bin_range,  x_cen[2, idx_f]  .+ bin_range)

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=probs, bins = bins, colorscale=log10, colorrange=colorrange)
	
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

# ╔═╡ 2e03e2e9-86aa-419b-a1c8-b05da3abce61
let
	c = cens[idx_f, :]
	
	global obs_c_galcen = lguys.Galactocentric(x=c.x, y=c.y, z=c.z, 
		v_x=c.vx, v_y=c.vy, v_z=c.vz)
	global obs_c = lguys.transform(lguys.Observation, obs_c_galcen)
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
obs = TOML.parsefile("../../mc_orbits/orbit1.toml")

# ╔═╡ 8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
let 
	fig = Figure()

	dy = 2
	dx = dy * 1/cosd(obs["dec"])
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=(obs["ra"] .+ (-dx, dx), obs["dec"] .+ (-dy, dy)),
		aspect = 1,
		
	)

	x = [o.ra for o in obs_pred]
	y = [o.dec for o in obs_pred]

	
	h = Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100, colorscale=log10, colorrange=(1e-10, nothing))
	errscatter!([obs["ra"]], [obs["dec"]], color=COLORS[3], ms=50)

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ ce1cc3f7-a421-4676-8531-496cb099c98e
Heatmap

# ╔═╡ 76d3b8ee-94f3-4d8e-adeb-b7f7346a3c5b
obs

# ╔═╡ 9d74ffbb-4c31-4062-9434-7755f53e4da0
let
	dr = 0.05
	limits = (obs["pm_ra"] .+ (-dr, dr), obs["pm_dec"] .+ (-dr, dr))
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\mu_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\mu_\delta / \textrm{mas\,yr^{-1}}",
	limits=limits,
	aspect=DataAspect())

	x = [o.pm_ra for o in obs_pred]
	y = [o.pm_dec for o in obs_pred]
	
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100, 
		colorscale=log10, colorrange=(1e-4, nothing))
	
	errscatter!([obs["pm_ra"]], [obs["pm_dec"]], xerr=[obs["pm_ra_err"]], yerr=[obs["pm_dec_err"]], color=:red)
	
	fig
end

# ╔═╡ 6fe6deb4-ae44-4ca0-9617-95841fdaf791
let
	fig = Figure()
	dx = 4
	dy = 20
	limits = (obs["distance"] .+ (-dx, dx), obs["radial_velocity"] .+ (-dy, dy))
	
	ax = Axis(fig[1,1],
		xlabel="distance / kpc",
		ylabel = "radial velocity / km/s",
		limits=limits
	)


	
	x = [o.distance for o in obs_pred]
	y = [o.radial_velocity for o in obs_pred]
	
	h = Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100,
		colorscale=log10, colorrange=(1e-5, nothing)
	)

	errscatter!([obs["distance"]], [obs["radial_velocity"]], xerr=[obs["distance_err"]], yerr=[obs["radial_velocity_err"]], color=:red)

	Colorbar(fig[1, 2], h)
	
	fig
end

# ╔═╡ 7d2bfa85-e364-4e60-8187-1ed63328427c
heatmap

# ╔═╡ a0a2e31d-7700-42ca-94f0-e6ae08da21ce
obs

# ╔═╡ 59202a69-fb7e-4b10-854d-5f37dbf1b1fa
let
	fig = Figure()
	ax = Axis(fig, limits=((1, 2), nothing))
	ax.limits.val
end

# ╔═╡ e23e85f9-7667-4f6e-8af6-2516fa292e2b
md"""
# Velocity profile
"""

# ╔═╡ 0dd09c1e-67c0-4f23-bd12-41cbef62e4de
import StatsBase: weights, std

# ╔═╡ f348c1b6-5196-4546-ba22-72677611a420
begin 
	ra0 = obs_c.ra
	dec0 = obs_c.dec

	xi, eta = lguys.to_tangent(df["ra"], df["dec"], ra0, dec0)
	
	r_ell = lguys.calc_r_ell(xi, eta, 1, 1, 0) * 60
end

# ╔═╡ f5e72ea1-0d69-4db6-921f-4d987f66d683
df

# ╔═╡ d32757d2-dc08-488c-9ddd-d3eefefa2db7
function calc_σv(r_ell, rv, mass;  r_max = 30)
	filt = r_ell .< r_max
	vel = rv[filt]

	return std(vel, weights(mass)[filt])
end

# ╔═╡ 194ff30a-31a7-44bc-ac54-722a629457fc
σv = calc_σv(r_ell, df["radial_velocity"], df["probability"], r_max=240)

# ╔═╡ fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
gaussian(x, μ, σ) = 1/sqrt(2π)* 1/σ * exp(-(x-μ)^2/(2σ^2))

# ╔═╡ ed097e29-d2dc-4fa7-94e5-483b380600cc
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="RV (km/s)",
		ylabel="density"
	)
	
	filt = r_ell .< 240
	rv = df["radial_velocity"][filt]
	mass = df["probability"][filt]

	μ = Arya.mean(rv, Arya.sb.weights(mass))
	
	h = Arya.histogram(rv, weights=mass, normalization=:pdf)

	scatter!(lguys.midpoint(h.bins), h.values)

	x_model = LinRange(80, 130, 100)
	y_model = gaussian.(x_model, μ, σv)
	lines!(x_model, y_model)
	fig
end

# ╔═╡ b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
@doc raw"""
	calc_rb(\sigma, delta_t)

Given a radial velocity dispersion (in km/s), and the time since the last pericentre (in Gyr), returns the break radius (in kpc).

See peñarrubia et al. (200X) for the original equation.
```math
r_{\rm break} = C\,\sigma_{v}\,\Delta t
```
where $C=0.55$ is an empirical fit
"""
function calc_rb(σ, delta_t)
	kpc_per_kms_Gyr = 0.9777
	return 0.55 * σ * (delta_t)
end

# ╔═╡ 750185c9-a317-4978-aa04-486e5bfb7a63
r_b_kpc = calc_rb(σv, orbit_props["t_last_peri"]) # kpc

# ╔═╡ 1866c280-89c3-4a71-9dbf-50b583360145
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel = L"\log \rho_\star", 
		limits=((-1.9, 1), (-8, 2)))

	plot_ρ_s!(snap_i, label="initial")
	plot_ρ_s!(snap_f, label="final")
	vlines!(log10.(r_b_kpc))
	
	axislegend(ax)
	fig
end

# ╔═╡ 02183d68-bf7a-4774-9b48-f50712eeb552
r_b_arcmin = 206265 / 60 * obs_c.distance / 1e3

# ╔═╡ 9b75409d-55f5-47c3-ab63-8168d31d3d54
md"""
# properties
"""

# ╔═╡ f293b8ce-286f-430b-a8c9-9bf4ba44528d
let
	props_f = Dict(
		"ra" => obs_c.ra,
		"dec" => obs_c.dec,
		"pm_ra" => obs_c.pm_ra,
		"pm_dec" => obs_c.pm_dec,
		"sigma_los" => σv,
		"r_b" => r_b_arcmin,
	)

	filename = splitext(basename(starsfile))[1] * ".yml"
	open(filename, "w") do f
		TOML.print(f, props_f)
	end

	props_f
end

# ╔═╡ 6408828d-d570-4585-8fa8-24661857fb35
function get_M_h(output::lguys.Output, radius; idxs=(1:10:length(output)))
	N = length(idxs)
	M = Vector{Float64}(undef, N)
	
	for i in eachindex(idxs)
		snap = output[idxs[i]]
		ϵ = lguys.calc_ϵ(snap)
		filt = ϵ .> 0
		
		rs = lguys.calc_r(snap[filt])
		filt2 = rs .< radius

		ps = probabilities[snap.index[filt][filt2]]
		M[i] = sum(ps)
	end

	return M
end

# ╔═╡ df2fc50f-3a03-424e-ae88-02ae207a04ee
probabilities[out[1].index]

# ╔═╡ b36b594c-23b4-4683-8a12-3fa1b4c8b0d9
r_h = 0.1

# ╔═╡ 03cdff7d-6e86-45ab-9837-9029293639fa
probabilities

# ╔═╡ cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="normalized stellar mass within $n_rh rh",
		# yscale=log10,
		#yticks=[1, 0.1],
	)
	

	M_dm_h = get_M_h(out, n_rh * r_h)
	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1])
	
	
	fig
end

# ╔═╡ f3db0a1d-78bc-4ac4-84b6-d53f27a28795
out[10].x_cen

# ╔═╡ a0391689-66a2-473f-9704-e12a3d033d13
import LinearAlgebra: dot

# ╔═╡ 44ab0a25-ab4c-4e90-8619-2a068a285755
function calc_v_rad(snap)
	x_vec = snap.positions .- snap.x_cen
	v_vec = snap.velocities .- snap.v_cen

	x_hat = x_vec / lguys.calc_r(x_vec)'

	v_rad = sum(x_hat .* v_vec, dims=1)
	return vec(v_rad)
end

# ╔═╡ c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
function calc_v_tan(snap)
	v = lguys.calc_v(snap)
	v_rad = calc_v_rad(snap)
	v_tan = @. sqrt(v^2 - v_rad^2)
	return vec(v_tan)
end

# ╔═╡ 89ea28dc-2677-4e6f-9066-d5ffff3b474c
function calc_v_los(snap)
	(snap.velocities[1, :] .- snap.v_cen[1])
end

# ╔═╡ b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
function calc_σ_v(snap; r_max=1)
	rs = lguys.calc_r(snap)
	filt = rs .< r_max
	vs = calc_v_los(snap[filt])
	masses = probabilities[snap[filt].index]
	return lguys.mean(vs .^ 2, weights(masses))^0.5 * lguys.V0
end

# ╔═╡ 76438e61-d227-41cf-b9ea-e658bc389772
let
	snap = snap_f
	
	mass = probabilities[snap.index]
	v_rad = calc_v_los(snap)  * lguys.V0
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< 1

	h = Arya.histogram(v_rad[filt], 
		weights=mass[filt], normalization=:pdf, limits=(-20, 20))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_{x} \quad [\textrm{km\,s^{-1}}]",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoint(h.bins), h.values)

	σ = calc_σ_v(snap)
	x_model = LinRange(-20, 20, 100)
	y_model = gaussian.(x_model, 0, σ)
	lines!(x_model, y_model, 
		color=COLORS[2], label="gaussian (σ = $(round(σ, digits=2)) km / s)"
	)

	axislegend()
	fig
	
end

# ╔═╡ 04629bcf-5c19-41eb-8903-72947c209cbf
let
	snap = snap_f
	
	mass = probabilities[snap.index]
	v_rad = calc_v_rad(snap) * lguys.V0
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< 1

	h = Arya.histogram(v_rad[filt], 
		weights=mass[filt], normalization=:pdf, limits=(-20, 20))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_\textrm{rad} = \vec{v} \cdot \hat{r} \  / \ \textrm{km\,s^{-1}}",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoint(h.bins), h.values)

	σ = calc_σ_v(snap)
	x_model = LinRange(-20, 20, 100)
	y_model = gaussian.(x_model, 0, σ*√3/2)
	lines!(x_model, y_model)
	fig
	
end

# ╔═╡ d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
calc_σ_v(snap_i) / sqrt(3)

# ╔═╡ d42795d0-bd69-4c2c-be5b-e27e85199ee3
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel=L"\langle v_{x,\,\star} \rangle_{r \leq 1\,\textrm{kpc}} \ / \ \textrm{km\,s^{-1}}",
		# yscale=log10,
		#yticks=[1, 0.1],
	)

	idx = 1:10:length(out)
	vs = [calc_σ_v(out[i]) for i in idx]

	scatter!(out.times[idx] * lguys.T0, vs)
	
	fig
end

# ╔═╡ 6e34b91c-c336-4538-a961-60833d37f070
function v_rad_hist(snap, bins=40)

	mass = probabilities[snap.index]
	v_rad = calc_v_rad(snap)
	logr = log10.(lguys.calc_r(snap))
	x_bins, v_bins = lguys.calc_histogram(logr, bins, weights=v_rad .* mass)
	_, counts = lguys.calc_histogram(logr, bins, weights=mass)

	return x_bins, v_bins ./ counts
end

# ╔═╡ 17c880ba-d6de-4d9c-b44b-d5ecdb4da226
calc_v_rad(snap_i)

# ╔═╡ 57d5decd-8259-4ec2-87ac-44d28625cd7b
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="log r",
		ylabel="mean radial velocity",
		limits=((nothing, 1), (-60, 60))
	)
	
	x, y= v_rad_hist(snap_f, 100)
	scatter!(lguys.midpoint(x), y * lguys.V0)

	vlines!(log10.(r_b_kpc))
	fig
end

# ╔═╡ 98bcfb2a-8eb5-4c7c-8dea-46d000a6d245
r_b_arcmin

# ╔═╡ 0c26f965-0381-4f79-a6ce-0772ae922b3f
let
	snap = snap_f
	
	mass = probabilities[snap.index]
	v_rad = calc_v_rad(snap) * lguys.V0
	logr = log10.(lguys.calc_r(snap))

	limits = (-2, 1, -100, 100)
	fig = Figure()
	ax = Axis(fig[1,1],
		limits=limits,
		xlabel=L"\log r / \textrm{kpc}",
		ylabel=L"v_\textrm{rad} / \textrm{km\,s^{-1}}"
	)
	
	h = Arya.hist2d!(ax, logr, v_rad, weights=mass, bins=200, colorrange=(1e-10, 1), colorscale=log10,
	)
	vlines!(log10.(r_b_kpc))

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ Cell order:
# ╠═377284f2-dcee-44d3-9a04-728605cea92a
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╠═8ed2e715-1095-4c99-9ecb-54fecbc27d7f
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═5e12d306-a430-4b15-b3a7-d4806a5856cd
# ╠═ef3481f8-2505-4c04-a04a-29bdf34e9abb
# ╠═172588cc-ae22-440e-8488-f508aaf7ce96
# ╠═973955ad-3214-42cf-831f-a782f1d2434a
# ╠═6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
# ╠═8f9ee2a7-de43-43e0-8257-93ecc630044f
# ╠═396a53a3-de0f-4d97-9693-40f3757d66f9
# ╠═f9fe37ef-de81-4d69-9308-cda968851ed2
# ╠═6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
# ╠═77479cd4-513c-4603-9aa0-1acd964c403a
# ╠═e33245ad-7562-47af-94b9-9708af45b2a4
# ╠═7bc2c15c-33f7-43f3-a47f-ca39ffc22071
# ╠═8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
# ╟─c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
# ╠═a1138e64-5fa5-4a0e-aeef-487ee78a7adc
# ╠═1866c280-89c3-4a71-9dbf-50b583360145
# ╟─f5a3ea2f-4f19-4b0e-af55-74902f2c6485
# ╠═1f722acb-f7b9-4d6c-900e-11eae85e0708
# ╠═2e03e2e9-86aa-419b-a1c8-b05da3abce61
# ╠═fe56d883-1e4b-400c-b5f2-03ac058bb0ad
# ╠═7e588ae3-89f3-4b91-8963-f6bf4391a859
# ╠═cf6a7cbb-1034-4026-a3f3-1e854d2929e2
# ╠═ead40db8-ab7c-4167-b2ca-b00b5cb3cde4
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
# ╠═ce1cc3f7-a421-4676-8531-496cb099c98e
# ╠═76d3b8ee-94f3-4d8e-adeb-b7f7346a3c5b
# ╠═9d74ffbb-4c31-4062-9434-7755f53e4da0
# ╠═6fe6deb4-ae44-4ca0-9617-95841fdaf791
# ╠═7d2bfa85-e364-4e60-8187-1ed63328427c
# ╠═a0a2e31d-7700-42ca-94f0-e6ae08da21ce
# ╠═59202a69-fb7e-4b10-854d-5f37dbf1b1fa
# ╟─e23e85f9-7667-4f6e-8af6-2516fa292e2b
# ╠═0dd09c1e-67c0-4f23-bd12-41cbef62e4de
# ╠═f348c1b6-5196-4546-ba22-72677611a420
# ╠═f5e72ea1-0d69-4db6-921f-4d987f66d683
# ╠═d32757d2-dc08-488c-9ddd-d3eefefa2db7
# ╠═194ff30a-31a7-44bc-ac54-722a629457fc
# ╠═fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
# ╠═ed097e29-d2dc-4fa7-94e5-483b380600cc
# ╟─76438e61-d227-41cf-b9ea-e658bc389772
# ╟─04629bcf-5c19-41eb-8903-72947c209cbf
# ╟─b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
# ╠═750185c9-a317-4978-aa04-486e5bfb7a63
# ╠═02183d68-bf7a-4774-9b48-f50712eeb552
# ╠═9b75409d-55f5-47c3-ab63-8168d31d3d54
# ╠═f293b8ce-286f-430b-a8c9-9bf4ba44528d
# ╠═6408828d-d570-4585-8fa8-24661857fb35
# ╠═df2fc50f-3a03-424e-ae88-02ae207a04ee
# ╠═b36b594c-23b4-4683-8a12-3fa1b4c8b0d9
# ╠═03cdff7d-6e86-45ab-9837-9029293639fa
# ╠═cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
# ╠═f3db0a1d-78bc-4ac4-84b6-d53f27a28795
# ╠═b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
# ╠═d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
# ╠═d42795d0-bd69-4c2c-be5b-e27e85199ee3
# ╠═a0391689-66a2-473f-9704-e12a3d033d13
# ╠═44ab0a25-ab4c-4e90-8619-2a068a285755
# ╠═c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
# ╠═89ea28dc-2677-4e6f-9066-d5ffff3b474c
# ╠═6e34b91c-c336-4538-a961-60833d37f070
# ╠═17c880ba-d6de-4d9c-b44b-d5ecdb4da226
# ╠═57d5decd-8259-4ec2-87ac-44d28625cd7b
# ╠═98bcfb2a-8eb5-4c7c-8dea-46d000a6d245
# ╠═0c26f965-0381-4f79-a6ce-0772ae922b3f
