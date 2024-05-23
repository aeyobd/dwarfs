### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ 60fe0995-3dd6-48ff-86e8-6fc5e099e39b
begin 
	using HDF5

	f = h5open("star_probabilities.hdf5")
	p_idx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
end

# ╔═╡ a227ef84-a874-4176-9253-36bc36a743fa
using LsqFit: curve_fit

# ╔═╡ fcb7ae4e-9e8b-4a64-8eda-edba3f4dac64
using SciPy

# ╔═╡ 8b25b3e6-70de-49f8-894a-d803353998d6
using Measurements

# ╔═╡ 7094bc54-deb4-48a5-bf09-9ee6c684ac3c
begin 
	out =  lguys.Output("out/combined.hdf5")
	
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
	out.x_cen .= x_cen
	out.v_cen .= v_cen

	out
end

# ╔═╡ f4d42909-5397-41e6-ac0c-07185bfece0e
r_h = 0.109

# ╔═╡ a71bb30e-4a4c-4e94-b318-a426e3ee3045
begin 
	orbit_expected = CSV.read("orbit.csv", DataFrame)
	x_cen_exp = transpose(hcat(orbit_expected.x, orbit_expected.y, orbit_expected.z))
	v_cen_exp = -transpose(hcat(orbit_expected.v_x, orbit_expected.v_y, orbit_expected.v_z))

end

# ╔═╡ 29887611-5a0b-4f3a-8a3f-2da94b1765d2
lguys.plot_xyz(x_cen, x_cen_exp, labels=["n body", "point particle"])

# ╔═╡ 46f6e4a2-6e2b-4f23-807c-b41be51a0960
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="R / kpc", ylabel="z / kpc"
	)
	x = x_cen[1, :]
	y = x_cen[2, :]
	z = x_cen[3, :]
	R = @. sqrt(x^2 + y^2)
	lines!(R, z)
	fig
end

# ╔═╡ a9e79439-16a4-4908-bfe0-f0770cdb26df
md"""
# Mass evolution
"""

# ╔═╡ dfa506e2-095a-47a2-9f14-baa04dd45d2d
function get_M_h(output::lguys.Output, radius; idxs=(1:10:length(output)))
	N = length(idxs)
	M = Vector{Float64}(undef, N)
	Ms = Vector{Float64}(undef, N)
	
	for i in eachindex(idxs)
		snap = output[idxs[i]]
		ϵ = lguys.calc_ϵ(snap)
		filt = ϵ .> 0
		
		rs = lguys.calc_r(snap[filt])
		filt2 = rs .< radius

		M[i] = sum(filt2)

		ps = probabilities[snap.index[filt]]

		Ms[i] = sum(ps[filt2])
		
	end

	return M, Ms
end

# ╔═╡ 807c46f6-56e0-4ab4-b31f-ac4a3fec9761
collect(1:4:10)

# ╔═╡ 3755ea87-bc79-4b2e-b84f-5b08df7aa086
let
	fig = Figure()
	ax = Axis(fig[1,1],
	xlabel="time / Gyr",
	ylabel="normalized mass within 1 rh"
	)

	M_dm_h, M_s_h = get_M_h(out, r_h)

	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="DM")
	scatter!(out.times[1:10:end] * lguys.T0, M_s_h ./ M_s_h[1], label="stars")

	Legend(fig[1, 2], ax)
	
	fig
end

# ╔═╡ 5bd14979-a02d-4c03-a6e9-090a63a25059
let
	fig = Figure()
	ax = Axis(fig[1,1],
	xlabel="time / Gyr",
	ylabel="normalized mass within 10 rh"
	)

	M_dm_h, M_s_h = get_M_h(out, 10r_h)

	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="DM")
	scatter!(out.times[1:10:end] * lguys.T0, M_s_h ./ M_s_h[1], label="stars")

	Legend(fig[1, 2], ax)
	println(M_dm_h[1], M_s_h[1])
	fig
end

# ╔═╡ 0fa11815-3ab0-4b19-9be7-186b7c2c1063
let
	fig = Figure()
	ax = Axis(fig[1,1],
	xlabel="time / Gyr",
	ylabel="normalized mass within 100 rh",
	yscale=log10
	)

	M_dm_h, M_s_h = get_M_h(out, 100r_h)

	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="DM")
	scatter!(out.times[1:10:end] * lguys.T0, M_s_h ./ M_s_h[1], label="stars")

	Legend(fig[1, 2], ax)
	println(M_dm_h[1], M_s_h[1])
	fig
end

# ╔═╡ df260566-9d3f-4f02-8081-30bc7f277d10
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="normalized bound mass",
		yscale=log10
	)

	M_dm_h, M_s_h = get_M_h(out, 1e10r_h)

	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="DM")
	scatter!(out.times[1:10:end] * lguys.T0, M_s_h ./ M_s_h[1], label="stars")

	Legend(fig[1, 2], ax)
	println(M_dm_h[1], M_s_h[1])
	fig
end

# ╔═╡ 0e50f90a-2024-433d-a5b9-b1352c8095cc
size(out[1].positions)

# ╔═╡ e14fa4a1-6175-4b9f-ad01-525c1617fe63
md"""

# Evolution of circular Velocity
"""

# ╔═╡ d2491f4d-f2b3-47c0-894c-94529b8cbe55
function get_Vc(snap, skip=10)
	ϵ = lguys.calc_ϵ(snap)
	filt = ϵ .> 0
	r = lguys.calc_r(snap[filt])
	ms = snap.masses[filt][sortperm(r)]
	r = sort(r)

	r = r[1:skip:end]
	M = cumsum(ms)[1:skip:end]

	Vc = @. sqrt(lguys.G * M / r) * lguys.V0
	r, Vc 
end

# ╔═╡ 3b044830-d631-4379-8401-cba29523e2f0
function plot_Vc!(snap; kwargs...)

	rc, Vc = get_Vc(snap)
	lines!(log10.(rc), Vc; kwargs...)
end

# ╔═╡ 86877090-cfc9-4621-9187-d3da81f0428f
import NaNMath as nm

# ╔═╡ 1ace0815-7dc6-486f-b4c4-95fb41bc4313
import StatsBase: percentile

# ╔═╡ b7492731-148e-4445-a35b-23c9a444ef48
function Vc_model(r, param)
	Rs,scale =param
	x = r ./ Rs
	inner = @. (nm.log(1+x) - x/(1+x)) / x
	return @. scale *nm.sqrt(inner)
end

# ╔═╡ 93e0c066-0229-4b49-95dc-f0a9a9f742e0
function fit_Vc(rc, Vc; percen=80)
	filt = Vc .> percentile(Vc, percen)
	fit = curve_fit(Vc_model, rc[filt], Vc[filt], [2., 30.])
	
	return (; r_c=fit.param[1], V_c=fit.param[2], fit=fit,
	r_min=minimum(rc[filt]),
	r_max=maximum(rc[filt])
	)
end

# ╔═╡ 8547a5a2-7170-47a3-8110-cc2a5d49c070
function V_r_max(Rs, Vs)

	r = 2.16258 * Rs
	return r, Vc_model(r, [Rs, Vs])
end

# ╔═╡ 086b4aa9-5671-4c31-9a01-f5585464bb8a
function get_Vh(snap, r_h)
	r = lguys.calc_r(snap)
	m = sum(snap.masses[r .< r_h])
	return sqrt(lguys.G * m / r_h) * lguys.V0
end

# ╔═╡ 9a8a5ad2-6aff-4d37-9748-7a415568f751
function get_V_r_max(output::lguys.Output; skip=1)
	V_maxs = Float64[]
	r_maxs = Float64[]
	V_hs = Float64[]

	for i in 1:skip:length(output)
		snap = output[i]
		rc, Vc = get_Vc(snap)
		fit = fit_Vc(rc, Vc)
		r_max, V_max = V_r_max(fit.r_c, fit.V_c)

		Vh = get_Vh(snap, r_h)
		
		push!(V_maxs, V_max)
		push!(r_maxs, r_max)
		push!(V_hs, Vh)
	end

	return r_maxs, V_maxs, V_hs
end

# ╔═╡ 9418347c-e49f-4483-be0f-98a3b9578bac
r_max, V_max, V_h = get_V_r_max(out, skip=10)

# ╔═╡ c068c177-e879-4b8e-b1af-18690af9b334
let 
	i = length(out)
	snap = out[i]

	
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"V_\textrm{circ}",
	title="time = $(out.times[i]  * lguys.T0) Gyr")

	plot_Vc!(snap, label="initial")
	r, V = get_Vc(snap)
	fit = fit_Vc(r, V)
	r_max, V_max = V_r_max(fit.r_c, fit.V_c)



	r_fit = LinRange(fit.r_min, fit.r_max, 100)
	V_fit = Vc_model(r_fit, [fit.r_c, fit.V_c])

	lines!(log10.(r_fit), V_fit, linewidth=3)
	
	scatter!(log10.(r_max), V_max, color=Arya.COLORS[3])


	ax2 = Axis(fig[2, 1], xlabel="r/kpc", ylabel="residual")

	filt = fit.r_min .<= r .<= fit.r_max
	dV = V[filt] .- Vc_model(r[filt], [fit.r_c, fit.V_c])
	scatter!(r[filt], dV)
	vlines!(r_max)
	
    rowsize!(fig.layout, 2, Relative(1/4))
	
	fig
end

# ╔═╡ 245721a6-01aa-43e7-922d-ed5da02207c1
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel=L"V_\text{circ} / \text{km\,s^{-1}}",
	limits=(nothing, (0, nothing)))
	x = out.times[1:10:end] * lguys.T0
	scatter!(x, V_max, label=L"maximum $V_\text{circ}$")
	scatter!(x, V_h, label=L"r=r_h")
	axislegend(ax)
	
	fig
end

# ╔═╡ c5796d82-013b-4cdc-a625-31249b51197d
md"""
# Density evolution
"""

# ╔═╡ f563cbc7-655b-46e3-8686-2e4561b2467a
function plot_ρ_dm!(snap, x_cen; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	pos .-= x_cen
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ abe6f9cd-2b49-4826-ba9c-56244a10bff0
let 
	i = 15

	snap = out[i]

	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"V_\textrm{circ}",
		title="circular velocity at snap$i"
	)


	plot_Vc!(snap)

	fig

	rc, Vc = get_Vc(snap)
	fit_Vc(rc, Vc)
end

# ╔═╡ 1ebb1ab6-c1a0-4d6b-8529-65155575b96e
mlog10 = Makie.pseudolog10

# ╔═╡ 743d20cd-9636-4384-9bf2-b2d7e259ae7d
# ╠═╡ disabled = true
#=╠═╡
r_h = 0.308
  ╠═╡ =#

# ╔═╡ 32cd2b3b-b756-4fbf-b348-5839fe1c0360
md"""
# Sky Projections
"""

# ╔═╡ 07c3bfc4-615e-4847-a66f-fb824f384ce8
begin 
	snap_cen = lguys.Snapshot(x_cen, v_cen, zeros(size(x_cen, 1)))

	obs_c = lguys.to_sky(snap_cen)
end

# ╔═╡ 4962c654-9f9f-44b9-b0ee-68791522562a
begin 
	snap_o = lguys.Snapshot(x_cen_exp, v_cen_exp, zeros(size(x_cen_exp, 1)))

	obs_o = lguys.to_sky(snap_o)
end

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	idx_i = 1
	idx_f = 1040
	snap_i = out[idx_i]
	snap_f = out[idx_f]
	
	"$(out.times[[idx_i, idx_f]] * lguys.T0) Gyr"
end
  ╠═╡ =#

# ╔═╡ d00f2ff5-a552-43ef-be40-b77b4dc33730
#=╠═╡
begin
	r = lguys.calc_r(x_cen)
	peri_filt = idx_f-200:idx_f
	t_last_peri_arg = argmin(r[peri_filt])
	t_last_peri = out.times[peri_filt[t_last_peri_arg]] * lguys.T0
	delta_t_peri = out.times[idx_f] * lguys.T0 - t_last_peri
end
  ╠═╡ =#

# ╔═╡ 256b5101-78db-4bcf-aa67-0dd51cce6222
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = lguys.calc_r(x_cen)
	lines!(out.times * lguys.T0, r)
	lines!(lguys.T0*(orbit_expected.t .- orbit_expected.t[begin]), lguys.calc_r(x_cen_exp))
	scatter!(out.times[idx_f] * lguys.T0, r[idx_f])
	vlines!(t_last_peri)
	fig
end
  ╠═╡ =#

# ╔═╡ db320665-f46d-4aed-a2b2-4b39bcb605c5
#=╠═╡
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"V_\textrm{circ}")

	plot_Vc!(snap_i, label="initial")

	plot_Vc!(snap_f, label="final")

	α = 0.4
	β = 0.65
	x = LinRange(1, 0.1, 100)

	y = @. 2^α * x^β * (1 + x^2)^(-α)
	lines!(log10.(x .* r_max[1]), y .* V_max[1], lw=10, label="EN21")

	scatter!(log10.(r_max), V_max, color=Arya.COLORS[4], label="Vmax")

	axislegend(ax)
	fig
end
  ╠═╡ =#

# ╔═╡ eda7c370-dc70-4b82-8648-45af67918298
#=╠═╡
begin
	rc, Vc = get_Vc(snap_f)
	t_cross = rc ./ Vc
	r_t_cross = lguys.lerp(t_cross, rc)
	
	r_break = r_t_cross(t_last_peri)
end
  ╠═╡ =#

# ╔═╡ c9a86f53-1909-4509-98a8-7122bf4e5e3b
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="r", ylabel="t cross")
	
	plot!(rc, t_cross)
	hlines!(t_last_peri)
	vlines!(r_break)
	fig
end
  ╠═╡ =#

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
#=╠═╡
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=(nothing, (-13, 0)))
	
	plot_ρ_dm!(snap_i, x_cen[:, idx_i], label="initial")
	plot_ρ_dm!(snap_f, x_cen[:, idx_f], label="final")
	
	axislegend(ax)
	lines!([0, 0] .+ log10.(r_break), [-11.5, -10],  color=:black)
	scatter!(log10.(r_break), -10, marker=:utriangle, color=:black)

	text!(L"r_\textrm{break}", position=(log10.(r_break),-11.5), space=:data, rotation=π/2, align=(:left, :baseline))
	fig
	# only include bound points in profile...
end
  ╠═╡ =#

# ╔═╡ 4801ff80-5761-490a-801a-b263b90d63fd
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bins = LinRange(-10, 10, 100)
	colorrange=(1e1, 1e3)
	
	bins = (x_cen[1, idx_i]  .+ bins,  x_cen[2, idx_i]  .+bins)
	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")

	bins = LinRange(-10, 10, 100)
	bins = (x_cen[1, idx_f]  .+ bins,  x_cen[2, idx_f]  .+bins)
	hm = Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], bins = bins, colorscale=log10, colorrange=colorrange)
	
	Colorbar(fig[:, end+1], hm)
	
	fig
end
  ╠═╡ =#

# ╔═╡ fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
#=╠═╡
let
	fig = Figure()
	r_max = 130
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc", title="dark matter",
	limits=(-r_max, r_max, -r_max, r_max))
	bins = LinRange(-r_max, r_max, 100)
	
	Arya.hist2d!(ax, snap_f.positions[2, :], snap_f.positions[3, :], bins = bins, colorscale=mlog10)

	fig
end
  ╠═╡ =#

# ╔═╡ 7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
#=╠═╡
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e2, 1e6)),
		xlabel=L"\epsilon", ylabel="count")
	stephist!(lguys.calc_ϵ(snap_i), )
	es = lguys.calc_ϵ(snap_f)
	es = es[es .> 0]
	stephist!(es)

	fig

end
  ╠═╡ =#

# ╔═╡ f2381a9c-96d7-4935-8967-32b8bfb3fb48
#=╠═╡
begin 
	p_min = 1e-2 * maximum(probabilities)
	snap_f_stars = snap_f[probabilities[snap_f.index] .> p_min]
	println(length(snap_f_stars))
	obs_pred = lguys.to_sky(snap_f_stars)
	m_star_f = probabilities[snap_f_stars.index]
end
  ╠═╡ =#

# ╔═╡ 68243ab7-0d82-48a2-8737-7cc105ae7325
#=╠═╡
let
	fig = Figure()

	dy = 2
	dx = dy * 1/cosd(obs.dec)
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=(obs.ra .+ (-dx, dx), obs.dec .+ (-dy, dy)),
		aspect = 1/cosd(obs.dec)
	)

	x = [o.ra for o in obs_pred]
	y = [o.dec for o in obs_pred]
	
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)

	scatter!(obs.ra, obs.dec, color=:red)

	idx = idx_f-10:idx_f+10
	x = [o.ra for o in obs_c[idx]]
	y = [o.dec for o in obs_c[idx]]
	
	l = lines!(x, y, color=out.times[idx])
	Colorbar(fig[1, 2], l)

	x = [o.ra for o in obs_o[end-30:end]]
	y = [o.dec for o in obs_o[end-30:end]]
	
	lines!(x, y, color=Arya.COLORS[2])
	
	fig
end
  ╠═╡ =#

# ╔═╡ 7ae8b02e-4c20-407f-bd00-b8866bb45bd6
#=╠═╡
let
	fig = Figure()

	dy = 1
	dx = dy * 1/cosd(obs.dec)
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=(obs.ra .+ (-dx, dx), obs.dec .+ (-dy, dy)),
		aspect = 1
	)

	x = [o.ra for o in obs_pred]
	y = [o.dec for o in obs_pred]
	
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)

	
	fig
end
  ╠═╡ =#

# ╔═╡ cc68d116-d91c-4774-ad6c-0520cee6b339
md"""
# 2D profile fits
"""

# ╔═╡ 815993ab-3527-4df2-99be-5f4b3e7163ba
function norm_hist(xs; bw=0.1, weights=nothing)
	bins = collect(minimum(xs):bw:(maximum(xs) + bw))
	x_h, y_h = lguys.calc_histogram(xs, bins, weights=weights)
	_, counts = lguys.calc_histogram(xs, bins)

    y_e = y_h ./ sqrt.(counts)
    
	area = length(xs)

	return x_h, y_h ./ area, y_e ./ area
end

# ╔═╡ 09038c96-974c-4b34-be01-d8356b88616d
function calc_Σ(log_r, hist)
	r = 10 .^ log_r
	Σ = hist ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
	return Σ
end

# ╔═╡ 8edd3e2b-fdf9-4124-bd6f-e5c7a5677a0d
import StatsBase: mean, weights

# ╔═╡ 98371ce4-a7e7-4091-9716-6a20cd9cb3f6
#=╠═╡
begin 
	ra_pred = [o.ra for o in obs_pred]
	dec_pred = [o.dec for o in obs_pred]
	ra_0 = mean(ra_pred, weights(m_star_f))
	dec_0 = mean(dec_pred, weights(m_star_f))

	xi, eta = lguys.to_tangent(ra_pred, dec_pred, ra_0, dec_0)
end
  ╠═╡ =#

# ╔═╡ c51fdd98-097c-4205-a494-87de0e662092
#=╠═╡
scatter(xi, eta, color=m_star_f)
  ╠═╡ =#

# ╔═╡ f8a4da63-b017-410e-8ac1-887b373b5578
#=╠═╡
radii = @. sqrt(xi^2 + eta^2) * 60
  ╠═╡ =#

# ╔═╡ d8eb16fe-8e16-45c7-9415-492530692adf
#=╠═╡
hist(radii)
  ╠═╡ =#

# ╔═╡ 6f783c76-e1d4-425a-8a6c-f3ff44457710
function calc_Γ(log_rs, ρs, step=1)
	dx = lguys.gradient(log_rs)
	dρ = lguys.gradient(log10.(ρs))

	return dρ ./ dx #lguys.gradient(log10.(ρs), log_rs)
end

# ╔═╡ c7d2f326-4d69-4647-b508-b0ac70269f1a
function calc_properties(rs; bw=0.1, weights=nothing)
    log_r_bin, counts, δ_counts = norm_hist(log10.(rs), bw=bw, weights=weights)
    counts = counts .± δ_counts
    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r .± δ_log_r
    r_bin = 10 .^ log_r_bin

    r = 10 .^ log_r

    As = π * diff((10 .^ log_r_bin) .^ 2)
    Σ = counts ./ As 
    # Σ_e = @. y_e / ys * Σ
    #Σ_m = calc_Σ_mean(xs, ys)

    M_in = cumsum(counts)
    A_in = @. π * (r_bin[2:end])^2
    Σ_m = M_in ./ A_in

    Γ = calc_Γ(log_r, Σ)
    Γ_max = @. 2*(1 - Σ / Σ_m)
    
    return (log_r=log_r, log_r_bins=log_r_bin, counts=counts, Σ=Σ, Σ_m=Σ_m, Γ=Γ, Γ_max=Γ_max, M_in=M_in, N=length(rs), M=sum(weights))
end

# ╔═╡ bb3284f4-2bdd-4c5b-8031-18af86dfa0f9
#=╠═╡
props_pred = calc_properties(radii, weights=m_star_f)
  ╠═╡ =#

# ╔═╡ 82043022-3ac8-451a-b8a1-236f03af1a98
#=╠═╡
props_pred
  ╠═╡ =#

# ╔═╡ d394bb02-4eca-4c13-858e-030a002f1108
function log_Σ_exp(r, A, r_s)
    return @. A + 1/log(10) * (-r / r_s)
end

# ╔═╡ e738405e-2b22-4768-a01d-ef167ed1dbbf
import LinearAlgebra: diag

# ╔═╡ b3171488-08c2-441a-b122-96c9f90a1cc9
import QuadGK: quadgk

# ╔═╡ 7a2e9210-084a-47f7-8dea-98bf48cc6712
function predict_properties(Σ_model; N=10_000, log_r_min=-2, log_r_max=2)
    log_r_bins = LinRange(log_r_min, log_r_max, 1000)
    log_r = lguys.midpoint(log_r_bins)
    r = 10 .^  log_r
    r_bins = 10 .^ log_r_bins
    
    Σ = Σ_model.(r)
    Γ = calc_Γ(log_r, Σ)
    M_in = [quadgk(rrr->2π*rrr*Σ_model(rrr), 0, rr)[1] for rr in r]
    Σ_m = M_in ./ (π * r .^ 2)
    Γ_max = 2*(1 .- Σ ./ Σ_m)
    counts =  Σ .* (2π *  r .* diff(r_bins) )
    
    return (log_r=log_r, log_r_bins=log_r_bins, counts=counts, Σ=Σ, Σ_m=Σ_m, Γ=Γ, Γ_max=Γ_max, M_in=M_in)

end

# ╔═╡ 225b2b60-d011-4340-810b-1f041b9552cd
log_r_label = "log r / arcmin"

# ╔═╡ c74d6489-b93a-4eb1-87a9-3b329ff6f527
COLORS=Arya.COLORS

# ╔═╡ 7bfb7485-5424-4c09-a0a9-4ed67b67fca9
begin 
	value(a::Measurement) = a.val
	value(a) = a
	err(a::Measurement) = a.err
	err(a) = 0
		function Makie.convert_single_argument(y::Array{Measurement{T}}) where T
		return value.(y)
	end
end

# ╔═╡ 4dbfea97-4151-4785-8584-9d7708dfd860
function plot_Σ_fit_res(obs, pred, res)
    fig = Figure()
    ax = Axis(fig[1, 1], 
        ylabel=L"\log \Sigma\ / \textrm{(stellar mass/arcmin^2)}")
    N = length(obs.log_r)
    y = log10.(obs.Σ )
    p = errorbars!(ax, value.(obs.log_r), value.(y), err.(y))
    scatter!(ax, value.(obs.log_r), value.(y))

    lines!(ax, pred.log_r, log10.(pred.Σ ), color=COLORS[2])
    
    ax2 = Axis(fig[2, 1],
        ylabel=L"\delta\log\Sigma", 
    	xlabel=log_r_label,
		limits = (nothing, (-1, 1))
	)

#     p2 = plot(ylabel=L"\Delta \log\Sigma", xlabel=log_r_label, ylim=(-2, 2))

	y = res
    scatter!(ax2, value.(obs.log_r), value.(y), err.(y), label="")
    errorbars!(ax2, value.(obs.log_r), value.(y), err.(y))


    hlines!(0, color=:black)
    
    rowsize!(fig.layout, 2, Relative(1/4))

    linkxaxes!(ax, ax2)
    hidexdecorations!(ax, grid=false)
#     return plot(p1, p2, layout=grid(2, 1, heights=(0.8, 0.2)), link=:x, bottom_margin=[-5Plots.mm 0Plots.mm])
    return Makie.FigureAxisPlot(fig, ax, p)
end

# ╔═╡ 8f8ad8f2-8e30-4f07-ba53-df4af34f7838
function fit_profile(obs; r_max=Inf, N=10_000)
    r_val = [10 ^ r.val for r in obs.log_r]
    log_Σ = log10.(obs.Σ)
    filt = r_val .< r_max
    filt .&= map(x->isfinite(x), log_Σ)
    filt .&= @. !isnan(log_Σ)
    
    r_val = r_val[filt]
    log_Σ = log_Σ[filt]


    log_Σ_val = [s.val for s in log_Σ]
    log_Σ_e = [s.err for s in log_Σ]

    popt, covt = SciPy.optimize.curve_fit(log_Σ_exp, r_val, log_Σ_val, 
        sigma=log_Σ_e, p0=[1, 0.1])
    
    popt_p = popt .± sqrt.(diag(covt))
    println("log_Σ_0 = $(popt_p[1])")
    println("r_s = $(popt_p[2])")
    props = predict_properties(r->10 .^ log_Σ_exp(r, popt...), 
        N=N, log_r_min=obs.log_r_bins[1], log_r_max=obs.log_r_bins[end])
    
    log_Σ_pred = log_Σ_exp.(10 .^ value.(obs.log_r), popt...)
    log_Σ_res = log10.(obs.Σ) .- log_Σ_pred
    return popt_p, props, log_Σ_res
end

# ╔═╡ e2e4fe18-aa1e-43f3-86fc-4967390eca4d
#=╠═╡
popt, pred, res = fit_profile(props_pred)	
  ╠═╡ =#

# ╔═╡ 7ad8d849-41e7-43c7-863b-f76e9e3869c8
#=╠═╡
pred
  ╠═╡ =#

# ╔═╡ 4fc2081c-0d5d-4fa6-88bf-9a26ef144249
#=╠═╡
let
	fig, ax, p = scatter(value.(props_pred.log_r), value.(props_pred.counts))

	ax.xlabel = "log r / arcmin"
	ax.ylabel = "count / bin"

	fig
end
  ╠═╡ =#

# ╔═╡ 8b4cb5dd-602b-416d-ae77-6928819bd0ce
#=╠═╡
let 
	fig, ax, p = plot_Σ_fit_res(props_pred, pred, res)

	log_x = LinRange(-1, 2, 1000)
	r_s = 6.59
	y = log_Σ_exp(10 .^ log_x, value.(popt[1]), r_s)
	lines!(ax,log_x, y, color=COLORS[3], linestyle=:dash, label="observed")
	fig
end
  ╠═╡ =#

# ╔═╡ adce4d28-ee26-4bcf-98d2-06360dc82004
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=(-0.8, 1.7, -8, 5))
	scatter!(value.(props_pred.log_r), value.(props_pred.Γ),)
	errorbars!(value.(props_pred.log_r), value.(props_pred.Γ), err.(props_pred.Γ) )
	
	lines!(pred.log_r, pred.Γ, label="exponential", color=COLORS[2])
	
	ax.xlabel = log_r_label
	ax.ylabel = L"\Gamma = d\,\log \Sigma / d\,\log r"
	
	fig
end
  ╠═╡ =#

# ╔═╡ d087bb85-0d8e-47f6-8519-f608db4df11f
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	scatter!(10 .^ value.(props_pred.log_r), value.(props_pred.Γ),)
	errorbars!(10 .^ value.(props_pred.log_r), value.(props_pred.Γ), err.(props_pred.Γ) )
	
	lines!(10 .^ pred.log_r, pred.Γ, label="exponential", color=COLORS[2])
	
	ax.xlabel = "r / arcmin"
	ax.ylabel = L"\Gamma = d\,\log \Sigma / d\,\log r"
	
	fig
end
  ╠═╡ =#

# ╔═╡ 91093e20-1b08-4ff6-a270-2db085c1999a
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=(-0.8, 2.2, -1, 2.2))
	scatter!(value.(props_pred.log_r), value.(props_pred.Γ_max),)
	errorbars!(value.(props_pred.log_r), value.(props_pred.Γ_max), err.(props_pred.Γ_max) )
	
	lines!(pred.log_r, pred.Γ_max, label="exponential", color=COLORS[2])
	
	ax.xlabel = log_r_label
	ax.ylabel = L"\Gamma_\mathrm{max} = 2(1 - \Sigma / \bar{\Sigma})"
	
	fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═7094bc54-deb4-48a5-bf09-9ee6c684ac3c
# ╠═60fe0995-3dd6-48ff-86e8-6fc5e099e39b
# ╠═f4d42909-5397-41e6-ac0c-07185bfece0e
# ╠═a71bb30e-4a4c-4e94-b318-a426e3ee3045
# ╠═29887611-5a0b-4f3a-8a3f-2da94b1765d2
# ╠═256b5101-78db-4bcf-aa67-0dd51cce6222
# ╠═d00f2ff5-a552-43ef-be40-b77b4dc33730
# ╠═46f6e4a2-6e2b-4f23-807c-b41be51a0960
# ╟─a9e79439-16a4-4908-bfe0-f0770cdb26df
# ╠═dfa506e2-095a-47a2-9f14-baa04dd45d2d
# ╠═807c46f6-56e0-4ab4-b31f-ac4a3fec9761
# ╠═3755ea87-bc79-4b2e-b84f-5b08df7aa086
# ╠═5bd14979-a02d-4c03-a6e9-090a63a25059
# ╠═0fa11815-3ab0-4b19-9be7-186b7c2c1063
# ╠═df260566-9d3f-4f02-8081-30bc7f277d10
# ╠═0e50f90a-2024-433d-a5b9-b1352c8095cc
# ╟─e14fa4a1-6175-4b9f-ad01-525c1617fe63
# ╠═d2491f4d-f2b3-47c0-894c-94529b8cbe55
# ╠═93e0c066-0229-4b49-95dc-f0a9a9f742e0
# ╠═3b044830-d631-4379-8401-cba29523e2f0
# ╠═a227ef84-a874-4176-9253-36bc36a743fa
# ╠═86877090-cfc9-4621-9187-d3da81f0428f
# ╠═1ace0815-7dc6-486f-b4c4-95fb41bc4313
# ╠═b7492731-148e-4445-a35b-23c9a444ef48
# ╠═8547a5a2-7170-47a3-8110-cc2a5d49c070
# ╠═086b4aa9-5671-4c31-9a01-f5585464bb8a
# ╠═9a8a5ad2-6aff-4d37-9748-7a415568f751
# ╠═9418347c-e49f-4483-be0f-98a3b9578bac
# ╠═db320665-f46d-4aed-a2b2-4b39bcb605c5
# ╠═c068c177-e879-4b8e-b1af-18690af9b334
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╟─c5796d82-013b-4cdc-a625-31249b51197d
# ╠═f563cbc7-655b-46e3-8686-2e4561b2467a
# ╠═abe6f9cd-2b49-4826-ba9c-56244a10bff0
# ╠═eda7c370-dc70-4b82-8648-45af67918298
# ╠═c9a86f53-1909-4509-98a8-7122bf4e5e3b
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═4801ff80-5761-490a-801a-b263b90d63fd
# ╠═1ebb1ab6-c1a0-4d6b-8529-65155575b96e
# ╠═fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
# ╠═743d20cd-9636-4384-9bf2-b2d7e259ae7d
# ╠═32cd2b3b-b756-4fbf-b348-5839fe1c0360
# ╠═f2381a9c-96d7-4935-8967-32b8bfb3fb48
# ╠═07c3bfc4-615e-4847-a66f-fb824f384ce8
# ╠═4962c654-9f9f-44b9-b0ee-68791522562a
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═68243ab7-0d82-48a2-8737-7cc105ae7325
# ╠═7ae8b02e-4c20-407f-bd00-b8866bb45bd6
# ╟─cc68d116-d91c-4774-ad6c-0520cee6b339
# ╠═815993ab-3527-4df2-99be-5f4b3e7163ba
# ╠═09038c96-974c-4b34-be01-d8356b88616d
# ╠═8edd3e2b-fdf9-4124-bd6f-e5c7a5677a0d
# ╠═98371ce4-a7e7-4091-9716-6a20cd9cb3f6
# ╠═c51fdd98-097c-4205-a494-87de0e662092
# ╠═f8a4da63-b017-410e-8ac1-887b373b5578
# ╠═d8eb16fe-8e16-45c7-9415-492530692adf
# ╠═fcb7ae4e-9e8b-4a64-8eda-edba3f4dac64
# ╠═6f783c76-e1d4-425a-8a6c-f3ff44457710
# ╠═8b25b3e6-70de-49f8-894a-d803353998d6
# ╠═c7d2f326-4d69-4647-b508-b0ac70269f1a
# ╠═4dbfea97-4151-4785-8584-9d7708dfd860
# ╠═82043022-3ac8-451a-b8a1-236f03af1a98
# ╠═8f8ad8f2-8e30-4f07-ba53-df4af34f7838
# ╠═bb3284f4-2bdd-4c5b-8031-18af86dfa0f9
# ╠═d394bb02-4eca-4c13-858e-030a002f1108
# ╠═e738405e-2b22-4768-a01d-ef167ed1dbbf
# ╠═7a2e9210-084a-47f7-8dea-98bf48cc6712
# ╠═b3171488-08c2-441a-b122-96c9f90a1cc9
# ╠═e2e4fe18-aa1e-43f3-86fc-4967390eca4d
# ╠═4fc2081c-0d5d-4fa6-88bf-9a26ef144249
# ╠═8b4cb5dd-602b-416d-ae77-6928819bd0ce
# ╠═225b2b60-d011-4340-810b-1f041b9552cd
# ╠═7ad8d849-41e7-43c7-863b-f76e9e3869c8
# ╠═c74d6489-b93a-4eb1-87a9-3b329ff6f527
# ╠═7bfb7485-5424-4c09-a0a9-4ed67b67fca9
# ╠═adce4d28-ee26-4bcf-98d2-06360dc82004
# ╠═d087bb85-0d8e-47f6-8519-f608db4df11f
# ╠═91093e20-1b08-4ff6-a270-2db085c1999a
