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

# ╔═╡ 4316d09d-43f6-4f6b-b830-0b8961743c61
if !@isdefined obs # read in sample (but only once)
	include("../../mc_orbits/sample.jl")
end

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

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	idx_i = 1
	idx_f = 1038
	snap_i = out[idx_i]
	snap_f = out[idx_f]
	
	"$(out.times[[idx_i, idx_f]] * lguys.T0) Gyr"
end
  ╠═╡ =#

# ╔═╡ a71bb30e-4a4c-4e94-b318-a426e3ee3045
begin 
	orbit_expected = CSV.read("orbit.csv", DataFrame)
	x_cen_exp = transpose(hcat(orbit_expected.x, orbit_expected.y, orbit_expected.z))
	v_cen_exp = -transpose(hcat(orbit_expected.v_x, orbit_expected.v_y, orbit_expected.v_z))

end

# ╔═╡ 29887611-5a0b-4f3a-8a3f-2da94b1765d2
lguys.plot_xyz(x_cen, x_cen_exp)

# ╔═╡ f563cbc7-655b-46e3-8686-2e4561b2467a
function plot_ρ_dm!(snap, x_cen; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	pos .-= x_cen
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

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
function fit_Vc(rc, Vc)
	filt = Vc .> percentile(Vc, 80)
	fit = curve_fit(Vc_model, rc[filt], Vc[filt], [2., 30.])
	return fit.param
end

# ╔═╡ 8547a5a2-7170-47a3-8110-cc2a5d49c070
function V_r_max(param)
	Rs, scale =param

	r = 2.16258 * Rs
	return r, Vc_model(r, param)
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
		param = fit_Vc(rc, Vc)
		r_max, V_max = V_r_max(param)

		Vh = get_Vh(snap, r_h)
		
		push!(V_maxs, V_max)
		push!(r_maxs, r_max)
		push!(V_hs, Vh)
	end

	return r_maxs, V_maxs, V_hs
end

# ╔═╡ 9418347c-e49f-4483-be0f-98a3b9578bac
r_max, V_max, V_h = get_V_r_max(out, skip=10)

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

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
#=╠═╡
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=(nothing, (-13, 0)))
	
	plot_ρ_dm!(snap_i, x_cen[:, idx_i], label="initial")
	plot_ρ_dm!(snap_f, x_cen[:, idx_f], label="final")
	
	axislegend(ax)
	fig
	# only include bound points in profile...
end
  ╠═╡ =#

# ╔═╡ d5316554-9f4d-45ea-93bf-f6d0f5eab4c2
function plot_ρ_s!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions, p_idx)
	rs = lguys.calc_r(pos, snap.x_cen)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=probabilities)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 4c331023-0777-43c7-90e5-3341aae0b141
#=╠═╡
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel = L"\log \rho_\star", 
		limits=((-1.9, 1), (-8, 2)))

	plot_ρ_s!(snap_i, label="initial")
	plot_ρ_s!(snap_f, label="final")
	
	axislegend(ax)
	fig
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

# ╔═╡ 1ebb1ab6-c1a0-4d6b-8529-65155575b96e
mlog10 = Makie.pseudolog10

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

# ╔═╡ 155222b2-f9df-4189-afbe-1af9d571d859
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 2cb67a4d-6941-4b9e-ae09-ad96f6fac51f
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 245721a6-01aa-43e7-922d-ed5da02207c1
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel="V_circ / km / s",
	limits=(nothing, (0, nothing)))
	x = out.times[1:10:end] * lguys.T0
	scatter!(x, V_max, label="max")
	scatter!(x, V_h, label=L"r=r_h")
	axislegend(ax)
	
	fig
end

# ╔═╡ 6497d973-d800-4052-a9b1-23f91dc3aa9f
begin 
	anim = @animate for i in 1:10:length(out)
		snap = out[i]
		plot(legend=false, grid=false, axis=false, dpi=100)
		scatter!(snap.positions[2, :], snap.positions[3, :], 
			ms=1, msw=0, ma=0.1, marker_z = m_star)
		plot!([0, 50], [-200, -200], color=:black)
		annotate!([25], [-200], ["50 kpc"], annotationfontsize=8, annotationhalign=:center, annotationvalign=:bottom)
		scatter!([0], [0], ms=2, msw=0)

		# scatter!([x_cen[2, i]], [x_cen[3,  i]], ms=1)
		xlims!(-200, 200)
		ylims!(-200, 200)
	end
	
	gif(anim, "sculptor.gif", fps = 12)
end

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

# ╔═╡ a1cc4e82-0b92-4fcf-9254-1ce788e408bb
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 2014a0f5-c99d-40de-9005-fa61ae9863b6
#=╠═╡
length(lguys.calc_ϵ(snap_i))
  ╠═╡ =#

# ╔═╡ aca25368-28c7-4f4a-bfcf-295902eaff9a
begin
	ρ_s(r, a, n) = exp( -(r/a)^(1/n))
end

# ╔═╡ 743d20cd-9636-4384-9bf2-b2d7e259ae7d
# ╠═╡ disabled = true
#=╠═╡
r_h = 0.308
  ╠═╡ =#

# ╔═╡ 7f168f93-849e-4586-a04f-6434165e6561
function gradient(x, y)
	s = sortperm(x)
	x = x[s]
	y = y[s]
	N = length(x)

	grad = Vector{Float64}(undef, N)

	grad[1] = (y[2] - y[1]) / (x[2] - x[1])
	grad[end] = (y[end] - y[end-1]) / (x[end] - x[end-1])
	for i in 2:(N-1)
		hs = x[i] - x[i-1]
		hd = x[i+1] - x[i]

		numerator = hs^2 * y[i+1] + (hd^2 - hs^2) * y[i] - hd^2*y[i-1]
		denom = hd*hs*(hd + hs)
		grad[i] = numerator/denom
	end
	return grad
end

# ╔═╡ b7629fac-2357-4e62-95ef-2ea12e92f335
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r_c = lguys.calc_r(x_cen)
	lines!(out.times * lguys.T0, r_c)
	fig
end

# ╔═╡ 42c677e3-8b25-4c8a-99ef-3abd0abcb891
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r_c = lguys.calc_r(x_cen)
	lines!(1:length(out), r_c)
	fig
end

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

# ╔═╡ f2381a9c-96d7-4935-8967-32b8bfb3fb48
#=╠═╡
begin 
	p_min = 1e-5
	snap_f_stars = snap_f[probabilities[snap_f.index] .> p_min]
	println(length(snap_f_stars))
	obs_pred = lguys.to_sky(snap_f_stars)
	m_star_f = probabilities[snap_f_stars.index]
end
  ╠═╡ =#

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

# ╔═╡ 68243ab7-0d82-48a2-8737-7cc105ae7325
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees"
	)

	x = [o.ra for o in obs_pred]
	y = [o.dec for o in obs_pred]
	
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)

	scatter!(obs.ra, obs.dec, color=:red)

	fig
end
  ╠═╡ =#

# ╔═╡ 3c25a197-b610-4bb2-b254-49d14534b343
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\mu_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\mu_\delta / \textrm{mas\,yr^{-1}}")

	x = [o.pm_ra for o in obs_pred]
	y = [o.pm_dec for o in obs_pred]
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)
	scatter!(obs.pm_ra, obs.pm_dec, color=:red)
	fig
end
  ╠═╡ =#

# ╔═╡ 270739f1-4dda-4c92-9fe1-371caf5dd322
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="distance / kpc",
		ylabel = "radial velocity / km/s")


	
	x = [o.distance for o in obs_pred]
	y = [o.radial_velocity for o in obs_pred]
	Arya.hist2d!(ax, x, y, weights=m_star_f, bins=100)

	scatter!(obs.distance, obs.radial_velocity, markersize=15, color=:red)

	fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═4316d09d-43f6-4f6b-b830-0b8961743c61
# ╠═7094bc54-deb4-48a5-bf09-9ee6c684ac3c
# ╠═60fe0995-3dd6-48ff-86e8-6fc5e099e39b
# ╠═f4d42909-5397-41e6-ac0c-07185bfece0e
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═29887611-5a0b-4f3a-8a3f-2da94b1765d2
# ╠═a71bb30e-4a4c-4e94-b318-a426e3ee3045
# ╠═f563cbc7-655b-46e3-8686-2e4561b2467a
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
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═d5316554-9f4d-45ea-93bf-f6d0f5eab4c2
# ╠═4c331023-0777-43c7-90e5-3341aae0b141
# ╠═4801ff80-5761-490a-801a-b263b90d63fd
# ╠═1ebb1ab6-c1a0-4d6b-8529-65155575b96e
# ╠═fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
# ╠═155222b2-f9df-4189-afbe-1af9d571d859
# ╠═2cb67a4d-6941-4b9e-ae09-ad96f6fac51f
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╠═6497d973-d800-4052-a9b1-23f91dc3aa9f
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
# ╠═a1cc4e82-0b92-4fcf-9254-1ce788e408bb
# ╠═2014a0f5-c99d-40de-9005-fa61ae9863b6
# ╠═aca25368-28c7-4f4a-bfcf-295902eaff9a
# ╠═743d20cd-9636-4384-9bf2-b2d7e259ae7d
# ╠═7f168f93-849e-4586-a04f-6434165e6561
# ╠═b7629fac-2357-4e62-95ef-2ea12e92f335
# ╠═42c677e3-8b25-4c8a-99ef-3abd0abcb891
# ╠═46f6e4a2-6e2b-4f23-807c-b41be51a0960
# ╠═f2381a9c-96d7-4935-8967-32b8bfb3fb48
# ╠═07c3bfc4-615e-4847-a66f-fb824f384ce8
# ╠═4962c654-9f9f-44b9-b0ee-68791522562a
# ╠═68243ab7-0d82-48a2-8737-7cc105ae7325
# ╠═3c25a197-b610-4bb2-b254-49d14534b343
# ╠═270739f1-4dda-4c92-9fe1-371caf5dd322
