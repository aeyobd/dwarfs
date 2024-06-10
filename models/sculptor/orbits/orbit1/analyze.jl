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

# ╔═╡ a227ef84-a874-4176-9253-36bc36a743fa
using LsqFit: curve_fit

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

# ╔═╡ a9e79439-16a4-4908-bfe0-f0770cdb26df
md"""
# Mass evolution
"""

# ╔═╡ dfa506e2-095a-47a2-9f14-baa04dd45d2d
function get_M_h(output::lguys.Output, radius; idxs=(1:10:length(output)))
	N = length(idxs)
	M = Vector{Float64}(undef, N)
	
	for i in eachindex(idxs)
		snap = output[idxs[i]]
		ϵ = lguys.calc_ϵ(snap)
		filt = ϵ .> 0
		
		rs = lguys.calc_r(snap[filt])
		filt2 = rs .< radius

		M[i] = sum(filt2)
	end

	return M
end

# ╔═╡ 807c46f6-56e0-4ab4-b31f-ac4a3fec9761
collect(1:4:10)

# ╔═╡ e4ffb45d-a7e2-4c57-b050-01483010fc4a
import YAML

# ╔═╡ 510706ac-ffbd-4996-af9e-67f1b910d51c
orbit_props = YAML.load_file("obital_properties.yml")

# ╔═╡ 53641449-c5b3-45ff-a692-a5cd717c8369
idx_f = orbit_props["idx_f"]

# ╔═╡ 7e3df305-9678-447e-a48e-f102cf6ebced
idx_i = 1

# ╔═╡ 9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
snap_i = out[idx_i]

# ╔═╡ 8d127679-401c-439d-913d-e2020df1c600
snap_f = out[idx_f]

# ╔═╡ 78d93a01-7048-4a97-b9ae-727be0e223d7
r_h = 0.11

# ╔═╡ 0fa11815-3ab0-4b19-9be7-186b7c2c1063
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="normalized mass within 100 rh",
		yscale=log10,
		yticks=[1, 0.1]
	)

	M_dm_h = get_M_h(out, Inf)
	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="all")

	
	M_dm_h = get_M_h(out, 100r_h)
	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="100")



	M_dm_h = get_M_h(out, 10r_h)
	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="10")

	M_dm_h = get_M_h(out, 1r_h)
	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1], label="1")
	
	Legend(fig[1, 2], ax)
	
	fig
end

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

# ╔═╡ db320665-f46d-4aed-a2b2-4b39bcb605c5
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

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log r / kpc", ylabel=L"\log\rho_\textrm{DM}",
	limits=(nothing, (-13, 0)))
	
	plot_ρ_dm!(snap_i, x_cen[:, idx_i], label="initial")
	plot_ρ_dm!(snap_f, x_cen[:, idx_f], label="final")
	
	axislegend(ax)
	#lines!([0, 0] .+ log10.(r_break), [-11.5, -10],  color=:black)
	#scatter!(log10.(r_break), -10, marker=:utriangle, color=:black)

	#text!(L"r_\textrm{break}", position=(log10.(r_break),-11.5), space=:data, rotation=π/2, align=(:left, :baseline))
	fig
	# only include bound points in profile...
end

# ╔═╡ 4801ff80-5761-490a-801a-b263b90d63fd
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

# ╔═╡ 1ebb1ab6-c1a0-4d6b-8529-65155575b96e
mlog10 = Makie.pseudolog10

# ╔═╡ fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
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

# ╔═╡ 7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
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

# ╔═╡ 743d20cd-9636-4384-9bf2-b2d7e259ae7d
# ╠═╡ disabled = true
#=╠═╡
r_h = 0.308
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═7094bc54-deb4-48a5-bf09-9ee6c684ac3c
# ╟─a9e79439-16a4-4908-bfe0-f0770cdb26df
# ╠═dfa506e2-095a-47a2-9f14-baa04dd45d2d
# ╠═807c46f6-56e0-4ab4-b31f-ac4a3fec9761
# ╠═e4ffb45d-a7e2-4c57-b050-01483010fc4a
# ╠═510706ac-ffbd-4996-af9e-67f1b910d51c
# ╠═53641449-c5b3-45ff-a692-a5cd717c8369
# ╠═7e3df305-9678-447e-a48e-f102cf6ebced
# ╠═9c3f79ee-89db-4fe1-aa62-4e706bdd73f8
# ╠═8d127679-401c-439d-913d-e2020df1c600
# ╠═78d93a01-7048-4a97-b9ae-727be0e223d7
# ╠═0fa11815-3ab0-4b19-9be7-186b7c2c1063
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
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═4801ff80-5761-490a-801a-b263b90d63fd
# ╠═1ebb1ab6-c1a0-4d6b-8529-65155575b96e
# ╠═fa9c08d6-98d1-46a4-a5d1-6cd79db77ace
# ╠═7c6f7fc7-e692-44a1-9ad0-a9377b0a5cdf
# ╠═743d20cd-9636-4384-9bf2-b2d7e259ae7d
