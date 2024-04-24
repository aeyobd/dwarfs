### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using Plots; gr()
	import Arya

	import LilGuys as lguys
end

# ╔═╡ 374489bc-627f-4fc9-9734-7c49456710ac
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 9e9d463b-feaa-41bd-96c7-9c320d933b71
using QuadGK

# ╔═╡ 9f75b286-b021-4fa1-a29d-7051c55c0a33
if !isdefined(Main, :PlutoRunner) # if running from file
	using ArgParse
	println("running as script")
	dirname2 = get_args()
end

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 82c76c56-e874-4eba-9367-569b656155a2
pwd()

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Inputs
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
dirname1 = "/cosma/home/durham/dc-boye1/data/dwarfs/models/sculptor/isolation/1e6"

# ╔═╡ 79b07d75-fb05-4833-ac2c-ea0e9c24e791
r_s_s = 0.109 # stellar scale radius

# ╔═╡ ef3d1d5f-0979-44a7-8f0f-bf4638ea5612
ρ_s(r) = 1/8π /r_s_s^3 * exp(-r / r_s_s)

# ╔═╡ 199c29c1-7325-41e2-91e8-854dba902930
quadgk(r->4π*r^2 * ρ_s(r), 0, Inf)

# ╔═╡ f645be76-6477-4970-b655-603d700a10e7
begin 
	if isdefined(Main, :PlutoRunner)
		dirname = dirname1
	else
		dirname = dirname2
	end
	cd(@__DIR__)
	cd(dirname)
	pwd()
end

# ╔═╡ 5717c18c-69cd-4740-a37a-d7ef80d68ae9
plots_dir = "$dirname/figures"

# ╔═╡ 5f7e3de9-a7fe-4217-a53c-0101d4b6314d
if dirname !== ""
	if !isdir(dirname)
		mkdir(plots_dir)
	end
	println("saving figures to $plots_dir")
else
	println("no directory specified")
end

# ╔═╡ 2fb86841-fc96-47e1-af57-898fa2690ff3
begin
	println("$dirname")


	f = h5open("star_probabilities.hdf5")
	pidx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(pidx)]
	pidx = sort(pidx)
	close(f)

end

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $dirname")
	out = lguys.Output("out/combined.hdf5")

	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ 97e98ab8-b60b-4b48-b465-a34a16858f88
md"""
# Initial/Final
"""

# ╔═╡ c5672da9-0dad-4d22-abe5-9e186ccde02d
begin
	snap_i = out[1]
	snap_i.x_cen = x_cen[:, 1]
	snap_i.v_cen = v_cen[:, 1]

	snap_f = out[end - 1]
	snap_f.x_cen = x_cen[:, end]
	snap_f.v_cen = v_cen[:, end]

end

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
# ╠═╡ disabled = true
#=╠═╡
begin 
	plot()

	rc, Vc = lguys.calc_V_circ(snap_i)
	scatter!(log10.(rc), Vc * lguys.V0, label="initial")

	rc, Vc = lguys.calc_V_circ(snap_f)
	scatter!(log10.(rc), Vc * lguys.V0, label="final")
	xlabel!("log r")
	ylabel!("V_circ")
end
  ╠═╡ =#

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
begin 
	histogram2d(log10.(lguys.calc_r(snap_i.positions)), lguys.calc_r(snap_i.velocities) * lguys.V0)
	xlabel!("log radius / kpc")
	ylabel!("velocity (km/s)")
end

# ╔═╡ 4e4c420e-6d35-42af-b859-b62861b1b6b9
begin 
	histogram2d(log10.(lguys.calc_r(snap_f.positions)), lguys.calc_r(snap_f.velocities) * lguys.V0)
	xlabel!("log radius / kpc")
	ylabel!("velocity (km/s)")
end

# ╔═╡ b9746093-0f2f-4478-82ba-00911c8fcceb
histogram2d(snap_i.positions[1, :], snap_i.positions[2, :], bins = LinRange(-10, 10, 100))

# ╔═╡ 528d0258-c218-4e7a-aa18-f096fdf45270
histogram2d(snap_f.positions[1, :] .- x_cen[1, end:end], snap_f.positions[2, :] .- x_cen[2, end:end], bins = LinRange(-10, 10, 100))

# ╔═╡ 29391db0-a45e-4edb-9c5c-7cfc8c141551
length(probabilities)

# ╔═╡ 9f6aa667-2478-4e96-9289-4d2c74bbf1c0
lguys.calc_ρ_hist(radii, weights=probabilities)

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40)
	plot!(log10.(lguys.midpoint(r)), log10.(ρ); xlabel="log r", ylabel="log ρ", kwargs...)
end

# ╔═╡ dd3a3e04-7b33-4f84-a935-0860283eca80
function plot_ρ_s!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions, pidx)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=probabilities)
	plot!(log10.(lguys.midpoint(r)), log10.(ρ); xlabel="log r", ylabel="log ρ", kwargs...)
end

# ╔═╡ 27f8deff-96ae-4d9a-a110-d146ac34965a
begin 
	function plot_ρ_dm(snap)
		plot()
		plot_ρ_dm!(snap)
	end

end

# ╔═╡ 60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
begin 
	plot()
	plot_ρ_dm!(snap_i, label="initial")
	plot_ρ_dm!(snap_f, label="final")

end

# ╔═╡ a8de1562-6afe-4002-b7bd-9d5541e8d354
begin 
	plot(xlim=(-2, 2), ylim=(-12, 3), title="stars")
	plot_ρ_s!(snap_i, label="initial")
	plot_ρ_s!(snap_f, label="final")
	vline!([log10(r_s_s)], label="r_s")

	log_r_pred = LinRange(-2, 2, 1000)
	ρ_s_pred = ρ_s.(10 .^ log_r_pred)

	plot!(log_r_pred, log10.(ρ_s_pred), label="expected")

end

# ╔═╡ e33d56a7-7a0e-4fa9-8f0d-041b43584d59
sum(lguys.calc_r(snap_i) .< r_s_s)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
histogram(log10.(lguys.calc_r(snap_i)), ylabel="count", xlabel="log r / kpc", label="")

# ╔═╡ fb0dec74-aaab-43a4-9b37-d13634c5dcad
md"""
# Time variation
"""

# ╔═╡ ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
begin 
	Ls = hcat([lguys.calc_L_tot(snap) for snap in out]...)
	plot(transpose(Ls))
	xlabel!("snapshot")
	ylabel!("L")
end

# ╔═╡ bfa7593c-4915-4e03-83e6-8f790de4c1a5
begin 
	Es = hcat([lguys.calc_E_tot(snap) for snap in out]...)
	plot(transpose(Es ./ Es[1]))
	xlabel!("snapshot")
	ylabel!("rel change energy")
end

# ╔═╡ fe476f91-5be3-4cf9-88df-55d9568ab88f
pos = lguys.extract(snap_i, :positions)

# ╔═╡ 70f9aa47-956c-4db3-ad3a-f18a18318d33
lguys.plot_xyz_layout(x_cen)

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
plot(transpose(v_cen) * lguys.V0)

# ╔═╡ dc221349-eb61-4ace-8de3-a6c50249aca0
begin 
	rs = Vector[]
	Ms = Vector[]
	rs_s = Vector[]

	percens = [0.003, 0.01, .05, .10, .20, .30, .50, .9]
	
	for i in 1:length(out)
		r = lguys.calc_r(out[i].positions .- x_cen[:, i])
		s_idx = sortperm(r)
		push!(rs, r[s_idx])
		
		
		ms = cumsum(probabilities[out[i].index[s_idx]])
		ms ./= ms[end]
		idxs = searchsortedfirst.([ms], percens)
		if i % 50 == 0
			println(idxs)
		end
		push!(Ms, ms)
		push!(rs_s, r[s_idx][idxs])
	end

	rs = hcat(rs...)
	rs_s = hcat(rs_s...)

end

# ╔═╡ f21cfe22-95f3-485d-902b-b022a41548c2
begin 
	plot(legend_position=:outertopright)
	for Npoints in round.(Int, length(out[1]) * percens)
		plot!(1:length(out), log10.(rs[Npoints, :]), label="f=$(round(Npoints/length(out[1]), digits=3))")
	end
	xlabel!("snapshot")
	ylabel!("log r contating fraction of particles")
end

# ╔═╡ 38aeb93b-8b79-4880-b0d1-cef180d13bc3
begin 
	plot(dpi=600, legend=false)
	for i in eachindex(percens)
		plot!(1:length(out), log10.(rs_s[i, :]), label="$(percens[i])")
	end

	xlabel!("snapshot")
	ylabel!("log r contating stellar mass")
end

# ╔═╡ 967136d3-8d58-4fdc-9537-aa3a85a92528
out.times * lguys.T0

# ╔═╡ 470466ff-3709-4cdb-93c0-b5b7848a940d
begin 
	anim = @animate for i in 1:10:length(out)
		snap = out[i]
		xr = 10
		plot(legend=false, grid=false, axis=false, dpi=100,
			xlim=(-xr, xr), ylim=(-xr, xr), fmt=:gif)
		scatter!(snap.positions[2, :], snap.positions[3, :],
		ms=1, msw=0, ma=1)
		scatter!([0], [0], ms=2, msw=0)
	
	end
	gif(anim, "isolation.gif", fps=12) #anim
end

# ╔═╡ 3b2bb553-0130-4c8a-80ad-6e1f7071a293
lguys.plot_xyz_layout(lguys.extract_vector(out, :positions, 4012))

# ╔═╡ ce8420db-b395-4db2-acb2-d900281add38
begin 
	w = 10
	plot(xlims=(-w, w), ylims=(-w, w))

	idx_p = 100
	snap_p = out[idx_p]

	filt_p = snap_p.positions[3, :] .- x_cen[3, idx_p] .< w
	scatter!(snap_p.positions[1, filt_p] .- x_cen[1, idx_p], snap_p.positions[2, filt_p] .- x_cen[2, idx_p], marker_z=asinh.(probabilities[snap_p.index][filt_p] ./ 0.01), ms=1, msw=0, label="")
end

# ╔═╡ a9fe4a34-f8ce-46bd-84fb-bdae17a508a4
histogram(probabilities, bins=20, yscale=:log10, xlabel="probability", ylabel="count")

# ╔═╡ 665ae91c-585c-4561-8f20-7541370cb837
# ╠═╡ disabled = true
#=╠═╡
begin 
	ps = lguys.extract(snap_i, :positions, p_filt)
	histogram2d(ps[1, :], ps[2, :], weights=probabilities)
end
  ╠═╡ =#

# ╔═╡ 106cbda4-57e0-459b-868b-b44339c944fc
begin 
	ps = lguys.extract_vector(snap_f, :positions)
	radii = lguys.calc_r(ps)

end

# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═374489bc-627f-4fc9-9734-7c49456710ac
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╠═82c76c56-e874-4eba-9367-569b656155a2
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═79b07d75-fb05-4833-ac2c-ea0e9c24e791
# ╠═ef3d1d5f-0979-44a7-8f0f-bf4638ea5612
# ╠═9e9d463b-feaa-41bd-96c7-9c320d933b71
# ╠═199c29c1-7325-41e2-91e8-854dba902930
# ╠═5717c18c-69cd-4740-a37a-d7ef80d68ae9
# ╠═f645be76-6477-4970-b655-603d700a10e7
# ╠═9f75b286-b021-4fa1-a29d-7051c55c0a33
# ╠═5f7e3de9-a7fe-4217-a53c-0101d4b6314d
# ╠═2fb86841-fc96-47e1-af57-898fa2690ff3
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╟─97e98ab8-b60b-4b48-b465-a34a16858f88
# ╠═c5672da9-0dad-4d22-abe5-9e186ccde02d
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╠═4e4c420e-6d35-42af-b859-b62861b1b6b9
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╠═528d0258-c218-4e7a-aa18-f096fdf45270
# ╠═665ae91c-585c-4561-8f20-7541370cb837
# ╠═106cbda4-57e0-459b-868b-b44339c944fc
# ╠═29391db0-a45e-4edb-9c5c-7cfc8c141551
# ╠═9f6aa667-2478-4e96-9289-4d2c74bbf1c0
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═dd3a3e04-7b33-4f84-a935-0860283eca80
# ╠═27f8deff-96ae-4d9a-a110-d146ac34965a
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═a8de1562-6afe-4002-b7bd-9d5541e8d354
# ╠═e33d56a7-7a0e-4fa9-8f0d-041b43584d59
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╟─fb0dec74-aaab-43a4-9b37-d13634c5dcad
# ╠═ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═fe476f91-5be3-4cf9-88df-55d9568ab88f
# ╠═70f9aa47-956c-4db3-ad3a-f18a18318d33
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═dc221349-eb61-4ace-8de3-a6c50249aca0
# ╠═f21cfe22-95f3-485d-902b-b022a41548c2
# ╠═38aeb93b-8b79-4880-b0d1-cef180d13bc3
# ╠═967136d3-8d58-4fdc-9537-aa3a85a92528
# ╠═470466ff-3709-4cdb-93c0-b5b7848a940d
# ╠═3b2bb553-0130-4c8a-80ad-6e1f7071a293
# ╠═ce8420db-b395-4db2-acb2-d900281add38
# ╠═a9fe4a34-f8ce-46bd-84fb-bdae17a508a4
