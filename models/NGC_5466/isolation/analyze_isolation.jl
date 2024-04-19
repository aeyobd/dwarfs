### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using Plots; gr()
	import LilGuys as lguys
end

# ╔═╡ 374489bc-627f-4fc9-9734-7c49456710ac
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 82c76c56-e874-4eba-9367-569b656155a2
pwd()

# ╔═╡ d9552869-9e46-4262-add5-1241710cbc8b
dirname = "fiducial"

# ╔═╡ f645be76-6477-4970-b655-603d700a10e7
begin 
	cd(@__DIR__)
	cd(dirname)
	pwd()
end

# ╔═╡ 5717c18c-69cd-4740-a37a-d7ef80d68ae9
plots_dir = "figures"

# ╔═╡ 5f7e3de9-a7fe-4217-a53c-0101d4b6314d
if plots_dir !== ""
	if !isdir(plots_dir)
		mkdir(plots_dir)
	end
	println("saving figures to $plots_dir")
else
	println("no directory specified")
end

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $dirname")
	out = lguys.Output("out")

end

# ╔═╡ 82517c9b-ff77-430e-bbed-239da814a851
begin
	println("loading centeres from $dirname")

	cens = CSV.read("data/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ 282484ff-a82c-413a-a391-cf45504a3c45
x_cen

# ╔═╡ c5672da9-0dad-4d22-abe5-9e186ccde02d
begin
	snap_i = out[1]
	idx_f = length(out)
	snap_f = out[idx_f]
end

# ╔═╡ 24937043-3513-4f29-ad90-4da70f761626
lguys.mean(snap_i.velocities .- v_cen[:, 1]) * lguys.V0 * 1e3 # m / s

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
begin 
	scatter(log10.(prof_i.rs), prof_i.Vs_circ * lguys.V0)
	scatter!(log10.(prof_f.rs), prof_f.Vs_circ * lguys.V0)

	xlabel!("log r")
	ylabel!("V_circ")
end

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
begin 
	histogram2d(log10.(lguys.calc_r(snap_i.positions, x_cen[:, 1])), lguys.calc_r(snap_i.velocities, v_cen[:, 1]) * lguys.V0)
	xlabel!("log radius / kpc")
	ylabel!("velocity (km/s)")
end

# ╔═╡ 4e4c420e-6d35-42af-b859-b62861b1b6b9
begin 
	histogram2d(log10.(lguys.calc_r(snap_f.positions, x_cen[:, idx_f])), lguys.calc_r(snap_f.velocities, v_cen[:, idx_f]) * lguys.V0)
	xlabel!("log radius / kpc")
	ylabel!("velocity (km/s)")
end

# ╔═╡ b9746093-0f2f-4478-82ba-00911c8fcceb
histogram2d(snap_i.positions[1, :], snap_i.positions[2, :], bins = LinRange(-0.1, 0.1, 100))

# ╔═╡ 528d0258-c218-4e7a-aa18-f096fdf45270
histogram2d(snap_f.positions[1, :], snap_f.positions[2, :], bins = LinRange(-0.1, 0.1, 100))

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap, x_c; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	rs = lguys.calc_r(pos, x_c)
	r, ρ = lguys.calc_ρ_hist(rs, 40)
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
	plot_ρ_dm!(snap_i, x_cen[:, 1],  label="initial")
	plot_ρ_dm!(snap_f, x_cen[:, idx_f], label="final")

end

# ╔═╡ 5915c34e-6562-4bb0-8411-231a9704aca0
r_h = 2.3 * 60 / 206265 * 16

# ╔═╡ a9095ec3-2fe3-4241-8bd0-fedc91aa0d4b
radii

# ╔═╡ ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
begin 
	Ls = hcat([lguys.calc_L_tot(snap) for snap in out]...)
	plot(transpose(Ls))
	xlabel!("snapshot")
	ylabel!("L/L_0")
end

# ╔═╡ bfa7593c-4915-4e03-83e6-8f790de4c1a5
begin 
	Es = hcat([lguys.calc_E_tot(out[i], v_cen[:, i]) for i in 1:length(out)]...)
	E_k = hcat([sum(lguys.calc_E_spec_kin(out[i], v_cen[:, i])) for i in 1:length(out)]...)

	E_phi = 1/2 * hcat([sum(out[i].Φs) for i in 1:length(out)]...)
	plot(transpose(Es ./ Es[1]), label="total")
	plot!(transpose(E_phi ./ E_phi[1]), label="grav")
	plot!(transpose(E_k ./ E_k[1]), label="potential")

	xlabel!("snapshot")
	ylabel!("rel change energy")
end

# ╔═╡ 6e98a5ce-ab2e-441b-b772-73df8e6655c2
E_phi ./ E_k

# ╔═╡ fe476f91-5be3-4cf9-88df-55d9568ab88f
pos = lguys.extract(snap_i, :positions)

# ╔═╡ d6fe3195-4d4d-4db8-a384-a195233dc159
plot(transpose(x_cen))

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
plot(transpose(v_cen) * lguys.V0)

# ╔═╡ 33fadf6d-23b7-40e7-b7ec-6863dc1980f6
begin 
	rs = Vector[]
	Ms = Vector[]
	rs_s = Vector[]

	percens = [0.01, .05, .10, .20, .30, .50]
	
	for i in 1:size(x_cen)[2]
		r = lguys.calc_r(out[i].positions .- x_cen[:, i])
		s_idx = sortperm(r)
		push!(rs, r[s_idx])
	end

	rs = hcat(rs...)
	rs_s = hcat(rs_s...)

end

# ╔═╡ e326570a-11d7-428f-90cf-520ddfce0fcf
lguys.percentile(radii, 50)

# ╔═╡ 4dff8e47-af0b-488e-90a0-1151a348d4fe
r_h

# ╔═╡ cf61f890-c512-432f-8dc3-160c1a00e793
Rs = @. sqrt(snap_i.positions[1, :]^2 + snap_i.positions[2, :]^2)

# ╔═╡ 802f1318-57ef-4870-a82e-c14aee9dd658
lguys.percentile(Rs, 50)

# ╔═╡ df627d8a-b2ca-49db-a22a-54c3c43c2aff
size(rs)

# ╔═╡ f21cfe22-95f3-485d-902b-b022a41548c2
begin 
	plot()
	for Npoints in round.(Int, length(out[1]) * percens)
		plot!(1:size(x_cen)[2], log10.(rs[Npoints, :]), label="N=$Npoints")
	end
	xlabel!("snapshot")
	ylabel!("log r contating N particles")
end

# ╔═╡ 470466ff-3709-4cdb-93c0-b5b7848a940d
begin 
	anim = @animate for i in 1:10:length(out)
		snap = out[i]
		xr = 10
		plot(legend=false, grid=false, axis=false, dpi=100,
			xlim=(-xr, xr), ylim=(-xr, xr))
		scatter!(snap.positions[2, :], snap.positions[3, :],
		ms=1, msw=0, ma=1)
		scatter!([0], [0], ms=2, msw=0)
	
	end
	gif(anim, "isolation.mp4", fps=12)
	anim
end

# ╔═╡ 3b2bb553-0130-4c8a-80ad-6e1f7071a293
lguys.plot_xyz_layout(1e3 * lguys.extract_vector(out, :positions, 3000)[:, 1:50])

# ╔═╡ 136692ab-39a6-4899-b95d-abe063beaed6
histogram(lguys.calc_r(snap_i.positions) * 1e3, xlim=(0, 30))

# ╔═╡ 6e14d250-eeca-4688-b212-2106f273c468
minimum(lguys.calc_r(snap_i.positions)) * 1e3

# ╔═╡ 08fdf8bc-4b63-4459-b0bd-6394ab9f13f8
plot([minimum(lguys.calc_r(out[i].positions, x_cen[:, i])) * 1e3 for i in 1:length(out)])

# ╔═╡ 19514763-a33b-48d1-8948-8a5ed5ea936a
out[1].positions[:, 2]

# ╔═╡ c60893bf-396e-4297-91e7-af726b73ec7b
lguys.percentile(lguys.calc_r(snap_i.positions), 0.1) * 1e3

# ╔═╡ 300791aa-ad1a-4640-8857-67f104c597a7
lguys.percentile(lguys.calc_r(snap_f.positions, x_cen[:, end]), 0.1) * 1e3

# ╔═╡ 665ae91c-585c-4561-8f20-7541370cb837
# ╠═╡ disabled = true
#=╠═╡
begin 
	ps = lguys.extract(snap_i, :positions, p_filt)
	histogram2d(ps[1, :], ps[2, :], weights=probabilities)
end
  ╠═╡ =#

# ╔═╡ 106cbda4-57e0-459b-868b-b44339c944fc
#=╠═╡
begin 
	ps = lguys.extract_vector(snap_f, :positions)
	radii = lguys.calc_r(ps)

end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═374489bc-627f-4fc9-9734-7c49456710ac
# ╠═82c76c56-e874-4eba-9367-569b656155a2
# ╠═d9552869-9e46-4262-add5-1241710cbc8b
# ╠═f645be76-6477-4970-b655-603d700a10e7
# ╠═5717c18c-69cd-4740-a37a-d7ef80d68ae9
# ╠═5f7e3de9-a7fe-4217-a53c-0101d4b6314d
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═82517c9b-ff77-430e-bbed-239da814a851
# ╠═282484ff-a82c-413a-a391-cf45504a3c45
# ╠═c5672da9-0dad-4d22-abe5-9e186ccde02d
# ╠═24937043-3513-4f29-ad90-4da70f761626
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╠═4e4c420e-6d35-42af-b859-b62861b1b6b9
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╠═528d0258-c218-4e7a-aa18-f096fdf45270
# ╠═665ae91c-585c-4561-8f20-7541370cb837
# ╠═106cbda4-57e0-459b-868b-b44339c944fc
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═27f8deff-96ae-4d9a-a110-d146ac34965a
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═5915c34e-6562-4bb0-8411-231a9704aca0
# ╠═a9095ec3-2fe3-4241-8bd0-fedc91aa0d4b
# ╠═ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═6e98a5ce-ab2e-441b-b772-73df8e6655c2
# ╠═fe476f91-5be3-4cf9-88df-55d9568ab88f
# ╠═d6fe3195-4d4d-4db8-a384-a195233dc159
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═33fadf6d-23b7-40e7-b7ec-6863dc1980f6
# ╠═e326570a-11d7-428f-90cf-520ddfce0fcf
# ╠═4dff8e47-af0b-488e-90a0-1151a348d4fe
# ╠═cf61f890-c512-432f-8dc3-160c1a00e793
# ╠═802f1318-57ef-4870-a82e-c14aee9dd658
# ╠═df627d8a-b2ca-49db-a22a-54c3c43c2aff
# ╠═f21cfe22-95f3-485d-902b-b022a41548c2
# ╠═470466ff-3709-4cdb-93c0-b5b7848a940d
# ╠═3b2bb553-0130-4c8a-80ad-6e1f7071a293
# ╠═136692ab-39a6-4899-b95d-abe063beaed6
# ╠═6e14d250-eeca-4688-b212-2106f273c468
# ╠═08fdf8bc-4b63-4459-b0bd-6394ab9f13f8
# ╠═19514763-a33b-48d1-8948-8a5ed5ea936a
# ╠═c60893bf-396e-4297-91e7-af726b73ec7b
# ╠═300791aa-ad1a-4640-8857-67f104c597a7
