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

# ╔═╡ 9f75b286-b021-4fa1-a29d-7051c55c0a33
if !isdefined(Main, :PlutoRunner) # if running from file
	using ArgParse
	println("running as script")
	dirname2 = get_args()
end

# ╔═╡ c0df53fd-35d1-48aa-b1ba-131087be3aed
using HDF5

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 82c76c56-e874-4eba-9367-569b656155a2
pwd()

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
dirname1 = "/cosma/home/durham/dc-boye1/sculptor/isolation/1e5/"

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

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $dirname")
	out = lguys.Output("out")
end

# ╔═╡ c5672da9-0dad-4d22-abe5-9e186ccde02d
begin
	snap_i = out[1]
	snap_f = out[end]
end

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
begin 
	scatter(log10.(prof_i.rs), prof_i.Vs_circ * lguys.V0)
	scatter!(log10.(prof_f.rs), prof_f.Vs_circ * lguys.V0)

	xlabel!("log r")
	ylabel!("V_circ")
end

# ╔═╡ 30f3a998-8c4b-4cc3-97d7-967537f9986f
begin 
	scatter(log10.(prof_i.rs), log10.(prof_i.ρs))
	scatter!(log10.(prof_f.rs), log10.(prof_f.ρs))
	xlabel!("log r / kpc")
	ylabel!("log rho / 1e10 Msun kpc^-3")

end

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
histogram2d(snap_f.positions[1, :], snap_f.positions[2, :], bins = LinRange(-10, 10, 100))

# ╔═╡ 665ae91c-585c-4561-8f20-7541370cb837
# ╠═╡ disabled = true
#=╠═╡
begin 
	ps = lguys.extract(snap_i, :positions, p_filt)
	histogram2d(ps[1, :], ps[2, :], weights=probabilities)
end
  ╠═╡ =#

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap)
	pos = lguys.extract_vector(snap, :positions)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40)
	plot!(log10.(lguys.midpoint(r)), log10.(ρ))
end

# ╔═╡ 27f8deff-96ae-4d9a-a110-d146ac34965a
begin 
	function plot_ρ_dm(snap)
		plot()
		plot_ρ_dm!(snap)
	end

end

# ╔═╡ 41d9b58e-5cc9-4ec3-b696-db6b734b6ff7
range(0, 10, 13)

# ╔═╡ 012f3aa8-6e48-432a-b470-2390df3b5ec0
searchsortedfirst.([[1,2,3,5,6]], 5)

# ╔═╡ 4cc6cd49-1c65-415b-aa03-3452de1ee359
issorted(lguys.extract(snap_f, :index))

# ╔═╡ 195e5809-b80f-4005-9faa-ca792eec384a
cat([[1,2,3], [1,2,3]], dims=4)

# ╔═╡ ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
begin 
	Ls = hcat([lguys.calc_L_tot(snap) for snap in out]...)
	plot(transpose(Ls ./ Ls[:, 1]))
	xlabel!("snapshot")
	ylabel!("L/L_0")
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

# ╔═╡ 4984a413-8fb7-4937-8c5c-f7e424573866
cens = lguys.ss_centre(out)

# ╔═╡ 3cf46f2f-33c1-4eb5-b2af-243a0e89d1be
begin 
	
	x_cen = [cen.x_c for cen in cens]
	v_cen = [cen.v_c for cen in cens]
end

# ╔═╡ d6fe3195-4d4d-4db8-a384-a195233dc159
plot(transpose(hcat(x_cen...)))

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
plot(transpose(hcat(v_cen...)))

# ╔═╡ dc221349-eb61-4ace-8de3-a6c50249aca0
begin 
	rs = Vector[]
	for i in 1:length(out)
		push!(rs, sort(lguys.calc_r(out[i].positions .- x_cen[i])))
	end

	rs = hcat(rs...)
end

# ╔═╡ f21cfe22-95f3-485d-902b-b022a41548c2
begin 
	plot()
	for Npoints in [100, 300, 1000, 3000]
		plot!(1:length(out), log10.(rs[Npoints, :]), label="N=$Npoints")
	end
	xlabel!("snapshot")
	ylabel!("log r contating N particles")
end

# ╔═╡ a9033f1f-82cb-4e21-b708-e378f682c20c


# ╔═╡ bbfcd821-016e-47ff-a265-63784fc004af
begin 
	f = h5open("star_probabilities.hdf5")
	pidx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(pidx)]
	pidx = sort(pidx)
	close(f)
	p_filt = probabilities .> 0
	p_filt .&= map(!, isnan.(probabilities))
	pidx=pidx[p_filt]
	probabilities = probabilities[p_filt]
end

# ╔═╡ 106cbda4-57e0-459b-868b-b44339c944fc
begin 
	ps = lguys.extract_vector(snap_f, :positions, p_filt)
	radii = lguys.calc_r(ps)

end

# ╔═╡ f187611a-7147-4437-8cff-522e2a58cb71
	bins = 10 .^ LinRange(minimum(log10.(radii)), maximum(log10.(radii)), 30)


# ╔═╡ 0aa997e4-b256-4723-adc8-5a4c83a81427
length(ps)

# ╔═╡ 739d6b6c-40a0-4d51-9cca-ad771ff6f00a
ps

# ╔═╡ 9f6aa667-2478-4e96-9289-4d2c74bbf1c0
lguys.calc_ρ_hist(radii, weights=probabilities)

# ╔═╡ dd3a3e04-7b33-4f84-a935-0860283eca80
function plot_ρ_s!(snap)
	pos = lguys.extract_vector(snap, :positions, pidx)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=probabilities)
	plot!(log10.(lguys.midpoint(r)), log10.(ρ))
end

# ╔═╡ a8de1562-6afe-4002-b7bd-9d5541e8d354
begin 
	plot(xlim=(-2, 2), ylim=(-15, 5))
	plot_ρ_s!(snap_i)
	plot_ρ_s!(snap_f)

end

# ╔═╡ e0f6d5b1-1d85-4493-baa7-bd26caa5a3d5
length(probabilities)

# ╔═╡ f2ad545e-6b79-473d-bf6a-5392b7885501
begin 
	#bins, counts = lguys.calc_histogram(radii, bins, weights=probabilities)
	areas = 4π/3 * diff(bins .^ 3)
	ρ = counts ./ areas
	ρ_s_i = lguys.calc_ρ_hist(radii, 20, weights=probabilities)
end

# ╔═╡ 4a593099-44cc-4de0-bea2-2926623006b8
begin
	plot()
	scatter!(log10.(lguys.midpoint(bins)), log10.(ρ))
	scatter!(log10.(lguys.midpoint(bins)), log10.(ρ_s_i))
	#plot!(x -> -10 .^ (x/0.3))
end

# ╔═╡ ae1ef414-e010-4c79-930b-ab9069b93155
ρ

# ╔═╡ 42fdea8e-26b3-47d2-9195-655eba559a75
sum((probabilities) .<= 0)

# ╔═╡ 4d5013fd-9650-4739-818e-1b3e37d21f34
probabilities

# ╔═╡ 80f68ab3-1c01-4781-adc2-22ed21d335eb
sortperm(snap_i.index) * sortperm(pidx)

# ╔═╡ 890f9284-99b5-44b9-8832-b1df73e99b00
snap_i.index

# ╔═╡ 470466ff-3709-4cdb-93c0-b5b7848a940d
# ╠═╡ disabled = true
#=╠═╡
begin 
	anim = @animate for i in 1:1:length(out)
		snap = out[i]
		xr = 10
		plot(legend=false, grid=false, axis=false, dpi=100,
			xlim=(-xr, xr), ylim=(-xr, xr))
		scatter!(snap.positions[2, :], snap.positions[3, :],
		ms=1, msw=0, ma=1)
		scatter!([0], [0], ms=2, msw=0)
	
	end
	gif(anim, "isolation.gif", fps=12)
end
  ╠═╡ =#

# ╔═╡ 3b2bb553-0130-4c8a-80ad-6e1f7071a293
lguys.plot_xyz(lguys.extract(out, :positions, 2217))[3]

# ╔═╡ f7ec3ec0-7d3f-4dcb-b508-40f389723cd1
out[1].index[sortperm(lguys.calc_r(out[1]))][3000]

# ╔═╡ 6ba9c052-b786-4721-8947-ffb9e8f015de
	snap_i.masses = probabilities[inverse(sortperm(snap_i.index))]


# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╠═82c76c56-e874-4eba-9367-569b656155a2
# ╠═9f75b286-b021-4fa1-a29d-7051c55c0a33
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═f645be76-6477-4970-b655-603d700a10e7
# ╠═5717c18c-69cd-4740-a37a-d7ef80d68ae9
# ╠═5f7e3de9-a7fe-4217-a53c-0101d4b6314d
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═c5672da9-0dad-4d22-abe5-9e186ccde02d
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═30f3a998-8c4b-4cc3-97d7-967537f9986f
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╠═4e4c420e-6d35-42af-b859-b62861b1b6b9
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╠═528d0258-c218-4e7a-aa18-f096fdf45270
# ╠═665ae91c-585c-4561-8f20-7541370cb837
# ╠═106cbda4-57e0-459b-868b-b44339c944fc
# ╠═9f6aa667-2478-4e96-9289-4d2c74bbf1c0
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═dd3a3e04-7b33-4f84-a935-0860283eca80
# ╠═27f8deff-96ae-4d9a-a110-d146ac34965a
# ╠═a8de1562-6afe-4002-b7bd-9d5541e8d354
# ╠═e0f6d5b1-1d85-4493-baa7-bd26caa5a3d5
# ╠═41d9b58e-5cc9-4ec3-b696-db6b734b6ff7
# ╠═f187611a-7147-4437-8cff-522e2a58cb71
# ╠═f2ad545e-6b79-473d-bf6a-5392b7885501
# ╠═4a593099-44cc-4de0-bea2-2926623006b8
# ╠═ae1ef414-e010-4c79-930b-ab9069b93155
# ╠═0aa997e4-b256-4723-adc8-5a4c83a81427
# ╠═739d6b6c-40a0-4d51-9cca-ad771ff6f00a
# ╠═42fdea8e-26b3-47d2-9195-655eba559a75
# ╠═4d5013fd-9650-4739-818e-1b3e37d21f34
# ╠═012f3aa8-6e48-432a-b470-2390df3b5ec0
# ╠═4cc6cd49-1c65-415b-aa03-3452de1ee359
# ╠═195e5809-b80f-4005-9faa-ca792eec384a
# ╠═80f68ab3-1c01-4781-adc2-22ed21d335eb
# ╠═ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═fe476f91-5be3-4cf9-88df-55d9568ab88f
# ╠═4984a413-8fb7-4937-8c5c-f7e424573866
# ╠═3cf46f2f-33c1-4eb5-b2af-243a0e89d1be
# ╠═d6fe3195-4d4d-4db8-a384-a195233dc159
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═dc221349-eb61-4ace-8de3-a6c50249aca0
# ╠═f21cfe22-95f3-485d-902b-b022a41548c2
# ╠═a9033f1f-82cb-4e21-b708-e378f682c20c
# ╠═c0df53fd-35d1-48aa-b1ba-131087be3aed
# ╠═bbfcd821-016e-47ff-a265-63784fc004af
# ╠═890f9284-99b5-44b9-8832-b1df73e99b00
# ╠═470466ff-3709-4cdb-93c0-b5b7848a940d
# ╠═3b2bb553-0130-4c8a-80ad-6e1f7071a293
# ╠═f7ec3ec0-7d3f-4dcb-b508-40f389723cd1
# ╠═6ba9c052-b786-4721-8947-ffb9e8f015de
