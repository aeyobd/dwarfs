### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ a893932c-f184-42bc-9a0e-0960f10520aa
begin
	import Pkg
	Pkg.activate()
	
	import LilGuys as lguys
end

# ╔═╡ 641946b3-e6f2-4d6d-8777-7698f353eb3d
begin 
	import QuadGK: quadgk
	using CairoMakie
	using DataFrames, CSV
	using NaNMath; nm = NaNMath
	using Arya
	using HDF5
	using Revise
end

# ╔═╡ 17ffde4b-5796-4915-9741-d594cf0c5ca7
md"""
# Daniel's script for painting stars
based on Rapha's codes. 
Eddington inversion

Given the directory, the script changes to there and works there
"""

# ╔═╡ 93045024-a91d-4b31-9a5a-7c999afdb9ec
md"""
# Inputs
"""

# ╔═╡ 7809e324-ba5f-4520-b6e4-c7727c227154
dirname1 = "/cosma/home/durham/dc-boye1/data/dwarfs/models/sculptor/isolation/1e6/"

# ╔═╡ 855fd729-22b6-4845-9d2b-e796d4a15811
begin 
	# parameters 
	r_s_s = 0.11
	ρ_s(r) = exp(-r/r_s_s)
	Nr = 150
	
	NE = 1000
	
	idx_i = 20
	#include(dirname1 * "/star_params.jl")
	idx_s = lpad(idx_i - 1,3,"0")
	snapname = "./out/snapshot_$idx_s.hdf5"
	#snapname = "../initial.hdf5"
	outname = "./star_probabilities.hdf5"


	overwrite = true
end

# ╔═╡ 9ea8bae1-bdd4-4987-94cc-599df21f23ef
NE

# ╔═╡ 92cc30ae-0b0a-4a72-aa0c-29150eeee5e0
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

# ╔═╡ 4cb09115-143d-456f-9c6a-19656f638677
begin 
	println("$dirname")
	snap = lguys.Snapshot(snapname)
	cens = CSV.read("out/centres.csv", DataFrame)
end

# ╔═╡ e5551ea6-7681-4544-a610-ddb5ee0add7b
cens[idx_i, :]

# ╔═╡ d7d83daf-6c70-4026-afa7-1f426ee805fa
cens[idx_i, :].t * lguys.T0

# ╔═╡ 0f5671e6-deb4-11ee-3178-4d3f920f23a2
begin
	# centre snapshot
	cen_i = cens[idx_i, :]
	snap_i = lguys.copy(snap)
	snap_i.positions #.-= 0 * [cen_i.x, cen_i.y, cen_i.z]
	snap_i.velocities #.-= 0 * [cen_i.vx, cen_i.vy, cen_i.vz]


	# sort by radii
	radii = lguys.calc_r(snap_i.positions)

	snap_i = snap_i[sortperm(radii)]
	radii = sort(radii)
	vels = lguys.calc_r(snap_i.velocities)

	# calculate energies
	Φs = lguys.calc_radial_discrete_Φ(snap_i.masses, radii)
	ke = @. 1/2 * vels^2
	ϵs = @. -Φs - ke

	# filter snapshot
	idx = ϵs .> 0
	idx_excluded = snap_i.index[map(!, idx)]
	
	snap_i = snap_i[idx]
	radii = radii[idx]
	vels = vels[idx]
	
	ke = ke[idx]
	ϵs = ϵs[idx]
	Φs = Φs[idx]
end

# ╔═╡ f79414b4-620e-440e-a421-7bc13d373546
Mins = cumsum(snap_i.masses) ./ sum(snap_i.masses)

# ╔═╡ 23158c79-d35c-411a-a35a-950f04214e19
begin 
	M_s_tot = 4π * quadgk(r-> r^2 * ρ_s(r), 0, Inf)[1]
	M_s(r) = 4π * quadgk(r-> r^2 * ρ_s(r) / M_s_tot, 0, r)[1]
end

# ╔═╡ b3663c19-c029-4a1e-ab82-a05177e3a5d0
import StatsBase: percentile

# ╔═╡ 36b4adbd-d706-4e72-a922-53080c67946c
begin
	log_radii = log10.(radii)
	r_min = radii[1]
	r_max = radii[end]
	r_e = 10 .^ LinRange(log10.(r_min), log10.(r_max), Nr+1)
	r_e = percentile(radii, LinRange(0, 100, Nr+1))
	r = lguys.midpoint(r_e)
end

# ╔═╡ e37bf6d7-9445-49bf-8333-f68ad25436b2
r_e

# ╔═╡ 7a39cd4f-9646-4969-9410-b093bca633cb
begin 
	# just a informative note
	N_s_out = 1 - M_s(r_max)
	N_s_in = M_s(r_min)
	println("missing $N_s_out stars outside")
	println("missing $N_s_in stars inside")
end

# ╔═╡ bd1bca1d-0982-47d8-823e-eadc05191b88
begin 
	m_dm = Vector{Float64}(undef, Nr)
	for i in 1:Nr
		filt = r_e[i] .<= radii .< r_e[i+1]
		m_dm[i] = sum(snap_i.masses[filt])
	end
end

# ╔═╡ dfa675d6-aa32-45c3-a16c-626e16f36083
begin 
	ψ = lguys.lerp(radii, -Φs).(r)
	
	ν_dm = lguys.calc_ρ_hist(radii, r_e)[2]
	ν_dm ./= length(snap)

	ν_s = ρ_s.(r) ./ M_s_tot

	M = lguys.lerp(radii, Mins).(r)
	Ms = M_s.(r)
end

# ╔═╡ 1fab14ec-6bfc-4243-bc48-915d1a129925
begin 
	ψ_i = lguys.lerp(radii, -Φs).(r)
	f_dm = lguys.calc_fϵ(ν_dm, ψ_i, r)
	f_s = lguys.calc_fϵ(ν_s, ψ_i, r)
end

# ╔═╡ 126c6825-723f-4d13-b5a3-64cba72fc867
md"""
The density distribution of the stars (analytic) and dark matter (calculated) from the snapsho
"""

# ╔═╡ 8bb8736d-a41b-4dac-a6cd-06d0d4704654
let
	fig = Figure()
	ax = Axis(fig[1,1], #limits=(-1.2, 3, -15, 2),
		xlabel="log r", ylabel="log density")
	lines!(log10.(r), log10.(ν_dm), label="DM")
	lines!(log10.(r), log10.(ν_s), label="stars (analytic)")

	fig
end

# ╔═╡ 7481f47a-2b8a-45a3-9b4e-31ea14e9d331
md"""
The potential $\psi = -\Phi$ as a function of log radii (for the spherically calculated & interpolated and actual snapshot from Gadget 4)
"""

# ╔═╡ b625d8a5-7265-4849-9bd6-ca8064d392eb
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Psi(r)")
	lines!(log10.(r), ψ, label="interpolated")
	lines!(log10.(radii), -snap_i.Φs, label="snapshot")

	fig
end

# ╔═╡ 78ce5a98-fd3f-4e39-981f-2bea58b117bf
begin 
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, NE + 1)
	f_dm_e = f_dm.(E)
	f_s_e = f_s.(E)
end

# ╔═╡ 8b66d00d-529b-4e8c-9386-b17117996579
md"""
The calculated distribution function as a function of log specific binding energy $\epsilon =  -\Phi - T$
"""

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="log ϵ", ylabel="log f", limits=(nothing, (-15, 7)) )
	lines!(log10.(E), nm.log10.(f_dm_e), label="DM")
	lines!(log10.(E), nm.log10.(f_s_e), label="stars")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 7409a024-cea4-47a5-84d2-846c96d88b7a
begin 
	probs = f_s_e ./ f_dm_e
	probs ./= sum(probs .* lguys.gradient(E)) # pdf, dN/dE
	prob = lguys.lerp(E, probs)
end

# ╔═╡ a389fde4-f05f-48ff-9a8c-8eb0643a849d
ϵs

# ╔═╡ 3b229c8e-9320-4c07-b948-c34a0c082341
begin 
	ps = prob.(ϵs)
	print(sum(ps .< 0), " negative probabilities")
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0
	ps ./= sum(ps)
end

# ╔═╡ 84fdc265-988c-40db-87e5-44ba55d0e412
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"$\epsilon$ (binding energy)", 
		ylabel="relative density", 
		)
	
	lines!((E), 20*prob.(E) ./maximum(prob.(E)), color=Arya.COLORS[3], label="f_s / f_dm")
	stephist!((ϵs), normalization=:pdf, label="dark matter")
	stephist!((ϵs), weights=100ps, label="stars (nbody)")
	vlines!([maximum(ϵs)], color="grey", ls=:dot, label=L"\epsilon_\textrm{max}")

	axislegend(ax, position=:lt)

	fig
end

# ╔═╡ daf54cb4-06a3-4a8e-a533-354ca8740aec
md"""
A histogram of the assigned (positive) stellar weights. See the number of negative probabilities above
"""

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="stellar weights", ylabel="frequency", yscale=log10)
	hist!(ps, bins=100)
	fig
end

# ╔═╡ 5b9c7242-5d4b-4f5a-83e5-3e4d11017aa5
sum(ps .> 0)

# ╔═╡ b1a34b9b-c6ca-4818-95b2-5b55fb32511e
sum(ps .== 0)

# ╔═╡ a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radii", ylabel="pdf (dn / dlog r)")
	stephist!(log10.(radii), normalization=:pdf, label="dark matter")
	stephist!(log10.(radii), weights=ps, normalization=:pdf, label="stars")
	fig
end

# ╔═╡ 33a26663-0c08-411b-902b-a509b2afa5ad
let
	fig = Figure()
	Axis(fig[1,1], xlabel="log radii", ylabel="pstar > 0.025")

	hist!(log10.(radii[ps .> 2e-4]))

	fig
end

# ╔═╡ a73db457-2cc2-4917-bb25-0429f4daecbd
length(radii)

# ╔═╡ 0e7a922f-4e61-4d02-ac28-3915f5d1c9da
length(ps)

# ╔═╡ 7d69e393-c4db-4ff0-ab5a-85ac50c785c2
maximum(ps)

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
begin 
	m_s = Vector{Float64}(undef, Nr)

	for i in 1:Nr
		filt = r_e[i] .<= radii .< r_e[i+1]
		m_s[i] = sum(ps[filt])
	end
	m_s ./= sum(m_s)
	ν_s_nbody = lguys.calc_ρ_hist(radii, r_e, weights=ps)[2]
end

# ╔═╡ b67a17dc-7a59-4b62-9120-4c2ada12298c
md"""
The main result. The reconstructed density profile
"""

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
let
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=(-1, 0.8, -15, 3)
		)
	lines!(log10.(r), nm.log10.(ν_s), label="stars")
	lines!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody")

	ax2 = Axis(fig[2,1], 
		xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Delta\log \nu ", 
		limits=(nothing, (-1, 1)))
	
	scatter!(log10.(r), nm.log10.(ν_s_nbody) .- nm.log10.(ν_s), label="")
	hlines!([0], color="black", label="", z_order=1)

	linkxaxes!(ax, ax2, )
	rowsize!(fig.layout, 2, Auto(0.3))
	hidexdecorations!(ax, grid=false)
	fig
end

# ╔═╡ 73752ab9-9a1d-4d30-bf95-376f894cafc3


# ╔═╡ b2452dfe-f3af-4a4f-a5da-8301cbf3fbf1
R = sqrt.(lguys.get_x(snap_i) .^2 .+  lguys.get_y(snap_i) .^2)

# ╔═╡ 131bd303-4f0f-4b24-bd90-5c65cf342f4c
radii

# ╔═╡ 2fe5b16c-c8fd-4d7a-899e-43ee06a98722
r_h_recon = round(R[findfirst(x->x>0.5, cumsum(ps))], digits=3)

# ╔═╡ ee3e2851-860f-450b-85ff-611262577373
md"""
reconstructed half light radius = $r_h_recon kpc
"""

# ╔═╡ 52bf618a-5494-4cfb-9fd0-ee82f9682116
radii[findfirst(x->x>0.5, cumsum(ps))]

# ╔═╡ ec46946a-0bf7-4374-af44-8d8b2ab6a3df
ps

# ╔═╡ 06e1b872-ce52-434f-a8d1-3b0a5055eed2
md"""
A histogram of the stellar density
"""

# ╔═╡ 90856551-fff8-4d15-be66-6353091b5e50
begin 
	r_hist = 5
	N_hist = 100
end

# ╔═╡ cc231e78-cfc0-4876-9f8b-980139f7d27f
let
	fig = Figure()
	filt = ps .> 0
	
	ax = Axis(fig[1,1], limits=(-r_hist, r_hist, -r_hist, r_hist), aspect=1)
	
	Arya.hist2d!(ax, 
		lguys.get_x(snap_i)[filt], lguys.get_y(snap_i)[filt], 
		weights=ps[filt], bins=LinRange(-r_hist, r_hist, N_hist))	
	fig
end

# ╔═╡ 9f179cf4-b9c2-4065-b1e2-16b389fd2e9d
eltype(ones(Int64, 10))

# ╔═╡ ee866164-c6f2-4f70-bde1-360abd5fd80e
md"""
Histogram of same region in dark matter only (over a smaller dynamic range). Dark matter is much more extended (as expected)
"""

# ╔═╡ a52e5e94-6068-4545-962f-e02a485b62f5
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(-r_hist, r_hist, -r_hist, r_hist), aspect=1)
	Arya.hist2d!(ax, lguys.get_x(snap_i), lguys.get_y(snap_i), bins=LinRange(-r_hist, r_hist, N_hist))	
	fig
end

# ╔═╡ a3071af2-beff-408d-b109-d4f289f8f7f4
md"""
## Last checks and saving profile
"""

# ╔═╡ a483545a-42ea-474b-a2ca-6f1bd4c7275b
p_idx = snap_i.index

# ╔═╡ 7bc02274-aa44-4609-b9a6-e409de5172af
begin
	idx_all = vcat(p_idx, idx_excluded)
	ps_all = vcat(ps, zeros(length(idx_excluded)))
end

# ╔═╡ f66a4a6d-b31c-481b-bf7a-4801c783ceb4
sort(idx_all) == sort(snap.index) # check we didn't lose any star particles

# ╔═╡ 92798e2e-2356-43c2-a4a9-82a70619a5f5
sum(ps_all[idx_excluded]) # should be zero

# ╔═╡ 4671b864-470d-4181-993b-4e64d5687460
sum(ps_all) # should be 1

# ╔═╡ bd8489da-ac53-46c6-979e-06d5dc6e25d1
function write_stars()
	if isfile(outname)
		if overwrite
			rm(outname)
		else
			println("file already exists")
			return
		end
	end


	h5write(outname, "/index", idx_all)
	h5write(outname, "/probabilities", ps_all)
	println("probabilities written to $outname")
end

# ╔═╡ f42cb7e1-64b7-47da-be05-dd50c2471fb3
write_stars()

# ╔═╡ Cell order:
# ╟─17ffde4b-5796-4915-9741-d594cf0c5ca7
# ╠═a893932c-f184-42bc-9a0e-0960f10520aa
# ╠═641946b3-e6f2-4d6d-8777-7698f353eb3d
# ╟─93045024-a91d-4b31-9a5a-7c999afdb9ec
# ╠═7809e324-ba5f-4520-b6e4-c7727c227154
# ╠═855fd729-22b6-4845-9d2b-e796d4a15811
# ╠═9ea8bae1-bdd4-4987-94cc-599df21f23ef
# ╟─92cc30ae-0b0a-4a72-aa0c-29150eeee5e0
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╠═e5551ea6-7681-4544-a610-ddb5ee0add7b
# ╠═d7d83daf-6c70-4026-afa7-1f426ee805fa
# ╠═e37bf6d7-9445-49bf-8333-f68ad25436b2
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╠═f79414b4-620e-440e-a421-7bc13d373546
# ╠═23158c79-d35c-411a-a35a-950f04214e19
# ╠═b3663c19-c029-4a1e-ab82-a05177e3a5d0
# ╠═36b4adbd-d706-4e72-a922-53080c67946c
# ╠═7a39cd4f-9646-4969-9410-b093bca633cb
# ╠═bd1bca1d-0982-47d8-823e-eadc05191b88
# ╠═1fab14ec-6bfc-4243-bc48-915d1a129925
# ╠═dfa675d6-aa32-45c3-a16c-626e16f36083
# ╟─126c6825-723f-4d13-b5a3-64cba72fc867
# ╠═8bb8736d-a41b-4dac-a6cd-06d0d4704654
# ╟─7481f47a-2b8a-45a3-9b4e-31ea14e9d331
# ╠═b625d8a5-7265-4849-9bd6-ca8064d392eb
# ╠═78ce5a98-fd3f-4e39-981f-2bea58b117bf
# ╟─8b66d00d-529b-4e8c-9386-b17117996579
# ╠═75d23b44-71e7-4e28-ad3e-c537f3d4422f
# ╠═7409a024-cea4-47a5-84d2-846c96d88b7a
# ╠═84fdc265-988c-40db-87e5-44ba55d0e412
# ╠═a389fde4-f05f-48ff-9a8c-8eb0643a849d
# ╠═3b229c8e-9320-4c07-b948-c34a0c082341
# ╟─daf54cb4-06a3-4a8e-a533-354ca8740aec
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═5b9c7242-5d4b-4f5a-83e5-3e4d11017aa5
# ╠═b1a34b9b-c6ca-4818-95b2-5b55fb32511e
# ╠═a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
# ╠═33a26663-0c08-411b-902b-a509b2afa5ad
# ╠═a73db457-2cc2-4917-bb25-0429f4daecbd
# ╠═0e7a922f-4e61-4d02-ac28-3915f5d1c9da
# ╠═7d69e393-c4db-4ff0-ab5a-85ac50c785c2
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╟─b67a17dc-7a59-4b62-9120-4c2ada12298c
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═73752ab9-9a1d-4d30-bf95-376f894cafc3
# ╠═ee3e2851-860f-450b-85ff-611262577373
# ╠═b2452dfe-f3af-4a4f-a5da-8301cbf3fbf1
# ╠═131bd303-4f0f-4b24-bd90-5c65cf342f4c
# ╠═2fe5b16c-c8fd-4d7a-899e-43ee06a98722
# ╠═52bf618a-5494-4cfb-9fd0-ee82f9682116
# ╠═ec46946a-0bf7-4374-af44-8d8b2ab6a3df
# ╟─06e1b872-ce52-434f-a8d1-3b0a5055eed2
# ╠═90856551-fff8-4d15-be66-6353091b5e50
# ╠═cc231e78-cfc0-4876-9f8b-980139f7d27f
# ╠═9f179cf4-b9c2-4065-b1e2-16b389fd2e9d
# ╟─ee866164-c6f2-4f70-bde1-360abd5fd80e
# ╠═a52e5e94-6068-4545-962f-e02a485b62f5
# ╟─a3071af2-beff-408d-b109-d4f289f8f7f4
# ╠═a483545a-42ea-474b-a2ca-6f1bd4c7275b
# ╠═7bc02274-aa44-4609-b9a6-e409de5172af
# ╠═f66a4a6d-b31c-481b-bf7a-4801c783ceb4
# ╠═92798e2e-2356-43c2-a4a9-82a70619a5f5
# ╠═4671b864-470d-4181-993b-4e64d5687460
# ╠═bd8489da-ac53-46c6-979e-06d5dc6e25d1
# ╠═f42cb7e1-64b7-47da-be05-dd50c2471fb3
