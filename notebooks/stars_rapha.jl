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
	using Plots
	using DataFrames, CSV
	using NaNMath; nm = NaNMath
	using Arya
	using HDF5
	using PlutoUI
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
dirname1 = "/cosma/home/durham/dc-boye1/data/dwarfs/models/sculptor/isolation/1e5/"

# ╔═╡ 855fd729-22b6-4845-9d2b-e796d4a15811
begin 
	# parameters 
	idx_i = 5
	idx_s = lpad(idx_i - 1,3,"0")
	snapname = "./out/snapshot_$idx_s.hdf5"
	outname = "./star_probabilities.hdf5"
	
	r_s_s = 0.109 # kpc
	
	ρ_s(r) = exp(-r/r_s_s)
	
	Nr = 99

	# number of energy bins, should not actually 
	# affect results, only changes the density f(E) is evaluated at
	NE = 100

	overwrite = true
end

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
	snap_i.positions .-= 0 * [cen_i.x, cen_i.y, cen_i.z]
	snap_i.velocities .-= 0 * [cen_i.vx, cen_i.vy, cen_i.vz]


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

# ╔═╡ 36b4adbd-d706-4e72-a922-53080c67946c
begin
	log_radii = log10.(radii)
	r_min = radii[1]
	r_max = radii[end]
	r_e = 10 .^ LinRange(log10.(r_min), log10.(r_max), Nr+1)
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
begin 
	plot(xlims=(-1.2, 3), ylims=(-15, 2),
	xlabel="log r", ylabel="log density")
	plot!(log10.(r), log10.(ν_dm), label="DM")
	plot!(log10.(r), log10.(ν_s), label="stars (analytic)")
end

# ╔═╡ 7481f47a-2b8a-45a3-9b4e-31ea14e9d331
md"""
The potential $\psi = -\Phi$ as a function of log radii (for the spherically calculated & interpolated and actual snapshot from Gadget 4)
"""

# ╔═╡ b625d8a5-7265-4849-9bd6-ca8064d392eb
begin 
	plot(xlabel="log r / kpc", ylabel="ψ(r)")
	plot!(log10.(r), ψ, label="interpolated")
	plot!(log10.(radii), -snap_i.Φs, label="snapshot")
end

# ╔═╡ 78ce5a98-fd3f-4e39-981f-2bea58b117bf
begin 
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, NE)
	f_dm_e = f_dm.(E)
	f_s_e = f_s.(E)
end

# ╔═╡ 8b66d00d-529b-4e8c-9386-b17117996579
md"""
The calculated distribution function as a function of log specific binding energy $\epsilon =  -\Phi - T$
"""

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
begin
	plot(xlabel="log ϵ", ylabel="log f", ylims=(-15, 7))
	plot!(log10.(E), nm.log10.(f_dm_e), label="DM")
	plot!(log10.(E), nm.log10.(f_s_e), label="stars")
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
begin 
	plot(xlabel="ϵ (binding energy)", ylabel="relative density", legend_position=:outertopright, size=(600, 300))
	plot!((E), 20*prob.(E) ./maximum(prob.(E)), label="f_s / f_dm")
	stephist!((ϵs), norm=:pdf, label="dark matter")
	stephist!((ϵs), weights=100ps, label="stars (nbody)")
	vline!([maximum(ϵs)], color="grey", ls=:dot, label=L"\epsilon_{\rm max}")
end

# ╔═╡ daf54cb4-06a3-4a8e-a533-354ca8740aec
md"""
A histogram of the assigned (positive) stellar weights. See the number of negative probabilities above
"""

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
begin 
	plot(xlabel="stellar weights", ylabel="frequency")
	histogram!(ps, nbins=100, yscale=:log10)
end

# ╔═╡ 5b9c7242-5d4b-4f5a-83e5-3e4d11017aa5
sum(ps .> 0)

# ╔═╡ b1a34b9b-c6ca-4818-95b2-5b55fb32511e
sum(ps .== 0)

# ╔═╡ a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
begin 
	plot(xlabel="log radii", ylabel="pdf (dn / dlog r)")
	stephist!(log10.(radii), norm=:pdf, label="dark matter")
	stephist!(log10.(radii), weights=ps, norm=:pdf, label="stars")
end

# ╔═╡ 33a26663-0c08-411b-902b-a509b2afa5ad
begin 
	plot(xlabel="log radii", ylabel="pstar > 0.025")

	histogram!(log10.(radii[ps .> 2e-4]))
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
begin
	p_nu_1 = plot(ylabel=L"\log \nu", ylim=(-4.8, 1.2), xlim=(-1.2, 0.3), fontfamily="Times")
	plot!(log10.(r), nm.log10.(ν_s), label="stars")
	plot!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody")

	p_nu_2 = plot(xlabel=L"\log  \rm r / kpc", ylabel=L"\Delta\log \nu ", ylim=(-1, 1), size=(600, 200), xlim=(-2.0, 0.7))
	scatter!(log10.(r), nm.log10.(ν_s_nbody) .- nm.log10.(ν_s), label="")
	hline!([0], color="black", label="", z_order=1)

	plot(p_nu_1, p_nu_2, layout=grid(2, 1, heights=[0.8, 0.2]))
end

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
begin 
	flt = ps .> 0
	histogram2d(lguys.get_x(snap_i)[flt], lguys.get_y(snap_i)[flt], weights=ps[flt], bins=LinRange(-r_hist, r_hist, N_hist), colorbar_scale=:log10, aspect_ratio=1, clim=(1e-8, 1), xlims=(-r_hist, r_hist))	
end

# ╔═╡ ee866164-c6f2-4f70-bde1-360abd5fd80e
md"""
Histogram of same region in dark matter only (over a smaller dynamic range). Dark matter is much more extended (as expected)
"""

# ╔═╡ a52e5e94-6068-4545-962f-e02a485b62f5
begin 
	histogram2d(lguys.get_x(snap_i), lguys.get_y(snap_i), bins=LinRange(-r_hist, r_hist, N_hist), colorbar_scale=:log10, aspect_ratio=1, xlims=(-r_hist, r_hist))	
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
# ╟─92cc30ae-0b0a-4a72-aa0c-29150eeee5e0
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╠═e5551ea6-7681-4544-a610-ddb5ee0add7b
# ╠═d7d83daf-6c70-4026-afa7-1f426ee805fa
# ╠═e37bf6d7-9445-49bf-8333-f68ad25436b2
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╠═f79414b4-620e-440e-a421-7bc13d373546
# ╠═23158c79-d35c-411a-a35a-950f04214e19
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
# ╠═ec46946a-0bf7-4374-af44-8d8b2ab6a3df
# ╟─06e1b872-ce52-434f-a8d1-3b0a5055eed2
# ╠═90856551-fff8-4d15-be66-6353091b5e50
# ╠═cc231e78-cfc0-4876-9f8b-980139f7d27f
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
