### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

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

# ╔═╡ 7809e324-ba5f-4520-b6e4-c7727c227154
dirname1 = "/cosma/home/durham/dc-boye1/data/dwarfs/models/crater_ii/isolation/1e4/fiducial"

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

# ╔═╡ 855fd729-22b6-4845-9d2b-e796d4a15811
begin 
	# parameters 
	snapname = "./out/snapshot_000.hdf5"
	outname = "./star_probabilities.hdf5"
	
	r_s_s = 0.9 # kpc
	
	ρ_s(r) = exp((-r/r_s_s))
	
	Nr = 50
	NE = 100
end

# ╔═╡ 10569544-637d-4ab3-b193-c9b7bf7b3f3e
md"""
overwrite $outname ? $(@bind overwrite CheckBox())
"""

# ╔═╡ 4cb09115-143d-456f-9c6a-19656f638677
begin 
	println("$dirname")
	snap = lguys.Snapshot(snapname)
	cens = CSV.read("data/centres.csv", DataFrame)
end

# ╔═╡ 0f5671e6-deb4-11ee-3178-4d3f920f23a2
begin
	cen_i = cens[1, :]
	snap_i = lguys.copy(snap)
	snap_i.positions .-= [cen_i.x, cen_i.y, cen_i.z]
	snap_i.velocities .-= [cen_i.vx, cen_i.vy, cen_i.vz]


	radii = lguys.calc_r(snap_i.positions)
	snap_i = snap_i[sortperm(radii)]
	radii = sort(radii)
	vels = lguys.calc_r(snap_i.velocities)
	
	Φs = lguys.calc_radial_discrete_Φ(snap_i.masses, radii)
	ke = @. 1/2 * vels^2
	ϵs = @. -Φs - ke
	idx = ϵs .> 0
	
	snap_i = snap_i[idx]
	ke = ke[idx]
	radii = radii[idx]
	vels = vels[idx]
	ϵs = ϵs[idx]
	Φs = Φs[idx]

	Mins = cumsum(snap_i.masses) ./ sum(snap_i.masses)
end

# ╔═╡ 23158c79-d35c-411a-a35a-950f04214e19
begin 
	M_s_tot = 4π * quadgk(r-> r^2 * ρ_s(r), 0, Inf)[1]
	M_s(r) = 4π * quadgk(r-> r^2 * ρ_s(r) / M_s_tot, 0, r)[1]
end

# ╔═╡ 7a39cd4f-9646-4969-9410-b093bca633cb
begin 
	log_radii = log10.(radii)
	r_min = radii[1]
	r_max = radii[end]
	r_e = 10 .^ LinRange(log10.(r_min), log10.(r_max), Nr+1)
	r = lguys.midpoint(r_e)
	N_s_out = 1 - M_s(r_max)
	N_s_in = M_s(r_min)
	print("missing $N_s_out stars outside")
	print("missing $N_s_in stars inside")

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

# ╔═╡ 8bb8736d-a41b-4dac-a6cd-06d0d4704654
begin 
	plot(xlims=(-1.2, 3), ylims=(-15, 0),
	xlabel="log r", ylabel="log density")
	plot!(log10.(r), log10.(ν_dm), label="DM")
	plot!(log10.(r), log10.(ν_s), label="stars (analytic)")
end

# ╔═╡ b625d8a5-7265-4849-9bd6-ca8064d392eb
begin 
	plot(xlabel="log r / kpc", ylabel="ψ(r)")
	plot!(log10.(r), ψ, label="interpolated")
	plot!(log10.(radii), -Φs, label="snapshot")
end

# ╔═╡ 78ce5a98-fd3f-4e39-981f-2bea58b117bf
begin 
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, NE)
	f_dm_e = f_dm.(E)
	f_s_e = f_s.(E)
end

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
begin
	plot(xlabel="log ϵ", ylabel="log f", ylims=(-15, 5))
	plot!(log10.(E), nm.log10.(f_dm_e), label="DM")
	plot!(log10.(E), nm.log10.(f_s_e), label="stars")
end

# ╔═╡ 7409a024-cea4-47a5-84d2-846c96d88b7a
begin 
	probs = f_s_e ./ f_dm_e
	probs ./= sum(probs) # arbitrary scaling...
	prob = lguys.lerp(E, probs)
end

# ╔═╡ 84fdc265-988c-40db-87e5-44ba55d0e412
begin 
	plot(xlabel="log ϵ", ylabel="f_s / f_dm")
	plot!(log10.(E), f_s_e ./ f_dm_e / 50, label="probability")
	stephist!(log10.(ϵs), norm=:pdf, label="dark matter")
end

# ╔═╡ 3b229c8e-9320-4c07-b948-c34a0c082341
begin 
	ps = prob.(ϵs)
	print(sum(ps .< 0), " negative probabilities")
	ps[ps .< 0] .= 0
	ps ./= sum(ps)
end

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
begin 
	plot(xlabel="stellar weights", ylabel="frequency")

	histogram!((ps), yscale=:log10)
end

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

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
begin
	plot(xlabel=L"\log  \rm r / kpc", ylabel=L"\log \nu", ylim=(-10, 0), xlim=(-1.5, 2.3), fontfamily="Times")
	plot!(log10.(r), nm.log10.(ν_s), label="stars")
	plot!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody")
end

# ╔═╡ 27071738-f294-4e48-99f0-a08f2e369137
lguys.plot_centre(snap_i.positions, marker_z = ps, cmap=cgrad(:greys, rev=true), ms=1.5, width=15)

# ╔═╡ ffec02d2-32ed-4396-b5ff-2021724476d9
lguys.plot_centre(snap_i.positions, ms=1, width=10)

# ╔═╡ 7348b63e-43fa-4ff0-98dd-9de7f3167243
begin 
	idx0 = snap.index
	p_idx = zeros(length(idx0))
	idx_s = snap_i.index
	for i in 1:length(idx_s)
		j = findfirst(isequal(idx_s[i]), idx0 )
		p_idx[j] = ps[i]
	end
end

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
		
		h5write(outname, "/index", idx0)
		h5write(outname, "/probabilities", p_idx)
		println("probabilities written to $outname")
end

# ╔═╡ f42cb7e1-64b7-47da-be05-dd50c2471fb3
write_stars()

# ╔═╡ Cell order:
# ╠═a893932c-f184-42bc-9a0e-0960f10520aa
# ╠═641946b3-e6f2-4d6d-8777-7698f353eb3d
# ╠═7809e324-ba5f-4520-b6e4-c7727c227154
# ╟─10569544-637d-4ab3-b193-c9b7bf7b3f3e
# ╠═92cc30ae-0b0a-4a72-aa0c-29150eeee5e0
# ╠═855fd729-22b6-4845-9d2b-e796d4a15811
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╠═23158c79-d35c-411a-a35a-950f04214e19
# ╠═7a39cd4f-9646-4969-9410-b093bca633cb
# ╠═bd1bca1d-0982-47d8-823e-eadc05191b88
# ╠═1fab14ec-6bfc-4243-bc48-915d1a129925
# ╠═dfa675d6-aa32-45c3-a16c-626e16f36083
# ╠═8bb8736d-a41b-4dac-a6cd-06d0d4704654
# ╠═b625d8a5-7265-4849-9bd6-ca8064d392eb
# ╠═78ce5a98-fd3f-4e39-981f-2bea58b117bf
# ╠═75d23b44-71e7-4e28-ad3e-c537f3d4422f
# ╠═7409a024-cea4-47a5-84d2-846c96d88b7a
# ╠═84fdc265-988c-40db-87e5-44ba55d0e412
# ╠═3b229c8e-9320-4c07-b948-c34a0c082341
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
# ╠═33a26663-0c08-411b-902b-a509b2afa5ad
# ╠═a73db457-2cc2-4917-bb25-0429f4daecbd
# ╠═0e7a922f-4e61-4d02-ac28-3915f5d1c9da
# ╠═7d69e393-c4db-4ff0-ab5a-85ac50c785c2
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═27071738-f294-4e48-99f0-a08f2e369137
# ╠═ffec02d2-32ed-4396-b5ff-2021724476d9
# ╠═7348b63e-43fa-4ff0-98dd-9de7f3167243
# ╠═f42cb7e1-64b7-47da-be05-dd50c2471fb3
# ╠═bd8489da-ac53-46c6-979e-06d5dc6e25d1
