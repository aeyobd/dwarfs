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
	import Interpolations: linear_interpolation, Line
	import QuadGK: quadgk
	using Plots
	using DataFrames, CSV
	using NaNMath; nm = NaNMath
	using Arya
	using HDF5
	using PlutoUI
end

# ╔═╡ 7809e324-ba5f-4520-b6e4-c7727c227154
dirname1 = "/cosma/home/durham/dc-boye1/sculptor/isolation/1e5/"

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
	snapname = "./out/snapshot_001.hdf5"
	outname = "./star_probabilities.hdf5"
	r_s_s = 0.308 # kpc
	
	ρ_s(r) = exp((-r/r_s_s))
	
	Nr = 50
	NE = 100
end

# ╔═╡ 10569544-637d-4ab3-b193-c9b7bf7b3f3e
md"""
overwrite $outname ? $(@bind overwrite CheckBox())
"""

# ╔═╡ 8113e42a-19d3-4244-a1b1-7ee0aa7015ec
stars_rapha = CSV.read("stars.csv", DataFrame)

# ╔═╡ 4cb09115-143d-456f-9c6a-19656f638677
snap = lguys.Snapshot(snapname)

# ╔═╡ 4bcc1596-9332-4d56-b7d2-7aa52ca6392c
begin 
	p_rapha = stars_rapha[invperm(sortperm(snap.index)), :probs]
	all(stars_rapha[invperm(sortperm(snap.index)), :PartIDs] .== snap.index)
end

# ╔═╡ da0b82b9-d0e8-4135-a793-efd121fca3c3
cen = lguys.ss_centre(snap)

# ╔═╡ 0f5671e6-deb4-11ee-3178-4d3f920f23a2
begin
	snap_i = lguys.copy(snap)
	snap_i.positions .-= cen.x_c
	snap_i.velocities .-= cen.v_c


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

# ╔═╡ 290bd73f-8065-4d97-b94e-9f18fdd5d7e4
prof = lguys.Profile(snap_i)

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
	ψ = linear_interpolation(radii, -Φs)(r)
	
	As = 4π/3 * diff(r_e .^ 3)
	ν_dm = m_dm ./ As

	ν_s = ρ_s.(r) ./ M_s_tot

	M = linear_interpolation(radii, Mins)(r)
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
	plot(xscale=:log10, yscale=:log10, xlims=(1e-1, 1e3), ylims=(1e-15, 2),
	xlabel="r", ylabel="ν")
	plot!(r, ν_dm, label="DM")
	plot!(r, ν_s, label="stars")
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
	probs ./= maximum(probs)
	prob = linear_interpolation(E, probs, extrapolation_bc=Line())
end

# ╔═╡ 3b229c8e-9320-4c07-b948-c34a0c082341
begin 
	ps = prob.(ϵs)
	print(sum(ps .< 0))
	ps[ps .< 0] .= 0
end

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
histogram(ps, yscale=:log10)

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
begin 
	m_s = Vector{Float64}(undef, Nr)
	m_s_r = Vector{Float64}(undef, Nr)

	for i in 1:Nr
		filt = r_e[i] .<= radii .< r_e[i+1]
		m_s[i] = sum(ps[filt])
		m_s_r[i] = sum(p_rapha[filt])
	end
	m_s ./= sum(m_s)
	m_s_r ./= sum(m_s_r)
	ν_s_nbody = m_s ./ As
	ν_s_nbody_r = m_s_r ./ As
end

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
begin
	plot(xlabel=L"\log  \rm r / kpc", ylabel=L"\log \nu", ylim=(-15, 2), xlim=(-2, 1.5))
	plot!(log10.(r), nm.log10.(ν_s), label="stars")
	plot!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody")
	#plot!(log10.(r) , nm.log10.(ν_s_nbody_r), label="nbody")

end

# ╔═╡ b7d81ab1-9bf3-4976-9a3c-784ddbae85f7
begin
	w = 5
	filt = abs.(snap_i.positions[3, :]) .< w
	plot(xlims=(-w, w), ylims=(-w, w), xlabel="X", ylabel="Y")
	scatter!(snap_i.positions[1, filt], snap_i.positions[2, filt], marker_z=(ps[filt]), msw=0, ms=1, label="", cmap=cgrad(:greys, rev=true))
end

# ╔═╡ 74e18b36-c914-4df7-8743-1e0a1bb0c391
begin
	plot(xlims=(-w, w), ylims=(-w, w), xlabel="X", ylabel="Y")
	scatter!(snap_i.positions[1, filt], snap_i.positions[2, filt], msw=0, ms=1, label="", color=:black)
end

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
begin 
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
# ╠═8113e42a-19d3-4244-a1b1-7ee0aa7015ec
# ╠═4bcc1596-9332-4d56-b7d2-7aa52ca6392c
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╠═da0b82b9-d0e8-4135-a793-efd121fca3c3
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╠═290bd73f-8065-4d97-b94e-9f18fdd5d7e4
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
# ╠═3b229c8e-9320-4c07-b948-c34a0c082341
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═b7d81ab1-9bf3-4976-9a3c-784ddbae85f7
# ╠═74e18b36-c914-4df7-8743-1e0a1bb0c391
# ╠═7348b63e-43fa-4ff0-98dd-9de7f3167243
# ╠═f42cb7e1-64b7-47da-be05-dd50c2471fb3
# ╠═bd8489da-ac53-46c6-979e-06d5dc6e25d1
