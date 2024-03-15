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
	import Arya
	using HDF5
	using PlutoUI
end

# ╔═╡ 855fd729-22b6-4845-9d2b-e796d4a15811
begin 
	# parameters 
	snapname = "./out/snapshot_001.hdf5"
	outname = "./star_probabilities.hdf5"
	r_s_s = 0.308 # kpc
	
	ρ_s(r) = exp((-r/r_s_s))
	
	Nr = 100
	NE = 100
end

# ╔═╡ 4cb09115-143d-456f-9c6a-19656f638677
snap = lguys.Snapshot(snapname)

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

# ╔═╡ 8bb8736d-a41b-4dac-a6cd-06d0d4704654
begin 
	plot(xscale=:log10, yscale=:log10, xlims=(1e-1, 1e3), ylims=(1e-10, 1),
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

# ╔═╡ 4d06de8b-70c8-4f3a-bd68-4d4cb569913e
begin 
	dν_dψ_dm = lguys.gradient(ψ, ν_dm)
	d2ν_dψ2_dm = lguys.gradient(ψ, dν_dψ_dm)
	dν_dψ_s = lguys.gradient(ψ, ν_s)
	d2ν_dψ2_s = lguys.gradient(ψ, dν_dψ_s)
	d2ν_dψ2_s = ifelse.(isnan.(d2ν_dψ2_s), 0, d2ν_dψ2_s)

	d2ν_dψ2_dm_i = linear_interpolation(reverse(ψ), reverse(d2ν_dψ2_dm), extrapolation_bc=Line())
	d2ν_dψ2_s_i = linear_interpolation(reverse(ψ), reverse(d2ν_dψ2_s), extrapolation_bc=Line())
end

# ╔═╡ 45d0f929-de54-4a62-b4bd-2fdcaf8594b3
plot(log10.(r), dν_dψ_s)

# ╔═╡ a537f12e-52f1-45b1-b6ac-0af8f392236a
plot(ψ, d2ν_dψ2_dm)

# ╔═╡ c12d34fb-4c4c-4974-8cf1-7ba212c4b42d
less_than_one = M .< 1 .- 1 ./ collect(1:Nr)

# ╔═╡ 41917a4f-453f-41da-bab9-9074021a5f7e
begin 
	f_s(ϵ) = 1/√8π^2 * quadgk(ψi -> d2ν_dψ2_s_i(ψi) / sqrt(ϵ - ψi), 0, ϵ)[1]
	f_dm(ϵ) = 1/√8π^2 * quadgk(ψi -> d2ν_dψ2_dm_i(ψi) / sqrt(ϵ - ψi), 0, ϵ)[1]
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
	plot(xlabel="log ϵ", ylabel="log f DM")
	plot!(log10.(E), nm.log10.(f_dm_e), label="")
end

# ╔═╡ 0d973031-5352-4d34-967a-c5beae114532
plot(log10.(E), nm.log10.(f_s_e), ylims=(-10, 5), xlabel="log E", ylabel="log f stars", legend=false)

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
stephist(ps, yscale=:log10)

# ╔═╡ 1bce7ae0-5081-4e29-bb72-cd0e937cee7e
Arya.init()

# ╔═╡ 14957bf0-4a3c-42fe-983c-7a71a07ad155
ps

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
begin 
	m_s = Vector{Float64}(undef, Nr)

	for i in 1:Nr
		filt = r_e[i] .<= radii .< r_e[i+1]
		m_s[i] = sum(ps[filt])
	end
	m_s ./= sum(m_s)
	ν_s_nbody = m_s ./ As
end

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
begin
	plot(xlabel="log r / kpc", ylabel="log ν", ylim=(-10, 0), xlim=(-1, 1))
	plot!(log10.(r), nm.log10.(ν_s), label="exponential")
	plot!(log10.(r), nm.log10.(ν_s_nbody), label="nbody reconstruction")

end

# ╔═╡ b7d81ab1-9bf3-4976-9a3c-784ddbae85f7
begin
	w = 0.5
	plot(xlims=(-w, w), ylims=(-w, w), xlabel="X", ylabel="Y")
	filt_in = abs.(snap_i.positions[3, :]) .< w
	scatter!(snap_i.positions[1, filt_in], snap_i.positions[2, filt_in], marker_z=asinh.(ps[filt_in] *1e3), msw=0, ms=3)
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

# ╔═╡ 10569544-637d-4ab3-b193-c9b7bf7b3f3e
md"""
overwrite $outname ? $(@bind overwrite CheckBox())
"""

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
# ╠═855fd729-22b6-4845-9d2b-e796d4a15811
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╠═da0b82b9-d0e8-4135-a793-efd121fca3c3
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╠═290bd73f-8065-4d97-b94e-9f18fdd5d7e4
# ╠═23158c79-d35c-411a-a35a-950f04214e19
# ╠═7a39cd4f-9646-4969-9410-b093bca633cb
# ╠═bd1bca1d-0982-47d8-823e-eadc05191b88
# ╠═dfa675d6-aa32-45c3-a16c-626e16f36083
# ╠═8bb8736d-a41b-4dac-a6cd-06d0d4704654
# ╠═b625d8a5-7265-4849-9bd6-ca8064d392eb
# ╠═4d06de8b-70c8-4f3a-bd68-4d4cb569913e
# ╠═45d0f929-de54-4a62-b4bd-2fdcaf8594b3
# ╠═a537f12e-52f1-45b1-b6ac-0af8f392236a
# ╠═c12d34fb-4c4c-4974-8cf1-7ba212c4b42d
# ╠═41917a4f-453f-41da-bab9-9074021a5f7e
# ╠═78ce5a98-fd3f-4e39-981f-2bea58b117bf
# ╠═75d23b44-71e7-4e28-ad3e-c537f3d4422f
# ╠═0d973031-5352-4d34-967a-c5beae114532
# ╠═7409a024-cea4-47a5-84d2-846c96d88b7a
# ╠═3b229c8e-9320-4c07-b948-c34a0c082341
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═1bce7ae0-5081-4e29-bb72-cd0e937cee7e
# ╠═14957bf0-4a3c-42fe-983c-7a71a07ad155
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═b7d81ab1-9bf3-4976-9a3c-784ddbae85f7
# ╠═7348b63e-43fa-4ff0-98dd-9de7f3167243
# ╟─10569544-637d-4ab3-b193-c9b7bf7b3f3e
# ╠═f42cb7e1-64b7-47da-be05-dd50c2471fb3
# ╠═bd8489da-ac53-46c6-979e-06d5dc6e25d1
