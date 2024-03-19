### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ fdf125d4-e5bd-11ee-2b47-15e50a226cc3
begin 
	import Pkg; Pkg.activate()

	import LilGuys as lguys
	using Plots
end

# ╔═╡ dc9b9af8-26b9-4fd4-a886-f935454d4443
begin 
	snap = lguys.Snapshot("isolation.hdf5")
	snapr = lguys.Snapshot("/cosma/home/durham/dc-boye1/unions1_rapha/model.hdf5")
	cen = lguys.ss_centre(snap)
	cenr = lguys.ss_centre(snapr)

	snap.positions .-= cen.x_c
	snap.velocities .-= cen.v_c
	snapr.positions .-= cenr.x_c
	snapr.velocities .-= cenr.v_c
end

# ╔═╡ 917ecda1-1699-408c-9ff5-f64856184e02
histogram2d(snap.positions[1, :], snap.positions[2, :], xlabel="x/kpc", ylabel="y/kpc")

# ╔═╡ 5f50df86-e64b-4935-8a74-1c3991d7c432
histogram2d(snap.velocities[1, :] * lguys.V0, snap.velocities[2, :] * lguys.V0, xlabel="vx/km s⁻¹", ylabel="vy/km s⁻¹", bins=40)

# ╔═╡ fe96b9e9-1dab-4b42-90eb-4f6bafd2e230
begin 
	prof = lguys.Profile(snap, zeros(3), zeros(3), Nr=100)
	prof_r = lguys.Profile(snapr, zeros(3), zeros(3), Nr=100)
end

# ╔═╡ 71eec295-a292-4cce-8269-fd067309155e
begin 
	rs = lguys.calc_r(snap)
	bins, cs = lguys.calc_histogram(log10.(rs))
	r_e = 10 .^ bins
	r_c = lguys.midpoint(r_e)
	areas = 4π/3 * diff(r_e .^ 3)
	ρs = cs ./ areas
end

# ╔═╡ c18f7856-a12c-4ed0-832e-75d2c47a98c2
function ρ_exp(r)
	return exp(-r)
end

# ╔═╡ f910e299-d85f-40a7-b787-bf327107816c
begin 
	plot()
	x = LinRange(-3.5, -1.8, 100)
	r_h = 0.003/2
	M_h = 16 / lguys.M0
	ρ0 = M_h / r_h ^ 3

	plot!(x, log10(ρ0) .+ log10.(ρ_exp.(10 .^ x / r_h)))
	scatter!(log10.(prof.rs), log10.(prof.ρs))
	scatter!(log10.(prof_r.rs),  log10.(prof_r.ρs))

end

# ╔═╡ d58eb47e-5c1e-481a-af9a-9d25f44b0916
sum(snap.masses) * lguys.M0

# ╔═╡ 488739be-5aec-4822-9bc0-b0627c0cd0ed
scatter(rs, lguys.calc_ϵ(snap), ms=1)

# ╔═╡ 0ac8f1e2-fd2e-4d47-9db6-9ab4fec7e3ef
begin 
	plot()
	scatter!(log10.(prof.rs), (prof.Vs_circ) * lguys.V0)
	scatter!(log10.(prof_r.rs), (prof_r.Vs_circ) * lguys.V0)
	xlabel!("log r / kpc")
	ylabel!(" circular velocity")

end

# ╔═╡ b50d5b35-c2df-4235-bf2a-2c71cd78565a
begin
	rs1 = lguys.calc_r(snap)
	vs1 = lguys.calc_r(snap.velocities)

	rs_r = lguys.calc_r(snapr)
	vs_r = lguys.calc_r(snapr.velocities)

	ϕ1 = lguys.calc_radial_discrete_Φ(snap)
	ϕ_r = lguys.calc_radial_discrete_Φ(snapr)

	ϵ1 = -ϕ1 .- 1/2 * vs1 .^ 2
	ϵ_r = -ϕ_r .- 1/2 * vs_r .^ 2

end

# ╔═╡ 5b3bc89a-a581-4c01-bdfd-db1512129888
begin 
	scatter(log10.(rs1), ϵ1, ms=1, msw=0)
	scatter!(log10.(rs_r), ϵ_r, ms=1, msw=0)
end

# ╔═╡ b3893bba-b5f7-4479-8573-86e4068e9014
begin 
	plot(xlabel="log r", ylabel="log v")
	scatter!(log10.(rs1), (vs1), ms=1, msw=0)
	scatter!(log10.(rs_r), (vs_r), ms=1, msw=0)
end

# ╔═╡ Cell order:
# ╠═fdf125d4-e5bd-11ee-2b47-15e50a226cc3
# ╠═dc9b9af8-26b9-4fd4-a886-f935454d4443
# ╠═917ecda1-1699-408c-9ff5-f64856184e02
# ╠═5f50df86-e64b-4935-8a74-1c3991d7c432
# ╠═fe96b9e9-1dab-4b42-90eb-4f6bafd2e230
# ╠═71eec295-a292-4cce-8269-fd067309155e
# ╠═c18f7856-a12c-4ed0-832e-75d2c47a98c2
# ╠═f910e299-d85f-40a7-b787-bf327107816c
# ╠═d58eb47e-5c1e-481a-af9a-9d25f44b0916
# ╠═488739be-5aec-4822-9bc0-b0627c0cd0ed
# ╠═0ac8f1e2-fd2e-4d47-9db6-9ab4fec7e3ef
# ╠═b50d5b35-c2df-4235-bf2a-2c71cd78565a
# ╠═5b3bc89a-a581-4c01-bdfd-db1512129888
# ╠═b3893bba-b5f7-4479-8573-86e4068e9014
