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
end

# ╔═╡ 917ecda1-1699-408c-9ff5-f64856184e02
lguys.scatter_xyz(snap.positions)

# ╔═╡ 5f50df86-e64b-4935-8a74-1c3991d7c432
lguys.scatter_xyz(snap.velocities)

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
	scatter!(log10.(r_c), log10.(ρs))
	plot!(log10.(r_c), log10.(1e11 * ρ_exp.(r_c/0.0015)))
end

# ╔═╡ d58eb47e-5c1e-481a-af9a-9d25f44b0916
sum(snap.masses) * lguys.M0

# ╔═╡ 488739be-5aec-4822-9bc0-b0627c0cd0ed
scatter(rs, lguys.calc_ϵ(snap), ms=1)

# ╔═╡ Cell order:
# ╠═fdf125d4-e5bd-11ee-2b47-15e50a226cc3
# ╠═dc9b9af8-26b9-4fd4-a886-f935454d4443
# ╠═917ecda1-1699-408c-9ff5-f64856184e02
# ╠═5f50df86-e64b-4935-8a74-1c3991d7c432
# ╠═71eec295-a292-4cce-8269-fd067309155e
# ╠═c18f7856-a12c-4ed0-832e-75d2c47a98c2
# ╠═f910e299-d85f-40a7-b787-bf327107816c
# ╠═d58eb47e-5c1e-481a-af9a-9d25f44b0916
# ╠═488739be-5aec-4822-9bc0-b0627c0cd0ed
