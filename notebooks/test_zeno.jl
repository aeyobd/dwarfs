### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 15d3d6a6-00ec-11ef-18de-2302746435df
begin
	import Pkg; Pkg.activate()

	import LilGuys as lguys
	using Plots, Arya
end

# ╔═╡ 4d4d3d2f-9d0a-4cbc-811b-4308616cdf33
begin 
	snap = lguys.Snapshot("../zeno/nfw_1e4.hdf5")
	snap = snap[sortperm(lguys.calc_r(snap))]
end

# ╔═╡ 8538b857-1ac8-4597-92ab-417f8dbd048d
radii = lguys.calc_r(snap)

# ╔═╡ a9245b79-f7e1-467a-b3a2-738355a682e4
M_s = 1; r_s=1

# ╔═╡ 14193d07-7acb-45ce-a536-185d50dbc54e
A(c) = log(1+c) - c/(1+c)

# ╔═╡ e8bdfcca-c4d0-422c-875a-8a225e22e1b3
M_in(r) = A(r/r_s) * M_s ./ A(1)

# ╔═╡ fc53c38e-be6f-45f2-acf9-cd9da1913522
A(1)

# ╔═╡ 2df25677-1747-4d8b-afc9-a05642af9701
begin 
	plot(log10.(radii), cumsum(snap.masses))
	plot!(log10.(radii), M_in.(radii))
end

# ╔═╡ Cell order:
# ╠═15d3d6a6-00ec-11ef-18de-2302746435df
# ╠═4d4d3d2f-9d0a-4cbc-811b-4308616cdf33
# ╠═8538b857-1ac8-4597-92ab-417f8dbd048d
# ╠═a9245b79-f7e1-467a-b3a2-738355a682e4
# ╠═14193d07-7acb-45ce-a536-185d50dbc54e
# ╠═e8bdfcca-c4d0-422c-875a-8a225e22e1b3
# ╠═fc53c38e-be6f-45f2-acf9-cd9da1913522
# ╠═2df25677-1747-4d8b-afc9-a05642af9701
