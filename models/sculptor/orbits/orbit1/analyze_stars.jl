### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 340ffbbe-17bd-11ef-35c6-63505bb128b7
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ f9fe37ef-de81-4d69-9308-cda968851ed2
begin 
	using HDF5

	f = h5open("star_probabilities.hdf5")
	p_idx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
end

# ╔═╡ c61ed7f6-593e-4d80-bcc8-645ecf70e213
Makie.save

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
begin 
	out =  lguys.Output("out/combined.hdf5")
	
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
	out.x_cen .= x_cen
	out.v_cen .= v_cen

	out
end

# ╔═╡ Cell order:
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═c61ed7f6-593e-4d80-bcc8-645ecf70e213
# ╠═f9fe37ef-de81-4d69-9308-cda968851ed2
# ╠═6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
