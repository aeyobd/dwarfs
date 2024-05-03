### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ ecfd2d72-08e0-11ef-23b5-55f91aef7d76
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ 043e6913-7b92-4d55-962f-33c1d6e6b524
cd("/cosma/home/durham/dc-boye1/sculptor/orbits/orbit1")

# ╔═╡ 88535209-6ff9-45ee-90ef-4939bd79789c
begin 
	out =  lguys.Output("out/combined.hdf5")
	
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
	out.x_cen .= x_cen
	out.v_cen .= v_cen

	out
end

# ╔═╡ 751c780d-0d52-446e-a53f-fcf535819f6c
snap_i = out[1]

# ╔═╡ 57497573-5d20-4d26-8b5a-7af72580c59b
let
	fig = Figure()
	r_max = 150
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc", title="dark matter",
	limits=(-r_max, r_max, -r_max, r_max))

	bins = LinRange(-r_max, r_max, 300)
	colorrange=(1, 1e3)

	framerate = 10
	idxs = vcat(1:91, 92:10:length(out))

	hm = Arya.hist2d!(ax, snap_i.positions[2, :], snap_i.positions[3, :], bins = bins, colorscale=log10, colorrange=colorrange)
	
	record(fig, "sculptor.mp4", idxs, framerate = framerate) do i
		snap = out[i]
		x = snap.positions[2, :]
		y = snap.positions[3, :]
		H, x_e, y_e = Arya.histogram2d(x, y, bins)

		println(maximum(H))
		println(minimum(H[H .> 0]))
		
		hm[3] = H
	end

	fig
end

# ╔═╡ Cell order:
# ╠═ecfd2d72-08e0-11ef-23b5-55f91aef7d76
# ╠═043e6913-7b92-4d55-962f-33c1d6e6b524
# ╠═88535209-6ff9-45ee-90ef-4939bd79789c
# ╠═751c780d-0d52-446e-a53f-fcf535819f6c
# ╠═57497573-5d20-4d26-8b5a-7af72580c59b
