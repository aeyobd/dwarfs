### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 1ecbfa5e-8439-11ee-3a96-fb3d2dc7a2d3
begin 
	import Pkg
	Pkg.activate()
	using Plots
	plotly()
	using DataFrames, CSV

	import LilGuys as lguys
end

# ╔═╡ ed6ec0e9-5734-4569-b214-2c86c22d3c55
plots = Ref{Vector{Plots.Plot}}() # so we can reuse this variable

# ╔═╡ cd9dea36-93f8-4461-85b1-87030d9367bb
pwd()

# ╔═╡ 9c7e5bb4-8db0-4527-a8ec-331e2aed958b
begin
	out = lguys.Output("out");
	out2 = lguys.Output("../mc_orbits/out")
end

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
begin
	function read_peris(filename)
		df = CSV.read(filename, DataFrame)
		df_idx = sortperm(df[!, :index])
		peris = df[df_idx, :pericenter]
		apos = df[df_idx, :apocenter];
		return peris, apos
	end

	peris, apos = read_peris("peris_apos.csv")
	peris2, apos2 = read_peris("../mc_orbits/peris_apos.csv")

end

# ╔═╡ fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
begin
	plots[] = []
	p1 = histogram(peris)
	histogram!(peris2)
	xlabel!("pericenter")
	push!(plots[], p1)

	p1 = histogram(apos)
	xlabel!("apocenter")
	histogram!(apos2)
	push!(plots[], p1)

	plots[]
end

# ╔═╡ 90b730b2-d64d-42c5-9aa6-31ee1c0927fe
idx = 1

# ╔═╡ a644a33e-57a4-4af4-98be-1e0e84948069
begin
	gp_orbit = CSV.read("/cosma/home/durham/dc-boye1/dwarfs/notebooks/sculptor_orbits.csv", DataFrame)
	positions3 = transpose(hcat(gp_orbit.x, gp_orbit.y, gp_orbit.z))
end

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = lguys.extract(out, :positions, idx)
	positions2 = lguys.extract(out2, :positions, idx)
	velocities = lguys.extract(out, :velocities, idx)
	velocities2 = lguys.extract(out2, :velocities, idx)
	accelerations = lguys.extract(out, :accelerations, idx)
	accelerations2 = lguys.extract(out2, :accelerations, idx)

	rs = lguys.calc_r(positions)
	rs2 = lguys.calc_r(positions2)
	rs3 = lguys.calc_r(Matrix(positions3))

	vs = lguys.calc_r(velocities)
	a_s = lguys.calc_r(accelerations)
end

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.plot_xyz(positions, positions2, positions3)

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities, velocities2)

# ╔═╡ 127c1123-29cf-45a7-b7ab-7f57b4d6e58f
lguys.plot_xyz(accelerations, accelerations2)

# ╔═╡ 2e9817a4-906d-4a0d-867c-d301cf3f2b1e
begin 
	plot(out.times, rs, label="gadget 4")
	plot!(out2.times, rs2, label="gadget 2")
	plot!(gp_orbit.t / lguys.T0 * 1e9, rs3, label="galpy")
	xlabel!("time (code units)")
	ylabel!("radius (kpc)")
end

# ╔═╡ Cell order:
# ╠═1ecbfa5e-8439-11ee-3a96-fb3d2dc7a2d3
# ╠═ed6ec0e9-5734-4569-b214-2c86c22d3c55
# ╠═cd9dea36-93f8-4461-85b1-87030d9367bb
# ╠═9c7e5bb4-8db0-4527-a8ec-331e2aed958b
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
# ╠═90b730b2-d64d-42c5-9aa6-31ee1c0927fe
# ╠═a644a33e-57a4-4af4-98be-1e0e84948069
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═127c1123-29cf-45a7-b7ab-7f57b4d6e58f
# ╠═2e9817a4-906d-4a0d-867c-d301cf3f2b1e
