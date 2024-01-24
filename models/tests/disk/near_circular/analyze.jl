### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 71e3f40f-a912-4bcb-aa30-313f8b5dae9e
begin 
	import Pkg
	Pkg.activate()
	import LilGuys as lguys
	using Plots; plotly()

	pwd()
end

# ╔═╡ 7fde36fd-8f27-45c3-b255-989e394998c2
out = lguys.Output("out")

# ╔═╡ 8b30ba3c-8441-432c-a893-8e38db3599f5
begin 
	id = 15
	positions = lguys.extract(out, :positions, id)
	velocities = lguys.extract(out, :velocities, id)
	Φs = lguys.extract(out, :Φs_ext, id)
end

# ╔═╡ 74088189-e6f3-4ac7-803b-78de2559c0de
begin 
	plots = Ref{Vector{Plots.Plot}}()
	p = Ref{Plots.Plot}()
end

# ╔═╡ f6b3f948-354e-4405-94da-8ec644030ecf
begin 
	a = 3.944
	b = 0.311
	M = 6
	calc_Φ(R, z) = -lguys.G * M / sqrt(R^2 + (a + sqrt(z^2 + b^2))^2 )
end

# ╔═╡ 05b392b8-43eb-4f97-81d5-985513497908
begin 
	x = positions[1, :]
	y = positions[2, :]
	
	R = @. sqrt(x^2 + y^2)
	z = positions[3, :]
	r = @. sqrt(R^2 + z^2)
	
	vx = velocities[1, :]
	vy = velocities[2, :]
	vz = velocities[3, :]
	vR = @. (vx*x + vy*y) / R
	vϕ = @.  (vx*y - vy*x) / R

	v = @. sqrt(vx^2 + vy^2 + vz^2)
	Lz = vϕ .* R
end

# ╔═╡ c60b4267-f042-400e-b14b-77b1205f23b4
begin 
	plots[] = []
	
	p[] = plot(R, z)
	xlabel!("R")
	ylabel!("z")
	push!(plots[], p[])

	p[] = plot(out.times, R)
	xlabel!("t")
	ylabel!("R")
	push!(plots[], p[])

	p[] = plot(out.times, z)
	xlabel!("t")
	ylabel!("z")
	push!(plots[], p[])

	plots[]
end

# ╔═╡ d1177d02-9899-4fe8-b9b8-06f7330ee98f
begin 
	plots[] = []
	
	p[] = plot(vR, vz)
	xlabel!("vR")
	ylabel!("vz")
	push!(plots[], p[])

	p[] = plot(out.times, vR)
	xlabel!("t")
	ylabel!("vR")
	push!(plots[], p[])

	p[] = plot(out.times, vz)
	xlabel!("t")
	ylabel!("vz")
	push!(plots[], p[])

	plots[]
end

# ╔═╡ f9dd10c3-40cd-492c-9d62-0fa505ee065b
begin 
	plots[] = []
	
	p[] = plot(calc_Φ.(R, z) ./ Φs .- 1)
	ylabel!("relerr in potential")
	push!(plots[], p[])

	p[] = plot(Lz)
	ylabel!("Lz")
	push!(plots[], p[])

	p[] = plot(out.times, lguys.calc_E_spec.(Φs, v))
	ylabel!("total energy")
	push!(plots[], p[])

end

# ╔═╡ Cell order:
# ╠═71e3f40f-a912-4bcb-aa30-313f8b5dae9e
# ╠═7fde36fd-8f27-45c3-b255-989e394998c2
# ╠═8b30ba3c-8441-432c-a893-8e38db3599f5
# ╠═74088189-e6f3-4ac7-803b-78de2559c0de
# ╠═f6b3f948-354e-4405-94da-8ec644030ecf
# ╠═05b392b8-43eb-4f97-81d5-985513497908
# ╠═c60b4267-f042-400e-b14b-77b1205f23b4
# ╠═d1177d02-9899-4fe8-b9b8-06f7330ee98f
# ╠═f9dd10c3-40cd-492c-9d62-0fa505ee065b
