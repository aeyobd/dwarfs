### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 71e3f40f-a912-4bcb-aa30-313f8b5dae9e
begin 
	import Pkg
	Pkg.activate()
	import LilGuys as lguys
	using CairoMakie
	using Arya
end

# ╔═╡ 7fde36fd-8f27-45c3-b255-989e394998c2
out = lguys.Output("out/combined.hdf5")

# ╔═╡ 8b30ba3c-8441-432c-a893-8e38db3599f5
begin 
	id = 2
	positions = lguys.extract_vector(out, :positions, id)
	velocities = lguys.extract_vector(out, :velocities, id)
	# Φs = lguys.extract(out, :Φs_ext, id)
end

# ╔═╡ d805a2b3-b180-4f59-930d-20ff079a3bcc
out[1].Φs

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

# ╔═╡ 250a2dc5-f40e-4a80-99df-429fba95df64
plot(x, y)

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

# ╔═╡ e798f63f-c5ad-46c7-a07e-478f95635367


# ╔═╡ f32c522f-a979-40c3-a64f-97a7053b95b9
begin
	plot(out.times, 1/2*v.^2, label="T")
	plot!(out.times, Φs, label="V")
	plot!(out.times, Φs + 1/2*v.^2, label="total")
	xlabel!("time")
	ylabel!("energy")
end

# ╔═╡ Cell order:
# ╠═71e3f40f-a912-4bcb-aa30-313f8b5dae9e
# ╠═7fde36fd-8f27-45c3-b255-989e394998c2
# ╠═8b30ba3c-8441-432c-a893-8e38db3599f5
# ╠═d805a2b3-b180-4f59-930d-20ff079a3bcc
# ╠═f6b3f948-354e-4405-94da-8ec644030ecf
# ╠═05b392b8-43eb-4f97-81d5-985513497908
# ╠═250a2dc5-f40e-4a80-99df-429fba95df64
# ╠═c60b4267-f042-400e-b14b-77b1205f23b4
# ╠═d1177d02-9899-4fe8-b9b8-06f7330ee98f
# ╠═f9dd10c3-40cd-492c-9d62-0fa505ee065b
# ╠═e798f63f-c5ad-46c7-a07e-478f95635367
# ╟─f32c522f-a979-40c3-a64f-97a7053b95b9
