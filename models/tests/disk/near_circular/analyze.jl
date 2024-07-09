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
	id = 1
	positions = lguys.extract_vector(out, :positions, id)
	velocities = lguys.extract_vector(out, :velocities, id)
	#Φs = lguys.extract(out, :Φs_ext, id)
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
	plot(R, z)
end

# ╔═╡ 8e596b73-7555-4eba-9065-d300f7f8228d
plot(out.times, R)

# ╔═╡ f4190584-2f7d-4ef8-87bb-b828af4d0a6e
plot(out.times, z)

# ╔═╡ 0fee7f88-6acb-4a21-b326-9a4f98247f6f
plot(x, y)

# ╔═╡ d1177d02-9899-4fe8-b9b8-06f7330ee98f
let 
	@info plot(vR, vz)


	@info plot(out.times, vR)


	@info plot(out.times, vz)

end

# ╔═╡ f9dd10c3-40cd-492c-9d62-0fa505ee065b
begin 
	plot(calc_Φ.(R, z) ./ Φs .- 1)
end

# ╔═╡ 0e511698-9702-408f-8013-f2aab5768fd1
plot(out.times, lguys.calc_E_spec.(Φs, v))

# ╔═╡ 04e83c2a-f807-4cec-940c-7beceb5fd2e6
 plot(Lz)

# ╔═╡ Cell order:
# ╠═71e3f40f-a912-4bcb-aa30-313f8b5dae9e
# ╠═7fde36fd-8f27-45c3-b255-989e394998c2
# ╠═8b30ba3c-8441-432c-a893-8e38db3599f5
# ╠═f6b3f948-354e-4405-94da-8ec644030ecf
# ╠═05b392b8-43eb-4f97-81d5-985513497908
# ╠═c60b4267-f042-400e-b14b-77b1205f23b4
# ╠═8e596b73-7555-4eba-9065-d300f7f8228d
# ╠═f4190584-2f7d-4ef8-87bb-b828af4d0a6e
# ╠═0fee7f88-6acb-4a21-b326-9a4f98247f6f
# ╠═d1177d02-9899-4fe8-b9b8-06f7330ee98f
# ╠═f9dd10c3-40cd-492c-9d62-0fa505ee065b
# ╠═0e511698-9702-408f-8013-f2aab5768fd1
# ╠═04e83c2a-f807-4cec-940c-7beceb5fd2e6
