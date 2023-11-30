### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 06484816-d7f4-4389-b0bf-e5a95214082e
begin
	using Pkg
	Pkg.activate()
	using Plots
	plotlyjs()

	using LilGuys
	using Optim
end

# ╔═╡ 5b956f10-b98d-465d-99fe-631d500228a3
using QuadGK

# ╔═╡ ddd5502b-5c4a-4608-93c8-511101dd5dfe
out = Output("/cosma/home/durham/dc-boye1/data/models/isolation_sculptor/out/")

# ╔═╡ 69b291b2-2589-489e-9a9b-19bdbc12114f
begin
snap_i = out[1, :]
snap_f = out[end, :]
end

# ╔═╡ 60817bd5-aff6-4981-8459-a392b978bd62
norm_p(x) = 4π * (1/(x + 1) + log(x + 1) - 1)

# ╔═╡ 373801b3-b4d4-4905-a14a-deca0ad0142f
norm_p(0)

# ╔═╡ 48db6a1d-0461-41af-b7c3-3134e507c97d
function ρ_nfw(r, params)
	x = r ./ params.r_s
	ρ = @. 1/params.r_s^3 * 1/norm_p(1e2/params.r_s) * 1/x* 1/ (1+x)^2
	return ρ
end

# ╔═╡ 49920be8-6057-434b-988f-37b398bb3b4c
params = LilGuys.NFWParams([0,0,0], 1)

# ╔═╡ 4ad3f249-9385-4e39-bb5c-5640846100b0
quadgk(x -> 4π*x^2 * ρ_nfw(x, params), 0, 1000)

# ╔═╡ 19f8aef2-9de7-47ed-a2cb-4350797598e3
function p_nfw(x::Matrix{F}, params)
	r = LilGuys.get_r(x .- params.x0)
	return ρ_nfw(r, params)
end

# ╔═╡ c3aec4c9-4be0-47bd-b1d7-5d4f228708d4
function nll(pos, params)
	if params.r_s < 0
		return Inf
	end
	return -sum(log.(p_nfw(pos, params)))
end

# ╔═╡ e8987520-8f46-405a-abff-5278fa91fde7
function f(x)
	x1, x2, x3, r = x
	params = LilGuys.NFWParams([x1, x2, x3], r)
	return nll(snap_i.pos, params)
end

# ╔═╡ 56d52bd5-b291-4df6-8902-985c67a1347e
begin
	result = optimize(f, [0.,0.,0.,10.], BFGS())
end

# ╔═╡ e4916530-c4c6-4731-8bfb-dc8fcdd9e156
result.minimizer

# ╔═╡ ee1e4fb7-2e64-4623-9fd6-4516a41dc6d9
nll(snap_i.pos, LilGuys.NFWParams([0.069, 0.46, 0.32], 6))

# ╔═╡ 1cdd3a9c-8d79-4452-9fc5-d0a7d445f2b4
out.softening^-1 * snap_i.m

# ╔═╡ 0283ad97-5f81-4daa-b522-e61cb6b86e39
r = LilGuys.get_r(snap_i.pos .- params.x0)

# ╔═╡ efaf3eda-6de3-4b36-bbea-e4002b6c3d6f
maximum(snap_i.Φ)

# ╔═╡ 41bec3db-f68e-480c-8eb2-8245a6d0f8ef
scatter(snap_i.pos[1, :], snap_i.pos[2, :])

# ╔═╡ d86be9ed-2010-4d56-a031-57acaba94921
scatter(snap_f.pos[1, :], snap_f.pos[2, :])

# ╔═╡ 943d7709-4082-4f99-b17e-701fcb82edab
snap_i.index

# ╔═╡ Cell order:
# ╠═06484816-d7f4-4389-b0bf-e5a95214082e
# ╠═ddd5502b-5c4a-4608-93c8-511101dd5dfe
# ╠═69b291b2-2589-489e-9a9b-19bdbc12114f
# ╠═373801b3-b4d4-4905-a14a-deca0ad0142f
# ╠═60817bd5-aff6-4981-8459-a392b978bd62
# ╠═48db6a1d-0461-41af-b7c3-3134e507c97d
# ╠═49920be8-6057-434b-988f-37b398bb3b4c
# ╠═4ad3f249-9385-4e39-bb5c-5640846100b0
# ╠═19f8aef2-9de7-47ed-a2cb-4350797598e3
# ╠═5b956f10-b98d-465d-99fe-631d500228a3
# ╠═c3aec4c9-4be0-47bd-b1d7-5d4f228708d4
# ╠═e8987520-8f46-405a-abff-5278fa91fde7
# ╠═56d52bd5-b291-4df6-8902-985c67a1347e
# ╠═e4916530-c4c6-4731-8bfb-dc8fcdd9e156
# ╠═ee1e4fb7-2e64-4623-9fd6-4516a41dc6d9
# ╠═1cdd3a9c-8d79-4452-9fc5-d0a7d445f2b4
# ╠═0283ad97-5f81-4daa-b522-e61cb6b86e39
# ╠═efaf3eda-6de3-4b36-bbea-e4002b6c3d6f
# ╠═41bec3db-f68e-480c-8eb2-8245a6d0f8ef
# ╠═d86be9ed-2010-4d56-a031-57acaba94921
# ╠═943d7709-4082-4f99-b17e-701fcb82edab
