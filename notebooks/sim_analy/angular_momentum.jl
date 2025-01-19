### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 64492622-d5d9-11ef-0503-8de99f2e38e1
begin
	using Pkg;Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 9a8824ed-e23f-42e1-bac8-f81c1dbcd79e
md"""
This notebook was a quick exploration to see if the angular momentum decay of a tidally stripped galaxy happened to be easily parameterizable. However, it seems that this decay does not directyl correspond with even radius but is also not well measured, so very challenging to determine.
"""

# ╔═╡ f93754f6-7612-4cca-bf77-c1288dcd58f4
modelname = "ursa_minor/1e6_v37_r5.0/orbit_mean"

# ╔═╡ ab029f52-8d09-486c-8c8d-af0ce4eb13d0
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", modelname) * "/"

# ╔═╡ c1ba4866-5432-4879-93d8-1c1108d01521
out = Output(modeldir)

# ╔═╡ 21d76586-ed8f-46cc-80fe-68081f9af5aa
x_cen = out.x_cen

# ╔═╡ 868cce89-3ce6-4f9d-af06-ff18e1338052
v_cen = out.v_cen

# ╔═╡ 5a1ee037-4980-461b-abc7-4a4bc34a7050
times = out.times

# ╔═╡ 3b25c4bc-307e-44e0-ad85-bba0f5764d3f
profs = LilGuys.read_structs_from_hdf5(modeldir * "profiles.hdf5", LilGuys.MassProfile3D)

# ╔═╡ a83ebd19-5453-4240-a0e8-ac68eccba04f
profs_idx = parse.(Int64, first.(profs))

# ╔═╡ 2dec79c4-e303-4b4d-bf0a-23d05fb0b640
masses = [prof.second.N_bound for prof in profs][sortperm(profs_idx)]

# ╔═╡ 884319be-a394-47dd-88f4-6a202e7af42f
L = LilGuys.calc_L_spec(x_cen, v_cen)

# ╔═╡ fc2019d5-641d-44b7-86dc-36f50a87c42a
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "L"
	)

	lines!(times * T2GYR, calc_r(L))

	fig
end

# ╔═╡ 306020b7-60bd-4bc1-9647-4642b55229b3
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "L"
	)

	lines!(times * T2GYR, calc_r(L) .* masses)

	fig
end

# ╔═╡ 1b46d5af-bc43-49c0-a47b-9ede83316a63
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "r / kpc",
		ylabel = "L"
	)

	lines!(calc_r(x_cen), calc_r(L))

	fig
end

# ╔═╡ 94307e69-68f5-47be-a365-ca57a8f48530
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "r / kpc",
		ylabel = L"\dot{L}"
	)

	scatter!(calc_r(x_cen), LilGuys.gradient(calc_r(L)))

	fig
end

# ╔═╡ c47785ca-a344-4dcf-bde4-55b86994244c
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "r / kpc",
		ylabel = L"\dot{L}"
	)

	scatter!(calc_r(x_cen), LilGuys.gradient(L[1, :]), color=times)

	fig
end

# ╔═╡ Cell order:
# ╟─9a8824ed-e23f-42e1-bac8-f81c1dbcd79e
# ╠═64492622-d5d9-11ef-0503-8de99f2e38e1
# ╠═f93754f6-7612-4cca-bf77-c1288dcd58f4
# ╠═ab029f52-8d09-486c-8c8d-af0ce4eb13d0
# ╠═c1ba4866-5432-4879-93d8-1c1108d01521
# ╠═21d76586-ed8f-46cc-80fe-68081f9af5aa
# ╠═868cce89-3ce6-4f9d-af06-ff18e1338052
# ╠═5a1ee037-4980-461b-abc7-4a4bc34a7050
# ╠═3b25c4bc-307e-44e0-ad85-bba0f5764d3f
# ╠═a83ebd19-5453-4240-a0e8-ac68eccba04f
# ╠═2dec79c4-e303-4b4d-bf0a-23d05fb0b640
# ╠═884319be-a394-47dd-88f4-6a202e7af42f
# ╠═fc2019d5-641d-44b7-86dc-36f50a87c42a
# ╠═306020b7-60bd-4bc1-9647-4642b55229b3
# ╠═1b46d5af-bc43-49c0-a47b-9ede83316a63
# ╠═94307e69-68f5-47be-a365-ca57a8f48530
# ╠═c47785ca-a344-4dcf-bde4-55b86994244c
