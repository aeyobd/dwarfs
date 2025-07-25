### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 6a2488b6-61b6-11f0-200b-5bb12181877f
begin
	import Pkg; Pkg.activate()

	import Agama
	using Arya, CairoMakie

	using LilGuys
end

# ╔═╡ ce9cc08a-a809-4070-a453-fdd00a31046b
import TOML

# ╔═╡ 10400b5e-ae6a-4d12-aaa6-ff540a9ec966
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ 8172501e-5d18-48b2-88c6-5d1af17e4b4a
coords_i = LilGuys.transform.(LilGuys.Galactocentric, LilGuys.rand_coords(obs_props, 100))

# ╔═╡ b63d16aa-3759-4308-bf9c-c6add56485e0
pos_i = hcat(LilGuys.position.(coords_i)...)

# ╔═╡ da57300c-4752-4d2f-906e-da7a0f632b0b
vel_i = hcat(LilGuys.velocity.(coords_i)...) ./ V2KMS

# ╔═╡ 841f94b7-3c63-4713-92f7-2491c42f4edf
LilGuys.transform.(LilGuys.Galactocentric, coords_i)

# ╔═╡ 546fa6a8-766e-4007-9b89-1ae8fd508295
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 39369256-56b3-4131-a3df-857c9782abee
tmax = -10/T2GYR

# ╔═╡ 70ba6747-e0f8-4c72-8bf4-8c1bb1a3fa7f
o = Agama.orbit(pot, pos_i, vel_i, timerange=(0, tmax), N=1001)

# ╔═╡ 952e39d8-2fc4-42ce-aaad-a18547d75b3a
o2 = LilGuys.agama_orbit(pot, coords_i, timerange=(0, tmax), N=1001)

# ╔═╡ 17fcce0f-ddc5-436e-83fd-55f3c4ddb1fd


# ╔═╡ 47dedc74-a0c7-46d9-a0f9-8bf4ec0ea898
o3 = [LilGuys.resample(LilGuys.leapfrog(pot, coord), o2[1].times) for coord in coords_i[1:100]]

# ╔═╡ 689dc319-df43-4111-8256-06fb2bf4be16
extra_steps_2 = 10

# ╔═╡ dca594ff-4846-4b40-af6c-ae40ed1d6706
o4 = [LilGuys.resample(LilGuys.leapfrog(pot, coord, method=LilGuys.step_dkd), o2[1].times) for coord in coords_i[1:100]]

# ╔═╡ 3db40811-316f-45b6-8f52-108918400615
LilGuys.resample(o4[1], o2[1].times)

# ╔═╡ eee8ac0b-615b-4627-a229-66a77396f84d
10/T2GYR

# ╔═╡ c98c94de-5317-43c2-961a-a42364122ca6
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())
		
	
	for i in eachindex(o)
		lines!(o[i].positions[2, :], o[i].positions[3, :], color=COLORS[1], alpha=0.2)

	end

	fig
end

# ╔═╡ 0a64dfb3-7768-4f00-9ade-b6d42422c23f
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())
		
	
	for i in eachindex(o2)
		lines!(o2[i].positions[2, :], o2[i].positions[3, :], color=COLORS[1], alpha=0.2)

	end

	fig
end

# ╔═╡ 2c59fdae-5799-4ba6-b6f0-cebf57f40579
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())
		
	
	for i in eachindex(o3)
		lines!(o3[i].positions[2, :], o3[i].positions[3, :], color=COLORS[1], alpha=0.2)

	end

	fig
end

# ╔═╡ b63755fa-9a35-487d-9cce-47d6de9c3ee3
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())
		
	
	for i in eachindex(o4)
		lines!(o4[i].positions[2, :], o4[i].positions[3, :], color=COLORS[1], alpha=0.2)

	end

	fig
end

# ╔═╡ c9cce424-536a-4b21-ba64-87197b0628fe
o3[1].times[1:end] .- o2[1].times

# ╔═╡ 49899135-5a86-4e79-947e-a20d2bc2c219
o2[1].times[end], o3[1].times[end]

# ╔═╡ f8229fc0-c287-439c-8e1c-9242384a5fe5
let
	fig = Figure()
	ax = Axis(fig[1,1])
	for i in eachindex(o4)
		lines!(radii(o4[i].positions[:, 1:end] ), color=COLORS[2], alpha=0.2)

	end
	fig
end

# ╔═╡ 0e5ef66f-fe3f-46ef-b2dd-eeb441b925fa
let
	fig = Figure()
	ax = Axis(fig[1,1])
	for i in eachindex(o4)
		lines!(radii(o4[i].positions[:, 1:end] .- o2[i].positions[:, :]), color=COLORS[2], alpha=0.2)

	end
	
	
	for i in eachindex(o3)
		lines!(radii(o3[i].positions[:, 1:end] .- o2[i].positions[:, :]), color=COLORS[1], alpha=0.2)

	end

	

	fig
end

# ╔═╡ 1853f1bb-86a9-41fc-9dee-8986704ee94b
hist(LilGuys.minimum.(LilGuys.radii.(getfield.(o, :positions))))

# ╔═╡ b2d1aa12-f853-40ab-a35b-331fdbd5b848
LilGuys.position.(coords_i)

# ╔═╡ 12dba36a-ff6d-45d1-913c-9c2aed43065a
pot._py + pot._py

# ╔═╡ 110a2a04-7ddd-4ab6-bf98-525a881f5d61
LilGuys.Orbit(o[1])

# ╔═╡ Cell order:
# ╠═6a2488b6-61b6-11f0-200b-5bb12181877f
# ╠═ce9cc08a-a809-4070-a453-fdd00a31046b
# ╠═10400b5e-ae6a-4d12-aaa6-ff540a9ec966
# ╠═8172501e-5d18-48b2-88c6-5d1af17e4b4a
# ╠═b63d16aa-3759-4308-bf9c-c6add56485e0
# ╠═da57300c-4752-4d2f-906e-da7a0f632b0b
# ╠═841f94b7-3c63-4713-92f7-2491c42f4edf
# ╠═546fa6a8-766e-4007-9b89-1ae8fd508295
# ╠═39369256-56b3-4131-a3df-857c9782abee
# ╠═70ba6747-e0f8-4c72-8bf4-8c1bb1a3fa7f
# ╠═952e39d8-2fc4-42ce-aaad-a18547d75b3a
# ╠═17fcce0f-ddc5-436e-83fd-55f3c4ddb1fd
# ╠═47dedc74-a0c7-46d9-a0f9-8bf4ec0ea898
# ╠═689dc319-df43-4111-8256-06fb2bf4be16
# ╠═dca594ff-4846-4b40-af6c-ae40ed1d6706
# ╠═3db40811-316f-45b6-8f52-108918400615
# ╠═eee8ac0b-615b-4627-a229-66a77396f84d
# ╠═c98c94de-5317-43c2-961a-a42364122ca6
# ╠═0a64dfb3-7768-4f00-9ade-b6d42422c23f
# ╠═2c59fdae-5799-4ba6-b6f0-cebf57f40579
# ╠═b63755fa-9a35-487d-9cce-47d6de9c3ee3
# ╠═c9cce424-536a-4b21-ba64-87197b0628fe
# ╠═49899135-5a86-4e79-947e-a20d2bc2c219
# ╠═f8229fc0-c287-439c-8e1c-9242384a5fe5
# ╠═0e5ef66f-fe3f-46ef-b2dd-eeb441b925fa
# ╠═1853f1bb-86a9-41fc-9dee-8986704ee94b
# ╠═b2d1aa12-f853-40ab-a35b-331fdbd5b848
# ╠═12dba36a-ff6d-45d1-913c-9c2aed43065a
# ╠═110a2a04-7ddd-4ab6-bf98-525a881f5d61
