### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 64492622-d5d9-11ef-0503-8de99f2e38e1
begin
	using Pkg;Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 0b72ae9c-7713-4733-bc65-75bcf45f5a1b
begin
	using PythonCall
	agama = pyimport("agama")
	np = pyimport("numpy")
end

# ╔═╡ 16b1c4f3-e5b5-4445-a3c3-3cdb5ef95fbc
using Measurements

# ╔═╡ 9a8824ed-e23f-42e1-bac8-f81c1dbcd79e
md"""
This notebook was a quick exploration to see if the angular momentum decay of a tidally stripped galaxy happened to be easily parameterizable. However, it seems that this decay does not directyl correspond with even radius but is also not well measured, so very challenging to determine.
"""

# ╔═╡ f93754f6-7612-4cca-bf77-c1288dcd58f4
modelname = "ursa_minor/1e6_v37_r5.0/orbit_mean"

# ╔═╡ ab029f52-8d09-486c-8c8d-af0ce4eb13d0
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", modelname) * "/"

# ╔═╡ e1e4a9da-2e81-4f63-ba22-2318312261c4
readdir(modeldir * "simulation")

# ╔═╡ 48c95a1d-4dd0-434b-a567-aee23a12150a
pot = agama.Potential(modeldir * "simulation/agama_potential.ini")

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

# ╔═╡ 57f349aa-de32-4e3d-bbc1-bccb6f30e9ae
K = 1/2 * calc_r(v_cen) .^ 2

# ╔═╡ 4bf63199-5bf5-48fa-a057-d983c8733827
af = agama.ActionFinder(pot)

# ╔═╡ 9fe8c08b-0043-4220-b222-26f5433d2863
act_py, ang_py, freq_py = af(np.array(vcat(x_cen, v_cen)'), angles=true, frequencies=true)

# ╔═╡ bfbaba5b-8fd4-423f-82b2-6f5a0bd632cb
py2mat(x) = pyconvert(Matrix{Float64}, x)'

# ╔═╡ e747b374-bf35-40c3-a6c6-c75a08f2264c
py2vec(x) = pyconvert(Vector{Float64}, x)

# ╔═╡ 08b1cfeb-18e2-4681-a3fc-ac24e7d6cc95
actions = py2mat(act_py)

# ╔═╡ 3f889717-adf5-47f2-a9be-96cee02aa188
py2mat(ang_py)

# ╔═╡ 0676ee62-5c6c-405f-8cbd-07702a5cfb49
Jr, Jz, Jϕ = eachrow(actions)

# ╔═╡ 528b9bd0-9fe8-4236-8691-f284a12cfca3
freq_py

# ╔═╡ 0e572199-cc99-44bc-8829-d0415a22d73a
Φ = pyconvert(Vector{Float64}, pot.potential(np.array(x_cen')))

# ╔═╡ 4d7db905-4f80-447a-86ea-0027f32e8ab1
E = Φ .+ K

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

# ╔═╡ e8c953ba-9a7d-4896-84bf-66559f214b74
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "K"
	)

	lines!(times * T2GYR, K)
	lines!(times * T2GYR, Φ)
	lines!(times * T2GYR, E)

	fig
end

# ╔═╡ 882fa476-2c4b-4f40-852d-7c8d55c3fb0a
R =  @. sqrt(x_cen[1, :]^2 + x_cen[2, :]^2)

# ╔═╡ fb9981a0-ef43-4b08-9f1a-e21bab92b3d0
z = x_cen[3, :]

# ╔═╡ e3431072-6206-40e0-9212-aeec0a2b007e
r = calc_r(x_cen)

# ╔═╡ 87a86e3e-ad3b-485e-8099-b2d2f19857c7
θ = atan.(R, z)

# ╔═╡ 8bf93541-5154-4a04-94ea-21df0328848c
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "K"
	)

	lines!(times * T2GYR, r)

	lines!(times * T2GYR, R)
	lines!(times * T2GYR, z)

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

# ╔═╡ 9f01bf9b-64a1-413c-89c8-e92fdfffb486
md"""
## Change in energy quantities
"""

# ╔═╡ 4a58cbc5-0b13-46eb-82a1-63289ee0d94e
obs_props_file = joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/observed_properties.toml")

# ╔═╡ cdb1102b-99f2-475a-b291-adc7f715f395
icrs = LilGuys.coord_from_file(obs_props_file)

# ╔═╡ d2e85adb-a278-4c9d-8c59-c6dd1f061aca
icrs_err = LilGuys.coord_err_from_file(obs_props_file)

# ╔═╡ f70ec6c6-2296-4746-99cd-e4d36ac87dc5
r_h = (12/60)

# ╔═╡ e4c21527-47ef-47d1-b011-baf50f5aeaeb
icrs_u = LilGuys.ICRS(
	ra = icrs.ra ± r_h,
	dec = icrs.dec ± r_h,
	pmra = icrs.pmra ± icrs_err.pmra,
	pmdec = icrs.pmdec ± icrs_err.pmdec,
	distance = icrs.distance ± icrs_err.distance,
	radial_velocity = icrs.radial_velocity ± icrs_err.radial_velocity
)

# ╔═╡ 802a338b-9465-4c11-98e9-d055047aadb9
gc = LilGuys.transform(Galactocentric, icrs_u)

# ╔═╡ 3400e6de-6b05-4376-86f6-744d8bc8c202
act_obs_py, ang_obs, _ = af(np.array(Measurements.value.([gc.x, gc.y, gc.z, gc.v_x/V2KMS, gc.v_y/V2KMS, gc.v_z/V2KMS])), angles=true, frequencies=true)

# ╔═╡ 8d14d9f0-ecf9-4130-94d0-d1df0fd83c1a
act_obs = py2vec(act_obs_py)

# ╔═╡ dcb399a2-0416-4b76-8322-f28b0f795406
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "L"
	)

	lines!(times * T2GYR, Jr, label = "Jr")
	lines!(times * T2GYR, Jz, label="Jz")
	lines!(times * T2GYR, Jϕ, label="Jϕ")

	lines!(times * T2GYR, L[3, :])

	hlines!(act_obs, linestyle=:dot)

	axislegend()
	fig
end

# ╔═╡ 4166ded2-826d-4c2a-8d01-098e9c4b9615
⊕(x::Real, y::Real) = sqrt(x^2 + y^2)

# ╔═╡ bd1adc5d-2564-4c96-84dc-341e39a0aed2
4 ⊕ 3

# ╔═╡ 117d309b-fb6c-4012-95cf-d2aa513e7c4e
calc_χ2(val, exp) = (val .- Measurements.value(exp)) .^2 / Measurements.uncertainty(exp)^2

# ╔═╡ 34308090-58be-49cb-9c43-7597d16a9351
vR = @. (v_cen[1, :] * x_cen[1, :] + v_cen[2, :] * x_cen[2, :]) ./ R


# ╔═╡ 391c0fb0-b184-4f00-9401-c5f0ecb2c0c5
vR_obs = (gc.v_x * gc.x + gc.v_y * gc.y ) / (gc.x ⊕ gc.y)

# ╔═╡ 8bc38fcc-a744-4fa4-972d-556b06be6434
vϕ_obs = (gc.v_x * gc.y - gc.v_y * gc.x ) / (gc.x ⊕ gc.y)

# ╔═╡ b6cef668-d06c-4f43-8052-a458f323f70a
R_obs = gc.x ⊕ gc.y

# ╔═╡ cde0b763-63b5-49c6-8e07-4d6eba1885b0
R_obs * vϕ_obs / V2KMS

# ╔═╡ bf921c3f-5f84-4c01-b400-80be780da089
LilGuys.calc_L_spec(LilGuys.position_of(gc), LilGuys.velocity_of(gc)) / V2KMS

# ╔═╡ 2d3a9b98-3faf-4493-9d54-b32c38e300cb
vϕ = @. (v_cen[1, :] * x_cen[2, :] - v_cen[2, :] * x_cen[1, :]) ./ R

# ╔═╡ 58ca6746-037e-432f-bc61-3dcb710bbbd9
begin 
	χ2 = calc_χ2(R, gc.x ⊕ gc.y)
	χ2 .+= calc_χ2(z, gc.z)
	χ2 .+= calc_χ2(vR * V2KMS, vR_obs)
	χ2 .+= calc_χ2(v_cen[3, :] * V2KMS, gc.v_z)
	χ2 .+= calc_χ2(vϕ * V2KMS, vϕ_obs)

end

# ╔═╡ d0a39048-7b8c-4077-b1ac-942d12d4c33b
lines(χ2)

# ╔═╡ 2347d72c-59a9-417a-9c07-a25c6703331c
χ2_filt = 100:212

# ╔═╡ e97afd43-394e-45c0-8d18-9e33629b9a82
idx_f = eachindex(times)[χ2_filt][argmin(χ2[χ2_filt])]

# ╔═╡ cb6cebfa-49aa-4e7a-967c-ee1cbdf0e211
let
	fig = Figure()
	ax = Axis(fig[1,1])


	idx = 140:idx_f
	lines!(R[idx], z[idx], color=times[idx])

	R1 = gc.x .⊕ gc.y
	z1 = gc.z
	errscatter!(Measurements.value(R1), Measurements.value(z1))
	fig
end

# ╔═╡ db9cd7db-895d-4b7e-a159-6847b6f32040
let
	fig = Figure()
	ax = Axis(fig[1,1])


	idx = 140:idx_f
	vz = v_cen[3, :]
	lines!(vR[idx] * V2KMS, vz[idx] * V2KMS, color=times[idx])

	errscatter!(Measurements.value(vR_obs), Measurements.value(gc.v_z))
	fig
end

# ╔═╡ Cell order:
# ╟─9a8824ed-e23f-42e1-bac8-f81c1dbcd79e
# ╠═64492622-d5d9-11ef-0503-8de99f2e38e1
# ╠═f93754f6-7612-4cca-bf77-c1288dcd58f4
# ╠═ab029f52-8d09-486c-8c8d-af0ce4eb13d0
# ╠═0b72ae9c-7713-4733-bc65-75bcf45f5a1b
# ╠═e1e4a9da-2e81-4f63-ba22-2318312261c4
# ╠═48c95a1d-4dd0-434b-a567-aee23a12150a
# ╠═c1ba4866-5432-4879-93d8-1c1108d01521
# ╠═21d76586-ed8f-46cc-80fe-68081f9af5aa
# ╠═868cce89-3ce6-4f9d-af06-ff18e1338052
# ╠═5a1ee037-4980-461b-abc7-4a4bc34a7050
# ╠═3b25c4bc-307e-44e0-ad85-bba0f5764d3f
# ╠═a83ebd19-5453-4240-a0e8-ac68eccba04f
# ╠═2dec79c4-e303-4b4d-bf0a-23d05fb0b640
# ╠═884319be-a394-47dd-88f4-6a202e7af42f
# ╠═57f349aa-de32-4e3d-bbc1-bccb6f30e9ae
# ╠═4bf63199-5bf5-48fa-a057-d983c8733827
# ╠═9fe8c08b-0043-4220-b222-26f5433d2863
# ╠═3400e6de-6b05-4376-86f6-744d8bc8c202
# ╠═8d14d9f0-ecf9-4130-94d0-d1df0fd83c1a
# ╠═bfbaba5b-8fd4-423f-82b2-6f5a0bd632cb
# ╠═e747b374-bf35-40c3-a6c6-c75a08f2264c
# ╠═08b1cfeb-18e2-4681-a3fc-ac24e7d6cc95
# ╠═3f889717-adf5-47f2-a9be-96cee02aa188
# ╠═0676ee62-5c6c-405f-8cbd-07702a5cfb49
# ╠═528b9bd0-9fe8-4236-8691-f284a12cfca3
# ╠═0e572199-cc99-44bc-8829-d0415a22d73a
# ╠═4d7db905-4f80-447a-86ea-0027f32e8ab1
# ╠═fc2019d5-641d-44b7-86dc-36f50a87c42a
# ╠═dcb399a2-0416-4b76-8322-f28b0f795406
# ╠═e8c953ba-9a7d-4896-84bf-66559f214b74
# ╠═882fa476-2c4b-4f40-852d-7c8d55c3fb0a
# ╠═fb9981a0-ef43-4b08-9f1a-e21bab92b3d0
# ╠═e3431072-6206-40e0-9212-aeec0a2b007e
# ╠═87a86e3e-ad3b-485e-8099-b2d2f19857c7
# ╠═8bf93541-5154-4a04-94ea-21df0328848c
# ╠═306020b7-60bd-4bc1-9647-4642b55229b3
# ╠═1b46d5af-bc43-49c0-a47b-9ede83316a63
# ╠═94307e69-68f5-47be-a365-ca57a8f48530
# ╠═c47785ca-a344-4dcf-bde4-55b86994244c
# ╠═9f01bf9b-64a1-413c-89c8-e92fdfffb486
# ╠═4a58cbc5-0b13-46eb-82a1-63289ee0d94e
# ╠═cdb1102b-99f2-475a-b291-adc7f715f395
# ╠═d2e85adb-a278-4c9d-8c59-c6dd1f061aca
# ╠═16b1c4f3-e5b5-4445-a3c3-3cdb5ef95fbc
# ╠═f70ec6c6-2296-4746-99cd-e4d36ac87dc5
# ╠═e4c21527-47ef-47d1-b011-baf50f5aeaeb
# ╠═802a338b-9465-4c11-98e9-d055047aadb9
# ╠═4166ded2-826d-4c2a-8d01-098e9c4b9615
# ╠═bd1adc5d-2564-4c96-84dc-341e39a0aed2
# ╠═cb6cebfa-49aa-4e7a-967c-ee1cbdf0e211
# ╠═db9cd7db-895d-4b7e-a159-6847b6f32040
# ╠═117d309b-fb6c-4012-95cf-d2aa513e7c4e
# ╠═34308090-58be-49cb-9c43-7597d16a9351
# ╠═391c0fb0-b184-4f00-9401-c5f0ecb2c0c5
# ╠═8bc38fcc-a744-4fa4-972d-556b06be6434
# ╠═b6cef668-d06c-4f43-8052-a458f323f70a
# ╠═cde0b763-63b5-49c6-8e07-4d6eba1885b0
# ╠═bf921c3f-5f84-4c01-b400-80be780da089
# ╠═2d3a9b98-3faf-4493-9d54-b32c38e300cb
# ╠═58ca6746-037e-432f-bc61-3dcb710bbbd9
# ╠═d0a39048-7b8c-4077-b1ac-942d12d4c33b
# ╠═2347d72c-59a9-417a-9c07-a25c6703331c
# ╠═e97afd43-394e-45c0-8d18-9e33629b9a82
