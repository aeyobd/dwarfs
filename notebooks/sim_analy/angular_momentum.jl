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

# ╔═╡ 86278c5b-dac2-4a93-a189-9720c362b9f3
using CSV, DataFrames

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

# ╔═╡ 26a86ae1-c042-4be0-90e6-4df515f5b2a4


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

# ╔═╡ 244bea42-fbb8-41ec-97f5-174b7661b713
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "r / kpc",
	)

	lines!(times * T2GYR, calc_r(x_cen))

	fig
end

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

# ╔═╡ ddd4c109-a467-46b1-b3dc-99af8d4aea96
function plot_meas!(meas)
	hlines!(Measurements.value(meas), color=:black, linestyle=:dot)

	y = Measurements.value(meas)
	ye = Measurements.uncertainty(meas)

	hspan!(y-ye, y+ye, color=(:black, 0.2))

end

# ╔═╡ df85719e-b62d-4d00-b592-78afca6dd7cb
t_scale = 0.94

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

# ╔═╡ 211bd090-9188-4fcb-9cd3-99f49a325a50
lines(times * T2GYR, E)


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

# ╔═╡ 569f3a8a-7134-4e74-a757-8fbe35f50cfd
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "mass"
	)

	lines!(times * T2GYR, log10.(masses))

	fig
end

# ╔═╡ 0a8949c6-3706-4fba-9255-015b72c01546
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!(log10.(masses), calc_r(L))

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

	lines!(times * T2GYR, -0.87L[1, :])
	lines!(times * T2GYR, -0.5L[2, :])
	lines!(times * T2GYR, L[3, :])

	hlines!(act_obs, linestyle=:dot)

	axislegend()
	fig
end

# ╔═╡ 66988f4b-aa3d-4106-b9dd-17cde5c9e529
L_obs = LilGuys.calc_L_spec(LilGuys.position_of(gc), LilGuys.velocity_of(gc) / V2KMS)

# ╔═╡ 615988d2-5ccd-4ac4-b72d-352eeb122125
let
	fig = Figure()

	for i in 1:3
		ax = Axis(fig[i,1],
			xlabel = "time / Gyr",
			ylabel = ["x", "y", "z"][i]
		)
	
		lines!(times * T2GYR, L[i, :])
		plot_meas!(L_obs[i])

	end

	fig
end

# ╔═╡ 4166ded2-826d-4c2a-8d01-098e9c4b9615
⊕(x::Real, y::Real) = sqrt(x^2 + y^2)

# ╔═╡ 4839120c-63ed-4200-833a-bd171f3c9237
begin 
	orbit_exp = CSV.read(modeldir * "simulation/orbit.csv", DataFrame)
	orbit_exp.t .-= orbit_exp.t[1]
	orbit_exp[!, :R] = orbit_exp.x .⊕ orbit_exp.y
	orbit_exp[!, :Vr] = @. 1/orbit_exp.R *( orbit_exp.v_x * orbit_exp.x + orbit_exp.v_y*orbit_exp.y)
	orbit_exp[!, :Vphi] = @. 1/orbit_exp.R *( orbit_exp.v_x * orbit_exp.y - orbit_exp.v_y*orbit_exp.x)
	
end

# ╔═╡ 94309696-9545-499e-990c-759c4488e1d2
r_f = gc.x ⊕ gc.y ⊕ gc.z

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

# ╔═╡ 662af96b-46a5-425a-85c0-42ddc038e304
let
	fig = Figure(
		size=(500, 700)
	)
	
	ax1 = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "R/kpc",
	)
	#vlines!(times[idx_f] * T2GYR)

	lines!(times * T2GYR, R)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.x .⊕ orbit_exp.y)
	plot_meas!(gc.x ⊕ gc.y)
	hidexdecorations!(ax1)
	ylims!(10, 60)


	ax2 = Axis(fig[2,1],
		xlabel = "time / Gyr",
		ylabel = "z/kpc"
	)

	lines!(times * T2GYR, z)
	plot_meas!(gc.z)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.z)
	hidexdecorations!(ax2)
	ylims!(30, 90)


	ax3 = Axis(fig[3,1],
		xlabel = "time / Gyr",
		ylabel = L"$v_R$ / km\,s$^{-1}$"
	)

	lines!(times * T2GYR, vR * V2KMS)
	plot_meas!(vR_obs)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.Vr * V2KMS)

	hidexdecorations!(ax3)
	ylims!(-150, 150)
	
	ax4 = Axis(fig[4,1],
		xlabel = "time / Gyr",
		ylabel = L"$v_\phi$ / km\,s$^{-1}$"
	)

	lines!(times * T2GYR, vϕ * V2KMS)
	plot_meas!(vϕ_obs)
	hidexdecorations!(ax4)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.Vphi * V2KMS)
	ylims!(20, 150)
	
	ax5 = Axis(fig[5,1],
		xlabel = "time / Gyr",
		ylabel = L"$v_z$ / km\,s$^{-1}$"
	)

	lines!(times * T2GYR, v_cen[3, :] * V2KMS)
	plot_meas!(gc.v_z)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.v_z * V2KMS)

	ylims!(-200, 100)
	linkxaxes!(ax1, ax2, ax3, ax4, ax5)
	xlims!(8.5, 9.5)
	fig
end

# ╔═╡ 58ca6746-037e-432f-bc61-3dcb710bbbd9
begin 
	χ2 = calc_χ2(R, gc.x ⊕ gc.y)
	χ2 .+= calc_χ2(z, gc.z)
	χ2 .+= calc_χ2(vR * V2KMS, vR_obs)
	χ2 .+= calc_χ2(v_cen[3, :] * V2KMS, gc.v_z)
	χ2 .+= calc_χ2(vϕ * V2KMS, vϕ_obs)

end

# ╔═╡ d0a39048-7b8c-4077-b1ac-942d12d4c33b
lines(log10.(χ2))

# ╔═╡ 2347d72c-59a9-417a-9c07-a25c6703331c
χ2_filt = 100:212

# ╔═╡ e97afd43-394e-45c0-8d18-9e33629b9a82
idx_f = eachindex(times)[χ2_filt][argmin(χ2[χ2_filt])]

# ╔═╡ eae7a157-a6fb-466a-a293-785139df0ef3
times[idx_f] * T2GYR

# ╔═╡ 4e8f8a44-1d3c-4666-8a64-6e708eb3f165
vR[idx_f] * V2KMS

# ╔═╡ 5f35b33a-855d-49f1-bb65-a0e5f60d61cb
let
	fig = Figure(
		size=(500, 700)
	)
	
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "R/kpc",
	)
	vlines!(times[idx_f] * T2GYR)

	lines!(times * T2GYR, R)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.x .⊕ orbit_exp.y)
	plot_meas!(gc.x ⊕ gc.y)
	hidexdecorations!(ax)


	ax = Axis(fig[2,1],
		xlabel = "time / Gyr",
		ylabel = "z/kpc"
	)

	lines!(times * T2GYR, z)
	plot_meas!(gc.z)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.z)
	hidexdecorations!(ax)


	ax = Axis(fig[3,1],
		xlabel = "time / Gyr",
		ylabel = "vR"
	)

	lines!(times * T2GYR, vR * V2KMS)
	plot_meas!(vR_obs)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.Vr * V2KMS)

	hidexdecorations!(ax)
	
	ax = Axis(fig[4,1],
		xlabel = "time / Gyr",
		ylabel = "vphi"
	)

	lines!(times * T2GYR, vϕ * V2KMS)
	plot_meas!(vϕ_obs)
	hidexdecorations!(ax)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.Vphi * V2KMS)

	
	ax = Axis(fig[5,1],
		xlabel = "time / Gyr",
		ylabel = "vz"
	)

	lines!(times * T2GYR, v_cen[3, :] * V2KMS)
	plot_meas!(gc.v_z)
	lines!(orbit_exp.t * T2GYR * t_scale, orbit_exp.v_z * V2KMS)


	fig
end

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

# ╔═╡ 658280af-5bf9-47b1-aee1-21550580c6a6
v_cen[3, idx_f] * V2KMS, gc.v_z

# ╔═╡ b80957dc-b24b-4e06-8aad-07018312b370
vϕ[idx_f] * V2KMS, vϕ_obs

# ╔═╡ 917d02ee-30f2-4338-8d8a-91f9fd5a5044
vR[idx_f] * V2KMS, vR_obs

# ╔═╡ 6e59f120-6f5c-4ca1-967d-1034ed4f3473
R[idx_f], gc.x ⊕ gc.y

# ╔═╡ 433308c7-ce37-4ace-9008-e1a456774d7e
z[idx_f], gc.z

# ╔═╡ a7265965-9092-451b-bf0b-8f32f8b28858
md"""
# Adjusting initial conditions
"""

# ╔═╡ f7a5d076-f97d-41bd-8ad6-a08b6d5ad1a8
import LinearAlgebra: ⋅, ×

# ╔═╡ 4dffdfaa-a1be-452d-9c8a-9479f0b7b767
pos_i = [orbit_exp.x[1], orbit_exp.y[1], orbit_exp.z[1]]

# ╔═╡ 05dc4779-28de-4a2e-a47f-231f3e909d9d
vel_i = [orbit_exp.v_x[1], orbit_exp.v_y[1], orbit_exp.v_z[1]]

# ╔═╡ aa48ab4a-d512-4568-9b46-394a7ef092c5
Φ_in(x) = pyconvert(Float64, pot.potential(x))

# ╔═╡ 7629ac1a-fc7c-4225-ac2b-752b569e04f1
function shift_initial_conditions(Φ, pos_i, vel_i; dE=0, dlogr=0, dLx=0, dLy=0, dLz=0)
	r_i = LilGuys.calc_r(pos_i)
	r_hat = pos_i / r_i
	
	L_i = LilGuys.calc_L_spec(pos_i, vel_i)

	E_i = Φ(pos_i) + 1/2 * calc_r(vel_i)^2


	# spherical approximation...;/
	r_apo_i = LilGuys.find_zero(r->Φ(r_hat * r) - E_i, r_i)


	E_new = E_i + dE
	r_apo_new = LilGuys.find_zero(r->Φ(r_hat * r) - E_new, r_i)
	
	r_new = r_i * r_apo_new / r_apo_i * 10^dlogr
	pos_new = r_hat * r_new
	L_new = [L_i[1] + dLx, L_i[2] + dLy, L_i[3] + dLz]
	
	v_new = sqrt(2E_new - 2*Φ(pos_new))		

	vel_new_sign = sign(pos_i ⋅ vel_i)
	vel_new = v_new * (L_new × r_hat / r_new/v_new + vel_new_sign*r_hat * sqrt(1-calc_r(L_new/r_new/v_new)^2))

	println(L_new)
	println(pos_new × vel_new)
	println(E_new, "\t", Φ(pos_new) + 1/2*calc_r(vel_new)^2)
	println("r", r_new, "\t", calc_r(pos_new), "\t", r_i)

	return pos_new, vel_new
end

# ╔═╡ a086ce7a-81b3-4bef-81a3-3e94a436b2e7
pos_i, vel_i

# ╔═╡ ba62386c-f7a2-4e15-a4bc-b6f57f06aa03
shift_initial_conditions(Φ_in, pos_i, vel_i, dE=0.000, dLx=0)

# ╔═╡ 79068d79-b896-45ca-b34f-9e2f6b835a1b
shift_initial_conditions(Φ_in, pos_i, vel_i, dE=0.06, dLx=-3)

# ╔═╡ 5e968467-0a07-4c59-8779-6a8977d707e7


# ╔═╡ 2da8ed33-bc9c-4ff1-b1c8-6bee99146ffa
E[end] - E[1]

# ╔═╡ 0c148bd9-2842-4fa0-9d7d-9eb5d5aeb680
L[:, end] .- L[:, 1]

# ╔═╡ Cell order:
# ╟─9a8824ed-e23f-42e1-bac8-f81c1dbcd79e
# ╠═64492622-d5d9-11ef-0503-8de99f2e38e1
# ╠═86278c5b-dac2-4a93-a189-9720c362b9f3
# ╠═f93754f6-7612-4cca-bf77-c1288dcd58f4
# ╠═ab029f52-8d09-486c-8c8d-af0ce4eb13d0
# ╠═0b72ae9c-7713-4733-bc65-75bcf45f5a1b
# ╠═e1e4a9da-2e81-4f63-ba22-2318312261c4
# ╠═48c95a1d-4dd0-434b-a567-aee23a12150a
# ╠═c1ba4866-5432-4879-93d8-1c1108d01521
# ╠═4839120c-63ed-4200-833a-bd171f3c9237
# ╠═26a86ae1-c042-4be0-90e6-4df515f5b2a4
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
# ╠═244bea42-fbb8-41ec-97f5-174b7661b713
# ╠═fc2019d5-641d-44b7-86dc-36f50a87c42a
# ╠═66988f4b-aa3d-4106-b9dd-17cde5c9e529
# ╠═615988d2-5ccd-4ac4-b72d-352eeb122125
# ╠═ddd4c109-a467-46b1-b3dc-99af8d4aea96
# ╠═df85719e-b62d-4d00-b592-78afca6dd7cb
# ╠═eae7a157-a6fb-466a-a293-785139df0ef3
# ╠═4e8f8a44-1d3c-4666-8a64-6e708eb3f165
# ╠═5f35b33a-855d-49f1-bb65-a0e5f60d61cb
# ╠═662af96b-46a5-425a-85c0-42ddc038e304
# ╠═dcb399a2-0416-4b76-8322-f28b0f795406
# ╠═e8c953ba-9a7d-4896-84bf-66559f214b74
# ╠═211bd090-9188-4fcb-9cd3-99f49a325a50
# ╠═882fa476-2c4b-4f40-852d-7c8d55c3fb0a
# ╠═fb9981a0-ef43-4b08-9f1a-e21bab92b3d0
# ╠═e3431072-6206-40e0-9212-aeec0a2b007e
# ╠═87a86e3e-ad3b-485e-8099-b2d2f19857c7
# ╠═8bf93541-5154-4a04-94ea-21df0328848c
# ╠═306020b7-60bd-4bc1-9647-4642b55229b3
# ╠═569f3a8a-7134-4e74-a757-8fbe35f50cfd
# ╠═0a8949c6-3706-4fba-9255-015b72c01546
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
# ╠═94309696-9545-499e-990c-759c4488e1d2
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
# ╠═658280af-5bf9-47b1-aee1-21550580c6a6
# ╠═b80957dc-b24b-4e06-8aad-07018312b370
# ╠═917d02ee-30f2-4338-8d8a-91f9fd5a5044
# ╠═6e59f120-6f5c-4ca1-967d-1034ed4f3473
# ╠═433308c7-ce37-4ace-9008-e1a456774d7e
# ╟─a7265965-9092-451b-bf0b-8f32f8b28858
# ╠═f7a5d076-f97d-41bd-8ad6-a08b6d5ad1a8
# ╠═4dffdfaa-a1be-452d-9c8a-9479f0b7b767
# ╠═05dc4779-28de-4a2e-a47f-231f3e909d9d
# ╠═aa48ab4a-d512-4568-9b46-394a7ef092c5
# ╠═7629ac1a-fc7c-4225-ac2b-752b569e04f1
# ╠═a086ce7a-81b3-4bef-81a3-3e94a436b2e7
# ╠═ba62386c-f7a2-4e15-a4bc-b6f57f06aa03
# ╠═79068d79-b896-45ca-b34f-9e2f6b835a1b
# ╠═5e968467-0a07-4c59-8779-6a8977d707e7
# ╠═2da8ed33-bc9c-4ff1-b1c8-6bee99146ffa
# ╠═0c148bd9-2842-4fa0-9d7d-9eb5d5aeb680
