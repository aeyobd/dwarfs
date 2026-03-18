### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ cc1546bc-2224-11f1-95bd-43f3581a0747
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Agama
	using Arya
	using CairoMakie
end

# ╔═╡ b7ece03b-5a27-4c79-be0e-3d31d9dfa1ca
dwarfs_root = ENV["DWARFS_ROOT"]

# ╔═╡ a80d207c-bd0c-49f8-9852-04408be7a31f
pot = Agama.Potential(file=joinpath(dwarfs_root, "agama/potentials/EP2020.ini"))

# ╔═╡ d8e1431d-686e-4690-8fa8-3b2ed3fcdf7d
action_finder = Agama.ActionFinder(pot)

# ╔═╡ c97c941d-ec95-43d0-ad94-54cd0189aaad
import TOML

# ╔═╡ c2a74848-fe0d-455d-9faa-046395c201d1
obs_props = TOML.parsefile(joinpath(dwarfs_root, "observations/ursa_minor/observed_properties.toml"))

# ╔═╡ 66c2765c-9978-4c77-8530-7de6c582f173
coord_0 = LilGuys.transform(Galactocentric, ICRS(obs_props))

# ╔═╡ 77889fbf-84c7-493b-9831-08f76b0e1d72
pos_i = LilGuys.position(coord_0)

# ╔═╡ b48aa428-ad81-4533-9981-5170e41e14df
vel_i = LilGuys.velocity(coord_0) / V2KMS

# ╔═╡ 6352a779-58c1-4c13-a8c2-b5d3a167186e
# acts, angs, freqs = Agama.actions_angles(action_finder, pos_i, vel_i)

# ╔═╡ 6ac3205b-0b03-4d55-987b-4076450f99a7
integrate = LilGuys.integrate

# ╔═╡ d4dae201-040c-4983-8532-afa9941f7dc1
orbit = LilGuys.agama_orbit(pot, coord_0, timerange=(0, 2_000), N=10_000)

# ╔═╡ a1943076-a14a-4725-b8ec-50fafacaf568
LilGuys.plot_xyz(orbit.positions)

# ╔═╡ 7d7c6202-71b1-449a-aced-7824905cdc4c
Rs = @. sqrt(orbit.positions[1, :]^2 + orbit.positions[2, :]^2)

# ╔═╡ e842269e-5c72-414d-95fd-ef05e82a0a63
zs = orbit.positions[3, :]

# ╔═╡ b0c6d416-eb06-496c-861e-0c6ced6d4867
phis = atan.(orbit.positions[2, :], orbit.positions[1, :])

# ╔═╡ 9db2a3ba-7465-46b4-8b9f-6f4f94aa09ed
vRs = @. (orbit.positions[1, :] * orbit.velocities[1, :] + orbit.positions[2, :] * orbit.velocities[2, :]) / Rs

# ╔═╡ 3e345975-5198-446a-b919-2d94ecc79d6f
vzs = orbit.velocities[3, :]

# ╔═╡ a314f9ee-ba56-44aa-940c-6c3626840816
vphis = @. (-orbit.positions[2, :] * orbit.velocities[1, :] + orbit.positions[1, :] * orbit.velocities[2, :]) ./ Rs

# ╔═╡ 270ec92e-bcf5-47a2-aea5-3fbbca1f376a
let 
	fig = Figure()

	ax = Axis(fig[1,1])
	lines!(orbit.times, Rs)

	vlines!(288 .* vec(1:5) .- 190)

	fig

end

# ╔═╡ 9e90e023-cadf-4fc5-a9fc-8f68d6c3c445
let
	f = lines(orbit.times, zs)
	vlines!(388 .* vec(1:5))

	f
end

# ╔═╡ f781f9ca-983d-4458-84da-01e52566b2a8
rs = LilGuys.radii(orbit)

# ╔═╡ 72d44232-efec-4658-8569-523ac55e59f2
r = LilGuys.lerp(orbit.times, rs)

# ╔═╡ 073fb890-a7ff-402a-87d7-87a811e96683
vrs =  dropdims(sum(orbit.positions .* orbit.velocities, dims=1), dims=1) ./ rs

# ╔═╡ 677fcc3f-a77e-42a4-980e-a6d0fb23a265
vr = LilGuys.lerp(orbit.times, vrs)

# ╔═╡ 6e6dc292-913c-4b07-91d4-2021d7b297b4
1/2π * integrate(t -> vr(t)^2, 50, 50+288)

# ╔═╡ cdfbf9a5-197d-41b1-978d-29a4712f9c26
lines(rs, vrs)

# ╔═╡ 369bce34-200a-4b82-8053-e3c56851b024
vzp = sqrt.(radii(orbit.velocities) .^2 + @. - vrs^2 - vphis^2)

# ╔═╡ cdc638c2-df78-47e2-bcae-4b9316c7befe
lines(zs, vzp)

# ╔═╡ d392ba01-5df2-41e9-b66c-8c9ac451c6f4
lines(zs, vzs)

# ╔═╡ a34a7a71-fe15-4da8-acdc-537474bf888f
vz = LilGuys.lerp(zs[1:600], vzs[1:600])

# ╔═╡ acb69673-5419-4dcd-9107-c0adae98b49c
2/π * integrate(vz, 45, -40)

# ╔═╡ 4434c960-6587-4b08-8f19-23573e4f2764
lines(phis, vphis)

# ╔═╡ d8c4743a-982c-4b83-bc0f-e99c3cd8a2a3
lines(phis)

# ╔═╡ 8cf236d9-8377-4765-bb22-fdb9375f0505
# 1/2π * integrate(t -> R vphi(t)^2, 0, 2π)

# ╔═╡ 0a65db1e-41b3-4673-abe0-395833c15e51
phi = LilGuys.lerp(orbit.times, phis)

# ╔═╡ a1ba500a-b999-4056-8aad-96a385b38bd8
vphi = LilGuys.lerp(phis, vphis)

# ╔═╡ 4d37bd00-7e60-4108-8c72-a6cf2730b7ff
md"""
# Stackel approximation
"""

# ╔═╡ ac7a589a-d407-4cea-b422-e0901abe3c71


# ╔═╡ 2f7ebeab-54ac-4aeb-b129-836c4d379c8f
pϕs = Rs .* vphis

# ╔═╡ 14591915-69e6-4e82-bc0b-a0dd05ec37c9
pϕs[1]

# ╔═╡ a80779cf-3515-4556-83c3-9ca7e94b9728
md"""
## action angles
"""

# ╔═╡ 459edda8-978f-4406-83e0-a85e2d50b252
Δ = 3

# ╔═╡ 7b1d8b1a-fd7e-4be0-89bf-5d57fd01ffc5
to_λ(R, z) = 1/2 * (R^2 + z^2 + Δ^2) + 1/2 * sqrt((R^2 + z^2 - Δ^2)^2 + 4R^2 * Δ^2)

# ╔═╡ b669732d-081c-4970-9608-4869aefdeddc
λs = to_λ.(Rs, zs)

# ╔═╡ 77faee53-808d-41ad-82d2-056a9eeb9040
idx_λ = argmin(λs[1:1500]):argmax(λs[1:1500])

# ╔═╡ 8735627a-5067-4a17-ba34-473ac9878521
let
	f = lines(λs)
	lines!(idx_λ, λs[idx_λ])
	f

end

# ╔═╡ 11aeb1fb-d905-4e0c-9683-04c1b0fe1307
λ0, λ1 = extrema(λs)

# ╔═╡ 986ac219-ab2e-4f47-ba0b-187809ac2b82
to_ν(R, z) = 1/2 * (R^2 + z^2 + Δ^2) - 1/2 * sqrt((R^2 + z^2 - Δ^2)^2 + 4R^2 * Δ^2)

# ╔═╡ a930ae2f-734a-456b-bb81-d2b5bff9eba5
νs = to_ν.(Rs, zs)

# ╔═╡ bd8dfe65-eef3-4ada-ae25-610291ebe6f7
idx_ν = (argmin(νs[800:2000]):argmax(νs[800:2000])) .+ 799

# ╔═╡ 14a12214-c96f-4cfd-9e5b-d622767793e9
let
	f = lines(νs)

	lines!(idx_ν, νs[idx_ν])
	f
end

# ╔═╡ 5cf9eecf-7007-48d7-955d-2f40c35dfec0
let
	f = lines(λs, νs)
	hlines!([minimum(νs), maximum(νs)])
	vlines!([minimum(λs), maximum(λs)])
	f
end

# ╔═╡ e5111573-a9a4-4c9e-b3bf-8ec85ad53b16
ν0, ν1 = 1e-12, maximum(νs)

# ╔═╡ affa10a6-e932-4be2-95e8-310e8b415d55
function Φλ(λ, ν=0)
	R = sqrt((λ - Δ^2)*(Δ^2 - ν)/Δ^2)
	z = sqrt(λ * ν / Δ^2)
	
	Agama.potential(pot, [R, 0, z])
end

# ╔═╡ 5cc07662-2a9c-42a2-8f20-219199b2caea
Φλ(10)

# ╔═╡ aa7281e4-68f7-4650-bd34-b26c0d91786f
Rs[100]^2 / (λs[100] - Δ^2) + zs[100]^2 / (λs[100])

# ╔═╡ a33ff3f8-3d52-4b33-9fd3-89c630c709ce
@. (λs - Δ^2) * (Δ^2 - νs) / Δ^2 / Rs^2

# ╔═╡ 8f521fba-e4b3-4d0a-a23f-07f8591083b4
@. (λs * νs) / Δ^2 / zs^2

# ╔═╡ 103105bd-081c-46fe-878c-5eb8e45cdfa2
function to_pλ(R, z, vR, vz)
	λ = to_λ(R, z)
	return R*vR / (2*(λ - Δ^2)) + z*vz / 2λ
end

# ╔═╡ d14e496f-de5f-4b29-9b57-bcc4226663d6
pλs = to_pλ.(Rs, zs, vRs, vzs)

# ╔═╡ 1735da39-cb04-4f3e-a37f-8d81da91dec6
lines(λs, pλs)

# ╔═╡ ef8e8e8a-6945-4d02-9049-80ed156282b7
pλ = LilGuys.lerp(λs[idx_λ], pλs[idx_λ])

# ╔═╡ 535a37f9-caf5-4818-81da-4656c3bab2b6
let
	f = lines(λs, pλs)
	lines!(λs, pλ.(λs))
	f

end

# ╔═╡ d1a4f07e-3671-4b00-b118-d4172ad50b00
integrate(pλ, minimum(λs), maximum(λs)) / π

# ╔═╡ bc233951-6c09-4125-9f9f-8f2474c35451
function to_pν(R, z, vR, vz)
	ν = to_ν(R, z)
	return R*vR / (2*(ν - Δ^2)) + z*vz / 2ν
end

# ╔═╡ 96ad9467-9ade-4f59-a893-7ac26a4d4bda
pνs = to_pν.(Rs, zs, vRs, vzs)

# ╔═╡ 3ac8e98e-6080-4fea-a83d-8e84b6893121
pν = LilGuys.lerp(νs[idx_ν], pνs[idx_ν])

# ╔═╡ 39465a61-b07c-4df9-90ad-fa812f09bed2
2integrate(pν, 0, maximum(νs)) / π

# ╔═╡ fd9022f3-3c61-4fe6-acd7-90def278222d
let
	f = lines(νs, pνs)
	ylims!(-1000Δ^-2, 1000Δ^-2)
	f
end

# ╔═╡ d891cd3a-82d5-4bf8-9cd4-8e1c9f3a2db3
Lz = pϕs[1]

# ╔═╡ f2d15eb8-d750-41d5-966e-d5322e4a5dbd
function integrate_action(f, λ1, λ2)
	λ(x) = (λ2 - λ1)/2 * sin(x) + (λ1+λ2)/2

	return integrate(x -> (λ2 - λ1)/2 * cos(x) * f(λ(x)), -π/2, π/2)
end

# ╔═╡ 87e1da82-4ebd-4cc3-b708-92e593b774b0
dJλ_dE = 1/(4π) * integrate_action(λ -> 1 / ((λ - Δ^2)*pλ(λ)), λ0, λ1)

# ╔═╡ cb0c9c01-9ce2-4fcd-a60e-c01933ccc614
1/dJλ_dE

# ╔═╡ 0c216ef2-f564-4f99-833d-ad629134f5bc
dJλ_dLz = -Lz/4π * integrate_action(λ -> 1/((λ - Δ^2)^2*pλ(λ)), λ0, λ1)

# ╔═╡ 712aac03-4e14-49f7-933e-a8457ff7744d
dJλ_dI3 = -1/4π * integrate_action(λ -> 1/((λ - Δ^2)*λ*pλ(λ)), λ0, λ1)

# ╔═╡ 8d9b01c9-ab40-4538-9353-b6486aad0d72
am = Agama.ActionMapper(pot)

# ╔═╡ 232462d9-f16b-4ed2-ba8b-4044efebefca
f(x) = -Φλ(x) * x

# ╔═╡ ddf92d71-f70e-40c8-b518-3c2049db6a94
1/dJλ_dE

# ╔═╡ 74425e17-4439-4400-bfe6-a6e8cd1b0e11
dJν_dE = 1/2π * integrate_action(λ -> 1/((λ - Δ^2)*pν(λ)), ν0, ν1)

# ╔═╡ 7b305db5-ae84-47ff-9f1a-3dc3c15a6266
1/dJν_dE

# ╔═╡ 7f31063b-42c6-4603-8ae9-9234a0f127dc
1/dJν_dE

# ╔═╡ 595a128c-cf06-4b54-8436-ef90d48e1d29
dJν_dLz = -Lz/2π * integrate_action(λ -> 1/((λ - Δ^2)^2*pν(λ)), ν0, ν1)

# ╔═╡ 88edee43-3a67-4885-9c9f-66259970230a
dJν_dI3 = -1/2π * integrate_action(λ -> 1/((λ - Δ^2)*λ*pν(λ)), ν0, ν1)

# ╔═╡ 7d86e4a6-decb-4a56-a265-e039d9d44e09
Jac = [
	dJλ_dE dJλ_dLz dJλ_dI3
	dJν_dE dJν_dLz dJν_dI3
	0 1 0
]

# ╔═╡ f2e46d1c-ff80-4df3-ac18-9513b65e287d
invJac = inv(Jac)

# ╔═╡ 2067020c-a4d6-4b7c-b89b-842c82fc1081
invJac

# ╔═╡ 65c8b16e-7ad5-4631-8b4b-b3bfde660277
invJac[1, :]

# ╔═╡ 8f5982d1-82c1-463b-af42-0695de1c3aa1
idx = 1 +16*5

# ╔═╡ dcfa6ede-5af2-4c55-a921-c46f646be2a9
E = Agama.potential(pot, orbit.positions[:, idx]) + radii(orbit.velocities[:, idx])^2/2

# ╔═╡ cfb02f12-23e5-4ff8-9d65-6301871f5fec
I3 = (Φλ(λs[idx], νs[idx]) - Φλ(λs[idx])) * λs[idx] + 1/2 * (
	zs[idx]^2 * vphis[idx]^2 
		+ (Rs[idx]*vzs[idx] - zs[idx] * vRs[idx])^2 
		+ vzs[idx]^2 *Δ^2)

# ╔═╡ a1e7aa30-fc20-428d-95d4-7e069a0b1897
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!(@. 2*(λs - Δ^2) * pλs^2)
	x = @. E - Lz^2 / (2*λs - Δ^2) - I3 / (λs) + f(λs) / (λs - Δ^2)
	lines!(x, linestyle=:dot)
	fig

end

# ╔═╡ 7872473f-44fc-4ce8-bf48-57652889e139
dS_dE = 1/4 * integrate(λ -> 1/((λ - Δ^2) * pλ(λ)), λ0, λs[idx]) +
	1/4 * integrate_action(λ -> 1/((λ - Δ^2) * pν(λ)), ν0, min(ν1, νs[idx]))

# ╔═╡ 3f7acdda-7702-4cb5-a543-5f16bbab6bac
dS_dLz = phis[idx] - Lz/4 * (integrate(λ -> 1/((λ - Δ^2)^2 * pλ(λ)), λ0, λs[idx]) +
	integrate_action(λ -> 1/((λ - Δ^2)^2 * pν(λ)), ν0, min(ν1, νs[idx]))
						  )

# ╔═╡ a3358660-8df8-4ee3-9fd2-b88cb02c3d43
dS_dI3 = 1/4 * (integrate(λ -> 1/((λ - Δ^2)*λ * pλ(λ)), λ0, λs[idx]) +
	integrate_action(λ -> 1/((λ - Δ^2)*λ * pν(λ)), ν0, min(ν1, νs[idx]))
						  )

# ╔═╡ 86e81fc8-811d-4c76-bb21-3f14b9ce6e98
invJac[:, 1]' * [dS_dE, dS_dLz, dS_dI3]

# ╔═╡ 3c1d1d6b-8204-4913-a238-91f81257a803
invJac[:, 2]' * [dS_dE, dS_dLz, dS_dI3]

# ╔═╡ 3a7c3fa3-d856-44c9-88b6-893a9c5a3679
invJac[:, 3]' * [dS_dE, dS_dLz, dS_dI3]

# ╔═╡ 8cc7ac6c-93eb-4cd2-9471-ce8290e8b2fa
diff(orbit.times)[1]

# ╔═╡ b25f11b8-1ced-431a-8adb-5d96ee7e0579
zs[1]

# ╔═╡ f6e60d7a-d5b3-4e12-9fda-4d2478c9dbe8
lines(phis)

# ╔═╡ a6256269-03b3-4fbe-80b0-2668b21098d4
acts, angs, freqs = Agama.actions_angles(action_finder, orbit.positions[:, idx], orbit.velocities[:, idx])

# ╔═╡ 9f0b2ca5-b86e-4a7d-b96e-497f0337d13c
acts

# ╔═╡ 68f91133-0e2e-43e4-9b11-c7ae1e7fcc4f
2π ./ freqs

# ╔═╡ 42395cc1-ff5d-47de-a346-2e6d5886bf41
acts[2]

# ╔═╡ 62691d87-aea2-4ed9-93e0-d946b79e6f11
acts[2]

# ╔═╡ 0522ee88-8aad-461e-9b6f-08e9ef38a64a
acts[3]

# ╔═╡ e5adfd56-a9f6-4498-bdb4-bf82445acd10
acts

# ╔═╡ 13c0de70-b3f8-4ec7-bb53-cfa126bbe00e
freqs

# ╔═╡ 8238d14a-f558-4d80-a689-3ba28ec0add5
freqs

# ╔═╡ 93dfce47-92bb-4050-ae58-4a285a2cf092
let 
	h = 1
	act_new = acts .+ [0, 0, h]
	pos, vel = Agama.from_actions(am, act_new, angs)

	E(pos, vel) = Agama.potential(pot, pos) + radii(vel)^2/2
	(E(pos, vel) - E(pos_i, vel_i) )/ h
end

# ╔═╡ f5fe0687-b2ab-40d5-8429-1de6189df8cf
freqs

# ╔═╡ 42a495ea-f2ee-499a-8d7a-260dad6ad031
freqs

# ╔═╡ af7183a8-5da8-4811-a0da-64d395abde68
angs

# ╔═╡ 232472ea-08d4-4881-8004-e5f235e2a371
angs

# ╔═╡ Cell order:
# ╠═cc1546bc-2224-11f1-95bd-43f3581a0747
# ╠═b7ece03b-5a27-4c79-be0e-3d31d9dfa1ca
# ╠═a80d207c-bd0c-49f8-9852-04408be7a31f
# ╠═d8e1431d-686e-4690-8fa8-3b2ed3fcdf7d
# ╠═c97c941d-ec95-43d0-ad94-54cd0189aaad
# ╠═c2a74848-fe0d-455d-9faa-046395c201d1
# ╠═66c2765c-9978-4c77-8530-7de6c582f173
# ╠═77889fbf-84c7-493b-9831-08f76b0e1d72
# ╠═b48aa428-ad81-4533-9981-5170e41e14df
# ╠═6352a779-58c1-4c13-a8c2-b5d3a167186e
# ╠═6ac3205b-0b03-4d55-987b-4076450f99a7
# ╠═d4dae201-040c-4983-8532-afa9941f7dc1
# ╠═a1943076-a14a-4725-b8ec-50fafacaf568
# ╠═7d7c6202-71b1-449a-aced-7824905cdc4c
# ╠═e842269e-5c72-414d-95fd-ef05e82a0a63
# ╠═b0c6d416-eb06-496c-861e-0c6ced6d4867
# ╠═9db2a3ba-7465-46b4-8b9f-6f4f94aa09ed
# ╠═3e345975-5198-446a-b919-2d94ecc79d6f
# ╠═a314f9ee-ba56-44aa-940c-6c3626840816
# ╠═270ec92e-bcf5-47a2-aea5-3fbbca1f376a
# ╠═6e6dc292-913c-4b07-91d4-2021d7b297b4
# ╠═9f0b2ca5-b86e-4a7d-b96e-497f0337d13c
# ╠═72d44232-efec-4658-8569-523ac55e59f2
# ╠═677fcc3f-a77e-42a4-980e-a6d0fb23a265
# ╠═68f91133-0e2e-43e4-9b11-c7ae1e7fcc4f
# ╠═9e90e023-cadf-4fc5-a9fc-8f68d6c3c445
# ╠═f781f9ca-983d-4458-84da-01e52566b2a8
# ╠═073fb890-a7ff-402a-87d7-87a811e96683
# ╠═cdfbf9a5-197d-41b1-978d-29a4712f9c26
# ╠═369bce34-200a-4b82-8053-e3c56851b024
# ╠═cdc638c2-df78-47e2-bcae-4b9316c7befe
# ╠═d392ba01-5df2-41e9-b66c-8c9ac451c6f4
# ╠═a34a7a71-fe15-4da8-acdc-537474bf888f
# ╠═acb69673-5419-4dcd-9107-c0adae98b49c
# ╠═42395cc1-ff5d-47de-a346-2e6d5886bf41
# ╠═4434c960-6587-4b08-8f19-23573e4f2764
# ╠═62691d87-aea2-4ed9-93e0-d946b79e6f11
# ╠═d8c4743a-982c-4b83-bc0f-e99c3cd8a2a3
# ╠═8cf236d9-8377-4765-bb22-fdb9375f0505
# ╠═0522ee88-8aad-461e-9b6f-08e9ef38a64a
# ╠═0a65db1e-41b3-4673-abe0-395833c15e51
# ╠═a1ba500a-b999-4056-8aad-96a385b38bd8
# ╠═4d37bd00-7e60-4108-8c72-a6cf2730b7ff
# ╠═cb0c9c01-9ce2-4fcd-a60e-c01933ccc614
# ╠═7b305db5-ae84-47ff-9f1a-3dc3c15a6266
# ╠═7b1d8b1a-fd7e-4be0-89bf-5d57fd01ffc5
# ╠═986ac219-ab2e-4f47-ba0b-187809ac2b82
# ╠═affa10a6-e932-4be2-95e8-310e8b415d55
# ╠═ac7a589a-d407-4cea-b422-e0901abe3c71
# ╠═aa7281e4-68f7-4650-bd34-b26c0d91786f
# ╠═5cc07662-2a9c-42a2-8f20-219199b2caea
# ╠═a33ff3f8-3d52-4b33-9fd3-89c630c709ce
# ╠═8f521fba-e4b3-4d0a-a23f-07f8591083b4
# ╠═103105bd-081c-46fe-878c-5eb8e45cdfa2
# ╠═bc233951-6c09-4125-9f9f-8f2474c35451
# ╠═b669732d-081c-4970-9608-4869aefdeddc
# ╠═a930ae2f-734a-456b-bb81-d2b5bff9eba5
# ╠═2f7ebeab-54ac-4aeb-b129-836c4d379c8f
# ╠═d14e496f-de5f-4b29-9b57-bcc4226663d6
# ╠═96ad9467-9ade-4f59-a893-7ac26a4d4bda
# ╠═1735da39-cb04-4f3e-a37f-8d81da91dec6
# ╠═fd9022f3-3c61-4fe6-acd7-90def278222d
# ╠═8735627a-5067-4a17-ba34-473ac9878521
# ╠═535a37f9-caf5-4818-81da-4656c3bab2b6
# ╠═77faee53-808d-41ad-82d2-056a9eeb9040
# ╠═14a12214-c96f-4cfd-9e5b-d622767793e9
# ╠═bd8dfe65-eef3-4ada-ae25-610291ebe6f7
# ╠═3ac8e98e-6080-4fea-a83d-8e84b6893121
# ╠═ef8e8e8a-6945-4d02-9049-80ed156282b7
# ╠═5cf9eecf-7007-48d7-955d-2f40c35dfec0
# ╠═11aeb1fb-d905-4e0c-9683-04c1b0fe1307
# ╠═e5111573-a9a4-4c9e-b3bf-8ec85ad53b16
# ╠═d1a4f07e-3671-4b00-b118-d4172ad50b00
# ╠═39465a61-b07c-4df9-90ad-fa812f09bed2
# ╠═14591915-69e6-4e82-bc0b-a0dd05ec37c9
# ╠═e5adfd56-a9f6-4498-bdb4-bf82445acd10
# ╟─a80779cf-3515-4556-83c3-9ca7e94b9728
# ╠═459edda8-978f-4406-83e0-a85e2d50b252
# ╠═d891cd3a-82d5-4bf8-9cd4-8e1c9f3a2db3
# ╠═7d86e4a6-decb-4a56-a265-e039d9d44e09
# ╠═13c0de70-b3f8-4ec7-bb53-cfa126bbe00e
# ╠═2067020c-a4d6-4b7c-b89b-842c82fc1081
# ╠═87e1da82-4ebd-4cc3-b708-92e593b774b0
# ╠═0c216ef2-f564-4f99-833d-ad629134f5bc
# ╠═712aac03-4e14-49f7-933e-a8457ff7744d
# ╠═f2d15eb8-d750-41d5-966e-d5322e4a5dbd
# ╠═8d9b01c9-ab40-4538-9353-b6486aad0d72
# ╠═8238d14a-f558-4d80-a689-3ba28ec0add5
# ╠═dcfa6ede-5af2-4c55-a921-c46f646be2a9
# ╠═cfb02f12-23e5-4ff8-9d65-6301871f5fec
# ╠═232462d9-f16b-4ed2-ba8b-4044efebefca
# ╠═a1e7aa30-fc20-428d-95d4-7e069a0b1897
# ╠═93dfce47-92bb-4050-ae58-4a285a2cf092
# ╠═ddf92d71-f70e-40c8-b518-3c2049db6a94
# ╠═7f31063b-42c6-4603-8ae9-9234a0f127dc
# ╠═f5fe0687-b2ab-40d5-8429-1de6189df8cf
# ╠═74425e17-4439-4400-bfe6-a6e8cd1b0e11
# ╠═595a128c-cf06-4b54-8436-ef90d48e1d29
# ╠═88edee43-3a67-4885-9c9f-66259970230a
# ╠═7872473f-44fc-4ce8-bf48-57652889e139
# ╠═3f7acdda-7702-4cb5-a543-5f16bbab6bac
# ╠═a3358660-8df8-4ee3-9fd2-b88cb02c3d43
# ╠═f2e46d1c-ff80-4df3-ac18-9513b65e287d
# ╠═65c8b16e-7ad5-4631-8b4b-b3bfde660277
# ╠═86e81fc8-811d-4c76-bb21-3f14b9ce6e98
# ╠═3c1d1d6b-8204-4913-a238-91f81257a803
# ╠═3a7c3fa3-d856-44c9-88b6-893a9c5a3679
# ╠═8f5982d1-82c1-463b-af42-0695de1c3aa1
# ╠═8cc7ac6c-93eb-4cd2-9471-ce8290e8b2fa
# ╠═42a495ea-f2ee-499a-8d7a-260dad6ad031
# ╠═af7183a8-5da8-4811-a0da-64d395abde68
# ╠═b25f11b8-1ced-431a-8adb-5d96ee7e0579
# ╠═f6e60d7a-d5b3-4e12-9fda-4d2478c9dbe8
# ╠═a6256269-03b3-4fbe-80b0-2668b21098d4
# ╠═232472ea-08d4-4881-8004-e5f235e2a371
