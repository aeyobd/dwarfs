### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ fa66d3a2-1ed6-11f0-187d-efbb138ffe4e
begin
	using Pkg; Pkg.activate()

	using LilGuys
	using Arya

	using CairoMakie
end

# ╔═╡ e228d9d1-ca40-41ca-bc00-76c7b77321d1
using DataFrames, Turing, PairPlots

# ╔═╡ c949368b-4526-4793-9406-99f7e4faab33
using Distributions

# ╔═╡ 9858575a-12f1-4c93-9bdd-7545300d947d
md"""
# Simulation RV study
"""

# ╔═╡ 9e497dad-ddb5-47ea-a3f8-e7c2f39638af
md"""
The goal of this notebook is to validate the RV methods used and determine what the resulting parameter derivation appears as using a mock sample.
"""

# ╔═╡ 6a34142b-64bd-4e7c-8914-63162699fd7a
import TOML

# ╔═╡ bab02077-c470-4716-a51f-5daea48a0ee8
import Measurements

# ╔═╡ 0e9a6e6e-3375-41e9-be6e-36f8d5263e6e
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ a784a9f5-9a97-4afd-bff9-69a0d20b7455
simulation_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "sculptor", "1e7_V31_r4.2")

# ╔═╡ 93d3e14e-114a-4bc1-8e13-ad1434e2ce86
⊕ = RVUtils.:⊕

# ╔═╡ fc1f2053-8d56-4c28-abea-13c02acd32ef
md"""
# Model creation
"""

# ╔═╡ 8a0f298d-fb85-45b4-93e3-c25d748ae742
pot = NFW(r_s=2.7, M_s=0.4)

# ╔═╡ d72cd1ad-230a-466b-b169-b1ce991785d7
prof = LilGuys.Exp2D(R_s=0.13)

# ╔═╡ febf1579-131f-4512-bb16-07f33a4fad95
Ψ(x) = -LilGuys.potential(pot, x)

# ╔═╡ ad80aad2-40de-495c-9a74-453e94bd076f
ρ(x) = LilGuys.density(prof, x)

# ╔═╡ 3d8dbe61-4d3d-4655-8247-25ef4c2bc9b0
r = 10 .^ LinRange(-5, 2, 1000)

# ╔═╡ ebedb175-9bc7-4929-a306-c5a9061d1ded
f_s = LilGuys.DistributionFunction(ρ.(r), Ψ.(r), r)

# ╔═╡ 9d953cb4-24d9-4756-8264-7762b67d9b9b
f_interp = LilGuys.lerp(Ψ.(r), f_s.(Ψ.(r)))

# ╔═╡ d4d90de3-3c5e-4bf8-88ae-7659a42d2355
Es = LinRange(Ψ(r[end]), Ψ(r[1]), 1000)

# ╔═╡ abc08e57-db68-4533-ba0a-e1d235e29c45
lines(Es, log10.(f_s.(Es)))

# ╔═╡ 8824a4d7-77c7-4605-b0ea-bea4bc532923
Npoints = 1_000

# ╔═╡ f0577c3e-a7f7-4455-9be2-9c62600819ac
rs = LilGuys.sample_density(ρ, Npoints)

# ╔═╡ d4d530c5-bab7-48c6-a863-e0cbaee77380
vs_model = LinRange(0, sqrt(2*Ψ(r[1])), 100)

# ╔═╡ dac56962-3235-487e-8b23-e6d267bb38ee
lines(vs_model, @. vs_model^2 * f_interp(Ψ(r[1]) - vs_model^2/2))

# ╔═╡ c3fe8ad8-31c7-447b-a26d-7db6d90d5ecc
function sample_velocity(r)
	psi = Ψ(r)
	v_max = sqrt(2psi)
	vs = LinRange(0, v_max, 100)
	f_vs = @. vs^2 * f_interp(psi - vs^2/2)
	cdf_vs = cumsum(f_vs) # constant stepsize
	cdf_vs ./= cdf_vs[end]

	v = LilGuys.lerp(cdf_vs, vs)(rand())
end

# ╔═╡ 55b72bd8-6ed5-4d5f-97fd-aca74b50c749
vs = sample_velocity.(rs)

# ╔═╡ 6cabf5bb-cde7-442f-a32e-252df3566042
sv_exp = sqrt(LilGuys.mean(vs .^ 2)) * V2KMS / √3

# ╔═╡ 3cb6335a-da69-404f-9336-1868929e5c5a
positions = rs' .* LilGuys.rand_unit(Npoints)

# ╔═╡ feb20ceb-9f0a-45ee-8535-1df7b0c796f2
velocities = vs' .* LilGuys.rand_unit(Npoints)

# ╔═╡ f218c6ce-1638-41dc-99ef-61dafa080eea
hist(velocities[1, :] * V2KMS)

# ╔═╡ 6630955e-b607-43c1-a677-6a5f562259bd
LilGuys.std(velocities[1, :] * V2KMS)

# ╔═╡ 2e48c25b-9de3-425c-9e8b-3ec68d28dc86
obs_props = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 423b1299-ccf8-466f-b8b4-347c735bd208
icrs = ICRS(obs_props)

# ╔═╡ 9c03ae01-39a1-47b3-aafb-1cd931c37d46
gc0 = LilGuys.transform(Galactocentric, icrs)

# ╔═╡ f43f151f-c662-43cc-acc2-9b42fbeb94a7
gcs = [Galactocentric(positions[:, i] .+ LilGuys.position(gc0), velocities[:, i] * V2KMS .+ LilGuys.velocity(gc0)) for i in 1:Npoints]

# ╔═╡ 40950a02-2767-48c8-97c3-35388c7acac7
icrss = LilGuys.transform.(ICRS, gcs)

# ╔═╡ ebca98d3-8f73-453c-a2bf-1b5219ac8f52
df = LilGuys.to_frame(icrss)

# ╔═╡ 4612abf8-9a80-45c3-8544-2a2dc70b8162
md"""
# Observing sample
"""

# ╔═╡ 7f5efe2c-6379-41be-b682-5734d2f5c1dc
dist_RV_err = LogUniform(0.3, 2)

# ╔═╡ 814dd267-bc58-4bee-8eef-aa35615dd855
icrs0 = RVUtils.icrs(obs_props)

# ╔═╡ ec4e41a4-51df-427a-b0e6-7532ec54d215
gsr0 = LilGuys.transform(GSR, icrs0)

# ╔═╡ 55a3f21e-e452-4c3e-8292-e15aab008fcf
df_obs = let
	df_out = copy(df)

	df_out[:, :RV_err] = rand(dist_RV_err, Npoints)
	df_out[:, :RV] = df_out.radial_velocity .+ randn(Npoints) .* df_out.RV_err

	
	
	Δrv = RVUtils.rv_correction(df_out.ra, df_out.dec, gsr0)
	df_out[:, :delta_rv] = Measurements.value.(Δrv)
	df_out[:, :delta_rv_err] = Measurements.uncertainty.(Δrv)
	
	RVUtils.add_gsr!(df_out, 
					 distance=obs_props["distance"],
					pmra=obs_props["pmra"], pmdec=obs_props["pmdec"])
	
	df_out[:, :vz] .= df_out.radial_velocity_gsr .+ df_out.delta_rv
	df_out[:, :vz_err] .= df_out.RV_err .⊕ df_out.delta_rv_err
	
	
	df_out

end

# ╔═╡ 669ab918-d3c3-488d-9c4f-d44dc10f9334
md"""
# Observational methods
"""

# ╔═╡ 956e4eff-ddd7-4573-8be3-c87a733e0d87
md"""
We expect

- vlos = $(gsr0.radial_velocity)
- sigma v = $(round(sv_exp, digits=2))
"""

# ╔═╡ d3eb02fa-d437-464d-ab6d-8e5f7769bf7c
samples = sample(RVUtils.model_vel_1c(df_obs.vz, df_obs.vz_err), NUTS(0.65), MCMCThreads(), 1000, 16)

# ╔═╡ ec033ffb-d52d-4d0b-8b4f-fe98fcad2f4b
pairplot(samples)

# ╔═╡ 9e7f681d-9551-4ede-870d-059b39c5318a
md"""
## Gradient model
we expect the exact same inference as the normal model and both gradient parameters to be zero.
"""

# ╔═╡ f51b9db7-c304-44eb-8ec7-af5e448bb321
xi, eta = [x * 60 for x in LilGuys.to_tangent(df_obs.ra, df_obs.dec, icrs.ra, icrs.dec)]

# ╔═╡ 533338cc-983f-4852-86e9-da68564f1b9c
model_gradient = RVUtils.model_vel_gradient(df_obs.vz, df_obs.vz_err, xi, eta)

# ╔═╡ fa9de589-25d3-4e0e-bbe3-9e48881f8246
samples_gradient = sample(model_gradient, NUTS(0.65), MCMCThreads(), 1000, 16)

# ╔═╡ d6ac7d8e-0628-4822-9a43-ec48cc929c51
pairplot(samples_gradient)

# ╔═╡ 2abfb3e3-31ff-43e9-bbb0-2d5d37c70f88
RVUtils.bayes_evidence(model_gradient, DataFrame(samples_gradient), ["A", "B"])

# ╔═╡ 994b6e87-5f3d-49ae-910b-48fc39ffe4fd
md"""
## Rell - sigma
"""

# ╔═╡ 14c22e16-6345-40dc-9a5b-984ea548aba6
LilGuys.centroid(velocities) *V2KMS

# ╔═╡ 40f43bff-49ad-4a2d-95dd-6583bb1f3175
R_ell = @. xi ⊕ eta

# ╔═╡ 6cf2126b-d9f3-41b1-b14f-138ad2bb6f90
model_Rell = RVUtils.model_vel_sigma_R(df_obs.vz, df_obs.vz_err, R_ell)

# ╔═╡ 8962ba6b-31d4-4059-9a32-57044aafacfb
samples_Rell = sample(model_Rell, NUTS(0.65), MCMCThreads(), 1_000, 16)

# ╔═╡ 07e5b191-178a-4e81-ab06-8f03fcefc6df
pairplot(samples_Rell)

# ╔═╡ 505fffb4-1754-43aa-8c26-a42bef3ad89b
RVUtils.bayes_evidence(model_Rell, DataFrame(samples_Rell), "dlσ_dlR")

# ╔═╡ a2e27b30-3c9f-4c60-b71b-5355e69d8483
md"""
# Density profile
"""

# ╔═╡ 41317a1e-e475-4d96-86ce-5ae969567bbe
scatter(xi, eta)

# ╔═╡ 793492fc-5cb6-4931-aa5f-05538ddb6bf6
scatter(xi, eta, color=df_obs.radial_velocity, colormap=:bluesreds, markersize=3)

# ╔═╡ d5cdacc5-7c4b-49af-98c8-2629602d61b6
scatter(xi, eta, color=df_obs.vz, colormap=:bluesreds, markersize=3)

# ╔═╡ aff8bffa-a647-46b2-a771-7b7d236fcf54
density_prof = LilGuys.StellarDensityProfile(R_ell, bins=LinRange(-0.1, 2, 20))

# ╔═╡ 4b43eeb7-2369-49cd-a825-5d4f0d80a2ce
R_sim = @. positions[1, :] ⊕ positions[2, :]

# ╔═╡ 94e4ca9c-eb97-4712-abe1-5c27c99259db
density_prof_1d = LilGuys.StellarDensityProfile(R_sim, bins=LinRange(-2.1, 0.0, 20))

# ╔═╡ de27fdff-6dfb-4ad9-9d71-2706181c1362
arcmin_to_kpc = LilGuys.arcmin2kpc(1, icrs.distance)

# ╔═╡ 053081de-8dd0-4986-8620-6b0b2fc59e35
let
	fig = Figure()
	ax = Axis(fig[1,1])

	errorscatter!(density_prof.log_R, density_prof.log_Sigma, yerror=LilGuys.error_interval.(density_prof.log_Sigma))

	x = LinRange(-1, 1.8, 1000)
	y = @. log10(LilGuys.surface_density(prof, 10 .^ x .* arcmin_to_kpc)) .+ log10(Npoints * arcmin_to_kpc^2)

	lines!(x, y)

	fig
end

# ╔═╡ 9665c806-9afd-42fe-a619-69f7f09abe6e
let
	fig = Figure()
	ax = Axis(fig[1,1])

	errorscatter!(density_prof_1d.log_R, density_prof_1d.log_Sigma, yerror=LilGuys.error_interval.(density_prof_1d.log_Sigma))

	x = LinRange(-2, 0, 1000)
	y = @. log10(LilGuys.surface_density(prof, 10 .^ x)) .+ log10(Npoints)

	lines!(x, y)

	fig
end

# ╔═╡ Cell order:
# ╠═9858575a-12f1-4c93-9bdd-7545300d947d
# ╠═9e497dad-ddb5-47ea-a3f8-e7c2f39638af
# ╠═fa66d3a2-1ed6-11f0-187d-efbb138ffe4e
# ╠═6a34142b-64bd-4e7c-8914-63162699fd7a
# ╠═bab02077-c470-4716-a51f-5daea48a0ee8
# ╠═e228d9d1-ca40-41ca-bc00-76c7b77321d1
# ╠═0e9a6e6e-3375-41e9-be6e-36f8d5263e6e
# ╠═a784a9f5-9a97-4afd-bff9-69a0d20b7455
# ╠═93d3e14e-114a-4bc1-8e13-ad1434e2ce86
# ╟─fc1f2053-8d56-4c28-abea-13c02acd32ef
# ╠═8a0f298d-fb85-45b4-93e3-c25d748ae742
# ╠═d72cd1ad-230a-466b-b169-b1ce991785d7
# ╠═febf1579-131f-4512-bb16-07f33a4fad95
# ╠═ad80aad2-40de-495c-9a74-453e94bd076f
# ╠═3d8dbe61-4d3d-4655-8247-25ef4c2bc9b0
# ╠═ebedb175-9bc7-4929-a306-c5a9061d1ded
# ╠═9d953cb4-24d9-4756-8264-7762b67d9b9b
# ╠═d4d90de3-3c5e-4bf8-88ae-7659a42d2355
# ╠═abc08e57-db68-4533-ba0a-e1d235e29c45
# ╠═8824a4d7-77c7-4605-b0ea-bea4bc532923
# ╠═f0577c3e-a7f7-4455-9be2-9c62600819ac
# ╠═d4d530c5-bab7-48c6-a863-e0cbaee77380
# ╠═dac56962-3235-487e-8b23-e6d267bb38ee
# ╠═c3fe8ad8-31c7-447b-a26d-7db6d90d5ecc
# ╠═55b72bd8-6ed5-4d5f-97fd-aca74b50c749
# ╠═6cabf5bb-cde7-442f-a32e-252df3566042
# ╠═3cb6335a-da69-404f-9336-1868929e5c5a
# ╠═feb20ceb-9f0a-45ee-8535-1df7b0c796f2
# ╠═f218c6ce-1638-41dc-99ef-61dafa080eea
# ╠═6630955e-b607-43c1-a677-6a5f562259bd
# ╠═2e48c25b-9de3-425c-9e8b-3ec68d28dc86
# ╠═423b1299-ccf8-466f-b8b4-347c735bd208
# ╠═9c03ae01-39a1-47b3-aafb-1cd931c37d46
# ╠═f43f151f-c662-43cc-acc2-9b42fbeb94a7
# ╠═40950a02-2767-48c8-97c3-35388c7acac7
# ╠═ebca98d3-8f73-453c-a2bf-1b5219ac8f52
# ╠═4612abf8-9a80-45c3-8544-2a2dc70b8162
# ╠═c949368b-4526-4793-9406-99f7e4faab33
# ╠═7f5efe2c-6379-41be-b682-5734d2f5c1dc
# ╠═814dd267-bc58-4bee-8eef-aa35615dd855
# ╠═ec4e41a4-51df-427a-b0e6-7532ec54d215
# ╠═55a3f21e-e452-4c3e-8292-e15aab008fcf
# ╟─669ab918-d3c3-488d-9c4f-d44dc10f9334
# ╠═956e4eff-ddd7-4573-8be3-c87a733e0d87
# ╠═d3eb02fa-d437-464d-ab6d-8e5f7769bf7c
# ╠═ec033ffb-d52d-4d0b-8b4f-fe98fcad2f4b
# ╟─9e7f681d-9551-4ede-870d-059b39c5318a
# ╠═f51b9db7-c304-44eb-8ec7-af5e448bb321
# ╠═533338cc-983f-4852-86e9-da68564f1b9c
# ╠═fa9de589-25d3-4e0e-bbe3-9e48881f8246
# ╠═d6ac7d8e-0628-4822-9a43-ec48cc929c51
# ╠═2abfb3e3-31ff-43e9-bbb0-2d5d37c70f88
# ╠═994b6e87-5f3d-49ae-910b-48fc39ffe4fd
# ╠═14c22e16-6345-40dc-9a5b-984ea548aba6
# ╠═40f43bff-49ad-4a2d-95dd-6583bb1f3175
# ╠═6cf2126b-d9f3-41b1-b14f-138ad2bb6f90
# ╠═8962ba6b-31d4-4059-9a32-57044aafacfb
# ╠═07e5b191-178a-4e81-ab06-8f03fcefc6df
# ╠═505fffb4-1754-43aa-8c26-a42bef3ad89b
# ╠═a2e27b30-3c9f-4c60-b71b-5355e69d8483
# ╠═41317a1e-e475-4d96-86ce-5ae969567bbe
# ╠═793492fc-5cb6-4931-aa5f-05538ddb6bf6
# ╠═d5cdacc5-7c4b-49af-98c8-2629602d61b6
# ╠═aff8bffa-a647-46b2-a771-7b7d236fcf54
# ╠═4b43eeb7-2369-49cd-a825-5d4f0d80a2ce
# ╠═94e4ca9c-eb97-4712-abe1-5c27c99259db
# ╠═de27fdff-6dfb-4ad9-9d71-2706181c1362
# ╠═053081de-8dd0-4986-8620-6b0b2fc59e35
# ╠═9665c806-9afd-42fe-a619-69f7f09abe6e
