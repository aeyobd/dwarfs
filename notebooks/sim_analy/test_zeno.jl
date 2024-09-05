### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 01165464-35b2-43e9-9d52-d06fa9edf3cf
using LilGuys

# ╔═╡ f979f2a8-3420-4ede-a739-7d727dfdf818
md"""
# Initial conditions
Check the initial DM halo. Loads a snapshot and checks if the halo has the expected density & velocity profile
"""

# ╔═╡ 920546bd-4838-413c-b687-f891a7f5e985
import TOML

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Input
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
model_dir = "/arc7/home/dboyea/dwarfs/zeno/halos"

# ╔═╡ d7a04cc7-369e-4687-b423-deda779f1c57
name = "nfw_1e7"

# ╔═╡ eb17e47b-b650-4362-ba29-77344e37bc48
md"""
# File loading
"""

# ╔═╡ 037a5670-ac79-4b82-9982-bec9c2495f91
Base.@kwdef struct NFWZeno <: lguys.SphericalProfile
	m_a = lguys.A_NFW(1)
	a = 1
	b = 64
	taper = :exp
end

# ╔═╡ 3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
begin 
	halo = NFW(M_s=1, r_s=1)
end

# ╔═╡ d3313f08-7e4e-43b8-b55d-ea099d031bfe
begin 
	using DataFrames, CSV

	zeno_prof = CSV.read("/astro/dboyea/dwarfs/zeno/profiles/nfw.csv", DataFrame)

	zeno_prof.Radius *= halo.r_s
	zeno_prof.Mass *= halo.M_s
	zeno_prof.Density *=  halo.M_s / halo.r_s^3

	zeno_prof[!, :V_circ] = sqrt.(zeno_prof.Mass ./ zeno_prof.Radius)

	zeno_prof
end

# ╔═╡ 2ff32520-d828-4c31-8005-f04ca4941177
halo_zeno = NFWZeno()

# ╔═╡ f0d5a936-dd2e-4a51-be5b-00217389e017
function calc_ρ_nfw(halo, r)
	return halo.m_a / (4π * lguys.A_NFW(1) * r * (halo.a + r)^2)
end

# ╔═╡ dd47be7c-99dc-4f66-8bb4-33f827438a26
function calc_ρ_e(halo, r)
	ρ_b = calc_ρ_nfw(halo, halo.b)

	b = halo.b
	γ = b/(b+halo.a) - 1/2
	return ρ_b * (b/r)^2 * exp(-2γ * (r/b - 1))
end

# ╔═╡ 4cae6065-6c19-4667-b5ae-87ed49958e70
function LilGuys.calc_ρ(halo::NFWZeno, r)
	if r < halo.b
		return calc_ρ_nfw(halo, r)
	else
		return calc_ρ_e(halo, r)
	end
end

# ╔═╡ ac9c4de3-ae76-43da-b89c-5b358c434c7c
calc_ρ(halo_zeno, halo_zeno.b)

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $model_dir")
	snap = lguys.Snapshot(joinpath(model_dir, "$name.hdf5"))

	if false
		snap.x_cen = lguys.centroid(snap.positions)
		snap.v_cen = lguys.centroid(snap.velocities)
		
		println(snap.x_cen)
		println(snap.v_cen)
		
		snap.positions .-= snap.x_cen
		snap.velocities .-= snap.v_cen
		
		snap.x_cen = zeros(3)
		snap.v_cen = zeros(3)
	end
	
	lguys.calc_r(snap)
end

# ╔═╡ 0cff7839-d23e-4ba6-bb44-9efc3aca5f65
maximum(calc_r(snap))

# ╔═╡ 0ccb9018-d88c-4cec-a8da-625be1289bfe
snap.x_cen

# ╔═╡ 5ebe92b8-602e-42be-8751-58898b7323b0
snap.v_cen

# ╔═╡ 14e3b593-17b9-4acd-a9cb-d5923662a02c
prof = lguys.calc_profile(snap, filt_bound=false)

# ╔═╡ 9fb58f1b-c98b-4a93-9683-ab478e44e2d7
prof.v_circ_max * V2KMS

# ╔═╡ 2e293959-9c05-4d9b-b889-a68584ca88f0
prof.r_circ_max 

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"$V_\textrm{circ}$ / km s$^{-1}$")

	log_r = LinRange(-2, 2.5, 1000)

	lines!(prof.log_r, prof.v_circ .* V2KMS, label="snapshot")

	lines!(log_r, lguys.V2KMS * lguys.calc_v_circ.(halo, 10 .^ log_r), label="analytic")

	# lines!(log10.(zeno_prof.Radius), zeno_prof.V_circ .* V2KMS, label="zeno profile")
	scatter!(log10.(lguys.calc_r_circ_max(halo)), lguys.calc_v_circ_max(halo) * lguys.V2KMS, label="vmax, rmax")

	axislegend()
	fig
end

# ╔═╡ a49d1735-203b-47dd-81e1-500ef42b054e
md"""
phase space distribution of star particles initial and final snapshot
"""

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radius / kpc", ylabel="velocity (km/s)" )
	Arya.hist2d!(ax, log10.(lguys.calc_r(snap.positions)), lguys.calc_r(snap.velocities) * lguys.V2KMS, bins=100)


	fig
end

# ╔═╡ a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
md"""
the dark matter distribution of  the snapshot (initial and final)
"""

# ╔═╡ b9746093-0f2f-4478-82ba-00911c8fcceb
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bins = LinRange(-10, 10, 30)
	Arya.hist2d!(ax, snap.positions[1, :], snap.positions[2, :], bins = bins, colorscale=log10)

	fig
end

# ╔═╡ 24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
md"""
stellar distribution of initial and final snapshot
"""

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 100, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	#limits=(-2, 5, -12, 2)
	)

	log_r = LinRange(-2, 3, 1000)
	y = log10.(calc_ρ.(halo, 10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)
	
	y = log10.(calc_ρ.(halo_zeno, 10 .^ log_r))
	lines!(log_r, y, label="analytic", linestyle=:dash)

	lines!(log10.(zeno_prof.Radius), 
		log10.(zeno_prof.Density ), 
		label="zeno", linestyle=:dash)
	
	lines!(prof.log_r, log10.(prof.rho))

	axislegend(ax)
	fig

end

# ╔═╡ a0b9d23b-05ca-4819-a9e6-7734cb8d1981
calc_ρ(halo_zeno, 1e3)

# ╔═╡ aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9
lguys.get_M_tot(halo)

# ╔═╡ 84b3759e-a598-4afc-a2b4-ce841e80ff96
let
	skip = 10
	snap_i = snap
	idx = 1:skip:length(snap_i)
	r = sort(lguys.calc_r(snap_i))[idx]
	
	M = cumsum(snap_i.masses)[idx]

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log M in",
		#limits=((-3, 3), nothing)
	)

	lines!(log10.(r), log10.(M))


	# lines!(log10.(zeno_prof.Radius), 
	# 	log10.(zeno_prof.Mass), 
	# 	label="zeno", linestyle=:dash)

	
	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_M.(halo, 10 .^ log_r) )
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	fig
end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log number of particles",
		#limits=((-3, 3), nothing)
	)

	lines!(prof.log_r_bins[2:end], log10.(cumsum(prof.counts)))

	fig
end

# ╔═╡ f01efc63-c4ac-45ae-8dac-209819a6249e
lguys.calc_M(halo, 2)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="log counts / bin", xlabel="log r / kpc",
		limits=(nothing, (-0.1, nothing))
	)
	
	scatter!(prof.log_r, log10.(prof.counts))

	fig
end

# ╔═╡ fddeb921-468b-4f00-b4cb-a6fc4faec555
R200 = lguys.calc_R200(halo)

# ╔═╡ da55fc43-69f7-4375-b8fc-c61dd606fb24
N200 = sum(lguys.calc_r(snap) .< R200)

# ╔═╡ d841539a-f755-460f-9994-16229aadca6a
grav_softening = 4R200 / sqrt(N200)

# ╔═╡ 5d5d72d6-8272-40c9-bce8-d7d90c670052
md"""
# Softening

From @power2003, we can estimate the ideal softening length with 

``
h = \frac{4R_{200}}{\sqrt{N_{200}}}
``

In our case, 
- R200 = $R200
- N200 = $N200
- so h= $grav_softening kpc
"""

# ╔═╡ ddea90f9-90f5-4624-a292-e7007da4247b
0.959 * 0.14

# ╔═╡ 98a01aea-7d8b-4978-b8cf-2865d6d04e28
0.42977 * 0.14

# ╔═╡ db034d78-f647-4382-b5e1-5e4623350d96
md"""
## Dynamical time
"""

# ╔═╡ 7032a304-1448-4182-b22b-6083a2efea5d
r_circs = 10 .^ LinRange(-3, 0, 1000)

# ╔═╡ ac22ac31-f4bf-497b-9971-c98a9900acfb
v_circs = lguys.calc_v_circ.(halo, r_circs)

# ╔═╡ d112b4d5-6a79-43ae-9ec1-23285e7c4a6e
t_dyn = 2π * r_circs ./ v_circs

# ╔═╡ c0fd9958-2ba9-45d4-87d0-0a401939811b
t_dyn_rho = @. 1 / sqrt(lguys.G * lguys.calc_ρ(halo, r_circs))

# ╔═╡ df08309c-8939-4abb-ac45-684c175c24f0
let
	fig, ax = FigAxis(xlabel="log r", ylabel = "t circ (code units)")

	lines!(log10.(r_circs), t_dyn)
	lines!(log10.(r_circs), t_dyn_rho)
	vlines!(log10(grav_softening))

	fig
end

# ╔═╡ 5798f8f8-32ce-4e9e-8489-4a165f6d240c
 @. 1 / sqrt(lguys.G * lguys.calc_ρ(halo, grav_softening))

# ╔═╡ 3e94b33f-35a9-4ccc-a90e-340b6beb310d
t_max = lguys.calc_r_circ_max(halo) / lguys.calc_v_circ_max(halo)

# ╔═╡ Cell order:
# ╟─f979f2a8-3420-4ede-a739-7d727dfdf818
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═01165464-35b2-43e9-9d52-d06fa9edf3cf
# ╠═920546bd-4838-413c-b687-f891a7f5e985
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═d7a04cc7-369e-4687-b423-deda779f1c57
# ╟─eb17e47b-b650-4362-ba29-77344e37bc48
# ╠═037a5670-ac79-4b82-9982-bec9c2495f91
# ╠═3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
# ╠═2ff32520-d828-4c31-8005-f04ca4941177
# ╠═f0d5a936-dd2e-4a51-be5b-00217389e017
# ╠═dd47be7c-99dc-4f66-8bb4-33f827438a26
# ╠═4cae6065-6c19-4667-b5ae-87ed49958e70
# ╠═ac9c4de3-ae76-43da-b89c-5b358c434c7c
# ╠═0cff7839-d23e-4ba6-bb44-9efc3aca5f65
# ╠═d3313f08-7e4e-43b8-b55d-ea099d031bfe
# ╠═0ccb9018-d88c-4cec-a8da-625be1289bfe
# ╠═5ebe92b8-602e-42be-8751-58898b7323b0
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═14e3b593-17b9-4acd-a9cb-d5923662a02c
# ╠═9fb58f1b-c98b-4a93-9683-ab478e44e2d7
# ╠═2e293959-9c05-4d9b-b889-a68584ca88f0
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╟─72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═a0b9d23b-05ca-4819-a9e6-7734cb8d1981
# ╠═aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9
# ╠═84b3759e-a598-4afc-a2b4-ce841e80ff96
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═f01efc63-c4ac-45ae-8dac-209819a6249e
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╟─5d5d72d6-8272-40c9-bce8-d7d90c670052
# ╠═fddeb921-468b-4f00-b4cb-a6fc4faec555
# ╠═da55fc43-69f7-4375-b8fc-c61dd606fb24
# ╠═d841539a-f755-460f-9994-16229aadca6a
# ╠═ddea90f9-90f5-4624-a292-e7007da4247b
# ╠═98a01aea-7d8b-4978-b8cf-2865d6d04e28
# ╟─db034d78-f647-4382-b5e1-5e4623350d96
# ╠═7032a304-1448-4182-b22b-6083a2efea5d
# ╠═ac22ac31-f4bf-497b-9971-c98a9900acfb
# ╠═d112b4d5-6a79-43ae-9ec1-23285e7c4a6e
# ╠═df08309c-8939-4abb-ac45-684c175c24f0
# ╠═c0fd9958-2ba9-45d4-87d0-0a401939811b
# ╠═5798f8f8-32ce-4e9e-8489-4a165f6d240c
# ╠═3e94b33f-35a9-4ccc-a90e-340b6beb310d
