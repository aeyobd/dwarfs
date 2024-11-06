### A Pluto.jl notebook ###
# v0.20.0

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
model_dir = ENV["DWARFS_ROOT"] * "/simulations/sculptor/1e6_V31_r3.2/vasiliev+21_mean/"

# ╔═╡ d7a04cc7-369e-4687-b423-deda779f1c57
name = "initial"

# ╔═╡ eb17e47b-b650-4362-ba29-77344e37bc48
md"""
# File loading
"""

# ╔═╡ 3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
halo = lguys.load_profile(joinpath(model_dir, "../halo.toml"))

# ╔═╡ 7f9db45f-38ea-4427-9af1-d5431429f612
halo.r_s

# ╔═╡ 54a2a708-d8ba-4c5c-9e67-ac656dd8e9f4
lguys.calc_M(halo, 1)

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $model_dir")
	snap = lguys.Snapshot(joinpath(model_dir, "$name.hdf5"))
	snap.Φs = - ones(length(snap))

	sss = lguys.Centres.SS_State(snap)
	lguys.calc_centre!(sss, snap)

	snap.x_cen = sss.centre.position
	snap.v_cen = sss.centre.velocity

	println(snap.x_cen)
	println(snap.v_cen)

	lguys.calc_r(snap)
end

# ╔═╡ 0ccb9018-d88c-4cec-a8da-625be1289bfe
snap.x_cen

# ╔═╡ 5ebe92b8-602e-42be-8751-58898b7323b0
snap.v_cen * V2KMS

# ╔═╡ 14e3b593-17b9-4acd-a9cb-d5923662a02c
prof = lguys.MassProfile3D(snap)

# ╔═╡ 9fb58f1b-c98b-4a93-9683-ab478e44e2d7
prof.v_circ_max * V2KMS

# ╔═╡ 2e293959-9c05-4d9b-b889-a68584ca88f0
prof.r_circ_max 

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"$V_\textrm{circ}$ / km s$^{-1}$")

	
	log_r = LinRange(-2, 4, 1000)

	lines!(log10.(prof.r_circ), prof.v_circ .* V2KMS, label="snapshot")

	lines!(log_r, lguys.V2KMS * lguys.calc_v_circ.(halo, 10 .^ log_r), label="analytic")

	scatter!(log10.(prof.r_circ_max), prof.v_circ_max * lguys.V2KMS, label="max; observed")
	
	scatter!(log10.(lguys.calc_r_circ_max(halo)), lguys.calc_v_circ_max(halo) * lguys.V2KMS, label="max; halo")
	
	axislegend()

	
	ax2 = Axis(fig[2,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"d v / v")
	ax2.limits=(nothing, nothing, -0.1, 0.1)

	y = prof.v_circ
	x = log10.(prof.r_circ)
	ye = calc_v_circ.(halo, 10 .^ x)

	res = (y .- ye) ./ ye

	err = prof.v_circ_err ./ ye

	errscatter!(x, res, yerr=err)

	rowsize!(fig.layout, 2, Relative(0.3))
	
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
	Arya.hist2d!(ax, log10.(lguys.calc_r(snap)), lguys.calc_v(snap) * lguys.V2KMS, bins=100)


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
	xlabel = "dx / kpc", ylabel="dy/kpc", title="initial")

	bins = LinRange(-10, 10, 30)
	Arya.hist2d!(ax, snap.positions[1, :] .- snap.x_cen[1], snap.positions[2, :] .- snap.x_cen[2], bins = bins, colorscale=log10)

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

# ╔═╡ 09139e38-fdb7-4754-9cc0-79e13a131b08
12*3600

# ╔═╡ 60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	#limits=(-2, 5, -12, 2)
	)

	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_ρ.(halo, 10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)

	
	lines!(prof.log_r, log10.(prof.rho))

	axislegend(ax)
	fig

end

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
		limits=(nothing, (-0.1, log10(length(snap)/10)))
	)
	
	scatter!(prof.log_r, log10.(prof.counts))

	fig
end

# ╔═╡ fddeb921-468b-4f00-b4cb-a6fc4faec555
R200 = lguys.calc_R200(halo)

# ╔═╡ da55fc43-69f7-4375-b8fc-c61dd606fb24
N200 = sum(lguys.calc_r(snap) .< R200)

# ╔═╡ d841539a-f755-460f-9994-16229aadca6a
grav_softening = 4R200 / sqrt(N200) / sqrt(10)

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
	fig, ax = FigAxis(xlabel="log r", ylabel = L"$t_\textrm{dyn}$ (code units)")

	lines!(log10.(r_circs), log10.(t_dyn), label="circ")
	lines!(log10.(r_circs), log10.(t_dyn_rho), label="grav")
	vlines!(log10(grav_softening), linestyle=:dot, color=:black, label="softening")

	axislegend()

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
# ╠═3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
# ╠═7f9db45f-38ea-4427-9af1-d5431429f612
# ╠═0ccb9018-d88c-4cec-a8da-625be1289bfe
# ╠═5ebe92b8-602e-42be-8751-58898b7323b0
# ╠═54a2a708-d8ba-4c5c-9e67-ac656dd8e9f4
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═14e3b593-17b9-4acd-a9cb-d5923662a02c
# ╠═9fb58f1b-c98b-4a93-9683-ab478e44e2d7
# ╠═2e293959-9c05-4d9b-b889-a68584ca88f0
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═09139e38-fdb7-4754-9cc0-79e13a131b08
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9
# ╠═84b3759e-a598-4afc-a2b4-ce841e80ff96
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═f01efc63-c4ac-45ae-8dac-209819a6249e
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╟─5d5d72d6-8272-40c9-bce8-d7d90c670052
# ╠═fddeb921-468b-4f00-b4cb-a6fc4faec555
# ╠═da55fc43-69f7-4375-b8fc-c61dd606fb24
# ╠═d841539a-f755-460f-9994-16229aadca6a
# ╟─db034d78-f647-4382-b5e1-5e4623350d96
# ╠═7032a304-1448-4182-b22b-6083a2efea5d
# ╠═ac22ac31-f4bf-497b-9971-c98a9900acfb
# ╠═d112b4d5-6a79-43ae-9ec1-23285e7c4a6e
# ╠═df08309c-8939-4abb-ac45-684c175c24f0
# ╠═c0fd9958-2ba9-45d4-87d0-0a401939811b
# ╠═5798f8f8-32ce-4e9e-8489-4a165f6d240c
# ╠═3e94b33f-35a9-4ccc-a90e-340b6beb310d
