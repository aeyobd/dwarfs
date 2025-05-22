### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 01165464-35b2-43e9-9d52-d06fa9edf3cf
using LilGuys

# ╔═╡ 9303ace1-bd7c-4391-bace-8a0a8eccd251
using PlutoUI

# ╔═╡ f979f2a8-3420-4ede-a739-7d727dfdf818
md"""
# Initial conditions
Check the initial DM halo. Loads a snapshot and checks if the halo has the expected density & velocity profile
"""

# ╔═╡ 920546bd-4838-413c-b687-f891a7f5e985
import TOML

# ╔═╡ d8998dd8-bac8-450d-9e9f-fa5d0b282e13
CairoMakie.activate!(type=:png)

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Input
"""

# ╔═╡ 80da94f9-6fdd-4591-93e3-c98ea1479c65
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 21080ba9-d6df-4ba3-b4f9-e97f8d850dc6
@bind inputs confirm(notebook_inputs(;
	modelname = TextField(70, default="ursa_minor/1e6_v31_r4.0/orbitname"),
))

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
model_dir = joinpath(ENV["DWARFS_ROOT"], "simulations", inputs.modelname)

# ╔═╡ d7a04cc7-369e-4687-b423-deda779f1c57
#name = "initial"
name = "initial"

# ╔═╡ 8b79ec3a-73d7-4dd6-8c91-8d3358f7896e
paramname = "../halo.toml"

# ╔═╡ eb17e47b-b650-4362-ba29-77344e37bc48
md"""
# File loading
"""

# ╔═╡ 3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
halo_in = lguys.load_profile(joinpath(model_dir,paramname ))

# ╔═╡ 7430e5ef-7cf6-4d0e-a1e8-3db3c35388e8
halo = NFW(r_s=halo_in.r_s, M_s=halo_in.M_s)

# ╔═╡ 7f9db45f-38ea-4427-9af1-d5431429f612
halo.r_s

# ╔═╡ 54a2a708-d8ba-4c5c-9e67-ac656dd8e9f4
lguys.mass(halo, (1))

# ╔═╡ 020bbb15-235d-468c-bc02-01f3b75708e9
LilGuys.expint(Complex(-0.01))

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $model_dir")
	snap = lguys.Snapshot(joinpath(model_dir, "$name.hdf5"))
	snap.potential = -ones(size(snap.positions, 2)) # all particles bound
	sss = lguys.Centres.SS_State(snap)
	lguys.calc_centre!(sss, snap)

	snap.x_cen = sss.centre.position
	snap.v_cen = sss.centre.velocity

	println(snap.x_cen)
	println(snap.v_cen)

	lguys.radii(snap)
end

# ╔═╡ 0ccb9018-d88c-4cec-a8da-625be1289bfe
snap.x_cen

# ╔═╡ 5ebe92b8-602e-42be-8751-58898b7323b0
snap.v_cen * V2KMS

# ╔═╡ c9cd2c2f-90b7-4d1c-8c48-586c1dac0257
snap.v_cen

# ╔═╡ 14e3b593-17b9-4acd-a9cb-d5923662a02c
prof = lguys.MassProfile(snap)

# ╔═╡ 0d9f39fa-acfb-4408-806b-924c7b177c05
density_prof = DensityProfile(snap)

# ╔═╡ ac08974f-8d28-4098-81a8-21884157e76a
props = LilGuys.MassScalars(snap, prof)

# ╔═╡ 9fb58f1b-c98b-4a93-9683-ab478e44e2d7
props.v_circ_max

# ╔═╡ 3d4f1a39-fbfc-4e07-8bd4-17e7813da5af
props.v_circ_max * V2KMS

# ╔═╡ 2e293959-9c05-4d9b-b889-a68584ca88f0
props.r_circ_max 

# ╔═╡ 3cc4a870-da17-4c0a-8ab8-27e134c71044
0.044 * 0.53612788

# ╔═╡ bc127b27-c573-462b-84f1-66890391b49b


# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure(
		size=(4*72, 4*72)
	)
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"$V_\textrm{circ}$ / km s$^{-1}$")

	
	log_r = LinRange(-2, 2, 1000)

	lines!(log10.(prof.radii), lguys.circular_velocity(prof) .* V2KMS, label="snapshot")

	lines!(log_r, lguys.V2KMS * lguys.v_circ.(halo, 10 .^ log_r), label="analytic")

	scatter!(log10.(props.r_circ_max), props.v_circ_max * lguys.V2KMS, label="max; observed")
	
	scatter!(log10.(lguys.r_circ_max(halo)), lguys.v_circ_max(halo) * lguys.V2KMS, label="max; halo")

	hidexdecorations!(ticks=false, minorticks=false)
	axislegend()

	
	ax2 = Axis(fig[2,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"d v / v")
	ax2.limits=(nothing, nothing, -0.2, 0.2)

	y = lguys.circular_velocity(prof)
	x = log10.(prof.radii)
	ye = v_circ.(halo, 10 .^ x)

	res = (y .- ye) ./ ye

	err = lguys.sym_error.(y) ./ ye

	errorscatter!(x, res, yerror=err)

	rowsize!(fig.layout, 2, Relative(0.3))

	linkxaxes!(ax, ax2)
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
	Arya.hist2d!(ax, log10.(lguys.radii(snap)), lguys.speeds(snap) * lguys.V2KMS, bins=100)


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
	Arya.hist2d!(ax, snap.positions[1, :] .- snap.x_cen[1], snap.positions[2, :] .- snap.x_cen[2], bins = bins, colorscale=log10, colorrange=(1, nothing))

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
	limits=(-2, 2, -12, 2)
	)

	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.density.(halo, 10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)

	
	lines!(density_prof.log_r, log10.(density_prof.rho))

	axislegend(ax)
	fig

end

# ╔═╡ 2b7c7517-4b28-4db0-860b-dc6e5e514e6c
LilGuys.mean(radii(snap) .> 2)

# ╔═╡ aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9
lguys.mass(halo)

# ╔═╡ 84b3759e-a598-4afc-a2b4-ce841e80ff96
let
	skip = 10
	snap_i = snap
	idx = 1:skip:length(snap_i)
	r = sort(lguys.radii(snap_i))[idx]
	
	M = cumsum(snap_i.masses)[idx]

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log M in",
		#limits=((-3, 3), nothing)
	)

	lines!(log10.(r), log10.(M))


	
	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.mass.(halo, 10 .^ log_r) )
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

	lines!(density_prof.log_r_bins[2:end], log10.(cumsum(density_prof.counts)))

	fig
end

# ╔═╡ f01efc63-c4ac-45ae-8dac-209819a6249e
lguys.mass(halo, 1)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="log counts / bin", xlabel="log r / kpc",
		limits=(nothing, (-0.1, log10(length(snap)/10)))
	)
	
	scatter!(density_prof.log_r, log10.(density_prof.counts))

	fig
end

# ╔═╡ fddeb921-468b-4f00-b4cb-a6fc4faec555
if halo isa lguys.ExpCusp
	R200 = 3halo.r_s
else
	R200 = lguys.R200(halo)
end

# ╔═╡ 233c5aca-4966-4ba5-b5ac-f5d0e0a727dc
EXTRA_SOFTENING = 1/sqrt(10)

# ╔═╡ da55fc43-69f7-4375-b8fc-c61dd606fb24
N200 = sum(lguys.radii(snap) .< R200)

# ╔═╡ d841539a-f755-460f-9994-16229aadca6a
grav_softening = 4R200 / sqrt(N200) * EXTRA_SOFTENING

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

# ╔═╡ 3fc4f98e-1652-49af-ad00-ec9b634b0715
4R200 / sqrt(length(snap)) * EXTRA_SOFTENING

# ╔═╡ b1a3577b-45dd-4a10-889d-1c05c2465433
halo.r_s

# ╔═╡ db034d78-f647-4382-b5e1-5e4623350d96
md"""
## Dynamical time
"""

# ╔═╡ 7032a304-1448-4182-b22b-6083a2efea5d
r_circs = 10 .^ LinRange(-3, log10(R200), 1000)

# ╔═╡ ac22ac31-f4bf-497b-9971-c98a9900acfb
v_circs = lguys.v_circ.(halo, r_circs)

# ╔═╡ d112b4d5-6a79-43ae-9ec1-23285e7c4a6e
t_dyn = 2π * r_circs ./ v_circs

# ╔═╡ c0fd9958-2ba9-45d4-87d0-0a401939811b
t_dyn_rho = @. 1 / sqrt(lguys.G * lguys.density(halo, r_circs))

# ╔═╡ df08309c-8939-4abb-ac45-684c175c24f0
let
	fig, ax = FigAxis(xlabel="log r", ylabel = L"$t_\textrm{dyn}$ (code units)")

	lines!(log10.(r_circs), log10.(t_dyn), label="circ")
	lines!(log10.(r_circs), log10.(t_dyn_rho), label="grav")
	vlines!(log10(grav_softening), linestyle=:dot, color=:black, label="softening")

	axislegend(position=:lt)

	fig
end

# ╔═╡ 5798f8f8-32ce-4e9e-8489-4a165f6d240c
 @. 1 / sqrt(lguys.G * lguys.density(halo, grav_softening))

# ╔═╡ 3e94b33f-35a9-4ccc-a90e-340b6beb310d
t_max = lguys.r_circ_max(halo) / lguys.v_circ_max(halo)

# ╔═╡ Cell order:
# ╟─f979f2a8-3420-4ede-a739-7d727dfdf818
# ╠═21080ba9-d6df-4ba3-b4f9-e97f8d850dc6
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═01165464-35b2-43e9-9d52-d06fa9edf3cf
# ╠═920546bd-4838-413c-b687-f891a7f5e985
# ╠═9303ace1-bd7c-4391-bace-8a0a8eccd251
# ╠═d8998dd8-bac8-450d-9e9f-fa5d0b282e13
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═80da94f9-6fdd-4591-93e3-c98ea1479c65
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═d7a04cc7-369e-4687-b423-deda779f1c57
# ╠═8b79ec3a-73d7-4dd6-8c91-8d3358f7896e
# ╟─eb17e47b-b650-4362-ba29-77344e37bc48
# ╠═3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
# ╠═7430e5ef-7cf6-4d0e-a1e8-3db3c35388e8
# ╠═7f9db45f-38ea-4427-9af1-d5431429f612
# ╠═0ccb9018-d88c-4cec-a8da-625be1289bfe
# ╠═5ebe92b8-602e-42be-8751-58898b7323b0
# ╠═c9cd2c2f-90b7-4d1c-8c48-586c1dac0257
# ╠═54a2a708-d8ba-4c5c-9e67-ac656dd8e9f4
# ╠═020bbb15-235d-468c-bc02-01f3b75708e9
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═14e3b593-17b9-4acd-a9cb-d5923662a02c
# ╠═0d9f39fa-acfb-4408-806b-924c7b177c05
# ╠═ac08974f-8d28-4098-81a8-21884157e76a
# ╠═9fb58f1b-c98b-4a93-9683-ab478e44e2d7
# ╠═3d4f1a39-fbfc-4e07-8bd4-17e7813da5af
# ╠═2e293959-9c05-4d9b-b889-a68584ca88f0
# ╠═3cc4a870-da17-4c0a-8ab8-27e134c71044
# ╠═bc127b27-c573-462b-84f1-66890391b49b
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═09139e38-fdb7-4754-9cc0-79e13a131b08
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═2b7c7517-4b28-4db0-860b-dc6e5e514e6c
# ╠═aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9
# ╠═84b3759e-a598-4afc-a2b4-ce841e80ff96
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═f01efc63-c4ac-45ae-8dac-209819a6249e
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╟─5d5d72d6-8272-40c9-bce8-d7d90c670052
# ╠═fddeb921-468b-4f00-b4cb-a6fc4faec555
# ╠═233c5aca-4966-4ba5-b5ac-f5d0e0a727dc
# ╠═da55fc43-69f7-4375-b8fc-c61dd606fb24
# ╠═d841539a-f755-460f-9994-16229aadca6a
# ╠═3fc4f98e-1652-49af-ad00-ec9b634b0715
# ╠═b1a3577b-45dd-4a10-889d-1c05c2465433
# ╟─db034d78-f647-4382-b5e1-5e4623350d96
# ╠═7032a304-1448-4182-b22b-6083a2efea5d
# ╠═ac22ac31-f4bf-497b-9971-c98a9900acfb
# ╠═d112b4d5-6a79-43ae-9ec1-23285e7c4a6e
# ╠═df08309c-8939-4abb-ac45-684c175c24f0
# ╠═c0fd9958-2ba9-45d4-87d0-0a401939811b
# ╠═5798f8f8-32ce-4e9e-8489-4a165f6d240c
# ╠═3e94b33f-35a9-4ccc-a90e-340b6beb310d
