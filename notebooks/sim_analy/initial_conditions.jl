### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 82c76c56-e874-4eba-9367-569b656155a2
pwd()

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Input
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
dirname = "/arc7/home/dboyea/sculptor/isolation/1e6/M0.5_c13.1"

# ╔═╡ 920546bd-4838-413c-b687-f891a7f5e985
import TOML

# ╔═╡ 900fda29-c87e-41e4-a46c-b7faea21c0ce
params = TOML.parsefile(joinpath(dirname, "halo.toml"))

# ╔═╡ 79b07d75-fb05-4833-ac2c-ea0e9c24e791
begin 
	R_s = params["r_s"]
	M_s = params["M_s"]
	c = 13.1
end

# ╔═╡ ef3d1d5f-0979-44a7-8f0f-bf4638ea5612
ρ_s(r) = 1/8π /r_s_s^3 * exp(-r / r_s_s)

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $dirname")
	snap = lguys.Snapshot("$dirname/initial.hdf5")

end

# ╔═╡ 97e98ab8-b60b-4b48-b465-a34a16858f88
md"""
# Initial/Final
"""

# ╔═╡ 367095eb-1821-4947-a00e-10b35fbc82d2
halo = lguys.NFW(M_s = M_s, r_s=R_s)

# ╔═╡ e47f633f-2e46-4ffd-82c9-b762980d7b07
lguys.calc_M(halo, R_s)

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"$V_\textrm{circ}$ / km s$^{-1}$")

	rc, Vc = lguys.calc_V_circ(snap)
	lines!(log10.(rc), Vc * lguys.V0, label="initial")


	A(x) = log(1+x)  - x/(1+x)
	V200 = sqrt(M_s * A(c) / (R_s * c))
	println(V200)
	V_nfw(x) = V200 * sqrt(A(x) / x / (A(c) / c))

	log_r = LinRange(-2, 2.5, 1000)
	y = V_nfw.(10 .^ log_r ./ R_s)
	lines!(log_r, y * lguys.V0)

	lines!(log_r, lguys.V0 * lguys.calc_V_circ.(halo, 10 .^ log_r))
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
	Arya.hist2d!(ax, log10.(lguys.calc_r(snap.positions)), lguys.calc_r(snap.velocities) * lguys.V0, bins=100)


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

	bins = LinRange(-10, 10, 100)
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
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 27f8deff-96ae-4d9a-a110-d146ac34965a
begin 
	function plot_ρ_dm(snap)
		plot()
		plot_ρ_dm!(snap)
	end

end

# ╔═╡ 60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	limits=(-3, 3, -12, 0))
	plot_ρ_dm!(snap, label="initial")

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	axislegend(ax)
	fig

end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	skip = 10
	snap_i = snap
	idx = 1:skip:length(snap_i)
	r = sort(lguys.calc_r(snap_i))[idx]
	
	M = idx

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log number of particles",
		limits=((-3, 3), nothing)
	)

	ax.yticks = 0:2:7
	lines!(log10.(r), log10.(M))

	fig
end

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="count", xlabel="log r / kpc", )
	
	hist!(log10.(lguys.calc_r(snap)), label="")

	fig
end

# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╠═82c76c56-e874-4eba-9367-569b656155a2
# ╠═7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═920546bd-4838-413c-b687-f891a7f5e985
# ╠═900fda29-c87e-41e4-a46c-b7faea21c0ce
# ╠═79b07d75-fb05-4833-ac2c-ea0e9c24e791
# ╠═ef3d1d5f-0979-44a7-8f0f-bf4638ea5612
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╟─97e98ab8-b60b-4b48-b465-a34a16858f88
# ╠═367095eb-1821-4947-a00e-10b35fbc82d2
# ╠═e47f633f-2e46-4ffd-82c9-b762980d7b07
# ╟─0e89851e-763f-495b-b677-b664501a17ef
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╟─72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═27f8deff-96ae-4d9a-a110-d146ac34965a
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
