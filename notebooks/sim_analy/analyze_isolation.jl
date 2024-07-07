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

# ╔═╡ 374489bc-627f-4fc9-9734-7c49456710ac
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 96c91860-f3cc-4531-a8cf-39c85887b394
import TOML

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 82c76c56-e874-4eba-9367-569b656155a2


# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Inputs
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
name = "/arc7/home/dboyea/sculptor/isolation/1e6_M0.8_c13/"

# ╔═╡ a29c993a-c7eb-4b57-a474-50bdbd0ce1ec
halo_params = TOML.parsefile(joinpath(name, "halo.toml"))

# ╔═╡ 79b07d75-fb05-4833-ac2c-ea0e9c24e791
halo = lguys.NFW(; lguys.dict_to_tuple(halo_params)...)

# ╔═╡ 46ce993c-50a1-43a5-8d55-132fac90de33
lguys.calc_M200(halo)

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $name")
	out = lguys.Output(joinpath(name, "out/combined.hdf5"))

	cens = CSV.read(joinpath(name, "out/centres.csv"), DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ 4710d1b5-fd6f-4369-8398-64e9123c8b41
idx_i = 10; idx_f = length(out)

# ╔═╡ 97f89831-00e6-49a2-a712-ac47fd2dee47
out.times[idx_f] * lguys.T0

# ╔═╡ 7a43ee65-0fce-4e4b-9226-714f4cb0e106
out.times[idx_i] * lguys.T0

# ╔═╡ 97e98ab8-b60b-4b48-b465-a34a16858f88
md"""
# Initial/Final
"""

# ╔═╡ c5672da9-0dad-4d22-abe5-9e186ccde02d
begin
	snap_i = out[idx_i]
	snap_i.x_cen = x_cen[:, idx_i]
	snap_i.v_cen = v_cen[:, idx_i]

	snap_f = out[idx_f]
	snap_f.x_cen = x_cen[:, idx_f]
	snap_f.v_cen = v_cen[:, idx_f]

end

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"V_\textrm{circ}")

	rc, Vc = lguys.calc_V_circ(snap_i)
	lines!(log10.(rc), Vc * lguys.V0, label="initial")

	rc, Vc = lguys.calc_V_circ(snap_f)
	lines!(log10.(rc), Vc * lguys.V0, label="final")

	V_nfw(x) = lguys.calc_V_circ(halo, x)

	log_r = LinRange(-2, 2.5, 1000)
	y = V_nfw.(10 .^ log_r)
	lines!(log_r, y * lguys.V0)
	fig
end

# ╔═╡ a49d1735-203b-47dd-81e1-500ef42b054e
md"""
phase space distribution of star particles initial and final snapshot
"""

# ╔═╡ 7d717638-1caf-4267-b9f5-c060c19e2849


# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1,1], xlabel="log radius / kpc", ylabel="velocity (km/s)" )

	hist_kwargs = (; bins=100, colorscale=log10, colorrange=(1, nothing))
	
	Arya.hist2d!(ax, log10.(lguys.calc_r(snap_i.positions)), lguys.calc_r(snap_i.velocities) * lguys.V0; hist_kwargs...)

	x = LinRange(0, 2, 1000)
	e = lguys.calc_Φ.(halo, 10 .^ x)
	y = sqrt.(-2e) .* lguys.V0
	lines!(x, y)

	y = lguys.calc_V_circ.(halo, 10 .^ x) .* lguys.V0
	lines!(x, y)

	ax2 = Axis(fig[1,2] )
	Arya.hist2d!(ax2, log10.(lguys.calc_r(snap_f.positions)), lguys.calc_r(snap_f.velocities) * lguys.V0;
	hist_kwargs...)

	linkaxes!(ax, ax2)
	hideydecorations!(ax2)

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
	hist_kwargs = (; bins=bins, colorscale=log10, colorrange=(1, nothing))
	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :]; hist_kwargs...)

	ax2 = Axis(fig[1,2], aspect=1,
	title="final")
	
	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :]; hist_kwargs...)
	hideydecorations!(ax2)
	linkaxes!(ax, ax2)
	fig
end

# ╔═╡ 82a8514e-c1de-4446-918b-9156734c213e
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "vx / kpc", ylabel="vy/kpc", title="initial")

	bins = LinRange(-60, 60, 100)
	hist_kwargs = (; bins=bins, colorscale=log10, colorrange=(1, nothing))
	Arya.hist2d!(ax, snap_i.velocities[1, :] .* lguys.V0, snap_i.velocities[2, :] .* lguys.V0; hist_kwargs...)

	ax2 = Axis(fig[1,2], aspect=1,
	title="final")
	
	Arya.hist2d!(ax2, snap_f.velocities[1, :] .* lguys.V0, snap_f.velocities[2, :] .* lguys.V0; hist_kwargs...)
	hideydecorations!(ax2)
	linkaxes!(ax, ax2)
	fig
end

# ╔═╡ 24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
md"""
stellar distribution of initial and final snapshot
"""

# ╔═╡ 106cbda4-57e0-459b-868b-b44339c944fc
begin 
	ps = lguys.extract_vector(snap_f, :positions)
	radii = lguys.calc_r(ps)

end

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
	plot_ρ_dm!(snap_i, label="initial")
	plot_ρ_dm!(snap_f, label="final")


	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_ρ.(halo, 10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	axislegend(ax)
	fig

end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	skip = 10
	snap = snap_i
	idx = 1:skip:length(snap)
	r = sort(lguys.calc_r(snap))[idx]
	
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

# ╔═╡ e33d56a7-7a0e-4fa9-8f0d-041b43584d59
sum(lguys.calc_r(snap_i) .< 0.1)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="count", xlabel="log r / kpc", )
	
	hist!(log10.(lguys.calc_r(snap_i)), label="")

	fig
end

# ╔═╡ 1ce0aa3c-953e-4f7b-bf7e-7c815c505b5e
println(sum(lguys.calc_r(snap_i) .< 0.05))

# ╔═╡ fb0dec74-aaab-43a4-9b37-d13634c5dcad
md"""
# Time variation
"""

# ╔═╡ 5b9a64a1-38c1-4f4d-aac3-d91663c368a3
lguys.calc_L_tot(snap_i)

# ╔═╡ ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
let 
	idx = 1:10:length(out)

	Ls = hcat([lguys.calc_L_tot(snap) for snap in out[idx]]...)

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="total angular momentum component"
	)

	t = out.times[idx] * lguys.T0
	
	for i in 1:3
		lines!(t, Ls[i, :], label=["x", "y", "z"][i])
	end

	axislegend(ax)
	fig
end

# ╔═╡ bfa7593c-4915-4e03-83e6-8f790de4c1a5
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="relative change in energy")
	idx = 1:10:length(out)
	Es = [lguys.calc_E_tot(snap) for snap in out[idx]]
	lines!((Es ./ Es[1]))
	fig
end

# ╔═╡ e3a45a8e-cc52-4e9d-9db3-97109b59fc77
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="-V / 2T (virial ratio approx 1)")

	idx = 1:10:length(out)
	
	E_kin = [sum(0.5 * lguys.calc_v(snap) .^ 2) for snap in out[idx]]
	E_pot = [0.5 * sum(snap.Φs) for snap in out[idx]]

	lines!(-E_pot ./ 2E_kin)
	fig
end

# ╔═╡ 91a44ed4-8466-4a58-b3ff-1e7630b8ac8c
lguys.plot_xyz(x_cen)

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
lguys.plot_xyz(v_cen * lguys.V0)

# ╔═╡ dc221349-eb61-4ace-8de3-a6c50249aca0
function find_radii_fracs(out, x_cen) 
	rs = Vector[]
	Ms = Vector[]
	rs_s = Vector[]

	percens = [0.003, 0.01, .03, .1, .3, .9]
	
	for i in 1:length(out)
		r = lguys.calc_r(out[i].positions .- x_cen[:, i])
		s_idx = sortperm(r)
		push!(rs, r[s_idx])
			end

	rs = hcat(rs...)
	rs_s = hcat(rs_s...)

	return percens, rs, rs_s

end

# ╔═╡ 34244a2e-9501-451c-bd77-bebfebde2a78
percens, rs, rs_s = find_radii_fracs(out, x_cen)

# ╔═╡ 967136d3-8d58-4fdc-9537-aa3a85a92528
times = out.times * lguys.T0

# ╔═╡ f21cfe22-95f3-485d-902b-b022a41548c2
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel="log r containing fraction of DM")
	
	for i in eachindex(percens)
		Npoints = round.(Int, length(out[1]) * percens[i])
		label = "$(round(Npoints/length(out[1]), digits=3))"
		y = log10.(rs[Npoints, :])
		lines!(times, y, 				
			label=label, color=i, colorrange=(1, length(percens))
		)
	end

	Legend(fig[1,2], ax, "f = ")
	fig
end

# ╔═╡ 3b2bb553-0130-4c8a-80ad-6e1f7071a293
lguys.plot_xyz(lguys.extract_vector(out, :positions, 100_000))

# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═374489bc-627f-4fc9-9734-7c49456710ac
# ╠═96c91860-f3cc-4531-a8cf-39c85887b394
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╠═82c76c56-e874-4eba-9367-569b656155a2
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═a29c993a-c7eb-4b57-a474-50bdbd0ce1ec
# ╠═4710d1b5-fd6f-4369-8398-64e9123c8b41
# ╠═79b07d75-fb05-4833-ac2c-ea0e9c24e791
# ╠═46ce993c-50a1-43a5-8d55-132fac90de33
# ╠═97f89831-00e6-49a2-a712-ac47fd2dee47
# ╠═7a43ee65-0fce-4e4b-9226-714f4cb0e106
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╟─97e98ab8-b60b-4b48-b465-a34a16858f88
# ╠═c5672da9-0dad-4d22-abe5-9e186ccde02d
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═a49d1735-203b-47dd-81e1-500ef42b054e
# ╠═7d717638-1caf-4267-b9f5-c060c19e2849
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╠═82a8514e-c1de-4446-918b-9156734c213e
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═106cbda4-57e0-459b-868b-b44339c944fc
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═27f8deff-96ae-4d9a-a110-d146ac34965a
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═e33d56a7-7a0e-4fa9-8f0d-041b43584d59
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╠═1ce0aa3c-953e-4f7b-bf7e-7c815c505b5e
# ╟─fb0dec74-aaab-43a4-9b37-d13634c5dcad
# ╠═5b9a64a1-38c1-4f4d-aac3-d91663c368a3
# ╠═ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═e3a45a8e-cc52-4e9d-9db3-97109b59fc77
# ╠═91a44ed4-8466-4a58-b3ff-1e7630b8ac8c
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═dc221349-eb61-4ace-8de3-a6c50249aca0
# ╠═34244a2e-9501-451c-bd77-bebfebde2a78
# ╠═f21cfe22-95f3-485d-902b-b022a41548c2
# ╠═967136d3-8d58-4fdc-9537-aa3a85a92528
# ╠═3b2bb553-0130-4c8a-80ad-6e1f7071a293
