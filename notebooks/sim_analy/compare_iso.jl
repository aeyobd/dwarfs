### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 84db4ab2-13d5-11ef-2687-b14402ce4daf
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ ee9629d6-8310-4aca-bef7-56df6e633223
using DataFrames, CSV

# ╔═╡ 78f3d3dd-c878-4a47-9cdf-43af0478b010
function plot_ρ_dm!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 543adda2-a9ac-407a-bb3b-62b3fdcf2d44
begin 
	R_s = 2.76
	M_s = 0.290
end

# ╔═╡ 598d69fe-81cc-42c7-9b6e-97a4c4b3224c
begin
	dirname = "/cosma/home/durham/dc-boye1/sculptor/isolation/"
	snapshots = Dict{String, Any}()
	snapshots["1e4"] = lguys.Snapshot(dirname * "1e4/initial.hdf5")
	snapshots["1e5"] = lguys.Snapshot(dirname * "1e5/initial.hdf5")
	snapshots["1e6"] = lguys.Snapshot(dirname * "1e6/initial.hdf5")
	snapshots["1e7"] = lguys.Snapshot(dirname * "1e7/initial.hdf5")
end

# ╔═╡ 26b85a2c-c57a-41dd-b295-972619ee9033
let
	fig = Figure(size=(600,600))
	ax_rho = Axis(fig[1,1],
		#xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel=L"\log\;\rho_\textrm{DM}\quad [10^{10}\textrm{M_\odot/kpc}]"
	)
	for (name, snap) in snapshots
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)


	axislegend(ax_rho)

	ax_n = Axis(fig[2,1],
		xlabel="log r / kpc",
		ylabel="log number of particles",
		yticks=0:7
	)

	for (name, snap) in snapshots
		if name == "1e7"
			skip = 10
		else
			skip = 1
		end

		r = lguys.calc_r(snap)
		sort!(r)
		N_r = 1:length(r)

		lines!(log10.(r), log10.(N_r))
	end


	linkxaxes!(ax_rho, ax_n)
	
	fig
end

# ╔═╡ 8f9e7dad-144e-40e8-969b-86a672eef923
expected = CSV.read("../zeno/profiles/nfw.csv", DataFrame)

# ╔═╡ 06ba05a5-6b8b-476e-9224-04de8ab21c6b
prof_default = CSV.read("../zeno/profiles/default.csv", DataFrame)

# ╔═╡ 13b18244-e94d-47a4-8ded-7a5c4959efd3
let

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/models/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["default"] = lguys.Snapshot(dirname_z * "default_1e6.hdf5")
	zeno_snaps["large range"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	
	R_s = 1
	M_s = 4
	
	fig = Figure()
	ax_rho = Axis(fig[1,1],
		xlabel=L"\log\;r/r_s",
		ylabel=L"\log\;\rho_\textrm{DM}",
		limits=(nothing, (-15, 5))
	)

	scatter!(log10.(prof_default.Radius), log10.(prof_default.Density) .- 1, 
		alpha=0.2, label="default profile")
	
	scatter!(log10.(expected.Radius), log10.(expected.Density) .- 0, 
		alpha=0.2, label="large range profile")
	for (name, snap) in zeno_snaps
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)

	axislegend(ax_rho)
	fig
end

# ╔═╡ f9209efb-3bb6-40e5-b1ec-b78565410ec5
let

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/models/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["1e4"] = lguys.Snapshot(dirname_z * "nfw_1e4.hdf5")
	zeno_snaps["1e5"] = lguys.Snapshot(dirname_z * "nfw_1e5.hdf5")
	zeno_snaps["1e6"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	zeno_snaps["1e7"] = lguys.Snapshot(dirname_z * "nfw_1e7.hdf5")
	
	R_s = 1
	M_s = 4
	
	fig = Figure()
	ax_rho = Axis(fig[1,1],
		xlabel=L"\log\;r/r_s",
		ylabel=L"\log\;\rho_\textrm{DM}",
		limits=(nothing, (-15, 5))
	)
	scatter!(log10.(expected.Radius), log10.(expected.Density), 
		alpha=0.1, label="zeno profile")

	for (name, snap) in zeno_snaps
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)

	axislegend(ax_rho)
	fig
end

# ╔═╡ c51a8c83-0317-49e1-a08c-0ce4859bfdfa
let

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["default"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	zeno_snaps["10"] = lguys.Snapshot(dirname_z * "nfw_1e6_np10.hdf5")
	zeno_snaps["100"] = lguys.Snapshot(dirname_z * "nfw_1e6_np1e2.hdf5")
	zeno_snaps["1e4"] = lguys.Snapshot(dirname_z * "nfw_1e6_np1e4.hdf5")
	zeno_snaps["big range"] = lguys.Snapshot(dirname_z * "nfw_1e6_br.hdf5")

	
	R_s = 1
	M_s = 4
	
	fig = Figure(size=(700, 700))
	ax_rho = Axis(fig[1,1],
		xlabel="",
		ylabel=L"\log\;\rho_\textrm{DM}"
	)
	for (name, snap) in zeno_snaps
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	axislegend(ax_rho)
	fig
end

# ╔═╡ 4a976529-1ad1-4649-a50a-031a51fa9508
let

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/models/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["default (r=4)"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	zeno_snaps["break r = 2"] = lguys.Snapshot(dirname_z * "nfw_1e6_b2.hdf5")
	zeno_snaps["break r = 64"] = lguys.Snapshot(dirname_z * "nfw_1e6_b64.hdf5")
	zeno_snaps["break r = 10"] = lguys.Snapshot(dirname_z * "nfw_1e6_b10.hdf5")
	
	R_s = 1
	M_s = 4
	
	fig = Figure()
	ax_rho = Axis(fig[1,1],
		xlabel="",
		ylabel=L"\log\;\rho_\textrm{DM}"
	)
	for (name, snap) in zeno_snaps
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	axislegend(ax_rho)
	fig
end

# ╔═╡ 9227d920-be3f-40e5-8821-02e2d93ca97e
let

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/models/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["default"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	zeno_snaps["gauss"] = lguys.Snapshot(dirname_z * "nfw_1e6_gauss.hdf5")
	zeno_snaps["sw"] = lguys.Snapshot(dirname_z * "nfw_1e6_sw.hdf5")
	
	R_s = 1
	M_s = 4
	
	fig = Figure()
	ax_rho = Axis(fig[1,1],
		xlabel="",
		ylabel=L"\log\;\rho_\textrm{DM}"
	)
	for (name, snap) in zeno_snaps
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	axislegend(ax_rho)
	fig
end

# ╔═╡ 5ecf013a-bfaf-41ac-bd76-2c6fa4cbddbc
4096 * 16

# ╔═╡ Cell order:
# ╠═84db4ab2-13d5-11ef-2687-b14402ce4daf
# ╠═78f3d3dd-c878-4a47-9cdf-43af0478b010
# ╠═543adda2-a9ac-407a-bb3b-62b3fdcf2d44
# ╠═598d69fe-81cc-42c7-9b6e-97a4c4b3224c
# ╠═ee9629d6-8310-4aca-bef7-56df6e633223
# ╠═26b85a2c-c57a-41dd-b295-972619ee9033
# ╠═8f9e7dad-144e-40e8-969b-86a672eef923
# ╠═06ba05a5-6b8b-476e-9224-04de8ab21c6b
# ╠═13b18244-e94d-47a4-8ded-7a5c4959efd3
# ╠═f9209efb-3bb6-40e5-b1ec-b78565410ec5
# ╠═c51a8c83-0317-49e1-a08c-0ce4859bfdfa
# ╠═4a976529-1ad1-4649-a50a-031a51fa9508
# ╠═9227d920-be3f-40e5-8821-02e2d93ca97e
# ╠═5ecf013a-bfaf-41ac-bd76-2c6fa4cbddbc
