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

# ╔═╡ 5e27a3da-1a97-46d9-aa31-3712f06b9a94
begin 
	dirname = "/cosma/home/durham/dc-boye1/sculptor/isolation/"
	snapshots = Dict{String, Any}()
	snapshots["1e4"] = lguys.Snapshot(dirname * "1e4/fiducial/out/snapshot_010.hdf5")
	snapshots["1e5"] = lguys.Snapshot(dirname * "1e5/out/snapshot_004.hdf5")
	snapshots["1e6"] = lguys.Snapshot(dirname * "1e6/out/snapshot_004.hdf5")
	snapshots["1e7"] = lguys.Snapshot(dirname * "1e7/out/snapshot_002.hdf5")

end

# ╔═╡ 26b85a2c-c57a-41dd-b295-972619ee9033
let
	fig = Figure(size=(700, 700))
	ax_rho = Axis(fig[1,1],
		xlabel="",
		ylabel=L"\log\;\rho_\textrm{DM}"
	)
	for (name, snap) in snapshots
		plot_ρ_dm!(snap, label=name)
	end

	ρ_0 = M_s / (4 * π * R_s^3)

	ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(ρ_nfw.(10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)


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

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["default (r=4)"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	zeno_snaps["break r = 2"] = lguys.Snapshot(dirname_z * "nfw_1e6_b2.hdf5")
	zeno_snaps["break r = 64"] = lguys.Snapshot(dirname_z * "nfw_1e6_b64.hdf5")
	zeno_snaps["break r = 10"] = lguys.Snapshot(dirname_z * "nfw_1e6_b10.hdf5")
	
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

# ╔═╡ 9227d920-be3f-40e5-8821-02e2d93ca97e
let

	dirname_z = "/cosma/home/durham/dc-boye1/data/dwarfs/zeno/"
	zeno_snaps = Dict{String, Any}()
	zeno_snaps["default"] = lguys.Snapshot(dirname_z * "nfw_1e6.hdf5")
	zeno_snaps["gauss"] = lguys.Snapshot(dirname_z * "nfw_1e6_gauss.hdf5")
	zeno_snaps["sw"] = lguys.Snapshot(dirname_z * "nfw_1e6_sw.hdf5")
	
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

# ╔═╡ 5ecf013a-bfaf-41ac-bd76-2c6fa4cbddbc
4096 * 16

# ╔═╡ Cell order:
# ╠═84db4ab2-13d5-11ef-2687-b14402ce4daf
# ╠═78f3d3dd-c878-4a47-9cdf-43af0478b010
# ╠═543adda2-a9ac-407a-bb3b-62b3fdcf2d44
# ╠═5e27a3da-1a97-46d9-aa31-3712f06b9a94
# ╠═26b85a2c-c57a-41dd-b295-972619ee9033
# ╠═c51a8c83-0317-49e1-a08c-0ce4859bfdfa
# ╠═4a976529-1ad1-4649-a50a-031a51fa9508
# ╠═9227d920-be3f-40e5-8821-02e2d93ca97e
# ╠═5ecf013a-bfaf-41ac-bd76-2c6fa4cbddbc
