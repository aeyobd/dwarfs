### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using CairoMakie
	using DataFrames, CSV

	using Printf
	using Revise # so we can overwrite file if anything happens
	import LilGuys as lguys

	using Arya
end


# ╔═╡ f823e80a-f6db-440f-8d25-56860618c82f
using LilGuys

# ╔═╡ 2bce531e-eaf1-4258-9ca7-9a05751cbd5b
include("paper_style.jl")

# ╔═╡ 49cca906-53e3-4531-a936-119d8b372c61
include(joinpath(ENV["DWARFS_ROOT"], "orbits/orbit_utils.jl"))

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
FIGDIR = "./figures"

# ╔═╡ 7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
md"""
## Setup
"""

# ╔═╡ 00ba3075-c3e2-4965-acf3-00cda0ef320f
import TOML

# ╔═╡ ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
module OrbitUtils
	include("orbit_utils.jl")
end

# ╔═╡ 4c70700e-a8dd-4585-b31d-598f05d615e8
import .OrbitUtils: axes_xyz_flat, plot_xyz!, plot_rt!, plot_rt_today!

# ╔═╡ a7111062-b025-43a9-bdb1-aee08deb60e9
CairoMakie.activate!(type=:png, px_per_unit=2)

# ╔═╡ 6b0b6c8f-eef9-4c45-8227-d1906fe6b80b
scale_theme_element!(:linewidth, 1/2)

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100 # number of orbits to plot

# ╔═╡ 3d417fe3-a75d-476e-a4f9-8ce25914c473
function read_orbits(modeldir)
	structs = LilGuys.read_ordered_structs(joinpath(modeldir, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	last.(structs)[filt]
end

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 15863916-6601-4f45-9f45-4cd303bbcc4d
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_L3M11_9Gyr")

# ╔═╡ 9efc0091-ea97-439a-bbf4-5c8b5f1127fc
modeldir_no = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_M11")

# ╔═╡ ff6522d0-84cb-4521-8400-61c02973d535
lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/vasiliev24_L3M11")) |> reverse

# ╔═╡ e201e22e-bab4-4da9-8273-a59ce73f83a3
lmc_orbit_M10 = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L3M10")) |> reverse

# ╔═╡ 9028c3ec-c8fd-4b2b-a76d-4c16682e24db
lmc_orbit_light = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L2M11")) |> reverse

# ╔═╡ cc852a14-63de-4094-821b-b5ed81fd9b7e
orbits = read_orbits(modeldir)

# ╔═╡ e96f758a-1cb9-436e-b351-cfa311520faa
orbits_no = read_orbits(modeldir_no)

# ╔═╡ 587e90ea-8597-445a-a1e9-8ec020469c35
pos_lmc_resampled = LilGuys.resample(lmc_orbit, orbits[1].times)

# ╔═╡ 46e1e951-55da-4c2f-9565-330b8853c7fc
orbits_light = read_orbits(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "L2M11"))

# ╔═╡ d6466c39-d703-49a5-a97b-26f51d0cdd85
pos_lmc_light_resampled = LilGuys.resample(lmc_orbit_light, orbits_light[1].times)

# ╔═╡ c6945657-7bf1-431f-a9e8-9b4dce2210f1
pos_lmc_M10_resampled = LilGuys.resample(lmc_orbit_M10, orbits_light[1].times)

# ╔═╡ 633275fd-fe32-483e-8612-fbb27e84bbb3
orbits_M10 = read_orbits(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "L3M10"))

# ╔═╡ 7beb5358-f1a9-4eb4-9ee1-0f946c9de388
orbits_v21 = read_orbits(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev2021"))

# ╔═╡ 0d342a60-f6b5-4466-b0a7-ea7d3fd01a21
orbits_umi = read_orbits(joinpath(ENV["DWARFS_ROOT"], "orbits/ursa_minor", "vasiliev24_L3M11"))

# ╔═╡ a24e6f3d-5c26-4d20-8f8f-9f7b328c64b1
orbits_umi_light = read_orbits(joinpath(ENV["DWARFS_ROOT"], "orbits/ursa_minor", "L2M11"))

# ╔═╡ 25ad3b62-3580-4b9f-81fb-df7217e955dd
orbits_umi_M10 = read_orbits(joinpath(ENV["DWARFS_ROOT"], "orbits/ursa_minor", "L3M10"))

# ╔═╡ c35bf067-73d7-4026-897c-5a737d0d01e7
md"""
# LMC mass effects
"""

# ╔═╡ e411b2d6-fc0d-4395-94cc-6442ff07fb52
let
	fig = Figure()

	# MW distance
	ax_rt = Axis(fig[1,1:2],
				 ylabel = L"$r_\textrm{sat-MW}$ / kpc",
				limits=(-9, 0.2, 0, 400))
	
	plot_rt!(ax_rt, orbits, color=COLORS[1], alpha=0.05, 
			 linestyle=:solid, label="L3M11" => (;alpha=0.5))
	plot_rt!(ax_rt, orbits_light, color=COLORS[2], alpha=0.05, 
			 label="L2M11" => (;alpha=0.5))
	plot_rt!(ax_rt, orbits_M10, color=COLORS[3], alpha=0.05,  label="L3M10" => (;alpha=0.5))
	
	plot_rt!(ax_rt, lmc_orbit, color=COLORS[1])
	plot_rt!(ax_rt, lmc_orbit_light, color=COLORS[2])
	plot_rt!(ax_rt, lmc_orbit_M10, color=COLORS[3])


	axislegend(position=:rt, merge=true, unique=true)
	text!(0.15, 0.9, align=(:left, :center), text="Sculptor", space=:relative, font=:bold)

	ax_rt.xticklabelsvisible=false

# MW distance
	ax_rt_umi = Axis(fig[2,1:2],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{sat-MW}$ / kpc",
				limits=(-9.0, 0.2, 0, 400))
	
	plot_rt!(ax_rt_umi, orbits_umi, color=COLORS[1], alpha=0.05, 
			 linestyle=:solid, label="L3M11" => (;alpha=0.5))
	plot_rt!(ax_rt_umi, orbits_umi_light, color=COLORS[2], alpha=0.05, 
			 label="L2M11" => (;alpha=0.5))
	plot_rt!(ax_rt_umi, orbits_umi_M10, color=COLORS[3], alpha=0.05,  label="L3M10" => (;alpha=0.5))
	# plot_rt!(ax_rt, orbits_v21, color=COLORS[4], alpha=0.05,  label="L3M10" => (;alpha=0.5))
	
	plot_rt!(ax_rt_umi, lmc_orbit, color=COLORS[1])
	plot_rt!(ax_rt_umi, lmc_orbit_light, color=COLORS[2])
	plot_rt!(ax_rt_umi, lmc_orbit_M10, color=COLORS[3])


	l1 = lines!([NaN], [NaN], color=:black, linestyle=:solid, alpha=0.2, linewidth=theme(:linewidth)[]/2)
	l2 = lines!([NaN], [NaN], color=:black, linestyle=:solid, alpha=1)

	axislegend(ax_rt_umi, [l1, l2], ["orbit sample", "LMC"], position=:rt)
	text!(0.15, 0.9, align=(:left, :center), text="Ursa Minor", space=:relative, font=:bold)


	@savefig "scl_lmc_orbits_mass_effect"
	fig
end

# ╔═╡ dfb3b17e-cbd0-459e-bb7d-2e0cb2708c5f
let
	fig = Figure()

	# MW distance
	ax_rt = Axis(fig[1,1:2],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{sat-LMC}$ / kpc",
				limits=(-9, 0, 0, nothing))
	
	plot_rt!(ax_rt, orbits .- [pos_lmc_resampled], color=COLORS[2], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, orbits_light .- [pos_lmc_light_resampled], color=COLORS[1], alpha=0.05)

	plot_rt!(ax_rt, orbits_M10 .- [pos_lmc_M10_resampled], color=COLORS[3], alpha=0.05)


	fig
end

# ╔═╡ 96cfb827-131b-4e3e-90db-7f91f34e9d31
let
	fig = Figure()

	# MW distance
	ax_rt = Axis(fig[1,1:2],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{sat-LMC}$ / kpc",
				limits=(-10, 0, 0, nothing))
	
	plot_rt!(ax_rt, orbits_umi .- [pos_lmc_resampled], color=COLORS[2], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, orbits_umi_light .- [pos_lmc_light_resampled], color=COLORS[1], alpha=0.05)
	plot_rt!(ax_rt, orbits_umi_M10 .- [pos_lmc_M10_resampled], color=COLORS[3], alpha=0.05)



	fig
end

# ╔═╡ Cell order:
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
# ╠═4c70700e-a8dd-4585-b31d-598f05d615e8
# ╠═2bce531e-eaf1-4258-9ca7-9a05751cbd5b
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═6b0b6c8f-eef9-4c45-8227-d1906fe6b80b
# ╠═3d417fe3-a75d-476e-a4f9-8ce25914c473
# ╠═35ce583b-0938-429e-af5d-b17b399f6690
# ╠═49cca906-53e3-4531-a936-119d8b372c61
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═15863916-6601-4f45-9f45-4cd303bbcc4d
# ╠═9efc0091-ea97-439a-bbf4-5c8b5f1127fc
# ╠═ff6522d0-84cb-4521-8400-61c02973d535
# ╠═e201e22e-bab4-4da9-8273-a59ce73f83a3
# ╠═9028c3ec-c8fd-4b2b-a76d-4c16682e24db
# ╠═cc852a14-63de-4094-821b-b5ed81fd9b7e
# ╠═e96f758a-1cb9-436e-b351-cfa311520faa
# ╠═587e90ea-8597-445a-a1e9-8ec020469c35
# ╠═d6466c39-d703-49a5-a97b-26f51d0cdd85
# ╠═c6945657-7bf1-431f-a9e8-9b4dce2210f1
# ╠═46e1e951-55da-4c2f-9565-330b8853c7fc
# ╠═633275fd-fe32-483e-8612-fbb27e84bbb3
# ╠═7beb5358-f1a9-4eb4-9ee1-0f946c9de388
# ╠═0d342a60-f6b5-4466-b0a7-ea7d3fd01a21
# ╠═a24e6f3d-5c26-4d20-8f8f-9f7b328c64b1
# ╠═25ad3b62-3580-4b9f-81fb-df7217e955dd
# ╟─c35bf067-73d7-4026-897c-5a737d0d01e7
# ╠═e411b2d6-fc0d-4395-94cc-6442ff07fb52
# ╠═dfb3b17e-cbd0-459e-bb7d-2e0cb2708c5f
# ╠═96cfb827-131b-4e3e-90db-7f91f34e9d31
