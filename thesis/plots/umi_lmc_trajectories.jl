### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 2b01d8f5-272e-4aa2-9825-58bb052acd10
import Agama

# ╔═╡ 00ba3075-c3e2-4965-acf3-00cda0ef320f
import TOML

# ╔═╡ ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
module OrbitUtils
	include("orbit_utils.jl")
end

# ╔═╡ a7111062-b025-43a9-bdb1-aee08deb60e9
CairoMakie.activate!(type=:png)

# ╔═╡ 6b0b6c8f-eef9-4c45-8227-d1906fe6b80b
scale_theme_element!(:linewidth, 1/2)

# ╔═╡ 4c70700e-a8dd-4585-b31d-598f05d615e8
import .OrbitUtils: axes_xyz_flat, plot_xyz!, plot_rt!, plot_rt_today!

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100 # number of orbits to plot

# ╔═╡ 15863916-6601-4f45-9f45-4cd303bbcc4d
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits/ursa_minor", "vasiliev24_L3M11")

# ╔═╡ 9efc0091-ea97-439a-bbf4-5c8b5f1127fc
modeldir_no = joinpath(ENV["DWARFS_ROOT"], "orbits/ursa_minor", "vasiliev24_M11")

# ╔═╡ ff6522d0-84cb-4521-8400-61c02973d535
lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/ursa_minor/vasiliev24_L3M11")) |> reverse

# ╔═╡ cc852a14-63de-4094-821b-b5ed81fd9b7e
idx, orbits = let
	structs = LilGuys.read_ordered_structs(joinpath(modeldir, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	first.(structs)[filt], last.(structs)[filt]
end

# ╔═╡ e96f758a-1cb9-436e-b351-cfa311520faa
idx_no, orbits_no = let
	structs = LilGuys.read_ordered_structs(joinpath(modeldir_no, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	first.(structs)[filt], last.(structs)[filt]
end

# ╔═╡ 5ec0129c-3075-44f1-bcdf-7090484bcd8d
md"""
# Plots
"""

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ 59bb1f11-987d-4e2f-bb07-6905cd09a3f2
t_min = -5 / T2GYR

# ╔═╡ 587e90ea-8597-445a-a1e9-8ec020469c35
pos_lmc_resampled = LilGuys.resample(lmc_orbit, orbits[1].times)

# ╔═╡ 123584a9-e9b4-4c48-acb5-655bd60aaeb6
sw = @lift $(theme(:linewidth)) / 2

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let
	fig = Figure(figsize=(388, 2*388))
	limits = tuple(fill((-150., 150.), 3)...)

	ax_xyz = axes_xyz_flat(fig, limits)
	for ax in ax_xyz
		ax.xticks = -100:100:100
		ax.yticks = -100:100:100
	end

	plot_xyz!(ax_xyz, orbits, color=COLORS[2], alpha=0.05, time_min=t_min, linestyle=:solid)

	plot_xyz!(ax_xyz, orbits_no, color=COLORS[3], alpha=0.05, time_min=t_min, linestyle=:solid)

	plot_xyz!(ax_xyz, lmc_orbit, time_min=t_min, linestyle=:solid, color=COLORS[1])

	OrbitUtils.plot_xyz_today!(ax_xyz, lmc_orbit, strokewidth=sw, color=COLORS[1])
	OrbitUtils.plot_xyz_today!(ax_xyz, orbits, strokewidth=sw, color=COLORS[2])
	
	OrbitUtils.plot_xyz_sun!(ax_xyz, strokewidth=sw)


	# MW distance
	ax_rt = Axis(fig[2,1:3 ],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{sat-MW}$ / kpc")
	
	plot_rt!(ax_rt, orbits, color=COLORS[2], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, orbits_no, color=COLORS[3], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, lmc_orbit, color=COLORS[1])

	plot_rt_today!(ax_rt, lmc_orbit, strokewidth=sw, color=COLORS[1])
	plot_rt_today!(ax_rt, orbits, strokewidth=sw, color=COLORS[2])

	hidexdecorations!(ticks=false, minorticks=false)
	

	
	xlims!(-9, 0.2)
	ylims!(0, 400)

	# LMC distance
	ax_rt2 = Axis(fig[3, 1:3 ],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{UMi-LMC}$ / kpc")
	
	plot_rt!(ax_rt2, orbits .- [pos_lmc_resampled], color=COLORS[2], alpha=0.05, linestyle=:solid)
	xlims!(-9, 0.2)
	plot_rt_today!(ax_rt2, orbits .- [pos_lmc_resampled], strokewidth=sw, color=COLORS[2])


	# labels
	lines!([NaN], [NaN], color=COLORS[3], alpha=0.5, label="UMi, MW only")
	lines!([NaN], [NaN], color=COLORS[2], alpha=0.5, label="UMi, MW+LMC")
	lines!([NaN], [NaN], color=COLORS[1], label="LMC")

	axislegend(position=:lt, backgroundcolor=(:white, 0.8))
	#Legend(fig[4, 2], ax_rt, tellwidth=false)

	
	rowsize!(fig.layout, 2, Aspect(1, 1.0))
	rowsize!(fig.layout, 3, Aspect(1, 1.0))

	resize_to_layout!(fig)
	@savefig "umi_lmc_xyzr_orbits"
	fig

end

# ╔═╡ Cell order:
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═2b01d8f5-272e-4aa2-9825-58bb052acd10
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
# ╠═2bce531e-eaf1-4258-9ca7-9a05751cbd5b
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═6b0b6c8f-eef9-4c45-8227-d1906fe6b80b
# ╠═4c70700e-a8dd-4585-b31d-598f05d615e8
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═35ce583b-0938-429e-af5d-b17b399f6690
# ╠═15863916-6601-4f45-9f45-4cd303bbcc4d
# ╠═9efc0091-ea97-439a-bbf4-5c8b5f1127fc
# ╠═49cca906-53e3-4531-a936-119d8b372c61
# ╠═ff6522d0-84cb-4521-8400-61c02973d535
# ╠═cc852a14-63de-4094-821b-b5ed81fd9b7e
# ╠═e96f758a-1cb9-436e-b351-cfa311520faa
# ╟─5ec0129c-3075-44f1-bcdf-7090484bcd8d
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═59bb1f11-987d-4e2f-bb07-6905cd09a3f2
# ╠═587e90ea-8597-445a-a1e9-8ec020469c35
# ╠═123584a9-e9b4-4c48-acb5-655bd60aaeb6
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
