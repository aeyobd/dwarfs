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

# ╔═╡ ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
include("orbit_utils.jl")

# ╔═╡ 49cca906-53e3-4531-a936-119d8b372c61
include(joinpath(ENV["DWARFS_ROOT"], "orbits/orbit_utils.jl"))

# ╔═╡ 7450144e-5464-4036-a215-b6e2cd270405
md"""
This notebook analyzes the result of the MC samples of orbits in the same potential to determine the plausable range of pericentres and apocentres.
In general, this notebook is meant to plot and calculate main properties.
If you would like to investigate a special case further, it is likely best to make a new notebook in the target analysis directory.
Additionally, see analyze_lmc.jl in this directory for a version which also plots additional plots for an MW-LMC potential.
"""

# ╔═╡ 2b9d49c6-74cc-4cce-b29e-04e94776863f
md"""
The most important variable is to set the modelname to the appropriate directory.
"""

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
FIGDIR = "./figures"

# ╔═╡ 43e80dd7-aa44-43e6-adcc-ebb7ae9e9eb8
t_max = -2/T2GYR

# ╔═╡ 7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
md"""
## Setup
"""

# ╔═╡ 2b01d8f5-272e-4aa2-9825-58bb052acd10
import Agama

# ╔═╡ 00ba3075-c3e2-4965-acf3-00cda0ef320f
import TOML

# ╔═╡ a7111062-b025-43a9-bdb1-aee08deb60e9
CairoMakie.activate!(type=:png)

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100 # number of orbits to plot

# ╔═╡ 15863916-6601-4f45-9f45-4cd303bbcc4d
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_L3M11")

# ╔═╡ 9efc0091-ea97-439a-bbf4-5c8b5f1127fc
modeldir_no = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_M11")

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

# ╔═╡ ff6522d0-84cb-4521-8400-61c02973d535
lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/vasiliev24_L3M11"))

# ╔═╡ 5ec0129c-3075-44f1-bcdf-7090484bcd8d
md"""
# Plots
"""

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ 59bb1f11-987d-4e2f-bb07-6905cd09a3f2
t_min = -2 / T2GYR

# ╔═╡ 587e90ea-8597-445a-a1e9-8ec020469c35
pos_lmc_resampled = LilGuys.resample(lmc_orbit, orbits[1].times)

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let
	fig = Figure()
	limits = LilGuys.limits_xyz(LilGuys.positions.(orbits)...)

	ax_xyz = axes_xyz_flat(fig, limits)

	plot_xyz!(ax_xyz, orbits, color=COLORS[1], alpha=0.05, time_min=t_min, linestyle=:solid)
	plot_xyz!(ax_xyz, orbits_no, color=COLORS[2], alpha=0.05, time_min=t_min, linestyle=:solid)
	plot_xyz!(ax_xyz, lmc_orbit, color=lmc_orbit.times[lmc_orbit.times .> t_min], time_min=t_min, linestyle=:solid, colorrange=(-10/T2GYR, 0), colormap=Reverse(:greys))
	# plot_xyz!(ax_xyz, best_orbit, color=:black, time_min=-5/T2GYR)
	# plot_xyz_today!(ax_xyz, best_orbit, length(best_orbit))

	gs = GridLayout(fig[2, 1:3])
	ax_rt = Axis(gs[1,1 ],
				xlabel = "time / Gyr", ylabel = "MW distance / kpc")
	
	plot_rt!(ax_rt, orbits, color=COLORS[1], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, orbits_no, color=COLORS[2], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, lmc_orbit, color=lmc_orbit.times, linestyle=:solid, colorrange=(-10/T2GYR, 0), colormap=Reverse(:greys))


	# labels
	lines!([NaN], [NaN], color=COLORS[2], alpha=0.5, label="Sculptor")
	lines!([NaN], [NaN], color=COLORS[1], alpha=0.5, label="+LMC")
	lines!([NaN], [NaN], color=:black, label="LMC")

	axislegend(position=:lt)
	
	# plot_rt!(ax_rt, best_orbit, color=:black)
	xlims!(-10, 0)
	ylims!(0, 400)
	# plot_rt_today!(ax_rt, best_orbit, length(best_orbit))


	ax_rt2 = Axis(gs[1,2 ],
				xlabel = "time / Gyr", ylabel = "LMC distance / kpc")
	
	plot_rt!(ax_rt2, orbits .- [pos_lmc_resampled], color=COLORS[3], alpha=0.05, linestyle=:solid)

	rowsize!(fig.layout, 2, Relative(0.5))

	@savefig "scl_lmc_xyzr_orbits"
	fig

end

# ╔═╡ Cell order:
# ╟─7450144e-5464-4036-a215-b6e2cd270405
# ╟─2b9d49c6-74cc-4cce-b29e-04e94776863f
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═43e80dd7-aa44-43e6-adcc-ebb7ae9e9eb8
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═2b01d8f5-272e-4aa2-9825-58bb052acd10
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═2bce531e-eaf1-4258-9ca7-9a05751cbd5b
# ╠═ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
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
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
