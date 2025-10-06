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

# ╔═╡ 6ca3fe17-3f13-43fe-967b-881078135ead
modelname = "EP2020"

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

# ╔═╡ 700dcfed-77bd-4a46-9389-e114d5d47221
import .OrbitUtils: axes_xyz_flat, plot_xyz_sun!, plot_xyz_today!, plot_xyz!, plot_rt_today!, plot_rt!

# ╔═╡ 75e9df93-6fea-4977-b4bc-535220e955c7
scale_theme_element!(:linewidth, 1/2)

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 180e128f-7c3f-4682-90fd-da8e752c1018
function load_best_orbit(galaxyname, modelname="EP2020_special_cases", orbitname="orbit_smallperi",)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)

	best_orbit = LilGuys.Orbit(joinpath(modeldir, orbitname * ".csv"))
end

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100 # number of orbits to plot

# ╔═╡ cc852a14-63de-4094-821b-b5ed81fd9b7e
function load_orbits(galaxyname, modelname="EP2020")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)
	
	idx, orbits = let
		structs = LilGuys.read_ordered_structs(joinpath(modeldir, "orbits.hdf5"), LilGuys.Orbit)
	
		filt = 1:min(Nmax, length(structs))
		first.(structs)[filt], last.(structs)[filt]
	end

	return idx, orbits
end

# ╔═╡ 15863916-6601-4f45-9f45-4cd303bbcc4d
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", modelname)

# ╔═╡ 10acf60c-8acf-4c53-9c18-fbd4692b28e3
_, orbits = load_orbits("sculptor")

# ╔═╡ 127a988e-74ab-41f4-8979-6800622f3ff4
scl_best = load_best_orbit("sculptor")

# ╔═╡ 5ec0129c-3075-44f1-bcdf-7090484bcd8d
md"""
# Plots
"""

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let
	fig = Figure()
	limits = LilGuys.limits_xyz(LilGuys.positions.(orbits)...)

	ax_xyz = axes_xyz_flat(fig, limits=limits)

	for ax in ax_xyz
		ax.xticks = -100:100:100
		ax.yticks = -100:100:100
	end

	plot_xyz!(ax_xyz, orbits, color=COLORS[3], alpha=0.05, linestyle=:solid)
	plot_xyz!(ax_xyz, scl_best, color=:black,)
	plot_xyz_today!(ax_xyz, scl_best, length(scl_best))
	plot_xyz_sun!(ax_xyz, strokewidth=OrbitUtils.strokewidth)

	ax_rt = Axis(fig[2, 1:3],
				xlabel = "time / Gyr", ylabel =   L"$r_\textrm{sat-MW}$ / kpc")
	
	plot_rt!(ax_rt, orbits, color=COLORS[3], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, scl_best, color=:black)
	xlims!(-10, 0.1)
	ylims!(0, 120)
	plot_rt_today!(ax_rt, scl_best, length(scl_best))

	rowsize!(fig.layout, 2, Aspect(1, 1.0))

	resize_to_layout!(fig)
	@savefig "scl_xyzr_orbits"
	fig

end

# ╔═╡ Cell order:
# ╠═6ca3fe17-3f13-43fe-967b-881078135ead
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═2b01d8f5-272e-4aa2-9825-58bb052acd10
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═2bce531e-eaf1-4258-9ca7-9a05751cbd5b
# ╠═ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═700dcfed-77bd-4a46-9389-e114d5d47221
# ╠═75e9df93-6fea-4977-b4bc-535220e955c7
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═cc852a14-63de-4094-821b-b5ed81fd9b7e
# ╠═180e128f-7c3f-4682-90fd-da8e752c1018
# ╠═35ce583b-0938-429e-af5d-b17b399f6690
# ╠═15863916-6601-4f45-9f45-4cd303bbcc4d
# ╠═10acf60c-8acf-4c53-9c18-fbd4692b28e3
# ╠═127a988e-74ab-41f4-8979-6800622f3ff4
# ╟─5ec0129c-3075-44f1-bcdf-7090484bcd8d
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
