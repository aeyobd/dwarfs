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

# ╔═╡ 7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
md"""
## Setup
"""

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
FIGDIR = "./figures"

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

# ╔═╡ 6aec88e1-68a9-41e9-93ce-655434904386
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100 # number of orbits to plot

# ╔═╡ 15863916-6601-4f45-9f45-4cd303bbcc4d
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits","sculptor", "vasiliev24_L3M11_9Gyr")

# ╔═╡ b919812e-2232-470d-9a42-f61f37fbc5ae
# modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", modelnames["scl_lmc_orbits"]...)

# ╔═╡ 9efc0091-ea97-439a-bbf4-5c8b5f1127fc
modeldir_no = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_M11")

# ╔═╡ fd7246ec-212c-4255-a658-61bd82b9b7e0
# modeldir_best = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_L3M11_2x_special_cases")

# ╔═╡ 453d7f2e-91a4-46b5-9697-a4548f145cec
# modeldir_best_long = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/1e6_new_v31_r3.2/", "L3M11_9Gyr_smallperi.3b")

# ╔═╡ ff6522d0-84cb-4521-8400-61c02973d535
lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/vasiliev24_L3M11")) |> reverse

# ╔═╡ e201e22e-bab4-4da9-8273-a59ce73f83a3
# lmc_orbit_L10 = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L3M10")) |> reverse

# ╔═╡ 9028c3ec-c8fd-4b2b-a76d-4c16682e24db
# lmc_orbit_light = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L2M11")) |> reverse

# ╔═╡ 4da40bb0-18d7-4e79-b031-1319a3296207
# lmc_orbit_M10 = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L3M10")) |> reverse

# ╔═╡ 3d417fe3-a75d-476e-a4f9-8ce25914c473
function read_orbits(modeldir)
	structs = LilGuys.read_ordered_structs(joinpath(modeldir, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	first.(structs)[filt], last.(structs)[filt]
end

# ╔═╡ cc852a14-63de-4094-821b-b5ed81fd9b7e
idx, orbits = read_orbits(modeldir)

# ╔═╡ e96f758a-1cb9-436e-b351-cfa311520faa
idx_no, orbits_no = read_orbits(modeldir_no)

# ╔═╡ c63a1c4c-6171-4853-906a-54ba74fb0766
best_orbit = OrbitUtils.load_best_orbit(modelnames["scl_lmc"][1:2]...) |> reverse

# ╔═╡ 782e9dae-b15f-40fd-9d75-9a334810aa81
best_m_lmc = best_orbit - LilGuys.resample(lmc_orbit, best_orbit.times)

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

# ╔═╡ d6466c39-d703-49a5-a97b-26f51d0cdd85
# pos_lmc_light_resampled = LilGuys.resample(lmc_orbit_light, orbits_light[1].times)

# ╔═╡ c6945657-7bf1-431f-a9e8-9b4dce2210f1
# pos_lmc_M10_resampled = LilGuys.resample(lmc_orbit_M10, orbits_light[1].times)

# ╔═╡ 123584a9-e9b4-4c48-acb5-655bd60aaeb6
sw = @lift $(theme(:linewidth)) / 2

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let
	fig = Figure(figsize=(6*72, 5*72))
	limits = LilGuys.limits_xyz(LilGuys.positions.(orbits)...)

	ax_xyz = axes_xyz_flat(fig, limits=limits)

	plot_xyz!(ax_xyz, orbits, color=COLORS[2], alpha=0.05, time_min=t_min, linestyle=:solid)

	plot_xyz!(ax_xyz, orbits_no, color=COLORS[3], alpha=0.05, time_min=t_min, linestyle=:solid)

	plot_xyz!(ax_xyz, lmc_orbit, time_min=t_min, linestyle=:solid, color=COLORS[1])


	plot_xyz!(ax_xyz, best_orbit, time_min=t_min, linestyle=:solid, color=:black)
	#plot_xyz!(ax_xyz, best_long_orbit, time_min=t_min, linestyle=:dash, color=:black)
	OrbitUtils.plot_xyz_today!(ax_xyz, best_orbit, strokewidth=sw, color=:black)
	OrbitUtils.plot_xyz_today!(ax_xyz, lmc_orbit, strokewidth=sw, color=COLORS[1])
	
	OrbitUtils.plot_xyz_sun!(ax_xyz, strokewidth=sw)


	# MW distance
	ax_rt = Axis(fig[2,1:3 ],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{sat-MW}$ / kpc")
	
	plot_rt!(ax_rt, orbits, color=COLORS[2], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, orbits_no, color=COLORS[3], alpha=0.05, linestyle=:solid)
	plot_rt!(ax_rt, lmc_orbit, color=COLORS[1])

	plot_rt!(ax_rt, best_orbit, linestyle=:solid, color=:black)
	#plot_rt!(ax_rt, best_long_orbit, linestyle=:dash, color=:black)
	plot_rt_today!(ax_rt, best_orbit, strokewidth=sw, color=:black)
	plot_rt_today!(ax_rt, lmc_orbit, strokewidth=sw, color=COLORS[1])
	hidexdecorations!(ticks=false, minorticks=false)

	

	
	xlims!(t_min*T2GYR, 0.2)
	ylims!(0, 300)

	

	# LMC distance
	ax_rt2 = Axis(fig[3, 1:3 ],
				xlabel = "time / Gyr", ylabel = L"$r_\textrm{sat-LMC}$ / kpc")
	
	plot_rt!(ax_rt2, orbits .- [pos_lmc_resampled], color=COLORS[2], alpha=0.05, linestyle=:solid)
	xlims!(t_min*T2GYR, 0.2)
	ylims!(0, 200)

	plot_rt!(ax_rt2, best_m_lmc, linestyle=:solid, color=:black)
	#plot_rt!(ax_rt2, best_long_m_lmc, linestyle=:dash, color=:black)
	plot_rt_today!(ax_rt2, best_m_lmc, strokewidth=sw, color=:black)


	# labels
	lines!([NaN], [NaN], color=COLORS[3], alpha=0.5, label="Scl, MW only")
	lines!([NaN], [NaN], color=COLORS[2], alpha=0.5, label="Scl, MW+LMC")
	lines!([NaN], [NaN], color=:black, label="Scl, selected")

	lines!([NaN], [NaN], color=COLORS[1], label="LMC")

	axislegend(position=:lt, backgroundcolor=(:white, 0.8))
	#Legend(fig[4, 2], ax_rt, tellwidth=false)

	Label(fig[0, :], "Sculptor (with LMC)", font=:bold)

	rowsize!(fig.layout, 2, Aspect(1, 1.0))
	rowsize!(fig.layout, 3, Aspect(1, 1.0))

	resize_to_layout!(fig)
	@savefig "scl_lmc_xyzr_orbits"
	fig

end

# ╔═╡ Cell order:
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
# ╠═4c70700e-a8dd-4585-b31d-598f05d615e8
# ╠═2bce531e-eaf1-4258-9ca7-9a05751cbd5b
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═6b0b6c8f-eef9-4c45-8227-d1906fe6b80b
# ╠═6aec88e1-68a9-41e9-93ce-655434904386
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═35ce583b-0938-429e-af5d-b17b399f6690
# ╠═15863916-6601-4f45-9f45-4cd303bbcc4d
# ╠═b919812e-2232-470d-9a42-f61f37fbc5ae
# ╠═9efc0091-ea97-439a-bbf4-5c8b5f1127fc
# ╠═fd7246ec-212c-4255-a658-61bd82b9b7e0
# ╠═453d7f2e-91a4-46b5-9697-a4548f145cec
# ╠═49cca906-53e3-4531-a936-119d8b372c61
# ╠═ff6522d0-84cb-4521-8400-61c02973d535
# ╠═e201e22e-bab4-4da9-8273-a59ce73f83a3
# ╠═9028c3ec-c8fd-4b2b-a76d-4c16682e24db
# ╠═4da40bb0-18d7-4e79-b031-1319a3296207
# ╠═3d417fe3-a75d-476e-a4f9-8ce25914c473
# ╠═cc852a14-63de-4094-821b-b5ed81fd9b7e
# ╠═e96f758a-1cb9-436e-b351-cfa311520faa
# ╠═c63a1c4c-6171-4853-906a-54ba74fb0766
# ╠═782e9dae-b15f-40fd-9d75-9a334810aa81
# ╟─5ec0129c-3075-44f1-bcdf-7090484bcd8d
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═59bb1f11-987d-4e2f-bb07-6905cd09a3f2
# ╠═587e90ea-8597-445a-a1e9-8ec020469c35
# ╠═d6466c39-d703-49a5-a97b-26f51d0cdd85
# ╠═c6945657-7bf1-431f-a9e8-9b4dce2210f1
# ╠═123584a9-e9b4-4c48-acb5-655bd60aaeb6
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
