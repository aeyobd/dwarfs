### A Pluto.jl notebook ###
# v0.20.24

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


# ╔═╡ 57db77a4-38a2-4f46-a89b-b48511664d5f
using OrderedCollections

# ╔═╡ f823e80a-f6db-440f-8d25-56860618c82f
using LilGuys

# ╔═╡ 2bce531e-eaf1-4258-9ca7-9a05751cbd5b
include("paper_style.jl")

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

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ cb4788e1-5b01-4387-9188-745890b48217
function plot_yz!(ax, positions::AbstractMatrix;
				  x_direction=2, y_direction=3,
        plot=lines!, kwargs...)

    y = positions[x_direction, :]
    z = positions[y_direction, :]

    plot(ax, y, z; kwargs...)
end

# ╔═╡ 185947eb-012f-4c7e-ae2a-fd04ef58e8f7
function plot_yz!(ax, orbit::Orbit; x_direction=2, y_direction=3, time_min=-10/T2GYR, plot_today=false, kwargs...)
    time_filt = orbit.times .> time_min
    plot_yz!(ax, orbit.positions[:, time_filt]; x_direction=x_direction, y_direction=y_direction, kwargs...)

	if plot_today

		scatter!(ax, orbit.positions[x_direction, 1], orbit.positions[y_direction, 1]; color=kwargs[:color])
	end
end

# ╔═╡ 6c1b2974-df12-4b7c-a504-8223f893c50c
function plot_rt_point!(ax, orbit::Orbit; time_min=-10/T2GYR, plot_today=false, kwargs...)
    time_filt = orbit.times .> time_min
    plot_rt!(ax, orbit.times[time_filt]*T2GYR, radii(orbit)[time_filt]; kwargs...)

	if plot_today
		scatter!(ax, orbit.times[1]*T2GYR, radii(orbit)[1]; color=kwargs[:color])
	end
end

# ╔═╡ 82c764ab-1194-48d2-a7aa-1272dbe4682d
rasterize = 4

# ╔═╡ 4da40bb0-18d7-4e79-b031-1319a3296207

function plot_yz!(axes, orbits::AbstractVector{<:Orbit}; time_min=-Inf, kwargs...)
    positions = zeros(3, 0)
    for orbit in orbits
        time_filt = orbit.times .> time_min
        positions = hcat(positions, orbit.positions[:, time_filt], [NaN, NaN, NaN])
    end

    plot_yz!(axes, positions; rasterize=rasterize, kwargs...)
end


# ╔═╡ 3d417fe3-a75d-476e-a4f9-8ce25914c473
function read_orbits(modeldir)
	structs = LilGuys.read_ordered_structs(joinpath(modeldir, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	
	last.(structs)[filt]
end

# ╔═╡ 3f3df4af-3f2b-4803-b007-4302a029f988
function read_orbits(galaxyname, modelname)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)
	return read_orbits(modeldir)
end

# ╔═╡ 7463aada-bf2d-41df-8968-03d027b499a7
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ eb102600-d9b2-42f5-87ce-c9d5861dea96
example_orbits = OrderedDict(
	peri => Orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020_special_cases/orbit_$(peri).csv")) |> reverse
	for peri in ["peri_1.5", "peri_4", "mean", "peri_12", "peri_18", "peri_26"]
)

# ╔═╡ 5ec0129c-3075-44f1-bcdf-7090484bcd8d
md"""
# Plots
"""

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ 59bb1f11-987d-4e2f-bb07-6905cd09a3f2
t_min = -3.5 / T2GYR

# ╔═╡ 123584a9-e9b4-4c48-acb5-655bd60aaeb6
sw = @lift $(theme(:linewidth)) / 2

# ╔═╡ 086494a3-a1d4-4558-bd1f-7ad3a3f21562
kwargs_best = (; colorrange=(1, length(example_orbits)), linestyle=:solid, time_min=t_min)

# ╔═╡ 833502cb-a187-4f14-8982-4ec9b215639f
kwargs_best_lmc = (; color=:black, label="N-body: MW+LMC", linestyle=:dot, rasterize=rasterize, plot_today=false, time_min=t_min)

# ╔═╡ 2a088b35-a974-4510-bd31-bdf25b0f781b
labels = OrderedDict(
	"peri_1.5" => "1.5",
	"peri_4" => "4",
	"mean" => "7 (mean)",
	"peri_12" => "12",
	"peri_18" => "18",
	"peri_26" => "26",
)

# ╔═╡ 9340b2c1-bf46-4959-8d48-9117926bbc1f
let
	fig = Figure(size=(6*72, 3*72))
	ax_xy = Axis(fig[1,1],
		limits = 125 .* (-1, 1, -1, 1),
		xlabel = "x / kpc",
		ylabel = "y / kpc",
)	
	
	for (i, (peri, orbit)) in enumerate(example_orbits)
		plot_yz!(ax_xy, orbit; label=labels[peri], color=i, x_direction=1, y_direction=2, kwargs_best...)
	end

	
	ax_yz = Axis(fig[1,2],
		limits = 125 .* (-1, 1, -1, 1),
		xlabel = "y / kpc",
		ylabel = "z / kpc",
)	
	for (i, (peri, orbit)) in enumerate(example_orbits)
		plot_yz!(ax_yz, orbit; color=i, kwargs_best...)
	end

	ax_rt = Axis(fig[1,3],
		xlabel = "time / Gyr",
		ylabel = L"$r_\textrm{Galcen}$ / kpc",
		limits=(-4, 0.5, 0, 150),
				 xticks=-4:1:0
				   )	

	for (i, (peri, orbit)) in enumerate(example_orbits)
		plot_rt!(ax_rt, orbit; color=i, kwargs_best...)
	end


	
	rowsize!(fig.layout, 1, Aspect(1, 1))

	Legend(fig[2, :], ax_xy, "pericentre / kpc", tellwidth=false, tellheight=true, nbanks=6)
	resize_to_layout!(fig)

	@savefig "boo3_orbits"
	fig

end

# ╔═╡ Cell order:
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═57db77a4-38a2-4f46-a89b-b48511664d5f
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═ff577282-1d04-4f6a-bb4e-74cf5a8d51e3
# ╠═4c70700e-a8dd-4585-b31d-598f05d615e8
# ╠═2bce531e-eaf1-4258-9ca7-9a05751cbd5b
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═6b0b6c8f-eef9-4c45-8227-d1906fe6b80b
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═cb4788e1-5b01-4387-9188-745890b48217
# ╠═185947eb-012f-4c7e-ae2a-fd04ef58e8f7
# ╠═6c1b2974-df12-4b7c-a504-8223f893c50c
# ╠═82c764ab-1194-48d2-a7aa-1272dbe4682d
# ╠═4da40bb0-18d7-4e79-b031-1319a3296207
# ╠═3d417fe3-a75d-476e-a4f9-8ce25914c473
# ╠═3f3df4af-3f2b-4803-b007-4302a029f988
# ╠═7463aada-bf2d-41df-8968-03d027b499a7
# ╠═eb102600-d9b2-42f5-87ce-c9d5861dea96
# ╟─5ec0129c-3075-44f1-bcdf-7090484bcd8d
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═59bb1f11-987d-4e2f-bb07-6905cd09a3f2
# ╠═123584a9-e9b4-4c48-acb5-655bd60aaeb6
# ╠═086494a3-a1d4-4558-bd1f-7ad3a3f21562
# ╠═833502cb-a187-4f14-8982-4ec9b215639f
# ╠═2a088b35-a974-4510-bd31-bdf25b0f781b
# ╠═9340b2c1-bf46-4959-8d48-9117926bbc1f
