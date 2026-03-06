### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100 # number of orbits to plot

# ╔═╡ 15863916-6601-4f45-9f45-4cd303bbcc4d
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits",)

# ╔═╡ ff6522d0-84cb-4521-8400-61c02973d535
# lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/vasiliev24_L3M11")) |> reverse

# ╔═╡ cb4788e1-5b01-4387-9188-745890b48217
function plot_yz!(ax, positions::AbstractMatrix;
        plot=lines!, kwargs...)

    # x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]

    plot(ax, y, z; kwargs...)
end

# ╔═╡ 185947eb-012f-4c7e-ae2a-fd04ef58e8f7
function plot_yz!(ax, orbit::Orbit; time_min=-Inf, plot_today=true, kwargs...)
    time_filt = orbit.times .> time_min
    plot_yz!(ax, orbit.positions[:, time_filt]; kwargs...)

	if plot_today

		scatter!(ax, orbit.positions[2, 1], orbit.positions[3, 1]; color=kwargs[:color])
	end
end

# ╔═╡ 6c1b2974-df12-4b7c-a504-8223f893c50c
function plot_rt_point!(ax, orbit::Orbit; time_min=-Inf, plot_today=true, kwargs...)
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

# ╔═╡ 3575de4d-797d-4822-a28a-a4424281c277
orbits_boo3 = read_orbits("bootes3", "EP2020")

# ╔═╡ 7463aada-bf2d-41df-8968-03d027b499a7
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ c63a1c4c-6171-4853-906a-54ba74fb0766
best_mw = OrbitUtils.load_best_orbit(modelnames["bootes3"][1:2]...) |> reverse

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

# ╔═╡ 123584a9-e9b4-4c48-acb5-655bd60aaeb6
sw = @lift $(theme(:linewidth)) / 2

# ╔═╡ 4ee69fe0-ab61-457e-86f3-e6b4b3228527
function plot_yz_sun!(;kwargs...)
    X_SUN = [-8.1219733661223, 0.0, 0.0208]

    scatter!(X_SUN[2], X_SUN[3]; marker=:star5, color=COLORS[9], label="the Sun", strokewidth=sw, kwargs...)
end

# ╔═╡ 51b371cc-6ba6-4cd5-8060-42261971cf91
kwargs_mw_only = (; color=COLORS[1], alpha=0.05, linestyle=:solid, label="Dwarf: MW-only" => (; alpha=0.5), time_min=t_min)

# ╔═╡ 7ef18abc-99ff-466a-8f68-88dc8cdefb31
kwargs_mw_lmc = (; color=COLORS[2], alpha=0.05, linestyle=:solid, label="Dwarf: MW+LMC" => (; alpha=0.5), time_min=t_min)

# ╔═╡ 4ca15028-782c-4480-afae-258d92d5e772
kwargs_lmc = (; color=COLORS[3],  label="LMC", time_min=t_min, linewidth=sw[]*3, linestyle=:dash, rasterize=rasterize)

# ╔═╡ 086494a3-a1d4-4558-bd1f-7ad3a3f21562
kwargs_best = (; color=:black, label="N-body: MW-only", linestyle=:solid)

# ╔═╡ 833502cb-a187-4f14-8982-4ec9b215639f
kwargs_best_lmc = (; color=:black, label="N-body: MW+LMC", linestyle=:dot, rasterize=rasterize, plot_today=false)

# ╔═╡ 9340b2c1-bf46-4959-8d48-9117926bbc1f
let
	fig = Figure(size=(3.5*72, 2*72))

	ax_yz = Axis(fig[1,1],
		limits = 150 .* (-1, 1, -1, 1),
		xlabel = "y / kpc",
		ylabel = "z / kpc",
)	
	plot_yz!(ax_yz, orbits_boo3; kwargs_mw_only...)
	plot_yz!(ax_yz, best_mw; kwargs_best...)
	plot_yz_sun!()

	ax_rt = Axis(fig[1,2],
		xlabel = "time / Gyr",
		ylabel = L"$r_\textrm{Galcen}$ / kpc",
		limits=(-5.5, 0.5, 0, 150),
				   )	

	plot_rt!(ax_rt, orbits_boo3; rasterize=rasterize, kwargs_mw_only...)
	plot_rt!(ax_rt, best_mw; kwargs_best...)


	# Legend(fig[1, 3], ax_scl, tellwidth=true, tellheight=false, nbanks=1)

	
	rowsize!(fig.layout, 1, Aspect(1, 1))
	resize_to_layout!(fig)

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
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═35ce583b-0938-429e-af5d-b17b399f6690
# ╠═15863916-6601-4f45-9f45-4cd303bbcc4d
# ╠═49cca906-53e3-4531-a936-119d8b372c61
# ╠═ff6522d0-84cb-4521-8400-61c02973d535
# ╠═cb4788e1-5b01-4387-9188-745890b48217
# ╠═185947eb-012f-4c7e-ae2a-fd04ef58e8f7
# ╠═6c1b2974-df12-4b7c-a504-8223f893c50c
# ╠═82c764ab-1194-48d2-a7aa-1272dbe4682d
# ╠═4da40bb0-18d7-4e79-b031-1319a3296207
# ╠═3d417fe3-a75d-476e-a4f9-8ce25914c473
# ╠═3f3df4af-3f2b-4803-b007-4302a029f988
# ╠═3575de4d-797d-4822-a28a-a4424281c277
# ╠═7463aada-bf2d-41df-8968-03d027b499a7
# ╠═c63a1c4c-6171-4853-906a-54ba74fb0766
# ╠═4ee69fe0-ab61-457e-86f3-e6b4b3228527
# ╟─5ec0129c-3075-44f1-bcdf-7090484bcd8d
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═59bb1f11-987d-4e2f-bb07-6905cd09a3f2
# ╠═123584a9-e9b4-4c48-acb5-655bd60aaeb6
# ╠═51b371cc-6ba6-4cd5-8060-42261971cf91
# ╠═7ef18abc-99ff-466a-8f68-88dc8cdefb31
# ╠═4ca15028-782c-4480-afae-258d92d5e772
# ╠═086494a3-a1d4-4558-bd1f-7ad3a3f21562
# ╠═833502cb-a187-4f14-8982-4ec9b215639f
# ╠═9340b2c1-bf46-4959-8d48-9117926bbc1f
