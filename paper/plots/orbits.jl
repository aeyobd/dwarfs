### A Pluto.jl notebook ###
# v0.20.20

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
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits",)

# ╔═╡ b919812e-2232-470d-9a42-f61f37fbc5ae


# ╔═╡ 9efc0091-ea97-439a-bbf4-5c8b5f1127fc
modeldir_no = joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor", "vasiliev24_M11")

# ╔═╡ fd7246ec-212c-4255-a658-61bd82b9b7e0


# ╔═╡ 453d7f2e-91a4-46b5-9697-a4548f145cec


# ╔═╡ ff6522d0-84cb-4521-8400-61c02973d535
lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/vasiliev24_L3M11")) |> reverse

# ╔═╡ e201e22e-bab4-4da9-8273-a59ce73f83a3
# lmc_orbit_L10 = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L3M10")) |> reverse

# ╔═╡ 9028c3ec-c8fd-4b2b-a76d-4c16682e24db
# lmc_orbit_light = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/L2M11")) |> reverse

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

# ╔═╡ 4da40bb0-18d7-4e79-b031-1319a3296207

function plot_yz!(axes, orbits::AbstractVector{<:Orbit}; time_min=-Inf, kwargs...)
    positions = zeros(3, 0)
    for orbit in orbits
        time_filt = orbit.times .> time_min
        positions = hcat(positions, orbit.positions[:, time_filt], [NaN, NaN, NaN])
    end

    plot_yz!(axes, positions; rasterize=true, kwargs...)
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

# ╔═╡ cc852a14-63de-4094-821b-b5ed81fd9b7e
orbits_scl_lmc = read_orbits("sculptor", "vasiliev24_L3M11_9Gyr")

# ╔═╡ 3575de4d-797d-4822-a28a-a4424281c277
orbits_scl = read_orbits("sculptor", "EP2020")

# ╔═╡ e96f758a-1cb9-436e-b351-cfa311520faa
orbits_umi_lmc = read_orbits("ursa_minor", "vasiliev24_L3M11")

# ╔═╡ b98a9e78-93f4-429c-beec-879bcabd06e1
orbits_umi = read_orbits("ursa_minor", "EP2020")

# ╔═╡ c63a1c4c-6171-4853-906a-54ba74fb0766
best_scl_lmc = OrbitUtils.load_best_orbit(modelnames["scl_lmc"][1:2]...) |> reverse

# ╔═╡ 538f7a29-6cf9-455a-ab56-de3c1aee8fc9
best_scl = OrbitUtils.load_best_orbit(modelnames["scl_smallperi"][1:2]...) |> reverse

# ╔═╡ 782e9dae-b15f-40fd-9d75-9a334810aa81
best_umi = OrbitUtils.load_best_orbit(modelnames["umi_smallperi"][1:2]...) |> reverse

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

# ╔═╡ 51b371cc-6ba6-4cd5-8060-42261971cf91
kwargs_mw_only = (; color=COLORS[1], alpha=0.05, linestyle=:solid, label="Dwarf: MW-only" => (; alpha=0.5), time_min=t_min)

# ╔═╡ 7ef18abc-99ff-466a-8f68-88dc8cdefb31
kwargs_mw_lmc = (; color=COLORS[2], alpha=0.05, linestyle=:solid, label="Dwarf: MW+LMC" => (; alpha=0.5), time_min=t_min)

# ╔═╡ 4ca15028-782c-4480-afae-258d92d5e772
kwargs_lmc = (; color=COLORS[3],  label="LMC", time_min=t_min, linewidth=sw[]*3, linestyle=:dash)

# ╔═╡ 086494a3-a1d4-4558-bd1f-7ad3a3f21562
kwargs_best = (; color=:black, label="N-body: MW-only", linestyle=:solid)

# ╔═╡ 4ee69fe0-ab61-457e-86f3-e6b4b3228527
function plot_yz_sun!(;kwargs...)
    X_SUN = [-8.1219733661223, 0.0, 0.0208]

    scatter!(X_SUN[2], X_SUN[3]; marker=:star5, color=COLORS[9], label="the Sun", strokewidth=sw, kwargs...)
end

# ╔═╡ 833502cb-a187-4f14-8982-4ec9b215639f
kwargs_best_lmc = (; color=:black, label="N-body: MW+LMC", linestyle=:dot, plot_today=false)

# ╔═╡ 9340b2c1-bf46-4959-8d48-9117926bbc1f
let
	fig = Figure(size=(3.5*72, 4*72))

	ax_scl = Axis(fig[1,1],
		title = "Sculptor",
		limits = 300 .* (-1, 1, -1, 1),
		xlabel = "y / kpc",
		ylabel = "z / kpc",
		xticks = -200:200:200,
		yticks = -200:200:200,

)	
	plot_yz!(ax_scl, orbits_scl; kwargs_mw_only...)
	plot_yz!(ax_scl, orbits_scl_lmc; kwargs_mw_lmc...)
	plot_yz!(ax_scl, best_scl; kwargs_best...)
	plot_yz!(ax_scl, best_scl_lmc; kwargs_best_lmc...)
	plot_yz!(ax_scl, lmc_orbit; kwargs_lmc...)
	plot_yz_sun!()

	ax_scl_t = Axis(fig[2,1],
		xlabel = "time / Gyr",
		ylabel = L"$r_\textrm{Galcen}$ / kpc",
		limits=(-9.5, 0.5, 0, 350),
		yticks=0:100:400,
				   )	

	plot_rt!(ax_scl_t, orbits_scl; kwargs_mw_only...)
	plot_rt!(ax_scl_t, orbits_scl_lmc; kwargs_mw_lmc...)
	plot_rt_point!(ax_scl_t, best_scl; kwargs_best...)
	plot_rt_point!(ax_scl_t, best_scl_lmc; kwargs_best_lmc...)
	plot_rt_point!(ax_scl_t, lmc_orbit; kwargs_lmc...)


	ax_umi = Axis(fig[1,2], 
		title = "Ursa Minor",
		limits = 150 .* (-1, 1, -1, 1),
		xticks = -100:100:100,
		yticks = -100:100:100,
		xlabel = "y / kpc",

				 )	

	plot_yz!(ax_umi, orbits_umi; kwargs_mw_only...)
	plot_yz!(ax_umi, orbits_umi_lmc; kwargs_mw_lmc...)
	plot_yz!(ax_umi, best_umi; kwargs_best...)
	plot_yz!(ax_umi, lmc_orbit; kwargs_lmc...)
	plot_yz_sun!()

	ax_umi_t = Axis(fig[2,2],
		limits=(-9.5, 0.5, 0, 150),
		xlabel = "time / Gyr",

	)	


	plot_rt!(ax_umi_t, orbits_umi; kwargs_mw_only...)
	plot_rt!(ax_umi_t, orbits_umi_lmc; kwargs_mw_lmc...)
	plot_rt_point!(ax_umi_t, best_umi; kwargs_best...)
	plot_rt_point!(ax_umi_t, lmc_orbit; kwargs_lmc...)


	Legend(fig[3, :], ax_scl, tellwidth=false, tellheight=true, nbanks=2)

	
	rowsize!(fig.layout, 1, Aspect(1, 1))
	resize_to_layout!(fig)

	@savefig "orbits"
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
# ╠═cb4788e1-5b01-4387-9188-745890b48217
# ╠═185947eb-012f-4c7e-ae2a-fd04ef58e8f7
# ╠═6c1b2974-df12-4b7c-a504-8223f893c50c
# ╠═4da40bb0-18d7-4e79-b031-1319a3296207
# ╠═3d417fe3-a75d-476e-a4f9-8ce25914c473
# ╠═3f3df4af-3f2b-4803-b007-4302a029f988
# ╠═cc852a14-63de-4094-821b-b5ed81fd9b7e
# ╠═3575de4d-797d-4822-a28a-a4424281c277
# ╠═e96f758a-1cb9-436e-b351-cfa311520faa
# ╠═b98a9e78-93f4-429c-beec-879bcabd06e1
# ╠═c63a1c4c-6171-4853-906a-54ba74fb0766
# ╠═538f7a29-6cf9-455a-ab56-de3c1aee8fc9
# ╠═782e9dae-b15f-40fd-9d75-9a334810aa81
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
