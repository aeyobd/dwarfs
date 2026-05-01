### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using CairoMakie

	import LilGuys as lguys
	using PyFITS
	
	using Arya
end


# ╔═╡ 48257649-efa8-448d-9296-d49517137277
using OrderedCollections

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

# ╔═╡ 4da40bb0-18d7-4e79-b031-1319a3296207

function plot_yz!(axes, orbits::AbstractVector{<:Orbit}; time_min=-Inf, kwargs...)
    positions = zeros(3, 0)
    for orbit in orbits
        time_filt = orbit.times .> time_min
        positions = hcat(positions, orbit.positions[:, time_filt], [NaN, NaN, NaN])
    end

    plot_yz!(axes, positions; rasterize=rasterize, kwargs...)
end


# ╔═╡ 3f3df4af-3f2b-4803-b007-4302a029f988
function read_orbits(galaxyname, modelname)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)
	return read_fits(joinpath(modeldir, "orbital_properties.fits"))
end

# ╔═╡ 3575de4d-797d-4822-a28a-a4424281c277
orbit_props = read_orbits("bootes3", "EP2020")

# ╔═╡ 0700fa00-3b38-476d-be4c-c5c61037ee3e
example_orbits = OrderedDict(
	peri => TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020_special_cases/orbit_$(peri).toml"))
	for peri in ["peri_1.5", "peri_4", "mean", "peri_12", "peri_18", "peri_26"]
)

# ╔═╡ b5ae57dc-c148-4dca-983e-a600e0c46f57
md"""
# Plots
"""

# ╔═╡ 7d28f877-b8b9-48b4-8c53-fbb7b7bdedc0
style_kwargs = (;
		  markersize=0.3, color=COLORS[3], alpha=0.3, rasterize=5)

# ╔═╡ 2ccbde9d-0475-4cca-9af1-a4cf697d2bea
kwargs_best = (; colorrange=(1, length(example_orbits)), marker=:circle, strokewidth=0)

# ╔═╡ 070a93af-b31d-4ea7-b763-4a1a36f0cf0c
function get_v_tan(df)
	icrs = ICRS(ra=df["ra"], dec=df["dec"], distance=df["distance"], 
			   pmra=df["pmra"], pmdec=df["pmdec"], radial_velocity=df["radial_velocity"])
	gsr = LilGuys.transform(GSR, icrs)
	vra = LilGuys.pm2kms(gsr.pmra, gsr.distance)
	vdec = LilGuys.pm2kms(gsr.pmdec, gsr.distance)
	return sqrt(vra^2 + vdec^2) * gsr.distance
end

# ╔═╡ eb9a4d6f-ba32-426e-8e76-2537d62ee25c
v_tan = get_v_tan.(eachrow(orbit_props))

# ╔═╡ c538b1ac-55d5-4dd5-ac1a-c812b096b33b
labels = OrderedDict(
	"peri_1.5" => "1.5",
	"peri_4" => "4",
	"mean" => "7 (mean)",
	"peri_12" => "12",
	"peri_18" => "18",
	"peri_26" => "26",
)

# ╔═╡ 0ba65849-4670-4fd7-b9a2-bdb512043fb4
[COLORS[4], COLORS[2], :black, COLORS[1], COLORS[3], COLORS[8]]

# ╔═╡ 2b9e8530-c4ad-4c46-8f1c-1413614c218a
let
	fig = Figure()
	ax = Axis(
		fig[1,1],
		xlabel = "pericentre / kpc",
		ylabel="density",
	)

	hist!(orbit_props.pericentre, bins=100, normalization=:pdf, )
		
	for (i, (peri, orbit)) in enumerate(example_orbits)
		scatter!(orbit["pericentre"], -0.005, color=i; kwargs_best...)
	end

	@savefig "peri_hist"
	fig
end

# ╔═╡ ca87bd37-b348-4535-a1ae-b016c59620a7
L = radii([orbit_props.Lx_f orbit_props.Ly_f orbit_props.Lz_f]')

# ╔═╡ 9340b2c1-bf46-4959-8d48-9117926bbc1f
let
	fig = Figure()
	ax = Axis(
		fig[1,1],
		xlabel = "helcen distance / kpc",
		ylabel = "pericentre / kpc"
	)

	scatter!(orbit_props.distance, orbit_props.pericentre; style_kwargs...)

	hidexdecorations!(ticks=false, minorticks=false)
		
	for (i, (peri, orbit)) in enumerate(example_orbits)
		scatter!(orbit["distance"], orbit["pericentre"], color=i; kwargs_best...)
	end

	
	ax_L = Axis(
		fig[2,1],
		xlabel = "helcen distance / kpc",
		ylabel = L"$|\mathbf{L}|$ / kpc km\,s$^{-1}$"
	)

	scatter!(orbit_props.distance, L * V2KMS; style_kwargs...)
	for (i, (peri, orbit)) in enumerate(example_orbits)
		
		scatter!(orbit["distance"], radii([orbit["Lx_f"], orbit["Ly_f"], orbit["Lz_f"]])*V2KMS,
				 color=i; label=labels[peri], kwargs_best...)
	end

	axislegend(position=:lt, nbanks=2)

	linkxaxes!(ax, ax_L)
	rowgap!(fig.layout, 1, 0)

	@savefig "peri_vs_dist"
	
	fig

end

# ╔═╡ Cell order:
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═48257649-efa8-448d-9296-d49517137277
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
# ╠═4da40bb0-18d7-4e79-b031-1319a3296207
# ╠═3f3df4af-3f2b-4803-b007-4302a029f988
# ╠═3575de4d-797d-4822-a28a-a4424281c277
# ╠═0700fa00-3b38-476d-be4c-c5c61037ee3e
# ╟─b5ae57dc-c148-4dca-983e-a600e0c46f57
# ╠═7d28f877-b8b9-48b4-8c53-fbb7b7bdedc0
# ╠═2ccbde9d-0475-4cca-9af1-a4cf697d2bea
# ╠═070a93af-b31d-4ea7-b763-4a1a36f0cf0c
# ╠═eb9a4d6f-ba32-426e-8e76-2537d62ee25c
# ╠═c538b1ac-55d5-4dd5-ac1a-c812b096b33b
# ╠═0ba65849-4670-4fd7-b9a2-bdb512043fb4
# ╠═9340b2c1-bf46-4959-8d48-9117926bbc1f
# ╠═2b9e8530-c4ad-4c46-8f1c-1413614c218a
# ╠═ca87bd37-b348-4535-a1ae-b016c59620a7
