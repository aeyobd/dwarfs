### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

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


# ╔═╡ fdea5667-7aa1-4d68-8d43-2742e9f5eb22
using PlutoUI

# ╔═╡ f823e80a-f6db-440f-8d25-56860618c82f
using LilGuys

# ╔═╡ a39d4241-ac9c-40bf-807e-9577cb2cc855
using PyFITS

# ╔═╡ d975d00c-fd69-4dd0-90d4-c4cbe73d9754
using Statistics, Distributions

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

# ╔═╡ 6ca3fe17-3f13-43fe-967b-881078135ead
@bind modelname TextField(24, default="example") |> confirm

# ╔═╡ 5ca2096b-6eb9-4325-9c74-421f3e0fdea2
module OrbitUtils
	include("orbit_utils.jl")
end

# ╔═╡ f311c4d6-88a4-48a4-a0a0-8a6a6013897c
agama_units = OrbitUtils.get_units(modelname)

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
FIGDIR = "./$modelname/figures"

# ╔═╡ c759d1cb-a412-4587-8aa3-f39e15c09ee2
-2/T2GYR

# ╔═╡ 7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
md"""
## Setup
"""

# ╔═╡ 2b01d8f5-272e-4aa2-9825-58bb052acd10
import Agama

# ╔═╡ a7111062-b025-43a9-bdb1-aee08deb60e9
CairoMakie.activate!(type=:png)

# ╔═╡ 00ba3075-c3e2-4965-acf3-00cda0ef320f
import TOML

# ╔═╡ 3b83205d-91c1-481e-9305-0d59bc692135
coord_labels = Dict(
	:ra => "ra / degrees",
	:dec => "dec / degrees",
	:pmra => L"$\mu_{\alpha *}$ / mas\,yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas\,yr$^{-1}$",
	:radial_velocity => L"$v_\textrm{los}$ / km\,s$^{-1}$",
	:distance => "distance / kpc",
)

# ╔═╡ 8b818798-69fb-481d-ade1-9fd436b1f281
kms_label = L" / km\,s$^{-1}$"

# ╔═╡ 18d6d521-1abf-4085-b9b0-6f45c0eb2feb
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ 8f70add4-effe-437d-a10a-4e15228f9fec
md"""
# Data loading
"""

# ╔═╡ b15fb3ac-2219-419a-854a-31a783acf891
md"""
The code below lets us import the obs variable from the sample.jl file in the simulation directory.
This is primarily used for tests later.
"""

# ╔═╡ fa790a4d-e74f-479b-8ff6-aa2f23cb573d
obs = lguys.ICRS(obs_props)

# ╔═╡ 88536e86-cf2a-4dff-ae64-514821957d40
md"""
warning, the simulation output and peris_apos.fits must analyze the same simulation. Checking the individual orbits below should confirm this and the loading code should fail if there are any issues.
"""

# ╔═╡ 26d616da-95ec-4fb9-b9a8-2f095d74c722
"""
	sort_snap(snap)

returns the snapshot sorted by the index
"""
function sort_snap(snap)
	return snap[sortperm(snap.index)]
end

# ╔═╡ da6e5566-f2df-4feb-9188-53eca9a1a0d5
df_props = read_fits(joinpath(modelname, "orbital_properties.fits"))

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ 35ce583b-0938-429e-af5d-b17b399f6690
Nmax = 100

# ╔═╡ 61c5e886-4c54-4080-8111-122765405ffe
pot = Agama.Potential(file=joinpath(modelname, "agama_potential.ini"))

# ╔═╡ 6fcf3685-3a60-4997-aa41-cc4cf4891797
#orbits = LilGuys.agama_orbit(pot, icrs0, agama_units=agama_units, timerange=(0, -10/T2GYR))

# ╔═╡ cc852a14-63de-4094-821b-b5ed81fd9b7e
idx, orbits = let
	structs = LilGuys.read_ordered_structs(joinpath(modelname, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	first.(structs)[filt], last.(structs)[filt]
end

# ╔═╡ 7cfec8fe-2b5a-4939-a04c-a0cdcc0292ca
icrs0 = LilGuys.coords_from_df(df_props[idx, :])

# ╔═╡ c4bd3b49-cb3b-45c3-8122-d623490e593a
Base.summarysize(orbits) / 2^30

# ╔═╡ 3eeb1784-bc35-4ffe-b02f-8ea738d41ac8
positions = LilGuys.positions.(orbits)

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
velocities = LilGuys.velocities.(orbits)

# ╔═╡ 1ce6663b-1435-4887-a6aa-7a5e9f6c5cde
times = orbits[1].times

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
begin
	rs = lguys.radii.(positions)
	vs = lguys.radii.(velocities)
end

# ╔═╡ 5ec0129c-3075-44f1-bcdf-7090484bcd8d
md"""
# Plots
"""

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="Galactocentric distance / kpc"
	)

	for i in eachindex(orbits)
		lines!(orbits[i].times * lguys.T2GYR, rs[i], alpha=0.1, color=COLORS[1])
	end

	@savefig "r_vs_time_samples"
	fig
end

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
@savefig "xyz_samples" lguys.plot_xyz(positions..., color=COLORS[1], alpha=0.1)

# ╔═╡ 57a8d1c8-3940-4430-8b46-375fb2bf1695
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "R / kpc",
		ylabel = "z / kpc",
		aspect = DataAspect(),
		limits = (0, nothing, nothing, nothing)
	)
	for i in eachindex(positions)
		x = positions[i][1, :]
		y = positions[i][2, :]
		z = positions[i][3, :]
		R = @. sqrt(x^2 + y^2)
	
		lines!(R, z, 
			alpha=0.1, color=COLORS[1]
		)
	end

	fig
end

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities..., units=" / km s⁻¹", alpha=0.1, color=COLORS[1])

# ╔═╡ ec00c319-284f-4b83-a670-fb49dd6dedc4
Φs_ext = [Agama.potential(pot, orbit.positions, agama_units) for orbit in orbits]

# ╔═╡ 0c90dc59-6641-4094-b283-fcef68271019
md"""
The two plots below show the 3D cloud of initial positions and velocities for visualization purposes.
"""

# ╔═╡ f6b27164-ee7c-428b-aefb-75e89d178f3e
let
	fig = Figure()
	ax = Axis3(fig[1, 1], 
		xlabel = L"$x$ / kpc",
		ylabel = L"$y$ / kpc",
		zlabel = L"$z$ / kpc",
	)
	
	scatter!(df_props.x, df_props.y, df_props.z,
)
	fig
end

# ╔═╡ 5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
let
	fig = Figure()
	ax = Axis3(fig[1, 1], 
		xlabel = L"$v_x$ %$kms_label",
		ylabel = L"$v_y$ %$kms_label",
		zlabel = L"$v_z$ %$kms_label",
	)
	
	scatter!(df_props.v_x * V2KMS, df_props.v_y * V2KMS, df_props.v_z * V2KMS,
)
	fig
end

# ╔═╡ c0151e67-0ca4-4f0c-8fbd-b397c4ff7de8
tides = OrbitUtils.scalar_tidal_forces.(pot, positions, agama_units=agama_units)

# ╔═╡ 0e2e3c40-7aca-4fd0-b0c8-0850a4e1c8b5
max_tides = maximum.(tides)

# ╔═╡ 809664c4-9a85-4f87-ab9e-eec4178f6f80
scatter(max_tides, df_props.pericentre[idx])

# ╔═╡ f0b92dbf-2226-4500-a232-f8f178e46b2e
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "galactocentric radius",
			  ylabel = "tidal strain"
			 )

	for i in eachindex(orbits)
		lines!(rs[i], tides[i], color=COLORS[1], alpha=0.1)
	end

	@savefig "tides_vs_radius"

	fig
end

# ╔═╡ 16dc40c7-dd4a-4490-ac0f-b2c36c0e984e
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "time / Gyr",
			 ylabel = "tidal strain")


	for i in idx
		lines!(times * T2GYR, tides[i], color=COLORS[1], alpha=0.1)
	end


	@savefig "tides_vs_time"
	fig
end

# ╔═╡ Cell order:
# ╟─7450144e-5464-4036-a215-b6e2cd270405
# ╟─2b9d49c6-74cc-4cce-b29e-04e94776863f
# ╠═6ca3fe17-3f13-43fe-967b-881078135ead
# ╠═5ca2096b-6eb9-4325-9c74-421f3e0fdea2
# ╠═f311c4d6-88a4-48a4-a0a0-8a6a6013897c
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═c759d1cb-a412-4587-8aa3-f39e15c09ee2
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═2b01d8f5-272e-4aa2-9825-58bb052acd10
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═fdea5667-7aa1-4d68-8d43-2742e9f5eb22
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═a39d4241-ac9c-40bf-807e-9577cb2cc855
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╠═18d6d521-1abf-4085-b9b0-6f45c0eb2feb
# ╟─8f70add4-effe-437d-a10a-4e15228f9fec
# ╠═b15fb3ac-2219-419a-854a-31a783acf891
# ╠═fa790a4d-e74f-479b-8ff6-aa2f23cb573d
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╟─26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╠═da6e5566-f2df-4feb-9188-53eca9a1a0d5
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═35ce583b-0938-429e-af5d-b17b399f6690
# ╠═7cfec8fe-2b5a-4939-a04c-a0cdcc0292ca
# ╠═61c5e886-4c54-4080-8111-122765405ffe
# ╠═6fcf3685-3a60-4997-aa41-cc4cf4891797
# ╠═cc852a14-63de-4094-821b-b5ed81fd9b7e
# ╠═c4bd3b49-cb3b-45c3-8122-d623490e593a
# ╠═3eeb1784-bc35-4ffe-b02f-8ea738d41ac8
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═1ce6663b-1435-4887-a6aa-7a5e9f6c5cde
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╟─5ec0129c-3075-44f1-bcdf-7090484bcd8d
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═ec00c319-284f-4b83-a670-fb49dd6dedc4
# ╟─0c90dc59-6641-4094-b283-fcef68271019
# ╟─f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╟─5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╠═c0151e67-0ca4-4f0c-8fbd-b397c4ff7de8
# ╠═0e2e3c40-7aca-4fd0-b0c8-0850a4e1c8b5
# ╠═809664c4-9a85-4f87-ab9e-eec4178f6f80
# ╠═f0b92dbf-2226-4500-a232-f8f178e46b2e
# ╠═16dc40c7-dd4a-4490-ac0f-b2c36c0e984e
