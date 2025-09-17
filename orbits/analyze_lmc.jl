### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 6fe44ded-579b-4d7c-8f5f-0cc5382b59a3
module PlotUtils 
	include("plot_utils.jl")
end

# ╔═╡ 5ca2096b-6eb9-4325-9c74-421f3e0fdea2
module OrbitUtils
	include("orbit_utils.jl")
end

# ╔═╡ c863152e-b67d-4e8a-a09f-29c60fc24ebe
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 425ca5e9-364f-437c-9986-3a09eb60affc
@bind inputs confirm(notebook_inputs(;
	modelname = TextField(24, default="example"),
	galaxyname = TextField(24, default="sculptor"),
	t_min = NumberField(-10.0:10, -10.0),
	Nmax = NumberField(0:1e6, 1e4),
	random_examples = CheckBox(),
))

# ╔═╡ c9eae5d6-ea64-4a81-80fb-950a33913824
modelname = inputs.modelname

# ╔═╡ 94f73107-6934-4e13-828e-b289bb0190ba
galaxyname = inputs.galaxyname

# ╔═╡ d6fab584-059b-41df-9981-530865ec48ae
Nmax = round(Int, inputs.Nmax)

# ╔═╡ 1ee918a2-0189-412c-af9a-0f435a9ecff9
t_min = inputs.t_min

# ╔═╡ d1409cf5-6faa-40a6-abf6-fa62d5fa6839
random_examples = inputs.random_examples

# ╔═╡ f311c4d6-88a4-48a4-a0a0-8a6a6013897c
agama_units = OrbitUtils.get_units(joinpath(galaxyname, modelname))

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
FIGDIR = "./$galaxyname/$modelname/figures"

# ╔═╡ b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
p_value =cdf(Normal(), -3)# 3sigma

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
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))

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

# ╔═╡ eff58c52-a32b-4faa-9b98-c8234d9b21fc
err = lguys.ICRS(;(prop => LilGuys.get_uncertainty(obs_props, string(prop)) for prop in propertynames(obs) if prop != :coord)...)

# ╔═╡ da6e5566-f2df-4feb-9188-53eca9a1a0d5
df_props = read_fits(joinpath(galaxyname, modelname, "orbital_properties.fits"))

# ╔═╡ 74469ce6-ae1f-4eea-bc02-f1d5b73648fd
df_props_special2 = read_fits(joinpath(galaxyname, modelname * "_special_cases", "orbital_properties.fits"))

# ╔═╡ 384be6a6-f9d9-47e0-9792-aef6689dcbdb
apos = df_props.apocentre

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_props.pericentre

# ╔═╡ 392315ee-72d2-4a14-9afc-5fd6424b3e83
md"""
## quantile properties & orbit selection
"""

# ╔═╡ 950e0210-e4fe-4bad-b82f-69247fd0edd8
md"""
We need to select orbits systematically from the probability cloud. 
Asya's method was to simply take the orbits from the MCMC samples with pericentres
closest to specified quantiles. 
I use this as an initial test, but prefer to take the 
median values in some percentile rage, such that the orbit tends to be closer to the middle of the cloud.
"""

# ╔═╡ 4a4b8b73-92c3-4e1f-93d8-e44369b8f148
quantiles = [0.5, p_value, 1-p_value]

# ╔═╡ 413d4e5d-c9cd-4aca-be1e-d132b2bd616d
peri_qs = lguys.quantile(peris, quantiles)

# ╔═╡ 3ddc8254-aebd-4fe1-a6cb-82e2a45fffc3
peri_lmc_qs = lguys.quantile(df_props.pericentre_lmc, quantiles)

# ╔═╡ 17a63cc8-84f4-4248-a7b0-c8378454b1f7
if random_examples
	idx = [argmin(abs.(p .- peris)) for p in peri_qs]
end

# ╔═╡ 61c5e886-4c54-4080-8111-122765405ffe
pot = Agama.Potential(file=joinpath(galaxyname, modelname, "agama_potential.ini"))

# ╔═╡ cbd5088d-0d00-488c-b174-1ee3b177d736
if random_examples
	icrs0 = LilGuys.coords_from_df(df_props[idx, :])
	
	orbits = LilGuys.agama_orbit(pot, icrs0, agama_units=agama_units, timerange=(0, -10/T2GYR))
	orbit_labels = ["mean", "smallperi", "largeperi"]
		
else
	@assert isdir(joinpath(galaxyname, modelname * "_special_cases"))
	orbit_ics = TOML.parsefile(joinpath(galaxyname, modelname * "_special_cases", "initial_conditions.toml"))
	orbit_labels = [o["name"] for o in orbit_ics["orbits"]]
	
	orbits = [Orbit(joinpath(galaxyname, modelname * "_special_cases", "orbit_" * name * ".csv")) for name in orbit_labels]

	icrs0 = [LilGuys.ICRS(o) for o in orbit_ics["orbits"]]
	#orbits = LilGuys.agama_orbit(pot, icrs0, agama_units=agama_units, timerange=(0, -10/T2GYR))
end

# ╔═╡ 2d7b23d0-600e-439f-b547-91df46802252
lmc_orbit = OrbitUtils.get_lmc_orbit(joinpath(galaxyname, modelname))

# ╔═╡ 4132f73b-e037-45f9-b94e-0e0dfe0be3c6
begin 
	props_special = hcat(OrbitUtils.orbital_properties(pot, reverse.(orbits); agama_units=agama_units), LilGuys.to_frame(icrs0), DataFrame("label" => orbit_labels))

	df_special_lmc = OrbitUtils.orbital_properties(pot, [
		reverse(orbit - LilGuys.resample(lmc_orbit, orbit.times)) for orbit in orbits]
												   , agama_units=agama_units)


	for key in [:pericentre, :apocentre, :time_last_peri]
		props_special[!, Symbol("$(key)_lmc")] = df_special_lmc[!, key]
	end
end

# ╔═╡ e4dac06c-aadc-475d-8066-31ced204b6d0
df_special_lmc.pericentre

# ╔═╡ 2b344d65-a142-4b83-9f46-0203367935ee
df_props.x

# ╔═╡ 1acef60e-60d6-47ba-85fd-f9780934788b
md"""
# plots
"""

# ╔═╡ 50baf5a6-fb5b-494e-95f3-53414a9f1cc0
md"""
## Histograms
I have some histograms of orbital properties below, to get a sense of the overal distribution.
"""

# ╔═╡ 049ef8a5-fe4d-4c18-95d0-a361e1abdf30
PlotUtils.plot_param_hist(df_props, props_special, :pericentre)

# ╔═╡ ca1c236e-795a-408b-845b-9c13bc838619
PlotUtils.plot_param_hist(df_props, props_special, :apocentre)

# ╔═╡ d54bfa19-850d-4061-a7f6-ec94f19cfd73
PlotUtils.plot_param_hist(df_props, props_special, :pericentre_lmc)

# ╔═╡ 46b4242b-8af7-4233-8ecf-d86740b4c884
PlotUtils.plot_param_hist(df_props, props_special, :apocentre_lmc)

# ╔═╡ 820b1ae7-9555-4ddc-9881-97c0e895d57b
PlotUtils.plot_param_hist(df_props, props_special, :time_last_peri)

# ╔═╡ dedb723f-db84-4f81-9ccf-159717e5e632
PlotUtils.plot_param_hist(df_props, props_special, :time_last_peri_lmc)

# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = df_props.distance

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
hist(lguys.radii(hcat(df_props.x, df_props.y, df_props.z)'),
	axis=(; xlabel="initial galactocentric radius / kpc",
	ylabel="count")
)

# ╔═╡ 68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
hist(lguys.radii(hcat(df_props.v_x, df_props.v_y, df_props.v_z)'),
	axis=(; xlabel="initial galactocentric velocity / km/s",
		ylabel="count"
	)
)

# ╔═╡ c48b4e73-480e-4a50-b5fc-db5f6c5b040e
PlotUtils.plot_correlations(df_props, props_special, :pericentre)

# ╔═╡ 8fed6eac-307e-4db0-8bc2-1d6b2391b5fc
PlotUtils.plot_correlations(df_props, props_special, :apocentre)

# ╔═╡ fcb9bd66-0c79-44a5-851c-b2daec357f42
PlotUtils.plot_correlations(df_props, props_special, :pericentre_lmc)

# ╔═╡ e59ceb7a-b8f5-42fc-97f0-d822701caa93
PlotUtils.plot_correlations(df_props, props_special, :apocentre_lmc)

# ╔═╡ 16b8d7f7-8cea-47f0-a653-bb39a4ef1b44
PlotUtils.plot_correlations(df_props, props_special, :time_last_peri)

# ╔═╡ ae93e5f6-bada-447e-b934-fe1eca5a2de3
PlotUtils.plot_correlations(df_props, props_special, :time_last_peri_lmc)

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

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

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="galaxy–MW distance / kpc"
	)

	for i in eachindex(orbits)
		lines!(orbits[i].times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		# hlines!([peris[idx[i]], apos[idx[i]]], linestyle=:solid, alpha=0.1, color=:black)
		scatter!(props_special.time_last_peri[i] .* lguys.T2GYR, props_special.pericentre[i], color=COLORS[i])

	end

	scatter!([NaN], [NaN], color=:black, label="pericentre")
	ylims!(0, nothing)
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 3def754c-46a3-43da-9b75-9ef5540a6022
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="galaxy–LMC distance / kpc"
	)

	for i in eachindex(orbits)
		lines!(orbits[i].times * lguys.T2GYR, radii(orbits[i] - LilGuys.resample(lmc_orbit, orbits[i].times)), label=orbit_labels[i])
	
		scatter!(props_special.time_last_peri_lmc[i] .* lguys.T2GYR, props_special.pericentre_lmc[i], color=COLORS[i])

	end

	scatter!([NaN], [NaN], color=:black, label="pericentre")
		ylims!(0, nothing)

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 02d66dc9-ac91-4512-afd4-b665abed7714


# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.plot_xyz(positions..., labels=orbit_labels)

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
	
		plot!(R, z, 
			
		)
	end

	fig
end

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities..., units=" / km / s", labels=orbit_labels)

# ╔═╡ 17522ba0-7e24-4dbe-9659-8f694defaaf9
orbits[3].times

# ╔═╡ c4a1a691-51bb-4c3e-89ab-398841b1d155
md"""
# Orbit Info
"""

# ╔═╡ 519a88f0-8e2d-4c09-83e0-3cc2ee147e35
function get_initial_t(j)
	i = idx[j]
	N = length(rs[j])
	t_ini = -1
	for t in N-2:-1:1
		if diff(rs[j])[t] >= 0 && diff(rs[j])[t+1] <= 0
			t_ini = t
			break
		end
	end
	if t_ini == -1
		t_ini = length(rs[j])
	end
	return t_ini
end
	

# ╔═╡ 025fb533-a0d6-4697-bdc7-821bfcf92153
orbits

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(orbits)
		t =1
		@printf "orbit: \t\t %i\n" i
		
		@printf "pericentre:\t %0.1f\n" props_special.pericentre[i]
		@printf "apocentre: \t %0.1f\n" props_special.apocentre[i]

		@printf "time of first apocentre: %f \n" times[end] - times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -1* velocities[i][:, t]...
		@printf "final position: [%f, %f, %f]\n" positions[i][:, 1]...
		@printf "final velocity: [%f, %f, %f]\n" -lguys.V2KMS * velocities[i][:, 1]...

		o = icrs0[i]
		@printf "ra: \t %f\n" o.ra
		@printf "dec: \t %f\n" o.dec
		@printf "pmra: \t %f\n" o.pmra
		@printf "pmdec: \t %f\n" o.pmdec
		@printf "dist: \t %f\n" o.distance
		@printf "rv: \t %f\n" o.radial_velocity

		println()
	end
end

# ╔═╡ 1c5d6ba8-fec2-4629-bcc3-e9b1501b04c0
md"""
# LMC Bound?
"""

# ╔═╡ 44719d22-2e49-4edc-b460-8f0529b74b0b
function is_bound_to_lmc(df_props, lmc_pot, lmc_orbit)
    units = agama_units

	time_0 = orbits[1].times[1]
	lmc_centre = LilGuys.resample(lmc_orbit, [time_0])

    is_bound = fill(false, size(df_props, 1))

    for i in eachindex(is_bound)
        pos = [df_props.x_i[i], df_props.y_i[i], df_props.z_i[i]]
        vel = [df_props.v_x_i[i], df_props.v_y_i[i], df_props.v_z_i[i]] / V2KMS
        phi = Agama.potential(lmc_pot, pos, units, t=time_0)

        ke = radii(vel, lmc_centre.velocities[:, 1]) ^ 2 / 2
        is_bound[i] = ke + phi < 0
    end

    return is_bound
end

# ╔═╡ 6bef3bf2-231d-4b07-9729-75a28fe03c1a
orbits[1].times[end]

# ╔═╡ 9a43f9b6-79cc-455f-b1a5-7403fa407eb6
lmc_pot = Agama.Potential(file=joinpath(galaxyname, modelname, "potential_lmc.ini"))

# ╔═╡ 07f31190-63ce-45a6-b77a-7a6141c22d59
lmc_bound = is_bound_to_lmc(df_props, lmc_pot, lmc_orbit)

# ╔═╡ 7b771b69-9472-4f31-acf7-212a768528d9
is_bound_to_lmc(df_props_special2, lmc_pot, lmc_orbit)

# ╔═╡ 6ebb0cde-357a-48e0-8723-72fabdbf8b3a


# ╔═╡ ded55ff2-5b96-484b-86df-1874ee2b4314
mean(lmc_bound)

# ╔═╡ 1dd8afb2-90f9-4cd6-b16a-22c3835e7037
is_bound_to_lmc(df_special_lmc, lmc_pot, lmc_orbit)

# ╔═╡ 696599b3-2fc7-4acf-8fd0-adceca5e9e01
PlotUtils.plot_correlations(df_props[lmc_bound, :], props_special, :pericentre_lmc)

# ╔═╡ Cell order:
# ╟─7450144e-5464-4036-a215-b6e2cd270405
# ╟─2b9d49c6-74cc-4cce-b29e-04e94776863f
# ╠═425ca5e9-364f-437c-9986-3a09eb60affc
# ╠═c9eae5d6-ea64-4a81-80fb-950a33913824
# ╠═94f73107-6934-4e13-828e-b289bb0190ba
# ╠═d6fab584-059b-41df-9981-530865ec48ae
# ╠═1ee918a2-0189-412c-af9a-0f435a9ecff9
# ╠═d1409cf5-6faa-40a6-abf6-fa62d5fa6839
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═2b01d8f5-272e-4aa2-9825-58bb052acd10
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═fdea5667-7aa1-4d68-8d43-2742e9f5eb22
# ╠═a7111062-b025-43a9-bdb1-aee08deb60e9
# ╠═f823e80a-f6db-440f-8d25-56860618c82f
# ╠═00ba3075-c3e2-4965-acf3-00cda0ef320f
# ╠═a39d4241-ac9c-40bf-807e-9577cb2cc855
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═6fe44ded-579b-4d7c-8f5f-0cc5382b59a3
# ╠═5ca2096b-6eb9-4325-9c74-421f3e0fdea2
# ╠═c863152e-b67d-4e8a-a09f-29c60fc24ebe
# ╠═f311c4d6-88a4-48a4-a0a0-8a6a6013897c
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╠═18d6d521-1abf-4085-b9b0-6f45c0eb2feb
# ╟─8f70add4-effe-437d-a10a-4e15228f9fec
# ╟─b15fb3ac-2219-419a-854a-31a783acf891
# ╠═fa790a4d-e74f-479b-8ff6-aa2f23cb573d
# ╠═eff58c52-a32b-4faa-9b98-c8234d9b21fc
# ╠═da6e5566-f2df-4feb-9188-53eca9a1a0d5
# ╠═74469ce6-ae1f-4eea-bc02-f1d5b73648fd
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╟─392315ee-72d2-4a14-9afc-5fd6424b3e83
# ╟─950e0210-e4fe-4bad-b82f-69247fd0edd8
# ╠═4a4b8b73-92c3-4e1f-93d8-e44369b8f148
# ╠═413d4e5d-c9cd-4aca-be1e-d132b2bd616d
# ╠═3ddc8254-aebd-4fe1-a6cb-82e2a45fffc3
# ╠═e4dac06c-aadc-475d-8066-31ced204b6d0
# ╠═17a63cc8-84f4-4248-a7b0-c8378454b1f7
# ╠═cbd5088d-0d00-488c-b174-1ee3b177d736
# ╠═61c5e886-4c54-4080-8111-122765405ffe
# ╠═2d7b23d0-600e-439f-b547-91df46802252
# ╠═4132f73b-e037-45f9-b94e-0e0dfe0be3c6
# ╠═2b344d65-a142-4b83-9f46-0203367935ee
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╟─50baf5a6-fb5b-494e-95f3-53414a9f1cc0
# ╠═049ef8a5-fe4d-4c18-95d0-a361e1abdf30
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═d54bfa19-850d-4061-a7f6-ec94f19cfd73
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═820b1ae7-9555-4ddc-9881-97c0e895d57b
# ╠═dedb723f-db84-4f81-9ccf-159717e5e632
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
# ╠═c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╠═8fed6eac-307e-4db0-8bc2-1d6b2391b5fc
# ╠═fcb9bd66-0c79-44a5-851c-b2daec357f42
# ╠═e59ceb7a-b8f5-42fc-97f0-d822701caa93
# ╠═16b8d7f7-8cea-47f0-a653-bb39a4ef1b44
# ╠═ae93e5f6-bada-447e-b934-fe1eca5a2de3
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═3eeb1784-bc35-4ffe-b02f-8ea738d41ac8
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═1ce6663b-1435-4887-a6aa-7a5e9f6c5cde
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═3def754c-46a3-43da-9b75-9ef5540a6022
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═02d66dc9-ac91-4512-afd4-b665abed7714
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╟─57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═17522ba0-7e24-4dbe-9659-8f694defaaf9
# ╟─c4a1a691-51bb-4c3e-89ab-398841b1d155
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═025fb533-a0d6-4697-bdc7-821bfcf92153
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
# ╠═1c5d6ba8-fec2-4629-bcc3-e9b1501b04c0
# ╠═44719d22-2e49-4edc-b460-8f0529b74b0b
# ╠═6bef3bf2-231d-4b07-9729-75a28fe03c1a
# ╠═9a43f9b6-79cc-455f-b1a5-7403fa407eb6
# ╠═07f31190-63ce-45a6-b77a-7a6141c22d59
# ╠═7b771b69-9472-4f31-acf7-212a768528d9
# ╠═6ebb0cde-357a-48e0-8723-72fabdbf8b3a
# ╠═ded55ff2-5b96-484b-86df-1874ee2b4314
# ╠═1dd8afb2-90f9-4cd6-b16a-22c3835e7037
# ╠═696599b3-2fc7-4acf-8fd0-adceca5e9e01
