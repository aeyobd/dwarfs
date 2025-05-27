### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ ca5f69e4-31c3-11f0-2b96-d3cead0ed832
begin
	import Pkg; Pkg.activate()

	using CairoMakie 
	using Arya
	
	using LilGuys
	
end

# ╔═╡ ad687c7e-a9d9-4a80-b7d1-dc24f7be8a76
using PythonCall

# ╔═╡ aff88239-9b72-44df-b895-3d3e0b6430dc
using OrderedCollections

# ╔═╡ 26886e5c-546f-45ca-a882-2b354a962f07
galaxy = "ursa_minor"

# ╔═╡ cc2c7b51-5483-4693-bda6-40563ffb7289
modelname = "1e6_v38_r4.0"

# ╔═╡ d9849d12-43fe-4788-aa26-145ed0d7db19
orbitname = "orbit_smallperi.3"

# ╔═╡ 17c701ff-2620-42b7-b2f2-02e5cbe1900d
vasiliev_units = false

# ╔═╡ 620d4b82-d277-46d1-b01c-329a65ca3626
r_break_obs_arcmin = 25

# ╔═╡ c73efb23-ba43-43ba-952a-18ed7fec9e91
lmc = false

# ╔═╡ 1374aa1e-ee19-4a17-b317-cf0be4196dfe


# ╔═╡ f7b10565-3af7-4fa3-8dbf-a34e82b044cb
md"""
## Background
"""

# ╔═╡ bea952a4-e48a-489e-8f31-62dd8827e869
md"""
The jacobi radius is the first order approximation of the zero velocity point, where stars stop becoming bound to a dwarf galaxy in orbit around the host. While this is a rough approximation it is still useful. The expression is 

``
r_J \approx D \left(\frac{m}{3M}\right)^{1/3}
``

Where $r_J$ is the distance from the satellite to the jacobi radius (BT89; eq 7-84]., $D$ is the distance between systems, $m$ is the mass within $r_J$ of the satellite, and $M$ is the mass within $D$ of the host. Note that $r_J$ also solves

``
3\bar \rho_{\rm host}(D) = \bar\rho_{\rm sat}(r_J)
``

i.e. $r_J$ is the point where the satellite mean density is 3 times that of the host.
"""

# ╔═╡ d0471521-d248-4281-a7d3-718ded56eb83
module AgamaUtils 
	include(joinpath(ENV["DWARFS_ROOT"], "utils/agama_utils.jl"))
end

# ╔═╡ 733b8ae6-681a-4362-a6aa-10fd5457c6c6
CairoMakie.activate!(type=:png)

# ╔═╡ 9e99b247-77db-4383-acf0-52e16df3bd97
import TOML

# ╔═╡ 44334862-c17b-4ed1-8239-1b0902748fa8
import Agama

# ╔═╡ 7ca75d0b-cd00-4fae-a03b-f88f740b6743
md"""
# Data loading
"""

# ╔═╡ 3b645724-08a1-4689-99c6-2ff45c5cf676
potential_dir = ENV["DWARFS_ROOT"] * "/agama/potentials/"

# ╔═╡ 8bb2a7f2-97de-4f09-afdd-6cbebe9b86f5
begin 
	potential_name = "simulation/agama_potential.ini"
	if lmc
		potential_name = joinpath(potential_dir, "vasiliev24", "L3M11", "potential_lmc_init.ini")
	end
end

# ╔═╡ 8db8ca90-4c26-4277-8697-5ffd370a7bec
readdir(potential_dir * "vasiliev24/L3M11")

# ╔═╡ f3d4a880-300f-450a-bd57-9f9639f5966e
readdir(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxy,modelname))

# ╔═╡ fa86bf35-2f3a-4830-bd1d-dd3f40085f3c
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxy, modelname, orbitname)

# ╔═╡ b9665a69-879a-485f-a5e0-ef16d55b09c0
readdir(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxy, modelname))

# ╔═╡ 23a1da00-2941-4911-8f14-ebed939dca19
dwarf_halo = LilGuys.load_profile(joinpath(modeldir, "../halo.toml"))

# ╔═╡ 0d8ae710-df73-45f7-910b-cf127f38e9d4
mw_halo = Agama.Potential(file=joinpath(modeldir, potential_name))

# ╔═╡ 29a0a1ee-d953-487f-8be3-1f9335fb098b
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))

# ╔═╡ d2ca0476-2641-4ab8-99a5-6647eb0cdeb0
if lmc
	orbital_props_lmc = TOML.parsefile(joinpath(modeldir, "orbital_properties_lmc.toml"))
else
end

# ╔═╡ 90559295-c061-4377-a65f-294ee576d521
orbital_props = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))

# ╔═╡ 944a99f4-f0a7-4206-bbbc-a12d6154234d
md"""
# Utils
"""

# ╔═╡ bc854ffd-64df-4037-8cbf-2811d8df6644
function density_axis(gs; kwargs...)
    return Axis(gs[1,1];
        xlabel=L"log $r$ / kpc",
        ylabel = L"$\log\,\bar\rho / 10^{10}M_\odot$\,kpc$^{-3}$",
				kwargs...
        )
end


# ╔═╡ 031f5ee5-1ce9-49f0-ad3d-5df668924e43
function calc_ρ_mean(prof::LilGuys.SphericalProfile, r)
    y = LilGuys.mass.(prof, r)
    rho_mean = @. y / (4π/3*r^3)
end

# ╔═╡ 2d61739c-4330-4cce-8203-e2dd92a99e46
function calc_ρ_mean(prof::Agama.Potential, r; vasiliev_units=vasiliev_units)
    M = Agama.enclosed_mass(prof, r)
    if vasiliev_units
        M *= AgamaUtils.V_M2MSUN ./ M2MSUN
    end

    y = @. M / (4π/3 * r^3)
end

# ╔═╡ 88df87cf-5b6c-4ad6-afe5-27dd2a42a674
function plot_ρ_mean!(ax, prof::LilGuys.SphericalProfile)
    x = LinRange(-2, 1.5, 1000)
    r = 10 .^ x
    y = calc_ρ_mean(prof, r,)

    lines!(x, log10.(y))
end

# ╔═╡ c1499c52-e66f-4946-abf0-39d14d5b2a09
function plot_ρ_mean!(ax, prof::Py; xlimits=(-2, 1.5), vasiliev_units=vasiliev_units, kwargs...)
    x = LinRange(xlimits[1], xlimits[2], 1000)
    r = 10 .^ x
 
    y = calc_ρ_mean(prof, r, vasiliev_units=vasiliev_units)
    lines!(x, log10.(y); kwargs...)
end

# ╔═╡ ac9a6b4f-50db-4915-ac03-2d72a84b3f4c
md"""
# Calculations
"""

# ╔═╡ bf082ceb-e772-4091-a178-a160f980a444
if lmc 
	peri = orbital_props_lmc["pericentre"]
else
	peri = orbital_props["pericentre"]
end

# ╔═╡ 7dcbcf23-7a19-4931-b63c-5308e5abed10
ρ_host = calc_ρ_mean(mw_halo, [peri, peri])[1]

# ╔═╡ 7784c958-e815-4f08-a255-ddd7cc2af8b0
r_J = LilGuys.find_zero(r -> calc_ρ_mean(dwarf_halo, r) - 3*ρ_host, 1.)

# ╔═╡ 4e072f88-f319-4ae7-8220-3f8d5b0766b9
r_J_arcmin = LilGuys.kpc2arcmin(r_J, orbital_props["distance_f"])

# ╔═╡ 40c80707-2d4d-449f-b4af-20e691d4d43d
md"""
# inverse
Given r_J what is required pericentre?

"""

# ╔═╡ 71cc56b5-8fb5-4d79-b8f9-4ec72963fc92
function calc_rperi_rj(host, satellite; limits=(0, 2), vasiliev_units=vasiliev_units)
    x = LinRange(limits[1], limits[2], 1000)

    r = 10 .^ x

    ρ_mean_host = LilGuys.lerp(r, calc_ρ_mean(host, r, vasiliev_units=vasiliev_units))

    r2 = 10 .^ LinRange(-2, 3, 1000)
    ρ_mean_sat = LilGuys.lerp(r2, calc_ρ_mean(satellite, r2))

    rj = zeros(length(r))

    rperi = r
    
    for i in eachindex(rj)
        rr = rperi[i]
        rj[i] = LilGuys.find_zero(r -> ρ_mean_sat(r) - 3*ρ_mean_host(rr), 1.)
    end

    return rperi, rj
end

    

# ╔═╡ 80ac2f40-678d-4de4-b77d-d3780e164f05
peris, rjs = calc_rperi_rj(mw_halo, dwarf_halo)

# ╔═╡ 7e2bde96-0232-47ff-aba0-6fbeab70ea50
r_break_obs = LilGuys.arcmin2kpc(r_break_obs_arcmin, obs_props["distance"])

# ╔═╡ 00e4bfa3-cfb5-4ce6-8384-37aaf1cb80b2
r_peri_req =LilGuys.find_zero(r -> calc_ρ_mean(dwarf_halo, r_break_obs) - 3*calc_ρ_mean(mw_halo, r), 1.)

# ╔═╡ eebdd830-b077-474e-b044-4ee6e14c9b84
let 
fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log peri", ylabel = "rJ / kpc")
			
	lines!(log10.(peris), log10.(rjs), linewidth=2)


	vlines!(log10(peri), color=COLORS[2], label="current orbit")
	hlines!(log10.(r_J), color=COLORS[2])

	hlines!(log10.(r_break_obs), color=COLORS[3], label="observed break")
	vlines!(log10.(r_peri_req), color=COLORS[3])

	axislegend(position=:rb)

	fig

end

# ╔═╡ afae5f49-accd-4869-9480-a9afbe4bcbfc
let
	fig = Figure(size=(5*72, 3*72))
	ax = density_axis(fig[1,1], title="dwarf")

	x = LinRange(-1, 1, 100)

	f_dwarf(x) = log10.(calc_ρ_mean(dwarf_halo, 10 .^x))
	f_mw(x) = log10.(calc_ρ_mean(mw_halo, 10 .^x))
					 
	lines!(x, f_dwarf(x))

	for i in 1:2
		x = [log10(r_J), log10(r_break_obs)][i]
		label = ["model", "required"][i]
		scatter!(x, f_dwarf(x), label=label)
	end

	axislegend(position=:rt)

	# MW
	ax_mw = density_axis(fig[1,2], title="MW", xreversed=true)

	x = LinRange(-1, 2, 100)

	lines!(x, log10.(calc_ρ_mean(mw_halo, 10 .^ x)))
	hideydecorations!(ax_mw, ticks=false, minorticks=false)

	for x in [log10(peri), log10(r_peri_req)]
		scatter!(x, f_mw(x))
	end


	linkyaxes!(ax, ax_mw)

	fig
end

# ╔═╡ df8d55f3-1629-415b-9461-b66d8ef2c35c
md"""
# Properties
"""

# ╔═╡ b77e000c-c66e-4bd6-8716-26216e6f7bb4
df_summary = OrderedDict(
	"r_J" => r_J_arcmin,
	"r_J_kpc" => r_J,
	"peri" => peri,
	"r_break" => r_break_obs_arcmin,
	"r_break_kpc" => r_break_obs,
	"peri_required" => r_peri_req,
	"rho_mean_host_peri" => ρ_host,
	"rho_mean_host_required" => calc_ρ_mean(mw_halo, r_peri_req),
)

# ╔═╡ 6021e71d-01fd-4a2e-bc3d-af45db173a5b
if lmc
	outname = "jacobi_lmc.toml"
else
	outname = "jacobi.toml"
end

# ╔═╡ 82f87cc8-6a2a-4cb3-ae79-941d6ed16fc5
open(joinpath(modeldir, outname), "w") do f
	TOML.print(f, df_summary)
end

# ╔═╡ Cell order:
# ╠═26886e5c-546f-45ca-a882-2b354a962f07
# ╠═cc2c7b51-5483-4693-bda6-40563ffb7289
# ╠═d9849d12-43fe-4788-aa26-145ed0d7db19
# ╠═17c701ff-2620-42b7-b2f2-02e5cbe1900d
# ╠═620d4b82-d277-46d1-b01c-329a65ca3626
# ╠═c73efb23-ba43-43ba-952a-18ed7fec9e91
# ╠═8bb2a7f2-97de-4f09-afdd-6cbebe9b86f5
# ╠═1374aa1e-ee19-4a17-b317-cf0be4196dfe
# ╠═8db8ca90-4c26-4277-8697-5ffd370a7bec
# ╟─f7b10565-3af7-4fa3-8dbf-a34e82b044cb
# ╟─bea952a4-e48a-489e-8f31-62dd8827e869
# ╠═ca5f69e4-31c3-11f0-2b96-d3cead0ed832
# ╠═d0471521-d248-4281-a7d3-718ded56eb83
# ╠═733b8ae6-681a-4362-a6aa-10fd5457c6c6
# ╠═9e99b247-77db-4383-acf0-52e16df3bd97
# ╠═ad687c7e-a9d9-4a80-b7d1-dc24f7be8a76
# ╠═44334862-c17b-4ed1-8239-1b0902748fa8
# ╟─7ca75d0b-cd00-4fae-a03b-f88f740b6743
# ╠═3b645724-08a1-4689-99c6-2ff45c5cf676
# ╠═f3d4a880-300f-450a-bd57-9f9639f5966e
# ╠═fa86bf35-2f3a-4830-bd1d-dd3f40085f3c
# ╠═b9665a69-879a-485f-a5e0-ef16d55b09c0
# ╠═23a1da00-2941-4911-8f14-ebed939dca19
# ╠═0d8ae710-df73-45f7-910b-cf127f38e9d4
# ╠═29a0a1ee-d953-487f-8be3-1f9335fb098b
# ╠═d2ca0476-2641-4ab8-99a5-6647eb0cdeb0
# ╠═90559295-c061-4377-a65f-294ee576d521
# ╟─944a99f4-f0a7-4206-bbbc-a12d6154234d
# ╠═bc854ffd-64df-4037-8cbf-2811d8df6644
# ╠═031f5ee5-1ce9-49f0-ad3d-5df668924e43
# ╠═2d61739c-4330-4cce-8203-e2dd92a99e46
# ╠═88df87cf-5b6c-4ad6-afe5-27dd2a42a674
# ╠═c1499c52-e66f-4946-abf0-39d14d5b2a09
# ╟─ac9a6b4f-50db-4915-ac03-2d72a84b3f4c
# ╠═bf082ceb-e772-4091-a178-a160f980a444
# ╠═7dcbcf23-7a19-4931-b63c-5308e5abed10
# ╠═7784c958-e815-4f08-a255-ddd7cc2af8b0
# ╠═4e072f88-f319-4ae7-8220-3f8d5b0766b9
# ╟─40c80707-2d4d-449f-b4af-20e691d4d43d
# ╠═71cc56b5-8fb5-4d79-b8f9-4ec72963fc92
# ╠═80ac2f40-678d-4de4-b77d-d3780e164f05
# ╠═7e2bde96-0232-47ff-aba0-6fbeab70ea50
# ╠═00e4bfa3-cfb5-4ce6-8384-37aaf1cb80b2
# ╠═eebdd830-b077-474e-b044-4ee6e14c9b84
# ╠═afae5f49-accd-4869-9480-a9afbe4bcbfc
# ╟─df8d55f3-1629-415b-9461-b66d8ef2c35c
# ╠═aff88239-9b72-44df-b895-3d3e0b6430dc
# ╠═b77e000c-c66e-4bd6-8716-26216e6f7bb4
# ╠═6021e71d-01fd-4a2e-bc3d-af45db173a5b
# ╠═82f87cc8-6a2a-4cb3-ae79-941d6ed16fc5
