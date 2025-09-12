### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 1985e943-655a-45f6-a3a8-d4125031ac8e
galaxy = "fornax"

# ╔═╡ 472e24b8-3b89-4731-9752-6333597f555b
import TOML

# ╔═╡ 20c338e3-3d10-41f4-b8ee-d0eda4e755bd
CairoMakie.activate!(type=:png)

# ╔═╡ d273f1e3-8cf7-4531-8757-6ef4735a01d8
md"""
# Data Loading
"""

# ╔═╡ 9732c616-9f53-4913-ab63-906c4fa5beee
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))

# ╔═╡ 5743880e-567f-4af5-acb5-e2ec420110aa
Mstar = 1.2*LilGuys.mag_to_L(obs_props["Mv"]) / M2MSUN

# ╔═╡ d5097843-1aad-49c3-a51a-5c580606315e
R_h =  LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])

# ╔═╡ dfa45730-fde1-4600-aa91-61a7d78ffd13
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 0e381a00-bb4d-4c88-b4c3-cbd6475c76df
R_s = R_h/ α_exp

# ╔═╡ a6fd5c3f-c9ae-40a0-a8b7-c09b9dd9d2af
prof_stars_ana = LilGuys.Exp2D(R_s=R_s, M=Mstar)

# ╔═╡ ce86303c-5c04-47e1-bf7d-c87edfdb7057
v_max = LilGuys.vel_from_M_s_fattahi(Mstar)

# ╔═╡ d6e0a0c7-6034-40ac-a20e-e70e889b021f
r_max = LilGuys.Ludlow.solve_rmax(v_max)

# ╔═╡ 932473fc-535f-4dab-bb03-53c8a71aff16
prof_dm_ana = NFW(r_circ_max=r_max, v_circ_max=v_max)

# ╔═╡ 1d6de9c4-c157-4a80-bf89-f2752419fe54
v_max * V2KMS

# ╔═╡ ee96c7dd-56a3-47c7-b6d8-fdff3f598f3d
LilGuys.mass(prof_stars_ana)

# ╔═╡ 36e17f1c-4372-4718-8d6c-ea721240f9f9
Mstar * M2MSUN / 1e7

# ╔═╡ 896de8b7-1023-499f-9162-f40e5ab39afe
mass(prof_dm_ana, 12)

# ╔═╡ b04db283-a222-4cc8-ad42-f1e1c60cef16
LilGuys.M200(prof_dm_ana), LilGuys.concentration(prof_dm_ana)

# ╔═╡ 899e958f-2ec1-4886-8fb5-053671f948b6
R_s

# ╔═╡ c1d328a6-a05f-4231-9494-58502c0f7942
@savefig "example_density_profiles" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "log radius / kpc",
		ylabel = L"log $\rho$ / $10^{10}\,\textrm{M}_\odot\,\textrm{kpc}^{-3}$",
		limits = (-1, 2, -10, 0.0)
	)

	x = LinRange(-2, 2, 1000)
	r = 10 .^ x
	
	lines!(x, log10.(LilGuys.density.(prof_dm_ana, r)), label="NFW dark matter", linestyle=:solid, color=COLORS[1])
	lines!(x, log10.(LilGuys.density.(prof_stars_ana, r)), label="exponential profile", color=COLORS[2], linestyle=:dot)


	axislegend(margin=(19,19,19,19), position=:rt)
	fig
end

# ╔═╡ Cell order:
# ╠═1985e943-655a-45f6-a3a8-d4125031ac8e
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═472e24b8-3b89-4731-9752-6333597f555b
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╟─d273f1e3-8cf7-4531-8757-6ef4735a01d8
# ╠═9732c616-9f53-4913-ab63-906c4fa5beee
# ╠═5743880e-567f-4af5-acb5-e2ec420110aa
# ╠═d5097843-1aad-49c3-a51a-5c580606315e
# ╠═0e381a00-bb4d-4c88-b4c3-cbd6475c76df
# ╠═dfa45730-fde1-4600-aa91-61a7d78ffd13
# ╠═a6fd5c3f-c9ae-40a0-a8b7-c09b9dd9d2af
# ╠═ce86303c-5c04-47e1-bf7d-c87edfdb7057
# ╠═d6e0a0c7-6034-40ac-a20e-e70e889b021f
# ╠═932473fc-535f-4dab-bb03-53c8a71aff16
# ╠═1d6de9c4-c157-4a80-bf89-f2752419fe54
# ╠═ee96c7dd-56a3-47c7-b6d8-fdff3f598f3d
# ╠═36e17f1c-4372-4718-8d6c-ea721240f9f9
# ╠═896de8b7-1023-499f-9162-f40e5ab39afe
# ╠═b04db283-a222-4cc8-ad42-f1e1c60cef16
# ╠═899e958f-2ec1-4886-8fb5-053671f948b6
# ╠═c1d328a6-a05f-4231-9494-58502c0f7942
