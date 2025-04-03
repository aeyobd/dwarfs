### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	import PythonCall
	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ 4c551542-92dc-4a61-b4aa-e272fd29d565
stars = LilGuys.read_fits(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/processed/rv_members.fits")

# ╔═╡ 0eac8d32-eb66-458d-a668-18492e658901
stars_good = filter(x->isfinite(x.VZ), stars)

# ╔═╡ 88025a55-f541-4aa1-8426-14e83bd790b7
CairoMakie.activate!(type=:png)

# ╔═╡ 63d5e3f6-6cec-4e50-8fa1-8907508aa47f
import TOML

# ╔═╡ b2f9ed4b-52be-4eeb-a419-d6406148e6c3
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "ursa_minor", "observed_properties.toml"))

# ╔═╡ 1f07b964-a443-4f70-abf5-35512f9d6b69
gsr = LilGuys.transform(GSR, ICRS(obs_props))

# ╔═╡ 23989a87-6aad-4ca1-bd49-f6f13977234b
v0 = gsr.radial_velocity

# ╔═╡ 87c47d8f-4639-49ca-a4e7-467f583c2e3e
δv = 5 * obs_props["sigma_v"]

# ╔═╡ 940549d0-666e-4238-bda8-ceb2fece356f
v0 + δv

# ╔═╡ 0e7ed733-f30b-4e5d-88b3-211b454125b3
extrema(stars.VZ[isfinite.(stars.VZ)])

# ╔═╡ ef53825e-879f-41e4-891e-11133ffbe80d
isfinite.(stars.VZ)

# ╔═╡ 90e396c1-8bae-475e-8638-3634bc271cdf
colorrange = (v0 - δv, v0 + δv)

# ╔═╡ 8ae4ddda-79eb-43f9-ae68-c6cb691f366f
@savefig "umi_v_z_obs_scatter" let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = L"$\xi$ / arcmin",
		ylabel = L"$\eta$ / arcmin",
		aspect = DataAspect(),
		xreversed=true,
		#limits = (-50, 50, -50, 50)
	)
			
	p = voronoiplot!(stars_good.xi, stars_good.eta, color = stars_good.VZ, colormap=:redsblues, markersize=2, strokewidth = 0, colorrange=colorrange, show_generators=false)
	
	scatter!(stars_good.xi, stars_good.eta, color = :black, markersize=1, alpha=0.5, strokewidth=0)

	Colorbar(fig[1,2], p, label = L"$v_z$ / km\,s$^{-1}$")
	fig
end

# ╔═╡ 64b5fd2b-e87e-4cbe-833d-900f6332a7ee
LilGuys.pm2kms(0.001, obs_props["distance"])

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═4c551542-92dc-4a61-b4aa-e272fd29d565
# ╠═0eac8d32-eb66-458d-a668-18492e658901
# ╠═88025a55-f541-4aa1-8426-14e83bd790b7
# ╠═63d5e3f6-6cec-4e50-8fa1-8907508aa47f
# ╠═b2f9ed4b-52be-4eeb-a419-d6406148e6c3
# ╠═1f07b964-a443-4f70-abf5-35512f9d6b69
# ╠═23989a87-6aad-4ca1-bd49-f6f13977234b
# ╠═87c47d8f-4639-49ca-a4e7-467f583c2e3e
# ╠═940549d0-666e-4238-bda8-ceb2fece356f
# ╠═0e7ed733-f30b-4e5d-88b3-211b454125b3
# ╠═ef53825e-879f-41e4-891e-11133ffbe80d
# ╠═90e396c1-8bae-475e-8638-3634bc271cdf
# ╠═8ae4ddda-79eb-43f9-ae68-c6cb691f366f
# ╠═64b5fd2b-e87e-4cbe-833d-900f6332a7ee
