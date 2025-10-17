### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

	using PyFITS
end

# ╔═╡ 1482481d-a5f3-48c2-a4d2-1353afe7fd72
using OrderedCollections

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 043fd044-e44c-4c9f-861e-22db402e1b10
import Agama

# ╔═╡ d992c565-3786-44dc-835f-9e6d7cf9a094
import TOML

# ╔═╡ 837df33b-02ee-4312-a1f8-db60b0baeba1
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ aeb0721a-f9e7-4a52-bc76-6fe0f0282f91
module Utils
	include("./utils.jl")
end

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ dcf00943-acab-4805-ad19-8a9d6b80f32b
md"""
# Utils
"""

# ╔═╡ 19a931d3-490b-4450-a1ed-b43ca337ba43
pot = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 8081e889-c1d4-4e50-b8e1-91aa8e82adf1
function modeldir(galaxyname, modelname)
	joinpath(ENV["DWARFS_ROOT"], "analysis",galaxyname, modelname,)
end

# ╔═╡ 6161cd69-594b-4eb9-83d9-12a682db7c88
function get_scalars(galaxyname, modelname)
	read_fits(joinpath(modeldir(galaxyname, modelname), "profiles_scalars.fits")
	)
end

# ╔═╡ 62ec2189-8fdd-4e7e-b679-05aade350dd0
function get_idx_f(galaxyname, modelname)
	return TOML.parsefile(joinpath(modeldir(galaxyname, modelname), "orbital_properties.toml"))["idx_f"]
end

# ╔═╡ 882e0efa-d365-40f2-a609-d2840b6a66ac
function get_peri(galaxyname, modelname, starsname="")
	return TOML.parsefile(joinpath(modeldir(galaxyname, modelname), "orbital_properties.toml"))["pericentre"]
end

# ╔═╡ 74857d22-6166-4c0c-88f0-961f42cc159f
function get_r_j(galaxyname, modelname, lmc=false)
	model_dir = modeldir(galaxyname, modelname)
	if lmc
		props = TOML.parsefile(model_dir * "/jacobi_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "/jacobi.toml")
	end

	return props["r_J_kpc"]
end

# ╔═╡ 2c799458-b133-4f34-ac59-e2f83dba5c68
function mass_profile(filename)
	prof_kwargs = TOML.parsefile(filename)

	prof_kwargs = LilGuys.collapse_errors(prof_kwargs) |> LilGuys.dict_to_tuple

	return LilGuys.MassProfile(;prof_kwargs...)
end

# ╔═╡ aa77e11f-65a0-40cd-af1c-43c43fb555ef
function get_stellar_profiles(galaxyname, modelname, starsname; Mstar=1)
	stars_dir = joinpath(modeldir(galaxyname, modelname), "stars", starsname)
	prof_i = mass_profile(joinpath(stars_dir, "stellar_profile_i_3d_mass.toml"))
	prof_f = mass_profile(joinpath(stars_dir, "stellar_profile_f_3d_mass.toml"))

	prof_i = LilGuys.scale(prof_i, 1, Mstar / prof_i.M_in[end].middle)
	prof_f = LilGuys.scale(prof_f, 1, Mstar / prof_f.M_in[end].middle)

	return prof_i, prof_f
end

# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(galaxyname, modelname, starsname)
	mass_profs = LilGuys.read_ordered_structs(
		joinpath(modeldir(galaxyname, modelname), "profiles.hdf5"),
		LilGuys.MassProfile
	)

	prof_stars_i, prof_stars_f = get_stellar_profiles(galaxyname, modelname, starsname, Mstar=Utils.get_M_star(galaxyname))

	idx_f = get_idx_f(galaxyname, modelname)
	
	return (dm_i=mass_profs[1].second, 
			dm_f=mass_profs[idx_f].second,
			stars_i = prof_stars_i,
			stars_f = prof_stars_f,
			scalars = get_scalars(galaxyname, modelname)
		   )
end

# ╔═╡ 2971bd21-552f-4dc4-a863-26a20e2b2cd8
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ 50a763c6-7db8-4e5e-a4e7-a3aceb010b6b
function i_f_legend(ax, profs)
	t_i = round(-profs.dm_f.time * T2GYR, digits=1)
	l1 = lines!([NaN], [NaN], linestyle=:dot, color=:black)
	l2 = lines!([NaN], [NaN], linestyle=:solid, color=:black)
	l3 = scatter!([NaN], [NaN], color=COLORS[1])
	axislegend(ax, [l1, l2, l3], [L"initial ($t=%$t_i$ Gyr)", L"final ($t=0$ Gyr)", "velocity max", ], position=:rb)
end

# ╔═╡ 07fe2997-559e-4bef-b4f4-934d850cdfe2
function plot_r_j!(r_j, y0)	
	y_mark = y0 * 10^-0.5
	arrows2d!([r_j], [y_mark], [0], [y_mark * (1 - 10^-0.30)], shaftwidth=theme(:linewidth)[]/2, tiplength=9, 
			  color=:grey, minshaftlength=0)
	text!(r_j, y_mark, text="Jacobi", rotation=-π/2, align=(:left, :center), color=:grey)

end

# ╔═╡ 59e9658c-7a2c-4964-8a24-fdf7b152ebcc
theme(:Annotation)[:style][]

# ╔═╡ 8b41afc4-a86a-4d1f-937d-3ff7a0e9cdd7
1/ sin(90)

# ╔═╡ 51c4cca9-3cc9-4242-9d41-53b37b3d8f04
function plot_prof(prof; text=nothing, color=nothing, kwargs...)
	x = radii(prof)
	y = middle.(LilGuys.circular_velocity(prof) * V2KMS)
	lines!(x, y; color=color, kwargs...)

	if !isnothing(text)
		Utils.text_along_line_log!(x, y, 0.05, text=text, color=color, align=(:left, :bottom))
	end

end

# ╔═╡ 59fee278-e9be-4733-9356-3d3bf9c08310
function plot_i_f(profs; kwargs...)
	plot_prof(profs.dm_i, text="dark matter", linestyle=:dot; kwargs...)
	plot_prof(profs.dm_f; kwargs...)
	
end

# ╔═╡ 20dbd304-e705-4dde-8ea4-558151e63f53
function plot_i_f_stars(profs; kwargs...)
	prof_i, prof_f = profs
	plot_prof(profs.stars_i, text="stars", linestyle=:dot; kwargs...)
	plot_prof(profs.stars_f; kwargs...)

end

# ╔═╡ 7f039e64-fbff-4091-924c-26523911074c


# ╔═╡ 70090f1a-0f6a-4082-a735-3b8168e26855
function plot_tidal_track!(df)
	scatter!(df.r_circ_max, df.v_circ_max * V2KMS, markersize=2/3 * theme(:markersize)[], color=COLORS[1])
end

# ╔═╡ 3d55cd39-5727-48f8-a8ae-de168e9e08c8
function get_mean_density(args...)
	peris = get_peri(args...)
	ρ =  Agama.enclosed_mass(pot, peris) ./ (4π/3 * peris .^ 3)

	return ρ
end

# ╔═╡ d732b00a-acbc-4772-b7b9-103dd90a9357
smallfontsize = 0.8 * theme(:fontsize)[]

# ╔═╡ 2d9e9ecf-32b8-4aa0-8fbb-a5bb24b7e9a1
function plot_tidal_track_expected!(df)
	x, y = LilGuys.EN21_tidal_track(df.r_circ_max[1], df.v_circ_max[1], x_min=df.r_circ_max[end] / df.r_circ_max[1])
	y *= V2KMS
	lines!(x, y, color=:black,linestyle=:dash, label="EN21 tidal track", linewidth=smalllinewidth)

	Utils.text_along_line_log!(x, y, x[end], text="EN2021", align=(:left, :bottom), fontsize=smallfontsize)
end

# ╔═╡ cc375e6c-4e7e-4f84-a55e-d72a8ee59ee3
get_mean_density(modelnames["scl_smallperi"]...)

# ╔═╡ b4de3539-9f82-45b0-b6a0-f51035a173f7
function ρ_to_v(ρ, r)
	v_scale = sqrt(ρ * 4π/3)
	v = v_scale .* r
end

# ╔═╡ add019df-a6de-46f5-8c24-74761251ff70
function plot_density_line!(ρ)
	
	r = 10 .^ LinRange(-2, 2, 10)
	v = ρ_to_v(3ρ, r)
	lines!(r, v*V2KMS, color=:grey, linewidth=smalllinewidth)
	Utils.text_along_line_log!(r, v*V2KMS, 0.3, text=L"3\bar\rho_\textrm{MW,\ peri}", align=(:center, :bottom), color=:grey)
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
function plot_tidal_track(galaxyname, modelname, starsname; title="")
	profs = get_profiles(galaxyname, modelname, starsname)
	fig = Figure(size=(4.5*72, 3.5*72))

	ax = Axis(fig[1,1], 
		xlabel = "radius / kpc",
		ylabel = L"cicular velocity / km\,s$^{-1}$",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
		yminorticks = [1:10; 15; 25; 35],
		title = title
	)
	
	plot_i_f(profs, color=COLORS[1])
	plot_tidal_track!(profs.scalars)
	plot_tidal_track_expected!(profs.scalars)


	# stars
	plot_i_f_stars(profs, color=COLORS[2], linewidth=smalllinewidth)

	# annotations
	plot_density_line!(get_mean_density(galaxyname, modelname, starsname)...)

	r_j = get_r_j(galaxyname, modelname)

	y0_j = LilGuys.lerp(log10.(profs[1].radii), middle.(LilGuys.circular_velocity(profs[1])))(log10(r_j)) * V2KMS
	plot_r_j!(r_j, 	y0_j)

	
	i_f_legend(ax, profs)

	v0 = profs.scalars.v_circ_max[1]*V2KMS
	xlims!(0.03, profs.scalars.r_circ_max[1] * 10^1)
	ylims!(v0* 10^-1.9, v0 * 10^0.1)

	fig
end

# ╔═╡ 415f249c-de75-454f-8ab8-be3bec7e4222
@savefig "scl_tidal_track"  plot_tidal_track(modelnames["scl_smallperi"]..., 
											 title="Sculptor: smallperi orbit")

# ╔═╡ f40dab86-d30b-4ba5-95b2-929b872dd887
R_h_3d_to_2d = LilGuys.R_h(LilGuys.Exp2D()) / LilGuys.r_h(LilGuys.Exp2D())

# ╔═╡ 119b9896-f53a-4d17-a7df-53e2408d9156
@savefig "umi_tidal_track"  plot_tidal_track(modelnames["umi_smallperi"]..., title="Ursa Minor: smallperi orbit")

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═043fd044-e44c-4c9f-861e-22db402e1b10
# ╠═d992c565-3786-44dc-835f-9e6d7cf9a094
# ╠═837df33b-02ee-4312-a1f8-db60b0baeba1
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═aeb0721a-f9e7-4a52-bc76-6fe0f0282f91
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╟─dcf00943-acab-4805-ad19-8a9d6b80f32b
# ╠═19a931d3-490b-4450-a1ed-b43ca337ba43
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═62ec2189-8fdd-4e7e-b679-05aade350dd0
# ╠═882e0efa-d365-40f2-a609-d2840b6a66ac
# ╠═74857d22-6166-4c0c-88f0-961f42cc159f
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═2c799458-b133-4f34-ac59-e2f83dba5c68
# ╠═aa77e11f-65a0-40cd-af1c-43c43fb555ef
# ╠═2971bd21-552f-4dc4-a863-26a20e2b2cd8
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═50a763c6-7db8-4e5e-a4e7-a3aceb010b6b
# ╠═07fe2997-559e-4bef-b4f4-934d850cdfe2
# ╠═59e9658c-7a2c-4964-8a24-fdf7b152ebcc
# ╠═8b41afc4-a86a-4d1f-937d-3ff7a0e9cdd7
# ╠═59fee278-e9be-4733-9356-3d3bf9c08310
# ╠═20dbd304-e705-4dde-8ea4-558151e63f53
# ╠═51c4cca9-3cc9-4242-9d41-53b37b3d8f04
# ╠═7f039e64-fbff-4091-924c-26523911074c
# ╠═2d9e9ecf-32b8-4aa0-8fbb-a5bb24b7e9a1
# ╠═70090f1a-0f6a-4082-a735-3b8168e26855
# ╠═3d55cd39-5727-48f8-a8ae-de168e9e08c8
# ╠═d732b00a-acbc-4772-b7b9-103dd90a9357
# ╠═add019df-a6de-46f5-8c24-74761251ff70
# ╠═cc375e6c-4e7e-4f84-a55e-d72a8ee59ee3
# ╠═b4de3539-9f82-45b0-b6a0-f51035a173f7
# ╠═415f249c-de75-454f-8ab8-be3bec7e4222
# ╠═f40dab86-d30b-4ba5-95b2-929b872dd887
# ╠═119b9896-f53a-4d17-a7df-53e2408d9156
