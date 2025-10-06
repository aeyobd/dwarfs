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

# ╔═╡ d992c565-3786-44dc-835f-9e6d7cf9a094
import TOML

# ╔═╡ aeb0721a-f9e7-4a52-bc76-6fe0f0282f91
module Utils
	include("./utils.jl")
end

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 8081e889-c1d4-4e50-b8e1-91aa8e82adf1
function modeldir(modelname)
	joinpath(ENV["DWARFS_ROOT"], "analysis",modelname,)
end

# ╔═╡ 6161cd69-594b-4eb9-83d9-12a682db7c88
function get_scalars(modelname)
	read_fits(joinpath(modeldir(modelname), "profiles_scalars.fits")
	)
end

# ╔═╡ 62ec2189-8fdd-4e7e-b679-05aade350dd0
function get_idx_f(modelname)
	return TOML.parsefile(joinpath(modeldir(modelname), "orbital_properties.toml"))["idx_f"]
end

# ╔═╡ 74857d22-6166-4c0c-88f0-961f42cc159f
function get_r_j(modelname, lmc=false)
	model_dir = modeldir(modelname)
	if lmc
		props = TOML.parsefile(model_dir * "/jacobi_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "/jacobi.toml")
	end

	return props["r_J_kpc"]
end

# ╔═╡ 2c799458-b133-4f34-ac59-e2f83dba5c68
function MassProfile(filename)
	prof_kwargs = TOML.parsefile(filename)

	prof_kwargs = LilGuys.collapse_errors(prof_kwargs) |> LilGuys.dict_to_tuple

	return LilGuys.MassProfile(;prof_kwargs...)
end

# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(modelname)
	mass_profs = LilGuys.read_ordered_structs(
		joinpath(modeldir(modelname), "profiles.hdf5"),
		LilGuys.MassProfile
	)

	prof_f_u = MassProfile(joinpath(modeldir(modelname),  "mass_profile_unbound.toml"))

	return mass_profs[1].second, mass_profs[get_idx_f(modelname)].second, prof_f_u
end

# ╔═╡ aa77e11f-65a0-40cd-af1c-43c43fb555ef
function get_stellar_profiles(modelname, starsname; Mstar=1)
	prof_i = MassProfile(joinpath(modeldir(modelname), "stars", starsname, "stellar_profile_i_3d_mass.toml"))
	prof_f = MassProfile(joinpath(modeldir(modelname), "stars", starsname, "stellar_profile_f_3d_mass.toml"))
	prof_f_u = MassProfile(joinpath(modeldir(modelname), "stars", starsname, "stellar_profile_unbound_f_3d_mass.toml"))

	prof_i = LilGuys.scale(prof_i, 1, Mstar / prof_i.M_in[end].middle)
	prof_f = LilGuys.scale(prof_f, 1, Mstar / prof_f.M_in[end].middle)

	return prof_i, prof_f
end

# ╔═╡ 11d210d4-3b54-42d9-89f0-b4677ec32348
function get_stellar_scalars(modelname, starsname; Mstar=1)
	scalars = read_fits(joinpath(modeldir(modelname), "stars", starsname, "stellar_profiles_3d_scalars.fits"))

	return scalars
end

# ╔═╡ 2971bd21-552f-4dc4-a863-26a20e2b2cd8
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ 50a763c6-7db8-4e5e-a4e7-a3aceb010b6b
function i_f_legend(ax, profs)
	t_i = round(-profs[end].time * T2GYR, digits=1)
	l1 = lines!([NaN], [NaN], linestyle=:dot, color=:black)
	l2 = lines!([NaN], [NaN], linestyle=:solid, color=:black)
	l3 = scatter!([NaN], [NaN], color=COLORS[1])
	axislegend(ax, [l1, l2, l3], [L"initial ($t=%$t_i$ Gyr)", L"final ($t=0$ Gyr)", "velocity max", ], position=:rb)
end

# ╔═╡ 07fe2997-559e-4bef-b4f4-934d850cdfe2
function plot_r_j!(r_j, y0)	
	arrows2d!([r_j], [y0], [0], [y0 * (1 - 10^0.10)], shaftwidth=theme(:linewidth)[]/2, tipwidth=6, tiplength=6)
	text!(r_j, y0, text="Jacobi", rotation=π/2, align=(:left, :center))

end

# ╔═╡ 7f039e64-fbff-4091-924c-26523911074c
function text_along_line!(x, y, x_0; text, h=0.03, kwargs...)
	f = LilGuys.lerp(x, y)
	y_0 = f(x_0)
	dy = Utils.log_derivative(f, x_0, h=h)
	rf = Utils.rotation_factor(Makie.current_axis(), true)
	θ = @lift atan($rf * dy)

	@info "x0 y0 text", x_0, y_0
	text!(x_0, y_0, text=text, rotation=θ; kwargs...)
end

# ╔═╡ 59fee278-e9be-4733-9356-3d3bf9c08310
function plot_i_f(profs; kwargs...)
	prof_i, prof_f = profs

	x_i, y_i = radii(prof_i), middle.(LilGuys.circular_velocity(prof_i)) * V2KMS
	lines!(x_i, y_i; linestyle=:dot, kwargs...)

	x, y = radii(prof_f), middle.(LilGuys.circular_velocity(prof_f)) * V2KMS

	lines!(x, y; kwargs...)


	text_along_line!(x_i, y_i, 0.05, text="dark matter", color=COLORS[1], align=(:left, :bottom), h=0.3)
	
end


# ╔═╡ 20dbd304-e705-4dde-8ea4-558151e63f53
function plot_i_f_stars(profs; kwargs...)
	prof_i, prof_f = profs
	x_i, y_i = radii(prof_i), middle.(LilGuys.circular_velocity(prof_i) * V2KMS)
	lines!(x_i, y_i, linestyle=:dot; kwargs...)

	x, y = radii(prof_f), middle.(LilGuys.circular_velocity(prof_f) * V2KMS)
	lines!(x, y; kwargs...)
	#lines!(radii(prof_f), LilGuys.circular_velocity(prof_f) * V2KMS, linestyle=:solid; kwargs...)

	text_along_line!(x_i, y_i, 0.05, text="stars", color=COLORS[2], align=(:left, :bottom))

end

# ╔═╡ d732b00a-acbc-4772-b7b9-103dd90a9357
smallfontsize = 0.8 * theme(:fontsize)[]

# ╔═╡ 2d9e9ecf-32b8-4aa0-8fbb-a5bb24b7e9a1
function plot_tidal_track!(df)
	x, y = LilGuys.EN21_tidal_track(df.r_circ_max[1], df.v_circ_max[1], x_min=df.r_circ_max[end] / df.r_circ_max[1])
	y *= V2KMS
	lines!(x, y, color=:black,linestyle=:dash, label="EN21 tidal track", linewidth=smalllinewidth)

	text_along_line!(x, y, x[end], text="EN2021", align=(:left, :bottom), fontsize=smallfontsize)
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
function plot_tidal_track(modelname, starsname; Mstar=1, title="", y0_j)
	df = get_scalars(modelname)
	profs = get_profiles(modelname)
	profs_stars = get_stellar_profiles(modelname, starsname, Mstar=Mstar)
	stellar_scalars = get_stellar_scalars(modelname, starsname,)
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "radius / kpc",
		ylabel = L"cicular velocity / km\,s$^{-1}$",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
		yminorticks = [1:10; 15; 25; 35],
		title = title
	)

	r_j = get_r_j(modelname)


	
	plot_i_f(profs, color=COLORS[1])
	scatter!(df.r_circ_max, df.v_circ_max * V2KMS, markersize=2/3 * theme(:markersize)[], color=COLORS[1])


	# stars
	plot_i_f_stars(profs_stars, color=COLORS[2], linewidth=smalllinewidth)

	plot_tidal_track!(df)
	
	# axislegend(position=:rb, patchsize=(24, 12), backgroundcolor=(:white, 2))

	v0 = df.v_circ_max[1]*V2KMS
	plot_r_j!(r_j, 	y0_j)
	
	i_f_legend(ax, profs)

	xlims!(0.03, df.r_circ_max[1] * 10^1)
	ylims!(v0* 10^-1.9, v0 * 10^0.1)

	fig
end

# ╔═╡ 415f249c-de75-454f-8ab8-be3bec7e4222
@savefig "scl_tidal_track"  plot_tidal_track("sculptor/1e7_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10", Mstar=3.1e-4, title="Sculptor: smallperi orbit", y0_j=3)

# ╔═╡ f40dab86-d30b-4ba5-95b2-929b872dd887
R_h_3d_to_2d = LilGuys.R_h(LilGuys.Exp2D()) / LilGuys.r_h(LilGuys.Exp2D())

# ╔═╡ 119b9896-f53a-4d17-a7df-53e2408d9156
@savefig "umi_tidal_track"  plot_tidal_track("ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5", "exp2d_rs0.10", title="Ursa Minor: smallperi orbit", Mstar=7e-5, y0_j=2)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═d992c565-3786-44dc-835f-9e6d7cf9a094
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═aeb0721a-f9e7-4a52-bc76-6fe0f0282f91
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═62ec2189-8fdd-4e7e-b679-05aade350dd0
# ╠═74857d22-6166-4c0c-88f0-961f42cc159f
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═2c799458-b133-4f34-ac59-e2f83dba5c68
# ╠═aa77e11f-65a0-40cd-af1c-43c43fb555ef
# ╠═11d210d4-3b54-42d9-89f0-b4677ec32348
# ╠═2971bd21-552f-4dc4-a863-26a20e2b2cd8
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═50a763c6-7db8-4e5e-a4e7-a3aceb010b6b
# ╠═07fe2997-559e-4bef-b4f4-934d850cdfe2
# ╠═59fee278-e9be-4733-9356-3d3bf9c08310
# ╠═20dbd304-e705-4dde-8ea4-558151e63f53
# ╠═7f039e64-fbff-4091-924c-26523911074c
# ╠═2d9e9ecf-32b8-4aa0-8fbb-a5bb24b7e9a1
# ╠═d732b00a-acbc-4772-b7b9-103dd90a9357
# ╠═415f249c-de75-454f-8ab8-be3bec7e4222
# ╠═f40dab86-d30b-4ba5-95b2-929b872dd887
# ╠═119b9896-f53a-4d17-a7df-53e2408d9156
