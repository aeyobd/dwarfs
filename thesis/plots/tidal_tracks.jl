### A Pluto.jl notebook ###
# v0.20.18

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
	prof_f_u = LilGuys.scale(prof_f_u, 1, Mstar / prof_f_u.M_in[end].middle)

	return prof_i, prof_f, prof_f_u
end

# ╔═╡ 72e106eb-8aee-445e-bf07-dc3c427aa229


# ╔═╡ 11d210d4-3b54-42d9-89f0-b4677ec32348
function get_stellar_scalars(modelname, starsname; Mstar=1)
	scalars = read_fits(joinpath(modeldir(modelname), "stars", starsname, "stellar_profiles_3d_scalars.fits"))

	return scalars
end

# ╔═╡ 0a4e4850-71c5-4b81-bc24-eb2daf78b341
function get_Jacobi_radius(modelname)
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
function plot_tidal_track(modelname, starsname; Mstar=1, title="")
	df = get_scalars(modelname)
	profs = get_profiles(modelname)
	profs_stars = get_stellar_profiles(modelname, starsname, Mstar=Mstar)
	stellar_scalars = get_stellar_scalars(modelname, starsname,)
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$r_\textrm{circ}$ / kpc",
		ylabel = L"$\textrm{v}_\textrm{circ}$ / km\,s$^{-1}$",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
		yminorticks = [1:10; 15; 25; 35],
		title = title
	)

	r_j = get_r_j(modelname)


	x, y = LilGuys.EN21_tidal_track(df.r_circ_max[1], df.v_circ_max[1], x_min=0.1)

	
	scatter!(df.r_circ_max, df.v_circ_max * V2KMS, label="max", markersize=theme(:markersize)[]/2, color=COLORS[1])

	lines!(radii(profs[1]), LilGuys.circular_velocity(profs[1]) * V2KMS, label=L"DM initial ($t=-9.2$\,Gyr)", color=COLORS[1], linestyle=:dot)
	lines!(radii(profs[2]), LilGuys.circular_velocity(profs[2]) * V2KMS, label=L"DM final ($t=0.0$\,Gyr)", color=COLORS[1], linestyle=:solid)
	lines!(radii(profs[3]), LilGuys.circular_velocity(profs[3]) * V2KMS, label="DM final + unbound", color=COLORS[1], linestyle=:dash)
	lines!(x, y * V2KMS, color=:black,linestyle=:dot, label="EN21 tidal track")

	# stars
	lines!(radii(profs_stars[1]), LilGuys.circular_velocity(profs_stars[1]) * V2KMS, label="stars initial", color=COLORS[2], linestyle=:dot)
	lines!(radii(profs_stars[end]), LilGuys.circular_velocity(profs_stars[end]) * V2KMS, label="stars final", color=COLORS[2], linestyle=:solid)
	
	scatter!(stellar_scalars.r_h, stellar_scalars.sigma_v * V2KMS, color=COLORS[3])

	
	text!(stellar_scalars.r_h[1], stellar_scalars.sigma_v[1] * V2KMS, text=L"r_h, \sigma_\textrm{v}", color=COLORS[3])

	
	axislegend(position=:rb, patchsize=(24, 12), backgroundcolor=(:white, 2))

	vlines!(r_j, color=:black, linewidth=theme(:linewidth)[]/2)
	text!(r_j, 8, text="Jacobi", rotation=π/2)

	xlims!(0.03, df.r_circ_max[1] * 10^1)
	ylims!(1, df.v_circ_max[1]*V2KMS * 10^0.1)

	fig
end

# ╔═╡ 415f249c-de75-454f-8ab8-be3bec7e4222
@savefig "scl_tidal_track"  plot_tidal_track("sculptor/1e7_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10", Mstar=3.1e-4, title="Sculptor Dark Matter: smallperi orbit")

# ╔═╡ f40dab86-d30b-4ba5-95b2-929b872dd887
R_h_3d_to_2d = LilGuys.R_h(LilGuys.Exp2D()) / LilGuys.r_h(LilGuys.Exp2D())

# ╔═╡ 119b9896-f53a-4d17-a7df-53e2408d9156
@savefig "umi_tidal_track"  plot_tidal_track("ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5", "exp2d_rs0.10", title="Ursa Minor: smallperi orbit", Mstar=7e-5)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═d992c565-3786-44dc-835f-9e6d7cf9a094
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═62ec2189-8fdd-4e7e-b679-05aade350dd0
# ╠═74857d22-6166-4c0c-88f0-961f42cc159f
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═2c799458-b133-4f34-ac59-e2f83dba5c68
# ╠═aa77e11f-65a0-40cd-af1c-43c43fb555ef
# ╠═72e106eb-8aee-445e-bf07-dc3c427aa229
# ╠═11d210d4-3b54-42d9-89f0-b4677ec32348
# ╠═0a4e4850-71c5-4b81-bc24-eb2daf78b341
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═415f249c-de75-454f-8ab8-be3bec7e4222
# ╠═f40dab86-d30b-4ba5-95b2-929b872dd887
# ╠═119b9896-f53a-4d17-a7df-53e2408d9156
