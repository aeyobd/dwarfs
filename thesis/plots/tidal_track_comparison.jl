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

# ╔═╡ c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 8081e889-c1d4-4e50-b8e1-91aa8e82adf1
function modeldir(modelname)
	joinpath(ENV["DWARFS_ROOT"], "analysis",modelname,)
end

# ╔═╡ 6161cd69-594b-4eb9-83d9-12a682db7c88
function get_scalars(modelname)
	df = read_fits(joinpath(modeldir(modelname), "profiles_scalars.fits")
	)
	idx_f = TOML.parsefile(joinpath(modeldir(modelname), "orbital_properties.toml"))["idx_f"]

	if idx_f < size(df, 1)
		df.time .-= df.time[idx_f]
	else
		@warn "idx_f > number calculated profiles for $modelname"
		df.time .-= df.time[end]
	end
	df
end

# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(modelname)
	mass_profs = LilGuys.read_ordered_structs(
		joinpath(modeldir(modelname), "profiles.hdf5"),
		LilGuys.MassProfile
	)
	return mass_profs[1], mass_profs[end]
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
function plot_tidal_track(modelname; title="")
	df = get_scalars(modelname)
	profs = get_profiles(modelname)
	
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$r_\textrm{circ}$ / kpc",
		ylabel = L"$\textrm{v}_\textrm{circ}$ / km\,s$^{-1}$",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
		title = title
	)


	x, y = LilGuys.EN21_tidal_track(df.r_circ_max[1], df.v_circ_max[1], x_min=0.3)

	
	scatter!(df.r_circ_max, df.v_circ_max * V2KMS, label="max", markersize=theme(:markersize)[]/2)

	lines!(radii(profs[1].second), LilGuys.circular_velocity(profs[1].second) * V2KMS, label=L"initial ($t=-9.2$\,Gyr)")
	lines!(radii(profs[end].second), LilGuys.circular_velocity(profs[end].second) * V2KMS, label=L"final ($t=0.0$\,Gyr)")
	lines!(x, y * V2KMS, color=:black,linestyle=:dot, label="EN21 tidal track")

	axislegend(position=:lt, patchsize=(24, 12), backgroundcolor=(:white, 0.5))
	xlims!(df.r_circ_max[1] * 10^-1.5, df.r_circ_max[1] * 10^1)
	ylims!(df.v_circ_max[1]*V2KMS * 10^-0.7, df.v_circ_max[1]*V2KMS * 10^0.1)

	fig
end

# ╔═╡ 6830deae-07e6-4668-b600-7790145328c7
function compare_final(modelnames; kwargs...)
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$r_\textrm{circ}$ / kpc",
		ylabel = L"$\textrm{v}_\textrm{circ}$ / km\,s$^{-1}$",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
	)

	for (label, modelname) in modelnames

		profs = get_profiles(modelname)

		lines!(radii(profs[end].second), LilGuys.circular_velocity(profs[end].second) * V2KMS; label=label)
	end

	axislegend(position=:rt)
	fig
end

# ╔═╡ cf57e5bb-8829-409c-8189-79a1fdb9fddb
function compare_tidal_tracks(modelnames)
	
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$r_\textrm{circ}$ / kpc",
		ylabel = L"$\textrm{v}_\textrm{circ}$ / km\,s$^{-1}$",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
	)

	r_max = 0
	v_max = 0
	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)
		r_max = max(r_max, df.r_circ_max[1])
		v_max = max(v_max, df.v_circ_max[1] * V2KMS)

		lines!(df.r_circ_max, df.v_circ_max * V2KMS, label=label)
	end


	axislegend(position=:lt, patchsize=(24, 12), backgroundcolor=(:white, 0.5))

	df = get_scalars(first(modelnames).second)
	xlims!(r_max * 10^-1.5, r_max * 10^1)
	ylims!(v_max * 10^-0.7, v_max * 10^0.1)

	fig
end

# ╔═╡ 7d8f9954-0bbc-4fc0-ba21-aad54b3aaee4
function compare_boundmass(modelnames; relative=false, legend_position=:rb)
	
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "simulation time / Gyr",
		ylabel = L"$\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
	)

	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)

		y = log10.(df.v_circ_max * V2KMS)
		if relative
			y .-=  log10(df.v_circ_max[1])
		end
		
		lines!(df.time * T2GYR, y, label=label)
	end




	ax2 = Axis(fig[1,2], 
		xlabel = L"$r_\textrm{max}$ / kpc",
		ylabel = L"$\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
	)

	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)

		y = log10.(df.v_circ_max * V2KMS)
		if relative
			y .-=  log10(df.v_circ_max[1])
		end
		
		lines!(df.r_circ_max, y, label=label)
	end
	axislegend(position=legend_position, patchsize=(24, 12), backgroundcolor=(:white, 0.5))


	linkyaxes!(ax, ax2)
	hideydecorations!(ax2, ticks=false, minorticks=false)

	rowsize!(fig.layout, 1, Aspect(1, 3/3))
	resize_to_layout!()
	fig
end

# ╔═╡ 415f249c-de75-454f-8ab8-be3bec7e4222
plot_tidal_track("sculptor/1e7_new_v31_r3.2/orbit_smallperi")

# ╔═╡ fc8596bc-237b-4b67-a747-f01c708a4f22
plot_tidal_track("sculptor/1e6_new_v43_r7/orbit_smallperi.3")

# ╔═╡ 450f3db7-7c66-40b9-990a-24329097858f
plot_tidal_track("sculptor/1e6_v43_r5_beta0.2_a4/orbit_smallperi")

# ╔═╡ 4f5b3992-021a-4a9c-9ccc-750669156f28
plot_tidal_track("sculptor/1e6_v48_r7_oblate_0.5/orbit_smallperi")

# ╔═╡ 254ef828-1a85-4c40-bc36-9c94345c21c0
plot_tidal_track("sculptor/1e6_Ms0.54_rs1.08_c1/orbit_smallperi")

# ╔═╡ faff05dd-31b7-47bf-b6da-32d94ca80514
@savefig "scl_mw_halo_boundmass" compare_boundmass(OrderedDict(
	"fiducial" => "sculptor/1e7_new_v31_r3.2/orbit_smallperi",
	"heavier halo" => "sculptor/1e6_new_v43_r7/orbit_smallperi",
	"heavier, new orbit" => "sculptor/1e6_new_v43_r7/orbit_smallperi.3",
	"lighter halo" => "sculptor/1e6_new_v25_r2.5/orbit_smallperi",
),
				 )

# ╔═╡ ffad6fd8-a767-459e-aa3d-e9fd6741a47f
@savefig "scl_orbits_boundmass"  compare_boundmass(OrderedDict(
	"fiducial" => "sculptor/1e7_new_v31_r3.2/orbit_smallperi",
	"mean orbit" => "sculptor/1e6_new_v31_r3.2/orbit_mean",
	"mw impact" => "sculptor/1e6_new_v31_r3.2/L3M11_9Gyr_smallperi.a4",

))

# ╔═╡ a4c31b0c-41d7-498f-9d78-74bfc441c6fd
@savefig "scl_mw_structure_boundmass" compare_boundmass(OrderedDict(
	"heavy NFW" => "sculptor/1e6_new_v43_r7/orbit_smallperi",
	"cored" => "sculptor/1e6_Ms0.54_rs1.08_c1/orbit_smallperi",
	"anisotropic" =>  "sculptor/1e6_v43_r5_beta0.2_a4/orbit_smallperi",
	"oblate" => "sculptor/1e6_v48_r7_oblate_0.5/orbit_smallperi",
),
				 )

# ╔═╡ e12645e3-5466-4dd3-b5d2-2ee9288d1ade
@savefig "scl_lmc_halo_boundmass" let 
fig = compare_boundmass(OrderedDict(
	"lighter" => "sculptor/1e7_new_v25_r2.5/smallperilmc",
	"fiducial" => "sculptor/1e6_new_v31_r4.2/smallperilmc",
	"cored" => "sculptor/1e6_v49_r7.1_c0.1/smallperilmc",
	#"expcusp" => "sculptor/1e6_expcusp/vasiliev24_L3M11_extremeperi",
),)

	ylims!(1.2, 1.7)

	fig

end

# ╔═╡ 39a6433a-281a-42b4-a231-0fc6e9379422
@savefig "umi_orbits_boundmass" compare_boundmass(OrderedDict(
	"smallperi" => "ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5",
	"mean" => "ursa_minor/1e6_v37_r5.0/orbit_mean.2",

	"LMC mean" => "ursa_minor/1e5_new_v38_r4.0/L3M11_mean",
	"LMC sat" => "ursa_minor/1e5_new_v38_r4.0/L3M11_lmcsat",
	"LMC smallperi" => "ursa_minor/1e5_new_v38_r4.0/L3M11_smallperi",
#	"LMC smallperilmc" => "ursa_minor/1e5_new_v38_r4.0/L3M11_smallperilmc",
), )

# ╔═╡ 6bc1364d-ab67-414b-9519-4cb8706100af
compare_boundmass(OrderedDict(
	"fiducial" => "ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5",
	"fiducial/1" => "ursa_minor/1e5_new_v38_r4.0/orbit_smallperi.1",
	"cored" => "ursa_minor/1e5_v51_r4_c0.1//orbit_smallperi",
	"big cored" => "ursa_minor/1e5_Ms0.5_r0.9_c1/orbit_smallperi",
	# "big cored2" => "ursa_minor/1e5_Ms0.25_r0.9_c1/orbit_smallperi",
	"anisotropic" =>  "ursa_minor/1e5_v38_r3_beta0.2_a4//orbit_smallperi",
),
)

# ╔═╡ b44f0e10-370d-42d7-8d6d-b3ce1dec9048


# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═6830deae-07e6-4668-b600-7790145328c7
# ╠═cf57e5bb-8829-409c-8189-79a1fdb9fddb
# ╠═7d8f9954-0bbc-4fc0-ba21-aad54b3aaee4
# ╠═415f249c-de75-454f-8ab8-be3bec7e4222
# ╠═fc8596bc-237b-4b67-a747-f01c708a4f22
# ╠═450f3db7-7c66-40b9-990a-24329097858f
# ╠═4f5b3992-021a-4a9c-9ccc-750669156f28
# ╠═254ef828-1a85-4c40-bc36-9c94345c21c0
# ╠═faff05dd-31b7-47bf-b6da-32d94ca80514
# ╠═ffad6fd8-a767-459e-aa3d-e9fd6741a47f
# ╠═a4c31b0c-41d7-498f-9d78-74bfc441c6fd
# ╠═e12645e3-5466-4dd3-b5d2-2ee9288d1ade
# ╠═39a6433a-281a-42b4-a231-0fc6e9379422
# ╠═6bc1364d-ab67-414b-9519-4cb8706100af
# ╠═b44f0e10-370d-42d7-8d6d-b3ce1dec9048
