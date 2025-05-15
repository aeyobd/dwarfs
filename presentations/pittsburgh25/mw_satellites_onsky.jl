### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 4c0da6d9-c5d9-4971-8b69-7ec9c107d981
begin
	using Pkg; Pkg.activate()

using CSV, DataFrames
	using CairoMakie, Arya

	using SkyCoords

	using LilGuys

	using GeoMakie
end

# ╔═╡ 726c1519-d147-4c7c-9901-55c0d4fab221
include("./style.jl")

# ╔═╡ 97993507-1ab6-4aae-8363-f439fb98c2e5
FIGDIR = "./figures"

# ╔═╡ 36021dc4-d970-4425-b39e-c3cbad8ce0cb
update_theme!(
	fonts = (; regular = "TeX Gyre Heros Makie",
			bold = "TeX Gyre Heros Makie Bold"
			)
)

# ╔═╡ f3881abc-9fb3-4b9a-89f6-3790c21fc6ed
abbreviations = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/iau_abbrev.tsv"), DataFrame, delim='\t')

# ╔═╡ 9179b4ed-eb7f-4a27-8719-fc028078a8c6
function shorten_name(name)
	for row in reverse(eachrow(abbreviations))
		if contains(name, row.nomative)
			name = replace(name, row.nomative => row.abbreviation)
		end
	end
	name
end

# ╔═╡ 11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
function read_pace(name)
	CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/pace/$name.csv"), DataFrame)
end

# ╔═╡ 2f76352c-fc02-4274-8350-e99ce3b7ccdc
alldwarfs = read_pace("pace_all")

# ╔═╡ a232715d-a27f-499b-8123-9adf5257acf5
ambiguous_table = read_pace("gc_ambiguous")

# ╔═╡ 88d507e8-98cb-42ca-be20-16580f5beb4b
read_pace("gc_other")

# ╔═╡ bce24033-3bca-42d4-9224-b40378cdedab
allcluster = vcat(read_pace("gc_harris"), read_pace("gc_mw_new"))

# ╔═╡ af7b55de-8ae6-4b19-80b8-dcfe5cae276f
baumgardt_columns = string.(split("Cluster         RA        DEC      R_Sun  DRSun    R_GC   DRGC   N_RV   N_PM    Mass         DM        V  Delta_V  M/L_V  DM/L    rc     rh,l     rh,m     rt    rho_c  rho_h,m   sig_c   sig_h,m lg(Trh) lg(Mini) T_Diss  M_Low  M_High    MF  Delta_MF  sig0    vesc   etac   etah  A_Rot Delta_AR P_Rot", r"\s+"))

# ╔═╡ 9f623db7-32db-43ce-9ddd-98d5ee4248ba
baumgardt = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/baumgardt_23.txt"), DataFrame, comment="#", ignorerepeated=true, delim=' ', header=baumgardt_columns)

# ╔═╡ 70477210-9891-4225-a544-b0030feb7af2
dwarfs = alldwarfs[.!ismissing.(alldwarfs.host) .& (alldwarfs.host .∈ [["mw", "lmc", "smc"]]), :]

# ╔═╡ 31d1d520-ee54-4d57-b0b4-73936f51c441
shorten_name.(dwarfs.name)[1]

# ╔═╡ 6142560a-9cd1-491f-a108-b3b45b85df2e
ambiguous = vcat(ambiguous_table, dwarfs[(dwarfs.confirmed_dwarf .=== 0), :])

# ╔═╡ ece9f710-a071-4b0e-a6df-667d12a4638f
dwarfs.ll

# ╔═╡ 1238d6ac-75eb-492d-add8-7132643609c1
function to_gal(ra, dec)
	icrs = @. ICRSCoords(deg2rad(ra), deg2rad(dec))
	gals = convert.(GalCoords, icrs)
	l = [rad2deg(g.l) for g in gals]
	b = [rad2deg(g.b) for g in gals]

	return l, b
end

# ╔═╡ d098bbec-0e99-4aaa-9096-717c785c17e9
icrs = [ICRSCoords(deg2rad(row.ra), deg2rad(row.dec)) for row in eachrow(dwarfs)]

# ╔═╡ aeb2540f-3832-42fe-a5f7-709bcf2932b8
gals = convert.(GalCoords, icrs)

# ╔═╡ b9a6058c-da01-44f5-965b-b2c55d5b124e
l = [rad2deg(g.l) for g in gals]

# ╔═╡ 39b09310-e4af-47db-ac4a-d4d7afe5302c
b = [rad2deg(g.b) for g in gals]

# ╔═╡ 6d8cfb4c-5602-4112-a691-92a319f5e224
classical_systems = ["sculptor_1", "ursa_minor_1", "sextans_1", "draco_1", 
		 "leo_1", "leo_2", "fornax_1", "carina_1", ]

# ╔═╡ ce9d55e5-c840-4f5c-a2de-7c3538f42bd2
classicals = dwarfs[dwarfs.key .∈ [classical_systems], :]

# ╔═╡ 5918ea56-83ad-48b7-8e7b-a4a5038e236b
magellanic = dwarfs[dwarfs.key .∈ [["smc", "lmc"]], :]

# ╔═╡ d8e0b701-166b-438c-9c36-5564e7e488c9
confirmed_dwarfs = dwarfs[(dwarfs.confirmed_dwarf .=== 1) .&
	(dwarfs.key .∉ [["smc", "lmc"]]), :]

# ╔═╡ 42c02a3d-5d7d-4559-8911-b197408343f1
key_dwarfs = dwarfs[dwarfs.key .∈ [["sculptor_1", "ursa_minor_1", "fornax_1"]], :]

# ╔═╡ 17236ffc-dbe5-4dbc-8b1d-fe490c3b50f4
magellanic[:, [:key, :ll, :bb]]

# ╔═╡ 89ae71c5-08eb-4c92-9301-ad2e630bf943
labels = vcat(confirmed_dwarfs.key, "ursa_major_3", "lmc", "smc")

# ╔═╡ 0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
CairoMakie.activate!(type=:png)

# ╔═╡ 594fbbef-40c8-4cb8-9d56-14c39c65a914
allsatalites = vcat(alldwarfs, allcluster, ambiguous_table)

# ╔═╡ a7b54b66-9d44-46ab-a4c2-d9326acccda5
fg_color = (:white, 0.5)

# ╔═╡ 7564bd74-d182-4d89-bf57-2d0fe286bc7a
grid_color = (:white, 0.3)

# ╔═╡ fa060f84-e39f-4d6a-a6d8-421f70a3c20e
theme(:fonts)

# ╔═╡ 3051cdfb-63fb-44e4-a75d-37a51fd40508
@savefig "mw_satellites_onsky" let
	fig = Figure(backgroundcolor=:transparent, size=(1920, 1080))
	
	ax = GeoAxis(fig[1,1];
		dest = "+proj=hammer",
		#limits = (0., 360, -90, 90),
		xgridcolor=grid_color,
		ygridcolor=grid_color,
		yticklabelsvisible=false,
		xticklabelsvisible=false,
		yticklabelsize=8,
		xgridwidth=0.5,
		ygridwidth=0.5,
	)
	xlims!(-180, 180)

	scatter!(ax, -allcluster.ll, allcluster.bb, markersize=2, label="cluster")

	scatter!(ax, -confirmed_dwarfs.ll, confirmed_dwarfs.bb, markersize=3, label="dwarf")
	scatter!(ax, -ambiguous.ll, ambiguous.bb, markersize=4, label="ambiguous")

	scatter!(ax, -classicals.ll, classicals.bb, markersize=6, label="classical")
	

	for row in eachrow(magellanic)
		color = ifelse(row.key == "lmc", COLORS[5], COLORS[5])
		
	
	end

	
	Legend(fig[1, 1], ax, tellwidth=false, tellheight=false, halign=:right, valign=:bottom, nbanks=2, backgroundcolor=:transparent, labelcolor=fg_color, framecolor=fg_color)

	for key in labels
		row = allsatalites[allsatalites.key .== key, :]
		if size(row, 1) != 1
			@info "not found $key"
		end
		
		name = row.name[1]
		offset = (3., 0.)

		if name ∈ ["Sculptor", "Ursa Minor", "Fornax"]
			color = COLORS[4]
			fontsize=48
		elseif name ∈ ["LMC", "SMC"]
			fontsize=36 
			color=COLORS[5]
			offset = (0., 0.)
		else
			color = fg_color
			fontsize=24
		end

		text!(ax, -row.ll, row.bb, text=shorten_name(name), fontsize=fontsize, align=(:left, :center), offset=offset,  color=color)
	end
	

	fig
end

# ╔═╡ 82e8db20-91cc-43b4-8183-7e302ce74041
@savefig "mw_satellites_onsky_all" let
	fig = Figure(backgroundcolor=:transparent, size=(1920, 1080))
	
	ax = GeoAxis(fig[1,1];
		dest = "+proj=hammer",
		#limits = (0., 360, -90, 90),
		xgridcolor=grid_color,
		ygridcolor=grid_color,
		yticklabelsvisible=false,
		xticklabelsvisible=false,
		yticklabelsize=8,
		xgridwidth=0.5,
		ygridwidth=0.5,
	)
	xlims!(-180, 180)

	scatter!(ax, -allcluster.ll, allcluster.bb, markersize=2, label="cluster")

	scatter!(ax, -confirmed_dwarfs.ll, confirmed_dwarfs.bb, markersize=3, label="dwarf")
	scatter!(ax, -ambiguous.ll, ambiguous.bb, markersize=4, label="ambiguous")

	scatter!(ax, -classicals.ll, classicals.bb, markersize=6, label="classical")
	

	for row in eachrow(magellanic)
		color = ifelse(row.key == "lmc", COLORS[5], COLORS[5])
		
	
	end

	
	Legend(fig[1, 1], ax, tellwidth=false, tellheight=false, halign=:right, valign=:bottom, nbanks=2, backgroundcolor=:transparent, labelcolor=fg_color, framecolor=fg_color)

	for key in labels
		row = allsatalites[allsatalites.key .== key, :]
		if size(row, 1) != 1
			@info "not found $key"
		end
		
		name = row.name[1]
		offset = (3., 0.)

		color = fg_color
		fontsize=24


		text!(ax, -row.ll, row.bb, text=shorten_name(name), fontsize=fontsize, align=(:left, :center), offset=offset,  color=color)
	end
	

	fig
end

# ╔═╡ 25514746-3762-4cb0-8b03-0e3e38d82a69


# ╔═╡ Cell order:
# ╠═4c0da6d9-c5d9-4971-8b69-7ec9c107d981
# ╠═97993507-1ab6-4aae-8363-f439fb98c2e5
# ╠═36021dc4-d970-4425-b39e-c3cbad8ce0cb
# ╠═726c1519-d147-4c7c-9901-55c0d4fab221
# ╠═f3881abc-9fb3-4b9a-89f6-3790c21fc6ed
# ╠═9179b4ed-eb7f-4a27-8719-fc028078a8c6
# ╠═31d1d520-ee54-4d57-b0b4-73936f51c441
# ╠═2f76352c-fc02-4274-8350-e99ce3b7ccdc
# ╠═a232715d-a27f-499b-8123-9adf5257acf5
# ╠═6142560a-9cd1-491f-a108-b3b45b85df2e
# ╠═11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
# ╠═88d507e8-98cb-42ca-be20-16580f5beb4b
# ╠═bce24033-3bca-42d4-9224-b40378cdedab
# ╠═af7b55de-8ae6-4b19-80b8-dcfe5cae276f
# ╠═9f623db7-32db-43ce-9ddd-98d5ee4248ba
# ╠═70477210-9891-4225-a544-b0030feb7af2
# ╠═ece9f710-a071-4b0e-a6df-667d12a4638f
# ╠═1238d6ac-75eb-492d-add8-7132643609c1
# ╠═d098bbec-0e99-4aaa-9096-717c785c17e9
# ╠═aeb2540f-3832-42fe-a5f7-709bcf2932b8
# ╠═b9a6058c-da01-44f5-965b-b2c55d5b124e
# ╠═39b09310-e4af-47db-ac4a-d4d7afe5302c
# ╠═6d8cfb4c-5602-4112-a691-92a319f5e224
# ╠═ce9d55e5-c840-4f5c-a2de-7c3538f42bd2
# ╠═5918ea56-83ad-48b7-8e7b-a4a5038e236b
# ╠═d8e0b701-166b-438c-9c36-5564e7e488c9
# ╠═42c02a3d-5d7d-4559-8911-b197408343f1
# ╠═17236ffc-dbe5-4dbc-8b1d-fe490c3b50f4
# ╠═89ae71c5-08eb-4c92-9301-ad2e630bf943
# ╠═0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
# ╠═594fbbef-40c8-4cb8-9d56-14c39c65a914
# ╠═a7b54b66-9d44-46ab-a4c2-d9326acccda5
# ╠═7564bd74-d182-4d89-bf57-2d0fe286bc7a
# ╠═fa060f84-e39f-4d6a-a6d8-421f70a3c20e
# ╠═3051cdfb-63fb-44e4-a75d-37a51fd40508
# ╠═82e8db20-91cc-43b4-8183-7e302ce74041
# ╠═25514746-3762-4cb0-8b03-0e3e38d82a69
