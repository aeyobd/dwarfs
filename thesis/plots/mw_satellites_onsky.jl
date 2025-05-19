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
include("./paper_style.jl")

# ╔═╡ 97993507-1ab6-4aae-8363-f439fb98c2e5
FIGDIR = "./figures"

# ╔═╡ 261a9e5c-93af-48d2-8a71-2c9e70ee8b1f
abbreviations = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/iau_abbrev.tsv"), DataFrame, delim='\t')

# ╔═╡ 38857534-663a-4764-b1d2-e9b861daf862
function shorten_name(name)
	for row in reverse(eachrow(abbreviations))
		if contains(name, row.nomative)
			name = replace(name, row.nomative => row.abbreviation)
		end
	end
	name
end

# ╔═╡ 3b149121-2b73-4bae-bd94-877e581ea5e6
shorten_name("Sagittarius")

# ╔═╡ 1a775735-4db9-4507-95e6-ea7ca052143c
function read_pace(name)
	CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/pace/$name.csv"), DataFrame)
end

# ╔═╡ f21a3642-d321-4d35-9472-df8049dab6c0
ambiguous_table = read_pace("gc_ambiguous")

# ╔═╡ 2f76352c-fc02-4274-8350-e99ce3b7ccdc
alldwarfs = read_pace("pace_all")

# ╔═╡ 915259bc-05ef-40c5-acf1-c49c1fba19f1
allcluster = vcat(read_pace("gc_harris"), read_pace("gc_mw_new"))

# ╔═╡ 0b63949c-24dc-40f6-8683-27a64a78459d


# ╔═╡ ece9f710-a071-4b0e-a6df-667d12a4638f
dwarfs.ll

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

# ╔═╡ 0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
CairoMakie.activate!(type=:png)

# ╔═╡ a7b54b66-9d44-46ab-a4c2-d9326acccda5
fg_color = (:white, 0.8)

# ╔═╡ 7564bd74-d182-4d89-bf57-2d0fe286bc7a
grid_color = (:white, 0.3)

# ╔═╡ 3051cdfb-63fb-44e4-a75d-37a51fd40508
@savefig "mw_satellites_onsky" let
	fig = Figure(backgroundcolor=:transparent)
	ax = GeoAxis(fig[1,1];
		dest = "+proj=aitoff",
		#yticks = [-75:15:75],
		#yticks = collect(-180:15.:180),
		#limits = [-180., 180, -90, 90],
		xgridcolor=grid_color,
		ygridcolor=grid_color,
		yticklabelsvisible=false,
		xticklabelsvisible=false,
		yticklabelsize=8,
		xreversed=false,
	)

	scatter!(ax, -confirmed_dwarfs.ll, confirmed_dwarfs.bb, markersize=3, label="dwarf")
	scatter!(ax, -classicals.ll, classicals.bb, markersize=5, label="classical")
	

	# scatter!(-magellanic.ll, magellanic.bb, markersize=10)
	for row in eachrow(magellanic)
		color = ifelse(row.key == "lmc", COLORS[4], fg_color)
		
		text!(ax, -row.ll, row.bb, text=row.name, fontsize=8, align=(:left, :center), offset=(4.0, 0.0),  color=color)
	end

	scatter!(-key_dwarfs.ll, key_dwarfs.bb, marker=:rect, color=COLORS[4], markersize=8, label="discussed here")

	# lines!(fill(180, 1000), LinRange(-90, 90, 1000), linewidth=1, color=:black)
	# lines!(fill(-180, 1000), LinRange(-90, 90, 1000), linewidth=1, color=:black)
	
	Legend(fig[1, 1], ax, tellwidth=false, tellheight=false, halign=:right, valign=:bottom)

	for row in eachrow(classicals)
		if row.name ∈ ["Sculptor", "Ursa Minor", "Fornax"]
			color = COLORS[4]
		else
			color = fg_color
		end
	
		text!(ax, -row.ll, row.bb, text=row.name, fontsize=8, align=(:left, :center), offset=(4.0, 0.0),  color=color)
	end


	#text!(-(150:-30:-150), zeros(11), text=string.(mod.(150:-30:-150, 360)), fontsize=8, color=grid_color, offset=(2, 0.))

	print(ax.yticks[])
	fig
end

# ╔═╡ 70477210-9891-4225-a544-b0030feb7af2
# ╠═╡ disabled = true
#=╠═╡
dwarfs = allsat[.!ismissing.(allsat.host) .& (allsat.host .∈ [["mw", "lmc", "smc"]]), :]
  ╠═╡ =#

# ╔═╡ 923f7988-c344-49e4-98c8-ee4b69527fcc
dwarfs = alldwarfs[.!ismissing.(alldwarfs.host) .& (alldwarfs.host .∈ [["mw", "lmc", "smc"]]) .& (alldwarfs.confirmed_dwarf .== 1), :]

# ╔═╡ Cell order:
# ╠═4c0da6d9-c5d9-4971-8b69-7ec9c107d981
# ╠═97993507-1ab6-4aae-8363-f439fb98c2e5
# ╠═726c1519-d147-4c7c-9901-55c0d4fab221
# ╠═261a9e5c-93af-48d2-8a71-2c9e70ee8b1f
# ╠═38857534-663a-4764-b1d2-e9b861daf862
# ╠═3b149121-2b73-4bae-bd94-877e581ea5e6
# ╠═1a775735-4db9-4507-95e6-ea7ca052143c
# ╠═f21a3642-d321-4d35-9472-df8049dab6c0
# ╠═2f76352c-fc02-4274-8350-e99ce3b7ccdc
# ╠═915259bc-05ef-40c5-acf1-c49c1fba19f1
# ╠═923f7988-c344-49e4-98c8-ee4b69527fcc
# ╟─70477210-9891-4225-a544-b0030feb7af2
# ╠═0b63949c-24dc-40f6-8683-27a64a78459d
# ╠═ece9f710-a071-4b0e-a6df-667d12a4638f
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
# ╠═0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
# ╠═a7b54b66-9d44-46ab-a4c2-d9326acccda5
# ╠═7564bd74-d182-4d89-bf57-2d0fe286bc7a
# ╠═3051cdfb-63fb-44e4-a75d-37a51fd40508
