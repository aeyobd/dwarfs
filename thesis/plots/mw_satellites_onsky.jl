### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 4c0da6d9-c5d9-4971-8b69-7ec9c107d981
begin
	using Pkg; Pkg.activate()

	using CSV, DataFrames
	using CairoMakie, GeoMakie
	
	using LilGuys, Arya
end

# ╔═╡ 222f3946-abd7-4e06-91fb-f83a4cb3239a
let
	include("./paper_style.jl")

	width  = 1920 / 72
    set_theme!(theme_arya(width=width, fontsize=40, px_per_unit=1))

    scale_theme_element!(:linewidth, 2)
    legend_attr = theme(:Legend)
    legend_attr.margin = legend_attr.padding
    update_theme!(Legend = legend_attr, figure_padding=theme(:fontsize)[])

	update_theme!(
	fonts = (; regular = "TeX Gyre Heros Makie",
			bold = "TeX Gyre Heros Makie Bold"
			)
)

	Arya.update_figsize!(width, 1080/1920 * width)

	CairoMakie.activate!(type=:png, px_per_unit=4)
end

# ╔═╡ 7216773a-bf24-49bb-a4ab-4afb35ce2a2b
md"""
# Setup
"""

# ╔═╡ 1ee09c26-16be-4919-a2c8-8a21612b9873
import SkyCoords

# ╔═╡ 75878bfc-6d21-4662-aa88-e2d501c06d95
import FileIO

# ╔═╡ 02df901a-cbae-454c-a907-32aaf07bf3e1
img = FileIO.load(("resources/Gaia_s_sky_in_colour.png"))

# ╔═╡ 97993507-1ab6-4aae-8363-f439fb98c2e5
FIGDIR = "./figures"

# ╔═╡ 5c124928-d05d-467a-93f9-dc4e9e3b1d04
md"""
We want to use the IAU abbreviations
"""

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

# ╔═╡ 483e5816-7f41-476f-9696-096ce83d4340
md"""
# Data loading
"""

# ╔═╡ 11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
function read_pace(name)
	CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/pace/$name.csv"), DataFrame)
end

# ╔═╡ 2f76352c-fc02-4274-8350-e99ce3b7ccdc
alldwarfs = read_pace("dwarf_all")

# ╔═╡ 516eb72e-44e6-4637-ba75-8b0c5ce6f12a
 dwarfs = read_pace("dwarf_mw")

# ╔═╡ 70477210-9891-4225-a544-b0030feb7af2
# dwarfs = alldwarfs[.!ismissing.(alldwarfs.host) .& (alldwarfs.host .∈ [["mw", "lmc", "smc"]]), :] # reproduces table

# ╔═╡ a232715d-a27f-499b-8123-9adf5257acf5
ambiguous_table = read_pace("gc_ambiguous")

# ╔═╡ 45499b5b-2c79-4a0e-8cc8-930581a8dec7
function to_galcoords(ra, dec)
	coords = SkyCoords.ICRSCoords.(deg2rad.(ra), deg2rad.(dec))
	galcoords = SkyCoords.convert.(SkyCoords.GalCoords, coords)

	l = SkyCoords.lon.(galcoords) .|> rad2deg
	b = SkyCoords.lat.(galcoords) .|> rad2deg

	return l, b
end

# ╔═╡ bce24033-3bca-42d4-9224-b40378cdedab
allcluster = vcat(read_pace("gc_harris"), read_pace("gc_mw_new"))

# ╔═╡ 88d507e8-98cb-42ca-be20-16580f5beb4b
read_pace("gc_other")

# ╔═╡ 0a1d1783-0c24-455e-84c1-027e168cfef7
md"""
# Filtering & classification
"""

# ╔═╡ 6d8cfb4c-5602-4112-a691-92a319f5e224
classical_systems = ["sculptor_1", "ursa_minor_1", "sextans_1", "draco_1", 
		 "leo_1", "leo_2", "fornax_1", "carina_1", "antlia_2", "sagittarius_1", "canes_venatici_1", "crater_2"]

# ╔═╡ fa47cae6-3879-4c62-8f7b-56aa12e8d105
md"""
- Sagittarius
- Fornax
- Leo I
- Sculptor
- Antlia II
- Leo II
- Carina
- Draco
- Ursa Minor
- Canes Venatici I
- Sextans I
- Crater II. 
"""

# ╔═╡ 6142560a-9cd1-491f-a108-b3b45b85df2e
ambiguous = vcat(ambiguous_table, dwarfs[(dwarfs.confirmed_dwarf .!== 1), :])

# ╔═╡ ce9d55e5-c840-4f5c-a2de-7c3538f42bd2
classicals = dwarfs[dwarfs.key .∈ [classical_systems], :]

# ╔═╡ 5918ea56-83ad-48b7-8e7b-a4a5038e236b
magellanic = alldwarfs[alldwarfs.key .∈ [["smc", "lmc", "m_033"]], :]

# ╔═╡ a09fe200-9239-4e74-a8c4-c548e67481f4
andromida = let
	df = DataFrame(
		name = ["M31", "M33"], 
		ra = [010.684708, 023.4620690621800], #
		dec = [41.268750, 30.6601751119800],
		ra_ref = ["2006AJ....131.1163S", "2020yCat.1350....0G"],
	)

	df[!, :ll], df[!, :bb] = to_galcoords(df.ra, df.dec)

	df[!, :key] = lowercase.(replace.(df.name, " "=>"_"))
	df
end

# ╔═╡ d16c6607-b4ea-4648-bf8d-d6a02feacbd2
alldwarfs.key[.!ismissing.(alldwarfs.host) .& (alldwarfs.host .== "m_031")]

# ╔═╡ b2d74e34-b6d9-4069-a398-9e33b5ea405c
unique(alldwarfs.host)

# ╔═╡ d8e0b701-166b-438c-9c36-5564e7e488c9
confirmed_dwarfs = dwarfs[(dwarfs.confirmed_dwarf .=== 1), :]

# ╔═╡ a401d47d-e064-429e-baf3-d38ba64b3c4a
confirmed_faint_dwarfs =  confirmed_dwarfs[confirmed_dwarfs.key .∉ [["smc", "lmc", classical_systems...]], :]

# ╔═╡ 42c02a3d-5d7d-4559-8911-b197408343f1
key_dwarfs = dwarfs[dwarfs.key .∈ [["sculptor_1", "ursa_minor_1", "fornax_1"]], :]

# ╔═╡ 17236ffc-dbe5-4dbc-8b1d-fe490c3b50f4
magellanic[:, [:key, :ll, :bb]]

# ╔═╡ 89ae71c5-08eb-4c92-9301-ad2e630bf943
labels = vcat(confirmed_faint_dwarfs.key, "ursa_major_3", "lmc", "smc", classical_systems,  andromida.key)

# ╔═╡ 594fbbef-40c8-4cb8-9d56-14c39c65a914
allsatalites = vcat(dwarfs, allcluster, ambiguous_table, andromida, cols=:union)

# ╔═╡ 095ed289-8086-4da3-aaeb-36dec2b4d349
md"""
# Plot utils
"""

# ╔═╡ a7b54b66-9d44-46ab-a4c2-d9326acccda5
fg_color = (:white, 0.8)

# ╔═╡ 7564bd74-d182-4d89-bf57-2d0fe286bc7a
grid_color = (:white, 0.3)

# ╔═╡ bddcd6bb-5dba-4b6e-ba07-e21df433f812
ms = theme(:markersize)[] * 1.5

# ╔═╡ b7fd5546-ecc9-444b-9862-1954c3c762e1
CairoMakie.activate!(type=:png)

# ╔═╡ c954772f-9e75-484e-bfdb-4455249b1f99
function clear_axis!(gs)	
	ax = GeoAxis(gs;
		dest = "+proj=hammer",
		limits = (0., 360, -90, 90),
		xgridcolor=grid_color,
		ygridcolor=grid_color,
		yticklabelsvisible=false,
		xticklabelsvisible=false,
		yticklabelsize=8,
		xgridwidth=0.5,
		ygridwidth=0.5,
				 valign=:center,
				 
	)
	xlims!(-180, 180)

	return ax
end

# ╔═╡ 21f0d8c8-76a2-41f0-a390-390c5b539ed3
function clear_figure(;backgroundcolor=:transparent, padding=0)
	fig = Figure(backgroundcolor=backgroundcolor, figure_padding = padding)
end

# ╔═╡ 8ea50ebf-87d3-4810-9699-89617463d9b4
function dark_axis()
	fig = Figure(backgroundcolor=:transparent, figure_padding = 0)
	
	ax = GeoAxis(fig[1,1];
		dest = "+proj=hammer",
		#limits = (0., 360, -90, 90),
		yticklabelsvisible=false,
		xticklabelsvisible=false,
		yticksvisible=false,
		xticksvisible=false,
		yticklabelsize=8,
		xgridwidth=0.5,
		ygridwidth=0.5,
		titlegap=0
	)
	xlims!(-180, 180)

	return fig, ax
end

# ╔═╡ 21373373-2dec-4bd8-a178-6af90df0835a
function clear_legend(fig, ax)
	
	Legend(fig[1, 1], ax, tellwidth=false, tellheight=false, halign=:right, valign=:bottom, nbanks=2, backgroundcolor=:transparent, labelcolor=fg_color, framecolor=grid_color, pad=0, margin=theme(:Legend).framewidth[] .* ones(4))
end

# ╔═╡ 5a537690-3c29-4a0b-9a31-15d19e456db6
function plot_points!(ax, x=:ll, y=:bb; xreverse=true, red=COLORS[4])

	if xreverse
		xfactor = -1
	else
		xfactor = 1
	end
	
	scatter!(ax, xfactor * allcluster[!, x], allcluster[!, y], label="cluster",
		markersize=2/3*ms, color=red, marker=:circle)

	scatter!(ax, xfactor * confirmed_faint_dwarfs[!, x], confirmed_faint_dwarfs[!, y], label="dwarf", 
		markersize=ms, color=COLORS[1], marker=:rect)
	
	scatter!(ax, xfactor * ambiguous[!, x], ambiguous[!, y], label="ambiguous",
		markersize=ms/3*2, color=:transparent, marker=:hexagon, 
		strokewidth=ms/8, strokecolor=COLORS[5])

	scatter!(ax, xfactor * classicals[!, x], classicals[!, y], label="classical",
		markersize=5/4*ms, color=COLORS[3], marker=:diamond, )


end

# ╔═╡ 3b31bf9e-1d07-460a-80d9-75625d06bf41
abbreviations

# ╔═╡ 756bb316-6afa-4217-b52b-eff39e5c02d6
function plot_labels!(ax; highlight=[])
	for key in labels
		row = allsatalites[allsatalites.key .== key, :]
		if size(row, 1) != 1
			@info "not found $key"
		end
		
		name = row.name[1]
		offset = (12., 0.)
		fontsize=24
		color = fg_color
		align = (:left, :center)

		text=shorten_name(name)
		fonts = :regular

		if key ∈ highlight
			color = COLORS[9]
			fontsize=48
			offset = (24., 0.)

		elseif key ∈ ["ursa_major_1", "leo_1", "leo_2", "pegasus_3", "pisces_2", "horologium_1"]
			align = (:left, :top)
		elseif key ∈ ["leo_5",]
			align = (:left, 0.25)
		elseif key ∈ ["tucana_2",  "columba_1"]
			align = (:left, :bottom)
			offset = (0, 8.)
		elseif key ∈ ["smc"]
		elseif key == "lmc"
			offset = (18., 0.)
		elseif text == "Car II"
			text = ""
		elseif text == "Car III"
			text = "Car II & III"
		elseif text == "Boo I"
			text = ""
		elseif text == "Boo II"
			text = "Boo I & II"
		elseif key ∈ ["grus_1"]
			align = (:right, :center)
			offset = (-offset[1], 0.)
		elseif key ∈ classical_systems
			offset = (18., 0.)
		end


		text!(ax, -row.ll, row.bb, text=text, fontsize=fontsize, align=align, offset=offset,  color=color, font=fonts)
	end

end

# ╔═╡ 3f44fb3f-865a-4acb-9101-283c7918e11b
function acknowledgment!(ax)
	text!(ax, 0, 0, text="Daniel Boyea, 2025\nBackground image: ESA/Gaia/DPAC\nData from Local Volume Database (A. Pace, 2024)", space=:relative, color=fg_color, fontsize=24, offset=(0, 0))
end

# ╔═╡ 2c2e54c6-4d2f-40b1-b012-5d21dc3bd4b9
function plot_gaia_image!(gs)
	ax2 = Axis(gs, aspect=DataAspect(), backgroundcolor=:transparent,)
	image!(ax2, rotr90(img), interpolate=true)
	hidespines!(ax2)
	hidedecorations!(ax2)
end

# ╔═╡ 89e31154-62d6-4233-9a5f-a0480efd8a18
md"""
# The plots
"""

# ╔═╡ 3051cdfb-63fb-44e4-a75d-37a51fd40508
@savefig "mw_satellites" let
	fig = clear_figure(backgroundcolor=:black, padding=20)
	plot_gaia_image!(fig[1,1])
	ax = Makie.current_axis()
	ax.valign = :top

	
	ax = clear_axis!(fig[1,1])

	
	plot_points!(ax, red=COLORS[4])
	clear_legend(fig, ax)


	#acknowledgment!(ax)

	plot_labels!(ax, highlight=["sculptor_1", "ursa_minor_1"])

	resize_to_layout!(fig)

	fig
end

# ╔═╡ 784b840d-08b8-4276-8dd4-6916a7f155ce
@savefig "mw_satellites_w_ack" let
	fig = clear_figure(backgroundcolor=:black, padding=20)
	plot_gaia_image!(fig[1,1])
	ax = Makie.current_axis()
	ax.valign = :top

	
	ax = clear_axis!(fig[1,1])

	
	plot_points!(ax, red=COLORS[4])
	clear_legend(fig, ax)

	ax3 = Axis(fig[1,1], backgroundcolor=:transparent)
	hidespines!(ax3)
	hidedecorations!(ax3)
	acknowledgment!(ax3)

	
	
	plot_labels!(ax, highlight=["sculptor_1", "ursa_minor_1"])

	resize_to_layout!(fig)

	fig
end

# ╔═╡ 4541aa9d-867f-4dca-a456-5dd8fc4fba9b
@savefig "mw_satellites_all" let
	fig = clear_figure(backgroundcolor=:black, padding=20)
	plot_gaia_image!(fig[1,1])
	ax = Makie.current_axis()
	ax.valign = :top

	
	ax = clear_axis!(fig[1,1])

	
	plot_points!(ax, red=COLORS[4])
	clear_legend(fig, ax)

	acknowledgment!(ax)


	plot_labels!(ax)

	resize_to_layout!(fig)

	fig
end

# ╔═╡ Cell order:
# ╟─7216773a-bf24-49bb-a4ab-4afb35ce2a2b
# ╠═4c0da6d9-c5d9-4971-8b69-7ec9c107d981
# ╠═1ee09c26-16be-4919-a2c8-8a21612b9873
# ╠═75878bfc-6d21-4662-aa88-e2d501c06d95
# ╠═222f3946-abd7-4e06-91fb-f83a4cb3239a
# ╠═02df901a-cbae-454c-a907-32aaf07bf3e1
# ╠═97993507-1ab6-4aae-8363-f439fb98c2e5
# ╟─5c124928-d05d-467a-93f9-dc4e9e3b1d04
# ╠═f3881abc-9fb3-4b9a-89f6-3790c21fc6ed
# ╠═9179b4ed-eb7f-4a27-8719-fc028078a8c6
# ╟─483e5816-7f41-476f-9696-096ce83d4340
# ╠═11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
# ╠═2f76352c-fc02-4274-8350-e99ce3b7ccdc
# ╠═516eb72e-44e6-4637-ba75-8b0c5ce6f12a
# ╠═70477210-9891-4225-a544-b0030feb7af2
# ╠═a232715d-a27f-499b-8123-9adf5257acf5
# ╠═45499b5b-2c79-4a0e-8cc8-930581a8dec7
# ╠═bce24033-3bca-42d4-9224-b40378cdedab
# ╠═88d507e8-98cb-42ca-be20-16580f5beb4b
# ╟─0a1d1783-0c24-455e-84c1-027e168cfef7
# ╠═6d8cfb4c-5602-4112-a691-92a319f5e224
# ╟─fa47cae6-3879-4c62-8f7b-56aa12e8d105
# ╠═6142560a-9cd1-491f-a108-b3b45b85df2e
# ╠═ce9d55e5-c840-4f5c-a2de-7c3538f42bd2
# ╠═5918ea56-83ad-48b7-8e7b-a4a5038e236b
# ╠═a09fe200-9239-4e74-a8c4-c548e67481f4
# ╠═d16c6607-b4ea-4648-bf8d-d6a02feacbd2
# ╠═b2d74e34-b6d9-4069-a398-9e33b5ea405c
# ╠═d8e0b701-166b-438c-9c36-5564e7e488c9
# ╠═a401d47d-e064-429e-baf3-d38ba64b3c4a
# ╠═42c02a3d-5d7d-4559-8911-b197408343f1
# ╠═17236ffc-dbe5-4dbc-8b1d-fe490c3b50f4
# ╠═89ae71c5-08eb-4c92-9301-ad2e630bf943
# ╠═594fbbef-40c8-4cb8-9d56-14c39c65a914
# ╟─095ed289-8086-4da3-aaeb-36dec2b4d349
# ╠═a7b54b66-9d44-46ab-a4c2-d9326acccda5
# ╠═7564bd74-d182-4d89-bf57-2d0fe286bc7a
# ╠═bddcd6bb-5dba-4b6e-ba07-e21df433f812
# ╠═b7fd5546-ecc9-444b-9862-1954c3c762e1
# ╠═c954772f-9e75-484e-bfdb-4455249b1f99
# ╠═21f0d8c8-76a2-41f0-a390-390c5b539ed3
# ╠═8ea50ebf-87d3-4810-9699-89617463d9b4
# ╠═21373373-2dec-4bd8-a178-6af90df0835a
# ╠═5a537690-3c29-4a0b-9a31-15d19e456db6
# ╠═3b31bf9e-1d07-460a-80d9-75625d06bf41
# ╠═756bb316-6afa-4217-b52b-eff39e5c02d6
# ╠═3f44fb3f-865a-4acb-9101-283c7918e11b
# ╠═2c2e54c6-4d2f-40b1-b012-5d21dc3bd4b9
# ╟─89e31154-62d6-4233-9a5f-a0480efd8a18
# ╠═3051cdfb-63fb-44e4-a75d-37a51fd40508
# ╠═784b840d-08b8-4276-8dd4-6916a7f155ce
# ╠═4541aa9d-867f-4dca-a456-5dd8fc4fba9b
