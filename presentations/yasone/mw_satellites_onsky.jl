### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 4c0da6d9-c5d9-4971-8b69-7ec9c107d981
begin
	using Pkg; Pkg.activate()

	using CSV, DataFrames
	using CairoMakie, GeoMakie
	
	using LilGuys, Arya

	import FileIO
end

# ╔═╡ 726c1519-d147-4c7c-9901-55c0d4fab221
include("./style.jl")

# ╔═╡ 7216773a-bf24-49bb-a4ab-4afb35ce2a2b
md"""
# Setup
"""

# ╔═╡ c04cab38-6e6e-4097-a82f-91cc8a3864f1
CairoMakie.activate!(type=:png, px_per_unit=4)

# ╔═╡ 02df901a-cbae-454c-a907-32aaf07bf3e1
img = FileIO.load("Gaia_s_sky_in_colour.png")

# ╔═╡ 97993507-1ab6-4aae-8363-f439fb98c2e5
FIGDIR = "./figures"

# ╔═╡ 36021dc4-d970-4425-b39e-c3cbad8ce0cb
update_theme!(
	fonts = (; regular = "TeX Gyre Heros Makie",
			bold = "TeX Gyre Heros Makie Bold"
			)
)

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

# ╔═╡ a8b06d5c-ff4e-4220-8176-834d0c53b05e
md"""
# Data loading and sample categorization
"""

# ╔═╡ 3ca0b64f-3688-4a14-ab25-4a7667858fd8
alldwarfs = CSV.read("dwarf_galaxies.csv", DataFrame)

# ╔═╡ acec9bde-4b5b-473b-80b1-eecf1c302ec5
yasone_halo = CSV.read("yasone_halo.csv", DataFrame)

# ╔═╡ 00d140f9-609e-4cbb-bdfa-842f6e0a40fe
yasone_disk = CSV.read("yasone_disk.csv", DataFrame)

# ╔═╡ 93ddb9df-3686-41ce-8bf6-862525961cf5
andromida = CSV.read("andromida.csv", DataFrame)

# ╔═╡ a4b99fc1-1076-4995-b35e-7df7937bb113
allcluster = CSV.read("globular_clusters.csv", DataFrame)

# ╔═╡ 59c7fdaa-ef6a-4891-8538-ec4d6b08db1b
ambiguous = CSV.read("ambiguous.csv", DataFrame)

# ╔═╡ 5fc42d73-8a44-4a97-9584-4dccb5911d0f
@assert size(allcluster, 1) == 155

# ╔═╡ c4148f34-468b-4934-bbbc-f30b73c63e6a
@assert size(alldwarfs, 1) == 65

# ╔═╡ bbe6dced-82eb-4cfc-a10c-69b79ef51606
@assert size(ambiguous, 1) == 25

# ╔═╡ 89ae71c5-08eb-4c92-9301-ad2e630bf943
labels = vcat(alldwarfs.key,yasone_halo.key, andromida.key, yasone_disk.key)

# ╔═╡ 594fbbef-40c8-4cb8-9d56-14c39c65a914
allsatalites = vcat(alldwarfs, allcluster, ambiguous, andromida, yasone_halo,yasone_disk, cols=:union)

# ╔═╡ faa1cd3d-7014-45ad-a14b-8b72ace707ca
classical_systems = ["sculptor_1", "ursa_minor_1", "sextans_1", "draco_1", 
		 "leo_1", "leo_2", "fornax_1", "carina_1", "antlia_2", "sagittarius_1", "canes_venatici_1", "crater_2"]

# ╔═╡ 6e337185-9fa9-4a9d-9576-40a7ede77ace
magellanic = alldwarfs[alldwarfs.key .∈ [["smc", "lmc", "m_033"]], :]

# ╔═╡ d511095f-9628-49c0-afc2-1daabb0d9a78
faintdwarfs =  alldwarfs[alldwarfs.key .∉ [["smc", "lmc"]], :]

# ╔═╡ 095ed289-8086-4da3-aaeb-36dec2b4d349
md"""
# Plot utils
"""

# ╔═╡ a7b54b66-9d44-46ab-a4c2-d9326acccda5
fg_color = (:white, 0.8)

# ╔═╡ 7564bd74-d182-4d89-bf57-2d0fe286bc7a
grid_color = (:white, 0.3)

# ╔═╡ bddcd6bb-5dba-4b6e-ba07-e21df433f812
ms = theme(:markersize)[] * 3/2

# ╔═╡ 8e76042d-6223-497d-bdc3-63e10db7de28
padding = 40

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
	fig = Figure(backgroundcolor=backgroundcolor, size=(1920, 1080), figure_padding = padding)
end

# ╔═╡ 8ea50ebf-87d3-4810-9699-89617463d9b4
function dark_axis()
	fig = Figure(backgroundcolor=:transparent, size=(1920, 1080), figure_padding = 0)
	
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
function clear_legend(fig, args...)
	
	Legend(fig[1, 1], args..., 
		   tellwidth=false, tellheight=false, 
		   halign=:center, valign=:center, 
		   nbanks = 6, 
		   backgroundcolor = :transparent, 
		   labelcolor = fg_color, 
		   framecolor = fg_color, 
		  # margin = padding .* ones(4),
	)
end

# ╔═╡ ae0935e4-4e9e-404e-9b91-1ca4564b1d95
styles = Dict(
	"cluster" => (;markersize=2/3*ms, color=COLORS[4], marker=:circle),
	"dwarf" => (markersize=ms, color=COLORS[1], marker=:rect),
	"ambiguous" => (
		markersize=ms/3*2, color=:transparent, marker=:hexagon, 
		strokewidth=ms/8, strokecolor=COLORS[5]),
	"classical" => (
			 markersize=5/4*ms, 
			 color=COLORS[3], 
			 marker=:diamond, 
				   ),
	"yasone_halo" => (markersize=ms*3/2, color=COLORS[5], marker=:hexagon, strokecolor=fg_color, strokewidth=ms/16),
	"yasone_disk" => (markersize=ms*3/2, color=COLORS[4], marker=:star5, strokecolor=fg_color, strokewidth=ms/16, ),
	
)

# ╔═╡ 5a537690-3c29-4a0b-9a31-15d19e456db6
function plot_points!(ax, x=:ll, y=:bb; xreverse=true, red=red, xfactor=1)
	if xreverse
		xfactor *= -1
	end
	
	scatter!(ax, xfactor * allcluster[!, x],  allcluster[!, y];
			 label="globular clusters", styles["cluster"]...)

	scatter!(ax, xfactor * faintdwarfs[!, x], faintdwarfs[!, y]; label="dwarf galaxies", 
			 styles["dwarf"]...)
	
	scatter!(ax, xfactor * ambiguous[!, x], ambiguous[!, y];
			label="ambiguous systems",  styles["ambiguous"]...
			)

end

# ╔═╡ 9b5b8c2f-916b-4d8f-81a2-c8746ae00bc3
function plot_yasone!(ax; x=:ll, y=:bb, xreverse=true, red=red, yfactor=1, xfactor=1)
	if xreverse
		xfactor *= -1
	end
	
	scatter!(ax, xfactor*yasone_halo[!, x], yfactor * yasone_halo[!, y]; 
			 label="Yasone-halo", styles["yasone_halo"]...)

	scatter!(ax, xfactor*yasone_disk[!, x],  yfactor * yasone_disk[!, y], 
			 label="Yasone-disk"; styles["yasone_disk"]...)
end

# ╔═╡ fa0e1f1d-fa13-4bad-b879-185c4d7a5e57
labels

# ╔═╡ 756bb316-6afa-4217-b52b-eff39e5c02d6
function plot_labels!(ax; highlight=[], yasone=false, yasone_disk = false)
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
		elseif startswith(key, "yasone") & !yasone
			continue
		elseif key ∈ ["ursa_major_1", "leo_1", "leo_2", "pegasus_3", "pisces_2", "pictor_2", "horologium_2", "draco_2"]
			align = (:left, :top)
		elseif key ∈ ["horologium_1", "virgo_1",]
			align = (:center, :top)
			offset = (0, -offset[1])
		elseif key ∈ ["tucana_2",  "columba_1", "yasone_1"]
			align = (:left, :bottom)
			offset = (0, 8.)
		elseif key ∈ ["SMC"]
		elseif key == "LMC"
			offset = (18., 0.)
		elseif text == "Car II"
			text = ""
		elseif text == "Car III"
			text = "Car II & III"
		elseif text == "Leo IV"
			text = "Leo IV & V"
		elseif text == "Leo V"
			text = ""
			align = (:left, 0.25)
		elseif text == "Boo I"
			text = ""
		elseif text == "Boo II"
			text = "Boo I & II"
		elseif text == "Tuc III"
			text = ""
		elseif text == "Tuc IV"
			text = "Tuc III & IV"
		elseif key ∈ ["grus_1","sextans_2", ]
			align = (:right, :center)
			offset = (-offset[1], 0.)
		elseif key == "yasone_1"
			
		end

		if startswith(key, "yasone")
			if !yasone_disk && (parse(Int, split(key, "_")[2]) > 4)
				continue
			end
			fonts = :bold
			fontsize=28
		end

		text!(ax, -row.ll, row.bb, text=text, fontsize=fontsize, align=align, offset=offset,  color=color, font=fonts)
	end

end

# ╔═╡ 726168cd-a133-41d3-8bb4-12404994781b
function yasone_legend(fig)

	yasone_elems = [MarkerElement(; styles["yasone_halo"]...), MarkerElement(; styles["yasone_disk"]...)]

	leg2 = clear_legend(fig, yasone_elems,
						["yasone halo", "yasone disk"])


		#leg2.tellheight = true
	leg2.halign=:left
	leg2.valign=:bottom
	leg2.nbanks = 1
	leg2.margin = zeros(4)

	leg2
end

# ╔═╡ d7533590-6c5e-4bbf-8d9a-3f40e946d1c2
sort(allcluster.name)

# ╔═╡ 5c0e1220-8f79-41b0-855e-f5623ea5e411
sort(ambiguous.name)

# ╔═╡ f6ad426a-4ad0-4e3c-a78f-f01ad3aff597
let
	fig = clear_figure(backgroundcolor=(0.2), padding=padding)
	
	ax_scatter = clear_axis!(fig[1,1])
	# ax_scatter.color = fg_color
	
	plot_points!(ax_scatter, red=COLORS[4])
	plot_labels!(ax_scatter, yasone=true)
	plot_yasone!(ax_scatter)

	Box(fig[1, 1], color = :transparent, strokecolor = fg_color)

	leg = clear_legend(fig[2,1], ax_scatter)
	leg.tellheight = true
	leg.halign=:center
	leg.valign=:bottom
	leg.nbanks = 5
	leg.margin = zeros(4)
	rowgap!(fig.layout, 0)


	# leg2 = yasone_legend(fig[2,1])
	# leg2.tellheight=false
	
	resize_to_layout!(fig)
	fig
end

# ╔═╡ 01ca3652-fbfc-406d-b5ec-d594519d6e39
let
	fig = clear_figure(backgroundcolor=:black, padding=padding)
	
	ax_scatter = Axis(fig[1,1], backgroundcolor=:black,  xticklabelcolor=grid_color, bottomspinecolor=grid_color, xtickcolor=grid_color, xminortickcolor=grid_color, leftspinecolor=grid_color, ytickcolor=grid_color, yminortickcolor=grid_color, yticklabelcolor=grid_color)
	# ax_scatter.color = fg_color
	
	plot_points!(ax_scatter, red=COLORS[4], xreverse=true)
	plot_labels!(ax_scatter, yasone=true, yasone_disk=true)
	plot_yasone!(ax_scatter, xreverse=true)

	Box(fig[1, 1], color = :transparent, strokecolor = fg_color)

	leg = clear_legend(fig[2,1], ax_scatter)
	leg.tellheight = true
	leg.halign=:center
	leg.valign=:bottom
	leg.nbanks = 5
	leg.margin = zeros(4)
	rowgap!(fig.layout, 0)

	ylims!(ax_scatter, -5, 5)

	# leg2 = yasone_legend(fig[2,1])
	# leg2.tellheight=false
	
	resize_to_layout!(fig)
	fig
end

# ╔═╡ 2c2e54c6-4d2f-40b1-b012-5d21dc3bd4b9
function plot_gaia_image!(gs)
	ax2 = Axis(gs, aspect=DataAspect(), backgroundcolor=:transparent,)
	image!(ax2, rotr90(img), interpolate=true)
	hidespines!(ax2)
	hidedecorations!(ax2)
	ax2
end

# ╔═╡ 3051cdfb-63fb-44e4-a75d-37a51fd40508
@savefig "mw_satellites_yasone" let
	fig = clear_figure(backgroundcolor=:black, padding=padding/2)
	
	ax_img = plot_gaia_image!(fig[1,1])
	ax_scatter = clear_axis!(fig[1,1])

	plot_points!(ax_scatter)
	plot_labels!(ax_scatter, yasone=true)
	plot_yasone!(ax_scatter)

	leg = clear_legend(fig[2,1], ax_scatter)
	leg.tellheight = true
	leg.halign=:center
	leg.valign=:bottom
	leg.nbanks = 5
	leg.margin = zeros(4)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╟─7216773a-bf24-49bb-a4ab-4afb35ce2a2b
# ╠═4c0da6d9-c5d9-4971-8b69-7ec9c107d981
# ╠═c04cab38-6e6e-4097-a82f-91cc8a3864f1
# ╠═02df901a-cbae-454c-a907-32aaf07bf3e1
# ╠═97993507-1ab6-4aae-8363-f439fb98c2e5
# ╠═726c1519-d147-4c7c-9901-55c0d4fab221
# ╠═36021dc4-d970-4425-b39e-c3cbad8ce0cb
# ╟─5c124928-d05d-467a-93f9-dc4e9e3b1d04
# ╠═f3881abc-9fb3-4b9a-89f6-3790c21fc6ed
# ╠═9179b4ed-eb7f-4a27-8719-fc028078a8c6
# ╠═a8b06d5c-ff4e-4220-8176-834d0c53b05e
# ╠═3ca0b64f-3688-4a14-ab25-4a7667858fd8
# ╠═acec9bde-4b5b-473b-80b1-eecf1c302ec5
# ╠═00d140f9-609e-4cbb-bdfa-842f6e0a40fe
# ╠═93ddb9df-3686-41ce-8bf6-862525961cf5
# ╠═a4b99fc1-1076-4995-b35e-7df7937bb113
# ╠═59c7fdaa-ef6a-4891-8538-ec4d6b08db1b
# ╠═5fc42d73-8a44-4a97-9584-4dccb5911d0f
# ╠═c4148f34-468b-4934-bbbc-f30b73c63e6a
# ╠═bbe6dced-82eb-4cfc-a10c-69b79ef51606
# ╠═89ae71c5-08eb-4c92-9301-ad2e630bf943
# ╠═594fbbef-40c8-4cb8-9d56-14c39c65a914
# ╠═faa1cd3d-7014-45ad-a14b-8b72ace707ca
# ╠═6e337185-9fa9-4a9d-9576-40a7ede77ace
# ╠═d511095f-9628-49c0-afc2-1daabb0d9a78
# ╟─095ed289-8086-4da3-aaeb-36dec2b4d349
# ╠═a7b54b66-9d44-46ab-a4c2-d9326acccda5
# ╠═7564bd74-d182-4d89-bf57-2d0fe286bc7a
# ╠═bddcd6bb-5dba-4b6e-ba07-e21df433f812
# ╠═8e76042d-6223-497d-bdc3-63e10db7de28
# ╠═c954772f-9e75-484e-bfdb-4455249b1f99
# ╠═21f0d8c8-76a2-41f0-a390-390c5b539ed3
# ╠═8ea50ebf-87d3-4810-9699-89617463d9b4
# ╠═21373373-2dec-4bd8-a178-6af90df0835a
# ╠═ae0935e4-4e9e-404e-9b91-1ca4564b1d95
# ╠═5a537690-3c29-4a0b-9a31-15d19e456db6
# ╠═9b5b8c2f-916b-4d8f-81a2-c8746ae00bc3
# ╠═fa0e1f1d-fa13-4bad-b879-185c4d7a5e57
# ╠═756bb316-6afa-4217-b52b-eff39e5c02d6
# ╠═726168cd-a133-41d3-8bb4-12404994781b
# ╠═d7533590-6c5e-4bbf-8d9a-3f40e946d1c2
# ╠═5c0e1220-8f79-41b0-855e-f5623ea5e411
# ╠═f6ad426a-4ad0-4e3c-a78f-f01ad3aff597
# ╠═01ca3652-fbfc-406d-b5ec-d594519d6e39
# ╠═3051cdfb-63fb-44e4-a75d-37a51fd40508
# ╠═2c2e54c6-4d2f-40b1-b012-5d21dc3bd4b9
