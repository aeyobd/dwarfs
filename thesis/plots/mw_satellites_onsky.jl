### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 4c0da6d9-c5d9-4971-8b69-7ec9c107d981
begin
	using Pkg; Pkg.activate()

	using CSV, DataFrames
	using CairoMakie, GeoMakie
	
	using LilGuys, Arya
end

# ╔═╡ 726c1519-d147-4c7c-9901-55c0d4fab221
include("./paper_style.jl")

# ╔═╡ 7216773a-bf24-49bb-a4ab-4afb35ce2a2b
md"""
# Setup
"""

# ╔═╡ 1ee09c26-16be-4919-a2c8-8a21612b9873
import SkyCoords

# ╔═╡ 75878bfc-6d21-4662-aa88-e2d501c06d95
import FileIO

# ╔═╡ 02df901a-cbae-454c-a907-32aaf07bf3e1
img = FileIO.load(("Gaia_s_sky_in_colour.png"))

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

# ╔═╡ c3442639-f5ad-4b60-b1a6-c9321c5a38c6
yasone_table = let
	df = CSV.read(
		Vector{UInt8}("""name,ra,dec
		Yasone 1,265.52019,13.17146
		Yasone 2,262.34921,6.42302
		Yasone 3,292.96258,-26.44994
		Yasone 4,225.68761,-0.91947
		"""),
			DataFrame
		)

	df[!, :key] = lowercase.(replace.(df.name, " "=>"_"))

	df[!, "ll"], df[!, "bb"] = to_galcoords(df.ra, df.dec)

	df

end

# ╔═╡ bce24033-3bca-42d4-9224-b40378cdedab
allcluster = vcat(read_pace("gc_harris"), read_pace("gc_mw_new"))

# ╔═╡ 71e21c9d-95c0-4397-929d-b415b0720dbf
hcat(to_galcoords(allcluster.ra, allcluster.dec)..., allcluster.ll, allcluster.bb)

# ╔═╡ 88d507e8-98cb-42ca-be20-16580f5beb4b
read_pace("gc_other")

# ╔═╡ 7258550c-c49a-4f67-83f7-58c157d4b5de
md"""
## Checking against baumgardt
"""

# ╔═╡ 357a7c6c-2d6a-40e5-ba85-7042ecc9bdda
md"""
The pace catalogue includes a few more globular clusters than Baumgardt, but reassuring to know that most are in both.
"""

# ╔═╡ af7b55de-8ae6-4b19-80b8-dcfe5cae276f
baumgardt_columns = string.(split("Cluster         RA        DEC      R_Sun  DRSun    R_GC   DRGC   N_RV   N_PM    Mass         DM        V  Delta_V  M/L_V  DM/L    rc     rh,l     rh,m     rt    rho_c  rho_h,m   sig_c   sig_h,m lg(Trh) lg(Mini) T_Diss  M_Low  M_High    MF  Delta_MF  sig0    vesc   etac   etah  A_Rot Delta_AR P_Rot", r"\s+"))

# ╔═╡ 9f623db7-32db-43ce-9ddd-98d5ee4248ba
baumgardt = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/baumgardt_23.txt"), DataFrame, comment="#", ignorerepeated=true, delim=' ', header=baumgardt_columns)

# ╔═╡ 58f03667-5903-48d7-a962-9fbbf5d72620
let
	fig = Figure(backgroundcolor=:transparent, size=(1920, 1080), pad = 0)
	
	ax = GeoAxis(fig[1,1];
		dest = "+proj=hammer",
		#limits = (0., 360, -90, 90),
		yticklabelsvisible=false,
		xticklabelsvisible=false,
		yticklabelsize=8,
		xgridwidth=0.5,
		ygridwidth=0.5,
				 valign=:top,
	
	)

	scatter!(baumgardt.RA, baumgardt.DEC, markersize=24)

	scatter!(allcluster.ra, allcluster.dec)

	fig
end

# ╔═╡ 0a1d1783-0c24-455e-84c1-027e168cfef7
md"""
# Filtering & classification
"""

# ╔═╡ 6d8cfb4c-5602-4112-a691-92a319f5e224
classical_systems = ["sculptor_1", "ursa_minor_1", "sextans_1", "draco_1", 
		 "leo_1", "leo_2", "fornax_1", "carina_1", ]

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
labels = vcat(confirmed_faint_dwarfs.key, "ursa_major_3", "lmc", "smc", classical_systems, yasone_table.key, andromida.key)

# ╔═╡ 0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
CairoMakie.activate!(type=:png)

# ╔═╡ 594fbbef-40c8-4cb8-9d56-14c39c65a914
allsatalites = vcat(dwarfs, allcluster, ambiguous_table, andromida, yasone_table, cols=:union)

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

# ╔═╡ d415fac3-d155-44fc-86af-35e6ead45734
COLORS

# ╔═╡ b27ae510-f567-4623-a2a8-d4f192d2f792
red = colorant"#e41a1c"

# ╔═╡ b7fd5546-ecc9-444b-9862-1954c3c762e1
CairoMakie.activate!(type=:png)

# ╔═╡ c954772f-9e75-484e-bfdb-4455249b1f99
function clear_axis!(gs)	
	ax = GeoAxis(gs;
		dest = "+proj=hammer",
		#limits = (0., 360, -90, 90),
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
function clear_legend(fig, ax)
	
	Legend(fig[1, 1], ax, tellwidth=false, tellheight=false, halign=:right, valign=:bottom, nbanks=2, backgroundcolor=:transparent, labelcolor=fg_color, framecolor=fg_color, pad=0, margin=theme(:Legend).framewidth[] .* ones(4))
end

# ╔═╡ 51eea365-1dd0-4ebf-8a89-30928c51ab5f
allcluster[!, :M_V]

# ╔═╡ 5a537690-3c29-4a0b-9a31-15d19e456db6
function plot_points!(ax, x=:ll, y=:bb; xreverse=true, red=red, yasone=false)

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

	if yasone
		scatter!(ax, xfactor*yasone_table[!, x], yasone_table[!, y], markersize=ms*3/2, color=COLORS[5], marker=:hexagon)
	end

end

# ╔═╡ 3b31bf9e-1d07-460a-80d9-75625d06bf41
abbreviations

# ╔═╡ 756bb316-6afa-4217-b52b-eff39e5c02d6
function plot_labels!(ax; highlight=[], yasone=false)
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
		elseif key ∈ ["ursa_major_1", "leo_1", "leo_2", "pegasus_3", "pisces_2", "horologium_1"]
			align = (:left, :top)
		elseif key ∈ ["leo_5",]
			align = (:left, 0.25)
		elseif key ∈ ["tucana_2",  "columba_1", "yasone_1"]
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
		elseif key ∈ ["grus_1", "yasone_4"]
			align = (:right, :center)
			offset = (-offset[1], 0.)
		elseif key ∈ classical_systems
			offset = (18., 0.)
		elseif key == "yasone_1"
			
		end

		if startswith(key, "yasone")
			fonts = :bold
			fontsize=28
		end

		text!(ax, -row.ll, row.bb, text=text, fontsize=fontsize, align=align, offset=offset,  color=color, font=fonts)
	end

end

# ╔═╡ c04cab38-6e6e-4097-a82f-91cc8a3864f1
CairoMakie.activate!(type=:png, px_per_unit=4)

# ╔═╡ 602f9f99-987b-4557-9709-b42025e4f54e
function scale_theme_element!(key, scale)
    @info "old $key value = $(theme(key)[])"
    update_theme!(; (; key => scale*theme(key)[])...)
    @info "new $key value = $(theme(key)[])"
end

# ╔═╡ 222f3946-abd7-4e06-91fb-f83a4cb3239a
let

    set_theme!(theme_arya(width=800/72 / Arya.HW_RATIO, fontsize=40, px_per_unit=1))
    CairoMakie.activate!(type=:svg, px_per_unit=1, pt_per_unit=1)

    scale_theme_element!(:linewidth, 2)
    legend_attr = theme(:Legend)
    legend_attr.margin = legend_attr.padding
    update_theme!(Legend = legend_attr, figure_padding=theme(:fontsize)[])

	update_theme!(
	fonts = (; regular = "TeX Gyre Heros Makie",
			bold = "TeX Gyre Heros Makie Bold"
			)
)
	
end

# ╔═╡ 15e5e568-560d-4a70-8ca3-b084d6d9c832
theme(:fontsize)

# ╔═╡ 3f44fb3f-865a-4acb-9101-283c7918e11b
function acknowledgment!(ax)
	text!(ax, 0, 0, text="Daniel Boyea, 2025\nBackground image: ESA/Gaia/DPAC\nData from Local Volume Database (A. Pace, 2024)", space=:relative, color=fg_color, fontsize=24, offset=(theme(:Legend).margin[][1], 0))
end

# ╔═╡ 58ea3ff0-40ab-4967-aa97-80d4aa903183


# ╔═╡ 2c2e54c6-4d2f-40b1-b012-5d21dc3bd4b9
function plot_gaia_image!(gs)
	ax2 = Axis(gs, aspect=DataAspect(), backgroundcolor=:transparent,)
	image!(ax2, rotr90(img), interpolate=true)
	hidespines!(ax2)
	hidedecorations!(ax2)
end

# ╔═╡ 3051cdfb-63fb-44e4-a75d-37a51fd40508
@savefig "mw_satellites" let
	fig = clear_figure(backgroundcolor=:black, padding=20)
	plot_gaia_image!(fig[1,1])
	ax = Makie.current_axis()
	ax.width = 1800
	ax.valign = :top

	
	ax = clear_axis!(fig[1,1])

	
	plot_points!(ax, yasone=false, red=COLORS[4])
	clear_legend(fig, ax)


	#acknowledgment!(ax)

	plot_labels!(ax, yasone=false, highlight=["sculptor_1", "ursa_minor_1"])

	resize_to_layout!(fig)

	fig
end

# ╔═╡ 82e8db20-91cc-43b4-8183-7e302ce74041


# ╔═╡ 25514746-3762-4cb0-8b03-0e3e38d82a69
let 
	fig = Figure( size=(1920, 1080))
	

	ax = Axis(fig[1,1], xscale = log10, xticks=Makie.automatic,
			xlabel = "galactocentric distance / kpc",
			  ylabel = "Mv",
			  yreversed=true
		)


	plot_points!(ax, :distance_gc, :M_V, xreverse=false)

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


		if key ∈ classical_systems
			offset = (18., 0.)
		elseif name ∈ ["LMC", "SMC", "M31", "M33"]
			continue 
			offset = (0., 0.)
		elseif startswith(key, "yasone")
			continue
		end

		@info row.distance_gc, row.M_V
		text!(ax, disallowmissing(row.distance_gc), disallowmissing(row.M_V), text=text, fontsize=fontsize, align=align, offset=offset,  color=:grey)
	end

	axislegend(position=:lb, margin=theme(:Legend).padding)

	@savefig "mw_satellites_distance_luminosity"
	fig
end

# ╔═╡ 6fb045c9-fce3-4d24-8edc-8ddde135972a
let 
	fig = Figure(size=(1920, 1080))
	

	ax = Axis(fig[1,1], xscale = log10, xticks=Makie.automatic,
			xlabel = "half light radius / parsecs",
			  ylabel = "Mv",
			  yreversed=true
		)


	plot_points!(ax, :rhalf_physical, :M_V, xreverse=false, yasone=false)

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


		if key ∈ classical_systems
			offset = (18., 0.)
		elseif name ∈ ["LMC", "SMC"]
			continue 
			offset = (0., 0.)
		elseif startswith(key, "yasone") | (key ∈ ["m31", "m33"])
			continue
		end

		text!(ax, disallowmissing(row.rhalf_physical), disallowmissing(row.M_V), text=text, fontsize=fontsize, align=align, offset=offset,  color=:grey)
	end

	axislegend(position=:rb, margin=theme(:Legend).padding)

	@savefig "mw_satellites_size_luminosity"
	fig
end

# ╔═╡ ebc0bd82-2b1e-4653-8aa2-46fbb56952ed
let 
	fig = Figure(size=(1920, 1080))
	

	ax = Axis(fig[1,1], xscale = log10, xticks=Makie.automatic,
			xlabel = "half light radius / arcmin",
			  ylabel = "Mv",
			  yreversed=true
		)


	plot_points!(ax, :rhalf, :M_V, xreverse=false, yasone=false)

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


		if key ∈ classical_systems
			offset = (18., 0.)
		elseif name ∈ ["LMC", "SMC"]
			continue 
			offset = (0., 0.)
		elseif startswith(key, "yasone") | (key ∈ ["m31", "m33"])
			continue
		end

		text!(ax, disallowmissing(row.rhalf), disallowmissing(row.M_V), text=text, fontsize=fontsize, align=align, offset=offset,  color=:grey)
	end

	axislegend(position=:rb, margin=theme(:Legend).padding)

	@savefig "mw_satellites_apparant_size_luminosity"
	fig
end

# ╔═╡ 7e544a72-62cf-4e0d-a23c-cba001a1e2fd
LilGuys.arcmin2kpc.(dwarfs.rhalf, dwarfs.distance)

# ╔═╡ d9201257-2e2f-4a21-8bb2-97ff349fe04b
dwarfs.rhalf_physical

# ╔═╡ 5f321550-4d2e-4f6c-8dbf-9ec1c72f1a44
icrs = [ICRS(ra=row.ra, dec=row.dec, distance=row.distance) for row in eachrow(confirmed_faint_dwarfs)]

# ╔═╡ 9326b0d6-847d-49fe-bfe3-380f328d5eb0
gc = LilGuys.transform.(Galactocentric, icrs)

# ╔═╡ c6140b22-0744-456e-a0fe-164d72c577af
confirmed_faint_dwarfs.distance_gc ./ radii.(LilGuys.position.(gc))

# ╔═╡ Cell order:
# ╟─7216773a-bf24-49bb-a4ab-4afb35ce2a2b
# ╠═4c0da6d9-c5d9-4971-8b69-7ec9c107d981
# ╠═1ee09c26-16be-4919-a2c8-8a21612b9873
# ╠═75878bfc-6d21-4662-aa88-e2d501c06d95
# ╠═222f3946-abd7-4e06-91fb-f83a4cb3239a
# ╠═02df901a-cbae-454c-a907-32aaf07bf3e1
# ╠═97993507-1ab6-4aae-8363-f439fb98c2e5
# ╠═726c1519-d147-4c7c-9901-55c0d4fab221
# ╟─5c124928-d05d-467a-93f9-dc4e9e3b1d04
# ╠═f3881abc-9fb3-4b9a-89f6-3790c21fc6ed
# ╠═9179b4ed-eb7f-4a27-8719-fc028078a8c6
# ╠═483e5816-7f41-476f-9696-096ce83d4340
# ╠═11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
# ╠═2f76352c-fc02-4274-8350-e99ce3b7ccdc
# ╠═516eb72e-44e6-4637-ba75-8b0c5ce6f12a
# ╠═70477210-9891-4225-a544-b0030feb7af2
# ╠═a232715d-a27f-499b-8123-9adf5257acf5
# ╠═45499b5b-2c79-4a0e-8cc8-930581a8dec7
# ╠═c3442639-f5ad-4b60-b1a6-c9321c5a38c6
# ╠═71e21c9d-95c0-4397-929d-b415b0720dbf
# ╠═bce24033-3bca-42d4-9224-b40378cdedab
# ╠═88d507e8-98cb-42ca-be20-16580f5beb4b
# ╟─7258550c-c49a-4f67-83f7-58c157d4b5de
# ╠═357a7c6c-2d6a-40e5-ba85-7042ecc9bdda
# ╠═58f03667-5903-48d7-a962-9fbbf5d72620
# ╠═af7b55de-8ae6-4b19-80b8-dcfe5cae276f
# ╠═9f623db7-32db-43ce-9ddd-98d5ee4248ba
# ╟─0a1d1783-0c24-455e-84c1-027e168cfef7
# ╠═6d8cfb4c-5602-4112-a691-92a319f5e224
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
# ╠═0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
# ╠═594fbbef-40c8-4cb8-9d56-14c39c65a914
# ╠═095ed289-8086-4da3-aaeb-36dec2b4d349
# ╠═a7b54b66-9d44-46ab-a4c2-d9326acccda5
# ╠═7564bd74-d182-4d89-bf57-2d0fe286bc7a
# ╠═bddcd6bb-5dba-4b6e-ba07-e21df433f812
# ╠═d415fac3-d155-44fc-86af-35e6ead45734
# ╠═b27ae510-f567-4623-a2a8-d4f192d2f792
# ╠═b7fd5546-ecc9-444b-9862-1954c3c762e1
# ╠═c954772f-9e75-484e-bfdb-4455249b1f99
# ╠═21f0d8c8-76a2-41f0-a390-390c5b539ed3
# ╠═8ea50ebf-87d3-4810-9699-89617463d9b4
# ╠═21373373-2dec-4bd8-a178-6af90df0835a
# ╠═51eea365-1dd0-4ebf-8a89-30928c51ab5f
# ╠═5a537690-3c29-4a0b-9a31-15d19e456db6
# ╠═3b31bf9e-1d07-460a-80d9-75625d06bf41
# ╠═756bb316-6afa-4217-b52b-eff39e5c02d6
# ╠═c04cab38-6e6e-4097-a82f-91cc8a3864f1
# ╠═602f9f99-987b-4557-9709-b42025e4f54e
# ╠═3051cdfb-63fb-44e4-a75d-37a51fd40508
# ╠═15e5e568-560d-4a70-8ca3-b084d6d9c832
# ╠═3f44fb3f-865a-4acb-9101-283c7918e11b
# ╠═58ea3ff0-40ab-4967-aa97-80d4aa903183
# ╠═2c2e54c6-4d2f-40b1-b012-5d21dc3bd4b9
# ╠═82e8db20-91cc-43b4-8183-7e302ce74041
# ╠═25514746-3762-4cb0-8b03-0e3e38d82a69
# ╠═6fb045c9-fce3-4d24-8edc-8ddde135972a
# ╠═ebc0bd82-2b1e-4653-8aa2-46fbb56952ed
# ╠═7e544a72-62cf-4e0d-a23c-cba001a1e2fd
# ╠═d9201257-2e2f-4a21-8bb2-97ff349fe04b
# ╠═5f321550-4d2e-4f6c-8dbf-9ec1c72f1a44
# ╠═9326b0d6-847d-49fe-bfe3-380f328d5eb0
# ╠═c6140b22-0744-456e-a0fe-164d72c577af
