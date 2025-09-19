### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 4c0da6d9-c5d9-4971-8b69-7ec9c107d981
begin
	using Pkg; Pkg.activate()

	using CSV, DataFrames
	using CairoMakie, GeoMakie
	
	using LilGuys, Arya

	import SkyCoords
	using OrderedCollections
	using PyFITS
end

# ╔═╡ 7216773a-bf24-49bb-a4ab-4afb35ce2a2b
md"""
# Setup
"""

# ╔═╡ 97993507-1ab6-4aae-8363-f439fb98c2e5
FIGDIR = "./figures"

# ╔═╡ 0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
CairoMakie.activate!(type=:png)

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

# ╔═╡ 45499b5b-2c79-4a0e-8cc8-930581a8dec7
function to_galcoords(ra, dec)
	coords = SkyCoords.ICRSCoords.(deg2rad.(ra), deg2rad.(dec))
	galcoords = SkyCoords.convert.(SkyCoords.GalCoords, coords)

	l = SkyCoords.lon.(galcoords) .|> rad2deg
	b = SkyCoords.lat.(galcoords) .|> rad2deg

	return l, b
end

# ╔═╡ 89b6a8ba-cb35-436d-a6b0-629505e46808
function to_galcoords(ra::Real, dec::Real)
	coords = SkyCoords.ICRSCoords(deg2rad(ra), deg2rad(dec))
	galcoords = SkyCoords.convert(SkyCoords.GalCoords, coords)

	l = SkyCoords.lon(galcoords) |> rad2deg
	b = SkyCoords.lat(galcoords) |> rad2deg

	return l, b
end

# ╔═╡ 11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
function read_pace(name)
	CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/pace/$name.csv"), DataFrame)
end

# ╔═╡ 516eb72e-44e6-4637-ba75-8b0c5ce6f12a
 dwarfs_pace = read_pace("dwarf_mw")

# ╔═╡ 17cc0a5f-6fc5-4d62-82fc-f4b207cb2b66
allcluster_pace = vcat(read_pace("gc_harris"), read_pace("gc_mw_new"))

# ╔═╡ edb5992a-9c96-4540-843d-9800a392a40e
harris = read_pace("gc_harris")

# ╔═╡ 22b55e07-9958-4515-b242-591338fa08e2
md"""
## Pace et al. (2024)
"""

# ╔═╡ b71bdb56-608d-4590-90fd-55dec7480077
alldwarf_pace = read_pace("dwarf_mw")

# ╔═╡ a232715d-a27f-499b-8123-9adf5257acf5
ambiguous_pace = read_pace("gc_ambiguous")

# ╔═╡ 88d507e8-98cb-42ca-be20-16580f5beb4b
read_pace("gc_other")

# ╔═╡ 83b25227-6b56-45f0-92e9-0dca1cf5eb2d
md"""
## Yasone et al. (2025)
"""

# ╔═╡ c3442639-f5ad-4b60-b1a6-c9321c5a38c6
yasone_halo = let
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

# ╔═╡ 20138d1b-7817-4732-9e0a-9880b5dc478e
yasone_disk = let
	df = CSV.read(
		Vector{UInt8}("""name,ra,dec
		Yasone 5,95.88220,20.76780
		Yasone 6,93.99560,10.87850
		Yasone 7,12.42938,70.43137
		Yasone 8,330.78870,60.75680
		Yasone 9,316.50120,45.98940
		Yasone 10,302.37750,40.73327
		Yasone 11,298.09783,29.94878
		Yasone 12,295.19450,30.15050
		Yasone 13,338.52182,57.71902
		Yasone 14,302.60948,40.92500
		Yasone 15,100.54790,11.73720
		Yasone 16,327.23950,59.01930
		Yasone 17,331.34584,55.14620
		Yasone 18,316.83027,45.91032
		Yasone 19,298.82535,35.49282
		Yasone 20,296.25171,16.28925
		"""), DataFrame
	)

		df[!, :key] = lowercase.(replace.(df.name, " "=>"_"))

	df[!, "ll"], df[!, "bb"] = to_galcoords(df.ra, df.dec)


	df
end

# ╔═╡ 0a1d1783-0c24-455e-84c1-027e168cfef7
md"""
# Filtering & classification
"""

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

# ╔═╡ 1f6137df-7514-4d5f-9724-5dbf6788a9f8
alldwarfs = dwarfs_pace# [dwarfs_pace.confirmed_dwarf .=== 1, :]

# ╔═╡ f6700ba9-7878-4a5b-bfa5-2f760dcaff5a
ambiguous = ambiguous_pace

# ╔═╡ 6df26b0d-12f1-41e2-9345-8d0f1225f11d
allcluster = harris

# ╔═╡ 6669bc73-efe4-49ba-999b-e1efcb3b6df8
CSV.write("dwarf_galaxies.csv", alldwarfs)

# ╔═╡ f48f3fed-f0ae-41ae-b3a0-83ca90eef2ad
CSV.write("globular_clusters.csv", allcluster)

# ╔═╡ fde48a6b-3311-484e-87df-6b5dff0a2c18
CSV.write("ambiguous.csv", ambiguous)

# ╔═╡ fa993dfa-07e3-4355-b29f-7d43bb6c42e9
CSV.write("andromida.csv", andromida)

# ╔═╡ 2a76ae4e-4268-40a8-b2cc-af98e4f0b8d9
CSV.write("yasone_halo.csv", yasone_halo)

# ╔═╡ 602a8e2d-575c-4341-a3f3-33a052a43981
CSV.write("yasone_disk.csv", yasone_disk)

# ╔═╡ 594fbbef-40c8-4cb8-9d56-14c39c65a914
allsatalites = vcat(alldwarfs, allcluster, ambiguous, andromida, yasone_halo,yasone_disk, cols=:union)

# ╔═╡ 2fd079d4-7fb1-4f2a-9f8a-e2c361d17137
md"""
# Extra datasets
"""

# ╔═╡ 7aa7fed2-2c6e-4db9-883e-a74d228f5e57
md"""
## McConnachie 2012
"""

# ╔═╡ 276edcfe-de8f-4f1f-9051-677b26bab2f2
m12_old_group = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", "all", "mcconnachie2012", "table1.dat"), DataFrame, stripwhitespace=true)

# ╔═╡ 8c0817e6-5a1e-4c00-a5b8-df9124baf0f1
m12_old_vel = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", "all", "mcconnachie2012", "table2.dat"), DataFrame, stripwhitespace=true, skipto=3, header=1, delim="|", types=Dict(:distance=>Float64, :distance_mw=>Float64, :velocity_mw=>Float64))

# ╔═╡ d4a142a7-d5e1-4ac6-acf8-86271ab3315d
m12_old_struct = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", "all", "mcconnachie2012", "table3.dat"), DataFrame, stripwhitespace=true,delim="|", types=Dict(:V=>Float64, :r_h => Float64, :M_V=>Float64, :ellipticity=>Float64))

# ╔═╡ 61cf738e-516d-4557-9296-787b5221ddbe
m12 = let 
	df = leftjoin(m12_old_group, m12_old_struct, on="name", makeunique=true)
	df = leftjoin(df, m12_old_vel, on="name"=>"Name", makeunique=true)
	df[!, :ra] = @. 360 / 24 * (df.ra_h + df.ra_m/60 + df.ra_s /60^2)
	df[!, :dec] = @. ifelse(df.dec_sign == Ref("+"), 1, -1) * (df.dec_deg + df.dec_min / 60 + df.dec_sec/60^2)

	df[!, :apparent_magnitude_v] = df.M_V .+ df.distance_modulus
	df[!, :rhalf] = ifelse.(ismissing.(df.r_h), NaN, df.r_h)
	df[!, :distance] = ifelse.(ismissing.(df.distance), NaN, df.distance)
	df[!, :rhalf_physical] = LilGuys.arcmin2kpc.(df.rhalf, df.distance) .* 1e3
	
	df
end

# ╔═╡ 58f03667-5903-48d7-a962-9fbbf5d72620
md"""
## Baumgart cluster catalogue
"""

# ╔═╡ b29048cd-f988-4bc7-8d71-3b659570a87c
md"""
To reproduce list in Baumgardt and Vasiliev 2021, we need to :
exclude the following 
-"Gran_1"
- "Gran_2"
- "Gran_3"
- "Gran_5"
- Patchick_126
- VVV-CL160

Include 
- Mercer 5?
"""

# ╔═╡ 4338b9a5-d6ed-4402-a4fb-f0d105c2c5ac
baumgardt_exclude = ["Gran_1", "Gran_2", "Gran_3", "Gran_5", "Patchick_126", "VVV-CL160"]

# ╔═╡ af7b55de-8ae6-4b19-80b8-dcfe5cae276f
baumgardt_columns = string.(split("Cluster         RA        DEC      R_Sun  DRSun    R_GC   DRGC   N_RV   N_PM    Mass         DM        V  Delta_V  M/L_V  DM/L    rc     rh_l     rh_m     rt    rho_c  rho_h_m   sig_c   sig_h_m lg(Trh) lg(Mini) T_Diss  M_Low  M_High    MF  Delta_MF  sig0    vesc   etac   etah  A_Rot Delta_AR P_Rot", r"\s+"))

# ╔═╡ c84d678c-f502-4e2f-a08b-3ccce27b17f4
mercer_5 = let
	df = Dict("Cluster"=> "Mercer_5", "RA"=>275.832, "DEC"=>−13.669)

	df["ll"], df["bb"] = to_galcoords(df["RA"], df["DEC"])

	df
end

# ╔═╡ 9f623db7-32db-43ce-9ddd-98d5ee4248ba
baumgardt = let
	df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/baumgardt_23.txt"), DataFrame, comment="#", ignorerepeated=true, delim=' ', header=baumgardt_columns)

	df[!, "ll"], df[!, "bb"] = to_galcoords(df.RA, df.DEC)

	df[!, :rhalf_physical] = df.rh_l
	df[!, :M_V] = 4.83 .- 5/2*log10.((df.Mass ./ df.:var"M/L_V" ))
	df[!, :apparent_magnitude_v] = df.V
	df[!, :rhalf] = LilGuys.kpc2arcmin.(df.rh_l./1e3, df.R_Sun)
	df

end;

# ╔═╡ d90f944d-f05e-4a7a-9450-5b235fe42d85
baumgard2021 = let
	df = baumgardt[baumgardt.Cluster .∉ Ref(baumgardt_exclude), :]
	append!(df, mercer_5, cols=:union)
	df[!, :key] = df.Cluster
	df[!, :name] = df.Cluster
	df
end

# ╔═╡ 216e636a-da55-48ca-ba37-e14b6c2509d8
md"""
## Ambiguous systems
"""

# ╔═╡ 7763389d-56c1-4fc8-9116-9e5025d1feaa
ambiguous_smith_keys = [
	"koposov_1", 
	"koposov_2",
	"segue_3",
	"munoz_1",
	"balbinot_1",
	"kim_1",
	"kim_2",
	"crater_1",
	"laevens_3",
	"draco_2",
	"eridanus_3",
	"pictor_1",
	"smash_1",
	"kim_3",
	"des_1",
	"des_J0111-1341",
	"des_J0225+0304",
	"des_3",
	"des_4",
	"des_5",
	"gaia_03",
	"ps1_1",
	"to_1",
	"bliss_1",
	"hsc_1",
	"delve_1",
	"delve_2",
	"ymca_1",
	"delve_3",
	"delve_4",
	"delve_5",
	"delve_6",
	"ursa_major_3"
]

# ╔═╡ 2694ffb9-f5d3-440d-980e-0e3ac1d67e61
ambiguous_smith = vcat(
	ambiguous_pace[ambiguous_pace.key .∈ Ref(ambiguous_smith_keys), :],
	allcluster_pace[allcluster_pace.key .∈ Ref(ambiguous_smith_keys), :],
	alldwarf_pace[alldwarf_pace.key .∈ Ref(ambiguous_smith_keys), :],
)

# ╔═╡ ff0a76c1-5845-43c8-b364-0986b126618b
ambiguous_smith_keys[ambiguous_smith_keys .∉ Ref(ambiguous_smith.key)]

# ╔═╡ d6769a11-de42-496a-8a4f-45f35238ac5b
md"""
## McConnachie 2021
"""

# ╔═╡ 6ecada0c-3b62-4965-a25f-62ebf3e0a80f
mcconnachie_v21 = let
	df = read_fits(joinpath(ENV["DWARFS_ROOT"], "observations/all/MV2020/NearbyGalaxies_Jan2021_PUBLIC.fits"))

	lb = [to_galcoords(LilGuys.sexadecimal_to_decimal(row.RA, row.Dec)...) for row in eachrow(df)]
	df[!, :ll] = first.(lb)
	df[!, :bb] = last.(lb)

	df[!, :r_gal] = map(eachrow(df)) do row	
		r = LilGuys.dm2kpc(row.dmod)
		ra, dec = LilGuys.sexadecimal_to_decimal(row.RA, row.Dec)
		gc = LilGuys.transform(Galactocentric, LilGuys.ICRS(ra=ra,dec=dec, distance=r),)
		d_gc = radii(LilGuys.position(gc))
	end

	df[:, :key] = [strip(replace(name, "(1)"=>"", "*"=>"")) for name in df.GalaxyName]

	df[:, :name] = [replace(name, "1"=>" I", "2"=>" II", "3"=>" III", "4"=>" IV", "5"=>" V", 
								 "dSph"=>"", 
								 "CanesVenatici"=>"Canes Venatici", "ComaBerenices"=>"Coma Berenices",
								"CanisMajor"=>"Canis Major", "UrsaMajor"=> "Ursa Major", "UrsaMinor"=>"Ursa Minor" ) for name in df.key]


	df[:, :distance] = LilGuys.dm2kpc.(df.dmod)
	df[:, :rhalf] = df.rh
	df[:, :rhalf_physical] = 1e3 * LilGuys.arcmin2kpc.(df.rhalf, df.distance)
	df[:, :M_V] = df.Vmag .- df.dmod
	df[:, :apparent_magnitude_v] = df.Vmag
	df
end

# ╔═╡ edac377b-0018-4645-9379-a7eedb2bb59b
function escape_velocity(R_gc)
	M = 100 # 10e12 solar mass
	return V2KMS * sqrt(2*M / R_gc)
end

# ╔═╡ 8be09a2b-8edc-4594-be37-66cc99409e46
mcconnachie_is_ambiguous = mcconnachie_v21.key .∈ Ref(["DESJ0225+0304", "Draco2"])

# ╔═╡ 33126032-0db6-4c60-bca2-f32eedd737b9
escape_velocity(180)

# ╔═╡ f27c08a6-01ac-4a8c-8d9e-854c12b2e80f
mcconnachie_is_galaxy = .!startswith.(mcconnachie_v21.GalaxyName, ["*"])

# ╔═╡ 2f9cdd0e-e7f1-4072-923a-548168e07ca9
mcconnachie_is_bound = map(eachrow(mcconnachie_v21)) do row
	d_gc = row.r_gal
	rv = row.vh

	if row.name ∈ ["Leo 1", "Antlia 2"]
		return true
	end
	if (rv == 999) .& (d_gc .< 300)
		return true
	end
	return (abs.(rv) .< escape_velocity(d_gc)) .& (d_gc .< 300)
end

# ╔═╡ aab3e114-9e7a-47fb-b93b-6fe0c34caf45
mcconnachie_v21_mw = mcconnachie_v21[mcconnachie_is_bound .& .! mcconnachie_is_ambiguous, :]

# ╔═╡ 095ed289-8086-4da3-aaeb-36dec2b4d349
md"""
# Plot utils
"""

# ╔═╡ 8f150143-8304-4729-81ed-b04115d03842
styles_mag = [
	(; markersize=5, alpha=0.5, color=COLORS[4]),
	(; color=COLORS[5], markersize=5, marker=:hexagon, alpha=0.5),
	(; color=COLORS[1], markersize=5, alpha=0.5, marker=:rect)
]

# ╔═╡ b9ad6544-ac60-4cc2-a86b-9ba78605ae9a
function plot_mag_size(samples; title="")

	fig = Figure(size=(3*72, 6*72))
	ax = Axis(fig[1,1], xscale = log10, xticks=Makie.automatic,
			xlabel = "rh / pc",
			  ylabel = "absolute magnitude Mv",
			  yreversed=true,
			  limits=(0.4, 500, -12.5, 5),
			  title=title
		)

	
	for (i, (label, sample)) in enumerate(samples)
		df = dropmissing(sample, [:rhalf_physical, :M_V])
		scatter!(df.rhalf_physical, df.M_V; label=label, styles_mag[i]...)

		@info "$label: $(size(df, 1)) rows"
	end
	axislegend(position=:rb, backgroundcolor=(:white, 0.8))

	ax = Axis(fig[2,1], xscale = log10, 
			xlabel = "rh / arcmin",
			  ylabel = "V",
			  yreversed=true,
			  limits=(0.15, 40, 3, 20),
			  xticks = [0.25, 0.5, 1, 2, 5, 10, 20, 40],
			  
		)

	
	
	for (i, (label, sample)) in enumerate(samples)
		df = dropmissing(sample, [:rhalf, :apparent_magnitude_v])

		@info "$label: $(size(df, 1)) rows"

		scatter!(df.rhalf, df.apparent_magnitude_v; label=label, styles_mag[i]...)

	end


	return fig
end

# ╔═╡ 7dda1ac1-a733-4b2f-b25e-9125f1737ff2
@savefig "magnitude_size_m12" plot_mag_size(OrderedDict(
	"baumgardt" => baumgardt,
	"ambiguous" => ambiguous,
	"M12 all" => m12
))

# ╔═╡ a45d2884-9a1f-42d3-933a-656e95855bf7
@savefig "magnitude_size_m12_mw" plot_mag_size(OrderedDict(
	"baumgardt" => baumgardt,
	"ambiguous" => ambiguous,
	"M12 MW only" => m12[m12.subgroup .== "MW", :]
))

# ╔═╡ 2ae52446-2741-46fb-8ed5-35b88253eeaf
@savefig "magnitude_size_pace" plot_mag_size(OrderedDict(
	"pace clusters" => allcluster,
	"ambiguous" => ambiguous,
	"pace dwarfs" => alldwarfs,
))

# ╔═╡ ff913a4f-b0ed-46ec-86fc-eecfb1106bdd
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xscale=log10, yreversed=true)
		
	
	df = alldwarfs[sortperm(alldwarfs.rhalf, rev=true), :]


	scatter!(df.rhalf[1:20], df.apparent_magnitude_v[1:20])
	text!(df.rhalf[1:20], df.apparent_magnitude_v[1:20], text=df.name[1:20], fontsize=5)

	fig
end

# ╔═╡ 9ee8283e-72a7-45e2-931b-7dd4f288fc7f
alldwarfs.apparent_magnitude_v

# ╔═╡ a7cdbf5d-fec0-4895-86be-102a42e2bbca
plot_mag_size(OrderedDict(
	"pace clusters" => allcluster,
	"ambiguous" => ambiguous,
	"pace dwarfs" => alldwarfs[alldwarfs.confirmed_dwarf .=== 1, :],
))

# ╔═╡ 35d385a4-b0f8-4476-a6e6-1b906927733d
@savefig "magnitude_size_m12_v21" plot_mag_size(OrderedDict(
	"pace clusters" => harris,
	"ambiguous" => ambiguous,
	"m21" => mcconnachie_v21,
))

# ╔═╡ Cell order:
# ╟─7216773a-bf24-49bb-a4ab-4afb35ce2a2b
# ╠═4c0da6d9-c5d9-4971-8b69-7ec9c107d981
# ╠═97993507-1ab6-4aae-8363-f439fb98c2e5
# ╠═0ee0aa61-3046-4dae-92e9-bd6b124d8f5f
# ╟─5c124928-d05d-467a-93f9-dc4e9e3b1d04
# ╠═f3881abc-9fb3-4b9a-89f6-3790c21fc6ed
# ╠═9179b4ed-eb7f-4a27-8719-fc028078a8c6
# ╟─483e5816-7f41-476f-9696-096ce83d4340
# ╠═45499b5b-2c79-4a0e-8cc8-930581a8dec7
# ╠═89b6a8ba-cb35-436d-a6b0-629505e46808
# ╠═11d8d92e-9ba7-48a4-8c38-69e0f4b5d14a
# ╠═516eb72e-44e6-4637-ba75-8b0c5ce6f12a
# ╠═17cc0a5f-6fc5-4d62-82fc-f4b207cb2b66
# ╠═edb5992a-9c96-4540-843d-9800a392a40e
# ╟─22b55e07-9958-4515-b242-591338fa08e2
# ╠═b71bdb56-608d-4590-90fd-55dec7480077
# ╠═a232715d-a27f-499b-8123-9adf5257acf5
# ╠═88d507e8-98cb-42ca-be20-16580f5beb4b
# ╟─83b25227-6b56-45f0-92e9-0dca1cf5eb2d
# ╠═c3442639-f5ad-4b60-b1a6-c9321c5a38c6
# ╠═20138d1b-7817-4732-9e0a-9880b5dc478e
# ╟─0a1d1783-0c24-455e-84c1-027e168cfef7
# ╠═a09fe200-9239-4e74-a8c4-c548e67481f4
# ╠═1f6137df-7514-4d5f-9724-5dbf6788a9f8
# ╠═f6700ba9-7878-4a5b-bfa5-2f760dcaff5a
# ╠═6df26b0d-12f1-41e2-9345-8d0f1225f11d
# ╠═6669bc73-efe4-49ba-999b-e1efcb3b6df8
# ╠═f48f3fed-f0ae-41ae-b3a0-83ca90eef2ad
# ╠═fde48a6b-3311-484e-87df-6b5dff0a2c18
# ╠═fa993dfa-07e3-4355-b29f-7d43bb6c42e9
# ╠═2a76ae4e-4268-40a8-b2cc-af98e4f0b8d9
# ╠═602a8e2d-575c-4341-a3f3-33a052a43981
# ╠═594fbbef-40c8-4cb8-9d56-14c39c65a914
# ╟─2fd079d4-7fb1-4f2a-9f8a-e2c361d17137
# ╟─7aa7fed2-2c6e-4db9-883e-a74d228f5e57
# ╠═276edcfe-de8f-4f1f-9051-677b26bab2f2
# ╠═8c0817e6-5a1e-4c00-a5b8-df9124baf0f1
# ╠═d4a142a7-d5e1-4ac6-acf8-86271ab3315d
# ╠═61cf738e-516d-4557-9296-787b5221ddbe
# ╟─58f03667-5903-48d7-a962-9fbbf5d72620
# ╟─b29048cd-f988-4bc7-8d71-3b659570a87c
# ╠═4338b9a5-d6ed-4402-a4fb-f0d105c2c5ac
# ╠═af7b55de-8ae6-4b19-80b8-dcfe5cae276f
# ╠═c84d678c-f502-4e2f-a08b-3ccce27b17f4
# ╠═9f623db7-32db-43ce-9ddd-98d5ee4248ba
# ╠═d90f944d-f05e-4a7a-9450-5b235fe42d85
# ╟─216e636a-da55-48ca-ba37-e14b6c2509d8
# ╟─7763389d-56c1-4fc8-9116-9e5025d1feaa
# ╠═2694ffb9-f5d3-440d-980e-0e3ac1d67e61
# ╠═ff0a76c1-5845-43c8-b364-0986b126618b
# ╟─d6769a11-de42-496a-8a4f-45f35238ac5b
# ╠═6ecada0c-3b62-4965-a25f-62ebf3e0a80f
# ╠═edac377b-0018-4645-9379-a7eedb2bb59b
# ╠═aab3e114-9e7a-47fb-b93b-6fe0c34caf45
# ╠═8be09a2b-8edc-4594-be37-66cc99409e46
# ╠═33126032-0db6-4c60-bca2-f32eedd737b9
# ╠═f27c08a6-01ac-4a8c-8d9e-854c12b2e80f
# ╠═2f9cdd0e-e7f1-4072-923a-548168e07ca9
# ╟─095ed289-8086-4da3-aaeb-36dec2b4d349
# ╠═8f150143-8304-4729-81ed-b04115d03842
# ╠═b9ad6544-ac60-4cc2-a86b-9ba78605ae9a
# ╠═7dda1ac1-a733-4b2f-b25e-9125f1737ff2
# ╠═a45d2884-9a1f-42d3-933a-656e95855bf7
# ╠═2ae52446-2741-46fb-8ed5-35b88253eeaf
# ╠═ff913a4f-b0ed-46ec-86fc-eecfb1106bdd
# ╠═9ee8283e-72a7-45e2-931b-7dd4f288fc7f
# ╠═a7cdbf5d-fec0-4895-86be-102a42e2bbca
# ╠═35d385a4-b0f8-4476-a6e6-1b906927733d
