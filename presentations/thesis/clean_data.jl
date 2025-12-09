### A Pluto.jl notebook ###
# v0.20.21

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
alldwarfs = dwarfs_pace[dwarfs_pace.confirmed_dwarf .=== 1, :]

# ╔═╡ f6700ba9-7878-4a5b-bfa5-2f760dcaff5a
ambiguous = ambiguous_pace

# ╔═╡ 6df26b0d-12f1-41e2-9345-8d0f1225f11d
allcluster = allcluster_pace

# ╔═╡ 6669bc73-efe4-49ba-999b-e1efcb3b6df8
CSV.write("dwarf_galaxies.csv", alldwarfs)

# ╔═╡ f48f3fed-f0ae-41ae-b3a0-83ca90eef2ad
CSV.write("globular_clusters.csv", allcluster)

# ╔═╡ fde48a6b-3311-484e-87df-6b5dff0a2c18
CSV.write("ambiguous.csv", ambiguous)

# ╔═╡ fa993dfa-07e3-4355-b29f-7d43bb6c42e9
CSV.write("andromida.csv", andromida)

# ╔═╡ 594fbbef-40c8-4cb8-9d56-14c39c65a914
allsatalites = vcat(alldwarfs, allcluster, ambiguous, andromida, cols=:union)

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
# ╟─0a1d1783-0c24-455e-84c1-027e168cfef7
# ╠═a09fe200-9239-4e74-a8c4-c548e67481f4
# ╠═1f6137df-7514-4d5f-9724-5dbf6788a9f8
# ╠═f6700ba9-7878-4a5b-bfa5-2f760dcaff5a
# ╠═6df26b0d-12f1-41e2-9345-8d0f1225f11d
# ╠═6669bc73-efe4-49ba-999b-e1efcb3b6df8
# ╠═f48f3fed-f0ae-41ae-b3a0-83ca90eef2ad
# ╠═fde48a6b-3311-484e-87df-6b5dff0a2c18
# ╠═fa993dfa-07e3-4355-b29f-7d43bb6c42e9
# ╠═594fbbef-40c8-4cb8-9d56-14c39c65a914
