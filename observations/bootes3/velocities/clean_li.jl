### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using PythonCall
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 23b3766e-15b7-4b1e-bb41-d6af36a59caf
using PyFITS

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Deimos

This notebook analyzes the RV stars in from Li+2026's new S5 observations.

"""

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ 36634dea-21bc-4823-8a15-7bce20b6fc17
import TOML

# ╔═╡ 9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../data/"

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Data loading
"""

# ╔═╡ d6a2ebb4-f125-4f28-9009-fc67d5dda5a1
raw_table = CSV.read("../data/li+2026.txt", DataFrame, ignorerepeated=true, delim=" ")

# ╔═╡ d5e1fc38-9b93-4024-9259-c9253914c72a
boo3_li = raw_table

# ╔═╡ 0137b308-f0c6-4e73-afa6-346797b6f304
j24 = read_fits("$data_dir/j24_2c.fits") 

# ╔═╡ 2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
xmatch_max = 2.0

# ╔═╡ 0ed1ac33-0537-4838-ae0b-200be73c501c
md"""
Depends on:
- `../data/j24_2c.fits`
- `../observed_properties.toml`
- `../data/li+2026.txt`
Creates:
- `processed/rv_li.fits`

xmatch by nearest star within $xmatch_max arcmin.
"""

# ╔═╡ d321c481-0ceb-4e3b-b5fc-12af735155e3
filt_xmatch, idx_xmatch = RVUtils.xmatch(boo3_li, j24, xmatch_max)

# ╔═╡ 96e1413d-e4fd-48af-8cde-73a5fcb4976a
begin 
	source_id = allowmissing(j24.source_id[idx_xmatch])
	source_id[.!filt_xmatch] .= missing
end

# ╔═╡ a6f7f1a1-2dad-4397-9565-7b2aaf1fa733
F_best = RVUtils.get_f_best.([j24], source_id)

# ╔═╡ 412148d2-3b14-400a-b93b-f4aa18965d8d
F_match = F_best .=== 1.

# ╔═╡ 28a22929-89f2-422d-9cf0-b06d7e45d9a4
df_out = let
	df = copy(boo3_li)
	df[!, :source_id_li] = df.source_id
	df[!, :source_id] .= source_id
	df[!, :F_match] = F_match
	df[!, :F_scatter] .= true
	df
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_li.fits", df_out, overwrite=true)

# ╔═╡ 4b8a45ce-4862-4842-8c60-bae9a64acd75
md"""
### Sanity checks with sourceid matching
"""

# ╔═╡ 423e446b-3cf6-4150-8b94-c4d7388fceb4
@assert length(unique(skipmissing(df_out.source_id))) == sum(filt_xmatch)

# ╔═╡ 011700b8-5dcc-409b-807c-952ea9726992
@assert size(df_out[nonunique(df_out, :source_id, keep=:noduplicates), :], 1) == 0 "duplicated source IDs"

# ╔═╡ 70729597-2e50-4773-89ad-87e97f7bdad6
@assert length(unique(df_out.source_id)) == sum(.!ismissing.(df_out.source_id)) == sum(filt_xmatch) == length(boo3_li.RV) == 21 "missing stars in final DF"

# ╔═╡ 35397f56-7bc5-48be-a715-1b20c67acaad
@assert all(df_out.source_id_li .== df_out.source_id)

# ╔═╡ 70dbd229-281a-4796-873d-4be2a4000e34
md"""
# Numbers
"""

# ╔═╡ 9ac3364e-5532-457e-b29b-bfb10d32e845
matched = j24[idx_xmatch, :]

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═0ed1ac33-0537-4838-ae0b-200be73c501c
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═23b3766e-15b7-4b1e-bb41-d6af36a59caf
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═36634dea-21bc-4823-8a15-7bce20b6fc17
# ╠═9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╠═d6a2ebb4-f125-4f28-9009-fc67d5dda5a1
# ╠═d5e1fc38-9b93-4024-9259-c9253914c72a
# ╠═0137b308-f0c6-4e73-afa6-346797b6f304
# ╠═2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
# ╠═d321c481-0ceb-4e3b-b5fc-12af735155e3
# ╠═96e1413d-e4fd-48af-8cde-73a5fcb4976a
# ╠═a6f7f1a1-2dad-4397-9565-7b2aaf1fa733
# ╠═412148d2-3b14-400a-b93b-f4aa18965d8d
# ╠═28a22929-89f2-422d-9cf0-b06d7e45d9a4
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╟─4b8a45ce-4862-4842-8c60-bae9a64acd75
# ╠═423e446b-3cf6-4150-8b94-c4d7388fceb4
# ╠═011700b8-5dcc-409b-807c-952ea9726992
# ╠═70729597-2e50-4773-89ad-87e97f7bdad6
# ╠═35397f56-7bc5-48be-a715-1b20c67acaad
# ╠═70dbd229-281a-4796-873d-4be2a4000e34
# ╠═9ac3364e-5532-457e-b29b-bfb10d32e845
