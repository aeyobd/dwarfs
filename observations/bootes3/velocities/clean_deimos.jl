### A Pluto.jl notebook ###
# v0.20.23

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

# ╔═╡ 2f5ca80c-d98d-45cc-a661-834002b620f6
using PairPlots

# ╔═╡ 23b3766e-15b7-4b1e-bb41-d6af36a59caf
using PyFITS

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
using Turing

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Deimos

This notebook analyzes the RV stars in from Geha+2026's reprocessing of legacy Keck/DEIMOS observations of Boo 3 (from Carlin et al. 2019) and creates a sample sutable to be combined with others for velocity analysis.

"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

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
raw_table = CSV.read("../data/deimos_velocities.txt", DataFrame, ignorerepeated=true, delim=" ", skipto=36, 
	header=["system", "ra", "dec", "r_mag", "g_r", "nmask", "texp", "SN", "RV",
			"RV_err", "WECaT", "EWCaT_err", "FE_H", "FE_H_err", "Var", "P_SAT_deimos"]
					)

# ╔═╡ 0272cf3b-4299-488a-9d63-4262813c617e
hist(raw_table.Var)

# ╔═╡ d5e1fc38-9b93-4024-9259-c9253914c72a
boo3_deimos = raw_table[raw_table[!, "system"] .== Ref("Boo3"), :]

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 0137b308-f0c6-4e73-afa6-346797b6f304
j24 = read_fits("$data_dir/j24_2c.fits") 

# ╔═╡ 2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
xmatch_max = 2.0

# ╔═╡ 0ed1ac33-0537-4838-ae0b-200be73c501c
md"""
Depends on:
- `../data/j24_2c.fits`
- `../observed_properties.toml`
- `../data/deimos_velocities.txt`
Creates:
- `processed/rv_deimos.fits`
- `processed/mcmc_summary_uncorrected_spencer+18.csv`

xmatch by nearest star within $xmatch_max arcmin.
"""

# ╔═╡ d321c481-0ceb-4e3b-b5fc-12af735155e3
filt_xmatch, idx_xmatch = RVUtils.xmatch(boo3_deimos, j24, xmatch_max)

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
	df = copy(boo3_deimos)
	df[!, :source_id] .= source_id
	df[!, :F_match] = F_match
	df[!, :F_scatter] .= true
	# df[:, :RV_err_orig] = copy(df.RV_err)
	# df[:, :RV_err] = @. sqrt(df.RV_err^2 *1.4^2 + 1.1^2 )

	df
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_deimos.fits", df_out, overwrite=true)

# ╔═╡ 4b8a45ce-4862-4842-8c60-bae9a64acd75
md"""
### Sanity checks with sourceid matching
"""

# ╔═╡ 423e446b-3cf6-4150-8b94-c4d7388fceb4
@assert length(unique(skipmissing(df_out.source_id))) == sum(filt_xmatch)

# ╔═╡ 011700b8-5dcc-409b-807c-952ea9726992
df_out[nonunique(df_out, :source_id, keep=:noduplicates), :]

# ╔═╡ 70729597-2e50-4773-89ad-87e97f7bdad6
length(unique(df_out.source_id)), sum(.!ismissing.(df_out.source_id)), sum(filt_xmatch)

# ╔═╡ 70dbd229-281a-4796-873d-4be2a4000e34
md"""
# Numbers
"""

# ╔═╡ 2612309e-f128-4ddf-96d8-91b86552de17
length(boo3_deimos.RV)

# ╔═╡ 8064bea8-ff54-4a71-807c-30f52f68f1f1
length(df_out.RV)

# ╔═╡ 7a48caaa-d953-4d1e-b7db-2d07d8773bf2
sum(F_match)

# ╔═╡ e11da237-4c8f-450f-bd41-b2c9cdf2f11d
sum(F_match)

# ╔═╡ 9ac3364e-5532-457e-b29b-bfb10d32e845
matched = j24[idx_xmatch, :]

# ╔═╡ 7962908c-4150-46a8-8388-2f5283b1d0c4
matched[matched.F_BEST .== 0.0, :]

# ╔═╡ 4290a51e-961b-4b38-b954-1a65e946080c
md"""
# Plots
"""

# ╔═╡ cee182e1-5545-4b86-9152-b288a9b11d86
hist(boo3_deimos.SN)

# ╔═╡ c0c02bce-e0f9-4fa0-bba6-a5360a825e23
xmatch_dists = lguys.angular_distance.(df_out.ra, df_out.dec, j24.ra[idx_xmatch], j24.dec[idx_xmatch]) * 60^2

# ╔═╡ d3e682ee-0dab-4aa5-9f3e-423a9d2f92b0
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "distance / arcsec",
			  yscale=log10,
			  yticks=Makie.automatic,
	)
	hist!(xmatch_dists, bins=0:0.05:3)
	vlines!(xmatch_max, color=:black )
	fig
end

# ╔═╡ c873ccf2-8ecf-4a85-a152-76acd1c1fc62
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "ra",
			  ylabel = "dec",
	)
	scatter!(j24.ra[idx_xmatch[filt_xmatch]], j24.dec[idx_xmatch[filt_xmatch]], markersize=2)

	scatter!(df_out.ra[filt_xmatch], df_out.dec[filt_xmatch], markersize=1,alpha=1, marker=:circ, color=COLORS[2])

	fig
end

# ╔═╡ 344f3c29-0873-4180-9025-51aeaeb2c681
boo3_members = let
	filt_memb = (boo3_deimos.P_SAT_deimos .> 0.5) .& (boo3_deimos.Var .<= 0)
	df = boo3_deimos[filt_memb, :]
	df
end

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis()
	bins = make_bins(boo3_deimos.RV, calc_limits(boo3_deimos.RV), bandwidth=1)

	h = histogram(boo3_deimos.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(boo3_members.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═0ed1ac33-0537-4838-ae0b-200be73c501c
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═2f5ca80c-d98d-45cc-a661-834002b620f6
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═23b3766e-15b7-4b1e-bb41-d6af36a59caf
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═36634dea-21bc-4823-8a15-7bce20b6fc17
# ╠═e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
# ╠═9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╠═d6a2ebb4-f125-4f28-9009-fc67d5dda5a1
# ╠═0272cf3b-4299-488a-9d63-4262813c617e
# ╠═d5e1fc38-9b93-4024-9259-c9253914c72a
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
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
# ╟─70dbd229-281a-4796-873d-4be2a4000e34
# ╠═2612309e-f128-4ddf-96d8-91b86552de17
# ╠═8064bea8-ff54-4a71-807c-30f52f68f1f1
# ╠═7a48caaa-d953-4d1e-b7db-2d07d8773bf2
# ╠═e11da237-4c8f-450f-bd41-b2c9cdf2f11d
# ╠═9ac3364e-5532-457e-b29b-bfb10d32e845
# ╠═7962908c-4150-46a8-8388-2f5283b1d0c4
# ╟─4290a51e-961b-4b38-b954-1a65e946080c
# ╠═cee182e1-5545-4b86-9152-b288a9b11d86
# ╠═c0c02bce-e0f9-4fa0-bba6-a5360a825e23
# ╠═d3e682ee-0dab-4aa5-9f3e-423a9d2f92b0
# ╠═c873ccf2-8ecf-4a85-a152-76acd1c1fc62
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
