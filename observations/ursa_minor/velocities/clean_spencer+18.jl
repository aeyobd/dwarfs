### A Pluto.jl notebook ###
# v0.20.6

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
# Spencer + 18

This notebook analyzes the RV stars in spencer+18 and creates a sample sutable to be combined with others for velocity analysis.

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

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 37e6d04b-0347-4f14-9b06-d657a03df225
spencer18_raw = read_fits("$data_dir/spencer+18_all_new.fits") |> x->rename(x, 
	:vlos => :RV,
	:e_vlos => :RV_err,
	:RAJ2000 => :ra,
	:DEJ2000 => :dec,
)


# ╔═╡ a3f87e29-6095-436f-9095-deb2d16b994f
md"""
Remove stars with significant skewness or kertosis, does not matter for provided table (prefiltered).
"""

# ╔═╡ aefa18bd-4203-4034-9ae9-9893fa8cc64d
filt_qual = abs.(spencer18_raw[:, "S-vlos"] .< 1) .&& (abs.(spencer18_raw[:, "K-vlos"] .- 3) .<= 1)

# ╔═╡ 4e39ac0f-60d0-4091-a5e8-8aa3977e9f39
sum(.!filt_qual)

# ╔═╡ 0137b308-f0c6-4e73-afa6-346797b6f304
j24 = read_fits("$data_dir/jensen+24_2c.fits") 

# ╔═╡ 153b1c3f-88a0-4b31-a2d7-65c084ca3a43
spencer18_all = let
	df = DataFrame()
	spencer18_raw[!, :ra] = round.(spencer18_raw.ra, digits=3)
	spencer18_raw[!, :dec] = round.(spencer18_raw.dec, digits=4)
	for group in groupby(spencer18_raw[filt_qual, :], [:ra, :dec])
		d = OrderedDict()
		d["ra"] = mean(group.ra)
		d["dec"] = mean(group.dec)

		xs = group.RV
		ws = @. 1/group.RV_err^2
		d["RV"] = lguys.mean(xs, ws)
		d["RV_err"] = 1/sqrt(sum(ws))
		d["RV_sigma"] = lguys.std(xs, ws)
		d["RV_count"] = length(xs)
		if length(xs) <= 1
			d["RV_sigma"] = NaN
		end
				
		append!(df, d, promote=true)

	end

	df
end

# ╔═╡ 2a6e8cda-a3cf-4824-9584-fa3d3c1da53a
sum(spencer18_all.RV_sigma .=== missing)

# ╔═╡ 2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
xmatch_max = 2.0

# ╔═╡ 0ed1ac33-0537-4838-ae0b-200be73c501c
md"""
Depends on:
- `../data/jensen+24_2c.fits`
- `../observed_properties.toml`
- `../data/spencer+18_all_new.fits`, `*_averaged.fits`, `*_repeated.fits`
Creates:
- `processed/rv_spencer+18.fits`
- `processed/mcmc_summary_uncorrected_spencer+18.csv`

xmatch by nearest star within $xmatch_max arcmin.
"""

# ╔═╡ d321c481-0ceb-4e3b-b5fc-12af735155e3
filt_xmatch, idx_xmatch = RVUtils.xmatch(spencer18_all, j24, xmatch_max)

# ╔═╡ 96e1413d-e4fd-48af-8cde-73a5fcb4976a
begin 
	source_id = allowmissing(j24.source_id[idx_xmatch])
	source_id[.!filt_xmatch] .= missing
end

# ╔═╡ 48ca64bd-fd0a-4e49-b281-37ac200ba3ff


# ╔═╡ f751b13c-51fc-47c5-af2f-e1188c5bba5f
p_chi2 = RVUtils.prob_chi2(spencer18_all)

# ╔═╡ a6f7f1a1-2dad-4397-9565-7b2aaf1fa733
F_scatter = @. (p_chi2 > 0.001) || isnan(p_chi2)

# ╔═╡ 761a7be6-420a-4428-addb-7947125596b3
(spencer18_all.RV_sigma ./ spencer18_all.RV_err)[.!F_scatter]

# ╔═╡ 25e927e1-773d-4d41-894b-1f4512f01233
function get_f_best(source_id)
	if ismissing(source_id)
		return missing
	end
	if source_id ∉ j24.source_id
		return missing
	end
	return j24.F_BEST[j24.source_id .== source_id] |> only
end

# ╔═╡ 412148d2-3b14-400a-b93b-f4aa18965d8d
F_match = get_f_best.(source_id) .=== 1.0

# ╔═╡ 28a22929-89f2-422d-9cf0-b06d7e45d9a4
df_out = let
	df = copy(spencer18_all)
	df[!, :source_id] .= source_id
	df[!, :F_match] = F_match
	df[!, :F_scatter] = F_scatter
	df
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_spencer+18.fits", df_out, overwrite=true)

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
sum(spencer18_all.RV_count), length(spencer18_raw.RV)

# ╔═╡ 8064bea8-ff54-4a71-807c-30f52f68f1f1
length(spencer18_all.RV)

# ╔═╡ 7a48caaa-d953-4d1e-b7db-2d07d8773bf2
sum(F_match)

# ╔═╡ 8b993569-173e-4326-a778-cf6239771447
sum(F_scatter)

# ╔═╡ e11da237-4c8f-450f-bd41-b2c9cdf2f11d
sum(F_match .& F_scatter)

# ╔═╡ 9ac3364e-5532-457e-b29b-bfb10d32e845
matched = j24[idx_xmatch, :]

# ╔═╡ 7962908c-4150-46a8-8388-2f5283b1d0c4
matched[matched.F_BEST .== 0.0, :]

# ╔═╡ 34b1fcc3-7355-4d0d-b0aa-8af0541e8651
mean(matched[matched.F_BEST .== 0.0, :F_CPAR])

# ╔═╡ 08967635-801f-44dc-bb7e-62abbe0f86e4
mean(matched[matched.F_BEST .== 0.0, :F_INGRID])

# ╔═╡ 1853d75f-65e9-4ad0-8a5e-68603d4ecddf
mean(matched[matched.F_BEST .== 0.0, :F_ASTROMETRIC])

# ╔═╡ 4290a51e-961b-4b38-b954-1a65e946080c
md"""
# Plots
"""

# ╔═╡ 8a176998-016d-4932-946f-eca78eb53d2b
md"""
## Histograms / exporatory
"""

# ╔═╡ cee182e1-5545-4b86-9152-b288a9b11d86
hist(spencer18_raw.S_N)

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

# ╔═╡ fe505b23-c610-4e44-98f4-fb4c9ffe0388
μ = obs_properties["radial_velocity"]

# ╔═╡ d3440c3e-be84-403f-a04d-68a237375b62
σ = obs_properties["sigma_v"]

# ╔═╡ 33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
filt_memb = μ - 3σ .< spencer18_all.RV .<  μ + 3σ

# ╔═╡ 344f3c29-0873-4180-9025-51aeaeb2c681
spencer18_members = spencer18_all[filt_memb, :]

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis()
	bins = make_bins(spencer18_all.RV, calc_limits(spencer18_all.RV), bandwidth=1)

	h = histogram(spencer18_all.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(spencer18_members.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ 1a75bca0-c305-47a0-959e-afa42bd75d11
md"""
## Probabilities
"""

# ╔═╡ 5774a0ef-8abf-4b6d-9b34-f80c29b4aab4
F_scatter_old = (spencer18_all.RV_count .< 2) .|| (spencer18_all.RV_sigma .< 3*spencer18_all.RV_err .* sqrt.(spencer18_all.RV_count))

# ╔═╡ 05e35f22-e682-4ec6-acf4-1d7bd9a4ca79
md"""
Histogram below should be uniform + excess at 0.
"""

# ╔═╡ d86097ad-2feb-46d8-9ceb-a2743a74a1a9
hist(filter(isfinite, p_chi2), 
	 axis = (; xlabel = L"P(<\chi^2)")
	)

# ╔═╡ 5de87579-1d6e-4dc0-860b-501816d75f8c
hist(filter(isfinite, spencer18_members.RV_sigma ./ spencer18_members.RV_err))

# ╔═╡ 214ed7f8-f338-4186-8f6d-870498486c0a
md"""
## MCMC
"""

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(spencer18_members.RV, spencer18_members.RV_err)

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = RVUtils.summarize(chain)

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(spencer18_members, samples, bins=30)

# ╔═╡ 60f45b81-8beb-4eb0-8f55-db04c5368eee
CSV.write("processed/mcmc_summary_uncorrected_spencer+18.csv", df_summary)

# ╔═╡ cb2923b2-dfda-4257-8d35-cd8f28188dbe
md"""
# other dataframe
"""

# ╔═╡ 28ff5e47-468d-4b54-bf24-d822818ec93d
md"""
Spencer+18 publishes a few other dataframes, including mixed data and averages across studies.
"""

# ╔═╡ 9847f6c4-05bd-481c-894d-4966d37be8c8
spencer_mixed_raw = read_fits(joinpath(data_dir, "spencer+18_repeated.fits")) 

# ╔═╡ fdfde937-39e5-4d51-8eb9-b2ec1ddd7dbd
spencer_mixed = let
	df = DataFrame()
	for group in groupby(spencer_mixed_raw, :ID)
		d = OrderedDict()
		d["ID"] = group.ID[1]
		d["ra"] = mean(group.RAJ2000)
		d["dec"] = mean(group.DEJ2000)

		xs = group.RV
		ws = @. 1/group.e_RV^2
		d["RV"] = lguys.mean(xs, ws)
		d["RV_err"] = 1/sqrt(sum(ws))
		d["RV_sigma"] = lguys.std(xs, ws)
		d["RV_count"] = length(xs)
		if length(xs) <= 1
			d["RV_sigma"] = NaN
		end
				
		append!(df, d, promote=true)

	end

	df
end

# ╔═╡ ae75e97c-b9c7-4ea9-8d91-7d5277d75a9e
spencer18_averaged = read_fits(joinpath(data_dir, "spencer+18_averaged.fits")) |> x-> rename!(x,
	:vlosmean => :RV,
	:e_vlosmean => :RV_err)

# ╔═╡ 05139b83-7d39-4011-b5e5-28f61186bef3
s18_ave_membs = filter(r->μ-3σ < r.RV < μ+3σ, spencer18_averaged)

# ╔═╡ d03045a6-6766-45e9-91cb-7830f93f99d7
model_averaged = RVUtils.model_vel_1c(s18_ave_membs.RV, s18_ave_membs.RV_err)

# ╔═╡ 6a72be02-4f86-4aee-ba4a-8e9413982f87
hist(s18_ave_membs.RV)

# ╔═╡ ff6046f5-2276-4dae-8a95-f000d6fbf3f9
hist(spencer_mixed.RV)

# ╔═╡ 2a43f4d1-ec91-40f5-a1f3-1750ffa06cfd
samples_averaged = sample(model_averaged, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ a46da768-07bb-4628-9f5b-abc708e30e91
pairplot(samples_averaged)

# ╔═╡ bc72b8b6-d302-459d-b07b-44dc909cb183
readdir(data_dir)

# ╔═╡ 1dbd4a17-f5eb-4e55-ac4c-4ed0065a1777
model_mixed = RVUtils.model_vel_1c(spencer_mixed.RV, spencer_mixed.RV_err .* sqrt.(spencer_mixed.RV_count))

# ╔═╡ b3bb14ee-3ced-436a-8067-f16052d0ea2d
samples_mixed = sample(model_mixed, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 35470ddd-1185-4d29-b1fa-ca3438bbfa3f
pairplot(samples_mixed)

# ╔═╡ 7538de40-970f-46ab-8e4c-ae1e6ed151fd
model_higherr = RVUtils.model_vel_1c(spencer18_members.RV[spencer18_members.RV_count .> 1], spencer18_members.RV_sigma[spencer18_members.RV_count .> 1])

# ╔═╡ 9437554c-2281-4602-a5a0-d5decc52af31
md"""
With multiple epoch observations, chi2 distribution looks quite good.
"""

# ╔═╡ 7a5f3817-9f49-43fe-a2cb-664feb83e81c
p_chi2_mixed = RVUtils.prob_chi2(spencer_mixed)

# ╔═╡ 48e544a1-b094-495f-be8e-655e2a32125a
hist(p_chi2_mixed)

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
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═37e6d04b-0347-4f14-9b06-d657a03df225
# ╟─a3f87e29-6095-436f-9095-deb2d16b994f
# ╠═aefa18bd-4203-4034-9ae9-9893fa8cc64d
# ╠═4e39ac0f-60d0-4091-a5e8-8aa3977e9f39
# ╠═0137b308-f0c6-4e73-afa6-346797b6f304
# ╠═153b1c3f-88a0-4b31-a2d7-65c084ca3a43
# ╠═2a6e8cda-a3cf-4824-9584-fa3d3c1da53a
# ╠═2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
# ╠═d321c481-0ceb-4e3b-b5fc-12af735155e3
# ╠═412148d2-3b14-400a-b93b-f4aa18965d8d
# ╠═96e1413d-e4fd-48af-8cde-73a5fcb4976a
# ╠═48ca64bd-fd0a-4e49-b281-37ac200ba3ff
# ╠═f751b13c-51fc-47c5-af2f-e1188c5bba5f
# ╠═a6f7f1a1-2dad-4397-9565-7b2aaf1fa733
# ╠═25e927e1-773d-4d41-894b-1f4512f01233
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
# ╠═8b993569-173e-4326-a778-cf6239771447
# ╠═e11da237-4c8f-450f-bd41-b2c9cdf2f11d
# ╠═9ac3364e-5532-457e-b29b-bfb10d32e845
# ╠═7962908c-4150-46a8-8388-2f5283b1d0c4
# ╠═34b1fcc3-7355-4d0d-b0aa-8af0541e8651
# ╠═08967635-801f-44dc-bb7e-62abbe0f86e4
# ╠═1853d75f-65e9-4ad0-8a5e-68603d4ecddf
# ╟─4290a51e-961b-4b38-b954-1a65e946080c
# ╟─8a176998-016d-4932-946f-eca78eb53d2b
# ╠═cee182e1-5545-4b86-9152-b288a9b11d86
# ╠═c0c02bce-e0f9-4fa0-bba6-a5360a825e23
# ╠═d3e682ee-0dab-4aa5-9f3e-423a9d2f92b0
# ╠═c873ccf2-8ecf-4a85-a152-76acd1c1fc62
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═fe505b23-c610-4e44-98f4-fb4c9ffe0388
# ╠═d3440c3e-be84-403f-a04d-68a237375b62
# ╠═33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╟─1a75bca0-c305-47a0-959e-afa42bd75d11
# ╠═761a7be6-420a-4428-addb-7947125596b3
# ╠═5774a0ef-8abf-4b6d-9b34-f80c29b4aab4
# ╟─05e35f22-e682-4ec6-acf4-1d7bd9a4ca79
# ╠═d86097ad-2feb-46d8-9ceb-a2743a74a1a9
# ╠═5de87579-1d6e-4dc0-860b-501816d75f8c
# ╟─214ed7f8-f338-4186-8f6d-870498486c0a
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═60f45b81-8beb-4eb0-8f55-db04c5368eee
# ╟─cb2923b2-dfda-4257-8d35-cd8f28188dbe
# ╟─28ff5e47-468d-4b54-bf24-d822818ec93d
# ╠═9847f6c4-05bd-481c-894d-4966d37be8c8
# ╠═fdfde937-39e5-4d51-8eb9-b2ec1ddd7dbd
# ╠═ae75e97c-b9c7-4ea9-8d91-7d5277d75a9e
# ╠═d03045a6-6766-45e9-91cb-7830f93f99d7
# ╠═05139b83-7d39-4011-b5e5-28f61186bef3
# ╠═6a72be02-4f86-4aee-ba4a-8e9413982f87
# ╠═ff6046f5-2276-4dae-8a95-f000d6fbf3f9
# ╠═2a43f4d1-ec91-40f5-a1f3-1750ffa06cfd
# ╠═a46da768-07bb-4628-9f5b-abc708e30e91
# ╠═bc72b8b6-d302-459d-b07b-44dc909cb183
# ╠═1dbd4a17-f5eb-4e55-ac4c-4ed0065a1777
# ╠═b3bb14ee-3ced-436a-8067-f16052d0ea2d
# ╠═35470ddd-1185-4d29-b1fa-ca3438bbfa3f
# ╠═7538de40-970f-46ab-8e4c-ae1e6ed151fd
# ╟─9437554c-2281-4602-a5a0-d5decc52af31
# ╠═7a5f3817-9f49-43fe-a2cb-664feb83e81c
# ╠═48e544a1-b094-495f-be8e-655e2a32125a
