### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe


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
No need to filter further...
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


# ╔═╡ aefa18bd-4203-4034-9ae9-9893fa8cc64d
filt_qual = abs.(spencer18_raw[:, "S-vlos"] .< 1) .&& (abs.(spencer18_raw[:, "K-vlos"] .- 3) .<= 1)

# ╔═╡ cee182e1-5545-4b86-9152-b288a9b11d86
hist(spencer18_raw.S_N)

# ╔═╡ 0137b308-f0c6-4e73-afa6-346797b6f304
j24 = read_fits("$data_dir/jensen+24_2c.fits") 

# ╔═╡ 153b1c3f-88a0-4b31-a2d7-65c084ca3a43
spencer18_all = let
	df = DataFrame()
	spencer18_raw[!, :ra] = round.(spencer18_raw.ra, digits=3)
	spencer18_raw[!, :dec] = round.(spencer18_raw.dec, digits=4)
	for group in groupby(spencer18_raw, [:ra, :dec])
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

# ╔═╡ 2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
xmatch_max = 2.0

# ╔═╡ d321c481-0ceb-4e3b-b5fc-12af735155e3
filt_xmatch, idx_xmatch = RVUtils.xmatch(spencer18_all, j24, xmatch_max)

# ╔═╡ 686ede3b-fb97-434a-abe1-337759e45ff9
sum(.!filt_xmatch)

# ╔═╡ 28a22929-89f2-422d-9cf0-b06d7e45d9a4
df_out = let
	df = copy(spencer18_all)
	source_id = allowmissing(j24.source_id[idx_xmatch])
	source_id[.!filt_xmatch] .= missing

	df[!, :source_id] .= source_id
	df
end

# ╔═╡ 423e446b-3cf6-4150-8b94-c4d7388fceb4
@assert length(unique(skipmissing(df_out.source_id))) == sum(filt_xmatch)

# ╔═╡ 011700b8-5dcc-409b-807c-952ea9726992
df_out[nonunique(df_out, :source_id, keep=:noduplicates), :]

# ╔═╡ 70729597-2e50-4773-89ad-87e97f7bdad6
length(unique(df_out.source_id)), sum(.!ismissing.(df_out.source_id)), sum(filt_xmatch)

# ╔═╡ 73132bc5-7908-4bf6-b2cf-cef8a3a8cf80
readdir("processed"), pwd()

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_spencer+18.fits", df_out, overwrite=true)

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

# ╔═╡ 4290a51e-961b-4b38-b954-1a65e946080c
md"""
# Plots
"""

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

# ╔═╡ 613d430d-b13e-492d-a21b-06350290f398
p_chi2 = 1 .- cdf.(Chisq.(max.(1., spencer18_members.RV_count .- 1)), spencer18_members.RV_sigma .^ 2 ./ spencer18_members.RV_err .^2)

# ╔═╡ d86097ad-2feb-46d8-9ceb-a2743a74a1a9
hist(p_chi2[spencer18_members.RV_count .> 1])

# ╔═╡ 5de87579-1d6e-4dc0-860b-501816d75f8c
hist(filter(isfinite, spencer18_members.RV_sigma ./ spencer18_members.RV_err))

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
CSV.write("processed/rv_spencer+18.csv", df_summary)

# ╔═╡ cb2923b2-dfda-4257-8d35-cd8f28188dbe
md"""
# other dataframe
"""

# ╔═╡ 9847f6c4-05bd-481c-894d-4966d37be8c8
spencer_mixed_raw = read_fits(joinpath(data_dir, "umi_rv_spencer+18.fits")) 

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

# ╔═╡ Cell order:
# ╠═811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═2f5ca80c-d98d-45cc-a661-834002b620f6
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
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
# ╠═aefa18bd-4203-4034-9ae9-9893fa8cc64d
# ╠═cee182e1-5545-4b86-9152-b288a9b11d86
# ╠═0137b308-f0c6-4e73-afa6-346797b6f304
# ╠═153b1c3f-88a0-4b31-a2d7-65c084ca3a43
# ╠═2c4eabb0-3af9-4b6f-b15d-9ed3b4290040
# ╠═d321c481-0ceb-4e3b-b5fc-12af735155e3
# ╠═686ede3b-fb97-434a-abe1-337759e45ff9
# ╠═28a22929-89f2-422d-9cf0-b06d7e45d9a4
# ╠═423e446b-3cf6-4150-8b94-c4d7388fceb4
# ╠═011700b8-5dcc-409b-807c-952ea9726992
# ╠═70729597-2e50-4773-89ad-87e97f7bdad6
# ╠═73132bc5-7908-4bf6-b2cf-cef8a3a8cf80
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╠═c0c02bce-e0f9-4fa0-bba6-a5360a825e23
# ╠═d3e682ee-0dab-4aa5-9f3e-423a9d2f92b0
# ╠═c873ccf2-8ecf-4a85-a152-76acd1c1fc62
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═fe505b23-c610-4e44-98f4-fb4c9ffe0388
# ╠═d3440c3e-be84-403f-a04d-68a237375b62
# ╠═33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
# ╟─4290a51e-961b-4b38-b954-1a65e946080c
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═d86097ad-2feb-46d8-9ceb-a2743a74a1a9
# ╠═613d430d-b13e-492d-a21b-06350290f398
# ╠═5de87579-1d6e-4dc0-860b-501816d75f8c
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═60f45b81-8beb-4eb0-8f55-db04c5368eee
# ╠═cb2923b2-dfda-4257-8d35-cd8f28188dbe
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
