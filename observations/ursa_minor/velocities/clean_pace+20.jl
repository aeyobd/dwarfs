### A Pluto.jl notebook ###
# v0.20.8

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
# Pace+22

This notebook analyzes the RV stars in tolstoy+23 and creates a sample sutable to be combined with others for velocity analysis.

"""

# ╔═╡ 89f5e91b-9970-4f8d-80ae-4766195e6a56
md"""
Depends on:
- `../data/jensen+24_2c.fits`
- `../data/pace+20.fits`
Creates:
- `processed/rv_pace+20.fits`
- `processed/mcmc_summary_uncorrected_pace+20.csv`

"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
not = !

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
pace20_raw = read_fits("$data_dir/pace+20.fits") |> x -> rename(x,
		:vLOS => :RV,
		:e_vLOS => :RV_err,
		:RAJ2000 => :ra,
		:DEJ2000 => :dec,
		:__Fe_H_ => :fe_h,
		:e__Fe_H_ => :fe_h_err,
	)

# ╔═╡ 0137b308-f0c6-4e73-afa6-346797b6f304
j24 = read_fits("$data_dir/jensen+24_2c.fits")

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
pace20_all = let
	df = DataFrame()

	for group in groupby(pace20_raw, :ID)
		d = OrderedDict()
		d["ID"] = group.ID[1]
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


		filt_fe = .!ismissing.(group.fe_h)
		xs = group.fe_h[filt_fe] |> disallowmissing
		ws = @.(1/group.fe_h_err[filt_fe]^2) |> disallowmissing
		d["fe_h"] = lguys.mean(xs, ws)
		d["fe_h_err"] = lguys.std(xs, ws)

		for col in ["GaiaDR2", "Megacam"]
			d[col] = group[1, col]
		end
		
		d["PdSph"] = mean(group.PdSph)
		d["PMR"] = mean(group.PMR)

				
		append!(df, d, promote=true)

	end

	df
end

# ╔═╡ d321c481-0ceb-4e3b-b5fc-12af735155e3
filt_xmatch, idx_xmatch = RVUtils.xmatch(pace20_all, j24, 1)

# ╔═╡ 686ede3b-fb97-434a-abe1-337759e45ff9
sum(.!filt_xmatch)

# ╔═╡ 032a02f9-1697-4f85-897c-a0f623a78da6
begin
	source_id = allowmissing(j24.source_id[idx_xmatch])
	source_id[.!filt_xmatch] .= missing
end

# ╔═╡ dce22a0d-a7dc-4642-abdc-68f83b714e94
pace20_all[.!filt_xmatch, :]

# ╔═╡ 3f6c0ba3-5f91-4c42-ba9a-833a8c3fa1d8
F_match = RVUtils.get_f_best.([j24], source_id) .=== 1.0

# ╔═╡ 29c38300-9a84-4f72-b6be-43fdc683d931
p_chi2 = RVUtils.prob_chi2(pace20_all)

# ╔═╡ cebb491c-88a2-47ae-b851-4a079b5066b7
F_scatter = RVUtils.filter_chi2(p_chi2)

# ╔═╡ 28a22929-89f2-422d-9cf0-b06d7e45d9a4
df_out = let
	df = copy(pace20_all)

	df[!, :F_match] = F_match
	df[!, :F_scatter] = F_scatter
	df[!, :source_id] = source_id
	df
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_pace+20.fits", df_out, overwrite=true)

# ╔═╡ 93d77566-540e-478f-88e1-76814f193175
md"""
# Numbers
"""

# ╔═╡ 02b17a2e-3bf9-4471-810c-78afa4a8403a
sum(pace20_all.RV_count), length(pace20_raw.RV)

# ╔═╡ 410b4a5c-0d64-4579-9a5e-312a6819be76
length(pace20_all.RV)

# ╔═╡ 5aef38da-cfea-4291-a021-34700ac874a7
sum(F_scatter)

# ╔═╡ 7ba28d8a-208e-473c-931b-a6ed33227db1
sum(.!ismissing.(F_match))

# ╔═╡ 11c198c8-1bcb-4112-b6ee-abc3b4cb4316
sum(skipmissing(F_match))

# ╔═╡ e6c02439-9043-452e-a4d5-306871b361ed
sum(skipmissing(F_scatter .& F_match))

# ╔═╡ c630d985-ee8a-4d75-8ebc-9e705712474a
matched = j24[idx_xmatch, :];

# ╔═╡ 8548b355-38f5-474c-b1e5-4b4fbfab4634
mean(matched[matched.F_BEST .== 0.0, :F_CPAR])

# ╔═╡ 55445dd0-6319-4307-8b8b-b76e9ae3b381
mean(matched[matched.F_BEST .== 0.0, :F_INGRID])

# ╔═╡ ee6a2453-33cb-47db-afcb-64db4e3ca0d8
mean(matched[matched.F_BEST .== 0.0, :F_ASTROMETRIC])

# ╔═╡ 4290a51e-961b-4b38-b954-1a65e946080c
md"""
# Plots
"""

# ╔═╡ 33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
filt_memb = RVUtils.sigma_clip(pace20_all.RV, 3) .& (pace20_all.PdSph .> 0.5) .& .!ismissing.(pace20_all.PdSph)

# ╔═╡ 344f3c29-0873-4180-9025-51aeaeb2c681
pace20_members = pace20_all[filt_memb, :]

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis()
	bins = make_bins(pace20_all.RV, calc_limits(pace20_all.RV), bandwidth=1)

	h = histogram(pace20_all.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(pace20_members.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ d86097ad-2feb-46d8-9ceb-a2743a74a1a9
hist(filter(isfinite, pace20_members.fe_h),
	 axis = (; xlabel = "[Fe/H]")
	)

# ╔═╡ 89d38852-1773-4c09-aa0a-3a5580816d57
md"""
## MCMC modeling
"""

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(pace20_members.RV, pace20_members.RV_err)

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = RVUtils.summarize(chain)

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(pace20_members, samples, bins=30)

# ╔═╡ 60f45b81-8beb-4eb0-8f55-db04c5368eee
CSV.write("processed/mcmc_summary_uncorrected_pace+20.csv", df_summary)

# ╔═╡ 32291b0f-1816-45a6-abcb-99c13512d8af
md"""
# Xmatch validation
"""

# ╔═╡ e856771f-100b-4229-aa51-13dc34860646
xmatch_max = 1

# ╔═╡ f8814fc8-8621-4e18-96f9-cf85b84ccfec
md"""
Uses normal xmatch and probability chi2 cuts, (max seperation = $xmatch_max arcsec)
"""

# ╔═╡ 1b45a429-59da-4156-ad75-d434a07dee0f
xmatch_dists = lguys.angular_distance.(df_out.ra, df_out.dec, j24.ra[idx_xmatch], j24.dec[idx_xmatch]) * 60^2

# ╔═╡ f671cf05-c097-4d22-8b73-49b75ccd26ac
sum(xmatch_dists .< xmatch_max), length(unique(df_out.source_id))

# ╔═╡ cb1b4342-2487-46fa-99ed-4cac8b74ee40
length(xmatch_dists)

# ╔═╡ a11d26b6-8416-4983-afc0-629433a45cbc
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

# ╔═╡ 88de547a-ec5d-446f-b092-0fa591d2709c
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

# ╔═╡ 94fbac70-e63b-4a24-8245-06ef2ac169c8
md"""
## Uncertainties
"""

# ╔═╡ 55954f5d-0f0b-4106-bc0a-57c1cd280cdb
hist(filter(isfinite, p_chi2))

# ╔═╡ a723770b-8a09-4e00-9766-910a5a91485a
md"""
## Failed xmatch
"""

# ╔═╡ 137cd49d-42cc-402d-a8bf-6883a9a071a4
let 
	fig, ax, p = hist(log10.(pace20_all.RV_err[filt_xmatch]))

	stephist!(log10.(pace20_all.RV_err[.!filt_xmatch]), color=COLORS[2])

	fig
end

# ╔═╡ 2c303dcd-ef0b-49fd-8b47-380675063070
scatter(1 ./ pace20_all.RV_err, filt_xmatch .+ rand(length(filt_xmatch)), alpha=0.1)

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═89f5e91b-9970-4f8d-80ae-4766195e6a56
# ╟─f8814fc8-8621-4e18-96f9-cf85b84ccfec
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
# ╠═0137b308-f0c6-4e73-afa6-346797b6f304
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═d321c481-0ceb-4e3b-b5fc-12af735155e3
# ╠═686ede3b-fb97-434a-abe1-337759e45ff9
# ╠═032a02f9-1697-4f85-897c-a0f623a78da6
# ╠═dce22a0d-a7dc-4642-abdc-68f83b714e94
# ╠═3f6c0ba3-5f91-4c42-ba9a-833a8c3fa1d8
# ╠═29c38300-9a84-4f72-b6be-43fdc683d931
# ╠═cebb491c-88a2-47ae-b851-4a079b5066b7
# ╠═28a22929-89f2-422d-9cf0-b06d7e45d9a4
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╟─93d77566-540e-478f-88e1-76814f193175
# ╠═02b17a2e-3bf9-4471-810c-78afa4a8403a
# ╠═410b4a5c-0d64-4579-9a5e-312a6819be76
# ╠═5aef38da-cfea-4291-a021-34700ac874a7
# ╠═7ba28d8a-208e-473c-931b-a6ed33227db1
# ╠═11c198c8-1bcb-4112-b6ee-abc3b4cb4316
# ╠═e6c02439-9043-452e-a4d5-306871b361ed
# ╠═c630d985-ee8a-4d75-8ebc-9e705712474a
# ╠═8548b355-38f5-474c-b1e5-4b4fbfab4634
# ╠═55445dd0-6319-4307-8b8b-b76e9ae3b381
# ╠═ee6a2453-33cb-47db-afcb-64db4e3ca0d8
# ╟─4290a51e-961b-4b38-b954-1a65e946080c
# ╠═33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═d86097ad-2feb-46d8-9ceb-a2743a74a1a9
# ╟─89d38852-1773-4c09-aa0a-3a5580816d57
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═60f45b81-8beb-4eb0-8f55-db04c5368eee
# ╟─32291b0f-1816-45a6-abcb-99c13512d8af
# ╠═e856771f-100b-4229-aa51-13dc34860646
# ╠═1b45a429-59da-4156-ad75-d434a07dee0f
# ╠═f671cf05-c097-4d22-8b73-49b75ccd26ac
# ╠═cb1b4342-2487-46fa-99ed-4cac8b74ee40
# ╠═a11d26b6-8416-4983-afc0-629433a45cbc
# ╠═88de547a-ec5d-446f-b092-0fa591d2709c
# ╟─94fbac70-e63b-4a24-8245-06ef2ac169c8
# ╠═55954f5d-0f0b-4106-bc0a-57c1cd280cdb
# ╟─a723770b-8a09-4e00-9766-910a5a91485a
# ╠═137cd49d-42cc-402d-a8bf-6883a9a071a4
# ╠═2c303dcd-ef0b-49fd-8b47-380675063070
