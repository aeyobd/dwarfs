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
# TOLSTOY analy

This notebook analyzes the RV stars in tolstoy+23 and creates a sample sutable to be combined with others for velocity analysis.

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
No need to filter further...
"""

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79

begin 
	tolstoy23_all = read_fits("$data_dir/tolstoy+23.fits")

	rename!(tolstoy23_all, 
		:Vlos => :RV,
		:e_Vlos => :RV_err_orig,
		:s_Vlos => :RV_sigma,
		:Nspec => :RV_count,
		:RAJ2000 => :ra,
		:DEJ2000 => :dec,
		:GaiaDR3 => :source_id,
		:__Fe_H_ => :fe_h,
		:e__Fe_H_ => :fe_h_err,
	)

	tolstoy23_all[:, :RV_err] = @. max((1.448 - 0.021*tolstoy23_all.S_N)*tolstoy23_all.RV_err_orig, 0.65) # correction from arroyo-polonio 2024
	
	tolstoy23_all = tolstoy23_all
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_tolstoy+23.fits", tolstoy23_all, overwrite=true)

# ╔═╡ 33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
filt_tolstoy = RVUtils.sigma_clip(tolstoy23_all.RV, 3) .& (tolstoy23_all.Mem .== "m")

# ╔═╡ 344f3c29-0873-4180-9025-51aeaeb2c681
tolstoy23 = tolstoy23_all[filt_tolstoy, :]

# ╔═╡ 4290a51e-961b-4b38-b954-1a65e946080c
md"""
# Plots
"""

# ╔═╡ 76dc33cd-110a-42e6-bd72-65ba6b7bf85a
scatter(tolstoy23.S_N, tolstoy23.RV_err .- tolstoy23.RV_err_orig)

# ╔═╡ c8e08e37-60d8-4f14-adea-f6f93017cbb9
scatter(tolstoy23.S_N, tolstoy23.RV_err_orig)

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis()
	bins = make_bins(tolstoy23_all.RV, calc_limits(tolstoy23_all.RV), bandwidth=1)

	h = histogram(tolstoy23_all.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(tolstoy23.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ d86097ad-2feb-46d8-9ceb-a2743a74a1a9
hist((collect ∘ skipmissing)(tolstoy23_all.fe_h),
	 axis = (; xlabel = "[Fe/H]")
	)

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(tolstoy23.RV, tolstoy23.RV_err)

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = RVUtils.summarize(chain)

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(tolstoy23, samples, bins=30)

# ╔═╡ 60f45b81-8beb-4eb0-8f55-db04c5368eee
CSV.write("processed/rv_tolstoy+23.csv", df_summary)

# ╔═╡ 8b909559-c7e9-4fac-9d29-bf3b654393fa
md"""
# Dart
"""

# ╔═╡ de449b47-2040-4514-aac3-9fe4729d6213
md"""
DART sample given from Federico. 

The sample is used in Battaglia + 2008 (Kinematic status of Scl). I cannot find the whole dataset publically available. 

Measurements taken using VLT/FLAMES. 

The paper reports $v_{\rm hel, sys} = 110.6\pm0.5$ and a velocity dispersion of $\sigma=10.1 \pm 0.3$
"""

# ╔═╡ ea398623-ebc2-420a-8377-65bef43c0154
begin 
	dart = CSV.read("$data_dir/Sculptor_DARTS_BEST_psat40.csv", DataFrame)
	rename!(dart, 
		"vel"=>"RV",
		"evel"=>"RV_err",
		"feh" => "Fe_H",		
		# "feh_lo" => "Fe_H_lo",		
		# "feh_hi" => "Fe_H_hi",		
	)

	dart = dart[:, [:ra, :dec, :source_id, :RV, :RV_err, ]]
end

# ╔═╡ e23e7d8c-e1ec-48db-a96d-83798da42e10
let
	filt, idx = RVUtils.xmatch(dart, tolstoy23)
	@info "matched $(sum(filt)) out of $(length(filt))"

	df_matched = copy(dart[filt, :])
	df_matched[:, :RV_tolstoy] = tolstoy23[idx[filt], :RV]
	df_matched[:, :RV_err_tolstoy] = tolstoy23[idx[filt], :RV_err]
	
	rv1  = df_matched.RV
	rv1_err  = df_matched.RV_err
	
	rv2  = df_matched.RV_tolstoy
	rv2_err  = df_matched.RV_err_tolstoy

	filt = RVUtils.filt_missing(rv1; low=75, high=150) .& RVUtils.filt_missing(rv2;  low=75, high=150)


	fig, ax = FigAxis(
		xlabel = "RV DART",
		ylabel = "RV Tolstoy",
		aspect=DataAspect()
	)

	w = @. 1 / (rv1_err^2 + rv2_err^2)
	μ = LilGuys.mean(rv1 .- rv2, w)
	δμ = sqrt(1/sum(w))
	@info "mean = $μ ± $δμ"
	errorscatter!(rv1[filt], rv2[filt], xerror=rv1_err[filt], yerror=rv2_err[filt])

	lines!([80, 150], [80, 150], color=:black)
	
	fig
end

# ╔═╡ ca774514-a1ea-4af2-8be4-cfef0f009574
model_dart = RVUtils.model_vel_1c(dart.RV, dart.RV_err)

# ╔═╡ 603b18e1-55d3-4a85-b281-0ae2a20222ea
chain_dart = sample(model_dart, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 1cca058a-57b8-4be2-9c81-e67b651a38fc
pairplot(chain_dart)

# ╔═╡ 3524c534-7d9d-41e2-8802-0526878ad498
df_summary_dart = RVUtils.summarize(chain_dart)

# ╔═╡ 155ac7dc-7ec6-43a1-bea4-3c079ab0f19c
samples_dart = DataFrame(chain_dart)

# ╔═╡ e9e30f25-76aa-4bfc-8f30-edad90efd47c
RVUtils.plot_samples(dart, samples_dart, bins=30)

# ╔═╡ db4256f6-23df-49e5-8a3b-d81c03257b5a
CSV.write("processed/rv_dart.csv", df_summary_dart)

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
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
# ╟─4290a51e-961b-4b38-b954-1a65e946080c
# ╠═76dc33cd-110a-42e6-bd72-65ba6b7bf85a
# ╠═c8e08e37-60d8-4f14-adea-f6f93017cbb9
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═d86097ad-2feb-46d8-9ceb-a2743a74a1a9
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═60f45b81-8beb-4eb0-8f55-db04c5368eee
# ╟─8b909559-c7e9-4fac-9d29-bf3b654393fa
# ╟─de449b47-2040-4514-aac3-9fe4729d6213
# ╠═ea398623-ebc2-420a-8377-65bef43c0154
# ╠═e23e7d8c-e1ec-48db-a96d-83798da42e10
# ╠═ca774514-a1ea-4af2-8be4-cfef0f009574
# ╠═603b18e1-55d3-4a85-b281-0ae2a20222ea
# ╠═1cca058a-57b8-4be2-9c81-e67b651a38fc
# ╠═3524c534-7d9d-41e2-8802-0526878ad498
# ╠═155ac7dc-7ec6-43a1-bea4-3c079ab0f19c
# ╠═e9e30f25-76aa-4bfc-8f30-edad90efd47c
# ╠═db4256f6-23df-49e5-8a3b-d81c03257b5a
