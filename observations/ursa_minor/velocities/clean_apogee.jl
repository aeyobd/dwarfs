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

# ╔═╡ 05c93bf4-0f94-44cc-9aab-7c6193831806
using PyFITS

# ╔═╡ 2f5ca80c-d98d-45cc-a661-834002b620f6
using PairPlots

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
using Turing

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Apogee RV analsysi.

This notebook analyzes the RV stars in APOGEE and creates a sample sutable to be combined with others for velocity analysis.

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

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
module RVUtils
	include("../../rv_utils.jl")
	not = !
end

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../data/"

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_fits("../data/jensen+24_2c.fits")

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Apogee
"""

# ╔═╡ 45d6b4aa-ca44-4a71-afb8-ba6b2e674c7a
md"""
APOGEE DR 17 sample from federico's paper is apogee_raw_2. I have also xmatched the field against apogee which I use for self-consistency.
"""

# ╔═╡ d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
apogee_raw = read_fits("processed/apogee_xmatch.fits")

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
begin 
	apogee_all = copy(apogee_raw)
	rename!(apogee_all, 
		"FE_H" => "fe_h",
		"FE_H_ERR" => "fe_h_err"
	)

	_apogee_all_filt = not.(ismissing.(apogee_all.RV))

	apogee_all = DataFrame(apogee_all[_apogee_all_filt, :])

	apogee_all[!, :RV] = float.(apogee_all.RV)

	apogee_all = apogee_all[:, [:source_id, :ra, :dec, :RV, :RV_err, :RV_sigma, :RV_count, :fe_h, :fe_h_err, :STARFLAG, :RV_FLAG, :SNR]]
	apogee_all[ismissing.(apogee_all.fe_h), :fe_h] .= NaN
	apogee_all[ismissing.(apogee_all.fe_h_err), :fe_h_err] .= NaN

	apogee_all
end

# ╔═╡ 961e230d-e549-49b5-a98c-aaa330fa42da
@assert 0 == sum(apogee_all.RV_FLAG .> 0) # no need to include more flags

# ╔═╡ bb57569e-0a4c-455e-b3c6-7bd2897bb843
sum(apogee_raw.RV_FLAG .> 0)

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_apogee.fits", apogee_all, overwrite=true)

# ╔═╡ 59a20996-7b90-4afd-818a-fffd8029350f
md"""
# Plots and analysis
"""

# ╔═╡ 4d63430e-f59c-4c68-97e1-7eaa5679e55f
apogee = apogee_all[-300 .< apogee_all.RV .< -230, :]

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis(
		xlabel = "RV",
		ylabel = "count"
	)
	bins = make_bins(apogee_all.RV, calc_limits(apogee_all.RV), bandwidth=1)

	h = histogram(apogee.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(apogee_all.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
filt_apogee = RVUtils.sigma_clip(apogee_all.RV, 3)

# ╔═╡ 1bc2541c-69ac-40aa-94f5-743707ca5a66
hist(apogee.RV)

# ╔═╡ 63792b12-dde2-4361-88b2-1f5a789631ff
hist(filter(isfinite, apogee.fe_h))

# ╔═╡ 5dafeaad-df9e-4558-99a7-721bdb91d28d
hist(apogee.RV_err)

# ╔═╡ 949bdb6a-4554-4659-b0d3-7e669f470c10
hist(apogee.RV_sigma)

# ╔═╡ fba487be-22a5-43b3-9ac5-8505462b7d56
hist(apogee.RV_sigma ./ (apogee.RV_err .* sqrt.(apogee.RV_count) ))

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(apogee.RV, apogee.RV_err)

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = RVUtils.summarize(chain)

# ╔═╡ 5db61bab-a71e-4bcb-a706-f150f1126f4c
CSV.write("processed/rv_apogee.csv", df_summary)

# ╔═╡ ba337b2c-9a0c-40b3-867f-c9ebf683e1e3
apogee.RV

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(apogee, samples, bins=30)

# ╔═╡ 0ba4323f-8ffe-4ac7-999e-981f8c282708


# ╔═╡ Cell order:
# ╠═811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═05c93bf4-0f94-44cc-9aab-7c6193831806
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═2f5ca80c-d98d-45cc-a661-834002b620f6
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═36634dea-21bc-4823-8a15-7bce20b6fc17
# ╠═e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╟─45d6b4aa-ca44-4a71-afb8-ba6b2e674c7a
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═961e230d-e549-49b5-a98c-aaa330fa42da
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╠═bb57569e-0a4c-455e-b3c6-7bd2897bb843
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╟─59a20996-7b90-4afd-818a-fffd8029350f
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═4d63430e-f59c-4c68-97e1-7eaa5679e55f
# ╠═bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
# ╠═1bc2541c-69ac-40aa-94f5-743707ca5a66
# ╠═63792b12-dde2-4361-88b2-1f5a789631ff
# ╠═5dafeaad-df9e-4558-99a7-721bdb91d28d
# ╠═949bdb6a-4554-4659-b0d3-7e669f470c10
# ╠═fba487be-22a5-43b3-9ac5-8505462b7d56
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═5db61bab-a71e-4bcb-a706-f150f1126f4c
# ╠═ba337b2c-9a0c-40b3-867f-c9ebf683e1e3
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═0ba4323f-8ffe-4ac7-999e-981f8c282708
