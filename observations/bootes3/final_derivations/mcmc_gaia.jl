### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using PyFITS
	using DataFrames, CSV
	import TOML
	
	using Turing
	import PairPlots
end

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	Nsamples = 3_000
	Nthreads = 16
	write_results = true
else
	Nsamples = 300
	Nthreads = 1
	write_results = false
end

# ╔═╡ 066c7b30-f818-4e4a-8db8-c8bac469f558
module MCMCUtils
	include("../mcmc_utils.jl")
end

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 508c40e3-030f-425d-ac79-e2f27e97ee37
md"""
## Inputs
"""

# ╔═╡ 9f043c2e-5357-41fb-a9ec-d5464e1a11dd
samplename = "j24_1c"

# ╔═╡ 6413e7ca-b614-4c59-9853-274444232c41
n_R_h_gc = 10

# ╔═╡ e10ec7e0-e71f-4290-ab52-189c8763f887
p_min = 0.0

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
outdir = "mcmc"

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	using LilGuys

	FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".$samplename.mcmc"
end

# ╔═╡ df3ab039-55ee-4431-bd68-1c77f843dd18
sampler = NUTS(0.65)

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ 7e8124ea-7bbe-465b-a9dc-4b14d268c39e
obs_props = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 398e5019-006b-4d4d-922c-7ce244470e8e
ra_0, dec_0 = obs_props["ra_original"], obs_props["dec_original"]

# ╔═╡ 01910592-9514-46aa-85ac-f46fa23262fc
allstars = read_fits("../data/$samplename.fits")

# ╔═╡ 4acaad09-7d87-444e-b5ec-7fcf79f895ef
best_stars = let
	df = allstars[allstars.F_BEST .== 1.0, :]
	df[!, :xi], df[!, :eta] = 60 .* LilGuys.to_tangent(df.ra, df.dec, ra_0, dec_0)
	df[!, :R_ell] = @. sqrt(df.xi^2 + df.eta^2)

	df[:, :L_sat] = df.L_CMD_SAT .* df.L_PM_SAT
	df[:, :L_bg] = df.L_CMD_BKD .* df.L_PM_BKD

	df
end

# ╔═╡ 5936850b-90b6-48a5-9a36-8ab4c4f05724
md"""
### removing GCs
"""

# ╔═╡ d17c50db-ce94-4aad-92e2-6355b3b3a91b
props_ngc5466 = TOML.parsefile("../observed_properties_ngc5466.toml")

# ╔═╡ 56724dcf-b46d-406a-83a9-44b773800e88
props_ngc5272 = TOML.parsefile("../observed_properties_ngc5272.toml")

# ╔═╡ 82a0c53a-989f-4e0d-81ad-1c49d4f7050e
R_ngc5466 = n_R_h_gc*props_ngc5466["R_h"]

# ╔═╡ f3fecad0-f214-4772-8da5-841772670592
R_ngc5272 = n_R_h_gc*props_ngc5272["R_h"]

# ╔═╡ 8eaa620f-a8db-4b5e-89c6-2f54105b5933
function radius_to_gc(stars, props_gc)
	xi, eta = LilGuys.to_tangent(stars.ra, stars.dec, props_gc["ra"], props_gc["dec"])
	R = @. sqrt(xi^2 + eta^2) .* 60
	return R
end

# ╔═╡ 61793e5e-6435-492a-85d4-383dc9cf4f01
function excise_gcs(stars)


	filt = radius_to_gc(stars, props_ngc5272) .> R_ngc5272

	filt .&= radius_to_gc(stars, props_ngc5466) .> R_ngc5466


	return filt
end

# ╔═╡ 2f9eca06-d0c9-405b-8112-dd6132df1273
md"""
### The final sample
"""

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = maximum(best_stars.R_ell)

# ╔═╡ a736178a-ae01-4940-afa3-e7dc3522ec47
stars = let
	df = copy(best_stars)
 
	filt = df.R_ell .< R_max
	filt .&= df.L_sat ./ df.L_bg .> p_min
	filt .&= excise_gcs(df)

	
	df[filt, :]
end

# ╔═╡ 36784213-6393-4e33-8fb9-c4ef5b1fbc48
area_tot = π * (R_max^2 - R_ngc5466^2 - R_ngc5272^2)  

# ╔═╡ e33170f3-8c57-4672-b8a7-f083d69437c1
md"""
# Removing gc's
"""

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Setup
"""

# ╔═╡ 626846b5-156a-4883-8542-463e3a90faaa
function write_samples_summary(samples, df_samples, model_suffix)
	if !write_results
		return
	end
	samplesout = joinpath(outdir, "samples.$samplename.mcmc_$model_suffix.csv")
	summaryout = joinpath(outdir, "summary.$samplename.mcmc_$model_suffix.csv")
	CSV.write(samplesout, df_samples)
	CSV.write(summaryout, MCMCUtils.summarize(samples))
end

# ╔═╡ fab117b3-53a7-48c8-a6a1-8f004d58a90d
N_stars = size(stars, 1)

# ╔═╡ 14ac9c2a-b218-43c5-a129-3eda0774608a


# ╔═╡ 8fa15363-191e-4875-95f1-3ce87762e582
function to_frame(samples)
	df = DataFrame(samples)
	df[!, :N_sat] = N_stars * df.f_sat
	df
end

# ╔═╡ e82534ab-257c-4d18-b08e-3c37a05ce623
md"""
# Models
"""

# ╔═╡ 641dbab7-224e-473f-8772-9539a00bdd69
mcmc_data = MCMCUtils.GaiaData(
	xi = stars.xi, 
	eta = stars.eta,
	L_bg = stars.L_bg,
	L_sat = stars.L_sat,
	source_id = stars.source_id,
)

# ╔═╡ 215cde3e-88a3-4273-8618-32919a53531f
md"""
## Elliptical Exp model
"""

# ╔═╡ 74c1c7dc-0d60-4da4-9cdb-a0c4d8e571f7
mcmc_model_exp = MCMCUtils.exp_ell_model(mcmc_data, area_tot=area_tot)

# ╔═╡ dd1aaa3d-e964-4655-9a16-41e23e78a21d
samples_exp = sample(mcmc_model_exp, sampler, MCMCThreads(), Nsamples, Nthreads)

# ╔═╡ 1f8fa1da-0328-40f6-80a5-dd2d17dbd356
df_exp = to_frame(samples_exp)

# ╔═╡ 8657ffb8-ccb1-4dec-95c2-3c6273639af5
@savefig "corner.exp_ell" PairPlots.pairplot(samples_exp)

# ╔═╡ 6d6bbafb-65f7-4beb-b737-ef56afc826d1
write_samples_summary(samples_exp, df_exp, "exp_ell")

# ╔═╡ 23c01f4f-7cce-4b15-bcc6-fa0601addf84
md"""
### Elliptical Plummer model
"""

# ╔═╡ 32e5278a-8674-45fc-af41-a93f009fc008
mcmc_model_ell = MCMCUtils.plummer_ell_model(mcmc_data, area_tot=area_tot)

# ╔═╡ 725a1902-2c48-4b14-a3c2-a87afcefa080
samples_ell = sample(mcmc_model_ell, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ 7819e763-cebe-476b-86ec-641f035bf910
df_ell = to_frame(samples_ell)

# ╔═╡ 09f26d6f-91de-45a0-b3c8-d6f8e0104b53
@savefig "corner.plummer_ell" PairPlots.pairplot(samples_ell)

# ╔═╡ 66a6a188-c5f4-4761-96ba-672320727dc9
write_samples_summary(samples_ell, df_ell, "ell")

# ╔═╡ a3d02bd8-f330-4f13-8a1c-59499537748c
md"""
## Elliptical Sérsic
"""

# ╔═╡ 1f789712-401a-4942-b499-073d2f25408c
mcmc_model_sersic_ell = MCMCUtils.sersic_ell_model(mcmc_data, area_tot=area_tot, R_max=R_max)

# ╔═╡ 43ccc903-7857-4207-8652-255339f5ddc2
samples_sersic_ell = sample(mcmc_model_sersic_ell, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ 0b7e6509-fd02-47cf-bd45-aa2cafdc18f1
df_sersic_ell = to_frame(samples_sersic_ell)

# ╔═╡ 496ede90-f72f-4d64-b27f-714dc27b0fcb
summary_sersic_ell = MCMCUtils.summarize(samples_sersic_ell)

# ╔═╡ 451850ad-2e06-4609-9a58-626984146b9e
@savefig "corner.sersic_ell" PairPlots.pairplot(samples_sersic_ell)

# ╔═╡ 4979f59b-5297-4596-98d7-b06e84ed3472
write_samples_summary(samples_sersic_ell, df_sersic_ell, "sersic_ell")

# ╔═╡ 94075c3a-ad48-4648-9d76-ec0fd3a875bb
scatter(stars.xi, stars.eta, markersize=1, alpha=0.1)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═066c7b30-f818-4e4a-8db8-c8bac469f558
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─508c40e3-030f-425d-ac79-e2f27e97ee37
# ╠═9f043c2e-5357-41fb-a9ec-d5464e1a11dd
# ╠═6413e7ca-b614-4c59-9853-274444232c41
# ╠═e10ec7e0-e71f-4290-ab52-189c8763f887
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═398e5019-006b-4d4d-922c-7ce244470e8e
# ╠═01910592-9514-46aa-85ac-f46fa23262fc
# ╠═4acaad09-7d87-444e-b5ec-7fcf79f895ef
# ╠═5936850b-90b6-48a5-9a36-8ab4c4f05724
# ╠═d17c50db-ce94-4aad-92e2-6355b3b3a91b
# ╠═56724dcf-b46d-406a-83a9-44b773800e88
# ╠═82a0c53a-989f-4e0d-81ad-1c49d4f7050e
# ╠═f3fecad0-f214-4772-8da5-841772670592
# ╠═8eaa620f-a8db-4b5e-89c6-2f54105b5933
# ╠═61793e5e-6435-492a-85d4-383dc9cf4f01
# ╟─2f9eca06-d0c9-405b-8112-dd6132df1273
# ╠═a736178a-ae01-4940-afa3-e7dc3522ec47
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═36784213-6393-4e33-8fb9-c4ef5b1fbc48
# ╟─e33170f3-8c57-4672-b8a7-f083d69437c1
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═626846b5-156a-4883-8542-463e3a90faaa
# ╠═fab117b3-53a7-48c8-a6a1-8f004d58a90d
# ╠═14ac9c2a-b218-43c5-a129-3eda0774608a
# ╠═8fa15363-191e-4875-95f1-3ce87762e582
# ╠═e82534ab-257c-4d18-b08e-3c37a05ce623
# ╠═641dbab7-224e-473f-8772-9539a00bdd69
# ╠═215cde3e-88a3-4273-8618-32919a53531f
# ╠═74c1c7dc-0d60-4da4-9cdb-a0c4d8e571f7
# ╠═dd1aaa3d-e964-4655-9a16-41e23e78a21d
# ╠═1f8fa1da-0328-40f6-80a5-dd2d17dbd356
# ╠═8657ffb8-ccb1-4dec-95c2-3c6273639af5
# ╠═6d6bbafb-65f7-4beb-b737-ef56afc826d1
# ╟─23c01f4f-7cce-4b15-bcc6-fa0601addf84
# ╠═32e5278a-8674-45fc-af41-a93f009fc008
# ╠═725a1902-2c48-4b14-a3c2-a87afcefa080
# ╠═7819e763-cebe-476b-86ec-641f035bf910
# ╠═09f26d6f-91de-45a0-b3c8-d6f8e0104b53
# ╠═66a6a188-c5f4-4761-96ba-672320727dc9
# ╠═a3d02bd8-f330-4f13-8a1c-59499537748c
# ╠═1f789712-401a-4942-b499-073d2f25408c
# ╠═43ccc903-7857-4207-8652-255339f5ddc2
# ╠═0b7e6509-fd02-47cf-bd45-aa2cafdc18f1
# ╠═496ede90-f72f-4d64-b27f-714dc27b0fcb
# ╠═451850ad-2e06-4609-9a58-626984146b9e
# ╠═4979f59b-5297-4596-98d7-b06e84ed3472
# ╠═94075c3a-ad48-4648-9d76-ec0fd3a875bb
