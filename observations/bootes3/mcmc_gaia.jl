### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using PyFITS
	import PairPlots
	
	using Turing
end

# ╔═╡ c1816838-9ab4-45bc-b9ae-98157c70a40c
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 2ef405a3-9bdf-4255-87d8-901487f2d27d
include(joinpath(ENV["DWARFS_ROOT"], "utils/pluto_utils.jl"))

# ╔═╡ 0a974ff8-7c22-49dd-98f4-738537175d20
@bind inputs confirm(notebook_inputs(;
		samplename = TextField(40, default="best_sample"),
		mock = CheckBox(),
		simple = CheckBox(),
		p_min = NumberField(default=0),
		n_R_h_gc = NumberField(default=10),
	))


# ╔═╡ dce2754b-b345-404f-a564-c26438cc835f
if !@isdefined(PlutoRunner)
	run_all = true
else
	md"run all models $(@bind run_all CheckBox())"
end

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	samplename = ARGS[1]
	Nsamples = 1_000
	Nthreads = 16
else

	Nsamples = 300
	Nthreads = 1

	samplename = inputs.samplename
end

# ╔═╡ 90c12d18-f1c5-4813-8289-8c030be49ab3
mock = inputs.mock

# ╔═╡ 17689578-56d6-44a2-955c-be934d41c7b8
simple = inputs.simple

# ╔═╡ e10ec7e0-e71f-4290-ab52-189c8763f887
p_min = inputs.p_min

# ╔═╡ a0274b6a-daba-48aa-a289-86063dba9393
n_R_h_gc = inputs.n_R_h_gc

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
outdir = "mcmc"

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	using LilGuys

	FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".$samplename.mcmc"
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 066c7b30-f818-4e4a-8db8-c8bac469f558
module MCMCUtils
	include("./mcmc_utils.jl")
end

# ╔═╡ 5426e4db-6810-443a-a98b-6e1917fee397
module GaiaFilters
	include("../../utils/gaia_filters.jl")
end

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 626846b5-156a-4883-8542-463e3a90faaa
function write_samples_summary(samples, df_samples, model_suffix)
	samplesout = joinpath(outdir, "samples.$samplename.mcmc_$model_suffix.csv")
	summaryout = joinpath(outdir, "summary.$samplename.mcmc_$model_suffix.csv")
	CSV.write(samplesout, df_samples)
	CSV.write(summaryout, MCMCUtils.summarize(samples))
end

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ 7e8124ea-7bbe-465b-a9dc-4b14d268c39e
obs_props = let
	df = TOML.parsefile("observed_properties.toml")
	if mock
		df["ra"] = 205 - 0.3
		df["dec"] = 22 + 0.2
	else
		
		df["ra"] = df["ra_original"]
		df["dec"] = df["dec_original"]
	end
	df["ellipticity"] = 0.0
	df
end

# ╔═╡ 398e5019-006b-4d4d-922c-7ce244470e8e
ra_0, dec_0 = obs_props["ra"], obs_props["dec"]

# ╔═╡ 4acaad09-7d87-444e-b5ec-7fcf79f895ef
allstars = let
	df = read_fits("samples/$samplename.fits")
	df[!, :xi], df[!, :eta] = 60 .* LilGuys.to_tangent(df.ra, df.dec, ra_0, dec_0)

	if "L_CMD_SAT" ∈ names(df)
		df[:, :L_sat] = df.L_CMD_SAT .* df.L_PM_SAT
		df[:, :L_bg] = df.L_CMD_BKD .* df.L_PM_BKD
	end
	df
end

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = GaiaFilters.calc_R_max(allstars.xi, allstars.eta, 0, 0)

# ╔═╡ 1b77bc6a-0512-4d5e-a817-ce0eab032129
@assert 60 < R_max < 12*60

# ╔═╡ e33170f3-8c57-4672-b8a7-f083d69437c1
md"""
# Removing gc's
"""

# ╔═╡ ea690a36-512e-4a80-a258-faece76c299c
# props_ngc5272 = TOML.parsefile("observed_properties_ngc5272.toml")

# ╔═╡ efea8cbe-dd33-419a-9a49-e07a3424dbe3
# R_ngc5272 = 0*n_R_h_gc*props_ngc5272["R_h"]

# ╔═╡ d17c50db-ce94-4aad-92e2-6355b3b3a91b
props_ngc5466 = TOML.parsefile("observed_properties_ngc5466.toml")

# ╔═╡ 82a0c53a-989f-4e0d-81ad-1c49d4f7050e
R_ngc5466 = n_R_h_gc*props_ngc5466["R_h"]

# ╔═╡ 36784213-6393-4e33-8fb9-c4ef5b1fbc48
area_tot = π * (R_max^2 - R_ngc5466^2)  

# ╔═╡ 61793e5e-6435-492a-85d4-383dc9cf4f01
function excise_gcs(stars)
	# xi, eta = LilGuys.to_tangent(stars.ra, stars.dec, props_ngc5272["ra"], props_ngc5272["dec"])
	# dist_ngc5272 = @. sqrt(xi^2 + eta^2) .* 60

	# filt = dist_ngc5272 .> R_ngc5272


	xi, eta = LilGuys.to_tangent(stars.ra, stars.dec, props_ngc5466["ra"], props_ngc5466["dec"])
	dist_ngc5466 = @. sqrt(xi^2 + eta^2) .* 60

	filt = dist_ngc5466 .> R_ngc5466

	return stars[filt, :]
end

# ╔═╡ a736178a-ae01-4940-afa3-e7dc3522ec47
stars = let
	df = copy(allstars)

 	
	df[!, :R_ell] = @. sqrt(df.xi^2 + df.eta^2)

	filt = (df.R_ell .< R_max) .& (df.L_sat ./ df.L_bg .> p_min)
	df = df[filt, :]

	
	df = excise_gcs(df)
	if simple
		df.L_sat .= 1
		df.L_bg .= 1
	end
	df
end

# ╔═╡ 4aac4483-6f0e-4b73-8c64-cfd3e9995b59
scatter(stars.xi, stars.eta, markersize=0.5, alpha=0.3, color=:black)

# ╔═╡ 0215620c-6193-4215-b018-6f775e31bb5c
let
	LLR_min = -2

	fig = Figure()
	ax = Axis(fig[1,1])
	LLR = log10.(stars.L_sat ./ stars.L_bg)
	p = scatter!(stars.xi, stars.eta, color=LLR, alpha=0.3, markersize=max.(LLR_min, LLR) ./ 3, colorrange=(LLR_min,2), colormap=Reverse(Arya.get_arya_cmap()))

	Colorbar(fig[1,2], p)

	fig
end

# ╔═╡ 23a250b8-a7aa-4d0b-b838-a5e4436b0bc1
hist(log10.(stars.L_sat ./ stars.L_bg), bins=LinRange(-30, 3, 100))

# ╔═╡ fb28f68a-4893-4feb-aa29-5cc6acd7d2a8
let
	filt = stars.L_sat ./ stars.L_bg .> 9

	fig = Figure()
	ax = Axis(fig[1,1])
	p = scatter!(stars.xi[filt], stars.eta[filt], color=( stars.phot_g_mean_mag[filt]), alpha=0.3)

	Colorbar(fig[1,2], p)

	fig
end

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Setup
"""

# ╔═╡ fab117b3-53a7-48c8-a6a1-8f004d58a90d
N_stars = size(stars, 1)

# ╔═╡ 14ac9c2a-b218-43c5-a129-3eda0774608a
function plot_corner(samples)
	df = DataFrame(samples)
	df[!, :log_R_h] = log10.(df.R_h)

	cols = [v for v in values(samples.info.varname_to_symbol)]
	cols[cols .== :R_h] .= :log_R_h

	PairPlots.pairplot(df[:, cols])
end

# ╔═╡ 8fa15363-191e-4875-95f1-3ce87762e582
function to_frame(samples)
	df = DataFrame(samples)
	df[!, :N_sat] = N_stars * df.f_sat
	df
end

# ╔═╡ e82534ab-257c-4d18-b08e-3c37a05ce623
md"""
## Run the model
"""

# ╔═╡ df3ab039-55ee-4431-bd68-1c77f843dd18
sampler = NUTS()

# ╔═╡ cd12f215-cf8c-4d15-93f8-ebb0d5762c18
md"""
### Spherical Plummer model
"""

# ╔═╡ b962db4f-f7b9-40a7-810e-d14c2f6bc39c
mcmc_data = MCMCUtils.GaiaData(
	source_id = 1:N_stars,
	xi = stars.xi, 
	eta = stars.eta,
	L_bg = stars.L_bg,
	L_sat = stars.L_sat
)

# ╔═╡ fe3d7c64-c3ff-431f-976a-9dab7e8e71b3
mcmc_model_plummer = MCMCUtils.plummer_model(mcmc_data, area_tot=area_tot)

# ╔═╡ 051f456b-5b6e-438f-8438-392456759ae6
samples_plummer = sample(mcmc_model_plummer, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ f721cfdc-e8cb-4fc6-870d-d50cfb671fc6
df_plummer = to_frame(samples_plummer)

# ╔═╡ 3f592c6f-f35a-4977-9135-2987060dfadb
@savefig "corner.plummer" PairPlots.pairplot(samples_plummer)

# ╔═╡ 0430ac07-77a7-40a8-86dd-bac3b3a4fe24
LilGuys.from_tangent(median(df_plummer.d_xi)/60, median(df_plummer.d_eta)/60, ra_0, dec_0)

# ╔═╡ 61e9d125-a0d6-4956-bc78-9ca507598e94
@info "number stars $(MCMCUtils.to_measurement(df_plummer.N_sat))"

# ╔═╡ 1eb86ece-0fc7-453a-aa4e-f8e50cc81229
mock && sum(.!ismissing.(stars.Mini)) # true number of stars

# ╔═╡ 4b1d12dc-ed81-46b2-bfc3-53feab3829d4
write_samples_summary(samples_plummer, df_plummer, "plummer")

# ╔═╡ 23c01f4f-7cce-4b15-bcc6-fa0601addf84
md"""
### Elliptical Plummer model
"""

# ╔═╡ 32e5278a-8674-45fc-af41-a93f009fc008
mcmc_model_ell = MCMCUtils.plummer_ell_model(mcmc_data, area_tot=area_tot)

# ╔═╡ 725a1902-2c48-4b14-a3c2-a87afcefa080
if run_all
	samples_ell = sample(mcmc_model_ell, sampler, MCMCThreads(), Nsamples, Nthreads) 
end

# ╔═╡ 7819e763-cebe-476b-86ec-641f035bf910
df_ell = to_frame(samples_ell)

# ╔═╡ 09f26d6f-91de-45a0-b3c8-d6f8e0104b53
@savefig "corner.plummer_ell" PairPlots.pairplot(samples_ell)

# ╔═╡ 0ddaa045-1b43-405e-af05-066e04bd3918
@info "number stars ell: $(MCMCUtils.to_measurement(df_ell.N_sat))"

# ╔═╡ 66a6a188-c5f4-4761-96ba-672320727dc9
write_samples_summary(samples_ell, df_ell, "ell")

# ╔═╡ a02fbdd1-b38d-4540-999c-2096bcc50c10
mock && sum(.!ismissing.(stars.Mini)) # true number of stars

# ╔═╡ c708b169-c593-4c38-a91f-40723249993f
md"""
### Sersic model
"""

# ╔═╡ 7b1b4c0f-aa49-4ee0-b860-0bd927db8768
mcmc_model_sersic = MCMCUtils.sersic_model(mcmc_data, area_tot=area_tot, R_max=R_max,
	prior_R_h = Uniform(0, 200))

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
if run_all
	samples_sersic = sample(mcmc_model_sersic, sampler, MCMCThreads(), Nsamples, Nthreads) 
end

# ╔═╡ aca9c6b7-b7be-4f79-ba47-438048285041
df_sersic = to_frame(samples_sersic)

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
summary_sersic = MCMCUtils.summarize(samples_sersic)

# ╔═╡ 428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
@savefig "corner_s" PairPlots.pairplot(samples_sersic)

# ╔═╡ 7ed82514-9484-47d9-a332-7501c5e5d84d
@info "number stars sersic: $(MCMCUtils.to_measurement(df_sersic.N_sat))"

# ╔═╡ 56cc5edd-d56f-434d-9f6e-47e6e00092ed
write_samples_summary(samples_sersic, df_sersic, "sersic")

# ╔═╡ d57d7605-3904-490b-b785-42320275b0c5
@info "acceptance rate (all, per-step)"  mean(df_sersic.acceptance_rate), mean(df_sersic.is_accept)

# ╔═╡ a3d02bd8-f330-4f13-8a1c-59499537748c
md"""
# Elliptical Sérsic
"""

# ╔═╡ 1f789712-401a-4942-b499-073d2f25408c
mcmc_model_sersic_ell = MCMCUtils.sersic_ell_model(mcmc_data, area_tot=area_tot, R_max=R_max)

# ╔═╡ 43ccc903-7857-4207-8652-255339f5ddc2
if run_all
	samples_sersic_ell = sample(mcmc_model_sersic_ell, sampler, MCMCThreads(), Nsamples, Nthreads) 
end

# ╔═╡ 0b7e6509-fd02-47cf-bd45-aa2cafdc18f1
df_sersic_ell = to_frame(samples_sersic_ell)

# ╔═╡ 496ede90-f72f-4d64-b27f-714dc27b0fcb
summary_sersic_ell = MCMCUtils.summarize(samples_sersic_ell)

# ╔═╡ 451850ad-2e06-4609-9a58-626984146b9e
@savefig "corner.sersic_ell" PairPlots.pairplot(samples_sersic_ell)

# ╔═╡ 2439d273-96f9-4ae8-b69e-3add855401ed
@info "number stars sersic: wll $(MCMCUtils.to_measurement(df_sersic_ell.N_sat))"

# ╔═╡ 4979f59b-5297-4596-98d7-b06e84ed3472
write_samples_summary(samples_sersic_ell, df_sersic_ell, "sersic_ell")

# ╔═╡ 2a16c20e-7fdd-4382-b41e-78aafc8e2311
md"""
# A few sanity_check plots
"""

# ╔═╡ 7857111a-9327-4a91-a3fb-f50cb26249b8
let
	f = scatter(stars.xi, stars.eta, markersize=0.5, alpha=0.3, color=:black)

	scatter!(0, 0)
	scatter!(median(df_plummer.d_xi), median(df_plummer.d_eta))

	if run_all
		scatter!(median(df_ell.d_xi), median(df_ell.d_eta))

		scatter!(median(df_sersic.d_xi), median(df_sersic.d_eta))
	end
	
	f
end

# ╔═╡ f82d9b9d-b690-449b-b5c8-2dac5d2526c6
π * R_max^2

# ╔═╡ 5c5f6d05-5e9b-4889-acdc-73bb87e971ee
let
	fig = Figure()
	ax = Axis(fig[1,1])
	prof = LilGuys.SurfaceDensityProfile(stars.R_ell)

	LilGuys.plot_log_Σ!(ax, prof)
	hlines!(log10(N_stars ./ area_tot))
	fig
end

# ╔═╡ 32f87db0-df52-4990-96f6-77a3c3f64f6a
let
	fig = Figure()
	ax = Axis(fig[1,1])
	xi, eta = LilGuys.to_tangent(stars.ra, stars.dec, props_ngc5466["ra"], props_ngc5466["dec"])
	dist_ngc5466 = @. sqrt(xi^2 + eta^2) .* 60

	prof = LilGuys.SurfaceDensityProfile(dist_ngc5466)

	LilGuys.plot_log_Σ!(ax, prof)
	hlines!(log10(N_stars ./ area_tot))

	prof_model = LilGuys.KingProfile(R_s=1.43, R_t=1.43 * 10^1.04, k=5e2)

	@info LilGuys.R_h(prof_model),  props_ngc5466["R_h"]
	x = LinRange(-0.5, 1.5, 100)
	y = LilGuys.surface_density.(prof_model, 10 .^ x)
	lines!(x, log10.(y  .+ 0.9*N_stars ./ area_tot), color=COLORS[2])
	vlines!(log10(R_ngc5466))
	fig
end

# ╔═╡ Cell order:
# ╠═0a974ff8-7c22-49dd-98f4-738537175d20
# ╟─dce2754b-b345-404f-a564-c26438cc835f
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═90c12d18-f1c5-4813-8289-8c030be49ab3
# ╠═17689578-56d6-44a2-955c-be934d41c7b8
# ╠═e10ec7e0-e71f-4290-ab52-189c8763f887
# ╠═a0274b6a-daba-48aa-a289-86063dba9393
# ╠═2ef405a3-9bdf-4255-87d8-901487f2d27d
# ╠═c1816838-9ab4-45bc-b9ae-98157c70a40c
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═066c7b30-f818-4e4a-8db8-c8bac469f558
# ╠═5426e4db-6810-443a-a98b-6e1917fee397
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═626846b5-156a-4883-8542-463e3a90faaa
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═398e5019-006b-4d4d-922c-7ce244470e8e
# ╠═4acaad09-7d87-444e-b5ec-7fcf79f895ef
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═1b77bc6a-0512-4d5e-a817-ce0eab032129
# ╠═a736178a-ae01-4940-afa3-e7dc3522ec47
# ╠═36784213-6393-4e33-8fb9-c4ef5b1fbc48
# ╟─e33170f3-8c57-4672-b8a7-f083d69437c1
# ╠═ea690a36-512e-4a80-a258-faece76c299c
# ╠═efea8cbe-dd33-419a-9a49-e07a3424dbe3
# ╠═82a0c53a-989f-4e0d-81ad-1c49d4f7050e
# ╠═61793e5e-6435-492a-85d4-383dc9cf4f01
# ╠═d17c50db-ce94-4aad-92e2-6355b3b3a91b
# ╠═4aac4483-6f0e-4b73-8c64-cfd3e9995b59
# ╠═0215620c-6193-4215-b018-6f775e31bb5c
# ╠═23a250b8-a7aa-4d0b-b838-a5e4436b0bc1
# ╠═fb28f68a-4893-4feb-aa29-5cc6acd7d2a8
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═fab117b3-53a7-48c8-a6a1-8f004d58a90d
# ╠═14ac9c2a-b218-43c5-a129-3eda0774608a
# ╠═8fa15363-191e-4875-95f1-3ce87762e582
# ╟─e82534ab-257c-4d18-b08e-3c37a05ce623
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╟─cd12f215-cf8c-4d15-93f8-ebb0d5762c18
# ╠═b962db4f-f7b9-40a7-810e-d14c2f6bc39c
# ╠═fe3d7c64-c3ff-431f-976a-9dab7e8e71b3
# ╠═051f456b-5b6e-438f-8438-392456759ae6
# ╠═f721cfdc-e8cb-4fc6-870d-d50cfb671fc6
# ╠═3f592c6f-f35a-4977-9135-2987060dfadb
# ╠═0430ac07-77a7-40a8-86dd-bac3b3a4fe24
# ╠═61e9d125-a0d6-4956-bc78-9ca507598e94
# ╠═1eb86ece-0fc7-453a-aa4e-f8e50cc81229
# ╠═4b1d12dc-ed81-46b2-bfc3-53feab3829d4
# ╠═23c01f4f-7cce-4b15-bcc6-fa0601addf84
# ╠═32e5278a-8674-45fc-af41-a93f009fc008
# ╠═725a1902-2c48-4b14-a3c2-a87afcefa080
# ╠═7819e763-cebe-476b-86ec-641f035bf910
# ╠═09f26d6f-91de-45a0-b3c8-d6f8e0104b53
# ╠═0ddaa045-1b43-405e-af05-066e04bd3918
# ╠═66a6a188-c5f4-4761-96ba-672320727dc9
# ╠═a02fbdd1-b38d-4540-999c-2096bcc50c10
# ╟─c708b169-c593-4c38-a91f-40723249993f
# ╠═7b1b4c0f-aa49-4ee0-b860-0bd927db8768
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═aca9c6b7-b7be-4f79-ba47-438048285041
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
# ╠═7ed82514-9484-47d9-a332-7501c5e5d84d
# ╠═56cc5edd-d56f-434d-9f6e-47e6e00092ed
# ╠═d57d7605-3904-490b-b785-42320275b0c5
# ╠═a3d02bd8-f330-4f13-8a1c-59499537748c
# ╠═1f789712-401a-4942-b499-073d2f25408c
# ╠═43ccc903-7857-4207-8652-255339f5ddc2
# ╠═0b7e6509-fd02-47cf-bd45-aa2cafdc18f1
# ╠═496ede90-f72f-4d64-b27f-714dc27b0fcb
# ╠═451850ad-2e06-4609-9a58-626984146b9e
# ╠═2439d273-96f9-4ae8-b69e-3add855401ed
# ╠═4979f59b-5297-4596-98d7-b06e84ed3472
# ╟─2a16c20e-7fdd-4382-b41e-78aafc8e2311
# ╠═7857111a-9327-4a91-a3fb-f50cb26249b8
# ╠═f82d9b9d-b690-449b-b5c8-2dac5d2526c6
# ╠═5c5f6d05-5e9b-4889-acdc-73bb87e971ee
# ╠═32f87db0-df52-4990-96f6-77a3c3f64f6a
