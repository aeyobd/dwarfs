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
	
	using Turing
end

# ╔═╡ adf5db84-58a5-4160-886a-1414087d1db7
using PyFITS

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
	samplename = ARGS[2]
	Nsamples = 1_000
	Nthreads = 16
else

	Nsamples = 1000
	Nthreads = 1

	md"""
	$(@bind galaxy confirm(TextField(default="bootes3")))
	
	$(@bind samplename confirm(TextField(default="delve_bhb_sample")))
	"""
end

# ╔═╡ 5c342791-a615-40da-8be8-8cb03d705645
readdir("samples")

# ╔═╡ 05517bcc-7967-4bc7-9396-c271e420665d
import PairPlots

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
outdir = joinpath("..", galaxy, "mcmc")

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	using LilGuys

	FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".$samplename.mcmc_sersic_simple"
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 066c7b30-f818-4e4a-8db8-c8bac469f558
module MCMCUtils
	include("../mcmc/mcmc_utils.jl")
end

# ╔═╡ 461ea0d0-9a33-4f82-8a99-d4bc936dee2d
module GaiaPlots
	include("../../thesis/plots/gaia_utils.jl")
end

# ╔═╡ 5426e4db-6810-443a-a98b-6e1917fee397
module GaiaFilters
	include("../../utils/gaia_filters.jl")
end

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ 7e8124ea-7bbe-465b-a9dc-4b14d268c39e
obs_props = let
	df = MCMCUtils.get_obs_props(galaxy)
	df["ra"] += 0.0
	df["dec"] -= 0.0
	df["ellipticity"] = 0.0
	df
end

# ╔═╡ 398e5019-006b-4d4d-922c-7ce244470e8e
obs_props["ra"]

# ╔═╡ 4acaad09-7d87-444e-b5ec-7fcf79f895ef
allstars = let
	df = read_fits("samples/$samplename.fits")
	df[!, :R_ell] = @. sqrt(df.xi^2 + df.eta^2)
	df
end

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = GaiaFilters.calc_R_max(allstars.xi, allstars.eta, 0, 0)

# ╔═╡ a736178a-ae01-4940-afa3-e7dc3522ec47
stars = allstars[allstars.R_ell .< R_max, :]

# ╔═╡ 67b865b5-6d7d-4a62-84cd-82983a76f8ba
data = stars

# ╔═╡ 93cba051-259f-46da-b99c-c837d3bba496


# ╔═╡ db29e430-3b71-47bd-a3a0-adbe9ab91b9a
scatter(stars.bp_rp, stars.G, markersize=1, alpha=0.2)

# ╔═╡ f35d18e3-27c7-423f-a6fa-440b841dbad8
hexbin(data.xi, data.eta, bins=50)

# ╔═╡ 4aac4483-6f0e-4b73-8c64-cfd3e9995b59
scatter(data.xi, data.eta, markersize=1, alpha=0.3, color=:black)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Sersic model
"""

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function ell_model(data::DataFrame, R_max::Real=R_max)
	d_xi  ~ Normal(0, 30)
	d_eta ~ Normal(0, 30)
	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	
	R_h ~ LogNormal(2, 2)
	# n ~ Uniform(0.4, 12)
	n = 1
	f_sat ~ Uniform(0, 1)
	

	R = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity, position_angle)

	prof = LilGuys.Plummer(r_s=R_h, M=1)
	mtot = LilGuys.mass_2D(prof, R_max)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R)
	L_bg_space = 1 / (π  * R_max^2)
	
	LL = sum(@. log10.(
		(1-f_sat) * L_bg_space
		+ f_sat *  L_sat_space
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ d6d5f794-fb61-4aa1-a181-b12a7f9e3ab0
@model function plummer_model(data::DataFrame, R_max::Real=R_max)
	d_xi  ~ Normal(0, 30)
	d_eta ~ Normal(0, 30)

	R_h ~ LogNormal(2, 2)
	# n ~ Uniform(0.4, 12)
	n = 1
	f_sat ~ Uniform(0, 1)
	

	R = @. sqrt((data.xi + d_xi)^2 + (data.eta + d_eta)^2)
	
	prof = LilGuys.Plummer(r_s=R_h, M=1)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R)
	L_bg_space = 1 / (π  * R_max^2)
	
	LL = sum(@. log10.(
		(1-f_sat) * L_bg_space
		+ f_sat *  L_sat_space
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ c7182ecf-d75e-44ec-b273-3e692cb5efa0
@model function sersic_model(data::DataFrame, R_max::Real=R_max)
	d_xi  ~ Normal(0, 30)
	d_eta ~ Normal(0, 30)
	R_h ~ LogNormal(2, 2)
	n ~ Uniform(0.4, 12)
	f_sat ~ Uniform(0, 1)
	
	R = @. sqrt((data.xi + d_xi)^2 + (data.eta + d_eta)^2)
	
	b_n = LilGuys.guess_b_n(n)
	prof = LilGuys.Sersic(n=n, R_h=R_h, _b_n=b_n)

	mtot = LilGuys.mass_2D(prof, R_max)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R) / mtot
	L_bg_space = 1 / (π  * R_max^2)
	
	LL = sum(@. log10.(
		(1-f_sat) * L_bg_space
		+ f_sat *  L_sat_space
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ a774454c-1c3b-4c35-b753-bfb03cc9245e
@model function ell_sersic_model(data::DataFrame, R_max::Real=R_max)
	d_xi  ~ Normal(0, 30)
	d_eta ~ Normal(0, 30)
	R_h ~ LogNormal(2, 2)
	n ~ Uniform(0.4, 12)
	f_sat ~ Uniform(0, 1)

	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	

	
	# R = @. sqrt((data.xi + d_xi)^2 + (data.eta + d_eta)^2)
	R = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity, position_angle)

	b_n = LilGuys.guess_b_n(n)
	prof = LilGuys.Sersic(n=n, R_h=R_h, _b_n=b_n)

	mtot = LilGuys.mass_2D(prof, R_max)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R) / mtot
	L_bg_space = 1 / (π  * R_max^2)
	
	LL = sum(@. log10.(
		(1-f_sat) * L_bg_space
		+ f_sat *  L_sat_space
	))
	
	Turing.@addlogprob!(LL)
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

# ╔═╡ fe3d7c64-c3ff-431f-976a-9dab7e8e71b3
mcmc_model_plummer = plummer_model(data)

# ╔═╡ 051f456b-5b6e-438f-8438-392456759ae6
samples_plummer = sample(mcmc_model_plummer, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ f721cfdc-e8cb-4fc6-870d-d50cfb671fc6
df_plummer = DataFrame(samples_plummer)

# ╔═╡ 3f592c6f-f35a-4977-9135-2987060dfadb
PairPlots.pairplot(samples_plummer)

# ╔═╡ 2b838aad-e5c8-431b-9b5d-5d4de82284b1
median(df_plummer.f_sat) * size(stars, 1)

# ╔═╡ 8671e6bd-43db-41be-96c4-1aae2dfd240d
quantile(df_plummer.f_sat, [0.16, 0.84]) * size(stars, 1)

# ╔═╡ 23c01f4f-7cce-4b15-bcc6-fa0601addf84
md"""
### Elliptical Plummer model
"""

# ╔═╡ 32e5278a-8674-45fc-af41-a93f009fc008
mcmc_model_ell = ell_model(data)

# ╔═╡ 725a1902-2c48-4b14-a3c2-a87afcefa080
samples_ell = sample(mcmc_model_ell, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ 7819e763-cebe-476b-86ec-641f035bf910
df_ell = DataFrame(samples_ell)

# ╔═╡ b07a4109-0395-404e-8fcd-59475bded7f9
stars.L_CMD_SAT_HB ./ stars.L_CMD_BKD

# ╔═╡ 09f26d6f-91de-45a0-b3c8-d6f8e0104b53
PairPlots.pairplot(samples_ell)

# ╔═╡ c708b169-c593-4c38-a91f-40723249993f
md"""
### Sersic model
"""

# ╔═╡ 7b1b4c0f-aa49-4ee0-b860-0bd927db8768
mcmc_model = sersic_model(data)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
samples = sample(mcmc_model, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ 0dea434b-8c46-401f-819f-6735fc318a3a


# ╔═╡ c6659d9d-f303-430d-818d-23499541011d
import NaNMath

# ╔═╡ aca9c6b7-b7be-4f79-ba47-438048285041
df_out = DataFrame(samples)

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
summary = MCMCUtils.summarize(samples)

# ╔═╡ 428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
@savefig "corner" PairPlots.pairplot(samples)

# ╔═╡ 14ac9c2a-b218-43c5-a129-3eda0774608a
function plot_corner(samples)
	df = DataFrame(samples)
	df[!, :log_R_h] = log10.(df.R_h)

	cols = [v for v in values(samples.info.varname_to_symbol)]
	cols[cols .== :R_h] .= :log_R_h

	PairPlots.pairplot(df[:, cols])
end

# ╔═╡ 2501b91a-1999-4284-bf5d-c898de2f271c
plot_corner(samples)

# ╔═╡ d57d7605-3904-490b-b785-42320275b0c5
@info "acceptance rate (all, per-step)"  mean(df_out.acceptance_rate), mean(df_out.is_accept)

# ╔═╡ 3416eb76-9c80-40d7-b49f-3ead0cf58f75
md"""
## Plots
"""

# ╔═╡ 01be3284-b52a-4f1c-acb4-2430d4085e37
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	scatter!(data.xi, data.eta, markersize=1, color=:black, 
	   )

	scatter!(0, 0)
	scatter!(-median(df_out.d_xi), -median(df_out.d_eta))
	scatter!(-median(df_plummer.d_xi), -median(df_plummer.d_eta))
	scatter!(-median(df_ell.d_xi), -median(df_ell.d_eta))
	
	fig
end

# ╔═╡ a16da68d-ecb7-40e3-9fad-5fe966990b23
prof_obs = LilGuys.SurfaceDensityProfile(stars.R_ell, bins=(1:0.1:2.5)) |> LilGuys.filter_empty_bins


# ╔═╡ ed5ab012-c855-4bea-93bc-44aa4fcc1bbf
prof_obs.log_Sigma

# ╔═╡ 2b780e99-54e3-4b90-beb9-85501449ff74
function plot_sersic!(df_out)


	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Sersic(R_h=df_out.R_h[i], n=df_out.n[i])
		M = size(stars, 1) * df_out.f_sat[i] / LilGuys.mass_2D(prof, R_max)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(x, y, color=COLORS[2], alpha=0.1)
	end

	ylims!(-6, 3)



end

# ╔═╡ 178a0f5f-a2ba-4834-8359-369eb215a64d
function plot_plummer!(df_out; color=COLORS[2])

	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Plummer(r_s=df_out.R_h[i], M=1)
		M = size(stars, 1) * df_out.f_sat[i]
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(x, y, color=color, alpha=0.1)
	end

end

# ╔═╡ 2314bd82-b6fd-4ee5-898b-3f0582074ca1
log_Σ_bg = median(prof_obs.log_Sigma[prof_obs.log_R .> 1.8])

# ╔═╡ 603b1716-8d78-43b3-bc14-c114cdae35d6
Sigma_sub = 10 .^ prof_obs.log_Sigma .- 10 .^ log_Σ_bg

# ╔═╡ 02976736-5020-4db6-941d-2f66e65eee9a
NaNMath.log10.(Sigma_sub .- LilGuys.lower_error.(Sigma_sub))

# ╔═╡ ff2534ca-c619-4503-84c0-e3567fdd788d
filt_bkg = NaNMath.log10.(Sigma_sub .- LilGuys.lower_error.(Sigma_sub)) .|> isfinite |> LilGuys.find_longest_consecutive_true

# ╔═╡ bafc69b0-a266-4789-82d5-6c9d99b4977e
begin
	log_Sigma_sub = log10.(Sigma_sub[filt_bkg])
	prof_sub = LilGuys.SurfaceDensityProfile(
		R_units=prof_obs.R_units,
		log_R=prof_obs.log_R[filt_bkg],
		log_R_bins=prof_obs.log_R_bins[filt_bkg],
		log_Sigma=log_Sigma_sub,
		log_R_scale=prof_obs.log_R_scale,
		log_m_scale=prof_obs.log_m_scale,
		annotations=prof_obs.annotations
	)
end

# ╔═╡ 9dfad57b-dbaa-4a46-8fa6-8115ef9c11cd
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=log_r_label, ylabel=log_Sigma_label)


	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))


	
	plot_sersic!(df_out)
	plot_plummer!(df_plummer, color=COLORS[3])
	# plot_plummer!(df_ell, color=COLORS[4])
	ylims!(-4, -0.5)
	hlines!(middle.((log_Σ_bg)))

	errorscatter!(prof_sub.log_R, (prof_sub.log_Sigma), yerror=error_interval.(prof_sub.log_Sigma), color=:black, label="bg-subtracted")

	lines!([NaN], [NaN], color=COLORS[2], label="Sersic")
	lines!([NaN], [NaN], color=COLORS[3], label="Plummer")
	axislegend(position=:lb)


	fig
end

# ╔═╡ 3d236a62-46a8-4ff6-a711-50840c17c32d
prof_sub.log_Sigma

# ╔═╡ 6e072f15-9699-47c7-b92e-473a3dce6ffd
function plot_plummer(df_out)
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=log_r_label, ylabel=log_Sigma_label)


	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))

	plot_plummer!(df_out)

	ylims!(-6, 3)

	fig

end

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fontsize=4
	Nc = size(samples, 3)
	Nvar = size(summary, 1)

	fig = Figure(size=(2*72, Nvar/4*72),
		yminorticksvisible=false
	)

	for i in 1:Nvar
		ax = Axis(fig[i, 1], 
			ylabelsize=fontsize, 
			xlabelsize=fontsize,
		    xticklabelsize=fontsize,
		    yticklabelsize=fontsize,
			ylabel=summary.parameters[i],
			xlabel="step",
			yminorticksvisible=false,
			ylabelrotation=0,
		)
		
		for c in 1:Nc
			y = samples[:, i, c]
			lines!((y), linewidth=0.1)
		end



		ax2 = Axis(fig[i, 2])
		hist!(vec(samples[:, i, :]), direction=:x)
		
		if i < Nvar
			hidexdecorations!(ax)
			
		end
		hidexdecorations!(ax2)
		hideydecorations!(ax2)
		linkyaxes!(ax, ax2)
	end
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)
	colsize!(fig.layout, 2, Relative(1/4))

	@savefig "chains"
	fig

end

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ 19b570dc-2339-459d-9ba9-ec0139ccd0be
median(DataFrame(samples).f_sat) * size(data, 1)

# ╔═╡ ef5eb016-4922-4935-ac86-a1b30bc9e727
median(DataFrame(samples_ell).f_sat) * size(data, 1)

# ╔═╡ f2e64722-2040-4308-9915-05f36fb43093
median(DataFrame(samples_plummer).f_sat) * size(data, 1)

# ╔═╡ d1009491-62de-42d7-89ad-b41b8385aaaf
samplesout = joinpath(outdir, "samples$FIGSUFFIX.csv")

# ╔═╡ a93579f2-d37d-492c-95e6-02bd7491dc8c
summaryout = joinpath(outdir, "summary$FIGSUFFIX.csv")

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
CSV.write(summaryout, summary)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═5c342791-a615-40da-8be8-8cb03d705645
# ╠═05517bcc-7967-4bc7-9396-c271e420665d
# ╠═adf5db84-58a5-4160-886a-1414087d1db7
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═066c7b30-f818-4e4a-8db8-c8bac469f558
# ╠═461ea0d0-9a33-4f82-8a99-d4bc936dee2d
# ╠═5426e4db-6810-443a-a98b-6e1917fee397
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═398e5019-006b-4d4d-922c-7ce244470e8e
# ╠═4acaad09-7d87-444e-b5ec-7fcf79f895ef
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═a736178a-ae01-4940-afa3-e7dc3522ec47
# ╠═67b865b5-6d7d-4a62-84cd-82983a76f8ba
# ╠═93cba051-259f-46da-b99c-c837d3bba496
# ╠═db29e430-3b71-47bd-a3a0-adbe9ab91b9a
# ╠═f35d18e3-27c7-423f-a6fa-440b841dbad8
# ╠═4aac4483-6f0e-4b73-8c64-cfd3e9995b59
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═d6d5f794-fb61-4aa1-a181-b12a7f9e3ab0
# ╠═c7182ecf-d75e-44ec-b273-3e692cb5efa0
# ╠═a774454c-1c3b-4c35-b753-bfb03cc9245e
# ╟─e82534ab-257c-4d18-b08e-3c37a05ce623
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╠═cd12f215-cf8c-4d15-93f8-ebb0d5762c18
# ╠═fe3d7c64-c3ff-431f-976a-9dab7e8e71b3
# ╠═051f456b-5b6e-438f-8438-392456759ae6
# ╠═f721cfdc-e8cb-4fc6-870d-d50cfb671fc6
# ╠═3f592c6f-f35a-4977-9135-2987060dfadb
# ╠═2b838aad-e5c8-431b-9b5d-5d4de82284b1
# ╠═8671e6bd-43db-41be-96c4-1aae2dfd240d
# ╠═23c01f4f-7cce-4b15-bcc6-fa0601addf84
# ╠═32e5278a-8674-45fc-af41-a93f009fc008
# ╠═725a1902-2c48-4b14-a3c2-a87afcefa080
# ╠═7819e763-cebe-476b-86ec-641f035bf910
# ╠═b07a4109-0395-404e-8fcd-59475bded7f9
# ╠═09f26d6f-91de-45a0-b3c8-d6f8e0104b53
# ╟─c708b169-c593-4c38-a91f-40723249993f
# ╠═7b1b4c0f-aa49-4ee0-b860-0bd927db8768
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═0dea434b-8c46-401f-819f-6735fc318a3a
# ╠═c6659d9d-f303-430d-818d-23499541011d
# ╠═ed5ab012-c855-4bea-93bc-44aa4fcc1bbf
# ╠═9dfad57b-dbaa-4a46-8fa6-8115ef9c11cd
# ╠═3d236a62-46a8-4ff6-a711-50840c17c32d
# ╠═aca9c6b7-b7be-4f79-ba47-438048285041
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
# ╠═2501b91a-1999-4284-bf5d-c898de2f271c
# ╠═14ac9c2a-b218-43c5-a129-3eda0774608a
# ╠═d57d7605-3904-490b-b785-42320275b0c5
# ╟─3416eb76-9c80-40d7-b49f-3ead0cf58f75
# ╠═01be3284-b52a-4f1c-acb4-2430d4085e37
# ╠═2b780e99-54e3-4b90-beb9-85501449ff74
# ╠═178a0f5f-a2ba-4834-8359-369eb215a64d
# ╠═2314bd82-b6fd-4ee5-898b-3f0582074ca1
# ╠═603b1716-8d78-43b3-bc14-c114cdae35d6
# ╠═02976736-5020-4db6-941d-2f66e65eee9a
# ╠═ff2534ca-c619-4503-84c0-e3567fdd788d
# ╠═bafc69b0-a266-4789-82d5-6c9d99b4977e
# ╠═a16da68d-ecb7-40e3-9fad-5fe966990b23
# ╠═6e072f15-9699-47c7-b92e-473a3dce6ffd
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═19b570dc-2339-459d-9ba9-ec0139ccd0be
# ╠═ef5eb016-4922-4935-ac86-a1b30bc9e727
# ╠═f2e64722-2040-4308-9915-05f36fb43093
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
