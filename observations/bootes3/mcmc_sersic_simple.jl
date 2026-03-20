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

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
	samplename = ARGS[2]
	Nsamples = 1_000
	Nthreads = 1
else

	Nsamples = 1000
	Nthreads = 16

	md"""
	$(@bind galaxy confirm(TextField(default="bootes3")))
	
	$(@bind samplename confirm(TextField(default="delve_loose")))
	"""
end

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

# ╔═╡ 14c7788e-0762-49f0-b120-69003cd6a321


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

# ╔═╡ d9d6cb26-d575-4740-b69b-b6503e8d6959
params = let
	p = GaiaFilters.read_paramfile(
	joinpath("../$galaxy/density_profiles/$samplename.toml")
) 
	pop!(p, "profile_kwargs", Dict())

	p["filename"] = joinpath("..", galaxy, "density_profiles", p["filename"])
	p |> LilGuys.dict_to_tuple
end

# ╔═╡ ab6286c2-7b8b-4db5-a8fc-2b24a7bec04f
filt_params = GaiaFilters.GaiaFilterParams(obs_props; params... )

# ╔═╡ 00380333-4c2f-4c32-9521-2764cffef265
all_stars = GaiaFilters.read_gaia_stars(filt_params)

# ╔═╡ 4acaad09-7d87-444e-b5ec-7fcf79f895ef
stars = GaiaFilters.select_members(all_stars, filt_params)

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = maximum(stars.R_ell)

# ╔═╡ 67b865b5-6d7d-4a62-84cd-82983a76f8ba
data = stars

# ╔═╡ db29e430-3b71-47bd-a3a0-adbe9ab91b9a
scatter(stars.bp_rp, stars.G, markersize=1, alpha=0.2)

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Sersic model
"""

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function sersic_model(data::DataFrame, R_max::Real=R_max)
	d_xi  ~ Normal(0, 30)
	d_eta ~ Normal(0, 30)
	# ellipticity ~ Uniform(0, 0.99)
	# position_angle ~ Uniform(0, 180)
	
	R_h ~ LogNormal(2, 2)
	# n ~ Uniform(0.4, 12)
	n = 1
	f_sat ~ Uniform(0, 1)
	

	# R = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity, position_angle)
	R = @. sqrt((data.xi + d_xi)^2 + (data.eta + d_eta)^2)
	
	# b_n = LilGuys.guess_b_n(n)
	# prof = LilGuys.Sersic(n=n, R_h=R_h, _b_n=b_n)

	prof = LilGuys.Plummer(r_s=R_h)
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

# ╔═╡ 7b1b4c0f-aa49-4ee0-b860-0bd927db8768
mcmc_model = sersic_model(data)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
samples = sample(mcmc_model, sampler, MCMCThreads(), Nsamples, Nthreads) 

# ╔═╡ f35d18e3-27c7-423f-a6fa-440b841dbad8
hexbin(data.xi[data.R_ell .< 150], data.eta[data.R_ell .< 150], bins=50)

# ╔═╡ afba6878-065d-4aa5-96fc-8ad44a097444
md"""
## Analyze
"""

# ╔═╡ aca9c6b7-b7be-4f79-ba47-438048285041
df_out = DataFrame(samples)

# ╔═╡ 01be3284-b52a-4f1c-acb4-2430d4085e37
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	scatter!(data.xi, data.eta, markersize=1, color=:black, 
	   )

	scatter!(0, 0)
	scatter!(-median(df_out.d_xi), -median(df_out.d_eta))
	
	fig
end

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
summary = MCMCUtils.summarize(samples)

# ╔═╡ d57d7605-3904-490b-b785-42320275b0c5
@info "acceptance rate (all, per-step)"  mean(df_out.acceptance_rate), mean(df_out.is_accept)

# ╔═╡ 3416eb76-9c80-40d7-b49f-3ead0cf58f75
md"""
## Plots
"""

# ╔═╡ 2b780e99-54e3-4b90-beb9-85501449ff74
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=log_r_label, ylabel=log_Sigma_label)

	prof_obs = LilGuys.SurfaceDensityProfile(stars.R_ell) |> LilGuys.filter_empty_bins

	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))

	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Plummer(r_s=df_out.R_h[i])
		M = size(stars, 1) * df_out.f_sat[i] / LilGuys.mass_2D(prof, R_max)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(x, y, color=COLORS[2], alpha=0.1)
	end

	ylims!(-6, 3)

	fig

end

# ╔═╡ 428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
@savefig "corner" PairPlots.pairplot(samples)

# ╔═╡ 19b570dc-2339-459d-9ba9-ec0139ccd0be
median(DataFrame(samples).f_sat) * size(data, 1)

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
# ╠═05517bcc-7967-4bc7-9396-c271e420665d
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
# ╠═14c7788e-0762-49f0-b120-69003cd6a321
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═398e5019-006b-4d4d-922c-7ce244470e8e
# ╠═d9d6cb26-d575-4740-b69b-b6503e8d6959
# ╠═ab6286c2-7b8b-4db5-a8fc-2b24a7bec04f
# ╠═00380333-4c2f-4c32-9521-2764cffef265
# ╠═4acaad09-7d87-444e-b5ec-7fcf79f895ef
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═67b865b5-6d7d-4a62-84cd-82983a76f8ba
# ╠═db29e430-3b71-47bd-a3a0-adbe9ab91b9a
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╟─e82534ab-257c-4d18-b08e-3c37a05ce623
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╠═7b1b4c0f-aa49-4ee0-b860-0bd927db8768
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═f35d18e3-27c7-423f-a6fa-440b841dbad8
# ╠═01be3284-b52a-4f1c-acb4-2430d4085e37
# ╟─afba6878-065d-4aa5-96fc-8ad44a097444
# ╠═aca9c6b7-b7be-4f79-ba47-438048285041
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═d57d7605-3904-490b-b785-42320275b0c5
# ╟─3416eb76-9c80-40d7-b49f-3ead0cf58f75
# ╠═2b780e99-54e3-4b90-beb9-85501449ff74
# ╠═428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
# ╠═19b570dc-2339-459d-9ba9-ec0139ccd0be
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
