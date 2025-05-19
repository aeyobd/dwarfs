### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 34728442-3039-11f0-1ae7-2fd35f259853
begin
	import Pkg; Pkg.activate()

	using Turing
	using CairoMakie
	using LilGuys
	using Arya
	using PairPlots
	using PyFITS
end

# ╔═╡ bb8ebb73-431c-4eec-82e6-303fbd1c3bfe
using DataFrames

# ╔═╡ 1327aa72-9e42-408f-8f09-b5a0ffb15a07
galaxy = "fornax"

# ╔═╡ bdd2586e-8cad-4962-81a3-979c273497bd
FIGDIR = joinpath(galaxy, "figures")

# ╔═╡ a703174b-dde4-4ebf-a438-ae03a126aeee
FIGSUFFIX = "_full_fit"

# ╔═╡ 3cfaffe4-4058-4b66-9980-8e7163ed2755
sampler = NUTS(0.65)

# ╔═╡ ce040b6b-834a-48d8-ad53-6eeacd07b2e4
num_samples = 1000

# ╔═╡ 060c6029-5da2-4bfa-b770-8085d50784e3
num_threads = 16

# ╔═╡ 42dec10b-aa8a-4478-bcf6-7ef0bcbc0159
import TOML

# ╔═╡ ba3c9630-72ea-4ea6-b796-6821912f294d
md"""
# Data Loading
"""

# ╔═╡ 9266bd23-ab72-4de0-a6bb-49ea6b4d4308
stars = read_fits(joinpath(galaxy, "samples", "fiducial_sample.fits"))

# ╔═╡ 3090c1e3-fbf6-46d0-be8e-896856abd9ba
scatter(stars.xi, stars.eta)

# ╔═╡ bf46f5ab-ef2b-4271-8809-4e5ce0460896
md"""
# Utils
"""

# ╔═╡ b8b80e8a-3c1b-469b-ae40-2ed57b54df48
function summarize_chain(chain; p=0.16)
	df = DataFrame(chain)
	summary = DataFrame(Turing.summarize(chain))

	Nr = size(summary, 1)
	meds = zeros(Nr)
	err_low = zeros(Nr)
	err_high = zeros(Nr)
	
	for i in 1:Nr
		sym = summary[i, :parameters]
		meds[i] = median(df[!, sym])
		err_low[i] = meds[i] - quantile(df[!, sym], p)
		err_high[i] = quantile(df[!, sym], 1-p) - meds[i]
	end

	summary[!, :median] = meds
	summary[!, :err_low] = err_low
	summary[!, :err_high] = err_high

	select!(summary, :parameters, :median, :err_low, :err_high, Not([:parameters, :median, :err_low, :err_high]))
	return summary
end

# ╔═╡ a301017e-ad94-4a4c-a125-de99f4b990a2
observed_properties = TOML.parsefile(joinpath(galaxy, "observed_properties.toml"))

# ╔═╡ b925dc34-71f8-42c6-a2d2-4ee1871352ab
mod(observed_properties["position_angle"], 180)

# ╔═╡ 8574dee7-38ab-4c44-8a1c-8de0c615881f
if 135 .> mod(observed_properties["position_angle"], 180) .> 45
	theta_prior = Uniform(0, 180)
else
	theta_prior = Uniform(-90, 90)
end

# ╔═╡ 62c3447b-0fd4-4959-ad74-167bc7d615a3
xi_prior = Normal(0, 0.1)

# ╔═╡ cf05bebc-c0eb-4c84-affe-81e0409fbde0
ell_prior = Uniform(0.01, 0.99)

# ╔═╡ 9eb949ce-0114-4506-bc57-bc76c6bdec2d
md"""
# Exp2D
"""

# ╔═╡ c4799b1f-385a-4a45-9b31-93a714ce751d
plot(LinRange(0.001, 60, 1000), LogNormal(2, 2), axis=(; xscale=log10, xticks = Makie.automatic))

# ╔═╡ 57a880c1-ac07-4377-b905-a239095f5f82
median(LogNormal(0, 2))

# ╔═╡ 0e079aba-0941-4b35-a9c8-dcbecc25f047
@model function turing_model_exp(xi, eta) 
	xi0 ~ xi_prior
	eta0 ~ xi_prior
	ell ~ ell_prior
	r_s ~ LogUniform(0.01, 30)
	θ ~ theta_prior

	R = LilGuys.calc_R_ell(xi .- xi0, eta .- eta0, ell, θ)
	prof = LilGuys.Exp2D(R_s=r_s)
	LL = sum(log.(LilGuys.surface_density.(prof, R)))
	Turing.Turing.@addlogprob! LL
end

# ╔═╡ 3a8a4fb0-3b29-445c-bb54-a7e4158f45af
model_exp = turing_model_exp(stars.xi, stars.eta)

# ╔═╡ c17551f7-2ab4-4063-a289-a843652f7e2d
samples_exp =  sample(model_exp, sampler, MCMCThreads(), num_samples, num_threads)

# ╔═╡ aa782702-083c-4a2c-a52f-eb22587e30c2
@savefig "exp2d" pairplot(samples_exp)

# ╔═╡ 239b64ba-a3ae-4866-8fff-a403698e45d6
summary_exp = summarize_chain(samples_exp)

# ╔═╡ 2ab6fde1-09c7-4cc3-b276-ee5c6b17cae0
md"""
# Inner exp
"""

# ╔═╡ 55d83cbb-00e7-477f-bc19-2106ae6f9c4b
@model function turing_model_inner_exp(xi, eta, R_inner=3) 
	xi0 ~ xi_prior
	eta0 ~ xi_prior
	ell ~ ell_prior
	r_s ~ LogUniform(0.01, 30)
	θ ~ theta_prior

	R = LilGuys.calc_R_ell(xi .- xi0, eta .- eta0, ell, θ)
	prof = LilGuys.Exp2D(R_s=r_s)
	LL = sum(log.(LilGuys.surface_density.(prof, R)) .* (R .< R_inner*r_s))
	Turing.Turing.@addlogprob! LL
end

# ╔═╡ 509c0bc0-493a-4fd6-a89c-0d41bd90a31f
model_inner_exp = turing_model_inner_exp(stars.xi, stars.eta)

# ╔═╡ 7180df2b-c97c-4858-ab50-b3432287c385
#samples_inner_exp =  sample(model_inner_exp, sampler, MCMCThreads(), num_samples, num_threads)

# ╔═╡ f30fcc17-bc61-48da-844f-772a1bd89bd4
@savefig "exp2d_inner" pairplot(samples_inner_exp)

# ╔═╡ faa2feb9-eded-4422-9881-834e951ee15d
summary_inner_exp = summarize_chain(samples_inner_exp)

# ╔═╡ aa9fab4f-2585-45ba-9819-da9d0d6cdbc5
md"""
# Plummer
"""

# ╔═╡ 357fb1f8-bf05-416a-8b4a-5bd912bcef85
@model function turing_model_plummer(xi, eta) 
	xi0 ~ xi_prior
	eta0 ~ xi_prior
	ell ~ ell_prior
	r_s ~ LogUniform(0.01, 30)
	θ ~ theta_prior

	R = LilGuys.calc_R_ell(xi .- xi0, eta .- eta0, ell, θ)
	prof = LilGuys.Plummer(r_s=r_s, M=1.0)
	LL = sum(log.(LilGuys.surface_density.(prof, R)))
	Turing.Turing.@addlogprob! LL
end

# ╔═╡ c33dfd66-9167-4602-b448-01120871ba53
model_plummer = turing_model_plummer(stars.xi, stars.eta)

# ╔═╡ f8a6de56-db13-403b-8795-bad01c360546
samples_plummer =  sample(model_plummer, sampler, MCMCThreads(), num_samples, num_threads)

# ╔═╡ 3a0defd6-9731-4c1f-809f-5d1cf54c1fb6
@savefig "plummer" pairplot(samples_plummer)

# ╔═╡ e9c0810f-4b96-4c7c-870d-c9380bbe174c
summary_plummer = summarize_chain(samples_plummer)

# ╔═╡ 0af2dab3-bc2f-42b9-97a9-dac38a328b3f
md"""
# Sersic
"""

# ╔═╡ 11745b10-80a6-40db-9d6c-eb61cb753630


# ╔═╡ 969b3ce6-5404-4b05-b0ff-4eb4d44ca161
@model function turing_model_sersic(xi, eta) 
	xi0 ~ xi_prior
	eta0 ~ xi_prior
	ell ~ ell_prior
	R_h ~ LogUniform(0.01, 30)
	n ~ Uniform(0, 10)
	θ ~ theta_prior

	R = LilGuys.calc_R_ell(xi .- xi0, eta .- eta0, ell, θ)
	prof = LilGuys.Sersic(R_h=R_h, n=n,_b_n=LilGuys.guess_b_n(n))
	M_scale = LilGuys.mass_2D(prof, Inf)
	prof = LilGuys.Sersic(R_h=R_h, n=n, _b_n=LilGuys.guess_b_n(n), Σ_h=1/M_scale)
	LL = sum(log.(LilGuys.surface_density.(prof, R)))
	Turing.Turing.@addlogprob! LL
end

# ╔═╡ a12339b8-e9ec-48d9-8256-9c84ddde6abe
model_sersic = turing_model_sersic(stars.xi, stars.eta)

# ╔═╡ 958f30e9-7a3a-4442-af54-132a4a8cf850
samples_sersic =  sample(model_sersic, sampler, MCMCThreads(), num_samples, num_threads)

# ╔═╡ 57b244f8-432a-40e8-8a86-db785b6ae026
@savefig "sersic" pairplot(samples_sersic)

# ╔═╡ 47bc2a43-7e62-4a81-ac43-65d6f5926eb8
summary_sersic = summarize_chain(samples_sersic)

# ╔═╡ 050b5981-3f91-4881-873a-97eb349f66d5
all_fits = OrderedDict(
	"sersic" => summary_sersic,
	"exp2d" => summary_exp,
	"plummer" => summary_plummer
)

# ╔═╡ 3235ccea-6db9-4412-98fa-a6de0c4a1ae3
minimum(DataFrame(samples_sersic).lp), mean(DataFrame(samples_sersic).lp)

# ╔═╡ e9da18a6-3f55-41a8-a7cf-e06010c984f2
minimum(DataFrame(samples_exp).lp), mean(DataFrame(samples_exp).lp)

# ╔═╡ 9e095080-dcc9-474b-8e03-c2abdf81c08b
minimum(DataFrame(samples_plummer).lp), mean(DataFrame(samples_plummer).lp)

# ╔═╡ ddb7b8a4-a10e-48ec-9645-405a0478bbf9
begin
	derived_props = OrderedDict()

	for (name, fit) in all_fits
		df = OrderedDict()

		for argname in fit.parameters
			i = findfirst(fit.parameters .== argname)
			df["$(argname)"] = fit.median[i]
			df["$(argname)_em"] = fit.err_low[i]
			df["$(argname)_ep"] = fit.err_high[i]
		end
		derived_props[name] = df

	end

	derived_props
end

# ╔═╡ 8f383aa2-1f76-4634-8996-db2f2396c748
outfile = joinpath(galaxy, "samples", "full_density_fits.toml")

# ╔═╡ 522bfa40-c204-4adc-bc58-de5de964367e
open(outfile, "w") do f
	TOML.print(f, derived_props)
end

# ╔═╡ cba7941f-cf14-46eb-9b65-1f07bac70a2d
TOML.parsefile(outfile)

# ╔═╡ Cell order:
# ╠═1327aa72-9e42-408f-8f09-b5a0ffb15a07
# ╠═bdd2586e-8cad-4962-81a3-979c273497bd
# ╠═a703174b-dde4-4ebf-a438-ae03a126aeee
# ╠═3cfaffe4-4058-4b66-9980-8e7163ed2755
# ╠═ce040b6b-834a-48d8-ad53-6eeacd07b2e4
# ╠═060c6029-5da2-4bfa-b770-8085d50784e3
# ╠═34728442-3039-11f0-1ae7-2fd35f259853
# ╠═42dec10b-aa8a-4478-bcf6-7ef0bcbc0159
# ╠═bb8ebb73-431c-4eec-82e6-303fbd1c3bfe
# ╟─ba3c9630-72ea-4ea6-b796-6821912f294d
# ╠═9266bd23-ab72-4de0-a6bb-49ea6b4d4308
# ╠═3090c1e3-fbf6-46d0-be8e-896856abd9ba
# ╟─bf46f5ab-ef2b-4271-8809-4e5ce0460896
# ╠═b8b80e8a-3c1b-469b-ae40-2ed57b54df48
# ╠═a301017e-ad94-4a4c-a125-de99f4b990a2
# ╠═b925dc34-71f8-42c6-a2d2-4ee1871352ab
# ╠═8574dee7-38ab-4c44-8a1c-8de0c615881f
# ╠═62c3447b-0fd4-4959-ad74-167bc7d615a3
# ╠═cf05bebc-c0eb-4c84-affe-81e0409fbde0
# ╟─9eb949ce-0114-4506-bc57-bc76c6bdec2d
# ╠═c4799b1f-385a-4a45-9b31-93a714ce751d
# ╠═57a880c1-ac07-4377-b905-a239095f5f82
# ╠═0e079aba-0941-4b35-a9c8-dcbecc25f047
# ╠═3a8a4fb0-3b29-445c-bb54-a7e4158f45af
# ╠═c17551f7-2ab4-4063-a289-a843652f7e2d
# ╠═aa782702-083c-4a2c-a52f-eb22587e30c2
# ╠═239b64ba-a3ae-4866-8fff-a403698e45d6
# ╠═2ab6fde1-09c7-4cc3-b276-ee5c6b17cae0
# ╠═55d83cbb-00e7-477f-bc19-2106ae6f9c4b
# ╠═509c0bc0-493a-4fd6-a89c-0d41bd90a31f
# ╠═7180df2b-c97c-4858-ab50-b3432287c385
# ╠═f30fcc17-bc61-48da-844f-772a1bd89bd4
# ╠═faa2feb9-eded-4422-9881-834e951ee15d
# ╟─aa9fab4f-2585-45ba-9819-da9d0d6cdbc5
# ╠═357fb1f8-bf05-416a-8b4a-5bd912bcef85
# ╠═c33dfd66-9167-4602-b448-01120871ba53
# ╠═f8a6de56-db13-403b-8795-bad01c360546
# ╠═3a0defd6-9731-4c1f-809f-5d1cf54c1fb6
# ╠═e9c0810f-4b96-4c7c-870d-c9380bbe174c
# ╟─0af2dab3-bc2f-42b9-97a9-dac38a328b3f
# ╠═11745b10-80a6-40db-9d6c-eb61cb753630
# ╠═969b3ce6-5404-4b05-b0ff-4eb4d44ca161
# ╠═a12339b8-e9ec-48d9-8256-9c84ddde6abe
# ╠═958f30e9-7a3a-4442-af54-132a4a8cf850
# ╠═57b244f8-432a-40e8-8a86-db785b6ae026
# ╠═47bc2a43-7e62-4a81-ac43-65d6f5926eb8
# ╠═050b5981-3f91-4881-873a-97eb349f66d5
# ╠═3235ccea-6db9-4412-98fa-a6de0c4a1ae3
# ╠═e9da18a6-3f55-41a8-a7cf-e06010c984f2
# ╠═9e095080-dcc9-474b-8e03-c2abdf81c08b
# ╠═ddb7b8a4-a10e-48ec-9645-405a0478bbf9
# ╠═8f383aa2-1f76-4634-8996-db2f2396c748
# ╠═522bfa40-c204-4adc-bc58-de5de964367e
# ╠═cba7941f-cf14-46eb-9b65-1f07bac70a2d
