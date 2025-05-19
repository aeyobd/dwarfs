### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using LilGuys
	using Turing
end

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	galaxy = "reticulum2"
end

# ╔═╡ af1d827c-bb7e-4875-936e-cfec7578f7db
md"""
Prior on log sigma relative. Can be fairly loose
"""

# ╔═╡ 72209cca-6101-4905-a286-0489e9247b94
prior = Uniform(-12, 6)

# ╔═╡ 8b3a3bd5-4126-4991-a27b-49939f90fecb
pvalue = 0.16

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 28550d54-ad88-4df5-ac04-f46a15588ff8
import DensityEstimators: bin_indices

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
module MCMCUtils
	include("mcmc_utils.jl")
end

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 54932cde-bc26-420f-9b36-8b926fa93f84
outdir = joinpath("..", galaxy, "mcmc")

# ╔═╡ 8c565a28-84bc-4bc7-8a0e-b2e0dff76665
if !isdir(joinpath(outdir))
	mkdir(joinpath(outdir))
end

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".mcmc_hist_fast"

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
obs_props = MCMCUtils.get_obs_props(galaxy)

# ╔═╡ b2adcf46-c538-43d0-b08b-193175a4c68a
struct_params = MCMCUtils.StructuralParams(; LilGuys.dict_to_tuple(TOML.parsefile(joinpath("..", galaxy, "mcmc", "default_bins.toml")))...)

# ╔═╡ c58ed939-07c1-4d34-90a4-3c5d9cc9910b
stars = MCMCUtils.get_fits(galaxy, obs_props)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
#  model
"""

# ╔═╡ 12d15434-e919-460c-ad1d-a1880ed6562f
bins = struct_params.bins

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function hist_model(Lsat::Vector{Float64}, Lbkd::Vector{Float64})
	log_Σ ~ prior

	Σ = 10 .^ log_Σ

	f = Σ / (1 + Σ)
	LL = sum(@. log10(f*Lsat + (1-f) * Lbkd) )
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 562fe2c7-3ed4-4adb-af7f-1eca2aa35d8b
Nbins = length(bins) - 1	

# ╔═╡ ba5ceef7-ac10-46f7-a89b-30f38b6ddec1
r_b = bin_indices(stars.R_ell, bins)

# ╔═╡ df84d941-7ddf-4c1c-96bb-1858b49bb710
Lsat = stars.L_CMD_SAT .* stars.L_PM_SAT

# ╔═╡ 97cbf4b6-025d-461b-b82e-044f86b713c1
Lbg = stars.L_CMD_BKD .* stars.L_PM_BKD

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
begin 
	chains = []
	for i in 1:Nbins
		filt = r_b .== i
		model = hist_model(Lsat[filt], Lbg[filt])
		chain = sample(model, NUTS(0.65), MCMCThreads(), 1000, 16) 
		push!(chains, chain)
	end
end

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fontsize=4
	Nvar = length(bins)-1

	fig = Figure(size=(2*72, Nvar/4*72),
		yminorticksvisible=false
	)

	for i in 1:Nvar
		ax = Axis(fig[i, 1], 
			ylabelsize=fontsize, 
			xlabelsize=fontsize,
		  xticklabelsize=fontsize,
		  yticklabelsize=fontsize,
			ylabel=L"\theta_{%$i}",
			xlabel="step",
			yminorticksvisible=false,
		)
		
		for c in 1:size(chains[1], 3)
			y = chains[i][:, 1, c]
			lines!((y), linewidth=0.1)
		end



		ax2 = Axis(fig[i, 2])
		hist!(vec(chains[i][:, 1, :]), direction=:x)
		
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

# ╔═╡ 370146c4-0414-467e-8a3b-ff3e0fdd09e0
md"""
# processing
"""

# ╔═╡ 22d1fd06-c689-41d9-bc16-aedde1c85444
begin 
	summaries = vcat([DataFrame(summarize(chain)) for chain in chains]...)
	summaries[!, "log_R"] = midpoints(log10.(bins))
	
	log_Sigma_median = zeros(Nbins)
	log_Sigma_low = zeros(Nbins)
	log_Sigma_high = zeros(Nbins)
	for i in 1:length(bins)-1
		log_Sigma_low[i], log_Sigma_median[i], log_Sigma_high[i] = quantile(chains[i][:, 1, :], [pvalue, 0.5, 1-pvalue])
	end

	summaries[!, :median] = log_Sigma_median
	summaries[!, :lower_error] = log_Sigma_median .- log_Sigma_low
	summaries[!, :upper_error] =  log_Sigma_high .- log_Sigma_median
	summaries[!, :parameters] = string.(summaries.parameters) .* string.(1:Nbins)

	summaries
end

# ╔═╡ 67b5b69a-2b2b-4431-955a-3829bd19f13c
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		
	)
		
	errorscatter!(summaries.log_R, summaries.median, yerror = collect(zip(summaries.lower_error, summaries.upper_error)))

	@savefig "derived parameters"

	fig
end

# ╔═╡ 36db7e1d-4b48-4510-99d5-7d567ac70d5d
begin 
	df_out = DataFrame()
	df_out[!, "iteration"] = DataFrame(chains[1])[!, "iteration"]
	df_out[!, "chain"] = DataFrame(chains[1])[!, "chain"]
	df_out[!, "lp"] .= 0

	for i in eachindex(chains)
		df = DataFrame(chains[i])
		df_out[!, "params[$i]"] = df[!, "log_Σ"]
		df_out[!, "lp"] .+= df.lp
		@assert df_out.chain == df.chain
		@assert df_out.iteration == df.iteration
	end

	df_out
end

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ d1009491-62de-42d7-89ad-b41b8385aaaf
samplesout = joinpath(outdir, "samples$FIGSUFFIX.csv")

# ╔═╡ 639c7d31-8693-45dd-88be-492b124804e9
summaryout = joinpath(outdir, "summary$FIGSUFFIX.csv")

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 7eee7173-6924-45bf-bc91-bf9625854a38
CSV.write(summaryout, summaries)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╟─af1d827c-bb7e-4875-936e-cfec7578f7db
# ╠═72209cca-6101-4905-a286-0489e9247b94
# ╠═8b3a3bd5-4126-4991-a27b-49939f90fecb
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═28550d54-ad88-4df5-ac04-f46a15588ff8
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═54932cde-bc26-420f-9b36-8b926fa93f84
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═f01b0b37-5fdf-4c72-92cb-7b1cbd6c88b9
# ╠═b2adcf46-c538-43d0-b08b-193175a4c68a
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═12d15434-e919-460c-ad1d-a1880ed6562f
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═562fe2c7-3ed4-4adb-af7f-1eca2aa35d8b
# ╠═ba5ceef7-ac10-46f7-a89b-30f38b6ddec1
# ╠═df84d941-7ddf-4c1c-96bb-1858b49bb710
# ╠═97cbf4b6-025d-461b-b82e-044f86b713c1
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═67b5b69a-2b2b-4431-955a-3829bd19f13c
# ╟─370146c4-0414-467e-8a3b-ff3e0fdd09e0
# ╠═22d1fd06-c689-41d9-bc16-aedde1c85444
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═639c7d31-8693-45dd-88be-492b124804e9
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═7eee7173-6924-45bf-bc91-bf9625854a38
