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
	galaxy = "draco"
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
FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".mcmc_hist_struct"

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


# ╔═╡ df84d941-7ddf-4c1c-96bb-1858b49bb710
Lsat = stars.L_CMD_SAT .* stars.L_PM_SAT

# ╔═╡ 97cbf4b6-025d-461b-b82e-044f86b713c1
Lbg = stars.L_CMD_BKD .* stars.L_PM_BKD

# ╔═╡ 8ad26d05-6551-4a18-a01a-4803225a2ff2
Nstruct = 1024

# ╔═╡ 17d57cb5-5b57-4c1c-9ee3-939d409a2536
Nchains = 1

# ╔═╡ ef076143-de9e-4768-a64e-7fb4f88e6ec5
Nsteps = 100

# ╔═╡ c3fe069d-c4dd-4b4e-a0b0-100f6f4d253d
import Logging

# ╔═╡ 5c1b8277-971e-41e7-b70c-20a14297ce60
logger = Logging.SimpleLogger(Logging.Error)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
begin 
	chain_dfs = Vector{DataFrame}(undef, Nstruct)
	Threads.@threads for i in 1:Nstruct

		params = struct_params
		dξ = rand(Normal(0, params.position_err))
		dη = Normal(0, params.position_err) |> rand
		ellipticity = truncated(Normal(params.ellipticity, params.ellipticity_err), 
			lower=0, upper=0.99) |> rand
		position_angle =  Normal(params.position_angle, params.position_angle_err) |> rand
		
		radii = LilGuys.calc_R_ell(stars.xi .+ dξ, stars.eta .+ dη, ellipticity, position_angle)
		r_b = bin_indices(radii, bins)
		dfs = DataFrame()
		@info "sample $i"
		
		for j in 1:Nbins
			#@info "bin $j"
			filt = r_b .== j
			
			model = hist_model(Lsat[filt], Lbg[filt])
			
			chains = Logging.with_logger(logger) do
				mapreduce(c -> sample(model, NUTS(0.65), Nsteps, progress=false, verbose=false) , chainscat, 1:Nchains)
			end
				
			df = DataFrame(chains)
			df[:, :bin] .= j
			append!(dfs, df)
		end

		dfs[:, :d_xi] .= dξ
		dfs[:, :d_eta] .= dη
		dfs[:, :ellipticity] .= ellipticity
		dfs[:, :position_angle] .= position_angle
		dfs[:, :sample_struct] .= i
		chain_dfs[i] = dfs
	end
end

# ╔═╡ 370146c4-0414-467e-8a3b-ff3e0fdd09e0
md"""
# processing
"""

# ╔═╡ 36db7e1d-4b48-4510-99d5-7d567ac70d5d
df_out = vcat(chain_dfs...)

# ╔═╡ 89cd85c9-d588-49cb-9f0f-19f1e66f1bec
df_out

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═╡ disabled = true
#=╠═╡
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
		
		for c in 1:Nchains
			for s in 1:Nstruct
				y = df_out[(df_out.chain .== c) .& (df_out.sample_struct .== s) .& (df_out.bin .== i), :log_Σ]
				lines!((y), linewidth=0.1)
			end
		end



		ax2 = Axis(fig[i, 2])
		
		y = df_out[df_out.bin .== i, :log_Σ]

		hist!(y, direction=:x)
		
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
  ╠═╡ =#

# ╔═╡ 5320bb1c-c866-44a0-a451-c159b93d8805
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
			xlabel="ell",
			yminorticksvisible=false,
		)

		df = df_out[df_out.bin .== i, :]
		y = df.log_Σ
		x = df.ellipticity

		hist2d!(x, y)

		ax2 = Axis(fig[i, 2], xlabel = "PA",
			xlabelsize=fontsize,
		  xticklabelsize=fontsize,
		  yticklabelsize=fontsize,
				  )
		x = df.position_angle
		hist2d!(x, y)


		ax3 = Axis(fig[i, 3], xlabel = "d xi",	xlabelsize=fontsize,
		  xticklabelsize=fontsize,
		  yticklabelsize=fontsize,)
		x = df.d_xi
		hist2d!(x, y)

	
		ax4 = Axis(fig[i, 4], xlabel = "d eta",	xlabelsize=fontsize,
		  xticklabelsize=fontsize,
		  yticklabelsize=fontsize,)
		x = df.d_eta
		hist2d!(x, y)
		
		if i < Nvar
			hidexdecorations!(ax)
			hidexdecorations!(ax2)
			hidexdecorations!(ax3)
			hidexdecorations!(ax4)

		end
		hideydecorations!(ax2)
		hideydecorations!(ax3)
		hideydecorations!(ax4)
		linkyaxes!(ax, ax2, ax3, ax4)
	end
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)

	@savefig "chains_coor"
	fig

end

# ╔═╡ 22d1fd06-c689-41d9-bc16-aedde1c85444
begin 
	summaries = DataFrame()
		#vcat([DataFrame(summarize(chain)) for chain in chains]...)
	summaries[!, "log_R"] = midpoints(log10.(bins))
	
	log_Sigma_median = zeros(Nbins)
	log_Sigma_low = zeros(Nbins)
	log_Sigma_high = zeros(Nbins)
	for i in 1:length(bins)-1
		y = df_out[df_out.bin .== i, :log_Σ]
		log_Sigma_low[i], log_Sigma_median[i], log_Sigma_high[i] = quantile(y, [pvalue, 0.5, 1-pvalue])
	end

	summaries[!, :median] = log_Sigma_median
	summaries[!, :lower_error] = log_Sigma_median .- log_Sigma_low
	summaries[!, :upper_error] =  log_Sigma_high .- log_Sigma_median
	#summaries[!, :parameters] = string.(summaries.parameters) .* string.(1:Nbins)

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
# ╠═8ad26d05-6551-4a18-a01a-4803225a2ff2
# ╠═17d57cb5-5b57-4c1c-9ee3-939d409a2536
# ╠═ef076143-de9e-4768-a64e-7fb4f88e6ec5
# ╠═c3fe069d-c4dd-4b4e-a0b0-100f6f4d253d
# ╠═5c1b8277-971e-41e7-b70c-20a14297ce60
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═89cd85c9-d588-49cb-9f0f-19f1e66f1bec
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═5320bb1c-c866-44a0-a451-c159b93d8805
# ╠═67b5b69a-2b2b-4431-955a-3829bd19f13c
# ╟─370146c4-0414-467e-8a3b-ff3e0fdd09e0
# ╠═22d1fd06-c689-41d9-bc16-aedde1c85444
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═639c7d31-8693-45dd-88be-492b124804e9
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═7eee7173-6924-45bf-bc91-bf9625854a38
