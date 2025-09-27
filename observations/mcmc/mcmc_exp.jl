### A Pluto.jl notebook ###
# v0.20.18

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
else
	@bind galaxy confirm(TextField(default="leo2"))
end

# ╔═╡ 05517bcc-7967-4bc7-9396-c271e420665d
import PairPlots

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
outdir = joinpath("..", galaxy, "mcmc")

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	# import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".mcmc_exp"
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 066c7b30-f818-4e4a-8db8-c8bac469f558
module MCMCUtils
	include("mcmc_utils.jl")
end

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ 7e8124ea-7bbe-465b-a9dc-4b14d268c39e
obs_props = MCMCUtils.get_obs_props(galaxy)

# ╔═╡ 00380333-4c2f-4c32-9521-2764cffef265
stars = MCMCUtils.get_fits(galaxy, obs_props)

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = maximum(stars.R_ell)

# ╔═╡ 67b865b5-6d7d-4a62-84cd-82983a76f8ba
data = MCMCUtils.GaiaData(stars)

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function exp_model(data::MCMCUtils.GaiaData)
	d_xi ~ Normal(0, 5)
	d_eta ~ Normal(0, 5)
	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	
	R_s ~ LogNormal(0, 5)
	prof = LilGuys.Exp2D(R_s=R_s, M=1)

	R = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity, position_angle)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R)
	L_bg_space = 1 / (π*R_max^2)

	f_sat ~ Uniform(0, 1)
	LL = sum(@. log10.(
		(1-f_sat) * data.L_bg * L_bg_space
		+ f_sat * data.L_sat * L_sat_space
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ df3ab039-55ee-4431-bd68-1c77f843dd18
sampler = NUTS()

# ╔═╡ 7b1b4c0f-aa49-4ee0-b860-0bd927db8768
mcmc_model = exp_model(data)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
samples = sample(mcmc_model, sampler, MCMCThreads(), 1000, 16) 

# ╔═╡ 36db7e1d-4b48-4510-99d5-7d567ac70d5d
df_out = DataFrame(samples)

# ╔═╡ d57d7605-3904-490b-b785-42320275b0c5
@info mean(df_out.acceptance_rate), mean(df_out.is_accept)

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
chain_summary = MCMCUtils.summarize(samples)

# ╔═╡ d1600a1c-1a73-4a0e-a36f-3008a5e9ed23
prof_obs = LilGuys.SurfaceDensityProfile(stars.R_ell[stars.PSAT .> 0.5]) |> LilGuys.filter_empty_bins

# ╔═╡ 51e22997-a5ec-48b3-9719-9ae0e26cc20c
prof_all = LilGuys.SurfaceDensityProfile(stars.R_ell)

# ╔═╡ 2b780e99-54e3-4b90-beb9-85501449ff74
let
	fig = Figure()
	ax = Axis(fig[1,1])


	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))

	errorscatter!(prof_all.log_R, (prof_all.log_Sigma), yerror=error_interval.(prof_all.log_Sigma))

	
	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	df = df_out
	for i in 1:100
		M = size(stars, 1) * df.f_sat[i]
		prof = LilGuys.Exp2D(R_s=df.R_s[i], M=M)
		y = @. log10(LilGuys.surface_density(prof, R))

		lines!(x, y, color=COLORS[3], alpha=0.1)
	end

	ylims!(-6, 3)

	fig

end

# ╔═╡ 428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
PairPlots.pairplot(samples)

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fontsize=4
	Nvar = size(chain_summary, 1)

	fig = Figure(size=(2*72, Nvar/4*72),
		yminorticksvisible=false
	)

	for i in 1:Nvar
		ax = Axis(fig[i, 1], 
			ylabelsize=fontsize, 
			xlabelsize=fontsize,
		  xticklabelsize=fontsize,
		  yticklabelsize=fontsize,
			ylabel=chain_summary.parameters[i],
			xlabel="step",
			yminorticksvisible=false,
			ylabelrotation=0,
		)
		
		for c in 1:size(samples.value, 3)
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
CSV.write(summaryout, chain_summary)

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
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═00380333-4c2f-4c32-9521-2764cffef265
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═67b865b5-6d7d-4a62-84cd-82983a76f8ba
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╠═7b1b4c0f-aa49-4ee0-b860-0bd927db8768
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═d57d7605-3904-490b-b785-42320275b0c5
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═d1600a1c-1a73-4a0e-a36f-3008a5e9ed23
# ╠═51e22997-a5ec-48b3-9719-9ae0e26cc20c
# ╠═2b780e99-54e3-4b90-beb9-85501449ff74
# ╠═428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
