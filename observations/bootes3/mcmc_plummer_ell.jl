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
	mag_max = parse(Float64, ARGS[1])
	Nt = 16
	N_samples = 3_000
	write = true
else
	Nt = 1
	N_samples = 1_000
	write = false

	@bind mag_max confirm(NumberField(17:0.1:22, default=20))
end

# ╔═╡ 55b15061-f400-4feb-b259-fddd5d25b34e
module MCMCUtils
	include("./mcmc_utils.jl")
end

# ╔═╡ 05517bcc-7967-4bc7-9396-c271e420665d
import PairPlots

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
outdir = joinpath("mcmc")

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	# import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".G_$mag_max.mcmc_plummer_ell"
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ 7e8124ea-7bbe-465b-a9dc-4b14d268c39e
obs_props = let
	df = MCMCUtils.get_obs_props("bootes3")
	df["ra"] = df["ra_original"]
	df["dec"] = df["dec_original"]
	df["ellipticity"] = 0
	df
end

# ╔═╡ 00380333-4c2f-4c32-9521-2764cffef265
stars = let
	df = MCMCUtils.get_fits("bootes3", obs_props)
	df[df.G .< mag_max, :]
end

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = maximum(stars.R_ell)

# ╔═╡ 67b865b5-6d7d-4a62-84cd-82983a76f8ba
data = MCMCUtils.GaiaData(stars)

# ╔═╡ 6e895965-0ddf-4044-8131-a5965e4b2e11
N_stars = size(stars, 1)

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ df3ab039-55ee-4431-bd68-1c77f843dd18
sampler = NUTS()

# ╔═╡ 7b1b4c0f-aa49-4ee0-b860-0bd927db8768
mcmc_model = plummer_model(data)

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
samples = sample(mcmc_model, sampler, MCMCThreads(), N_samples, Nt) 

# ╔═╡ d57d7605-3904-490b-b785-42320275b0c5
@info mean(df_out.acceptance_rate), mean(df_out.is_accept)

# ╔═╡ 36db7e1d-4b48-4510-99d5-7d567ac70d5d
df_out = let
	df = DataFrame(samples)
	df[!, :N_memb] = N_stars * df.f_sat
	df
end

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
chain_summary = MCMCUtils.summarize(samples)

# ╔═╡ 51e22997-a5ec-48b3-9719-9ae0e26cc20c
prof_obs = LilGuys.SurfaceDensityProfile(stars.R_ell[stars.LLR_nospace .> 2.0], bins=0.5:0.1:2.7) |> LilGuys.filter_empty_bins

# ╔═╡ 91bc2e1b-6337-41b2-97d3-20f3719f066f
log_Σ_bg = median(prof_obs.log_Sigma[end-2:end])

# ╔═╡ 2b780e99-54e3-4b90-beb9-85501449ff74
let
	fig = Figure()
	ax = Axis(fig[1,1])


	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))


	
	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	df = df_out
	for i in 1:100
		M = size(stars, 1) * df.f_sat[i]
		prof = LilGuys.Plummer(r_s=df.R_h[i], M=M)
		y = @. log10(LilGuys.surface_density(prof, R) .+ 10^log_Σ_bg.middle)

		lines!(x, y, color=COLORS[3], alpha=0.1)
	end

	ylims!(-4, -1)

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

	if write
		@savefig "chains"
	end
	fig

end

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ bb9dc499-8ff5-4dae-8e17-dc2346276822
median(df_out.N_memb)

# ╔═╡ d1009491-62de-42d7-89ad-b41b8385aaaf
samplesout = joinpath(outdir, "samples$FIGSUFFIX.csv")

# ╔═╡ a93579f2-d37d-492c-95e6-02bd7491dc8c
summaryout = joinpath(outdir, "summary$FIGSUFFIX.csv")

# ╔═╡ f0611371-d69d-4683-b24b-36c607da2857
if write
	CSV.write(samplesout, df_out)
	CSV.write(summaryout, chain_summary)
end

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45


# ╔═╡ 4c0d5b06-d99b-41de-8a42-e29dbd6e0e53


# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═55b15061-f400-4feb-b259-fddd5d25b34e
# ╠═05517bcc-7967-4bc7-9396-c271e420665d
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═00380333-4c2f-4c32-9521-2764cffef265
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═67b865b5-6d7d-4a62-84cd-82983a76f8ba
# ╠═6e895965-0ddf-4044-8131-a5965e4b2e11
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╠═7b1b4c0f-aa49-4ee0-b860-0bd927db8768
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═d57d7605-3904-490b-b785-42320275b0c5
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═51e22997-a5ec-48b3-9719-9ae0e26cc20c
# ╠═91bc2e1b-6337-41b2-97d3-20f3719f066f
# ╠═2b780e99-54e3-4b90-beb9-85501449ff74
# ╠═428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═bb9dc499-8ff5-4dae-8e17-dc2346276822
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═f0611371-d69d-4683-b24b-36c607da2857
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
