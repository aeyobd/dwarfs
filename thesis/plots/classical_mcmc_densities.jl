### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ cf655228-1528-4e77-9629-07e99984951f
using OrderedCollections

# ╔═╡ 6d2bff6d-b347-49c4-87df-ad58d8a27ff3
include("./paper_style.jl")

# ╔═╡ 67dd3488-79e1-4ba2-b9ac-f417765d55de
import TOML

# ╔═╡ 4ee54995-a31f-4ae1-8204-e55884d786bb
α = LilGuys.R_h(LilGuys.Exp2D(R_s=1))

# ╔═╡ 40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
function load_profile(galaxyname; algname="../mcmc/hist_fast")
	# if galaxyname ∈ ["crater2", "antlia2"]
	# 	algname = "mcmc_hist_fast"
	# else
	# 	algname = "mcmc_hist"
	# end
	# filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
	# 	"processed/profile.$(algname).toml")
	@info galaxyname
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_profile.toml")
	
    prof = LilGuys.StellarDensityProfile(filename) 


	density_fit = TOML.parsefile(filename * "_inner_fits.toml")
	R_h = 10 .^ density_fit["log_R_s_exp2d_inner"] * α
	R_h_u = α * 10 ^ LilGuys.Measurement(density_fit["log_R_s_exp2d_inner"], density_fit["log_R_s_exp2d_inner_em"], density_fit["log_R_s_exp2d_inner_ep"])
	@info "R_h = $R_h_u arcmin"
	@info "counts = $(sum(prof.counts))"

	Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(R_h))
	M_s = Σ_h * R_h .^ 2

	
    prof = LilGuys.scale(prof, 1/R_h, 1/M_s)

	filt = maximum.(LilGuys.error_interval.(prof.log_Sigma)).< 1

	prof.log_R = prof.log_R[filt]
	prof.log_Sigma = prof.log_Sigma[filt]

	prof
end

# ╔═╡ 552e9438-862f-4710-a7d4-c8d798b5f1aa
galaxynames = [
	"fornax",
	"leo1",
	"leo2",
	"carina",
	"sextans1",
	"draco",	
]

# ╔═╡ f4e8b66c-5f18-45fd-8859-32479d7227bc
profiles = OrderedDict(name => load_profile(name) for name in galaxynames)

# ╔═╡ 8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
prof_scl = load_profile("sculptor")

# ╔═╡ c3b52d5a-a8b4-4207-a8c7-9d32914aca93
prof_umi = load_profile("ursa_minor")

# ╔═╡ b31bfb7e-8550-4f97-8b42-91b92edaa255
CairoMakie.activate!(type="svg", pt_per_unit=2)

# ╔═╡ dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
begin 
	plummer = LilGuys.Plummer()
	Σ_0 = surface_density(plummer, 1)
	plummer = LilGuys.Plummer(M = 1/surface_density(plummer, 1))
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
				xticks = -2:1:1,
		limits=(nothing, nothing, -3.1, 1)

	)

    x = LinRange(-1.4, 1, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y, color=:black, label="Plummer", linestyle=:dash)


	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="classical dwarfs")
	end


	scatterlines!(prof_scl.log_R, prof_scl.log_Sigma, color=COLORS[2], label="Sculptor", linewidth=1, markersize=3)
	scatterlines!(prof_umi.log_R, prof_umi.log_Sigma, color=COLORS[3], label="Ursa Minor", linewidth=1, markersize=3)

	axislegend(position=:lb, merge=true, unique=true)



	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta \log\,\Sigma",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.4, 1, -1, 1)
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)

	
	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5)
	end

	scatterlines!(prof_scl.log_R, prof_scl.log_Sigma .- f(prof_scl.log_R), color=COLORS[2], linewidth=1, markersize=3)
	scatterlines!(prof_umi.log_R, prof_umi.log_Sigma .- f(prof_umi.log_R), color=COLORS[3], linewidth=1, markersize=3)


	linkxaxes!(ax, ax_res)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/4))
	hidexdecorations!(ax, ticks=false, minorticks=false)

	#@savefig "classical_dwarf_profiles"
	fig
end

# ╔═╡ f706fc01-fade-4d4f-97d3-ce08264ec680
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
				xticks = -2:1:1,
		limits=(nothing, nothing, -4.2, 1.2)

	)


	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.2, label="classical dwarfs")
	end


    x = LinRange(-1.4, 1, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Sérsic (n=1)")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

	axislegend(position=:lb, merge=true, unique=true)



	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta \log\,\Sigma",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.4, 1, -1, 2)
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.2)
	end


    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Sérsic (n=1)")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


	linkxaxes!(ax, ax_res)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/4))
	hidexdecorations!(ax, ticks=false, minorticks=false)

	# @savefig "classical_dwarf_profiles"
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═67dd3488-79e1-4ba2-b9ac-f417765d55de
# ╠═4ee54995-a31f-4ae1-8204-e55884d786bb
# ╠═40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
# ╠═552e9438-862f-4710-a7d4-c8d798b5f1aa
# ╠═cf655228-1528-4e77-9629-07e99984951f
# ╠═f4e8b66c-5f18-45fd-8859-32479d7227bc
# ╠═8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
# ╠═c3b52d5a-a8b4-4207-a8c7-9d32914aca93
# ╠═b31bfb7e-8550-4f97-8b42-91b92edaa255
# ╠═6d2bff6d-b347-49c4-87df-ad58d8a27ff3
# ╠═dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═f706fc01-fade-4d4f-97d3-ce08264ec680
