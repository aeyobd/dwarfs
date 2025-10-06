### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ b71fe082-4ef3-41db-ac41-69dee16981a6
md"""
# Utilities
"""

# ╔═╡ 9f06a5c4-a8ae-49d2-991e-1c1576361700
function get_R_h(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors
	R_h = obs_props["R_h"]
	return R_h
end

# ╔═╡ 40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
function load_profile(galaxyname; algname="fiducial", inner=false)
	@info galaxyname
	
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_profile.toml")
	
    prof = LilGuys.SurfaceDensityProfile(filename) |> LilGuys.filter_empty_bins
	# remove large errors
	filt = maximum.(LilGuys.error_interval.(prof.log_Sigma)).< 1
	prof.log_Sigma[.!filt] .= Measurement(NaN, NaN)
	prof = prof |> LilGuys.filter_empty_bins

	
	R_h = get_R_h(galaxyname)
	
	# rescale to R_h and \Sigma_R_h
	Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(middle(R_h)))
	M_s = Σ_h * middle(R_h) .^ 2

    prof = LilGuys.scale(prof, 1/middle(R_h), 1/M_s)

	prof
end

# ╔═╡ 8254a9ac-b612-4df1-8ab1-59c499b9006e
md"""
# Data loading
"""

# ╔═╡ 8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
prof_scl = load_profile("sculptor")

# ╔═╡ c3b52d5a-a8b4-4207-a8c7-9d32914aca93
prof_umi = load_profile("ursa_minor")

# ╔═╡ 740d04a4-4327-47cd-b5aa-b71ed94a0610
prof_fornax = load_profile("fornax")

# ╔═╡ 552e9438-862f-4710-a7d4-c8d798b5f1aa
galaxynames = [
	"fornax",
	"leo1",
	# "antlia2", this one is unreliable...
	 "leo2",
	 "carina",
	 "draco",	
	"canes_venatici1",
	 "sextans1",
	"crater2",
]

# ╔═╡ f4e8b66c-5f18-45fd-8859-32479d7227bc
profiles = let
	profiles = OrderedDict(name => load_profile(name) for name in galaxynames)
	profiles["sculptor"] = prof_scl
	profiles["ursa_minor"] = prof_umi
	profiles
end

# ╔═╡ 971b2120-dfdf-4820-95ab-9b264b566bcf
md"""
# Plots
"""

# ╔═╡ 07e1390f-2789-4e20-80d6-b2cb4631be5c
galaxies_scatter = OrderedDict(
	"fornax" => (
		color = COLORS[3],
		marker = :circle,
		label = "Fornax",
	),
	"sculptor" => (
		color = COLORS[2], 
		marker=:rect,
		label = "Sculptor",
	),
	"ursa_minor" => (
		color = COLORS[4],
		marker = :utriangle,
		label = "Ursa Minor",
		markersize=8,
	),
)

# ╔═╡ dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
plummer = let 
	prof = LilGuys.Plummer()
	Σ_0 = surface_density(prof, 1)
	prof = LilGuys.Plummer(M = 1/Σ_0)
end

# ╔═╡ 5b5c426e-ce27-48bd-a3ac-7f23e7663ed7
let
	fig = Figure(size=(3.5*72, 2.5*72))

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
		xlabel = L"\log\,R\ / \ R_h",
		#xticks = -2:1:1,
		limits=(-1.0, 0.9, -3.5, 1.0)

	)

	# exponential reference
    x = LinRange(-1.5, 1.2, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="2D exponential", 
		   linewidth=0.8*theme(:linewidth)[])
    

	for (galaxy, kwargs) in galaxies_scatter
		prof = profiles[galaxy]

		yerr = LilGuys.error_interval.(prof.log_Sigma)
		errorscatter!(ax, prof.log_R, prof.log_Sigma, yerror=yerr,
					  markersize=theme(:markersize)[]*0.8; kwargs...)
	end

	axislegend(position=:lb, merge=true, unique=true, patchsize=(24, 6))



	@savefig "scl_umi_fornax_exp"
	fig
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
				xticks = -2:1:1,
		limits=(-1.5, nothing, -4.2, 1.1)

	)

    x = LinRange(-1.4, 1, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y, color=:black, label="Plummer", linestyle=:dash)


	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="classical dwarfs")
	end


	scatterlines!(prof_scl.log_R, prof_scl.log_Sigma, color=COLORS[2], label="Sculptor", )
	scatterlines!(prof_umi.log_R, prof_umi.log_Sigma, color=COLORS[3], label="Ursa Minor", )

	axislegend(position=:lb, merge=true, unique=true)



	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta \log\,\Sigma",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.5, 1, -1.2, 2.2)
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)

	
	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5)
	end

	scatterlines!(prof_scl.log_R, prof_scl.log_Sigma .- f(prof_scl.log_R), color=COLORS[2])
	scatterlines!(prof_umi.log_R, prof_umi.log_Sigma .- f(prof_umi.log_R), color=COLORS[3])


	linkxaxes!(ax, ax_res)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/4))
	hidexdecorations!(ax, ticks=false, minorticks=false)

	@savefig "classical_dwarf_profiles"
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═67dd3488-79e1-4ba2-b9ac-f417765d55de
# ╠═cf655228-1528-4e77-9629-07e99984951f
# ╠═6d2bff6d-b347-49c4-87df-ad58d8a27ff3
# ╠═b71fe082-4ef3-41db-ac41-69dee16981a6
# ╠═9f06a5c4-a8ae-49d2-991e-1c1576361700
# ╠═40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
# ╠═8254a9ac-b612-4df1-8ab1-59c499b9006e
# ╠═8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
# ╠═c3b52d5a-a8b4-4207-a8c7-9d32914aca93
# ╠═740d04a4-4327-47cd-b5aa-b71ed94a0610
# ╠═552e9438-862f-4710-a7d4-c8d798b5f1aa
# ╠═f4e8b66c-5f18-45fd-8859-32479d7227bc
# ╟─971b2120-dfdf-4820-95ab-9b264b566bcf
# ╠═07e1390f-2789-4e20-80d6-b2cb4631be5c
# ╠═dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
# ╠═5b5c426e-ce27-48bd-a3ac-7f23e7663ed7
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
