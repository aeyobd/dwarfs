### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	import PyPlot as plt
	import PyCall

end

# ╔═╡ cf655228-1528-4e77-9629-07e99984951f
using OrderedCollections

# ╔═╡ 67dd3488-79e1-4ba2-b9ac-f417765d55de
import TOML

# ╔═╡ 444aa05a-a838-44fc-b1d0-1ea8fdfc19f8
PyCall.@pyimport arya

# ╔═╡ 7bcd6483-72b3-4148-b22d-964e599c91e9
arya.style.init("apj")

# ╔═╡ d2310706-3d20-4de8-848b-f76192a3e40d
COLORS = arya.COLORS

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
		color = COLORS[3],
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

# ╔═╡ 76d51ea5-8b21-4a4f-b2e5-c0adcc673df3
mlw = 2 * 2/3

# ╔═╡ 17113383-254e-44e2-b7c3-2900dec2cd7d
smallfontsize =  8

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
let
	fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(3.35, 3),
	                              gridspec_kw=Dict("height_ratios" => [3, 1], 
												  "hspace" => 0))
		

	plt.sca(ax)
	for (galaxy, prof) in profiles
		plt.plot((prof.log_R), middle.(prof.log_Sigma), color=COLORS[1], alpha=0.5)
	end

	plt.plot([], [], color=COLORS[1], alpha=0.5, label="classical dwarfs")

	

	plt.plot(prof_scl.log_R, middle.(prof_scl.log_Sigma), color=COLORS[2], label="Sculptor", marker="s")
	plt.plot(prof_umi.log_R, middle.(prof_umi.log_Sigma), color=COLORS[3], label="Ursa Minor", marker="^")


    x = LinRange(-1.4, 1, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    plt.plot(x, y, color="k", label="2D exponential", lw=mlw)
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    plt.plot(x, y, ls="--", color="k", label="Plummer", lw=mlw)

	
	plt.legend()



	plt.sca(ax_res)

	ax_res.annotate("residual against 2D exponential", 
					(0.0, 1.0), 
					xytext=(smallfontsize, -smallfontsize/2),
					textcoords = "offset points",
					xycoords = "axes fraction",
	            
	            va="top", ha="left",
	            fontsize=smallfontsize,
				
	            )

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	

	
	for (galaxy, prof) in profiles
		plt.plot(prof.log_R, middle.(prof.log_Sigma) .- f(prof.log_R), color=COLORS[1], alpha=0.5,)
	end

	plt.plot(prof_scl.log_R, middle.(prof_scl.log_Sigma) .- f(prof_scl.log_R), color=COLORS[2], marker="s")
	plt.plot(prof_umi.log_R, middle.(prof_umi.log_Sigma) .- f(prof_umi.log_R), color=COLORS[3], marker="^")

    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    plt.plot(x, y .- f(x), color=:black, label="Exp2D",linewidth=mlw)
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    plt.plot(x, y .- f(x), color=:black, label="Plummer", linestyle="--", linewidth=mlw)



	# ax_res setup
	ax_res.set_ylabel(raw"$\Delta\log\,\Sigma$")
	ax_res.set_xlabel(raw"$\log\,R_\mathrm{ell}\ /\ R_h$")
	ax_res.set_xlim(-1.5, 1.0)
	ax_res.set_ylim(-1.2, 2.2)
	ax_res.set_xticks(-1:1:1)
	ax_res.tick_params(right=false)	

	
	# Primary (bottom) x-axis — linear in log10 space
	ax.set_xlabel(raw"$\log\,R_\mathrm{ell}\ /\ R_h$")
	ax.set_ylabel(raw"$\log\,\Sigma\ /\ \Sigma_h$")
	ax.set_xlim(-1.5, 1.0)
	ax.set_ylim(-4.2, 1.1)
	ax.set_xticks(-1:1:1)
	ax.tick_params(top=false, right=false)
	
	# # Secondary (top) x-axis — shows the actual R_ell/R_h values
	ax2 = ax.twiny()
	ax2.set_xscale("log")
	ax2.set_xlim(10^-1.5, 10^1)
	ax2.set_xticks([0.1, 1, 10])
	ax2.set_xticklabels(["0.1", "1", "10"])
	ax2.set_xlabel(raw"$R_\mathrm{ell}\ /\ R_h$")
	ax2.tick_params(bottom=false, right=false)
	
	plt.savefig("figures/classical_dwarf_profiles.pdf")
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═67dd3488-79e1-4ba2-b9ac-f417765d55de
# ╠═cf655228-1528-4e77-9629-07e99984951f
# ╠═444aa05a-a838-44fc-b1d0-1ea8fdfc19f8
# ╠═7bcd6483-72b3-4148-b22d-964e599c91e9
# ╠═d2310706-3d20-4de8-848b-f76192a3e40d
# ╟─b71fe082-4ef3-41db-ac41-69dee16981a6
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
# ╠═76d51ea5-8b21-4a4f-b2e5-c0adcc673df3
# ╠═17113383-254e-44e2-b7c3-2900dec2cd7d
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
