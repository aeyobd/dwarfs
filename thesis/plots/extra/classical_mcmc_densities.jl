### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ e8c8a05c-ccf8-416f-b16d-860d883c49e7
using CSV, DataFrames

# ╔═╡ cf655228-1528-4e77-9629-07e99984951f
using OrderedCollections

# ╔═╡ 6d2bff6d-b347-49c4-87df-ad58d8a27ff3
include("./paper_style.jl")

# ╔═╡ 67dd3488-79e1-4ba2-b9ac-f417765d55de
import TOML

# ╔═╡ 4ee54995-a31f-4ae1-8204-e55884d786bb
α = LilGuys.R_h(LilGuys.Exp2D(R_s=1))

# ╔═╡ a65be29f-bf1c-410d-9f62-cf2c2646f40a
function get_R_h(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors
	R_h = obs_props["R_h"]
	return middle(R_h)
end

# ╔═╡ 25b6d546-9edb-479f-bfd3-8a1097911694
function get_distance(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors
	d = obs_props["distance"]
	return middle(d)
end

# ╔═╡ b3c13c19-9022-4a7c-9355-104d930ae0e1
function get_R_h_fit(galaxyname)

	density_fit = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "mcmc/summary.mcmc_sersic.csv"), DataFrame)
	row = density_fit[density_fit.parameters .== ["R_h"], :] |> only

	return row.median
end

# ╔═╡ e826e3e3-dc01-4d56-9a35-0d69a2af3bbb
function get_fit_param(galaxyname, param)

	density_fit = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "mcmc/summary.mcmc_sersic.csv"), DataFrame)
	row = density_fit[density_fit.parameters .== [param], :] |> only

	return Measurement(row.median, row.lower_error, row.upper_error)
end

# ╔═╡ c425a56f-bbd1-4a30-a8a8-0d6e83cbfc99
get_R_h_fit("sculptor")

# ╔═╡ 40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
function load_profile(galaxyname; algname="../mcmc/hist_fast")
	@info galaxyname
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_profile.toml")
	
    prof = LilGuys.SurfaceDensityProfile(filename) 

	R_h = get_R_h_fit(galaxyname)
	@info "counts = $(sum(prof.counts))"

	Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(R_h))
	M_s = Σ_h * R_h .^ 2

	
    prof = LilGuys.scale(prof, 1/R_h, 1/M_s)

	filt = maximum.(LilGuys.error_interval.(prof.log_Sigma)) .< 1
	filt = LilGuys.find_longest_consecutive_true(filt)

	LilGuys.filter_by_bin(prof, filt)
end

# ╔═╡ 552e9438-862f-4710-a7d4-c8d798b5f1aa
galaxynames = [
	"fornax",
	"leo1",
	#"antlia2", # this one is unreliable...
	 "leo2",
	 "carina",
	 "draco",	
	"canes_venatici1",
	 "sextans1",
	"crater2",
]

# ╔═╡ e2b8987f-dc62-4cb1-9e96-1236e87bb096
n_values = get_fit_param.(galaxynames, "n")

# ╔═╡ cf3c2a09-9403-4ef8-964f-91957a57f3cf
let
	x = [get_fit_param(gal, "ellipticity") for gal in galaxynames]
	fig = errorscatter(x, n_values, yerror=error_interval.(n_values))


	n = [get_fit_param("sculptor", "n")]
	x = get_fit_param("sculptor", "ellipticity")
	errorscatter!([x], n, yerror=error_interval.(n))
	
	n = [get_fit_param("ursa_minor", "n")]
	x = get_fit_param("ursa_minor", "ellipticity")

	errorscatter!([x], n, yerror=error_interval.(n))

	hlines!(1)
	fig

end

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
		limits=(-1.5, 1.1, -1.2, 2)
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

	@savefig "classical_dwarf_mcmc_profiles"
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

# ╔═╡ 24ad252e-4823-407d-9f13-34a6c2dc4c73
md"""
# Some data for the table
"""

# ╔═╡ 5bd8d9b0-c662-464d-afbf-55509478f6d8
all_profiles = merge(profiles, Dict("sculptor"=>prof_scl, "ursa_minor"=>prof_umi))

# ╔═╡ 0e295946-eab3-45dd-8a09-e1528f638999
for (galaxy, prof) in all_profiles
	idx_max = length(prof.log_R)
	log_R_max = prof.log_R_bins[idx_max+1]
	log_Sigma_min = prof.log_Sigma[idx_max]
	@info galaxy, 10^log_R_max * get_R_h_fit(galaxy), log_Sigma_min, 10^prof.log_R[idx_max]

end

# ╔═╡ f1b597a7-6761-478a-9fb6-59c1a241e64d
"""
    latex_string(m::Measurement)

Return a LaTeX string like `1.23^{+0.04}_{-0.05}` with correct rounding for asymmetric uncertainties.
"""
function latex_string(m::Measurement)
    lo = m.lower
    hi = m.upper
    mid = m.middle

    # find order of magnitude of uncertainties
    sig_lo = floor(Int, log10(lo))
    sig_hi = floor(Int, log10(hi))
    sig = min(sig_lo, sig_hi)

    # we want 1 or 2 sig figs for uncertainty
    round_digits = -sig + 1

    lo_r = round(lo; digits=round_digits)
    hi_r = round(hi; digits=round_digits)
    mid_r = round(mid; digits=round_digits)

    return "\$$(mid_r)^{+$(hi_r)}_{-$(lo_r)}\$"
end


# ╔═╡ 56ce76cd-30e6-4a71-b4ea-c36bb0f25212
for (galaxy, prof) in all_profiles
	idx_max = length(prof.log_R)
	log_R_max = prof.log_R_bins[idx_max+1]
	log_Sigma_min = prof.log_Sigma[idx_max]
	log_distance = log10(get_distance(galaxy))
	println(galaxy, "\t", latex_string(log_Sigma_min .+ 2*log_distance))

end

# ╔═╡ c7999ba5-bb78-4bc4-8676-40aea5e7b50b
for (galaxy, prof) in all_profiles
	idx_max = length(prof.log_R)
	log_R_max = prof.log_R_bins[idx_max+1]
	log_Sigma_min = prof.log_Sigma[idx_max]
	println(galaxy, "\t", latex_string(get_fit_param(galaxy, "n")), "\t", latex_string(get_fit_param(galaxy, "R_h")), "\t", get_R_h(galaxy))

end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═67dd3488-79e1-4ba2-b9ac-f417765d55de
# ╠═e8c8a05c-ccf8-416f-b16d-860d883c49e7
# ╠═4ee54995-a31f-4ae1-8204-e55884d786bb
# ╠═a65be29f-bf1c-410d-9f62-cf2c2646f40a
# ╠═25b6d546-9edb-479f-bfd3-8a1097911694
# ╠═b3c13c19-9022-4a7c-9355-104d930ae0e1
# ╠═e826e3e3-dc01-4d56-9a35-0d69a2af3bbb
# ╠═e2b8987f-dc62-4cb1-9e96-1236e87bb096
# ╠═cf3c2a09-9403-4ef8-964f-91957a57f3cf
# ╠═c425a56f-bbd1-4a30-a8a8-0d6e83cbfc99
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
# ╠═24ad252e-4823-407d-9f13-34a6c2dc4c73
# ╠═5bd8d9b0-c662-464d-afbf-55509478f6d8
# ╠═0e295946-eab3-45dd-8a09-e1528f638999
# ╠═56ce76cd-30e6-4a71-b4ea-c36bb0f25212
# ╠═c7999ba5-bb78-4bc4-8676-40aea5e7b50b
# ╠═f1b597a7-6761-478a-9fb6-59c1a241e64d
