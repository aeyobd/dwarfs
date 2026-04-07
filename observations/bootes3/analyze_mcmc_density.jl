### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 61558ffd-4ccc-47b6-89e8-20861e7827d7
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using PyFITS
	import PairPlots
	
	using Turing
end

# ╔═╡ a060d1e1-1e5a-49d6-ac68-2a6bd1d1a455
using CSV, DataFrames

# ╔═╡ 8442396b-50ac-4614-becb-c7eb8bf041ff
samplename = "delve_matched_filter_6deg"

# ╔═╡ 2f85f664-7498-4960-bf8a-163af435a192
begin 
	using LilGuys

	FIGDIR = joinpath("mcmc", "figures"); FIGSUFFIX=".$samplename.mcmc"
end

# ╔═╡ 22574221-b8ff-4e76-ac1a-f0c9b316499b
import NaNMath

# ╔═╡ 813eca57-4dc6-4999-892b-cf910b02127d
module GaiaFilters
	include("../../utils/gaia_filters.jl")
end

# ╔═╡ 2a6238fa-bc93-46a1-8091-0d19b1b0266b


# ╔═╡ c2835710-a51a-410e-8cfb-192e3e879236
md"""
## Plots
"""

# ╔═╡ 322b527d-c444-48fc-bf89-396fd0f9dfce
log_r_label = "log r / arcmin"

# ╔═╡ 633e6c10-3609-4979-824f-d82436937677
log_Sigma_label = L"log \Sigma"

# ╔═╡ d47c1295-84b5-4673-b367-041dbbbffd0e
allstars = read_fits("samples/$samplename.fits")

# ╔═╡ 2ccd02ea-dc0e-4605-9d50-a305dd7fe06f
stars = allstars[allstars.L_sat ./ allstars.L_bg .> 0, :]

# ╔═╡ 90902e53-7d7c-41f2-9d71-158b1cf80e9d
df_plummer = CSV.read("mcmc/samples.$samplename.mcmc_plummer.csv", DataFrame)

# ╔═╡ b76fcd10-3cbb-4d58-b66f-d942eb02f8ea
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	scatter!(stars.xi, stars.eta, markersize=0.3, color=:black, 
	   )

	scatter!(0, 0)
	# scatter!(-median(df_out.d_xi), -median(df_out.d_eta))
	scatter!(median(df_plummer.d_xi), median(df_plummer.d_eta))
	# scatter!(-median(df_ell.d_xi), -median(df_ell.d_eta))
	
	fig
end

# ╔═╡ 199f497a-40b0-4fc1-9e40-397307b11965
Measurement(median(df_plummer.f_sat))

# ╔═╡ 3902daf7-c819-4159-8b89-17abea4a7fcf
R_max = GaiaFilters.calc_R_max(allstars.xi, allstars.eta, 0, 0)

# ╔═╡ 7ace8961-8ed4-4387-8518-07f4643da622
N_stars = median(df_plummer.N_sat ./ df_plummer.f_sat)

# ╔═╡ 1b470198-af01-4118-9d81-bcb4febcbee6
N_stars / (π*R_max^2)

# ╔═╡ 0427aed9-eb83-44f8-8469-0c00d5532cd6
LR = (stars.L_sat ./ stars.L_bg) ./ LilGuys.mean((stars.L_sat ./ stars.L_bg))

# ╔═╡ ebf4ddbd-c1ed-4073-bcc5-3e7fd12c4c99
sum(LR)

# ╔═╡ c2942ebf-157f-45fc-9dcd-911b73725b16
length(stars.R)

# ╔═╡ 050a5206-f0d9-4811-9248-91f5e937bfaf
log10(maximum(stars.R))

# ╔═╡ 88980f89-b8f9-4e72-8b9e-7f5d668a2776
# ╠═╡ disabled = true
#=╠═╡
prof_obs = LilGuys.SurfaceDensityProfile(stars.R, bins=(0:0.1:2.541), weights= LR) |> LilGuys.filter_empty_bins
  ╠═╡ =#

# ╔═╡ 9421bde9-41b4-4e67-8512-f094a75e5e1b
prof_obs = LilGuys.SurfaceDensityProfile(stars.R, bins=(0:0.1:2.541), weights= LR) |> LilGuys.filter_empty_bins

# ╔═╡ cd9df222-ac47-4dce-b74d-77d11486f31c
function plot_sersic!(df_out)


	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Sersic(R_h=df_out.R_h[i], n=df_out.n[i])
		M = N_stars * df_out.f_sat[i] / LilGuys.mass_2D(prof, R_max)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(x, y, color=COLORS[2], alpha=0.1)
	end

	ylims!(-6, 3)



end

# ╔═╡ a281ea08-f49f-4fa4-b22a-3b830d29981d
function plot_plummer_bkg!(df_out, N_stars; color=COLORS[2])

	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	for i in 1:100
		prof =  LilGuys.Plummer(r_s=df_out.R_h[i], M=1)
		M = N_stars * df_out.f_sat[i]
		Σ_bkg = N_stars * (1- df_out.f_sat[i]) / (π*R_max^2)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M .+ Σ_bkg)

		lines!(x, y, color=color, alpha=0.1)
	end

end

# ╔═╡ b3505169-98cf-4b25-a8da-ea36bbcc4a23
function plot_plummer!(df_out; color=COLORS[2])

	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Plummer(r_s=df_out.R_h[i], M=1)
		M = N_stars * df_out.f_sat[i]
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(x, y, color=color, alpha=0.1)
	end

end

# ╔═╡ 35a82d05-8207-475d-bd4f-ad2266735df8
log_Σ_bg = median(prof_obs.log_Sigma[prof_obs.log_R .> 1.8])

# ╔═╡ c6f80088-c6e5-43ed-95db-6d93873aebf8
Sigma_sub = 10 .^ prof_obs.log_Sigma .- 10 .^ log_Σ_bg

# ╔═╡ f648277f-8eb5-4527-8cfb-8e1a99fe16f9
filt_bkg = NaNMath.log10.(Sigma_sub .- LilGuys.lower_error.(Sigma_sub)) .|> isfinite |> LilGuys.find_longest_consecutive_true

# ╔═╡ 389dc523-5188-4548-9881-e17d761fbb59
begin
	log_Sigma_sub = log10.(Sigma_sub[filt_bkg])
	prof_sub = LilGuys.SurfaceDensityProfile(
		R_units=prof_obs.R_units,
		log_R=prof_obs.log_R[filt_bkg],
		log_R_bins=prof_obs.log_R_bins[filt_bkg[1]:filt_bkg[end]+1],
		log_Sigma=log_Sigma_sub,
		log_R_scale=prof_obs.log_R_scale,
		log_m_scale=prof_obs.log_m_scale,
		annotations=prof_obs.annotations
	)
end

# ╔═╡ ddba7645-45d8-4fc1-9cd9-d33c486910ee
N_prof = sum(π * diff(LilGuys.radii_bins(prof_sub) .^ 2) .* 10 .^ prof_sub.log_Sigma)

# ╔═╡ 500057bf-c451-414b-a9bb-4b0138c6486b
prof_shift = log10.((median(df_plummer.f_sat) .* N_stars ./ N_prof ))

# ╔═╡ ae519c10-c59e-4011-a4e0-798d1346137f
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=log_r_label, ylabel=log_Sigma_label)




	plot_plummer_bkg!(df_plummer, N_stars * 10^(2prof_shift.middle); color=COLORS[2])

	# plot_sersic!(df_sersic)
	plot_plummer!(df_plummer, color=COLORS[3])
	# plot_plummer!(df_ell, color=COLORS[4])
	ylims!(-4, -0.5)
	hlines!(middle.((log_Σ_bg)) .+ prof_shift.middle)

	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma .+ prof_shift.middle), yerror=error_interval.(prof_obs.log_Sigma), label="observed", color=COLORS[1])

	errorscatter!(prof_sub.log_R, (prof_sub.log_Sigma .+ prof_shift.middle), yerror=error_interval.(prof_sub.log_Sigma), color=:black, label="bg-subtracted")

	lines!([NaN], [NaN], color=COLORS[2], label="Sersic")
	lines!([NaN], [NaN], color=COLORS[3], label="Plummer")
	ylims!(-3, -0.5)
	axislegend(position=:lb)


	fig
end

# ╔═╡ 2cf3dc2c-14ec-4483-884e-2d859453ce50
hist(LR)

# ╔═╡ b80ead46-7202-495e-865b-1ed6b6a00844
# prof_obs = LilGuys.SurfaceDensityProfile("density_profiles/delve_mf_sub_profile.toml") |> LilGuys.filter_empty_bins

# ╔═╡ 64b08c1d-8fc8-44ce-8e2c-0962180f46d0
function plot_plummer(df_out)
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=log_r_label, ylabel=log_Sigma_label)


	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))

	plot_plummer!(df_out)

	ylims!(-6, 3)

	fig

end

# ╔═╡ Cell order:
# ╠═8442396b-50ac-4614-becb-c7eb8bf041ff
# ╠═61558ffd-4ccc-47b6-89e8-20861e7827d7
# ╠═22574221-b8ff-4e76-ac1a-f0c9b316499b
# ╠═a060d1e1-1e5a-49d6-ac68-2a6bd1d1a455
# ╠═2f85f664-7498-4960-bf8a-163af435a192
# ╠═813eca57-4dc6-4999-892b-cf910b02127d
# ╠═2a6238fa-bc93-46a1-8091-0d19b1b0266b
# ╟─c2835710-a51a-410e-8cfb-192e3e879236
# ╠═322b527d-c444-48fc-bf89-396fd0f9dfce
# ╠═633e6c10-3609-4979-824f-d82436937677
# ╠═d47c1295-84b5-4673-b367-041dbbbffd0e
# ╠═2ccd02ea-dc0e-4605-9d50-a305dd7fe06f
# ╠═b76fcd10-3cbb-4d58-b66f-d942eb02f8ea
# ╠═90902e53-7d7c-41f2-9d71-158b1cf80e9d
# ╠═199f497a-40b0-4fc1-9e40-397307b11965
# ╠═ddba7645-45d8-4fc1-9cd9-d33c486910ee
# ╠═500057bf-c451-414b-a9bb-4b0138c6486b
# ╠═3902daf7-c819-4159-8b89-17abea4a7fcf
# ╠═7ace8961-8ed4-4387-8518-07f4643da622
# ╠═ae519c10-c59e-4011-a4e0-798d1346137f
# ╠═cd9df222-ac47-4dce-b74d-77d11486f31c
# ╠═1b470198-af01-4118-9d81-bcb4febcbee6
# ╠═a281ea08-f49f-4fa4-b22a-3b830d29981d
# ╠═b3505169-98cf-4b25-a8da-ea36bbcc4a23
# ╠═35a82d05-8207-475d-bd4f-ad2266735df8
# ╠═c6f80088-c6e5-43ed-95db-6d93873aebf8
# ╠═f648277f-8eb5-4527-8cfb-8e1a99fe16f9
# ╠═389dc523-5188-4548-9881-e17d761fbb59
# ╠═0427aed9-eb83-44f8-8469-0c00d5532cd6
# ╠═ebf4ddbd-c1ed-4073-bcc5-3e7fd12c4c99
# ╠═c2942ebf-157f-45fc-9dcd-911b73725b16
# ╠═050a5206-f0d9-4811-9248-91f5e937bfaf
# ╠═88980f89-b8f9-4e72-8b9e-7f5d668a2776
# ╠═9421bde9-41b4-4e67-8512-f094a75e5e1b
# ╠═2cf3dc2c-14ec-4483-884e-2d859453ce50
# ╠═b80ead46-7202-495e-865b-1ed6b6a00844
# ╠═64b08c1d-8fc8-44ce-8e2c-0962180f46d0
