### A Pluto.jl notebook ###
# v0.20.20

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

# ╔═╡ af2976f5-b269-465b-9e1d-3a6f62230128
using CSV, DataFrames

# ╔═╡ f633c12c-5351-4169-911e-2afff2403fea
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 7564bff6-1ba9-4255-a94a-606450085eba
import TOML

# ╔═╡ f4b51722-b44b-44d7-a293-1e186a4c37cb
module Utils
	include("./model_utils.jl")
end

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ a15d7c08-4b28-4290-9a3c-96f0d35dd0ae
md"""
# Data Loading
"""

# ╔═╡ ea5a2bc1-9664-402c-a827-de3b96174605
function get_param(exp_fit, param)
	idx = only(findall(exp_fit.parameters .== [param]))
	return exp_fit.median[idx]
end

# ╔═╡ 04686bc6-3888-4df3-8f0c-502cc20ce86a
fit_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/mcmc/summary.mcmc_2exp.csv"), DataFrame)

# ╔═╡ 1ca18923-c3cd-4ea0-9e3d-adec98ccfe6e
fit_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/mcmc/summary.mcmc_2exp.csv"), DataFrame)

# ╔═╡ 4f256bd9-b96c-40c6-98b0-512ae2376ba3
multipop_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", "multipop", "processed/sculptor.mcmc_2pop_vel_fe.summary.csv"), DataFrame)

# ╔═╡ 661caeeb-75ac-4af2-af29-daba92f30970


# ╔═╡ 21688d34-0e1c-4d37-a182-ca100a0623e8
multipop_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations", "multipop", "processed/ursa_minor.mcmc_2pop_vel_fe.summary.csv"), DataFrame)

# ╔═╡ fdb0e872-0ba0-4471-9306-8e9b8760bc29
rv_scl = read_fits(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/velocities/processed/rv_combined_x_wide_2c_psat_0.2.fits"))

# ╔═╡ 1a95a435-acdc-45d9-8eb4-5d8edaaa98cf
rv_umi = read_fits(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/velocities/processed/rv_combined_x_2c_psat_0.2.fits"))

# ╔═╡ 8366afbb-c528-40cc-8035-18927f77ca3b
metals_scl = let
	df = copy(rv_scl)

	df[!, "fe_h"] = df.fe_h_t23
	df[!, "fe_h_err"] = df.fe_h_err_t23
	df[!, "study"] .= "T+23"

	filt = ismissing.(df.fe_h)
	df[filt, "fe_h"] .= df.fe_h_apogee[filt]
	df[filt, "fe_h_err"] = df.fe_h_err_apogee[filt]
	df[filt, "study"] .= "apogee"

	filt = ismissing.(df.fe_h)
	df[filt, "fe_h"] .= df.fe_h_gmos[filt]
	df[filt, "fe_h_err"] = df.fe_h_err_gmos[filt]
	df[filt, "study"] .= "S+23"

	df[.!ismissing.(df.fe_h), :]
end

# ╔═╡ 689bb4f2-c1f7-4ec8-a11f-4b4078160bc5
metals_umi = let
	df = copy(rv_umi)

	df[!, "fe_h"] = df.fe_h_p20
	df[!, "fe_h_err"] = df.fe_h_err_p20
	df[!, "study"] .= "p+20"

	filt = ismissing.(df.fe_h)
	df[filt, "fe_h"] .= df.fe_h_apogee[filt]
	df[filt, "fe_h_err"] = df.fe_h_err_apogee[filt]
	df[filt, "study"] .= "apogee"

	filt = ismissing.(df.fe_h)
	df[filt, "fe_h"] .= df.fe_h_graces[filt]
	df[filt, "fe_h_err"] .= df.fe_h_err_graces[filt]
	df[filt, "study"] .= "S+23"
	
	df[.!ismissing.(df.fe_h), :]
end

# ╔═╡ 9b6781e0-7490-45dd-b807-4f49309eb8ca
md"""
# Plots
"""

# ╔═╡ ca1d0e3b-5b47-44a4-9ee0-f7583615908c
import StatsBase: median, mean, sem

# ╔═╡ fc5aaacc-72b1-41b1-9816-bd0ac47fff55
function median_plot!(x, y, y_err; num_bins=20, kwargs...)
	num_per_bin = length(x) ÷ num_bins
	bins = LilGuys.bins_equal_number(x, nothing,  num_per_bin=num_per_bin)
	Nb = length(bins) - 1

	x_m = midpoints(bins)
	y_m = zeros(Nb)
	y_m_err = zeros(Nb)
	for i in 1:Nb
		filt = bins[i] .<= x .< bins[i+1]
		y_m[i] = mean(y[filt])
		y_m_err[i] = sem(y[filt]) + sqrt(1 / sum(y_err[filt] .^-2))
	end

	errorscatter!(x_m, y_m, yerror=y_m_err; color=:black, markersize=4, kwargs...)
end

# ╔═╡ 74705326-077f-47c8-82ae-2d3c23fe853b
function get_profiles(df, dist=nothing)
	R_1 = get_param(df, "R_s")
	f_outer = get_param(df, "f_outer")
	R_2 = get_param(df, "R_s_outer")

	if !isnothing(dist)
		R_1 = LilGuys.arcmin2kpc(R_1, dist)
		R_2 = LilGuys.arcmin2kpc(R_2, dist)
	end
	
	prof = LilGuys.Exp2D(R_s=R_1, M=(1-f_outer))
	prof_outer = LilGuys.Exp2D(R_s=R_2, M=f_outer)

	return prof, prof_outer
end

# ╔═╡ 8bbf2aa6-0804-4f18-b5e6-4b919635663e
function plot_median_profile(gs, df, prof_obs)
	ax = Axis(gs, xlabel = "log R / arcmin", ylabel = "log surface density")
	
	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma), color=:black, markersize=4)

	
	x = LinRange(-1, 2.5, 1000)

	R = 10 .^ x
	
	prof, prof_outer = get_profiles(df)
	prof_double = LilGuys.DoubleExp2D(prof.M, prof.R_s, prof_outer.M, prof_outer.R_s)


	scale = sum(prof_obs.counts)
	y = @. log10(LilGuys.surface_density(prof, R) * scale)
	y2 = @. log10(LilGuys.surface_density(prof_outer, R) * scale)
	y3 = @. log10(LilGuys.surface_density(prof_double, R) * scale)
	lw = theme(:linewidth)[]/2
	lines!(x, y, color=COLORS[1], linewidth=lw)
	lines!(x, y2, color=COLORS[4], linewidth=lw)
	lines!(x, y3, color=COLORS[5], linewidth=lw)

	
	ylims!(-6, 3)

	ax

end

# ╔═╡ 4f9e9a52-9846-4307-aa7b-9cfc9583737f
function plot_mixture!(df, y_low, y_high)
	
	x = LinRange(-1, 2.5, 1000)
	R = 10 .^ x
	prof, prof_outer = get_profiles(df)
	
	y1 = @. (LilGuys.surface_density(prof, R))
	y2 = @. (LilGuys.surface_density(prof_outer, R))

	y = (y1 .* y_high .+ y2 * y_low) ./ (y1 .+ y2)
	lines!(x, y, color=COLORS[5], linewidth=1)

	ylims!(-6, 3)

end

# ╔═╡ 23f21352-f175-414c-8f21-e6a27f73d634
styles = Dict(
	"S+23" => (;markersize=5, color=COLORS[2], marker=:star5),
	"T+23" => (;),
	"apogee" => (;color=COLORS[3], marker=:rect),
	"p+20" => (;)
)

# ╔═╡ f1bc7e02-0ee8-46e1-bd7a-0d833ab306b4
df_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/velocities/processed/vz_r_ell_binned.rv_combined_x_wide_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ e4076fb8-6d7a-4bfb-883d-25f452d935e0
df_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/velocities/processed/vz_r_ell_binned.rv_combined_x_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ 23bcdc62-1c21-47b3-b8a3-f616ed779bd3
function plot_σv_obs!(galaxyname)
	if galaxyname == "sculptor"
		df = df_scl
	elseif galaxyname == "ursa_minor"
		df = df_umi
	end

	errorscatter!(log10.(df.x), (df.σ), yerror=df.σ_em, color=:black, markersize=6)
end

# ╔═╡ bd656734-acb8-4653-9fb8-684272bc2ce5
get_obs_props(galaxyname) = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))

# ╔═╡ a7441217-1826-4d15-b579-766493009c4f
rmax_0, vmax_0 = LilGuys.r_circ_max(LilGuys.ExpCusp(M=1, r_s=1)), LilGuys.v_circ_max(LilGuys.ExpCusp(M=1, r_s=1)) * V2KMS

# ╔═╡ ba50a481-42ff-4ddf-b635-69452b180e7b
halo_scl = NFW(r_circ_max = 2.5, v_circ_max=25/V2KMS)

# ╔═╡ 31a9c31d-a466-4c61-aecd-5d673424fe54
halo_umi = NFW(r_circ_max=5, v_circ_max=27/V2KMS) #LilGuys.ExpCusp(r_s = 1 / rmax_0, M=(19.4/vmax_0)^2 * (1 / rmax_0))

# ╔═╡ 6ab650e4-174e-4f2e-af27-af922056d68e
LilGuys.r_circ_max(halo_umi), LilGuys.v_circ_max(halo_umi) * V2KMS

# ╔═╡ 054cc22b-4595-47fb-b4b7-e04638bf77ed
function rho_sigma2_r(halo, prof, r)
	integrand(r) = LilGuys.density(prof, r) * LilGuys.mass(halo, r) / r^2

	return LilGuys.integrate(integrand, r, Inf)
end

# ╔═╡ d7807df6-55a2-43dd-821d-1e476be084b3
function lerp_rho_sigma2(halo, prof, rs)
	ys = rho_sigma2_r.(halo, prof, rs)
	return LilGuys.lerp(rs, ys)
end

# ╔═╡ caf4a15e-6005-4d28-b76a-b0c519ec6f80
function sigma_los(halo::LilGuys.SphericalProfile, prof, R)
	# f = lerp_rho_sigma2(halo, prof, logrange(0.0001, 1000, 3000))
	integrand(r) = rho_sigma2_r(halo, prof, r) * r / sqrt(r^2 - R^2)

	Sigma_sigma2 = 2*LilGuys.integrate(integrand, R* (1+1e-10), Inf)

	return sqrt(Sigma_sigma2 / LilGuys.surface_density(prof, R))
end

# ╔═╡ 8a757f44-0ef3-4ed2-bf90-f9edbda8d102
function plot_mixture_sigma!(df, halo, galaxyname; plot_components=false, color=COLORS[5])
	
	x = LinRange(-0.5, 2.5, 100)
	R = 10 .^ x
	
	dist = get_obs_props(galaxyname)["distance"]

	R_kpc = LilGuys.arcmin2kpc.(R, dist)
	x_kpc = log10.(R_kpc)
	prof, prof_outer = get_profiles(df, dist)
	prof_double = LilGuys.DoubleExp2D(prof.M, prof.R_s, prof_outer.M, prof_outer.R_s)


	σ = sigma_los.(halo, prof_double, R_kpc)


	lines!(x, σ * V2KMS, color=color, linewidth=1)
	
	if plot_components
		σ1 = sigma_los.(halo, prof, R_kpc)
		σ2 = sigma_los.(halo, prof_outer, R_kpc)

		lines!(x, σ1*V2KMS, color=COLORS[1], linewidth=1)
		lines!(x, σ2*V2KMS, color=COLORS[4], linewidth=1)
	end

	ylims!(-6, 3)

end

# ╔═╡ 0b722f4c-368c-4308-abcb-2584252f595f
sigma_los(halo_scl, get_profiles(fit_scl, get_obs_props("sculptor")["distance"])[1], 0.2) * V2KMS

# ╔═╡ d8bb89ad-c76c-4a28-b818-b12b782b9966
scatter(log10.(df_umi.x), df_umi.σ)

# ╔═╡ c2b44dcf-cc5e-4630-b4c7-46fcb1725b66
let
	fig = Figure(size=(3.5, 4) .* 72)

	ax_scl_dens = plot_median_profile(fig[1, 1], fit_scl,Utils.load_expected_density_profile("sculptor") )

	ylims!(-3, 2)
	ax_scl_dens.title = "Sculptor"
	
	ax_scl = Axis(fig[2,1], ylabel="[Fe/H]", xlabel=L"$\log\, R_\textrm{ell}$ / arcmin")
	for study in ["T+23", "apogee", "S+23"]
		df = metals_scl[metals_scl.study .== study, :]
		label_key = Dict("T+23" => "T+23/P+20", "apogee"=>"APOGEE", "S+23"=>"S+23a/b")
		scatter!(log10.(df.R_ell), df.fe_h, markersize=2; label=label_key[study], styles[study]...)
	end

	lines!([NaN], [NaN], color=COLORS[1], label="inner model", linewidth=1)
	lines!([NaN], [NaN], color=COLORS[4], label="outer",  linewidth=1)
	lines!([NaN], [NaN], color=COLORS[5], label="combined",  linewidth=1)
	
	median_plot!(log10.(metals_scl.R_ell), metals_scl.fe_h, metals_scl.fe_h_err)
	# vlines!(r_trans_scl, color=COLORS[2], linewidth=theme(:linewidth)[]/2)
	plot_mixture!(fit_scl, get_param(multipop_scl, "mu_fe_b"), get_param(multipop_scl, "mu_fe_a"))

	# axislegend(position=:rt, patchsize=(5, 5), backgroundcolor=(:white, 0.8))
	ylims!(-4, -0.5)
	xlims!(-0.2, 2.2)


	ax_umi_dens = plot_median_profile(fig[1, 2], fit_umi,Utils.load_expected_density_profile("ursa_minor") )
	ax_umi_dens.title = "Ursa Minor"
	ax_umi_dens.ylabel = ""
	ax_umi_dens.yticklabelsvisible = false
	ylims!(-3.5, 2)
	# vlines!(r_trans_umi, color=COLORS[2], linewidth=theme(:linewidth)[]/2)

	
	ax_umi = Axis(fig[2,2], ylabel="[Fe/H]", xlabel=L"$\log\, R_\textrm{ell}$ / arcmin")
	for study in ["p+20", "apogee", "S+23"]
		df = metals_umi[metals_umi.study .== study, :]
		scatter!(log10.(df.R_ell), df.fe_h, markersize=2; label=study, styles[study]...)
	end
	ax_umi.ylabel = ""
	ax_umi.yticklabelsvisible = false
	# axislegend(position=:lt)

	median_plot!(log10.(metals_umi.R_ell), metals_umi.fe_h, metals_umi.fe_h_err)
	# vlines!(r_trans_umi, color=COLORS[2], linewidth=theme(:linewidth)[]/2)
	plot_mixture!(fit_umi, get_param(multipop_umi, "mu_fe_b"), get_param(multipop_umi, "mu_fe_a"))

	ylims!(-3.5, -0.5)
	xlims!(-0.5, 2.2)

	ax_scl_sigma = Axis(fig[3,1], xlabel=L"$\log\, R_\textrm{ell}$ / arcmin", ylabel=L"$\sigma_\textrm{v, los}$ / km s$^{-1}$")

	plot_σv_obs!("sculptor")
	plot_mixture_sigma!(fit_scl, halo_scl, "sculptor", plot_components=true)

	ylims!(4, 15)
	xlims!(-0.5, 2.2)



	ax_umi_sigma = Axis(fig[3,2] ,xlabel=L"$\log\, R_\textrm{ell}$ / arcmin")
	ax_umi_sigma.ylabel = ""
	ax_umi_sigma.yticklabelsvisible = false


	plot_σv_obs!("ursa_minor")
	plot_mixture_sigma!(fit_umi, halo_umi, "ursa_minor", plot_components=true)
	ylims!(4, 15)
	xlims!(-0.5, 2.2)

	linkxaxes!(ax_scl, ax_scl_dens)
	linkxaxes!(ax_umi, ax_umi_dens)

	linkyaxes!(ax_scl, ax_umi)
	linkyaxes!(ax_scl_dens, ax_umi_dens)

	# hidexdecorations!(ax_scl, ticks=false, minorticks=false)
	hidexdecorations!(ax_scl_dens, ticks=false, minorticks=false)
	hidexdecorations!(ax_scl, ticks=false, minorticks=false)

	hidexdecorations!(ax_umi_dens, ticks=false, minorticks=false)
	hidexdecorations!(ax_umi, ticks=false, minorticks=false)

	rowgap!(fig.layout, 5)
	colgap!(fig.layout, 5)

	# rowsize!(fig.layout, 1, Aspect(1, 1))
	# rowsize!(fig.layout, 2, Aspect(1, 1))


	Legend(fig[4, :], ax_scl, nbanks=3, tellheight=true)
	resize_to_layout!()


	@savefig "scl_umi_fe_h_gradient"
	fig
end

# ╔═╡ 22598e4a-b7bb-4ed3-a71b-e3e29b37e58a
md"""
# Naive transition radius fits
"""

# ╔═╡ b70ef4cc-4ea4-46cf-b701-2d6a3e92374c


# ╔═╡ 5838512b-a957-444c-8e15-0ce6304daeb5
 get_param(multipop_scl, "mu_fe_b"),  get_param(multipop_scl, "sigma_fe_b")

# ╔═╡ e491768b-8092-4d0f-b8d2-af9afe0eb108
 get_param(multipop_scl, "mu_fe_a"),  get_param(multipop_scl, "sigma_fe_a")

# ╔═╡ 18a4d0e6-f56b-4b4f-923e-169161e7ccb2


# ╔═╡ 4e5653d5-5723-4763-ab4f-18849dcabd85
get_param(multipop_umi, "mu_fe_b"),  get_param(multipop_umi, "sigma_fe_b")

# ╔═╡ 3f077048-62d6-4355-9adf-50b6e09140ad
get_param(multipop_umi, "mu_fe_a"),  get_param(multipop_umi, "sigma_fe_a")

# ╔═╡ 5e45fcdb-1812-4853-89c7-7f8a3ec3e478
function weighted_mean(x, x_err)
	w = @.  1/x_err^2
	return LilGuys.mean(x, w), LilGuys.std(x, w), sqrt(1 / sum(w))
end

# ╔═╡ 60e28a3c-afc8-4e46-82f1-8b424ca836cb
function weighted_mean(filt, x, x_err)
	filt = filt .& isfinite.(1 ./ x_err)
	return weighted_mean(x[filt], x_err[filt])
end

# ╔═╡ e5e5e8e7-0ff3-4bee-b6d3-4f425dc23c02
function get_r_trans(df)
	N = size(df, 1)
	r_trans = zeros(N)

	prof, prof_outer = get_profiles(df)
	r_trans = LilGuys.find_zero(r -> LilGuys.surface_density(prof, r) - LilGuys.surface_density(prof_outer, r), 10)
	
	r_trans
end
		

# ╔═╡ cfa85373-0eed-482b-a68a-1f3516145fb5
r_trans_umi = (get_r_trans(fit_umi))

# ╔═╡ e00cb812-94da-41b2-bedd-38c88ba2c92c
weighted_mean(metals_umi.R_ell .> r_trans_umi, disallowmissing(metals_umi.fe_h), metals_umi.fe_h_err)

# ╔═╡ 4138393e-e923-4f2f-a569-cf9e5bfb3a5e
weighted_mean(metals_umi.R_ell .< r_trans_umi, disallowmissing(metals_umi.fe_h), metals_umi.fe_h_err)

# ╔═╡ e6eb753c-7989-4bf1-bc90-6027096762fe
r_trans_scl = get_r_trans(fit_scl)

# ╔═╡ e5c280b0-bba3-4f38-9781-00daca41d508
weighted_mean(metals_scl.R_ell .> r_trans_scl, disallowmissing(metals_scl.fe_h), metals_scl.fe_h_err)

# ╔═╡ 9ccb9aa0-9f0b-4dd4-9d6c-3a6174a2bab6
weighted_mean(metals_scl.R_ell .< r_trans_scl, disallowmissing(metals_scl.fe_h), metals_scl.fe_h_err)

# ╔═╡ 39866a0d-ff2d-4963-b61b-93bd3b6a27c7
LilGuys.mean(metals_scl.fe_h[metals_scl.R_ell .< r_trans_scl], 1 ./ metals_scl.fe_h_err[metals_scl.R_ell .< r_trans_scl] .^ 2), get_param(multipop_scl, "mu_fe_a")

# ╔═╡ 33b41ea2-b1b6-4edb-80fa-73583db11d7e
LilGuys.mean(metals_scl.fe_h[metals_scl.R_ell .< r_trans_scl], 1 ./ metals_scl.fe_h_err[metals_scl.R_ell .< r_trans_scl] .^ 2), get_param(multipop_scl, "mu_fe_a")

# ╔═╡ d3a907cc-241e-49d8-8d91-2fd9061e520b
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═7564bff6-1ba9-4255-a94a-606450085eba
# ╠═af2976f5-b269-465b-9e1d-3a6f62230128
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═f4b51722-b44b-44d7-a293-1e186a4c37cb
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f633c12c-5351-4169-911e-2afff2403fea
# ╟─a15d7c08-4b28-4290-9a3c-96f0d35dd0ae
# ╠═ea5a2bc1-9664-402c-a827-de3b96174605
# ╠═04686bc6-3888-4df3-8f0c-502cc20ce86a
# ╠═1ca18923-c3cd-4ea0-9e3d-adec98ccfe6e
# ╠═4f256bd9-b96c-40c6-98b0-512ae2376ba3
# ╠═661caeeb-75ac-4af2-af29-daba92f30970
# ╠═21688d34-0e1c-4d37-a182-ca100a0623e8
# ╠═fdb0e872-0ba0-4471-9306-8e9b8760bc29
# ╠═1a95a435-acdc-45d9-8eb4-5d8edaaa98cf
# ╠═8366afbb-c528-40cc-8035-18927f77ca3b
# ╠═689bb4f2-c1f7-4ec8-a11f-4b4078160bc5
# ╠═9b6781e0-7490-45dd-b807-4f49309eb8ca
# ╠═ca1d0e3b-5b47-44a4-9ee0-f7583615908c
# ╠═fc5aaacc-72b1-41b1-9816-bd0ac47fff55
# ╠═8bbf2aa6-0804-4f18-b5e6-4b919635663e
# ╠═74705326-077f-47c8-82ae-2d3c23fe853b
# ╠═4f9e9a52-9846-4307-aa7b-9cfc9583737f
# ╠═23f21352-f175-414c-8f21-e6a27f73d634
# ╠═f1bc7e02-0ee8-46e1-bd7a-0d833ab306b4
# ╠═e4076fb8-6d7a-4bfb-883d-25f452d935e0
# ╠═23bcdc62-1c21-47b3-b8a3-f616ed779bd3
# ╠═bd656734-acb8-4653-9fb8-684272bc2ce5
# ╠═a7441217-1826-4d15-b579-766493009c4f
# ╠═ba50a481-42ff-4ddf-b635-69452b180e7b
# ╠═31a9c31d-a466-4c61-aecd-5d673424fe54
# ╠═6ab650e4-174e-4f2e-af27-af922056d68e
# ╠═054cc22b-4595-47fb-b4b7-e04638bf77ed
# ╠═d7807df6-55a2-43dd-821d-1e476be084b3
# ╠═caf4a15e-6005-4d28-b76a-b0c519ec6f80
# ╠═8a757f44-0ef3-4ed2-bf90-f9edbda8d102
# ╠═0b722f4c-368c-4308-abcb-2584252f595f
# ╠═d8bb89ad-c76c-4a28-b818-b12b782b9966
# ╠═c2b44dcf-cc5e-4630-b4c7-46fcb1725b66
# ╟─22598e4a-b7bb-4ed3-a71b-e3e29b37e58a
# ╠═b70ef4cc-4ea4-46cf-b701-2d6a3e92374c
# ╠═cfa85373-0eed-482b-a68a-1f3516145fb5
# ╠═e6eb753c-7989-4bf1-bc90-6027096762fe
# ╠═e5c280b0-bba3-4f38-9781-00daca41d508
# ╠═5838512b-a957-444c-8e15-0ce6304daeb5
# ╠═9ccb9aa0-9f0b-4dd4-9d6c-3a6174a2bab6
# ╠═e491768b-8092-4d0f-b8d2-af9afe0eb108
# ╠═18a4d0e6-f56b-4b4f-923e-169161e7ccb2
# ╠═e00cb812-94da-41b2-bedd-38c88ba2c92c
# ╠═4e5653d5-5723-4763-ab4f-18849dcabd85
# ╠═4138393e-e923-4f2f-a569-cf9e5bfb3a5e
# ╠═3f077048-62d6-4355-9adf-50b6e09140ad
# ╠═39866a0d-ff2d-4963-b61b-93bd3b6a27c7
# ╠═5e45fcdb-1812-4853-89c7-7f8a3ec3e478
# ╠═60e28a3c-afc8-4e46-82f1-8b424ca836cb
# ╠═33b41ea2-b1b6-4edb-80fa-73583db11d7e
# ╠═e5e5e8e7-0ff3-4bee-b6d3-4f425dc23c02
# ╠═d3a907cc-241e-49d8-8d91-2fd9061e520b
