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

# ╔═╡ af2976f5-b269-465b-9e1d-3a6f62230128
using CSV, DataFrames

# ╔═╡ f633c12c-5351-4169-911e-2afff2403fea
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ f4b51722-b44b-44d7-a293-1e186a4c37cb
module Utils
	include("./model_utils.jl")
end

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 04686bc6-3888-4df3-8f0c-502cc20ce86a
samples_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/mcmc/samples.mcmc_2exp.csv"), DataFrame)

# ╔═╡ 0bd742a6-e683-4337-b195-e24d0e9db65c


# ╔═╡ 1ca18923-c3cd-4ea0-9e3d-adec98ccfe6e
samples_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/mcmc/samples.mcmc_2exp.csv"), DataFrame)

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

	errorscatter!(x_m, y_m, yerror=y_m_err; color=:black, markersize=6, kwargs...)
end

# ╔═╡ e5e5e8e7-0ff3-4bee-b6d3-4f425dc23c02
function get_r_trans(df)
	N = size(df, 1)
	r_trans = zeros(N)
	for i in 1:N
		prof = LilGuys.Exp2D(R_s=df.R_s[i], M=(1-df.f_outer[i]))
		prof_outer = LilGuys.Exp2D(R_s=df.R_s_outer[i], M=df.f_outer[i])

		r_trans[i] = LilGuys.find_zero(r -> LilGuys.surface_density(prof, r) - LilGuys.surface_density(prof_outer, r), 10)
	end
	r_trans
end
		

# ╔═╡ e6eb753c-7989-4bf1-bc90-6027096762fe
r_trans_scl = median(log10.(get_r_trans(samples_scl)))

# ╔═╡ cfa85373-0eed-482b-a68a-1f3516145fb5
r_trans_umi = median(log10.(get_r_trans(samples_umi)))

# ╔═╡ 3288393c-a1b0-4772-815f-965be98d5a29
10^r_trans_scl

# ╔═╡ 52e12e1b-740d-4ae4-886f-1f3f849b5d7e
10^r_trans_umi

# ╔═╡ 8bbf2aa6-0804-4f18-b5e6-4b919635663e
function plot_samples(gs, df, prof_obs)
	ax = Axis(gs, xlabel = "log R / arcmin", ylabel = "log surface density")
	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma), color=:black)

	
	x = LinRange(-1, 2.5, 1000)

	R = 10 .^ x

	for i in 1:100
		M = sum(prof_obs.counts)
		prof = LilGuys.Exp2D(R_s=df.R_s[i], M=M* (1-df.f_outer[i]))
		prof_outer = LilGuys.Exp2D(R_s=df.R_s_outer[i], M=M * df.f_outer[i])
		y = @. log10(LilGuys.surface_density(prof, R))
		y2 = @. log10(LilGuys.surface_density(prof_outer, R))
		lw = theme(:linewidth)[]/2
		lines!(x, y, color=COLORS[1], alpha=0.03, linewidth=lw)
		lines!(x, y2, color=COLORS[2], alpha=0.03, linewidth=lw)
		# lines!(x, log10.(10 .^ y2 .+ 10 .^ y), color=:black, alpha=0.1)
	end


	
	ylims!(-6, 3)

	ax

end

# ╔═╡ 23f21352-f175-414c-8f21-e6a27f73d634
styles = Dict(
	"S+23" => (;markersize=8, color=COLORS[4], marker=:star5),
	"T+23" => (;),
	"apogee" => (;color=COLORS[3], marker=:rect),
	"p+20" => (;)
)

# ╔═╡ c2b44dcf-cc5e-4630-b4c7-46fcb1725b66
let
	fig = Figure()

	ax_scl_dens = plot_samples(fig[1, 1], samples_scl,Utils.load_expected_density_profile("sculptor") )
	vlines!(r_trans_scl, color=COLORS[2], linewidth=theme(:linewidth)[]/2)

	ylims!(-3, 2)
	ax_scl_dens.title = "Sculptor"
	
	ax_scl = Axis(fig[2,1], ylabel="[Fe/H]", xlabel=L"$\log\, R_\textrm{ell}$ / arcmin")
	for study in ["T+23", "apogee", "S+23"]
		df = metals_scl[metals_scl.study .== study, :]
		scatter!(log10.(df.R_ell), df.fe_h, markersize=2; styles[study]...)
	end
	
	median_plot!(log10.(metals_scl.R_ell), metals_scl.fe_h, metals_scl.fe_h_err)
	vlines!(r_trans_scl, color=COLORS[2], linewidth=theme(:linewidth)[]/2)
	ylims!(-3.5, -0.5)
	xlims!(-0.2, 2.2)



	ax_umi_dens = plot_samples(fig[1, 2], samples_umi,Utils.load_expected_density_profile("ursa_minor") )
	ax_umi_dens.title = "Ursa Minor"
	ylims!(-3.5, 1.5)
	vlines!(r_trans_umi, color=COLORS[2], linewidth=theme(:linewidth)[]/2)

	
	ax_umi = Axis(fig[2,2], ylabel="[Fe/H]", xlabel=L"$\log\, R_\textrm{ell}$ / arcmin")
	for study in ["p+20", "apogee", "S+23"]
		df = metals_umi[metals_umi.study .== study, :]
		scatter!(log10.(df.R_ell), df.fe_h, markersize=2; styles[study]...)
	end
	
	median_plot!(log10.(metals_umi.R_ell), metals_umi.fe_h, metals_umi.fe_h_err)
	vlines!(r_trans_umi, color=COLORS[2], linewidth=theme(:linewidth)[]/2)

	xlims!(-0.5, 2.2)

	linkxaxes!(ax_scl, ax_scl_dens)
	linkxaxes!(ax_umi, ax_umi_dens)
	#hidexdecorations!(ax2, ticks=false, minorticks=false)
	hidexdecorations!(ax_scl_dens, ticks=false, minorticks=false)
	hidexdecorations!(ax_umi_dens, ticks=false, minorticks=false)

	#rowsize!(fig.layout, 0, Relative(1/2 / (2+1/2)))

	@savefig "scl_umi_fe_h_gradient"
	fig
end

# ╔═╡ d3a907cc-241e-49d8-8d91-2fd9061e520b
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 7564bff6-1ba9-4255-a94a-606450085eba
import TOML

# ╔═╡ 661caeeb-75ac-4af2-af29-daba92f30970
get_R_h(galaxyname) = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))["R_h"]

# ╔═╡ 157c6a2c-1078-4a88-8ea5-3dc4f0802916
get_R_h("sculptor")

# ╔═╡ b149555c-f132-43c6-b8b0-e89336407178
function jax_profiles(galaxyname, B, r_s, ell)
	R_s = get_R_h(galaxyname) / α_exp
	
	M = 1
	prof_inner = LilGuys.Exp2D(R_s=R_s, M=M)
	
	R_s_outer =  r_s * sqrt(1 - ell) * 60
	M_outer = M * B * (R_s_outer/R_s)^2
	prof_outer = LilGuys.Exp2D(R_s=R_s_outer, M=M_outer)

	return prof_inner, prof_outer
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═af2976f5-b269-465b-9e1d-3a6f62230128
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═f4b51722-b44b-44d7-a293-1e186a4c37cb
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═157c6a2c-1078-4a88-8ea5-3dc4f0802916
# ╠═f633c12c-5351-4169-911e-2afff2403fea
# ╠═04686bc6-3888-4df3-8f0c-502cc20ce86a
# ╠═0bd742a6-e683-4337-b195-e24d0e9db65c
# ╠═1ca18923-c3cd-4ea0-9e3d-adec98ccfe6e
# ╠═fdb0e872-0ba0-4471-9306-8e9b8760bc29
# ╠═1a95a435-acdc-45d9-8eb4-5d8edaaa98cf
# ╠═8366afbb-c528-40cc-8035-18927f77ca3b
# ╠═689bb4f2-c1f7-4ec8-a11f-4b4078160bc5
# ╠═ca1d0e3b-5b47-44a4-9ee0-f7583615908c
# ╠═fc5aaacc-72b1-41b1-9816-bd0ac47fff55
# ╠═e5e5e8e7-0ff3-4bee-b6d3-4f425dc23c02
# ╠═e6eb753c-7989-4bf1-bc90-6027096762fe
# ╠═cfa85373-0eed-482b-a68a-1f3516145fb5
# ╠═3288393c-a1b0-4772-815f-965be98d5a29
# ╠═52e12e1b-740d-4ae4-886f-1f3f849b5d7e
# ╠═8bbf2aa6-0804-4f18-b5e6-4b919635663e
# ╠═23f21352-f175-414c-8f21-e6a27f73d634
# ╠═c2b44dcf-cc5e-4630-b4c7-46fcb1725b66
# ╠═661caeeb-75ac-4af2-af29-daba92f30970
# ╠═d3a907cc-241e-49d8-8d91-2fd9061e520b
# ╠═b149555c-f132-43c6-b8b0-e89336407178
# ╠═7564bff6-1ba9-4255-a94a-606450085eba
