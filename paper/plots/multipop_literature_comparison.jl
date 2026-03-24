### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ e252e94a-2704-11f1-a85f-65d8ce46059d
begin
	import Pkg;
	Pkg.activate()
	using LilGuys
	using Arya, CairoMakie

	import CSV
	import DataFrames: DataFrame
	import TOML

	FIGDIR = "./figures"
end

# ╔═╡ 32ca3ab5-1093-4822-9539-64290b93c657
include("./paper_style.jl")

# ╔═╡ ce81066a-e9ed-405a-a423-dc65f3b7e70e
color_inner = COLORS[1]

# ╔═╡ aee42c04-caff-4f16-b1b0-f0c0197b464a
color_outer = COLORS[4]

# ╔═╡ 2698ead1-29f8-4891-8fd3-221d4afc7893
color_combined = COLORS[5]

# ╔═╡ d57b7e40-32ee-46c9-84b4-3690e1b73965
CairoMakie.activate!(type=:png)

# ╔═╡ aaeb6eea-5e92-494a-8d4f-aa3b9efb66d6
fit_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/mcmc/summary.mcmc_2exp.csv"), DataFrame)

# ╔═╡ c9ffa9d2-82e1-46dd-b04b-5320af3cc336
fit_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/mcmc/summary.mcmc_2exp.csv"), DataFrame)

# ╔═╡ 9f26b06b-d3c6-47e0-9596-a67c9fd8faba
module Utils
	include("./model_utils.jl")
end

# ╔═╡ 73535d55-aac5-423f-b14c-c42bc412897d
function plot_composite_density!(profs, scale, labels=["inner", "outer", "third"])
	x = LinRange(-1, 2.5, 1000)
	R = 10 .^ x

	y_tot = zeros(length(R))

	ys = []
	for prof in profs
		y = @. LilGuys.surface_density(prof, R) * scale
		push!(ys, y)
		y_tot .+= y
	end

	# lines!(log10.(R), log10.(y_tot), color=color_combined, depth_shift=-1)
	r_trans = LilGuys.find_zero(log_r -> LilGuys.surface_density(profs[1], 10^log_r) - LilGuys.surface_density(profs[2], 10^log_r), 0.5, linewidth=theme(:linewidth)[]/2)

	@info "r_trans = $(exp10(r_trans))"

	vlines!(r_trans, color=:black, alpha=0.2)


	for (i, y) in enumerate(ys)
		lines!(log10.(R), log10.(y), linewidth=theme(:linewidth)[]/2, 
			   linestyle = [:solid, :dash, :dot][i], label=labels[i], 
			   color=[color_inner, color_outer, COLORS[3]][i])
	end

end

# ╔═╡ d29b3b14-347c-41b5-89f4-e65af22a438c
function plot_obs_profile(gs, prof_obs)
	ax = Axis(gs, xlabel = L"$\log\ R_\textrm{ell}$ / arcmin", ylabel = L"log $\Sigma_\star$ / arcmin$^{-2}$",
			 yticks = -6:2:2,)
	
	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma), color=:black, markersize=4, marker=:rect)

	
	ylims!(-6, 3)
	ax
end

# ╔═╡ e0adf5d7-9a9b-48f9-acca-a48d0856edfa
function get_param(exp_fit, param)
	idx = only(findall(exp_fit.parameters .== [param]))
	return exp_fit.median[idx]
end

# ╔═╡ ce7a184c-488e-4034-8c1e-96a8871b34b6
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

# ╔═╡ b9cbee73-550f-49a7-af88-7029ca22369b
function panel_label!(s)
	fs = 0.8 * theme(:fontsize)[]
	text!(1.0, 1.0, text=s, space=:relative, align=(:right, :top), fontsize=fs, offset=(-fs/2, -fs/2))
end

# ╔═╡ 7e01de64-437b-4415-9e53-82474f8e42f7
function plot_scl_this_work(gs)
	prof_obs = Utils.load_expected_density_profile("sculptor")
	M = sum(prof_obs.counts)
	ax = plot_obs_profile(gs, prof_obs)

	profs = get_profiles(fit_scl)
	plot_composite_density!(profs, M)

	ax
end


# ╔═╡ f7d851c5-001b-4832-957d-51024927f067
function plot_umi_this_work(gs)
	prof_obs = Utils.load_expected_density_profile("ursa_minor")
	M = sum(prof_obs.counts)
	ax = plot_obs_profile(gs, prof_obs)

	profs = get_profiles(fit_umi)
	plot_composite_density!(profs, M)

	ax
end

# ╔═╡ 5dd4dc5b-ff07-4549-bc29-625da7ccf3ed
r_h_ap24 = [0.129, 0.27] * 60 * sqrt(1 - 0.36)# spherical plummer models...

# ╔═╡ fc1a5729-21b6-46ef-a50c-c829bd77aa90
r_h_ap24_3 = [0.126, 0.26, 1.32] * 60 * sqrt(1 - 0.36)# spherical plummer models...

# ╔═╡ 7ef3f0b0-31bf-45b4-aa3b-21c13261c081
f_outer_ap24 = 0.66

# ╔═╡ 342ead8a-200f-442f-aa6f-88dcccc65824
function plot_scl_ap(gs)

	prof_obs = Utils.load_expected_density_profile("sculptor")
	M = sum(prof_obs.counts)
	ax = plot_obs_profile(gs, prof_obs)

	prof_ap_inner = LilGuys.Plummer(r_s=r_h_ap24[1], M=(1 - f_outer_ap24))
	prof_ap_outer = LilGuys.Plummer(r_s=r_h_ap24[2], M =f_outer_ap24)
	prof_ap_3 = LilGuys.Plummer(r_s=r_h_ap24_3[3], M = 0.017)
	plot_composite_density!([prof_ap_inner, prof_ap_outer, prof_ap_3], M, [nothing, nothing, "third"])

	ax
end


# ╔═╡ a59ba0e9-ab71-4351-a4e4-483dbfc3beb2
r_h_pace20 = [11.444078065902497, 23.159557192656436] # their MMT models converted to arcmin given their adopted distance

# ╔═╡ cec58838-89db-42aa-a2f6-478e9cb5d2a5
f_outer_pace20 = 0.24

# ╔═╡ 7ac5cf51-6b7f-49b3-ae8e-df17499f8766
function plot_umi_pace(gs)
	prof_obs = Utils.load_expected_density_profile("ursa_minor")
	M = sum(prof_obs.counts)
	ax = plot_obs_profile(gs, prof_obs)
	prof_inner = LilGuys.Plummer(r_s=r_h_pace20[1], M=(1 - f_outer_pace20))
	prof_outer = LilGuys.Plummer(r_s=r_h_pace20[2], M =f_outer_pace20)


	plot_composite_density!([prof_inner, prof_outer], M)

end

# ╔═╡ 2c30e0e2-dd4e-49e4-80de-d9d08bd8819e
let
	fig = Figure()

	ax = plot_scl_this_work(fig[1,1])
	ax.title[] = "Sculptor"
	ax.ylabel[] = ""
	hidexdecorations!(ticks=false, minorticks=false)
	axislegend(position=:lb)
	panel_label!("this work")


	ax = plot_scl_ap(fig[2, 1], )
	ax.ylabel[] = ""
	axislegend(position=:lb)
	panel_label!("AP+24")

	ax = plot_umi_this_work(fig[1, 2])
	ax.title[] = "Ursa Minor"
	panel_label!("this work")

	hidexdecorations!(ticks=false, minorticks=false)
	hideydecorations!(ticks=false, minorticks=false)

	plot_umi_pace(fig[2,2])
	hideydecorations!(ticks=false, minorticks=false)
	panel_label!("Pace+20")

	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)

	Label(fig[:, 0],  L"log $\Sigma_\star$ / arcmin$^{-2}$", rotation=π/2)

	@savefig "decomposition_comparison"
	fig
end

# ╔═╡ 1c7df7e0-bf03-4768-8673-bbbee299f0fc
let
	fig = Figure()
	prof_obs = Utils.load_expected_density_profile("ursa_minor")
	
	M = sum(prof_obs.counts)
	
	ax_left = plot_obs_profile(fig[1, 1], prof_obs)

	profs = get_profiles(fit_umi)
	plot_composite_density!(profs, M)
	


	ax_right = plot_obs_profile(fig[1, 2], prof_obs)
	prof_inner = LilGuys.Plummer(r_s=r_h_pace20[1], M=(1 - f_outer_pace20))
	prof_outer = LilGuys.Plummer(r_s=r_h_pace20[2], M =f_outer_pace20)


	plot_composite_density!([prof_inner, prof_outer], M)
	fig
end

# ╔═╡ Cell order:
# ╠═e252e94a-2704-11f1-a85f-65d8ce46059d
# ╠═ce81066a-e9ed-405a-a423-dc65f3b7e70e
# ╠═aee42c04-caff-4f16-b1b0-f0c0197b464a
# ╠═2698ead1-29f8-4891-8fd3-221d4afc7893
# ╠═32ca3ab5-1093-4822-9539-64290b93c657
# ╠═d57b7e40-32ee-46c9-84b4-3690e1b73965
# ╠═aaeb6eea-5e92-494a-8d4f-aa3b9efb66d6
# ╠═c9ffa9d2-82e1-46dd-b04b-5320af3cc336
# ╠═9f26b06b-d3c6-47e0-9596-a67c9fd8faba
# ╠═73535d55-aac5-423f-b14c-c42bc412897d
# ╠═d29b3b14-347c-41b5-89f4-e65af22a438c
# ╠═ce7a184c-488e-4034-8c1e-96a8871b34b6
# ╠═e0adf5d7-9a9b-48f9-acca-a48d0856edfa
# ╠═b9cbee73-550f-49a7-af88-7029ca22369b
# ╠═2c30e0e2-dd4e-49e4-80de-d9d08bd8819e
# ╠═7e01de64-437b-4415-9e53-82474f8e42f7
# ╠═f7d851c5-001b-4832-957d-51024927f067
# ╠═7ac5cf51-6b7f-49b3-ae8e-df17499f8766
# ╠═342ead8a-200f-442f-aa6f-88dcccc65824
# ╠═5dd4dc5b-ff07-4549-bc29-625da7ccf3ed
# ╠═fc1a5729-21b6-46ef-a50c-c829bd77aa90
# ╠═7ef3f0b0-31bf-45b4-aa3b-21c13261c081
# ╠═a59ba0e9-ab71-4351-a4e4-483dbfc3beb2
# ╠═cec58838-89db-42aa-a2f6-478e9cb5d2a5
# ╠═1c7df7e0-bf03-4768-8673-bbbee299f0fc
