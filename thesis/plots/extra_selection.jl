### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	import PythonCall
	using DataFrames 
	using CSV
	using CairoMakie
end

# ╔═╡ c3af8de5-611a-457e-8ba0-636803c5a76c
using PyFITS

# ╔═╡ 2d5297cd-6a01-4b26-ac77-995b878d765d
using Arya

# ╔═╡ ae29bed0-6700-47f1-8952-35e867ce126b
using OrderedCollections

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using LilGuys
	FIGDIR = "figures"
end

# ╔═╡ 69c98029-165c-407b-9a63-a27e06e30e45
include("paper_style.jl")

# ╔═╡ d3bd7158-ea70-47a0-9800-9bfc08f3557c
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sammple

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ 0107ca79-9c8a-43cf-86af-b99c2560743f
CairoMakie.activate!(type=:png, px_per_unit=2)

# ╔═╡ a4848423-9be3-48d7-98c2-8d1096eb2560
module Utils
	include("gaia_utils.jl")
end

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ 1a0f25ad-2f7b-4c47-9003-3d57beb0d8ac
function load_samples(galaxyname, profname)
	observed_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/" * galaxyname * "/observed_properties.toml")
	
	filepath = joinpath(ENV["DWARFS_ROOT"], "observations", "$galaxyname/density_profiles/$profname.toml")
	filt_kwargs = read_paramfile(filepath)
	filt_kwargs["filename"] = joinpath(dirname(filepath), filt_kwargs["filename"])
	
	profile_kwargs = pop!(filt_kwargs, "profile_kwargs")


	filt_params = GaiaFilterParams(observed_properties; LilGuys.dict_to_tuple(filt_kwargs)...)

	all_stars = read_gaia_stars(filt_params)

	members = select_members(all_stars, filt_params)

	if :F_BEST ∈ names(all_stars)
		best_stars = all_stars[all_stars.F_BEST .== 1.0, :]
	else
		best_stars = all_stars
	end


	return OrderedDict(
		:best => best_stars,
		:best_inner => best_stars[best_stars.R_ell .< observed_properties["R_h"], :],
		:members => members,
	), filt_params

end

# ╔═╡ 55def0bc-9eac-44da-ab24-95efdd4c0cd7
scl_simple, params_scl_simple = load_samples("sculptor", "simple")

# ╔═╡ 123bf199-adbc-44ab-9197-e2cd70a05bde
umi_simple, params_umi_simple = load_samples("ursa_minor", "simple")

# ╔═╡ 7cefea73-5a23-4b8a-b450-a613eeaec708
scl_delve, params_scl_delve = load_samples("sculptor", "delve_new")

# ╔═╡ 084503a0-3783-47f6-ba23-d5f888fc8e5d
params_scl_delve.cmd_y

# ╔═╡ 75078e58-3731-4df2-b9e6-d70fe1a9420d
umi_unions, params_umi_unions = load_samples("ursa_minor", "unions")

# ╔═╡ dc00b7ce-ab13-4c66-816b-89c0daa83ece
obs_props_umi = Utils.get_obs_props("ursa_minor")

# ╔═╡ 4632600f-bb97-454b-a7af-52c7c307d787
obs_props_scl = Utils.get_obs_props("sculptor")

# ╔═╡ 1a317c97-8584-4a71-9836-ffa7dfc2786c
plot_kwargs = Dict(
	:best => (; markersize=0.5, alpha=1, color=RGBf(0.8, 0.8, 0.8)),
	:best_inner => (; markersize=0.5, alpha=1, color=RGBf(0.0, 0.0, 0.0)),
	:members => (; markersize=0.5, alpha=1, color=COLORS[2])
)

# ╔═╡ b30602a3-d3a2-4258-80e0-dbef5be3ee21
bp_rp_label = L"G_\textrm{BP} - G_\textrm{RP}"

# ╔═╡ 57dca48b-bdca-4bc7-8a05-625ea14bb808
G_label = L"G"

# ╔═╡ 6fb63758-c638-427a-8cb3-86b45a054282
label_survey!(surveyname; x=0, y=1, offset=(5, -5)) = text!(x, y, offset=offset, text=surveyname, fontsize=theme(:fontsize)[] * 0.8, align=(:left, :top), space=:relative)

# ╔═╡ ce2d1f6f-8ff3-4137-a979-cf86e57b1ef4
R_h_scl = obs_props_scl["R_h"]

# ╔═╡ 527cd57b-d90f-4230-b771-da2c96e1fbb3
R_h_umi = obs_props_umi["R_h"]

# ╔═╡ bb80dedf-415b-4b65-9c9e-bd11df4f8f9d
scl_delve[:best]

# ╔═╡ 08a5354b-148c-4b2a-8a22-4c298fe01458
umi_unions[:best].R_ell

# ╔═╡ 40911b4f-0a33-4d1a-b6a7-3e9019db5a01
import KernelDensity

# ╔═╡ c6cf0e57-d586-44e6-8ddf-089b59b4ef37
function plot_kde(gs, df, bandwidth=2; 
				  styles=plot_kwargs, x1=:phot_bp_mean_mag, x2=:phot_rp_mean_mag, y=:phot_g_mean_mag, R_max=90, xlabel=L"r-i", ylabel= L"i", 
				  nlevels=10, cmin=-4, kwargs...)

	ax = Axis(gs; 
		xlabel = L"$\xi$ / arcmin",
		ylabel = L"$\eta$ / arcmin",
		xreversed=true,
		limits = ((-R_max, R_max), (-R_max, R_max)),
			  kwargs...
			 )


	scatter!(df.xi, df.eta, alpha=0.1, color=:black, markersize=1, rasterize=true)

	k = KernelDensity.kde((df.xi, df.eta), bandwidth=(bandwidth,bandwidth), boundary=((-R_max, R_max), (-R_max, R_max)))
	cmax = log10(maximum(k.density))
	contour!(k, levels = 10 .^ LinRange(cmax+cmin, cmax, nlevels))


end

# ╔═╡ aaed8b9a-b7ab-4f4d-bbc9-05e698052457
let
	fig = Figure()

	

	plot_kde(fig[1,1], scl_delve[:members], cmin=-2, title="Sculptor")
	label_survey!("DELVE")

	
	plot_kde(fig[1,2], umi_unions[:members], cmin=-1.2, title="Ursa Minor")
	label_survey!("UNIONS")


	rowsize!(fig.layout, 1, Aspect(1, 1))
	resize_to_layout!()

	@savefig "delve_unions_tangent"
	fig
end

# ╔═╡ 75ace0e8-2141-4c15-9754-62bec0456ce6
function plot_tangent(gs, datasets; 
				  styles=plot_kwargs, x1=:phot_bp_mean_mag, x2=:phot_rp_mean_mag, y=:phot_g_mean_mag, R_max=Inf, xlabel=L"r-i", ylabel= L"i", kwargs...)

	ax = Axis(gs; 
		xlabel = L"\xi",
		ylabel = L"\eta",
		xreversed=true,
			  kwargs...
			 )

	for (label, df) in datasets
		xs = df.xi
		ys = df.eta
		scatter!(xs, ys, label=label, rasterize=true, alpha=0.03; styles[label]...)
	end


end

# ╔═╡ b5636d2a-24eb-40a2-9df3-695953245265


# ╔═╡ 25fec3d5-d80d-43ad-b64b-6c889c3b5e41
function plot_cmd(gs, datasets; 
				  styles=plot_kwargs, x1=:phot_bp_mean_mag, x2=:phot_rp_mean_mag, y=:phot_g_mean_mag, R_max=Inf, xlabel=L"r-i", ylabel= L"i", kwargs...)

	ax = Axis(gs; 
		xlabel = xlabel,
		ylabel = ylabel,
		yreversed=true,
			  kwargs...
			 )

	for (label, df) in datasets
		filt = df.R_ell .< R_max
		xs = df[filt, x1] .- df[filt, x2]
		ys = df[filt, y]
		scatter!(xs, ys, label=label, rasterize=true; styles[label]...)
	end


end

# ╔═╡ fc60f85e-10fb-4674-b538-972c273a9ef2
let
	fig = Figure(size=(4.5*72, 4.5*72))



	plot_cmd(fig[1,1], scl_simple, title="Sculptor", xlabel=bp_rp_label, ylabel=G_label)
	xlims!(-0.5, 2.5)
	ylims!(21, 15.5)
	label_survey!(L"Gaia")



	plot_cmd(fig[1,2], umi_simple, title="Ursa Minor", xlabel=bp_rp_label, ylabel=G_label)

	xlims!(-0.5, 2.5)
	ylims!(21, 15.5)

	label_survey!(L"Gaia")

	
	# scl - delve
	plot_cmd(fig[2,1], scl_delve, x1=:mag_psf_g, x2=:mag_psf_r, y=:mag_psf_g, 
			 xlabel = L"g-r", ylabel=L"g",
			)
	xlims!(-1, 2)
	ylims!(24, 17)

	label_survey!("DELVE")


	plot_cmd(fig[2,2], umi_unions, x1=:RMAG, x2=:IMAG, y=:IMAG,)
	xlims!(-1, 1)
	ylims!(24, 17)
	label_survey!("UNIONS")



	@savefig "extra_cmd_selection"
	fig
end

# ╔═╡ 00a96fce-0e9d-4e93-9c3b-5d6e1c86fa9d
function plot_cmd_cut!(filt_params)
	cmd = filt_params.cmd_cut
	N = length(cmd) ÷ 2

	x = [cmd[1 + 2(i-1)] for i in 1:N]
	y = [cmd[2i] for i in 1:N]


	x = [x; x[1]]
	y = [y; y[1]]
	lines!(x, y)
end

# ╔═╡ 871aafc3-e71c-475a-b651-e83a1badf470
let
	fig = Figure()

	plot_cmd(fig[1,1], scl_simple)
	xlims!(-0.5, 2.5)
	ylims!(21, 15.5)
	plot_cmd_cut!(params_scl_simple)

	fig
end

# ╔═╡ 590da3e2-3c3c-445a-a5aa-27a536cbb686
let
	fig = Figure()

	plot_cmd(fig[1,1], umi_simple)

	xlims!(-0.5, 2.5)
	ylims!(21, 15.5)
	plot_cmd_cut!(params_umi_simple)

	fig
end

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═c3af8de5-611a-457e-8ba0-636803c5a76c
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═0107ca79-9c8a-43cf-86af-b99c2560743f
# ╠═2d5297cd-6a01-4b26-ac77-995b878d765d
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═a4848423-9be3-48d7-98c2-8d1096eb2560
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═d3bd7158-ea70-47a0-9800-9bfc08f3557c
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═1a0f25ad-2f7b-4c47-9003-3d57beb0d8ac
# ╠═55def0bc-9eac-44da-ab24-95efdd4c0cd7
# ╠═123bf199-adbc-44ab-9197-e2cd70a05bde
# ╠═7cefea73-5a23-4b8a-b450-a613eeaec708
# ╠═084503a0-3783-47f6-ba23-d5f888fc8e5d
# ╠═75078e58-3731-4df2-b9e6-d70fe1a9420d
# ╠═dc00b7ce-ab13-4c66-816b-89c0daa83ece
# ╠═4632600f-bb97-454b-a7af-52c7c307d787
# ╠═1a317c97-8584-4a71-9836-ffa7dfc2786c
# ╠═b30602a3-d3a2-4258-80e0-dbef5be3ee21
# ╠═57dca48b-bdca-4bc7-8a05-625ea14bb808
# ╠═6fb63758-c638-427a-8cb3-86b45a054282
# ╠═ce2d1f6f-8ff3-4137-a979-cf86e57b1ef4
# ╠═527cd57b-d90f-4230-b771-da2c96e1fbb3
# ╠═fc60f85e-10fb-4674-b538-972c273a9ef2
# ╠═aaed8b9a-b7ab-4f4d-bbc9-05e698052457
# ╠═bb80dedf-415b-4b65-9c9e-bd11df4f8f9d
# ╠═08a5354b-148c-4b2a-8a22-4c298fe01458
# ╠═40911b4f-0a33-4d1a-b6a7-3e9019db5a01
# ╠═c6cf0e57-d586-44e6-8ddf-089b59b4ef37
# ╠═75ace0e8-2141-4c15-9754-62bec0456ce6
# ╠═b5636d2a-24eb-40a2-9df3-695953245265
# ╠═25fec3d5-d80d-43ad-b64b-6c889c3b5e41
# ╠═871aafc3-e71c-475a-b651-e83a1badf470
# ╠═590da3e2-3c3c-445a-a5aa-27a536cbb686
# ╠═00a96fce-0e9d-4e93-9c3b-5d6e1c86fa9d
