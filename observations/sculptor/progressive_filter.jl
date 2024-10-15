### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d5bec398-03e3-11ef-0930-f3bd4f3c64fd
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames
	using Measurements
	
	using CairoMakie
	import TOML
		
	import LilGuys as lguys
	using Arya

end

# ╔═╡ 91d87251-c9c6-467f-9cae-4f30bfea8acc
include("../../utils/gaia_filters.jl")

# ╔═╡ 48caecb2-180c-4ce4-a57b-6fed82328b01
md"""
See filter notebook as well. This notebook is used to understand the background density
"""

# ╔═╡ 47b0d3e6-a79b-4f49-b847-e708e5d6aabf
md"""
 # setup
"""

# ╔═╡ 25066fa5-a811-47e9-8579-5b1e167e2a86
import DensityEstimators: histogram

# ╔═╡ ff92927e-b078-45fd-9c13-1ce5a009d0bb
red = COLORS[6];

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ c15da236-6c47-4144-b17a-a28578ba619a
galaxy_dir = "processed"

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "centre"

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
begin 
	param_file = joinpath(galaxy_dir, "$name.toml")

	params_json = read_paramfile(param_file)
	params = GaiaFilterParams(params_json)
end

# ╔═╡ 4093a7d6-2f74-4c37-a4a8-270934ede924
md"""
# loading stars
"""

# ╔═╡ 2b9cb3d6-e6ec-4c85-9c27-0f5e2090a0ff
all_stars_unfiltered = read_gaia_stars(joinpath(galaxy_dir, params.filename), params)

# ╔═╡ 9371bcb3-0e26-4162-8a40-dc1bf1dacdda
r_ell_max = 60lguys.calc_r_max(all_stars_unfiltered.ra, all_stars_unfiltered.dec,
	params.ellipticity, params.PA, centre=(params.ra, params.dec)
)

# ╔═╡ 5e80e0ad-9ed9-467d-a5d6-fb585bfecf18
maximum(all_stars_unfiltered.r_ell)

# ╔═╡ f14e8857-340f-4f30-85da-28560361360f
filt_r_ell = apply_filter(all_stars_unfiltered, r_ell_filter, params) .& ang_dist_filter(all_stars_unfiltered, params)

# ╔═╡ 6c092147-c295-4ee5-9ee3-6e04c2aaaf98
all_stars = all_stars_unfiltered[filt_r_ell, :]

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ 9fd9954c-f670-4b69-b66b-5d392825fd97
import LilGuys.Plots as LP

# ╔═╡ d7d81ed8-0427-4ee5-8213-320ce5a6711f
function plot_tangent!(grid, all_stars, members=nothing; markersize=2, kwargs...)
	ax = Axis(grid, 
		xlabel = "xi / degrees",
		ylabel = "eta / degrees",
	)
		
	scatter!(ax, all_stars.xi, all_stars.eta, markersize=2,
        color=(:black, 0.2))

	if !isnothing(members )
		scatter!(ax, members.xi, members.eta; 
		markersize=markersize, color=red, kwargs...)
	end
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	
	return ax
end

# ╔═╡ d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
function plot_tangent(all_stars, members=nothing; markersize=2, kwargs...) 
	fig = Figure()

	ax = plot_tangent!(fig[1,1], all_stars, markersize=2,
        color=(:black, 0.2))

	if !isnothing(members )
		scatter!(ax, members.xi, members.eta; markersize=markersize, color=red, kwargs...)
	end
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	
	return fig
end

# ╔═╡ 75701b7d-e1f1-47f2-800e-c07eec01a4ff
function plot_pms!(grid, all_stars, members=nothing; da=10, markersize=2, kwargs...)
	ax = Axis(grid, limits=(-da, da, -da, da), aspect=1,
	    xlabel=L"\mu_{\alpha *} / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}")
	scatter!(ax, all_stars.pmra, all_stars.pmdec, 
	    color=(:black, 0.2), markersize=1)

	if !isnothing(members)
		scatter!(ax, members.pmra, members.pmdec, 
		    color=red, markersize=markersize; kwargs...)
	end

	ax
end

# ╔═╡ bffe11bd-4233-4a7c-9411-0dfb1ac79077
function plot_pms(all_stars, members=nothing; kwargs...)
	fig = Figure()
	plot_pms!(fig[1, 1], all_stars, members; kwargs...)
	fig
end

# ╔═╡ 0f002b56-8b8f-4025-8d7b-fb51423e8da0
function plot_cmd!(grid, all_stars, members=nothing; markersize=2, kwargs...)
	ax = Axis(grid, aspect=1,
	    limits=(-0.5, 3, 15, 22), yreversed=true,
	    xlabel="bp - rp", 
	    ylabel="G",)
	scatter!(ax, all_stars.bp_rp, all_stars.phot_g_mean_mag, 
	    color=(:black, 0.2), markersize=1)

	if !isnothing(members)
		scatter!(ax, members.bp_rp, members.phot_g_mean_mag;
		    color=red, markersize=markersize, kwargs...)
	end
	
	ax.xgridvisible = false
	ax.ygridvisible = false

	ax
end

# ╔═╡ ca345f40-6743-4d43-b1d4-f509b120512a
function plot_cmd(all_stars, members=nothing; kwargs...)
	fig = Figure()
	ax = plot_cmd!(fig[1,1], all_stars, members; kwargs...)
	return fig
end

# ╔═╡ 049ff11e-c04c-41d9-abf1-ec040b799649
function plot_parallax!(grid, all_stars, members=nothing; markersize=2, kwargs...)
	da = 15
	ax = Axis(grid, aspect=1, limits=(-10, 10, 0, 4),
	    xlabel=L"\varpi / \mathrm{ mas }", 
	    ylabel=L"\delta\varpi / \mathrm{ mas }")
	scatter!(ax, all_stars.parallax, all_stars.parallax_error, 
	    color=(:grey, 0.2), markersize=1)

	if !isnothing(members)
		scatter!(ax, members.parallax, members.parallax_error; 
		    color=red, markersize=markersize, kwargs...)
	end
	ax
end

# ╔═╡ 2ec0ac71-8bb4-4459-a5a7-049ab711075e
function plot_parallax(all_stars, members=nothing; kwargs...)
	fig = Figure()
	plot_parallax!(fig[1,1], all_stars, members; kwargs...)
	fig
end

# ╔═╡ b57a31e4-5394-4c32-9020-5ade8a1018c1
function plot_sample(all_stars, members=nothing; kwargs...)
	@info plot_tangent(all_stars, members; kwargs...)
	@info plot_parallax(all_stars, members; kwargs...)
	@info plot_cmd(all_stars, members; kwargs...)
	@info plot_pms(all_stars, members; kwargs...)
end

# ╔═╡ c5ca7506-5952-4973-879e-e9848bb72d03
plot_sample(all_stars)

# ╔═╡ 933be42e-c33d-4afe-a40d-f016d1839519
params

# ╔═╡ f3eac4e6-0b4b-4d54-83ac-213b78fcb522
md"""
# filters
"""

# ╔═╡ f74f597a-908b-4569-8de8-219abba31afd
filt_parallax = parallax_filter(all_stars, params)

# ╔═╡ eab8e0d4-0a40-48b6-9689-0c102d883a96
filt_pm = pm_filter(all_stars, params)

# ╔═╡ 20fe601c-874c-4a37-b3fa-b229232a0f4f
params.ruwe_max

# ╔═╡ b0e0a366-8bac-4cf6-8b40-33d517b41e47
filt_cmd = cmd_filter(all_stars, params)

# ╔═╡ 90c347ff-74fe-4032-93ec-931b181d3829
filt_ruwe = ruwe_filter(all_stars, params)

# ╔═╡ 6ee77988-5d0b-49db-87cb-a041c1ed4a28
plot_sample(all_stars, all_stars[filt_pm .& filt_parallax .& filt_ruwe, :]; markersize=5, alpha=0.3)

# ╔═╡ 29dba8f1-b4a6-41da-8323-437447c9d888
filt = filt_cmd .& filt_pm .& filt_parallax .& filt_ruwe

# ╔═╡ 99726828-62ad-465a-a961-4340b143cf6a
plot_sample(all_stars, all_stars[filt, :]; markersize=5, alpha=0.3)

# ╔═╡ cf7b3e70-2a92-4fe3-928f-6f842899893c
md"""
# Filter verification
"""

# ╔═╡ 8775835f-3e50-4a88-a3fe-6439a1bcbd22
let 
	fig = plot_tangent(all_stars_unfiltered, all_stars)
	if params.max_ang_dist !== nothing
		arc!(Point2f(0), params.max_ang_dist, 0, 2π)
	end
	fig
end

# ╔═╡ 740a8f04-5f19-466b-9fcd-767de3ad06b9
let 
	fig, ax = FigAxis(
		xlabel = "ruwe",
		ylabel = "count",
		xscale=log10,
		limits=((0.5, 100), nothing)
	)

	nan_filt = (!).(isnan.(all_stars.ruwe))
	println(sum( @. filt_ruwe & ! nan_filt))
	
	h = histogram(all_stars.ruwe[nan_filt], normalization=:none)
	h1 = histogram(all_stars.ruwe[filt_ruwe .& nan_filt], normalization=:none)
	lines!(midpoints(h.bins), h.values)
	lines!(midpoints(h1.bins), h1.values)
	vlines!(params.ruwe_max, color=:black)
	fig
end

# ╔═╡ 1823b424-5e65-49d7-b1a1-75035b2006df
let
	fig = plot_parallax(all_stars, all_stars[filt_parallax, :])

	x = [-5, 0, 5]
	y = abs.(1/params.n_sigma_dist .* x)
	lines!(x, y)

	fig
end

# ╔═╡ daad6357-9231-4025-a084-267ecc28eac2
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=1, #limits=(-10, 10, 0, 4),
	    xlabel=L"\mu_{\alpha*} / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\delta\mu_{\alpha*} / \mathrm{ mas\, yr^{-1} }",
		limits=(-5, 5, nothing, nothing),
		xgridvisible=false,
		ygridvisible=false,
	)
	
	scatter!(ax, all_stars.pmra, all_stars.pmra_error, 
	    color=(:grey, 0.2), markersize=1)


	scatter!(ax, all_stars[filt, :pmra], all_stars[filt, :pmra_error]; 
		color=red, markersize=2)
	
	if params.n_sigma_pm !== nothing
		x = [-5, 0, 5]
		y = abs.(1/params.n_sigma_pm .* x)
		lines!(x, y)
	end
	
	fig
end

# ╔═╡ 03636b4c-6d14-4ae3-8f72-411eddc60eaa
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=1, #limits=(-10, 10, 0, 4),
	    xlabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\delta\mu_\delta / \mathrm{ mas\, yr^{-1} }",
		limits=(-5, 5, nothing, nothing),
		xgridvisible=false,
		ygridvisible=false,
	)

	
	scatter!(ax, all_stars.pmdec, all_stars.pmdec_error, 
	    color=(:grey, 0.2), markersize=1)


	scatter!(ax, all_stars[filt, :pmdec], all_stars[filt, :pmdec_error]; 
		color=red, markersize=2)

	if params.n_sigma_pm !== nothing
		x = [-5, 0, 5]
		y = abs.(1/params.n_sigma_pm .* x)
		lines!(x, y)
	end
	
	fig
end

# ╔═╡ dbc48359-f4a2-41a3-a26d-57dadf0bb3e4
import StatsBase: median

# ╔═╡ 70aad4cf-0a09-46c2-999e-a9bec307a4f4
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=1, #limits=(-10, 10, 0, 4),
	    xlabel=L"\delta\mu_{\alpha*} / \mathrm{ mas\, yr^{-1}}", 
		title="PM, CMD, and parallax filters",
		ylabel="count"
	)

	hist!(all_stars[filt, :].pmdec_error, bins=40)
	println(median(all_stars[filt, :].pmdec_error))
	fig
end

# ╔═╡ 49ae0572-5d6b-4935-bc95-0a845bb3df2f
md"""
# Background density
"""

# ╔═╡ 65a10161-cbeb-49fe-b8d1-075ffe346e43
function get_density(df)
	r = df.r_ell
	props = lguys.StellarProfile(r, bins=40, normalization=:none)

	println("stars left ", length(r))
	println("counts in last bin ", props.counts[end-2: end])
	println("densities ", (props.log_Sigma .± props.log_Sigma_err)[end-5:end])
	
	return props.log_r, props.log_Sigma, props.log_Sigma_err
end

# ╔═╡ 08b5251d-bbe4-48f6-9cb3-01e0a8364c1d
get_density(all_stars[filt, :])

# ╔═╡ 198ba5a6-9d04-49de-9726-4e1d9be9780f
value = lguys.value

# ╔═╡ b55d72b9-e957-4480-b6db-97c9798b4d68
function plot_density!(grid::GridPosition, all_stars, sample=nothing)
	ax = Axis(grid,
		xlabel="log radius / arcmin",
		ylabel=L"\log\;\Sigma \, / \, \textrm{stars arcmin^{-2}}",
		limits=((-1., 2.3), (-4, 2.2))
	)
	
	plot_density!(ax, all_stars, yerr=false, color=:black)
	if sample !== nothing
		y_end = plot_density!(ax, sample, color=red)
		N = length(sample.r_ell)

		hlines!(value.(y_end))	
		
		text!(ax, 0.1, 0.1, space=:relative, 
			text="$N stars ")
	
		text!(ax, -1, value.(y_end), text=		
			L"\log\Sigma_\textrm{bg} = %$y_end", fontsize=14)
		
	end





	return ax
end

# ╔═╡ bd4212fe-3b70-44ba-85fb-b199840a6f71
function plot_density!(ax::Axis, sample; yerr=true, kwargs...)
	x, y, y_err = get_density(sample)

	
	if !yerr
		y_err = nothing
	end
	errscatter!(ax, x, y, yerr=y_err; kwargs...)

	y_end = (y .± yerr)[end-3:end]
	y_end = minimum(y_end)
end

# ╔═╡ d0c00f3e-f3d3-4cd1-8541-9ae239420174
function plot_density(all_stars, sample=nothing)
	fig = Figure()
	plot_density!(fig[1,1], all_stars, sample)

	return fig
end

# ╔═╡ 1c06e13b-8f5d-414b-8f2a-3d6f051ad495
function plot_sample_w_dens(all_stars, sample=nothing)
	fig = Figure(size=(800, 800))
	plot_tangent!(fig[1,1], all_stars, sample)
	plot_pms!(fig[1,2], all_stars, sample)
	plot_cmd!(fig[2,1], all_stars, sample)
	plot_density!(fig[2,2], all_stars, sample)

	fig
end

# ╔═╡ 289b5daa-a0a3-41c8-9d9c-3a3280139f2a
plot_sample_w_dens(all_stars_unfiltered)

# ╔═╡ 498399c5-43f7-44db-bddd-03ec0f0fdd6b
plot_sample_w_dens(all_stars)

# ╔═╡ fab13f99-662c-4385-b266-029ce2eea233
plot_sample_w_dens(all_stars, 
	all_stars[filt_ruwe, :]
)

# ╔═╡ aac80d15-72c9-4091-befa-1b7fa7127c63
plot_sample_w_dens(all_stars, 
	all_stars[filt_parallax, :]
)

# ╔═╡ 4bfc9046-a5b2-4a72-82a7-3547066d7064
plot_sample_w_dens(all_stars, 
	all_stars[filt_cmd, :]
)

# ╔═╡ d74e062b-d7b0-422e-9f1b-102afee26d7b
plot_sample_w_dens(all_stars, 
	all_stars[filt_pm, :]
)

# ╔═╡ 8e52c3d5-aef1-4b13-9316-dd13f76e2ab3
plot_sample_w_dens(all_stars, 
	all_stars[filt_pm .& filt_parallax, :]
)

# ╔═╡ 16bd1702-b471-44da-96df-aca63c46a604
plot_sample_w_dens(all_stars, 
	all_stars[filt, :]
)

# ╔═╡ e2dc5d97-0c8f-49dd-bb1a-320b45132cca
sum(all_stars.r_ell .< 4) / π /4^2

# ╔═╡ b76fbfe9-8a0d-4816-ae81-6fb8597aaf80
all_stars[filt_cmd, :]

# ╔═╡ 81de6949-c021-49ea-82f5-b300c03a2cc0
log10(60 * 4)

# ╔═╡ 12220e1c-9a98-447f-b3a4-b579b2989b68
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log radius / arcmin",
		ylabel=L"\log\;\Sigma \, / \, \textrm{stars arcmin^{-2}}",
		#limits=(nothing, (-2, 3))
	)
	
	plot_density!(ax, all_stars, label="no filter")
	plot_density!(ax, all_stars[filt_ruwe, :], label="RUWE only")

	plot_density!(ax, all_stars[filt_parallax, :], label="parallax only")
	plot_density!(ax, all_stars[filt_cmd, :], label="CMD only")

	plot_density!(ax, all_stars[filt_pm, :], label="proper motion only")
	plot_density!(ax, all_stars[filt, :], label="all cuts")

	axislegend(position=:lb)

	fig
end

# ╔═╡ Cell order:
# ╟─48caecb2-180c-4ce4-a57b-6fed82328b01
# ╟─47b0d3e6-a79b-4f49-b847-e708e5d6aabf
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═25066fa5-a811-47e9-8579-5b1e167e2a86
# ╠═91d87251-c9c6-467f-9cae-4f30bfea8acc
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═c15da236-6c47-4144-b17a-a28578ba619a
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╠═4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═2b9cb3d6-e6ec-4c85-9c27-0f5e2090a0ff
# ╠═9371bcb3-0e26-4162-8a40-dc1bf1dacdda
# ╠═5e80e0ad-9ed9-467d-a5d6-fb585bfecf18
# ╠═f14e8857-340f-4f30-85da-28560361360f
# ╠═6c092147-c295-4ee5-9ee3-6e04c2aaaf98
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
# ╠═9fd9954c-f670-4b69-b66b-5d392825fd97
# ╠═d7d81ed8-0427-4ee5-8213-320ce5a6711f
# ╠═bffe11bd-4233-4a7c-9411-0dfb1ac79077
# ╠═75701b7d-e1f1-47f2-800e-c07eec01a4ff
# ╠═0f002b56-8b8f-4025-8d7b-fb51423e8da0
# ╠═ca345f40-6743-4d43-b1d4-f509b120512a
# ╠═049ff11e-c04c-41d9-abf1-ec040b799649
# ╠═2ec0ac71-8bb4-4459-a5a7-049ab711075e
# ╠═b57a31e4-5394-4c32-9020-5ade8a1018c1
# ╠═c5ca7506-5952-4973-879e-e9848bb72d03
# ╠═6ee77988-5d0b-49db-87cb-a041c1ed4a28
# ╠═99726828-62ad-465a-a961-4340b143cf6a
# ╠═933be42e-c33d-4afe-a40d-f016d1839519
# ╟─f3eac4e6-0b4b-4d54-83ac-213b78fcb522
# ╠═f74f597a-908b-4569-8de8-219abba31afd
# ╠═eab8e0d4-0a40-48b6-9689-0c102d883a96
# ╠═20fe601c-874c-4a37-b3fa-b229232a0f4f
# ╠═b0e0a366-8bac-4cf6-8b40-33d517b41e47
# ╠═90c347ff-74fe-4032-93ec-931b181d3829
# ╠═29dba8f1-b4a6-41da-8323-437447c9d888
# ╟─cf7b3e70-2a92-4fe3-928f-6f842899893c
# ╠═8775835f-3e50-4a88-a3fe-6439a1bcbd22
# ╠═740a8f04-5f19-466b-9fcd-767de3ad06b9
# ╠═1823b424-5e65-49d7-b1a1-75035b2006df
# ╠═daad6357-9231-4025-a084-267ecc28eac2
# ╠═03636b4c-6d14-4ae3-8f72-411eddc60eaa
# ╠═70aad4cf-0a09-46c2-999e-a9bec307a4f4
# ╠═dbc48359-f4a2-41a3-a26d-57dadf0bb3e4
# ╟─49ae0572-5d6b-4935-bc95-0a845bb3df2f
# ╠═65a10161-cbeb-49fe-b8d1-075ffe346e43
# ╠═08b5251d-bbe4-48f6-9cb3-01e0a8364c1d
# ╠═d0c00f3e-f3d3-4cd1-8541-9ae239420174
# ╠═198ba5a6-9d04-49de-9726-4e1d9be9780f
# ╠═b55d72b9-e957-4480-b6db-97c9798b4d68
# ╠═bd4212fe-3b70-44ba-85fb-b199840a6f71
# ╠═1c06e13b-8f5d-414b-8f2a-3d6f051ad495
# ╠═289b5daa-a0a3-41c8-9d9c-3a3280139f2a
# ╠═498399c5-43f7-44db-bddd-03ec0f0fdd6b
# ╠═fab13f99-662c-4385-b266-029ce2eea233
# ╠═aac80d15-72c9-4091-befa-1b7fa7127c63
# ╠═4bfc9046-a5b2-4a72-82a7-3547066d7064
# ╠═d74e062b-d7b0-422e-9f1b-102afee26d7b
# ╠═8e52c3d5-aef1-4b13-9316-dd13f76e2ab3
# ╠═16bd1702-b471-44da-96df-aca63c46a604
# ╠═e2dc5d97-0c8f-49dd-bb1a-320b45132cca
# ╠═b76fbfe9-8a0d-4816-ae81-6fb8597aaf80
# ╠═81de6949-c021-49ea-82f5-b300c03a2cc0
# ╠═12220e1c-9a98-447f-b3a4-b579b2989b68
