### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ d5bec398-03e3-11ef-0930-f3bd4f3c64fd
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames
	using Measurements
	
	using GLMakie
	import TOML
		
	import LilGuys as lguys
	using Arya

end

# ╔═╡ 91d87251-c9c6-467f-9cae-4f30bfea8acc
include("filter_utils.jl")

# ╔═╡ 48caecb2-180c-4ce4-a57b-6fed82328b01
md"""
See filter notebook as well. This notebook is used to understand the background density
"""

# ╔═╡ 47b0d3e6-a79b-4f49-b847-e708e5d6aabf
md"""
 # setup
"""

# ╔═╡ ff92927e-b078-45fd-9c13-1ce5a009d0bb
red = COLORS[6];

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "sculptor/simple"

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
begin 
	param_file = "$name.toml"

	params_json = read_file(param_file)
	params = DensityParams(params_json)
end

# ╔═╡ 4093a7d6-2f74-4c37-a4a8-270934ede924
md"""
# functions
"""

# ╔═╡ 2b9cb3d6-e6ec-4c85-9c27-0f5e2090a0ff
begin
	all_stars_unfiltered = DataFrame(FITS(params.filename)[2])
	xi, eta = lguys.to_tangent(all_stars_unfiltered.ra, all_stars_unfiltered.dec, params.ra, params.dec)
	
	r_ell = 60lguys.calc_r_ell(xi, eta, params.ecc, params.PA)

	all_stars_unfiltered[!, :r_ell] = r_ell
	all_stars_unfiltered[!, :xi] = xi
	all_stars_unfiltered[!, :eta] = eta


end

# ╔═╡ 9371bcb3-0e26-4162-8a40-dc1bf1dacdda
r_max = 60lguys.calc_r_max(all_stars_unfiltered.ra, all_stars_unfiltered.dec,
	params.ecc, params.PA
)

# ╔═╡ 0307085b-816f-469f-8811-60af02cfcb67
# ╠═╡ disabled = true
#=╠═╡
r_max = 60*sqrt(maximum(xi .^ 2 .+ eta .^ 2))
  ╠═╡ =#

# ╔═╡ 29859854-0d17-42ba-8b1d-8788511840c9
r_ell_max = r_max * (1 - params.ecc^2)^(1/4)

# ╔═╡ 6c092147-c295-4ee5-9ee3-6e04c2aaaf98
begin 
	filt_r_ell = r_ell .< r_ell_max
	all_stars = all_stars_unfiltered[filt_r_ell, :]
end

# ╔═╡ 0c498087-0184-4da2-a079-e972dd987712
md"""
The next three plots compare how different the r_ell and xi and eta calculated here are from what is (presumably) in the given catalogue.
"""

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
function plot_tangent(all_stars, members=nothing) 
	fig = Figure()

	ax = plot_all_tangent!(fig[1,1], all_stars, markersize=2,
        color=(:black, 0.2))

	if !isnothing(members )
		plot_all_tangent!(ax, members, markersize=2, color=red)
	end
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	
	return fig
end

# ╔═╡ d7d81ed8-0427-4ee5-8213-320ce5a6711f
function plot_tangent!(grid, all_stars, members=nothing)
	ax = plot_all_tangent!(grid, all_stars, markersize=2,
        color=(:black, 0.2))

	if !isnothing(members )
		plot_all_tangent!(ax, members, markersize=2, color=red)
	end
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	
	return ax
end

# ╔═╡ 75701b7d-e1f1-47f2-800e-c07eec01a4ff
function plot_pms!(grid, all_stars, members=nothing; da=10)
	ax = Axis(grid, limits=(-da, da, -da, da), aspect=1,
	    xlabel=L"\mu_{\alpha *} / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}")
	scatter!(ax, all_stars.pmra, all_stars.pmdec, 
	    color=(:black, 0.2), markersize=1)

	if !isnothing(members)
		scatter!(ax, members.pmra, members.pmdec, 
		    color=(red, 1), markersize=1)
	end

	ax
end

# ╔═╡ bffe11bd-4233-4a7c-9411-0dfb1ac79077
function plot_pms(all_stars, members=nothing)
	fig = Figure()
	plot_pms!(fig[1, 1], all_stars, members)
	fig
end

# ╔═╡ 0f002b56-8b8f-4025-8d7b-fb51423e8da0
function plot_cmd!(grid, all_stars, members=nothing)
	ax = Axis(grid, aspect=1,
	    limits=(-0.5, 3, 10, 22), yreversed=true,
	    xlabel="bp - rp", 
	    ylabel="G",)
	scatter!(ax, all_stars.bp_rp, all_stars.phot_g_mean_mag, 
	    color=(:black, 0.2), markersize=1)

	if !isnothing(members)
		scatter!(ax, members.bp_rp, members.phot_g_mean_mag, 
		    color=(red, 1), markersize=1)
	end
	
	ax.xgridvisible = false
	ax.ygridvisible = false

	ax
end

# ╔═╡ ca345f40-6743-4d43-b1d4-f509b120512a
function plot_cmd(all_stars, members=nothing)
	fig = Figure()
	ax = plot_cmd!(fig[1,1], all_stars, members)
	return fig
end

# ╔═╡ 049ff11e-c04c-41d9-abf1-ec040b799649
function plot_parallax!(grid, all_stars, members=nothing)
	da = 15
	ax = Axis(grid, aspect=1, limits=(-10, 10, 0, 4),
	    xlabel=L"\varpi / \mathrm{ mas }", 
	    ylabel=L"\delta\varpi / \mathrm{ mas }")
	scatter!(ax, all_stars.parallax, all_stars.parallax_error, 
	    color=(:grey, 0.2), markersize=1)

	if !isnothing(members)
		scatter!(ax, members.parallax, members.parallax_error, 
		    color=(red, 1), markersize=1)
	end
	ax
end

# ╔═╡ 2ec0ac71-8bb4-4459-a5a7-049ab711075e
function plot_parallax(all_stars, members=nothing)
	fig = Figure()
	plot_parallax!(fig[1,1], all_stars, members)
	fig
end

# ╔═╡ b57a31e4-5394-4c32-9020-5ade8a1018c1
function plot_sample(all_stars, members=nothing)
	@info plot_tangent(all_stars, members)
	@info plot_parallax(all_stars, members)
	@info plot_cmd(all_stars, members)
	@info plot_pms(all_stars, members)
end

# ╔═╡ ab387f41-8d0f-4a1f-8dbf-b90396914c75
maximum(all_stars.parallax_over_error)

# ╔═╡ c5ca7506-5952-4973-879e-e9848bb72d03
plot_sample(all_stars)

# ╔═╡ f74f597a-908b-4569-8de8-219abba31afd
filt_parallax = apply_filter(all_stars, parallax_filter, params.dist, params.dist_err, params.n_sigma_dist)

# ╔═╡ eab8e0d4-0a40-48b6-9689-0c102d883a96
filt_pm = apply_filter(all_stars, pm_filter, params.pmra, params.pmdec, params.dpm)

# ╔═╡ c15d2045-f6c2-42d1-947e-fb4b31daeb99
filt_ruwe = apply_filter(all_stars, max_filter, :ruwe, 1.3)

# ╔═╡ b0e0a366-8bac-4cf6-8b40-33d517b41e47
filt_cmd = apply_filter(all_stars, cmd_filter, params.cmd_cut)

# ╔═╡ 29dba8f1-b4a6-41da-8323-437447c9d888
filt = filt_cmd .& filt_pm .& filt_parallax .& filt_ruwe

# ╔═╡ 620aa7d7-0b02-475d-80a9-1e721d6144bf
mymemb = all_stars[filt, :]

# ╔═╡ 49ae0572-5d6b-4935-bc95-0a845bb3df2f
md"""
# Background density
"""

# ╔═╡ 31096853-9eaa-40e2-90aa-b248df77f73f
lguys.calc_properties

# ╔═╡ 65a10161-cbeb-49fe-b8d1-075ffe346e43
function get_density(df)
	r = df.r_ell
	props = lguys.calc_properties(r, bins=40, normalization=false)

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

# ╔═╡ 16bd1702-b471-44da-96df-aca63c46a604
plot_sample_w_dens(all_stars, 
	all_stars[filt, :]
)

# ╔═╡ e2dc5d97-0c8f-49dd-bb1a-320b45132cca
sum(all_stars.r_ell .< 4) / π /4^2

# ╔═╡ fb95bc32-9cd2-485e-96f1-e8bebfbc8b59
let
	fig = Figure()
	
	ax = plot_density!(fig[1,1], all_stars, all_stars[filt_cmd, :])

	hlines!(-0.761)
	text!(0, -0.761, text="hi")
	text!(ax, 0.2, 0.2, space=:relative, text="100\nhi")
	fig
end

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
# ╠═91d87251-c9c6-467f-9cae-4f30bfea8acc
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╟─4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═2b9cb3d6-e6ec-4c85-9c27-0f5e2090a0ff
# ╠═9371bcb3-0e26-4162-8a40-dc1bf1dacdda
# ╠═0307085b-816f-469f-8811-60af02cfcb67
# ╠═29859854-0d17-42ba-8b1d-8788511840c9
# ╠═6c092147-c295-4ee5-9ee3-6e04c2aaaf98
# ╟─0c498087-0184-4da2-a079-e972dd987712
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
# ╠═d7d81ed8-0427-4ee5-8213-320ce5a6711f
# ╠═bffe11bd-4233-4a7c-9411-0dfb1ac79077
# ╠═75701b7d-e1f1-47f2-800e-c07eec01a4ff
# ╠═0f002b56-8b8f-4025-8d7b-fb51423e8da0
# ╠═ca345f40-6743-4d43-b1d4-f509b120512a
# ╠═049ff11e-c04c-41d9-abf1-ec040b799649
# ╠═2ec0ac71-8bb4-4459-a5a7-049ab711075e
# ╠═b57a31e4-5394-4c32-9020-5ade8a1018c1
# ╠═ab387f41-8d0f-4a1f-8dbf-b90396914c75
# ╠═c5ca7506-5952-4973-879e-e9848bb72d03
# ╠═f74f597a-908b-4569-8de8-219abba31afd
# ╠═eab8e0d4-0a40-48b6-9689-0c102d883a96
# ╠═c15d2045-f6c2-42d1-947e-fb4b31daeb99
# ╠═b0e0a366-8bac-4cf6-8b40-33d517b41e47
# ╠═29dba8f1-b4a6-41da-8323-437447c9d888
# ╠═620aa7d7-0b02-475d-80a9-1e721d6144bf
# ╟─49ae0572-5d6b-4935-bc95-0a845bb3df2f
# ╠═31096853-9eaa-40e2-90aa-b248df77f73f
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
# ╠═16bd1702-b471-44da-96df-aca63c46a604
# ╠═e2dc5d97-0c8f-49dd-bb1a-320b45132cca
# ╠═fb95bc32-9cd2-485e-96f1-e8bebfbc8b59
# ╠═b76fbfe9-8a0d-4816-ae81-6fb8597aaf80
# ╠═81de6949-c021-49ea-82f5-b300c03a2cc0
# ╠═12220e1c-9a98-447f-b3a4-b579b2989b68
