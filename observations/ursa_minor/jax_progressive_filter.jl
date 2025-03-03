### A Pluto.jl notebook ###
# v0.20.0

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
See filter notebook as well. This notebook is used to understand the background density using the columns in the J+24 sample.
Also compares the simple method to the J+24 method of filtering.
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
param_file = "processed/fiducial.toml"

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
params = GaiaFilterParams(read_paramfile(param_file))

# ╔═╡ a86b5875-c30a-446a-8897-763b07f76c68
begin 
	simple_file = "processed/simple.toml"
	params_simple = GaiaFilterParams(read_paramfile(simple_file))
end

# ╔═╡ 4093a7d6-2f74-4c37-a4a8-270934ede924
md"""
# functions
"""

# ╔═╡ 9371bcb3-0e26-4162-8a40-dc1bf1dacdda
all_stars_unfiltered = read_gaia_stars("processed/" * params.filename, params)

# ╔═╡ 0307085b-816f-469f-8811-60af02cfcb67
# ╠═╡ disabled = true
#=╠═╡
r_max = 60*sqrt(maximum(xi .^ 2 .+ eta .^ 2))
  ╠═╡ =#

# ╔═╡ 29859854-0d17-42ba-8b1d-8788511840c9
all_stars = all_stars_unfiltered[r_ell_filter(all_stars_unfiltered, params), :]

# ╔═╡ c26a1042-1c5f-4539-992f-022cc7441a1a
maximum(all_stars_unfiltered.r_ell)

# ╔═╡ 9c2da11e-79fb-406e-9c5a-83eefb620814
params.filename

# ╔═╡ 6c092147-c295-4ee5-9ee3-6e04c2aaaf98
begin 

	gaia_all = read_gaia_stars("data/gaia_6deg.fits", params)
	gaia_all = gaia_all[r_ell_filter(gaia_all, params), :]
	
end

# ╔═╡ fecf285c-1d7b-4e91-8108-8c906197cd20
memb_simple = select_members(gaia_all, params_simple)

# ╔═╡ c403ee08-852b-4bbe-a2ee-05b52be35210
let 
	fig, ax = FigAxis(aspect=DataAspect())

	xi_p, eta_p = lguys.shear_points_to_ellipse(all_stars.xi, all_stars.eta, params.ellipticity, params.PA)
	scatter!(xi_p, eta_p, markersize=3, alpha=0.05, color=:black)

	fig
end

# ╔═╡ 6895c52d-9b93-48a3-b05f-87a6d66349b3
let 
	fig, ax = FigAxis(aspect=DataAspect())

	scatter!(all_stars.xi, all_stars.eta, markersize=3, alpha=0.05, color=:black)

	fig
end

# ╔═╡ a53f9db4-df89-4717-8535-c06c989307bd
filt_best = all_stars.F_BEST .== 1

# ╔═╡ 75fadb00-d768-4dcc-a619-f0f71de30263
names(all_stars)

# ╔═╡ 44a44f97-9115-4610-9706-33acf065d0e7
members = select_members(all_stars, params)

# ╔═╡ 0c498087-0184-4da2-a079-e972dd987712
md"""
The next three plots compare how different the r_ell and xi and eta calculated here are from what is (presumably) in the given catalogue.
"""

# ╔═╡ 52bb6b36-736a-45a8-b1e1-7f174b366ec8
let
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel="PSAT", 
		ylabel="number of stars",
		yscale=log10)

	hist!(ax, all_stars.PSAT[all_stars.PSAT .>= 0], 
		bins=20, label="j24")
	f
end

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ d7d81ed8-0427-4ee5-8213-320ce5a6711f
function plot_tangent!(grid, all_stars, members=nothing; markersize=2, kwargs...)
	ax = Axis(grid, 
		xlabel = "xi / degrees",
		ylabel = "eta / degrees",
		aspect=DataAspect(),
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
function plot_pms!(grid, all_stars, members=nothing; da=10)
	ax = Axis(grid, limits=(-da, da, -da, da), aspect=1,
	    xlabel=L"\mu_{\alpha *} / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}",
		xgridvisible=false,
		ygridvisible=false,
	)
	
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

# ╔═╡ 082dcd7f-6b4b-42a1-b9f5-a65380787cee
all_stars_unfiltered.F_ASTROMETRIC

# ╔═╡ 31d4e643-67eb-45bf-9c3b-1e031ed0fed8
sum(all_stars_unfiltered.F_ASTROMETRIC .!== 1.0 * (all_stars_unfiltered.ruwe .< 1.3))

# ╔═╡ 494befb5-007c-40a7-8a8d-3e9908a8fb19
sum(all_stars_unfiltered.F_NOTANAN .!== @. 1.0 * !(isnan(all_stars_unfiltered.pmra)
	# || isnan(all_stars_unfiltered.parallax)
	# || isnan(all_stars_unfiltered.pmdec)
	|| isnan(all_stars_unfiltered.phot_bp_mean_mag)
	|| isnan(all_stars_unfiltered.phot_rp_mean_mag)

))

# ╔═╡ 1219b4de-5e92-4951-8a21-a2b8a9e462af
sum(all_stars_unfiltered.F_INGRID .!==
@. 1.0 * !(
		abs(all_stars_unfiltered.pmra) > 10
		|| abs(all_stars_unfiltered.pmdec) > 10
		|| (all_stars_unfiltered.phot_g_mean_mag) > 21
		|| (all_stars_unfiltered.phot_g_mean_mag) < 16
		|| (all_stars_unfiltered.bp_rp > 2.5) 

)
)

# ╔═╡ e09b031c-a1d3-4e10-bf80-776f0cee69f4
sum(all_stars_unfiltered.F_CPAR .!==
	@. 1.0 * !(
		abs.(all_stars_unfiltered.parallax_over_error) > 3
	)
)

# ╔═╡ ab387f41-8d0f-4a1f-8dbf-b90396914c75
maximum(all_stars.parallax_over_error)

# ╔═╡ bbc4400c-66cc-42fb-95c3-f21c489ec316
sum(all_stars_unfiltered.F_BEST .!==
	@. (
		all_stars_unfiltered.F_ASTROMETRIC
		* all_stars_unfiltered.F_CPAR
		* all_stars_unfiltered.F_INGRID
		* all_stars_unfiltered.F_NOTANAN
		* all_stars_unfiltered.F_FLUXEXCESS
		#* (all_stars_unfiltered.phot_g_mean_mag < 21)
	))

# ╔═╡ 18fd8d15-f64d-47fa-88e5-a449ca62f156
ingrid = all_stars_unfiltered[all_stars_unfiltered.F_INGRID .== 1, :]

# ╔═╡ 4d2cc8cc-8d8b-4582-b80c-a69a72bac9c1
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="G")
	hist!(all_stars.phot_g_mean_mag[filt_best], bins=100)
	fig
end

# ╔═╡ c5ca7506-5952-4973-879e-e9848bb72d03
plot_sample(all_stars)

# ╔═╡ f74f597a-908b-4569-8de8-219abba31afd
filt_parallax = all_stars.F_CPAR .== 1; sum(filt_parallax)

# ╔═╡ eab8e0d4-0a40-48b6-9689-0c102d883a96
filt_pm = all_stars.PSAT_PM .> params.PSAT_min; sum(filt_pm)

# ╔═╡ c15d2045-f6c2-42d1-947e-fb4b31daeb99
filt_ruwe = all_stars.F_ASTROMETRIC .== 1; sum(filt_parallax)

# ╔═╡ b0e0a366-8bac-4cf6-8b40-33d517b41e47
filt_cmd = all_stars.PSAT_CMD .> 0.1

# ╔═╡ 7e1a766d-3ee1-45b9-a9d6-8bf1bd583454
all_stars.PSAT_CMD

# ╔═╡ 13809b36-b2b4-43df-a36c-f3f30b80accc
all_stars.PSAT_PM

# ╔═╡ 9fbef5d9-1530-4fe4-9190-6671fb9e81eb
filt_qso = (.! all_stars.in_qso_candidates) .&& (.! all_stars.in_galaxy_candidates)

# ╔═╡ 461e4520-5cce-4f45-b0f6-4819bb53b0af
sum(all_stars.in_qso_candidates)

# ╔═╡ 25c847ee-b0f0-4058-a0c0-2cfba2847732
sum(.! filt_qso)

# ╔═╡ 29dba8f1-b4a6-41da-8323-437447c9d888
filt = filt_cmd .& filt_pm .& filt_parallax .& filt_ruwe

# ╔═╡ 1a7953dc-2cd9-4d0b-9e06-5a476f48ac3b
filt_good = all_stars.F_BEST .== 1

# ╔═╡ e7847510-e9dc-4ca4-aa37-662d975b74dc
sum(filt .& filt_good)

# ╔═╡ 620aa7d7-0b02-475d-80a9-1e721d6144bf
mymemb = all_stars[filt, :]

# ╔═╡ 686bd76e-893c-461e-bb50-a9d12492d123
plot_sample(all_stars[filt, :], members, )

# ╔═╡ 39f615b4-b62c-493e-8d11-be45910d79a8
md"""
# density calculation
"""

# ╔═╡ 9b3288b4-3c17-4325-9e2c-94f96328f3c3
function plot_rh!()
    vline!([log10.(r_h)], color="grey", s=:dash, z_order=1, label=L"r_h")
end

# ╔═╡ 62388099-b80b-4272-9bde-0c7315b15c19
r = 60 * members.r_ell # arcminutes

# ╔═╡ a254083b-61f2-45ad-b678-c4df6f16964b
let
	fig = Figure()
	ax = PolarAxis(fig[1,1])
	ϕ = atan.(members.eta, members.xi)
	scatter!(ax, ϕ, members.r_ell, markersize=5, alpha=0.1)

	fig
end

# ╔═╡ b20058a8-a9b6-49ff-b8ff-a2d45c76f645
R_s_over_R_h = 1.6783

# ╔═╡ a7209445-84e9-435d-9144-90f8eb5e70cb
md"""
# Membership selection effects
"""

# ╔═╡ 13fb3ebc-50c0-43aa-88e9-1a7543e4e202
hist(members.phot_g_mean_mag)

# ╔═╡ 80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
size(members)

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

# ╔═╡ 1fc6160b-17d1-4d98-ad63-d4894c5816c1
lguys.histogram(randn(100), normalization=:none)

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

# ╔═╡ 52d80029-24fa-43da-9d43-b03babc65079
plot_tangent(all_stars, all_stars[.!filt_qso, :])


# ╔═╡ 6648b46f-b026-4998-8b27-c38ee6032b58
plot_cmd(all_stars, all_stars[.!filt_qso, :])


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

# ╔═╡ a2cf1532-f8c9-40cf-9e65-cc446ca6db34
plot_sample_w_dens(all_stars, 
	all_stars[filt .& filt_good, :]
)

# ╔═╡ e2dc5d97-0c8f-49dd-bb1a-320b45132cca
sum(all_stars.r_ell .< 4) / π /4^2

# ╔═╡ 62f00571-a1d0-4f77-9f09-a4b87c5aa63f
plot_sample_w_dens(all_stars, 
	all_stars[all_stars.PSAT .> 0.2, :]
)

# ╔═╡ d96709fc-8c70-4b08-91ae-0aeedf4ca54a
sum(.!filt_qso)

# ╔═╡ f85b8e62-5179-44e7-ba34-9253df297d91
plot_sample_w_dens(all_stars, 
	all_stars[.!filt_qso, :]
)

# ╔═╡ 9f25a484-5e75-431a-9720-680cafdbbf8e
plot_sample_w_dens(all_stars, 
	all_stars[(all_stars.PSAT .> 0.2) .& filt_qso, :]
)

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
	plot_density!(ax, all_stars[filt_good, :], label="good only")

	plot_density!(ax, all_stars[filt_parallax, :], label="parallax only")
	plot_density!(ax, all_stars[filt_cmd, :], label="CMD only")

	plot_density!(ax, all_stars[filt_pm, :], label="proper motion only")

	axislegend(position=:lb)

	fig
end

# ╔═╡ fd11c145-2743-446f-854e-0f45a344a404
maximum(all_stars_unfiltered.r_ell)

# ╔═╡ e01ed571-4661-4ce1-8013-ba91237fa856
maximum(all_stars.r_ell)

# ╔═╡ 691d6285-6bd3-47a1-a9e8-024874fb6bbd
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log radius / arcmin",
		ylabel=L"\log\;\Sigma \, / \, \textrm{stars arcmin^{-2}}",
		limits=(nothing, (-5, 3))
	)
	
	plot_density!(ax, gaia_all, label="no filter")

	plot_density!(ax, memb_simple, label="simple cuts")

	plot_density!(ax, all_stars[all_stars.PSAT_NOSPACE .> 0.2, :], label="PSAT_nospace > 0.2")
	plot_density!(ax, all_stars[all_stars.PSAT .> 0.2, :], label="PSAT > 0.2")

	axislegend()

	fig
end

# ╔═╡ 699fa15a-c88d-4f88-a31c-c2b2d0ff3cd3
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log radius / arcmin",
		ylabel=L"\log\;\Sigma \, / \, \textrm{stars arcmin^{-2}}",
		limits=(nothing, (-5, 3))
	)
	
	plot_density!(ax, all_stars[all_stars.PSAT .> 0.2, :], label="PSAT > 0.2")
	plot_density!(ax, all_stars[(all_stars.PSAT .> 0.2) .& filt_qso, :], label="+ not in QSO/Galaxy cand.")

	axislegend()

	fig
end

# ╔═╡ ba6ac07c-94f1-4024-8838-9c4d9b1c6cb1
let
	fig, ax = FigAxis(
		limits=(-1, 1, -1, 1)
	)

	errscatter!(members.pmra, members.pmdec, xerr=members.pmra_error, yerr=members.pmdec_error, alpha=0.1)
	fig
end

# ╔═╡ 6d69d5d8-2e5b-4bcb-aadf-fb0dac2a49bb
let 
	pm_disp = 0.02
	n_sig_pmra = (members.pmra .- params.pmra) ./ (members.pmra_error .+ pm_disp)
	n_sig_pmdec = (members.pmdec .- params.pmdec) ./ (members.pmdec_error .+ pm_disp)
	n_sig_pm = @. sqrt(n_sig_pmra^2 + n_sig_pmdec^2)

	fig, ax = FigAxis(xlabel="z-score proper motion consistent with UMi",
		limits=((0,6), nothing)
	)

	hist!(n_sig_pm)
	fig
end

# ╔═╡ 040dab05-b690-4e7b-839d-e65196b1856a
let 
	pm_disp = 0.02
	n_sig_pmra = (members.pmra ) ./ (members.pmra_error .+ pm_disp)
	n_sig_pmdec = (members.pmdec) ./ (members.pmdec_error .+ pm_disp)
	n_sig_pm = @. sqrt(n_sig_pmra^2 + n_sig_pmdec^2)

	fig, ax = FigAxis(xlabel="z-score proper motions consistent with 0", limits=((0,6), nothing))

	hist!(n_sig_pm)
	fig
end

# ╔═╡ Cell order:
# ╠═48caecb2-180c-4ce4-a57b-6fed82328b01
# ╟─47b0d3e6-a79b-4f49-b847-e708e5d6aabf
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═91d87251-c9c6-467f-9cae-4f30bfea8acc
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╟─8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╠═a86b5875-c30a-446a-8897-763b07f76c68
# ╠═fecf285c-1d7b-4e91-8108-8c906197cd20
# ╟─4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═9371bcb3-0e26-4162-8a40-dc1bf1dacdda
# ╠═0307085b-816f-469f-8811-60af02cfcb67
# ╠═29859854-0d17-42ba-8b1d-8788511840c9
# ╠═c26a1042-1c5f-4539-992f-022cc7441a1a
# ╠═9c2da11e-79fb-406e-9c5a-83eefb620814
# ╠═6c092147-c295-4ee5-9ee3-6e04c2aaaf98
# ╠═c403ee08-852b-4bbe-a2ee-05b52be35210
# ╠═6895c52d-9b93-48a3-b05f-87a6d66349b3
# ╠═a53f9db4-df89-4717-8535-c06c989307bd
# ╠═75fadb00-d768-4dcc-a619-f0f71de30263
# ╠═44a44f97-9115-4610-9706-33acf065d0e7
# ╟─0c498087-0184-4da2-a079-e972dd987712
# ╟─52bb6b36-736a-45a8-b1e1-7f174b366ec8
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
# ╠═082dcd7f-6b4b-42a1-b9f5-a65380787cee
# ╠═31d4e643-67eb-45bf-9c3b-1e031ed0fed8
# ╠═494befb5-007c-40a7-8a8d-3e9908a8fb19
# ╠═1219b4de-5e92-4951-8a21-a2b8a9e462af
# ╠═e09b031c-a1d3-4e10-bf80-776f0cee69f4
# ╠═ab387f41-8d0f-4a1f-8dbf-b90396914c75
# ╠═bbc4400c-66cc-42fb-95c3-f21c489ec316
# ╠═18fd8d15-f64d-47fa-88e5-a449ca62f156
# ╠═4d2cc8cc-8d8b-4582-b80c-a69a72bac9c1
# ╠═c5ca7506-5952-4973-879e-e9848bb72d03
# ╠═f74f597a-908b-4569-8de8-219abba31afd
# ╠═eab8e0d4-0a40-48b6-9689-0c102d883a96
# ╠═c15d2045-f6c2-42d1-947e-fb4b31daeb99
# ╠═b0e0a366-8bac-4cf6-8b40-33d517b41e47
# ╠═7e1a766d-3ee1-45b9-a9d6-8bf1bd583454
# ╠═13809b36-b2b4-43df-a36c-f3f30b80accc
# ╠═9fbef5d9-1530-4fe4-9190-6671fb9e81eb
# ╠═461e4520-5cce-4f45-b0f6-4819bb53b0af
# ╠═25c847ee-b0f0-4058-a0c0-2cfba2847732
# ╠═29dba8f1-b4a6-41da-8323-437447c9d888
# ╠═1a7953dc-2cd9-4d0b-9e06-5a476f48ac3b
# ╠═e7847510-e9dc-4ca4-aa37-662d975b74dc
# ╠═620aa7d7-0b02-475d-80a9-1e721d6144bf
# ╠═686bd76e-893c-461e-bb50-a9d12492d123
# ╟─39f615b4-b62c-493e-8d11-be45910d79a8
# ╠═9b3288b4-3c17-4325-9e2c-94f96328f3c3
# ╠═62388099-b80b-4272-9bde-0c7315b15c19
# ╠═a254083b-61f2-45ad-b678-c4df6f16964b
# ╠═b20058a8-a9b6-49ff-b8ff-a2d45c76f645
# ╟─a7209445-84e9-435d-9144-90f8eb5e70cb
# ╠═13fb3ebc-50c0-43aa-88e9-1a7543e4e202
# ╠═80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
# ╟─49ae0572-5d6b-4935-bc95-0a845bb3df2f
# ╠═65a10161-cbeb-49fe-b8d1-075ffe346e43
# ╠═08b5251d-bbe4-48f6-9cb3-01e0a8364c1d
# ╠═1fc6160b-17d1-4d98-ad63-d4894c5816c1
# ╠═d0c00f3e-f3d3-4cd1-8541-9ae239420174
# ╠═198ba5a6-9d04-49de-9726-4e1d9be9780f
# ╠═b55d72b9-e957-4480-b6db-97c9798b4d68
# ╠═bd4212fe-3b70-44ba-85fb-b199840a6f71
# ╠═52d80029-24fa-43da-9d43-b03babc65079
# ╠═6648b46f-b026-4998-8b27-c38ee6032b58
# ╠═1c06e13b-8f5d-414b-8f2a-3d6f051ad495
# ╠═289b5daa-a0a3-41c8-9d9c-3a3280139f2a
# ╠═498399c5-43f7-44db-bddd-03ec0f0fdd6b
# ╠═fab13f99-662c-4385-b266-029ce2eea233
# ╠═aac80d15-72c9-4091-befa-1b7fa7127c63
# ╠═4bfc9046-a5b2-4a72-82a7-3547066d7064
# ╠═d74e062b-d7b0-422e-9f1b-102afee26d7b
# ╠═16bd1702-b471-44da-96df-aca63c46a604
# ╠═a2cf1532-f8c9-40cf-9e65-cc446ca6db34
# ╠═e2dc5d97-0c8f-49dd-bb1a-320b45132cca
# ╠═62f00571-a1d0-4f77-9f09-a4b87c5aa63f
# ╠═d96709fc-8c70-4b08-91ae-0aeedf4ca54a
# ╠═f85b8e62-5179-44e7-ba34-9253df297d91
# ╠═9f25a484-5e75-431a-9720-680cafdbbf8e
# ╠═b76fbfe9-8a0d-4816-ae81-6fb8597aaf80
# ╠═81de6949-c021-49ea-82f5-b300c03a2cc0
# ╠═12220e1c-9a98-447f-b3a4-b579b2989b68
# ╠═fd11c145-2743-446f-854e-0f45a344a404
# ╠═e01ed571-4661-4ce1-8013-ba91237fa856
# ╠═691d6285-6bd3-47a1-a9e8-024874fb6bbd
# ╠═699fa15a-c88d-4f88-a31c-c2b2d0ff3cd3
# ╠═ba6ac07c-94f1-4024-8838-9c4d9b1c6cb1
# ╠═6d69d5d8-2e5b-4bcb-aadf-fb0dac2a49bb
# ╠═040dab05-b690-4e7b-839d-e65196b1856a
