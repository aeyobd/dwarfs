### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ a893932c-f184-42bc-9a0e-0960f10520aa
begin
	import Pkg
	Pkg.activate()
	using CairoMakie
	using Arya
	
	using LilGuys
end

# ╔═╡ 17ffde4b-5796-4915-9741-d594cf0c5ca7
md"""
# Distribution functions

This notebook is a scratch notebook which calculates the distribution function analytically from a potential and calculates some moments of the distribution function.

"""

# ╔═╡ 93045024-a91d-4b31-9a5a-7c999afdb9ec
md"""
# Inputs
"""

# ╔═╡ 73bb13e4-4a36-42b3-a04d-f235a8b11362
R_s_s = 0.10

# ╔═╡ 64578172-da38-49ed-8777-51b538aa9b18
prof_halo = LilGuys.TruncNFW(v_circ_max=31/V2KMS, r_circ_max=3.2, trunc=100)

# ╔═╡ 5a42e848-caee-4aea-a541-7812c4fbeb13
profile = LilGuys.Exp2D(R_s = R_s_s)

# ╔═╡ aa69bde5-ab93-4105-9d48-ad0ace88e3f0
r_h = LilGuys.calc_r_h(profile)

# ╔═╡ 1b9f3101-4891-4a8c-8c73-41a222b6c26f
ρ_s(r) = calc_ρ(profile, r)

# ╔═╡ e935a80a-7f4f-4143-ae6c-bee62b82c30e
md"""
# Core calculation
"""

# ╔═╡ 23158c79-d35c-411a-a35a-950f04214e19
begin 
	M_s_tot = calc_M(profile, Inf)
	M_s(r) = calc_M(profile, r)
end

# ╔═╡ 29619cc3-1be3-4b24-92e0-ceccfd4a3f59
"""
given the radii and parameters for stellar profile, returns the
radial bins used in the calculation. 

"""
function make_radius_bins(radii::AbstractVector, params::Dict)
	log_radii = log10.(radii)

	
	lr_min = minimum(log_radii)
	lr_max = maximum(log_radii)
	Nr = params["num_radial_bins"]
	
	r_e = 10 .^ LinRange(lr_min, lr_max, Nr+1)

	return r_e
	
end

# ╔═╡ deb46b0b-3504-4231-9f85-99f27caf924b
r_e = 10 .^ LinRange(-5, 4, 10000)

# ╔═╡ 36b4adbd-d706-4e72-a922-53080c67946c
r = midpoints(r_e)

# ╔═╡ f3e95cfc-0087-4afa-8e88-a4e1628c50a0
"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(r_bins, M_s)
	N_s_out = 1 - M_s(r_e[end])
	N_s_in = M_s(r_e[1])
	println("missing $N_s_out stars outside")
	println("missing $N_s_in stars inside")
end

# ╔═╡ 41a08be6-d822-4a19-afe5-62c7dc9ff118
print_missing(r_e, M_s)

# ╔═╡ dfa675d6-aa32-45c3-a16c-626e16f36083
ψ = - LilGuys.calc_Φ.(prof_halo, r)

# ╔═╡ ab4d4458-fd1d-417a-8a74-5e06f41af166
ν_dm = calc_ρ.(prof_halo, r)

# ╔═╡ 20f858d6-a9f5-4880-a431-60b395cc7e50
ν_s = max.(ρ_s.(r) ./ M_s_tot, 0)

# ╔═╡ cf2f2756-5775-4302-b0eb-50ed530cec1a
f_dm = LilGuys.DistributionFunction(ν_dm, ψ, r)

# ╔═╡ 1fab14ec-6bfc-4243-bc48-915d1a129925
f_s = LilGuys.DistributionFunction(ν_s, ψ, r)

# ╔═╡ 126c6825-723f-4d13-b5a3-64cba72fc867
md"""
The density distribution of the stars (analytic) and dark matter (calculated) from the snapsho
"""

# ╔═╡ 7481f47a-2b8a-45a3-9b4e-31ea14e9d331
md"""
The potential $\psi = -\Phi$ as a function of log radii (for the spherically calculated & interpolated and actual snapshot from Gadget 4)
"""

# ╔═╡ 1c0899c6-2692-45ab-b9e6-668dc576d679
function make_energy_bins(ψ, N; k=5)
	E_max = ψ[1]
	E_min = ψ[end]

	x = LinRange(-k, k, N+1)
	E = E_min .+ (1 .+ tanh.(x)) ./ 2 * (E_max - E_min)

	E
end

# ╔═╡ 3025546e-3ce1-4a78-824d-24f644238e32
E = ψ

# ╔═╡ 500c67b2-8c4a-4d03-b4e4-39beff43a46c
f_dm_e = f_dm.(E)

# ╔═╡ 9e492a55-7b20-4eca-aead-c7afeee63f11
f_s_e = f_s.(E)

# ╔═╡ e9a12600-bf3c-43a6-8ab6-a06b463229e8
function calc_g(ϵ, radii=r)
	integrand(r) = ((-LilGuys.calc_Φ(prof_halo, r) - ϵ) > 0) && sqrt( 2 * (-LilGuys.calc_Φ(prof_halo, r) - ϵ)) * r^2
	
	return 16π^2 * sum(integrand.(radii) .* LilGuys.gradient(radii))
end

# ╔═╡ 5992781a-8041-405c-ae3d-c60b1f8ef951
r_m = LilGuys.lerp(ψ, r)

# ╔═╡ 96554ebf-bf80-4fac-b9bc-3a7d6c524d6d
function calc_g_int(ϵ)
	integrand(r) = ((-LilGuys.calc_Φ(prof_halo, r) - ϵ) > 0) && sqrt( 2 * (-LilGuys.calc_Φ(prof_halo, r) - ϵ)) * r^2
	
	return 16π^2 * LilGuys.integrate(integrand, 0, r_m(ϵ))
end

# ╔═╡ 94c9c9cb-762d-44a9-b0c5-3634dffa2694
calc_g(E[50], r)

# ╔═╡ 2e0cf0f5-0ada-43be-9045-207bf10e0a97
calc_g_int(E[50])

# ╔═╡ 8b66d00d-529b-4e8c-9386-b17117996579
md"""
The calculated distribution function as a function of log specific binding energy $\epsilon =  -\Phi - T$
"""

# ╔═╡ 9f55e818-0ab8-4c25-987b-386ebb7f8e31
f_dm_i = LilGuys.lerp(E, f_dm_e)

# ╔═╡ 81e027df-311f-4036-a8ec-5e03ddbfc20d
f_s_i = LilGuys.lerp(E, f_s_e)

# ╔═╡ 7409a024-cea4-47a5-84d2-846c96d88b7a
begin 
	probs = f_s_e ./ f_dm_e
	probs ./= sum(probs .* LilGuys.gradient(E)) # pdf, dN/dE
	prob = LilGuys.lerp(E, probs)
end

# ╔═╡ 8bb8736d-a41b-4dac-a6cd-06d0d4704654
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(nothing, nothing, -15, 3),
		xlabel=L"\log\ r", 
		ylabel=L"\log\ \rho"
	)
	lines!(log10.(r), log10.(ν_dm), label="DM")

	lines!(log10.(r), log10.(ν_s), label="stars")


	
	axislegend()
	fig
end

# ╔═╡ b625d8a5-7265-4849-9bd6-ca8064d392eb
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log\ r / \textrm{kpc}", ylabel=L"\Psi(r)")
	lines!(log10.(r), ψ, label="theoretical")

	fig
end

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="ϵ [code units]", ylabel="log f(ϵ)", limits=(nothing, (-15, 10)) )
	lines!((E), log10.(f_s_e), label="stars")
	lines!((E), log10.(f_dm_e), label="DM")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 29a07d14-5de3-44ee-894e-4a5e6969f2dd
g = calc_g.(E)

# ╔═╡ 8b810ccf-5ff3-405a-be99-82e753dad894
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="ϵ [code units]", ylabel="log dN/dϵ", limits=(nothing, (-15, 3)) )
	lines!((E), log10.(f_s_e .* g), label="stars")
	lines!((E), log10.(f_dm_e .* g), label="DM")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 652ab3a4-8ac4-44df-be5b-555447c63522
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel="ϵ [code units]", 
		ylabel="log probability(ϵ)", 
		limits=(nothing, (-15, 4)) 
	)
	lines!(E, log10.(f_s_e ./ f_dm_e), label="stars")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ a5ccf909-79b7-4ac2-b3d0-0b7796e14354
md"""
## New, calculating the velocity dispersion
"""

# ╔═╡ c10d65b8-daa6-4a1f-b07b-06c507dc6a83
E

# ╔═╡ 0d217f82-8316-40b8-a344-c4a409e9c70d
function phase_density(r, vel, f)
	V = LilGuys.calc_Φ(prof_halo, r)
	T = @. 1/2 * vel^2
	E = -T .- V

	if E < 0
		return 0
	end

	return f(E)
end

# ╔═╡ df6e75a8-a6ca-409c-b2b1-25a2be74af45
begin 
	xs = 10 .^ LinRange(log10.(r_e[1]), log10(r_e[end]), 500)
	vs = LinRange(0.00, (2ψ[1])^0.5, 500)
	dx = LilGuys.gradient(xs)
	dv = LilGuys.gradient(vs)
	
	Z = hcat([[phase_density(x, v, f_s_i) for x in xs] for v in vs]...)
	Z_dm = hcat([[phase_density(x, v, f_dm_i) for x in xs] for v in vs]...)

end

# ╔═╡ aceb7c70-de8d-424c-bcda-20d7874d21f2
let
	fig, ax = FigAxis(
		xlabel=L"$\log\ |r\,|$ / kpc",
		ylabel = L"$|v\,|$ / km s$^{-1}$",
		limits=(nothing, nothing)
	)	
	
	h = heatmap!(log10.(xs), V2KMS * vs, Z,
		colorscale=log10, colorrange=(1e-12, maximum(Z))
	)
	
	y = sqrt.( - 2LilGuys.calc_Φ.(prof_halo, xs)) .* V2KMS
	lines!(log10.(xs), y)

	Colorbar(fig[1, 2], h, label="stellar distribution function")

	fig
end

# ╔═╡ f31d102e-946e-404f-9d2f-30b48aad4b04
vs'

# ╔═╡ 7a85cdc2-0001-4d87-817a-9c71a2df9971
let
	fig, ax = FigAxis(
		xlabel=L"$\log\ |r\,|$ / kpc",
		ylabel = L"$|v\,|$ / km s$^{-1}$",
		limits=(nothing, nothing)
	)
	
	h = heatmap!(log10.(xs), V2KMS * vs, Z_dm,
	colorscale=log10, colorrange=( 1e-15 * maximum(Z), maximum(Z))
	)
	y = sqrt.( - 2LilGuys.calc_Φ.(prof_halo, xs)) .* V2KMS
	lines!(log10.(xs), y)

	Colorbar(fig[1, 2], h, label="DM distribution function")

	fig
end

# ╔═╡ 81f2a43e-65b5-4578-858f-76c33088e65a
sum(
	@. Z * 16π^2 * xs^2 * vs'^2 * dx * dv'
)

# ╔═╡ 2dbc5266-68ff-4d4b-83dd-f144d5554f6b
sum(
	@. Z_dm * 16π^2 * xs^2 * vs'^2 * dx * dv'
)

# ╔═╡ 35db741c-5f30-476e-9a65-5d86313ad290
LilGuys.calc_M(prof_halo, xs[end])

# ╔═╡ 5be4f410-e39c-429b-bef1-2188ccb9e556
md"""
# Moments of the distribution function
"""

# ╔═╡ 68d0d571-99f1-4c39-bf17-7c481686cfaf
md"""
Formally, f(...) is the phase space density, ie number of particles within some box in dx, dy, dz, dvx, dvy, dvz. So, many more directly understandable quantities arise from moments of the distribution function.

First of all, the integral of the distribution function should be equal to the system mass,

``
M = \iint f d{\bf v}^3 d{\bf x}^3
``

Similarly, the density is given by

``
\rho = \int f d{\bf v}^3
``

.
The velocity dispersion is the integral

``
\sigma_v^2 = \iint v^2 f d{\bf x}^3 d{\bf v}^3
``

"""

# ╔═╡ 6c55bde8-df17-4c9c-85f7-119fd1a1c290
nu_calc = dropdims(sum(
		@. Z * (vs')^2 * dv' * 4π
	, dims=2), dims=2) 

# ╔═╡ eae166b1-7884-46b1-a863-76ee5c3817e1
sum(nu_calc .* dx .* xs .^ 2 .* 4π)

# ╔═╡ b32d21fe-5aa5-4fab-b716-b7f9c8dc3b1e
nu_dm = dropdims(sum(
		@. Z_dm * (vs')^2 * dv' * 4π
	, dims=2), dims=2) 

# ╔═╡ cc8a9d43-9613-4b1f-bf00-b3118efb0e94
let
	fig, ax = FigAxis(
		limits=(nothing, (-15, 5)),
		xlabel="log r / kpc",
		ylabel = L"$\rho_\star$ / M$_\odot$r$^{-3}$",
	)
	
	scatter!(log10.(xs), log10.(nu_calc), 
		label="integrated")

	y_model = LilGuys.calc_ρ.(profile, xs)

	lines!(log10.(xs), log10.(y_model), label="analytic")

	axislegend()
	fig
end

# ╔═╡ 6da67070-1cf6-4ffa-b62d-bc87677056cc
let
	fig, ax = FigAxis(
		limits=(nothing, (-5, 14)),
		xlabel="log r / kpc",
		ylabel = L"$\rho_\textrm{DM}$ / M$_\odot$r$^{-3}$",
	)
	
	y_model = LilGuys.calc_ρ.(prof_halo, xs) * M2MSUN
	lines!(log10.(xs), log10.(y_model), label="analytic")

	scatter!(log10.(xs), log10.(nu_dm * M2MSUN), 
		markersize=3, color=COLORS[2], label="integrated")

	axislegend()


	fig
end

# ╔═╡ 94a09909-2af8-481f-a46c-d5a0dc8cbe50
v_dens = dropdims(sum(
		@. Z * 4π * xs^2 * dx
	, dims=1), dims=1) 

# ╔═╡ 72348a89-ad25-4e75-bb7e-11623d23240b
v_dens_dm = dropdims(sum(
		@. Z_dm * 4π * xs^2 * dx
	, dims=1), dims=1) 

# ╔═╡ 4b88d424-bf0d-45eb-a6ae-594b1ffd61b6
σ_v = sqrt(
	sum(
	@. v_dens * vs^2 * 4π * vs^2 * dv
)
)

# ╔═╡ 93d17b6b-f909-41d9-9547-31a4ccdcc3e8
σ_v_dm = sqrt(
	sum(
	@. v_dens_dm * vs^2 * 4π * vs^2 * dv
)
)

# ╔═╡ 22291d60-6f67-4960-a3f9-ca3cc2220674
σ_v * V2KMS

# ╔═╡ e7ee66a4-f783-4f72-b6d2-02dfec7f6977
σ_v_dm * V2KMS

# ╔═╡ c9a5a155-d588-4d74-9bbc-fd8103eecc8b
σ_v * V2KMS / sqrt(3) # average x-velocity dispersion

# ╔═╡ 04d5ec93-ec4c-4c83-bf03-5fb7456bb0f4
let
	fig, ax = FigAxis(
		xlabel="|v| / km / s",
		ylabel = "density / dv^3",
		limits=(nothing, (-5, 3))
	)

	lines!(vs * V2KMS, log10.(v_dens / V2KMS))
	fig
end

# ╔═╡ Cell order:
# ╟─17ffde4b-5796-4915-9741-d594cf0c5ca7
# ╠═a893932c-f184-42bc-9a0e-0960f10520aa
# ╟─93045024-a91d-4b31-9a5a-7c999afdb9ec
# ╠═73bb13e4-4a36-42b3-a04d-f235a8b11362
# ╠═64578172-da38-49ed-8777-51b538aa9b18
# ╠═5a42e848-caee-4aea-a541-7812c4fbeb13
# ╠═aa69bde5-ab93-4105-9d48-ad0ace88e3f0
# ╠═1b9f3101-4891-4a8c-8c73-41a222b6c26f
# ╟─e935a80a-7f4f-4143-ae6c-bee62b82c30e
# ╠═23158c79-d35c-411a-a35a-950f04214e19
# ╟─29619cc3-1be3-4b24-92e0-ceccfd4a3f59
# ╠═deb46b0b-3504-4231-9f85-99f27caf924b
# ╠═36b4adbd-d706-4e72-a922-53080c67946c
# ╟─f3e95cfc-0087-4afa-8e88-a4e1628c50a0
# ╠═41a08be6-d822-4a19-afe5-62c7dc9ff118
# ╠═dfa675d6-aa32-45c3-a16c-626e16f36083
# ╠═ab4d4458-fd1d-417a-8a74-5e06f41af166
# ╠═20f858d6-a9f5-4880-a431-60b395cc7e50
# ╠═cf2f2756-5775-4302-b0eb-50ed530cec1a
# ╠═1fab14ec-6bfc-4243-bc48-915d1a129925
# ╟─126c6825-723f-4d13-b5a3-64cba72fc867
# ╟─7481f47a-2b8a-45a3-9b4e-31ea14e9d331
# ╠═1c0899c6-2692-45ab-b9e6-668dc576d679
# ╠═3025546e-3ce1-4a78-824d-24f644238e32
# ╠═500c67b2-8c4a-4d03-b4e4-39beff43a46c
# ╠═9e492a55-7b20-4eca-aead-c7afeee63f11
# ╠═96554ebf-bf80-4fac-b9bc-3a7d6c524d6d
# ╠═e9a12600-bf3c-43a6-8ab6-a06b463229e8
# ╠═5992781a-8041-405c-ae3d-c60b1f8ef951
# ╠═94c9c9cb-762d-44a9-b0c5-3634dffa2694
# ╠═2e0cf0f5-0ada-43be-9045-207bf10e0a97
# ╟─8b66d00d-529b-4e8c-9386-b17117996579
# ╠═9f55e818-0ab8-4c25-987b-386ebb7f8e31
# ╠═81e027df-311f-4036-a8ec-5e03ddbfc20d
# ╠═7409a024-cea4-47a5-84d2-846c96d88b7a
# ╠═8bb8736d-a41b-4dac-a6cd-06d0d4704654
# ╠═b625d8a5-7265-4849-9bd6-ca8064d392eb
# ╠═75d23b44-71e7-4e28-ad3e-c537f3d4422f
# ╠═29a07d14-5de3-44ee-894e-4a5e6969f2dd
# ╠═8b810ccf-5ff3-405a-be99-82e753dad894
# ╠═652ab3a4-8ac4-44df-be5b-555447c63522
# ╟─a5ccf909-79b7-4ac2-b3d0-0b7796e14354
# ╠═c10d65b8-daa6-4a1f-b07b-06c507dc6a83
# ╠═0d217f82-8316-40b8-a344-c4a409e9c70d
# ╠═df6e75a8-a6ca-409c-b2b1-25a2be74af45
# ╟─aceb7c70-de8d-424c-bcda-20d7874d21f2
# ╠═f31d102e-946e-404f-9d2f-30b48aad4b04
# ╟─7a85cdc2-0001-4d87-817a-9c71a2df9971
# ╠═81f2a43e-65b5-4578-858f-76c33088e65a
# ╠═2dbc5266-68ff-4d4b-83dd-f144d5554f6b
# ╠═35db741c-5f30-476e-9a65-5d86313ad290
# ╠═eae166b1-7884-46b1-a863-76ee5c3817e1
# ╟─5be4f410-e39c-429b-bef1-2188ccb9e556
# ╟─68d0d571-99f1-4c39-bf17-7c481686cfaf
# ╠═6c55bde8-df17-4c9c-85f7-119fd1a1c290
# ╠═b32d21fe-5aa5-4fab-b716-b7f9c8dc3b1e
# ╟─cc8a9d43-9613-4b1f-bf00-b3118efb0e94
# ╟─6da67070-1cf6-4ffa-b62d-bc87677056cc
# ╠═94a09909-2af8-481f-a46c-d5a0dc8cbe50
# ╠═72348a89-ad25-4e75-bb7e-11623d23240b
# ╠═4b88d424-bf0d-45eb-a6ae-594b1ffd61b6
# ╠═93d17b6b-f909-41d9-9547-31a4ccdcc3e8
# ╠═22291d60-6f67-4960-a3f9-ca3cc2220674
# ╠═e7ee66a4-f783-4f72-b6d2-02dfec7f6977
# ╠═c9a5a155-d588-4d74-9bbc-fd8103eecc8b
# ╠═04d5ec93-ec4c-4c83-bf03-5fb7456bb0f4
