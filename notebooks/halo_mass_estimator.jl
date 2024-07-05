### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 2c90aad4-87e7-11ee-2fbd-05f67c7e0343
begin 
	using Pkg; Pkg.activate()

	import LilGuys as lguys
	
	using CairoMakie, Roots, Arya
end

# ╔═╡ 07a710d8-0ae4-4d9f-9759-002750730010
md"""
# Estimation of halo mass
The first step of running any n-body simulation is to create an initial halo to realize. However our only direct observational constraints are 1. the stellar mass, and 2. the velocity dispersion. So, motivated by cosmological simulations, we use a relationship between the stellar mass and the total halo properties to estimate the initial NFW profile we need to match observations.

"""

# ╔═╡ 33f433ec-e94e-41fd-a808-6254bf4d34ce
md"""
## Procedure overview

First, we calculate the total stellar mass of the halo. 

Next, we use the emperical fit of the maximum circular velocity of a halo to the stellar mass from @fattahi2018; 
``
M_\star = m_0\,\nu^\alpha \exp({-\nu^\gamma})
``
where $m_0 = 3\times10^8 M_\odot$, $\alpha=3.36$, $\gamma=-2.4$, and $\nu = V_{\rm max}/$50km/s. 

Finally, the characteristic radius is from @ludlow2016 (see their appendix C). The full equations are complex, but briefly, from a given M and z, you can calculate their fit to the expected concentration. We then need to solve for a value of $M_s$ which gives a value of $c$ which together predict the correct $V_{\rm circ, max}$. 

Using that the maximum circular velocity occurs at $\alpha r_s$ where $\alpha = 2.163$ (is a numerical solution), and the equation for circular velocity of a NFW profile
``
\left(V_c/V_{200}\right)^2 = \frac{c}{r/r_s}\frac{A(r/r_s)}{A(c)}
``

where $V_{200} = \sqrt{G\ M_{200} / r_{200}}$. 

See my notes in `dwarfs/notes/gravitational_potentials.md` for more information on the equations related to an NFW halo
"""

# ╔═╡ 45ac95e4-a9a1-4f8f-9523-8f670f076ae5
md"""
## Example

Asya's crater II paper provides one example of this procedure.

1. Crater II has an absolute magnitude of $M_V = -11.1$ [@mcconnachie2012], and the sun has an absolute magnitude of 4.83, so with a stellar mass to light ratio of 1.6, this gives a total stellar mass of $M_\star=3.8\times10^6\,$M$_\odot$ .
2. Solving for circular velocity, we find $\nu=0.643$ or $V_{\rm max}=32.2$km/s.
3. The solution here is $M_{200} = 5.4e9$ where $c=13.0$, and $r_s = 2.86$ and $M_s = 6.12e8$. 


Or for Fornax, Mstr = 2.39e7, leaving Vmax = 39.6km / s, r_max = 8
"""

# ╔═╡ 67e892f1-3969-484f-9671-e0ac018babf2
md"""
## Input: stellar mass
This is done from the observed absolute magnitude and an expected mass to light for stars ratio. 
"""

# ╔═╡ fc490d05-d494-4771-9dec-7d1837485ff0
""" magnitude ratio to relative luminosity"""
mag_to_L(m1, m2) = 10^((m2-m1)/2.5)

# ╔═╡ e338d0b5-fd4f-43ed-9379-27d50787d387
MV_sol = 4.83 # solar magnitude

# ╔═╡ aec3d9a7-d447-4e9f-aa83-fa5b20541b5c
begin
	# mass to lig
	mass_to_light_stars = 1.7
	MV = -11.1 
	
	M_s_0 = mag_to_L(MV, MV_sol) * mass_to_light_stars / 1e10
	# Y_s = 1.6
	# M_s_0 = 2.56e5 # crater ii
	# M_s_0 = 2.39e-3 # fornax
end

# ╔═╡ bdd43353-da4c-48fb-bb67-b2aec628dd71
md"""
## Solving for Vcirc_max
"""

# ╔═╡ ca65b2b1-8414-48dc-badb-09c2d50879ad
md"""
Fattahi2018 analyzes a sample of dwarfs from apostle simulations. Their fit is described in the caption of Figure 1 and is without uncertainties. The fit is to the maximum values of circular velocity and stellar mass during the evolution of the satalite. There is maybe a 20% uncertainty in the circular velocity at a given stellar mass
"""

# ╔═╡ 18a20162-1b0f-4da8-931d-5ff95da22f54
"""
Given the maximum circular velocity (in km/s), returns the stellar mass (in Msun).
Uses the fit from @fattahi2018.
"""
function M_s_from_vel(v_max)
	ν = v_max * lguys.V0 / 50 # km/s
	α = 3.36
	γ = -2.4
	m_0 = 3e-2 # 1e10 Msun
	return m_0 * ν^α * exp(-ν^γ)
end

# ╔═╡ 3fc2fab3-2e37-4927-94f0-eeecaa89d575
V_circ_max_0 = find_zero(x -> M_s_from_vel(x) - M_s_0, 30)

# ╔═╡ d9cb3699-6a0f-4d32-b3f8-410cf2c92019
let
	v_1 = LinRange(18, 100, 1000) ./ lguys.V0

	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "v_circ_max / km s",
		ylabel = "M_s",
		yscale=log10,
		xscale=log10,
		limits=((1, 100), (1e-8, 1e0)),
		xticks=[1, 10, 100],
		xminorticks=IntervalsBetween(9),
		aspect=1
	)
	lines!(v_1 * lguys.V0, M_s_from_vel.(v_1))

	hlines!([M_s_0])
	vlines!([V_circ_max_0 * lguys.V0])

	fig

end

# ╔═╡ fd225a65-b833-4971-a308-09adc0f0b10f
md"""
## Concentration from mass & redshift
"""

# ╔═╡ 41b442e9-6691-4389-9517-7b0b55a25118
md"""
Ludlow2016 derive an analytic relationship between the halo mass and concentration with redshift. The y also compare results from cosmological simulations. The cosmological simulations have an error of 0.1 dex (log c), see their Figure 10.
"""

# ╔═╡ ef4ad7fc-fdf8-4e60-b2a6-e583ae226749
begin
	# from plank cosmology 2018
	h = 0.674 #pm 0.05
	Ω_m0 = 0.315
	Ω_Λ0 = 1 - Ω_m0
	σ8 = 0.811

	n_s = 0.965
	
	ρ_crit = 277.5366*h^2 # M⊙/kpc = (3H^2 / 8πG)

	const G = 4.30091e-6 # km^2/s^2 kpc / M⊙

end

# ╔═╡ 92def240-121d-467e-8026-f30e99e97b06
begin 
	# intermediate equations from @ludlow2016
	
	Ω_Λ(z) = Ω_Λ0 / (Ω_Λ0 + Ω_m0 * (1+z)^3)
	Ω_m(z) = 1 - Ω_Λ(z)

	Ψ(z) = Ω_m(z)^(4/7) - Ω_Λ(z) + (1 + Ω_m(z)/2) * (1 + Ω_Λ(z)/70)
	D(z) = Ω_m(z) / Ω_m0 * Ψ(0) / Ψ(z) * (1+z)^-1

	ξ(M) = 1/(h * M/1e10) # with M in M⊙
	σ(M, z) = D(z) * 22.26*ξ(M)^0.292 / (1 + 1.53*ξ(M)^0.275 + 3.36*ξ(M)^0.198)
end

# ╔═╡ c185ed67-e042-422f-9582-6a815c292fd2
begin
	# more intermediate equations from @ludlow2016

	c_0(z) = 3.395 * (1+z)^-0.215
	β(z) = 0.307 * (1+z)^0.540
	γ_1(z) = 0.628  * (1+z)^-0.047
	γ_2(z) = 0.317 * (1+z)^-0.893
	a(z) = 1/(1+z)
	ν_0(z) = 1/D(z) * (4.135 - 0.564/a(z) - 0.210/a(z)^2 + 0.0557/a(z)^3 - 0.00348/a(z)^4)
	δ_sc = 1.686
	ν(M, z) = δ_sc / σ(M, z)
end

# ╔═╡ 57071589-7f70-4739-bb6f-5ce184a8a09a
"""
Calculates the approximate concentration of a halo given the initial mass and redshift using the n-body fit from @ludlow2016
"""
function calc_c(M, z)
	x = ν(M * lguys.M0, z)/ν_0(z)
	result = c_0(z)
	result *= x^(-γ_1(z))
	result *= (1 + x^(1/β(z)))^(-β(z) * (γ_2(z) - γ_1(z)))
	return result
end

# ╔═╡ effb14d6-56ad-4baa-a484-83820aa72bb4
let
	M = 10 .^ LinRange(-8, 17, 1000) / 1e10
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log M_{200}\ / \ M_\odot",
		ylabel=L"\log c",
	)

	for z in [0, 1, 2, 3, 5, 7]
		lines!(log10.(M / h * lguys.M0), log10.(calc_c.(M, z)), label="z=$z", color=z, colorrange=(0, 8), colormap=:blues)
	end

	axislegend()
	fig
end

# ╔═╡ 312e16fd-8a35-4067-bd02-3260c213a9eb
md"""
# Solving for virial mass
"""

# ╔═╡ 56d46bf2-a4a3-4d93-97d7-44c63a35b1db
function NFW(M200::Real; z::Real=0)
	c = calc_c(M200, z)

	M_s = M200 / lguys.A_NFW(c)
	R200 = lguys.calc_R200(M200)
	r_s = R200 / c

	halo = lguys.NFW(M_s=M_s, r_s=r_s)
	
	if abs(halo.c - c) > 1e-6
		error("halo does not agree with c: $c != $(halo.c)")
	end
	return halo
end

# ╔═╡ a4c594f9-3de9-43d8-810f-b0f5d4cb693a
function NFW_from_M200_c(M200::Real, c::Real)
	M_s = M200 / lguys.A_NFW(c)
	R200 = lguys.calc_R200(M200)
	r_s = R200 / c

	halo = lguys.NFW(M_s=M_s, r_s=r_s)

	return halo
end

# ╔═╡ af38fc99-0822-4edb-b88b-20cf35a5ac21
NFW(10)

# ╔═╡ f6294a94-81e6-47a1-8131-23e33044a529
lguys.ρ_crit

# ╔═╡ a04c4445-9107-4546-8d73-fa043b68036e
function calc_V_circ(halo, r)
	x = r/halo.r_s

	v_v200_sq = halo._c/x * lguys.A_NFW(x) / lguys.A_NFW(halo._c)
	V200 = lguys.calc_V200(halo)
	return sqrt(v_v200_sq ) * V200
end

# ╔═╡ ca454ae0-4f9b-4bf5-8ac1-fd8d77a0159d
lguys.calc_V_circ_max(NFW(1.04)) * lguys.V0

# ╔═╡ 1ed56766-966d-4199-a969-f7c7c74ae5a0
p = NFW(1.04)

# ╔═╡ 17374f2c-9513-4231-a1bc-52d5532ef77c
lguys.calc_M200(p) / (4π / 3 * lguys.calc_R200(p) ^ 3)

# ╔═╡ f5a5a91c-0f3f-4904-b11b-1849ced8d83f
lguys.ρ_crit * 200

# ╔═╡ d64bd724-bcb2-4e44-baaf-073c8057ef32
lguys.calc_c(p)

# ╔═╡ e51b4948-b774-481b-9f77-d0771c9780d6
h

# ╔═╡ 0cd889e2-a6d7-4be4-bb31-c72563c8368e
lguys.calc_M200(p)

# ╔═╡ 3cbe361e-7448-4d67-94ca-13aa8e60e4ed
ρ_crit

# ╔═╡ 5c8e1f2d-747f-45f4-953d-fb02462c87e6
begin 
	M200 = find_zero(x->lguys.calc_V_circ_max(NFW(x)) - V_circ_max_0,  1)
	halo = NFW(M200)
end

# ╔═╡ 9a9e2e8c-af8c-4e00-98f7-c3c3759e4cd0
r_max = lguys.calc_r_max(halo)

# ╔═╡ 437d7d57-ca66-4e53-b7bb-21ac7bab5d69
V_circ_max = lguys.calc_V_circ_max(halo)

# ╔═╡ 21c7cecb-4a5f-4aba-94f1-47a1268a7595
md"""
## Results
- M_200 = $(round(M200, digits=4)) x 10^10 Msun
- r200 = $(round(halo.c*halo.r_s, digits=1)) kpc
- Ms = $(round(halo.M_s, digits=4)) x 10^10 Msun
- rs = $(round(halo.r_s, digits=2)) kpc
- M(r<rs) = $(round(halo.M_s *lguys.A_NFW(1), digits=4))
- c = $(round(halo.c, digits=3))
- Vmax = $(round(V_circ_max * lguys.V0, digits=1)) km / s
- rmax = $(round(r_max, digits=2)) kpc

"""

# ╔═╡ 832f5932-1950-4ea1-8b44-e6d40742752a
V_circ_max * lguys.V0

# ╔═╡ 2dbb0181-3018-4d68-a232-4d6d83350f7e
lguys.V0 * V_circ_max_0

# ╔═╡ 5e68334b-246d-4f17-9862-4a3cf1743f4f
begin
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="radius", ylabel="circular velocity")
	
	r = LinRange(0, 10, 1000)
	v = map(rr->lguys.calc_V_circ(halo, rr), r)
	plot!(r, v)
	hlines!([lguys.calc_V_circ_max(halo), V_circ_max_0])
	vlines!(lguys.calc_r_max(halo), )

	fig
end

# ╔═╡ b0197f89-92a5-4ac3-9ab2-b0d05aa368d7
lguys.calc_V_circ(halo, r_max)

# ╔═╡ d07c09ea-d1e3-4837-9969-46d3e4760235
lguys.calc_M(halo, lguys.calc_R200(halo))

# ╔═╡ bc0e8db9-ec7b-4c28-a326-8ffae92990e3
M200

# ╔═╡ 6bbeae54-9ff6-4935-9f12-8b0980dace52
lguys.calc_M(halo, 0.5)

# ╔═╡ 080e2e46-9eea-4010-abe2-c2286a00733e
md"""
## sanity checks
"""

# ╔═╡ 2232ac95-2873-4ef9-be90-96f1d0fe8049
function calc_V_max2(M200)
	r_max = calc_r_max(M200)
	c = calc_c(M200, 0)
	r = calc_r_s(M200)
	r200 = c*r
	G = 4.30091e-6 # km^2/s^2 kpc / M⊙

	V200 = sqrt(G * M200 / r200)

	x = r_max/r
	v_rel_sq = c/x * A_NFW(x)/A_NFW(c)
	return V200 * sqrt(v_rel_sq)
end

# ╔═╡ 52bafbb3-2424-4c6a-8f90-edc419351c5a
function calc_V_circ2(halo, r)
	M = calc_M(halo, r)
	return sqrt(G * M / r)
end

# ╔═╡ 7c154d6a-dbfc-415a-81d1-351257f8cdbb
function calc_V_max3(halo) # this one is wrong somewhere
	return 1.64 * halo.r_s * √(G * calc_ρ_s(halo))
end

# ╔═╡ f2c8daac-63a4-488d-af2a-4a5d3c7733ce
lguys.calc_M200(halo)

# ╔═╡ f677644c-ff10-4f0e-b6e3-26622ab92e43
md"""
# The full curve
"""

# ╔═╡ 8787c646-68d7-4277-b044-8a8aee0497cd
begin 
	M200_mean = 10 .^ LinRange(-1, 1, 1000)
	halo_mean = NFW.(M200_mean)
	Vc_mean = lguys.calc_V_circ_max.(halo_mean)
	Ms_mean = M_s_from_vel.(Vc_mean)
end

# ╔═╡ 8681ae51-7b55-4c6a-9fac-4f62f01af755
N_samples = 10_000

# ╔═╡ aeebdad1-552e-44f3-9b25-37b752ee5e69
begin 
	M200_samples = 10 .^ ( -1 .+ 2 * rand(N_samples))
	
	σ_c = 0.1
	c_samples = calc_c.(M200_samples, 0) .+ σ_c .* randn(N_samples)

	halo_samples = NFW_from_M200_c.(M200_samples, c_samples)

	Vc_samples = lguys.calc_V_circ_max.(halo_samples)
end

# ╔═╡ d8c978f3-57e9-48cb-b20a-6f441a241f7e
begin 
	Ms_samples = M_s_from_vel.(Vc_samples) .* 10 .^ (0 .+ (0.2 .- 0.5.* Vc_samples) .* randn(N_samples))
end

# ╔═╡ 51f743f3-001f-408d-8e11-4880a17d670e
hist(Vc_samples)

# ╔═╡ 90d563f3-d9a4-4e50-82b9-e37eecff87d9
ax_M200_kwargs = (xscale=log10, 
xlabel="M200 [code units]",
xticks = [0.1, 1, 10])

# ╔═╡ a620cbcc-5ac1-4f2a-ae93-6c345ee5f5ce
let
	fig = Figure()
	ax_ms = Axis(fig[1, 1];
		ylabel = "stellar mass [code units]",
		yscale = log10, 
		ax_M200_kwargs...
	)


	scatter!(M200_samples, Ms_samples, color=:black, alpha=0.05)
	lines!(M200_mean, Ms_mean)
	
	fig
end

# ╔═╡ 736c8755-6ae5-44c9-aa17-44a783c819d5
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = "concentration",
		ax_M200_kwargs...
	)

	y = [h.c for h in halo_samples]
	scatter!(M200_samples, y, color=:black, alpha=0.01)
	
	y = [h.c for h in halo_mean]
	lines!(M200_mean, y)

	fig
end

# ╔═╡ d0a5efe1-b0d3-4804-93a0-a98124452a1a
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = L"M_s",
		yscale=log10,
	yticks=[0.1, 1],
		ax_M200_kwargs...
	)
	y = [h.M_s for h in halo_mean]
	lines!(M200_mean, y)

	fig
end

# ╔═╡ 77315496-903c-491a-becc-c31c0bf4c88f
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = L"r_s",
		ax_M200_kwargs...
	)

	y = [h.r_s for h in halo_samples]

	scatter!(M200_samples, y,  color=:black, alpha=0.1, markersize=3)

	y = [h.r_s for h in halo_mean]
	lines!(M200_mean, y)
	fig
end

# ╔═╡ 4f3b2cb5-8855-4521-895c-accea0f9c115
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = L"V_c",
		yscale=log10, 
		yticks=[20, 30, 40, 50, 60, 70, 80],
		ax_M200_kwargs...
	)

	scatter!(M200_samples, Vc_samples * lguys.V0,  color=:black, alpha=0.1, markersize=3)
	lines!(M200_mean, Vc_mean * lguys.V0)

	fig
end

# ╔═╡ 1edd1d75-3cbf-49be-a80b-425feff7c73e
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		yscale=log10,
		xscale=log10,
		xlabel="Vc / km/s",
		ylabel = "stellar mass [code units]"
	)

	scatter!(Vc_samples * lguys.V0, Ms_samples, alpha=0.05, color=:black, markersize=5)
	lines!(Vc_mean * lguys.V0, Ms_mean)

	fig
end

# ╔═╡ 023f12d0-317a-46db-b5c5-81a7d2216c55
let
	v_1 = LinRange(18, 100, 1000) ./ lguys.V0

	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "v_circ_max / km s",
		ylabel = "M_s",
		yscale=log10,
		xscale=log10,
		limits=((1, 100), (1e-8, 1e0)),
		xticks=[1, 10, 100],
		xminorticks=IntervalsBetween(9),
		aspect=1
	)

	scatter!(Vc_samples * lguys.V0, Ms_samples, alpha=0.05, color=:black, markersize=5)
	lines!(v_1 * lguys.V0, M_s_from_vel.(v_1))

	lines!(Vc_mean * lguys.V0, Ms_mean)
	fig

end

# ╔═╡ e6995206-ebf2-40aa-ae7e-cc101e07f8ac
lllerp_ms_vc = lguys.lerp(log10.(Vc_mean), log10.(Ms_mean))

# ╔═╡ 22b0843a-a90c-40ab-8edc-d82a982abd43
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		xscale=log10,
		xlabel="Vc / km/s",
		ylabel = "log difference in Ms"
	)

	scatter!(Vc_samples * lguys.V0, log10.(Ms_samples) .- lllerp_ms_vc.(log10.(Vc_samples)), alpha=0.05, color=:black, markersize=5)

	fig
end

# ╔═╡ 122794fb-82ab-4f3c-a458-29462fbab1fd
md"""
## Scratch space
"""

# ╔═╡ Cell order:
# ╟─07a710d8-0ae4-4d9f-9759-002750730010
# ╟─33f433ec-e94e-41fd-a808-6254bf4d34ce
# ╠═45ac95e4-a9a1-4f8f-9523-8f670f076ae5
# ╠═2c90aad4-87e7-11ee-2fbd-05f67c7e0343
# ╟─67e892f1-3969-484f-9671-e0ac018babf2
# ╟─fc490d05-d494-4771-9dec-7d1837485ff0
# ╠═e338d0b5-fd4f-43ed-9379-27d50787d387
# ╠═aec3d9a7-d447-4e9f-aa83-fa5b20541b5c
# ╠═21c7cecb-4a5f-4aba-94f1-47a1268a7595
# ╟─bdd43353-da4c-48fb-bb67-b2aec628dd71
# ╟─ca65b2b1-8414-48dc-badb-09c2d50879ad
# ╠═18a20162-1b0f-4da8-931d-5ff95da22f54
# ╠═3fc2fab3-2e37-4927-94f0-eeecaa89d575
# ╠═d9cb3699-6a0f-4d32-b3f8-410cf2c92019
# ╟─fd225a65-b833-4971-a308-09adc0f0b10f
# ╠═41b442e9-6691-4389-9517-7b0b55a25118
# ╠═ef4ad7fc-fdf8-4e60-b2a6-e583ae226749
# ╠═92def240-121d-467e-8026-f30e99e97b06
# ╠═c185ed67-e042-422f-9582-6a815c292fd2
# ╠═57071589-7f70-4739-bb6f-5ce184a8a09a
# ╠═effb14d6-56ad-4baa-a484-83820aa72bb4
# ╟─312e16fd-8a35-4067-bd02-3260c213a9eb
# ╠═56d46bf2-a4a3-4d93-97d7-44c63a35b1db
# ╠═a4c594f9-3de9-43d8-810f-b0f5d4cb693a
# ╠═af38fc99-0822-4edb-b88b-20cf35a5ac21
# ╠═f6294a94-81e6-47a1-8131-23e33044a529
# ╠═a04c4445-9107-4546-8d73-fa043b68036e
# ╠═ca454ae0-4f9b-4bf5-8ac1-fd8d77a0159d
# ╠═1ed56766-966d-4199-a969-f7c7c74ae5a0
# ╠═17374f2c-9513-4231-a1bc-52d5532ef77c
# ╠═f5a5a91c-0f3f-4904-b11b-1849ced8d83f
# ╠═d64bd724-bcb2-4e44-baaf-073c8057ef32
# ╠═e51b4948-b774-481b-9f77-d0771c9780d6
# ╠═0cd889e2-a6d7-4be4-bb31-c72563c8368e
# ╠═3cbe361e-7448-4d67-94ca-13aa8e60e4ed
# ╠═5c8e1f2d-747f-45f4-953d-fb02462c87e6
# ╠═9a9e2e8c-af8c-4e00-98f7-c3c3759e4cd0
# ╠═437d7d57-ca66-4e53-b7bb-21ac7bab5d69
# ╠═832f5932-1950-4ea1-8b44-e6d40742752a
# ╠═2dbb0181-3018-4d68-a232-4d6d83350f7e
# ╠═5e68334b-246d-4f17-9862-4a3cf1743f4f
# ╠═b0197f89-92a5-4ac3-9ab2-b0d05aa368d7
# ╠═d07c09ea-d1e3-4837-9969-46d3e4760235
# ╠═bc0e8db9-ec7b-4c28-a326-8ffae92990e3
# ╠═6bbeae54-9ff6-4935-9f12-8b0980dace52
# ╟─080e2e46-9eea-4010-abe2-c2286a00733e
# ╠═2232ac95-2873-4ef9-be90-96f1d0fe8049
# ╠═52bafbb3-2424-4c6a-8f90-edc419351c5a
# ╠═7c154d6a-dbfc-415a-81d1-351257f8cdbb
# ╠═f2c8daac-63a4-488d-af2a-4a5d3c7733ce
# ╠═f677644c-ff10-4f0e-b6e3-26622ab92e43
# ╠═8787c646-68d7-4277-b044-8a8aee0497cd
# ╠═8681ae51-7b55-4c6a-9fac-4f62f01af755
# ╠═aeebdad1-552e-44f3-9b25-37b752ee5e69
# ╠═d8c978f3-57e9-48cb-b20a-6f441a241f7e
# ╠═51f743f3-001f-408d-8e11-4880a17d670e
# ╠═90d563f3-d9a4-4e50-82b9-e37eecff87d9
# ╠═a620cbcc-5ac1-4f2a-ae93-6c345ee5f5ce
# ╠═736c8755-6ae5-44c9-aa17-44a783c819d5
# ╠═d0a5efe1-b0d3-4804-93a0-a98124452a1a
# ╠═77315496-903c-491a-becc-c31c0bf4c88f
# ╠═4f3b2cb5-8855-4521-895c-accea0f9c115
# ╠═1edd1d75-3cbf-49be-a80b-425feff7c73e
# ╠═023f12d0-317a-46db-b5c5-81a7d2216c55
# ╠═e6995206-ebf2-40aa-ae7e-cc101e07f8ac
# ╠═22b0843a-a90c-40ab-8edc-d82a982abd43
# ╠═122794fb-82ab-4f3c-a458-29462fbab1fd
