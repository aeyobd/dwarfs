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

# ╔═╡ f4269d1a-940d-41da-9183-5130978f7d7b
using PairPlots

# ╔═╡ ef991557-7bf7-4add-b29f-f16187297e46
using DataFrames

# ╔═╡ e7ab194c-63a4-4274-aaba-43c3d369ce0d
using KernelDensity

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

# ╔═╡ 21619de9-8960-46fe-8a0c-808776b33c6d
import StatsBase: quantile, median

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
	M_L_star = 1.5
	M_L_star_err = 0.3
	MV = -10.82
	MV_err = 0.14
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

# ╔═╡ 2ec0b911-6338-4364-a250-90664e90f1a6
lMs_to_lVc_err = 0.1

# ╔═╡ 40694bea-0b6e-4235-ac1b-6b25dd91ac38
md"""
Compare to the figure 1 in the paper.
"""

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

# ╔═╡ 50ab4a6c-c40c-45e7-90aa-759e00a415be
md"""
compare to appendix figure
"""

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

	halo = lguys.NFW(M_s=M_s, r_s=r_s, c=c)

	return halo
end

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

# ╔═╡ a04c4445-9107-4546-8d73-fa043b68036e
function calc_V_circ(halo, r)
	x = r/halo.r_s

	v_v200_sq = halo._c/x * lguys.A_NFW(x) / lguys.A_NFW(halo._c)
	V200 = lguys.calc_V200(halo)
	return sqrt(v_v200_sq ) * V200
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

# ╔═╡ f677644c-ff10-4f0e-b6e3-26622ab92e43
md"""
# The full curve
"""

# ╔═╡ 8787c646-68d7-4277-b044-8a8aee0497cd
begin 
	M200_mean = 10 .^ LinRange(-2, 2, 1000)
	c_mean2 = calc_c.(M200_mean, 0)
	halo_mean = NFW.(M200_mean)
	Vc_mean = lguys.calc_V_circ_max.(halo_mean)
	Rc_mean = lguys.calc_r_circ_max.(halo_mean)
	Ms_mean = M_s_from_vel.(Vc_mean)
end

# ╔═╡ 6213f5b6-0465-482d-a530-727f49220d79
lMs_to_lVc = lguys.lerp(log10.(Ms_mean), log10.(Vc_mean))

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

	x = 10 .^ LinRange(0, -8, 1000) 
	y = 10 .^ lMs_to_lVc.(log10.(x))
	lines!(y * lguys.V0, x, linestyle=:dash, label="interpolated")

	axislegend(position=:lt)

	fig

end

# ╔═╡ bab1a32c-8d35-4958-af84-4b239175f602
lVc_to_lM200 = lguys.lerp(log10.(Vc_mean), log10.(M200_mean))

# ╔═╡ 8681ae51-7b55-4c6a-9fac-4f62f01af755
N_samples = 10000

# ╔═╡ aeebdad1-552e-44f3-9b25-37b752ee5e69
begin 
	M200_samples = 10 .^ ( -1 .+ 2 * rand(N_samples))
	
	σ_c = 0.1 #uncertainty from paper, approx 0.09 scatter + 0.03 model -> fitting function
	c_samples = calc_c.(M200_samples, 0) .* 10 .^ (0 .+ σ_c .* randn(N_samples))
end

# ╔═╡ c9e713c8-74b9-4b3d-bc7d-071955dc4615
halo_samples = NFW_from_M200_c.(M200_samples, c_samples)

# ╔═╡ d9db5c56-b8d6-4e08-aa47-979c191496b0
Vc_samples = lguys.calc_V_circ_max.(halo_samples)

# ╔═╡ 9695b55c-ff53-4369-9d12-324dd35ab021
Rc_samples = lguys.calc_r_circ_max.(halo_samples)

# ╔═╡ d8c978f3-57e9-48cb-b20a-6f441a241f7e
begin 
	Ms_samples = M_s_from_vel.(
		Vc_samples .* 10 .^ (lMs_to_lVc_err * randn(N_samples))
	)
end

# ╔═╡ 90d563f3-d9a4-4e50-82b9-e37eecff87d9
ax_M200_kwargs = (xscale=log10, 
xlabel="M200 [code units]",
xticks = [0.1, 1, 10])

# ╔═╡ a2ca4b2b-b198-4026-ac67-4aad70cdf1fa
0.15*lguys.V0

# ╔═╡ 4f971234-565a-43a4-a4bf-94a9415e6d7e
Vc_label = L"V_\textrm{circ,\,max}\  / \textrm{km\,s^{-1}}";

# ╔═╡ 4bb1dcef-04e2-4aa5-aa8f-d76a21502e01
rc_label = L"r_\textrm{circ,\,max}\  / \textrm{kpc}"

# ╔═╡ a5ed13a2-41f1-40dc-8eff-27fcb7e5aca9
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = Vc_label,
		xlabel = rc_label,
		yscale=log10, 
		yticks=[10, 30, 100],
		xscale=log10,
	)

	scatter!(Rc_samples, Vc_samples * lguys.V0,  color=:black, alpha=0.1, markersize=3)
	lines!(Rc_mean, Vc_mean * lguys.V0)

	fig
end

# ╔═╡ 023f12d0-317a-46db-b5c5-81a7d2216c55
let
	v_1 = LinRange(18, 100, 1000) ./ lguys.V0

	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = Vc_label,
		ylabel = L"M_\star",
		yscale=log10,
		xscale=log10,
		limits=((10, 100), (1e-8, 1e0)),
		xticks=[10, 100],
		xminorticks=IntervalsBetween(9),
		aspect=1
	)

	scatter!(Vc_samples * lguys.V0, Ms_samples, alpha=0.05, color=:black, markersize=5)
	lines!(v_1 * lguys.V0, M_s_from_vel.(v_1))

	lines!(Vc_mean * lguys.V0, Ms_mean)
	fig

end

# ╔═╡ be32e183-a0c4-4261-a3a8-3486db0734b7
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = L"M_\star",
		xlabel = L"r_\textrm{max}",
		yscale=log10, 
		xscale=log10,
		limits=(nothing, (1e-8, 1))
	)

	scatter!(Rc_samples, Ms_samples,  color=:black, alpha=0.1, markersize=3)
	lines!(Rc_mean, Ms_mean)

	fig
end

# ╔═╡ e6995206-ebf2-40aa-ae7e-cc101e07f8ac
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ 122794fb-82ab-4f3c-a458-29462fbab1fd
md"""
## Monte Carlo mass estimates
"""

# ╔═╡ fa6b2468-164f-42b7-88f5-6eff0a9c5ce1
function solve_M200_c(Vcmax, σ_c=σ_c)
	dc = 10 ^ (0 + σ_c * randn())

	f(M200) = lguys.calc_V_circ_max(NFW_from_M200_c(M200, dc * calc_c(M200, 0.))) - Vcmax

	M200 = find_zero(f, [0.001, 100])
	return M200, calc_c(M200, 0.) * dc
end

# ╔═╡ 197cbbf0-beb9-45e4-99dc-5d9baa58c696
function solve_rmax(Vcmax, σ_c=0)
	M200, c = solve_M200_c(Vcmax, σ_c)
	return lguys.calc_r_circ_max(lguys.NFW(M200=M200, c=c))
end

# ╔═╡ bbd5e178-eddc-4e75-b4c5-252048b75fb4
begin 
	samples = DataFrame(
		MV = MV .+ MV_err * randn(N_samples),
	)

	samples[!, :L] = mag_to_L.(samples.MV, MV_sol)
	samples[!, :Y] = M_L_star .+ M_L_star_err .* randn(N_samples)
	samples[!, :Ms] = samples.L .* samples.Y ./lguys.M0
	samples[!, :log_Vc] = (
		lMs_to_lVc.(log10.(samples.Ms)) .+ lMs_to_lVc_err * randn(N_samples)
	)

	samples[!, :Vc] = 10 .^ samples.log_Vc
	Mcs = [solve_M200_c(x) for x in samples.Vc]
	samples[!, :M200] = first.(Mcs)
	samples[!, :log_M200] = log10.(samples.M200)

	samples[!, :c] = last.(Mcs)
	halos = NFW_from_M200_c.(samples.M200, samples.c)

	halos = [lguys.NFW(M200=m, c=c) for (m, c) in zip(samples.M200, samples.c)]
	samples[!, :r_max] = lguys.calc_r_circ_max.(halos)
	samples[!, :V_max] = lguys.calc_V_circ_max.(halos)
	samples[!, :Ms2] = samples.Ms * 1e4

	samples
end

# ╔═╡ 5cdc46ec-6d1b-403f-a3a2-a15b2c82ee16
pairplot(samples[:, [:MV, :L, :Ms]])

# ╔═╡ 7fd48721-749d-4e99-91cc-5ffde830487d
pairplot(samples[:, [:Ms2, :V_max, :r_max]], labels=Dict(
	:Ms2 => L"$M_s$", 
	:log_M200 => L"\log M_{200} / 10^{10}\,\textrm{M}_\odot ", 
	:c => L"c", ))

# ╔═╡ 4e1290b5-171e-4715-a87b-28a2cfcb4325
samples.Ms * 1e10

# ╔═╡ 9fe20cff-b7a0-4e1d-ac7c-f9a94bba9a0a
lc_err = 0.01

# ╔═╡ 7d1ba57d-4332-4d6e-872b-42925c14ecfb
c_mean = [h.c for h in halo_mean]

# ╔═╡ d2616939-fba8-42ff-921a-8b0bcaf7acb9
function describe(x::Array; p=0.16)
	p1 = quantile(x, p)
	p2 = quantile(x, 1-p)
	m = median(x)

	return p1-m, m, p2-m
end

# ╔═╡ 4c0c93ad-5471-45d4-b44a-2d15a782491b
describe(samples.L)

# ╔═╡ 10aac004-8c52-4c91-944a-0002bdffa99d
describe(samples.Ms) .* lguys.M0 ./ 1e6

# ╔═╡ 0ba3f10c-46ad-48d5-9b15-eb67e9e505ed
describe(samples.M200)

# ╔═╡ ecade01b-f703-4780-bc5d-ccc4c448b676
describe(samples.c)

# ╔═╡ b85256a1-786f-4dee-a6f1-f55406c3b18e
md"""
# Self-consistency calculator
"""

# ╔═╡ c49ac57e-8e8d-4ed6-ad35-be400863f6b4
begin 
	V_max_in = 50 / lguys.V0
	# fiducial is 31.3
	r_max_in = 10.783

	V_max_in
end

# ╔═╡ 5e71b022-bd95-4c3c-8930-51100fb9ab1c
r_max_exp = solve_rmax(V_max_in, 0)

# ╔═╡ bbee444e-079b-4208-9faf-0a7fe5f81455
let
	fig, ax = FigAxis(
		xlabel="V max",
		ylabel="r max / kpc",
		limits=((00, 70), (0, 20))
	)

	x = samples.Vc * lguys.V0
	y = samples.r_max
	k = kde([x y])

	h = hexbin!(x, y, bins=50, colorscale=log10)

	c = contour!(k, colormap=:bluesgreens)

	scatter!(V_max_in * lguys.V0, r_max_in)

	lines!(Vc_mean * lguys.V0, Rc_mean)
	
	fig
end

# ╔═╡ 6baff9d8-a96d-4d1b-898f-089003459c19
let
	v_1 = LinRange(18, 100, 1000) ./ lguys.V0

	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = Vc_label,
		ylabel = L"M_\star",
		yscale=log10,
		xscale=log10,
		limits=((10, 100), (1e-5, 1e-3)),
		xticks=[10, 100],
		xminorticks=IntervalsBetween(9),
		aspect=1
	)

	x = samples.Vc * lguys.V0
	y = samples.Ms
	k = kde([x y])

	c = contour!(k, colormap=:bluesgreens)
	
	
	lines!(v_1 * lguys.V0, M_s_from_vel.(v_1))

	lines!(Vc_mean * lguys.V0, Ms_mean)

	vlines!(V_max_in * lguys.V0)
	fig

end

# ╔═╡ 2b1d7ec9-b3f7-4230-9f5a-80ece86f709d
halo_in = lguys.NFW(V_circ_max=V_max_in, r_circ_max=r_max_in)

# ╔═╡ ecda5f20-27cd-41a8-8545-9f3a6b91a80e
lguys.calc_M200(halo_in)

# ╔═╡ ab5bc36d-d203-45be-8e7d-d8c5c9473388
lguys.calc_R200(halo_in)

# ╔═╡ 84892a5a-f816-450b-9757-e4135a40aebc
lguys.calc_V_circ_max(halo_in) * lguys.V0

# ╔═╡ cb42028b-eeb0-48a7-b23d-5abe792eb4f0
md"""
# calculating H
"""

# ╔═╡ e582e62c-d851-4dd0-95d4-82fd1d97c26c
function calc_h(halo, N=1e6)
	return 4 * lguys.calc_R200(halo) / sqrt(N)
end

# ╔═╡ c3e1326c-772b-4d8c-aabe-a997b77bede4
calc_h(halo_in)

# ╔═╡ ba292975-8eea-4c3a-8b9c-23ff1b60c7ab
Np = 1e6

# ╔═╡ a75e163e-3489-4f82-901b-c511d1e44395
m = lguys.calc_M200(halo_in) / Np

# ╔═╡ 8b7f61d5-dcef-4889-8aab-a85b0ccdaf02
let
	fig, ax = FigAxis(
		xlabel = "log r / kpc",
		ylabel = "N ( r < r)"
	)

	lx = LinRange(-2, 2, 100)
	
	y = lguys.calc_M.(halo_in, 10 .^ lx) ./ m
	

	lines!(lx, log10.(y))
	vlines!(log10.(h), label="softening")
	vlines!(log10.(r_max_in), label="rmax")
	axislegend()
	fig
end

# ╔═╡ 93715b38-1790-4b1e-a551-17d69b68876f
lguys.calc_M(halo_in, lguys.calc_R200(halo_in))

# ╔═╡ 896ccd5d-8e92-4c26-944d-04d42550caef
lguys.calc_R200(halo_in)

# ╔═╡ 2ae19d4f-659f-4985-8fc2-98d94cc52730
4^2 * lguys.G * m / calc_h(halo_in)^2

# ╔═╡ b97bdec3-2dbe-4a1d-bad9-4116b2d2e614
lguys.G * lguys.calc_M200(halo_in) / lguys.calc_R200(halo_in)^2

# ╔═╡ Cell order:
# ╟─07a710d8-0ae4-4d9f-9759-002750730010
# ╟─33f433ec-e94e-41fd-a808-6254bf4d34ce
# ╟─45ac95e4-a9a1-4f8f-9523-8f670f076ae5
# ╠═2c90aad4-87e7-11ee-2fbd-05f67c7e0343
# ╠═f4269d1a-940d-41da-9183-5130978f7d7b
# ╠═ef991557-7bf7-4add-b29f-f16187297e46
# ╠═21619de9-8960-46fe-8a0c-808776b33c6d
# ╟─67e892f1-3969-484f-9671-e0ac018babf2
# ╟─fc490d05-d494-4771-9dec-7d1837485ff0
# ╠═e338d0b5-fd4f-43ed-9379-27d50787d387
# ╠═aec3d9a7-d447-4e9f-aa83-fa5b20541b5c
# ╟─bdd43353-da4c-48fb-bb67-b2aec628dd71
# ╟─ca65b2b1-8414-48dc-badb-09c2d50879ad
# ╠═18a20162-1b0f-4da8-931d-5ff95da22f54
# ╠═6213f5b6-0465-482d-a530-727f49220d79
# ╠═2ec0b911-6338-4364-a250-90664e90f1a6
# ╟─d9cb3699-6a0f-4d32-b3f8-410cf2c92019
# ╠═40694bea-0b6e-4235-ac1b-6b25dd91ac38
# ╟─fd225a65-b833-4971-a308-09adc0f0b10f
# ╟─41b442e9-6691-4389-9517-7b0b55a25118
# ╠═ef4ad7fc-fdf8-4e60-b2a6-e583ae226749
# ╠═92def240-121d-467e-8026-f30e99e97b06
# ╠═c185ed67-e042-422f-9582-6a815c292fd2
# ╠═57071589-7f70-4739-bb6f-5ce184a8a09a
# ╟─effb14d6-56ad-4baa-a484-83820aa72bb4
# ╠═50ab4a6c-c40c-45e7-90aa-759e00a415be
# ╠═bab1a32c-8d35-4958-af84-4b239175f602
# ╟─312e16fd-8a35-4067-bd02-3260c213a9eb
# ╠═56d46bf2-a4a3-4d93-97d7-44c63a35b1db
# ╠═a4c594f9-3de9-43d8-810f-b0f5d4cb693a
# ╟─080e2e46-9eea-4010-abe2-c2286a00733e
# ╠═2232ac95-2873-4ef9-be90-96f1d0fe8049
# ╠═a04c4445-9107-4546-8d73-fa043b68036e
# ╠═52bafbb3-2424-4c6a-8f90-edc419351c5a
# ╠═7c154d6a-dbfc-415a-81d1-351257f8cdbb
# ╠═f677644c-ff10-4f0e-b6e3-26622ab92e43
# ╠═8787c646-68d7-4277-b044-8a8aee0497cd
# ╠═8681ae51-7b55-4c6a-9fac-4f62f01af755
# ╠═aeebdad1-552e-44f3-9b25-37b752ee5e69
# ╠═c9e713c8-74b9-4b3d-bc7d-071955dc4615
# ╠═d9db5c56-b8d6-4e08-aa47-979c191496b0
# ╠═9695b55c-ff53-4369-9d12-324dd35ab021
# ╠═d8c978f3-57e9-48cb-b20a-6f441a241f7e
# ╠═90d563f3-d9a4-4e50-82b9-e37eecff87d9
# ╠═a2ca4b2b-b198-4026-ac67-4aad70cdf1fa
# ╠═a5ed13a2-41f1-40dc-8eff-27fcb7e5aca9
# ╠═4f971234-565a-43a4-a4bf-94a9415e6d7e
# ╠═4bb1dcef-04e2-4aa5-aa8f-d76a21502e01
# ╠═023f12d0-317a-46db-b5c5-81a7d2216c55
# ╠═be32e183-a0c4-4261-a3a8-3486db0734b7
# ╠═e6995206-ebf2-40aa-ae7e-cc101e07f8ac
# ╟─122794fb-82ab-4f3c-a458-29462fbab1fd
# ╠═fa6b2468-164f-42b7-88f5-6eff0a9c5ce1
# ╠═197cbbf0-beb9-45e4-99dc-5d9baa58c696
# ╠═bbd5e178-eddc-4e75-b4c5-252048b75fb4
# ╠═5cdc46ec-6d1b-403f-a3a2-a15b2c82ee16
# ╠═7fd48721-749d-4e99-91cc-5ffde830487d
# ╠═4e1290b5-171e-4715-a87b-28a2cfcb4325
# ╠═9fe20cff-b7a0-4e1d-ac7c-f9a94bba9a0a
# ╠═7d1ba57d-4332-4d6e-872b-42925c14ecfb
# ╠═d2616939-fba8-42ff-921a-8b0bcaf7acb9
# ╠═4c0c93ad-5471-45d4-b44a-2d15a782491b
# ╠═10aac004-8c52-4c91-944a-0002bdffa99d
# ╠═0ba3f10c-46ad-48d5-9b15-eb67e9e505ed
# ╠═ecade01b-f703-4780-bc5d-ccc4c448b676
# ╟─b85256a1-786f-4dee-a6f1-f55406c3b18e
# ╠═e7ab194c-63a4-4274-aaba-43c3d369ce0d
# ╠═c49ac57e-8e8d-4ed6-ad35-be400863f6b4
# ╠═5e71b022-bd95-4c3c-8930-51100fb9ab1c
# ╠═bbee444e-079b-4208-9faf-0a7fe5f81455
# ╠═6baff9d8-a96d-4d1b-898f-089003459c19
# ╠═2b1d7ec9-b3f7-4230-9f5a-80ece86f709d
# ╠═ecda5f20-27cd-41a8-8545-9f3a6b91a80e
# ╠═ab5bc36d-d203-45be-8e7d-d8c5c9473388
# ╠═84892a5a-f816-450b-9757-e4135a40aebc
# ╟─cb42028b-eeb0-48a7-b23d-5abe792eb4f0
# ╠═e582e62c-d851-4dd0-95d4-82fd1d97c26c
# ╠═c3e1326c-772b-4d8c-aabe-a997b77bede4
# ╠═ba292975-8eea-4c3a-8b9c-23ff1b60c7ab
# ╠═a75e163e-3489-4f82-901b-c511d1e44395
# ╠═8b7f61d5-dcef-4889-8aab-a85b0ccdaf02
# ╠═93715b38-1790-4b1e-a551-17d69b68876f
# ╠═896ccd5d-8e92-4c26-944d-04d42550caef
# ╠═2ae19d4f-659f-4985-8fc2-98d94cc52730
# ╠═b97bdec3-2dbe-4a1d-bad9-4116b2d2e614
