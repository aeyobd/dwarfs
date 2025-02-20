### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 2c90aad4-87e7-11ee-2fbd-05f67c7e0343
begin 
	using Pkg; Pkg.activate()

	using LilGuys
	
	using CairoMakie, Arya
end

# ╔═╡ f4269d1a-940d-41da-9183-5130978f7d7b
using PairPlots

# ╔═╡ ef991557-7bf7-4add-b29f-f16187297e46
using DataFrames

# ╔═╡ ff9ea455-d618-4ebf-b9e9-5d92e23b4f37
using Measurements

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

# ╔═╡ 60eba328-c047-4b59-8b8c-8b24a4884541
import CairoMakie: save

# ╔═╡ 21619de9-8960-46fe-8a0c-808776b33c6d
import StatsBase: quantile, median

# ╔═╡ a334b0a2-b408-4c81-a074-67cb1bcb58be
import TOML

# ╔═╡ 24b4184d-7f69-47fc-a758-50f21840283a
obs_props = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml")

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
	M_L_star = obs_props["M_L_s"]
	M_L_star_err = obs_props["M_L_s_err"]
	MV = obs_props["Mv"]
	MV_err = obs_props["Mv_err"]
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
	ν = v_max * V2KMS / 50 # km/s
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

Compare the below plot to their appendix figure.
"""

# ╔═╡ effb14d6-56ad-4baa-a484-83820aa72bb4
let
	M = 10 .^ LinRange(-8, 17, 1000) / 1e10
	h = 1 # don't remember the purpose of this?
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log M_{200}\ / \ M_\odot",
		ylabel=L"\log c",
	)

	for z in [0, 1, 2, 3, 5, 7]
		lines!(log10.(M / h * M2MSUN), log10.(LilGuys.Ludlow.c_ludlow.(M, z)), label="z=$z", color=z, colorrange=(0, 8), colormap=:blues)
	end

	axislegend()
	fig
end

# ╔═╡ f677644c-ff10-4f0e-b6e3-26622ab92e43
md"""
# The full curve
"""

# ╔═╡ 8787c646-68d7-4277-b044-8a8aee0497cd
begin 
	M200_mean = 10 .^ LinRange(-2, 1.5, 1000)
	c_mean2 = LilGuys.Ludlow.c_ludlow.(M200_mean, 0)
	halo_mean = [NFW(M200=M200_mean[i], c=c_mean2[i]) for i in eachindex(M200_mean)]
	Vc_mean = calc_v_circ_max.(halo_mean)
	Rc_mean = calc_r_circ_max.(halo_mean)
	Ms_mean = M_s_from_vel.(Vc_mean)
end

# ╔═╡ 6213f5b6-0465-482d-a530-727f49220d79
lMs_to_lVc = LilGuys.lerp(log10.(Ms_mean), log10.(Vc_mean))

# ╔═╡ d9cb3699-6a0f-4d32-b3f8-410cf2c92019
let
	v_1 = LinRange(18, 100, 1000) ./V2KMS

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
	lines!(v_1 * V2KMS,  M_s_from_vel.(v_1))

	x = 10 .^ LinRange(0, -8, 1000) 
	y = 10 .^ lMs_to_lVc.(log10.(x))
	lines!(y * V2KMS, x, linestyle=:dash, label="interpolated")

	axislegend(position=:lt)

	fig

end

# ╔═╡ bab1a32c-8d35-4958-af84-4b239175f602
lVc_to_lM200 = LilGuys.lerp(log10.(Vc_mean), log10.(M200_mean))

# ╔═╡ 8681ae51-7b55-4c6a-9fac-4f62f01af755
N_samples = 10000

# ╔═╡ aeebdad1-552e-44f3-9b25-37b752ee5e69
begin 
	M200_samples = 10 .^ ( -1 .+ 3 * rand(N_samples))
	
	σ_c = 0.1 #uncertainty from paper, approx 0.09 scatter + 0.03 model -> fitting function
	c_samples = LilGuys.Ludlow.c_ludlow.(M200_samples, 0) .* 10 .^ (0 .+ σ_c .* randn(N_samples))
end

# ╔═╡ c9e713c8-74b9-4b3d-bc7d-071955dc4615
halo_samples = [NFW(M200=M200, c=c) for (M200, c) in zip(M200_samples, c_samples)]

# ╔═╡ d9db5c56-b8d6-4e08-aa47-979c191496b0
Vc_samples = LilGuys.calc_v_circ_max.(halo_samples)

# ╔═╡ 9695b55c-ff53-4369-9d12-324dd35ab021
Rc_samples = LilGuys.calc_r_circ_max.(halo_samples)

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
0.15*V2KMS

# ╔═╡ a5ed13a2-41f1-40dc-8eff-27fcb7e5aca9
let
	fig = Figure()

	ax_c = Axis(fig[1, 1];
		ylabel = L"$\log\, v_\textrm{circ\,max}$ / km\,s$^{-1}$",
		xlabel = L"$\log \,r_\textrm{circ\,max}$ / kpc",
	)



	scatter!(log10.(Rc_samples), log10.(Vc_samples * V2KMS),
		color=:black, alpha=0.1, markersize=3)
	lines!(log10.(Rc_mean), log10.(Vc_mean * V2KMS))
	v = log10.(Vc_mean * V2KMS)

	y = LilGuys.Ludlow.solve_rmax.(Vc_mean, 0.1)
	lines!(log10.(y), v, color=COLORS[1], linestyle=:dash)

	y = LilGuys.Ludlow.solve_rmax.(Vc_mean, -0.1)
	lines!(log10.(y), v, color=COLORS[1], linestyle=:dash)
	hlines!(log10(50))
	fig
end

# ╔═╡ 60c1f90f-6604-464e-bc77-d182a8cc8785
argmin(abs.(Vc_mean .- 50 / V2KMS))

# ╔═╡ 4f971234-565a-43a4-a4bf-94a9415e6d7e
Vc_label = L"V_\textrm{circ,\,max}\  / \textrm{km\,s^{-1}}";

# ╔═╡ 4bb1dcef-04e2-4aa5-aa8f-d76a21502e01
rc_label = L"r_\textrm{circ,\,max}\  / \textrm{kpc}"

# ╔═╡ 023f12d0-317a-46db-b5c5-81a7d2216c55
let
	v_1 = LinRange(18, 100, 1000) ./ V2KMS

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

	scatter!(Vc_samples * V2KMS, Ms_samples, alpha=0.05, color=:black, markersize=5)
	lines!(v_1 * V2KMS, M_s_from_vel.(v_1))

	lines!(Vc_mean * V2KMS, Ms_mean)
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

# ╔═╡ 8109003b-60cf-458a-a43d-9cd084d27344
md"""
# Mean (by hand)
"""

# ╔═╡ db6f7825-4f9a-48ef-8700-6f4434032751
MV

# ╔═╡ 0541d9c1-c6fb-4bb6-8541-7529a65e446d
mag_to_L(MV ± MV_err, MV_sol)

# ╔═╡ b9d918b4-3c03-474c-8659-2ccc601c487d
Lm, Lm_err = mag_to_L(MV, MV_sol) .* (1, MV_err .* 0.2 / log(10))

# ╔═╡ da3746bb-db05-4981-9f77-fc92238f61ed
Ms_m = Lm * M_L_star / 1e10

# ╔═╡ 07cb7647-be33-4a16-8382-d44302f9e6d3
 lMs_to_lVc(log10(Ms_m)) + log10(V2KMS)

# ╔═╡ 02c22085-8e4d-4ae8-a1d7-ac746dc5b42f
vc_m_err = 0.1

# ╔═╡ c0f9d7f0-7983-4a44-b85c-15534623d534
lMs_to_lVc(log10(Ms_m*0.2))

# ╔═╡ 122794fb-82ab-4f3c-a458-29462fbab1fd
md"""
## Monte Carlo mass estimates
"""

# ╔═╡ 667937b2-006f-4429-b134-91934adb53c8
stellar_profile = LilGuys.Exp2D(R_s=0.13) # best fit to present day

# ╔═╡ 5dc54372-14c4-424b-8419-55d639f00baa
function calc_σv_star_mean_old(p; stellar_profile=stellar_profile, R_max=Inf, R = 10 .^ LinRange(-3, 3, 1000))
	integrand(r) = calc_ρ(stellar_profile, r) * calc_M(p, r) * LilGuys.G / r^2
	
	weighted_sigma(r) = LilGuys.integrate(integrand, r, R_max) * 4π * r^2
	mass(r) = calc_ρ(stellar_profile, r) * 4π * r^2
	sigmas = weighted_sigma.(R)
	
	sqrt(sum(sigmas) / sum(mass.(R)))
end

# ╔═╡ 28e09f56-7a1c-4fc9-95da-36ef0ed75a74
function calc_σv_star_mean(p; stellar_profile=stellar_profile, R_min=0, R_max=Inf)

	ρ(r) = calc_ρ(stellar_profile, r)
	M(r) = calc_M(p, r)
	integrand_1(r) = ρ(r) * M(r) * LilGuys.G / r^2
	integrand_2(r) = ρ(r) * M(r)  * LilGuys.G * 4π/3 * r

	weighted_σ2 = (
		#4π * R_max * LilGuys.integrate(integrand_1, R_min, R_max) -
		LilGuys.integrate(integrand_2, R_min, R_max)
	)

	Mtot = LilGuys.calc_M(stellar_profile, R_max)
	sqrt(weighted_σ2 / Mtot)
end

# ╔═╡ 0ad2bc38-9145-4e3a-8a07-2db6dad4404d
0.04 * LilGuys.V2KMS

# ╔═╡ 9fe20cff-b7a0-4e1d-ac7c-f9a94bba9a0a
lc_err = 0.10

# ╔═╡ e65a951d-21ce-45fe-837c-825c782e76d2
function sample_halo(;
		MV=MV, MV_err=MV_err, lc_err=lc_err, lMs_to_lVc_err=lMs_to_lVc_err

	)
    # Sample a single halo and return a dictionary of the results
    MV = MV + MV_err * randn()
    L = mag_to_L(MV, MV_sol)
    Y = M_L_star + M_L_star_err * randn()
    Ms = L * Y / M2MSUN
	
    log_v_circ_max = lMs_to_lVc(log10(Ms)) + lMs_to_lVc_err * randn()
	v_circ_max = 10 .^ log_v_circ_max
    r_circ_max = LilGuys.Ludlow.solve_rmax(v_circ_max, lc_err * randn())

    halo = NFW(r_circ_max=r_circ_max, v_circ_max=v_circ_max)

	M200 = LilGuys.calc_M200(halo) .* M2MSUN

	σv_kms = calc_σv_star_mean(halo) * V2KMS
	
    return Dict(
        :MV => MV,
        :L => L,
        :Y => Y,
        :Ms => Ms * M2MSUN,
        :log_v_circ_max => log_v_circ_max .+ log10(V2KMS),
        :r_circ_max => r_circ_max,
        :v_circ_max => v_circ_max * V2KMS,
        :M200 => M200,
        :log_M200 => log10.(M200),
        :c => halo.c,
        :σv => σv_kms
    )
end

# ╔═╡ bbd5e178-eddc-4e75-b4c5-252048b75fb4
function generate_samples(N_samples)
    # Create a dataframe to hold the samples
    samples = DataFrame(sample_halo())
    
    # Populate the dataframe using the sampling function in a loop
    for i in 1:N_samples
        halo_data = sample_halo()
        push!(samples, halo_data)
    end
    
    return samples
end

# ╔═╡ b665cd0f-8669-4520-974b-9ab9fba8ad15
samples = generate_samples(N_samples)

# ╔═╡ 5cdc46ec-6d1b-403f-a3a2-a15b2c82ee16
pairplot(samples[:, [:MV, :L, :Ms]])

# ╔═╡ 0af29495-0d78-4e32-bbf9-69587e8ae222
pairplot(samples[:, [:Ms, :log_M200, :c]])

# ╔═╡ 4e1290b5-171e-4715-a87b-28a2cfcb4325
samples.Ms * 1e4

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
describe(samples.Ms) .* M2MSUN ./ 1e6

# ╔═╡ 0ba3f10c-46ad-48d5-9b15-eb67e9e505ed
describe(samples.M200)

# ╔═╡ ecade01b-f703-4780-bc5d-ccc4c448b676
describe(samples.c)

# ╔═╡ b85256a1-786f-4dee-a6f1-f55406c3b18e
md"""
# Self-consistency calculator
"""

# ╔═╡ 76003206-050b-4913-9ab3-d4c0c2dd05f8
fig_dir = "./figures/"

# ╔═╡ 7fd48721-749d-4e99-91cc-5ffde830487d
let
	fig = pairplot(samples[:, [:Ms, :log_v_circ_max, :r_circ_max]], labels=Dict(
	:Ms2 => L"$M_\star$", 
	:V_max => L"$v_\textrm{circ, max}$ / km\,s$^{-1}$", 
	:r_max => L"$r_\textrm{circ, max}$ / kpc", ))

	save(joinpath(fig_dir, "v_max_r_max_mcmc.pdf"), fig)

	fig
end

# ╔═╡ 82c2844a-87ad-4a37-b8e0-6c11d30ec7c4
halos_ex = Dict(
	:mean => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 5.9),
	:heavy => NFW(v_circ_max = 42 / V2KMS, r_circ_max = 5.9),
	:compact => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 3.2),
	:middle => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 4.2),
)

# ╔═╡ 9ad6c9af-6922-46b9-919e-1c9a5efb2b3e


# ╔═╡ 8a06c0ad-2b48-48ba-a542-85926e164712
labels_ex = [:mean, :heavy, :compact, :middle]

# ╔═╡ c49ac57e-8e8d-4ed6-ad35-be400863f6b4
begin 
	V_max_in = 0.10
	# fiducial is 31.3
	r_max_in = 4.2

	V_max_in * V2KMS
end

# ╔═╡ afb4aff2-8ecc-40e1-a43a-e42cbb9536f6
31 / V2KMS

# ╔═╡ 0df6a87c-74b0-4713-a5e8-d73efbcb0b26
LilGuys.Ludlow.solve_rmax.(20 / V2KMS, 0.2)

# ╔═╡ d0173c19-507f-4a68-8a2a-a4f5d084a7e1
NFW_small = LilGuys.TruncNFW(r_circ_max=1.85, v_circ_max = 20 / V2KMS, trunc=10)

# ╔═╡ eee22cf7-79c6-427e-9542-ef7a042e9787
LilGuys.calc_M(NFW_small, 10000)

# ╔═╡ a0be8f75-a69b-447c-b616-442cf423fffd
LilGuys.Ludlow.solve_rmax.(31 / V2KMS, 0.2)

# ╔═╡ 6c50c1c9-13e0-4837-af21-f4d5343af609
LilGuys.Ludlow.solve_rmax.(47 / V2KMS, 0.2)

# ╔═╡ 5e71b022-bd95-4c3c-8930-51100fb9ab1c
r_max_exp = LilGuys.Ludlow.solve_rmax(V_max_in, 0)

# ╔═╡ f335d286-04d3-4248-b9bd-4bb6d8e82e33
r_max_exp, (LilGuys.Ludlow.solve_rmax.(V_max_in, 0.1), LilGuys.Ludlow.solve_rmax.(V_max_in, -0.1)) .- r_max_exp

# ╔═╡ b5a53e42-ef26-47b0-82f9-d404a4d3a544
log10(1 + 3.726 / r_max_exp)

# ╔═╡ 21acda7b-a116-4bc1-af7d-e715e2b32b86


# ╔═╡ 8196e6b2-8355-438c-abc1-ffce2e29b8f2
let
	fig, ax = FigAxis(
		ylabel=L"$\log\,v_\textrm{circ, max}$ / km\,s$^{-1}$",
		xlabel=L"$\log\,r_\textrm{circ, max}$ / kpc",
		limits=(0.2, 1.5, 1.2, 1.9)
	)

	for label in labels_ex
		y = LilGuys.calc_v_circ_max(halos_ex[label]) * V2KMS
		x = LilGuys.calc_r_circ_max(halos_ex[label]) 
		scatter!(log10(x), log10(y), label=string(label))
	end


	
	v = log10.(Vc_mean * V2KMS)
	lines!(log10.(Rc_mean), v, color=:grey, label="Ludlow+16")
	
	xl = log10.(LilGuys.Ludlow.solve_rmax.(Vc_mean, 0.1))
	xh =  log10.(LilGuys.Ludlow.solve_rmax.(Vc_mean, -0.1))
	x = [xl; reverse(xh)]
	y = [v; reverse(v)]
	poly!(x, y, color=(:grey, 0.2))

	hspan!(1.49 - 0.1, 1.49+0.1, color=(COLORS[2], 0.1))
	hlines!(1.49, color=(COLORS[2], 0.5), label="Fattahi+18")

	x = log10.(samples.r_circ_max )
	y = log10.(samples.v_circ_max)
	k = kde([x y])

	contour!(k)

	
	axislegend(position=:lt, title="halo")

	fig
end

# ╔═╡ bbee444e-079b-4208-9faf-0a7fe5f81455
let
	fig, ax = FigAxis(
		xlabel=L"$v_\textrm{circ, max}$ / km\,s$^{-1}$",
		ylabel=L"$r_\textrm{circ, max}$ / kpc",
		limits=((00, 70), (0, 20)),
	)

	x = samples.v_circ_max 
	y = samples.r_circ_max
	k = kde([x y])

	h = hexbin!(x, y, bins=50)

	#c = contour!(k, linewidth=6)

	scatter!(V_max_in * V2KMS, r_max_in)

	lines!(Vc_mean * V2KMS, Rc_mean)


	Colorbar(fig[1, 2], h, label="counts")

	save(joinpath(fig_dir, "r_max_v_max_in.pdf"), fig)
	
	fig
end

# ╔═╡ 41283b0b-6563-4b2b-b978-4e65f32c8240
calc_σv_star_mean(LilGuys.TruncNFW(r_circ_max=4.9, v_circ_max= 30/V2KMS, trunc=100)) * V2KMS

# ╔═╡ fb835e20-957e-4866-8ba6-2e64df38f68e
(median(samples.r_circ_max))

# ╔═╡ c33fd86a-ec88-4dc0-a883-519755f675c3
(median(samples.v_circ_max))

# ╔═╡ 34ef3a12-c2f9-4220-944d-f348efab86df
10 ^ 0.5

# ╔═╡ 7b3777cc-e2f7-4032-86a9-774b03a8a141
let
	fig, ax = FigAxis(
		xlabel=L"$\log r_\textrm{circ, max}$ / km\,s$^{-1}$",
		ylabel=L"$\log v_\textrm{circ, max}$ / kpc",
	)

	x = log10.(samples.r_circ_max )
	y = log10.(samples.v_circ_max )
	z = samples.σv

	σ_thresh = 9.0

	#c = contour!(k, linewidth=6)
	filt = z .>= σ_thresh
	h = scatter!(x[filt], y[filt], color=z[filt], colorrange=extrema(z[filt]))
	# scatter!(x[.!filt], y[.!filt], color=z[.!filt], colorrange=extrema(z[filt]), markersize=2)

	lines!(log10.(Rc_mean), log10.(Vc_mean * V2KMS))


	Colorbar(fig[1, 2], h, label=L"$\sigma_v$ / km\,s$^{-1}$")

	save(joinpath(fig_dir, "r_max_v_max_in.pdf"), fig)
	
	fig
end

# ╔═╡ 6baff9d8-a96d-4d1b-898f-089003459c19
let
	v_1 = LinRange(18, 100, 1000) ./V2KMS

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

	x = samples.v_circ_max
	y = samples.Ms
	k = kde([x y])

	c = contour!(k)
	
	
	lines!(v_1 * V2KMS, M_s_from_vel.(v_1))

	lines!(Vc_mean * V2KMS, Ms_mean)

	vlines!(V_max_in * V2KMS)
	fig

end

# ╔═╡ 2b1d7ec9-b3f7-4230-9f5a-80ece86f709d
halo_in = NFW(v_circ_max=V_max_in, r_circ_max=r_max_in)

# ╔═╡ db89cc22-f770-4f10-b3c0-3bc3298fb6a3
calc_σv_star_mean(halo_in) * V2KMS

# ╔═╡ ecda5f20-27cd-41a8-8545-9f3a6b91a80e
LilGuys.calc_M200(halo_in)

# ╔═╡ ab5bc36d-d203-45be-8e7d-d8c5c9473388
LilGuys.calc_R200(halo_in)

# ╔═╡ 84892a5a-f816-450b-9757-e4135a40aebc
LilGuys.calc_v_circ_max(halo_in) * V2KMS

# ╔═╡ cb42028b-eeb0-48a7-b23d-5abe792eb4f0
md"""
# calculating H
"""

# ╔═╡ ba292975-8eea-4c3a-8b9c-23ff1b60c7ab
Np = 1e7

# ╔═╡ e582e62c-d851-4dd0-95d4-82fd1d97c26c
function calc_h(halo, N=Np)
	return 4 * LilGuys.calc_R200(halo) / sqrt(N) / sqrt(10)
end

# ╔═╡ c3e1326c-772b-4d8c-aabe-a997b77bede4
calc_h(halo_in)

# ╔═╡ a75e163e-3489-4f82-901b-c511d1e44395
m = LilGuys.calc_M200(halo_in) / Np

# ╔═╡ 8b7f61d5-dcef-4889-8aab-a85b0ccdaf02
let
	fig, ax = FigAxis(
		xlabel = "log r / kpc",
		ylabel = "N ( r < r)"
	)

	lx = LinRange(-2, 2, 100)
	
	y = LilGuys.calc_M.(halo_in, 10 .^ lx) ./ m
	

	lines!(lx, log10.(y))
	vlines!(log10.(r_max_in), label="rmax")
	vlines!(log10(calc_h(halo_in)))
	axislegend()
	fig
end

# ╔═╡ 93715b38-1790-4b1e-a551-17d69b68876f
LilGuys.calc_M(halo_in, LilGuys.calc_R200(halo_in))

# ╔═╡ 896ccd5d-8e92-4c26-944d-04d42550caef
LilGuys.calc_R200(halo_in)

# ╔═╡ 2ae19d4f-659f-4985-8fc2-98d94cc52730
4^2 * LilGuys.G * m / calc_h(halo_in)^2

# ╔═╡ b97bdec3-2dbe-4a1d-bad9-4116b2d2e614
LilGuys.G * LilGuys.calc_M200(halo_in) / LilGuys.calc_R200(halo_in)^2

# ╔═╡ Cell order:
# ╟─07a710d8-0ae4-4d9f-9759-002750730010
# ╟─33f433ec-e94e-41fd-a808-6254bf4d34ce
# ╟─45ac95e4-a9a1-4f8f-9523-8f670f076ae5
# ╠═2c90aad4-87e7-11ee-2fbd-05f67c7e0343
# ╠═60eba328-c047-4b59-8b8c-8b24a4884541
# ╠═f4269d1a-940d-41da-9183-5130978f7d7b
# ╠═ef991557-7bf7-4add-b29f-f16187297e46
# ╠═21619de9-8960-46fe-8a0c-808776b33c6d
# ╠═a334b0a2-b408-4c81-a074-67cb1bcb58be
# ╠═24b4184d-7f69-47fc-a758-50f21840283a
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
# ╟─effb14d6-56ad-4baa-a484-83820aa72bb4
# ╠═bab1a32c-8d35-4958-af84-4b239175f602
# ╟─f677644c-ff10-4f0e-b6e3-26622ab92e43
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
# ╠═60c1f90f-6604-464e-bc77-d182a8cc8785
# ╠═4f971234-565a-43a4-a4bf-94a9415e6d7e
# ╠═4bb1dcef-04e2-4aa5-aa8f-d76a21502e01
# ╠═023f12d0-317a-46db-b5c5-81a7d2216c55
# ╠═be32e183-a0c4-4261-a3a8-3486db0734b7
# ╠═e6995206-ebf2-40aa-ae7e-cc101e07f8ac
# ╠═8109003b-60cf-458a-a43d-9cd084d27344
# ╠═db6f7825-4f9a-48ef-8700-6f4434032751
# ╠═ff9ea455-d618-4ebf-b9e9-5d92e23b4f37
# ╠═0541d9c1-c6fb-4bb6-8541-7529a65e446d
# ╠═b9d918b4-3c03-474c-8659-2ccc601c487d
# ╠═da3746bb-db05-4981-9f77-fc92238f61ed
# ╠═07cb7647-be33-4a16-8382-d44302f9e6d3
# ╠═02c22085-8e4d-4ae8-a1d7-ac746dc5b42f
# ╠═c0f9d7f0-7983-4a44-b85c-15534623d534
# ╟─122794fb-82ab-4f3c-a458-29462fbab1fd
# ╠═5dc54372-14c4-424b-8419-55d639f00baa
# ╠═28e09f56-7a1c-4fc9-95da-36ef0ed75a74
# ╠═667937b2-006f-4429-b134-91934adb53c8
# ╠═e65a951d-21ce-45fe-837c-825c782e76d2
# ╠═bbd5e178-eddc-4e75-b4c5-252048b75fb4
# ╠═b665cd0f-8669-4520-974b-9ab9fba8ad15
# ╠═5cdc46ec-6d1b-403f-a3a2-a15b2c82ee16
# ╠═0af29495-0d78-4e32-bbf9-69587e8ae222
# ╠═0ad2bc38-9145-4e3a-8a07-2db6dad4404d
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
# ╠═76003206-050b-4913-9ab3-d4c0c2dd05f8
# ╠═82c2844a-87ad-4a37-b8e0-6c11d30ec7c4
# ╠═9ad6c9af-6922-46b9-919e-1c9a5efb2b3e
# ╠═8a06c0ad-2b48-48ba-a542-85926e164712
# ╠═c49ac57e-8e8d-4ed6-ad35-be400863f6b4
# ╠═afb4aff2-8ecc-40e1-a43a-e42cbb9536f6
# ╠═f335d286-04d3-4248-b9bd-4bb6d8e82e33
# ╠═0df6a87c-74b0-4713-a5e8-d73efbcb0b26
# ╠═d0173c19-507f-4a68-8a2a-a4f5d084a7e1
# ╠═eee22cf7-79c6-427e-9542-ef7a042e9787
# ╠═a0be8f75-a69b-447c-b616-442cf423fffd
# ╠═6c50c1c9-13e0-4837-af21-f4d5343af609
# ╠═b5a53e42-ef26-47b0-82f9-d404a4d3a544
# ╠═5e71b022-bd95-4c3c-8930-51100fb9ab1c
# ╠═db89cc22-f770-4f10-b3c0-3bc3298fb6a3
# ╠═21acda7b-a116-4bc1-af7d-e715e2b32b86
# ╠═8196e6b2-8355-438c-abc1-ffce2e29b8f2
# ╠═bbee444e-079b-4208-9faf-0a7fe5f81455
# ╠═41283b0b-6563-4b2b-b978-4e65f32c8240
# ╠═fb835e20-957e-4866-8ba6-2e64df38f68e
# ╠═c33fd86a-ec88-4dc0-a883-519755f675c3
# ╠═34ef3a12-c2f9-4220-944d-f348efab86df
# ╠═7b3777cc-e2f7-4032-86a9-774b03a8a141
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
