### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ ed4b328a-59a1-11ef-2a75-092fdb7659b8
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya
	using LilGuys
end

# ╔═╡ 9d548f5a-223c-4827-b2d3-8361d0ced243
using DataFrames: rename!

# ╔═╡ 73f523e3-39eb-49ea-92d1-39e9d5e91d16
include("../load_mist.jl")

# ╔═╡ b653409b-fa89-4399-98d2-09c794aa88dc
md"""
# A comparison of J+24 samples

In this notebook, I compare the different samples from Jensen et al. (2024), crossmatch the probabilities to estimate the uncertanties with different assumptions, and validate all columns Jax added.

"""

# ╔═╡ 0a497022-a1c0-497c-ae34-aee52cba4679
md"""
## Short summary of Jaclyn et al. 2024
 Bayesian probability estmates using Gaia data with components for the CMD, proper motions, and spatial distribution


Assumed: literature position angle, ellipticity, and scale radius.

The likelyhood is the product of the likelihood for all three components. The reported parameters are the central proper motions, scale radius, satalite fraction, and relative aplitude of the second exponential component (if fit). If using elliptical radii, 


"""

# ╔═╡ 7db0c811-a50f-434b-9d89-05a82a6d265c
md"""
### Spatial likelihood
The spatial likelihood function is 

``
{\cal L}_S \propto \exp(-r/r_e) + B\,\exp(-r/r_s)
``

The MW component is assumed to be constant.

"""

# ╔═╡ 6e7743d1-1399-483b-9d3a-a8a4c6afd721
md"""
### CMD
The CMD likelihood is based on a likelihood map created by
- 12Gyr Padova isochrone
- convolved with FWHM 0.1mag + mean gaia uncertainty for color uncertainty
- convolved with uncertainty in distance modulus
- A constant likelihood value in G mag ranging from B-R=-0.25 to the RGB 
- HB uncertainty of 0.1 mag in G plus the mean error in B-R (in quadrature)

MW likelihood is based on emperical distribution. Unclear where


### Proper motion
5$\sigma$ consistency with dwarf proper motion derived from single component case.

"""

# ╔═╡ 951dc3e9-496b-450e-9b0a-1862af2e8626
md"""
## Probabilities
The probability a star belongs to the satallite population (and not the background) is given b 

``
P_{\rm sat} = \frac{f_{\rm sat}\,{\cal L}_{\rm sat}}{f_{\rm sat}\,{\cal L}_{\rm sat} + (1-f_{\rm sat}) {\mathcal L}_{\rm MW}}
``


"""

# ╔═╡ 4d967c49-4074-4e07-b175-841dde78f609
md"""
# Setup :)
"""

# ╔═╡ 2dacbd3e-b027-44cd-bb39-d3bfcf99366f
import Statistics: mean

# ╔═╡ 27d01941-3d40-4447-847f-e21191a20dd1
md"""
## Added columnes
- `xi`: ra-wise tangent coordinate (deg)
- `eta`: dec-wise tangent coordinate (deg)
- `r_ell` elliptical radius
- `dG` uncertainty on gaia_phot_g_mean_mag in dex ?
- `dBp`
- `dRP`
- `PSAT` The satalite probability considering spatial, CMD, and PM information
- `PSAT_S` The satalite probability only considering spatial likelihoods
- `PSAT_CMD` The satalite probability only considering CMD likelihoods
- `PSAT_PM` The satalite probability only considering proper motion likelihoods
- `PSAT_NOSPACE` The satalite probability only considering PM and CMD (no spatial) likelihoods
- `F_CPAR` flag for consistent parallax with galaxy
- `F_NOTANAN` flag for no NANs in COLUMNS?
- `F_ASTROMETRIC` flag for pmra/pmdec less than 10?
- `F_FLUXEXCESS` formula in ...
- `F_INGRID` inside CMD grid (i.e.)
- `F_BEST` stars which satisfy the above flags
- `L_SAT` The total likelihood for satalite membership (CMD * PM * S)
- `L_BKD` The total likelihood for MW membership
- `L_S_SAT` The spatial likelihood for the satalite
- `L_S_SAT_Inner`
- `L_S_SAT_Outer`
- `L_S_BKD`
- `L_PM_SAT`
- `L_PM_BKD`
- `L_CMD_SAT`
- `L_CMD_SAT_HB`
- `L_CMD_SAT_RGB`
- `L_CMD_BKD`

"""

# ╔═╡ 0d9da520-3d03-4d9a-8bb0-6d698e819abc
md"""
## Load samples
"""

# ╔═╡ d452a47a-22a2-4523-a7cf-f7b668cc7dda
data_dir = "../../data"

# ╔═╡ 330a4e01-59e0-4eb6-9900-d23db1159dd5
scl_ell = LilGuys.load_fits(joinpath(data_dir, "Sculptor.GAIASOURCE.RUWE.VELS.PROB.2-comp-ell.fits"))

# ╔═╡ b37ad796-1695-4665-b592-b68a99356907
scl_circ = LilGuys.load_fits(joinpath(data_dir, "Sculptor.GAIASOURCE.RUWE.VELS.PROB.2-comp-circ.fits"))

# ╔═╡ 39be2272-202a-4bb9-95cb-a9e584589d97
scl_1comp = LilGuys.load_fits(joinpath(data_dir, "Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits"))

# ╔═╡ 08cf3854-66e9-4358-8925-e069f59f65e6
scl_gaia = LilGuys.load_fits("sculptor/gaia_4deg_cen.fits")

# ╔═╡ 8a8aadb5-b84f-4ed0-8942-2c407ccc2517
setdiff(scl_1comp.source_id, scl_ell.source_id)

# ╔═╡ 9ec35690-9048-4f7f-954f-33a474d33bf3
setdiff(scl_1comp.source_id, scl_circ.source_id)

# ╔═╡ f66e116c-ad7a-4804-a3e9-099b7743fafa
md"""
# Double checking calculated properties
"""

# ╔═╡ 8ac66e16-932c-466f-bb56-0befe5681c04
df = scl_ell;

# ╔═╡ e8dc3faf-e3e4-44e6-988d-bbfdc9098d07
f_sat = 0.46

# ╔═╡ 11a6ee78-eeb7-49a6-a0d4-91245f855d2d
ra = 15.039

# ╔═╡ 263af8ce-56df-44bb-bf61-0cf268be178a
dec = -33.709

# ╔═╡ 2ce486a8-25d6-4688-b734-f61ebda7faa6
ell = 0.37

# ╔═╡ bf765255-c016-40f1-a409-2b0b9e5e94a4
PA = 94

# ╔═╡ 2bcff33a-eecb-4739-a3ff-3f06a38eeb81
xi, eta = LilGuys.to_tangent(df.ra, df.dec, ra, dec)

# ╔═╡ 69510524-dd91-4199-ab9e-b9252ec1e1b5
md"""
## Dataset
"""

# ╔═╡ b28733d4-7335-4157-b0ab-173aa6c56505
scl_gaia[:, "r_ell"] = LilGuys.calc_r_ell_sky(scl_gaia.ra, scl_gaia.dec, 0, PA; centre=(ra, dec))

# ╔═╡ 1e1f45d7-2b42-4ab3-b895-2e9fca8f4bc1
scl_gaia_filt = scl_gaia[scl_gaia.r_ell .< 120, :]

# ╔═╡ 8a9f4bcc-1fda-4204-8d80-918825bdcb49
missing_ids = setdiff(scl_ell.source_id, scl_gaia_filt.source_id)

# ╔═╡ 2f34b5ac-8fd6-4079-ad3f-e5de14187129
sqrt.(xi .^ 2 .+ eta .^ 2)[scl_ell.source_id .∈ [missing_ids], :]

# ╔═╡ cbac00f2-9c80-4f45-a692-1d37fe0f29d0
md"""
the difference in the sets is due to rounding errors (all at about 2 degrees)
"""

# ╔═╡ 29f3aff7-cf7c-4507-9069-78e29aa1ea9d
md"""
### Tangent plane coordinates and position angle...
"""

# ╔═╡ dd25169a-6391-44c0-b309-5f2a9a60fee8
scatter(xi, df.xi .- xi;
	axis=(; xlabel="xi (degrees)", ylabel="xi jax - xi me")
)

# ╔═╡ 6a58d6b1-5320-4515-8dad-18481cdca74f
scatter(eta, df.eta .- eta,
	axis=(; xlabel="eta (degrees)", ylabel="eta jax - xi me")
)

# ╔═╡ 52bd977c-cec6-4483-bd0c-470a3a8f071f
r_ell = LilGuys.calc_r_ell(xi, eta, ell, PA) *60

# ╔═╡ f5547f0a-b4ad-4d98-97a6-953a7be11950
scatter(r_ell, df.r_ell .- r_ell / 12.33 / sqrt(1-ell);
	axis=(;
		xlabel = "elliptical radius (rh)",
		ylabel= "rell me - jax (rh)",
	)
)

# ╔═╡ 737a2243-f7da-4441-9c2f-2472c6b65faf
md"""
### Likelihood functions
"""

# ╔═╡ df988526-ffa2-4c92-8a90-89f13204a115
md"""
I do not think that it is worth rederiving the CMD likelihood, so here is a figure which matches with the heatmap in the paper.
"""

# ╔═╡ 4243c523-edfd-40a1-9e91-a2563c0ea3a0
isocmd = ISOCMD("../../data/MIST_v1.2_vvcrit0.4_UBVRIplus/MIST_v1.2_feh_m2.00_afe_p0.0_vvcrit0.4_UBVRIplus.iso.cmd")

# ╔═╡ 93d499d2-4744-4613-878c-de4e425f7bfd
isocmd[9]

# ╔═╡ 502ad3a7-2589-4346-80bb-c0184356de9f
begin
	iso = isocmd[10.]
	iso[!, :G] = iso.Gaia_G_EDR3
	iso[!, :b_r] = iso.Gaia_BP_EDR3 .- iso.Gaia_RP_EDR3
	iso = iso[iso.phase .<= 3., :]
end

# ╔═╡ bfc412eb-ccae-4ae3-942a-ef8f8ea22fca
A_v = 0.0

# ╔═╡ b2ac5282-d503-4c97-b7f3-8f92e043bcaf
DM = 20.2

# ╔═╡ a66757cf-a936-436c-8315-045fffe4fca3
let
	fig, ax = FigAxis(
		xlabel = "bp - rp",
		ylabel = "G",
		yreversed=true,
		limits=(-0.5, 2.5, 15, 22)
	)
	
	p = scatter!(df.bp_rp, df.phot_g_mean_mag, color=(df.L_CMD_SAT), markersize=3)

	Colorbar(fig[1, 2], p, label="CMD likelihood")


	lines!(iso.b_r .- A_v, iso.G .+ DM)

	fig
end

# ╔═╡ 30c7bf30-2d50-4527-981e-d28acde0df5e
let
	fig, ax = FigAxis(
		xlabel = "delta bp - rp",
		ylabel = "G",
		yreversed=true,
		limits=(-0.01, 1, 15, 22)
	)
	
	p = scatter!(df.dBP, df.phot_g_mean_mag, color=(df.L_CMD_SAT), markersize=3)

	Colorbar(fig[1, 2], p, label="CMD likelihood")


	fig
end

# ╔═╡ 04903a00-c480-4ba5-8852-9b7c100e285a
lsat = scl_ell.L_CMD_SAT .* scl_ell.L_PM_SAT

# ╔═╡ 8550a1f6-8f54-43ee-9c9b-51e4ddd093b5
lbkd = scl_ell.L_CMD_BKD .* scl_ell.L_PM_BKD

# ╔═╡ da5b85b2-6406-4ef0-9c53-15555a795fd7
@. lsat * f_sat / (lsat * f_sat + (1-f_sat) * lbkd)

# ╔═╡ 3dc4d608-30b3-4c11-83e9-f5ae31df4786
let
	fig, ax = FigAxis(
		xlabel = "xi",
		ylabel = "eta",
	)

	df = scl_circ
	p = scatter!(df.xi, df.eta, color=df.L_S_BKD, 
		markersize=3, colorscale=log10, colorrange=(1e-1, 1e1)
	)

	Colorbar(fig[1, 2], p, label="spatial likelihood for MW", )

	fig
end

# ╔═╡ 59038c14-c565-4be4-baf5-021ce89fea1a
let
	fig, ax = FigAxis(
		xlabel = "xi",
		ylabel = "eta",
	)

	df = scl_circ
	p = scatter!(df.xi, df.eta, color=df.L_S_SAT, 
		markersize=3, colorscale=log10, colorrange=(1e-1, 1e1)
	)

	Colorbar(fig[1, 2], p, label="spatial likelihood for satalite", )

	fig
end

# ╔═╡ 326619b2-01ca-4f60-8619-516dadc58b98
rs = 0.28 * 60 

# ╔═╡ 0e336f85-16ec-44f5-a664-35dc569a3af8
re = 10.09/ 1.68

# ╔═╡ e9516503-2ef9-4ac3-97f5-b1c705bdd579
B = 0.015

# ╔═╡ b0ec46b5-8fbe-4514-9505-df638c534354
1/re

# ╔═╡ ffc1a8a8-3851-43d7-b017-0c17ee991b6f
l_s_sat_model(r) =  exp(-r / re) + B * exp(-r / rs) 

# ╔═╡ cc1e4080-93a4-4616-bdfe-74f23a8b3bf2
let
	fig, ax = FigAxis(
		yscale=log10,
		limits=(0, 150, 0.1, 10)
	)
	
	scatter!(r_ell, l_s_sat_model.(r_ell) ./ df.L_S_SAT .* maximum(df.L_S_SAT))
	hlines!(1)

	fig
end

# ╔═╡ bb030f67-ab58-40e8-83b8-e2dd7dbba3b3
let
	fig, ax = FigAxis(
		xlabel = "r ell",
		ylabel = "L sat",
		yscale=log10,
	)
	
	scatter!(r_ell, df.L_S_SAT ./ maximum(df.L_S_SAT))
	scatter!(r_ell, l_s_sat_model.(r_ell))

	fig
end

# ╔═╡ 61d1658a-57be-4ff1-a8a3-4fcfa98924d4
let
	fig, ax = FigAxis(
		limits=(nothing, nothing, 0, 1)
	)
	
	scatter!(r_ell, (maximum(log.(df.L_S_SAT)) .- log.(df.L_S_SAT)) ./ r_ell)
	fig
end


# ╔═╡ 4df3119c-b230-479d-bcaa-d38109e8536c
let
	fig, ax = FigAxis(
		limits=(nothing, nothing, -1, 1)
	)
	
	scatter!(r_ell, scl_1comp.L_S_SAT ./ maximum(scl_1comp.L_S_SAT) .- exp.(-0.1664 * r_ell))

	scatter!(r_ell, scl_1comp.L_S_SAT ./ maximum(scl_1comp.L_S_SAT) .- l_s_sat_model.(r_ell))
	
	fig
end


# ╔═╡ 73dbaf23-156b-4d2f-a4e0-e4ed80913306
1/0.166428 * 1.68 

# ╔═╡ a585a4c1-40dc-49ec-92cc-217484d82ee1
(maximum(log.(scl_1comp.L_S_SAT)) .- log.(scl_1comp.L_S_SAT)) ./ r_ell

# ╔═╡ 7e945f83-36f6-47d8-ae2b-7f0e2c3b1f57
md"""
### Proper motions
"""

# ╔═╡ 75512420-ee17-4786-8adc-c239b1df1e7b
pmra =  0.099 #.+ 0.004

# ╔═╡ 29a3794e-0fe2-4613-8b41-c3ff02229926
pmdec = -0.15

# ╔═╡ 53517031-31bc-488e-96e4-cd55069320ff
dpm = 0.000

# ╔═╡ d0ce4a42-2253-479f-a221-493ce2b16185
pmra_err = dpm .+ df.pmra_error 

# ╔═╡ da80d829-28ec-4cd8-93b6-831dc10a38a7
pmdec_err = dpm .+ df.pmdec_error 

# ╔═╡ a43555aa-ce37-4a98-8d29-e154c41f6f1a
mean(df.pmra_pmdec_corr[isfinite.(df.pmra_pmdec_corr)])

# ╔═╡ 1753360f-d8ea-4056-99f9-0b629a371fab
pm_ρ = df.pmra_pmdec_corr

# ╔═╡ 7c8854af-73dd-4a42-aa95-ba907ca8cf61
pm_dx = @. (pmra - df.pmra) ./ pmra_err

# ╔═╡ bf247e8c-39f6-4355-a3ad-0870ea1a9803
pm_dy = @. (pmdec - df.pmdec) ./ pmdec_err

# ╔═╡ 828c7416-4e43-4fc9-975c-fc04cb0ecc77
pm_num = @. (pm_dx^2 + pm_dy^2 - 2*pm_ρ * pm_dx * pm_dy
)

# ╔═╡ 41958113-1ba5-4018-9386-dfb42a1f2dd8
n_sigma_pm = @. pm_num / (2*(1-pm_ρ^2))

# ╔═╡ 41a60b57-e722-4315-9caf-d823da4d45c9
df.L_PM_SAT

# ╔═╡ 0e006ab4-e87c-41c7-a43c-79022fbe9df2
pm_factor = 0.5

# ╔═╡ 29386139-c2a6-4ecc-8d95-35a5f2fd78f9
@doc raw"""
	bivariate_normal(x, y, μx, μy, σx, σy, ρ)

A bivariate normal distribution on input vectors x, y assuming a mean (μx, μy), standard deviations of (σx, σy), and a correlation ρ ∈ [-1, 1].

``
\frac{1}{2\pi \sigma_x \sigma_y \sqrt{1 - \rho^2}} \exp\left[-\frac{1}{2(1-\rho^2)}\left(\frac{(x-\mu_x)^2}{\sigma_x^2} + \frac{(y-\mu_y)^2}{\sigma_y^2} - 2\rho \frac{(x-\mu_x)(y-\mu_y)}{\sigma_x\sigma_y}\right)\right]
``

"""
function bivariate_normal(x, y, μx, μy, σx, σy, ρ)
	A = 1 / (2π * σx * σy * √(1 - ρ^2))

	zx = (x - μx) / σx
	zy = (y - μy) / σy

	return A * exp(-1/(2*(1 - ρ^2)) * (
		zx^2 + zy^2 - 2ρ * zx * zy
	))
end

# ╔═╡ 8c90b632-8c17-4468-9d9c-e19d8f86fa58
"""
	⊕(x, y)

Add x and y in quadrature
"""
function ⊕(x, y)
	return sqrt(x^2 + y^2)
end

# ╔═╡ e965e22c-ef55-410f-8cf9-28403f01239e
L_pm = bivariate_normal.(df.pmra, df.pmdec, pmra, pmdec, df.pmra_error .⊕ dpm, df.pmdec_error .⊕ dpm, df.pmra_pmdec_corr)

# ╔═╡ 81aaaaef-8fff-41b9-981b-cd4b08dc62cd
scatter(df.L_PM_SAT, L_pm;
	alpha=0.1,
	markersize=5, 
	axis=(; limits=(-0.1, 4, -0.1, 4))
)

# ╔═╡ 4accc157-bf98-42b9-b8bf-8395bff5f49f
sum(.! isnan.(df.L_PM_SAT))

# ╔═╡ 54e97bc1-593b-4b74-8f77-433f5cad7b8b
let
	fig, ax = FigAxis(
		xlabel = "pmra",
		ylabel = "pmdec",
		limits=(-5, 5, -5, 5)
	)

	df = scl_ell
	filt = df.L_PM_SAT .> 0e-10
	# filt .&= isfinite.(scl_ell.L_PM_SAT)
	p = scatter!(df.pmra[filt], df.pmdec[filt], color=df.L_PM_SAT[filt], 
		markersize=3, colorscale=log10, colorrange=(1e-2, 1e1)
	)

	Colorbar(fig[1, 2], p, label="pm likelihood", )

	fig
end

# ╔═╡ 34dc5d88-16ce-49e7-912f-0e6f36c81bd3
pmra, pmdec

# ╔═╡ c6dda910-7ade-4668-9592-899264248339
let
	fig, ax = FigAxis(
		xlabel = "pmra",
		ylabel = "pmdec",
		limits=(-5, 5, -5, 5)
	)
	
	p = scatter!(scl_ell.pmra, scl_ell.pmdec, color=L_pm, 
		markersize=3, colorscale=log10, colorrange=(1e-2, 1e1)
	)

	Colorbar(fig[1, 2], p, label="PM likelihood", )

	fig
end

# ╔═╡ ee79b752-fcdd-4742-af23-e40eb86edbe9
md"""
### Likelihood products
"""

# ╔═╡ 60deb437-7803-4ce3-b2eb-c2f3c78f4bc2
let
	fig, ax = FigAxis(limits=(1e-10, 1e3, nothing, nothing), xscale=log10)
	
	scatter!(scl_1comp.L_SAT, scl_1comp.L_CMD_SAT .*scl_1comp.L_PM_SAT .* scl_1comp.L_S_SAT ./ scl_1comp.L_SAT, alpha=0.1)

	fig
end

# ╔═╡ 43271a93-c879-425d-8896-dd3ab856a98b
let
	fig, ax = FigAxis(limits=(1e-10, 1e3, nothing, nothing), xscale=log10)
	
	scatter!(scl_1comp.L_BKD, scl_1comp.L_CMD_BKD .*scl_1comp.L_PM_BKD .* scl_1comp.L_S_BKD ./ scl_1comp.L_BKD, alpha=0.1)

	fig
end

# ╔═╡ 8cacb247-4a07-43d7-9344-156c283a9cf5
md"""
### Likelihood comparisons
"""

# ╔═╡ 620d34a4-b96d-46eb-a213-10527344bd15
for i in 1:3
	for j in ((1:3)[1:3 .> i])
		println("comparing $i, $j")

		dfs = [scl_circ, scl_ell, scl_1comp]
		names = ["2c. circ", "2c ell", "1c"]
		df1 = dfs[i]
		df2 = dfs[j]

		println("1 = ", names[i])
		println("2 = ", names[j])
		xl = names[i]
		yl = names[j]

		@info scatter(log10.(df1.L_CMD_SAT), df2.L_CMD_SAT ./ df1.L_CMD_SAT, alpha=0.1;
			axis=(; xlabel="L cmd $xl", ylabel="L cmd ratio $yl / $xl",
				limits=(-20, nothing, nothing, nothing)
			)
		)
	end
end

# ╔═╡ 88673e74-ac12-468f-84b2-6205d447089e
for i in 1:3
	for j in ((1:3)[1:3 .> i])
		println("comparing $i, $j")

		dfs = [scl_circ, scl_ell, scl_1comp]
		names = ["2c. circ", "2c ell", "1c"]
		df1 = dfs[i]
		df2 = dfs[j]

		println("1 = ", names[i])
		println("2 = ", names[j])
		xl = names[i]
		yl = names[j]

		@info scatter(log10.(df1.L_PM_SAT), df2.L_PM_SAT ./ df1.L_PM_SAT, alpha=0.1;
			axis=(; xlabel="L pm $xl", ylabel="L pm ratio $yl / $xl",
			limits=(-10, nothing, nothing, nothing)
			)
		)
	end
end

# ╔═╡ bfe46773-5b74-4bcf-a65e-57912d8c2326
for i in 1:3
	for j in ((1:3)[1:3 .> i])
		println("comparing $i, $j")

		dfs = [scl_circ, scl_ell, scl_1comp]
		names = ["2c. circ", "2c ell", "1c"]
		df1 = dfs[i]
		df2 = dfs[j]

		println("1 = ", names[i])
		println("2 = ", names[j])
		xl = names[i]
		yl = names[j]

		@info scatter(log10.(df1.L_S_SAT), df2.L_S_SAT ./ df1.L_S_SAT, alpha=0.1;
			axis=(; xlabel="L spatial $xl", ylabel="L spatial ratio $yl / $xl", limits=(-10, nothing, nothing, nothing))
		)
	end
end

# ╔═╡ 2a632338-79b7-4a81-aa8a-67888844a388
md"""
# Satalite probabilities
"""

# ╔═╡ 66ee7b4e-fbd7-4fb5-86ea-b20a339a8a7a
psat = @. scl_ell.L_SAT * f_sat / (scl_ell.L_SAT * f_sat + (1-f_sat) * scl_ell.L_BKD)

# ╔═╡ 6ce28a2b-2f66-494c-8419-bc78eab3ec9a
scl_ell.PSAT

# ╔═╡ 793c845e-92b6-4aba-8fdc-4bd88f65bc0a
mean(scl_ell.PSAT[.! isnan.(scl_ell.PSAT)])

# ╔═╡ a79814b3-3b9d-4e26-9001-9e84cbe67bec
for i in 1:3
	for j in ((1:3)[1:3 .> i])
		println("comparing $i, $j")

		dfs = [scl_circ, scl_ell, scl_1comp]
		names = ["2c. circ", "2c ell", "1c"]
		df1 = dfs[i]
		df2 = dfs[j]

		println("1 = ", names[i])
		println("2 = ", names[j])
		xl = names[i]
		yl = names[j]
		@info scatter(df1.PSAT, df2.PSAT, alpha=0.1;
			axis=(; aspect=DataAspect(), xlabel="PSAT $xl", ylabel="PSAT $yl")
		)

	end
end

# ╔═╡ ad6cf0cc-42c5-414d-bf07-27311d7716c2
scatter(scl_ell.r_ell, scl_circ.PSAT_S .- scl_ell.PSAT_S, color=scl_ell.PSAT)

# ╔═╡ 8486d28e-5ae8-40f2-9a93-ce81606f1496
md"""
# Combining properties
"""

# ╔═╡ 0a48ae68-348b-4494-9c81-d11c9bbe87a0
begin 
	df_out = copy(df)

	df_out[!, :L_S_SAT_1C] .= scl_1comp.L_S_SAT
	df_out[!, :PSAT_S_1C] .= scl_1comp.PSAT_S
	df_out[!, :PSAT_1C] .= scl_1comp.PSAT

	df_out[!, :L_S_SAT_CIRC] .= scl_circ.L_S_SAT
	df_out[!, :PSAT_S_CIRC] .= scl_circ.PSAT_S
	df_out[!, :PSAT_CIRC] .= scl_circ.PSAT

	rename!(df_out, :PSAT => :PSAT_ELL, 
		:L_S_SAT => :L_S_SAT_ELL,
		:PSAT_S => :PSAT_S_ELL,
	)

	df_out[!, :PSAT] = max.(df_out.PSAT_CIRC, df_out.PSAT_1C, df_out.PSAT_ELL)

end

# ╔═╡ dd1bbf56-4d9a-445c-9722-5541b63d227d
scatter(df_out.PSAT_CIRC, df_out.PSAT_1C)

# ╔═╡ d39d9fba-288f-48fc-8db5-eb10854cb4a7
let
	fig, ax = FigAxis()
	
	scatter!(df_out.PSAT_CIRC, df_out.PSAT_ELL)

	fig
end

# ╔═╡ b546c73a-823b-4f19-a9a3-d01946c4982f
@assert all(scl_1comp.source_id .== df.source_id)

# ╔═╡ ff29e7c1-5aa4-4295-b114-e38fc3a208c7
@assert all(scl_circ.source_id .== df.source_id)

# ╔═╡ 1d0ef71c-e584-4612-8685-09e6bb79423b
convert(Float64, 1)

# ╔═╡ 9b57518d-98b0-4cb6-a290-9b9ea22cefac
LilGuys.write_fits("$data_dir/j24_sculptor_all.fits", df_out, verbose=true)

# ╔═╡ Cell order:
# ╟─b653409b-fa89-4399-98d2-09c794aa88dc
# ╟─0a497022-a1c0-497c-ae34-aee52cba4679
# ╟─7db0c811-a50f-434b-9d89-05a82a6d265c
# ╟─6e7743d1-1399-483b-9d3a-a8a4c6afd721
# ╟─951dc3e9-496b-450e-9b0a-1862af2e8626
# ╠═4d967c49-4074-4e07-b175-841dde78f609
# ╠═ed4b328a-59a1-11ef-2a75-092fdb7659b8
# ╠═2dacbd3e-b027-44cd-bb39-d3bfcf99366f
# ╠═9d548f5a-223c-4827-b2d3-8361d0ced243
# ╟─27d01941-3d40-4447-847f-e21191a20dd1
# ╟─0d9da520-3d03-4d9a-8bb0-6d698e819abc
# ╠═d452a47a-22a2-4523-a7cf-f7b668cc7dda
# ╠═330a4e01-59e0-4eb6-9900-d23db1159dd5
# ╠═b37ad796-1695-4665-b592-b68a99356907
# ╠═39be2272-202a-4bb9-95cb-a9e584589d97
# ╠═08cf3854-66e9-4358-8925-e069f59f65e6
# ╠═8a8aadb5-b84f-4ed0-8942-2c407ccc2517
# ╠═9ec35690-9048-4f7f-954f-33a474d33bf3
# ╟─f66e116c-ad7a-4804-a3e9-099b7743fafa
# ╠═8ac66e16-932c-466f-bb56-0befe5681c04
# ╠═e8dc3faf-e3e4-44e6-988d-bbfdc9098d07
# ╠═11a6ee78-eeb7-49a6-a0d4-91245f855d2d
# ╠═263af8ce-56df-44bb-bf61-0cf268be178a
# ╠═2ce486a8-25d6-4688-b734-f61ebda7faa6
# ╠═bf765255-c016-40f1-a409-2b0b9e5e94a4
# ╠═2bcff33a-eecb-4739-a3ff-3f06a38eeb81
# ╟─69510524-dd91-4199-ab9e-b9252ec1e1b5
# ╠═b28733d4-7335-4157-b0ab-173aa6c56505
# ╠═1e1f45d7-2b42-4ab3-b895-2e9fca8f4bc1
# ╠═8a9f4bcc-1fda-4204-8d80-918825bdcb49
# ╠═2f34b5ac-8fd6-4079-ad3f-e5de14187129
# ╟─cbac00f2-9c80-4f45-a692-1d37fe0f29d0
# ╟─29f3aff7-cf7c-4507-9069-78e29aa1ea9d
# ╠═dd25169a-6391-44c0-b309-5f2a9a60fee8
# ╠═6a58d6b1-5320-4515-8dad-18481cdca74f
# ╠═52bd977c-cec6-4483-bd0c-470a3a8f071f
# ╠═f5547f0a-b4ad-4d98-97a6-953a7be11950
# ╟─737a2243-f7da-4441-9c2f-2472c6b65faf
# ╟─df988526-ffa2-4c92-8a90-89f13204a115
# ╠═73f523e3-39eb-49ea-92d1-39e9d5e91d16
# ╠═4243c523-edfd-40a1-9e91-a2563c0ea3a0
# ╠═93d499d2-4744-4613-878c-de4e425f7bfd
# ╠═502ad3a7-2589-4346-80bb-c0184356de9f
# ╠═bfc412eb-ccae-4ae3-942a-ef8f8ea22fca
# ╠═b2ac5282-d503-4c97-b7f3-8f92e043bcaf
# ╠═a66757cf-a936-436c-8315-045fffe4fca3
# ╠═30c7bf30-2d50-4527-981e-d28acde0df5e
# ╠═04903a00-c480-4ba5-8852-9b7c100e285a
# ╠═8550a1f6-8f54-43ee-9c9b-51e4ddd093b5
# ╠═da5b85b2-6406-4ef0-9c53-15555a795fd7
# ╠═3dc4d608-30b3-4c11-83e9-f5ae31df4786
# ╠═59038c14-c565-4be4-baf5-021ce89fea1a
# ╠═326619b2-01ca-4f60-8619-516dadc58b98
# ╠═0e336f85-16ec-44f5-a664-35dc569a3af8
# ╠═e9516503-2ef9-4ac3-97f5-b1c705bdd579
# ╠═b0ec46b5-8fbe-4514-9505-df638c534354
# ╠═ffc1a8a8-3851-43d7-b017-0c17ee991b6f
# ╠═cc1e4080-93a4-4616-bdfe-74f23a8b3bf2
# ╠═bb030f67-ab58-40e8-83b8-e2dd7dbba3b3
# ╠═61d1658a-57be-4ff1-a8a3-4fcfa98924d4
# ╠═4df3119c-b230-479d-bcaa-d38109e8536c
# ╠═73dbaf23-156b-4d2f-a4e0-e4ed80913306
# ╠═a585a4c1-40dc-49ec-92cc-217484d82ee1
# ╟─7e945f83-36f6-47d8-ae2b-7f0e2c3b1f57
# ╠═75512420-ee17-4786-8adc-c239b1df1e7b
# ╠═29a3794e-0fe2-4613-8b41-c3ff02229926
# ╠═53517031-31bc-488e-96e4-cd55069320ff
# ╠═d0ce4a42-2253-479f-a221-493ce2b16185
# ╠═da80d829-28ec-4cd8-93b6-831dc10a38a7
# ╠═a43555aa-ce37-4a98-8d29-e154c41f6f1a
# ╠═1753360f-d8ea-4056-99f9-0b629a371fab
# ╠═7c8854af-73dd-4a42-aa95-ba907ca8cf61
# ╠═bf247e8c-39f6-4355-a3ad-0870ea1a9803
# ╠═828c7416-4e43-4fc9-975c-fc04cb0ecc77
# ╠═41958113-1ba5-4018-9386-dfb42a1f2dd8
# ╠═41a60b57-e722-4315-9caf-d823da4d45c9
# ╠═0e006ab4-e87c-41c7-a43c-79022fbe9df2
# ╟─29386139-c2a6-4ecc-8d95-35a5f2fd78f9
# ╟─8c90b632-8c17-4468-9d9c-e19d8f86fa58
# ╠═e965e22c-ef55-410f-8cf9-28403f01239e
# ╠═81aaaaef-8fff-41b9-981b-cd4b08dc62cd
# ╠═4accc157-bf98-42b9-b8bf-8395bff5f49f
# ╠═54e97bc1-593b-4b74-8f77-433f5cad7b8b
# ╠═34dc5d88-16ce-49e7-912f-0e6f36c81bd3
# ╠═c6dda910-7ade-4668-9592-899264248339
# ╟─ee79b752-fcdd-4742-af23-e40eb86edbe9
# ╠═60deb437-7803-4ce3-b2eb-c2f3c78f4bc2
# ╠═43271a93-c879-425d-8896-dd3ab856a98b
# ╠═8cacb247-4a07-43d7-9344-156c283a9cf5
# ╠═620d34a4-b96d-46eb-a213-10527344bd15
# ╠═88673e74-ac12-468f-84b2-6205d447089e
# ╠═bfe46773-5b74-4bcf-a65e-57912d8c2326
# ╠═2a632338-79b7-4a81-aa8a-67888844a388
# ╠═66ee7b4e-fbd7-4fb5-86ea-b20a339a8a7a
# ╠═6ce28a2b-2f66-494c-8419-bc78eab3ec9a
# ╠═793c845e-92b6-4aba-8fdc-4bd88f65bc0a
# ╠═a79814b3-3b9d-4e26-9001-9e84cbe67bec
# ╠═ad6cf0cc-42c5-414d-bf07-27311d7716c2
# ╟─8486d28e-5ae8-40f2-9a93-ce81606f1496
# ╠═0a48ae68-348b-4494-9c81-d11c9bbe87a0
# ╠═dd1bbf56-4d9a-445c-9722-5541b63d227d
# ╠═d39d9fba-288f-48fc-8db5-eb10854cb4a7
# ╠═b546c73a-823b-4f19-a9a3-d01946c4982f
# ╠═ff29e7c1-5aa4-4295-b114-e38fc3a208c7
# ╠═1d0ef71c-e584-4612-8685-09e6bb79423b
# ╠═9b57518d-98b0-4cb6-a290-9b9ea22cefac
