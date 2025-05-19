### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 93838644-cad6-4df3-b554-208b7afeb3b8
using PyFITS

# ╔═╡ 72f1febc-c6ea-449a-8cec-cd0e49c4e20c
using DataFrames

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
using PairPlots

# ╔═╡ 02b266e0-20ee-4992-a83a-9cec72f0d975
using ArviZ

# ╔═╡ bd6dfd17-02ee-4855-be37-fecfdab6776f
using LilGuys; FIGDIR = "figures"

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ 34e43f4a-bcff-41cb-92c4-0c8d600fd053
import CSV

# ╔═╡ 2a422e88-fc0d-4a89-a841-42f3c5c8dace
import KernelDensity

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 2ab0018f-1628-4f5e-b7da-370eb20c00d0
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ 5ec475a1-14bb-40f6-856a-69fa9efe087a
⊕ = RVUtils.:⊕

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
σv = obs_properties["sigma_v"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership
"""

# ╔═╡ a162219f-df5a-41d8-bf54-927a355f6431
memb_stars = read_fits(joinpath(data_dir, "rv_members_all.fits"))

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ 55ce0f69-8a96-4bbb-a59f-ee6503624ea6
md"""
# Numbers
"""

# ╔═╡ f9d4eade-c648-4f20-8403-07be993fb8c1
sum(.!ismissing.(memb_stars.RV_gmos))

# ╔═╡ 3377f632-713d-4fec-84a9-b0211b02cb43
median(memb_stars.RV_err)

# ╔═╡ 31a6c2e4-538c-4adc-bbda-5043680b17f7
extrema(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f
md"""
## MCMC Priors
"""

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(RVUtils.model_vel_1c(memb_stars.RV, memb_stars.RV_err), Prior(), 1000)
)

# ╔═╡ ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
pairplot(prior_samples[:, [:μ, :σ]])

# ╔═╡ 74ad07df-f15c-436f-b390-ce95b27f7fab
let
	fig, ax = FigAxis(
		xgridvisible=false,
		ygridvisible=false,
		xlabel=xlabel,
		ylabel="density",
		title="priors"
	)
	
	RVUtils.plot_samples!(prior_samples, LinRange(-130, 130, 100))

	fig
end

# ╔═╡ c897719b-4af9-4c3e-824d-e8ebcd0080d7
md"""
# Simple
"""

# ╔═╡ f8d81628-eea6-457e-b000-92a6ed952f88
model_vz = RVUtils.model_vel_1c(memb_stars.vz, memb_stars.vz_err)

# ╔═╡ b7574fe6-c365-4c2b-879e-69c9932abb7a
samples_vz = sample(model_vz, NUTS(0.65), 1000)

# ╔═╡ fb882aa1-2722-471a-95e4-1ef486801461
ll_vz = pointwise_loglikelihoods(model_vz, samples_vz)

# ╔═╡ b93defe6-6034-4ab4-b83c-10eb289ca280
waic(ll_vz["x"]).estimates

# ╔═╡ 433d1d51-4b2f-4bc0-807a-813eb215231a
md"""
# Rell sigma
"""

# ╔═╡ b32a6521-dae5-4d2a-89ab-6e8b3bc43be3
model_Rell = RVUtils.model_vel_sigma_R(memb_stars.vz, memb_stars.vz_err, memb_stars.R_ell)

# ╔═╡ 93ee52ad-2c6d-4eac-85cd-ad2a88a70428
samples_Rell = sample(model_Rell,NUTS(0.65), 1000)

# ╔═╡ 3f8c8a37-7381-40c7-8826-482aa52d0f6c
ll_Rell = Turing.pointwise_loglikelihoods(model_Rell, samples_Rell)

# ╔═╡ a8096350-34cf-4cdd-b7fd-f9d1ac3d2d2f
ll_arr = reshape(reduce(hcat, values(ll_Rell)), 1, 1000)

# ╔═╡ 15e3b145-604e-47f8-a0c3-9466c34277df
idata_Rell = from_mcmcchains(samples_Rell, library="Turing", log_likelihood=Dict("vz" => ll_arr))

# ╔═╡ 3adfa6eb-02f4-4c62-a10f-fedf108a0cb6
missing_obs = similar(memb_stars.RV, Missing)

# ╔═╡ 54ae48c6-c2ac-409b-aaad-b31f73707b3d
model_Rell_predict = RVUtils.model_vel_sigma_R(missing_obs, missing_obs, missing_obs)

# ╔═╡ 464af51d-9287-474a-8bcd-f542e7fbe9d2
waic(ll_arr).estimates

# ╔═╡ d072218f-9da8-4770-9b00-b22c638869bc
ArviZ.summary(idata_Rell)

# ╔═╡ b983cbd3-a90f-4ed1-b0a7-fd21e903f81d
summary_Rell = RVUtils.summarize(samples_Rell)

# ╔═╡ 527ccde2-64db-40f5-9b77-8898a92b3b6a
df_Rell = DataFrame(samples_Rell)

# ╔═╡ 820b82a1-63a4-4a83-8d6e-5659ef75ce2d
pairplot(samples_Rell)

# ╔═╡ 924369b0-aff6-46e4-a0bf-c6766ff93cbf
median(df_Rell.μ) + Δv_gsr

# ╔═╡ 4f00f714-6b16-41c7-a3c8-f2743d6023bf


# ╔═╡ 137bbc9c-1e02-4f78-acd1-0e30cf982615
#bf_sigma_Rell = RVUtils.bayes_evidence(model_Rell, df_Rell, "dlσ_dlR")

# ╔═╡ 3b7cb3ea-3684-4df9-b9e6-919fa11e9ffa
kde_Rell = KernelDensity.kde(df_Rell.dlσ_dlR)

# ╔═╡ f8cc1ebb-617a-46a0-ab22-2f84d75da2ea
md"""
In the plot below, we just want to make sure that the KDE density estimate looks reasonable at zero
"""

# ╔═╡ 9158e64f-5a8d-48bf-8043-8e93d41bd626
let
	fig = Figure()
	ax = Axis(fig[1,1],
		yscale=log10, 
		yticks = Makie.automatic,
		ylabel = "density",
		xlabel = "dlσ/dlR"
	)
		
	
	lines!(kde_Rell)
	vlines!(0, color=:black)

	x = df_Rell.dlσ_dlR[df_Rell.dlσ_dlR .< quantile(df_Rell.dlσ_dlR, 0.001)]
	scatter!(x, 10^-6 .* (1 .+ rand(length(x))), color=COLORS[2], markersize=1)
	fig
end

# ╔═╡ 4a83224a-1f2a-4f7c-b670-6a694fc8826a
KernelDensity.default_bandwidth(df_Rell.dlσ_dlR)

# ╔═╡ d051635b-123c-45e2-94da-6511f179a2bb
CSV.write("processed/mcmc_samples_Rell_sigma.csv", df_Rell)

# ╔═╡ 062994bc-fb0c-4f08-b7f7-7cc6714bad1e
md"""
## Gradient
"""

# ╔═╡ 25d9799f-5e89-4ec3-995c-a7309efb798f
model_gradient = RVUtils.model_vel_gradient(memb_stars.vz, memb_stars.vz_err, memb_stars.xi, memb_stars.eta)

# ╔═╡ cb440be7-1f0e-4a20-b15a-82b1820d1ced
samples_gradient = sample(model_gradient, NUTS(0.65), 1000)

# ╔═╡ e88403cc-e502-4e83-a947-43a9bae89af5
ll_gradient = pointwise_loglikelihoods(model_gradient, samples_gradient)

# ╔═╡ c8949f29-0957-4ca7-bad4-b52eff44fb72
waic(ll_gradient["x"]).estimates

# ╔═╡ e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
@savefig "gradient_corner" pairplot(samples_gradient)

# ╔═╡ 184b4a5d-cbab-44b2-9620-bf928ad81d0e
df_gradient = let
	df = DataFrame(samples_gradient)
	df[:, :A]
	df[:, :B]
	df[!, :r_grad] = @. 60 * (df.A ⊕ df.B )
	df[!, :Θ_grad] = @. atand(df.A, df.B) 

	df
end

# ╔═╡ af21597b-b54a-4445-aa09-22dbad18b42a
waic(df_gradient.lp')

# ╔═╡ e79d9d82-84cd-4bc5-9ea8-be6b07cacf6d
 pairplot(df_gradient[:, [:μ, :σ, :r_grad, :Θ_grad]])

# ╔═╡ f2aa15e2-824b-4f0e-8ad0-42611abe06da
median(df_gradient.μ) + Δv_gsr

# ╔═╡ ad607c95-0aea-4d84-8f18-e5b929b9bfca
bf_gradient = RVUtils.bayes_evidence(model_gradient, df_gradient, ["A", "B"])

# ╔═╡ 4d7e7b96-74a0-4066-88cf-739c043c7f47
summary_gradient = RVUtils.summarize(samples_gradient)

# ╔═╡ 8555e608-53c1-40d3-b21e-413af8953c30
md"""
code below validates induced PM gradient (should be approx 2).
"""

# ╔═╡ 7a9d2e0f-4dcb-4b69-9f68-d43d6dde8bf2
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ 3195286d-d85f-43a3-aa25-dae2134f570b
xi_rot, eta_rot = lguys.to_orbit_coords(memb_stars.ra, memb_stars.dec, obs_properties["ra"], obs_properties["dec"], θ_m

# ╔═╡ 88f2918e-e126-420a-96a2-5746a8010f73
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ c48bb30d-1186-4940-b061-91f53e8335e1
vec_pm = LilGuys.pm2kms.([icrs0.pmra, icrs0.pmdec], icrs0.distance) /(180/π)

# ╔═╡ 0ca7dc1b-3b41-4089-9c89-20c6e48213ea
@savefig "v_gradient_derived" let
	fig = Figure()

	ax=Axis(fig[1,1];
		  limits=(-10, 10, -10, 10), 
		  aspect=DataAspect(),
		  xlabel = L"$\partial \,v_z / \partial\, \xi$ / km\,s$^{-1}$\,degree$^{-1}$",
		  ylabel = L"$\partial \,v_z / \partial\, \eta$ / km\,s$^{-1}$\,degree$^{-1}$",
		xreversed=true
		 )
	
	scatter!(60df_gradient.A, 60df_gradient.B, alpha=0.1, markersize=1, 

	   )

	scatter!(0, 0, color=:black)
	arrows!([0], [0], [vec_pm[1]], [vec_pm[2]])

	fig
end

# ╔═╡ af0d2050-b42e-4a7f-aabb-5b08d23381e9
lguys.transform(ICRS, lguys.GSR(ra=icrs0.ra + 2/vec_pm[1] *cos(icrs0.dec), dec=icrs0.dec + 2/vec_pm[2], distance=icrs0.distance, radial_velocity=0)).radial_velocity  .- Δv_gsr

# ╔═╡ 4a473039-79f0-4d77-aa0c-681e2fba4f4c
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ ed35eb68-74f7-4009-9b68-dfca2ea547af
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

# ╔═╡ f523cf48-82bf-4b20-9d7c-215bbe10a193
kde = KernelDensity.kde((df_gradient.A, df_gradient.B))

# ╔═╡ 2b88915c-7222-44be-a483-b967ea131b80
log(pdf(kde, 0, 0) ./ lguys.gaussian(0, 0., 0.1)^2)

# ╔═╡ 51c7c22e-405b-4a47-beb1-b0fcf3efa391
KernelDensity.default_bandwidth((df_gradient.A, df_gradient.B))

# ╔═╡ 39bb4bb9-8445-4ee1-a3af-639d8fa96f65
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	
	contour!(kde.x, kde.y, asinh.(kde.density ./ 1e-5), levels=100,)

	scatter!(0, 0, color=:black)


	fig
end

# ╔═╡ 3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 60, normalization=:pdf)
	
	RVUtils.plot_samples!(DataFrame(samples_gradient), LinRange(40, 110, 100), thin=15)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ Cell order:
# ╠═6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═34e43f4a-bcff-41cb-92c4-0c8d600fd053
# ╠═2a422e88-fc0d-4a89-a841-42f3c5c8dace
# ╠═93838644-cad6-4df3-b554-208b7afeb3b8
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═02b266e0-20ee-4992-a83a-9cec72f0d975
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═bd6dfd17-02ee-4855-be37-fecfdab6776f
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╠═5ec475a1-14bb-40f6-856a-69fa9efe087a
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═a162219f-df5a-41d8-bf54-927a355f6431
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╟─55ce0f69-8a96-4bbb-a59f-ee6503624ea6
# ╠═f9d4eade-c648-4f20-8403-07be993fb8c1
# ╠═3377f632-713d-4fec-84a9-b0211b02cb43
# ╠═31a6c2e4-538c-4adc-bbda-5043680b17f7
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╟─6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╠═c897719b-4af9-4c3e-824d-e8ebcd0080d7
# ╠═f8d81628-eea6-457e-b000-92a6ed952f88
# ╠═b7574fe6-c365-4c2b-879e-69c9932abb7a
# ╠═fb882aa1-2722-471a-95e4-1ef486801461
# ╠═b93defe6-6034-4ab4-b83c-10eb289ca280
# ╠═433d1d51-4b2f-4bc0-807a-813eb215231a
# ╠═b32a6521-dae5-4d2a-89ab-6e8b3bc43be3
# ╠═93ee52ad-2c6d-4eac-85cd-ad2a88a70428
# ╠═e88403cc-e502-4e83-a947-43a9bae89af5
# ╠═af21597b-b54a-4445-aa09-22dbad18b42a
# ╠═3f8c8a37-7381-40c7-8826-482aa52d0f6c
# ╠═15e3b145-604e-47f8-a0c3-9466c34277df
# ╠═a8096350-34cf-4cdd-b7fd-f9d1ac3d2d2f
# ╠═3adfa6eb-02f4-4c62-a10f-fedf108a0cb6
# ╠═54ae48c6-c2ac-409b-aaad-b31f73707b3d
# ╠═464af51d-9287-474a-8bcd-f542e7fbe9d2
# ╠═c8949f29-0957-4ca7-bad4-b52eff44fb72
# ╠═d072218f-9da8-4770-9b00-b22c638869bc
# ╠═b983cbd3-a90f-4ed1-b0a7-fd21e903f81d
# ╠═527ccde2-64db-40f5-9b77-8898a92b3b6a
# ╠═820b82a1-63a4-4a83-8d6e-5659ef75ce2d
# ╠═924369b0-aff6-46e4-a0bf-c6766ff93cbf
# ╠═4f00f714-6b16-41c7-a3c8-f2743d6023bf
# ╠═137bbc9c-1e02-4f78-acd1-0e30cf982615
# ╠═3b7cb3ea-3684-4df9-b9e6-919fa11e9ffa
# ╠═f8cc1ebb-617a-46a0-ab22-2f84d75da2ea
# ╠═9158e64f-5a8d-48bf-8043-8e93d41bd626
# ╠═4a83224a-1f2a-4f7c-b670-6a694fc8826a
# ╠═d051635b-123c-45e2-94da-6511f179a2bb
# ╟─062994bc-fb0c-4f08-b7f7-7cc6714bad1e
# ╠═25d9799f-5e89-4ec3-995c-a7309efb798f
# ╠═cb440be7-1f0e-4a20-b15a-82b1820d1ced
# ╠═e16274c0-3b5a-4dc5-9330-f4f1fa06fa87
# ╠═184b4a5d-cbab-44b2-9620-bf928ad81d0e
# ╠═e79d9d82-84cd-4bc5-9ea8-be6b07cacf6d
# ╠═f2aa15e2-824b-4f0e-8ad0-42611abe06da
# ╠═ad607c95-0aea-4d84-8f18-e5b929b9bfca
# ╠═4d7e7b96-74a0-4066-88cf-739c043c7f47
# ╠═2b88915c-7222-44be-a483-b967ea131b80
# ╠═0ca7dc1b-3b41-4089-9c89-20c6e48213ea
# ╠═c48bb30d-1186-4940-b061-91f53e8335e1
# ╟─8555e608-53c1-40d3-b21e-413af8953c30
# ╠═af0d2050-b42e-4a7f-aabb-5b08d23381e9
# ╠═7a9d2e0f-4dcb-4b69-9f68-d43d6dde8bf2
# ╠═3195286d-d85f-43a3-aa25-dae2134f570b
# ╠═4a473039-79f0-4d77-aa0c-681e2fba4f4c
# ╠═88f2918e-e126-420a-96a2-5746a8010f73
# ╠═ed35eb68-74f7-4009-9b68-dfca2ea547af
# ╠═f523cf48-82bf-4b20-9d7c-215bbe10a193
# ╠═51c7c22e-405b-4a47-beb1-b0fcf3efa391
# ╠═39bb4bb9-8445-4ee1-a3af-639d8fa96f65
# ╠═3a9fee80-3ba2-4dc7-9c2a-c57cc11678e9
