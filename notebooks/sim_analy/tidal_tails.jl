### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ fb8bb8ba-34ad-11ef-23e6-1d890b60e0b9
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ 784169c7-0ef2-423a-be8a-5ab62fc5cbb6
using PyFITS

# ╔═╡ b918a965-9e54-42a4-9126-4dd614024ff5
using StatsBase: median, weights, mad, std, sample

# ╔═╡ a4fa1e76-8c2d-4402-b612-2f454bd06b8b
models_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/bootes3")

# ╔═╡ d0d1ecad-4a8d-4c1a-af2b-49f0d3d16bf2
model_dir = "$models_dir/1e6_v35_r3.0/orbit_5Gyr_largeperi"

# ╔═╡ cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
starsfile = "$model_dir/stars/exp2d_rs0.20/final.fits"

# ╔═╡ a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
obs_today_filename = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml")

# ╔═╡ 9c7035e7-c1e7-40d5-8ab6-38f0bb682111
md"""
# Tidal Tails
A detailed analysis of the stars in sculptor
"""

# ╔═╡ c4008b83-b61b-4baa-9fd5-9cced2dc6be8
import TOML

# ╔═╡ 0e9e2085-bf96-4f9d-a758-8af36edf02da
import DensityEstimators as DE

# ╔═╡ e539cb98-fcf2-456c-8ef9-a7061fce46ae
CairoMakie.activate!(type=:png)

# ╔═╡ 217527cb-7f25-4fd9-a4a8-78cb2c744c2b
FIGDIR = joinpath(dirname(starsfile), "figures")

# ╔═╡ 3a913d3f-db30-4ace-b9d0-b16889f7aa2c
using LilGuys; FIGDIR

# ╔═╡ 3ac6fd20-be9b-4d39-a055-5a51e9d4c876
md"""
## Data Loading
"""

# ╔═╡ 0ba05fdb-f859-4381-b2d0-145aa04f7bbf
obs_today_file = TOML.parsefile(obs_today_filename)

# ╔═╡ 4ef955fb-a813-46ad-8f71-7c8a3d371eee
obs_today = lguys.ICRS(;
	ra=obs_today_file["ra"], dec=obs_today_file["dec"],
	distance=obs_today_file["distance"],
	pmra=obs_today_file["pmra"],
	pmdec=obs_today_file["pmdec"],
	radial_velocity=obs_today_file["radial_velocity"],
)

# ╔═╡ 910267ee-aa39-4f04-961b-70f0862d27e2
obs_today_gsr = lguys.transform(lguys.GSR, obs_today)

# ╔═╡ 0d2dafb6-9b64-42c0-ba71-c662ff699aa2
obs_today_df = Dict(
	:ra => obs_today_file["ra"], 
	:dec => obs_today_file["dec"],
	:distance => obs_today_file["distance"],
	:pmra => obs_today_file["pmra"],
	:pmdec => obs_today_file["pmdec"],
	:radial_velocity => obs_today_file["radial_velocity"],
	:pmra_gsr => obs_today_gsr.pmra,
	:pmdec_gsr => obs_today_gsr.pmdec,
	:radial_velocity_gsr => obs_today_gsr.radial_velocity,

)

# ╔═╡ e37559b2-229c-4a37-b516-c6cb7c022b71
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 6d0cff99-3efb-4405-a11d-f13200fa5334
idx_f = orbit_props["idx_f"]

# ╔═╡ 69f4a705-b298-4a15-9a07-0fea39473841
orbit = Orbit(joinpath(model_dir, "centres.hdf5"))

# ╔═╡ 89a34cdf-c724-4c7e-ae42-11da6077feb2
function add_xi_eta_p!(df)
	xi_p, eta_p = lguys.to_orbit_coords(df.ra, df.dec, orbit_props["ra_f"], orbit_props["dec_f"], orbit_props["theta0"])

	df[!, :xi_p] = xi_p
	df[!, :eta_p] = eta_p

	return df
end

# ╔═╡ 7a92c896-7552-4f35-9761-5709d23e9adf
stars = read_fits(starsfile) |> add_xi_eta_p!

# ╔═╡ 6c76cfae-928b-47b3-babe-b0b9a5d68e65
obs_cen = stars[1, :]

# ╔═╡ 075ae901-bbbd-4d10-9d91-c393fc86a8e7
sky_orbit = read_fits(joinpath(model_dir, "skyorbit.fits")) |> add_xi_eta_p!

# ╔═╡ adca7da4-c7e5-4f24-9f50-0413c91fad1d
σ_v = std(stars.radial_velocity, weights(stars.weights))

# ╔═╡ ccc714c4-b2dd-4c83-939b-f0ca3d04c255
r_b_kpc = lguys.break_radius(orbit_props["t_last_peri"]/lguys.T2GYR, σ_v/lguys.V2KMS)

# ╔═╡ c351926c-1e3d-44d7-8f39-6fb5893b9b0d
r_b_arcmin = lguys.kpc2arcmin(r_b_kpc, orbit_props["distance_f"])

# ╔═╡ 3884f295-9806-4c03-8d5f-ea9ac680c2b8
md"""
# Plots
"""

# ╔═╡ 4c1d6f6f-e257-4126-9b4a-8e5aa8470295
function ra_dec_axis(ddeg=5; kwargs...)
	fig = Figure(;kwargs...)
	
	dy = ddeg
	dx = dy * 1/cosd(obs_cen.dec)
	limits = (obs_cen.ra .+ (-dx, dx), obs_cen.dec .+ (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=limits,
		aspect = 1,
		xgridvisible=false,
		ygridvisible=false,
		xreversed=true,
		
	)

	return fig, ax
end

# ╔═╡ 816b9db9-26c6-4ac8-9a46-82209d2cdc85
idx_orbit = idx_f - 5: idx_f

# ╔═╡ e9e35643-168e-4e87-a880-6831b46145c7
let 
	fig, ax = ra_dec_axis()

	bins = 100
	limits = ax.limits.val
	x = stars.ra
	y = stars.dec


	hi = Arya.histogram2d(x, y, bins, weights=stars.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-10, maximum(hi.values)))

	lines!(sky_orbit.ra[idx_orbit], sky_orbit.dec[idx_orbit])
	
	Colorbar(fig[1, 2], h,
		label="stellar density",
			 ticks = Makie.automatic
	)

	fig
end

# ╔═╡ 2c9e8ad0-bbff-4086-974e-89269868e324
pmra_f, pmdec_f = sky_orbit[idx_f, :].pmra, sky_orbit[idx_f, :].pmdec

# ╔═╡ a89adc86-67a0-453f-b382-96f721f74d39
diff(sky_orbit.ra)[idx_f-1], diff(sky_orbit.dec)[idx_f-1]

# ╔═╡ 93a266da-bbfa-4c83-9ed0-54a6283383f5
sky_orbit.ra[idx_f], orbit_props["ra_f"]

# ╔═╡ 2fefefa3-ae22-4979-9cab-6d121306c739
sky_orbit.dec[idx_f], orbit_props["dec_f"]

# ╔═╡ 37ef9662-e253-4d99-b26f-4141a58d738b
lguys.mean(stars.ra, stars.weights .^ 3), lguys.mean(stars.dec, stars.weights .^ 3)

# ╔═╡ 9419722f-9c65-45a7-b9c5-c666b22bf70e
sind(orbit_props["theta0"] ), cosd(orbit_props["theta0"])

# ╔═╡ 12c8d1f6-30fb-4616-b3dd-eb0aefd7d450
md"""
# Transform to stream coordinates
"""

# ╔═╡ 237c25a1-5537-400c-826d-5827c52c38a3
function xi_eta_p_axis(dx=10, dy=5; kwargs...)
	fig = Figure(;kwargs...)
	

	limits = ((-dx, dx), (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel=L"$\xi'$ / degrees",
		ylabel=L"$\eta'$ / degrees",
		limits=limits,
		aspect = DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
		xreversed=true,
		
	)

	return fig, ax
end

# ╔═╡ 7688c0ac-be70-470c-9db5-f2ed937e2fb5
function xi_eta_axis(dx=10, dy=dx; kwargs...)
	fig = Figure(;kwargs...)
	

	limits = ((-dx, dx), (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel=L"$\xi$ / degrees",
		ylabel=L"$\eta$ / degrees",
		limits=limits,
		aspect = DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
		xreversed=true,
		
	)

	return fig, ax
end

# ╔═╡ 3868afc5-a406-4ced-ada0-96336e5fbe96
let 
	fig, ax = xi_eta_axis()

	bins = 100
	limits = ax.limits.val
	x = stars.xi
	y = stars.eta


	hi = Arya.histogram2d(x, y, bins, weights=stars.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-5, 1) .* maximum(hi.values))

	lines!(sky_orbit.xi[idx_orbit], sky_orbit.eta[idx_orbit])
	arc!(Point2f(0, 0), r_b_arcmin / 60, -π, π, color=COLORS[3])
	text!(0, -r_b_arcmin/60, text=L"r_b", color=COLORS[3])

		
	pmra, pmdec = 30 .* (pmra_f, pmdec_f) ./ radii([pmra_f, pmdec_f])

	# annotation!([0], [0], [pmra], [pmdec], )

	Colorbar(fig[1, 2], h,
		label="stellar density",
		ticks = Makie.automatic
	)

	@savefig "tidal_stream_xi_eta"
	fig
end

# ╔═╡ d1564bf3-a3c4-410b-b61b-69ac56b618e5
let 
	fig, ax = xi_eta_p_axis(20, 10)
	
	bins = 100
	limits = ax.limits.val
	x = stars.xi_p
	y = stars.eta_p


	hi = Arya.histogram2d(x, y, bins, weights=stars.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-5, 1) .* maximum(hi.values))

	
	Colorbar(fig[1, 2], h,
		label="stellar density", ticks=Makie.automatic, 
			 tellheight=false
	)

	resize_to_layout!()
	@savefig "tidal_stream_xi_eta_prime"

	fig
end

# ╔═╡ 49705acd-d623-4860-a5e7-759f112553ab
md"""
# Seperation into stream components
"""

# ╔═╡ 9241f6a2-d829-4370-a497-35f6fb3358dd
r_centre = 60 # arcminutes

# ╔═╡ b09e9622-b522-4308-ad95-c939ed0a5315
r_max = 600 # arcminutes

# ╔═╡ 95b7627e-c17a-4be6-ac66-c3212a8feb23
not = !

# ╔═╡ 82b9c535-4991-486b-9fb7-d98159bfda8f
filt_cen = stars.r_ell .< r_centre

# ╔═╡ 2cb76b38-b768-4184-b73c-e8fb350351d8
filt_leading = not.(filt_cen) .& (stars.xi_p .> 0)

# ╔═╡ c0166054-7c1f-474f-a4bc-3b122034e923
filt_trailing = not.(filt_cen) .& (stars.xi_p .< 0)

# ╔═╡ aa157085-0705-44ed-9d65-5626658d71e7
filt_dist = stars.r_ell .< r_max

# ╔═╡ dd6faf9c-8206-46f8-bd6f-4c2e7a4b3887
filt_excl = not.(filt_dist)

# ╔═╡ d240246d-f79a-4f95-a5b9-b355e7bc092f
sum(filt_leading .& filt_dist), sum(stars.weights[filt_leading .& filt_dist])

# ╔═╡ cb8c2e79-a5c2-4a59-af40-393b83f64bc7
sum(filt_trailing .& filt_dist), sum(stars.weights[filt_trailing .& filt_dist])

# ╔═╡ 5ea9a308-13b4-4d64-84b2-2143be6e9a23
sum(stars.weights[stars.r_ell .> r_b_arcmin])

# ╔═╡ 51d391f8-130a-4c85-a165-f1cebaea3812
sum(stars.weights[filt_dist])

# ╔═╡ 268bf2e9-d6a0-4fca-b77b-005fc79e9497
r_b_arcmin / 60

# ╔═╡ 9c1839d6-1076-4598-8da6-49c02ec10580
sum(stars.weights[filt_cen])

# ╔═╡ e58263aa-48d5-4a30-8c2d-f6ed19ada10f
1 - sum(stars.weights[filt_cen])

# ╔═╡ c57edba1-09f9-4fe9-bc73-5746884ee7c2
md"""
Below, we plot our defined filters.
"""

# ╔═╡ 119040e2-ef45-4fc2-809d-aec333e276d0
let 
	fig, ax = ra_dec_axis(1.3r_max / 60)

	bins = 100
	limits = ax.limits.val
	x = stars.ra
	y = stars.dec

	filt = filt_cen .& filt_dist
	scatter!(x[filt], y[filt], label="centre", alpha=0.03, markersize=1)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=0.01, markersize=2)

	filt = filt_leading .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=0.01, markersize=2)

	filt = filt_excl
	scatter!(x[filt], y[filt], label="excluded", alpha=0.01, markersize=1)
	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ 443ac755-70fc-4193-945d-0e98622c5919
let 
	fig, ax = ra_dec_axis(r_centre / 60)
	ax.title = "central stars cut"

	x = stars.ra
	y = stars.dec
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=0.05, markersize=3)
	fig
end

# ╔═╡ 1d35a894-1eca-4ee0-9d1a-6ac4704b912c
md"""
# Plots of each stream arm
"""

# ╔═╡ 6b0a0ea3-f6f7-4cee-be11-bce6872ab870
let
	dr = 0.2
		
	fig = Figure()
	limits = (obs_cen.pmra - dr, obs_cen.pmra + dr, obs_cen.pmdec - dr, obs_cen.pmdec + dr)
	
	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="stars within $r_max arcmin",
		limits=limits,
	aspect=DataAspect())

	filt = stars.r_ell .< r_max

	x = stars.pmra[filt]
	y = stars.pmdec[filt]
	w = weights=stars.weights[filt]
	h = Arya.hist2d!(ax, x, y, bins=100, weights=w, limits=limits,
		colorscale=log10, colorrange=(1e-15, nothing))
	

	Colorbar(fig[1, 2], h, label="stellar density")
	
	fig
end

# ╔═╡ 54205139-3c7b-4beb-b9cb-c87b272df58a
let
	dr = 1
		
	fig = Figure()
	limits = (obs_cen.pmra - dr, obs_cen.pmra + dr, obs_cen.pmdec - dr, obs_cen.pmdec + dr)
	
	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="all stars",
		limits=limits,
	aspect=DataAspect())


	x = stars.pmra
	y = stars.pmdec
	w = weights=stars.weights
	h = Arya.hist2d!(ax, x, y, bins=100, weights=w, limits=limits,
		colorscale=log10, colorrange=(1e-20, nothing))
	

	Colorbar(fig[1, 2], h, label="stellar density")
	
	fig
end

# ╔═╡ 91955a7b-009e-4afa-aca1-312f4a779bb9
let
	dr = 1
		
	fig = Figure()
	limits = (obs_cen.pmra_gsr - dr, obs_cen.pmra_gsr + dr, obs_cen.pmdec_gsr - dr, obs_cen.pmdec_gsr + dr)
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="all stars",
		limits=limits,
	aspect=DataAspect())


	x = stars.pmra_gsr
	y = stars.pmdec_gsr
	w = weights=stars.weights
	h = Arya.hist2d!(ax, x, y, bins=100, weights=w, limits=limits,
		colorscale=log10, colorrange=(1e-20, nothing))
	

	Colorbar(fig[1, 2], h, label="stellar density")
	
	fig
end

# ╔═╡ 19bdc540-63a8-46ad-8d6f-a425d64cdd81
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="stars within $r_max arcmin",
	aspect=DataAspect())

	bins = 100
	limits = ax.limits.val
	x = stars.pmra
	y = stars.pmdec
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=0.03, markersize=3)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=0.01, markersize=5)

	filt = filt_leading .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=0.01, markersize=5)

	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ 51e4628a-a799-48a3-9ad0-b5c2bd7a8d87
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="stars within $r_max arcmin",
	aspect=DataAspect())

	bins = 100
	limits = ax.limits.val
	x = stars.pmra_gsr
	y = stars.pmdec_gsr
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=0.03, markersize=3)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=1, markersize=5)

	filt = filt_leading .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=1, markersize=5)

	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ e098260a-718a-4f43-b12c-cab0a1302b90
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{\alpha*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="",)

	bins = 100
	limits = ax.limits.val
	x = stars.distance
	y = stars.radial_velocity
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=0.03, markersize=3)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=0.1, markersize=5)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=0.1, markersize=5)

	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ be0158cf-c322-4147-8f01-8b6a715bc0dc
let
	fig = Figure()

	dr = 0.1
	limits = (obs_today.pmra - dr, obs_today.pmra + dr, obs_today.pmdec - dr, obs_today.pmdec + dr)
	
	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="Centre",
		aspect=DataAspect()
	)

	filt = filt_cen
	hist2d!(stars.pmra[filt], stars.pmdec[filt], weights=stars.weights[filt], limits=limits, colorscale=log10, colorrange=(1e-15, nothing), bins=100)

	errorscatter!(obs_today.pmra, obs_today.pmdec)
	

	fig
end

# ╔═╡ 5eb22695-57c4-4ebf-b412-588f5e366bf5
let
		fig = Figure()

	dr = 0.1
	limits = (obs_cen.pmra - dr, obs_cen.pmra + dr, obs_cen.pmdec - dr, obs_cen.pmdec + dr)

	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="leading arm",
		aspect=DataAspect()
	)

	filt = filt_leading .& filt_dist
	hist2d!(stars.pmra[filt], stars.pmdec[filt], weights=stars.weights[filt], limits=limits, colorscale=log10, colorrange=(1e-15, nothing), bins=100)

	scatter!(obs_cen.pmra, obs_cen.pmdec)
	

	fig
end

# ╔═╡ 79dd7e3d-8564-4389-960b-05512b0143c0
let
	fig = Figure()
	
	dr = 0.1
	limits = (obs_cen.pmra - dr, obs_cen.pmra + dr, obs_cen.pmdec - dr, obs_cen.pmdec + dr)
	
	ax = Axis(fig[1,1],
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="trailing arm",
		aspect=DataAspect()
	)

	filt = filt_trailing .& filt_dist
	hist2d!(stars.pmra[filt], stars.pmdec[filt], weights=stars.weights[filt], colorscale=log10, colorrange=(1e-15, nothing), bins=100, limits=limits)

	scatter!(obs_cen.pmra, obs_cen.pmdec)
	

	fig
end

# ╔═╡ d4d50209-e9a8-40dc-9f45-4e8000f70b39
md"""
# Functions of orbital coordinates
"""

# ╔═╡ bf83963a-5291-4de2-9cd0-5fc318c6fd53
md"""
### Utility functions
"""

# ╔═╡ 1d4a7103-3d37-4597-a180-88a4c3733497
bins = -0.4r_max/60:0.5:0.4r_max/60

# ╔═╡ 45b7b28f-9125-48bb-8895-a19d60b825d0
eta_max = 1

# ╔═╡ e12ac0ef-7903-4cbc-84f9-6052a701771c
function plot_rb!(ax, r_b; y0=nothing, dy=nothing)
	ylim = ax.limits[][3:4]
	println(ylim)
	if y0 === nothing
		
		y0 =  ylim[1] + 0.1 * (ylim[2] - ylim[1])
		dy =  0.07 * (ylim[2] - ylim[1])
	end
	
	annotation!(0, -40, r_b, y0, text="r_b")
	# text!(r_b, y0, text=L"r_b")
end

# ╔═╡ 404195ec-1a84-4c64-97ec-fe4a583a11a5


# ╔═╡ f1cb5d80-57a8-43e8-8819-4e1bc134edea
function weighted_mad(samples::Vector{T}, weights::Vector{T}) where T
    # Ensure that the number of samples matches the number of weights
    if length(samples) != length(weights)
        throw(ArgumentError("Length of samples and weights must be equal."))
    end

    # Compute the weighted median
    weighted_med = lguys.quantile(samples, weights, 0.5)

    # Compute the absolute deviations from the weighted median
    abs_deviations = abs.(samples .- weighted_med)

    # Compute the weighted median of the absolute deviations
    mad = lguys.quantile(abs_deviations, weights, 0.5)

    return mad
end

# ╔═╡ 4b0515a9-e55d-4bed-b6d3-0a2f12d0a238
function binned_median(x, y, bins; w)
	idxs = DE.bin_indices(x, bins)

	N =  length(bins) - 1
	y_m = Vector{Float64}(undef, N)
	y_l = Vector{Float64}(undef, N)
	y_h = Vector{Float64}(undef, N)
	
	for i in eachindex(y_m)
		filt = idxs .== i
		y_m[i] = median(y[filt], weights(w[filt]))

		s = weighted_mad(y[filt], (w[filt]))
		y_l[i], y_h[i] = lguys.quantile(y[filt], (w[filt]), [0.16, 0.84])

		# y_l[i] = -y_l[i] + y_m[i]
		# y_h[i] =  y_h[i] - y_m[i]
		
	end

	return y_m, y_l, y_h
end

# ╔═╡ 57528a9e-fbcf-4926-a3a1-574935df69b4
function binned_mean(x, y, bins; w)
	idxs = DE.bin_indices(x, bins)

	N =  length(bins) - 1
	y_m = Vector{Float64}(undef, N)
	y_l = Vector{Float64}(undef, N)
	y_h = Vector{Float64}(undef, N)
	
	for i in eachindex(y_m)
		filt = idxs .== i
		y_m[i] = lguys.mean(y[filt], (w[filt]))

		s =  lguys.std(y[filt], (w[filt]))

		y_l[i] = y_m[i] - s
		y_h[i] = y_m[i] + s
		
	end

	return y_m, y_l, y_h
end

# ╔═╡ 20e6efaa-95ab-4163-a584-f1fea47f02c7
function hist_eff(x, bins; w)
	idxs = DE.bin_indices(x, bins)

	N =  length(bins) - 1
	y = Vector{Float64}(undef, N)
	
	for i in eachindex(y)
		filt = idxs .== i

		z = w[filt]
		y[i] = sum(z) .^ 2 / sum(z .^ 2)
	end

	return y
end

# ╔═╡ 51555d9f-65fc-4c5c-86ca-a1a0153772ed


# ╔═╡ 35f1ac55-eff5-4847-ad40-b6e040ec168d
md"""
### Radial Velocity
"""

# ╔═╡ 2c636a11-0908-4e5e-b31b-af71833eb6f5
obs_labels = Dict(
	:radial_velocity => L"{v}_\textrm{los} / \textrm{km\,s^{-1}}",
	:radial_velocity_gsr => L"\tilde{v}_\textrm{los} / \textrm{km\,s^{-1}}",
	:distance => "distance / kpc",
	:pmra => L"{\mu}_{\alpha*}/ \textrm{mas\,yr^{-1}}",
	:pmra_gsr => L"\tilde{\mu}_{\alpha*}/ \textrm{mas\,yr^{-1}}",
	:pmdec => L"{\mu}_{\delta}/ \textrm{mas\,yr^{-1}}",
	:pmdec_gsr => L"\tilde{\mu}_{\delta}/ \textrm{mas\,yr^{-1}}",
	:xi_p => L"\xi' / \textrm{degrees}",
	:eta_p => L"\eta' / \textrm{degrees}",
)

# ╔═╡ f171df4d-2a52-4132-9a80-66eec6663179
obs_dy = Dict(
	:radial_velocity => 40,
	:radial_velocity_gsr => 40,
	:distance => 20,
	:pmra => 0.2,
	:pmra_gsr =>  0.2,
	:pmdec =>  0.2,
	:pmdec_gsr =>  0.2,
)

# ╔═╡ 499f3cd7-b7fc-4fe2-acc0-5e2d2697c3d6
function plot_hist2d_coord(sym; xsym = :xi_p)
	mass = stars.weights
	x = stars[:, xsym]
	
	r_max = 5
	
	y0 = obs_today_df[sym]
	dy0 = obs_dy[sym]
	y = stars[:, sym]

	limits = (-r_max, r_max,  y0 - dy0, y0 + dy0)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=obs_labels[xsym],
		ylabel=obs_labels[sym]
	)
	
	h = Arya.hist2d!(ax, x, y, weights=mass, bins=100, limits=limits,
		colorrange=(1e-10, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density", ticks=Makie.automatic)

	return fig
end

# ╔═╡ 52d96bbd-6ce5-42eb-8cef-9270666793f2
import Random

# ╔═╡ 5bb11976-1528-441a-9d10-9468e204dfb2


# ╔═╡ 594e0f49-6d4a-467f-a357-b58da4e7d261
function scatter_coord(sym; xsym = :xi_p)
	mass = stars.weights
	x = stars[:, xsym]
	
	r_max = 5
	
	y0 = obs_today_df[sym]
	dy0 = obs_dy[sym]
	y = stars[:, sym]

	limits = (-r_max, r_max,  y0 - dy0, y0 + dy0)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=obs_labels[xsym],
		ylabel=obs_labels[sym]
	)

	idx = sample(eachindex(x), weights(mass), 10_000)
	scatter!(x[idx], y[idx], )

	return fig
end

# ╔═╡ d81cf267-4020-4398-9cdd-0da783f122c5
scatter_coord(:pmra)

# ╔═╡ 88186b14-13ff-48a4-ba6a-9cd988a7e156
scatter_coord(:pmdec)

# ╔═╡ 445fa03b-c540-44c7-b596-366162523f8c
scatter_coord(:radial_velocity)

# ╔═╡ cc2be8c5-bdfc-4c55-815c-1ccc1d261678
for sym in [:pmra, :pmra_gsr, :pmdec, :pmdec_gsr, :radial_velocity, :distance, ]
	@info plot_hist2d_coord(sym)
end

# ╔═╡ a52f14bf-eb6d-4523-b6ee-38af9b822a49
for sym in [:pmra, :pmra_gsr, :pmdec, :radial_velocity, :distance, ]
	@info plot_hist2d_coord(sym, xsym=:eta_p)
end

# ╔═╡ c86fddf0-d293-4735-8397-5a1889f2d1be
function plot_obs_bin_means!(gs, ysym; xsym = :xi_p, eta_max=1, bins=bins)
	if xsym == :xi_p
		filt = abs.(stars.eta_p) .< eta_max
	elseif xsym == :eta_p
		filt = abs.(stars.xi_p) .< eta_max
	end

	w = stars.weights[filt]
	y = stars[filt, ysym]
	x = stars[filt, xsym]

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	
	yl0 = minimum(y_l)
	yh0 = maximum(y_h)
	dy0 = yh0 - yl0
	limits=(minimum(bins), maximum(bins), yl0 - 0.05dy0, yh0 + 0.05dy0)

	ax = Axis(gs,
		xlabel=obs_labels[xsym],
		ylabel=obs_labels[ysym],
		limits=limits
	)


	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)
	
	plot_rb!(ax, r_b_arcmin / 60)
	plot_rb!(ax, -r_b_arcmin / 60)

	ax
end

# ╔═╡ 6ac742ef-ff1f-4a1c-a812-e2962bd6a50d
function plot_obs_bin_means(ysym; xsym = :xi_p, eta_max=1, r_b_arcmin=r_b_arcmin, bins=bins, stars=stars)
	if xsym == :xi_p
		filt = abs.(stars.eta_p) .< eta_max
	elseif xsym == :eta_p
		filt = abs.(stars.xi_p) .< eta_max
	end

	w = stars.weights[filt]
	y = stars[filt, ysym]
	x = stars[filt, xsym]

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	
	fig = Figure()

	yl0 = minimum(y_l)
	yh0 = maximum(y_h)
	dy0 = yh0 - yl0
	limits=(minimum(bins), maximum(bins), yl0 - 0.05dy0, yh0 + 0.05dy0)

	ax = Axis(fig[1,1],
		xlabel=obs_labels[xsym],
		ylabel=obs_labels[ysym],
		limits=limits
	)


	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)
	
	plot_rb!(ax, r_b_arcmin / 60)

	fig
end

# ╔═╡ 32ce8cf5-94fa-4bb4-83eb-d3c6329d497a
1 ÷ 2

# ╔═╡ 5eaaa70a-3210-4ab2-9318-b56ad4962c22
let
	fig = Figure(size=(600,500))

	for i in 1:4
		ysym = [:pmra, :pmdec, :distance, :radial_velocity][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		plot_obs_bin_means!(gs, ysym)	|> lguys.hide_grid!


	end

	fig
end

# ╔═╡ 6eb0507a-c7dd-4009-b699-ae6423623f63
bins_cen = -1.5:0.1:1.5

# ╔═╡ 041ba6c2-0a5c-4a59-b536-95db6e68716f
let
	fig = Figure(size=(600,500))

	for i in 1:4
		ysym = [:pmra_gsr, :pmdec_gsr, :distance, :radial_velocity_gsr][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		plot_obs_bin_means!(gs, ysym)	|> lguys.hide_grid!

	end
	@savefig "velocities_along_orbit"
	fig
end

# ╔═╡ c6224653-9dcd-4e17-9999-ccb69da6c222
let
	fig = Figure(size=(600,500))

	for i in 1:4
		ysym = [:pmra_gsr, :pmdec_gsr, :distance, :radial_velocity_gsr][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		plot_obs_bin_means!(gs, ysym, xsym=:eta_p)	|> lguys.hide_grid!

	end

	fig
end

# ╔═╡ 13c53be9-4922-4872-8ebf-8ed147a72fff


# ╔═╡ 4be2b3c0-98df-4d25-8200-a051a5dba4e9
let
	fig = Figure(size=(600,500))

	for i in 1:4
		ysym = [:pmra, :pmdec, :distance, :radial_velocity][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		plot_obs_bin_means!(gs, ysym; xsym=:eta_p)	|> lguys.hide_grid!

	end

	fig
end

# ╔═╡ e3ecf2c5-67e9-4571-8743-04c67755a23b
for ysym in [:pmra, :pmra_gsr, :pmdec, :pmdec_gsr,  :distance, :radial_velocity, :radial_velocity_gsr]
	@info plot_obs_bin_means(ysym)
end

# ╔═╡ 4b4938c2-bf21-455e-88d5-7d3843a6e005
let
	xsym = :r_ell
	ysym = :radial_velocity_gsr
	bins = 60 * (0:0.1:1.5)

	filt = stars.r_ell .< 1.5 * 60
	w = stars.weights[filt]
	y = stars[filt, ysym]
	x = stars[filt, xsym]

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_mean(x, y, bins, w=w)

	sy = y_m - y_l
	limits=(minimum(bins), maximum(bins), 0.9minimum(sy), 1.1maximum(sy))

	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log r / arcmin",
		ylabel="velocity dispersion",
		limits=limits
	)


	scatter!(x_m, sy)
	
	fig
end

# ╔═╡ cae30107-1aee-48c0-8104-720cf5535b02
for ysym in [:pmra, :pmra_gsr, :pmdec, :distance, :radial_velocity]
	@info plot_obs_bin_means(ysym, xsym=:eta_p)
end

# ╔═╡ 925be0a1-0afa-4b5a-b9c9-003833e80a28
let	
	w = stars.weights
	x = stars.xi_p

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel="effective num of observations",
		yscale=log10
	)
	
	scatter!(midpoints(bins), hist_eff(x, bins, w=w))

	fig
end

# ╔═╡ ae9d8d21-1269-4053-89a3-3a8287a4ca70
md"""
# Validation
"""

# ╔═╡ 4c8b9b51-5e72-42a7-a919-ca2c2d1c5294
md"""
Below is a scatter plot of all stars in polar coordinates wrt ``\xi'`` and distance. The highlighted region (orange) is all stars within the specified distance of the centre for further analysis. This is just to make sure nothing crazy is left in the filter
"""

# ╔═╡ 7178c98a-15d1-4f10-9850-2079a6fdaa27
# ╠═╡ disabled = true
#=╠═╡
let
	fig = Figure()
	ax = PolarAxis(fig[1, 1])
	
	scatter!(deg2rad.(stars.xi_p), stars.distance, alpha=0.1)
	scatter!(deg2rad.(stars.xi_p[filt_dist]), stars.distance[filt_dist], alpha=0.1)

	fig
end
  ╠═╡ =#

# ╔═╡ 90e0d677-a837-4a17-8c64-68ed350dd3c8
# ╠═╡ disabled = true
#=╠═╡
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel=L"$\eta'$ / degrees",
		ylabel="distance / kpc"
	)
	
	scatter!(stars.eta_p, stars.distance, alpha=0.1)
	scatter!(stars.eta_p[filt_dist], stars.distance[filt_dist], alpha=0.1)

	fig
end
  ╠═╡ =#

# ╔═╡ 26d120cb-0bd4-432d-9ce1-fc9a23a5067d
md"""
The next two plots validate the rotated coordinate system.
"""

# ╔═╡ 65ca1318-8f62-4133-ab14-f3d3e1190029
# ╠═╡ disabled = true
#=╠═╡
let 
	fig, ax = ra_dec_axis()

	bins = 100
	limits = ax.limits.val
	x = stars.ra
	y = stars.dec
	color = stars.eta_p
	h = scatter!(x, y, color=color, colorrange=(-5, 5), colormap=:redsblues)
	
	Colorbar(fig[1, 2], h,
		label=L"\eta'"
	)
	fig
end
  ╠═╡ =#

# ╔═╡ 1b31eabc-3b6f-4d82-9cbb-62f1c53408de
# ╠═╡ disabled = true
#=╠═╡
let 
	fig, ax = ra_dec_axis()

	bins = 100
	limits = ax.limits.val
	x = stars.ra[filt_dist]
	y = stars.dec[filt_dist]
	color = stars.xi_p[filt_dist]
	h = scatter!(x, y, color=color, #colorrange=(-5, 5),
	colormap=:redsblues)
	
	Colorbar(fig[1, 2], h,
		label=L"\xi'"
	)
	fig
end
  ╠═╡ =#

# ╔═╡ 7ec71b3a-4c34-4803-923f-2effeff7cfcf
md"""
The orbit in the rotated coordinate system. We expect the orbit to be ~ a horizontal line and increasing (lighter color) going to positive ``\xi'``
"""

# ╔═╡ 84e4a23c-c926-4c86-ab11-07a08ef7d1b3
# ╠═╡ disabled = true
#=╠═╡
let
	fig, ax = xi_eta_axis()

	scatter!(stars.xi_p, stars.eta_p)
	p = lines!(sky_orbit.xi_p[idx_orbit], sky_orbit.eta_p[idx_orbit], 
		color=idx_orbit, linewidth=4)

	Colorbar(fig[1, 2], p, label="time / Gyr")
	fig
end
  ╠═╡ =#

# ╔═╡ e013dcb6-ea17-4fd5-af4d-318188600107
# ╠═╡ disabled = true
#=╠═╡
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.eta_p[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\eta'",
		aspect=DataAspect(),
		limits=(-5, 5, nothing, nothing)
	)

	bin_idxs = DE.bin_indices(x, bins)

	N = length(bins) - 1

	for i in 1:N
		filt = bin_idxs .== i
		scatter!(x[filt], y[filt], 
			color=i, colorrange=(1, N), markersize=3)
	end

	lines!(sky_orbit.xi_p[idx_orbit], sky_orbit.eta_p[idx_orbit], color=COLORS[3])

	vlines!(bins)
	fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═a4fa1e76-8c2d-4402-b612-2f454bd06b8b
# ╠═d0d1ecad-4a8d-4c1a-af2b-49f0d3d16bf2
# ╠═cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
# ╠═a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
# ╟─9c7035e7-c1e7-40d5-8ab6-38f0bb682111
# ╠═fb8bb8ba-34ad-11ef-23e6-1d890b60e0b9
# ╠═c4008b83-b61b-4baa-9fd5-9cced2dc6be8
# ╠═3a913d3f-db30-4ace-b9d0-b16889f7aa2c
# ╠═784169c7-0ef2-423a-be8a-5ab62fc5cbb6
# ╠═0e9e2085-bf96-4f9d-a758-8af36edf02da
# ╠═b918a965-9e54-42a4-9126-4dd614024ff5
# ╠═e539cb98-fcf2-456c-8ef9-a7061fce46ae
# ╠═217527cb-7f25-4fd9-a4a8-78cb2c744c2b
# ╠═3ac6fd20-be9b-4d39-a055-5a51e9d4c876
# ╠═0ba05fdb-f859-4381-b2d0-145aa04f7bbf
# ╠═910267ee-aa39-4f04-961b-70f0862d27e2
# ╠═0d2dafb6-9b64-42c0-ba71-c662ff699aa2
# ╠═4ef955fb-a813-46ad-8f71-7c8a3d371eee
# ╠═6d0cff99-3efb-4405-a11d-f13200fa5334
# ╠═e37559b2-229c-4a37-b516-c6cb7c022b71
# ╠═7a92c896-7552-4f35-9761-5709d23e9adf
# ╠═6c76cfae-928b-47b3-babe-b0b9a5d68e65
# ╠═69f4a705-b298-4a15-9a07-0fea39473841
# ╠═89a34cdf-c724-4c7e-ae42-11da6077feb2
# ╠═075ae901-bbbd-4d10-9d91-c393fc86a8e7
# ╠═adca7da4-c7e5-4f24-9f50-0413c91fad1d
# ╠═ccc714c4-b2dd-4c83-939b-f0ca3d04c255
# ╠═c351926c-1e3d-44d7-8f39-6fb5893b9b0d
# ╠═3884f295-9806-4c03-8d5f-ea9ac680c2b8
# ╠═4c1d6f6f-e257-4126-9b4a-8e5aa8470295
# ╠═e9e35643-168e-4e87-a880-6831b46145c7
# ╠═3868afc5-a406-4ced-ada0-96336e5fbe96
# ╠═816b9db9-26c6-4ac8-9a46-82209d2cdc85
# ╠═2c9e8ad0-bbff-4086-974e-89269868e324
# ╠═a89adc86-67a0-453f-b382-96f721f74d39
# ╠═93a266da-bbfa-4c83-9ed0-54a6283383f5
# ╠═2fefefa3-ae22-4979-9cab-6d121306c739
# ╠═37ef9662-e253-4d99-b26f-4141a58d738b
# ╠═9419722f-9c65-45a7-b9c5-c666b22bf70e
# ╟─12c8d1f6-30fb-4616-b3dd-eb0aefd7d450
# ╠═237c25a1-5537-400c-826d-5827c52c38a3
# ╠═7688c0ac-be70-470c-9db5-f2ed937e2fb5
# ╠═d1564bf3-a3c4-410b-b61b-69ac56b618e5
# ╟─49705acd-d623-4860-a5e7-759f112553ab
# ╠═9241f6a2-d829-4370-a497-35f6fb3358dd
# ╠═b09e9622-b522-4308-ad95-c939ed0a5315
# ╠═95b7627e-c17a-4be6-ac66-c3212a8feb23
# ╠═82b9c535-4991-486b-9fb7-d98159bfda8f
# ╠═2cb76b38-b768-4184-b73c-e8fb350351d8
# ╠═c0166054-7c1f-474f-a4bc-3b122034e923
# ╠═aa157085-0705-44ed-9d65-5626658d71e7
# ╠═dd6faf9c-8206-46f8-bd6f-4c2e7a4b3887
# ╠═d240246d-f79a-4f95-a5b9-b355e7bc092f
# ╠═cb8c2e79-a5c2-4a59-af40-393b83f64bc7
# ╠═5ea9a308-13b4-4d64-84b2-2143be6e9a23
# ╠═51d391f8-130a-4c85-a165-f1cebaea3812
# ╠═268bf2e9-d6a0-4fca-b77b-005fc79e9497
# ╠═9c1839d6-1076-4598-8da6-49c02ec10580
# ╠═e58263aa-48d5-4a30-8c2d-f6ed19ada10f
# ╟─c57edba1-09f9-4fe9-bc73-5746884ee7c2
# ╠═119040e2-ef45-4fc2-809d-aec333e276d0
# ╟─443ac755-70fc-4193-945d-0e98622c5919
# ╟─1d35a894-1eca-4ee0-9d1a-6ac4704b912c
# ╟─6b0a0ea3-f6f7-4cee-be11-bce6872ab870
# ╟─54205139-3c7b-4beb-b9cb-c87b272df58a
# ╟─91955a7b-009e-4afa-aca1-312f4a779bb9
# ╠═19bdc540-63a8-46ad-8d6f-a425d64cdd81
# ╟─51e4628a-a799-48a3-9ad0-b5c2bd7a8d87
# ╟─e098260a-718a-4f43-b12c-cab0a1302b90
# ╟─be0158cf-c322-4147-8f01-8b6a715bc0dc
# ╟─5eb22695-57c4-4ebf-b412-588f5e366bf5
# ╟─79dd7e3d-8564-4389-960b-05512b0143c0
# ╟─d4d50209-e9a8-40dc-9f45-4e8000f70b39
# ╟─bf83963a-5291-4de2-9cd0-5fc318c6fd53
# ╠═1d4a7103-3d37-4597-a180-88a4c3733497
# ╠═45b7b28f-9125-48bb-8895-a19d60b825d0
# ╠═e12ac0ef-7903-4cbc-84f9-6052a701771c
# ╠═404195ec-1a84-4c64-97ec-fe4a583a11a5
# ╠═f1cb5d80-57a8-43e8-8819-4e1bc134edea
# ╠═4b0515a9-e55d-4bed-b6d3-0a2f12d0a238
# ╠═57528a9e-fbcf-4926-a3a1-574935df69b4
# ╠═20e6efaa-95ab-4163-a584-f1fea47f02c7
# ╠═51555d9f-65fc-4c5c-86ca-a1a0153772ed
# ╟─35f1ac55-eff5-4847-ad40-b6e040ec168d
# ╠═2c636a11-0908-4e5e-b31b-af71833eb6f5
# ╠═f171df4d-2a52-4132-9a80-66eec6663179
# ╠═499f3cd7-b7fc-4fe2-acc0-5e2d2697c3d6
# ╠═52d96bbd-6ce5-42eb-8cef-9270666793f2
# ╠═5bb11976-1528-441a-9d10-9468e204dfb2
# ╠═594e0f49-6d4a-467f-a357-b58da4e7d261
# ╠═d81cf267-4020-4398-9cdd-0da783f122c5
# ╠═88186b14-13ff-48a4-ba6a-9cd988a7e156
# ╠═445fa03b-c540-44c7-b596-366162523f8c
# ╠═cc2be8c5-bdfc-4c55-815c-1ccc1d261678
# ╠═a52f14bf-eb6d-4523-b6ee-38af9b822a49
# ╠═c86fddf0-d293-4735-8397-5a1889f2d1be
# ╠═6ac742ef-ff1f-4a1c-a812-e2962bd6a50d
# ╠═32ce8cf5-94fa-4bb4-83eb-d3c6329d497a
# ╠═5eaaa70a-3210-4ab2-9318-b56ad4962c22
# ╠═6eb0507a-c7dd-4009-b699-ae6423623f63
# ╠═041ba6c2-0a5c-4a59-b536-95db6e68716f
# ╠═c6224653-9dcd-4e17-9999-ccb69da6c222
# ╠═13c53be9-4922-4872-8ebf-8ed147a72fff
# ╠═4be2b3c0-98df-4d25-8200-a051a5dba4e9
# ╠═e3ecf2c5-67e9-4571-8743-04c67755a23b
# ╠═4b4938c2-bf21-455e-88d5-7d3843a6e005
# ╠═cae30107-1aee-48c0-8104-720cf5535b02
# ╟─925be0a1-0afa-4b5a-b9c9-003833e80a28
# ╟─ae9d8d21-1269-4053-89a3-3a8287a4ca70
# ╟─4c8b9b51-5e72-42a7-a919-ca2c2d1c5294
# ╟─7178c98a-15d1-4f10-9850-2079a6fdaa27
# ╟─90e0d677-a837-4a17-8c64-68ed350dd3c8
# ╟─26d120cb-0bd4-432d-9ce1-fc9a23a5067d
# ╟─65ca1318-8f62-4133-ab14-f3d3e1190029
# ╟─1b31eabc-3b6f-4d82-9cbb-62f1c53408de
# ╟─7ec71b3a-4c34-4803-923f-2effeff7cfcf
# ╟─84e4a23c-c926-4c86-ab11-07a08ef7d1b3
# ╠═e013dcb6-ea17-4fd5-af4d-318188600107
