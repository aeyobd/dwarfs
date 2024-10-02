### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ b918a965-9e54-42a4-9126-4dd614024ff5
using StatsBase: median, weights, mad

# ╔═╡ 9c7035e7-c1e7-40d5-8ab6-38f0bb682111
md"""
# Tidal Tails
A detailed analysis of the stars in sculptor
"""

# ╔═╡ 0e9e2085-bf96-4f9d-a758-8af36edf02da
import DensityEstimators as DE

# ╔═╡ a4fa1e76-8c2d-4402-b612-2f454bd06b8b
models_dir = "/arc7/home/dboyea/dwarfs/analysis/sculptor"

# ╔═╡ 82e8f2e4-d3ea-43c5-8813-aaebbca71cda
r_b_arcmin = 120

# ╔═╡ d0d1ecad-4a8d-4c1a-af2b-49f0d3d16bf2
model_dir = "$models_dir/1e6_V31_r3.2/orbit1/"

# ╔═╡ cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
starsfile = "$model_dir/stars/exp2d_rs0.10/final.fits"

# ╔═╡ a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
obs_today_filename = "/astro/dboyea/dwarfs/observations/sculptor/observed_properties.toml"

# ╔═╡ c4008b83-b61b-4baa-9fd5-9cced2dc6be8
import TOML

# ╔═╡ e37559b2-229c-4a37-b516-c6cb7c022b71
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 89a34cdf-c724-4c7e-ae42-11da6077feb2
function add_xi_eta_p!(df)
	xi_p, eta_p = lguys.to_orbit_coords(df.ra, df.dec, orbit_props["ra0"], orbit_props["dec0"], orbit_props["theta0"])

	df[!, :xi_p] = xi_p
	df[!, :eta_p] = eta_p

	return df
end

# ╔═╡ 7a92c896-7552-4f35-9761-5709d23e9adf
stars = lguys.read_fits(starsfile) |> add_xi_eta_p!

# ╔═╡ 6c76cfae-928b-47b3-babe-b0b9a5d68e65
obs_cen = stars[1, :]

# ╔═╡ 075ae901-bbbd-4d10-9d91-c393fc86a8e7
sky_orbit = lguys.read_fits(joinpath(model_dir, "skyorbit.fits")) |> add_xi_eta_p!

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

# ╔═╡ 6d0cff99-3efb-4405-a11d-f13200fa5334
idx_f = orbit_props["idx_f"]

# ╔═╡ 816b9db9-26c6-4ac8-9a46-82209d2cdc85
idx_orbit = idx_f - 5: idx_f + 5

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
		label="stellar density"
	)

	fig
end

# ╔═╡ 2c9e8ad0-bbff-4086-974e-89269868e324
sky_orbit[idx_f, :].pmra, sky_orbit[idx_f, :].pmdec

# ╔═╡ a89adc86-67a0-453f-b382-96f721f74d39
diff(sky_orbit.ra)[idx_f], diff(sky_orbit.dec)[idx_f]

# ╔═╡ 9419722f-9c65-45a7-b9c5-c666b22bf70e
sind(orbit_props["theta0"] ), cosd(orbit_props["theta0"])

# ╔═╡ 12c8d1f6-30fb-4616-b3dd-eb0aefd7d450
md"""
# Transform to stream coordinates
"""

# ╔═╡ 8dd15727-d9b0-47b2-a031-d7765ebd3fba
orbit_props

# ╔═╡ 7688c0ac-be70-470c-9db5-f2ed937e2fb5
function xi_eta_axis(dx=10, dy=5; kwargs...)
	fig = Figure(;kwargs...)
	

	limits = ((-dx, dx), (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel=L"$\xi'$ / degrees",
		ylabel=L"$\eta'$ / degrees",
		limits=limits,
		aspect = DataAspect(),
		xgridvisible=false,
		ygridvisible=false
		
	)

	return fig, ax
end

# ╔═╡ d1564bf3-a3c4-410b-b61b-69ac56b618e5
let 
	fig, ax = xi_eta_axis()
	
	bins = 100
	limits = ax.limits.val
	x = stars.xi_p
	y = stars.eta_p


	hi = Arya.histogram2d(x, y, bins, weights=stars.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-10, maximum(hi.values)))

	
	Colorbar(fig[1, 2], h,
		label="stellar density"
	)
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

# ╔═╡ 9c1839d6-1076-4598-8da6-49c02ec10580
sum(stars.weights[filt_cen])

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
	scatter!(x[filt], y[filt], label="centre", alpha=0.03, markersize=3)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=0.1, markersize=5)

	filt = filt_leading .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=0.1, markersize=5)

	filt = filt_excl
	scatter!(x[filt], y[filt], label="excluded", alpha=0.1, markersize=5)
	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ 443ac755-70fc-4193-945d-0e98622c5919
let 
	fig, ax = ra_dec_axis(r_centre / 60)
	ax.title = "central stars cut"

	x = stars.ra
	y = stars.dec
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=1, markersize=3)
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
		title="",
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
		title="",
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
		title="",
	aspect=DataAspect())

	bins = 100
	limits = ax.limits.val
	x = stars.pmra
	y = stars.pmdec
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=0.03, markersize=3)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=1, markersize=5)

	filt = filt_leading .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=1, markersize=5)

	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ 51e4628a-a799-48a3-9ad0-b5c2bd7a8d87
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="",
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
		xlabel=L"{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
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
	hist2d!(stars.pmra[filt], stars.pmdec[filt], weights=stars.weights[filt], limits=limits, colorscale=log10, colorrange=(1e-7, nothing), bins=100)

	errscatter!(obs_today.pmra, obs_today.pmdec)
	

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

# ╔═╡ 92b761ba-dc5d-4b99-9462-2c9b6faf680d
let	
	mass = stars.weights
	v_rad = stars.radial_velocity
	x = stars.xi_p
	r_max = 120

	limits = (-r_max / 60, r_max/60, obs_today.radial_velocity - 40, obs_today.radial_velocity + 40)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"{v}_\textrm{los} / \textrm{km\,s^{-1}}"
	)
	
	h = Arya.hist2d!(ax, x, v_rad, weights=mass, bins=100, limits=limits,
		colorrange=(1e-10, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ 1d4a7103-3d37-4597-a180-88a4c3733497
bins = -0.4r_max/60:0.5:0.4r_max/60

# ╔═╡ f1af9a92-a80d-4b7c-ba66-4f8cd43e0157
let	
	mass = stars.weights
	d = stars.distance
	x = stars.xi_p

	limits = (-r_max / 60, r_max/60, obs_cen.distance - 20, obs_cen.distance + 20)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel="distance / kpc"
	)
	
	h = Arya.hist2d!(ax, x, d, weights=mass, bins=100, limits=limits,
		colorrange=(1e-10, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ 302b3ac4-94cc-44d9-a5e6-314460e76a9f
let	
	mass = stars.weights
	v_rad = stars.pmra
	x = stars.xi_p

	dy = 0.1

	limits = (-r_max / 60, r_max/60, obs_cen.pmra - dy, obs_cen.pmra + dy)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\alpha}/ \textrm{mas\,yr^{-1}}"
	)
	
	h = Arya.hist2d!(ax, x, v_rad, weights=mass, bins=100, limits=limits,
		colorrange=(1e-15, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

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

		y_l[i] = y_m[i] - s
		y_h[i] = y_m[i] + s
		
	end

	return y_m, y_l, y_h
end

# ╔═╡ e3ecf2c5-67e9-4571-8743-04c67755a23b
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.radial_velocity[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"{v}_\textrm{los} / \textrm{km\,s^{-1}}"
	)


	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ 1f58f320-7bd6-40d6-b66b-0ec0db29de4a
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.radial_velocity_gsr[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{v}_\textrm{los} / \textrm{km\,s^{-1}}"
	)


	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ 0b88d438-0666-4875-b2f6-66480b652617
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.distance[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel="distance / kpc"
	)

	bins = LinRange(-0.5r_max/60, 0.5r_max/60, 20)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ 90f9d044-d03d-4e0f-a89e-c67b97ce9fb0
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.pmra_gsr[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\alpha*}/ \textrm{mas\,yr^{-1}}"
	)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ 380cef17-05a1-4111-8479-d88c5757a62e
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.pmra[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"{\mu}_{\alpha*}/ \textrm{mas\,yr^{-1}}"
	)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ a293ddc2-461c-42d6-9d02-7077c66163e5
let	
	w = stars.weights
	y = stars.pmra_gsr
	x = stars.eta_p

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\eta' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\alpha*}/ \textrm{mas\,yr^{-1}}"
	)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
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

# ╔═╡ b79ea320-c652-4850-a617-58bfb0c09be8
let	
	mass = stars.weights
	y = stars.pmdec_gsr
	x = stars.xi_p
	y_m = median(y, weights(mass))

	
	limits = (-r_max / 60, r_max/60, y_m - 0.1, y_m + 0.1)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\delta*}/ \textrm{mas\,yr^{-1}}"
	)
	
	h = Arya.hist2d!(ax, x, y, weights=mass, bins=100, limits=limits,
		colorrange=(1e-15, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ 3ebaa8ee-4ddf-49e8-b9ac-f0df9b45876b
let	
	mass = stars.weights
	y = stars.pmdec
	x = stars.xi_p
	y_m = median(y, weights(mass))

	
	limits = (-r_max / 60, r_max/60, y_m - 0.1, y_m + 0.1)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"{\mu}_{\delta*}/ \textrm{mas\,yr^{-1}}"
	)
	
	h = Arya.hist2d!(ax, x, y, weights=mass, bins=100, limits=limits,
		colorrange=(1e-15, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ de36ca44-5c34-4673-a75f-e705e5cd83a8
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.pmdec_gsr[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\delta}/ \textrm{mas\,yr^{-1}}"
	)

	bins = LinRange(-0.5r_max/60, 0.5r_max/60, 20)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ 14b0f654-81e2-4f29-95d3-30918293d852
let	
	filt = abs.(stars.eta_p) .< 2
	w = stars.weights[filt]
	y = stars.pmdec[filt]
	x = stars.xi_p[filt]

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"{\mu}_{\delta}/ \textrm{mas\,yr^{-1}}"
	)

	bins = LinRange(-0.5r_max/60, 0.5r_max/60, 20)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

	fig
end

# ╔═╡ d3d1e03d-d99a-4f79-818d-5bdc5b4a005a
let	
	w = stars.weights
	y = stars.pmdec_gsr
	x = stars.eta_p

	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\eta' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\delta}/ \textrm{mas\,yr^{-1}}"
	)

	x_m  = midpoints(bins)
	y_m, y_l, y_h = binned_median(x, y, bins, w=w)

	scatter!(x_m, y_m)
	errorbars!(x_m, y_m, y_m-y_l, y_h - y_m)

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
let
	fig = Figure()
	ax = PolarAxis(fig[1, 1])
	
	scatter!(deg2rad.(stars.xi_p), stars.distance, alpha=0.1)
	scatter!(deg2rad.(stars.xi_p[filt_dist]), stars.distance[filt_dist], alpha=0.1)

	fig
end

# ╔═╡ 90e0d677-a837-4a17-8c64-68ed350dd3c8
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

# ╔═╡ 26d120cb-0bd4-432d-9ce1-fc9a23a5067d
md"""
The next two plots validate the rotated coordinate system.
"""

# ╔═╡ 65ca1318-8f62-4133-ab14-f3d3e1190029
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

# ╔═╡ 1b31eabc-3b6f-4d82-9cbb-62f1c53408de
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

# ╔═╡ 7ec71b3a-4c34-4803-923f-2effeff7cfcf
md"""
The orbit in the rotated coordinate system. We expect the orbit to be ~ a horizontal line and increasing (lighter color) going to positive ``\xi'``
"""

# ╔═╡ 84e4a23c-c926-4c86-ab11-07a08ef7d1b3
let
	fig, ax = xi_eta_axis()

	scatter!(stars.xi_p, stars.eta_p)
	p = lines!(sky_orbit.xi_p[idx_orbit], sky_orbit.eta_p[idx_orbit], 
		color=idx_orbit, linewidth=4)

	Colorbar(fig[1, 2], p, label="time / Gyr")
	fig
end

# ╔═╡ e013dcb6-ea17-4fd5-af4d-318188600107
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

# ╔═╡ Cell order:
# ╟─9c7035e7-c1e7-40d5-8ab6-38f0bb682111
# ╠═fb8bb8ba-34ad-11ef-23e6-1d890b60e0b9
# ╠═0e9e2085-bf96-4f9d-a758-8af36edf02da
# ╠═b918a965-9e54-42a4-9126-4dd614024ff5
# ╠═a4fa1e76-8c2d-4402-b612-2f454bd06b8b
# ╠═82e8f2e4-d3ea-43c5-8813-aaebbca71cda
# ╠═d0d1ecad-4a8d-4c1a-af2b-49f0d3d16bf2
# ╠═cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
# ╠═7a92c896-7552-4f35-9761-5709d23e9adf
# ╠═6c76cfae-928b-47b3-babe-b0b9a5d68e65
# ╠═89a34cdf-c724-4c7e-ae42-11da6077feb2
# ╠═075ae901-bbbd-4d10-9d91-c393fc86a8e7
# ╠═a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
# ╠═c4008b83-b61b-4baa-9fd5-9cced2dc6be8
# ╠═e37559b2-229c-4a37-b516-c6cb7c022b71
# ╠═0ba05fdb-f859-4381-b2d0-145aa04f7bbf
# ╠═4ef955fb-a813-46ad-8f71-7c8a3d371eee
# ╠═910267ee-aa39-4f04-961b-70f0862d27e2
# ╠═4c1d6f6f-e257-4126-9b4a-8e5aa8470295
# ╠═6d0cff99-3efb-4405-a11d-f13200fa5334
# ╠═e9e35643-168e-4e87-a880-6831b46145c7
# ╠═816b9db9-26c6-4ac8-9a46-82209d2cdc85
# ╠═2c9e8ad0-bbff-4086-974e-89269868e324
# ╠═a89adc86-67a0-453f-b382-96f721f74d39
# ╠═9419722f-9c65-45a7-b9c5-c666b22bf70e
# ╟─12c8d1f6-30fb-4616-b3dd-eb0aefd7d450
# ╠═8dd15727-d9b0-47b2-a031-d7765ebd3fba
# ╠═7688c0ac-be70-470c-9db5-f2ed937e2fb5
# ╠═d1564bf3-a3c4-410b-b61b-69ac56b618e5
# ╠═49705acd-d623-4860-a5e7-759f112553ab
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
# ╠═9c1839d6-1076-4598-8da6-49c02ec10580
# ╟─c57edba1-09f9-4fe9-bc73-5746884ee7c2
# ╟─119040e2-ef45-4fc2-809d-aec333e276d0
# ╟─443ac755-70fc-4193-945d-0e98622c5919
# ╟─1d35a894-1eca-4ee0-9d1a-6ac4704b912c
# ╠═6b0a0ea3-f6f7-4cee-be11-bce6872ab870
# ╠═54205139-3c7b-4beb-b9cb-c87b272df58a
# ╠═91955a7b-009e-4afa-aca1-312f4a779bb9
# ╠═19bdc540-63a8-46ad-8d6f-a425d64cdd81
# ╠═51e4628a-a799-48a3-9ad0-b5c2bd7a8d87
# ╠═e098260a-718a-4f43-b12c-cab0a1302b90
# ╠═be0158cf-c322-4147-8f01-8b6a715bc0dc
# ╠═5eb22695-57c4-4ebf-b412-588f5e366bf5
# ╠═79dd7e3d-8564-4389-960b-05512b0143c0
# ╠═d4d50209-e9a8-40dc-9f45-4e8000f70b39
# ╠═92b761ba-dc5d-4b99-9462-2c9b6faf680d
# ╠═1d4a7103-3d37-4597-a180-88a4c3733497
# ╠═e3ecf2c5-67e9-4571-8743-04c67755a23b
# ╠═1f58f320-7bd6-40d6-b66b-0ec0db29de4a
# ╠═f1af9a92-a80d-4b7c-ba66-4f8cd43e0157
# ╠═0b88d438-0666-4875-b2f6-66480b652617
# ╟─302b3ac4-94cc-44d9-a5e6-314460e76a9f
# ╠═90f9d044-d03d-4e0f-a89e-c67b97ce9fb0
# ╠═380cef17-05a1-4111-8479-d88c5757a62e
# ╠═a293ddc2-461c-42d6-9d02-7077c66163e5
# ╠═925be0a1-0afa-4b5a-b9c9-003833e80a28
# ╠═f1cb5d80-57a8-43e8-8819-4e1bc134edea
# ╠═4b0515a9-e55d-4bed-b6d3-0a2f12d0a238
# ╠═20e6efaa-95ab-4163-a584-f1fea47f02c7
# ╠═b79ea320-c652-4850-a617-58bfb0c09be8
# ╠═3ebaa8ee-4ddf-49e8-b9ac-f0df9b45876b
# ╠═de36ca44-5c34-4673-a75f-e705e5cd83a8
# ╠═14b0f654-81e2-4f29-95d3-30918293d852
# ╠═d3d1e03d-d99a-4f79-818d-5bdc5b4a005a
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
