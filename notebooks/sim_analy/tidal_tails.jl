### A Pluto.jl notebook ###
# v0.19.43

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

# ╔═╡ 20ac4d7c-d835-4a7f-9600-261e43f4b290
using FITSIO

# ╔═╡ 9c7035e7-c1e7-40d5-8ab6-38f0bb682111
md"""
# Tidal Tails
A detailed analysis of the stars in sculptor
"""

# ╔═╡ a4fa1e76-8c2d-4402-b612-2f454bd06b8b
models_dir = "/arc7/home/dboyea/sculptor"

# ╔═╡ 82e8f2e4-d3ea-43c5-8813-aaebbca71cda
r_b_arcmin = 64

# ╔═╡ d0d1ecad-4a8d-4c1a-af2b-49f0d3d16bf2
model_dir = "$models_dir/orbits/V50_r0.5/"

# ╔═╡ cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
starsfile = "$model_dir/stars/exp2d_rs0.08_today.fits"

# ╔═╡ 7a92c896-7552-4f35-9761-5709d23e9adf
stars = lguys.load_fits(starsfile)

# ╔═╡ 6c76cfae-928b-47b3-babe-b0b9a5d68e65
obs_cen = stars[1, :]

# ╔═╡ 075ae901-bbbd-4d10-9d91-c393fc86a8e7
sky_orbit = lguys.load_fits(joinpath(model_dir, "skyorbit.fits"))

# ╔═╡ a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
obs_today_filename = "/astro/dboyea/dwarfs/sculptor_obs_properties.toml"

# ╔═╡ c4008b83-b61b-4baa-9fd5-9cced2dc6be8
import TOML

# ╔═╡ e37559b2-229c-4a37-b516-c6cb7c022b71
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 0ba05fdb-f859-4381-b2d0-145aa04f7bbf
obs_today_file = TOML.parsefile(obs_today_filename)

# ╔═╡ 4ef955fb-a813-46ad-8f71-7c8a3d371eee
obs_today_icrs = lguys.ICRS(;
	ra=obs_today_file["ra"], dec=obs_today_file["dec"],
	distance=obs_today_file["distance"],
	pmra=obs_today_file["pmra"],
	pmdec=obs_today_file["pmdec"],
	radial_velocity=obs_today_file["radial_velocity"],
)

# ╔═╡ c3a3129e-18c6-4348-81c0-b03c5836b785
frame = lguys.GSR

# ╔═╡ 910267ee-aa39-4f04-961b-70f0862d27e2
obs_today = lguys.transform(frame, obs_today_icrs)

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
		ygridvisible=false
		
	)

	return fig, ax
end

# ╔═╡ 6d0cff99-3efb-4405-a11d-f13200fa5334
idx_f = orbit_props["idx_f"]

# ╔═╡ 816b9db9-26c6-4ac8-9a46-82209d2cdc85
idx_orbit = idx_f - 20: idx_f + 20

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

# ╔═╡ 12c8d1f6-30fb-4616-b3dd-eb0aefd7d450
md"""
# Transform to stream coordinates
"""

# ╔═╡ 6b516488-80d7-4904-9331-6161e5829638
"""
Given the slope of an orbit in ddec/dra, calculates a rotated
sky frame centred on RA, DEC and with the x-axis aligned with the orbit.
"""
function to_orbit_coords(ra, dec, ra0, dec0, PA)

	# want to rotate to dec, ra of centre, then rotate 

	α = deg2rad(ra0)
	δ = deg2rad(dec0)
	ϖ = deg2rad(90 - PA)
	Rmat = lguys.Rx_mat(ϖ) * lguys.Ry_mat(δ) * lguys.Rz_mat(-α) 

	coords = lguys.unit_vector(ra, dec)
	coords =  Rmat * coords' 
	skycoords = lguys.cartesian_to_sky(coords[1, :], coords[2, :], coords[3, :])[:, 1:2]

	ra = skycoords[:, 1] 
	ra .-= 360 * (ra .> 180)
	dec = skycoords[:, 2]

	ra, dec
end

# ╔═╡ 8dd15727-d9b0-47b2-a031-d7765ebd3fba
orbit_props

# ╔═╡ a8a4442e-b180-4ffe-8f91-522e45cac47e
xi_p, eta_p = to_orbit_coords(stars.ra, stars.dec, orbit_props["ra0"], 
orbit_props["dec0"], orbit_props["theta0"])

# ╔═╡ 22afcc2d-95d1-4999-9b42-0a7c75d054c4
xi_p_orbit, eta_p_orbit = to_orbit_coords(sky_orbit.ra, sky_orbit.dec, orbit_props["ra0"], 
orbit_props["dec0"], orbit_props["theta0"])

# ╔═╡ b9998808-90e3-4bdf-a060-5c4c22df745f
begin
	sky_orbit[:, :xi_p] = xi_p_orbit
	sky_orbit[:, :eta_p] = eta_p_orbit

	stars[:, :xi_p] = xi_p
	stars[:, :eta_p] = eta_p
end

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
filt_cen = stars.r_ell .< 120

# ╔═╡ 2cb76b38-b768-4184-b73c-e8fb350351d8
filt_leading = not.(filt_cen) .& (stars.xi_p .> 0)

# ╔═╡ c0166054-7c1f-474f-a4bc-3b122034e923
filt_trailing = not.(filt_cen) .& (stars.xi_p .< 0)

# ╔═╡ aa157085-0705-44ed-9d65-5626658d71e7
filt_dist = stars.r_ell .< r_max

# ╔═╡ d4d5c328-b135-4ddb-8d12-f8f13ba1da3d
scatter(stars.xi_p[filt_dist], stars.eta_p[filt_dist])

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
	fig, ax = ra_dec_axis(2r_centre / 60)

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
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
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
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
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

# ╔═╡ 19bdc540-63a8-46ad-8d6f-a425d64cdd81
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="",
	aspect=DataAspect())

	bins = 100
	limits = ax.limits.val
	x = stars.pmra
	y = stars.pmdec
	scatter!(x[filt_cen], y[filt_cen], label="centre", alpha=0.03, markersize=3)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="trailing arm", alpha=1, markersize=5)

	filt = filt_trailing .& filt_dist
	scatter!(x[filt], y[filt], label="leading arm", alpha=1, markersize=5)

	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ e098260a-718a-4f43-b12c-cab0a1302b90
let 
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
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
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
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
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="leading arm",
		aspect=DataAspect()
	)

	filt = filt_leading .& filt_dist
	hist2d!(stars.pmra[filt], stars.pmdec[filt], weights=stars.weights[filt], limits=limits, colorscale=log10, colorrange=(1e-15, nothing), bins=100)

	scatter!(obs_cen.pmra, obs_cen.pmdec)
	

	fig
end

# ╔═╡ fc3a1a85-c317-4b45-b9c4-50d15e0eb9da
obs_cen.pmra

# ╔═╡ 79dd7e3d-8564-4389-960b-05512b0143c0
let
	fig = Figure()
	
	dr = 0.1
	limits = (obs_cen.pmra - dr, obs_cen.pmra + dr, obs_cen.pmdec - dr, obs_cen.pmdec + dr)
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
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

	limits = (-r_max / 60, r_max/60, 30, 110)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{v}_\textrm{los} / \textrm{km\,s^{-1}}"
	)
	
	h = Arya.hist2d!(ax, x, v_rad, weights=mass, bins=100, limits=limits,
		colorrange=(1e-10, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ f1af9a92-a80d-4b7c-ba66-4f8cd43e0157
let	
	mass = stars.weights
	d = stars.distance
	x = stars.xi_p

	limits = (-r_max / 60, r_max/60, 70, 100)
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

	limits = (-r_max / 60, r_max/60, -0.35, -0.15)
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

# ╔═╡ b79ea320-c652-4850-a617-58bfb0c09be8
let	
	mass = stars.weights
	v_rad = stars.pmdec
	x = stars.xi_p

	limits = (-r_max / 60, r_max/60, 0.2, 0.4)
	fig = Figure()
	ax = Axis(fig[1,1],
		#limits=limits,
		xlabel=L"\xi' / \textrm{degrees}",
		ylabel=L"\tilde{\mu}_{\delta*}/ \textrm{mas\,yr^{-1}}"
	)
	
	h = Arya.hist2d!(ax, x, v_rad, weights=mass, bins=100, limits=limits,
		colorrange=(1e-15, nothing), colorscale=log10, normalization=:density
	)

	Colorbar(fig[1, 2], h, label="stellar mass density")
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
	scatter!(eta_p[filt_dist], stars.distance[filt_dist], alpha=0.1)

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

	scatter!(xi_p, eta_p)
	p = scatter!(xi_p_orbit[idx_orbit], eta_p_orbit[idx_orbit], color=idx_orbit)

	Colorbar(fig[1, 2], p, label="time / Gyr")
	fig
end

# ╔═╡ 6c611b16-a7db-44ad-8656-2fc76f7c989c
scatter(stars.xi_p, stars.eta_p)

# ╔═╡ Cell order:
# ╟─9c7035e7-c1e7-40d5-8ab6-38f0bb682111
# ╠═fb8bb8ba-34ad-11ef-23e6-1d890b60e0b9
# ╠═a4fa1e76-8c2d-4402-b612-2f454bd06b8b
# ╠═82e8f2e4-d3ea-43c5-8813-aaebbca71cda
# ╠═d0d1ecad-4a8d-4c1a-af2b-49f0d3d16bf2
# ╠═cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
# ╠═20ac4d7c-d835-4a7f-9600-261e43f4b290
# ╠═7a92c896-7552-4f35-9761-5709d23e9adf
# ╠═6c76cfae-928b-47b3-babe-b0b9a5d68e65
# ╠═075ae901-bbbd-4d10-9d91-c393fc86a8e7
# ╠═a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
# ╠═c4008b83-b61b-4baa-9fd5-9cced2dc6be8
# ╠═e37559b2-229c-4a37-b516-c6cb7c022b71
# ╠═0ba05fdb-f859-4381-b2d0-145aa04f7bbf
# ╠═4ef955fb-a813-46ad-8f71-7c8a3d371eee
# ╠═c3a3129e-18c6-4348-81c0-b03c5836b785
# ╠═910267ee-aa39-4f04-961b-70f0862d27e2
# ╠═4c1d6f6f-e257-4126-9b4a-8e5aa8470295
# ╠═6d0cff99-3efb-4405-a11d-f13200fa5334
# ╠═e9e35643-168e-4e87-a880-6831b46145c7
# ╠═816b9db9-26c6-4ac8-9a46-82209d2cdc85
# ╟─12c8d1f6-30fb-4616-b3dd-eb0aefd7d450
# ╠═6b516488-80d7-4904-9331-6161e5829638
# ╠═8dd15727-d9b0-47b2-a031-d7765ebd3fba
# ╠═a8a4442e-b180-4ffe-8f91-522e45cac47e
# ╠═22afcc2d-95d1-4999-9b42-0a7c75d054c4
# ╠═b9998808-90e3-4bdf-a060-5c4c22df745f
# ╠═7688c0ac-be70-470c-9db5-f2ed937e2fb5
# ╠═d1564bf3-a3c4-410b-b61b-69ac56b618e5
# ╠═49705acd-d623-4860-a5e7-759f112553ab
# ╠═9241f6a2-d829-4370-a497-35f6fb3358dd
# ╠═b09e9622-b522-4308-ad95-c939ed0a5315
# ╠═95b7627e-c17a-4be6-ac66-c3212a8feb23
# ╠═d4d5c328-b135-4ddb-8d12-f8f13ba1da3d
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
# ╠═443ac755-70fc-4193-945d-0e98622c5919
# ╟─1d35a894-1eca-4ee0-9d1a-6ac4704b912c
# ╠═6b0a0ea3-f6f7-4cee-be11-bce6872ab870
# ╠═54205139-3c7b-4beb-b9cb-c87b272df58a
# ╠═19bdc540-63a8-46ad-8d6f-a425d64cdd81
# ╠═e098260a-718a-4f43-b12c-cab0a1302b90
# ╠═be0158cf-c322-4147-8f01-8b6a715bc0dc
# ╠═5eb22695-57c4-4ebf-b412-588f5e366bf5
# ╠═fc3a1a85-c317-4b45-b9c4-50d15e0eb9da
# ╠═79dd7e3d-8564-4389-960b-05512b0143c0
# ╠═d4d50209-e9a8-40dc-9f45-4e8000f70b39
# ╠═92b761ba-dc5d-4b99-9462-2c9b6faf680d
# ╠═f1af9a92-a80d-4b7c-ba66-4f8cd43e0157
# ╠═302b3ac4-94cc-44d9-a5e6-314460e76a9f
# ╠═b79ea320-c652-4850-a617-58bfb0c09be8
# ╟─ae9d8d21-1269-4053-89a3-3a8287a4ca70
# ╟─4c8b9b51-5e72-42a7-a919-ca2c2d1c5294
# ╟─7178c98a-15d1-4f10-9850-2079a6fdaa27
# ╟─90e0d677-a837-4a17-8c64-68ed350dd3c8
# ╟─26d120cb-0bd4-432d-9ce1-fc9a23a5067d
# ╟─65ca1318-8f62-4133-ab14-f3d3e1190029
# ╟─1b31eabc-3b6f-4d82-9cbb-62f1c53408de
# ╟─7ec71b3a-4c34-4803-923f-2effeff7cfcf
# ╟─84e4a23c-c926-4c86-ab11-07a08ef7d1b3
# ╠═6c611b16-a7db-44ad-8656-2fc76f7c989c
