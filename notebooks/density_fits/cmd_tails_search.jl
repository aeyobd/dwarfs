### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ a27c968e-5fde-11ef-27c0-f5f22481829a
begin
	using Pkg; Pkg.activate()

	using Arya, CairoMakie
	
	using LilGuys
end

# ╔═╡ f0bdd44b-0bcb-41ab-966f-afd514142629
include("read_iso.jl")

# ╔═╡ 0e59c6c9-f053-439f-a6b5-5081e7a7a1a2
include("filter_utils.jl")

# ╔═╡ 003e78ba-37b7-4650-ab1e-ba7db2a15e70
import LinearAlgebra: normalize, ×

# ╔═╡ a5329eeb-9ec4-49ed-9c36-6307a67e3e83
iso_dir = "/astro/dboyea/dwarfs/data/MIST_v1.2_vvcrit0.0_UBVRIplus/"

# ╔═╡ fd364bd0-3740-4cb7-941f-10fe47143778
readdir(iso_dir)

# ╔═╡ 506510de-af00-4668-871e-ea00ac6810ef
isochrone = ISOCMD(joinpath(iso_dir, "MIST_v1.2_feh_m1.75_afe_p0.0_vvcrit0.0_UBVRIplus.iso.cmd"))

# ╔═╡ 5438d9d9-0256-4779-8c64-56b560c2f730
filt_params = DensityParams(read_file("sculptor/simple.toml"))

# ╔═╡ 1629ba90-c0ba-477f-9bc4-897bf61426fe
function cmd_axis(gs)
	return Axis(gs,
		xlabel = "Bp-Rp",
		ylabel = "G",
		yreversed=true,
		limits = (-0.5, 2, 16, 21),
		xgridvisible=false,
		ygridvisible=false,
	)
end

# ╔═╡ 784654ff-1596-49eb-8c38-a28200a8bd0c
begin
	iso = copy(isochrone[log10(12e9)])

	iso[!, :G] = iso.Gaia_G_EDR3
	iso[!, :bp_rp] = iso.Gaia_BP_EDR3 .- iso.Gaia_RP_EDR3

	iso = iso[iso.phase .<= 2, :]
end

# ╔═╡ d4e6fe00-e886-4807-8c92-0ce50fb124be
obs_props = TOML.parsefile("/astro/dboyea/dwarfs/sculptor_obs_properties.toml")

# ╔═╡ c4a62177-c50a-4800-8fbc-d443216c571e
orbit_props = TOML.parsefile("/astro/dboyea/sculptor/orbits/orbit1/1e6/V32_r5/orbital_properties.toml")

# ╔═╡ f33adaef-228d-483e-938d-0d2ccd828ee8
idx_f = orbit_props["idx_f"]

# ╔═╡ 8e5e7421-042c-4876-a992-9f847f1e7ba8
begin 
	orbit = lguys.load_fits("/astro/dboyea/sculptor/orbits/orbit1/1e6/V32_r5/skyorbit.fits")

	orbit[!, :xi], orbit[!, :eta] = lguys.to_tangent(orbit.ra, orbit.dec, orbit.ra[idx_f], orbit.dec[idx_f])

	orbit
end

# ╔═╡ d5ef30f4-7a58-41bd-9c42-62689c257bcd
idx_orbit=1775:1800

# ╔═╡ 5b791110-87d8-4139-abeb-91f6e2688304
dm = 5*log10(83.2 * 1e3) - 5

# ╔═╡ a3e53d78-40ac-43f2-b5a9-b3b7d5990a1d
dm_err = 0.05

# ╔═╡ 9438a961-f66f-4be3-885c-74fb0f27661d
lguys

# ╔═╡ 8ddacf0c-91f1-4543-ac2e-e37871c977f8
begin 
	all_stars = LilGuys.load_fits("sculptor/gaia_4deg_cen.fits")
	all_stars[:, :G] = all_stars.phot_g_mean_mag

	all_stars[:, :xi], all_stars[:, :eta] = lguys.to_tangent(all_stars.ra, all_stars.dec, filt_params.ra, filt_params.dec)

	all_stars[:, :xi_p], all_stars[:, :eta_p] = lguys.to_orbit_coords(all_stars.ra, all_stars.dec, filt_params.ra, filt_params.dec, orbit_props["theta0"])
	all_stars
end

# ╔═╡ 33ca51be-05cc-4d27-92ba-d7979aa44f7e
sum(isnan.(all_stars.bp_rp))

# ╔═╡ 5708825c-36c6-4937-b673-21ae1cd0449c
cmd_x = [filt_params.cmd_cut[1:2:end]; filt_params.cmd_cut[1]]

# ╔═╡ 340d6d44-2204-42f5-847f-0f9645e4df2c
cmd_y = [filt_params.cmd_cut[2:2:end]; filt_params.cmd_cut[2]]

# ╔═╡ 917516ad-0807-4fb1-9afb-d28044249c72
begin
	filt_basic = ruwe_filter(all_stars, filt_params)
	filt_basic .&= parallax_filter(all_stars, filt_params)
end

# ╔═╡ 3b0992aa-1d1a-42b9-9db4-ef32ea8d7a4c
let
	fig = Figure()

	ax = cmd_axis(fig[1, 1])

	filt = filt_basic

	scatter!(all_stars.bp_rp[filt], all_stars.G[filt], markersize=1, alpha=0.3, color=:black)

	#scatter!(all_stars.bp_rp[filt_cmd], all_stars.G[filt_cmd])
	#lines!(iso.bp_rp, iso.G .+ dm)

	poly!(cmd_x, cmd_y, color=:transparent, strokecolor=COLORS[1], strokewidth=2)
	fig
end

# ╔═╡ c8830c99-f0a2-4f3a-9080-e4ae539565a8
filt_pm = bivariate_z.(all_stars.pmra, all_stars.pmdec, filt_params.pmra, filt_params.pmdec, all_stars.pmra_error, all_stars.pmdec_error, all_stars.pmra_pmdec_corr) .< 6

# ╔═╡ 06f103a5-05ae-419d-af6f-21149d162af7
filt_cmd = cmd_filter(all_stars, filt_params)

# ╔═╡ 7b8487e3-b85c-433e-b8e8-47ef399caab9
let
	fig, ax = FigAxis(
		limits=(-2, 2, -2, 2)
	)
	scatter!(all_stars.pmra, all_stars.pmdec)
	scatter!(all_stars.pmra[filt_pm], all_stars.pmdec[filt_pm], markersize=6)

	fig
end

# ╔═╡ 32731e29-c988-437a-9e5f-ade7e6e9ffe9
filt_params

# ╔═╡ 0ab3e71f-8220-4293-bf74-95be5d443338
function get_mean_pm(df)
	pmra_cen = lguys.mean(df.pmra, lguys.weights(df.pmra_error .^ -2))
	pmra_cen_err = sqrt(1 / sum(df.pmra_error .^ -2))
	pmdec_cen = lguys.mean(df.pmdec, lguys.weights(df.pmdec_error .^ -2))
	pmdec_cen_err = sqrt(1 / sum(df.pmdec_error .^ -2))

	return pmra_cen, pmdec_cen, pmra_cen_err, pmdec_cen_err
end

# ╔═╡ fb8b12bf-a7f2-4f02-b1c8-edca5f2980ab
rs_test = [0.5, 1, 1.25, 1.5, 3.5]

# ╔═╡ d01537bf-6a71-47dc-b818-58693aac34fb
r_circ_cut = 0.25

# ╔═╡ ec205739-94e5-4ad1-a006-d20f89f160ec
function plot_cmd!(df; kwargs...)
	scatter!(df.bp_rp, df.G, alpha=1; kwargs...)
end

# ╔═╡ dfb0e6da-1871-4dc7-a2b9-b77f08982111
function select_region(allstars, centre; radius=0.5)
	x_cen, y_cen = centre

	r = @. sqrt((allstars.xi - x_cen)^2 + (allstars.eta - y_cen)^2)

	filt = r .< radius

	return allstars[filt, :]
end

# ╔═╡ 845dd99c-1534-47cb-a84d-42968bfbcb8d
function plot_cmd_members!(centre, radius)
	df = select_region(all_stars, centre, radius=radius)
	plot_cmd!(df, color=:grey, markersize=5, alpha=0.5)
	
	df = select_region(all_stars[filt_basic, :], centre, radius=radius)
	plot_cmd!(df)

	df = select_region(all_stars[filt_basic .& filt_cmd, :], centre, radius=radius)
	
	plot_cmd!(df, color=COLORS[2])

	N = size(df, 1)
	text!(0.05, 0.9, text="$N stars in CMD", color=:black, space=:relative)
end

# ╔═╡ d6ca26cb-e4f7-4d17-ba46-0e082fb02d14
function plot_pm_members!(centre, radius)
	df = select_region(all_stars[filt_basic, :], centre, radius=radius)
	scatter!(df.pmra, df.pmdec, color=:grey, markersize=5, alpha=0.5, label="parallax & quality cut")
	
	#df = select_region(all_stars[filt_basic, :], centre, radius=radius)
	#scatter!(df.pmra, df.pmdec, markersize=8, label="and pm cut")

	df = select_region(all_stars[filt_basic .& filt_cmd, :], centre, radius=radius)
	scatter!(df.pmra, df.pmdec, markersize=10, label="and cmd cut")

	pmra_cen = lguys.mean(df.pmra, lguys.weights(df.pmra_error .^ -2))
	pmra_cen_err = sqrt(1 / sum(df.pmra_error .^ -2))
	pmdec_cen = lguys.mean(df.pmdec, lguys.weights(df.pmdec_error .^ -2))
	pmdec_cen_err = sqrt(1 / sum(df.pmdec_error .^ -2))

	errscatter!([pmra_cen], [pmdec_cen], xerr=[pmra_cen_err], yerr=[pmdec_cen_err], color=:black)

end

# ╔═╡ 28bd7605-91d1-4c8c-a5f4-91af51a2d484
function pm_axis(gp; dpm=5, kwargs...)
	return Axis(gp, 
		xlabel=L"$\mu_{\alpha*}$ / mas\,yr$^{-1}$",
		ylabel=L"$\mu_\delta$ / mas\,yr$^{-1}$",
		aspect=DataAspect(),
		limits= dpm .* (-1, 1, -1, 1),
		xgridvisible=false,
		ygridvisible=false,
	)
end

# ╔═╡ b35026bd-8ca4-4552-8d68-2b8811f7b3d5
orbit_vector = normalize([
	orbit.xi[idx_orbit][end] - orbit.xi[idx_orbit][1],
	orbit.eta[idx_orbit][end] - orbit.eta[idx_orbit][1]
])

# ╔═╡ 6736916a-3484-4580-a909-6ba7e32e174a
let
	fig = Figure()
	ax = cmd_axis(fig[1, 1])
	
	plot_cmd_members!([0, 0], r_circ_cut)[1]
	fig
end

# ╔═╡ 4cb80651-8b11-47ed-bf0f-c858b420978c
let
	fig = Figure()
	ax = pm_axis(fig[1, 1], dpm=30)
	
	plot_pm_members!(0*orbit_vector, 4r_circ_cut)[1]

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 72af6da6-22d1-4bbd-b009-bd801133b271
let
	rs = rs_test

	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = cmd_axis(fig[i, 1])
		plot_cmd_members!(orbit_vector * r, r_circ_cut)
		
		text!(0.05, 0.8, text="$r degrees along orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ 89c96e28-14ba-4492-8f05-5474da8eef15
bg_vector = ([orbit_vector; 0] × [0, 0, 1])[1:2]

# ╔═╡ 01717a66-0e21-49e1-ba8e-4e44a0826bb2
let
	rs = rs_test

	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = cmd_axis(fig[i, 1])
		plot_cmd_members!(bg_vector * r, r_circ_cut)
		
		text!(0.05, 0.5, text="$r degrees ⟂ to orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ 7b1a2c03-9df2-4ec2-8eeb-64c9c9c2cdd4
let
	rs = rs_test

	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = pm_axis(fig[i, 1])
		plot_pm_members!(orbit_vector * r, r_circ_cut)
		
		text!(0.05, 0.8, text="$r degrees along orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ f826de98-01bb-4778-9a17-d2674a6cfae6
begin
	r_bins = -2:0.25:2

	pmra_means = Float64[]
	pmdec_means = Float64[]
	pmra_means_err = Float64[]
	pmdec_means_err = Float64[]
	pm_counts = Int[]
	for i in eachindex(r_bins)
		r = r_bins[i]
		df = select_region(all_stars[filt_cmd .& filt_pm .& filt_pm, :], orbit_vector * r, radius=r_circ_cut)
		x, y, xe, ye = get_mean_pm(df)
		push!(pmra_means, x)
		push!(pmdec_means, y)
		push!(pmra_means_err, xe)
		push!(pmdec_means_err, ye)
		push!(pm_counts, size(df, 1))
	end
end

# ╔═╡ 3bd33e91-89af-4e65-b25e-df33709823fd
members = all_stars[filt_cmd .& filt_pm .& filt_pm, :]

# ╔═╡ 0bea286f-b58c-4e91-8bfd-faf546f7ac10
let
	fig, ax = FigAxis(
		xlabel="distance along orbit / degrees", 
		ylabel="# of members in 0.25 deg radius circle",
		yscale=log10
	)

	scatter!(r_bins,pm_counts)

    fig
	

end

# ╔═╡ 086e9e5b-e5cd-4e4c-a48c-954792754115
eta_p_max = 0.25

# ╔═╡ 92f350c8-06bb-4625-a1ea-d27af8dc2864
begin
	r_bins_2 = -2.125:0.25:2.125
	
	pmra_means2 = Float64[]
	pmdec_means2 = Float64[]
	pmra_means_err2 = Float64[]
	pmdec_means_err2 = Float64[]
	pm_counts2 = Int[]
	for i in eachindex(r_bins_2)[1:end-1]
		filt = members.xi_p .> r_bins_2[i]
		filt .&= members.xi_p .< r_bins_2[i+1]
		filt .&= abs.(members.eta_p) .< eta_p_max
		
		df = members[filt, :]
		x, y, xe, ye = get_mean_pm(df)
		push!(pmra_means2, x)
		push!(pmdec_means2, y)
		push!(pmra_means_err2, xe)
		push!(pmdec_means_err2, ye)
		push!(pm_counts2, size(df, 1))
	end
end

# ╔═╡ 90b17eba-df5e-465f-a9d9-c8e2ba41b48c


# ╔═╡ 60e6c0c9-e382-4751-b0c9-e5062884b3c8
let

	fig, ax = FigAxis(
		xlabel="distance along orbit / degrees", 
		ylabel=L"$\mu_{\alpha*}$ / mas\,yr$^{-1}$",
		limits=(-1.1, 1.1, 0, 0.3)
	)

	errscatter!(r_bins, pmra_means, yerr=pmra_means_err)
	hlines!(filt_params.pmra)

    fig
	

end

# ╔═╡ ab81099c-5ffc-49d9-b46c-57cb6e6711c3
let

	fig, ax = FigAxis(
		xlabel=L"$\xi'$ / degrees", 
		ylabel=L"$\mu_{\alpha*}$ / mas\,yr$^{-1}$",
		limits=(-1.1, 1.1, 0, 0.3)
	)

	errscatter!(midpoints(r_bins_2), pmra_means2, yerr=pmra_means_err2)
	hlines!(filt_params.pmra)

    fig
	

end

# ╔═╡ 0d56d343-adef-488a-b494-2b9aa3f0bfbf
let

	fig, ax = FigAxis(
		xlabel=L"$\xi'$ / degrees", 
		ylabel=L"$\mu_{\alpha*}$ / mas\,yr$^{-1}$",
		limits=(-1.1, 1.1, -0.3, -0.1)
	)

	errscatter!(midpoints(r_bins_2), pmdec_means2, yerr=pmdec_means_err2)
	hlines!(filt_params.pmdec)

    fig
	

end

# ╔═╡ e898391c-2690-4a25-8d44-2b49ddf1a097
let

	fig, ax = FigAxis(
		xlabel="distance along orbit / degrees", 
		ylabel=L"$\mu_{\delta}$ / mas\,yr$^{-1}$",

	)

	errscatter!(r_bins, pmdec_means, yerr=pmdec_means_err)
	hlines!(filt_params.pmdec, color=:black)

    fig
	

end

# ╔═╡ e1a8309c-c626-4604-b2a6-7d11ef0472a4
let

	fig, ax = FigAxis(
		xlabel="distance along orbit / degrees", 
		ylabel=L"$\mu_{\delta}$ / mas\,yr$^{-1}$",
		limits=(-1.1, 1.1, -0.35, 0)

	)

	errscatter!(r_bins, pmdec_means, yerr=pmdec_means_err)
	hlines!(filt_params.pmdec, color=:black)

    fig
	

end

# ╔═╡ 3d981975-a2f7-45a7-92df-c464e56370b3
pmra_means_err

# ╔═╡ 6cd0161f-657f-446d-a046-f84fdb8d792f
let
	fig = Figure()
	ax = pm_axis(fig[1, 1])
	ax.limits=(-1, 1, -1, 1)
	
	errscatter!(pmra_means, pmdec_means, xerr=pmra_means_err, yerr=pmdec_means_err, color=:black, alpha=0.2)
	errscatter!(pmra_means, pmdec_means, color=r_bins, colormap=:bluesreds)
	fig
end

# ╔═╡ 7fd40134-633c-4cb9-87cd-2254c485e615
let
	fig = Figure()
	ax = pm_axis(fig[1, 1])
	ax.limits = (nothing, nothing)

	filt = abs.(r_bins) .< 0.8
	errscatter!(pmra_means[filt], pmdec_means[filt], xerr=pmra_means_err[filt], yerr=pmdec_means_err[filt], color=:black, alpha=0.2)
	p = errscatter!(pmra_means[filt], pmdec_means[filt], color=r_bins[filt], colormap=:bluesreds)

	Colorbar(fig[1, 2], p, label="distance along orbit / degrees")
	fig
end

# ╔═╡ 204c31a5-4270-4194-bf1a-be079ba04ed4
function xieta_axis()
	return FigAxis(
		xlabel=L"$\xi$ / degrees",
		ylabel=L"$\eta$ / degrees",
		aspect=DataAspect(),
	)
end

# ╔═╡ 9ea84bd5-c381-49e2-99a1-1168f6568994
function plot_region!(df; kwargs...)
	scatter!(df.xi, df.eta; kwargs...)
end

# ╔═╡ 44e4b452-4fc1-4452-a8f0-9d71cdb5346e
let
	fig, ax = xieta_axis()

	rs = rs_test

	println("background, +")
	for i in eachindex(rs)
		r = rs[i]
		df = select_region(all_stars, r * bg_vector, radius=r_circ_cut)
		plot_region!(df, color=r, colorrange=extrema(rs), colormap=:blues, label="⟂, +; r=$r")

		println(size(df, 1))
	end

	println("background, -")
	for i in eachindex(rs)
		r = -rs[i]
		df = select_region(all_stars, r * bg_vector, radius=r_circ_cut)
		plot_region!(df, color=r, 
			colorrange=extrema(-rs), colormap=Reverse(:greens), label="⟂, -; r=$r")

		println(size(df, 1))
	end

	println("along orbit, +")
	for i in eachindex(rs)
		r = rs[i]
		df = select_region(all_stars, r * orbit_vector, radius=r_circ_cut)
		plot_region!(df, color=r, 
			colorrange=extrema(rs), colormap=:reds, label="along orbit, +; r=$r")
		println(size(df, 1))
	end

	println("along orbit, -")
	for i in eachindex(rs)
		r = -rs[i]
		df = select_region(all_stars, r * orbit_vector, radius=r_circ_cut)
		plot_region!(df, color=r, colorrange=(-rs[end], 0), colormap=(:greys),
		label = "along orbit, -; r=$r")
		println(size(df, 1))
	end
	
	lines!(orbit.xi[idx_orbit], orbit.eta[idx_orbit])

	#Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ Cell order:
# ╠═a27c968e-5fde-11ef-27c0-f5f22481829a
# ╠═003e78ba-37b7-4650-ab1e-ba7db2a15e70
# ╠═f0bdd44b-0bcb-41ab-966f-afd514142629
# ╠═0e59c6c9-f053-439f-a6b5-5081e7a7a1a2
# ╠═a5329eeb-9ec4-49ed-9c36-6307a67e3e83
# ╠═fd364bd0-3740-4cb7-941f-10fe47143778
# ╠═506510de-af00-4668-871e-ea00ac6810ef
# ╠═5438d9d9-0256-4779-8c64-56b560c2f730
# ╠═1629ba90-c0ba-477f-9bc4-897bf61426fe
# ╠═784654ff-1596-49eb-8c38-a28200a8bd0c
# ╠═d4e6fe00-e886-4807-8c92-0ce50fb124be
# ╠═c4a62177-c50a-4800-8fbc-d443216c571e
# ╠═f33adaef-228d-483e-938d-0d2ccd828ee8
# ╠═8e5e7421-042c-4876-a992-9f847f1e7ba8
# ╠═d5ef30f4-7a58-41bd-9c42-62689c257bcd
# ╠═5b791110-87d8-4139-abeb-91f6e2688304
# ╠═a3e53d78-40ac-43f2-b5a9-b3b7d5990a1d
# ╠═9438a961-f66f-4be3-885c-74fb0f27661d
# ╠═8ddacf0c-91f1-4543-ac2e-e37871c977f8
# ╠═33ca51be-05cc-4d27-92ba-d7979aa44f7e
# ╠═3b0992aa-1d1a-42b9-9db4-ef32ea8d7a4c
# ╠═5708825c-36c6-4937-b673-21ae1cd0449c
# ╠═340d6d44-2204-42f5-847f-0f9645e4df2c
# ╠═917516ad-0807-4fb1-9afb-d28044249c72
# ╠═c8830c99-f0a2-4f3a-9080-e4ae539565a8
# ╠═06f103a5-05ae-419d-af6f-21149d162af7
# ╠═7b8487e3-b85c-433e-b8e8-47ef399caab9
# ╠═32731e29-c988-437a-9e5f-ade7e6e9ffe9
# ╠═845dd99c-1534-47cb-a84d-42968bfbcb8d
# ╠═d6ca26cb-e4f7-4d17-ba46-0e082fb02d14
# ╠═0ab3e71f-8220-4293-bf74-95be5d443338
# ╠═fb8b12bf-a7f2-4f02-b1c8-edca5f2980ab
# ╠═d01537bf-6a71-47dc-b818-58693aac34fb
# ╠═ec205739-94e5-4ad1-a006-d20f89f160ec
# ╠═dfb0e6da-1871-4dc7-a2b9-b77f08982111
# ╠═28bd7605-91d1-4c8c-a5f4-91af51a2d484
# ╠═b35026bd-8ca4-4552-8d68-2b8811f7b3d5
# ╠═6736916a-3484-4580-a909-6ba7e32e174a
# ╠═4cb80651-8b11-47ed-bf0f-c858b420978c
# ╠═72af6da6-22d1-4bbd-b009-bd801133b271
# ╠═89c96e28-14ba-4492-8f05-5474da8eef15
# ╠═01717a66-0e21-49e1-ba8e-4e44a0826bb2
# ╠═7b1a2c03-9df2-4ec2-8eeb-64c9c9c2cdd4
# ╠═f826de98-01bb-4778-9a17-d2674a6cfae6
# ╠═3bd33e91-89af-4e65-b25e-df33709823fd
# ╠═92f350c8-06bb-4625-a1ea-d27af8dc2864
# ╠═0bea286f-b58c-4e91-8bfd-faf546f7ac10
# ╠═086e9e5b-e5cd-4e4c-a48c-954792754115
# ╠═90b17eba-df5e-465f-a9d9-c8e2ba41b48c
# ╠═60e6c0c9-e382-4751-b0c9-e5062884b3c8
# ╠═ab81099c-5ffc-49d9-b46c-57cb6e6711c3
# ╠═0d56d343-adef-488a-b494-2b9aa3f0bfbf
# ╠═e898391c-2690-4a25-8d44-2b49ddf1a097
# ╠═e1a8309c-c626-4604-b2a6-7d11ef0472a4
# ╠═3d981975-a2f7-45a7-92df-c464e56370b3
# ╠═6cd0161f-657f-446d-a046-f84fdb8d792f
# ╠═7fd40134-633c-4cb9-87cd-2254c485e615
# ╠═204c31a5-4270-4194-bf1a-be079ba04ed4
# ╠═9ea84bd5-c381-49e2-99a1-1168f6568994
# ╟─44e4b452-4fc1-4452-a8f0-9d71cdb5346e
