### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 3a02cce0-3430-11f1-abe3-ef9766dea525
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	
	using PyFITS
end

# ╔═╡ 84d23982-3393-4319-95b8-ded7e66be1c8
using OrderedCollections

# ╔═╡ f912dcd3-51fe-456f-8189-b960dda0cdc8
import TOML

# ╔═╡ 84e4e3a4-63c0-4567-869a-8aa889913a6f
CairoMakie.activate!(type=:png)

# ╔═╡ 4f5c53c0-3552-4ab3-8bc9-86042b738622
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "bootes3")

# ╔═╡ a2daf333-f3ab-45bf-9160-a37d27bcf0b7
function get_time_f(name)
	props = TOML.parsefile(joinpath(model_dir, name, "orbital_properties.toml"))

	props["t_f_gyr"] / T2GYR
end


# ╔═╡ d8e2a511-8209-4a14-9692-0cd90d714d76
function read_scalars(name)
	df= read_fits(joinpath(model_dir, name, "profiles_scalars.fits"))
	df.time .-= get_time_f(name)
	df
end

# ╔═╡ 507c3285-04c0-4ce8-966a-3ad684d921ac
function read_stars(name, starsname)
	df = read_fits(joinpath(model_dir, name, "stars", starsname, "stellar_profiles_3d_scalars.fits"))

	df.time .-= get_time_f(name)
	df

end

# ╔═╡ 7af7cf5a-d149-4fe0-9566-d23b293af3ec
scalars = OrderedDict(
	"compact: 1x1.5kpc" => read_scalars("1e6_v30_r3.0/1_peri_1.5kpc"),
	"compact: 2x7kpc" => read_scalars("1e6_v30_r3.0/2_peri_7kpc"),
	"compact: 5x18kpc" => read_scalars("1e6_v30_r3.0/5_peri_18kpc"),
	"mean: 1x12kpc" => read_scalars("1e6_v22_r3.9/1_peri_12kpc"),
	"mean: 3x26kpc" => read_scalars("1e6_v22_r3.9/3_peri_26kpc"),

)

# ╔═╡ 76789b1e-022a-425c-ab8f-0e42f0e422f2
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = L"r_\text{max}", 
			 ylabel = L"\text{v}_\text{max}",
			 xscale=log10,
			 yscale=log10)

	for (label,df) in scalars
		lines!(df.r_circ_max ./ df.r_circ_max[1], df.v_circ_max / df.v_circ_max[1], label=label)
	end

	x, y = LilGuys.EN21_tidal_track(1, 1, x_min=0.05)
	lines!(x, y, color=:black)
	

	axislegend(position=:rb)

	fig

end

# ╔═╡ 7c7cf8cf-66e3-4a90-ba31-d41c993b2d81
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = L"\text{v}_\text{max}",
			 )

	for (label,df) in scalars
		lines!(df.time * T2GYR, df.v_circ_max*V2KMS, label=label, alpha=1)
	end

	fig

end

# ╔═╡ 63472478-d8be-422d-8c91-ceb2f8a3ed9b
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ b152eb97-6e1a-4886-87e7-88e80663051f
σv = obs_props["sigma_v"]

# ╔═╡ 693cc2b8-02b1-4830-8dea-303d195f1e8c
σv_range = obs_props["sigma_v"] - obs_props["sigma_v_em"], obs_props["sigma_v"] + obs_props["sigma_v_ep"]

# ╔═╡ a27bea66-cdcb-4aee-ad48-a309351c0a6e
R_h = LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])

# ╔═╡ a8be9a3e-6f5f-4371-88ec-ce0306ec9975
α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())

# ╔═╡ b5d50f54-1bfb-49e6-9c3d-473a36630586
prof_stars = LilGuys.Exp2D(R_s=R_h / LilGuys.R_h(LilGuys.Exp2D()))

# ╔═╡ 2794a7c3-cbf8-4333-82e9-eb74ff9380fe
function calc_σv(r_circ_max, v_circ_max)
	h = LilGuys.ExpCusp(r_s=r_circ_max / α_exp_cusp, M=1)
	v_scale = v_circ_max / LilGuys.v_circ_max(h) 
	h = LilGuys.ExpCusp(r_s = r_circ_max/α_exp_cusp, M=v_scale^2)

	return LilGuys.σv_star_mean(h, prof_stars)
end

# ╔═╡ 6cbce1fb-2abe-440e-b2b5-4defdfc6e9ab
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = L"\sigma_\text{v}",
			  limits=(-10, 1, 5, 14)
			 )


	for (label,df) in scalars
		filt = .!ismissing.(df.r_circ_max)
		σv = calc_σv.(df.r_circ_max[filt], df.v_circ_max[filt])
		lines!(df.time[filt] * T2GYR, σv*V2KMS, label=label, alpha=1)
	end

	hlines!(σv, color=:black)
	hspan!(σv_range..., alpha=0.2, color=:black)

	Legend(fig[1, 2], ax)
	fig

end

# ╔═╡ 39104ba2-0b62-46e3-b8e1-f56413d9c020
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = L"\sigma_\text{v}",
			  limits=(-0.5, 0.5, 5, 14)
			 )


	for (label,df) in scalars
		filt = .!ismissing.(df.r_circ_max)
		σv = calc_σv.(df.r_circ_max[filt], df.v_circ_max[filt])
		lines!(df.time[filt] * T2GYR, σv*V2KMS, label=label, alpha=1)
	end

	hlines!(σv, color=:black)
	hspan!(σv_range..., alpha=0.2, color=:black)

	vlines!(0)
	Legend(fig[1, 2], ax)
	fig

end

# ╔═╡ 5103affe-1820-4d4e-82a9-c540011f7175
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 03287f9d-158c-4369-8ddc-8285bbb34d5a
md"""
# Stellar tracks
"""

# ╔═╡ 69ad681c-bf5a-436c-b1a4-7148ec996666
stars_tracks = OrderedDict(
	"compact: 1x1.5kpc" => read_stars("1e6_v30_r3.0/1_peri_1.5kpc", "exp2d_rs0.20"),
	"compact: 2x7kpc" => read_stars("1e6_v30_r3.0/2_peri_7kpc", "exp2d_rs0.20"),
	"compact: 5x18kpc" => read_stars("1e6_v30_r3.0/5_peri_18kpc", "exp2d_rs0.20"),
	"mean: 1x12kpc" => read_stars("1e6_v22_r3.9/1_peri_12kpc", "exp2d_rs0.20"),
	"mean: 3x26kpc" => read_stars("1e6_v22_r3.9/3_peri_26kpc", "exp2d_rs0.20"),

)

# ╔═╡ 8ef648fb-9e38-4892-b0ba-401a1d7bc86b
function plot_stars(modelname)
	dm_track = scalars[modelname]
	stars_track = stars_tracks[modelname]

	fig = plot_stars(dm_track[dm_track.time .< 0, :], 
					 stars_track[stars_track.time .< 0, :]
					)
	fig.content[1].title = modelname
	fig
end

# ╔═╡ 58c4be38-b1ea-4d4b-9640-fbeb47943ce3
function plot_stellar_tracks(stars_tracks; xlims=(nothing, nothing))
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1], xlabel = "time / Gyr", ylabel = "bound mass fraction")

	for (label, track) in stars_tracks
		lines!(track.time * T2GYR, track.bound_mass, label=label)
	end
	vlines!(0, linestyle=:dash, color=:black)
	axislegend(position=:lb)
	xlims!(xlims...)


	ax = Axis(fig[1,2], xlabel = "time / Gyr", ylabel = L"$\sigma_\text{v}$ / km\,s$^{-1}$")

	for (label, track) in stars_tracks
		lines!(track.time * T2GYR, track.sigma_v * V2KMS)
	end
	vlines!(0, linestyle=:dash, color=:black)
	hlines!(σv, color=:black)
	hspan!(σv_range..., alpha=0.2, color=:black)
	xlims!(xlims...)

	fig
end

# ╔═╡ e9e9acf2-3f3d-4ce5-931e-b9ff78982b8f
plot_stellar_tracks(stars_tracks)

# ╔═╡ 66b9d281-3182-45d7-92c9-5e5757128772
plot_stellar_tracks(stars_tracks, xlims=(-0.5, 0.5))

# ╔═╡ 19fdf1e4-9131-4ece-9271-231f056b541d
function plot_stars(dm_track, stars_track)

	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale = log10, 
			 yscale = log10,
			 xticks = [0.01, 0.1, 1, 10, 100], 
			 yticks = collect(0.2:0.2:1),
			  yminorticks = 0.1:0.1:1,
			  xlabel = L"r_\text{max} / r_\text{max, 0}",
			  ylabel = L"\text{v}_\text{max} / \text{v}_\text{max, 0}",
			 )
	
	h_i = LilGuys.TruncNFW(r_circ_max=1, v_circ_max=1, trunc=20, xi=3)
	r = logrange(0.01, 10, 100)
	v = LilGuys.v_circ.(h_i, r)
	lines!(r, v, color=:black)


	r_scale  = dm_track.r_circ_max[1]
	v_scale = dm_track.v_circ_max[1]


	scatterlines!(dm_track.r_circ_max ./ r_scale, dm_track.v_circ_max ./ v_scale, label="DM")
	scatterlines!(stars_track.r_h ./ r_scale, stars_track.sigma_v .* sqrt(3) ./ v_scale, label="stars")

	r, v = LilGuys.EN21_tidal_track(1, 1, x_min=0.03)

	lines!(r, v, color=:black, linestyle=:dot)

	axislegend(position=:lt)
	
	fig
end

# ╔═╡ bfdc875a-43e5-4336-a31c-835399c2707b
plot_stars("mean: 1x12kpc")

# ╔═╡ 7bb69ca0-e379-4699-888a-48b363d877f7
plot_stars("compact: 1x1.5kpc")

# ╔═╡ aa39b874-1fdf-42c4-b615-3de554296fc4
plot_stars("compact: 2x7kpc")

# ╔═╡ 9c711002-3c5f-4943-a6b8-88d4728c3258
plot_stars("mean: 3x26kpc")

# ╔═╡ b3deb14d-83e8-481f-a35f-ded18fdcd8a7
plot_stars("compact: 5x18kpc")

# ╔═╡ 6e1f98e2-a53a-448b-bcc1-8a3a73d94a73
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = L"r_\text{max}", 
			 ylabel = L"\text{v}_\text{max}",
			 xscale=log10,
			 yscale=log10)

	for (label,df) in scalars
		lines!(df.r_circ_max ./ df.r_circ_max[1], df.v_circ_max / df.v_circ_max[1], label=label)
	end

	x, y = LilGuys.EN21_tidal_track(1, 1, x_min=0.05)
	lines!(x, y, color=:black)
	

	axislegend(position=:rb)

	fig

end

# ╔═╡ Cell order:
# ╠═3a02cce0-3430-11f1-abe3-ef9766dea525
# ╠═f912dcd3-51fe-456f-8189-b960dda0cdc8
# ╠═84d23982-3393-4319-95b8-ded7e66be1c8
# ╠═84e4e3a4-63c0-4567-869a-8aa889913a6f
# ╠═4f5c53c0-3552-4ab3-8bc9-86042b738622
# ╠═a2daf333-f3ab-45bf-9160-a37d27bcf0b7
# ╠═d8e2a511-8209-4a14-9692-0cd90d714d76
# ╠═507c3285-04c0-4ce8-966a-3ad684d921ac
# ╠═7af7cf5a-d149-4fe0-9566-d23b293af3ec
# ╠═76789b1e-022a-425c-ab8f-0e42f0e422f2
# ╠═7c7cf8cf-66e3-4a90-ba31-d41c993b2d81
# ╠═2794a7c3-cbf8-4333-82e9-eb74ff9380fe
# ╠═6cbce1fb-2abe-440e-b2b5-4defdfc6e9ab
# ╠═39104ba2-0b62-46e3-b8e1-f56413d9c020
# ╠═b152eb97-6e1a-4886-87e7-88e80663051f
# ╠═693cc2b8-02b1-4830-8dea-303d195f1e8c
# ╠═63472478-d8be-422d-8c91-ceb2f8a3ed9b
# ╠═a27bea66-cdcb-4aee-ad48-a309351c0a6e
# ╠═a8be9a3e-6f5f-4371-88ec-ce0306ec9975
# ╠═b5d50f54-1bfb-49e6-9c3d-473a36630586
# ╠═5103affe-1820-4d4e-82a9-c540011f7175
# ╟─03287f9d-158c-4369-8ddc-8285bbb34d5a
# ╠═69ad681c-bf5a-436c-b1a4-7148ec996666
# ╠═8ef648fb-9e38-4892-b0ba-401a1d7bc86b
# ╠═58c4be38-b1ea-4d4b-9640-fbeb47943ce3
# ╠═e9e9acf2-3f3d-4ce5-931e-b9ff78982b8f
# ╠═66b9d281-3182-45d7-92c9-5e5757128772
# ╠═19fdf1e4-9131-4ece-9271-231f056b541d
# ╠═bfdc875a-43e5-4336-a31c-835399c2707b
# ╠═7bb69ca0-e379-4699-888a-48b363d877f7
# ╠═aa39b874-1fdf-42c4-b615-3de554296fc4
# ╠═9c711002-3c5f-4943-a6b8-88d4728c3258
# ╠═b3deb14d-83e8-481f-a35f-ded18fdcd8a7
# ╠═6e1f98e2-a53a-448b-bcc1-8a3a73d94a73
