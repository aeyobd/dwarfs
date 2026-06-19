### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 998add30-6048-11f1-b3d8-71f16959ead9
begin
	import Pkg; Pkg.activate()

	using FITSIO
	using LilGuys
	using CairoMakie, Arya
end

# ╔═╡ c26bb9fb-4136-4e03-85da-065980bd5db5
using Distributions

# ╔═╡ 75d2c5d3-518c-433b-ad8f-74f9e6b9cb12
import TOML

# ╔═╡ b41ecde4-9ff8-4be9-a0d3-1a46bb2ca5ec
md"""
# matched filter loading
"""

# ╔═╡ 25031e14-432b-4c75-8a10-5b6a4ecd64cb
file = "DELVE/matchedfilter_density_G09_v2.fits"

# ╔═╡ c9285d36-093c-4b87-9fc8-5771a526208d
f = FITS(file, "r")

# ╔═╡ 696a39cb-f5e2-42ba-887e-e315ada9e154
mfm = read(f[1])

# ╔═╡ 47627af3-35ca-4e21-801a-fa930456eab5
read_header(f[1])

# ╔═╡ 0d2886e9-9244-44fd-b256-4b960e089166
crval1 = 22.820584885662086

# ╔═╡ dad0b1d1-f395-4f3b-b0c4-77d9e29ab2b3
crdelt1 = -0.1

# ╔═╡ dc5dff4c-40f2-4677-ac4a-e33a7c630dcf
crval2 = -20.62692894083525

# ╔═╡ 8c67d8ae-b478-4457-8195-1ff9824835fd
crdelt2 = 0.1

# ╔═╡ 63a3de10-42f5-43fe-a87d-f46b02287a84
xi = crval1 .+ crdelt1*(0:size(mfm, 1) .- 0.5)

# ╔═╡ 97b39683-045d-4a7b-b302-09fef5b49187
eta = crval2 .+ crdelt2*(0:size(mfm, 2) .- 0.5)

# ╔═╡ e544c971-18d7-4251-af12-f0a1752d0582
md"""
## sanity checking coordinate system
"""

# ╔═╡ 213c04fc-847a-4b21-95cb-43275362ad6d
import PythonCall

# ╔═╡ 15b18332-4bec-45c4-a520-29916718eadf
nddata = PythonCall.pyimport("astropy.nddata")

# ╔═╡ c9d5a1d5-d2c3-48f1-887f-709b4e119ee6
py_f = nddata.CCDData.read(file, unit="adu")

# ╔═╡ 57b162ac-e72b-45a7-9142-73c15597a8fc
wcs = py_f.wcs

# ╔═╡ 22d5787e-d7e3-4c5b-81dd-41215fcef29d
wcs.pixel_to_world(5, 9), xi[5], eta[9]

# ╔═╡ d3c5d1c5-a495-4807-8b4e-99a06cef5990
md"""
# selection
"""

# ╔═╡ 14775e8a-7315-408d-bc23-75febc0c08b4
let
	fig = heatmap(xi, eta, mfm, colorscale=log10, colorrange=(0.1, 3), axis=(; 
																  xreversed=true, aspect=DataAspect()))

	arc!((0,0), 10, 0, 2π)

	arc!((-20, -10), 10, 0, 2π)

	lines!([18, -40, -40, 18, 18], [-18, -18, 7, 7, -18])

	fig
end

# ╔═╡ eaaeeee1-48c7-45a1-9320-417baa4a20f6
x = [x for x in xi, y in eta]

# ╔═╡ a226754c-84f5-40a9-8129-74b7288b2e58
y = [y for x in xi, y in eta]

# ╔═╡ f6fc8c2c-60fc-43bd-93c0-04d67489e532
mfm

# ╔═╡ 04c40bb3-a3b8-44c2-aeda-92b8b2aed5d6
let
	f = heatmap(xi, eta, mfm, colorscale=log10, colorrange=(0.1, 3), axis=(; 
																  xreversed=true))
	xlims!(-10, 10)
	ylims!(-10, 10)


	arc!((0,0), 10, 0, 2π)
	f
end

# ╔═╡ 6093b60e-9f01-4e73-9479-d121ccabd4d0
sum(filter(isfinite, mfm))

# ╔═╡ eb4b561c-9068-4a05-814e-bd0dda79937e
let
	fig = Figure()

	ax = Axis(fig[1,1])

	
	filt = @. sqrt(x^2 + y^2) < 10
	
	stephist!(filter(isfinite, mfm[filt]), bins=LinRange(-0.1, 3, 300), normalization=:pdf)


	
	filt = @. sqrt((x .- 10)^2 + (y .- 10)^2) < 10
	
	stephist!(filter(isfinite, mfm[filt]), bins=LinRange(-0.1, 3, 100), normalization=:pdf)



	# stephist!(vec(mfm[-20 .< xi .< 18, -7 .< eta .< 0]), bins=LinRange(-0.1, 3, 100), normalization=:pdf)

	plot!(Normal(0.55, 0.14))
	xlims!(-0.1, 3)
	fig
end

# ╔═╡ cee7ea09-a507-4b17-aeab-02dde6f4f599
let
	filt = @. sqrt(x^2 + y^2) < 1.5
	sum(mfm[filt]) .- sum(filt) * 0.55
end

# ╔═╡ 266a33ca-587a-46f9-8142-0f10139b52e3
λ = median(filter(isfinite, mfm))

# ╔═╡ 856e97b9-a1b6-4305-a47f-ad928294e5d4
md"""
# MC sampling
"""

# ╔═╡ a8e650a9-dada-4259-9d22-90898584964b
box_length = 3

# ╔═╡ 3f01a40d-53f9-42a5-a9b5-0c6e7e5d1935
box_width = 0.5

# ╔═╡ 4a3ed690-fbfa-45ec-b13a-bdfabc769eab
function draw_box(; xrange=(-20, 15), yrange=(-7, -1), width=box_width, length=box_length)
	xlims = xrange[1] + length/2, xrange[2] - length/2
	ylims = yrange[1] + length/2, yrange[2] - length/2

	x0 = rand(Uniform(xlims[1], xlims[2]))
	y0 = rand(Uniform(ylims[1], ylims[2]))
	θ = rand(Uniform(0, 2π))

	dx = x .- x0
	dy = y .- y0
	x_p = @. dx * cos(θ) + dy*sin(θ)
	y_p = @. -dx * sin(θ) + dy*cos(θ)

	filt = @. (abs(x_p) < length/2 ) && (abs(y_p) < width/2)
end	

# ╔═╡ ae58ad0f-1ee1-49a5-85ab-2616738915ec
sum(draw_box())

# ╔═╡ b5a2f637-9a95-45ce-af8b-059c6e3966df
heatmap(xi, eta, draw_box(), 
	   axis=(; aspect=DataAspect()))

# ╔═╡ 5e9a972f-3e7f-4eea-a076-1da41c83a4b1
let
	box = Int.(draw_box())

	for i in 1:100
		box .+= draw_box()
	end

	box

	heatmap(box)
end

# ╔═╡ c0e7a3fa-cd36-4360-b593-15d8311d36cc
heatmap(mfm)

# ╔═╡ 83092123-5480-4549-8734-072b4cf0dc93
let
	f = heatmap(xi, eta, mfm, colorscale=log10, colorrange=(0.1, 3), axis=(; 
																  xreversed=true))
	xlims!(-20, 18)
	ylims!(-7, -0)


	f
end

# ╔═╡ e2a81ffe-dc4a-488c-9f9f-33bdd436332b
heatmap(mfm[-20 .< xi .< 15, -7 .< eta .< -1], colorscale=log10)

# ╔═╡ dd03b16f-d736-443d-8573-e9019f3c647c
N_box_samples = [sum(filter(isfinite, mfm[draw_box()])) for _ in 1:1000]

# ╔═╡ b9d2e4ab-37a7-491b-9170-0b65b515881e
std(N_box_samples), quantile(N_box_samples, 0.84) - median(N_box_samples)

# ╔═╡ 366e6995-a0b1-4b98-b623-294cb5dcdda3
σ_px = std(mfm[-20 .< xi .< 18, -7 .< eta .< 0])

# ╔═╡ 83e7578a-d400-4932-8ac0-f602bf0c29e9
sum(draw_box())

# ╔═╡ b93b7f88-779a-41c2-b471-2952553e3bda
Npx = box_width * box_length * 1/ abs(crdelt1*crdelt2)

# ╔═╡ b7260462-4e00-4040-90e4-e619377fc1ab
σ_px * sqrt(Npx)

# ╔═╡ 9defcef2-0748-4196-a3f5-e3decfadebf5
md"""
# Calculation of surface density
"""

# ╔═╡ 604abdbc-782b-483c-90bb-dc45c976a58e
box_area = box_length * box_width * 60^2

# ╔═╡ 909c65f5-1d1a-49c1-bd69-be492a5c3319
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = L"$\Sigma_\star - \bar\Sigma_\star$ / arcmin$^{-2}$ (random box)", ylabel="density")

	μ = median(N_box_samples / box_area) 

	hist!(N_box_samples / box_area .- μ, bins=50, normalization=:pdf)

	# = quantile(samples, 0.84) - median(samples)
	σ = std(N_box_samples / box_area)
	plot!(Normal(0, σ), color=COLORS[2])

	fig
end

# ╔═╡ 71306426-d6e6-42f1-85c2-4ce290f4a545
Σ_bkg = std(N_box_samples) / box_area

# ╔═╡ 5df10ad0-90b0-4410-8268-0f0717d9e28d
Σ_bkg_quantile = (quantile(N_box_samples, 0.84) - median(N_box_samples)) / box_area

# ╔═╡ a81801ac-ffa2-4bf9-a0a8-a032bdd3ffbf
Mstar_to_N = 8 # need to triple check this??

# ╔═╡ 4b7a9a3b-5179-4ffd-9975-588f6f79e22e
std(N_box_samples) * Mstar_to_N

# ╔═╡ d97bb93b-82f2-4411-b28b-5a614480547e
std(N_box_samples) * Mstar_to_N / box_area

# ╔═╡ 3520dbc5-94c6-43f4-8adb-c6638f682f31
Σ_bkg_quantile * Mstar_to_N

# ╔═╡ 4801e481-2c5a-4c7e-b9d8-f34daab335ed
md"""
# Proper calculation of surface brightness
"""

# ╔═╡ a0319c0e-d7a5-4231-8d8e-2507958ec51c
md"""
## Luminosity samplers
"""

# ╔═╡ 8e9058c4-fcb7-48f6-a067-b38c3a32050f
module Utils # only used for obs_dir and kroupa_imf
	obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")
	include(joinpath(obs_dir, "delve_utils.jl"))
end

# ╔═╡ 33025d08-6a6d-4c94-a888-a75029597eab
obs_dir = Utils.obs_dir

# ╔═╡ 6ac31fad-a5b8-4040-a8d6-0e13ef4aa584
read_mist_file = Utils.read_mist_file

# ╔═╡ 14ea35fc-b573-4e9e-85f5-d3dd50b5840d
kroupa_imf = Utils.kroupa_imf

# ╔═╡ 6b2a2036-ddbf-4516-a840-e00c2c7f33ec
iso_age = 12

# ╔═╡ a2815479-3a58-465d-acc8-eb183ba7663d
isochrones_delve = read_mist_file("isochrone.decam.$(iso_age)Gyrs.dat")

# ╔═╡ 417659b1-175d-404a-a754-e22f227d0d05
function get_isochrone(isochrones, M_H; stage_max=5)

	isos = isochrones

	if (M_H < minimum(isos.MH)) || (M_H > maximum(isos.MH))
		throw(DomainError(M_H, "metallicity out of isochrone range: $(extrema(isos.MH))"))
	end

	M_Hs = unique(isos.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]

	filt = isapprox.(isos.MH, M_H_adopted)
	filt .&= isos.label .<= stage_max # only keep through HB

	return isos[filt, :]
end

# ╔═╡ 313dc43b-c831-4fd3-ad3f-d9da0d265d59
obs_props = TOML.parsefile(joinpath(Utils.obs_dir, "observed_properties.toml"))

# ╔═╡ ab9a7df3-eb9c-4e17-9247-7f045123ffb7
dm_0 = obs_props["distance_modulus"]

# ╔═╡ 4db095ec-b56e-4e56-9787-02986c5a28bd
function find_tot_number_stars(N_obs; mag_limit, N_max=300_000, distance_modulus=dm_0, mag_interp, mass_sampler)

	count = 0
	mass_tot = 0
	L_tot = 0
	N_tot = 0
	
	for i in 1:N_max
		mass = mass_sampler()
		mag = mag_interp(mass) + distance_modulus
		count += mag < mag_limit

		mag_v = mag #mag_v_interp(mass) + distance_modulus 

		N_tot += 1
		L_tot += 10^(-0.4*mag_v)
		mass_tot += mass

		if count >= N_obs
			break
		end
	end

	if N_tot == N_max
		@warn "failed to converge"
		@info count / N_obs
		return NaN, NaN, NaN
	end

	return N_tot, mass_tot, -2.5*log10(L_tot) - distance_modulus

end

# ╔═╡ 329f398a-aa97-4cbf-a9f3-e7fe5c445cfd
dm_err = obs_props["distance_modulus_em"]

# ╔═╡ a7202e2a-a70f-412b-961d-cd18e6d42795
isochrones_delve

# ╔═╡ b60fd709-806f-4b20-9df0-bf2165a7b414
function interpolate_magnitude(iso, mag_col)
	return LilGuys.lerp(iso.Mini, iso[!, mag_col])
end

# ╔═╡ 57e99986-e81f-42b9-9954-b06e7cc8943c
let
	fig = Figure()
	ax = Axis(fig[1,1])

	iso = get_isochrone(isochrones_delve, -2.19, )
	gmag = interpolate_magnitude(iso, "gmag")
	rmag = interpolate_magnitude(iso, "rmag")


	M = LinRange(0.5, maximum(iso.Mini), 1000)
	mag = gmag.(M)
	colour = mag .- rmag.(M)

	scatter!(colour, mag .+ dm_0)
	ylims!(17, 24)

	ax.yreversed[] = true
	ax.aspect[] = 0.5
	fig
end

# ╔═╡ b440a540-9466-4b85-8713-3c0d43f4d185
function sample_tot_stars(obs_samples, N_samples; 
						  M_H, all_isochrones=isochrones_delve,
						  gcol = "gmag",
						  kwargs...)
	N_tots = Vector{Int}(undef, N_samples)
	MVs = Vector{Float64}(undef, N_samples)
	M_tots = Vector{Float64}(undef, N_samples)


	iso = get_isochrone(all_isochrones, M_H, )
	mag_interp = interpolate_magnitude(iso, gcol)
	# iso_v = get_isochrone(all_isochrones_ubvri, M_H)
	# mag_v_interp = interpolate_magnitude(iso_v, "Vmag")

	
	M_max = maximum(iso.Mini)
	mass_sampler = Utils.create_kroupa_sampler(M_max, mass_min=0.08)

	
	for i in 1:N_samples
		distance_modulus = dm_0 + dm_err*randn()

		N_tots[i], M_tots[i], MVs[i] = find_tot_number_stars(obs_samples[i]; 
			distance_modulus=distance_modulus,
			mass_sampler=mass_sampler, 
			mag_interp=mag_interp,
			kwargs...)
		
		if (i % 100) == 0
			@info "completed $i / $N_samples"
		end
	end

	return N_tots, M_tots, MVs
end

# ╔═╡ ea9d0046-df60-4d9f-a6cb-0d11a43aa2f1
N_tots, M_tots, M_Vsmples = sample_tot_stars(N_box_samples .- median(N_box_samples), 1000, M_H=-2.19, mag_limit=24)

# ╔═╡ 43f4a295-f5e0-48fd-af52-723f0b13a977
hist(M_tots ./ box_area, bins=50)

# ╔═╡ 69497e7d-cbea-4482-b822-f80e7aac5ee4
hist(N_tots ./ box_area, bins=50)

# ╔═╡ 8dadb348-8f96-4294-8a11-efeff41009e4
quantile(M_tots ./ N_box_samples, 0.5)

# ╔═╡ 37ea4296-3dce-4588-9ad3-76746af11d35
quantile(M_tots ./ box_area, 0.84) - quantile(M_tots ./box_area, 0.5)

# ╔═╡ 88068f3f-dbba-4e34-8c37-54ef916710a1
(quantile(M_tots ./ box_area, 0.84) - quantile(M_tots ./box_area, 0.5)) * 60^2

# ╔═╡ 5bfc4499-4fb8-4a44-99db-5504ea253ba8
quantile(M_tots ./box_area, 0.5)

# ╔═╡ 3189334a-55ef-404b-b07f-c339f19ac8a9
quantile(N_tots ./ box_area, 0.84) - quantile(N_tots ./box_area, 0.5)

# ╔═╡ 9ad36a57-f389-43d6-94e6-776f4d05d499
Mstar_to_N

# ╔═╡ 32b8533a-2d92-4cd8-b8bb-c17b7c8ccb8c
md"""
## Integration
"""

# ╔═╡ f8a23af2-47d1-4197-8071-fca36c0524f5
let
	iso = get_isochrone(isochrones_delve, -2.19)

	gmag = iso.gmag .+ dm_0 .+ 0dm_err

	M = iso.Mini
	w = kroupa_imf.(iso.Mini) .* LilGuys.gradient(iso.Mini)
	f = lines(LilGuys.gradient(iso.Mini))
	lines!(diff(iso.Mini))
	f
end

# ╔═╡ 09b5d2f7-b1a2-40f9-8433-2867cc861c4d
Mstar_to_N

# ╔═╡ 3e8ba495-9d22-4d9e-984c-d809382e755d
let
	iso = get_isochrone(isochrones_delve, -2.19)

	gmag = iso.gmag .+ dm_0 .+ 0dm_err

	M = iso.Mini
	w = kroupa_imf.(iso.Mini) .* LilGuys.gradient(iso.Mini)
	@info LilGuys.gradient(iso.Mini)
	@info diff(iso.Mini)

	M_tot = sum(iso.Mass .* w) ./ sum(w)
	N_obs = sum(w[gmag .< 24]) / sum(w)

	M_tot, N_obs
	M_tot / N_obs
end

# ╔═╡ 2d48f5a8-0e9a-4192-bb49-5d0ea92b121e
md"""
# Reverse sampling
"""

# ╔═╡ 66248a8b-ad2f-43ca-ad65-b21306f58169
quantile(M_tots, 0.84)

# ╔═╡ ebfd219f-d7bd-4a85-aaee-8fd25ecbe347


# ╔═╡ 9ef7d0f6-ddf6-4a96-9eb7-cc8e46ca0b52
let 
	M_tot = quantile(M_tots, 0.84)
	stars = sample_tot_stars


	iso = get_isochrone(isochrones_delve, -2.19)

	gmag = iso.gmag .+ dm_0 .+ 0dm_err

	M_max = maximum(iso.Mini)
	mass_sampler = Utils.create_kroupa_sampler(M_max, mass_min=0.08)

	mag_interp = interpolate_magnitude(iso, "gmag")
	M = 0
	N = 0
	m_min = 100
	Ntot = 0
	while M < M_tot
		m = mass_sampler()
		if mag_interp(m) .< 24 .- dm_0
			N += 1
			m_min = min(m, m_min)
		end
		M += m
		Ntot += 1
	end

	@info "m_min: $m_min, $M_max"

	@info M_tot ./ N
	@info Ntot ./ N
	xs = rand(Uniform(-1.5, 1.5), N)
	ys = rand(Uniform(-0.25, 0.25), N)

	@info N ./ box_area
	@info M ./ box_area
	@info N
	scatter(xs, ys)
end

# ╔═╡ 6472eee5-08f8-4d84-adc4-3488e6497bdd
1 / (1 - LilGuys.integrate(kroupa_imf, 0.08, 0.7) ./ LilGuys.integrate(kroupa_imf, 0.08, 0.8))

# ╔═╡ 45f6ae62-379e-4fb1-b5fc-d7f5b61598a7
(kroupa_imf(0.8) .- kroupa_imf(0.7))

# ╔═╡ d2fd1f7e-6e2e-4c93-af7c-19e62d090408
Σ_bkg_quantile

# ╔═╡ Cell order:
# ╠═998add30-6048-11f1-b3d8-71f16959ead9
# ╠═75d2c5d3-518c-433b-ad8f-74f9e6b9cb12
# ╟─b41ecde4-9ff8-4be9-a0d3-1a46bb2ca5ec
# ╠═25031e14-432b-4c75-8a10-5b6a4ecd64cb
# ╠═c9285d36-093c-4b87-9fc8-5771a526208d
# ╠═696a39cb-f5e2-42ba-887e-e315ada9e154
# ╠═47627af3-35ca-4e21-801a-fa930456eab5
# ╠═0d2886e9-9244-44fd-b256-4b960e089166
# ╠═dad0b1d1-f395-4f3b-b0c4-77d9e29ab2b3
# ╠═dc5dff4c-40f2-4677-ac4a-e33a7c630dcf
# ╠═8c67d8ae-b478-4457-8195-1ff9824835fd
# ╠═63a3de10-42f5-43fe-a87d-f46b02287a84
# ╠═97b39683-045d-4a7b-b302-09fef5b49187
# ╠═e544c971-18d7-4251-af12-f0a1752d0582
# ╠═213c04fc-847a-4b21-95cb-43275362ad6d
# ╠═15b18332-4bec-45c4-a520-29916718eadf
# ╠═c9d5a1d5-d2c3-48f1-887f-709b4e119ee6
# ╠═57b162ac-e72b-45a7-9142-73c15597a8fc
# ╠═22d5787e-d7e3-4c5b-81dd-41215fcef29d
# ╠═d3c5d1c5-a495-4807-8b4e-99a06cef5990
# ╠═14775e8a-7315-408d-bc23-75febc0c08b4
# ╠═eaaeeee1-48c7-45a1-9320-417baa4a20f6
# ╠═a226754c-84f5-40a9-8129-74b7288b2e58
# ╠═f6fc8c2c-60fc-43bd-93c0-04d67489e532
# ╠═04c40bb3-a3b8-44c2-aeda-92b8b2aed5d6
# ╠═6093b60e-9f01-4e73-9479-d121ccabd4d0
# ╠═eb4b561c-9068-4a05-814e-bd0dda79937e
# ╠═cee7ea09-a507-4b17-aeab-02dde6f4f599
# ╠═266a33ca-587a-46f9-8142-0f10139b52e3
# ╠═c26bb9fb-4136-4e03-85da-065980bd5db5
# ╟─856e97b9-a1b6-4305-a47f-ad928294e5d4
# ╠═a8e650a9-dada-4259-9d22-90898584964b
# ╠═3f01a40d-53f9-42a5-a9b5-0c6e7e5d1935
# ╠═4a3ed690-fbfa-45ec-b13a-bdfabc769eab
# ╠═ae58ad0f-1ee1-49a5-85ab-2616738915ec
# ╠═b5a2f637-9a95-45ce-af8b-059c6e3966df
# ╠═5e9a972f-3e7f-4eea-a076-1da41c83a4b1
# ╠═c0e7a3fa-cd36-4360-b593-15d8311d36cc
# ╠═83092123-5480-4549-8734-072b4cf0dc93
# ╠═e2a81ffe-dc4a-488c-9f9f-33bdd436332b
# ╠═dd03b16f-d736-443d-8573-e9019f3c647c
# ╠═909c65f5-1d1a-49c1-bd69-be492a5c3319
# ╠═b9d2e4ab-37a7-491b-9170-0b65b515881e
# ╠═366e6995-a0b1-4b98-b623-294cb5dcdda3
# ╠═83e7578a-d400-4932-8ac0-f602bf0c29e9
# ╠═b93b7f88-779a-41c2-b471-2952553e3bda
# ╠═b7260462-4e00-4040-90e4-e619377fc1ab
# ╟─9defcef2-0748-4196-a3f5-e3decfadebf5
# ╠═604abdbc-782b-483c-90bb-dc45c976a58e
# ╠═71306426-d6e6-42f1-85c2-4ce290f4a545
# ╠═5df10ad0-90b0-4410-8268-0f0717d9e28d
# ╠═a81801ac-ffa2-4bf9-a0a8-a032bdd3ffbf
# ╠═4b7a9a3b-5179-4ffd-9975-588f6f79e22e
# ╠═d97bb93b-82f2-4411-b28b-5a614480547e
# ╠═3520dbc5-94c6-43f4-8adb-c6638f682f31
# ╟─4801e481-2c5a-4c7e-b9d8-f34daab335ed
# ╟─a0319c0e-d7a5-4231-8d8e-2507958ec51c
# ╠═8e9058c4-fcb7-48f6-a067-b38c3a32050f
# ╠═33025d08-6a6d-4c94-a888-a75029597eab
# ╠═6ac31fad-a5b8-4040-a8d6-0e13ef4aa584
# ╠═14ea35fc-b573-4e9e-85f5-d3dd50b5840d
# ╠═6b2a2036-ddbf-4516-a840-e00c2c7f33ec
# ╠═a2815479-3a58-465d-acc8-eb183ba7663d
# ╠═417659b1-175d-404a-a754-e22f227d0d05
# ╠═4db095ec-b56e-4e56-9787-02986c5a28bd
# ╠═313dc43b-c831-4fd3-ad3f-d9da0d265d59
# ╠═ab9a7df3-eb9c-4e17-9247-7f045123ffb7
# ╠═329f398a-aa97-4cbf-a9f3-e7fe5c445cfd
# ╠═a7202e2a-a70f-412b-961d-cd18e6d42795
# ╠═b60fd709-806f-4b20-9df0-bf2165a7b414
# ╠═57e99986-e81f-42b9-9954-b06e7cc8943c
# ╠═b440a540-9466-4b85-8713-3c0d43f4d185
# ╠═ea9d0046-df60-4d9f-a6cb-0d11a43aa2f1
# ╠═43f4a295-f5e0-48fd-af52-723f0b13a977
# ╠═69497e7d-cbea-4482-b822-f80e7aac5ee4
# ╠═8dadb348-8f96-4294-8a11-efeff41009e4
# ╠═37ea4296-3dce-4588-9ad3-76746af11d35
# ╠═88068f3f-dbba-4e34-8c37-54ef916710a1
# ╠═5bfc4499-4fb8-4a44-99db-5504ea253ba8
# ╠═3189334a-55ef-404b-b07f-c339f19ac8a9
# ╠═9ad36a57-f389-43d6-94e6-776f4d05d499
# ╠═32b8533a-2d92-4cd8-b8bb-c17b7c8ccb8c
# ╠═f8a23af2-47d1-4197-8071-fca36c0524f5
# ╠═09b5d2f7-b1a2-40f9-8433-2867cc861c4d
# ╠═3e8ba495-9d22-4d9e-984c-d809382e755d
# ╟─2d48f5a8-0e9a-4192-bb49-5d0ea92b121e
# ╠═66248a8b-ad2f-43ca-ad65-b21306f58169
# ╠═ebfd219f-d7bd-4a85-aaee-8fd25ecbe347
# ╠═9ef7d0f6-ddf6-4a96-9eb7-cc8e46ca0b52
# ╠═6472eee5-08f8-4d84-adc4-3488e6497bdd
# ╠═45f6ae62-379e-4fb1-b5fc-d7f5b61598a7
# ╠═d2fd1f7e-6e2e-4c93-af7c-19e62d090408
