### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 43307c80-298c-11ef-3cfe-afb7123432ca
begin
	import Pkg; Pkg.activate()

	import LilGuys as lguys
	using GLMakie

	using Arya
end

# ╔═╡ 624a165e-73f2-40a5-86e7-852a4ec405ca
using FITSIO

# ╔═╡ 1c55d5eb-1065-4ef0-bb97-ada81c195432
N = 10000

# ╔═╡ 1093880a-8973-4aa9-bed2-25376166d705
function sample_radii(profile, N=10000; radii = 10 .^ LinRange(-4, 4, 10000))

	ms = lguys.calc_M.(profile, radii)
	ms ./= ms[end]

	icdf = lguys.lerp(ms, radii)

	p = rand(N)
	return icdf.(p)
end

# ╔═╡ d586226b-37dc-4f56-b632-d8d411848401
R_s=π^2 / 100

# ╔═╡ 97ae40b8-5b48-4828-a139-549d8cc9a124
profile = lguys.Exp2D(R_s=R_s)

# ╔═╡ 1cc8aaaa-b20b-401c-be78-6162da324d7e
radii = sample_radii(profile, N)

# ╔═╡ 194f9836-58cf-43ec-9689-f30bcc1feba1
pos = radii' .* lguys.rand_unit(N)

# ╔═╡ 18e0fe22-54d0-4f94-bb87-2f4c5a30a145
vel = zeros(size(pos));

# ╔═╡ 444ba26c-1f99-467e-a304-3879190dc03a
hist(radii)

# ╔═╡ 07e7f1a4-8cea-4832-8d61-5a0e6ef5f48d
snap = lguys.Snapshot(pos, vel, ones(N) / N)

# ╔═╡ 73661ad2-caa1-478a-984b-d1061c031dc4
r_prof, rho_prof = lguys.calc_ρ_hist(radii, 50, weights=snap.masses)

# ╔═╡ b2ea83f1-375c-41c5-8eab-4a436e4c46bb
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xscale=log10,
		yscale=log10,
		xlabel="radius / kpc",
		ylabel=L"\rho",
	)

	scatter!(Arya.midpoint(r_prof), rho_prof)

	x = 10 .^ LinRange(-2, 0, 1000)
	y = lguys.calc_ρ.(profile, x)
	
	vlines!(R_s)
	vlines!(lguys.get_r_h(profile))
	r_h1 = lguys.percentile(radii, 50)
	
	vlines!(r_h1)

	
	lines!(x, y)
	fig
end

# ╔═╡ 1befbfe1-65d4-4473-95aa-6cc456714a11
R = @. sqrt(pos[1, :]^2 + pos[2, :]^2)

# ╔═╡ 6cce5fb9-c12f-4af6-b188-a79c53c37fa5
props = lguys.calc_properties(R, bins=100)

# ╔═╡ 32e87484-508c-4280-8d38-c9275ba70c74
let
	fig = Figure()
	ax = Axis(fig[1, 1],
	)

	scatter!(props.log_r, props.log_Sigma)

	x = 10 .^ LinRange(-3, 2.2, 1000)
	y = lguys.calc_Σ.(profile, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[2])
	fig
end

# ╔═╡ 3637933e-7454-44da-a396-5d9e9e47dffd
skypos = lguys.to_sky(snap)

# ╔═╡ 4b698f55-932e-4c5a-b03b-f633f5abdd30
ra = [o.ra for o in skypos]

# ╔═╡ 313c1e22-470b-40c4-adc3-25bd76689e79
dec = [o.dec for o in skypos]

# ╔═╡ 6a2b4ab4-65f2-4e87-b688-031230580736
scatter(ra, dec)

# ╔═╡ b8a32ed2-c7ea-4d0d-a9b4-59e97a31db61
md"""
Should be centred at Sag A*
"""

# ╔═╡ f8f890a1-2dc6-49f5-8756-bac581af478c
ra0, dec0 = lguys.calc_centre2D(ra, dec, snap.masses, "mean")

# ╔═╡ f5ff4ee3-c437-4152-82ce-b8747dd4a6a3
xi, eta = lguys.to_tangent(ra, dec, ra0, dec0)

# ╔═╡ af89d81c-cd8c-4306-ab93-790ce45f1631
let
	global r_ell
	ecc = PA = 0
	
	b = sqrt(1-ecc)
	a = 1/b
	
	r_ell = lguys.calc_r_ell(xi, eta, a, b, PA-90)

	r_ell .*= 60 # to arcmin
end

# ╔═╡ a75a7819-3d25-4455-bd05-f8482b79d2e5
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel="xi", ylabel="eta",
		aspect=DataAspect()
	)

	
	scatter!(xi, eta)

	fig
end

# ╔═╡ 81871062-7ac7-41d4-a321-2d83424c88f5
π/(10800)

# ╔═╡ d9a0d436-2f4b-4efa-877d-bbbe4a5b556a
π/180/60

# ╔═╡ 55a9ac8c-c153-4354-9f95-e34b91ae1612
60 / 206265

# ╔═╡ 24fd05e6-1d01-4909-b348-21192b26861c
props_sky = lguys.calc_properties(r_ell)

# ╔═╡ 27aa2cd1-b365-4701-b698-fb2fb7d0f3eb
distance = 8.122

# ╔═╡ 3ab11ec8-1899-486c-97c5-0133607309c8
R_s_arcmin = R_s *40 / (distance / 86)

# ╔═╡ b0bfab5e-8ca7-4ffa-9b43-07fa9f1d11f4
lguys.kpc_to_arcmin(R_s, distance)

# ╔═╡ 94d54949-1118-40fe-9e68-62b716bbc1ae
prof_sky = lguys.Exp2D(R_s=R_s_arcmin)

# ╔═╡ 132b0c72-27fd-48e1-86b6-6a5291a2c145
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log r / arcmin",
		ylabel=L"\log \Sigma"
	)

	scatter!(props_sky.log_r, props_sky.log_Sigma)

	x = 10 .^ LinRange(-1, 2.7, 1000)
	y = lguys.calc_Σ.(prof_sky, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[2])
	fig
end

# ╔═╡ ede57bf6-c009-4c77-a439-e41121c3beca
filename = "test_sky_recon.fits"

# ╔═╡ 9965c814-bf78-4b0b-9818-778029786afa
let
	f = FITS(filename, "w")

	df = Dict(
		"ra" => ra,
		"dec" => dec
	)
	write(f, df)

	println("wrote to $filename")
	close(f)
	df
end

# ╔═╡ 4661b78d-2ea0-45c3-8e2e-c1d5bcf3d7f7
md"""
run calc_2D_profile on above and then fit_profile
"""

# ╔═╡ d3e8bd8d-47ce-41b9-b0bb-d6fc0b5d51ac
R_s_fit = 43.04

# ╔═╡ 5a3357dc-4794-4824-a319-43f568d9c0c0
lguys.arcmin_to_kpc(R_s_fit, distance)

# ╔═╡ 73947b81-a32b-4f3c-be80-f62e7d620cb4
R_s

# ╔═╡ 772bcaa1-9e99-471e-bbd0-7316a4ac65b9
prof_recon = lguys.ObsProfile("test_sky_recon_profile.toml")

# ╔═╡ be5b28f0-c5d3-4e32-8244-5c100c880eeb
prof_sky_recon = lguys.Exp2D(R_s=R_s_fit)

# ╔═╡ edfb1e71-90b7-407f-bfdc-8bb3d0dba95e
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log r / arcmin",
		ylabel=L"\log \Sigma"
	)

	errscatter!(prof_recon.log_r, prof_recon.log_Sigma, yerr=prof_recon.log_Sigma_err)

	x = 10 .^ LinRange(-1, 2.7, 1000)
	y = lguys.calc_Σ.(prof_sky, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[2])


	x = 10 .^ LinRange(-1, 2.7, 1000)
	y = lguys.calc_Σ.(prof_sky_recon, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[3])
	fig
end

# ╔═╡ Cell order:
# ╠═43307c80-298c-11ef-3cfe-afb7123432ca
# ╠═1c55d5eb-1065-4ef0-bb97-ada81c195432
# ╠═1093880a-8973-4aa9-bed2-25376166d705
# ╠═d586226b-37dc-4f56-b632-d8d411848401
# ╠═97ae40b8-5b48-4828-a139-549d8cc9a124
# ╠═1cc8aaaa-b20b-401c-be78-6162da324d7e
# ╠═194f9836-58cf-43ec-9689-f30bcc1feba1
# ╠═18e0fe22-54d0-4f94-bb87-2f4c5a30a145
# ╠═444ba26c-1f99-467e-a304-3879190dc03a
# ╠═07e7f1a4-8cea-4832-8d61-5a0e6ef5f48d
# ╠═73661ad2-caa1-478a-984b-d1061c031dc4
# ╠═b2ea83f1-375c-41c5-8eab-4a436e4c46bb
# ╠═1befbfe1-65d4-4473-95aa-6cc456714a11
# ╠═6cce5fb9-c12f-4af6-b188-a79c53c37fa5
# ╠═32e87484-508c-4280-8d38-c9275ba70c74
# ╠═3637933e-7454-44da-a396-5d9e9e47dffd
# ╠═4b698f55-932e-4c5a-b03b-f633f5abdd30
# ╠═313c1e22-470b-40c4-adc3-25bd76689e79
# ╠═6a2b4ab4-65f2-4e87-b688-031230580736
# ╟─b8a32ed2-c7ea-4d0d-a9b4-59e97a31db61
# ╠═f8f890a1-2dc6-49f5-8756-bac581af478c
# ╠═f5ff4ee3-c437-4152-82ce-b8747dd4a6a3
# ╠═af89d81c-cd8c-4306-ab93-790ce45f1631
# ╠═a75a7819-3d25-4455-bd05-f8482b79d2e5
# ╠═81871062-7ac7-41d4-a321-2d83424c88f5
# ╠═d9a0d436-2f4b-4efa-877d-bbbe4a5b556a
# ╠═55a9ac8c-c153-4354-9f95-e34b91ae1612
# ╠═24fd05e6-1d01-4909-b348-21192b26861c
# ╠═27aa2cd1-b365-4701-b698-fb2fb7d0f3eb
# ╠═3ab11ec8-1899-486c-97c5-0133607309c8
# ╠═b0bfab5e-8ca7-4ffa-9b43-07fa9f1d11f4
# ╠═94d54949-1118-40fe-9e68-62b716bbc1ae
# ╠═132b0c72-27fd-48e1-86b6-6a5291a2c145
# ╠═624a165e-73f2-40a5-86e7-852a4ec405ca
# ╠═ede57bf6-c009-4c77-a439-e41121c3beca
# ╠═9965c814-bf78-4b0b-9818-778029786afa
# ╠═4661b78d-2ea0-45c3-8e2e-c1d5bcf3d7f7
# ╠═d3e8bd8d-47ce-41b9-b0bb-d6fc0b5d51ac
# ╠═5a3357dc-4794-4824-a319-43f568d9c0c0
# ╠═73947b81-a32b-4f3c-be80-f62e7d620cb4
# ╠═772bcaa1-9e99-471e-bbd0-7316a4ac65b9
# ╠═be5b28f0-c5d3-4e32-8244-5c100c880eeb
# ╠═edfb1e71-90b7-407f-bfdc-8bb3d0dba95e
