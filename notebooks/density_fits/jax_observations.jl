### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames 
	using CSV
	using CairoMakie

	import NaNMath as nm
	using KernelDensity
	using Measurements

end

# ╔═╡ 69c98029-165c-407b-9a63-a27e06e30e45
using Arya

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
Do observations from different sources ll agree?
"""

# ╔═╡ 489f6d21-9e9a-4b1e-b05c-c63a44ba1951
import StatsBase: percentile, mean

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
end

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	f = FITS("./data/Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits")
	all_stars = DataFrame(f[2])
end

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = all_stars[all_stars.PSAT .> 0.2, :]

# ╔═╡ f3606421-b7a1-40cc-a743-8bd28169779f
xi0 = mean(members.xi)

# ╔═╡ 89b03ca3-c25c-4dc8-9710-7137e02e8d2c
eta0 = mean(members.eta)

# ╔═╡ 731ea468-5003-44e9-95b8-7fa7ef4b475b
Nmemb = size(members, 1)

# ╔═╡ bf1d1b84-1d11-4997-b083-32e94a01ea68
snap = lguys.Snapshot(hcat(members.xi, members.eta, zeros(Nmemb))', zeros(3, Nmemb), ones(Nmemb))

# ╔═╡ 223abc41-5158-49c2-96bf-df55b7be1114
cen = lguys.Centres.calc_centre(lguys.Centres.SS_State, snap)

# ╔═╡ 3e40a3f8-efbc-4807-b185-22fbb2e99adf
begin 
	xi = members.xi .- cen.position[1]
	eta = members.eta .- cen.position[2]
end

# ╔═╡ 02d19c07-411f-40d6-9ac3-010ebfd4bdfe
let 
	dθ = 10
	f = Figure()
	ax = Axis(f[1, 1], xlabel=L"\xi / \textrm{arcmin}", ylabel=L"\eta / \textrm{arcmin}", aspect=1, limits=(-dθ, dθ, -dθ, dθ), xgridvisible=false, ygridvisible=false)
	
	scatter!(60*all_stars.xi, 60*all_stars.eta, color=all_stars.phot_g_mean_mag, colormap=:greys, markersize=3)

	f
end

# ╔═╡ c54224d3-c528-41bb-bce6-c910fc680b55
hexbin(all_stars[:, "ra"], all_stars[:, "dec"], bins=300)

# ╔═╡ 4f538de6-8827-454f-a97f-8f6c2cd7ea3f
hexbin(members.ra, members.dec, bins=300)

# ╔═╡ 1942cac0-2da8-4bec-80d0-a9302ddc9ea1
begin 
	kde_d = kde((xi, eta))
	ik = InterpKDE(kde_d)
	Σ_kde = @. pdf([ik], kde_d.x', kde_d.y)
	Σ_kde = Σ_kde'
end

# ╔═╡ c6115ac4-57de-43f1-abfe-adfde1294bb5
function plot_circle!(radius, x=0, y=0; kwargs...)
	t = LinRange(0, 2π, 10000)
	lines!(x .+ radius*cos.(t), y .+ radius*sin.(t); kwargs...)
end

# ╔═╡ cdd246ae-e6bb-49a9-9512-085bd37827bd
let 
	dθ = 1

	fig, ax, p = heatmap(kde_d.x, kde_d.y, Σ_kde, colormap=:greys, colorscale=x -> asinh(x/1e-1))
	ax.aspect = 1
	contour!(kde_d.x, kde_d.y, asinh.(Σ_kde ./ 1), levels=10 )

	ax.limits = (-dθ, dθ,-dθ, dθ)
	ax.aspect=1

	filt = members.r_ell .> 1.5
	scatter!(ax, xi[filt], eta[filt], markersize=3, color=Arya.COLORS[3])

	Colorbar(fig[1, 2], p, ticks=[10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0])
	fig
end

# ╔═╡ 9c5a7fd6-fb24-4db7-8c42-21191ecfdc16
minimum(asinh.(Σ_kde[Σ_kde .> 01e-22]))

# ╔═╡ d49c5f63-e2a4-4ff5-af4a-d4239b50abae
begin
	kde_y = reshape(kde_d.y, (1, :))
	kde_r = @. sqrt((kde_d.x')^2 + (kde_d.y)^2)
	kde_theta = atan.(kde_d.y, kde_d.x')
end

# ╔═╡ 042f8215-b6fe-49db-a166-055567da4bb0
begin 
	θ_bins = LinRange(-π, π, 15)
	θ_bm = lguys.midpoint(θ_bins)
end

# ╔═╡ d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
lguys.calc_histogram(vec(kde_theta), θ_bins , weights=vec(Σ_kde))

# ╔═╡ 81374c1a-cf4a-4504-b669-7a2bbe5a6b5c
let
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
	title = "sculptor members")
	
	hexbin!(ax, members.bp_rp, members.phot_g_mean_mag, bins=100)
	f
end

# ╔═╡ 83b506bb-f385-464e-8f5c-7adfff45105a
begin 
	scatter(members.bp_rp, members.phot_g_mean_mag, color=log10.(members.r_ell), 
	)

end

# ╔═╡ b057629d-cc8a-48bd-ba5b-6d2203f988ba
begin 
	Arya.hist2d(members.parallax, members.parallax_error, bins=100)
end

# ╔═╡ a4e857c0-39d9-4fa0-871f-ccc66cb17c25
begin 
	#plot(colorbar_scale=:log10, xlabel=L"\varpi / \rm mas", ylabel=L"\delta\varpi/\rm mas"
	#)
	Arya.hist2d(all_stars.parallax, all_stars.parallax_error, bins=100)
end

# ╔═╡ c70fce56-9367-4f6a-92bb-f0d0d7a616e0
function calc_Σ(rs)
	#r_bins = lguys.make_equal_number_bins(rs, 500)
	r_bins = 10 .^ LinRange(log10(minimum(rs)), log10(maximum(rs)), 50)
	
	N_bins = length(r_bins) - 1
	ν = zeros(N_bins)
	r_mid = zeros(N_bins)
	for i in 1:(length(r_bins) - 1)
		fl = r_bins[i] .< rs .< r_bins[i+1]
		ν[i] = sum(fl)
		r_mid[i] = lguys.mean(rs[fl])
	end
	
	Areas = π  * diff(r_bins.^2 )
	Σs = ν ./ Areas
	return r_mid, Σs
end

# ╔═╡ cd8f63a6-bd53-4cf4-a3f4-a8943e0a1caf
scatter(members[:, :phot_g_mean_mag], members.r_ell)

# ╔═╡ 9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
let
	dθ = 1
	#plot(xlabel="ξ", ylabel="η", xlim=(-dθ, dθ), ylim=(-dθ, dθ), aspect_ratio=:equal)

	f = Figure()
	ax = Axis(f[1, 1], limits=(-dθ,dθ,-dθ,dθ), aspect=1,
		xlabel=L"\xi",
		ylabel=L"\eta")
	for a in LinRange(0.05dθ, dθ, 8)
		t = LinRange(0, 360, 100)
		x = @. a*cosd(t)
		y = @. a*sind(t)
		lines!(x, y, label="", color="black")
	end
	scatter!(members.xi, members.eta, markersize=2,label="", alpha=0.3)
	f
end

# ╔═╡ ec1c1408-5c3c-4bfa-b8bd-f32dbc3379a9
ϕ = atan.(eta, xi)

# ╔═╡ 9282a867-b9c7-4a50-b58a-da3ceb068cad
radii = @. sqrt(xi^2 + eta^2)

# ╔═╡ 6bb8fe50-d7c0-4d76-9d05-a403b7e3239f
rh = percentile(radii, 50)

# ╔═╡ 6ddd0b7b-6c51-442e-97c4-ddc321899da1
60rh

# ╔═╡ 468588eb-73ac-4be5-bfcd-5dcb3c79a8aa
let
	f = Figure()
	dθ = 0.5*rh

	ax = Axis(f[1,1], 
		xlabel = L"\xi / \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		limits = (-dθ, dθ, -dθ, dθ)
	)
	scatter!(ax, all_stars.xi, all_stars.eta)
	f
end

# ╔═╡ 74327a32-0084-4a3c-986e-0bd0b81996f9
let
	dθ = 0.7
	fig = Figure()
	ax = Axis(fig[1,1],
		aspect=1, limits=(-dθ, dθ, -dθ, dθ),
		xlabel = L"\xi/ \textrm{degree}",
		ylabel = L"\eta / \textrm{degree}",
		xgridvisible=false,
		ygridvisible=false
	)

	
	#scatter!(ax, xi, eta, markersize=4, color=Arya.COLORS[3])
	
	contourf!(kde_d.x, kde_d.y, asinh.(Σ_kde), levels=130 , colormap=:magma,)
	
	#plot_circle!(rh, 
	#	lw=2, ls=:dot, color=:black)

	θ_text = 45
	r_text = 2rh/60

	# annotate!(r_text * cosd(θ_text), r_text * sind(θ_text), Plots.text(L"$2 r_{h}$", 12,  rotation = -90+θ_text, valign=:bottom))
	fig
end

# ╔═╡ 4ff675ef-8609-49df-bfac-4070f14e3c25
maximum(radii)

# ╔═╡ 52d3db7f-8bbf-43cf-89f2-16134247d713
let
	fig = Figure()
	ax = PolarAxis(fig[1,1], rlimits=(0, 0.5))

	scatter!(ax, ϕ, radii, markersize=2)
	fig
end

# ╔═╡ 2771e601-5dc9-4115-bba1-9096aa99ee28
mean(sqrt.(all_stars.xi .^ 2 + all_stars.eta .^ 2) ./ all_stars.r_ell) * 60 

# ╔═╡ b3334ad6-f7bf-4ca2-b895-a582472f2e90
ecc = 0.37

# ╔═╡ 0090b5d8-cbd9-4353-af3a-f793bc71ec37
let
	fig = Figure()
	ax = PolarAxis(fig[1,1])
	ϕ = atan.(all_stars.eta, all_stars.xi)
	scatter!(ax, ϕ, all_stars.r_ell, markersize=2)

	bins = LinRange(0, 2π, 20)
	for i in 1:length(bins)-1
		filt = ϕ .> bins[i]
		filt .&= ϕ .< bins[i+1]
		if sum(filt) > 0
			println(sum(filt))
			println(bins[i], "\t", maximum(all_stars.r_ell[filt]))
		end
	end
	
	fig
end

# ╔═╡ 0f36be65-b713-4dcd-a8a4-6c9898e9d36d
15.379 / 9.786

# ╔═╡ 9ac3b2c6-1e4c-4223-9a31-84aee5ef9a54
1/(1-ecc)

# ╔═╡ 9427aa86-2d46-4890-b2c5-36de3ca83264
r_cuts = [percentile(radii, p) for p in LinRange(0, 100, 8)]

# ╔═╡ a563a2c5-bbca-4785-9024-41741119bf03
percentile(members.r_ell, 0)

# ╔═╡ ac3c0cf9-b00b-4f57-9e36-837806ab733c
maximum(members.r_ell)

# ╔═╡ 75b8cbb1-519e-4174-a7b9-22e661cbde5c
N_cuts = length(r_cuts) - 1

# ╔═╡ 4a6f9b97-3e14-4683-b338-e1bd67a07ea2
hist(log10.(members.r_ell))

# ╔═╡ a37b2c44-9976-4bf8-8392-891f352a5da0
ecdfplot(members.r_ell)

# ╔═╡ ebb24519-ee2d-4781-a2bd-16ccbc957060
let
	f = Figure(size=(1200, 500))
	ax = Axis(f[1, 1], xlabel=L"\theta / \textrm{radians}", ylabel="normalized density")
	bins = 10
	ps = []
	labels = []
	for i in 1:N_cuts
		filt = r_cuts[i] .<= radii .<= r_cuts[i+1] 
		println(sum(filt))
		if sum(filt) > 0
			bins, counts = lguys.calc_histogram(ϕ[filt], bins)
			counts .*= 1/sum(diff(bins) .* counts) * sum(diff(bins))
 			p = lines!(ax, lguys.midpoint(bins), counts, color=i, colorrange=(1, N_cuts))
			push!(ps, p)

			r_l = round(r_cuts[i], digits=3)
			r_h = round(r_cuts[i+1], digits=3)
			push!(labels, L"r \in [%$r_l, %$r_h ]")
		end
	end

	Legend(f[1,2], ps, labels)
	f
end

# ╔═╡ 9d58a063-a995-44eb-877a-6e8c62b377eb


# ╔═╡ 7af6fd4f-a7f4-4e7f-a1e3-8899717dd2cd
let
	f = Figure(size=(1200, 500))
	ax = Axis(f[1, 1], xlabel=L"\theta / \textrm{radians}", ylabel="normalized density for rell")
	bins = 10
	ps = []
	labels = []
	for i in 1:N_cuts
		filt = r_cuts[i] .<= members.r_ell * rh .<= r_cuts[i+1] 
		println(sum(filt))
		if sum(filt) > 0
			bins, counts = lguys.calc_histogram(ϕ[filt], bins)
			counts .*= 1/sum(diff(bins) .* counts) * sum(diff(bins))
 			p = lines!(ax, lguys.midpoint(bins), counts, color=i, colorrange=(1, N_cuts))
			push!(ps, p)

			r_l = round(r_cuts[i], digits=3)
			r_h = round(r_cuts[i+1], digits=3)
			push!(labels, L"rell \in [%$r_l, %$r_h ]")
		end
	end

	Legend(f[1,2], ps, labels)
	f
end

# ╔═╡ 49a5de24-907c-49e1-9747-8be5bb28b64d
let 
	fig, ax, h = hist(ϕ, bins=50, label="")
	ax.xlabel = "angle"
	ax.ylabel = "count"
	fig
end

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═489f6d21-9e9a-4b1e-b05c-c63a44ba1951
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═6bb8fe50-d7c0-4d76-9d05-a403b7e3239f
# ╠═6ddd0b7b-6c51-442e-97c4-ddc321899da1
# ╠═f3606421-b7a1-40cc-a743-8bd28169779f
# ╠═89b03ca3-c25c-4dc8-9710-7137e02e8d2c
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═bf1d1b84-1d11-4997-b083-32e94a01ea68
# ╠═223abc41-5158-49c2-96bf-df55b7be1114
# ╠═3e40a3f8-efbc-4807-b185-22fbb2e99adf
# ╠═4ff675ef-8609-49df-bfac-4070f14e3c25
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═02d19c07-411f-40d6-9ac3-010ebfd4bdfe
# ╠═c54224d3-c528-41bb-bce6-c910fc680b55
# ╠═4f538de6-8827-454f-a97f-8f6c2cd7ea3f
# ╠═468588eb-73ac-4be5-bfcd-5dcb3c79a8aa
# ╠═1942cac0-2da8-4bec-80d0-a9302ddc9ea1
# ╠═c6115ac4-57de-43f1-abfe-adfde1294bb5
# ╠═cdd246ae-e6bb-49a9-9512-085bd37827bd
# ╠═9c5a7fd6-fb24-4db7-8c42-21191ecfdc16
# ╠═74327a32-0084-4a3c-986e-0bd0b81996f9
# ╠═d49c5f63-e2a4-4ff5-af4a-d4239b50abae
# ╠═042f8215-b6fe-49db-a166-055567da4bb0
# ╠═d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
# ╠═81374c1a-cf4a-4504-b669-7a2bbe5a6b5c
# ╠═83b506bb-f385-464e-8f5c-7adfff45105a
# ╠═b057629d-cc8a-48bd-ba5b-6d2203f988ba
# ╠═a4e857c0-39d9-4fa0-871f-ccc66cb17c25
# ╠═c70fce56-9367-4f6a-92bb-f0d0d7a616e0
# ╠═cd8f63a6-bd53-4cf4-a3f4-a8943e0a1caf
# ╠═9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
# ╠═ec1c1408-5c3c-4bfa-b8bd-f32dbc3379a9
# ╠═9282a867-b9c7-4a50-b58a-da3ceb068cad
# ╠═52d3db7f-8bbf-43cf-89f2-16134247d713
# ╠═2771e601-5dc9-4115-bba1-9096aa99ee28
# ╠═b3334ad6-f7bf-4ca2-b895-a582472f2e90
# ╠═0090b5d8-cbd9-4353-af3a-f793bc71ec37
# ╠═0f36be65-b713-4dcd-a8a4-6c9898e9d36d
# ╠═9ac3b2c6-1e4c-4223-9a31-84aee5ef9a54
# ╠═9427aa86-2d46-4890-b2c5-36de3ca83264
# ╠═a563a2c5-bbca-4785-9024-41741119bf03
# ╠═ac3c0cf9-b00b-4f57-9e36-837806ab733c
# ╠═75b8cbb1-519e-4174-a7b9-22e661cbde5c
# ╠═4a6f9b97-3e14-4683-b338-e1bd67a07ea2
# ╠═a37b2c44-9976-4bf8-8392-891f352a5da0
# ╠═ebb24519-ee2d-4781-a2bd-16ccbc957060
# ╠═9d58a063-a995-44eb-877a-6e8c62b377eb
# ╠═7af6fd4f-a7f4-4e7f-a1e3-8899717dd2cd
# ╠═49a5de24-907c-49e1-9747-8be5bb28b64d
