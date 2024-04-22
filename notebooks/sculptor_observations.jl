### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames 
	using CSV
	using Plots; gr()

	import NaNMath as nm
	using KernelDensity
	using Measurements

end

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using Arya
end

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
Do observations from different sources ll agree?
"""

# ╔═╡ ebfe3e7b-c791-4699-b039-09c4be31ea0d
begin 
	P08 = CSV.read("first_apo.csv", DataFrame)
	P08_no_tides = CSV.read("no_tides.csv", DataFrame)
end

# ╔═╡ 32293baf-0c39-488d-b495-a86f42eed178
begin 
	ra0 = 15.03916666
	dec0 = -33.7091666
	ecc = 0.37
	rh = 12.33 # arcmin +- 0.05 # Fedrico's papper
	a = rh / 60 # deg
	b = (1-ecc) * a
	PA = 94 #position angle

end

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	f = FITS("../Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits")
	j24 = DataFrame(f[2])
	f2 = FITS("B22_sculptor.fits")
	b22 = DataFrame(f2[2])

	f3 = FITS("../data/sculptor_gaia_all.fits")
	scl_all = DataFrame(f3[2])
	scl_all[:, "xi"], scl_all[:, "eta"] = lguys.to_tangent(scl_all.ra, scl_all.dec, ra0, dec0, )
end

# ╔═╡ d1a2aa9e-13d3-4cd1-9326-3039991eebe3
begin 
	plot(xlabel=L"\xi / \textrm{degrees}", ylabel=L"\eta / \textrm{degrees}", aspect_ratio=1)
	scatter!(scl_all.xi, scl_all.eta, ms=1, alpha=0.1, label="gaia all")
end

# ╔═╡ 1d7bbfed-bb40-469b-86e4-411c3fbfdec8
begin 
	plot(xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}", aspect_ratio=1, xlim=(-10, 10), ylim=(-10, 10), dpi=400, fontfamily="Computer Modern", title="Gaia All")
	
	scatter!(60*scl_all.xi, 60*scl_all.eta, ms=2, alpha=1, label="", marker_z=scl_all.phot_g_mean_mag, cmap=cgrad(:greys), cb_title="G mag", clims=(8, 22))
	
end

# ╔═╡ 02d19c07-411f-40d6-9ac3-010ebfd4bdfe
begin 
	plot(xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}", aspect_ratio=1, xlim=(-10, 10), ylim=(-10, 10), dpi=400, fontfamily="Computer Modern", title="J24 All")
	
	scatter!(60*j24.xi, 60*j24.eta, ms=2, alpha=1, label="", marker_z=j24.phot_g_mean_mag, cmap=cgrad(:greys), cb_title="G mag", clims=(8, 22))
	
end

# ╔═╡ f6a5c138-10bc-48e2-aec7-45fd61b9b17f
begin 
	filt = j24.PSAT .> 0.2
	j24_filtered = j24[filt, :]

	filt_b22 = b22.Pmemb .> 0.2
	b22_filtered = b22[filt_b22, :]
end

# ╔═╡ c54224d3-c528-41bb-bce6-c910fc680b55
histogram2d(j24[:, "ra"], j24[:, "dec"], bins=300)

# ╔═╡ 5e5afcb8-9344-4213-a5b0-33f437ea0263
histogram2d(b22[:, "RA_ICRS"], b22[:, "DE_ICRS"], bins=300)

# ╔═╡ ba71616e-dadf-4025-9afd-66dc40d4e65b
begin 
	xi, eta = lguys.to_tangent(j24_filtered.ra, j24_filtered.dec, ra0, dec0,)
	r_ell = lguys.calc_r_ell(xi, eta, a, b, PA-90) * sqrt(a*b)

	xi_b, eta_b = lguys.to_tangent(b22_filtered.RA_ICRS, b22_filtered.DE_ICRS, ra0, dec0, )

	b22_filtered[:, "r_ell"] = lguys.calc_r_ell(xi_b, eta_b, a, b, PA-90)
end

# ╔═╡ eba107b4-c322-42a7-ad68-fb3ea1330f61
plot(log10.(sort(r_ell)))

# ╔═╡ 01b26f1d-f6bd-4e20-a344-0b55a2f91bf4
begin 
	stephist(j24.PSAT, bins=20)
	# stephist!(b22.Pmemb)
end

# ╔═╡ 37a0afd4-73b9-4d8b-b65a-eabae4d01bba
begin 
	plot(
		aspect_ratio=:equal, xlims=(-1.2, 1.2), ylims=(-1, 1),
		xlabel = L"\xi / \textrm{degree}", ylabel = L"\eta / \textrm{degree}"
	)
	histogram2d!(xi, eta,)
end

# ╔═╡ 930c5c65-46d9-4b33-af6f-0727fac63c5f


# ╔═╡ 468588eb-73ac-4be5-bfcd-5dcb3c79a8aa
begin 
	dθ = 0.5*rh/60
	nb = 20
	histogram2d(xi, eta,
		aspect_ratio=:equal, bins=(LinRange(-dθ, dθ, nb), LinRange(-dθ, dθ, nb)),
	xlim=(-dθ, dθ))
	xlabel!(L"\xi / \textrm{degree}")
	ylabel!(L"\eta / \textrm{degree}")
end

# ╔═╡ 1942cac0-2da8-4bec-80d0-a9302ddc9ea1
begin 
	kde_d = kde((xi, eta))
	ik = InterpKDE(kde_d)
	Σ_kde = @. pdf([ik], kde_d.x', kde_d.y)
end

# ╔═╡ 5786e566-fc49-42d6-9032-0ce4a7d4497a
# ╠═╡ disabled = true
#=╠═╡
scatter(xi, eta, marker_z= pdf.([ik], xi, eta), xlims=(-1, 1), ylims=(-1, 1), ms=1)
  ╠═╡ =#

# ╔═╡ c6115ac4-57de-43f1-abfe-adfde1294bb5
function plot_circle!(radius, x=0, y=0; kwargs...)
	t = LinRange(0, 2π, 10000)
	plot!(x .+ radius*cos.(t), y .+ radius*sin.(t); kwargs...)
end

# ╔═╡ 08135ca7-78c7-4f40-a014-cdb9a0833188
begin 
	scatter(xi, eta, aspect_ratio=:equal, xlim=(-1, 1), ylim=(-1, 1), ms=1, label="", )
	contour!(kde_d.x, kde_d.y, Σ_kde, lw=1, levels=10)
	plot_circle!(2rh/60, lw=2, ls=:dot, color=:black, label=L"$r_{h}$")

	xlabel!(L"\xi / \rm degree")
	ylabel!(L"\eta / \rm degree")
end

# ╔═╡ d49c5f63-e2a4-4ff5-af4a-d4239b50abae
begin
	kde_y = reshape(kde_d.y, (1, :))
	kde_r = @. sqrt((kde_d.x')^2 + (kde_d.y)^2)
	kde_theta = atan.(kde_d.y, kde_d.x')
end

# ╔═╡ 3fd48123-0ba9-4300-8406-76c36da982c5
begin 
	plot(proj=:polar)
	scatter!(kde_theta, kde_r, proj=:polar, label="", marker_z = Σ_kde, ms=3, ylims=(0, 0.5))
end

# ╔═╡ 042f8215-b6fe-49db-a166-055567da4bb0
begin 
	θ_bins = LinRange(-π, π, 15)
	θ_bm = lguys.midpoint(θ_bins)
end

# ╔═╡ d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
lguys.calc_histogram(vec(kde_theta), θ_bins , weights=vec(Σ_kde))

# ╔═╡ 51c415b2-c5c8-46db-81f6-13d242c79484
begin 
	r_cuts_2 = [0, 0.1, 0.2, 0.3, 0.5, 3]
	p1 = plot(xlabel="θ", ylabel="relative mean surface density in annulus")
	
	for i in 1:length(r_cuts_2)-1
		filt_r = r_cuts_2[i] .< kde_r .< r_cuts_2[i+1]
		y = lguys.calc_histogram(vec(kde_theta[filt_r]), θ_bins, weights=vec(Σ_kde[filt_r]))[2]
		
		plot!(θ_bm,y ./ lguys.mean(y), label="r = $(r_cuts_2[i]) - $(r_cuts_2[i+1])˚")
	end
	p1
end

# ╔═╡ e49c409b-188f-4d65-a193-288a5fc8e656
kde_d.x[argmax(Σ_kde)[2]]

# ╔═╡ e696721c-f3a6-4812-a17a-2062b3f957e7
kde_d.y[argmax(Σ_kde)[1]]

# ╔═╡ 40ad7745-df35-4898-af0f-8b24ee752aae
# ╠═╡ disabled = true
#=╠═╡
begin 
	#plot()
	dpm = 10
	histogram2d(data.pmra, data.pmdec, bins=30, xlims=(-10, 10), ylims=(-10, 10))
		#xlims=(-dpm, dpm), ylims=(-dpm, dpm))
	
	scatter!(data[filt, "pmra"], data[filt, "pmdec"], alpha=0.1)

end
  ╠═╡ =#

# ╔═╡ 47370a21-5127-4b53-84d5-b6307c7cd585
begin 
	plot(xlabel="rell", ylabel="count / bin")
	scatterhist!(log10.(j24_filtered.r_ell), label="J+24")
	scatterhist!(log10.(b22_filtered.r_ell), label="B+22")
end

# ╔═╡ 0e6d1577-cc27-40a5-a483-aacc6f735a2f
begin
	plot(xlabel="datapoint", ylabel="rell(jax) / rell(me)")
	scatter!(j24_filtered.r_ell ./ r_ell .* sqrt(a*b), alpha=0.1, label="") 
end

# ╔═╡ 1d5f749f-8ed3-4f92-b8c2-4d4bee8f46f1
gr()

# ╔═╡ 9e95d81c-5305-4b67-9867-b0e04efd63aa
histogram2d(j24_filtered.pmra, j24_filtered.pmdec, xlabel=L"\mu_{\alpha *} / \rm mas\, yr^{-1}", ylabel=L"\mu_\delta", aspect_ratio=:equal,)

# ╔═╡ 3dccd820-de5c-422c-8585-7c74e6a09b07
histogram2d(j24_filtered.pmra, j24_filtered.pmdec, bins=LinRange(-2, 2, 100), xlabel=L"\mu_{\alpha *} / \rm mas\, yr^{-1}", ylabel=L"\mu_\delta", colorbar_scale=:log10, aspect_ratio=1, xlims=(-2, 2))

# ╔═╡ 856a5528-7941-4691-b723-31ad3da7786e
begin 
	plot(xlabel=L"\mu_{\alpha *} / \rm mas\, yr^{-1}", ylabel=L"\mu_\delta  / \rm mas\, yr^{-1}", colorbar_scale=:log10, aspect_ratio=:equal)
	
	histogram2d!(j24.pmra, j24.pmdec, bins=LinRange(-20, 20, 100), xlim=(-20, 20))
end

# ╔═╡ 81374c1a-cf4a-4504-b669-7a2bbe5a6b5c
begin 
	plot(colorbar_scale=:log10, yflip=true,
	xlabel="BP - RP", ylabel="G")
	
	histogram2d!(j24_filtered.bp_rp, j24_filtered.phot_g_mean_mag, bins=100)
end

# ╔═╡ eb071752-34ac-4e97-9444-9f9bca62dc3b
begin 
	plot(colorbar_scale=:log10,yflip=true,
	xlabel="BP - RP", ylabel="G")
	histogram2d!(j24.bp_rp, j24.phot_g_mean_mag, bins=200)
end

# ╔═╡ 83b506bb-f385-464e-8f5c-7adfff45105a
begin 
	scatter(j24_filtered.bp_rp, j24_filtered.phot_g_mean_mag, marker_z=log10.(r_ell), 
	yflip=true, xlabel="bp rp", ylabel="G")

end

# ╔═╡ b057629d-cc8a-48bd-ba5b-6d2203f988ba
begin 
	plot(colorbar_scale=:log10, xlabel=L"\varpi / \rm mas", ylabel=L"\delta\varpi/\rm mas"
	)
	histogram2d!(j24_filtered.parallax, j24_filtered.parallax_error)
end

# ╔═╡ a4e857c0-39d9-4fa0-871f-ccc66cb17c25
begin 
	plot(colorbar_scale=:log10, xlabel=L"\varpi / \rm mas", ylabel=L"\delta\varpi/\rm mas"
	)
	histogram2d!(j24.parallax, j24.parallax_error, colorbar_scale=:log10)
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

# ╔═╡ 24b30dd5-a2a8-47b2-9695-8b48f6bf2667
begin 
	r_mid, Σs = calc_Σ(r_ell)
	#r_mid_b, Σs_b = calc_Σ(filtered_b22.r_ell)
end

# ╔═╡ 596d3898-d7b5-42f6-bcda-fe2e91c64a05
begin 
	plot(  )
	dx = -log10(rh)
	dy = 5.5
	scatter!(log10.(r_mid), log10.(Σs), label="sculptor", msa=0)
	#scatter!(log10.(r_mid_b * rh), log10.(Σs_b/rh^2), label="b22", msa=0)

	plot!(P08.x .+ dx, P08[!, " y"] .+ dy,  label="P08 tides")
	plot!(P08_no_tides.x .+ dx, P08_no_tides[!, " y"] .+ dy, ls=:dash, label="P08 no tides")
	vline!([log10(rh/60)], label="", color="black", ls=:dot)
	annotate!(log10(rh/60),  2, Plots.text(L"$r_{h}$", 12,  rotation = 90, valign=:top, halign=:right), )

	xlabel!(L"$\log r_{\rm ell}$ / degrees")
	ylabel!(L"surface density (stars deg$^{-2}$)")
	savefig("sculptor_profile.png")
end

# ╔═╡ cd8f63a6-bd53-4cf4-a3f4-a8943e0a1caf
scatter(j24_filtered[:, :phot_g_mean_mag], r_ell, ma=0.1, ms=3)

# ╔═╡ 802348e2-8e11-42f5-b2e8-4a66a39936a2
md"""
# Cumulative plots
"""

# ╔═╡ 54477bfa-bfc0-4d44-9f03-72531ec3de7d
function plot_cdf!(rs, weights=ones(length(rs)); norm=:density, kwargs...)
	idx = sortperm(rs)
	if norm == :density
		ys = cumsum(weights[idx]) / sum(weights[idx])
	elseif norm == :count
		ys = cumsum(weights[idx])
	else
		error("Norm not known : $norm")
	end
	
	plot!(rs[idx], ys; kwargs... )
	
	return rs[idx],ys
end

# ╔═╡ cd9ed83d-6827-4809-8929-c0e58ae5da57
begin 
	p2 = plot()
	r_cuts = [0, 0.5, 1, 1.5, 2, 3, 5]
	for i in 1:(length(r_cuts) - 1)
		
		r_filt = r_cuts[i] .< r_ell .< r_cuts[i + 1]
		plot_cdf!(j24_filtered[r_filt, :phot_g_mean_mag], ones(sum(r_filt)) / sum(r_filt), 
			label="$(r_cuts[i]) < r/r_h < $(r_cuts[i+1])")
	end

	xlabel!("G mag")
	ylabel!("CDF")
	p2
end

# ╔═╡ 786f8e2c-3b12-47a9-b431-18e30cdf2a1e
idxs = sortperm(r_ell)

# ╔═╡ 38553483-6028-4d78-9805-80581e55e214
M = cumsum(ones(length(idxs)))

# ╔═╡ 20d997f1-72de-4c6b-9aa3-836c85221c21
p = Ref{Plots.Plot}()

# ╔═╡ b16a4815-22c6-44cd-b0d9-67a349500637
begin 
	p[] = plot(xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}", aspect_ratio=1, xlim=(-10, 10), ylim=(-10, 10), dpi=400, fontfamily="Computer Modern", title="J24 members")
	
	scatter!(60*j24_filtered.xi, 60*j24_filtered.eta, ms=2, alpha=1, label="", marker_z=j24_filtered.phot_g_mean_mag, cmap=cgrad(:greys), cb_title="G mag", clims=(8, 22), dpi=1600, show=true)
	p[]
	savefig("j24_memb.png")
	p[]
end

# ╔═╡ 74327a32-0084-4a3c-986e-0bd0b81996f9
begin 
	p[] = plot(aspect_ratio = 1, ticks=:naitive, extra_kwargs=:subplot, xlim=(-dθ, dθ), ylim=(-dθ, dθ), legend=false, dpi=300
)
	
	scatter!(xi, eta, 
		ms=1, label="")
	
	contour!(kde_d.x, kde_d.y, Σ_kde, 
		lw=2, colorbar=false)
	
	plot_circle!(2rh/60, 
		lw=2, ls=:dot, color=:black)

	θ_text = 45
	r_text = 2rh/60
	annotate!(r_text * cosd(θ_text), r_text * sind(θ_text), Plots.text(L"$2 r_{h}$", 12,  rotation = -90+θ_text, valign=:bottom))

	xlabel!(L"\xi/ \rm degree")
	ylabel!(L"\eta / \rm degree")
	#savefig("sculptor_tangent_plane.tex")
	p[]
end

# ╔═╡ 71b76e09-565d-48e9-bfb3-2d457fc0c8af
r_circ = @. sqrt(j24_filtered.xi ^2 + j24_filtered.eta ^ 2)

# ╔═╡ 2de3acc9-0af4-408f-a3fb-4f3e7fc21514
begin 
	p[] = plot( xscale=:log10, yscale=:log10, xlabel="radius (degree)", ylabel="CDF")
	rs, cdfs = plot_cdf!(r_ell[idxs], label="elliptical", norm=:count)
	plot_cdf!(r_circ, label="circular", norm=:count)
	p[]
end

# ╔═╡ 523f6558-ff6b-4543-9612-554b97d78029
r2 = @. sqrt((xi + 0.0028) ^2 + (eta -0.00195)  ^ 2)

# ╔═╡ c1a82bee-8ac0-4038-8a72-76ebfd85e2d6
begin
	μ_xi = lguys.mean(xi)
	μ_eta = lguys.mean(eta)
	r_mean =  @.sqrt.((xi .- μ_xi).^2 + (eta - μ_eta).^2)
end

# ╔═╡ c463f463-30ac-4a3d-91ad-65ed428b8a2f
begin 
	p[] = plot( xscale=:log10, yscale=:log10, xlabel="radius (degree)", ylabel="cumulative star count")

	plot_cdf!(r_circ, label="circular", norm=:count)
	plot_cdf!(r_mean, label="mean centre",  norm=:count)

	p[]

end

# ╔═╡ e7ac218f-a88d-4775-84b2-6ca0d8a6c74f
lguys.std(xi) / sqrt(length(xi))

# ╔═╡ 06ee46e0-ab72-40dc-8aad-1f7f1cebc031
lguys.std(eta) / sqrt(length(eta))

# ╔═╡ 86ed9a86-e4d5-4f63-9642-fb28638af89b
begin
	plot(xlabel="log r / kpc", ylabel="log Σ")

	bins, ys = lguys.calc_histogram(log10.(r_ell), 50)
	areas = π * diff((10 .^ bins) .^ 2)

	plot!(lguys.midpoint(bins), log10.(ys ./ areas), label="elliptical")
	
	bins, ys = lguys.calc_histogram(log10.(r_circ), 50)
	plot!(lguys.midpoint(bins), log10.(ys ./areas), label="circular")
end

# ╔═╡ 9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
begin 
	plot(xlabel="ξ", ylabel="η", xlim=(-dθ, dθ), ylim=(-dθ, dθ), aspect_ratio=:equal)

	for a in LinRange(0.1dθ, 2dθ, 8)
		plot!(x->a*cosd(x + PA), x->a*(1-ecc)*sind(x + PA), LinRange(0, 360, 100), label="", color="black")
	end
	scatter!(xi, eta, ms=1,label="", alpha=0.3)
end

# ╔═╡ 03bd2cca-dc14-41a7-b1d5-02b31b41a6b1
ϕ = atan.(eta .- lguys.mean(eta), xi .- lguys.mean(xi))

# ╔═╡ b797849c-3617-4f85-a977-8493f0a694ac
scatter(ϕ, log10.(r_ell),  proj = :polar, ms=0.5, ylims=(-2.3, 0.2), label="")

# ╔═╡ 0f7fe2ac-94f5-43c5-9bae-9a49e7f832ef
scatter(ϕ .+ 1*(π .- 2π * rand(length(r_ell))), (r_ell),  proj = :polar, ms=0.5, ylims=(0, 1), label="", yticks=:none)

# ╔═╡ 877801fa-1062-47ae-851c-f5c12f4c49ed
scatter(xi, eta, ms=3, xlim=(-0.15, 0.15), ylim=(-0.15, 0.15))

# ╔═╡ 36f7623b-a1bb-4971-880d-45f99bcfd302
rand(10)

# ╔═╡ ebb24519-ee2d-4781-a2bd-16ccbc957060
begin
	h = plot()
	for r in 10 .^ LinRange(-1.3, 0, 10)
		if sum(0.9r .< r_ell .< 1.1r) > 0
			scatterhist!(ϕ[0.9r .< r_ell .< 1.1r], bins=20, label="$r", norm=:pdf)
		end
	end
	h
end

# ╔═╡ 49a5de24-907c-49e1-9747-8be5bb28b64d
begin 
	plot(xlabel="angle", ylabel="count")
	histogram!(ϕ, bins=50, label="")
end

# ╔═╡ 4a493f03-594a-4a4b-811a-6001926d44ae
md"""
# KDE density
"""

# ╔═╡ 48369c0e-f9ae-4a75-b884-4d3aeef0f87e
function gaussian_kernel(x, x0=0, bw=1, w=1; nsd=5)
	dx = x .- x0
	return @. ifelse(abs(dx) > nsd*bw, 0,
		w*1/√(2π) * 1/bw * exp(-dx^2/2bw^2)
	)
end

# ╔═╡ 7b4d350a-aca2-4c7c-a1c7-0a25dbea3eeb
function uniform_kernel(x, x0=0, bw=1, w=1)
	dx = x .- x0
	return @. ifelse(abs(dx) > bw, 0,
		w*1/2bw
	)
end

# ╔═╡ 11590585-48a0-4ec7-a5d4-f266cb4c90f5
import NearestNeighbors as nn

# ╔═╡ 179414e1-f51b-4e86-a606-f5757c6b4da4
"""
	calc_akde(positions; kernel, bw, N)

Assuming positions are sorted, calculates the adaptive KDE 
"""
function calc_akde(positions; kernel=gaussian_kernel, bw=1, k=5, N=10_000, cutoff=5)
	x = positions[(!).(isnan.(positions))]
	xmin = minimum(x)
	xmax = maximum(x)
	x_sample = LinRange(xmin, xmax, N)

	tree = nn.BallTree(reshape(x, (1, :)))
	hs = last.(nn.knn(tree, reshape(x, (1, :)), k, true)[2]) 
	hs .*= bw / sqrt(k)

	pdf = zeros(N)

	for i in 1:length(x)
		xl = searchsortedlast(x_sample, x[i] - cutoff*hs[i])
		xh = searchsortedfirst(x_sample, x[i] + cutoff*hs[i])

		xl = max(1, xl)
		xh = min(xh, length(x_sample))
		xl = 1
		xh = length(x_sample)

		pdf[xl:xh] .+= 1/hs[i]^2 .* [kernel((x_sample[j] - x[i])/ hs[i]) for j in xl:xh]

	end
		
	dx = (xmax - xmin) / N
	pdf ./= sum(pdf * dx)
	return x_sample, pdf, hs
end

# ╔═╡ 149d9308-e6e5-40c8-bc57-d01b44db2ad1
"""
	calc_kde(positions; kernel, bw, N)

Assuming positions are sorted, can calculate the kde
"""
function calc_kde(positions; kernel=gaussian_kernel, bw=1, N=1_000, cutoff=5)
	x = positions[(!).(isnan.(positions))]
	xmin = minimum(x)
	xmax = maximum(x)
	x_sample = LinRange(xmin, xmax, N)


	pdf = zeros(N)

	for i in 1:length(positions)
		xl = searchsortedlast(x_sample, positions[i] - cutoff*bw)
		xh = searchsortedfirst(x_sample, positions[i] + cutoff*bw)

		xl = max(1, xl)
		xh = min(xh, length(x_sample))

		pdf[xl:xh] .+= 1/bw^2* [kernel((x_sample[j] - positions[i])/ bw) for j in xl:xh]

	end
		
	dx = (xmax - xmin) / N
	pdf ./= sum(pdf * dx)
	return x_sample, pdf
end

# ╔═╡ e399f3f1-d623-45a3-a809-86e7db9a71ce
log10.(r_ell)

# ╔═╡ 45264d8f-2c81-4d9a-a6b4-737e384f01a7
histogram(log10.(r_ell))

# ╔═╡ 1bd6e99f-2a85-4301-a23b-e0f24ba75219
log10.(r_ell)

# ╔═╡ 7f6c0ac3-7111-45ee-9945-6c98c491208c
length(r_ell)

# ╔═╡ c71e1ee7-714a-4455-89b2-d479e7384be3
x, y, hs = calc_akde(sort((r_ell)), bw=4, k=500, kernel=gaussian_kernel)

# ╔═╡ 1e38e2f1-38e8-4813-a56e-b4fce91ea635
begin 
	scatter(sort(log10.(r_ell)), log10.(hs))
	scatter!(log10.(lguys.midpoint(sort(r_ell))), log10.(250*lguys.diff(sort(r_ell))), alpha=0.1)
	plot!([-3, 0.5], [-3, 0.5])
end

# ╔═╡ ad904c0e-cfe9-46cf-96b8-f4cb0093dff5
begin 
	plot(log10.(x), 5*y .* x .* log(10) ./ lguys.gradient(log10.(x)) /length(r_ell) )
	plot!(calc_kde(sort(log10.(r_ell)), bw=0.01)...)
	scatterhist!(log10.(r_ell), norm=:pdf)
end

# ╔═╡ 50138654-cfd9-4dc4-b35d-80f3f375960c
begin 
	plot()
	plot!(x, y )
	scatterhist!(r_ell, norm=:pdf)
	plot!(kde(r_ell))
end

# ╔═╡ 8c1ca1f4-70d5-486e-b8c0-5eaeea1aa096
begin 
	p[] = plot(xlabel="log r / degree", ylabel="counts / bin", yscale=:log10, ylims=(1, 1e3))
	#plot!(kde_prof.x, kde_prof.density .* length(r_ell)/9, label="kde")
	bins2, ys2 = lguys.calc_histogram(log10.(r_ell), 29)
	ys2_err = sqrt.(ys2)
	print(ys2_err)
	
	scatter!(lguys.midpoint(bins2), ys2, yerr=ys2_err, label="", lw=2)
	ys2_err = ys2_err ./ sum(ys2 .* diff(bins2))

	ys2 ./= sum(ys2 .* diff(bins2))
	p[]
end

# ╔═╡ 0c9d5bfc-a218-45a8-b342-58d9b7affef7


# ╔═╡ 36e7c3f1-9dc7-4233-a348-082cc7f623b1
begin 
	plot(xlabel="log r / degrees", ylabel="log density", ylims=(-10, 5))
	#plot!(log10.(lguys.midpoint(r_kde)), log10.(σ_kde), label="kde")
	rs2 = 10 .^ bins2
	rm2 = lguys.midpoint(rs2)
	areas2 = π * diff(rs2 .^ 2)
	σs = ys2 .* diff(bins2) ./ areas2
	σ_errs = ys2_err .* diff(bins2) ./ areas2
	logσ_e = hcat([log10.(1 .- σ_errs ./ σs), log10.(1 .+ σ_errs ./ σs)])
	scatter!(lguys.midpoint(bins2), log10.(σs), yerr=transpose(logσ_e), label="histogram", lw=1)
	scatter!(log10.(r_ell), -8 .+ 0.3*randn(length(r_ell)), label="stars", ms=1, alpha=1)
end

# ╔═╡ bf3e6ca0-d3dd-43da-89cc-5bb6d15ba0ce
maximum(r_ell)

# ╔═╡ 8863b312-eb5c-40a5-a638-138b8be214e5
ys2

# ╔═╡ f30253cd-8b8d-414f-8565-75f68c6688a9
ys2_err

# ╔═╡ 6a5ff6f9-fdc6-40d4-98b7-7066a005c54e
logσ_e

# ╔═╡ 697a3f87-143b-44f1-9edc-dafbeb220f4e
begin 
	plot(xlabel="log r / degrees", ylabel="log count")
	bins3, ys3 = lguys.calc_histogram(log10.(r_ell), 29)

	ys_un = ys3 * length(ys3)
	ys_un_err = sqrt.(ys_un)
	yerr = hcat(ys_un .- ys_un_err, ys_un .+ ys_un_err)
	scatter!(lguys.midpoint(bins3), log10.(ys_un), yerr=log10.(yerr ./ ys_un), lw=1, label="members")

	bins3, ys3 = lguys.calc_histogram(log10.(j24.r_ell[j24.PSAT .< 0.2] * √(a*b)), 29)
	ys_un = ys3 * length(ys3)
	ys_un_err = sqrt.(ys_un)
	yerr = hcat(ys_un .- ys_un_err, ys_un .+ ys_un_err)
	scatter!(lguys.midpoint(bins3), log10.(ys_un), yerr=log10.(yerr ./ ys_un), lw=1, label="nonmembers")

end

# ╔═╡ e00404f0-68bb-4cff-89f0-fc7ae9b78e99
begin
	p3 = plot(xlabel="log r / degrees", ylabel="log density")
	for pcut in [0, 0.01, 0.3, 0.5, 0.9, 0.99]
		bins_nm, ys_nm = lguys.calc_histogram(log10.(j24.r_ell[pcut .< j24.PSAT]) .+ log10(sqrt(a*b)), 50)
	
		ys_nm .*= diff(bins_nm)
		areas = π * diff((10 .^bins_nm) .^ 2)
		plot!(lguys.midpoint(bins_nm) , nm.log10.(ys_nm ./ areas), label="p > $pcut", line_z=pcut, cmap=cgrad(:blues, rev = true), colorbar=false)
	end
	p3
end 

# ╔═╡ 37e0e1b9-2e7b-48e2-9245-0060581b13c2
begin
	p[] = plot(xlabel="log r / degrees", ylabel="log density")
	for Gcut in [18, 19, 20, 21]
		bins_nm, ys_nm = lguys.calc_histogram(log10.(j24_filtered.r_ell[Gcut .> j24_filtered.phot_g_mean_mag]) .+ log10(sqrt(a*b)), 50)
	
		ys_nm .*= diff(bins_nm)
		areas = π * diff((10 .^bins_nm) .^ 2)
		plot!(lguys.midpoint(bins_nm) , nm.log10.(ys_nm ./ areas), label="G < $Gcut", line_z=Gcut, cmap=cgrad(:blues, rev = true), colorbar=false)
	end
	p[]
end 

# ╔═╡ Cell order:
# ╠═47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═ebfe3e7b-c791-4699-b039-09c4be31ea0d
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═d1a2aa9e-13d3-4cd1-9326-3039991eebe3
# ╠═1d7bbfed-bb40-469b-86e4-411c3fbfdec8
# ╠═02d19c07-411f-40d6-9ac3-010ebfd4bdfe
# ╠═b16a4815-22c6-44cd-b0d9-67a349500637
# ╠═f6a5c138-10bc-48e2-aec7-45fd61b9b17f
# ╠═c54224d3-c528-41bb-bce6-c910fc680b55
# ╠═5e5afcb8-9344-4213-a5b0-33f437ea0263
# ╠═eba107b4-c322-42a7-ad68-fb3ea1330f61
# ╠═32293baf-0c39-488d-b495-a86f42eed178
# ╠═ba71616e-dadf-4025-9afd-66dc40d4e65b
# ╠═01b26f1d-f6bd-4e20-a344-0b55a2f91bf4
# ╠═37a0afd4-73b9-4d8b-b65a-eabae4d01bba
# ╠═930c5c65-46d9-4b33-af6f-0727fac63c5f
# ╠═468588eb-73ac-4be5-bfcd-5dcb3c79a8aa
# ╠═1942cac0-2da8-4bec-80d0-a9302ddc9ea1
# ╠═5786e566-fc49-42d6-9032-0ce4a7d4497a
# ╠═c6115ac4-57de-43f1-abfe-adfde1294bb5
# ╠═74327a32-0084-4a3c-986e-0bd0b81996f9
# ╠═08135ca7-78c7-4f40-a014-cdb9a0833188
# ╠═d49c5f63-e2a4-4ff5-af4a-d4239b50abae
# ╠═3fd48123-0ba9-4300-8406-76c36da982c5
# ╠═042f8215-b6fe-49db-a166-055567da4bb0
# ╠═d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
# ╠═51c415b2-c5c8-46db-81f6-13d242c79484
# ╠═e49c409b-188f-4d65-a193-288a5fc8e656
# ╠═e696721c-f3a6-4812-a17a-2062b3f957e7
# ╠═40ad7745-df35-4898-af0f-8b24ee752aae
# ╠═47370a21-5127-4b53-84d5-b6307c7cd585
# ╠═0e6d1577-cc27-40a5-a483-aacc6f735a2f
# ╠═24b30dd5-a2a8-47b2-9695-8b48f6bf2667
# ╠═596d3898-d7b5-42f6-bcda-fe2e91c64a05
# ╠═1d5f749f-8ed3-4f92-b8c2-4d4bee8f46f1
# ╠═9e95d81c-5305-4b67-9867-b0e04efd63aa
# ╠═3dccd820-de5c-422c-8585-7c74e6a09b07
# ╠═856a5528-7941-4691-b723-31ad3da7786e
# ╠═81374c1a-cf4a-4504-b669-7a2bbe5a6b5c
# ╠═eb071752-34ac-4e97-9444-9f9bca62dc3b
# ╠═83b506bb-f385-464e-8f5c-7adfff45105a
# ╠═cd9ed83d-6827-4809-8929-c0e58ae5da57
# ╠═b057629d-cc8a-48bd-ba5b-6d2203f988ba
# ╠═a4e857c0-39d9-4fa0-871f-ccc66cb17c25
# ╠═c70fce56-9367-4f6a-92bb-f0d0d7a616e0
# ╠═cd8f63a6-bd53-4cf4-a3f4-a8943e0a1caf
# ╠═802348e2-8e11-42f5-b2e8-4a66a39936a2
# ╠═54477bfa-bfc0-4d44-9f03-72531ec3de7d
# ╠═786f8e2c-3b12-47a9-b431-18e30cdf2a1e
# ╠═38553483-6028-4d78-9805-80581e55e214
# ╠═20d997f1-72de-4c6b-9aa3-836c85221c21
# ╠═2de3acc9-0af4-408f-a3fb-4f3e7fc21514
# ╠═c463f463-30ac-4a3d-91ad-65ed428b8a2f
# ╠═71b76e09-565d-48e9-bfb3-2d457fc0c8af
# ╠═523f6558-ff6b-4543-9612-554b97d78029
# ╠═c1a82bee-8ac0-4038-8a72-76ebfd85e2d6
# ╠═e7ac218f-a88d-4775-84b2-6ca0d8a6c74f
# ╠═06ee46e0-ab72-40dc-8aad-1f7f1cebc031
# ╠═86ed9a86-e4d5-4f63-9642-fb28638af89b
# ╠═9002fcb7-ffb2-4e2b-a4e9-57599d6e233d
# ╠═03bd2cca-dc14-41a7-b1d5-02b31b41a6b1
# ╠═b797849c-3617-4f85-a977-8493f0a694ac
# ╠═0f7fe2ac-94f5-43c5-9bae-9a49e7f832ef
# ╠═877801fa-1062-47ae-851c-f5c12f4c49ed
# ╠═36f7623b-a1bb-4971-880d-45f99bcfd302
# ╠═ebb24519-ee2d-4781-a2bd-16ccbc957060
# ╠═49a5de24-907c-49e1-9747-8be5bb28b64d
# ╠═4a493f03-594a-4a4b-811a-6001926d44ae
# ╠═48369c0e-f9ae-4a75-b884-4d3aeef0f87e
# ╠═7b4d350a-aca2-4c7c-a1c7-0a25dbea3eeb
# ╠═11590585-48a0-4ec7-a5d4-f266cb4c90f5
# ╠═179414e1-f51b-4e86-a606-f5757c6b4da4
# ╠═149d9308-e6e5-40c8-bc57-d01b44db2ad1
# ╠═e399f3f1-d623-45a3-a809-86e7db9a71ce
# ╠═45264d8f-2c81-4d9a-a6b4-737e384f01a7
# ╠═1bd6e99f-2a85-4301-a23b-e0f24ba75219
# ╠═7f6c0ac3-7111-45ee-9945-6c98c491208c
# ╠═c71e1ee7-714a-4455-89b2-d479e7384be3
# ╠═1e38e2f1-38e8-4813-a56e-b4fce91ea635
# ╠═ad904c0e-cfe9-46cf-96b8-f4cb0093dff5
# ╠═50138654-cfd9-4dc4-b35d-80f3f375960c
# ╠═8c1ca1f4-70d5-486e-b8c0-5eaeea1aa096
# ╠═0c9d5bfc-a218-45a8-b342-58d9b7affef7
# ╠═36e7c3f1-9dc7-4233-a348-082cc7f623b1
# ╠═bf3e6ca0-d3dd-43da-89cc-5bb6d15ba0ce
# ╠═8863b312-eb5c-40a5-a638-138b8be214e5
# ╠═f30253cd-8b8d-414f-8565-75f68c6688a9
# ╠═6a5ff6f9-fdc6-40d4-98b7-7066a005c54e
# ╠═697a3f87-143b-44f1-9edc-dafbeb220f4e
# ╠═e00404f0-68bb-4cff-89f0-fc7ae9b78e99
# ╠═37e0e1b9-2e7b-48e2-9245-0060581b13c2
