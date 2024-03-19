### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()
	
	using FITSIO
	using DataFrames 
	using CSV
	
	import LilGuys as lguys
end

# ╔═╡ d31ccfbf-6628-4745-a0ca-fd2061821416
begin
	using KernelDensity
	using Turing, Distributions, StatsPlots
end

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	using Plots; 
	using Arya
 	Arya.set_default()
	gr()
# 	pythonplot()
	using LaTeXStrings

end

# ╔═╡ ebfe3e7b-c791-4699-b039-09c4be31ea0d
begin 
	P08 = CSV.read("first_apo.csv", DataFrame)
	P08_no_tides = CSV.read("no_tides.csv", DataFrame)
end

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	f = FITS("../Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits")
	j24 = DataFrame(f[2])
	f2 = FITS("B22_sculptor.fits")
	b22 = DataFrame(f2[2])
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
		xlabel = "ξ / ˚", ylabel = "η / ˚"
	)
	histogram2d!(xi, eta,)
end

# ╔═╡ 468588eb-73ac-4be5-bfcd-5dcb3c79a8aa
begin 
	dθ = 0.5*rh/60
	nb = 20
	histogram2d(xi, eta,
		aspect_ratio=:equal, bins=(LinRange(-dθ, dθ, nb), LinRange(-dθ, dθ, nb)),
	xlim=(-dθ, dθ))
	xlabel!("ξ / ˚")
	ylabel!("η / ˚")
end

# ╔═╡ 1942cac0-2da8-4bec-80d0-a9302ddc9ea1
begin 
	kde_d = kde((xi, eta))
	ik = InterpKDE(kde_d)
	Σ_kde = @. pdf([ik], kde_d.x', kde_d.y)
end

# ╔═╡ 5786e566-fc49-42d6-9032-0ce4a7d4497a
scatter(xi, eta, marker_z= pdf.([ik], xi, eta), xlims=(-1, 1), ylims=(-1, 1), ms=1)

# ╔═╡ c6115ac4-57de-43f1-abfe-adfde1294bb5
function plot_circle!(radius, x=0, y=0; kwargs...)
	t = LinRange(0, 2π, 10000)
	plot!(x .+ radius*cos.(t), y .+ radius*sin.(t); kwargs...)
end

# ╔═╡ 74327a32-0084-4a3c-986e-0bd0b81996f9
begin 
	plot(aspect_ratio=:equal, xlim=(-dθ, dθ), ylim=(-dθ, dθ), legend=false, dpi=300)
	
	scatter!(xi, eta, 
		ms=1, label="")
	
	contour!(kde_d.x, kde_d.y, Σ_kde, 
		lw=2, colorbar=false)
	
	plot_circle!(2rh/60, 
		lw=2, ls=:dot, color=:black)

	θ_text = 45
	r_text = 2rh/60
	annotate!(r_text * cosd(θ_text), r_text * sind(θ_text), Plots.text(L"$2 r_{h}$", 12,  rotation = -90+θ_text, valign=:bottom))

	xlabel!(L"$\xi$ / degree")
	ylabel!(L"$\eta$ / degree")
	# savefig("sculptor_tangent_plane.png")
end

# ╔═╡ 08135ca7-78c7-4f40-a014-cdb9a0833188
begin 
	scatter(xi, eta, aspect_ratio=:equal, xlim=(-1, 1), ylim=(-1, 1), ms=1, label="", )
	contour!(kde_d.x, kde_d.y, Σ_kde, lw=1, levels=10)
	plot_circle!(2rh/60, lw=2, ls=:dot, color=:black, label="r_h")

	xlabel!("ξ / ˚")
	ylabel!("η / ˚")
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
	θ_bins = LinRange(-π, π, 30)
	θ_bm = lguys.midpoint(θ_bins)
end

# ╔═╡ d84e9494-bf96-4dd4-b6c9-6317d6ac24a7
lguys.calc_histogram(vec(kde_theta), θ_bins , weights=vec(Σ_kde))

# ╔═╡ 51c415b2-c5c8-46db-81f6-13d242c79484
begin 
	r_cuts_2 = [0, 0.3, 0.6, 1, 3, 10]
	p1 = plot(xlabel="θ", ylabel="relative mean surface density in annulus")
	
	for i in 1:length(r_cuts_2)-1
		filt_r = r_cuts_2[i] .< kde_r .< r_cuts_2[i+1]
		y = lguys.calc_histogram(vec(kde_theta[filt_r]), θ_bins, weights=vec(Σ_kde[filt_r]))[2]
		
		plot!(θ_bm,y ./ maximum(y), label="r = $(r_cuts_2[i]) - $(r_cuts_2[i+1])˚")
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
histogram2d(j24_filtered.pmra, j24_filtered.pmdec, xlabel=L"$\mu_{\alpha *}$ / mas yr$^{-1}$", ylabel=L"$\mu_\delta$")

# ╔═╡ 83b506bb-f385-464e-8f5c-7adfff45105a
begin 
	scatter(j24_filtered.bp_rp, j24_filtered.phot_g_mean_mag, marker_z=log10.(r_ell), 
	yflip=true, xlabel="bp rp", ylabel="G")

end

# ╔═╡ b057629d-cc8a-48bd-ba5b-6d2203f988ba
histogram2d(j24_filtered.parallax, j24_filtered.parallax_error)

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
	μ_xi = mean(xi)
	μ_eta = mean(eta)
	r_mean =  @.sqrt.((xi .- μ_xi).^2 + (eta - μ_eta).^2)
end

# ╔═╡ c463f463-30ac-4a3d-91ad-65ed428b8a2f
begin 
	p[] = plot( xscale=:log10, yscale=:log10, xlabel="radius (degree)", ylabel="cumulative star count")

	plot_cdf!(r_circ, label="circular", norm=:count)
	plot_cdf!(r_mean, label="mean centre",  norm=:count)

	p[]

end

# ╔═╡ 64be0b75-9aa5-4b1a-a505-78b6bb43f261
mean(xi)

# ╔═╡ e7ac218f-a88d-4775-84b2-6ca0d8a6c74f
std(xi) / sqrt(length(xi))

# ╔═╡ 06ee46e0-ab72-40dc-8aad-1f7f1cebc031
std(eta) / sqrt(length(eta))

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
ϕ = atan.(eta .- mean(eta), xi .- mean(xi))

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

# ╔═╡ d1519ebd-f240-43cf-827b-e0960a7d3448
kde_prof = kde(log10.(r_ell), bandwidth=0.05)

# ╔═╡ 8c1ca1f4-70d5-486e-b8c0-5eaeea1aa096
begin 
	plot(xlabel="log r / degree", ylabel="density")
	plot!(kde_prof.x, kde_prof.density, label="kde")
	bins2, ys2 = lguys.calc_histogram(log10.(r_ell), 29)

	ys2 ./= sum(ys2 .* diff(bins2))
	scatter!(lguys.midpoint(bins2), ys2, label="histogram")
end

# ╔═╡ b94bffa3-7bee-401e-a88e-b04867c2e52c
begin
	d_kde = lguys.gradient(kde_prof.x, kde_prof.density)
	r_kde = 10 .^ kde_prof.x
	rm_kde = lguys.midpoint(r_kde)
	σ_kde =  1 ./ (π* diff(r_kde .^ 2) ) .* lguys.midpoint(kde_prof.density) .* diff(kde_prof.x)
end

# ╔═╡ 36e7c3f1-9dc7-4233-a348-082cc7f623b1
begin 
	plot(xlabel="log r / degrees", ylabel="log density", ylims=(-8, 5))
	plot!(log10.(lguys.midpoint(r_kde)), log10.(σ_kde), label="kde")
	rs2 = 10 .^ bins2
	rm2 = lguys.midpoint(rs2)
	areas2 = π * diff(rs2 .^ 2)
	σs = ys2 .* diff(bins2) ./ areas2
	σ_errs = sqrt.(ys2 * length(ys2)) ./ length(ys2) .* diff(bins2) ./ areas2
	scatter!(lguys.midpoint(bins2), log10.(σs), yerr=σ_errs, label="histogram")
	scatter!(log10.(r_ell), -8 .+ 0.3*randn(length(r_ell)), label="stars", ms=1, alpha=1)
end

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

# ╔═╡ b237411c-236d-4aab-a5c4-167575a2f74b


# ╔═╡ 2f7fbf40-c8b0-417c-bf8a-96fe375d9ede
import NaNMath as nm

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

# ╔═╡ 351c4370-78ff-43e0-9062-a8597329d1b1
function calc_Γ(rs, σs)
	xs = nm.log10.(rs)
	ys = nm.log10.(σs)
	return lguys.gradient(xs, ys)
end

# ╔═╡ 8ad632cd-28ee-4613-a806-deeec4ffb03f
begin 
	plot(xlabel="log r / degrees", ylabel="Γ")
	plot!(log10.(rm2), calc_Γ(rm2, σs), label="histogram" )
	plot!(log10.(rm_kde), calc_Γ(rm_kde, σ_kde), label="kde" )
	vline!([log10.(rh/60)], label="")
end

# ╔═╡ Cell order:
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═d31ccfbf-6628-4745-a0ca-fd2061821416
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═ebfe3e7b-c791-4699-b039-09c4be31ea0d
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═f6a5c138-10bc-48e2-aec7-45fd61b9b17f
# ╠═c54224d3-c528-41bb-bce6-c910fc680b55
# ╠═5e5afcb8-9344-4213-a5b0-33f437ea0263
# ╠═eba107b4-c322-42a7-ad68-fb3ea1330f61
# ╠═32293baf-0c39-488d-b495-a86f42eed178
# ╠═ba71616e-dadf-4025-9afd-66dc40d4e65b
# ╠═01b26f1d-f6bd-4e20-a344-0b55a2f91bf4
# ╠═37a0afd4-73b9-4d8b-b65a-eabae4d01bba
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
# ╠═cd9ed83d-6827-4809-8929-c0e58ae5da57
# ╠═83b506bb-f385-464e-8f5c-7adfff45105a
# ╠═b057629d-cc8a-48bd-ba5b-6d2203f988ba
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
# ╠═64be0b75-9aa5-4b1a-a505-78b6bb43f261
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
# ╠═d1519ebd-f240-43cf-827b-e0960a7d3448
# ╠═8c1ca1f4-70d5-486e-b8c0-5eaeea1aa096
# ╠═b94bffa3-7bee-401e-a88e-b04867c2e52c
# ╠═36e7c3f1-9dc7-4233-a348-082cc7f623b1
# ╠═697a3f87-143b-44f1-9edc-dafbeb220f4e
# ╠═e00404f0-68bb-4cff-89f0-fc7ae9b78e99
# ╠═351c4370-78ff-43e0-9062-a8597329d1b1
# ╠═8ad632cd-28ee-4613-a806-deeec4ffb03f
# ╠═b237411c-236d-4aab-a5c4-167575a2f74b
# ╠═2f7fbf40-c8b0-417c-bf8a-96fe375d9ede
