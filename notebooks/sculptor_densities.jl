### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()
	
	using FITSIO
	using Plots; 
	using Arya
	
	using DataFrames 
	using CSV
	using LaTeXStrings
	
	import LilGuys as lguys
end

# ╔═╡ 40eb3191-3c32-4928-9344-05b48da4bdd0
using Optim

# ╔═╡ d31ccfbf-6628-4745-a0ca-fd2061821416
begin
	using Turing, Distributions, StatsPlots
end

# ╔═╡ 02b0d62c-2a47-4308-ab96-ede30029b8c8
begin 
	using AbstractGPs, TemporalGPs
	using BSplines
end

# ╔═╡ 4cc65ca0-7a21-4dc0-b8f2-df1cdd55f915
using LogExpFunctions

# ╔═╡ 98c47808-7b07-4438-a855-d3f0167f3f34
md"""
what are some reasonable ways of calculating densities given points
"""

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
Arya.set_default()

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	f = FITS("../Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits")
	data = DataFrame(f[2])
	#f2 = FITS("b19.fits")
	#b22 = DataFrame(f2[2])
end

# ╔═╡ f6a5c138-10bc-48e2-aec7-45fd61b9b17f
begin 
	filt = data.PSAT .> 0.2
	filtered = data[filt, :]
end

# ╔═╡ 6ca5439b-b4d7-44d5-b7bd-303f14ee7c15
function to_tangent(α_0, δ_0, α, δ)
	denom = @. (sind(δ) * sind(δ_0) 
		+ cosd(δ) * cosd(δ_0) * cosd(α - α_0)
	)
	
	eta_num = @. (sind(δ_0) * cosd(δ) * cosd(α-α_0)
		-cosd(δ_0) * sind(δ) 
	)
	
	xi_num = @. cosd(δ) * sind(α - α_0)
	
	xi = @. rad2deg(xi_num/denom)
	eta = @. rad2deg(eta_num / denom)

	return xi, -eta
end

# ╔═╡ d04f8836-bf9a-4a7b-8ec9-a1ba783c7295
function calc_r_ell(x, y, a, b, PA)
	θ = @. deg2rad(PA)
	x_p = @. x * cos(θ) + -y * sin(θ)
	y_p = @. x * sin(θ) + y * cos(θ)

	r_sq = @. (x_p / a)^2 + (y_p / b)^2
	return sqrt.(r_sq)
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

# ╔═╡ ba71616e-dadf-4025-9afd-66dc40d4e65b
begin 
	xi, eta = to_tangent(ra0, dec0, filtered.ra, filtered.dec)
	r_ell = calc_r_ell(xi, eta, a, b, PA-90) * sqrt(a*b)
end

# ╔═╡ efe133bc-670a-49a5-9f22-c937456a11b5
md"""
# Regressions to the profile
"""

# ╔═╡ 19344fb8-7eb4-4f56-af6a-14441226cabd
not = !

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

# ╔═╡ 0bf61650-adae-43e0-8caf-02951aa6082c
r_mid, Σs = calc_Σ(r_ell)

# ╔═╡ 695ac014-86ee-48cc-a6e6-6d8a0c444be1
begin 
	calc_Γ(rs, Σs) = @. (log10(Σs[2:end] / Σs[1:end-1] )) / log10(rs[2:end] / rs[1:end-1])
	
	Γs = calc_Γ(r_mid, Σs)
	#Γs_b = calc_Γ(r_mid_b, Σs_b)

end

# ╔═╡ cd8f63a6-bd53-4cf4-a3f4-a8943e0a1caf
scatter(filtered[:, :phot_g_mean_mag], r_ell, ma=0.1, ms=3)

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

# ╔═╡ 786f8e2c-3b12-47a9-b431-18e30cdf2a1e
idxs = sortperm(r_ell)

# ╔═╡ 38553483-6028-4d78-9805-80581e55e214
M = cumsum(ones(length(idxs)))

# ╔═╡ 20d997f1-72de-4c6b-9aa3-836c85221c21
p = Ref{Plots.Plot}()

# ╔═╡ b0de0957-acc9-49e8-a988-f531844cb9c0
md"""
# Gaussian Process regression
"""

# ╔═╡ abb50575-f5b9-43cd-8d73-92a7986fd9f0


# ╔═╡ 299c30ce-1f9c-48a6-b432-06c20aef5e64
@model function gp_model(log_r, y)
	jitter = 0.01

	v ~ Gamma(2, 1)
    l ~ Gamma(4, 1)
    f = GP(v * with_lengthscale(SEKernel(), l))
	f1 = to_sde(f, SArrayStorage(Float64))
    f_latent = f1(log_r, jitter)
	f_posterior = posterior(f_latent, y)

	σ ~ LogNormal(0, 0.3)
	μ = mean(f_posterior(log_r))
    y ~ product_distribution(Normal.(μ, σ))
end

# ╔═╡ a91c47f3-ca0c-4169-85a5-20aeb6dfb7dd
#optimize!(gp; method=ConjugateGradient())   # Optimise the hyperparameters

# ╔═╡ 2ba5de18-7155-4b88-9e99-a0646d97ebb4
begin 
	plot()
	rs, cdfs = plot_cdf!(r_ell[idxs], label="elliptical", norm=:count)
end

# ╔═╡ f3eb1aa2-4fa1-4917-9212-2dfc99a8033e
log_r = log10.(rs)

# ╔═╡ 3d589596-7f7d-4629-8508-ee15aec84ebf
y =  logit.(cdfs[1:end-1] ./ cdfs[end])

# ╔═╡ b7d0c9e4-e040-49e4-a1d0-b483277437d8
m = gp_model(log_r[1:end-1], y)

# ╔═╡ 0d34d28e-d924-4580-b788-6a29b3f68752
m()

# ╔═╡ acc84796-a9cc-436c-b399-8a99af3bae14
chn = sample(m, NUTS(), 100)

# ╔═╡ 1293ba28-c3f0-4d5c-bc55-9e5e35153acb
begin 
	# f1 = GP( 2*with_lengthscale(SEKernel(), 0.3))
	f_n = GP( 5with_lengthscale(Matern52Kernel(), 0.1))
	f1 = to_sde(f_n, SArrayStorage(Float64))
	fx = f1(log_r[1:end-1], 0.1)
	f_posterior = posterior(fx, y)
end

# ╔═╡ 1a869cdf-3a2c-4168-9e48-863614270617
rand(fx)

# ╔═╡ a7803020-e9ac-42cf-9e9e-622928624ce9
logpdf(fx, y)

# ╔═╡ 4b76891f-2896-4621-a5d3-d1c7cb7de6dd
begin 
	plot(f_posterior(log_r))
	scatter!(log_r, logit.(cdfs ./cdfs[end]))
end

# ╔═╡ bbaf8376-7fbc-4b3f-aef0-2c43fe693c75
md"""
# Spline thyme
"""

# ╔═╡ 89c145a0-3b75-4e7b-a8e8-a84ad339ddca
import Dierckx

# ╔═╡ 2ec6615a-116f-4e65-8630-e09cc21dcfec
begin 
	r_pred = 10 .^ LinRange(log_r[1] * 0.9, log_r[end] * 1.1, 10000)
	cdf_smoothed = logistic.(mean(f_posterior(log10.(r_pred))))
	cdf_s = Dierckx.Spline1D(log10.(r_pred), cdf_smoothed, k=5)
	#cdf_s = Spline1D(log10.(rs), cdfs, s=1e5)
end

# ╔═╡ c43fd6ac-a960-4bc0-8615-833392972852
# Example function to fit a spline
# x and y are the observed data
# n_knots is the number of knots you want to use for the spline
Turing.@model function spline_model(x, y, n_knots, k)
	# Define priors for spline parameters
	σ ~ Exponential(0.01)
	β ~ MvNormal(0.5 .+ zeros(n_knots), 0.5 * I)

	μ = spline_predict(x, α, β)
	# Likelihood of the observations
	for i in eachindex(y)
		y[i] ~ Normal(μ[i], σ)
	end
end


# ╔═╡ b1ee02d3-692c-4b4a-9c3b-499494829346
import Distributions: I

# ╔═╡ fd28d670-a7d3-488b-beea-3dbecb2f38f0
n_knots = 40

# ╔═╡ 59858f57-5389-4805-b2a4-c611fac6fa6d
k=4

# ╔═╡ e1151e36-0552-4f93-825d-70bb2889409c
α = LinRange(-2.5, 0.2, n_knots-2+2k)

# ╔═╡ c8bef74a-0402-4466-9e04-36c580f3437c
α_m = (α[k-1:end-k])

# ╔═╡ e9178188-0475-4dfc-8b75-d994506286ff
function spline_predict(x, α, β)
	basis = BSplineBasis(k, α)
	coeffs = vcat(fill(β[1], k), β, fill(β[end], k))
	spl = Spline(basis, coeffs)

	y = spl.(x)
	y = @. ifelse(x < α[1], β[1], y)
	y = @. ifelse(x > α[end], β[end], y) # continue BC's

	return y
end

# ╔═╡ 4adee558-cabb-45bf-8fb9-756d176d2058
s_model = spline_model(log10.(rs), cdfs ./ cdfs[end], n_knots, k)

# ╔═╡ ddbc552a-ff52-47c6-a168-d3f34e08fae4
s_map = optimize(s_model, MAP())

# ╔═╡ 3281b681-0c20-467e-a8da-9bffd9730a6d
s_chain = sample(s_model, NUTS(0.65), 100, initial_params=s_map.values.array)

# ╔═╡ 0b0e899b-acf2-4094-a8b7-61db3bb28e03
function N(i, k, x, T)
    if k == 0
        return T[i] ≤ x < T[i+1] ? 1.0 : 0.0
    else
        c1 = (x - T[i]) / (T[i+k] - T[i] + eps()) # Add eps() to avoid division by zero
        c2 = (T[i+k+1] - x) / (T[i+k+1] - T[i+1] + eps())
        return c1 * N(i, k-1, x, T) + c2 * N(i+1, k-1, x, T)
    end
end

# ╔═╡ e9ca6f15-06e1-4468-8769-090cb2f8a3eb
begin 
	p[] = plot()
	
	scatter!(log10.(rs), cdfs ./ cdfs[end])
	
	for i in 1:10:length(s_chain)
		#ac = [s_chain[Symbol("α[$j]")][i] for j in 1:n_knots]
		β1 = [s_chain[Symbol("β[$j]")][i] for j in 1:(n_knots)]

		ys = spline_predict(log10.(rs), α, β1)
		
		plot!(log10.(rs), ys, alpha=0.3, legend=false, color="black", label="" )
		scatter!(α_m, β1, alpha=0.2, mc=4)
	end
	bc1 = [s_map.values[Symbol("β[$j]")] for j in 1:(n_knots)]
	ys = spline_predict(log10.(rs), α, bc1)
	plot!(log10.(rs), ys, color=2, lw=0.5)
	p[]
end

# ╔═╡ 086e4d78-2cff-4c41-8e1a-00de3b4820da
begin 
	basis = BSplineBasis(k, α)
	coeffs = vcat(fill(bc1[1], k), bc1, fill(bc1[end], k))
	spl = Spline(basis, coeffs)
	
	cdf0(r) = spl(log10(r))
	cdf1(r) = spl(log10(r), Derivative(1))
	cdf2(r) = spl(log10(r), Derivative(2))

	cdf0(r) = cdf_s(log10(r))
	cdf1(r) = Dierckx.derivative(cdf_s, log10.(r))
	cdf2(r) = Dierckx.derivative(cdf_s, log10(r), nu=2)

end

# ╔═╡ 5dbc63c3-0ac3-4541-b4d7-0ff5724d311e
begin
	dM(r) = 1/r * cdf1(r) / log(10)
	d2M(r) = 1/r^2 * (cdf2(r) / log(10)^2 - cdf1(r) / log(10) )
	dΣ(r) = d2M(r) / (2π*r) - dM(r) / (2π*r ^ 2)
	Σ(r) = dM(r) / (2π * r)
end

# ╔═╡ 2d90b0fc-ece7-4073-96f5-c268c3ba8265
begin 
	plot(xscale=:log10, yscale=:log10, ylims=(1e-6, 1e3))
	scatter!(rs, Σ,)
	scatter!(r_mid, Σs ./cdfs[end], label="sculptor", msa=0)

end

# ╔═╡ 3a8f1bba-0fcc-4a45-92fd-c02d8fec4a1e
begin 
	plot(ylims=(-8, 2))

	scatter!(log10.(r_mid[2:end]), Γs, label="sculptor")
	plot!(log10.(rs), log10.(rs) * -0.7 .- 2)
	scatter!(log10.(rs), dΣ.(rs) .* rs ./ Σ.(rs) )

	xlabel!(L"$r_{\rm ell}$")
	ylabel!(L"$\Gamma = d\ \log(\Sigma) / d\ \log (r_{\rm ell})$")
end

# ╔═╡ da8aa554-c521-4d15-a59e-31082d90463f
begin 
	plot(ylabel="CDF residuals")
	scatter!(log10.(rs), (cdf0.(rs) .- cdfs ./ cdfs[end]), ms=5)
end

# ╔═╡ 4a5fb4c2-015f-401d-be4e-65eff647e009
scatter(α_m[2:end], bc1, color="green")

# ╔═╡ 3db5cd7e-d039-4f10-a466-7da84e1d2c5b
plot(bc1)

# ╔═╡ 775e37c5-e1fc-45bf-9337-8bba4d94352a
bc1

# ╔═╡ 21ab38bd-05bb-4b75-9bf2-cd416202d391
spline_predict(365, α, bc1)

# ╔═╡ 4ef09052-4547-4688-a8a7-afd12c76fbdb
[α[end]]

# ╔═╡ Cell order:
# ╠═98c47808-7b07-4438-a855-d3f0167f3f34
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═40eb3191-3c32-4928-9344-05b48da4bdd0
# ╠═d31ccfbf-6628-4745-a0ca-fd2061821416
# ╠═02b0d62c-2a47-4308-ab96-ede30029b8c8
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═f6a5c138-10bc-48e2-aec7-45fd61b9b17f
# ╠═6ca5439b-b4d7-44d5-b7bd-303f14ee7c15
# ╠═d04f8836-bf9a-4a7b-8ec9-a1ba783c7295
# ╠═32293baf-0c39-488d-b495-a86f42eed178
# ╠═ba71616e-dadf-4025-9afd-66dc40d4e65b
# ╟─efe133bc-670a-49a5-9f22-c937456a11b5
# ╠═19344fb8-7eb4-4f56-af6a-14441226cabd
# ╠═c70fce56-9367-4f6a-92bb-f0d0d7a616e0
# ╠═0bf61650-adae-43e0-8caf-02951aa6082c
# ╠═695ac014-86ee-48cc-a6e6-6d8a0c444be1
# ╠═cd8f63a6-bd53-4cf4-a3f4-a8943e0a1caf
# ╟─802348e2-8e11-42f5-b2e8-4a66a39936a2
# ╠═54477bfa-bfc0-4d44-9f03-72531ec3de7d
# ╠═786f8e2c-3b12-47a9-b431-18e30cdf2a1e
# ╠═38553483-6028-4d78-9805-80581e55e214
# ╠═20d997f1-72de-4c6b-9aa3-836c85221c21
# ╠═b0de0957-acc9-49e8-a988-f531844cb9c0
# ╠═4cc65ca0-7a21-4dc0-b8f2-df1cdd55f915
# ╠═abb50575-f5b9-43cd-8d73-92a7986fd9f0
# ╠═299c30ce-1f9c-48a6-b432-06c20aef5e64
# ╠═b7d0c9e4-e040-49e4-a1d0-b483277437d8
# ╠═0d34d28e-d924-4580-b788-6a29b3f68752
# ╠═acc84796-a9cc-436c-b399-8a99af3bae14
# ╠═f3eb1aa2-4fa1-4917-9212-2dfc99a8033e
# ╠═1293ba28-c3f0-4d5c-bc55-9e5e35153acb
# ╠═3d589596-7f7d-4629-8508-ee15aec84ebf
# ╠═a7803020-e9ac-42cf-9e9e-622928624ce9
# ╠═4b76891f-2896-4621-a5d3-d1c7cb7de6dd
# ╠═2ec6615a-116f-4e65-8630-e09cc21dcfec
# ╠═1a869cdf-3a2c-4168-9e48-863614270617
# ╠═086e4d78-2cff-4c41-8e1a-00de3b4820da
# ╠═5dbc63c3-0ac3-4541-b4d7-0ff5724d311e
# ╠═a91c47f3-ca0c-4169-85a5-20aeb6dfb7dd
# ╠═da8aa554-c521-4d15-a59e-31082d90463f
# ╠═2d90b0fc-ece7-4073-96f5-c268c3ba8265
# ╠═2ba5de18-7155-4b88-9e99-a0646d97ebb4
# ╠═3a8f1bba-0fcc-4a45-92fd-c02d8fec4a1e
# ╠═bbaf8376-7fbc-4b3f-aef0-2c43fe693c75
# ╠═89c145a0-3b75-4e7b-a8e8-a84ad339ddca
# ╠═e1151e36-0552-4f93-825d-70bb2889409c
# ╠═c8bef74a-0402-4466-9e04-36c580f3437c
# ╠═e9178188-0475-4dfc-8b75-d994506286ff
# ╠═c43fd6ac-a960-4bc0-8615-833392972852
# ╠═b1ee02d3-692c-4b4a-9c3b-499494829346
# ╠═fd28d670-a7d3-488b-beea-3dbecb2f38f0
# ╠═59858f57-5389-4805-b2a4-c611fac6fa6d
# ╠═4adee558-cabb-45bf-8fb9-756d176d2058
# ╠═ddbc552a-ff52-47c6-a168-d3f34e08fae4
# ╠═3281b681-0c20-467e-a8da-9bffd9730a6d
# ╠═0b0e899b-acf2-4094-a8b7-61db3bb28e03
# ╠═e9ca6f15-06e1-4468-8769-090cb2f8a3eb
# ╠═4a5fb4c2-015f-401d-be4e-65eff647e009
# ╠═3db5cd7e-d039-4f10-a466-7da84e1d2c5b
# ╠═775e37c5-e1fc-45bf-9337-8bba4d94352a
# ╠═21ab38bd-05bb-4b75-9bf2-cd416202d391
# ╠═4ef09052-4547-4688-a8a7-afd12c76fbdb
