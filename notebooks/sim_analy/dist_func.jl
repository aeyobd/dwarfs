### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 963a0312-6eec-11ef-3070-53f3c18efb6f
begin 
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya
	using LilGuys
end

# ╔═╡ 9e3600af-5954-4ec1-ad0d-b86e83e63d79
LP = LilGuys.Plots

# ╔═╡ 0febc738-cfe8-47e8-8acb-7e09a6f4191c
import DensityEstimators as DE

# ╔═╡ 879566d4-c24c-415a-9781-4ec0a4413680
haloname = "/astro/dboyea/sculptor/isolation/1e6/fiducial/stars"

# ╔═╡ ee946594-16da-4edf-a884-967f8ae9ff53
snapname = joinpath(haloname, "../halos/fiducial.hdf5")

# ╔═╡ f31b6dab-bec1-4716-bea1-79507cf6fa03
snap = Snapshot(snapname)

# ╔═╡ ffaf29c3-8861-452e-8a42-d548dbcb4a54
halo = LilGuys.load_profile(joinpath(haloname, "../halo.toml"))

# ╔═╡ fbef83a8-8224-4a66-9dd1-c04ca37b3001
radii = LilGuys.calc_r(snap)

# ╔═╡ 9a5d4c4d-7128-4ad3-ad4d-f41ad51669c9
begin 
	r_bins = 10 .^ LinRange(minimum(log10.(radii)), maximum(log10.(radii)), 100)
	r_bins = r_bins[r_bins .< 1]
end

# ╔═╡ 77abd400-ea8c-4728-bdb0-9e86435608f7
minimum(radii)

# ╔═╡ 23c77c0c-c2b0-4245-b48d-638358eac6e9
r_m = midpoints(r_bins)

# ╔═╡ cc75d352-83b8-4fe2-9b9a-173912e97086
logr_m = log10.(r_m)

# ╔═╡ b7c9a69f-6135-41ec-8c57-527a4be335ab
h = DE.histogram(radii, r_bins)

# ╔═╡ b86eaab9-610c-4fd3-82d1-d9f1460d211a
plot(logr_m, log10.(h.values))

# ╔═╡ 992c9e66-c53e-441c-ade8-dcade5e61f52
h.values

# ╔═╡ a7d46b4a-7a79-4320-9968-9e89365b59a2
ρ = LilGuys.calc_ρ_from_hist(r_bins, h.values)

# ╔═╡ ea8ecc31-1e44-4cff-84f8-61dbdcea687e
Φ = LilGuys.calc_radial_discrete_Φ(snap)

# ╔═╡ 12d93b96-135c-4ee5-a226-fdfc64b9e916
ϕ = LilGuys.lerp(sort(radii), -Φ[sortperm(radii)]).(r_m)

# ╔═╡ e91afedb-9fcf-4303-ab27-a6ab53c2ce17
let
	fig, ax = FigAxis(
		xlabel = LP.log_r_label,
		ylabel = L"\psi"
	)


	lines!(logr_m, ϕ)

	fig
end

# ╔═╡ b180ed7f-ef91-4d62-91de-e38b50153b83
DF = LilGuys.DistributionFunction(ρ, ϕ, r_m)

# ╔═╡ 912db6da-822e-4b63-82bc-4f3329e3f89c
let
	fig = Figure(size=(400, 800))
	
	ax = Axis(fig[1,1],
		ylabel = LP.log_rho_label,
	)

	errscatter!(logr_m, log10.(ρ), yerr=1 / log(10) ./ sqrt.(h.values))

	x = r_m 
	y = ρ
	ax2 = Axis(fig[2, 1],
		ylabel = L"\textrm{asinh}\; d^2\rho/d\,r^2"
	)
	
	scatter!(logr_m, asinh.(LilGuys.gradient(LilGuys.gradient(y, x), x)))


	ax3 = Axis(fig[3, 1],
		xlabel = LP.log_r_label,
		ylabel = L"\textrm{asinh}\;d^2\rho/d\,\psi^2"
	)

	
	scatter!(logr_m, asinh.(DF.d2ρ_dψ2.(ϕ)))

	linkxaxes!(ax, ax2, ax3)
	hidexdecorations!(ax, grid=false)
	hidexdecorations!(ax2, grid=false)

	fig
end

# ╔═╡ b6f9bb10-105b-48b1-9f5d-447112a86ed5
ϵ = LinRange(ϕ[1], ϕ[end], 1000)

# ╔═╡ 30fe63ee-8478-45e9-9789-b73902c4074d
f_ϵ = DF.(ϵ)

# ╔═╡ ae225d08-7333-4a01-86d5-a228785ff51f
let
	fig, ax = FigAxis(
		xlabel = "energy",
	)

	scatter!(ϵ, asinh.(DF.d2ρ_dψ2.(ϵ)))

	dρ_dψ = LilGuys.gradient(ρ, ϕ)
	scatter!(ϕ, asinh.(dρ_dψ))
	scatter!(ϕ, asinh.(ρ))

	d2ρ_dψ2 = LilGuys.gradient(LilGuys.gradient(ρ, ϕ), ϕ)

	scatter!(ϕ, asinh.(d2ρ_dψ2))
	fig
end


# ╔═╡ cf3320c2-1192-4e76-9dcd-8782563f412e
let
	fig, ax = FigAxis(
		xlabel = "ψ",
		limits=(0.08, 0.099, 8, 12),
		ylabel="ρ"
	)

	scatter!(ϕ, asinh.(ρ))

	fig
end


# ╔═╡ 194a76f7-4e47-48a7-82da-aed603ba923d
let
	x = r_m 
	y = ϕ

	fig, ax = FigAxis(
		xlabel = LP.log_r_label,
		title="spherical approx"
	)
	scatter!(logr_m, y, label = L"\psi")
	scatter!(logr_m, LilGuys.gradient(y, x), label = L"\psi'(r)")
	scatter!(logr_m, LilGuys.gradient(LilGuys.gradient(y, x), x), label = L"\psi''(r)")


	logx = LinRange(logr_m[1], logr_m[end], 1000)
	x = 10 .^ logx 
	y = -LilGuys.calc_Φ.(halo, x)

	lines!(logx, y, label="analytic")
	lines!(logx, LilGuys.gradient(y, x))
	lines!(logx, LilGuys.gradient(LilGuys.gradient(y, x), x))

	axislegend()
	fig
end

# ╔═╡ 3ab29510-392e-447a-a742-797e219d408f
let
	x = r_m 
	y = ρ

	fig, ax = FigAxis(
		xlabel = LP.log_r_label,
		title="density and derivatives"
	)
	scatter!(logr_m, log10.(y), label = L"\log\rho")
	scatter!(logr_m, asinh.(LilGuys.gradient(y, x)), label = L"\textrm{arsinh}\;\rho'(r)")
	scatter!(logr_m, asinh.(LilGuys.gradient(LilGuys.gradient(y, x), x)), label = L"\textrm{arsinh}\;\rho''(r)")


	logx = LinRange(logr_m[1], logr_m[end], 1000)
	x = 10 .^ logx 
	y = LilGuys.calc_ρ.(halo, x) * length(snap.masses)

	lines!(logx, log10.(y), label="analytic")
	lines!(logx, asinh.(LilGuys.gradient(y, x)))
	lines!(logx, asinh.(LilGuys.gradient(LilGuys.gradient(y, x), x)))

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ f0eae949-5b4f-4a23-adf5-2a583d2c9e75
let
	x = r_m 
	
	y = LilGuys.lerp(sort(radii), snap.Φs[sortperm(radii)]).(x)

	fig, ax = FigAxis(
		xlabel = LP.log_r_label,
		title="interpolating snapshot potential"
	)
	scatter!(logr_m, y, label = L"\psi")
	scatter!(logr_m, LilGuys.gradient(y, x), label = L"\psi'(r)")
	scatter!(logr_m, LilGuys.gradient(LilGuys.gradient(y, x), x), label = L"\psi''(r)")


	logx = LinRange(logr_m[1], logr_m[end], 1000)
	x = 10 .^ logx 
	y = -LilGuys.calc_Φ.(halo, x)

	lines!(logx, y, label="analytic")
	lines!(logx, LilGuys.gradient(y, x))
	lines!(logx, LilGuys.gradient(LilGuys.gradient(y, x), x))

	axislegend()
	fig
end

# ╔═╡ Cell order:
# ╠═963a0312-6eec-11ef-3070-53f3c18efb6f
# ╠═9e3600af-5954-4ec1-ad0d-b86e83e63d79
# ╠═0febc738-cfe8-47e8-8acb-7e09a6f4191c
# ╠═879566d4-c24c-415a-9781-4ec0a4413680
# ╠═ee946594-16da-4edf-a884-967f8ae9ff53
# ╠═f31b6dab-bec1-4716-bea1-79507cf6fa03
# ╠═ffaf29c3-8861-452e-8a42-d548dbcb4a54
# ╠═fbef83a8-8224-4a66-9dd1-c04ca37b3001
# ╠═9a5d4c4d-7128-4ad3-ad4d-f41ad51669c9
# ╠═77abd400-ea8c-4728-bdb0-9e86435608f7
# ╠═23c77c0c-c2b0-4245-b48d-638358eac6e9
# ╠═cc75d352-83b8-4fe2-9b9a-173912e97086
# ╠═b7c9a69f-6135-41ec-8c57-527a4be335ab
# ╠═b86eaab9-610c-4fd3-82d1-d9f1460d211a
# ╠═992c9e66-c53e-441c-ade8-dcade5e61f52
# ╠═a7d46b4a-7a79-4320-9968-9e89365b59a2
# ╠═912db6da-822e-4b63-82bc-4f3329e3f89c
# ╠═ea8ecc31-1e44-4cff-84f8-61dbdcea687e
# ╠═12d93b96-135c-4ee5-a226-fdfc64b9e916
# ╠═e91afedb-9fcf-4303-ab27-a6ab53c2ce17
# ╠═b180ed7f-ef91-4d62-91de-e38b50153b83
# ╠═b6f9bb10-105b-48b1-9f5d-447112a86ed5
# ╠═30fe63ee-8478-45e9-9789-b73902c4074d
# ╠═ae225d08-7333-4a01-86d5-a228785ff51f
# ╠═cf3320c2-1192-4e76-9dcd-8782563f412e
# ╠═194a76f7-4e47-48a7-82da-aed603ba923d
# ╠═3ab29510-392e-447a-a742-797e219d408f
# ╠═f0eae949-5b4f-4a23-adf5-2a583d2c9e75
