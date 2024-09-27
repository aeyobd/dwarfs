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

# ╔═╡ 8ec7fdd1-6805-48aa-b700-9e6f4da765a5
md"""
some plots and rough calculations for the distribution function
"""

# ╔═╡ 9e3600af-5954-4ec1-ad0d-b86e83e63d79
LP = LilGuys.Plots

# ╔═╡ 0febc738-cfe8-47e8-8acb-7e09a6f4191c
import DensityEstimators as DE

# ╔═╡ 879566d4-c24c-415a-9781-4ec0a4413680
haloname = "/astro/dboyea/sculptor/isolation/1e6/s0.14/stars"

# ╔═╡ ee946594-16da-4edf-a884-967f8ae9ff53
snapname = joinpath(haloname, "../out/1")

# ╔═╡ 1455ca37-525e-46bf-9fd5-b99a058db1f3
md"""
# File loading
"""

# ╔═╡ ffaf29c3-8861-452e-8a42-d548dbcb4a54
halo = LilGuys.load_profile(joinpath(haloname, "../halo.toml"))

# ╔═╡ f31b6dab-bec1-4716-bea1-79507cf6fa03
snap = Snapshot(snapname)

# ╔═╡ de1519ce-22ae-4f3e-bd50-66b6b891d2c3
df_emp = LilGuys.read_hdf5_table(haloname * "/distribution_function.hdf5")

# ╔═╡ 617c927a-8204-4a84-8ff9-531fab8b9150
df_df = LilGuys.read_hdf5_table(haloname * "/../stars/distribution_function.hdf5")

# ╔═╡ 14100a85-74d9-46cf-9f8c-4cba33b9a6b0
df_E = LilGuys.read_hdf5_table(haloname * "/energies.hdf5")

# ╔═╡ 2a025f1a-ac58-4ec6-8066-069030b0a7f3
md"""
# Calculations
"""

# ╔═╡ fbef83a8-8224-4a66-9dd1-c04ca37b3001
radii = calc_r(snap)

# ╔═╡ 9a5d4c4d-7128-4ad3-ad4d-f41ad51669c9
begin 
	r_bins = 10 .^ LinRange(minimum(log10.(radii)), maximum(log10.(radii)), 80)
	#r_bins = LilGuys.quantile(radii, LinRange(0, 1, 200))
end

# ╔═╡ 23c77c0c-c2b0-4245-b48d-638358eac6e9
r_m = midpoints(r_bins)

# ╔═╡ cc75d352-83b8-4fe2-9b9a-173912e97086
logr_m = log10.(r_m)

# ╔═╡ b7c9a69f-6135-41ec-8c57-527a4be335ab
h = DE.histogram(radii, r_bins, weights=snap.masses)

# ╔═╡ a7d46b4a-7a79-4320-9968-9e89365b59a2
ρ = LilGuys.calc_ρ_from_hist(r_bins, h.values)

# ╔═╡ ea8ecc31-1e44-4cff-84f8-61dbdcea687e
Φ = LilGuys.calc_radial_discrete_Φ(snap)

# ╔═╡ 12d93b96-135c-4ee5-a226-fdfc64b9e916
ψ = LilGuys.lerp(sort(radii), -Φ[sortperm(radii)]).(r_m)

# ╔═╡ fa8d3a08-162e-4693-8274-4bdb0cb8964a
issorted(Φ[sortperm(LilGuys.calc_r(snap))])

# ╔═╡ 8993a8fc-ebc0-4d9b-80dc-f8db05c46269
Φ2 = LilGuys.calc_radial_discrete_Φ(calc_r(snap), snap.masses)

# ╔═╡ e91afedb-9fcf-4303-ab27-a6ab53c2ce17
let
	fig, ax = FigAxis(
		xlabel = LP.log_r_label,
		ylabel = L"\psi"
	)


	# lines!(logr_m, ϕ)
	scatter!(log10.(df_emp.radii), df_emp.psi, label="emperical", color=COLORS[3], markersize=3)

	skip = 100
	phi = -snap.Φs[1:skip:end]
	radii = calc_r(snap)[1:skip:end]

	scatter!(log10.(radii), phi, label="snapshot", markersize=3, color=COLORS[2])
	lines!(log10.(df_df.radii), df_df.psi, label="analytic")

	axislegend(position=:rb)

	fig
end

# ╔═╡ b180ed7f-ef91-4d62-91de-e38b50153b83
DF = LilGuys.DistributionFunction(ρ, ψ, r_m, force_positive=true)

# ╔═╡ 912db6da-822e-4b63-82bc-4f3329e3f89c
let
	fig = Figure(size=(400, 800))
	
	ax = Axis(fig[1,1],
		ylabel = LP.log_rho_label,
		limits=((-2, 2), (-15, 10))
	)

	#errscatter!(logr_m, log10.(ρ), yerr=1 / log(10) ./ sqrt.(h.values))
	lines!(log10.(df_df.radii), log10.(df_df.rho))
	scatter!(log10.(df_emp.radii), log10.(df_emp.rho))
	scatter!(logr_m, log10.(ρ))
	
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

	
	scatter!(logr_m, asinh.(DF.d2ρ_dψ2.(ψ)))

	linkxaxes!(ax, ax2, ax3)
	hidexdecorations!(ax, grid=false)
	hidexdecorations!(ax2, grid=false)

	fig
end

# ╔═╡ 354acb89-2e84-46f7-8ac0-d1283abd2061
DF_ana = LilGuys.DistributionFunction(calc_ρ.(halo, r_m), -LilGuys.calc_Φ.(halo, r_m), r_m)

# ╔═╡ b6f9bb10-105b-48b1-9f5d-447112a86ed5
ϵ = LinRange(minimum(ψ), maximum(ψ), 1000)

# ╔═╡ 30fe63ee-8478-45e9-9789-b73902c4074d
f_ϵ = DF.(ϵ)

# ╔═╡ 25257acd-7aac-4ed1-b5af-865fa295d97b
let
	fig, ax = FigAxis(
		xlabel="energy", ylabel="asinh distribution function",
		limits=((0.08, 0.102), nothing)
	)

	
	scatter!(ϵ, asinh.(f_ϵ), label="this notebook")

	e = -LilGuys.calc_Φ.(halo, r_m)
	lines!(e, asinh.(DF_ana.(e)), label="analytic")
	scatter!(df_emp.psi, asinh.(df_emp.f), label="emperical")

	lines!(df_df.psi, asinh.(df_df.f), color=COLORS[3], label="analytic 2")

	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ ae225d08-7333-4a01-86d5-a228785ff51f
let
	fig, ax = FigAxis(
		xlabel = "energy",
	)

	scatter!(ϵ, asinh.(DF.d2ρ_dψ2.(ϵ)))

	dρ_dψ = LilGuys.gradient(ρ, ψ)
	scatter!(ψ, asinh.(dρ_dψ))
	scatter!(ψ, asinh.(ρ))

	d2ρ_dψ2 = LilGuys.gradient(LilGuys.gradient(ρ, ψ), ψ)

	scatter!(ψ, asinh.(d2ρ_dψ2))
	fig
end


# ╔═╡ cf3320c2-1192-4e76-9dcd-8782563f412e
let
	fig, ax = FigAxis(
		xlabel = "ψ",
		limits=(0.08, 0.099, 8, 12),
		ylabel="ρ"
	)

	scatter!(ψ, asinh.(ρ))

	fig
end


# ╔═╡ 3ab29510-392e-447a-a742-797e219d408f
let
	h = 1e1
	x = r_m 
	y = ρ ./ h
	logx = LinRange(logr_m[1], logr_m[end], 1000)
	x_m = 10 .^ logx 
	y_m = LilGuys.calc_ρ.(halo, x_m) ./ h

	fig = Figure(size=(400, 800))
	
	ax = Axis(fig[1, 1],
		title="density and derivatives",
		limits=(nothing, (-10, 3))
	)
	scatter!(logr_m, log10.(y), label = L"\log\rho")
	lines!(logx, log10.(y_m), label="analytic")



	ax_2 = Axis(fig[2, 1], ylabel=L"\textrm{arsinh}\;\rho'(r)")
	scatter!(logr_m, asinh.(LilGuys.gradient(y, x)))
	lines!(logx, asinh.(LilGuys.gradient(y_m, x_m)))

	ax_3 = Axis(fig[3, 1],
			xlabel = LP.log_r_label,
		ylabel=L"\textrm{arsinh}\;\rho''(r)",
		limits=((-2, 3), nothing)
	)

	
	scatter!(logr_m, asinh.(LilGuys.gradient(LilGuys.gradient(y, x), x)))
	lines!(logx, asinh.(LilGuys.gradient(LilGuys.gradient(y_m, x_m), x_m)))

	linkxaxes!(ax, ax_2, ax_3)
	
	fig
end

# ╔═╡ 194a76f7-4e47-48a7-82da-aed603ba923d
let
	x = r_m 
	y = ψ

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

# ╔═╡ f0eae949-5b4f-4a23-adf5-2a583d2c9e75
let
	x = r_m 
	
	y = -LilGuys.lerp(sort(radii), snap.Φs[sortperm(radii)]).(x)

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
# ╠═8ec7fdd1-6805-48aa-b700-9e6f4da765a5
# ╠═963a0312-6eec-11ef-3070-53f3c18efb6f
# ╠═9e3600af-5954-4ec1-ad0d-b86e83e63d79
# ╠═0febc738-cfe8-47e8-8acb-7e09a6f4191c
# ╠═879566d4-c24c-415a-9781-4ec0a4413680
# ╠═ee946594-16da-4edf-a884-967f8ae9ff53
# ╟─1455ca37-525e-46bf-9fd5-b99a058db1f3
# ╠═ffaf29c3-8861-452e-8a42-d548dbcb4a54
# ╠═f31b6dab-bec1-4716-bea1-79507cf6fa03
# ╠═de1519ce-22ae-4f3e-bd50-66b6b891d2c3
# ╠═617c927a-8204-4a84-8ff9-531fab8b9150
# ╠═14100a85-74d9-46cf-9f8c-4cba33b9a6b0
# ╠═2a025f1a-ac58-4ec6-8066-069030b0a7f3
# ╠═fbef83a8-8224-4a66-9dd1-c04ca37b3001
# ╠═9a5d4c4d-7128-4ad3-ad4d-f41ad51669c9
# ╠═23c77c0c-c2b0-4245-b48d-638358eac6e9
# ╠═cc75d352-83b8-4fe2-9b9a-173912e97086
# ╠═b7c9a69f-6135-41ec-8c57-527a4be335ab
# ╠═a7d46b4a-7a79-4320-9968-9e89365b59a2
# ╠═ea8ecc31-1e44-4cff-84f8-61dbdcea687e
# ╠═12d93b96-135c-4ee5-a226-fdfc64b9e916
# ╠═912db6da-822e-4b63-82bc-4f3329e3f89c
# ╠═fa8d3a08-162e-4693-8274-4bdb0cb8964a
# ╠═8993a8fc-ebc0-4d9b-80dc-f8db05c46269
# ╠═e91afedb-9fcf-4303-ab27-a6ab53c2ce17
# ╠═b180ed7f-ef91-4d62-91de-e38b50153b83
# ╠═354acb89-2e84-46f7-8ac0-d1283abd2061
# ╠═b6f9bb10-105b-48b1-9f5d-447112a86ed5
# ╠═25257acd-7aac-4ed1-b5af-865fa295d97b
# ╠═30fe63ee-8478-45e9-9789-b73902c4074d
# ╠═ae225d08-7333-4a01-86d5-a228785ff51f
# ╠═cf3320c2-1192-4e76-9dcd-8782563f412e
# ╠═3ab29510-392e-447a-a742-797e219d408f
# ╠═194a76f7-4e47-48a7-82da-aed603ba923d
# ╠═f0eae949-5b4f-4a23-adf5-2a583d2c9e75
