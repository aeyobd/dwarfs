### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	import LilGuys as lguys
	using Plots; #plotlyjs()
	using CSV, DataFrames
	using LaTeXStrings
end

# ╔═╡ 07cf49a5-909a-48dc-982b-6c6953e5e414
using NaNMath; nm=NaNMath

# ╔═╡ 055e698b-9092-47e4-a70f-a03d5db71bb4
begin 
	import Dierckx
	import QuadGK: quadgk
end

# ╔═╡ 6964cbd2-9ef9-4682-9ac1-80c10518c374
snap = lguys.Snapshot("out/snapshot_001.hdf5")

# ╔═╡ 79ce82c5-47d8-4157-942a-ff8e524edfd4
cen = lguys.ss_centre(snap)

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
begin
	snap_i = lguys.copy(snap)
	snap_i.positions .-= cen.x_c
	snap_i.velocities .-= cen.v_c


	radii = lguys.calc_r(snap_i.positions)
	vels = lguys.calc_r(snap_i.velocities)
	Φ_f = lguys.calc_radial_Φ(radii, snap_i.masses)

	ke = @. 1/2 * vels^2
	ϵ = @. -Φ_f(radii) - ke
	filt = ϵ .> 0
	
	p = sortperm(radii)
	idx = p[filt[p]]

	snap_i = snap_i[idx]
	ke = ke[idx]
	radii = radii[idx]
	vels = vels[idx]
	ϵ = ϵ[idx]

	Min = cumsum(snap_i.masses)
end

# ╔═╡ 4ae5506f-0e5e-44ac-b9e4-161ca7b448af
lguys.scatter_xyz(snap_i.positions, ms=1)[1]

# ╔═╡ 5214efbc-4445-401f-9ef0-3d4c777fd1bd
md"""
# Density Profile methods
"""

# ╔═╡ c7aa3122-4bdd-4889-895a-d143f8919405
prof = lguys.Profile(snap_i, zeros(3), zeros(3), Nr=200)

# ╔═╡ d8c66c89-8e79-4214-82d1-05cfe248ab8a
begin 
	plot(xscale=:log10)
	scatter!(radii, Min)
	sp = Dierckx.Spline1D(radii, Min, s=0.0001)
	plot!(radii, sp.(radii))
end

# ╔═╡ 9b9f0ee6-0856-494d-8135-98ce77995d42
begin
	plot(xscale=:log10)
	#scatter!(rs_mid, Ψ, label="gravity (Ψ)")
	#plot!(rs_mid, Ψs, label="gravity (Ψ)")
	scatter!(radii, -snap_i.Φs, label="Ψ")

	scatter!(radii, 1/2 * vels .^ 2, label="kinetic")
	#scatter!(radii, ϵ, label="binding (ϵ)")
	xlabel!("r")
	ylabel!("specific energy")
end

# ╔═╡ 32411e1f-da6f-405c-9b6c-20c4a2dbf4a5
rapha = CSV.read("DFDM.csv", DataFrame)

# ╔═╡ 80b32d08-405f-4e1a-a6a0-3325afc1e283
md"""
# Distribution functions
The smoothed CDF (cumulative mass function), which everything follows from
"""

# ╔═╡ 8938d5ce-adc9-446c-8cc3-4b637d81e480
length(prof.Ms)

# ╔═╡ bd21d3e8-5a5c-40ee-9f16-3e429492d4d0
begin 
	log_radii = log10.(radii)
	rs = 10 .^ LinRange(log_radii[1], log_radii[end], 1000)
	rs = prof.rs
	r_end = rs[end]
	rs_mid = lguys.midpoint(rs)
	log_rs = log10.(rs)

	Mr = Dierckx.Spline1D(rs, prof.Ms)
	M1(r) = Dierckx.derivative(Mr, r)

	log_ρ = Dierckx.Spline1D(log10.(rs), log10.(prof.ρs), s=0.0)
	log_ρ1(log_r) = Dierckx.derivative(log_ρ, log_r)
	log_ρ2(log_r) = Dierckx.derivative(log_ρ, log_r, nu=2)
	ρ = Dierckx.Spline1D((rs), (prof.ρs), s=0.0)
	ρ1(r) =  Dierckx.derivative(ρ, r, nu=1)
	ρ2(r) = Dierckx.derivative(ρ, r, nu=2)

	ms = diff(prof.Ms)
end

# ╔═╡ 2e8de7c6-83e8-4a84-a7a4-4910401b0f4e
begin
	plot()
	#plot(log10.(rs_mid), log10.(ρ.(rs_mid)), label="GP")
	plot!(log10.(prof.rs), log10.(prof.ρs), label="histogram method")
	plot!(log10.(rs), log10.(ρ.(rs)))
end

# ╔═╡ 8db8fdc3-5b85-44ba-a8d9-99d197117536
begin 
	plot(xlabel="log radius", ylabel="log M(r)")
	scatter!(log10.(rs), log10.(Mr(rs)))
end

# ╔═╡ ecd35011-2f81-4a10-b931-250a2e0b9b2c
begin
	# takes sorted rs and cumulative masses Ms
	#function calc_f_E(M, r_end = M.t[end])
	# M1(r) = 
	# M2(r) = Dierckx.derivative(Mr, r, nu=2)
	# M3(r) = Dierckx.derivative(Mr, r, nu=3)
	# ρ2(r) = @. M3(r) / (4π * r^2) - M2(r) / (π * r^3) + M1(r) * 3/(2π * r^4) 
	# ρ(r) = @. M1(r) / (4π * r^2)
	# ρ1(r) = @. M2(r) / (4π * r^2) - M1(r) / (2π * r^3)
	
	end

# ╔═╡ aad56d80-3e06-4e98-b7c9-cd013482dfdf
begin 
	Ψs = -lguys.calc_radial_discrete_Φ(ms, rs_mid)
	Ψ = Dierckx.Spline1D(rs_mid, Ψs)
	Ψ_inv = Dierckx.Spline1D(reverse(Ψs), reverse(rs_mid))# reverse so xs increase

	Ψ1(r) = @. -lguys.G * Mr(r) / r^2
	Ψ2(r) = @. -lguys.G / r^3 * (r*M1(r) - 2*Mr(r))


	d2ρ_dΨ2(r) = ρ2(r) / Ψ1(r)^2 - ρ1(r) * Ψ2(r) / Ψ1(r)^3 

	function f_integrand(r, E)
		if E - Ψ(r) > 0
			return d2ρ_dΨ2(r) * Ψ1(r) / sqrt(E - Ψ(r)) 
		else
			return 0
		end
	end
	
	f(E) =  -1/(√8 * π^2) * (
		quadgk(r->f_integrand(r, E), Ψ_inv(E) * (1+1e-5), r_end)[1]
		+ 0#+ 1/sqrt(E) * ρ1(0) / Ψ1(0))
	)
end

# ╔═╡ a7646485-1dc9-41c0-96d3-5e0a740be11e
begin 
	# this is rapha's method more exactly
	ρ0 = ms ./ (4π/3 * diff(rs .^ 3))
	Φ0 = lguys.calc_radial_Φ(rs_mid, ms)
	Ψ0 = -Φ0.(rs_mid)

	dρ1 = lguys.gradient(Ψ0, ρ0)
	dρ2 = lguys.gradient(Ψ0, dρ1)

	dρ2_i = Dierckx.Spline1D(reverse(Ψ0), reverse(dρ2))

	f_2(E) =  1/(√8 * π^2) * quadgk(p -> dρ2_i(p)/ sqrt(E - p), Ψs[end], E*(1-1e-6))[1]

end

# ╔═╡ 6c35dc16-643d-47bc-9b7e-8f3069e0188d
f_2(100)

# ╔═╡ b1e5f70f-d724-40ec-a30b-df169ffd1e8b
Es = LinRange(0.01, (0.99maximum(Ψs)), 300)

# ╔═╡ 26cec8e1-bca0-4642-b98a-79ad2d97b076
scatter(Es, f_2.(Es))

# ╔═╡ 1ce021f2-e732-4acf-ae1d-b605513a019f
plot(rs, d2ρ_dΨ2.(rs), xscale=:log10)

# ╔═╡ df571e5d-2b8c-4476-b459-1c725062ac8e
begin
	plot(xlabel="ϵ", ylabel="DF")
	#plot!(Es, f.(Es))
	plot!(Es,f_2.(Es))
	#plot!(rapha.E / rapha.E[end-10] * Es[end], 0.7 * rapha.DFDM / maximum(rapha.DFDM[1:end-1]))
	#vline!([Es[end]])
	

end

# ╔═╡ 709b2c42-550c-4b82-88ae-6d940055e13f
begin 
	plot(xlabel=L"$\epsilon$", ylabel="counts")
	histogram!(ϵ)
end

# ╔═╡ 2677f891-b5a9-4c29-9aa0-0f305af61f95
begin 
	plot()
	plot!(Ψ.(rs), d2ρ_dΨ2.(rs),)
	plot!(Ψs, lguys.gradient(Ψs, lguys.gradient(Ψs, ρ.(rs_mid))), label="gradient")
	plot!(Ψs, dρ2_i.(Ψs), label="rapha method")

end

# ╔═╡ 2028a467-6a81-46e9-a438-b73f884c0d20
begin 
	plot(ylims=(-2e3, 1e3), xscale=:log10)
	plot!((rs), d2ρ_dΨ2.(rs),)
	plot!(rs_mid, lguys.gradient(Ψ0, lguys.gradient(Ψ0, ρ.(rs_mid))), label="gradient")

end

# ╔═╡ aba99571-2013-49e3-af86-eb59a0a38109
begin 
	plot(ylim=(-1, 1))
	plot!(Ψ.(rs), ρ1(rs) ./ Ψ1.(rs) )
	plot!(Ψs, lguys.gradient(Ψs, ρ.(rs_mid)), label="gradient")

end

# ╔═╡ 191fd0b0-f8d3-400b-afdb-6f1022854f35
Ψs

# ╔═╡ 3ef7e24e-ed44-445a-870b-e91e6ce07c65
Ψs[end]

# ╔═╡ 97baf791-8f97-4dfb-94ca-d1496162f688
begin 
	plot(radii, -snap_i.Φs)
	plot!(radii, Ψ.(radii))
	plot!( Ψ_inv.(Ψ0), Ψ0)
	plot!(rs_mid, Ψ0)
	plot!(xscale=:log10, ylabel="Ψ")

end

# ╔═╡ 2306053e-dadf-4e1c-bc6f-f13f4e43fe2a
begin 
	plot( ylabel="d2ρ/dr2", xlabel="r", xscale=:log)
	plot!(rs, ρ2.(rs))
	plot!(rs, lguys.gradient(rs, lguys.gradient(rs, ρ(rs))), label="gradient")
end

# ╔═╡ 7c52b903-00c4-450b-941a-640b43a696a8
begin 
	plot(rs, ρ1.(rs), xscale=:log)
	plot!(rs, lguys.gradient(rs, ρ(rs)), label="gradient")
	plot!(ylabel="d2ρ/dr2", xlabel="r")
end

# ╔═╡ 840cda9d-a632-4f2c-8770-e1001d6d975e
begin 
	plot(rs, 0*Ψ.(rs), xscale=:log)
	Φ2 = lguys.calc_radial_Φ(radii, snap_i.masses)
	plot!(rs, Ψ.(rs) .+ Φ2.(rs))
end

# ╔═╡ 5f29df4c-4e27-44eb-87f1-e1e8090c8a18
begin 
	plot(rs, Ψ1.(rs), xscale=:log)
	plot!(rs_mid, lguys.gradient(rs_mid, Ψs), label="gradient")
end

# ╔═╡ a23f1d63-4a48-41de-b454-c9ec4246b3c7
begin 
	plot(rs, Ψ2.(rs), xscale=:log, ylabel="d2Ψ/dr2")
	plot!(rs_mid, lguys.gradient(rs_mid, lguys.gradient(rs_mid, Ψs)), label="gradient")
end

# ╔═╡ 8b081079-a8cf-4275-a4c5-ccaa544b7d45
begin 
	plot(rs, M1.(rs), xscale=:log)
	plot!(rs, lguys.gradient(rs, Mr.(rs)), label="gradient")
	plot!()
end

# ╔═╡ 6faa3b9a-cc27-46cc-98f5-43016161e9af
begin 
	plot(rs, M2.(rs), xscale=:log)
	scatter!(rs, lguys.gradient(rs, lguys.gradient(rs, Mr.(rs))), label="gradient")
	plot!()
end

# ╔═╡ dd85599f-2d47-47bc-97d1-a952b9a1ee3a
begin 
	plot(ylims=[-0.07, 0.2], xlims=[0.08, 1], xscale=:log, ylabel="d3M(r)/dr3")
	plot!(rs, M3.(rs))
	plot!(rs, lguys.gradient(rs, lguys.gradient(rs, lguys.gradient(rs, Mr.(rs)))), label="gradient")
end

# ╔═╡ 0dccfbcf-06ca-4f94-a0d0-0e0d1ddabfce
md"""
# Adding stars
"""

# ╔═╡ aca25368-28c7-4f4a-bfcf-295902eaff9a
begin
	ρ_s(r, Γ=-0.72) = exp( Γ*r )
end

# ╔═╡ a17cff15-a6c1-4bbc-ac98-7a5bd8445f75
plot(log10.(rs), log10.(ρ_s.(rs)))

# ╔═╡ facb59d3-7dd9-403a-a3e1-4f441bc58ba4
begin 
	ρss = ρ_s.(rs_mid)
	Mss = cumsum(4π * diff(rs) .* ρss .* diff(rs .^3))
end

# ╔═╡ 691e94b6-15e3-4e78-99a6-8218d902b4e4


# ╔═╡ 3a12cfa5-61c2-4098-be34-cfa58771e5d5


# ╔═╡ b3d9c5d2-1607-498f-8838-3fdc8e814aa0
function calc_f_E(rs, Ms, Ψs)
	Mr = Dierckx.Spline1D(rs, Ms, k=5)
	rs_mid = lguys.midpoint(rs)

	r_end = rs[end]
	M1(r) = Dierckx.derivative(Mr, r)
	M2(r) = Dierckx.derivative(Mr, r, nu=2)
	M3(r) = Dierckx.derivative(Mr, r, nu=3)
	ρ2(r) = @. M3(r) / (4π * r^2) - M2(r) / (π * r^3) + M1(r) * 3/(2π * r^4) 
	ρ(r) = @. M1(r) / (4π * r^2)
	ρ1(r) = @. M2(r) / (4π * r^2) - M1(r) / (2π * r^3)
	
	Ψ = Dierckx.Spline1D(rs, Ψs)
	Ψ_inv = Dierckx.Spline1D(reverse(Ψs), reverse(rs)) # reverse so xs increase

	Ψ1(r) = @. -lguys.G * Mr(r) / r^2
	Ψ2(r) = @. -lguys.G / r^3 * (r*M1(r) - 2*Mr(r))


	d2ρ_dΨ2(r) = ρ2(r) / Ψ1(r)^2 - ρ1(r) * Ψ2(r) / Ψ1(r)^3 

	function f_integrand(r, E) 
		if E - Ψ(r) > 0
			d2ρ_dΨ2(r) * Ψ1(r) / sqrt(E - Ψ(r))
		else
			return 0
		end
	end
	
	f(E) =  -1/(√8 * π^2) * (
		quadgk(r->f_integrand(r, E) , Ψ_inv(E), r_end)[1]
		+ 0#+ 1/sqrt(E) * ρ1(0) / Ψ1(0))
	)
	return f
end

# ╔═╡ 89ab3f09-3e71-42f2-8f06-c34eb3810790
fs = calc_f_E(rs_mid, Mss, Ψ0)

# ╔═╡ 5cbb0912-9012-4cbc-911a-8ea5fda80746
begin 
	plot(ylim=(-0.01, 0.03))
	plot!(Es, fs.(Es))
	plot!(Es, f.(Es))
	vline!([Es[end]])
end

# ╔═╡ dab099c5-466d-4a97-b09d-85ca35c6b652
begin 
	plot(ylim=(0, 20))
	plot!(Es, fs.(Es) ./ f.(Es))
	vline!([Es[end]])
end

# ╔═╡ 8afbc6b8-a397-4b06-aaa3-aba90a692349
begin 
	fdm_i = Dierckx.Spline1D(Es, f.(Es))
	fs_i = Dierckx.Spline1D(Es, fs.(Es))

	probs = fs_i.(ϵ) ./ fdm_i.(ϵ)
	probs ./= maximum(probs)
	probs = ifelse.(probs .< 0, 0, probs)
end

# ╔═╡ afd043d3-000a-4d24-9654-65d4dfd042b0
begin 
	plot(log10.(radii), cumsum(probs) / sum(probs))
	plot!(log10.(radii), cumsum(ρ_s.(radii) * 4π .* radii.^2) / sum(ρ_s.(radii) * 4π .* radii.^2))
end

# ╔═╡ 776a0d95-81a0-44f5-b943-9bc9a5cb620c
begin 
	histogram(log10.(radii), weights=probs ./ radii .^ 2)
	plot!(log10.(radii), ρ_s.(radii))
end

# ╔═╡ c9f75a8a-895e-4ad1-8f02-f1744c6609ed
histogram(probs)

# ╔═╡ 34de8939-0b66-4931-ae51-17f6887da473
histogram(nm.log10.(fdm_i.(ϵ)))

# ╔═╡ 96b9aca9-d890-41fb-b1da-e74892656163
scatter(ϵ, fs_i.(ϵ) )

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═055e698b-9092-47e4-a70f-a03d5db71bb4
# ╠═07cf49a5-909a-48dc-982b-6c6953e5e414
# ╠═6964cbd2-9ef9-4682-9ac1-80c10518c374
# ╠═79ce82c5-47d8-4157-942a-ff8e524edfd4
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═4ae5506f-0e5e-44ac-b9e4-161ca7b448af
# ╟─5214efbc-4445-401f-9ef0-3d4c777fd1bd
# ╠═c7aa3122-4bdd-4889-895a-d143f8919405
# ╠═d8c66c89-8e79-4214-82d1-05cfe248ab8a
# ╠═9b9f0ee6-0856-494d-8135-98ce77995d42
# ╠═2e8de7c6-83e8-4a84-a7a4-4910401b0f4e
# ╠═32411e1f-da6f-405c-9b6c-20c4a2dbf4a5
# ╠═8db8fdc3-5b85-44ba-a8d9-99d197117536
# ╟─80b32d08-405f-4e1a-a6a0-3325afc1e283
# ╠═8938d5ce-adc9-446c-8cc3-4b637d81e480
# ╠═bd21d3e8-5a5c-40ee-9f16-3e429492d4d0
# ╠═ecd35011-2f81-4a10-b931-250a2e0b9b2c
# ╠═aad56d80-3e06-4e98-b7c9-cd013482dfdf
# ╠═a7646485-1dc9-41c0-96d3-5e0a740be11e
# ╠═6c35dc16-643d-47bc-9b7e-8f3069e0188d
# ╠═b1e5f70f-d724-40ec-a30b-df169ffd1e8b
# ╠═26cec8e1-bca0-4642-b98a-79ad2d97b076
# ╠═1ce021f2-e732-4acf-ae1d-b605513a019f
# ╠═df571e5d-2b8c-4476-b459-1c725062ac8e
# ╠═709b2c42-550c-4b82-88ae-6d940055e13f
# ╠═2677f891-b5a9-4c29-9aa0-0f305af61f95
# ╠═2028a467-6a81-46e9-a438-b73f884c0d20
# ╠═aba99571-2013-49e3-af86-eb59a0a38109
# ╠═191fd0b0-f8d3-400b-afdb-6f1022854f35
# ╠═3ef7e24e-ed44-445a-870b-e91e6ce07c65
# ╠═97baf791-8f97-4dfb-94ca-d1496162f688
# ╠═2306053e-dadf-4e1c-bc6f-f13f4e43fe2a
# ╠═7c52b903-00c4-450b-941a-640b43a696a8
# ╠═840cda9d-a632-4f2c-8770-e1001d6d975e
# ╠═5f29df4c-4e27-44eb-87f1-e1e8090c8a18
# ╠═a23f1d63-4a48-41de-b454-c9ec4246b3c7
# ╠═8b081079-a8cf-4275-a4c5-ccaa544b7d45
# ╠═6faa3b9a-cc27-46cc-98f5-43016161e9af
# ╠═dd85599f-2d47-47bc-97d1-a952b9a1ee3a
# ╟─0dccfbcf-06ca-4f94-a0d0-0e0d1ddabfce
# ╠═aca25368-28c7-4f4a-bfcf-295902eaff9a
# ╠═a17cff15-a6c1-4bbc-ac98-7a5bd8445f75
# ╠═facb59d3-7dd9-403a-a3e1-4f441bc58ba4
# ╠═89ab3f09-3e71-42f2-8f06-c34eb3810790
# ╠═5cbb0912-9012-4cbc-911a-8ea5fda80746
# ╠═dab099c5-466d-4a97-b09d-85ca35c6b652
# ╠═afd043d3-000a-4d24-9654-65d4dfd042b0
# ╠═776a0d95-81a0-44f5-b943-9bc9a5cb620c
# ╠═691e94b6-15e3-4e78-99a6-8218d902b4e4
# ╠═8afbc6b8-a397-4b06-aaa3-aba90a692349
# ╠═c9f75a8a-895e-4ad1-8f02-f1744c6609ed
# ╠═34de8939-0b66-4931-ae51-17f6887da473
# ╠═3a12cfa5-61c2-4098-be34-cfa58771e5d5
# ╠═96b9aca9-d890-41fb-b1da-e74892656163
# ╠═b3d9c5d2-1607-498f-8838-3fdc8e814aa0
