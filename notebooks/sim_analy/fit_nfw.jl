### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 6eebcbba-d378-11ef-155b-39f2b4447aa2
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya
	using LilGuys

end

# ╔═╡ e832a69c-bc33-4b4e-a27b-cced818519cc
md"""
Given a model, the goal of this notebook is to read in the dark matter profile and fit it with an analytic model.
"""

# ╔═╡ 97a7d911-8a1c-4d8d-b44f-8f55df7453f5
import LinearAlgebra: diag

# ╔═╡ 897ffe20-ea2a-4e0f-b0d1-1b3e37619561
modelname = "sculptor/1e7_V31_r3.2/orbit_smallperi"

# ╔═╡ 481636bd-1cdb-4640-8425-8cefa76b0279
idx_p = -1

# ╔═╡ 4d57b5d0-b045-46c9-ae7a-6b9e7e4ac842
let
	global profs, prof
	
	profs = LilGuys.read_structs_from_hdf5(joinpath(
	ENV["DWARFS_ROOT"], "analysis", modelname, "profiles.hdf5"
), LilGuys.MassProfile3D)

	idx = parse.(Int, first.(profs))

	profs = last.(profs[sortperm(idx)])

	if idx_p > 0
		prof = profs[idx_p]
	else
		prof = profs[end + 1 + idx_p]
	end
end

# ╔═╡ e9c11c63-3357-450d-a9a5-397fb0863465
log_rc = log10.(prof.r_circ)

# ╔═╡ 82c6793e-fae6-4737-880b-39f03036ab77
vc = prof.v_circ

# ╔═╡ cf790d58-0c3e-4c77-a124-bac7482bd3bc
vc_err = prof.v_circ_err

# ╔═╡ ef411010-e176-46d6-b2d0-f25de8483e7b
weights = 1 ./ vc_err .^ 2

# ╔═╡ b688e34d-15c1-490e-9ee1-0c96784b87a8
halo_type = LilGuys.TruncNFW

# ╔═╡ 8fb36116-0b32-4bfe-a972-780a939a5044
halo_type

# ╔═╡ 76039468-5b7e-4669-aaf9-64f675c55d38
function make_halo(args)
	if any(args .< 0)
		return NFW(M_s=0, r_s=1, c=1)
	end
	
	LilGuys.TruncNFW(M200=args[1], r_s=args[2], r_t=args[3])
end

# ╔═╡ 5b2dca6f-dd92-4e52-b2ea-5c10c479528d
function make_halo2(args)
	if any(args .< 0)
		return NFW(M_s=0, r_s=1, c=1)
	end
	
	LilGuys.ExpCusp(args[1], args[2])
end

# ╔═╡ f6ab76e4-6c09-4067-a3fb-c4df79d5236b
model(x, args) = LilGuys.calc_v_circ.(make_halo(args), 10 .^ x)

# ╔═╡ f0bfbcc0-87bb-41e4-ad70-f9cfbcf9fe76
model2(x, args) = LilGuys.calc_v_circ.(make_halo2(args), 10 .^ x)

# ╔═╡ d5b03cc0-006d-4213-9461-a7be97648d79
p0 = [0.10, 1, 3] # M0, r_s r_t

# ╔═╡ 24f16177-cc3b-4454-8fd5-693c5d87dc1d
md"""
## Density
"""

# ╔═╡ 41b6899f-920f-481d-8007-c97acc5060b6
model_ρ(x, args) = log10.(calc_ρ.(make_halo(args), 10 .^ x))

# ╔═╡ 30e8a4c2-7b52-4c7a-a552-d348802d4149
model_ρ2(x, args) = log10.(calc_ρ.(make_halo2(args), 10 .^ x))

# ╔═╡ 83c20412-7f38-432a-a4d6-49fe14b1df94
begin
	_ρ_filt = prof.counts .> 10
	log_ρ = log10.(prof.rho[_ρ_filt])
	log_r = prof.log_r[_ρ_filt]
	log_ρ_err = (prof.rho_err ./ prof.rho)[_ρ_filt]
end

# ╔═╡ f7b56544-873b-4619-b1ad-347a56aa0ffb
weights_ρ = 1 ./ log_ρ_err .^ 2

# ╔═╡ 371c8adc-3c24-45f9-92b1-4c488c8ac79b
log_r_max = 0.5

# ╔═╡ ec1ce63a-6c53-4ca4-b054-dd4d2a9c435c
popt, covt = LilGuys.curve_fit(model, log_rc[log_rc .< log_r_max], vc[log_rc .< log_r_max], weights[log_rc .< log_r_max], p0)

# ╔═╡ bc17d2c3-9a21-4633-a132-a397cb1374d9
sqrt.(diag(covt))

# ╔═╡ 4b6a3782-62d6-4dcc-844f-9bdfcaef6be5
popt_2, covt_2 = LilGuys.curve_fit(model2, log_rc[log_rc .< log_r_max], vc[log_rc .< log_r_max], weights[log_rc .< log_r_max], p0)

# ╔═╡ 7fd35c7c-83bc-4fa4-9fd4-caf90ca4ef04
r_filt = log_r .< log_r_max

# ╔═╡ 35723d03-3c41-488e-9d58-c498720d5423
popt2, covt2 = LilGuys.curve_fit(model_ρ, log_r[r_filt], log_ρ[r_filt], [2, 1.2, 2.5])

# ╔═╡ bb2e20a2-4278-4401-b8f6-efdfd4a1d86d
popt_rho2, covt_rho2 = LilGuys.curve_fit(model_ρ2, log_r[r_filt], log_ρ[r_filt], [0.01, 0.5])

# ╔═╡ 097f32f0-bd6f-4ffe-9481-48656d1420f0
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = L"$\log\,r$ / kpc",
		ylabel = L"$v_\textrm{circ}$",
	)

	lines!(log_rc, vc)

	# y_m = model(log_rc, popt)

	# lines!(log_rc, y_m)

	# y_m = model(log_rc, popt2)

	# lines!(log_rc, y_m)


	y_m = model2(log_rc, popt_rho2)

	lines!(log_rc, y_m)

	fig
end

# ╔═╡ 7e5e580b-eb06-4920-91e2-41297335d4a0
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = L"$\log\,r$ / kpc",
		ylabel = L"$\log \rho$",
		limits=(nothing, nothing, -10,1)
	)

	lines!(log_r, log_ρ)

	# y_m = log10.(calc_ρ.(best_halo, 10 .^ prof.log_r))

	# lines!(prof.log_r, y_m)

	# y_m = log10.(calc_ρ.(best_halo_2, 10 .^ prof.log_r))
	# lines!(prof.log_r, y_m)

	# y_m = log10.(calc_ρ.(LilGuys.ExpCusp(0.01, 0.6), 10 .^ prof.log_r))
	# lines!(prof.log_r, y_m)

	y_m = model_ρ2(prof.log_r, popt_rho2)
	lines!(prof.log_r, y_m)

	fig
end

# ╔═╡ b4184145-9668-4b40-a585-e78184c486a7
best_halo = LilGuys.TruncNFW(M_s=popt[1], r_s=popt[2], r_t=popt[3])

# ╔═╡ 9009d2b9-c9d1-4683-99f2-7790be1fae14
best_halo_2 = make_halo(popt2)

# ╔═╡ fbb0f59f-e1d3-406d-b524-1794455a68bb
make_halo2(popt_rho2)

# ╔═╡ Cell order:
# ╠═e832a69c-bc33-4b4e-a27b-cced818519cc
# ╠═6eebcbba-d378-11ef-155b-39f2b4447aa2
# ╠═97a7d911-8a1c-4d8d-b44f-8f55df7453f5
# ╠═897ffe20-ea2a-4e0f-b0d1-1b3e37619561
# ╠═481636bd-1cdb-4640-8425-8cefa76b0279
# ╠═4d57b5d0-b045-46c9-ae7a-6b9e7e4ac842
# ╠═e9c11c63-3357-450d-a9a5-397fb0863465
# ╠═82c6793e-fae6-4737-880b-39f03036ab77
# ╠═cf790d58-0c3e-4c77-a124-bac7482bd3bc
# ╠═ef411010-e176-46d6-b2d0-f25de8483e7b
# ╠═b688e34d-15c1-490e-9ee1-0c96784b87a8
# ╠═8fb36116-0b32-4bfe-a972-780a939a5044
# ╠═76039468-5b7e-4669-aaf9-64f675c55d38
# ╠═5b2dca6f-dd92-4e52-b2ea-5c10c479528d
# ╠═f6ab76e4-6c09-4067-a3fb-c4df79d5236b
# ╠═f0bfbcc0-87bb-41e4-ad70-f9cfbcf9fe76
# ╠═d5b03cc0-006d-4213-9461-a7be97648d79
# ╠═ec1ce63a-6c53-4ca4-b054-dd4d2a9c435c
# ╠═4b6a3782-62d6-4dcc-844f-9bdfcaef6be5
# ╠═35723d03-3c41-488e-9d58-c498720d5423
# ╠═bb2e20a2-4278-4401-b8f6-efdfd4a1d86d
# ╠═bc17d2c3-9a21-4633-a132-a397cb1374d9
# ╠═097f32f0-bd6f-4ffe-9481-48656d1420f0
# ╠═24f16177-cc3b-4454-8fd5-693c5d87dc1d
# ╠═41b6899f-920f-481d-8007-c97acc5060b6
# ╠═30e8a4c2-7b52-4c7a-a552-d348802d4149
# ╠═83c20412-7f38-432a-a4d6-49fe14b1df94
# ╠═f7b56544-873b-4619-b1ad-347a56aa0ffb
# ╠═371c8adc-3c24-45f9-92b1-4c488c8ac79b
# ╠═7fd35c7c-83bc-4fa4-9fd4-caf90ca4ef04
# ╠═7e5e580b-eb06-4920-91e2-41297335d4a0
# ╠═b4184145-9668-4b40-a585-e78184c486a7
# ╠═9009d2b9-c9d1-4683-99f2-7790be1fae14
# ╠═fbb0f59f-e1d3-406d-b524-1794455a68bb
