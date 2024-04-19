### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 949f4b7e-f76f-11ee-31ed-69a587dc7973
begin
	import Pkg; Pkg.activate()

	using Plots; #pythonplot()
	using Arya

	import LilGuys as lguys

	import CSV
	using DataFrames
end

# ╔═╡ 0c85a109-f75b-45b2-be77-1664bfdd4567
path = "fiducial"

# ╔═╡ 7bb77640-17c4-44b1-832f-9258fb790d8b
out = lguys.Output("$path/out")

# ╔═╡ 390bc099-3434-4c19-924c-bde742f9f46c
begin
	cens = CSV.read("$path/data/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ e82adf25-15c3-44dc-ba1f-46e962ef9b14
lguys.plot_xyz_layout(x_cen)

# ╔═╡ 2c426184-6e1f-4256-bd5b-fc70d20e0481
begin 
	plot(xlabel="time / Gyr", ylabel="galactocentric distance / kpc")
	plot!(cens.t * lguys.T0, lguys.calc_r(x_cen), label="")
end

# ╔═╡ 277b7a1c-f9e7-4645-97e2-e3bc70b145c3
begin 
	r_ini = lguys.calc_r(out[1].positions, x_cen[:, 1])[sortperm(out[1].index)]
end

# ╔═╡ 0bf41856-2db7-4f4a-b888-f0e9fc5deeaa
begin 
	pos = lguys.extract_vector(out[1], :positions)
	scatter(pos[1, :], pos[2, :], marker_z=r_ini, ms=1, colorbar_scale=:log10, label="", clim=(1e-3, 0.3e-1))
end

# ╔═╡ 8f0bb0fc-c1df-4885-ad2a-00b5ec70743e
# ╠═╡ disabled = true
#=╠═╡
begin 
	anim2 =  @animate for i in 1:20:length(out)
		pos2 = lguys.extract_vector(out[i], :positions)
		scatter(pos2[1, :] .- x_cen[1, i], pos2[2, :] .- x_cen[2, i], marker_z=r_ini, ms=1, colorbar_scale=:log10, label="", xlim=(-0.1, 0.1), ylim=(-0.1, 0.1), aspect_ratio=1, clim=(1e-3, 0.3e-1))
	end
	gif(anim2, "centre.gif", fps=6)
end
  ╠═╡ =#

# ╔═╡ 40d16524-39ea-49c9-8c0e-69708c74c512
sort(out[1].index)

# ╔═╡ 4b456983-017f-47e1-a76a-0688fa53b6fb
begin 
	plot(xlabel="R / kpc", ylabel="z / kpc", aspect_ratio=1, xlim=[0, 50])
	R_cen = @. sqrt(x_cen[1, :]^2 + x_cen[2, :]^2)
	plot!(R_cen, x_cen[3, :], label="")
end

# ╔═╡ 4f0a2096-0e67-43f9-99da-6fdb00fe2f19
# ╠═╡ disabled = true
#=╠═╡
begin 
	anim = @animate for i in 1:20:length(out)
		snap = out[i]
		lguys.plot_centre(snap.positions, 
			ms=1, msw=0, ma=0.1, width=50, z_filter=:none)
		# scatter!(x_cen[1, i:i], x_cen[2, i:i])
	end
	
	gif(anim, "orbit1.gif", fps = 12)
end
  ╠═╡ =#

# ╔═╡ 6f933ea2-0fbe-409a-912d-01ac9b853160
begin 
	lguys.plot_centre(out[end].positions .- x_cen[:, end], z_filter=:none, ms=2, alpha=0.1)
	scatter!([0], [0])
end

# ╔═╡ 66085044-36ed-4198-a80b-f4e733960a5c
dpi=200

# ╔═╡ fa96effc-6e63-48ae-a988-32fa8265043a
# pythonplot(); Arya.set_mpl_theme(); theme(:arya, ticks=:naitive)

# ╔═╡ 9e7ec2b0-2f11-459c-8ee2-8294b6c7e5f0
pythonplot()

# ╔═╡ 2fe3863a-9b3a-4186-9158-89f557051460
Arya.set_mpl_theme()

# ╔═╡ f8fe0a8f-c88e-4a62-b181-87f362817856
begin 
	jax = CSV.read("../jax_membs.tsv", DataFrame)
end

# ╔═╡ 3b95c66e-0a45-4e69-86d7-4208c64ffaa6
idx_f = length(out) - 38

# ╔═╡ f83658c4-dc0d-41e8-9707-26bcd1344d1b
begin 
	plot(xlabel="x / kpc", ylabel="y / kpc")
	scatter!(out[idx_f].positions[1, :], out[idx_f].positions[2, :], ms=1, aspect_ratio=1, alpha=0.1, label="")
end

# ╔═╡ 59382ef6-ecc6-4915-a3d5-02827729ab25
x_cen[:, idx_f]

# ╔═╡ fa20f05a-951f-44e1-ae8b-6e24393304f6
snap_f = out[idx_f]

# ╔═╡ 10456dbb-1aa4-436b-a80f-39ec29fa8137
particles_gc = [lguys.Galactocentric(snap_f.positions[:, i]..., (lguys.V0 * snap_f.velocities[:, i])...) for i in 1:length(snap_f)]

# ╔═╡ 1b027221-000b-4929-b893-7f2aa829b3ed
obs_pred = lguys.transform.(lguys.Observation, particles_gc)

# ╔═╡ 44353e99-0968-4619-9c7d-bea32be287fe
begin 
	plot(xlabel="ra", ylabel="dec")
	scatter!([o.ra for o in obs_pred], [o.dec for o in obs_pred], colorbar_scale=:log10, ms=1)
	scatter!(jax.RA, jax.DEC, ms=2, label="Jensen+2021")

end

# ╔═╡ 104d4ce4-3324-43be-8978-9c1650b15f51
histogram([o.radial_velocity for o in obs_pred], bins=LinRange(90, 120, 100))

# ╔═╡ ef36fcd0-abef-45d2-ab13-00c3fea80ff3
begin 
	plot()
	scatter!([o.distance for o in obs_pred], [o.radial_velocity for o in obs_pred], ms=2, alpha=0.3)
	scatter!([16], [106])
end

# ╔═╡ 6fe7539f-fb34-45e4-abc8-d1000b6ad1ab
begin 
	plot()
	scatter!([o.pm_ra for o in obs_pred], [o.pm_dec for o in obs_pred], alpha=0.5, label="", ms=1)
	scatter!([-5.41], [-0.79])
end

# ╔═╡ 620e20ca-1d7a-42be-bd4d-b62f5573f48b
obs_m = lguys.transform(lguys.Observation, lguys.Galactocentric(x_cen[:, idx_f]..., (v_cen[:, idx_f] * lguys.V0)...))

# ╔═╡ a80ec84b-6bf9-4130-9e87-b00ee6d340cf
begin 
	dra = 30
	ddec = 10
	plot(xlabel=L"\rm RA\, /\,^{\circ}", ylabel=L"\rm Dec\, /\,^{\circ} ", 
		xlims=(obs_m.ra - dra, obs_m.ra + dra), 
	ylims=(obs_m.dec - ddec, obs_m.dec + ddec), ticks=:native,
	size=dpi/1.5 .* (3.35, 2), dpi=dpi, font_size=22)
	
	scatter!([o.ra for o in obs_pred], [o.dec for o in obs_pred], alpha=0.5, label="simulation", ms=1)
	scatter!(jax.RA, jax.DEC, ms=2, label="Jensen+2021")
	#scatter!([obs_m.ra], [obs_m.dec], label="centre")

	#savefig("stream.pdf")
end

# ╔═╡ 04fc7b34-f672-46fd-8d8b-7c640874a2c4
xi, eta = lguys.to_tangent([o.ra for o in obs_pred], [o.dec for o in obs_pred], obs_m.ra, obs_m.dec)

# ╔═╡ d7e56838-7881-4b7c-af50-9b0f3d259a92
begin 
	plot(aspect_ratio=1,
	xlabel="\\xi", ylabel="\\eta")
	histogram2d!(xi*60, eta*60, bins=LinRange(-120, 120, 300), label="", colorbar_scale=:log10)
	#plot!(t->r_h_amin*sin(t), t->r_h_amin*cos(t), LinRange(0, 2π, 1000))
end

# ╔═╡ ec6636df-98d4-41a9-b343-8ec81224b620
begin 
	r_f, rho_f = lguys.calc_ρ_hist(lguys.calc_r(snap_f.positions, x_cen[:, idx_f]),20,  weights=snap_f.masses)
	r_f = lguys.midpoint(r_f)

	r_i, rho_i = lguys.calc_ρ_hist(lguys.calc_r(out[1].positions, x_cen[:, 1]), 20, weights=out[1].masses)
	r_i = lguys.midpoint(r_i)

	idx_m = 500
	r_m, rho_m = lguys.calc_ρ_hist(lguys.calc_r(out[idx_m].positions, x_cen[:, idx_m]), weights=out[idx_m].masses)
	r_m = lguys.midpoint(r_m)
end

# ╔═╡ d9a28fdc-d702-4cbd-9aaf-5b214a0d0a23
lguys.calc_r(out[1].positions, x_cen[:, 1])

# ╔═╡ e85ea859-4998-42c8-90ff-ccde7aeb070b
pythonplot

# ╔═╡ b8d2b64c-d908-4f2e-af5d-bd2e8bf15471
begin
	df_rho = DataFrame()
	df_rho[!, "r"] = r_f
	df_rho[!, "rho"] = rho_f
end

# ╔═╡ 9279f480-1e0e-4a64-b959-161939bfb38e
vs = lguys.calc_r(out[idx_f].velocities, v_cen[:, idx_f])

# ╔═╡ f0e72387-3247-4a1d-a582-9e73c436ed5a
rs = lguys.calc_r(out[idx_f].positions, x_cen[:, idx_f])

# ╔═╡ f3d672e0-6cc0-48c2-958d-d635c67cf97c
rs_s = sort(rs)

# ╔═╡ b01fba62-cfcd-4c21-886b-119837eafc2e
pot = lguys.calc_radial_discrete_Φ(out[idx_f].masses, rs_s)

# ╔═╡ eeffe348-46f7-475a-a2bd-02bd291b0b31
gr()

# ╔═╡ aadb3a9a-f54d-4333-b82e-35768bef7a86
ϵs = 1/2 * vs[sortperm(rs)].^2 .+ pot

# ╔═╡ ddf3c730-882d-4f86-b61d-8f936c52e06a
histogram(ϵs, bins=LinRange(-0.01, 0.01, 100))

# ╔═╡ b405662f-2708-4437-b00d-6bb9ea4b03a9
sum(ϵs .< 0)

# ╔═╡ 955abecb-54d2-4917-ac4b-003d04a305e9
rho_i

# ╔═╡ a5620a88-0afd-4bea-8a93-bae7f0d12086
r_h = 2.3 * 60 / 206265 * 16

# ╔═╡ d76d03d9-7080-4fc4-b7ea-9a3a69516d47
begin 
	plot(xlabel="log r / kpc", ylabel=L"\log \rho", xlims=(-2, 2), ylims=(-7, 2))
	plot!(log10.(r_i ./ r_h), log10.(rho_i), label="initial")
	#plot!(log10.(r_m ./ r_h), log10.(rho_m), label="intermediate")

	plot!(log10.(r_f ./ r_h), log10.(rho_f), label="present day")
end

# ╔═╡ 7b48fcf4-a66b-400a-8d23-64d8693d254e
histogram(vs[rs .< 3*r_h])

# ╔═╡ dcd893e8-84ca-4e54-8ee4-c4363a228947
histogram(rs ./ r_h, bins=LinRange(0, 10_000, 100), yscale=:log10)

# ╔═╡ 8814836e-1089-48bf-b349-912a0e532d6e
histogram(log10.(rs ./ r_h), bins=60)

# ╔═╡ 579be7d8-0cac-45ce-acc4-7b7866339dfd
r_h_amin = 2.3

# ╔═╡ 96f791cf-63ad-4cc5-bead-aebd794446b6
begin 
	plot(xlims=(-5, 5), ylims=(-5, 5), aspect_ratio=1, size=(330, 330),
	xlabel=L"\xi\, /\, \rm arcmin", ylabel=L"\eta\, /\, \rm arcmin")
	scatter!(xi*60, eta*60, ms=2, alpha=0.1, label="")
	plot!(t->r_h_amin*sin(t), t->r_h_amin*cos(t), LinRange(0, 2π, 1000), color=:black, label="")
	savefig("tangent_centre.pdf")
end

# ╔═╡ f08d508c-2e53-4840-832b-ec71a1fe3f67
m_p = out[1].masses[1]

# ╔═╡ acdefcc1-3fd0-47c1-955c-622f8aafe60b
v_circ = [sqrt(sum(rs .< r) * m_p / r) for r in rs]

# ╔═╡ 372f081d-a039-4747-9e30-0ed748d87608
scatter(log10.(rs), v_circ * lguys.V0)

# ╔═╡ 92162199-8716-4606-ad1b-876bc702c175
histogram(r_ini)

# ╔═╡ 3d69bb2a-b479-47f6-a015-21ac5b82a4f4
lguys.percentile(r_ini, 25), lguys.percentile(r_ini, 75) # inspiration for cuts below

# ╔═╡ 757c6f95-bfdd-4435-9964-f6567c9d10d1
begin 
	M_hs = []
	M_3hs = []
	f_bound = []

	# bound fractions inside, middle, and out
	f_bound_i = []
	f_bound_m = []
	f_bound_o = []
	r_in = 0.008 # how to split the particles
	r_out = 0.02

	v_max = []
	
	for i in 1:length(out)
		radii = lguys.calc_r(out[i].positions, x_cen[:, i])
		velocities = lguys.calc_r(out[i].velocities, v_cen[:, i])
		pots = out[i].Φs
		es = 1/2 * velocities .^ 2 .+ pots
		m = sum(radii .< r_h)
		push!(M_hs, m)
		
		m = sum(radii .< 3r_h)
		push!(M_3hs, m)

		push!(f_bound, lguys.mean(es .< 0))
		
		index = out[i].index
		ri =  r_ini[index]
		
		filt = ri .< r_in
		push!(f_bound_i, lguys.mean(es[filt] .< 0))
		filt = r_in .< ri .< r_out
		push!(f_bound_m, lguys.mean(es[filt] .< 0))
		filt = r_out .< ri
		push!(f_bound_o, lguys.mean(es[filt] .< 0))


		radii = sort(radii)
		M_r = m_p .* collect(1:length(radii))
		v_circ = @. sqrt(lguys.G * M_r ./ radii) * lguys.V0 # in kms
		push!(v_max, maximum(v_circ))
	end
	
	M_hs .*= m_p * lguys.M0
	M_3hs .*= m_p * lguys.M0
end

# ╔═╡ 6554d590-7269-473a-bb3f-da7f741e53f6
begin
	plot(xlabel="time / Gyr", ylabel="M(r_h)")
	scatter!(out.times, M_hs, s=1)
end

# ╔═╡ 1f9a3ec3-8c65-419b-abd2-4921b19aa092
begin 
	plot(xlabel="time / Gyr", ylabel="bound fraction")
	plot!(out.times, f_bound, s=1, label="")
end

# ╔═╡ e757ce6c-ca9a-4eda-a986-aa7f00451679
begin 
	plot(xlabel="time / Gyr", ylabel="max circular velocity")
	plot!(out.times, v_max, s=1, label="")
end

# ╔═╡ a92a2800-dc22-4461-848e-4c6cc0669c16
begin 
	plot(xlabel="time / Gyr", ylabel="bound fraction")
	plot!(out.times, f_bound_i, s=1, label="innermost")
	plot!(out.times, f_bound_m, s=1, label="middle")
	plot!(out.times, f_bound_o, s=1, label="outer")
end

# ╔═╡ c35607a8-ed18-4491-8425-ff76160d7ead
begin
	plot(xlabel="time / Gyr", ylabel="M(r_h)")
	scatter!(out.times, M_3hs, ms=1)
end

# ╔═╡ 7ca8f7cf-14c5-4d2b-aac6-093a77156183
pwd()

# ╔═╡ cc6dec8f-04e1-4794-a061-124e99cbd0d2
begin 
	time_evol = DataFrame()
	time_evol[!, "times"] = out.times .- out.times[idx_f]
	time_evol[!, "M_3r"] = M_3hs
	time_evol[!, "f"] = f_bound
	time_evol[!, "f_i"] = f_bound_i
	time_evol[!, "f_m"] = f_bound_m
	time_evol[!, "f_o"] = f_bound_o
	time_evol[!, "rs"] = lguys.calc_r(x_cen)
	time_evol[!, "xs"] = x_cen[1, :]
	time_evol[!, "ys"] = x_cen[2, :]
	time_evol[!, "zs"] = x_cen[3, :]

	time_evol
end

# ╔═╡ 6e5ec3e8-eb48-417a-abd0-b346601df64b
CSV.write("$path/nbody_orbit.csv", time_evol )

# ╔═╡ 38ffe923-2092-481c-b3ac-f5c15d43109f
CSV.write("$path/density_prof.csv", df_rho)

# ╔═╡ Cell order:
# ╠═949f4b7e-f76f-11ee-31ed-69a587dc7973
# ╠═0c85a109-f75b-45b2-be77-1664bfdd4567
# ╠═7bb77640-17c4-44b1-832f-9258fb790d8b
# ╠═390bc099-3434-4c19-924c-bde742f9f46c
# ╠═e82adf25-15c3-44dc-ba1f-46e962ef9b14
# ╠═2c426184-6e1f-4256-bd5b-fc70d20e0481
# ╠═277b7a1c-f9e7-4645-97e2-e3bc70b145c3
# ╠═0bf41856-2db7-4f4a-b888-f0e9fc5deeaa
# ╠═8f0bb0fc-c1df-4885-ad2a-00b5ec70743e
# ╠═40d16524-39ea-49c9-8c0e-69708c74c512
# ╠═4b456983-017f-47e1-a76a-0688fa53b6fb
# ╠═4f0a2096-0e67-43f9-99da-6fdb00fe2f19
# ╠═f83658c4-dc0d-41e8-9707-26bcd1344d1b
# ╠═6f933ea2-0fbe-409a-912d-01ac9b853160
# ╠═59382ef6-ecc6-4915-a3d5-02827729ab25
# ╠═fa20f05a-951f-44e1-ae8b-6e24393304f6
# ╠═10456dbb-1aa4-436b-a80f-39ec29fa8137
# ╠═1b027221-000b-4929-b893-7f2aa829b3ed
# ╠═620e20ca-1d7a-42be-bd4d-b62f5573f48b
# ╠═44353e99-0968-4619-9c7d-bea32be287fe
# ╠═66085044-36ed-4198-a80b-f4e733960a5c
# ╠═fa96effc-6e63-48ae-a988-32fa8265043a
# ╠═a80ec84b-6bf9-4130-9e87-b00ee6d340cf
# ╠═04fc7b34-f672-46fd-8d8b-7c640874a2c4
# ╠═d7e56838-7881-4b7c-af50-9b0f3d259a92
# ╠═9e7ec2b0-2f11-459c-8ee2-8294b6c7e5f0
# ╠═2fe3863a-9b3a-4186-9158-89f557051460
# ╠═96f791cf-63ad-4cc5-bead-aebd794446b6
# ╠═f8fe0a8f-c88e-4a62-b181-87f362817856
# ╠═104d4ce4-3324-43be-8978-9c1650b15f51
# ╠═ef36fcd0-abef-45d2-ab13-00c3fea80ff3
# ╠═6fe7539f-fb34-45e4-abc8-d1000b6ad1ab
# ╠═3b95c66e-0a45-4e69-86d7-4208c64ffaa6
# ╠═ec6636df-98d4-41a9-b343-8ec81224b620
# ╠═d9a28fdc-d702-4cbd-9aaf-5b214a0d0a23
# ╠═e85ea859-4998-42c8-90ff-ccde7aeb070b
# ╠═d76d03d9-7080-4fc4-b7ea-9a3a69516d47
# ╠═b8d2b64c-d908-4f2e-af5d-bd2e8bf15471
# ╠═9279f480-1e0e-4a64-b959-161939bfb38e
# ╠═f0e72387-3247-4a1d-a582-9e73c436ed5a
# ╠═7b48fcf4-a66b-400a-8d23-64d8693d254e
# ╠═dcd893e8-84ca-4e54-8ee4-c4363a228947
# ╠═8814836e-1089-48bf-b349-912a0e532d6e
# ╠═acdefcc1-3fd0-47c1-955c-622f8aafe60b
# ╠═f3d672e0-6cc0-48c2-958d-d635c67cf97c
# ╠═b01fba62-cfcd-4c21-886b-119837eafc2e
# ╠═eeffe348-46f7-475a-a2bd-02bd291b0b31
# ╠═aadb3a9a-f54d-4333-b82e-35768bef7a86
# ╠═ddf3c730-882d-4f86-b61d-8f936c52e06a
# ╠═b405662f-2708-4437-b00d-6bb9ea4b03a9
# ╠═372f081d-a039-4747-9e30-0ed748d87608
# ╠═955abecb-54d2-4917-ac4b-003d04a305e9
# ╠═a5620a88-0afd-4bea-8a93-bae7f0d12086
# ╠═579be7d8-0cac-45ce-acc4-7b7866339dfd
# ╠═f08d508c-2e53-4840-832b-ec71a1fe3f67
# ╠═92162199-8716-4606-ad1b-876bc702c175
# ╠═3d69bb2a-b479-47f6-a015-21ac5b82a4f4
# ╠═757c6f95-bfdd-4435-9964-f6567c9d10d1
# ╠═6554d590-7269-473a-bb3f-da7f741e53f6
# ╠═1f9a3ec3-8c65-419b-abd2-4921b19aa092
# ╠═e757ce6c-ca9a-4eda-a986-aa7f00451679
# ╠═a92a2800-dc22-4461-848e-4c6cc0669c16
# ╠═c35607a8-ed18-4491-8425-ff76160d7ead
# ╠═7ca8f7cf-14c5-4d2b-aac6-093a77156183
# ╠═cc6dec8f-04e1-4794-a061-124e99cbd0d2
# ╠═6e5ec3e8-eb48-417a-abd0-b346601df64b
# ╠═38ffe923-2092-481c-b3ac-f5c15d43109f
