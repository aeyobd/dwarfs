### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using Plots
	using DataFrames, CSV

	using Printf
	import LilGuys as lguys
end


# ╔═╡ d975d00c-fd69-4dd0-90d4-c4cbe73d9754
using Statistics

# ╔═╡ a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
if !@isdefined obs # read in sample (but only once)
	include("sample.jl")
	obs
end

# ╔═╡ ed6ec0e9-5734-4569-b214-2c86c22d3c55
plots = Ref{Vector{Plots.Plot}}() # so we can reuse this variable

# ╔═╡ cd9dea36-93f8-4461-85b1-87030d9367bb
pwd()

# ╔═╡ 9c7e5bb4-8db0-4527-a8ec-331e2aed958b
begin
	filename = "out/combined.hdf5"
	out = lguys.Output(filename);
end

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
begin
	df = CSV.read("peris_apos.csv", DataFrame)
	df_idx = sortperm(df[!, :index])
	peris = df[df_idx, :pericenter]
	apos = df[df_idx, :apocenter];
end

# ╔═╡ 4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
begin 
	snap = out[1]
	snap = snap[sortperm(snap.index)]
end

# ╔═╡ fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
# ╠═╡ disabled = true
#=╠═╡
begin
	plots[] = []
	p1 = histogram(peris)
	xlabel!("pericenter")
	push!(plots[], p1)

	p1 = histogram(apos, bins=20)
	xlabel!("apocenter")
	push!(plots[], p1)

	plots[]
end
  ╠═╡ =#

# ╔═╡ e10615ab-0ece-4ebb-81d0-47ffdd6ea6c5
begin
	percies = [16, 84]
	quantiles = lguys.percentile(peris, [16, 84])
	idx_qs = [findfirst(isequal(p), peris) for p in quantiles]
	idx = vcat([1], idx_qs)
end

# ╔═╡ 9fad353f-751e-4577-b704-aa2eadeb0969
quantiles

# ╔═╡ 92b28dd8-c353-4959-b181-a843367b3223
histogram(lguys.calc_r(snap.velocities))

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
histogram(lguys.calc_r(snap.positions))

# ╔═╡ 5a5ca70c-286f-4325-ac4a-143713a844c4
histogram(snap.Φs_ext)

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = lguys.calc_E_spec_kin(out[1]) + out[1].Φs_ext

# ╔═╡ 8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
maximum(ϵ)

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	points = [lguys.Galactocentric(snap.positions[:, i]*lguys.R0, -snap.velocities[:, i]*lguys.V0) for i in 1:length(snap)]
	observations = lguys.transform.(lguys.Observation, points)
end


# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ ac81acd8-4a78-4230-bc70-3b78a861b618
begin
plots[] = []
	
for sym in [:distance, :pm_ra, :pm_dec, :radial_velocity]
    x = getproperty.(observations, sym)
    p = histogram(x, normalize=:pdf)

    μ = getproperty(obs, sym)
    σ = getproperty(err, sym)
    
    x_mod = LinRange(μ - 3σ, μ + 3σ, 1000)
    y_mod = normal_dist.(x_mod, μ, σ)
    plot!(p, x_mod, y_mod, z_order=2, lw=2)
    title!(p, string(sym))

    push!(plots[], p)
end

plots[]
end

# ╔═╡ d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╠═╡ disabled = true
#=╠═╡
begin 

plots_scat = []
for sym in [:distance, :pm_ra, :pm_dec, :radial_velocity, :ra, :dec]
    x = [getproperty(o, sym) for o in observations]
    y = peris
    p = scatter(x, y, ms=0.3)
	xlabel!(p, string(sym))
	ylabel!(p, "pericenter")

    push!(plots_scat, p)
end

plots_scat
end
  ╠═╡ =#

# ╔═╡ 4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
begin 
	peri1 = lguys.percentile(peris, 10)
	peri2 = lguys.percentile(peris, 20)

	idx_s = @. peri1 < peris < peri2

	for sym in [:distance, :pm_ra, :pm_dec, :radial_velocity, :ra, :dec]
		md = median(getproperty.(observations[idx_s], sym))
		println("d $sym = $(md - getproperty(obs, sym))")
	end
end

# ╔═╡ e5825c4a-b446-44a3-8fd5-d94664965aca
for f in fieldnames(lguys.Observation)
	a = getproperty(observations[1], f)
	b = getproperty(obs, f)
	println("error of $f:\t", (a-b)/b)
end

# ╔═╡ a644a33e-57a4-4af4-98be-1e0e84948069


# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = lguys.extract_vector(out, :positions, idx)
	velocities = lguys.extract_vector(out, :velocities, idx)
	accelerations = lguys.extract_vector(out, :accelerations, idx)
	positions = [positions[:, i, :] for i in 1:length(idx)]
	velocities = [velocities[:, i, :] for i in 1:length(idx)]
	accelerations = [accelerations[:, i, :] for i in 1:length(idx)]
	Φs_ext = lguys.extract(out, :Φs_ext, idx)
	Φs = lguys.extract(out, :Φs, idx)


end

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
begin
	rs = lguys.calc_r.(positions)
	vs = lguys.calc_r.(velocities)
	accs = lguys.calc_r.(accelerations)
end

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
begin
	p = plot(out.times, rs)
	hline!([peris[idx], apos[idx]])
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.plot_xyz(positions...)

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities...)

# ╔═╡ 127c1123-29cf-45a7-b7ab-7f57b4d6e58f
lguys.plot_xyz(accelerations...)

# ╔═╡ ad078920-225d-436e-835b-d87a9db53c49
scatter(rs, vs, marker_z=out.times)

# ╔═╡ 8f8ea990-6ae8-45c2-a4f5-b7fc7a425b1d
scatter(rs, accs, marker_z=out.times)

# ╔═╡ 5d5f719c-e5b7-4e05-94b0-71523da46b66
begin 
	plots[] = [plot()]
	for i in 1:length(idx)
		scatter!(positions[i][1, :], accs[i])
	end
	plots[]
end 

# ╔═╡ 09bbae0d-ca3e-426d-b77c-69dd68ca42cc
begin
	M_b = 115
	r_b = 20
	
	a_exp(r) =  M_b / (r+r_b)^2
	phi_exp(r) = -M_b / (r+r_b)
end

# ╔═╡ fa4e7992-1ac6-4d71-a923-8b3cf81d0030
begin 
	plots[] = [plot()]
	for i in 1:length(idx)
		scatter!(out.times, accs[i] .- a_exp.(rs[i]))
	end
	plots[]
end 

# ╔═╡ 35f0ea14-a945-4745-910c-365b730676c5
begin 
	plots[] = [plot()]
	
	for i in 1:length(idx)
		scatter!(out.times, Φs_ext[i, :] .- phi_exp.(rs[i]))
	end
	plots[]
end

# ╔═╡ 2e7c1798-4066-4c46-b5ed-732263728ac0
out[14]

# ╔═╡ f6b27164-ee7c-428b-aefb-75e89d178f3e
lguys.scatter_xyz(snap.positions)

# ╔═╡ 5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
lguys.scatter_xyz(-snap.velocities .* lguys.V0)

# ╔═╡ 2dfe9a85-6553-4632-81e0-33c148fd1102
reverse(out.times)

# ╔═╡ ca334fc0-3840-4182-8bbf-b78375bb7ed9
positions

# ╔═╡ 519a88f0-8e2d-4c09-83e0-3cc2ee147e35
function get_initial_t(j)
	i = idx[j]
	N = length(rs[j])
	t_ini = -1
	for t in N-2:-1:1
		if diff(rs[j])[t] >= 0 && diff(rs[j])[t+1] <= 0
			t_ini = t
			break
		end
	end
	return t_ini
end
	

# ╔═╡ 5f45e7c7-e447-48bf-ade4-38f516df2dad
for i in 1:length(idx)
	fname = "orbit$i.csv"
	vel = -reverse(velocities[i], dims=2)
	pos = reverse(positions[i], dims=2)
	t0 = size(pos, 2) - get_initial_t(i)
	println(t0)
	t = out.times[end] .- reverse(out.times)

	df = DataFrame(
		t = t[t0:end],
		x = pos[1, t0:end],
		y = pos[2, t0:end],
		z = pos[3, t0:end],
		v_x = vel[1, t0:end],
		v_y = vel[2, t0:end],
		v_z = vel[3, t0:end],
	)

	println("saving to $fname")
	CSV.write(fname, df)
end

# ╔═╡ 5316884b-3971-4ca7-9106-f638241d3388
get_initial_t(1)

# ╔═╡ 34b812d2-21c0-4983-9d14-7dbef08ab670
begin 
	plot()
	plot!(rs[1])
	hline!(apos[1:1])
end

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(idx)
		t = get_initial_t(i)
		@printf "orbit: \t\t %i\n" i
		@printf "index: \t\t %i\n" idx[i]
		
		@printf "pericentre:\t %0.1f\n" peris[idx[i]]
		@printf "apocentre: \t %0.1f\n" apos[idx[i]]

		@printf "time of first apocentre: %f \n" out.times[end] - out.times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -velocities[i][:, t]...
		@printf "final position: [%f, %f, %f]\n" positions[i][:, 1]...
		@printf "final velocity: [%f, %f, %f]\n" -velocities[i][:, 1]...

		o = observations[idx[i]]
		@printf "ra: \t %f\n" o.ra
		@printf "dec: \t %f\n" o.dec
		@printf "pmra: \t %f\n" o.pm_ra
		@printf "pmdec: \t %f\n" o.pm_dec
		@printf "dist: \t %f\n" o.distance
		@printf "rv: \t %f\n" o.radial_velocity

		println()
	end
end

# ╔═╡ Cell order:
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
# ╠═ed6ec0e9-5734-4569-b214-2c86c22d3c55
# ╠═cd9dea36-93f8-4461-85b1-87030d9367bb
# ╠═9c7e5bb4-8db0-4527-a8ec-331e2aed958b
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
# ╠═fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
# ╠═e10615ab-0ece-4ebb-81d0-47ffdd6ea6c5
# ╠═9fad353f-751e-4577-b704-aa2eadeb0969
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═92b28dd8-c353-4959-b181-a843367b3223
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═5a5ca70c-286f-4325-ac4a-143713a844c4
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
# ╟─a644a33e-57a4-4af4-98be-1e0e84948069
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═127c1123-29cf-45a7-b7ab-7f57b4d6e58f
# ╠═ad078920-225d-436e-835b-d87a9db53c49
# ╠═8f8ea990-6ae8-45c2-a4f5-b7fc7a425b1d
# ╠═5d5f719c-e5b7-4e05-94b0-71523da46b66
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╠═fa4e7992-1ac6-4d71-a923-8b3cf81d0030
# ╠═35f0ea14-a945-4745-910c-365b730676c5
# ╟─2e7c1798-4066-4c46-b5ed-732263728ac0
# ╠═f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╠═5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╠═2dfe9a85-6553-4632-81e0-33c148fd1102
# ╠═5f45e7c7-e447-48bf-ade4-38f516df2dad
# ╠═ca334fc0-3840-4182-8bbf-b78375bb7ed9
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═5316884b-3971-4ca7-9106-f638241d3388
# ╠═34b812d2-21c0-4983-9d14-7dbef08ab670
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
                                    