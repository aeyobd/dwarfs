### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using Plots
	plotly()
	using DataFrames, CSV

	import LilGuys as lguys
end


# ╔═╡ a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
if !@isdefined obs # read in sample (but only once)
	include("sample.jl")
end

# ╔═╡ ed6ec0e9-5734-4569-b214-2c86c22d3c55
plots = Ref{Vector{Plots.Plot}}() # so we can reuse this variable

# ╔═╡ cd9dea36-93f8-4461-85b1-87030d9367bb
pwd()

# ╔═╡ 9c7e5bb4-8db0-4527-a8ec-331e2aed958b
begin
	filename = "out"
	out = lguys.Output(filename);
end

# ╔═╡ 4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
begin 
	snap = out[1]
	snap = snap[sortperm(snap.index)]
end

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
begin
	df = CSV.read("peris_apos.csv", DataFrame)
	df_idx = sortperm(df[!, :index])
	peris = df[df_idx, :pericenter]
	apos = df[df_idx, :apocenter];
end

# ╔═╡ fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
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

# ╔═╡ a449025f-647c-4b68-afad-cefc2941a4da
findfirst(isequal(lguys.percentile(peris, 16)), peris)

# ╔═╡ d75d4521-f4d3-4162-a01f-830d6a4b85b8
findfirst(isequal(lguys.percentile(peris, 84)), peris)

# ╔═╡ 795cfff6-7951-4f2b-a956-a8f8c8c40113
idx = [1, 2]

# ╔═╡ 0b9cb669-ef83-44a3-b4a8-f1c2e07cba0b
out.times

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

# ╔═╡ f26e923c-7b92-410d-bd5b-9ab3baed4562
obs[end]

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

# ╔═╡ e5825c4a-b446-44a3-8fd5-d94664965aca
for f in fieldnames(lguys.Observation)
	a = getproperty(observations[1], f)
	b = getproperty(obs, f)
	println("error of $f:\t", (a-b)/b)
end

# ╔═╡ a644a33e-57a4-4af4-98be-1e0e84948069
begin
	gp_orbit = CSV.read("/cosma/home/durham/dc-boye1/dwarfs/notebooks/sculptor_orbits.csv", DataFrame)
end

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = lguys.extract(out, :positions, idx)
	velocities = lguys.extract(out, :velocities, idx)
	accelerations = lguys.extract(out, :accelerations, idx)
	positions = [positions[:, i, :] for i in idx]
	velocities = [velocities[:, i, :] for i in idx]
	accelerations = [accelerations[:, i, :] for i in idx]
	# Φs_ext = lguys.extract(out, :Φs_ext, idx)
	Φs = lguys.extract(out, :Φs, idx)


end

# ╔═╡ fb2f5133-0249-48cd-a4fb-e9e46002f557
plot(transpose(Φs))

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
	for i in idx
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
	for i in idx
		scatter!(out.times, accs[i] .- a_exp.(rs[i]))
	end
	plots[]
end 

# ╔═╡ 35f0ea14-a945-4745-910c-365b730676c5
begin 
	plots[] = [plot()]
	
	for i in idx
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

# ╔═╡ Cell order:
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
# ╠═ed6ec0e9-5734-4569-b214-2c86c22d3c55
# ╠═cd9dea36-93f8-4461-85b1-87030d9367bb
# ╠═9c7e5bb4-8db0-4527-a8ec-331e2aed958b
# ╠═4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
# ╠═a449025f-647c-4b68-afad-cefc2941a4da
# ╠═d75d4521-f4d3-4162-a01f-830d6a4b85b8
# ╠═795cfff6-7951-4f2b-a956-a8f8c8c40113
# ╠═0b9cb669-ef83-44a3-b4a8-f1c2e07cba0b
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═92b28dd8-c353-4959-b181-a843367b3223
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═5a5ca70c-286f-4325-ac4a-143713a844c4
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═f26e923c-7b92-410d-bd5b-9ab3baed4562
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
# ╟─a644a33e-57a4-4af4-98be-1e0e84948069
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═fb2f5133-0249-48cd-a4fb-e9e46002f557
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
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
