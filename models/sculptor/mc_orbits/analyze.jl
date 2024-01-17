### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 1ecbfa5e-8439-11ee-3a96-fb3d2dc7a2d3
begin 
	import Pkg
	Pkg.activate()
	using Plots
	plotly()
	using DataFrames, CSV

	using LilGuys

end

# ╔═╡ a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
if !@isdefined err
	include("sample.jl")
end

# ╔═╡ cd9dea36-93f8-4461-85b1-87030d9367bb
pwd()

# ╔═╡ 7371d47e-b406-414c-8738-596f9c57a945
begin
	x = Ref{Vector{Int}}()
	y = Ref{Vector{Int}}()
	z = Ref{Vector{Int}}()
	pos = Ref{Vector{Int}}()
	idx = Ref{Vector{Int}}()

end

# ╔═╡ 9c7e5bb4-8db0-4527-a8ec-331e2aed958b
begin
	filename = "out"
	out = Output(filename);
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
	hist_plots = []
	p1 = histogram(peris)
	xlabel!("pericenter")
	push!(hist_plots, p1)

	p1 = histogram(apos, bins=20)
	xlabel!("apocenter")
	push!(hist_plots, p1)

	hist_plots
end

# ╔═╡ b2004011-d651-4201-a973-3c417f354e2b
begin 
	ps = LilGuys.extract(out, :positions)
	pos1 = hcat([p[:, 84] for p in ps]...)
	
	rs = hcat([LilGuys.calc_r(p) for p in ps]...)

end

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
begin
	idx[] = [2]
	plot(out.times, rs[idx[][1], :])
	hline!([peris[idx[]], apos[idx[]]])
end

# ╔═╡ 45a1a7b4-a9e8-4f24-8100-0dc010ab7262
snap.masses

# ╔═╡ b8e526b1-e6c5-40e1-a8ce-bf8e7e846182
out.times

# ╔═╡ 567f2899-9766-4e55-bc71-3f6f17a5a127
peris

# ╔═╡ 0793fe69-11ce-4009-bf8e-4022bfbe5b1a
out[10].accelerations

# ╔═╡ 3158f1ea-e68c-4ec3-bedd-a8d6f2b5814d
pos1

# ╔═╡ 2db2ce78-d599-4a2a-a3af-54a5215263e5
begin
	xy_plots = []

	p = plot()
	for i in idx[]
		xs = hcat([p[:, i] for p in ps]...)
		
		plot!(xs[1, :], xs[2, :])
	end
	xlabel!(p, "x")
	ylabel!(p, "y")
	push!(xy_plots, p)

	p = plot()
	for i in idx[]
		xs = hcat([p[:, i] for p in ps]...)
		plot!(xs[2, :], xs[3, :])
	end	
	xlabel!(p, "y")
	ylabel!(p, "z")
	push!(xy_plots, p)

	
	p = plot()
	for i in idx[]
		xs = hcat([p[:, i] for p in ps]...)
		plot!(xs[1, :], xs[3, :])
	end
	xlabel!(p, "x")
	
	ylabel!(p, "z")
	push!(xy_plots, p)

	xy_plots
end

# ╔═╡ 92b28dd8-c353-4959-b181-a843367b3223
histogram(LilGuys.calc_r(snap.velocities))

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
histogram(LilGuys.calc_r(snap.positions))

# ╔═╡ 5a5ca70c-286f-4325-ac4a-143713a844c4
histogram(snap.Φs_ext)

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = LilGuys.calc_E_spec_kin(out[1]) + out[1].Φs_ext

# ╔═╡ 8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
maximum(ϵ)

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
observations = to_sky(snap, invert_velocity=true)

# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

# ╔═╡ 0eb5b370-56b4-4546-9849-82d80e707f0f
idx[] = [argmax(abs.(dists .- d) .< 0.1) for d in [77, 80, 85, 86, 90, 92]]

# ╔═╡ a449025f-647c-4b68-afad-cefc2941a4da
argmax(abs.(dists .- LilGuys.percentile(dists, 16)) .< 0.1)

# ╔═╡ 0460a804-a14a-40cc-9490-0a1ee7551fbb
argmax(abs.(dists .- LilGuys.percentile(dists, 60)) .< 0.1)

# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ ac81acd8-4a78-4230-bc70-3b78a861b618
begin
plots = []
	
for sym in [:distance, :pm_ra, :pm_dec, :radial_velocity]
    x = [getproperty(o, sym) for o in observations]
    p = histogram(x, normalize=:pdf)

    μ = getproperty(obs, sym)
    σ = getproperty(err, sym)
    
    x_mod = LinRange(μ - 3σ, μ + 3σ, 1000)
    y_mod = normal_dist.(x_mod, μ, σ)
    plot!(p, x_mod, y_mod, z_order=2, lw=2)
    title!(p, string(sym))

    push!(plots, p)
end

plots
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

# ╔═╡ Cell order:
# ╠═1ecbfa5e-8439-11ee-3a96-fb3d2dc7a2d3
# ╠═a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
# ╠═cd9dea36-93f8-4461-85b1-87030d9367bb
# ╠═7371d47e-b406-414c-8738-596f9c57a945
# ╠═9c7e5bb4-8db0-4527-a8ec-331e2aed958b
# ╠═4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
# ╠═0eb5b370-56b4-4546-9849-82d80e707f0f
# ╠═a449025f-647c-4b68-afad-cefc2941a4da
# ╠═0460a804-a14a-40cc-9490-0a1ee7551fbb
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═b2004011-d651-4201-a973-3c417f354e2b
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═45a1a7b4-a9e8-4f24-8100-0dc010ab7262
# ╠═b8e526b1-e6c5-40e1-a8ce-bf8e7e846182
# ╠═567f2899-9766-4e55-bc71-3f6f17a5a127
# ╠═0793fe69-11ce-4009-bf8e-4022bfbe5b1a
# ╠═3158f1ea-e68c-4ec3-bedd-a8d6f2b5814d
# ╠═2db2ce78-d599-4a2a-a3af-54a5215263e5
# ╠═92b28dd8-c353-4959-b181-a843367b3223
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═5a5ca70c-286f-4325-ac4a-143713a844c4
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
