### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using Plots; gr()
	import LilGuys as lguys
end

# ╔═╡ 9f75b286-b021-4fa1-a29d-7051c55c0a33
if !isdefined(Main, :PlutoRunner) # if running from file
	using ArgParse
	println("running as script")
	dirname2 = get_args()
end

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 82c76c56-e874-4eba-9367-569b656155a2
pwd()

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
dirname1 = "/cosma/home/durham/dc-boye1/sculptor/isolation/1e4_s0.2/"

# ╔═╡ f645be76-6477-4970-b655-603d700a10e7
begin 
	if isdefined(Main, :PlutoRunner)
		dirname = dirname1
	else
		dirname = dirname2
	end
	cd(@__DIR__)
	cd(dirname)
	pwd()
end

# ╔═╡ 5717c18c-69cd-4740-a37a-d7ef80d68ae9
plots_dir = "$dirname/figures"

# ╔═╡ 5f7e3de9-a7fe-4217-a53c-0101d4b6314d
if dirname !== ""
	if !isdir(dirname)
		mkdir(plots_dir)
	end
	println("saving figures to $plots_dir")
else
	println("no directory specified")
end

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $dirname")
	out = lguys.Output("out")
end

# ╔═╡ ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
begin 
	Ls = hcat([lguys.calc_L_tot(snap) for snap in out]...)
	plot(transpose(Ls ./ Ls[:, 1]))
	xlabel!("snapshot")
	ylabel!("L/L_0")
end

# ╔═╡ bfa7593c-4915-4e03-83e6-8f790de4c1a5
begin 
	Es = hcat([lguys.calc_E_tot(snap) for snap in out]...)
	plot(transpose(Es ./ Es[1]))
	xlabel!("snapshot")
	ylabel!("rel change energy")
end

# ╔═╡ fe476f91-5be3-4cf9-88df-55d9568ab88f
pos = lguys.extract(out, :positions)

# ╔═╡ 4984a413-8fb7-4937-8c5c-f7e424573866
cens = lguys.ss_centre(out)

# ╔═╡ c5672da9-0dad-4d22-abe5-9e186ccde02d
begin
	snap_i = out[1]
	snap_f = out[end]
	prof_i = lguys.Profile(snap_i, cens[1].x_c, cens[1].v_c)
	prof_f = lguys.Profile(snap_f, cens[end].x_c, cens[end].v_c)
end

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
begin 
	scatter(log10.(prof_i.rs), prof_i.Vs_circ * lguys.V0)
	scatter!(log10.(prof_f.rs), prof_f.Vs_circ * lguys.V0)

	xlabel!("log r")
	ylabel!("V_circ")
end

# ╔═╡ 30f3a998-8c4b-4cc3-97d7-967537f9986f
begin 
	scatter(log10.(prof_i.rs), log10.(prof_i.ρs))
	scatter!(log10.(prof_f.rs), log10.(prof_f.ρs))
	xlabel!("log r / kpc")
	ylabel!("log rho / 1e10 Msun kpc^-3")

end

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
begin 
	histogram2d(log10.(lguys.calc_r(snap_i.positions)), lguys.calc_r(snap_i.velocities) * lguys.V0)
	xlabel!("log radius / kpc")
	ylabel!("velocity (km/s)")
end

# ╔═╡ 4e4c420e-6d35-42af-b859-b62861b1b6b9
begin 
	histogram2d(log10.(lguys.calc_r(snap_f.positions)), lguys.calc_r(snap_f.velocities) * lguys.V0)
	xlabel!("log radius / kpc")
	ylabel!("velocity (km/s)")
end

# ╔═╡ 3cf46f2f-33c1-4eb5-b2af-243a0e89d1be
begin 
	
	x_cen = [cen.x_c for cen in cens]
	v_cen = [cen.v_c for cen in cens]
end

# ╔═╡ d6fe3195-4d4d-4db8-a384-a195233dc159
plot(transpose(hcat(x_cen...)))

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
plot(transpose(hcat(v_cen...)))

# ╔═╡ dc221349-eb61-4ace-8de3-a6c50249aca0
begin 
	rs = Vector[]
	for i in 1:length(out)
		push!(rs, sort(lguys.calc_r(out[i].positions .- x_cen[i])))
	end

	rs = hcat(rs...)
end

# ╔═╡ f21cfe22-95f3-485d-902b-b022a41548c2
begin 
	plot()
	for Npoints in [100, 300, 1000, 3000]
		plot!(1:length(out), log10.(rs[Npoints, :]), label="N=$Npoints")
	end
	xlabel!("snapshot")
	ylabel!("log r contating N particles")
end

# ╔═╡ 470466ff-3709-4cdb-93c0-b5b7848a940d
# ╠═╡ disabled = true
#=╠═╡
begin 
	anim = @animate for i in 1:1:length(out)
		snap = out[i]
		xr = 10
		plot(legend=false, grid=false, axis=false, dpi=100,
			xlim=(-xr, xr), ylim=(-xr, xr))
		scatter!(snap.positions[2, :], snap.positions[3, :],
		ms=1, msw=0, ma=1)
		scatter!([0], [0], ms=2, msw=0)
	
	end
	gif(anim, "isolation.gif", fps=12)
end
  ╠═╡ =#

# ╔═╡ 3b2bb553-0130-4c8a-80ad-6e1f7071a293
lguys.plot_xyz(lguys.extract(out, :positions, 2217))[3]

# ╔═╡ f7ec3ec0-7d3f-4dcb-b508-40f389723cd1
out[1].index[sortperm(lguys.calc_r(out[1]))][3000]

# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╠═82c76c56-e874-4eba-9367-569b656155a2
# ╠═9f75b286-b021-4fa1-a29d-7051c55c0a33
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═f645be76-6477-4970-b655-603d700a10e7
# ╠═5717c18c-69cd-4740-a37a-d7ef80d68ae9
# ╠═5f7e3de9-a7fe-4217-a53c-0101d4b6314d
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═c5672da9-0dad-4d22-abe5-9e186ccde02d
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═30f3a998-8c4b-4cc3-97d7-967537f9986f
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╠═4e4c420e-6d35-42af-b859-b62861b1b6b9
# ╠═ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═fe476f91-5be3-4cf9-88df-55d9568ab88f
# ╠═4984a413-8fb7-4937-8c5c-f7e424573866
# ╠═3cf46f2f-33c1-4eb5-b2af-243a0e89d1be
# ╠═d6fe3195-4d4d-4db8-a384-a195233dc159
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═dc221349-eb61-4ace-8de3-a6c50249aca0
# ╠═f21cfe22-95f3-485d-902b-b022a41548c2
# ╠═470466ff-3709-4cdb-93c0-b5b7848a940d
# ╠═3b2bb553-0130-4c8a-80ad-6e1f7071a293
# ╠═f7ec3ec0-7d3f-4dcb-b508-40f389723cd1
