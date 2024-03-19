### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	import LilGuys as lguys # this is my codes. see 
	using Plots;
	using CSV, DataFrames
	using LaTeXStrings
	import Arya
end

# ╔═╡ 6af8db7f-016e-4a48-9df0-1a4721db048e
md"""
# Assignment 6

Sorry -- this assignment is very much related to my masters project, so I have extensive code to already help with the analysis (similar to what Rapha has provided). My complete set of codes and analysis is in a [github repository](https://github.com/aeyobd/dwarfs)
I have also found it much easier to work in julia with the simulation outputs
"""

# ╔═╡ 6964cbd2-9ef9-4682-9ac1-80c10518c374
begin 
	out1 = lguys.Output("orbit1/out")
	out2 = lguys.Output("orbit2/out")
	out3 = lguys.Output("orbit3/out")

	out1c = lguys.Output("orbitc/out")

	out_r = lguys.Output("/cosma/home/durham/dc-boye1/unions1_rapha/snapshots/")
end

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	cens1 = lguys.ss_centre(out1)
	cens2 = lguys.ss_centre(out2)
	cens3 = lguys.ss_centre(out3)
	cens1c = lguys.ss_centre(out1c)
	cens_r = lguys.ss_centre(out_r)

	prof_1_i = lguys.Profile(out1[1], cens1[1].x_c, cens1[1].v_c, Nr=200)
	prof_r_i = lguys.Profile(out_r[1], cens_r[1].x_c, cens_r[1].v_c, Nr=200)

end
  ╠═╡ =#

# ╔═╡ b7629fac-2357-4e62-95ef-2ea12e92f335
#=╠═╡
begin 
	plot(xlabel="time / Myr", ylabel="radius / kpc")
	r_c_1 = [lguys.calc_r(cen.x_c) for cen in cens1]
	r_c_1c = [lguys.calc_r(cen.x_c) for cen in cens1c]
	r_c_r = [lguys.calc_r(cen.x_c) for cen in cens_r]

	plot!(out1.times * lguys.T0 / 1e6, r_c_1, label="me")
	plot!(out1c.times * lguys.T0 / 1e6, r_c_1c, label="circular orbit")
	plot!(out_r.times * lguys.T0 / 1e6, r_c_r, label="rapha", ls=:dash)

end
  ╠═╡ =#

# ╔═╡ e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
#=╠═╡
begin 
	x_cen = hcat([cen.x_c for cen in cens1]...)
	x_cenr = hcat([cen.x_c for cen in cens_r]...)
	x_cenc = hcat([cen.x_c for cen in cens1c]...)

	v_cen = hcat([cen.v_c for cen in cens1]...)
end
  ╠═╡ =#

# ╔═╡ ceaaac14-1552-45f1-a7a6-a9ccdd66533c
#=╠═╡
lguys.plot_xyz(x_cen, x_cenc, legend=false)[1]
  ╠═╡ =#

# ╔═╡ 3be7f9ab-8217-4258-8dcf-d4c084b56be8
#=╠═╡
lguys.plot_xyz(x_cen, x_cenc, legend=false)[2]
  ╠═╡ =#

# ╔═╡ 7a22c06e-514e-414c-850d-d8de1b3eb1f9
#=╠═╡
lguys.plot_xyz(x_cen,x_cenc, legend=false)[3]
  ╠═╡ =#

# ╔═╡ 29c303ee-bd50-43d3-83af-11ca8a28e171
function calc_M_10(out, cens)
	M_10 = Vector{Float64}(undef, length(out))
	r0 = 10/1e3 # 10 pc

	for i in 1:length(out)
		M_10[i] = sum(out[i].masses .* (lguys.calc_r(out[i].positions .- cens[i].x_c) .< r0))
	end
	return M_10
end


# ╔═╡ 71f2dde1-f487-45cd-a80e-c1e60557abf4
#=╠═╡
begin 
	plot(xlabel="time / Myr", ylabel="mass (solar masses)", yscale=:log10)
	scatter!(out1.times * lguys.T0 / 1e6, calc_M_10(out1, cens1) * lguys.M0, label="simulation")
	
	scatter!(out2.times * lguys.T0 / 1e6, calc_M_10(out2, cens2) * lguys.M0, label="simulation 2")
	
	scatter!(out3.times * lguys.T0 / 1e6, calc_M_10(out3, cens3) * lguys.M0, label="simulation 2")

	scatter!(out1c.times * lguys.T0 / 1e6, calc_M_10(out1c, cens1c) * lguys.M0, label="simulation 2")

	scatter!(out_r.times * lguys.T0 / 1e6, calc_M_10(out_r, cens_r) * lguys.M0, label="rapha's version")

end
  ╠═╡ =#

# ╔═╡ 693bcb86-29b8-4dd6-989b-641af3ff4b42
p = Ref{Plots.Plot}()

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
#=╠═╡
begin 
	p[] = plot()
	scatter!(log10.(prof_1_i.rs), log10.(prof_1_i.ρs))

	scatter!(log10.(prof_r_i.rs), log10.(prof_r_i.ρs))
	xlabel!("log r / kpc")
	ylabel!("log density")


	# only include bound points in profile...
end
  ╠═╡ =#

# ╔═╡ 28c6cb7f-4434-4def-b4ba-a1e0d9ca2da9
#=╠═╡
begin 
	p[] = plot()
	scatter!(log10.(prof_1_i.rs), (prof_1_i.Vs_circ) * lguys.V0)
	scatter!(log10.(prof_r_i.rs), (prof_r_i.Vs_circ) * lguys.V0)
	xlabel!("log r / kpc")
	ylabel!(" circular velocity")

end
  ╠═╡ =#

# ╔═╡ 92babde9-d830-4440-9839-8d079b28f274
#=╠═╡
begin
	rs1 = lguys.calc_r(out1[1].positions, cens1[1].x_c)
	vs1 = lguys.calc_r(out1[1].velocities, cens1[1].v_c)

	rs_r = lguys.calc_r(out_r[1].positions, cens_r[1].x_c)
	vs_r = lguys.calc_r(out_r[1].velocities, cens_r[1].v_c)

	ϕ1 = lguys.calc_radial_discrete_Φ(out1[1].masses, rs1, )
	ϕ_r = lguys.calc_radial_discrete_Φ(out_r[1].masses, rs_r, )

	ϵ1 = -ϕ1 .- 1/2 * vs1 .^ 2
	ϵ_r = -ϕ_r .- 1/2 * vs_r .^ 2
	scatter(log10.(rs1), ϵ1, ms=1)
	scatter!(log10.(rs_r), ϵ_r, ms=1)
end
  ╠═╡ =#

# ╔═╡ f5edcc58-8d85-4fbd-a032-e36cadeddbc1
#=╠═╡
begin 
	plot(xlabel="log r", ylabel="log v")
	scatter!(log10.(rs1), (vs1), ms=1)
	scatter!(log10.(rs_r), (vs_r), ms=1)
end
  ╠═╡ =#

# ╔═╡ a5f23b9f-0e47-4ab5-a4ad-3d1ae079090a
begin 
	plot()
	scatter!(out1[end].positions[1,:], out1[end].positions[2, :], msw=0, ms=2, alpha=0.1,  label="")
end

# ╔═╡ 6497d973-d800-4052-a9b1-23f91dc3aa9f
#=╠═╡
begin 
	anim = @animate for i in 1:10:length(out1)
		snap = out1[i]
		plot(legend=false, grid=false, axis=false, dpi=100, xlims=(-0.3, 0.3), ylims=(-0.3, 0.3))
		scatter!(snap.positions[2, :] .- x_cen[2, i], snap.positions[3, :] .- x_cen[3, i], 
			ms=1, msw=1, ma=0.2)

		scatter!([x_cen[2, i]], [x_cen[3,  i]], ms=1)
	end
	
	gif(anim, "unions.gif", fps = 12)
end
  ╠═╡ =#

# ╔═╡ d6fff211-18ae-4e29-af3b-5599ca36c3fa
md"""
can calculate circular velocity curve
"""

# ╔═╡ 9b666249-bc63-4515-9677-ab2c57f197f4
#=╠═╡
begin 
	import LinearAlgebra: norm, ⋅, ×

	a0 = lguys.mean(lguys.calc_r(out1[1].accelerations))
	v_circ0 = sqrt(lguys.calc_r(cens1[1].x_c) * a0) 
	x0 = cens1[1].x_c
	v0  = cens1[1].v_c
	# the tangental velocity by subtracting the component in the direction we need
	vt = v0 .- v0 ⋅ x0  * x0 / norm(x0)^2

	v_ini = vt / norm(vt) * v_circ0
	println(v_ini * lguys.V0)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─6af8db7f-016e-4a48-9df0-1a4721db048e
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═6964cbd2-9ef9-4682-9ac1-80c10518c374
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═b7629fac-2357-4e62-95ef-2ea12e92f335
# ╠═e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
# ╠═ceaaac14-1552-45f1-a7a6-a9ccdd66533c
# ╠═3be7f9ab-8217-4258-8dcf-d4c084b56be8
# ╠═7a22c06e-514e-414c-850d-d8de1b3eb1f9
# ╠═29c303ee-bd50-43d3-83af-11ca8a28e171
# ╠═71f2dde1-f487-45cd-a80e-c1e60557abf4
# ╠═693bcb86-29b8-4dd6-989b-641af3ff4b42
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═28c6cb7f-4434-4def-b4ba-a1e0d9ca2da9
# ╠═92babde9-d830-4440-9839-8d079b28f274
# ╠═f5edcc58-8d85-4fbd-a032-e36cadeddbc1
# ╠═a5f23b9f-0e47-4ab5-a4ad-3d1ae079090a
# ╠═6497d973-d800-4052-a9b1-23f91dc3aa9f
# ╠═d6fff211-18ae-4e29-af3b-5599ca36c3fa
# ╠═9b666249-bc63-4515-9677-ab2c57f197f4
