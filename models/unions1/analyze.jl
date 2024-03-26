### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	using Revise
	using Plots; gr()
	using CSV, DataFrames

	# this is my codes. see my notes
	import LilGuys as lguys 
end

# ╔═╡ d4a521fb-b065-4717-9a87-1fe43ceeba2a
using Arya #only sets styles to be nice...

# ╔═╡ 6af8db7f-016e-4a48-9df0-1a4721db048e
md"""
# Assignment 6
Daniel Boyea | Astronomy 507 | March 22

Note that this assignment 

Sorry -- this assignment is very much related to my masters project, so I have extensive code to already help with the analysis (similar to what Rapha has provided). My complete set of codes and analysis is in a [github repository](https://github.com/aeyobd/dwarfs). I have liked using Julia instead of python for this work as it performs much better at certain tasks. I understand if this doesn't work as well, but I figured that it would be better to just use the tools I already had setup so I could run it more easily on Cosma and not use Sangenovese so other people had more access. 


I am more than happy to answer any questions or clarify anything as needed.
I have also found it much easier to work in julia with the simulation outputs
"""

# ╔═╡ 6964cbd2-9ef9-4682-9ac1-80c10518c374
begin 
	# load in some of the simulation outputs (arrays of snapshots)
	filenames = ["orbit1", "orbit2", "orbit4", "orbit5", "orbitc"]
	outputs = lguys.Output[]
	labels = []
	for file in filenames
		out = lguys.Output("$file/out")
		push!(labels, sum(out[1].masses) * lguys.M0)
		push!(outputs, out)
	end
	labels = round.(labels)

	out1c = lguys.Output("orbitc/out")
	out_r = lguys.Output("rapha/snapshots/")
	out_m = outputs[1]
end

# ╔═╡ 9f74dc37-f2e2-40ad-a104-97417bbae487
cen5 = lguys.ss_centre(outputs[end], percen=90)

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# calculate the centres for each snapshot
	centres = [lguys.ss_centre(out) for out in outputs]
	cens1c = lguys.ss_centre(out1c)
	cens_r = lguys.ss_centre(out_r)

	cen_m = centres[1]

end
  ╠═╡ =#

# ╔═╡ 9b666249-bc63-4515-9677-ab2c57f197f4
#=╠═╡
begin 
	import LinearAlgebra: norm, ⋅, ×

	a0 = lguys.mean(lguys.calc_r(outputs[1][1].accelerations))
	v_circ0 = sqrt(lguys.calc_r(centres[1][1].x_c) * a0) 
	x0 = centres[1][1].x_c
	v0  = centres[1][1].v_c
	# the tangental velocity by subtracting the component in the direction we need
	vt = v0 .- v0 ⋅ x0  * x0 / norm(x0)^2

	v_ini = vt / norm(vt) * v_circ0
	println(v_ini * lguys.V0)
end
  ╠═╡ =#

# ╔═╡ b7629fac-2357-4e62-95ef-2ea12e92f335
#=╠═╡
begin 
	plot(xlabel="time / Myr", ylabel="radius / kpc")
	r_c_1 = [lguys.calc_r(cen.x_c) for cen in centres[1]]
	r_c_1c = [lguys.calc_r(cen.x_c) for cen in cens1c]
	r_c_r = [lguys.calc_r(cen.x_c) for cen in cens_r]

	plot!(outputs[1].times * lguys.T0 / 1e6, r_c_1, label="given orbit")
	plot!(out1c.times * lguys.T0 / 1e6, r_c_1c, label="circular orbit")
	#plot!(out_r.times * lguys.T0 / 1e6, r_c_r, label="rapha", ls=:dash)

end
  ╠═╡ =#

# ╔═╡ 4c3e06cb-eddf-4c87-a66b-7d2aca8319c9
lguys.T0

# ╔═╡ d5a4ccf5-2e71-4101-b64f-02ef2052578a
#=╠═╡
begin 
	fs = ["orbit1/", "orbit2/", "orbit4/", "orbit5/", "orbitc/"]
	
	for i in 1:length(centres)
		xs = hcat([cen.x_c for cen in centres[i]]...)
		vs = hcat([cen.v_c for cen in centres[i]]...)
		ts = outputs[i].times
		df = DataFrame(
			(t=ts, x=xs[1, :], y=xs[2, :], z=xs[3, :], 
			vx=vs[1, :], vy=vs[2, :], vz=vs[3, :])
		)
		CSV.write("$(fs[i])centres.csv", df)
	end
end
  ╠═╡ =#

# ╔═╡ e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
#=╠═╡
begin 
	x_cen = hcat([cen.x_c for cen in cen_m]...)
	x_cenr = hcat([cen.x_c for cen in cens_r]...)
	x_cenc = hcat([cen.x_c for cen in cens1c]...)
end
  ╠═╡ =#

# ╔═╡ 9232168c-d8f9-4bae-ba54-20fe727d98c8
#=╠═╡
begin 
	R_gcc =@. sqrt(x_cenc[1, :]^2 + x_cenc[2, :]^2)
	z_gcc = x_cenc[3, :]

	R_gc =@. sqrt(x_cen[1, :]^2 + x_cen[2, :]^2)
	z_gc = x_cen[3, :]
end
  ╠═╡ =#

# ╔═╡ 9bb83ca0-be38-4b4b-996e-dede9e3fd5c7
#=╠═╡
begin 
	plot(xlabel="R / kpc", ylabel="z / kpc", aspect_ratio=:equal, dpi=1000, xlims=(0, 40), legend_position=:outertopright)
	plot!(R_gc, z_gc, label="fiducial")
	plot!(R_gcc, z_gcc, label="circular")

	scatter!([R_gc[1], R_gcc[1]], [z_gc[1], z_gcc[1]], label="today")
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


# ╔═╡ b6f8065a-dfd4-4cd4-a0bc-66f182928acf
1/exp(1)

# ╔═╡ c468d36c-b610-4e0c-93b5-990a240614aa
#=╠═╡
histogram(lguys.calc_r(out_m[1].positions .- 0*centres[1][1].x_c))
  ╠═╡ =#

# ╔═╡ 4f37a7d5-6ff1-4030-b926-6ba38c3bb7a1
out_m[1].positions

# ╔═╡ 71f2dde1-f487-45cd-a80e-c1e60557abf4
#=╠═╡
begin 
	plot(xlabel="time / Myr", ylabel="mass (solar masses)", yscale=:log10, ylims=(1e-1, 150), colorbar=false, legend_position=:outertopright)
	for i in 1:length(outputs)
		plot!(outputs[i].times * lguys.T0 / 1e6, calc_M_10(outputs[i], centres[i]) * lguys.M0, label="M0=$(labels[i])", line_z=labels[i], cmap=cgrad(:blues))
	end
	


	# plot!(out_m.times * lguys.T0 / 1e6, calc_M_10(out_m, cen_m) * lguys.M0, label="fiducial",)
	plot!(out1c.times * lguys.T0 / 1e6, calc_M_10(out1c, cens1c) * lguys.M0, label="circular", ls=:dash)

	#scatter!(out_r.times * lguys.T0 / 1e6, calc_M_10(out_r, cens_r) * lguys.M0, label="rapha's version")

end
  ╠═╡ =#

# ╔═╡ 4b93392b-8c87-4a4c-913a-b06ccbb42d83
#=╠═╡
begin 
	plot(xlabel="time / Gyr", ylabel="fraction left < 10 pc", colorbar=false, legend_position=:outertopright)
	for i in 1:length(outputs)
		plot!(outputs[i].times * lguys.T0 / 1e9, calc_M_10(outputs[i], centres[i]) * lguys.M0 ./ labels[i], label="M0=$(labels[i])", line_z=labels[i], cmap=cgrad(:blues))
	end
	hline!([1/exp(1)], label="", color=:black, ls=:dot)
	plot!(out1c.times * lguys.T0 / 1e9, calc_M_10(out1c, cens1c) * lguys.M0 / 16, label="circular", ls=:dash)
end
  ╠═╡ =#

# ╔═╡ 2f310239-681e-4574-9fc4-7ca9d06b5187
1e9/lguys.T0 * 13

# ╔═╡ 693bcb86-29b8-4dd6-989b-641af3ff4b42
p = Ref{Plots.Plot}()

# ╔═╡ ffb50d71-082d-4e27-9b89-1d80bace6050
#=╠═╡
begin 
	# compute the initial density profiles for two models
	prof_1_i = lguys.Profile(out_m[1], cen_m[1].x_c, cen_m[1].v_c, Nr=200)
	prof_r_i = lguys.Profile(out_r[1], cens_r[1].x_c, cens_r[1].v_c, Nr=200)
end
  ╠═╡ =#

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
#=╠═╡
begin 
	p[] = plot()
	scatter!(log10.(prof_1_i.rs), log10.(prof_1_i.ρs))

	scatter!(log10.(prof_r_i.rs), log10.(prof_r_i.ρs))
	xlabel!("log r / kpc")
	ylabel!("log density")
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

# ╔═╡ a5f23b9f-0e47-4ab5-a4ad-3d1ae079090a
begin 
	plot()
	scatter!(out_m[1].positions[2,:],out_m[1].positions[3, :], msw=0, ms=2, alpha=0.1,  label="")
end

# ╔═╡ 6497d973-d800-4052-a9b1-23f91dc3aa9f
#=╠═╡
begin 
	anim = @animate for i in 1:10:length(outputs[4])

		snap = outputs[4][i]
		c = cen5[i].x_c
		r =  c[:]
		r_hat = r ./ norm(r)
		
		plot(legend=false, grid=false, axis=false, dpi=100, xlims=(-0.3, 0.3), ylims=(-0.3, 0.3))
		scatter!(snap.positions[2, :] .- c[2], snap.positions[3, :] .- c[3], 
			ms=1, msw=0, ma=0.1)

		x_ref = -[0.2, 0.3] * r_hat[2]
		y_ref = -[0.2, 0.3] * r_hat[3]

		plot!(x_ref, y_ref, color=:black)
		

	end
	
	gif(anim, "unions.gif", fps = 5)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═6af8db7f-016e-4a48-9df0-1a4721db048e
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═d4a521fb-b065-4717-9a87-1fe43ceeba2a
# ╠═6964cbd2-9ef9-4682-9ac1-80c10518c374
# ╠═9b666249-bc63-4515-9677-ab2c57f197f4
# ╠═9f74dc37-f2e2-40ad-a104-97417bbae487
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═b7629fac-2357-4e62-95ef-2ea12e92f335
# ╠═4c3e06cb-eddf-4c87-a66b-7d2aca8319c9
# ╠═d5a4ccf5-2e71-4101-b64f-02ef2052578a
# ╠═e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
# ╠═9232168c-d8f9-4bae-ba54-20fe727d98c8
# ╠═9bb83ca0-be38-4b4b-996e-dede9e3fd5c7
# ╠═ceaaac14-1552-45f1-a7a6-a9ccdd66533c
# ╠═3be7f9ab-8217-4258-8dcf-d4c084b56be8
# ╠═7a22c06e-514e-414c-850d-d8de1b3eb1f9
# ╠═29c303ee-bd50-43d3-83af-11ca8a28e171
# ╠═b6f8065a-dfd4-4cd4-a0bc-66f182928acf
# ╠═c468d36c-b610-4e0c-93b5-990a240614aa
# ╠═4f37a7d5-6ff1-4030-b926-6ba38c3bb7a1
# ╠═71f2dde1-f487-45cd-a80e-c1e60557abf4
# ╠═4b93392b-8c87-4a4c-913a-b06ccbb42d83
# ╠═2f310239-681e-4574-9fc4-7ca9d06b5187
# ╠═693bcb86-29b8-4dd6-989b-641af3ff4b42
# ╠═ffb50d71-082d-4e27-9b89-1d80bace6050
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═28c6cb7f-4434-4def-b4ba-a1e0d9ca2da9
# ╠═a5f23b9f-0e47-4ab5-a4ad-3d1ae079090a
# ╠═6497d973-d800-4052-a9b1-23f91dc3aa9f
