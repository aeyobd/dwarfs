### A Pluto.jl notebook ###
# v0.19.36

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

# ╔═╡ 6964cbd2-9ef9-4682-9ac1-80c10518c374
out = lguys.Output("out")

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
begin
	cens = lguys.ss_centre(out)
	
	snap_i = out[1]
	snap_f = out[end]

	prof_i = lguys.Profile(snap_i, cens[1].x_c, cens[1].v_c)
	prof_f = lguys.Profile(snap_f,  cens[end].x_c, cens[end].v_c)
end

# ╔═╡ c9d4096f-abbd-4789-b762-4c3156c57232
profs = [lguys.Profile(out[i], cens[i].x_c, cens[i].v_c) for i in 1:length(out)]

# ╔═╡ db320665-f46d-4aed-a2b2-4b39bcb605c5
begin 
	Nt = length(out)
	V_maxs = Vector{Float64}(undef, Nt)
	r_maxs = Vector{Float64}(undef, Nt)

	for i in 1:Nt
		idx = argmax(profs[i].Vs_circ)
		V_maxs[i] = profs[i].Vs_circ[idx]
		r_maxs[i] = profs[i].rs[idx]
	end
end

# ╔═╡ 78e2da2d-126b-442f-91d0-4cf9c1a4e82c
prof_i

# ╔═╡ 693bcb86-29b8-4dd6-989b-641af3ff4b42
p = Ref{Plots.Plot}()

# ╔═╡ dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
begin 
	p[] = plot()
	for prof in [prof_i, prof_f]
		scatter!(log10.(prof.rs), log10.(prof.ρs))
	end
	p[]
	xlabel!("log r / kpc")
	ylabel!("log density")


	# only include bound points in profile...
end

# ╔═╡ 6619ffaa-2e99-4be4-ad96-b22afda64b88
begin 
	p[] = plot()
	for prof in [prof_i, prof_f]
		scatter!(log10.(prof.rs), lguys.V0 * prof.Vs_circ)
	end
	p[]
	xlabel!("log r / kpc")
	ylabel!("circular velocity")

	# only include bound points in profile...
end

# ╔═╡ 245721a6-01aa-43e7-922d-ed5da02207c1
begin 
	plot(out.times * lguys.T0/1e9, V_maxs * lguys.V0)
	xlabel!("time / Gyr")
	ylabel!("V_circ max / km / s")
end

# ╔═╡ 9d587b4e-5670-4854-9c17-f5c7f7c6a57f
md"""
from EN 21, 
"""

# ╔═╡ 04df30cf-0931-4a63-a4ac-2d8409bc5f7d
begin 
	scatter(r_maxs, V_maxs*lguys.V0, label="")
	x =  LinRange(1, 0.1, 1000)
	α = 0.4
	β = 0.65
	y = @. 2^α * x^β * (1 + x^2)^(-α)
	plot!(x .* r_maxs[1], y .* V_maxs[1] * lguys.V0, lw=10, label="EN21")
	plot!(xlabel=L"$R_{\rm max}$", ylabel=L"$V_{\rm max}\ / \ {\rm km / s^{-1}}$")
end

# ╔═╡ e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
begin 
	x_cen = hcat([cen.x_c for cen in cens]...)
	v_cen = hcat([cen.v_c for cen in cens]...)
end

# ╔═╡ 365542d6-ed5d-4c9d-9c00-5bdf5a3bfe75
floor(Int, 1/2)

# ╔═╡ 6497d973-d800-4052-a9b1-23f91dc3aa9f
begin 
	anim = @animate for i in 1:10:length(out)
		snap = out[i]
		scatter(snap.positions[2, :], snap.positions[3, :], 
			ms=1, msw=0, ma=0.1, legend=false)
		
		scatter!([x_cen[2, i]], 
			[x_cen[3,  i]], ms=1)
		xlims!(-200, 200)
		ylims!(-200, 200)
	end
	
	gif(anim, "sculptor.gif", fps = 15)
end

# ╔═╡ a71bb30e-4a4c-4e94-b318-a426e3ee3045
df = CSV.read("orbit.csv", DataFrame)

# ╔═╡ 052b814a-e830-458b-aa79-533a63a89285
lguys.plot_xyz(x_cen, transpose(hcat(df.x, df.y, df.z)))


# ╔═╡ 8e3e6607-b7de-4de6-a616-927488588579
begin
	plot(x_cen[2, :], x_cen[3, :], label="n body")
	plot!(df.y, df.z, label="single particle")
	scatter!(x_cen[2, 1:1], x_cen[3,1:1], label="start")
	xlabel!("x / kpc")
	ylabel!("y / kpc")
end

# ╔═╡ f25bdc49-6b88-4bc1-906a-baff1ca7d1f7
cens[1].v_c

# ╔═╡ 51aad59f-96ff-4e88-bff2-f3ce378dd73b
prof_i.rs

# ╔═╡ 59227684-c081-4465-bd39-c0f1698c8988
begin
	radii = lguys.calc_r(snap_i.positions, cens[1].x_c)
	vels = lguys.calc_r(snap_i.velocities, cens[1].v_c)
	Φ_f = lguys.calc_radial_Φ(radii, snap_i.masses)

	ke = @. 1/2 * vels^2
	ϵ = @. -Φ_f(radii) - ke
	filt = @. ϵ > 0
end

# ╔═╡ aca25368-28c7-4f4a-bfcf-295902eaff9a
begin
	ρ_s(r, a, n) = exp( -(r/a)^(1/n))
end

# ╔═╡ 743d20cd-9636-4384-9bf2-b2d7e259ae7d
r_h = 0.308

# ╔═╡ c7aa3122-4bdd-4889-895a-d143f8919405
prof = lguys.Profile(snap_i[filt], cens[1].x_c, cens[1].v_c, Nr=200)

# ╔═╡ 7f168f93-849e-4586-a04f-6434165e6561
function gradient(x, y)
	s = sortperm(x)
	x = x[s]
	y = y[s]
	N = length(x)

	grad = Vector{Float64}(undef, N)

	grad[1] = (y[2] - y[1]) / (x[2] - x[1])
	grad[end] = (y[end] - y[end-1]) / (x[end] - x[end-1])
	for i in 2:(N-1)
		hs = x[i] - x[i-1]
		hd = x[i+1] - x[i]

		numerator = hs^2 * y[i+1] + (hd^2 - hs^2) * y[i] - hd^2*y[i-1]
		denom = hd*hs*(hd + hs)
		grad[i] = numerator/denom
	end
	return grad
end

# ╔═╡ ea62a7e0-5fda-4d0c-aedb-61f6f97c1b8f
histogram(ϵ)

# ╔═╡ 33e0615a-2d5a-4b07-808f-98e3cd9d22b2
begin 
	rs = lguys.calc_r(snap_i[filt].positions, cens[1].x_c)
	vs = lguys.calc_r(snap_i[filt].positions, cens[1].x_c)
	idx = sortperm(rs)
	rs = rs[idx]
	vs = vs[idx]

end

# ╔═╡ 821793b5-514a-4dd8-95a3-1955fb65c855


# ╔═╡ b7629fac-2357-4e62-95ef-2ea12e92f335
begin 
	plot()
	r_c = [lguys.calc_r(cen.x_c) for cen in cens]
	plot!(out.times, r_c)
end

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═6964cbd2-9ef9-4682-9ac1-80c10518c374
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═c9d4096f-abbd-4789-b762-4c3156c57232
# ╠═db320665-f46d-4aed-a2b2-4b39bcb605c5
# ╠═78e2da2d-126b-442f-91d0-4cf9c1a4e82c
# ╠═693bcb86-29b8-4dd6-989b-641af3ff4b42
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═6619ffaa-2e99-4be4-ad96-b22afda64b88
# ╠═245721a6-01aa-43e7-922d-ed5da02207c1
# ╠═9d587b4e-5670-4854-9c17-f5c7f7c6a57f
# ╠═04df30cf-0931-4a63-a4ac-2d8409bc5f7d
# ╠═e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
# ╠═365542d6-ed5d-4c9d-9c00-5bdf5a3bfe75
# ╠═6497d973-d800-4052-a9b1-23f91dc3aa9f
# ╠═052b814a-e830-458b-aa79-533a63a89285
# ╠═8e3e6607-b7de-4de6-a616-927488588579
# ╠═a71bb30e-4a4c-4e94-b318-a426e3ee3045
# ╠═f25bdc49-6b88-4bc1-906a-baff1ca7d1f7
# ╠═51aad59f-96ff-4e88-bff2-f3ce378dd73b
# ╠═59227684-c081-4465-bd39-c0f1698c8988
# ╠═aca25368-28c7-4f4a-bfcf-295902eaff9a
# ╠═743d20cd-9636-4384-9bf2-b2d7e259ae7d
# ╠═c7aa3122-4bdd-4889-895a-d143f8919405
# ╠═7f168f93-849e-4586-a04f-6434165e6561
# ╠═ea62a7e0-5fda-4d0c-aedb-61f6f97c1b8f
# ╠═33e0615a-2d5a-4b07-808f-98e3cd9d22b2
# ╠═821793b5-514a-4dd8-95a3-1955fb65c855
# ╠═b7629fac-2357-4e62-95ef-2ea12e92f335
