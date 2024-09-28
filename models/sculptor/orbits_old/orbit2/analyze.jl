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
V_maxs = [maximum(prof.Vs_circ) for prof in profs]

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
	savefig("v_max.pdf")
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

# ╔═╡ 8e3e6607-b7de-4de6-a616-927488588579
lguys.plot_xyz(x_cen, transpose(hcat(df.x, df.y, df.z)))[3]

# ╔═╡ 79c0ae44-53c1-41e8-99d1-c20a5b4ae8cf
df.x

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
# ╠═e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
# ╠═365542d6-ed5d-4c9d-9c00-5bdf5a3bfe75
# ╠═6497d973-d800-4052-a9b1-23f91dc3aa9f
# ╠═8e3e6607-b7de-4de6-a616-927488588579
# ╠═a71bb30e-4a4c-4e94-b318-a426e3ee3045
# ╠═79c0ae44-53c1-41e8-99d1-c20a5b4ae8cf
