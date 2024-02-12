### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()

	import LilGuys as lguys
	using Plots
end

# ╔═╡ 6964cbd2-9ef9-4682-9ac1-80c10518c374
out = lguys.Output("out")

# ╔═╡ dc6bef2c-e4eb-41ab-8f37-3569c068573a
begin
	cens = lguys.ss_centre(out)
	
	snap_i = out[1]
	snap_f = out[end]

	prof_i = lguys.Profile(snap_i, cens[1].x_c, cens[1].v_c)
	prof_f = lguys.Profile(snap_f, cens[end].x_c, cens[end].v_c)
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

# ╔═╡ e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
begin 
	x_cen = hcat([cen.x_c for cen in cens]...)
	v_cen = [cen.v_c for cen in cens]
end

# ╔═╡ 6497d973-d800-4052-a9b1-23f91dc3aa9f
begin 
	anim = @animate for i in 1:10:(length(out))
		snap = out[i]
		scatter(snap.positions[2, :], snap.positions[3, :], 
			ms=1, msw=0, ma=0.1, legend=false)
		xlims!(-200, 200)
		ylims!(-200, 200)
		scatter!([x_cen[2, i]], [x_cen[3, i]])
	end
	
	gif(anim, "sculptor.gif", fps = 15)
end

# ╔═╡ 5835c87b-ded1-4dd3-8ee8-333658f866ca
scatter(x_cen[2, :], x_cen[3, :])

# ╔═╡ 5899c047-b370-4a26-88bd-678bef4b0d5c
Ns = [sum(c.filt) for c in cens]

# ╔═╡ 10b92861-639e-4b35-b356-3baf9fc36b45
plot(Ns)

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═6964cbd2-9ef9-4682-9ac1-80c10518c374
# ╠═dc6bef2c-e4eb-41ab-8f37-3569c068573a
# ╠═78e2da2d-126b-442f-91d0-4cf9c1a4e82c
# ╠═693bcb86-29b8-4dd6-989b-641af3ff4b42
# ╠═dfa6a5aa-e7ff-4e8b-b249-600ca7a02bc3
# ╠═6619ffaa-2e99-4be4-ad96-b22afda64b88
# ╠═e6cbe208-2d1e-42c6-aadc-db3c2bd7c737
# ╠═6497d973-d800-4052-a9b1-23f91dc3aa9f
# ╠═5835c87b-ded1-4dd3-8ee8-333658f866ca
# ╠═5899c047-b370-4a26-88bd-678bef4b0d5c
# ╠═10b92861-639e-4b35-b356-3baf9fc36b45
