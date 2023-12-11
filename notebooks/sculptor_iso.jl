### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 06484816-d7f4-4389-b0bf-e5a95214082e
begin
	using Pkg
	Pkg.activate()
	using Plots
	plotlyjs()

	using LilGuys
	using Optim
	using StatsBase: mean, std
end

# ╔═╡ ddd5502b-5c4a-4608-93c8-511101dd5dfe
out = Output("/cosma/home/durham/dc-boye1/data/models/sculptor/isolation/out/")

# ╔═╡ 69b291b2-2589-489e-9a9b-19bdbc12114f
begin
snap_i = out[1, :]
snap_f = out[end, :]
end

# ╔═╡ efaf3eda-6de3-4b36-bbea-e4002b6c3d6f
maximum(snap_i.Φ)

# ╔═╡ 41bec3db-f68e-480c-8eb2-8245a6d0f8ef
scatter(snap_i.pos[1, :], snap_i.pos[2, :])

# ╔═╡ d86be9ed-2010-4d56-a031-57acaba94921
scatter(snap_f.pos[1, :], snap_f.pos[2, :])

# ╔═╡ 943d7709-4082-4f99-b17e-701fcb82edab
snap_i.index

# ╔═╡ 13c83390-a33a-4d92-9cc8-57f552b70446
begin 
	p0 = mean(snap_f.pos, dims=2)
	v0 = mean(snap_f.vel, dims=2)
	
	rs = LilGuys.r(snap_f.pos .- p0)
	vs = LilGuys.r(snap_f.vel .- v0)
	# vs_a = LilGuys.get_r(Matrix(transpose(reshape(snap_i.vel, (10_000, 3)))))
	# rs_a = LilGuys.get_r(Matrix(transpose(reshape(snap_i.pos, (10_000, 3)))))

	scatter(log10.(rs), vs)
end

# ╔═╡ 8bebe29d-0157-490e-bf4c-bc847dc1c442
idx = argmax(rs)

# ╔═╡ 7d09ea26-c5fb-432f-8b21-b4983db0b13a
vs1 = LilGuys.r(snap_f.vel .- snap_f.vel[idx])

# ╔═╡ e80ffd4b-8a4d-4e33-9f0e-e0f1f41a99ae
histogram(vs1)

# ╔═╡ 781dd39a-aca8-4f46-bc15-95a10d766edf
rs1 = LilGuys.r(snap_f.pos .- snap_f.pos[idx])

# ╔═╡ 9d688a4c-f825-47f9-9a0f-f1a233c09835
histogram(log10.(rs1))

# ╔═╡ 6aaa45fb-1610-421e-932f-e0a97613581f
std(LilGuys.r(snap_f.pos .- p0)) / length(rs)^0.5

# ╔═╡ ba1182dd-7b42-4483-9044-66c7c4fe36e7
std(LilGuys.r(snap_f.vel .- v0)) / length(vs)^0.5

# ╔═╡ e8846565-0ae8-4e0d-a77e-8176225713e5
snap_f.pos[:, idx]

# ╔═╡ ecff81d3-38ff-4635-bab2-e159e5bbb5bf
snap_f.vel[:, idx]

# ╔═╡ b983472b-3310-447f-b017-ff490a8eb156
snap_f.acc

# ╔═╡ a8c405cc-a504-46ca-b137-aa203bd9f7fd
begin
	v = LilGuys.r(snap_f.vel[idx:idx] .- v0)[1]
	pot = snap_f.Φ[idx]
	e = pot + 0.5v^2
	e_l = pot + 0.5*(v-2*0.075)^2
	e_h = pot + 0.5*(v+2*0.075)^2
end

# ╔═╡ b48bdbe5-968e-40f5-a899-699c98764eec
histogram(snap_f.Φ)

# ╔═╡ 63c0ddfb-e988-4621-8091-de20bc583f8b
scatter(rs, snap_f.Φ)

# ╔═╡ 14a4c719-4820-45e7-8c5b-0d15eee3cdd7
ϵs = LilGuys.E_spec(snap_f)

# ╔═╡ f054630a-0cd2-429e-991b-455c0e06a0df
histogram(ϵs)

# ╔═╡ 16a6cabe-5008-4ce1-96bf-f1271e422a1d
dxs, dvs = LilGuys.knn_err(snap_f)

# ╔═╡ f4588f02-8abd-4ebc-8eaf-1924fdc47669
function circle(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

# ╔═╡ 6324702c-7ad0-408f-a64b-b3a574b07a89
begin
	circles = Plots.Shape[]
	snap = snap_f
	k = 100
	rs2 = LilGuys.r(snap.pos)
	r_cut = sort(rs2)[k]
	filt = rs .< r_cut
	snap_filt = snap_f[filt]
	for i in 1:k
		x = snap_filt.pos[1, i]
		y = snap_filt.pos[2, i]
		z = snap_filt.pos[3, i]
		if abs(z) > 10
			continue
		end
		r = dxs[i]
		push!(circles, circle(x, y, r))
	end

	kwargs = (fa=0.05, color=:black, lw=0, legend=false,  grid=false, aspect_ratio=:equal)
	plot(circles; kwargs...)
end

# ╔═╡ 52646f53-f055-430d-a53c-9362cac35b86
begin
	circles2 = Plots.Shape[]
	for i in 1:k
		x = snap_filt.vel[1, i]
		y = snap_filt.vel[2, i]
		z = snap_filt.vel[3, i]
		r = LilGuys.r(snap_filt.pos[:, i])

		r = dvs[i]
		push!(circles2, circle(x, y, r))
	end

	kwargs2 = (fa=0.05, color=:black, lw=0, legend=false,  grid=false, aspect_ratio=:equal)
	plot(circles2; kwargs2...)
end

# ╔═╡ 445924a7-fd07-43b4-9939-5f2bb34fc51a
snap2 = Snapshot("/cosma/home/durham/dc-boye1/data/models/fornax/mw/out/snapshot_090.hdf5")

# ╔═╡ 1fff02a4-8746-4f8f-a1c4-4fb00e450dd5
ps = LilGuys.bound_probabilities(snap2)

# ╔═╡ bb0ecbcd-2530-45c3-93e4-7b7925cce647
histogram(ps)

# ╔═╡ 6aa3a0e6-751b-4a87-97c9-2165a2ea4306
es = LilGuys.E_spec(snap2)

# ╔═╡ f42c5bee-1858-4e20-8e6e-ce4a7bb105e7
scatter(asinh.(es), ps)

# ╔═╡ 2720212c-55e0-49f6-825a-49f9e6ea0eef
length(snap2)

# ╔═╡ 0f006653-9f09-49b7-b19b-d35d804c75ab
sum((0.1 .< ps) .&& (0.9 .> ps))

# ╔═╡ 34250926-90dc-459f-9148-ea4043a172d1
scatter(snap2.pos[2, :], snap2.pos[1, :], lw=0, ms=1, ma=0.05)

# ╔═╡ 6ee1d6e7-4879-46b9-b05a-b50c2cc3f868
begin 
	xc, vc = LilGuys.energetic_center(snap2)
	snap2.pos .-= xc
	snap2.vel .-= vc
end

# ╔═╡ Cell order:
# ╠═06484816-d7f4-4389-b0bf-e5a95214082e
# ╠═ddd5502b-5c4a-4608-93c8-511101dd5dfe
# ╠═69b291b2-2589-489e-9a9b-19bdbc12114f
# ╠═efaf3eda-6de3-4b36-bbea-e4002b6c3d6f
# ╠═41bec3db-f68e-480c-8eb2-8245a6d0f8ef
# ╠═d86be9ed-2010-4d56-a031-57acaba94921
# ╠═943d7709-4082-4f99-b17e-701fcb82edab
# ╠═13c83390-a33a-4d92-9cc8-57f552b70446
# ╠═8bebe29d-0157-490e-bf4c-bc847dc1c442
# ╠═7d09ea26-c5fb-432f-8b21-b4983db0b13a
# ╠═e80ffd4b-8a4d-4e33-9f0e-e0f1f41a99ae
# ╠═781dd39a-aca8-4f46-bc15-95a10d766edf
# ╠═9d688a4c-f825-47f9-9a0f-f1a233c09835
# ╠═6aaa45fb-1610-421e-932f-e0a97613581f
# ╠═ba1182dd-7b42-4483-9044-66c7c4fe36e7
# ╠═e8846565-0ae8-4e0d-a77e-8176225713e5
# ╠═ecff81d3-38ff-4635-bab2-e159e5bbb5bf
# ╠═b983472b-3310-447f-b017-ff490a8eb156
# ╠═a8c405cc-a504-46ca-b137-aa203bd9f7fd
# ╠═b48bdbe5-968e-40f5-a899-699c98764eec
# ╠═63c0ddfb-e988-4621-8091-de20bc583f8b
# ╠═14a4c719-4820-45e7-8c5b-0d15eee3cdd7
# ╠═f054630a-0cd2-429e-991b-455c0e06a0df
# ╠═16a6cabe-5008-4ce1-96bf-f1271e422a1d
# ╠═f4588f02-8abd-4ebc-8eaf-1924fdc47669
# ╠═6324702c-7ad0-408f-a64b-b3a574b07a89
# ╠═52646f53-f055-430d-a53c-9362cac35b86
# ╠═445924a7-fd07-43b4-9939-5f2bb34fc51a
# ╠═1fff02a4-8746-4f8f-a1c4-4fb00e450dd5
# ╠═bb0ecbcd-2530-45c3-93e4-7b7925cce647
# ╠═6aa3a0e6-751b-4a87-97c9-2165a2ea4306
# ╠═f42c5bee-1858-4e20-8e6e-ce4a7bb105e7
# ╠═2720212c-55e0-49f6-825a-49f9e6ea0eef
# ╠═0f006653-9f09-49b7-b19b-d35d804c75ab
# ╠═34250926-90dc-459f-9148-ea4043a172d1
# ╠═6ee1d6e7-4879-46b9-b05a-b50c2cc3f868
