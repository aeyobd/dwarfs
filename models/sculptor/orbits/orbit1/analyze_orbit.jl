### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 061b1886-1878-11ef-3806-b91643300982
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ ab57edae-2292-4aef-9c1f-53802dbc0600
import YAML

# ╔═╡ 96a57df5-a7b7-447a-a4a6-2b05e391a5c6
obs_today = YAML.load_file("properties.yml")

# ╔═╡ b75f0fb1-be59-416c-a61f-4109bada9ae9
r_h = 0.11

# ╔═╡ a609f221-0721-4b4b-a393-49b386393c66
obs_today["ra_err"] = r_h 

# ╔═╡ 68d805a4-c5eb-4f2c-ba10-48c2a53f2874
obs_today["dec_err"] = r_h 

# ╔═╡ b250bf10-c228-4b14-938a-35561ae871d7
begin 
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
begin 
	orbit_expected = CSV.read("orbit.csv", DataFrame)
	x_cen_exp = transpose(hcat(orbit_expected.x, orbit_expected.y, orbit_expected.z))
	v_cen_exp = -transpose(hcat(orbit_expected.v_x, orbit_expected.v_y, orbit_expected.v_z))

end

# ╔═╡ a1c992c6-ad12-4968-b105-adfa1f327e76
lguys.plot_xyz(x_cen, x_cen_exp, labels=["n body", "point particle"])

# ╔═╡ 5255c605-56ea-4eb3-bd20-5134e3a96705
lguys.plot_xyz(v_cen, v_cen_exp, units=" / km/ s")

# ╔═╡ aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="R / kpc", ylabel="z / kpc",
		aspect=DataAspect()
	)
	x = x_cen[1, :]
	y = x_cen[2, :]
	z = x_cen[3, :]
	R = @. sqrt(x^2 + y^2)
	lines!(R, z)
	fig
end

# ╔═╡ f88b909f-c3dc-41e0-bdb1-25e229964d27
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = lguys.calc_r(x_cen)
	lines!(cens.t * lguys.T0, r)
	lines!(lguys.T0*(orbit_expected.t .- orbit_expected.t[begin]), lguys.calc_r(x_cen_exp))
	fig
end

# ╔═╡ 5e99f047-7eaa-4bd5-852a-a4beaeecec86
begin
		snap_cen = lguys.Snapshot(x_cen, v_cen, zeros(size(x_cen, 1)))
	
		obs_c = lguys.to_sky(snap_cen)
end

# ╔═╡ 89a1523b-d0c4-4b2f-b7c5-f1c4bafe01f5
getfield(obs_c[1],Symbol("ra"))

# ╔═╡ a179323f-4878-4021-b8d4-69ca733658cb
function calc_χ2s(obs_c, obs_today)
	χ2 = zeros(length(obs_c))
	for name in ["ra", "dec", "pm_ra", "pm_dec", "distance", "radial_velocity"]
		μ = obs_today[name]
		x = [getfield(o, Symbol(name)) for o in obs_c]
		σ = obs_today[name * "_err"]
		χ2 .+= @. (x -μ)^2/(σ)^2
	end
	return χ2
end

# ╔═╡ ecf7c820-81a4-4cb7-a794-b7835c77811e
χ2 = calc_χ2s(obs_c, obs_today)

# ╔═╡ cdde517a-1b3e-4d96-9156-4a8f72b795e9
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time",
		ylabel="chi2",
		yscale=log10
	)
	
	lines!(cens.t, χ2,)

	fig
end

# ╔═╡ ddd5821b-5e65-4220-831f-886b6713d026
idx_f = argmin(χ2)

# ╔═╡ 9530e936-1225-4cfc-aa9a-bf7644d612f5
r = lguys.calc_r(x_cen)

# ╔═╡ 882d4fc5-07ae-4b06-8da5-67f0894595db
import LinearAlgebra: dot

# ╔═╡ 7a30bd90-946e-418c-8339-be64c37cda76
vr = [dot(x_cen[:, i], v_cen[:, i]) / r[i] for i in 1:size(x_cen, 2)]

# ╔═╡ 0c69519c-9650-46b9-89f9-cc37227f5b1a
v_cen

# ╔═╡ d95c457b-c9be-4570-bc90-b4bbb7de56e2
plot(cens.t, vr)

# ╔═╡ 5ee77968-567a-49c4-8e19-4d88bb930558
function find_last_peri(vr, idx_f; minima=true)
	v_sign = sign(vr[idx_f])
	decreasing = vr[idx_f] < 0

	for i in idx_f-1:-1:1
		decreasing_new = vr[i] < 0
		
		if !decreasing && decreasing_new
			decreasing = decreasing_new
			if minima
				return i
			end
		elseif decreasing && !decreasing_new
			decreasing = decreasing_new
			if !minima
				return i
			end
		end
	end
	return -1
end

# ╔═╡ efcbae60-cf7c-4e74-aae4-39d19b74b6fa
idx_peri = find_last_peri(vr, idx_f)

# ╔═╡ 809095ae-0aee-47ec-9d54-9d314e2cc11d
idx_apo = find_last_peri(vr, idx_f, minima=false)

# ╔═╡ c3d6b68e-1b2a-4f84-a620-fb8dbe02867a
idx_anteperi = find_last_peri(vr, idx_apo, minima=true)

# ╔═╡ 7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
begin
	peri_filt = idx_f-200:idx_f
	t_last_peri_arg = argmin(r[peri_filt])
	t_last_peri = cens.t[peri_filt[t_last_peri_arg]] * lguys.T0
	delta_t_peri = cens.t[idx_f] * lguys.T0 - t_last_peri
end

# ╔═╡ 04d29fcb-70a0-414b-a487-7a18c44b9d58
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = lguys.calc_r(x_cen)
	lines!(cens.t * lguys.T0, r)
	scatter!(cens.t[idx_f] * lguys.T0, r[idx_f], 
		label="adpoted end", marker=:rect
	)
	scatter!(cens.t[idx_f] * lguys.T0, lguys.calc_r(x_cen_exp)[end], 
		marker=:+, markersize=10, label="expected"
	)
	
	scatter!(cens.t[idx_peri] * lguys.T0, r[idx_peri], zorder=10,
		label="last pericentre"
	)

	scatter!(cens.t[idx_apo] * lguys.T0, r[idx_apo], zorder=10,
		label="last apocentre"
	)
	
	scatter!(cens.t[idx_anteperi] * lguys.T0, r[idx_anteperi], zorder=10,
		label="last last pericentre"
	)
	
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ f5539446-eefc-498d-b422-19c47f94a6f1
available_marker_symbols()

# ╔═╡ af8a50bd-e761-4439-9fc9-80048c264d5b
begin 
	t_peri = lguys.T0 * cens.t[idx_peri]
	t_apo = lguys.T0 * cens.t[idx_apo]
	t_anteperi = lguys.T0 * cens.t[idx_anteperi]

	r_peri = r[idx_peri]
	r_apo = r[idx_apo]
end

# ╔═╡ 73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
t_f = cens.t[idx_f] * lguys.T0

# ╔═╡ 14eebce8-04f7-493b-824a-7808c7fa35dd
md"""
# validating today
"""

# ╔═╡ cdabdc7d-76a1-45f5-b83a-2454576d3964
let
	for (x, y) in [("ra", "dec"), ("pm_ra", "pm_dec"), ("distance", "radial_velocity")]
		fig = Figure()
		ax = Axis(fig[1,1],
			xlabel=x,
			ylabel=y
		)
	
		idx = idx_f-10:idx_f+10
	
		xs = [getfield(o, Symbol(x)) for o in obs_c]
		ys = [getfield(o, Symbol(y)) for o in obs_c]
		scatter!(xs[idx], ys[idx], color=log10.(χ2[idx]))
		
		errscatter!([obs_today[x]], [obs_today[y]],
			xerr=[obs_today[x * "_err"]], yerr=[obs_today[y * "_err"]]
		)
	
	
		@info fig
	end
end

# ╔═╡ 7f2305c1-3431-4e0d-801f-a89afbc3b55b
md"""
note, these plots are not actually the motion across the sky, but the position on the sky by adopting different $t_0$ values. TODO: implement a static galactocentric frame which correctly transforms into a stationary solar observer, but not super necessary
"""

# ╔═╡ 76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
begin
	
	orbital_properties = Dict(
		"pericentre" => r_peri,
		"apocentre" => r_apo,
		"period" => t_peri - t_anteperi,
		"idx_f" => idx_f,
		"idx_peri" => idx_peri,
		"t_last_peri" => t_f - t_peri
	)

	YAML.write_file("obital_properties.yml", orbital_properties)

	println("saved")
	orbital_properties
end

# ╔═╡ Cell order:
# ╠═061b1886-1878-11ef-3806-b91643300982
# ╠═ab57edae-2292-4aef-9c1f-53802dbc0600
# ╠═96a57df5-a7b7-447a-a4a6-2b05e391a5c6
# ╠═a609f221-0721-4b4b-a393-49b386393c66
# ╠═68d805a4-c5eb-4f2c-ba10-48c2a53f2874
# ╠═b75f0fb1-be59-416c-a61f-4109bada9ae9
# ╠═b250bf10-c228-4b14-938a-35561ae871d7
# ╠═bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
# ╠═a1c992c6-ad12-4968-b105-adfa1f327e76
# ╠═5255c605-56ea-4eb3-bd20-5134e3a96705
# ╠═aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
# ╠═f88b909f-c3dc-41e0-bdb1-25e229964d27
# ╠═5e99f047-7eaa-4bd5-852a-a4beaeecec86
# ╠═89a1523b-d0c4-4b2f-b7c5-f1c4bafe01f5
# ╠═a179323f-4878-4021-b8d4-69ca733658cb
# ╠═ecf7c820-81a4-4cb7-a794-b7835c77811e
# ╠═cdde517a-1b3e-4d96-9156-4a8f72b795e9
# ╠═ddd5821b-5e65-4220-831f-886b6713d026
# ╠═9530e936-1225-4cfc-aa9a-bf7644d612f5
# ╠═882d4fc5-07ae-4b06-8da5-67f0894595db
# ╠═7a30bd90-946e-418c-8339-be64c37cda76
# ╠═0c69519c-9650-46b9-89f9-cc37227f5b1a
# ╠═d95c457b-c9be-4570-bc90-b4bbb7de56e2
# ╠═5ee77968-567a-49c4-8e19-4d88bb930558
# ╠═efcbae60-cf7c-4e74-aae4-39d19b74b6fa
# ╠═809095ae-0aee-47ec-9d54-9d314e2cc11d
# ╠═c3d6b68e-1b2a-4f84-a620-fb8dbe02867a
# ╠═7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
# ╠═04d29fcb-70a0-414b-a487-7a18c44b9d58
# ╠═f5539446-eefc-498d-b422-19c47f94a6f1
# ╠═af8a50bd-e761-4439-9fc9-80048c264d5b
# ╠═73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
# ╠═14eebce8-04f7-493b-824a-7808c7fa35dd
# ╠═cdabdc7d-76a1-45f5-b83a-2454576d3964
# ╠═7f2305c1-3431-4e0d-801f-a89afbc3b55b
# ╠═76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
