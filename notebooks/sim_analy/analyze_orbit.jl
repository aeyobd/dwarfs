### A Pluto.jl notebook ###
# v0.19.43

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

# ╔═╡ 8b41af50-9ae0-475b-bacc-3799e2949b30
md"""
Analyzes the orbit of a n-body halo in a gravitational potential.
Requires the centres to be calculated prior.
"""

# ╔═╡ ab57edae-2292-4aef-9c1f-53802dbc0600
import TOML

# ╔═╡ 643cd0bf-77b3-4201-9ff7-09dd5aee277c
md"""
# inputs
"""

# ╔═╡ ac2c7484-9acd-4fda-9699-fdf17da507c2
dir = "/arc7/home/dboyea/sculptor/orbit1/"

# ╔═╡ d142b7bd-3002-4331-a725-577873c42f28
properties_file =  joinpath(dir, "properties.toml")

# ╔═╡ 0dd476fd-be53-4e9b-a686-a4462485c64c
orbit_file = joinpath(dir, "orbit.csv")

# ╔═╡ 2bc762ad-e590-443e-b3c2-91dc42a8a4d9
outfile = joinpath(dir, "orbital_properties.toml")

# ╔═╡ b75f0fb1-be59-416c-a61f-4109bada9ae9
r_h = 0.11 # order of mag, for chi sq fit

# ╔═╡ 30969f77-667e-4ae4-9897-82c1c1182652
md"""
# File loading
"""

# ╔═╡ 96a57df5-a7b7-447a-a4a6-2b05e391a5c6
obs_today = TOML.parsefile(properties_file)

# ╔═╡ a609f221-0721-4b4b-a393-49b386393c66
obs_today["ra_err"] = r_h 

# ╔═╡ 68d805a4-c5eb-4f2c-ba10-48c2a53f2874
obs_today["dec_err"] = r_h 

# ╔═╡ b250bf10-c228-4b14-938a-35561ae871d7
begin 
	cens = CSV.read(dir * "out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
begin 
	orbit_expected = CSV.read(orbit_file, DataFrame)
	x_cen_exp = transpose(hcat(orbit_expected.x, orbit_expected.y, orbit_expected.z))
	v_cen_exp = -transpose(hcat(orbit_expected.v_x, orbit_expected.v_y, orbit_expected.v_z))

end

# ╔═╡ 08c3df42-738b-47c4-aa6b-fc39a9cfc02f
md"""
# plots
"""

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
	lines!(cens.t * lguys.T0, r, label="model")
	lines!(lguys.T0*(orbit_expected.t .- orbit_expected.t[begin]), lguys.calc_r(x_cen_exp),
		label="expected"
	)

	axislegend(ax)
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
		ylabel="chi2 property fit",
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

# ╔═╡ a5ce5442-73ca-4aaf-915a-72fe9936e791
d_idx = 20

# ╔═╡ 7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
begin
	peri_filt = idx_f-d_idx:idx_f
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
	
	scatter!(cens.t[idx_peri] * lguys.T0, r[idx_peri], 
		label="last pericentre"
	)

	scatter!(cens.t[idx_apo] * lguys.T0, r[idx_apo], 
		label="last apocentre"
	)
	
	scatter!(cens.t[idx_anteperi] * lguys.T0, r[idx_anteperi], 
		label="last last pericentre"
	)
	
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ af8a50bd-e761-4439-9fc9-80048c264d5b
begin 
	if idx_peri > 0
		t_peri = lguys.T0 * cens.t[idx_peri]
		r_peri = r[idx_peri]

	else 
		t_peri = NaN
		r_peri = NaN
	end

	if idx_apo > 0
		t_apo = lguys.T0 * cens.t[idx_apo]
		r_apo = r[idx_apo]
	else
		t_apo = NaN
		r_apo = NaN
	end

	if idx_anteperi > 0
		t_anteperi = lguys.T0 * cens.t[idx_anteperi]
	else
		t_anteperi = NaN
	end

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


	open(outfile, "w") do f
		TOML.print(f, orbital_properties)
	end

	println("saved properties to $outfile")
	orbital_properties
end

# ╔═╡ Cell order:
# ╠═8b41af50-9ae0-475b-bacc-3799e2949b30
# ╠═061b1886-1878-11ef-3806-b91643300982
# ╠═ab57edae-2292-4aef-9c1f-53802dbc0600
# ╟─643cd0bf-77b3-4201-9ff7-09dd5aee277c
# ╠═ac2c7484-9acd-4fda-9699-fdf17da507c2
# ╠═d142b7bd-3002-4331-a725-577873c42f28
# ╠═0dd476fd-be53-4e9b-a686-a4462485c64c
# ╠═2bc762ad-e590-443e-b3c2-91dc42a8a4d9
# ╠═b75f0fb1-be59-416c-a61f-4109bada9ae9
# ╟─30969f77-667e-4ae4-9897-82c1c1182652
# ╠═96a57df5-a7b7-447a-a4a6-2b05e391a5c6
# ╠═a609f221-0721-4b4b-a393-49b386393c66
# ╠═68d805a4-c5eb-4f2c-ba10-48c2a53f2874
# ╠═b250bf10-c228-4b14-938a-35561ae871d7
# ╠═bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
# ╠═08c3df42-738b-47c4-aa6b-fc39a9cfc02f
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
# ╠═a5ce5442-73ca-4aaf-915a-72fe9936e791
# ╠═7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
# ╠═04d29fcb-70a0-414b-a487-7a18c44b9d58
# ╠═af8a50bd-e761-4439-9fc9-80048c264d5b
# ╠═73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
# ╠═14eebce8-04f7-493b-824a-7808c7fa35dd
# ╠═cdabdc7d-76a1-45f5-b83a-2454576d3964
# ╠═7f2305c1-3431-4e0d-801f-a89afbc3b55b
# ╠═76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
