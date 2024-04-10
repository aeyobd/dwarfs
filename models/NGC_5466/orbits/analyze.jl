### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 949f4b7e-f76f-11ee-31ed-69a587dc7973
begin
	import Pkg; Pkg.activate()

	using Plots
	import Arya

	import LilGuys as lguys

	import CSV
	using DataFrames
end

# ╔═╡ 0c85a109-f75b-45b2-be77-1664bfdd4567
path = "orbit1"

# ╔═╡ 7bb77640-17c4-44b1-832f-9258fb790d8b
out = lguys.Output("$path/out")

# ╔═╡ 390bc099-3434-4c19-924c-bde742f9f46c
begin
	cens = CSV.read("$path/data/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ e82adf25-15c3-44dc-ba1f-46e962ef9b14
lguys.plot_xyz_layout(x_cen)

# ╔═╡ 2c426184-6e1f-4256-bd5b-fc70d20e0481
begin 
	plot(xlabel="time / Gyr", ylabel="galactocentric distance / kpc")
	plot!(cens.t * lguys.T0, lguys.calc_r(x_cen), label="")
end

# ╔═╡ 4b456983-017f-47e1-a76a-0688fa53b6fb
begin 
	plot(xlabel="R / kpc", ylabel="z / kpc", aspect_ratio=1, xlim=[0, 50])
	R_cen = @. sqrt(x_cen[1, :]^2 + x_cen[2, :]^2)
	plot!(R_cen, x_cen[3, :], label="")
end

# ╔═╡ 0deee36d-3f91-40ba-9a46-e07a0e813574
pwd()

# ╔═╡ 4f0a2096-0e67-43f9-99da-6fdb00fe2f19
begin 
	anim = @animate for i in 1:10:length(out)
		snap = out[i]
		lguys.plot_centre(snap.positions, 
			ms=1, msw=0, ma=0.1, width=50, z_filter=:none)
	end
	
	gif(anim, "orbit1.gif", fps = 12)
end

# ╔═╡ f83658c4-dc0d-41e8-9707-26bcd1344d1b
begin 
	plot(xlabel="x / kpc", ylabel="y / kpc")
	scatter!(out[end].positions[1, :], out[end].positions[2, :], ms=3, aspect_ratio=1, alpha=0.1, label="")
end

# ╔═╡ fa20f05a-951f-44e1-ae8b-6e24393304f6
snap_f = out[end]

# ╔═╡ 10456dbb-1aa4-436b-a80f-39ec29fa8137
particles_gc = [lguys.Galactocentric(snap_f.positions[:, i]..., (lguys.V0 * snap_f.velocities[:, i])...) for i in 1:length(snap_f)]

# ╔═╡ 1b027221-000b-4929-b893-7f2aa829b3ed
obs_pred = lguys.transform.(lguys.Observation, particles_gc)

# ╔═╡ 44353e99-0968-4619-9c7d-bea32be287fe
begin 
	plot(xlabel="ra", ylabel="dec")
	histogram2d!([o.ra for o in obs_pred], [o.dec for o in obs_pred], colorbar_scale=:log10)
end

# ╔═╡ 620e20ca-1d7a-42be-bd4d-b62f5573f48b
obs_m = lguys.transform(lguys.Observation, lguys.Galactocentric(x_cen[:, end]..., (v_cen[:, end] * lguys.V0)...))

# ╔═╡ a80ec84b-6bf9-4130-9e87-b00ee6d340cf
begin 
	dra = 3
	plot(xlabel="ra", ylabel="dec", xlims=(obs_m.ra - dra, obs_m.ra + dra), 
	ylims=(obs_m.dec - dra, obs_m.dec + dra))
	scatter!([o.ra for o in obs_pred], [o.dec for o in obs_pred], alpha=0.5, label="")
end

# ╔═╡ Cell order:
# ╠═949f4b7e-f76f-11ee-31ed-69a587dc7973
# ╠═0c85a109-f75b-45b2-be77-1664bfdd4567
# ╠═7bb77640-17c4-44b1-832f-9258fb790d8b
# ╠═390bc099-3434-4c19-924c-bde742f9f46c
# ╠═e82adf25-15c3-44dc-ba1f-46e962ef9b14
# ╠═2c426184-6e1f-4256-bd5b-fc70d20e0481
# ╠═4b456983-017f-47e1-a76a-0688fa53b6fb
# ╠═0deee36d-3f91-40ba-9a46-e07a0e813574
# ╠═4f0a2096-0e67-43f9-99da-6fdb00fe2f19
# ╠═f83658c4-dc0d-41e8-9707-26bcd1344d1b
# ╠═fa20f05a-951f-44e1-ae8b-6e24393304f6
# ╠═10456dbb-1aa4-436b-a80f-39ec29fa8137
# ╠═1b027221-000b-4929-b893-7f2aa829b3ed
# ╠═44353e99-0968-4619-9c7d-bea32be287fe
# ╠═620e20ca-1d7a-42be-bd4d-b62f5573f48b
# ╠═a80ec84b-6bf9-4130-9e87-b00ee6d340cf
