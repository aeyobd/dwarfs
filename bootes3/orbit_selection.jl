### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 3a02cce0-3430-11f1-abe3-ef9766dea525
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	
end

# ╔═╡ 84d23982-3393-4319-95b8-ded7e66be1c8
using OrderedCollections

# ╔═╡ 84e4e3a4-63c0-4567-869a-8aa889913a6f
CairoMakie.activate!(type=:png)

# ╔═╡ 4f5c53c0-3552-4ab3-8bc9-86042b738622
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "bootes3", "1e5_v30_r2.2")

# ╔═╡ d8e2a511-8209-4a14-9692-0cd90d714d76
function read_orbit(name)
	Orbit(joinpath(model_dir, name, "simulation/orbit.csv"))
end

# ╔═╡ 7af7cf5a-d149-4fe0-9566-d23b293af3ec
scalars = OrderedDict(
	"largeperi+" => read_orbit("orbit_largeperi_long.1"),
	"mean" => read_orbit("orbit_mean.1"),
	"smallperi" => read_orbit("orbit_smallperi.1"),
	# "largeperi" => read_scalars("orbit_largeperi.1"),

	# "one mean" => read_scalars("one_peri.1"),
	# "one smallperi" => read_scalars("one_smallperi.1"),

)

# ╔═╡ 1e8d2a70-53ec-47eb-8e7b-4090a81aafd4
LilGuys.plot_xyz(LilGuys.positions.(values(scalars))..., )

# ╔═╡ c9706f79-839f-4447-9bc4-1bbcc59c8b8d
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = "galcen radius / kpc")

	for (label, orbit) in scalars
		lines!(T2GYR * LilGuys.times(orbit), radii(orbit), label=label)
	end

	axislegend(position=:lb)
	ylims!(0, nothing)

	fig
end

# ╔═╡ 75941aea-04ba-4c15-b9de-3f93163069af
LilGuys.positions.(values(scalars))

# ╔═╡ 5be62fd5-2451-4bd7-9344-0f56eb9325a8
import Agama

# ╔═╡ dd63e2ec-efe4-4a86-aa57-668d6a3f7c31
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 23d44f9f-9d4b-4c6a-a1c5-bd6ba2a8acaa
L0 = last.(radii.(LilGuys.angular_momenta.(values(scalars))))

# ╔═╡ 784bd0a5-cb82-4905-8215-b1a07c841ab3
r0 = last.(radii.((values(scalars))))

# ╔═╡ 44eaec8d-f315-421f-b291-e68bc781c42f
x0 = [orbit.positions[:, end] for orbit in values(scalars)]

# ╔═╡ 4f4fa650-5c49-4fd9-bd45-4b31dee3d659
v0 = last.((LilGuys.speeds.(values(scalars))))

# ╔═╡ 6ce13075-6c10-43a7-846c-81d5cde85dbe
E0s = @. 1/2 * v0^2 + Agama.potential(pot, x0)

# ╔═╡ d7f6faac-5874-4166-b10c-c5d4c6d09739

let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="L initial", ylabel = "R^3")
		
	scatter!(L0, r0 .^2)

	t = LinRange(0.2, 3, 1000)
	lines!(L0[2] .* t, r0[2] ^ 2 .* ((1 .+ 0.35(t .- 1))))

	fig
end


# ╔═╡ 2c41e8bd-6f33-476a-9459-58dae3ed52c1
r0 = 

# ╔═╡ b634719e-1160-47a2-aa78-4d72149ac0b4
first(values(scalars)).times

# ╔═╡ 8bb5bee1-d78d-449f-a421-0ee24f817f4c
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = "L")

	for (label, orbit) in scalars
		lines!(T2GYR * LilGuys.times(orbit), radii(LilGuys.angular_momenta(orbit)), label=label)
	end

	axislegend(position=:lb)
	ylims!(0, nothing)

	fig
end

# ╔═╡ b14727bc-ca1c-449d-8093-6c3ea875a32e
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = "L")

	for (label, orbit) in scalars
		r = LilGuys
		lines!(T2GYR * LilGuys.times(orbit), radii(LilGuys.angular_momenta(orbit)), label=label)
	end

	axislegend(position=:lb)
	ylims!(0, nothing)

	fig
end

# ╔═╡ 29f75e8f-ab85-4655-bcc8-c418f2a092ec


# ╔═╡ Cell order:
# ╠═3a02cce0-3430-11f1-abe3-ef9766dea525
# ╠═84d23982-3393-4319-95b8-ded7e66be1c8
# ╠═84e4e3a4-63c0-4567-869a-8aa889913a6f
# ╠═4f5c53c0-3552-4ab3-8bc9-86042b738622
# ╠═d8e2a511-8209-4a14-9692-0cd90d714d76
# ╠═7af7cf5a-d149-4fe0-9566-d23b293af3ec
# ╠═1e8d2a70-53ec-47eb-8e7b-4090a81aafd4
# ╠═c9706f79-839f-4447-9bc4-1bbcc59c8b8d
# ╠═75941aea-04ba-4c15-b9de-3f93163069af
# ╠═5be62fd5-2451-4bd7-9344-0f56eb9325a8
# ╠═dd63e2ec-efe4-4a86-aa57-668d6a3f7c31
# ╠═23d44f9f-9d4b-4c6a-a1c5-bd6ba2a8acaa
# ╠═784bd0a5-cb82-4905-8215-b1a07c841ab3
# ╠═44eaec8d-f315-421f-b291-e68bc781c42f
# ╠═4f4fa650-5c49-4fd9-bd45-4b31dee3d659
# ╠═6ce13075-6c10-43a7-846c-81d5cde85dbe
# ╠═d7f6faac-5874-4166-b10c-c5d4c6d09739
# ╠═2c41e8bd-6f33-476a-9459-58dae3ed52c1
# ╠═b634719e-1160-47a2-aa78-4d72149ac0b4
# ╠═8bb5bee1-d78d-449f-a421-0ee24f817f4c
# ╠═b14727bc-ca1c-449d-8093-6c3ea875a32e
# ╠═29f75e8f-ab85-4655-bcc8-c418f2a092ec
