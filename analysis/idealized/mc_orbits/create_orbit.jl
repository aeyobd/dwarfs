### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8f8f7fea-f463-11ef-0daf-2b30c3565e25
begin
	using Pkg; Pkg.activate()

	using LilGuys

	using PythonCall
end

# ╔═╡ dc8c001c-1af8-4a7b-8cdb-4678f5a24c8f
using CairoMakie

# ╔═╡ 6dfd08be-6833-4901-b8e3-70f854697728
using PlutoUI

# ╔═╡ c8d08554-186e-4e3d-8299-330b168d1f86
using CSV, DataFrames

# ╔═╡ c12a5fe5-59dd-41db-9355-ca3e303c2c9e
using OrderedCollections

# ╔═╡ 8335702e-2cc5-44df-a9c1-8df713ced453
include(ENV["DWARFS_ROOT"] * "/utils/agama_utils.jl")

# ╔═╡ 38fbfd8d-27f6-47d9-9fb5-167ee002b085
md"""
This notebook creates an idealized orbit in the specified plane with a given peri and apocentre.
"""

# ╔═╡ 910f0e20-f077-4d83-aed7-58245469b6e0
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 802de584-2a6b-4bab-8502-70ff50489246
@bind inputs confirm(notebook_inputs(;
	apocentre = NumberField(10:0.1:300, default=100),
	pericentre = NumberField(0.1:0.1:100, default=5),
))

# ╔═╡ 539105c2-c885-4a81-8be2-d25d1a7ba0f3
Phi_py = agama.Potential("../agama_potential.ini")

# ╔═╡ c6f999ec-0e57-4c50-8afd-ec9d060181f0
py2f(x) = pyconvert(Float64, x)

# ╔═╡ f911eb19-ac18-4778-8f51-0d8e22dc9c09
r_hat = [1, 0, 0]

# ╔═╡ 93165b16-53b4-4bd1-bdd7-fd6d4fd129cc
v_hat = [0, 1, 0]

# ╔═╡ a4d3681b-dcb3-48ed-8f29-5301e70e7a07
Φ(r) = py2f(Phi_py.potential(r_hat * r))

# ╔═╡ 0c110304-3fce-426f-8447-1aa12b78ff89
r_peri = inputs.pericentre

# ╔═╡ c07f60eb-e3a8-4207-8d51-158cdfd2c607
r_apo = inputs.apocentre

# ╔═╡ cf476bb9-2511-4276-a2c2-718aaf3f7470
x0 = r_hat * r_apo

# ╔═╡ a626708d-dbba-477e-b72a-c4e879ede7e9
L = sqrt(2 * (Φ(r_apo) - Φ(r_peri)) / (1/r_peri^2 - 1 / r_apo^2))

# ╔═╡ 5801fe76-0112-46d3-aff8-4b5f54932044
v0 = L / r_apo * v_hat

# ╔═╡ d3d9ca2f-6aa3-46ec-b144-491327e3688f
gc = Galactocentric(x0, v0 * V2KMS)

# ╔═╡ d8dfa384-d19a-43b1-bb8c-f30c7de05f34
orbit = calc_orbit(gc, Phi_py, time = 10 / T2GYR)

# ╔═╡ 70b934f1-34a7-42ea-8d3f-41f3492f0a8d
LilGuys.plot_xyz(orbit.positions)

# ╔═╡ e4289648-e661-4209-ace5-92c4b736b811
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "galactocentric radius / kpc"
	)

	lines!(orbit.times * T2GYR, calc_r(orbit.positions))

	fig

end

# ╔═╡ 0b840b8b-91b7-4fe1-9b6a-5926f02a7d2f
minimum(calc_r(orbit.positions))

# ╔═╡ 7d681d34-2bf6-4aae-b2b6-3ef6759c19b8
r_peri; @bind write_orbit CheckBox()

# ╔═╡ d78f8d6b-1a96-48c9-a67a-4b4f958eca44
filename = "orbit_$(r_peri)_$(r_apo)"

# ╔═╡ 1fe08ed1-77bd-4efb-9193-1f0e2202733a
import TOML

# ╔═╡ c9fc5439-8c5c-452e-a041-573fa82e1da2
pwd()

# ╔═╡ 24cda804-99f6-456c-9495-c37d6511e72a
orbit

# ╔═╡ d3b89db4-bdfb-4274-a2e0-b359c650bd20
df_props = OrderedDict(
	:position_i => x0,
	:velocity_i => v0,
)

# ╔═╡ 6c4dc92c-7f0c-4e0d-871e-494fb207cc4d
if write_orbit
	write(filename * ".csv", orbit)

	open(filename * ".toml", "w") do f
		TOML.print(f, df_props)
	end

end

# ╔═╡ Cell order:
# ╟─38fbfd8d-27f6-47d9-9fb5-167ee002b085
# ╠═802de584-2a6b-4bab-8502-70ff50489246
# ╠═8f8f7fea-f463-11ef-0daf-2b30c3565e25
# ╠═dc8c001c-1af8-4a7b-8cdb-4678f5a24c8f
# ╠═6dfd08be-6833-4901-b8e3-70f854697728
# ╠═910f0e20-f077-4d83-aed7-58245469b6e0
# ╠═8335702e-2cc5-44df-a9c1-8df713ced453
# ╠═539105c2-c885-4a81-8be2-d25d1a7ba0f3
# ╠═c6f999ec-0e57-4c50-8afd-ec9d060181f0
# ╠═f911eb19-ac18-4778-8f51-0d8e22dc9c09
# ╠═93165b16-53b4-4bd1-bdd7-fd6d4fd129cc
# ╠═a4d3681b-dcb3-48ed-8f29-5301e70e7a07
# ╠═0c110304-3fce-426f-8447-1aa12b78ff89
# ╠═c07f60eb-e3a8-4207-8d51-158cdfd2c607
# ╠═cf476bb9-2511-4276-a2c2-718aaf3f7470
# ╠═a626708d-dbba-477e-b72a-c4e879ede7e9
# ╠═5801fe76-0112-46d3-aff8-4b5f54932044
# ╠═d3d9ca2f-6aa3-46ec-b144-491327e3688f
# ╠═d8dfa384-d19a-43b1-bb8c-f30c7de05f34
# ╠═70b934f1-34a7-42ea-8d3f-41f3492f0a8d
# ╠═e4289648-e661-4209-ace5-92c4b736b811
# ╠═0b840b8b-91b7-4fe1-9b6a-5926f02a7d2f
# ╠═c8d08554-186e-4e3d-8299-330b168d1f86
# ╠═7d681d34-2bf6-4aae-b2b6-3ef6759c19b8
# ╠═d78f8d6b-1a96-48c9-a67a-4b4f958eca44
# ╠═1fe08ed1-77bd-4efb-9193-1f0e2202733a
# ╠═6c4dc92c-7f0c-4e0d-871e-494fb207cc4d
# ╠═c9fc5439-8c5c-452e-a041-573fa82e1da2
# ╠═24cda804-99f6-456c-9495-c37d6511e72a
# ╠═c12a5fe5-59dd-41db-9355-ca3e303c2c9e
# ╠═d3b89db4-bdfb-4274-a2e0-b359c650bd20
