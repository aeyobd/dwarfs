### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 6a2488b6-61b6-11f0-200b-5bb12181877f
begin
	import Pkg; Pkg.activate()

	import Agama
	using Arya, CairoMakie

	using LilGuys
end

# ╔═╡ 31aec302-7ac6-46da-8c82-bf5e4514753e
galaxyname = "tucana4"

# ╔═╡ 5061aa07-d028-4b3e-b85e-4e4d91888db5
potentialname = "vasiliev24/L2M11/potential.ini"

# ╔═╡ 546fa6a8-766e-4007-9b89-1ae8fd508295
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/$potentialname"))

# ╔═╡ 28757126-e326-4935-a2ec-8a2480b7f4f4
units = :vasiliev

# ╔═╡ 39369256-56b3-4131-a3df-857c9782abee
tmax = -10/T2GYR

# ╔═╡ ce9cc08a-a809-4070-a453-fdd00a31046b
import TOML

# ╔═╡ 6452886d-8706-4bde-b792-ab3238cbc85a
CairoMakie.activate!(type=:png)

# ╔═╡ 14fcc92e-9a19-4a69-9e33-3393e67a5953
units_dict = Dict(
	:code => Agama.AgamaUnits(),
	:vasiliev => Agama.VASILIEV_UNITS,
	:portail => Agama.AgamaUnits(),
)

# ╔═╡ 776c6c7b-1809-463b-944b-79c0dbfad2f8
if units ∈ keys(units_dict)
	agama_units = units_dict[units]
else
	error("units $units not known, availabe options: $(keys(units_dict))")
end

# ╔═╡ 10400b5e-ae6a-4d12-aaa6-ff540a9ec966
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))

# ╔═╡ 8172501e-5d18-48b2-88c6-5d1af17e4b4a
coords_i = LilGuys.rand_coords(obs_props, 100)

# ╔═╡ 952e39d8-2fc4-42ce-aaad-a18547d75b3a
o = LilGuys.agama_orbit(pot, coords_i, timerange=(0, tmax), N=1001, agama_units=agama_units)

# ╔═╡ f0c70301-d0cb-4ecf-983b-d248bd8bb389
let
	fig = LilGuys.plot_xyz(LilGuys.positions.(o)..., color=COLORS[1], alpha=0.1)

	LilGuys.plot_xyz!(fig.content, (oo.positions[:, 1] for oo in o)..., color=COLORS[2], plot=:scatter, markersize=3)

	
	fig
end

# ╔═╡ de1c1da1-5a90-4f09-a817-2435c5c7f723
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "radius / kpc",
		limits = (nothing, nothing, 0, nothing)
	)


	for oo in o
		lines!(oo.times * T2GYR, radii(oo), alpha=0.1, color=COLORS[1])
	end


	fig

end

# ╔═╡ 1853f1bb-86a9-41fc-9dee-8986704ee94b
hist(LilGuys.minimum.(LilGuys.radii.(getfield.(o, :positions))))

# ╔═╡ Cell order:
# ╠═31aec302-7ac6-46da-8c82-bf5e4514753e
# ╠═5061aa07-d028-4b3e-b85e-4e4d91888db5
# ╟─546fa6a8-766e-4007-9b89-1ae8fd508295
# ╠═28757126-e326-4935-a2ec-8a2480b7f4f4
# ╟─776c6c7b-1809-463b-944b-79c0dbfad2f8
# ╠═39369256-56b3-4131-a3df-857c9782abee
# ╠═6a2488b6-61b6-11f0-200b-5bb12181877f
# ╠═ce9cc08a-a809-4070-a453-fdd00a31046b
# ╠═6452886d-8706-4bde-b792-ab3238cbc85a
# ╠═14fcc92e-9a19-4a69-9e33-3393e67a5953
# ╠═10400b5e-ae6a-4d12-aaa6-ff540a9ec966
# ╠═8172501e-5d18-48b2-88c6-5d1af17e4b4a
# ╠═952e39d8-2fc4-42ce-aaad-a18547d75b3a
# ╠═f0c70301-d0cb-4ecf-983b-d248bd8bb389
# ╠═de1c1da1-5a90-4f09-a817-2435c5c7f723
# ╠═1853f1bb-86a9-41fc-9dee-8986704ee94b
