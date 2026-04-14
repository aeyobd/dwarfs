### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 3a02cce0-3430-11f1-abe3-ef9766dea525
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	
	using PyFITS
end

# ╔═╡ 84d23982-3393-4319-95b8-ded7e66be1c8
using OrderedCollections

# ╔═╡ 84e4e3a4-63c0-4567-869a-8aa889913a6f
CairoMakie.activate!(type=:png)

# ╔═╡ 4f5c53c0-3552-4ab3-8bc9-86042b738622
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "bootes3", "1e5_v30_r2.2")

# ╔═╡ d8e2a511-8209-4a14-9692-0cd90d714d76
function read_scalars(name)
	read_fits(joinpath(model_dir, name, "profiles_scalars.fits"))
end

# ╔═╡ 7af7cf5a-d149-4fe0-9566-d23b293af3ec
scalars = OrderedDict(
	"largeperi+" => read_scalars("orbit_largeperi_long.1"),

	"mean" => read_scalars("orbit_mean.1"),
	"smallperi" => read_scalars("orbit_smallperi.1"),
	# "largeperi" => read_scalars("orbit_largeperi.1"),

	# "one mean" => read_scalars("one_peri.1"),
	# "one smallperi" => read_scalars("one_smallperi.1"),

)

# ╔═╡ 76789b1e-022a-425c-ab8f-0e42f0e422f2
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = L"r_\text{max}", 
			 ylabel = L"\text{v}_\text{max}",
			 xscale=log10,
			 yscale=log10)

	for (label,df) in scalars
		lines!(df.r_circ_max, df.v_circ_max*V2KMS, label=label)
	end

	axislegend(position=:rb)

	fig

end

# ╔═╡ 7c7cf8cf-66e3-4a90-ba31-d41c993b2d81
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = L"\text{v}_\text{max}",
			 yscale=log10)

	for (label,df) in scalars
		lines!(df.time * T2GYR, df.v_circ_max*V2KMS, label=label, alpha=1)
	end

	fig

end

# ╔═╡ 63472478-d8be-422d-8c91-ceb2f8a3ed9b


# ╔═╡ Cell order:
# ╠═3a02cce0-3430-11f1-abe3-ef9766dea525
# ╠═84d23982-3393-4319-95b8-ded7e66be1c8
# ╠═84e4e3a4-63c0-4567-869a-8aa889913a6f
# ╠═4f5c53c0-3552-4ab3-8bc9-86042b738622
# ╠═d8e2a511-8209-4a14-9692-0cd90d714d76
# ╠═7af7cf5a-d149-4fe0-9566-d23b293af3ec
# ╠═76789b1e-022a-425c-ab8f-0e42f0e422f2
# ╠═7c7cf8cf-66e3-4a90-ba31-d41c993b2d81
# ╠═63472478-d8be-422d-8c91-ceb2f8a3ed9b
