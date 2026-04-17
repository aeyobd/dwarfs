### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ f912dcd3-51fe-456f-8189-b960dda0cdc8
import TOML

# ╔═╡ 84e4e3a4-63c0-4567-869a-8aa889913a6f
CairoMakie.activate!(type=:png)

# ╔═╡ 4f5c53c0-3552-4ab3-8bc9-86042b738622
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "bootes3")

# ╔═╡ d8e2a511-8209-4a14-9692-0cd90d714d76
function read_scalars(name)
	read_fits(joinpath(model_dir, name, "profiles_scalars.fits"))
end

# ╔═╡ 507c3285-04c0-4ce8-966a-3ad684d921ac
function read_stars(name, starsname)
	read_fits(joinpath(model_dir, name, "stars", starsname, "stellar_profiles_3d_scalars.fits"))
end

# ╔═╡ 6af2157f-9c6c-4773-8cb8-f8f5ada00e82
stars_1e6 = read_stars("1e6_v35_r3.0/orbit_5Gyr_largeperi/", "exp2d_rs0.20")

# ╔═╡ 7af7cf5a-d149-4fe0-9566-d23b293af3ec
scalars = OrderedDict(
	"18kpc" => read_scalars("1e5_v30_r2.2/orbit_largeperi_long.1"),
	"7kpc" => read_scalars("1e5_v30_r2.2/orbit_mean.1"),
	"1.5kpc" => read_scalars("1e5_v30_r2.2/orbit_smallperi.1"),
	"1e6" => read_scalars("1e6_v35_r3.0/orbit_5Gyr_largeperi/"),
	# "largeperi-light" => read_scalars("1e5_v22_r3.9/orbit_largeperi.1"),
	# "largeperi" => read_scalars("orbit_largeperi.1"),

	# "one mean" => read_scalars("one_peri.1"),
	# "one smallperi" => read_scalars("one_smallperi.1"),

)

# ╔═╡ 21be881c-a87c-4604-a4fb-e7ecb51897ab


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
			 )

	for (label,df) in scalars
		lines!(df.time * T2GYR, df.v_circ_max*V2KMS, label=label, alpha=1)
	end

	fig

end

# ╔═╡ 9eb86b52-33ed-4e01-b5dc-3b8e54a14608
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!(stars_1e6.time, stars_1e6.r_h)

	ylims!(0, 3)
	fig

end

# ╔═╡ fedfa580-bf42-4440-8b82-1e5435e9e85b
let
	fig = Figure()
	ax = Axis(fig[1,1])

	M_h = []

	df = scalars["1e6"]
	
	for row in eachrow(stars_1e6)
		idx = argmin(abs.(row.time .- df.time))
		@info idx
		rmax, vmax = df.r_circ_max[idx], df.v_circ_max[idx]

		halo = ExpCusp(rmax, vmax)
		M = LilGuys.mass(halo, row.r_h)

		push!(M_h, M)
		
	end


	y = stars_1e6.r_h .* M_h
	lines!(stars_1e6.time, y ./ y[1])
	y =M_h
	lines!(stars_1e6.time, y ./ y[1])
	ylims!(0, 2)

	fig

end


# ╔═╡ 63472478-d8be-422d-8c91-ceb2f8a3ed9b
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ b152eb97-6e1a-4886-87e7-88e80663051f
σv = obs_props["sigma_v"]

# ╔═╡ 693cc2b8-02b1-4830-8dea-303d195f1e8c
σv_range = obs_props["sigma_v"] - obs_props["sigma_v_em"], obs_props["sigma_v"] + obs_props["sigma_v_ep"]

# ╔═╡ 6cbce1fb-2abe-440e-b2b5-4defdfc6e9ab
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = L"\sigma_\text{v}",
			 )


	for (label,df) in scalars
		filt = .!ismissing.(df.r_circ_max)
		σv = calc_σv.(df.r_circ_max[filt], df.v_circ_max[filt])
		lines!(df.time[filt] * T2GYR, σv*V2KMS, label=label, alpha=1)
	end

	hlines!(σv, color=:black)
	hspan!(σv_range..., alpha=0.2, color=:black)

	Legend(fig[1, 2], ax)
	fig

end

# ╔═╡ 9e0156c6-55b7-464b-9696-1b0982f16a47
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr", 
			 ylabel = L"\sigma_\text{v}",
			 )

	df = scalars["1e6"]
	filt = .!ismissing.(df.r_circ_max)
	sigmas = calc_σv.(df.r_circ_max[filt], df.v_circ_max[filt], 0.2)
	lines!(df.time[filt] * T2GYR, sigmas*V2KMS, alpha=1, label="approximate")

	lines!(stars_1e6.time * T2GYR, stars_1e6.sigma_v * V2KMS, label="actual")

	hlines!(σv, color=:black)
	hspan!(σv_range..., alpha=0.2, color=:black)

	Legend(fig[1, 2], ax)
	fig

end

# ╔═╡ a27bea66-cdcb-4aee-ad48-a309351c0a6e
R_h = LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])

# ╔═╡ a8be9a3e-6f5f-4371-88ec-ce0306ec9975
α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())

# ╔═╡ b5d50f54-1bfb-49e6-9c3d-473a36630586
prof_stars = LilGuys.Exp2D(R_s=R_h / LilGuys.R_h(LilGuys.Exp2D()))

# ╔═╡ 5103affe-1820-4d4e-82a9-c540011f7175
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 49249f24-9ef6-4a30-911e-de948c04a19f



# ╔═╡ 06ebda0e-7f63-4428-a7ad-81506589ee7c
calc_σv(2, 0.12) * V2KMS

# ╔═╡ Cell order:
# ╠═3a02cce0-3430-11f1-abe3-ef9766dea525
# ╠═f912dcd3-51fe-456f-8189-b960dda0cdc8
# ╠═84d23982-3393-4319-95b8-ded7e66be1c8
# ╠═84e4e3a4-63c0-4567-869a-8aa889913a6f
# ╠═4f5c53c0-3552-4ab3-8bc9-86042b738622
# ╠═d8e2a511-8209-4a14-9692-0cd90d714d76
# ╠═507c3285-04c0-4ce8-966a-3ad684d921ac
# ╠═6af2157f-9c6c-4773-8cb8-f8f5ada00e82
# ╠═7af7cf5a-d149-4fe0-9566-d23b293af3ec
# ╠═21be881c-a87c-4604-a4fb-e7ecb51897ab
# ╠═76789b1e-022a-425c-ab8f-0e42f0e422f2
# ╠═7c7cf8cf-66e3-4a90-ba31-d41c993b2d81
# ╠═6cbce1fb-2abe-440e-b2b5-4defdfc6e9ab
# ╠═9eb86b52-33ed-4e01-b5dc-3b8e54a14608
# ╠═fedfa580-bf42-4440-8b82-1e5435e9e85b
# ╠═9e0156c6-55b7-464b-9696-1b0982f16a47
# ╠═b152eb97-6e1a-4886-87e7-88e80663051f
# ╠═693cc2b8-02b1-4830-8dea-303d195f1e8c
# ╠═63472478-d8be-422d-8c91-ceb2f8a3ed9b
# ╠═a27bea66-cdcb-4aee-ad48-a309351c0a6e
# ╠═a8be9a3e-6f5f-4371-88ec-ce0306ec9975
# ╠═b5d50f54-1bfb-49e6-9c3d-473a36630586
# ╠═5103affe-1820-4d4e-82a9-c540011f7175
# ╠═49249f24-9ef6-4a30-911e-de948c04a19f
# ╠═06ebda0e-7f63-4428-a7ad-81506589ee7c
