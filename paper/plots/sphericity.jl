### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 217b3ca2-65a9-11f1-bd0e-b504553e8444
begin
	import Pkg; Pkg.activate()

	using LilGuys, Arya
	using CairoMakie
	import TOML
end

# ╔═╡ 899c9b60-aa63-44e1-b0d6-0f846fe50858
model_keys = TOML.parsefile("model_key.toml")

# ╔═╡ 0efd64c3-bf38-4ee2-9d29-54cf3ec4a5d1
galaxyname = "ursa_minor"

# ╔═╡ adce841e-9d3f-4d24-abde-9f544ff276f8
if galaxyname == "sculptor"
	modelname, starsname = model_keys["scl_smallperi"][2:3]
elseif galaxyname == "ursa_minor"
	modelname, starsname = model_keys["umi_smallperi"][2:3]
end

# ╔═╡ aa58ea65-9af3-445e-b627-6a3e04c6c362
md"""
# Data Loading
"""

# ╔═╡ cf2e9266-40f3-4cc1-b989-4c52dd84183f
stars = LilGuys.read_hdf5_table(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "../stars", starsname, "probabilities_stars.hdf5"))

# ╔═╡ 93ff0843-3773-4f54-9c05-c970bdb3a35a
out = Output(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname), weights=stars.probability);

# ╔═╡ b2d90997-054d-4caa-9bb1-623f513b31c5
orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "orbital_properties.toml"))

# ╔═╡ 991c5ec4-8c3b-449e-9f37-a25ffd875730
idx_f = orbit_props["idx_f"]

# ╔═╡ 12386b3f-e659-400b-bc95-088a53cf9e8f
snap_f = out[idx_f]

# ╔═╡ 49d87ef9-269a-4df1-8220-0ad1a924bacf
snap_i = out[1]

# ╔═╡ bde2ac5f-1982-44b5-8ebd-054d1b414aaa
md"""
# Calculations
"""

# ╔═╡ 07751d6c-b615-4b82-b268-c13b3b6f8cf8
import LinearAlgebra: eigvals, eigen, norm, normalize, ⋅, ×

# ╔═╡ e6aa77ec-06a5-4489-bad7-0fb861a14064
function moments(snap)
	m = snap.weights
	x = snap.positions .- snap.x_cen
	return moments(x, m)
end

# ╔═╡ 5400dea0-9cd2-4475-847a-e16ad72e8352
function moments(x, m)
	I = (m' .* x ./ (radii(x) .^ 2)') * x'

	return I ./ sum(m)

end

# ╔═╡ 9684c24f-fc47-4e11-9f85-4b678f23bdff
function moments(x, m, r)
	I = (m' .* x ./ (r .^ 2)') * x'

	return I ./ sum(m)

end

# ╔═╡ 2d048d8a-6cd6-437d-9683-cf33ac643193


# ╔═╡ be8a7bd6-4b44-4d4b-9498-6f4c8f898c29
function iterate_eigen(snap)
	I = moments(snap)
	λs_new, A = eigen(I)

	λs_old = zeros(3)
	
	for i in 1:100
		# transform to rotated coordinates
		A_inv = inv(A)
		xp =  A_inv * (snap.positions .- snap.x_cen)
		q, p, _ = sqrt.(λs_new ./ λs_new[3])
		r = radii(xp ./ [q, p, 1])
		
		I = moments((snap.positions .- snap.x_cen), snap.weights, r)

		# printing diags
		println("iter $i")
		println("basis = ", round.(A, digits=2))
		println("p, q = ", p, "\t", q)
		println("λ = ", round.(λs_new, digits=3))

		
		λs_old = copy(λs_new)
		λs_new, A = eigen(I)

		if radii(log10.(λs_old ./ λs_new)) < 0.0001
			break
		end
	end

	return λs_new
end

# ╔═╡ 4e1c61c4-a3ec-4ea7-b6d6-18116f9a57f5
function ellipticity(args...)
	λs = iterate_eigen(args...)
	return 1 - sqrt(minimum(λs) / maximum(λs))
end

# ╔═╡ ce738872-b99f-4536-9b3e-9e641295a4ca
ellipticity(snap_i)

# ╔═╡ 48a4b97f-72a0-4034-9593-103d83e6097d
ellipticity(snap_f[LilGuys.bound_particles(snap_f)])

# ╔═╡ 83a32c70-fd87-4ca6-9255-25e2b8f936bd
ellipticity(snap_f)

# ╔═╡ 77c36219-d2aa-431e-80be-b9381bdefd1c
md"""
# Test
"""

# ╔═╡ c92f8c4c-9a66-4ed8-895d-b1fb40c5c792
q = 0.95; p = 0.98

# ╔═╡ 2b158c7e-c1aa-4b31-8927-6a5d18c8412f
x_hat = vec(LilGuys.rand_unit())

# ╔═╡ d1ff9e6b-1a47-4b17-af0b-898d4327c73c
y_hat = let
	y = vec(LilGuys.rand_unit())
	y .-= y ⋅ x_hat * x_hat
	normalize(y)
end

# ╔═╡ a18480e6-4c6b-417b-bf99-8cb31ffe39d7
z_hat = x_hat × y_hat

# ╔═╡ ad8771ed-7f96-44b3-a397-d8a3951d6e16
norm(z_hat), norm(x_hat), norm(y_hat)

# ╔═╡ 05bf6a56-af54-4ea5-8409-f128174fd901
x_hat ⋅ y_hat, x_hat ⋅ z_hat, y_hat ⋅ z_hat

# ╔═╡ 4d399693-7d44-4445-ac32-ab88e36bfb82
snap_test = let
	snap = LilGuys.sample_potential(LilGuys.TruncNFW(r_s=1, M_s=1, trunc=20, xi=3), 100_000)

	a = snap.positions[1, :]
	b = snap.positions[2, :] * p
	c = snap.positions[3, :] * q

	pos_new = a' .* x_hat .+ b' .* y_hat .+ c' .* z_hat
	
	snap.positions = pos_new
	snap
end
	

# ╔═╡ 174c7549-cbed-437a-a152-f2fa52292f37
scatter(snap_test.positions[1, :], snap_test.positions[2, :], alpha=0.1, markersize=1,  axis=(;
																   limits=(-2, 2, -2, 2), aspect=DataAspect()))

# ╔═╡ 00259a94-5619-409a-94d8-5f28fe728be3
scatter(snap_test.positions[1, :], snap_test.positions[3, :], alpha=0.1, markersize=1,  axis=(;
																   limits=(-2, 2, -2, 2), aspect=DataAspect()))

# ╔═╡ aae17938-10e8-4464-ac21-1670ffc83f20
I_test = moments(snap_test)

# ╔═╡ 14b5a3b9-3427-49bd-904e-88246acb698c
λs, A = eigen(I_test)

# ╔═╡ b5b2220f-23e6-4b6e-baa1-9e2776af2481
(I_test * A[:, 3]) ./ A[:, 3] ./ λs[3]

# ╔═╡ c6281a48-847e-4338-8289-7573f229d3d8
Ainv = A'

# ╔═╡ b4881515-323f-4ddd-9b1b-48b55c5bbd62
Ainv * I_test

# ╔═╡ c51bd889-00da-4546-8336-78dbea138e0c
A * A'

# ╔═╡ 306f4202-aa8b-4140-9d69-b60d51179460
xp = Ainv * (snap_test.positions )

# ╔═╡ f7642a73-93dd-4b6b-a6b5-ea82984523d1
ap = radii(xp ./ [sqrt(λs[3] / λs[1]), sqrt(λs[2] / λs[1]), 1])

# ╔═╡ c1c9ef98-1fad-4919-b179-41b8cc1cfffe
moments(xp ./ sqrt.(λs ./ λs[3]), 
		snap_test.weights, ap)

# ╔═╡ 45384b0c-68db-4fc4-98ec-6b0b14d567bc
Ip = moments(xp, snap_test.weights, ap)

# ╔═╡ b3faa0fe-41d2-4de2-81e1-0b9c345fae3b
λs

# ╔═╡ 43323750-1b59-45ed-8e63-62fffbd5d61e
1 - q

# ╔═╡ 93167c39-5c59-4e56-8c97-6ab36d22f3ac
ellipticity(snap_test)

# ╔═╡ Cell order:
# ╠═217b3ca2-65a9-11f1-bd0e-b504553e8444
# ╠═899c9b60-aa63-44e1-b0d6-0f846fe50858
# ╠═0efd64c3-bf38-4ee2-9d29-54cf3ec4a5d1
# ╠═adce841e-9d3f-4d24-abde-9f544ff276f8
# ╟─aa58ea65-9af3-445e-b627-6a3e04c6c362
# ╠═cf2e9266-40f3-4cc1-b989-4c52dd84183f
# ╠═93ff0843-3773-4f54-9c05-c970bdb3a35a
# ╠═b2d90997-054d-4caa-9bb1-623f513b31c5
# ╠═991c5ec4-8c3b-449e-9f37-a25ffd875730
# ╠═12386b3f-e659-400b-bc95-088a53cf9e8f
# ╠═49d87ef9-269a-4df1-8220-0ad1a924bacf
# ╠═bde2ac5f-1982-44b5-8ebd-054d1b414aaa
# ╠═07751d6c-b615-4b82-b268-c13b3b6f8cf8
# ╠═e6aa77ec-06a5-4489-bad7-0fb861a14064
# ╠═5400dea0-9cd2-4475-847a-e16ad72e8352
# ╠═9684c24f-fc47-4e11-9f85-4b678f23bdff
# ╠═2d048d8a-6cd6-437d-9683-cf33ac643193
# ╠═4e1c61c4-a3ec-4ea7-b6d6-18116f9a57f5
# ╠═be8a7bd6-4b44-4d4b-9498-6f4c8f898c29
# ╠═ce738872-b99f-4536-9b3e-9e641295a4ca
# ╠═48a4b97f-72a0-4034-9593-103d83e6097d
# ╠═83a32c70-fd87-4ca6-9255-25e2b8f936bd
# ╟─77c36219-d2aa-431e-80be-b9381bdefd1c
# ╠═c92f8c4c-9a66-4ed8-895d-b1fb40c5c792
# ╠═2b158c7e-c1aa-4b31-8927-6a5d18c8412f
# ╠═d1ff9e6b-1a47-4b17-af0b-898d4327c73c
# ╠═a18480e6-4c6b-417b-bf99-8cb31ffe39d7
# ╠═ad8771ed-7f96-44b3-a397-d8a3951d6e16
# ╠═05bf6a56-af54-4ea5-8409-f128174fd901
# ╠═4d399693-7d44-4445-ac32-ab88e36bfb82
# ╠═174c7549-cbed-437a-a152-f2fa52292f37
# ╠═00259a94-5619-409a-94d8-5f28fe728be3
# ╠═aae17938-10e8-4464-ac21-1670ffc83f20
# ╠═14b5a3b9-3427-49bd-904e-88246acb698c
# ╠═b5b2220f-23e6-4b6e-baa1-9e2776af2481
# ╠═c6281a48-847e-4338-8289-7573f229d3d8
# ╠═b4881515-323f-4ddd-9b1b-48b55c5bbd62
# ╠═c51bd889-00da-4546-8336-78dbea138e0c
# ╠═306f4202-aa8b-4140-9d69-b60d51179460
# ╠═c1c9ef98-1fad-4919-b179-41b8cc1cfffe
# ╠═f7642a73-93dd-4b6b-a6b5-ea82984523d1
# ╠═45384b0c-68db-4fc4-98ec-6b0b14d567bc
# ╠═b3faa0fe-41d2-4de2-81e1-0b9c345fae3b
# ╠═43323750-1b59-45ed-8e63-62fffbd5d61e
# ╠═93167c39-5c59-4e56-8c97-6ab36d22f3ac
