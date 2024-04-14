### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 21c8ee88-f68f-11ee-31eb-cb16b9d55d6c
begin 
	import Pkg; Pkg.activate()
	import LilGuys as lguys
	import Roots: find_zero
	import QuadGK: quadgk
	using Plots
	import Arya
	using Revise
end

# ╔═╡ e1be7316-59ac-48aa-b6e9-8745e2acc87e
include("/home/dboyea/project/dwarfs/scripts/nbody_imf.jl")

# ╔═╡ ca1421a7-a41e-4180-aba6-53b9a9d912d9
begin 
	ρ_s(x) =  1/π * exp(-2x) # 3d exponential profile in terms of r/r_h_2d
	σ_s(R) = quadgk(z->ρ_s(sqrt(R^2 + z^2)), -Inf, Inf)[1]
end

# ╔═╡ d8191a4a-75b3-4ff1-ad60-ab5e4843cbb2
begin
	# x is in parsecs
	# y is magnitude
x_bg = [0.15522171200877458, 0.4713977312156778, 0.7759112456012794, 1.0884431679356195, 1.3636089579100605, 1.6268708306315338, 2.0632259219899836, 2.331067439691089, 2.556319699692587, 2.9910857045868773, 3.3767720991259496, 3.6472433216507416, 3.971349188806047, 5.157013158385001, 5.344459082411294, 6.044640098190248, 42.12695838759916, 37.171613669721246, 39.34181812216651, 33.992717321125255, 27.56760362207403, 24.069326575980654, 20.20720794683603, 18.1053899624672, 16.394411932589335, 13.174677321205444, 11.11383314963997, 9.64809482567988, 8.304863605934587, 7.389243960384498, 6.8251767261546155]
y_bg = [21.276920833940807, 21.330339699664673, 21.204995869174322, 21.500704670262916, 21.168294697963745, 21.232055207270253, 21.02156776984011, 21.323030568110028, 21.323030568110028, 21.650852894007873, 21.58079409048938, 21.485386596685622, 21.45607231374836, 21.653496622442532, 21.822772998979445, 21.905894931233902, 28.090431063809106, 27.801953637556494, 27.59729795402634, 26.953861107061282, 26.073810565194147, 25.441570685717064, 24.66291490499101, 24.451183360062206, 24.097157019973757, 23.386616124799534, 23.000554016620498, 22.67405355493998, 22.45897847110852, 22.200981678573164, 22.07470476745881]
end

# ╔═╡ 955a027b-b223-4129-8cec-6a6af618dcaf
quadgk(r-> 4π * r.^2 * ρ_s(r), 0, 100)

# ╔═╡ b1994153-8b07-4b1f-bb7b-140a3c791781
begin 
	Rs = LinRange(0, 10, 10000)
	dR = lguys.gradient(Rs)
	σs = σ_s.(Rs)
	Ms_in = cumsum(@. 2π * Rs * σs * dR)
	Ms_in_3 = cumsum(@. 4π * Rs^2 * ρ_s.(Rs) * dR)

end

# ╔═╡ f9d3d42b-a42e-4686-ae02-3838e3558650
begin 
	plot(xlabel="R", ylabel="M inside (2D)")
	plot!(Rs, Ms_in)
	scatter!([1], [0.5], label="r_H")
	plot!(Rs, Ms_in_3, label="3d")

end

# ╔═╡ 710f94c8-a1c0-4aec-a5ad-8fd3872a1243
Rs[argmin(abs.(Ms_in .- 0.5))]

# ╔═╡ 01180cd1-5288-44b5-862c-67e8a36c4111
σ_exp(r) = 1/2.3361*exp(-1.64r)

# ╔═╡ 671bf9d8-cef3-4108-a7dd-75430b2d44a5
quadgk(r->2π*r*σ_exp(r), 0, 100)

# ╔═╡ ba26b7ef-0dd1-4dbf-b6d3-867c1d01a936
begin 
	plot(xlabel="log R", ylabel=raw"log \sigma")
	plot!(log10.(Rs), log10.(σs), label="3d exponential")
	plot!(log10.(Rs), log10.(σ_exp.(Rs)), label="2d exponential")
	vline!([0], label="r_H")
end

# ╔═╡ 5cd47f0f-5c1e-4f88-bfa7-51d22f2e5f4e
md"""
The snapshot r_h = 2.3 arcmins at a distance of 16 kpc
"""

# ╔═╡ f18631fa-ec02-42a1-a5a1-a68997c67b42
r_h = 2.3 * 60 / 206265 * 16

# ╔═╡ 14d3abe1-740c-4fc9-ba9e-750cde00f41d
5e4/lguys.M0

# ╔═╡ eb96f80e-fd4c-4f62-a77d-0c39a27a3861
snap_exp = lguys.Snapshot("isolation/exponential.hdf5")

# ╔═╡ a1b9da13-6cce-4ad8-96aa-b409480b60cf
r_bins, ρ_n = lguys.calc_ρ_hist(lguys.calc_r(snap_exp), 60, weights=snap_exp.masses)

# ╔═╡ 98729a7c-a413-423f-bb44-ec3425ad71db
begin 
	plot(log10.(lguys.midpoint(r_bins)), log10.(ρ_n), label="nbody model")
	plot!(log10.(Rs), log10.(σs), label="3d exponential")
end

# ╔═╡ 67dcfd2c-efcf-4d24-aaaa-40425801bebc
snap_i = lguys.Snapshot("isolation/varmass/initial.hdf5")

# ╔═╡ c0737916-605e-495b-aac1-f008438cc937
r_bins_i, ρ_n_i = lguys.calc_ρ_hist(lguys.calc_r(snap_i), 60, weights=snap_i.masses)

# ╔═╡ 3eec7494-95cd-4158-b52a-42d60ed8c5a5
begin 
	plot(xlabel="log r / kpc", ylabel="log density")
	plot!(log10.(lguys.midpoint(r_bins_i)), log10.(ρ_n_i), label="nbody model")
	plot!(log10.(Rs * r_h), log10.(σs * 4e-6/r_h^3), label="3d exponential in 2d")
	plot!(log10.(Rs * r_h), log10.(ρ_s.(Rs) * 4e-6/r_h^3), label="3d exponential")

	scatter!(log10.(x_bg / 1e3), -1/2.5*y_bg .+ 8.6, label="Baumgardt observations")
end

# ╔═╡ 79db0d71-f5ff-4098-89c1-c83cedf7af2c
sum(1/2 * lguys.calc_radial_discrete_Φ(snap_i))

# ╔═╡ 4d78f141-7a95-4877-9b21-d0547e6c7942
sum(1/2 * lguys.calc_v(snap_i) .^ 2) * 2

# ╔═╡ 58ae2346-b169-4651-8f13-3908728ee7e0
minimum(lguys.calc_r(snap_i))

# ╔═╡ 4bf3da00-7c91-4107-92f9-942ab2372b78
begin 
	plot(xlabel="x / pc", ylabel="y / pc", aspect_ratio=1, legend=false, xlims=[-70, 70], ylims=[-70, 70])
	scatter!(snap_i.positions[1, :]*1e3, snap_i.positions[2, :]*1e3, ms=1, alpha=0.3)
end

# ╔═╡ d54945f0-a06c-41eb-99c9-8998c759d5c6
md"""
The IMF Sampling method
"""

# ╔═╡ b1ed4a7c-e841-4d37-b832-20958ab0e04a
ms = LinRange(0.08, 100, 1000)

# ╔═╡ c16979e2-4915-418e-ba17-fd80c5a4ecf5
mm = sample_ms(m_tot=5.0e-6)

# ╔═╡ bb4d8d6e-238f-4a42-875f-36e9ffb07344
ps = LinRange(0, 1, 100)

# ╔═╡ 57651083-e99f-4aa7-9420-e3e4e1ddcf16
sum(mm)

# ╔═╡ c5df6ebc-1084-466c-a8a0-7a818e0cfc97
length(mm)

# ╔═╡ ce619bc1-1687-446a-8bd9-147cb711006a
iimf = calc_inv_imf()

# ╔═╡ f4b9bc26-a154-489b-88b6-ebee08973e02
plot(ps, iimf.(ps), yscale=:log10)

# ╔═╡ 4ff91067-ed4c-41e6-9193-27b58c548cbc
begin 
	plot()
	bins, h = lguys.calc_histogram(mm, 1000)
	h ./= length(mm)
	h ./= bins[2] - bins[1]
	scatter!(log10.(lguys.midpoint(bins)), log10.(h))
	plot!(log10.(ms), log10.(kroupa_imf.(ms)))
end

# ╔═╡ 70541f54-d8d3-425a-b3e9-f71f5d1103b4
length(h)

# ╔═╡ 826336e3-359e-45be-9996-6317cba988df
ms

# ╔═╡ 29c87835-1b93-47da-b315-f7a2560656df
quadgk(m->kroupa_imf(m) * m, 0.08, 100)

# ╔═╡ 0779e921-6d25-4306-a6cc-538a49e6948f
quadgk(m->kroupa_imf(m), 0.08, 100)

# ╔═╡ 63525479-85fe-4f7e-bd8a-e86699e27508
lguys.randu(0, 1, 100)

# ╔═╡ 8cc4b88a-392c-4131-b586-b345af0227f0
quadgk(m->kroupa_imf(m), 0.08, 100)

# ╔═╡ 258c8953-2250-4df5-a47c-c628ff546d4e
lguys.M0

# ╔═╡ Cell order:
# ╠═21c8ee88-f68f-11ee-31eb-cb16b9d55d6c
# ╠═ca1421a7-a41e-4180-aba6-53b9a9d912d9
# ╠═d8191a4a-75b3-4ff1-ad60-ab5e4843cbb2
# ╠═955a027b-b223-4129-8cec-6a6af618dcaf
# ╠═b1994153-8b07-4b1f-bb7b-140a3c791781
# ╠═f9d3d42b-a42e-4686-ae02-3838e3558650
# ╠═710f94c8-a1c0-4aec-a5ad-8fd3872a1243
# ╠═01180cd1-5288-44b5-862c-67e8a36c4111
# ╠═671bf9d8-cef3-4108-a7dd-75430b2d44a5
# ╠═ba26b7ef-0dd1-4dbf-b6d3-867c1d01a936
# ╠═5cd47f0f-5c1e-4f88-bfa7-51d22f2e5f4e
# ╠═f18631fa-ec02-42a1-a5a1-a68997c67b42
# ╠═14d3abe1-740c-4fc9-ba9e-750cde00f41d
# ╠═eb96f80e-fd4c-4f62-a77d-0c39a27a3861
# ╠═a1b9da13-6cce-4ad8-96aa-b409480b60cf
# ╠═98729a7c-a413-423f-bb44-ec3425ad71db
# ╠═67dcfd2c-efcf-4d24-aaaa-40425801bebc
# ╠═c0737916-605e-495b-aac1-f008438cc937
# ╠═3eec7494-95cd-4158-b52a-42d60ed8c5a5
# ╠═79db0d71-f5ff-4098-89c1-c83cedf7af2c
# ╠═4d78f141-7a95-4877-9b21-d0547e6c7942
# ╠═58ae2346-b169-4651-8f13-3908728ee7e0
# ╠═4bf3da00-7c91-4107-92f9-942ab2372b78
# ╠═d54945f0-a06c-41eb-99c9-8998c759d5c6
# ╠═e1be7316-59ac-48aa-b6e9-8745e2acc87e
# ╠═b1ed4a7c-e841-4d37-b832-20958ab0e04a
# ╠═c16979e2-4915-418e-ba17-fd80c5a4ecf5
# ╠═bb4d8d6e-238f-4a42-875f-36e9ffb07344
# ╠═57651083-e99f-4aa7-9420-e3e4e1ddcf16
# ╠═c5df6ebc-1084-466c-a8a0-7a818e0cfc97
# ╠═ce619bc1-1687-446a-8bd9-147cb711006a
# ╠═f4b9bc26-a154-489b-88b6-ebee08973e02
# ╠═4ff91067-ed4c-41e6-9193-27b58c548cbc
# ╠═70541f54-d8d3-425a-b3e9-f71f5d1103b4
# ╠═826336e3-359e-45be-9996-6317cba988df
# ╠═29c87835-1b93-47da-b315-f7a2560656df
# ╠═0779e921-6d25-4306-a6cc-538a49e6948f
# ╠═63525479-85fe-4f7e-bd8a-e86699e27508
# ╠═8cc4b88a-392c-4131-b586-b345af0227f0
# ╠═258c8953-2250-4df5-a47c-c628ff546d4e
