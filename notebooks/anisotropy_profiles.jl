### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ ef4b8bdc-957a-11f0-1b7c-eb9a330482dd
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

end

# ╔═╡ 3c5dd571-2889-414b-b4fb-395328851b0b
using LilGuys

# ╔═╡ 18a8ef55-b296-4cc7-b62e-b96457f4dad3
md"""
# Anisotropy profiles.
"""

# ╔═╡ 22d43864-dd05-43fb-8ecc-96b1b2948688
md"""
 taken as maxima from Navarro+ 2010
"""

# ╔═╡ c75e74e0-c69d-4c75-a0d1-f30f89c3fb9a
log_r_n10 = [-1.2931616403333048, -0.793368201978315, -0.3144102080426532, 0.05202235542633549, 0.35473260010444996, 0.501375778503558, 0.6574662291197355, 0.9893367422500408, 1.123305609902487, 1.285499372520286]


# ╔═╡ 845c5238-8ebc-412e-9be5-3ea1b8f10058
beta_n10 = [0.053596055858816924, 0.10058551007093838, 0.15739815510649016, 0.23895252330891645, 0.3234963787564724, 0.40908077396338294, 0.4751591778084251, 0.4732721672130877, 0.29203243841408527, 0.17021496230107935]

# ╔═╡ 959c69ff-c08d-4cf5-bd9b-fd2707dbab4d
function β_om(r, a, β0)
	return @. (β0 + (r/a)^2) / (1 + (r/a)^2)
end

# ╔═╡ 7ea5e47e-af87-4fbb-8888-a9ac4e7f5aad
snap = LilGuys.Snapshot(joinpath(ENV["DWARFS_ROOT"], "agama/halos/nfw_1e6_beta0.2_ra4.0.hdf5"))

# ╔═╡ 0b05671e-0ee4-4d6a-944b-b918cfd51222
bins = LinRange(-2, 1, 100)

# ╔═╡ a33c8581-59a7-46c5-a53c-a2332dd008cb
σ, β = LilGuys.β_prof(snap, r_bins=10 .^bins)

# ╔═╡ 1c79460f-a414-43e1-83e6-5eb6e52b7d6f
halo = LilGuys.NFW(M_s=1, r_s=1)

# ╔═╡ d777349d-d7cd-45e8-bfa1-e95e86eff14e
hist(log10.(radii(snap)))

# ╔═╡ f4ce80f9-f9d5-47d8-9eec-91ff12fd0717
prof = LilGuys.DensityProfile(snap)

# ╔═╡ c35a2a39-7c4f-4844-8ede-31fda3282938
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!(prof.log_r, log10.(prof.rho))
	lines!(bins, log10.(LilGuys.density.(halo, 10 .^ bins)))

	ylims!(-12, 3)
	fig
end
	

# ╔═╡ f9e9cb48-589b-4091-a963-4ef839ca29ac
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	lines!(midpoints(bins), σ)

	r = 10 .^ bins
	ρ = LilGuys.density.(halo, r)
	y = @. (r^1.875 * ρ) ^ (2/3)
	lines!(bins, y .* maximum(filter(isfinite, σ)) / maximum(y))
	@info r[argmax(y)]

	fig
end

# ╔═╡ 55245916-1740-414d-92fe-c0513773eb41
r_σmax = 0.7564633275546291 * halo.r_s

# ╔═╡ 2f003396-9395-4cdb-95a7-4b7e928a160a
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "log radius / rmax",
			 ylabel = L"\beta")


	scatter!(log_r_n10, beta_n10)

	x = LinRange(-2, 1, 1000)
	y = β_om(10 .^ x, 5, 0.2)

	lines!(x, y)
	ylims!(-0.5, 1)
	lines!(midpoints(bins) .- log10(r_σmax), β)
	fig

end

# ╔═╡ 192c6c2c-5007-48d1-a687-bef17dc76581
5 * r_σmax

# ╔═╡ 72d2a8bd-9256-4a9d-a31d-6e680f856a94


# ╔═╡ Cell order:
# ╟─18a8ef55-b296-4cc7-b62e-b96457f4dad3
# ╠═ef4b8bdc-957a-11f0-1b7c-eb9a330482dd
# ╠═3c5dd571-2889-414b-b4fb-395328851b0b
# ╠═22d43864-dd05-43fb-8ecc-96b1b2948688
# ╠═c75e74e0-c69d-4c75-a0d1-f30f89c3fb9a
# ╠═845c5238-8ebc-412e-9be5-3ea1b8f10058
# ╠═959c69ff-c08d-4cf5-bd9b-fd2707dbab4d
# ╠═2f003396-9395-4cdb-95a7-4b7e928a160a
# ╠═7ea5e47e-af87-4fbb-8888-a9ac4e7f5aad
# ╠═0b05671e-0ee4-4d6a-944b-b918cfd51222
# ╠═a33c8581-59a7-46c5-a53c-a2332dd008cb
# ╠═1c79460f-a414-43e1-83e6-5eb6e52b7d6f
# ╠═d777349d-d7cd-45e8-bfa1-e95e86eff14e
# ╠═f4ce80f9-f9d5-47d8-9eec-91ff12fd0717
# ╠═c35a2a39-7c4f-4844-8ede-31fda3282938
# ╠═f9e9cb48-589b-4091-a963-4ef839ca29ac
# ╠═55245916-1740-414d-92fe-c0513773eb41
# ╠═192c6c2c-5007-48d1-a687-bef17dc76581
# ╠═72d2a8bd-9256-4a9d-a31d-6e680f856a94
