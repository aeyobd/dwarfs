### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 0e1bdc30-2921-11f1-8312-7f57617def29
begin
	import Pkg;
	Pkg.activate()

	using CSV, DataFrames
	using LilGuys
	using Arya, CairoMakie
end

# ╔═╡ 25d2a280-6946-4359-b997-6ffd8c02d0ca
delve = CSV.read("data/delve_dr3_boo3.csv", DataFrame)

# ╔═╡ 1b0121d0-8016-4323-9e1d-efaa8e27c5e4
import TOML

# ╔═╡ dadc562e-7fe0-4b4c-8db1-be029a96abd9
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ 462f3ff1-a2ca-4031-9f7c-be5d61c96894


# ╔═╡ fba8af77-3db3-48e5-88cd-603478d1aeab
xieta = LilGuys.to_tangent.((delve.ra), vec(delve.dec), obs_props["ra"], obs_props["dec"])

# ╔═╡ 9cd5dff7-20a1-41b5-b89b-5c4e45db5545
xi = first.(xieta)

# ╔═╡ abf04b52-22e4-4c12-bf14-efdb3be8e3d9
eta = last.(xieta)

# ╔═╡ 1d5af445-58b6-43f2-ad38-f6722d429df5
let
	fig = Figure()

	ax = Axis(fig[1,1], xreversed=true, 
			  xlabel = "RA / deg",
			  ylabel = "dec / deg",
			  autolimitaspect = 1/ cosd(obs_props["dec"])
			 )


	scatter!(delve.ra, delve.dec)
	scatter!(obs_props["ra"], obs_props["dec"])

	fig
end

# ╔═╡ 5fe8cd32-7379-46f7-b814-84a4b8612ae6
let
	fig = Figure()

	ax = Axis(fig[1,1], xreversed=true, 
			  xlabel = "RA / deg",
			  ylabel = "dec / deg",
			  autolimitaspect = 1/ cosd(obs_props["dec"])
			 )

	scatter!(xi, eta)

	fig
end

# ╔═╡ Cell order:
# ╠═0e1bdc30-2921-11f1-8312-7f57617def29
# ╠═25d2a280-6946-4359-b997-6ffd8c02d0ca
# ╠═1b0121d0-8016-4323-9e1d-efaa8e27c5e4
# ╠═dadc562e-7fe0-4b4c-8db1-be029a96abd9
# ╠═462f3ff1-a2ca-4031-9f7c-be5d61c96894
# ╠═fba8af77-3db3-48e5-88cd-603478d1aeab
# ╠═9cd5dff7-20a1-41b5-b89b-5c4e45db5545
# ╠═abf04b52-22e4-4c12-bf14-efdb3be8e3d9
# ╠═1d5af445-58b6-43f2-ad38-f6722d429df5
# ╠═5fe8cd32-7379-46f7-b814-84a4b8612ae6
