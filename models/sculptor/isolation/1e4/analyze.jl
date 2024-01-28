### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using Plots; plotly()
	import LilGuys as lguys
end

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
out = lguys.Output("out")

# ╔═╡ c5672da9-0dad-4d22-abe5-9e186ccde02d
begin
	snap_i = lguys.Snapshot("initial.hdf5")
	snap_f = out[2]
	prof_i = lguys.Profile(snap_i)
	prof_f = lguys.Profile(snap_f)
end

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
begin 
	scatter(log10.(prof_i.rs), prof_i.Vs_circ * lguys.V0)
	scatter!(log10.(prof_f.rs), prof_f.Vs_circ * lguys.V0)

	xlabel!("log r")
	ylabel!("V_circ")
end

# ╔═╡ 30f3a998-8c4b-4cc3-97d7-967537f9986f
begin 
	scatter(log10.(prof_i.rs), log10.(prof_i.ρs))
	scatter!(log10.(prof_f.rs), log10.(prof_f.ρs))
	xlabel!("log r / kpc")
	ylabel!("log rho / 1e10 Msun kpc^-3")

end

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
begin 
	histogram2d(log10.(lguys.calc_r(snap_i.positions)), lguys.calc_r(snap_i.velocities) * lguys.V0)
end

# ╔═╡ Cell order:
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═c5672da9-0dad-4d22-abe5-9e186ccde02d
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═30f3a998-8c4b-4cc3-97d7-967537f9986f
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
