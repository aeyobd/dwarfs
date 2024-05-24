### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 6238c3fa-1974-11ef-1539-2f50fb31fe9a
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 141d3b97-e344-4760-aa12-a48b1afb125c
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 5b1dd353-a437-47cd-94be-7da9684581da
modeldir = "/cosma/home/durham/dc-boye1/data/dwarfs/models/sculptor/isolation/1e6"

# ╔═╡ 28ba4f0f-6cc6-44e2-a7bc-4eee460d91b0
starsname = "exp2d_small"

# ╔═╡ 21adbbe7-c8cc-4094-9e75-b68d97fa211a
starsfile = "stars/$(starsname)_stars.hdf5"

# ╔═╡ 312576e2-16da-4285-9c19-a7a8005acf25
paramname = "stars/$(starsname).yml"

# ╔═╡ feb6cc17-a25d-4ba8-a152-78412f422b80
import YAML

# ╔═╡ 4445c2a1-3f00-42fa-8259-5d2b33e61d75
params = YAML.load_file(joinpath(modeldir, paramname)); 


# ╔═╡ 051b314e-6cbc-405e-8873-b0288a93378d
profile_class = getproperty(lguys, Symbol(params["profile"]))

# ╔═╡ 7dfc4066-eee5-4118-9d12-48861aa66e03
NamedTuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

# ╔═╡ c76c8acc-ea88-4ce1-81cb-be0b67ef23fd
profile = profile_class(;NamedTuple(params["profile_kwargs"])...)

# ╔═╡ 0671aa68-3d02-4b6c-a238-adffede13bd8
ρ_s(r) = lguys.calc_ρ(profile, r)

# ╔═╡ fc7e58ec-ac03-40a4-99ae-2a16166402d3


# ╔═╡ e32e2d66-2dad-4f09-84f1-a3081e3891a7
begin 
	out = lguys.Output("$modeldir/out/combined.hdf5")

	cens = CSV.read("$modeldir/out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ 3150cdfd-7573-4db9-86b7-ef614150a7b9
times = out.times * lguys.T0

# ╔═╡ 2cc047db-ae05-42c5-898f-702ae3b83bd6
idx_i = 30; idx_f = length(out)

# ╔═╡ 052d229a-0362-42eb-b03c-fdab0f0bc6b4
begin
	snap_i = out[idx_i]
	snap_i.x_cen = x_cen[:, idx_i]
	snap_i.v_cen = v_cen[:, idx_i]

	snap_f = out[idx_f]
	snap_f.x_cen = x_cen[:, idx_f]
	snap_f.v_cen = v_cen[:, idx_f]

end

# ╔═╡ a888cc8f-e27c-4989-a535-6a2862c60c91
begin
	f = h5open(joinpath(modeldir, starsfile))
	pidx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(pidx)]
	pidx = sort(pidx)
	close(f)
end

# ╔═╡ addf19a4-088b-4aff-89a9-df73e8049f2c
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bins = LinRange(-1, 1, 100)
	probs = probabilities[snap_i.index]
	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=probs, bins = bins)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")

	probs = probabilities[snap_f.index]

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=probs, bins = bins)
	hideydecorations!(ax2)
	linkaxes!(ax, ax2)
	fig
end

# ╔═╡ de8ebcd0-d6e1-4b16-aaec-5bcd47cad1bd
function plot_ρ_s!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions, pidx)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 40, weights=probabilities)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ e76583dc-eea9-43c8-8051-a58a5c68a942
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel =  L"\log\, \rho_\star\; [10^{10} M_\odot / \textrm{kpc}^3]", 
		limits=((-1.9, 0.5), (-7, 2)))

	#vlines!(log10(r_s_s), label="r_s")

	plot_ρ_s!(snap_i, label="initial")
	plot_ρ_s!(snap_f, label="final")

	log_r_pred = LinRange(-2, 2, 1000)
	ρ_s_pred = ρ_s.(10 .^ log_r_pred)

	lines!(log_r_pred, log10.(ρ_s_pred), label="expected", color="black", linestyle=:dot)
	
	axislegend(ax)
	fig
end

# ╔═╡ 341440a0-9567-4ebf-8acb-cf327edfa4fb
function find_radii_fracs(out, x_cen, probabilities) 
	rs = Vector[]
	Ms = Vector[]
	rs_s = Vector[]

	percens = [0.003, 0.01, .03, .1, .3, .9]
	
	for i in 1:length(out)
		r = lguys.calc_r(out[i].positions .- x_cen[:, i])
		s_idx = sortperm(r)
		push!(rs, r[s_idx])
		
		ps = probabilities[out[i].index[s_idx]]
		ms = cumsum(ps)
		ms ./= ms[end]
		idxs = searchsortedfirst.([ms], percens)
		if i % 50 == 0
			println(idxs)
		end
		push!(Ms, ms)
		push!(rs_s, r[s_idx][idxs])
	end

	rs = hcat(rs...)
	rs_s = hcat(rs_s...)

	return percens, rs, rs_s

end

# ╔═╡ ab0b1a03-325c-489a-ac2f-309560541085
percens, rs, rs_s = find_radii_fracs(out, x_cen, probabilities)

# ╔═╡ d8f546d3-9e2e-4703-b652-5bea7bbbbd26
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel="log r containing stellar mass")
	for i in eachindex(percens)
		lines!(times, log10.(rs_s[i, :]), 
			color=i, colorrange=(1, length(percens)),
			label="$(percens[i])")
	end
	Legend(fig[1,2], ax, "fraction")
	fig
end

# ╔═╡ 1c83bf4c-12b2-4413-b4a9-87b75da1dcd8
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="log r containing stellar mass")
	for i in eachindex(percens)
		lines!(1:length(times), log10.(rs_s[i, :]), 
			color=i, colorrange=(1, length(percens)),
			label="$(percens[i])")
	end
	Legend(fig[1,2], ax, "fraction")
	fig
end

# ╔═╡ 193273c9-5b13-4af6-a345-4326cdebcf04
let
	fig = Figure()
	ax = Axis(fig[1,1],
	yscale=log10, xlabel="probability", ylabel="count")
	
	hist!(probabilities, bins=20, )
	fig
end

# ╔═╡ Cell order:
# ╠═6238c3fa-1974-11ef-1539-2f50fb31fe9a
# ╠═141d3b97-e344-4760-aa12-a48b1afb125c
# ╠═5b1dd353-a437-47cd-94be-7da9684581da
# ╠═28ba4f0f-6cc6-44e2-a7bc-4eee460d91b0
# ╠═21adbbe7-c8cc-4094-9e75-b68d97fa211a
# ╠═312576e2-16da-4285-9c19-a7a8005acf25
# ╠═feb6cc17-a25d-4ba8-a152-78412f422b80
# ╠═4445c2a1-3f00-42fa-8259-5d2b33e61d75
# ╠═051b314e-6cbc-405e-8873-b0288a93378d
# ╠═c76c8acc-ea88-4ce1-81cb-be0b67ef23fd
# ╠═7dfc4066-eee5-4118-9d12-48861aa66e03
# ╠═0671aa68-3d02-4b6c-a238-adffede13bd8
# ╠═fc7e58ec-ac03-40a4-99ae-2a16166402d3
# ╠═e32e2d66-2dad-4f09-84f1-a3081e3891a7
# ╠═3150cdfd-7573-4db9-86b7-ef614150a7b9
# ╠═2cc047db-ae05-42c5-898f-702ae3b83bd6
# ╠═052d229a-0362-42eb-b03c-fdab0f0bc6b4
# ╠═a888cc8f-e27c-4989-a535-6a2862c60c91
# ╠═addf19a4-088b-4aff-89a9-df73e8049f2c
# ╠═de8ebcd0-d6e1-4b16-aaec-5bcd47cad1bd
# ╠═e76583dc-eea9-43c8-8051-a58a5c68a942
# ╠═341440a0-9567-4ebf-8acb-cf327edfa4fb
# ╠═ab0b1a03-325c-489a-ac2f-309560541085
# ╠═d8f546d3-9e2e-4703-b652-5bea7bbbbd26
# ╠═1c83bf4c-12b2-4413-b4a9-87b75da1dcd8
# ╠═193273c9-5b13-4af6-a345-4326cdebcf04
