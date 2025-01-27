### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ a71257da-a168-11ef-02e4-e19693f816fa
begin
	using Pkg; Pkg.activate()

	using LilGuys
	using CSV, DataFrames
	using Arya
	using CairoMakie, GeoMakie
end

# ╔═╡ ebde81c2-3d5d-4736-8ffc-78b8fa925b51
df_pms = CSV.read("mv22_pm.csv", DataFrame, delim = "\t", ignorerepeated=true)

# ╔═╡ 35ea37d8-e012-428f-af7d-d6c592b4d92c
df_b22 = CSV.read("battaglia+22_galaxies.tsv", DataFrame, delim="\t", ignorerepeated=true)

# ╔═╡ eaaa6d20-84a0-4482-be2c-3559456bb651
begin
	df_all = innerjoin(df_pms, df_b22, on=:Galaxy)
	df_all[!, :distance] = 0.010 * 10 .^ (df_all.dm/5)

	rename!(df_all, :vlos => :radial_velocity)
	df_all
end

# ╔═╡ 34c184e2-3a11-4f29-a63c-b242f2a27d31
scatter(df_all.ra .- 360 .* (df_all.ra .> 180), df_all.dec)

# ╔═╡ a9d3769e-8e93-4941-bae9-ee19e46654f6
let
	fig = Figure()
	ax = GeoAxis(fig[1,1], dest="+proj=aitoff")
	
	p = arrows!(df_all.ra, df_all.dec, df_all.pmra, df_all.pmdec, color=df_all.radial_velocity, colormap=Reverse(:redblue))

	Colorbar(fig[1,2], p)
	fig
end

# ╔═╡ 7f37a48b-1169-4d53-8521-72c4ba0594d8
CSV.write("all_galaxies.csv", df_all)

# ╔═╡ Cell order:
# ╠═a71257da-a168-11ef-02e4-e19693f816fa
# ╠═ebde81c2-3d5d-4736-8ffc-78b8fa925b51
# ╠═35ea37d8-e012-428f-af7d-d6c592b4d92c
# ╠═eaaa6d20-84a0-4482-be2c-3559456bb651
# ╠═34c184e2-3a11-4f29-a63c-b242f2a27d31
# ╠═a9d3769e-8e93-4941-bae9-ee19e46654f6
# ╠═7f37a48b-1169-4d53-8521-72c4ba0594d8
