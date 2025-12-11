### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 8210d650-c4a3-11f0-ba24-efea8c8c5414
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
	using DataFrames, CSV
end

# ╔═╡ 478f92ef-0e95-4d0f-a9ac-8104280902d6
CairoMakie.activate!(type=:png)

# ╔═╡ e320a3f6-49e5-47ad-8698-117c7b3a505a
saga = CSV.read("saga.txt", DataFrame)

# ╔═╡ c72b1f54-b60f-4d51-8e65-fa51d33cf296
names(saga)

# ╔═╡ 681962ca-bf54-4e53-8c77-dc3c64224a13
saga.Object

# ╔═╡ fe297631-d40c-434c-94f0-bc554c4c8de7
scatter(saga[!, "log g"], saga.Teff)

# ╔═╡ da41bcb1-e867-426f-8b31-c092af016417
function get_abund_ratio(df, col)
	ele1, ele2 = split(col, "/")
	if ele2 == "H"
		xs = df[!, "[$ele1/H]"]
	else
		xs = df[!, "[$ele1/H]"] .- df[!, "[$ele2/H]"]
	end

	return xs
end

# ╔═╡ baf0d831-80ed-44e1-8a56-09cb016e557d
function plot_abund!(df, x, y; color=nothing, kwargs...)
	xs = get_abund_ratio(df, x)
	ys = get_abund_ratio(df, y)


	filt = .!ismissing.(xs)
	filt .&= .!ismissing.(ys)
	if (color isa String) && (length(split(color, "/")) == 2)
		color = get_abund_ratio(df, color)
		filt .&= .!ismissing.(color)
		color = color[filt]
	end
	
	scatter!(xs[filt], ys[filt]; color=color, kwargs...)

end

# ╔═╡ 5eabead9-3978-4dcb-82ad-e58b75355d63
function plot_abund(df, x, y; kwargs...)

	fig = Figure()
	ax = Axis(fig[1,1], xlabel="[$x]", ylabel="[$y]")
	p = plot_abund!(df, x, y; kwargs...)

	if length(split(get(kwargs, :color, ""), "/")) == 2
		Colorbar(fig[1,2], p, label="[$(kwargs[:color])]")
	end
	fig
end

# ╔═╡ 41a00c0b-8165-4a56-b73c-39dafa788416
plot_abund(saga, "Mg/H", "Mg/Fe")

# ╔═╡ 713b2fed-6615-4b21-a514-de291760bb5d
plot_abund(saga, "Mg/H", "Ba/Sc")

# ╔═╡ 91c29d8e-2903-4efa-9bf6-22caae0fe9c4
plot_abund(saga, "Ba/Mg", "Ba/Sr")

# ╔═╡ b100e134-0393-40cb-8660-c2af75e4b804
plot_abund(saga, "Ba/Fe", "Sr/Ba", color="Fe/H", markersize=3)

# ╔═╡ 3fa5d21d-2ede-48ab-81e8-0d103b033fc5
plot_abund(saga, "Fe/H", "Sr/Fe", color="Mg/Fe", markersize=3)

# ╔═╡ 18ab55da-32c9-4839-91db-8ab9b447601c
plot_abund(saga, "Fe/H", "Sr/Fe", markersize=2,)

# ╔═╡ 0f250408-78b9-486a-9397-4b5134a0e1f4
plot_abund(saga, "Si/H", "Sr/Si", markersize=2,)

# ╔═╡ 1142ffda-98ee-45a6-b170-5344302eddf6
plot_abund(saga, "Si/H", "Sr/Ba", color="Mg/Si")

# ╔═╡ add4ed9a-91f4-42fc-9f34-3dbab7e7bd77
plot_abund(saga, "Ba/H", "Sr/Ba", color="Fe/H")

# ╔═╡ e5d716d0-194c-422b-aa4a-7c52b3a1dbde
plot_abund(saga, "Mg/H", "Al/Mg")

# ╔═╡ 11cb0780-5060-41fc-9476-3dce2debc4cb
plot_abund(saga, "Mg/H", "Ba/Sr")

# ╔═╡ 0c3af7ec-f40e-4280-9dbc-a2372804b4d4
plot_abund(saga, "Eu/Mg", "Ba/Eu", color="Mg/H")

# ╔═╡ 5f9e9594-d71f-4fce-9f9c-f6abc5e8818b
plot_abund(saga, "Mg/H", "Ba/Eu", color="Mg/Fe")

# ╔═╡ 80d12e0a-319e-43f6-be85-58dea00c246d
plot_abund(saga, "Mg/H", "Ba/Sr", color="Mg/Ni")

# ╔═╡ 2c3f0092-0b22-4d33-a081-57ef58462df5
plot_abund(saga, "Mg/H", "Al/Mg", color="Mg/Fe")

# ╔═╡ 2a06fa2f-5d3a-40ba-8b10-ee206dbea9cc
plot_abund(saga, "Na/Mg", "K/Mg", color="Mg/H")

# ╔═╡ 7c3a76e0-0e72-4ea2-b4b5-18029a24d9ce
plot_abund(saga, "Mg/H", "C/H", color="Mg/Fe")

# ╔═╡ abd96e4e-33ee-46fc-a159-1f351d9bab8f
plot_abund(saga, "O/H", "C/O", color="O/Fe")

# ╔═╡ e4b89911-dedb-4a0f-addb-037e402c1664
plot_abund(saga, "Mg/H", "C/Mg", color="Mg/Fe")

# ╔═╡ 9974f5a2-896b-44fc-96d0-6e9d7904dd9d
plot_abund(saga, "Mg/H", "C/Mg")

# ╔═╡ ad84d695-e765-480a-9f90-904e4c657d8b
plot_abund(saga, "Mg/H", "C/Mg", color="Ba/Mg")

# ╔═╡ 2de27e2d-1377-49db-be0b-a8ccd81256e4
plot_abund(saga, "Mg/H", "C/Mg", color="Eu/Mg")

# ╔═╡ Cell order:
# ╠═8210d650-c4a3-11f0-ba24-efea8c8c5414
# ╠═478f92ef-0e95-4d0f-a9ac-8104280902d6
# ╠═e320a3f6-49e5-47ad-8698-117c7b3a505a
# ╠═c72b1f54-b60f-4d51-8e65-fa51d33cf296
# ╠═681962ca-bf54-4e53-8c77-dc3c64224a13
# ╠═fe297631-d40c-434c-94f0-bc554c4c8de7
# ╠═da41bcb1-e867-426f-8b31-c092af016417
# ╠═baf0d831-80ed-44e1-8a56-09cb016e557d
# ╠═5eabead9-3978-4dcb-82ad-e58b75355d63
# ╠═41a00c0b-8165-4a56-b73c-39dafa788416
# ╠═713b2fed-6615-4b21-a514-de291760bb5d
# ╠═91c29d8e-2903-4efa-9bf6-22caae0fe9c4
# ╠═b100e134-0393-40cb-8660-c2af75e4b804
# ╠═3fa5d21d-2ede-48ab-81e8-0d103b033fc5
# ╠═18ab55da-32c9-4839-91db-8ab9b447601c
# ╠═0f250408-78b9-486a-9397-4b5134a0e1f4
# ╠═1142ffda-98ee-45a6-b170-5344302eddf6
# ╠═add4ed9a-91f4-42fc-9f34-3dbab7e7bd77
# ╠═e5d716d0-194c-422b-aa4a-7c52b3a1dbde
# ╠═11cb0780-5060-41fc-9476-3dce2debc4cb
# ╠═0c3af7ec-f40e-4280-9dbc-a2372804b4d4
# ╠═5f9e9594-d71f-4fce-9f9c-f6abc5e8818b
# ╠═80d12e0a-319e-43f6-be85-58dea00c246d
# ╠═2c3f0092-0b22-4d33-a081-57ef58462df5
# ╠═2a06fa2f-5d3a-40ba-8b10-ee206dbea9cc
# ╠═7c3a76e0-0e72-4ea2-b4b5-18029a24d9ce
# ╠═abd96e4e-33ee-46fc-a159-1f351d9bab8f
# ╠═e4b89911-dedb-4a0f-addb-037e402c1664
# ╠═9974f5a2-896b-44fc-96d0-6e9d7904dd9d
# ╠═ad84d695-e765-480a-9f90-904e4c657d8b
# ╠═2de27e2d-1377-49db-be0b-a8ccd81256e4
