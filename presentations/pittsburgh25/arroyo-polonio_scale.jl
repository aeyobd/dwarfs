### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ cdaf2bae-afaf-4858-b07a-eb00107e64c3
import Pkg; Pkg.activate()

# ╔═╡ f115979c-3f07-11f0-2ff6-7ddc61640d3e
using Arya, CairoMakie

# ╔═╡ 6aae8cc9-1682-4233-8196-14f2dbba3fa8
using LilGuys

# ╔═╡ bdca1255-292a-434e-9be9-914cdc0fcc95
R_h = 9.46 / sqrt(1 - 0.36) / 60

# ╔═╡ 428559de-3172-438d-91ed-18a62f0212c5
range_deg = (0.01786759961324086,  1.27)

# ╔═╡ a0b04e51-fbb7-48a5-9ccd-06cca6241752
range_R_h = range_deg ./ R_h

# ╔═╡ a0c3397d-5b38-4489-9510-c57e5764429e
let
	fig = Figure(figure_padding=20)
	ax1 = Axis(fig[1,1], limits=(range_deg, (0.0,1e-5)), 
			   xticksmirrored=false, aspect=1e4,xaxisposition=:top, 
			   xscale=log10, xticks=[0.1, 1.0], xtickalign=0.0, xminortickalign=0.0)

	
	ax2 = Axis(fig[1,1], limits=(log10.(range_R_h), (0,1e-5)), xticksmirrored=false, aspect=1e4, xlabel=L"R / R_h", xtickalign=0.0, xminortickalign=0.0)

	
	hidespines!(ax2)
	hideydecorations!(ax1)
	hideydecorations!(ax2)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 164f32b6-e6e6-4112-8cf5-0ccca1d49861
let
	fig = Figure(figure_padding=20)
	ax1 = Axis(fig[1,1], limits=(log10.(range_deg .* 60), (0.0,1e-5)), 
			   xticksmirrored=false, aspect=1e4,xaxisposition=:top, 
			   xlabel="Radius / arcmin", xtickalign=0.0, xminortickalign=0.0)

	
	ax2 = Axis(fig[1,1], limits=(log10.(range_R_h), (0,1e-5)), xticksmirrored=false, aspect=1e4, xlabel=L"R / R_h", xtickalign=0.0, xminortickalign=0.0)

	
	hidespines!(ax2)
	hideydecorations!(ax1)
	hideydecorations!(ax2)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ f75f27a1-d791-45de-9504-1f0ae306fe60
prof_inner = LilGuys.Plummer(r_s=0.129, M=0.34)

# ╔═╡ 703eaed6-39f1-4393-9509-403d83430e9f
prof_outer = LilGuys.Plummer(r_s=0.27, M=0.66)

# ╔═╡ fef20079-1315-40c3-a94f-2dd76623101f
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(-2, 0.3, 1000)
	r = 10 .^ x

	sigma_in = LilGuys.surface_density.(prof_inner, r)
	sigma_out = LilGuys.surface_density.(prof_outer, r)


	lines!(x, log10.(sigma_in))
	lines!(x, log10.(sigma_out))
	lines!(x, log10.(sigma_out .+ sigma_in))

	fig
end

# ╔═╡ 2ddb452d-7782-4113-8c2a-9d3df09a19a3


# ╔═╡ e0f11f1e-d2c6-42ab-ad02-d2419e5f2fe4
let
	fig = Figure(figure_padding=20)
	ax1 = Axis(fig[1,1], limits=(log10.(range_deg .* 60), (0.0,1e-5)), 
			   xticksmirrored=false, aspect=1e4,xaxisposition=:top, 
			   xlabel="Radius / arcmin", xtickalign=0.0, xminortickalign=0.0)

	
	ax2 = Axis(fig[1,1], limits=(log10.(range_R_h), (0,1e-5)), xticksmirrored=false, aspect=1e4, xlabel=L"R / R_h", xtickalign=0.0, xminortickalign=0.0)

	
	hidespines!(ax2)
	hideydecorations!(ax1)
	hideydecorations!(ax2)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ af12f8e6-49cd-45b1-9610-9838bd632bc2
[0.01:0.01:0.1; 0.1:0.1:1]

# ╔═╡ Cell order:
# ╠═cdaf2bae-afaf-4858-b07a-eb00107e64c3
# ╠═f115979c-3f07-11f0-2ff6-7ddc61640d3e
# ╠═6aae8cc9-1682-4233-8196-14f2dbba3fa8
# ╠═bdca1255-292a-434e-9be9-914cdc0fcc95
# ╠═428559de-3172-438d-91ed-18a62f0212c5
# ╠═a0b04e51-fbb7-48a5-9ccd-06cca6241752
# ╠═a0c3397d-5b38-4489-9510-c57e5764429e
# ╠═164f32b6-e6e6-4112-8cf5-0ccca1d49861
# ╠═f75f27a1-d791-45de-9504-1f0ae306fe60
# ╠═703eaed6-39f1-4393-9509-403d83430e9f
# ╠═fef20079-1315-40c3-a94f-2dd76623101f
# ╠═2ddb452d-7782-4113-8c2a-9d3df09a19a3
# ╠═e0f11f1e-d2c6-42ab-ad02-d2419e5f2fe4
# ╠═af12f8e6-49cd-45b1-9610-9838bd632bc2
